#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my $version="v1.0.2";
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
my $notename = `hostname`; chomp $notename;
&log_current_time("$Script start...");
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($mi_rna,$lnc_rna,$lnc_tar,$lnc_fpkm, $cfg, $od, $step, $step_by_step,$gtf,$mi_m_tar);

GetOptions(
    "lnc_fpkm:s" =>\$lnc_fpkm,
    "cfg:s" =>\$cfg,
    "od:s"   =>\$od,
    "mi_rna:s" =>\$mi_rna,
    "lnc_rna:s"=>\$lnc_rna,
    "lnc_tar:s"=>\$lnc_tar,
    "mi_m_tar:s"=>\$mi_m_tar,
    "step:s" =>\$step,
    "sbs"  =>\$step_by_step,
    "gtf:s"   =>\$gtf,
    #"CSE:s"  =>\$customer_service_exe,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($lnc_rna and $lnc_tar and $cfg and $od );
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------
my %para;
&LOAD_PARA($cfg,\%para);
$para{'Known_db'}= "noncode" unless (defined $para{'Known_db'});
if (!exists $para{'Project_key'}){
	warn "the Project_key must be given !!";
}

########get step analysis 
my %step;
$step ||= ($step_by_step) ? '1' : join ',',(1..4);
&steps_process($step,$step_by_step,\%step);

if (defined $para{'sRNA_target'} && $para{'sRNA_target'} eq "off"){
	delete $step{3};
	delete $step{4};
}

if ($para{'Key'}=~m/Chromalveolata|Metazoa|Mycetozoa|Viridiplantae|Viruses/){
    print "analysis the shu !!";
}
########get mirna database
my $MIRNA=$config{miRBase};
########
mkdir  $od unless (-d $od);
mkdir "$od/work_sh" unless (-d "$od/work_sh");
if ($step{1}) {
	print "Known_lncRNA_Analysis start\n\n";
    print "Known_lncRNA_Analysis work_sh dir:$od/work_sh/Known_lncRNA_Analysis.sh\n\n";
    &MKDIR("$od/Known_lncRNA/");
	open SH1,">$od/work_sh/Known_lncRNA_Analysis.sh" or die $!;
	print SH1 "perl $Bin/bin/known_lncRNA/known_lncRNA.pl -fa $lnc_rna -key $para{'Known_db'} -out $od/Known_lncRNA/ ";
	&Cut_shell_qsub("$od/work_sh/Known_lncRNA_Analysis.sh",18,$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/Known_lncRNA_Analysis.sh");
	
}
if ($step{2}){
    print "lncRNA_Precursor_analysis start \n\n";
    print "lncRNA_Precursor_analysis work_sh dir:$od/work_sh/precursor_analyses.sh\n\n";
    #######################
    &MKDIR("$od/precursor");
    open OUT,">$od/work_sh/precursor_analyses.sh";
    if (exists $para{'Key'}){
    print "get the key for $para{'Project_key'} :$para{'Key'}";
    print OUT "perl $Bin/bin/precursor_analysis/Precursor_analysis.pl -fa $lnc_rna -od $od/precursor -key $para{'Key'} ";
    &Cut_shell_qsub("$od/work_sh/precursor_analyses.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
    &Check_qsub_error("$od/work_sh/precursor_analyses.sh");
    
    }else{
    print "the key for $para{'Project_key'} missing !!\n\nplease check the $cfg and make sure to skip the step!!";
    }
    
   close OUT; 
}
if ($step{3}){
    print "lncRNA miRNA target analysis start ..........\n\n";
    print "lncRNA miRNA target analysis work_sh dir:$od/work_sh/lncRNA_target2miRNA.sh\n\n";
    &MKDIR("$od/miRNA_Target2LncRNA");
    open OUT,">$od/work_sh/lncRNA_target2miRNA.sh";
    if (exists $para{'miRNA_fa'} and -e $para{'miRNA_fa'}) {
            print OUT "perl $Bin/bin/target_prediction/mir2target.cut.pl -miR $para{'miRNA_fa'} -lnc $lnc_rna -type $para{'SPECIES_TYPE'} -cfg $cfg -key $para{'Project_key'} -od $od/miRNA_Target2LncRNA  &&  ";
    }
    if (!exists $para{'miRNA_fa'}) {
        print OUT "perl $Bin/bin/target_prediction/select_fa.pl -fa $MIRNA -i $para{'Key'} -o $od/miRNA_Target2LncRNA/$para{'Key'}_sRNA.fa &&  ";
        print OUT "perl $Bin/bin/target_prediction/mir2target.cut.pl -miR $od/miRNA_Target2LncRNA/$para{'Key'}_sRNA.fa -lnc $lnc_rna -type $para{'SPECIES_TYPE'} -cfg $cfg -key $para{'Project_key'} -od $od/miRNA_Target2LncRNA  &&  ";
    }
    
    print OUT "perl $Bin/bin/target_prediction/result_tranver.pl -mi_tar $od/miRNA_Target2LncRNA/$para{'Project_key'}.mir2target.list -out $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list ";
    &Cut_shell_qsub("$od/work_sh/lncRNA_target2miRNA.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
    &Check_qsub_error("$od/work_sh/lncRNA_target2miRNA.sh");
}
if ($step{4}){
    print "get lncRNA_miRNA_mRNA interactiong start\n\n";
    print "lncRNA_miRNA_mRNA interactiong analysis work_sh dir:$od/work_sh/mRNA_lncRNA_miRNA.sh\n\n";
    &MKDIR("$od/mRNA_lncRNA_miRNA");
    open OUT,">$od/work_sh/mRNA_lncRNA_miRNA.sh";
    if (exists $para{'miRNA2mRNA'}){
        if (-e $para{'miRNA2mRNA'}){
                print OUT "perl $Bin/bin/miRNA_lnc_mRNA/result.pl -miRNA2mRNA $mi_m_tar -miRNA2lncRNA $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lncRNA2mRNAtar $lnc_tar -od $od/mRNA_lncRNA_miRNA  ";
        }
    }else{
            print OUT "perl $Bin/bin/miRNA_lnc_mRNA/result.pl -miRNA2lncRNA $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lncRNA2mRNAtar $lnc_tar -od $od/mRNA_lncRNA_miRNA  ";

    }
    #print OUT "perl $Bin/bin/miRNA_lnc_mRNA/result.pl -miRNA2mRNA $mi_m_tar -miRNA2lncRNA $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lncRNA2mRNAtar $lnc_tar -od $od/mRNA_lncRNA_miRNA  ";
    &Cut_shell_qsub("$od/work_sh/mRNA_lncRNA_miRNA.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
    &Check_qsub_error("$od/work_sh/mRNA_lncRNA_miRNA.sh");
}
##############################################################################################
###############################################################################################



sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
        '1' => 'Known_lncRNA_Analysis',
        '2' => 'Precursor_Analysis',
        '3' => 'LncRNA2miRNA_target',
        '4' => 'LncRNA_miRNA_mRNA_analysis',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ($step_str =~/^\d+$/) {
            for my $s ($step_str..4) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            die "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
        }
    } else {
        print "step-by-step: OFF\n";
        for my $s (split /,/,$step_str) {
            if ($s =~/^\d+$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                die "ERROR: illegal steps specified by --step.\n";
            }
        }
    }

    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &log_current_time("step check done.\n");
}

#############################################################################################################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
        # &Cut_shell_qsub($shell,$cpu,$vf,$queue);
        my $shell = shift;
        my $cpu = shift;
        my $vf = shift;
        my $queue = shift;

        my $line = `less -S $shell |wc -l `;
        if ($line<=1000) {
                if ($notename=~/login\-0\-4/) {
                        system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
                }
                else
                {
                        system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
                }
        }
        if ($line>1000) {
                my @div=glob("$shell.div*");
                foreach (@div) {
                        if (-e $_) {
                                system "rm $_";
			}
                }
                @div=();
                my $div_index=1;
                my $line_num=1;
                open IN,"$shell" || die;
                while (<IN>) {
                        chomp;
                        open OUT,">>$shell.div.$div_index.sh" || die;
                        if ($line_num<1000) {
                                print OUT "$_\n";
                                $line_num++;
                        }
                        else {
                                print OUT "$_\n";
                                $div_index++;
                                $line_num=1;
                                close OUT;
                        }
                }
                if ($line_num!=1) {
                        close OUT;
                }
                @div=glob "$shell.div*";
                foreach my $div_file (@div) {
                        if ($notename=~/login\-0\-4/) {
                                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
                        }
                        else
                        {
                                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
                        }
                }
        }
}



sub Check_qsub_error {#
        # Check The qsub process if error happend 
        my $sh=shift;
        my @Check_file=glob "$sh*.qsub/*.Check";
        my @sh_file=glob "$sh*.qsub/*.sh";

        if ($#sh_file!=$#Check_file) {
                print "Their Some Error Happend in $sh qsub, Please Check..\n";
                die;
        }
        else {
                print "$sh qsub is Done!\n";
        }
}
######################load para
sub LOAD_PARA
{
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}


#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


########################################################################################################

sub MKDIR{
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}


########################################################################################################
sub ABSOLUTE_DIR{ 
        my $cur_dir=`pwd`;chomp($cur_dir);
        my ($in)=@_;
        my $return="";
        if(-f $in){
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }elsif(-d $in){
                chdir $in;$return=`pwd`;chomp $return;
        }else{
                warn "Warning just for file and dir\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}
#############################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: niulg <renhd\@biomarker.com.cn> 
      Date: 

     Usage:
            --cfg     <FILE>  detail config, analysis parameters ,must be give;
            --od        <DIR>   ./lncRNA_com_miRNA  ,must be give;
	    --lnc_rna	<FILE>   ./lncRNA.fa,must be give;
	    --lnc_tar       <FILE>  ./novel_lncRNA_target.xls,(for mRNA_lncRNA_miRNA regulation ),choise step4 needed; 
	   ###################combine analysis project 
		--mi_rna	<FILE>  ./*miRNA.fa (for miRNA_lncRNA target gene predict),choise step3 needed;
		--mi_m_tar	<FILE>	Target_Predict/Human.mir2target.list,(for mRNA_lncRNA_miRNA regulation ),choise step4 needed;
					
            --step      <INT>   step to start from or steps to run (split by comma) [1]
                        1	Known_lncRNA_Analysis,
                        2	Precursor_Analysis,
                        3	LncRNA2miRNA_target,
                        4	LncRNA_miRNA_mRNA_analysis
            --sbs               step-by-step analysis, start from the step specified by -step
                                or not, only run the steps specified by -step

            
            --h                 help documents

   Example:
            perl $Script --cfg personality.cfg --od ./Personality  --lnc_rna ./lncRNA.fa

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

