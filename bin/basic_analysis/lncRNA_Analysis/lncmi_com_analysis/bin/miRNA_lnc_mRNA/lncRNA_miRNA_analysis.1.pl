#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Data::Dumper;
use File::Basename qw(basename dirname);
my ($mi_rna,$lnc_rna,$od,$cfg,$step,$oneStepOnly,$lnc_tar);
GetOptions(
    "mi_rna:s" =>\$mi_rna,
    "lnc_rna:s"=>\$lnc_rna,
    "lnc_tar:s"=>\$lnc_tar,
    "od:s"=>\$od,
    "cfg:s"=>\$cfg,
    "step:s"=>\$step,
    "oneStepOnly:s"=>\$oneStepOnly,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($mi_rna and $lnc_rna and $od and $cfg and $lnc_tar);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
&MKDIR("$od/work_sh");
$mi_rna=&ABSOLUTE_DIR($mi_rna);
$lnc_rna=&ABSOLUTE_DIR($lnc_rna);
$lnc_tar=&ABSOLUTE_DIR($lnc_tar);
$cfg=&ABSOLUTE_DIR($cfg);
#$m_tar=&ABSOLUTE_DIR($m_tar);
my $notename = `hostname`; chomp $notename;
$step = $step || 1;
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my %para;
open (IN,"$cfg") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	
	my @tmp=split /\s+/,$_;
	$para{$tmp[0]}=$tmp[1];
}
close IN;

############################################################################################
if ($step==1){
    print "lncRNA_Precursor_analysis start";
    &MKDIR("$od/precursor");
    #&MKDIR("$od/precursor/work_sh");
    
    open OUT,">$od/work_sh/precursor_analyses.sh";
    if (exists $para{'Key'}){
	print "get the key for $para{'SPECIES_NAME'} :$para{'Key'}";
	print OUT "perl $Bin/precursor_analysis/Precursor_analysis.pl -fa $lnc_rna -od $od/precursor -key $para{'Key'} ";
	&Cut_shell_qsub("$od/work_sh/precursor_analyses.sh",18,"15G","general.q");
	&Check_qsub_error("$od/work_sh/precursor_analyses.sh");
	$step++ unless ($oneStepOnly) ;
    }else{
	print "the key for $para{'SPECIES_NAME'} missing !!\n\nplease check the $cfg and make sure to skip the step!!";
	$step++ unless ($oneStepOnly) ;
    }
    
   close OUT; 
    
}
##############################################################################################
if ($step==2){
	print "lncRNA miRNA target analysis start ..........";
	&MKDIR("$od/miRNA_Target2LncRNA");
	open OUT,">$od/work_sh/lncRNA_target2miRNA.sh";
	print OUT "perl $Bin/target_prediction/mir2target.pl -miR $mi_rna -gene $lnc_rna -type $para{'SPECIES_TYPE'} -cfg $cfg -key $para{'SPECIES_NAME'} -od $od/miRNA_Target2LncRNA";
	print OUT "  &&  perl $Bin/target_prediction/result_tranver.pl -mi_tar $od/miRNA_Target2LncRNA/$para{'SPECIES_NAME'}.mir2target.list -od $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list ";
	&Cut_shell_qsub("$od/work_sh/lncRNA_target2miRNA.sh",18,"15G","general.q");
	&Check_qsub_error("$od/work_sh/lncRNA_target2miRNA.sh");
	$step++ unless ($oneStepOnly) ;
	
}
###############################################################################################
my ($lnc_diff,$miRNA_diff);
if ($step==3){
	print "get lncRNA_miRNA_mRNA interactiong start";
	&MKDIR("$od/mRNA_lncRNA_miRNA");
	open OUT,">$od/work_sh/mRNA_lncRNA_miRNA.sh";
	
	if (exists $para{'Lnc_Diff'} && exists $para{'miRNA_Diff'} ){
		$lnc_diff=&ABSOLUTE_DIR($para{'Lnc_Diff'});
		$miRNA_diff=&ABSOLUTE_DIR($para{'miRNA_Diff'});
		if (-d $lnc_diff && -d $miRNA_diff){
			print OUT "perl $Bin/miRNA_lnc_mRNA/get_interaction_result.pl -mi_lnc $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lnc_m $lnc_tar -od $od/mRNA_lncRNA_miRNA -in1 $lnc_diff -in2 $miRNA_diff ";
		}
	}
	if (exists $para{'Lnc_Diff'} && !exists $para{'miRNA_Diff'} ){
		$lnc_diff=&ABSOLUTE_DIR($para{'Lnc_Diff'});
		if (-d $lnc_diff ){
			print OUT "perl $Bin/miRNA_lnc_mRNA/nondeg_result.pl -mi_lnc $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lnc_m $lnc_tar -od $od/mRNA_lncRNA_miRNA -in1 $lnc_diff  ";
		}
	}
	if (!exists $para{'Lnc_Diff'} && !exists $para{'miRNA_Diff'} ){
		print OUT "perl $Bin/miRNA_lnc_mRNA/result.pl -mi_lnc $od/miRNA_Target2LncRNA/lncRNA_target2mirna.list -lnc_m $lnc_tar -od $od/mRNA_lncRNA_miRNA  ";
			
	}
	if (exists $para{'mi2mRNA_tar'}){
		my $mi_m_tar=&ABSOLUTE_DIR($para{'mi2mRNA_tar'});
		print OUT "-mi_m $mi_m_tar ";
	}
	&Cut_shell_qsub("$od/work_sh/mRNA_lncRNA_miRNA.sh",18,"15G","general.q");
	&Check_qsub_error("$od/work_sh/mRNA_lncRNA_miRNA.sh");
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
##########################################################################################################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		if ($notename=~/login\-0\-4/) {
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
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
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
}
###############################################################################################
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
################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
###############################################################################################
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: 1.0
     Usage:
            -mi_rna      <FILE>  ./*miRNA.fa
            -lnc_rna	<FILE>   ./lncRNA.fa
	    -lnc_tar	<FILE>	./novel_lncRNA_target.xls
	    -cfg	*cfg
            -od          dir
	    -step	defaut 1
			1:lncRNA_Precursor_analysis
			2:lncRNA miRNA target analysis
			3:lncRNA miRNA mRNA regulation network
	    -oneStepOnly 
	    -h 		help documents

   Example:
            perl $Script -mi_rna ./miRNA.fa -lnc_rna lncRNA.fa -lnc_tar ./novel_lncRNA_target.xls  -cfg free_pipeline.cfg -od  ./
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
