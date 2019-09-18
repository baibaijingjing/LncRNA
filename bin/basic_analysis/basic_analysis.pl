#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
my $version="2.4.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $odir, $step,$oneStepOnly,$test);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "od:s"   =>\$odir,
    "step:i" =>\$step,
    "test" =>\$test,
    "oneStepOnly:s"=>\$oneStepOnly,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir);
my $notename = `hostname`; chomp $notename;
system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$data_cfg=abs_path($data_cfg);
$detail_cfg=abs_path($detail_cfg);

&log_current_time("$Script start...");
# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
my (%data_cfg, %detail_cfg);

# read data config
&data_cfg_read($data_cfg,\%data_cfg);

# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);

$step ||= 1;
$detail_cfg{'Queue_type'}||=  "general.q";
$detail_cfg{'CPU'} ||= "30";
$detail_cfg{'Memory'} ||="30G";
# make primary dir
my $sh_dir = "$odir/work_sh";
my $tophat_cufflinks_dir = "$odir/Tophat_Cufflinks";
my $gene_expression_dir  = "$odir/geneExpression";
my $iso_expression_dir = "$odir/LncExpression";
mkdir $sh_dir unless (-d $sh_dir);
mkdir $tophat_cufflinks_dir unless (-d $tophat_cufflinks_dir);
mkdir $gene_expression_dir unless (-d $gene_expression_dir);
mkdir $iso_expression_dir unless (-d $iso_expression_dir);


#createLog($Title,$version,$$,"$odir",$test);
my $cmd;
# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
######################### Tophat & Cufflinks
if ($step==1) {
	open SH1,">$sh_dir/s2.1.tophat_process.sh";
	
     $cmd = "perl $Bin/tophat_cufflinks/tophat_prosess.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $tophat_cufflinks_dir";
       
        print SH1 $cmd;
        runOrDie("$sh_dir/s2.1.tophat_process.sh");
        $step++ unless ($oneStepOnly);
}
if ($step==2) {
	open SH1,">$sh_dir/s2.2.cufflinks_statistic.sh";
	if (exists $detail_cfg{'Assem'}){
            $cmd = "perl $Bin/tophat_cufflinks/tophat_cufflinks.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $tophat_cufflinks_dir \n";
        }else{
                $cmd = "perl $Bin/tophat_cufflinks/tophat_cut_filter.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $tophat_cufflinks_dir \n";
        }
        $cmd.="perl $Bin/tophat_cufflinks/Statistics_draw.pl  -cfg1 $data_cfg -cfg2 $detail_cfg -od $tophat_cufflinks_dir  ";
        print SH1 $cmd;
        qsubOrDie("$sh_dir/s2.2.cufflinks_statistic.sh",$detail_cfg{'Queue_type'},$detail_cfg{'CPU'},$detail_cfg{'Memory'});
        $step++ unless ($oneStepOnly);
}

##################################
my $genome = $detail_cfg{'Ref_seq'};
my $index = $detail_cfg{'Project_key'};
######################### Gene Expression Abstract

if ($step==3) {#需要分开为两步，转录本表达水平和lnc基因表达水平
    $cmd = "perl $Bin/gene_expression/gene_expression.pl -gtf $tophat_cufflinks_dir/Compare/Compare.gtf -genome $genome -cuffnorm $tophat_cufflinks_dir/Cuffnorm -cufflinks $tophat_cufflinks_dir/Cuffliks -index $index -od $gene_expression_dir ";#niulg,20150824
    $cmd.="-singleexon " if (exists $detail_cfg{'single_exon'});
    $cmd .= "  \n  perl $Bin/iso_expression/bin/isoExpression.pl -f $tophat_cufflinks_dir/Compare/Compare.gtf -c $tophat_cufflinks_dir/Cuffnorm -o  $iso_expression_dir\n";
    &step_cmd_process($cmd,"s2.3.gene_expression.sh",$sh_dir);
    $step++ unless ($oneStepOnly);
}
############################
######################### Subsequence
if ($step==4) {
    $cmd = "";
    for my $sample (sort keys %{$data_cfg{rawdata}}) {
        my $fq1 = $data_cfg{rawdata}{$sample}{fq1};
        my $gene_expression_list = "$gene_expression_dir/$sample.geneExpression.xls";
	$cmd .="$config{Rscript} $Bin/gene_expression/bin/fpkm_saturation_v3.R --i $gene_expression_list --o $gene_expression_dir/$sample.Saturation.png \n";#2015/08/11,modified by niulg
    }
    $cmd .= "cp $tophat_cufflinks_dir/Map_Stat/*png $gene_expression_dir\n";
    $cmd .= "cp $tophat_cufflinks_dir/Map_Stat/*mappedStat.xls $gene_expression_dir\n";
	$cmd .= "cp $tophat_cufflinks_dir/Map_Stat/*png $iso_expression_dir\n";
	$cmd .= "cp $tophat_cufflinks_dir/Map_Stat/*mappedStat.xls $iso_expression_dir\n";
    &step_cmd_process($cmd,"s2.4.Subsequence.sh",$sh_dir);
    #$step++ unless ($oneStepOnly);
}




#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);

        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
    }
    close CFG;
    &log_current_time("data config done.");
}

#############################################################################################################
sub detail_cfg_read {
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/ or /^$/);
        my ($key, $value) = (split /\s+/)[0,1];

        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_anno') {
            die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg'  or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
        }
	if ($key eq 'Assem') {
            $detail_cfg->{$key} = $value;
        }

    }
    close CFG;
    &log_current_time("detail config done.");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
#    &log_current_time("CMD: $cmd");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
}
########################
sub step_cmd_process_1 {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

#    $flag = system("sh $sh_file > $log_file");
#    if ($flag != 0){
#        log_current_time("Error: command failed: $cmd");
#        exit(1);
#    } else {
#        my $escaped_time = (time()-$start_time)."s";
#        &log_current_time("$sh_name done, escaped time: $escaped_time.");
#    }
}
#############################################################################################################
sub Cut_shell_qsub {
	my $shell = shift;
        my $cpu = shift;
        my $vf = shift;
        my $queue = shift;

        my $line = `less -S $shell |wc -l `;
        chomp $line;
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
                        open OUT,">>$shell.div.$div_index" || die;
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
#############################
sub Check_qsub_error {
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
#######################
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

#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            --cfg1      <FILE>  data config, rawdata & refseq path
            --cfg2      <FILE>  detail config, analysis parameters
            --od        <DIR>   analysis output directory

            --step      <INT>   step to start from                  [1]
                          1     Tophat & Cufflinks
                          2     Gene Expression & Iso Expression
			  			  3	    lncRNA_analysis
                          4     Subsequence
                          5     Lnc Subsequence
	    -oneStepOnly
            --h                 help documents

   Example:
            perl $Script --cfg1 data.cfg --cfg2 detail.cfg --od Analysis/

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
