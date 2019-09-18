#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($step,$mRNA_gtf,$lnc_gtf,$mRNA_fa,$miRNA_fa,$lncRNA_fa,$inter,$exp,$odir,$oneStepOnly,$cfg);
GetOptions(
	"help|?"=>\&USAGE,
	"step:s"=>\$step,
	"m:s"=>\$mRNA_gtf,
	"l:s"=>\$lnc_gtf,
	"mfa:s"=>\$mRNA_fa,
	"mifa:s"=>\$miRNA_fa,
	"lncfa:s"=>\$lncRNA_fa,
	"oneStepOnly:s"=>\$oneStepOnly,
	"i:s"=>\$inter,
	"e:s"=>\$exp,
	"od:s"=>\$odir,
	"cfg:s"=>\$cfg,
)or &USAGE;
&USAGE unless ($mRNA_gtf and $lnc_gtf and $mRNA_fa and $lncRNA_fa and $exp and $odir );

system "mkdir -p $odir" unless (-d $odir);
$odir = abs_path($odir);
my $sh_dir="$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
$cfg=abs_path($cfg);
$step=$step||1;
$inter= $inter||0;
&log_current_time("$Script start...");
#------------------------------------------------------------------------------
my $cmd;
my $list="$odir/../Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list";

if ($step==1){
	#len exon stat
	$cmd = "perl $Bin/mRNA_lnc_stat.pl  -i $mRNA_gtf -o $odir -k mRNA \n";
	$cmd .= "perl $Bin/mRNA_lnc_stat.pl -i $lnc_gtf  -o $odir -k lncRNA\n";
	&step_cmd_process($cmd,"1.mRNA_lnc_stat.sh",$sh_dir);
	$step++ unless ($oneStepOnly);
}
if($step==2){
	#get orf seq
	$cmd = "perl $Bin/Getorf_Extract.pl -fa $mRNA_fa -od $odir \n";
	$cmd .= "perl $Bin/get_lnc_orf_Extract.pl -fa $lncRNA_fa -od $odir \n";
	&step_cmd_process($cmd,"2.orf_stat.sh",$sh_dir);
	$step ++ unless ($oneStepOnly);#unless $inter==1;
}
if($step==3){
	$cmd = "perl $Bin/isoform-20.pl -merged $mRNA_gtf -lnc_gtf $lnc_gtf -od $odir \n";
	&step_cmd_process($cmd,"3.isoform_stat.sh",$sh_dir);
	$step ++ unless ($oneStepOnly);
}
if ($step==4){
	$cmd ="perl $Bin/Compare_Expression.pl $list $exp  $odir";
	&step_cmd_process($cmd,"4.expression_stat.sh",$sh_dir);
	$step ++ unless ($oneStepOnly);
}
if ($step==5){
	$cmd = "perl $Bin/proportion_lncRNA_gene.pl --in1 $mRNA_gtf --in2 $lnc_gtf --outdir $odir/chr_num_stat/";
	&step_cmd_process($cmd,"5.proportion_lncRNA_gene.sh",$sh_dir);
	$step ++ unless ($oneStepOnly);
}
if ($step==6){
	#inter action analysis of miRNA and lncRNA
	$cmd = "perl $Bin/lncRNA_mRNA_DEG.pl -i $odir/../Basic_Analysis -deg  $odir/../DEG_Analysis -lncdeg $odir/../Lnc_Diff_Analysis -cfg $cfg -od $odir/lncRNA_mRNA_DEG \n";
	print "perl $Bin/lncRNA_mRNA_DEG.pl -i $odir/../Basic_Analysis -deg  $odir/../DEG_Analysis -lncdeg $odir/../Lnc_Diff_Analysis -cfg $cfg -od $odir/lncRNA_mRNA_DEG \n";
	&step_cmd_process($cmd,"6.lncRNA_mRNA_DEG.sh",$sh_dir);
}
if ($step==7 && defined $inter){
	#inter action analysis of miRNA and lncRNA
	$cmd = "perl $Bin/miRNA_lncRNA_interaction.pl -mi $miRNA_fa -lnc $lncRNA_fa -od $odir/inter_action \n";
	&step_cmd_process($cmd,"4.inter_action.sh",$sh_dir);
}
#################################################################################################
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_; 
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
#    &log_current_time("CMD: $cmd");
#
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

sub log_current_time {
	my ($info) = @_;
	my $curr_time = date_time_format(localtime(time()));
	print "[$curr_time] $info\n";
}
sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
-----------------------------------------------------------------------------------------------
	Program: $Script
	Contact: renhd\@biomarker.com.cn
	   Date:2015/3/13
	  Usage:
	  		-m	mRNA_gtf,   gtf file
			-l	lncRNA_gtf, gtf file
			-mfa 	mRNA_fa,    fasta file
			-lncfa  lncTNA_fa,  fasta file
			-e      all_trans_fpkm.xls
			-cfg	detail.cfg
			-od	output dir	diretory 
			-h      help documents

			 ---------------------------------------------------------------------
			|-step   anslysis step [option]
			|	1 exon and len stat
			|	2 orf stat
			|	3 expression stat
			|	4 miRNA and lncRNA cor_analysis [optional]
			|-i      do miRNA lncRNA inter action analysis [option] value: 1
			|	-mifa   miRNA_fa,   fasta file [optional :miRNAºÍlncRNA»¥×÷·ÖÎö]
-----------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

