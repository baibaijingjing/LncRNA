#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

my ($miRNA,$lncRNA,$p,$ts,$gn,$step,$od);
GetOptions(
		"help|?" =>\&USAGE,
		"mi:s"    =>\$miRNA,
		"lnc:s"   =>\$lncRNA,
		"cpu:s"   =>\$p,
		"score:s" =>\$ts,
		"gap:s"   =>\$gn,
		"step:s"  =>\$step,
		"od:s" =>\$od,
		) or &USAGE;
&USAGE unless ($miRNA and $lncRNA and $od) ;
#####main
my $notename = `hostname`;
chomp $notename;
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");
$p=$p||1;
$ts=$ts||2.5;
$gn=$gn||1;
$step = $step || 1;

my $startTime = GetTime();
my $RNAdecoys = "/share/bioCloud/renhd/testing/tools/bin/scripts/RNAdecoys";

if ($step ==1){
	open OUT,"$od/work_sh/interaction.sh";
	print OUT "/share/bioCloud/renhd/testing/tools/bin/scripts/RNAdecoys -s $miRNA -t -o $od/interaction.txt -ts $ts -fp 2 -tp 17 -gl 9 -gr 12 -p $p -gn $gn \n";
	close OUT;

&Cut_shell_qsub("$od/work_sh/interaction.sh",6,"5G","general.q");
&Check_qsub_error("$od/work_sh/interaction.sh");
$step++;
}

if ($step==2){
	open IN ,"$od/interaction.txt";
	open OUT,">$od/result.txt";
	$/='>';
	<IN>;
	while (<IN>){
	chomp;
	next if ( /^#/|| /^$/ );
	my ($inter,$Query,$dot,$Subject)= split /\n+/,$_,4;
	$inter =~s/\t/ /g;
	my $Q=(split (/\s+/,$Query))[2];
	$dot=~s/ //g;
	my $S=(split (/\s+/,$Subject))[2];
	#if ($Q= ~/\-\-/){print "$Q\n$S\n$dot\n" ;}
	#print "\n$Q\n$S\n$dot\n" ;
	if ($Q =~/--/ ){
			print OUT ">$inter\n$Q\n$dot\n$S\n";
			}
	}
	close IN;
	close OUT;
}
####sub
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
	if ($line>=1000) {
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
sub Check_qsub_error {## Check The qsub process if error happend
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";
	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else{
		print "$sh qsub is Done!\n";
	}
}
sub MKDIR
{ # &MKDIR($out_dir);
	    my ($dir)=@_;
	    rmdir($dir) if(-d $dir);
	    mkdir($dir) if(!-d $dir);
}
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	  
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}   
	elsif(-d $in)
	{
	chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
	    exit;
	}
    chdir $cur_dir;
    return $return;
}
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
##########################################################################
sub USAGE{
my $usage=<<"USAGE";
Program: miRNA_lncRNA_Interaction Procedure
Version: $version
Contact:Ren Hengda <renhd\@biomarker.com.cn>

Usage:
	-mi   miRNA form sRNA analysis  fasta format
	-lnc  lnc RNA  sequence         fasta format
	-od   output dir                diretory
	
	Optional:
		-cpu    number of processors use ,default 1
		-score  target penalty score, lower is better (0-5), default 2.5
		-gap    number of gaps/bulges permitted (0-5) ,default 1
		-step   analysis step
				1 RNAdecoys analysis
				2 get interaction result  for RNAdecoys output

USAGE
	print $usage;
	exit;
}
