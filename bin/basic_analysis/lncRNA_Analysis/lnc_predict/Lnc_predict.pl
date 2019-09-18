#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($fa,$od, $step, $Type,$key,$log, $oneStepOnly,$db);
GetOptions(
				"help|?" =>\&USAGE,
				"-fa:s"  =>\$fa,
				"-type:s"  =>\$Type,
				"-db:s"  =>\$db,
				"od:s"   =>\$od,
				"key:s" =>\$key,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($fa and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$fa = &ABSOLUTE_DIR($fa);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);
&MKDIR("$od/work_sh");

$step = $step || 1;
$db = $db || "PLN";

#==================================================================
# bins 
#==================================================================
my $CPC_BIN     = "/share/bioCloud/renhd/testing/software/cpc-0.9-r2/bin/run_predict_local_v2.sh";
my $CNCI_BIN       ="/share/bioCloud/renhd/testing/software/CNCI_V2/CNCI.py";
my $PFAM_BIN   = "/share/bioCloud/renhd/testing/software/PfamScan/pfam_scan.pl";
my $PYTHON = " python ";
my $PFAM_DATA = "/share/nas2/database/pfam/27.0";


#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/LncPredict.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data  file:  $fa\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# load config file 
#==================================================================

my $total_read = `grep -c '>' $fa ` ;
my %result;

print $log "There are $total_read for Long Noncoding RNA Prediction\n";

#==================================================================
# pipeline 
#==================================================================

#######################################
#
# step 1: Run CPC Prediction
#
######
if ($step == 1) {
	print STDOUT "=== Run CPC Prediction:\n===\n";
	print STDOUT "CPC Prediction shell file: $od/work_sh/Predict.sh\n";
	print $log "=== Run CPC Prediction:\n ===\n";
	print $log "CPC Prediction shell file: $od/work_sh/Predict.sh\n";
	open OUT1,">$od/work_sh/Predict.sh" || die;
	my $CPC_dir = "$od/CPC";
	&MKDIR($CPC_dir) unless -d "$od/CPC" ;
	print OUT1 "cd $CPC_dir && ";
	print OUT1 "$CPC_BIN  $fa  $CPC_dir/$key.result.txt  $CPC_dir $CPC_dir/$key.evidence.txt  $db \n";
	close OUT1;	
	$step++ unless ($oneStepOnly) ;
}

if ($step == 2) {
	print STDOUT "=== Run CNCI Prediction:\n===\n";
	print STDOUT "CNCI Prediction shell file: $od/work_sh/Predict.sh\n";
	print $log "=== Run CNCI Prediction:\n ===\n";
	print $log "CNCI Prediction shell file: $od/work_sh/Predict.sh\n";
	open OUT2,">>$od/work_sh/Predict.sh" || die;
	&MKDIR("$od/CNCI") unless -d "$od/CNCI" ;
	print OUT2 "cd $od  && ";
	print OUT2 " $PYTHON $CNCI_BIN  -f  $fa  -o $od/CNCI -p 15 -m $Type  \n";
	close OUT2;
	$step++ unless ($oneStepOnly) ;
}

if ($step == 3) {
	print STDOUT "=== Run Pfam Prediction: \n===\n";
	print STDOUT "Pfam Prediction shell file: $od/work_sh/Predict.sh\n";
	print $log "=== Run Pfam Prediction:\n ===\n";
	print $log "Pfam Prediction shell file: $od/work_sh/Predict.sh\n";
	open OUT3,">>$od/work_sh/Predict.sh" || die;
	&MKDIR("$od/Pfam") unless -d "$od/Pfam" ;
	print OUT3 "cd $od/Pfam && ";
	print OUT3 " perl  $PFAM_BIN   -translate orf  -fasta $fa  -dir $PFAM_DATA  -outfile  $od/Pfam/Pfam_result.txt  -cpu 15  \n";
	close OUT3;
#	$step++ unless ($oneStepOnly) ;
	&Cut_shell_qsub("$od/work_sh/Predict.sh",30, "15G",  "general.q");
	&Check_qsub_error("$od/work_sh/Predict.sh");
	$step++ unless ($oneStepOnly) ;
}
if ($step==4){
	$result{'-cpc'} = "$od/CPC/$key.result.txt";
	$result{'-cnci'} = "$od/CNCI/CNCI.index";
	$result{'-pfam'} = "$od/Pfam/Pfam_result.txt";
	my @cmd ;
	foreach my $k (keys(%result)){
		push (@cmd ,$k," $result{$k} ");
	}
	print @cmd ;
	` perl $Bin/lnc_predict_veen.pl @cmd -fa $fa -od $od `;
	
}

#&Cut_shell_qsub("$od/work_sh/Predict.sh",6, "15G",  "general.q");
#&Check_qsub_error("$od/work_sh/Predict.sh");

print STDOUT "All Finished\n";
print $log "All Finished\n";

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================
sub fa_size{
	my $fa = shift;
	my %fa_stat;
	open IN, $fa || die;
	$/='>';
	<IN>;
	while (<IN>) {
    chomp;
    my ($id,$seq)=split /\n+/,$_,2;
    my $seq_id=(split /\s+/,$id)[0];
    $seq=~s/\s+//g;
	$fa_stat{$seq_id} = length($seq)
	}
	close IN ;
	return %fa_stat;

}
sub LOAD_PARA {
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
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
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


sub USAGE {#
my $usage=<<"USAGE";
Program: Tophat&Scripture_Analysis Procedure
Version: $version
Contact: Zhang QiuXue <zhangqx\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat + Scripture Combination£ºDesigned for RNA Analysis with a Reference Genome

Usage:
	-fa		fasta file ,must be given;
	-type	DataBase For CNCI predict:"ve" for vertebrate species, "pl" for plat species;
	-db 	DataBase For CPC predict:"PLN,BCT,INV,MAM,PHG,PRI,ROD,SYN,UNA,VRL,VRT,ENV"
	-od		output dir, must be given;
	-key		keyword  of outfile,must be given;;
	-s		step of the program   option,default 1;
				1  Run CPC 
				2  Run CNCI
				3  Run Pfam  
				4  Venn 
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
