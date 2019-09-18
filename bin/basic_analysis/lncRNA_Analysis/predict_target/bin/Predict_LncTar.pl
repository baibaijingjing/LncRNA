#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")}; 
my $BEGIN_TIME=time();


# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($gene,$lncSeq,$step,$sbs,$odir,$div,$q,$m);

GetOptions(
    "gene:s"=>\$gene,
	"lncSeq:s"=>\$lncSeq,
	"odir:s"=>\$odir,
	"div:s"=>\$div,
	"step:i"=>\$step,
	"sbs:s"=>\$sbs,
    "q:s"=>\$q,
    "m:s"=>\$m,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($gene and $lncSeq and $odir);
$div ||=1000;
$step ||=1;
$q ||= "general.q";
$m ||="15G";
my ($sTime,$eTime,$stime,$etime,$sTimeG, $eTimeG, $stimeG, $etimeG);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings);
($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$second = "0$second" if($second =~ /^\d$/);
$sTimeG = "$hour:$minute:$second";
$stimeG = time;

my $notename=`hostname`;
chomp $notename;


$odir=abs_path($odir);
mkdir $odir unless(-d $odir);
my $database="$odir/database";
mkdir $database unless (-d $database);
my $sh_dir="$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
my $Cut_dir="$odir/CUT";
mkdir $Cut_dir unless (-d $Cut_dir);
my $INPUT="$odir/INPUT";
mkdir $INPUT unless (-d $INPUT);
my $OUTPUT="$odir/OUTPUT";
mkdir $OUTPUT unless (-d $OUTPUT);
####################################
my $formatdb=$config{formatdb} || "formatdb";
my $blastall=$config{blastall} || "blastall";
my $cmd;
######### blast
if ($step==1){
	&start();
	open (IN,$gene) or die $!;
	open (OUT,">$odir/database/mRNA.fa") or die $!;
	$/=">";
	<IN>;
	while (<IN>){
		chomp;
		my ($one,$seq)=(split /\n/,$_)[0,1];
		my ($id)=(split /\s+/,$one)[0];
		$seq=~s/N//g;
		my @tem=split //,$seq;
		my $len=@tem;
		next if ($len<10);
		print OUT ">$id\n$seq\n";
	}
	close IN;
	close OUT;
	$/="\n";
	
	$cmd = "cd $database\n";
	$cmd .="$formatdb -i mRNA.fa -p F\n";
	$cmd .="$blastall -b 100 -v 100 -p blastn -d $database/mRNA.fa -i $lncSeq -S 2 -e 1e-5 -F F -m 8 -a 4 -o $odir/blast.pair\n";
	&step_cmd_process($cmd,"1.blast.sh",$sh_dir);
	&end();
	$step++ unless($sbs);
}

###### Cut blast.pair
if ($step==2){
	&log_current_time("Cut blast pairs start..");
	`less $odir/blast.pair |awk '{print \$1\"\\t\"\$2}'|uniq >$odir/tmp.pair`;
	#`rm $odir/blast.pair`;
	#`mv $odir/tmp.pair $odir/blast.pair`;
	my $blast="$odir/tmp.pair";
	&CUT($blast,$Cut_dir,$div);
	&log_current_time("Cut blast pairs end..");
	$step++ unless($sbs);
}
	
if($step==3){
	&start();
	print "===makepair start===\n";
	my @files=glob("$Cut_dir/*");
	open (SH,">$sh_dir/makepair.sh") or die $!;
	foreach my $file(@files){
		my $name=basename($file);
		print SH "perl $Bin/makepair.pl -lncSeq $lncSeq -mRNA $gene -blastPair $file -out $INPUT/$name.pair ";
		print SH "&& perl $Bin/LncTar.pl -p 2 -f $INPUT/$name.pair -d -0.10 -s F -o $OUTPUT/$name.test2\n";
	}
	close SH;
	&Cut_shell_qsub("$sh_dir/makepair.sh",18,$m,$q);
	&Check_qsub_error("$sh_dir/makepair.sh");
	&end();
	print "===makepair end===\n";
	$step++ unless($sbs);
}

if ($step==4){
	print "===Result extract start===\n";
	&start();
	my @outfiles=glob("$odir/OUTPUT/*.test2");
	`cat @outfiles >$odir/total.test2`;
	open (IN,"$odir/total.test2") or die $!;
	open (OUT,">$odir/LncTar_basepair.target") or die $!;
	print OUT "#LncRNA_ID\tTarget_mRNA_ID\n";
	my %lnc_mRNA;
	while (<IN>){
		chomp;
		my ($lnc,$mRNA)=(split /\s+/,$_)[0,2];
		next if (/^Query/);
		push @{$lnc_mRNA{$lnc}},$mRNA;
	}
	foreach my $k(keys %lnc_mRNA){
		my $num=@{$lnc_mRNA{$k}};
		print OUT "$k\t",join(";",@{$lnc_mRNA{$k}}),"\n";
	}
	close IN;
	close OUT;
	&end();
	print "===Result extract end===\n";
}
	



sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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


sub start{
	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$second = "0$second" if($second =~ /^\d$/);
    $sTime = "$hour:$minute:$second";
    $stime = time;
    print STDERR "started: $sTime\n"
}

sub end{
    $etime = time - $stime;
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $second = "0$second" if($second =~ /^\d$/);
    $eTime = "$hour:$minute:$second";

    print STDERR "
	ended: $eTime
	total:", int($etime / 3600),"h:",int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
}


sub CUT{
	my ($blast,$od,$div) =@_;
	my $out=basename($blast);
	my $num=1;
	my $cout=1;
	open (IN,$blast) or die $!;
	open (OUT,">$od/$out.cut.$num") or die $!;
	while(<IN>){
		chomp;
		if($cout>$div*$num){
			close OUT;
			$num++;
			open(OUT,">$od/$out.cut.$num") or die $!;
		}
		print OUT "$_\n";
		$cout++;
	}
	close OUT;
	close IN;
}


=c
sub CUT_FILE{
	my $in = shift;
	my $od = shift;
	my $div=shift;
	
	my $line = `less -S $in |wc -l `; chomp $line;
	if ($line>$div) {
		my @files=glob "$od/$in.div*";
		foreach (@files) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@files=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$in" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$od/$in.div.$div_index" || die;
			if ($line_num<=$div) {
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
		
	}
}
=cut


sub step_cmd_process {
	my ($cmd,$sh_name,$sh_dir) = @_;
	my $sh_file = "$sh_dir/$sh_name";
	my $log_file = "$sh_file.log";
	my $flag = 0;
	my $start_time = time();
	&log_current_time("$sh_name start...");
	if (-e $sh_file) {
		system "cat $sh_file >>$sh_file.bak";
		open (SH,">$sh_file") or die "$!:$sh_file\n";
		print SH "$cmd\n";
		close SH;
	}else {
		open (SH,">$sh_file") or die "$!:$sh_file\n";
		print SH "$cmd\n";
		close SH;
	}
	$flag = system ("sh $sh_file >$log_file");
	if ($flag !=0) {
		&log_current_time("Error: command failed: $cmd");
		exit(1);
	}else {
		my $escaped_time = (time()-$start_time)."s";
		&log_current_time("$sh_name done, escaped time: $escaped_time.");
	}
}

sub log_current_time{
	my ($info) = @_;   # get parameter
	my $curr_time = date_time_format(localtime(time()));      # get current time with string
	print "[$curr_time] $info\n";  # print info with time
}

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
   Contact: Yao Bei <yaob\@biomarker.com.cn> 
      Date: 2015-10-16

     Usage:
		-gene		<FILE> mRNA.fa file
		-lncSeq		<FILE>	LncRNA fa file
		-odir		<directory> output directory
		-div		<blast result division>	(default 1000)
		-step		
			1: blast
			2: cut blast pairs
			3: makepair and LncTar
			4: result extract
		-sbs		step by step
		-h		help documents

   Example:
            perl $Script -lncSeq lnc.fa -gene mRNA.fa -odir LncTar 

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

