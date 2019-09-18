#!/usr/bin/perl -w
#use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($cfg1, $cfg2, $od, $step, $log, $oneStepOnly);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg1:s"  =>\$cfg1,
				"cfg2:s"  =>\$cfg2,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($cfg1 and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

$step = $step || 1;
my $thread ||= 10;
#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/Bwa.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $cfg1\n";
print $log "detail config file:  $cfg2\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins 
#==================================================================


my $GFFREAD_BIN     = $config{GFFREAD_BIN};     # 2014-12-17 ~ 
my $bwa_dir = $config{bwa};
my $CUFFQUANT_BIN   = $config{CUFFQUANT_BIN};   # 2015-08-04 
my $CUFFNORM_BIN    = $config{CUFFNORM_BIN};    # 2015-08-04  
my $samtools = $config{samtools};
#==================================================================
# load config file 
#==================================================================

my %total_read;
my %para;
my %sample;

open (IN,"cat $cfg2 $cfg1|") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	
	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/Sample/) {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;


#==================================================================
# pipeline 
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome;
my $idx_prefix;
my $gtf;
my $gff;

if ($step!=1) {
    $genome = basename($para{Ref_seq});
    $idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
}
if ($step==1) {
	print STDOUT "=== check and build bwa index ===\n";
	print $log  "=== check and build bwa index ===\n";

#	&MKDIR("$od/Ref_Genome")
	mkdir("$od/Ref_Genome") if(!-d "$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome = basename($para{Ref_seq});
    $idx_prefix = $para{Ref_seq};
#   $idx_prefix =~s/.fa$//;

 #   if (!-f "$idx_prefix.1.bt2" and !-f "$idx_prefix.1.bt2l" ) {
	if (-e "$idx_prefix.sa" and -e "$idx_prefix.pac" and -e "$idx_prefix.ann" and -e "$idx_prefix.amb" and -e "$idx_prefix.bwt" ){
        system "ln -s $idx_prefix* ./";
	}
	elsif(! -e "$od/Ref_Genome/$genome.sa" or !-e "$od/Ref_Genome/$genome.pac" or !-e "$od/Ref_Genome/$genome.ann" or !-e "$od/Ref_Genome/$genome.amb" or !-e "$od/Ref_Genome/$genome.bwt") {	
		system "ln -s $idx_prefix ./ ";
        system "$bwa_dir index -a bwtsw $genome";
    }
	################################## gff2gtf
    $gff = basename($para{Ref_ann});
    $gtf = basename($para{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gtf";
	chdir "../";

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";

}

#################################### 
#
# step 2: Align the RNA-seq read to genome using Bwa
#
#########

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using Bwa  ===\n";
	print STDOUT "Bwa mapping shell file: $od/work_sh/Bwa.sh\n";
	print $log "=== Align the RNA-seq read to genome using Bwa  ===\n";
	print $log "Bwa mapping shell file: $od/work_sh/Bwa.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/Bwa.sh" || die;
	&MKDIR("$od/Bwa");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Bwa/$sam");
		print OUT "cd $od/Bwa/$sam && ";
		print OUT "$bwa_dir mem -T 19 -t $thread $od/Ref_Genome/$genome $sample{$sam}{FQ1} $sample{$sam}{FQ2} > $sam.sam";
		print OUT " && $samtools view -bS $od/Bwa/$sam/$sam.sam -o $od/Bwa/$sam/$sam.bam &&\n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type}||="middle.q";
	#qsubOrDie("$od/work_sh/Bwa.sh","$para{Queue_type}",18,"$para{Memory}");

	&Cut_shell_qsub("$od/work_sh/Bwa.sh",18,"$para{Memory}","$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/Bwa.sh");
	

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}


#################################### 
#
# step 4: Cat All Samples SAM file
#
########

if ($step==3) {
	my $catCMD = "cat ";
	my $flag = 0;
	foreach my $sam (sort keys %sample) {
		open IN , "$od/Bwa/$sam/$sam.sam" or die $!;
		open OUT ,">$od/Bwa/$sam/$sam.update.sam";
		while(<IN>)
		{
			if($flag==0 && $_=~/^@/)
			{
				print OUT $_;
			}
			elsif($flag!=0 && $_=~/^@/)
			{
				next;
			}
			else
			{
				my $line = "$sam:".$_;
				print OUT $line;
			}
		}
		$catCMD .= "$od/Bwa/$sam/$sam.update.sam ";
		close IN;
		close OUT;
		$flag++;
	}
	`$catCMD > $od/Bwa/All.sam && touch $od/Bwa/All.sam.Check` if(!-e "$od/Bwa/All.sam.Check");
	$step++ unless ($oneStepOnly);
}
if ($step==4) {
	open SH,">$od/work_sh/circRNA_identify.sh" || die;
	#$step++ unless ($oneStepOnly);
	#&MKDIR("$od/CircRNA_identify");
	print SH "perl $Bin/circRNA_analysis.pl --data_cfg $cfg1 --detail_cfg $cfg2 --gtf $od/Ref_Genome/$gtf --sam $od/Bwa/All.sam --od $od \n";
	close SH;
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type}||="middle.q";
	#qsubOrDie("$od/work_sh/Bwa.sh","$para{Queue_type}",18,"$para{Memory}");

	&Cut_shell_qsub("$od/work_sh/circRNA_identify.sh","$para{CPU}","$para{Memory}","$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/circRNA_identify.sh");
	

	#$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
	
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================
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
		system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $shell";
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
			system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $div_file";
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
Program: Bwa_Analysis Procedure
Version: $version
Contact: Mengmeng Song <songmm\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	The program will Generate SAM file

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Bwa analysis
                            3  Cat All Samples SAM file
				4  circRNA identify analysis 
							
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
