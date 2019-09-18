#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my @Original_ARGV=@ARGV;

# ==============================================================
# Get Options
# ==============================================================
my ($cfg,$dOut,$tophat,$step,$gff);

GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"tophat:s"=>\$tophat,
				"gff:s"=>\$gff,
				"od:s"=>\$dOut,
				"step:s"=>\$step,
				) or &USAGE;
&USAGE unless ($cfg and $tophat);
#===============================================================
# Default optional value 
#===============================================================
###### software dir ###########
my $_SNP_		=	"$Bin/bin/GATK_calling_Trans_Ref.pl";
$dOut||="./SNP_Analysis ";
$step||=1;
###### global value ###########
my %para;
my %sample;
my $cmd;
my $notename = `hostname`; chomp $notename;
$dOut=abs_path($dOut);
$tophat=abs_path($tophat);
#===============================================================
# pipeline
#===============================================================
log_current_time("$Script start...");
&read_config($cfg);
&MKDIR("$dOut");
my $genome=basename($para{'Ref_seq'});

################# bioinfor_pipeline.log ###################
&show_log2("step_6: SNP analysis start.");
###########################################################

my $Rscript = $config{Rscript};

if ($step==1) {

	&MKDIR("$dOut/aligndir");
	&MKDIR("$dOut/work_sh");
	my @bam=glob("$tophat/*/accepted_hits.bam");

	open (OUT ,">$dOut/work_sh/bam_process.sh") or die $!;
	foreach my $bam (@bam) {
		print OUT "perl $Bin/bin/bin/bam_process.pl -od $dOut/aligndir -bam  $bam \n";
	}
	close (OUT);

	&Cut_shell_qsub("$dOut/work_sh/bam_process.sh",$para{'CPU'},"20G",$para{'Queue_type'});
	&Check_qsub_error("$dOut/work_sh/bam_process.sh");
	system "perl $Bin/bin/bin/sort_fa.pl -fa $para{'Ref_seq'} -od $dOut/aligndir ";
	$step++;
}
if ($step==2) {
	$cmd = "perl $_SNP_ -ref $dOut/aligndir/$genome.fa -aligndir $dOut/aligndir -ploidy $para{'ploidy'} -win $para{'window'} -clu $para{'cluster'} -QD $para{'QD'} -FS $para{'FS'} -od $dOut/SNP -doRecal $para{'Recal'} -doIndel $para{'ReAlignIndel'}   ";
	$cmd.="-qphred $para{'qphred'}" if (exists $para{'qphred'});
	$cmd.=">/dev/null";
	log_current_time("GATK calling start...");
	log_current_time("CMD: $cmd");
	system $cmd;
	log_current_time("GATK calling done.");

	##SNP sites stat
	log_current_time("pairwised SNP abstrct and SNP density plot start...");
	if (defined $gff){
		$cmd=`cat $para{'Ref_ann'} $gff > $dOut/Integrated.gff` ;
		$cmd="perl $Bin/util/get_gene_fa.pl $dOut/SNP/Genome/$genome.fa $gff $dOut/New_gene";
		print "$cmd\n";
		system $cmd;
	}
	else {
		$cmd=`cp $para{'Ref_ann'}  $dOut/Integrated.gff` ;
	}

	$cmd="perl $Bin/util/get_gene_fa.pl $dOut/SNP/Genome/$genome.fa $dOut/Integrated.gff $dOut/All_gene";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/SNP_indel_anno/snp_indel_anno.pl -id $dOut -r $dOut/SNP/Genome/$genome.fa -s $para{'Project_key'} >/dev/null ";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/util/SNP_stat.pl -snp $dOut/stat/final.snp.anno.gatk.all.list  -gff $dOut/Integrated.gff -od $dOut/stat/ >/dev/null ";
	print "$cmd\n";
	system $cmd;

	$cmd="$Rscript $Bin/bin/bin/final_SNP_type.r  $dOut/stat/All.snp_type.stat  $dOut/stat/ >/dev/null ";
	print "$cmd\n";
	system $cmd;

	$cmd="perl $Bin/util/pairwised_SNP_abstrct_density_plot.pl --ref $dOut/aligndir/$genome.fa --snp $dOut/stat/final.snp.anno.gatk.all.list --od $dOut/ >/dev/null ";
	print "$cmd\n";
	system $cmd;
	log_current_time("pairwised SNP abstrct and SNP density plot done.");
}

#######################################################################################
my $elapsed_time = (time()-$BEGIN_TIME).'s';
log_current_time("$Script done. elapsed time: $elapsed_time.");
####################################################################################################

##############################
&show_log2("step_6: SNP analysis finished.");
#close LOG;
##############################

sub read_config {#
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	while (<IN>) {
		chomp;
		s/\r$//;
		s/\s+$//;
		next if (/\#/ || /^$/);
		my @tmp=split /\s+/,$_;
		if ($tmp[0]=~m/Sample/) {
			my $fq1=<IN>;chomp $fq1;
			my $fq2=<IN>;chomp $fq2;
			my @fq_1=split /\s+/,$fq1;
			$sample{$tmp[1]}{FQ1}=$fq_1[1];
			my @fq_2=split /\s+/,$fq2;
			$sample{$tmp[1]}{FQ2}=$fq_2[1];
		}
		$para{$tmp[0]}=$tmp[1];
	}
	close IN;
}

####################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		if ($notename=~/cluster/) {
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
			if ($notename=~/cluster/) {
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
}
sub check_file {#
	my ($file) = @_;
	if (-e $file) {
		return 1;
	}else{
		return 0;
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

####################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}
#############################################################################################################

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log1("cmd")
sub show_log1()
{
    open LOG, ">$dOut/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$dOut/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
    Program: $0
    Version: $version
    Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
Discription:
      Usage:
        Options:
        -cfg    <file>  required, ref_trans.detail.cfg
        -od     <path>  optional, directory where output file produced, default [./SNP_Trans]
        -tophat     <path>  required, directory where the tophat result						(template: ./Basic_Analysis/Tophat_Cufflinks/Tophat)
        -step		<num>	optional	run from this step,default [1]
						1 :  Preparation bam and genome fa;
						2 :  SNP with GATK
        -gff  <file>  new_gene.gff ,optinal
        -help           help

USAGE
	print $usage;
	exit;
}
