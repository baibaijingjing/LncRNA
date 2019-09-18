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
my ($fa, $gtf, $od, $step, $key, $type, $db, $log,$cuffnorm);
GetOptions(
				"help|?" =>\&USAGE,
				"fa:s"  =>\$fa,
				"gtf:s"  =>\$gtf,
				"od:s"   =>\$od,
				"key:s"    =>\$key,
				"type:s"=>\$type,
				"db:s"=>\$db,
				"cuffnorm:s"=>\$cuffnorm,
				) or &USAGE;
&USAGE unless ($fa and $gtf and $type and $db and $od and $cuffnorm) ;
################################

my $notename = `hostname`; chomp $notename;
#$cfg = &ABSOLUTE_DIR($cfg);
$fa = &ABSOLUTE_DIR($fa);
$gtf = &ABSOLUTE_DIR($gtf);
$cuffnorm=&ABSOLUTE_DIR($cuffnorm);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");
$type = $type || "MAM";
$db = $db || "ve";

#$step = $step || 1;

#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/lncRNA_identfy.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $fa\n";
print $log "detail config file:  $gtf\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins 
#==================================================================

my $GTF2FA_BIN = "/share/nas2/genome/biosoft/tophat/2.0.7/gtf_to_fasta";


#my %para;
#open (IN,"cat $cfg1 $cfg2|") || die "$!\n";
#while (<IN>) {
#        chomp;
#        s/\r$//;s/^\s+//;s/\s+$//;
#        next if (/^\#/ || /^$/);
#
#        my @tmp=split /\s+/,$_;
#        $para{$tmp[0]}=$tmp[1];
#}


#==================================================================
open OUT1,">$od/work_sh/Lnc_identify.sh" || die;	
&MKDIR("$od/code_filter");
#open OUT1,">$od/work_sh/Lnc_filter.sh" || die;
#	open OUT1,">$od/work_sh/1.Lnc_identify.sh" || die;
	print OUT1 "perl $Bin/stat.class_code.pl -i $gtf -o $od \n";

	print OUT1 "perl $Bin/lnc_merger_filter.pl -gtf  $gtf  -out $od/merged_filter.gtf && ";
	######################????????????????????????????????????????????????????????????#####################################
	
	print OUT1 "perl $Bin/gtf2fa.pl  -gtf $od/merged_filter.gtf  -fa $fa  -o $od/merged_filter.fa && ";
	my $lnc_fa="$od/merged_filter.fa";
	#&MKDIR("$od/Lnc_filter/code_filter");
	print OUT1 "perl $Bin/Lnc_predict.v2.pl -fa $lnc_fa -type $type -db $db -od $od/code_filter -key lnc_code_filter && ";
	print OUT1 "grep -v 'NA' $od/code_filter/list.txt |cut -f 1 |awk -F ':' '{print \$1}' >$od/lnc_id.list  && ";
	print OUT1 "perl $Bin/lncRNA_fpkm_filter.pl --indir $cuffnorm --lnc_list $od/lnc_id.list --od $od/lnc_filter_id.list";
	print OUT1 "perl $Bin/abstract_gtf_seq_by_transid.pl -i $od/lnc_filter_id.list -gtf $od/merged_filter.gtf -o $od/filter_final.gtf && ";
	print OUT1 "perl $Bin/gtf_to_gff.pl -i $od/filter_final.gtf -o $od/filter_final.gff && ";
	print OUT1 "perl $Bin/lnc_gtf2fa.pl -fa $fa -gtf $od/filter_final.gtf -o $od/lnc_filter_final.fa -gff $od/filter_final.gff &&  " ;
	#print OUT1 "perl $Bin/../predict_target/bin/get_mRNA_by_gff.pl -genome $para{Ref_seq} -gff $od/Lnc_filter/filter_final.gff -o $od/Lnc_filter/lnc_filter_final.fa &&" ;
	if ($key){
		print OUT1 "perl $Bin/Precursor_analysis.pl -fa $od/lnc_filter_final.fa -od $od/ -key $key && ";
	}	
	#print OUT1 "perl $Bin/../lnc_predict/Precursor_analysis.pl -fa $od/Lnc_filter/lnc_filter_final.fa -od $od/Lnc_filter -key $para{Key} && ";
	print OUT1 "perl $Bin/lnc.stat.class_code.pl -i $od/filter_final.gtf -o $od ";
	close OUT1;
	
	&Cut_shell_qsub("$od/work_sh/Lnc_identify.sh",18,"15G","general.q");
	&Check_qsub_error("$od/work_sh/Lnc_identify.sh");


	


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		if ($notename=~/login\-0\-4/) {
			system "perl /share/nas31/niulg/bin/submit_sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "perl /share/nas31/niulg/bin/submit_sge.pl $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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
				system "perl /share/nas31/niulg/bin/submit_sge.pl $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "perl /share/nas31/niulg/bin/submit_sge.pl $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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
#	rmdir($dir) if(-d $dir);
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

Description:
	This program is a Procedure deal with lncRNA prediction 


Usage:
	-fa      the genome fasta sequence,must be given;
	-gtf	the compare gtf file ,must be given;
	-type	DataBase For CNCI predict:"ve" for vertebrate species, "pl" for plat species;,must be given;
	-db 	DataBase For CPC predict:"PLN,BCT,INV,MAM,PHG,PRI,ROD,SYN,UNA,VRL,VRT,ENV",must be given;
	-od 	output dir ,must be given;
	-key	lncRNA Precursor analysis
	-cuffnorm	cufffnorm dir  
USAGE
	print $usage;
	exit;
}
