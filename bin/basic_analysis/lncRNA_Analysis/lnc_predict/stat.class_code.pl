#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my ($in,$odir);
GetOptions(
	        "help|?"=>\&USAGE,
			"o:s"=>\$odir,
			"i:s"=>\$in,
)or &USAGE;
&USAGE unless ($odir and $in );

system "mkdir -p $odir" unless (-d $odir);
$odir = abs_path($odir);
&log_current_time("$Script start...");

#stat class code
my %class;
open IN ,"$in";
while (<IN>){
	chomp;
	next if /^$/ || /^#/;
	my @tmp=split /\t/,$_;
	$tmp[8]=~/transcript_id \"(.+?)\".*?class_code \"(.+?)\"/;
#	$ID=~/gene_id \"(.+?)\".*?transcript_id \"(.+?)\"/;
	my $transid=$1;
	my $classcode=$2;
	if (exists $class{$transid}){
		next;
	}else {
		$class{$transid}=$classcode;
	}
}
close IN;
open OUT ,">$odir/all_class_code.xls";
for my $id (sort keys %class){
	print OUT "$id\t$class{$id}\n";
}
close OUT;
my %code;
#id list
open IN ,"$odir/all_class_code.xls";
while(<IN>){
	chomp;
	next if /^$/ || /^#/;
	my $c=(split /\t/,$_)[1];
	if (exists $code{$c}){
		$code{$c}++;
	}else{
		$code{$c}=1
	}
}
close IN;
#stat
open OUT,">$odir/plot.class_code.stat";
for my $key (sort keys %code){
	print OUT "$key\t$code{$key}\n";
}
close OUT;
`$config{Rscript} $Bin/simpleBar.r --infile $odir/plot.class_code.stat --outfile $odir/class_code.stat.png --x.col 1 --y.col 2 --x.lab type --y.lab Number `;

#################################
sub log_current_time{
	my ($info) = @_;
	my $curr_time = date_time_format(localtime(time()));
	print "[$curr_time] $info\n";
}
sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

###############################
sub USAGE {
	my $usage=<<"USAGE";
------------------------------------------------------------------------------------------------
	Program: $Script
	Contact: renhd\@biomarker.com.cn
	   Date:
	  Usage:


		-i input file ,mergerd.gtf 
		-o output diretory
  		-h      help documents
------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit();
}
