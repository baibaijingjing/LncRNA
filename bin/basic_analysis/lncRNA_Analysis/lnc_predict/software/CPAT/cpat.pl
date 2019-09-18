#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$hexamer,$logicmodel,$od,$cutoff);

GetOptions(
    "fa:s"=>\$fa,
	"od:s"=>\$od,
	"hexamer:s"=>\$hexamer,
	"logicmodel:s"=>\$logicmodel,
	"cutoff:s"=>\$cutoff,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($fa and $od);

$cutoff ||=0.38;

$od=abs_path($od);
mkdir $od unless (-d $od);

my $hexamer ||= "$Bin/Hexamer.tsv";
my $logicmodel ||="$Bin/logitModel.RData";
#my $CPAT_PY="/share/nas2/genome/biosoft/Python/2.7.8/bin/cpat.py";
my $CPAT_PY="/share/nas1/linhj/software/cpat.py";

my $cmd;

#$cmd="/share/nas2/genome/biosoft/Python/2.7.8/bin/python /share/nas2/genome/biosoft/Python/2.7.8/bin/cpat.py -g $fa -x $hexamer -d $logicmodel -o $od/tmp.txt \n";
$cmd="/share/nas2/genome/biosoft/Python/2.7.8/bin/python /share/nas1/linhj/software/cpat.py -g $fa -x $hexamer -d $logicmodel -o $od/tmp.txt \n";
$cmd .="perl $Bin/extract.pl -i $od/tmp.txt -cutoff $cutoff -o $od/cpat.txt\n";

print "$cmd\n";
system ($cmd);
	
#############################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei<yaob\@biomarker.com.cn> 
      Date: 

     Usage:
            -fa	merged.fa
			-od	output directory
			-hexamer	hexamer.csv
			-logicmodel	logitModel.RData
			-cutoff	the threshold of noncoding RNA (default : 0.38)

   Example:
            perl $Script -fa merged.fa -od CPAT -cutoff cutoff 

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

