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
my ($in,$out);

GetOptions(
    "in:s" => \$in,
	"out:s" => \$out,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in and $out);

########### fomrat the input lncRNA file
open (IN, $in) or die $!;
open (OUT,">$out") or die $!;
$/=">";
<IN>;
while (<IN>) {
	chomp;
	my @lines = (split /\n/,$_,2);
	my $id = (split /\s+/,$lines[0])[0];
	$lines[1] =~ s/\s+//g;
	print OUT ">$id\t$lines[1]\n";
}
$/="\n";

#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao bei <yaob\@biomarker.com.cn> 
      Date: 

     Usage:
            

   Example:
            perl $Script -in lnc_filter_final.fa -out format_lnc.fa

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
