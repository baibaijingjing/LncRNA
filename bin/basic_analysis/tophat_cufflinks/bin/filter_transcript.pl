#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my ($in1,$in2,$out);
GetOptions(
    "in1:s" =>\$in1,
    "in2:s" =>\$in2,
    "out:s"=>\$out,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in1 and $in2 and $out);
open (IN, "$in1") or die "$!";
open (OUT, ">$out");

my %hash;
while (<IN>){
    chomp;
    my @tem=split('\s+',$_);
    my $sep=@tem[3..$#tem];
    
    
    my $num=$sep=s/(q\d+\:)/$1/g;
    #print "$num\n";
    if ($num >= 2){
        #print OUT "$_\n";
        $hash{$tem[0]}=1;
    }
#print Dumper(\$hash);
    
}
open IN1, "$in2" or die "$!";
while (<IN1>){
    chomp;
    my @tem=split("\"",$_);
    if (exists $hash{$tem[3]}){
        print OUT "$_\n";
    }
    
}
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: Script
   Version: version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            -in1      <FILE>  file.tracking
            -in2      <FILE>   compare.gtf 
            -out     <file>  output compare.tracking
            -h                 help documents

   Example:
            perl Script -in1 all.tracking -in2 file.gtf -out compare.tracking
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
