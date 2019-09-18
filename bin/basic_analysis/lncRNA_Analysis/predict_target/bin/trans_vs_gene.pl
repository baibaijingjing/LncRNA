#!/usr/bin/perl -w
#use waining;
use strict;


open  IN ,"$ARGV[0]";
open OUT ,">$ARGV[1]";
while (<IN>){
	chomp;
	next if /^$/ || /^#/;
#	next unless /transcript/;
	my @tmp=split /\t/,$_;
	#next unless $tmp[2]=~/transcript/;
	$tmp[8]=~/gene_name \"(.+?)\".*?oId \"(.+?)\"/;
#	$ID=~/gene_id \"(.+?)\".*?transcript_id \"(.+?)\"/;
	my $gene=$1;
	my $trans_id=$2;
	print OUT "$gene\t$trans_id\n" ;#unless $gene=~/EN/;
}
