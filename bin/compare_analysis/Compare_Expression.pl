#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin $Script);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};

#lnc id list
if (@ARGV==0){
	print "\nperl Compare_Expression.pl lnc_target_id.list all_trans_fpkm.xls ./  \n\n";
	exit;
}

my %lnc;
open IN,"$ARGV[0]";
while (<IN>){
	chomp;
	my $id=(split (/\s+/,$_))[0];
	$lnc{$id}=1;
}
close IN;

##fpkm
my @lncRNA;
my @mRNA;
open IN,"$ARGV[1]";
open OUT1,">$ARGV[2]/lnc.FPKM.xls";
open OUT2,">$ARGV[2]/mRNA.FPKM.xls";
open OUT3,">$ARGV[2]/plot.txt";
while (<IN>){
	chomp;
	next if /^$/;
	if (/^#/){
		print OUT1 "$_\n";
		print OUT2 "$_\n";
	}else{
		my @tmp=(split (/\t/,$_));
		if (exists $lnc{$tmp[0]}){
			print OUT1 "$_\n";
			#for my $tmp(1..$#tmp){
			#	push @lncRNA,$tmp;
			#	print OUT3 "$tmp\t";
			#}
			for my $i(1..$#tmp){
				print OUT3 "lncRNA\t$tmp[$i]\n" unless ($tmp[$i]==0);
			}
		}else{
			print OUT2 "$_\n";
			for my $i(1..$#tmp){
				print OUT3 "mRNA\t$tmp[$i]\n" unless ($tmp[$i]==0);
			}
		}
	}
}
`$config{Rscript} $Bin/oneFactorBox.r --infile $ARGV[2]/plot.txt --outfile $ARGV[2]/lnc_vs_mRNA.fpkm.png --value.col 2 --x.col 1 --x.lab type --y.lab log10\\\(FPKM+1\\\) --title.lab lnc_vs_mRNA`;

close IN;
close OUT1;
close OUT2;
close OUT3;


