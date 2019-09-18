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
my ($lncSeq,$mRNA,$blastPair,$pair);

GetOptions(
    "lncSeq:s" =>\$lncSeq,
    "mRNA:s" =>\$mRNA,
	"blastPair:s"=>\$blastPair,
	"out:s"=>\$pair,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($lncSeq and $mRNA and $blastPair and $pair);



my %LncRNASeq;
my %mRNASeq;
open (IN,$lncSeq) or die $!;
$/=">";
<IN>;
while(<IN>){
	chomp;
	my ($one,$seq)=(split /\n/,$_)[0,1];
	my $id=(split /\s+/,$one)[0];
	$seq=~s/\s*//g;
	$seq=~s/N//g;
	$LncRNASeq{$id}=$seq;
}
close IN;

open (IN ,$mRNA) or die $!;
<IN>;
while(<IN>){
	chomp;
	my ($one,$seq)=(split /\n/,$_)[0,1];
	my ($id)=(split /\s+/,$one)[0];
	$seq =~ s/\s*//g;
	$seq=~s/N//g;
	$mRNASeq{$id}=$seq;
}
close IN;
$/="\n";

open(IN,$blastPair) or die $!;
open (OUT,">$pair") or die $!;
while(<IN>){
	chomp;
	my ($lnc,$mRNA)=(split /\s+/,$_)[0,1];
	print OUT "$lnc\t$LncRNASeq{$lnc}\t$mRNA\t$mRNASeq{$mRNA}\n";
}
close IN;
close OUT;





#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei <yaob\@biomarker.com.cn> 
      Date: 2015-10-16

     Usage:
		-lncSeq		<FILE>	LncRNA fa file
		-mRNA	<FILE>	mRNA fa file
		-blastPair	<FILE>	blast output
		-out		<FILE>	LncRNA and mRNA pair
		-h		help documents

   Example:
            perl $Script -lncSeq lnc.fa -mRNASeq mRNA.fa -blastPair blast.out -out Lnc_mRNA.pair

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
