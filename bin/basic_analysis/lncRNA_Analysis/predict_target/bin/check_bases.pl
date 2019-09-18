#!#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$out);

GetOptions(
	"i:s" => \$in,
	"o:s" => \$out,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in and $out);


open (IN, $in) or die $!;
open (OUT,">$out") or die $!;
while (<IN>) {
	chomp;
	my ($id,$seq) = split /\t/,$_;
	#print OUT "$seq\n";
	my @bases = split //,$seq;
	#print OUT "@bases";
	foreach my $base(@bases) {
		if ($base =~ /[^atcgATCG]/) {
			$base =~ s/$base//;
		}		
	}
	$seq = join ("",@bases);
	print OUT "$id\t$seq\n";
}
close IN;
close OUT;




#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei <yaob\@biomarker.com.cn> 
      Date: 

     Usage:
		-indir	<input directory>
		-od	<output directory>
		-h	help documents

   Example:
            perl $Script  -i in -o out

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
