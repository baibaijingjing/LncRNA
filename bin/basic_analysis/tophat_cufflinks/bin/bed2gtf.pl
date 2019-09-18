#!/usr/bin/perl

#
# Author : Ahfyth
#

use strict;
use warnings;
use	Getopt::Long;
use	Data::Dumper;
use List::Util qw(max min sum);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
my ($gtf,$bed);#,$isoform_len,$isoform_num,$isoform_cov);
GetOptions(
				"help|?" =>\&USAGE,
				"bed:s"=>\$bed,
				"gtf:s"=>\$gtf,
				) or &USAGE;
&USAGE unless ($gtf and $bed );

#$isoform_len||=200;
#$isoform_num||=	2;
#$isoform_cov||= 2;

open OUT, ">$gtf" or die( "$!" );

open IN, "$bed" or die( "$!" );
print "Loading file: $bed\n";
while( <IN> ){
	chomp;
	my $line = $_;
	# format explanation http://www.broadinstitute.org/software/scripture/book/export/html/12
	
	my @l = split /\t/,$line;	#chr start end name score strand thickStart thickEnd Rgb blockCount blockSizes blockStarts FWER Enrichment READS coverage RPKM lamda
	my $fpkm = $l[16];
	my $s = $l[1] + 1;	# adapt to 1-base, as GFF format used this method to record positions.
	my $strand = $l[5];
#	my @len = split /,/, $l[10];
#	my $exon_num = @len;
#	my $trans_len = sum(@len);
#	my $cov = $l[15];
	my $reads = $l[14];
	$strand = '.' unless $strand=~/^[+-]$/;
#	print "$exon_num\t$isoform_cov\t$isoform_len\n";
#	if ( $isoform_num <= $exon_num  and $isoform_len <= $trans_len and $isoform_cov <= $cov){
#		print OUT "$l[0]\tScripture\ttranscript\t$l[1]\t$l[2]\t.\t$strand\t.\ttranscript_id \"Scripture.$l[0].$.\"; FPKM \"$fpkm\";Exon_num \"$exon_num\";Length \"$trans_len\";Counts \"$reads\";Coverage \"$cov\" \n";
		my @len = split /,/, $l[10];
		my @starts = split /,/, $l[11];
		for(my $i=0; $i<$l[9]; ++$i )
		{
			my $exonstart = $starts[$i] + $s;
			my $exonend = $exonstart + $len[$i] - 1;
#			print OUT "$l[0]\tScripture\texon\t$exonstart\t$exonend\t.\t$strand\t.\ttranscript_id \"Scripture.$l[0].$.\"; FPKM \"$fpkm\" \n";
			print OUT "$l[0]\tScripture\texon\t$exonstart\t$exonend\t.\t$strand\t.\ttranscript_id \"Scripture.$l[0].$.\"; FPKM \"$fpkm\"\n";	
			#	}
	}
}
close IN;
close OUT;

sub USAGE{
	print << "	Usage End.";
Description:The process trans gtf file from bed file.
version:$ver
Usage:
	-bed	<STR>	bed file                         must be given;
	-gtf	<STR>	gtf Out file;

	Usage End.
		exit;
}
