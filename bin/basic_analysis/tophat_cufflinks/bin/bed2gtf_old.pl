#!/usr/bin/perl

#
# Author : Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <output.gtf> <scripture.bed> [scripture.bed ...]\n";
	print STDERR "\nThis program is designed to translate the scripture's bed file into gtf file.";
	print STDERR "\nYou can add as many bed files as you want.";
	print STDERR "\n\nIMPORTANT NOTE: The BED format defines that chromStart (column 2) is 0-based while chromEnd (column 3) is 1-based.";
	print STDERR "\nIt is your responsibility to make sure that your files satisfy this criterion.\n\n";
	exit 1;
}

my $out = shift;
open OUT, ">$out" or die( "$!" );

foreach my $bed ( @ARGV )
{
	open IN, "$bed" or die( "$!" );
	print "Loading file: $bed\n";
	while( <IN> )
	{
		chomp;
		# format explanation http://www.broadinstitute.org/software/scripture/book/export/html/12
		my @l = split /\t/;	#chr start end name score strand thickStart thickEnd Rgb blockCount blockSizes blockStarts FWER Enrichment READS coverage RPKM lamda
		my $fpkm = $l[16];
		my $s = $l[1] + 1;	# adapt to 1-base, as GFF format used this method to record positions.
		my $strand = $l[5];
		$strand = '.' unless $strand=~/^[+-]$/;
		print OUT "$l[0]\tScripture\ttranscript\t$l[1]\t$l[2]\t.\t$strand\t.\tttranscript_id \"Scripture.$l[0].$.\"; FPKM \"$fpkm\"\n";
		my @len = split /,/, $l[10];
		my @starts = split /,/, $l[11];
		for(my $i=0; $i<$l[9]; ++$i )
		{
			my $exonstart = $starts[$i] + $s;
			my $exonend = $exonstart + $len[$i] - 1;
			print OUT "$l[0]\tScripture\texon\t$exonstart\t$exonend\t.\t$strand\t.\ttranscript_id \"Scripture.$l[0].$.\"; FPKM \"$fpkm\"\n";
		}
	}
	close IN;
}
close OUT;

