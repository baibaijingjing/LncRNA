#!/usr/bin/perl
use strict;
use warnings;


open( IN,   "$ARGV[0]" )     || die "open file failed\n";
open( FILE, "$ARGV[1]" )   || die "open file failed\n";
open( OUT,  ">$ARGV[2]" ) || die "open or create file failed\n";
my %hash;
while (<IN>) {
	chomp;
	my @line = split/\t/,$_;
	$hash{$line[0]}=$line[1];
}
close IN;
my %ex;
while(<FILE>)
{
	next if(/#/);
	my @line = split/\t/,$_;
	if(!exists $ex{$line[0]})
	{
		$ex{$line[0]}=1;
		if(exists $hash{$line[1]})
		{
			my $length = $hash{$line[1]};
			my $ratio = (($line[3]-$line[4]-$line[5])/$length)*100;
			next if($ratio<90);
			printf OUT ( "%s\t%s\t%.2f\n", $line[0], $line[1], $ratio);
		}
	}
	
}