#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($k,$l,$m,$tg,$out);
GetOptions(
    "l:s" =>\$l,
    "m:s" =>\$m,
    "tg:s" =>\$tg,
    "out:s"=>\$out,
    "k:s"=>\$k,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($l and $m and $tg and $out);
mkdir $out unless -d $out;
my %lnc;
my %m;
my %tar;
open IN,$l or die $!;
while (<IN>){
	chomp;
	next if (/^#/ or /^$/);
	my @tem=split /\s+/,$_;
	my $lncid=$tem[0];
	$lnc{$tem[0]}=1;

}
close IN;
open IN1,$tg or die $!;
while (<IN1>){
        chomp;
        next if (/^#/ or /^$/);
        my @te=split /\s+/,$_;
        #my $lncid=$tem[0];
	if (exists $lnc{$te[0]}){
		if ($te[1]=~/;/){
			my @targ=split /;/,$te[1];
			foreach my $t (@targ){
				#print "$t";
				if (!exists $tar{$t}){
					$tar{$t}="lncRNA_target";
				}
			}
		}else{
			$tar{$te[1]}="lncRNA_target";
		}
	}
        foreach my $g (keys%tar){
		#print "$g\n";
	}

}
close IN1;
open OUT,">$out/$k.lncRNA_DEG_targene.xls" or die $!;
my $head_line;
open IN2,$m or die $!;
while (<IN2>){
        chomp;
        if (/^#/){
                $head_line=$_ ;
                my @head=split /\s+/,$head_line;
                my $h=join "\t",@head[1..$#head];
                print OUT "#ID\tfeature\t$h\n";
        }
	next if (/^$/);
        my @ge=split /\s+/,$_;
        if (exists $tar{$ge[0]}){
		my $str=join "\t",@ge[1..$#ge];
		print OUT "$ge[0]\t$tar{$ge[0]}\t$str\n";
	}

}
close OUT;
close IN2;



sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script

     Usage:
            -l      <FILE>  L01_vs_L02.DEG_final.xls
            -m      <FILE>  L01_vs_L02.final.xls
            -tg     <file>  novel_lncRNA_target.xls	
            -out    <dir> 
            -k      L01_vs_L02
            -h      help documents

   Example:
            perl Script -l L01_vs_L02.DEG_final.xls -m L01_vs_L02.final.xls -tg novel_lncRNA_target.xls -out input
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
