#!/usr/bin/perl -w
use autodie qw(:all);
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw (abs_path getcwd);
use File::Basename qw (basename dirname);
my $version="1.0";

my ($gtf,$od,$cutoff);
GetOptions(
	"h|?"=>\&USAGE,
	"gtf:s"=>\$gtf,
	"od:s"=>\$od,
	"cut:s"=>\$cutoff,
) or &USAGE;
&USAGE unless ($gtf and $od );
$cutoff ||= 2000;

$od=abs_path($od);
mkdir $od unless (-d $od);
my $CUTGTF="$od/CUTGTF";
mkdir $CUTGTF unless (-d $CUTGTF);

##################Compare.gtf cut
&CUTGTF($gtf,"$CUTGTF",$cutoff,'cut');

sub CUTGTF{
	my ($gtf,$od,$cut,$name)=@_;
	my %IDS;
	open (IN,$gtf) or die $!;
	while(my $line=<IN>){
		chomp $line;
		my $attribute=(split /\t+/,$line)[8];
		my $gene=(split /;/,$attribute)[0];
		push @{$IDS{$gene}},$line;
	}
	close IN;
	my @cus=glob "$od/*.gtf";
	foreach my $cu(@cus){
		system "rm $_";
	}

	my @aa=keys %IDS;
	my $index=0;
	LAB: for (my $i=1;;){
		my $num=0;
		open (OUT,">$od/$name.$i.gtf") or die $!;
		for ($index..$#aa){
			#print "$index\n";
			if ($num<$cut){
				print OUT join ("\n",@{$IDS{$aa[$index]}});
				print OUT "\n";
				$num++;
				$index++;
			}if ($num>=$cut){
				$num=0;
				$i++;
				close OUT;
				if ($index==$#aa+1){
					last;
				}else {
					next LAB;
				}
			}
		}
		if ($num){
			close OUT;
		}
		last;
	}
}



sub USAGE {
	my $usage =<< "USAGE";
Program: $Script
Version: $version
Usage:
	Options:
	-gtf	Compare.gtf
	-od	output dir
	-cut	cutoff,default(2000)
	-h	help
USAGE
	print $usage;
	exit;
}