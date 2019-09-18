#!/usr/bin/perl -w
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$Bin";
use predict_one_duplex;

use strict;
use warnings;

if (@ARGV!=2) {
	print "Usage: <RNAhybrid_aln.file> <svm_out> \n";
	exit;
}

my ($aln,$out)=@ARGV;
my $name_aln=basename $aln;
#>cel-let-7      22      F13D11.2.1|F13D11.2.1   1458    1265    -28.2   0.052460
#miRNA  3' UUGAUAUGUUGGAU-GAUGGAGU- 5'
#           |:|||||||||   | ||||||    
#target 5' UAUUAUACAACCGUUCCACCUCAA 3'
open IN,"$aln" || die $!;
open OUT,">$out" || die $!;
while(<IN>) {
  chomp;
  if(/^>/) {
    my $info=$_;
    if(<IN> =~ /^miRNA\s+3'\s+(\S+)\s+5'/) {
      my $miRNA = $1;
      my $match = <IN>;
      chomp $match;
      if(&GU_mismatch($match)) {next;}
      if(<IN> =~ /^target\s+5'\s+(\S+)\s+3'/) {
        my $target = $1;
        if((my $prediction = predict_one_duplex($target, $miRNA,$name_aln)) >= 0) {
          print OUT join("\n", $info."\t".$prediction, 
                           "miRNA  3' ".$miRNA. " 5'",
                           $match,
                           "target 5' ".$target." 3'")."\n\n";
        }
      }
    }
  }
}
close IN;
close OUT;

#seed region is 1 to 9 nt next to 5' end
#two consecutive GU matches are not allowed
sub GU_mismatch {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
  my $i=length($string)-1;
  my $count=0;
  while($i>=length($string)-9) {
    if(substr($string,$i,1) eq ':') {$count++;}
    else {$count=$count+0;}
    if($count>=2) {return $count;}
    $i--;
  }
	return 0;
}
