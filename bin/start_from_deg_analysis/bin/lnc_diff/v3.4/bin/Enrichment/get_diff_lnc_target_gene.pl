#!/usr/bin/perl -w
##get target gene from deg_final.xls
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
###########################################
# 										  #
#  get target gene from deg_final.xls     #
# 										  #
###########################################

# -----------GetOptions-------------------#

my ($infile,$tar,$odir);
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$infile,
	"t:s"=>\$tar,
	"od:s"=>\$odir,
)or &USAGE;
&USAGE unless ($infile and $tar and $odir);
mkdir $odir unless -d $odir;
$odir=&ABSOLUTE_DIR($odir);
my $y=basename($infile);
#$y=(split /\./,$y)[0];
#print $y;
system `perl $Bin/data_extract_by_ids.pl -idfile $infile -destfile $tar -out $odir/temp.txt`;
#
my %T;
open IN ,"$odir/temp.txt" || die $!;
open OUT,">$odir/target_gene.list" ;
print OUT "#Target_gene_ID\tfeature\n";
<IN>;
while (<IN>){
	chomp;
	my $tar_gene_ids=(split /\t/,$_)[1];
	my @ids=split /\;/,$tar_gene_ids;
	for my $id(sort values @ids){
		#print OUT  "$id\tlncRNA_target\n";
		$T{$id}="lncRNA_target";
	}
}
for my $key (sort keys %T){
	print OUT "$key\t$T{$key}\n";
}
close OUT;
close IN;
`rm $odir/temp.txt`;

#################################################################
sub ABSOLUTE_DIR{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE{
	my $usage=<<"USAGE";
ProgramName:$Script
Version:    $version
Contact:
Program Date:
Usage:
	Example:	perl get_diff_lnc_target_gene.pl -i lncRNA.DEG_final.xls -t novel_lncRNA_target.xls -od ./
	Options:
	-i <file> lncRNA.DEG_final.xls
	-t <file> novel_target_gene.xls
	-od <file> output dir,forced
	-h Help
USAGE
	print $usage;
	exit;
}
