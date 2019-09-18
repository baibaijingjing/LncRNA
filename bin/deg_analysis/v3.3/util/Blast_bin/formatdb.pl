#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
my $version="v1.0";
my ($input,$odir,$proornuc,$suoyinname);
GetOptions(
			"help|?"=>\&USAGE,
			"i:s"=>\$input,
			"odir:s"=>\$odir,
			"o:s"=>\$suoyinname,
			)or &USAGE;
&USAGE unless ($input and $odir);
unless (-d $odir) {
	` mkdir $odir -p `;
}
$odir = ABSOLUTE_DIR($odir);
#print $odir;die;
my $basename = basename($input);
$suoyinname ||= "F";
`cp $input $odir `;
#print "cd $odir ";die;

$/=">";
open (IN,$input)||die $!;
while (<IN>) {
	chomp;
	next if ($_=~/^\s*$/);
	my($id,$seq)=split(/\n/,$_,2);
	$seq=~s/\n//g;
	if ($seq=~/[^ATGCUN]/i) {
		$proornuc="T";
	}
	else{$proornuc="F";};
	last;
}
$/="\n";
close(IN);

chdir $odir ;
#print "formatdb -i $odir/$basename -p $proornuc -o $suoyinname -odir \n";die;
`$config{formatdb}  -i "$odir/$basename" -p $proornuc -o $suoyinname  `;

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
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

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  cloud_supply
Version:	$version
Contact:	He yang <hey\@biomarker.com.cn> 
Program Date:   2014/7/9
Usage:
Options:
-i    input file
-o    Parse options
      T - True: Parse SeqId and create indexes.
      F - False: Do not parse SeqId. Do not create indexes.
-odir Output directory
-h Help

USAGE
	print $usage;
	exit;
}
