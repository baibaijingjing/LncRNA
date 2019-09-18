#!/usr/local/bin/perl -w
# Copyright (c) BMK 2013
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2013.
# Modifier:       lium <lium@biomarker.com.cn>
# Last Modified:  2012-7-28
my $ver="1.0.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作


# get option
our %opts;
GetOptions(\%opts,"indir=s","out=s","h" );
# check option
if(!defined($opts{indir}) || !defined($opts{out}) || defined($opts{h})){
	&help();
	exit;
}
# get indir
my $indir=&ABSOLUTE_DIR($opts{indir});


###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################









############### Load geneExpression.xls files
my @geneExpress = glob "$indir/*.geneExpression.xls" ;
my %All_expression;
my %GeneLength;
my @sample;

foreach my $file (@geneExpress) {
	my $name=basename($file);
	$name=~/(.*)\.geneExpression\.xls/;
	my $key=$1;
	push @sample,$key;
	open (IN,"$file") or die $!;
	while (<IN>) {
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/ || /^Gene/);
		my @tmp=split/\s+/,$_;
		$All_expression{$tmp[0]}{$key}=int($tmp[5]);
		$GeneLength{$tmp[0]} = $tmp[1];
	}
#	print Dumper \%All_expression;
	close (IN) ;
}



############## combinate read count and output them into one file ############
open (OUT, ">$opts{out}") or die $!;
print OUT "#ID\t";
print OUT join("\t",(sort @sample)),"\tgeneLength\n";
foreach my $id (sort keys %All_expression) {
	print OUT "$id";
	my $str;
	foreach my $sam (sort @sample) {
		if (!defined $All_expression{$id}{$sam}) {
			$str.="\t"."0";
			next;
		}
		$str.="\t"."$All_expression{$id}{$sam}";
	}
	print OUT "$str\t$GeneLength{$id}\n";
}
close (OUT) ;







###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);








sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub help {
print <<"Usage End.";
Description: abstract read count data from the directory;
Version: $ver
Usage:
-indir		isoExpress files dir [ forced ]
-out		the filename of outfile [ forced ]
-h		for help
Usage End.
exit;
}

