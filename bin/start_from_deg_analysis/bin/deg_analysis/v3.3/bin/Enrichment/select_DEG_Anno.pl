#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $o);


my $head_line;
my $limit=0;
my %DEG;
open (IN,"$deg") or die $!;
while (<IN>) {
	$_=~s/\s+$//;
	$head_line=$_ if /^\#/;
	next if /^\#/;
	my ($name,$val)=split/\t/,$_,2;
	$DEG{$name}=$val;
}
close IN;
open (OUT,">$o") or die $!;

my %Data;

open I,"$fIn/Integrated_Function.annotation.xls" or die $!;
while(<I>){
        chomp;
        if( $.==1){
		my @tm=split/\t/;
		shift @tm;
                print OUT "$head_line\t".join("\t",@tm[0..$#tm])."\n";
                next;
        }
        next if (/^\s*$/);
        my @tmp=split/\t/;
        if(exists $DEG{$tmp[0]}){
                print OUT "$tmp[0]\t$DEG{$tmp[0]}\t".join("\t",@tmp[1..$#tmp])."\n";
        }
}

close(OUT);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

################################################################################################################

sub max{#&max(lists or arry);
	#���б��е����ֵ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#���б��е���Сֵ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.17
Usage:
  Options:
  -i    <file>  input file,All_Database_annotation.xls,forced 
  
  -deg  <file>  deg file,forced

  -o    <file>  output file,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
