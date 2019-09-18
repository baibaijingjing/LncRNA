#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my $BEGIN_TIME=time();
my $version="1.0.0";
my ($od,$subject,$align,$col,$key);
GetOptions(
			"help|?"=>\&USAGE,
			"od:s"=>\$od,
			"subject:s"=>\$subject,
			"align:s"=>\$align,
			"col:s"=>\$col,
			"key:s"=>\$key,
			) or &USAGE;
&USAGE unless ($subject and $align and $col and $key);
$od||="./";

`mkdir $od ` unless (-d $od);

$od=abs_path($od);
$subject=abs_path($subject);
$align=abs_path($align);


#################################################

my %seq_len;
my $total_len;

open (IN,$subject)||die $!;
$/=">";
while (<IN>) {
	chomp;
	next if $_=~/^\s*$/;
	my($subject_id,$seq)=split(/\n/,$_,2);
	$subject_id=(split(/\s+/,$subject_id))[0];
	$seq=~s/\s+//g;
	my $len=length($seq);
	$seq_len{$subject_id}=$len;
	$total_len+=$len;
}
$/="\n";
close(IN);

##################################################

my %map_site;
$col=~s/c//g;
my($one,$two,$three)=split(/,/,$col);
my ($subject_id,$subject_start,$subject_end);
open (L,$align)||die $!;
while (<L>) {
	chomp;
	next if (/^\s*$/);
	next if (/^\#/) ;
	next if (/^Query/) ;
	my @alignment=split(/\s+/,);
	
	 $subject_id=$alignment[$one-1];
	 $subject_start=$alignment[$two-1];
	 $subject_end=$alignment[$three-1];

	if ($subject_start>$subject_end) {
		my $tmp=$subject_start;
		$subject_start=$subject_end;
		$subject_end=$tmp;
	}

	if (exists $map_site{$subject_id}{$subject_start}) {
		if ($subject_end>$map_site{$subject_id}{$subject_start}) {
			$map_site{$subject_id}{$subject_start}=$subject_end;
		}
	}
	else {
		$map_site{$subject_id}{$subject_start}=$subject_end;
	
	}
}
close (L);

=c
my ($num,@tmp);
$num=@tmp=sort keys %map_site;
print $num,"\n";
=cut

my %map_gap;
my ($start,$end);
foreach $subject_id (sort keys %map_site) {
	
	$start=(sort keys %{$map_site{$subject_id}})[0];
	$end=$map_site{$subject_id}{$start};
	$map_gap{$subject_id}{$start}=$map_site{$subject_id}{$start};

	foreach $subject_start (sort {$a<=>$b}keys %{$map_site{$subject_id}}) {

		if ($end>=$subject_start ) {
		if ($map_gap{$subject_id}{$start}<$map_site{$subject_id}{$subject_start}) {
		$map_gap{$subject_id}{$start}=$map_site{$subject_id}{$subject_start};
		$end=$map_site{$subject_id}{$subject_start};
			} 
		}
		else {
		$start=$subject_start;
		$map_gap{$subject_id}{$start}=$map_site{$subject_id}{$subject_start};	
		$end=$map_site{$subject_id}{$subject_start};
		}
		
	}
}



my $length=0;
my $total_map_length=0;

open OUT,">$od/$key.xls"||die $!;
print OUT "#Subject_id\tMapped_len\tTotal_len\tCoverage\n";
foreach $subject_id (sort keys %map_gap) {
	print OUT $subject_id,"\t";
	$length=0;
	foreach $start (sort {$a<=>$b} keys %{$map_gap{$subject_id}}) {
		$length+=$map_gap{$subject_id}{$start}-$start+1;
		$total_map_length+=$map_gap{$subject_id}{$start}-$start+1;
	}
	print OUT "$length\t$seq_len{$subject_id}\t";
	printf OUT "%.2f%%\n",$length/$seq_len{$subject_id}*100;
}

print OUT "Total\t$total_map_length\t$total_len\t";
printf OUT "%.2f%%\n",$total_map_length/$total_len*100;

close (OUT);



################################################################################
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

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  coverage
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2015/1/8
Usage:
Options:
	-od     output dir
	-subject      <file> fa format,forced 
	-align       <file>   alignment ,including Subject_id	Subject_start	Subject_end,forced
	-col         <str>   Column number :Subject_id,Subject_start,Subject_end,
	                     e.g. c1,c2,c3,forced
	-key         <str>  Key of output file,forced
	-h      Help

USAGE
	print $usage;
	exit;
}
