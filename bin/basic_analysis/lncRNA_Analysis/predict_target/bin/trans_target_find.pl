#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($q,$t,$c,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"q:s"=>\$q,
				"t:s"=>\$t,
				"c:s"=>\$c,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($q and $t and $od);

 $c||=30;
 mkdir "$od" unless -d "$od";
 $od=&ABSOLUTE_DIR($od);
 $q=&ABSOLUTE_DIR($q);
 mkdir "$od/Mid_Dir" unless -d "$od/Mid_Dir";
 mkdir "$od/work_sh" unless -d "$od/work_sh";
 mkdir "$od/Mid_Dir" unless -d "$od/Mid_Dir";
 my $limit=0;
 my @Mid_File=("Mid_File_0");
 my %T_N;
 $/=">";
 open (IN,$t) or die $!;
 open (OUT1,">$od/Mid_Dir/$Mid_File[-1]") or die $!;
 <IN>;
 while (<IN>) {
 	chomp;
 	$limit++;
 	print OUT1 ">$_";
 	$T_N{(split/\s/,$_)[0]}=(split/\n/,$_)[0];
# 	unless ($limit%100) {
	unless ($limit%5) {
 		my $new_file=@Mid_File;
 		$new_file="Mid_File_".$new_file;
 		push @Mid_File,$new_file;
 		close OUT1;
 		open (OUT1,">$od/Mid_Dir/$Mid_File[-1]") or die;
 	}
 }
 close IN;
 close OUT1;
 $/="\n";

 open (SH,">$od/work_sh/trans_target.sh") or die $!;
 foreach my $file_name (@Mid_File) {
 	print SH "/share/nas1/zhangxc/bin/Process/RNAplex/bin/RNAplex -t $od/Mid_Dir/$file_name -f 1 -q $q -c $c -e -20 > $od/Mid_Dir/$file_name.out \n";
 }
 close SH;
 if (`hostname`=~/login\-0\-4/) {
 	`sh $config{qsub_sh} --independent --maxproc 100 --reqsub $od/work_sh/trans_target.sh`;
 }else{
 	`sh $config{qsub_sh} --independent --maxproc 100 --reqsub $od/work_sh/trans_target.sh"`;
 }


`cat $od/Mid_Dir/*.out >$od/original.out`;

my %T;
$/="\n\n";
open (IN,"$od/original.out") or die $!;
while (<IN>) {
	chomp;
	my @A=split/\n/,$_;
	my $target=shift @A;
	$target=~s/^>//;
	my $target_chr=(split / /,$T_N{$target})[-2];
	my $lncRNA=shift @A;
	$lncRNA=~s/^>//;
	my @B=split/;/,$lncRNA;
	my $lncRNA_chr=($#B==0)?"THISISNOTLIMIT":$B[-2];
	next if $target_chr eq $lncRNA_chr;
	my $value=0;
	foreach my $line (@A) {
		$line=~/\((.*)\)/;
		$value=1 if $1<-20;
	}
	$T{$lncRNA}{$T_N{$target}}=1 if $value==1;
}

close IN;
$/="\n";

open (OUT,">$od/lncRNA_target.xls") or die $!;
print OUT "#lncRNA\tGenes\n";
foreach my $lncRNA (keys %T) {
	my $out_line="$lncRNA";
	my $limit=keys %{$T{$lncRNA}};
	if ($limit==0) {
		print OUT "$out_line\n";
	}
	else{
		$out_line.="\t";
		foreach my $target (keys %{$T{$lncRNA}}) {
			$out_line.="$target;";
		}
		$out_line=~s/;$//;
		$out_line.="\n";
		print OUT $out_line;
	}
}
close OUT;

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
	#求列表中的最大值
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
	#求列表中的最小值
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
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
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
Program Date:   2013.7.12
Usage:
  Options:

  -q  <file>  lncRNA file,fasta format,forced 
  
  -t  <file>  gene file,fasta format,forced 
  
  -c  <num>   Set the extension penalty in [dacal/mol]. default 30
  
  -od <dir>   output dir,forced
 
  -h         Help

USAGE
	print $usage;
	exit;
}
