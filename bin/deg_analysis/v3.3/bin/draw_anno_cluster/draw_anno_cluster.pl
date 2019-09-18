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
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"k:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $od);
mkdir $od unless -d $od;
my $sh_dir="$od/work_sh";
mkdir $sh_dir unless -d $sh_dir;

$od=&ABSOLUTE_DIR($od);
$deg=&ABSOLUTE_DIR($deg);
$sh_dir=&ABSOLUTE_DIR($sh_dir);

my ($mid_file_tree,$mid_file_pathway);
runOrDie("ln -sf $fIn/*GO_tree.stat.xls $od");
$mid_file_tree=(glob("$od/*GO_tree.stat.xls"))[0];
die"Err:$mid_file_tree does not exis\n" unless -e $mid_file_tree;

&runOrDie("ln -sf $fIn/*Kegg.pathway $od/") ;
$mid_file_pathway=(glob("$od/*Kegg.pathway"))[0];
die "Err:file *.Kegg.ko or *.Kegg.pathway does not exist\n" if(!-e $mid_file_pathway);



open(OUT,">$sh_dir/draw_anno_cluster.sh") or die $!;
if (-e $mid_file_pathway) {
	print OUT "perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k kegg  -od $od/KEGG_Cluster \n";
#    `perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k kegg  -od $od/KEGG_Cluster`;
	system "perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k kegg  -od $od/KEGG_Cluster";
	`rm $mid_file_pathway`;
}
if (-e $mid_file_tree) {
	print OUT "perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k go	-od $od/GO_Cluster \n";
    #`perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k go	-od $od/GO_Cluster`;
	system "perl $Bin/kegggo_clustre.pl -i $od -deg $deg -k go -od $od/GO_Cluster";
	`rm $mid_file_tree`;
}

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
Program Date:   2013.10.17
Usage:
  Options:
  -i     <dir>  input dir,Allgene_Anno/Result forced 
  
  -deg   <file>  All_DEG/All.DEG_final.xls ,forced 
  
  -od    <file>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
