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
my $Title="LncRNA";   
#my %config=%{selectconf("/share/nas1/niulg/pipline/v2.2.5/Config")};
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$lncdeg,$od,$cfg);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"lncdeg:s"=>\$lncdeg,
				"od:s"=>\$od,
				"cfg:s"=>\$cfg,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $lncdeg and $cfg and $od);
$fIn=&ABSOLUTE_DIR($fIn);
$deg=&ABSOLUTE_DIR($deg);
$lncdeg=&ABSOLUTE_DIR($lncdeg);
mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);
mkdir "$od/work_sh" unless -d "$od/work_sh";
#####################read cfg file
my %para=%{readconf("$cfg")};
my ($FC,$FDR);
if (exists $para{'fold'}) {
	$FC=$para{'fold'};
	if ($FC<=0) {
		print "Fold Change Value Should Greater than zero!";die;
	}
}
if (!exists $para{'fold'}){
	$FC=2;
}

if (exists $para{'FDR'}) {
	$FDR=$para{'FDR'};
}else{
	$FDR=0.01;
}
my ($all_gene,$all_lnc,$genmoe_file);
$genmoe_file="$fIn/Tophat_Cufflinks/Ref_Genome/genome_size.txt" if (-f "$fIn/Tophat_Cufflinks/Ref_Genome/genome_size.txt");
$all_gene="$fIn/geneExpression/AllSample.genes_expression.xls" if (-f "$fIn/geneExpression/AllSample.genes_expression.xls");
$all_lnc="$fIn/LncExpression/AllSample.isoforms_expression.xls" if (-f "$fIn/LncExpression/AllSample.isoforms_expression.xls");
my @DIR=glob("$deg/*vs*");
open(OUT,">$od/work_sh/draw.sh") or die $!;
foreach my $dir (@DIR){
	my $group=basename($dir);
	print "$group\n";
	my $gene_deg="$dir/$group.DEG_final.xls";
	my $lnc_deg="$lncdeg/$group/$group.DEG_final.xls";
	my $gene_final="$dir/$group.final.xls";
	my $lnc_final="$lncdeg/$group/$group.final.xls";
	my ($control,$treated)=split/_vs_/,$group;
	print OUT "perl $Bin/draw_FC_FDR_volcano.pl -all $gene_final $lnc_final -de $gene_deg -lnc $lnc_deg -th $FC,$FDR -o $od/$group.FC_FDR.png \n";
	print OUT "perl $Bin/draw_FC_count_MA.pl -treated $treated -control $control -all $gene_final $lnc_final -de $gene_deg -lnc $lnc_deg -o $od/$group.FC_count.png \n";
	print OUT "perl $Bin/compare_deg_circos.pl --chr $genmoe_file --od $od --deg $gene_deg --all $all_gene --lncdeg $lnc_deg --lncall $all_lnc --outfile $group.circos \n";
}

close OUT;
system "sh $od/work_sh/draw.sh";
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
Program Date:   2013.10.17
Usage:
  Options:
  -i     <dir>  input dir Basic_Analysis,forced 
  
  -deg   <dir>  DEG_Analysis,forced 
  
  -lncdeg    <dir>   Lnc_Diff_Analysis,forced 
  
  -od    <file>  output dir,forced 
  -cfg	<file> 		detail.cfg ,forced
  -h         Help

USAGE
	print $usage;
	exit;
}
