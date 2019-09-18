#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir,$odir);

GetOptions(
	"indir:s" => \$indir,
	"od:s" => \$odir,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($indir and $odir);
my $dir_abs = abs_path($indir);
mkdir $odir unless (-d $odir);

my %table_info;
my $base_num=0;
my $min_data=100;
my $min_Q30 = 1000;
my ($min_map_ratio,$max_map_ratio) = (100,0);


############ Get $sample_num, $base_num, $min_data, $min_Q30
my $sample = glob(&get_table('样品测序数据评估统计表')) ;
chomp (my $sample_num = `less $sample |wc -l`);
$sample_num = $sample_num -1;
open (IN,"$sample") or die $!;
while (<IN>) {
	chomp;
	next if (/^#/);
	my ($sampleid,$BMK_ID,$readsum,$basesum,$GC,$N,$Q30) = split /\t/,$_;
	$basesum = sprintf "%0.2f",$basesum/1000000000;			### convert to Gb
	$base_num += $basesum;
	$min_data = $basesum if ($basesum<$min_data);
	$min_Q30 = $Q30 if ($Q30<$min_Q30);
}
close IN;

############ Get $min_map_ratio, $max_map_ratio
my @mapped_datas = glob (&get_table('Clean Data与参考基因组比对结果统计表'));
foreach (my $i=0;$i<=$#mapped_datas;$i++){
	open (IN,$mapped_datas[$i]) or die $!;
	while (<IN>) {
		chomp;
		next if (/^\s+/);
		if (/^mapped Reads/){
			my ($id,$bases,$ratio) = split /\t+/,$_;
			$ratio =~ s/%$//;
			$min_map_ratio = $ratio if ($ratio < $min_map_ratio);
			$max_map_ratio = $ratio if ($ratio > $max_map_ratio);
		}		
	}
	close IN;
}
$min_map_ratio = $min_map_ratio.'%';
$max_map_ratio = $max_map_ratio.'%';

########### Get $new_gene_num
chomp (my $new_gene_num = `grep ">" $dir_abs/mRNA/NewGene/*.newGene.longest_transcript.fa |wc -l`)	;

########## Get $new_gene_ann_num
my $new_gene_ann_num;
my $newGene_anno = glob (&get_table('新基因功能注释结果统计'));
open (IN,$newGene_anno) or die $!;
while (<IN>) {
	chomp;
	next if (/^#/);
	if (/^All_Annotated/) {
		$new_gene_ann_num = (split /\s+/,$_)[1];
	}
}
close IN;

########## Get $deg_num, $lnc_num, $diff_lnc, $optimized_gene_num
chomp (my $deg_num = `less $dir_abs/mRNA/DEG/All_DEG/All.DEG_final.xls |wc -l`);
$deg_num = $deg_num-1;

chomp (my $diff_lnc = `less $dir_abs/LncRNA/DEG/All_DEG/All.DEG_final.xls |wc -l`);
$diff_lnc = $diff_lnc-1;

chomp (my $lnc_num = `less $dir_abs/LncRNA/Identify/LncRNA.id |wc -l`);
$lnc_num -= 1;

my $optimized_gene_path = glob("$dir_abs/mRNA/Gene_Structure_Optimize/*.xls");
#print "$optimized_gene_path\n";
my $optimized_gene_num = `less $optimized_gene_path |awk '{print \$1}' |uniq |wc -l`;
chomp ($optimized_gene_num);
$optimized_gene_num -= 1;

open (OUT, ">$odir/data.txt") or die $!;
print OUT "Sample_num\t$sample_num\n";
print OUT "base_num\t$base_num\n";
print OUT "min_data\t$min_data\n";
print OUT "min_Q30\t$min_Q30"."%\n";
print OUT "min_map_ratio\t$min_map_ratio\n";
print OUT "max_map_ratio\t$max_map_ratio\n";
print OUT "new_gene_num\t$new_gene_num\n";
print OUT "new_gene_ann_num\t$new_gene_ann_num\n";
print OUT "deg_num\t$deg_num\n";
print OUT "diff_lnc\t$diff_lnc\n";
print OUT "lnc_num\t$lnc_num\n";
print OUT "optimized_gene_num\t$optimized_gene_num\n";





sub get_table {
	chomp (my $id = shift);
	%table_info=(
		"样品测序数据评估统计表" => "$dir_abs/QC/Clean_Data/AllSample_GC_Q.stat",
		"Clean Data与参考基因组比对结果统计表" => "$dir_abs/mRNA/Expression/*.mappedStat.xls",
		"新基因功能注释结果统计" => "$dir_abs/mRNA/NewGene/Annotation/Function_Annotation.stat.xls",
	);
	return $table_info{$id};
}



sub data_cfg_read {
	my $cfg = shift;
	my $num=0;
	open (CFG,$cfg) or die $!;
	while (<CFG>) {
			chomp;
			next if (/^\s+/ or /^#/);
			if (/^Sample/) {
				$num++;
			}			
	}
	return $num;
	close CFG;
}




sub detail_cfg_read {
	my $cfg=shift;
	my %detail;
	
	open (CFG,$cfg) or die $!;
	while (<CFG>) {
		chomp;
		next if (/^\s+/ or /^#/);
		my ($key,$value) = split /\s+/,$_;
		if ($key eq 'Project_id') {
			$detail{$key} = $value;
		}
	}
	return %detail;
	close CFG;
}





#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei <yaob\@biomarker.com.cn> 
      Date: 

     Usage:
		-indir	<input directory>
		-od	<output directory>
		-h	help documents

   Example:
            perl $Script  -indir indir -od odir

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
