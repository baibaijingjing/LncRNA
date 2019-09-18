#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};
my ($ingtf,$odir,$key);
GetOptions(
	"help|?"=>\&USAGE,
	"i:s"   =>\$ingtf,
	"o:s"   =>\$odir,
	"k:s"   =>\$key,
	)or &USAGE;
&USAGE unless ($ingtf and $odir and $key);
#system "mkdir -p $odir" unless (-d $odir);
$odir = abs_path($odir);
#main 
#contact :renhd\@biiomarker.com.cn
#if (@ARGV != 2) {
#	print "\n\tFunction: Extract len and exon_number from filter_final.gtf & merged.gtf\n\n";
#	print "\tInfiles: merge.gtf  filter_final.gtf \n\n";
#	print "\tperl mRNA_lnc_stat.pl <gtf>  <out_dir>\n\n";
#	exit;
#}

my %info;
open IN,"$ingtf" || die $!;
open OUT ,">$odir/$key.stat.xls";
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/);
	#my $line = $_;
	my @tmp=split /\t+/,$_;
	my ($trans_id) = $_ =~/transcript_id "(\S+)";/;
	#$tmp[8]=~/transcript_id\s\"([^\"]+)\";\sexon_number\s\"(\d+)\";/;
	#($trans_id,$info{$trans_id}{'exon_number'})=($1,$2);
	$info{$trans_id}{'len'}=$tmp[4]-$tmp[3]+1;
	$info{$trans_id}{'exon_number'}=0;
	($info{$trans_id}{'exon_number'}) = $_ =~/exon_number\s\"(\d+)\";/;
#	$info{$trans_id}{exon_number}= $line =~/exon_number "(\d+)";/;
	if (exists $info{$trans_id}){
		$info{$trans_id}{'length'}+=$info{$trans_id}{'len'};
		$info{$trans_id}{'exon_number'}=&max($info{$trans_id}{'exon_number'},0);
	}else{
		$info{$trans_id}{'length'}=$info{$trans_id}{'len'};
		#	$info{$trans_id}{'exon_number'}=$exon;
	}
	print OUT "$trans_id\t$info{$trans_id}{'exon_number'}\t$info{$trans_id}{'length'}\n";
}
close IN;
close OUT;
my %trans;
my %tra;
open IN,"$odir/$key.stat.xls";
open OUT,">$odir/$key.len.xls";
open OUT1,">$odir/$key.exon.xls";
while (<IN>){
	chomp;
	my @line=split /\t/,$_;
	if  (exists $trans{$line[0]}){
#		$trans{$line[0]}{len}=&max($line[2],$trans{$line[0]}{len});
		$tra{$line[0]}=&max($line[1],$tra{$line[0]});
		$trans{$line[0]}=&max($line[2],$trans{$line[0]});
	}else{
		$tra{$line[0]}=$line[1];
#		$trans{$line[0]}{len}=$line[2];
		$trans{$line[0]}=$line[2];
	}
#	print Dumper \%trans;

}
for my $k1 (sort keys %trans){
#	my %array = %{$trans{$k1}}; 
		print OUT "$k1\t$trans{$k1}\n";
		print OUT1 "$k1\t$tra{$k1}\n";
	}


close IN;
close OUT;
`rm $odir/$key.stat.xls`;
#`less  $odir/$key.len.xls |awk '{print \$2}'|sort |uniq -c > $odir/$key.LEN`;
#`less $odir/$key.exon.xls |awk '{print \$2}'|sort |uniq -c > $odir/$key.EXON`;
#`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/simpleBar.r --infile $odir/$key.LEN --outfile $key.len --x.col 1 --y.col 2 --x.lab length --y.lab count `;
#`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/simpleBar.r --infile $odir/$key.EXON --outfile $key.exon --x.col 2 --y.col 1 --x.lab exon --y.lab count `;

my %len;
open IN ,"$odir/$key.len.xls";
$len{200}=$len{400}=$len{600}=$len{800}=$len{1000}=$len{1200}=$len{1400}=$len{1600}=$len{1800}=$len{2000}=$len{2200}=$len{2400}=$len{2600}=$len{2800}=$len{3000}=$len{'>=3000'}=0;
while (<IN>){
	my @tmp= (split /\s+/,$_);
	if ($tmp[1]>=0 && $tmp[1]<200){$len{200}++;}	if ($tmp[1]>=200 && $tmp[1]<400){$len{400}++;}
	if ($tmp[1]>=400 && $tmp[1]<600){$len{600}++;}	if ($tmp[1]>=600 && $tmp[1]<800){$len{800}++;}
	if ($tmp[1]>=800 && $tmp[1]<1000){$len{1000}++;}	if ($tmp[1]>=1000 && $tmp[1]<1200){$len{1200}++;}
	if ($tmp[1]>=1200 && $tmp[1]<1400){$len{1400}++;}	if ($tmp[1]>=1400 && $tmp[1]<1600){$len{1600}++;}
	if ($tmp[1]>=1600 && $tmp[1]<1800){$len{1800}++;}	if ($tmp[1]>=1800 && $tmp[1]<2000){$len{2000}++;}
	if ($tmp[1]>=2000 && $tmp[1]<2200){$len{2200}++;}	if ($tmp[1]>=2200 && $tmp[1]<2400){$len{2400}++;}
	if ($tmp[1]>=2400 && $tmp[1]<2600){$len{2600}++;}	if ($tmp[1]>=2600 && $tmp[1]<2800){$len{2800}++;}
	if ($tmp[1]>=2800 && $tmp[1]<3000){$len{3000}++;}	if ($tmp[1]>=3000){$len{'>=3000'}++;}

}
my $all_len;
close IN;
open OUT ,">$odir/$key.stat.len";
print OUT "200\t$len{200}\n400\t$len{400}\n600\t$len{600}\n800\t$len{800}\n1000\t$len{1000}\n1200\t$len{1200}\n1400\t$len{1400}\n1600\t$len{1600}\n1800\t$len{1800}\n2000\t$len{2000}\n2200\t$len{2200}\n2400\t$len{2400}\n2600\t$len{2600}\n2800\t$len{2800}\n3000\t$len{3000}\n>=3000\t$len{'>=3000'}";
close OUT;
#plot
`$config{Rscript} $Bin/simpleBar.r --infile $odir/$key.stat.len --outfile $odir/$key.len.png  --x.col 1 --y.col 2 --x.lab $key.len --y.lab count --axis.size 10`;                              #输出的是原来的长度个数的柱状图             
#`mv $odir/$key.len $odir/$key.len.png`;
#########################################################修改的地方################################################################################################################################################
open OUT ,">$odir/$key.stat.len.proportation";
$all_len=$len{200}+$len{400}+$len{600}+$len{800}+$len{1000}+$len{1200}+$len{1400}+$len{1600}+$len{1800}+$len{2000}+$len{2200}+$len{2400}+$len{2600}+$len{2800}+$len{3000}+$len{'>=3000'};     #增加了一个总的长度
#计算每个lncRNA长度所占的比例大小，改变的原来的值
foreach $key(keys %len){
$len{$key}=$len{$key}/$all_len;
}
print OUT "200\t$len{200}\n400\t$len{400}\n600\t$len{600}\n800\t$len{800}\n1000\t$len{1000}\n1200\t$len{1200}\n1400\t$len{1400}\n1600\t$len{1600}\n1800\t$len{1800}\n2000\t$len{2000}\n2200\t$len{2200}\n2400\t$len{2400}\n2600\t$len{2600}\n2800\t$len{2800}\n3000\t$len{3000}\n>=3000\t$len{'>=3000'}";
close OUT;
`$config{Rscript} $Bin/simpleBar_proportation.r --infile $odir/$key.stat.len.proportation --outfile $odir/$key.len.proportion.png  --x.col 1 --y.col 2 --x.lab $key.len --y.lab Proportion --axis.size 10`;                              #将Y轴名称变为Proportion

#####################################################################################################################################################################################################################

my %exon;
my $i;
$exon{1}=$exon{2}=$exon{3}=$exon{4}=$exon{5}=$exon{6}=$exon{7}=$exon{8}=$exon{9}=$exon{10}=$exon{11}=$exon{12}=$exon{13}=$exon{14}=$exon{15}=$exon{16}=$exon{17}=$exon{18}=$exon{19}=$exon{20}=$exon{21}=$exon{22}=$exon{23}=$exon{24}=$exon{25}=$exon{26}=$exon{27}=$exon{28}=$exon{29}=$exon{30}=$i=0;
open IN,"$odir/$key.exon.xls";
while (<IN>){	
	my @tmp =(split /\s+/,$_);	#if (exists $exon{$tmp[1]}){$exon{$tmp[1]}++;}else{	$exon{$tmp[1]}=1}
	if ($tmp[1]==1 ){
		if ($key=~/lncRNA/){ next;}
		if ($key=~/mRNA/){$exon{1}++;}
	}
	if ($tmp[1]==2){ $exon{2}++;}  if ($tmp[1]==3){ $exon{3}++;}  if ($tmp[1]==4){ $exon{4}++;}
	if ($tmp[1]==5){ $exon{5}++;}  if ($tmp[1]==6){ $exon{6}++;}  if ($tmp[1]==7){ $exon{7}++;}  if ($tmp[1]==8){ $exon{8}++;}
	if ($tmp[1]==9){ $exon{9}++;}  if ($tmp[1]==10){ $exon{10}++;}    if ($tmp[1]==11){ $exon{11}++;}    if ($tmp[1]==12){ $exon{12}++;}
	if ($tmp[1]==13){ $exon{13}++;}   if ($tmp[1]==14){ $exon{14}++;}    if ($tmp[1]==15){ $exon{15}++;}    if ($tmp[1]==16){ $exon{16}++;}
	if ($tmp[1]==17){ $exon{17}++;}   if ($tmp[1]==18){ $exon{18}++;}    if ($tmp[1]==19){ $exon{19}++;}    if ($tmp[1]==20){ $exon{20}++;}
	if ($tmp[1]==21){ $exon{21}++;}   if ($tmp[1]==22){ $exon{22}++;}    if ($tmp[1]==23){ $exon{23}++;}    if ($tmp[1]==24){ $exon{24}++;}
	if ($tmp[1]==25){ $exon{25}++;}   if ($tmp[1]==26){ $exon{26}++;}    if ($tmp[1]==27){ $exon{27}++;}    if ($tmp[1]==28){ $exon{28}++;}
	if ($tmp[1]==29){ $exon{29}++;}		if ($tmp[1] >= 30){ $i++;  }
}
close IN;
open OUT ,">$odir/$key.stat.exon";
#$exon{1}=$exon{2}=$exon{3}=$exon{4}=$exon{5}=$exon{6}=$exon{7}=$exon{8}=$exon{9}=$exon{10}=$exon{11}=$exon{12}=$exon{13}=$exon{14}=$exon{15}=$exon{16}=$exon{17}=$exon{18}=$exon{19}=$exon{20}=$exon{21}=$exon{22}=$exon{23}=$exon{24}=$exon{25}=$exon{26}=$exon{27}=$exon{28}=$exon{29}=$exon{30}=$i=0;
#for my $key (sort keys %exon){
#	if ($key==1){ $exon{1}++;}	if ($key==2){ $exon{2}++;}	if ($key==3){ $exon{3}++;}	if ($key==4){ $exon{4}++;}
#	if ($key==5){ $exon{5}++;}	if ($key==6){ $exon{6}++;}	if ($key==7){ $exon{7}++;}	if ($key==8){ $exon{8}++;}
#	if ($key==9){ $exon{9}++;}	if ($key==10){ $exon{10}++;}	if ($key==11){ $exon{11}++;}	if ($key==12){ $exon{12}++;}
#	if ($key==13){ $exon{13}++;}	if ($key==14){ $exon{14}++;}	if ($key==15){ $exon{15}++;}	if ($key==16){ $exon{16}++;}
#	if ($key==17){ $exon{17}++;}	if ($key==18){ $exon{18}++;}	if ($key==19){ $exon{19}++;}	if ($key==20){ $exon{20}++;}
#	if ($key==21){ $exon{21}++;}	if ($key==22){ $exon{22}++;}	if ($key==23){ $exon{23}++;}	if ($key==24){ $exon{24}++;}
#	if ($key==25){ $exon{25}++;}	if ($key==26){ $exon{26}++;}	if ($key==27){ $exon{27}++;}	if ($key==28){ $exon{28}++;}
#	if ($key==29){ $exon{29}++;} if ($key==30){ $exon{30}++;}
#	if ($key >= 30){ $i++;	}
#} 
print OUT "1\t$exon{1}\n2\t$exon{2}\n3\t$exon{3}\n4\t$exon{4}\n5\t$exon{5}\n6\t$exon{6}\n7\t$exon{7}\n8\t$exon{8}\n9\t$exon{9}\n10\t$exon{10}\n11\t$exon{11}\n12\t$exon{12}\n13\t$exon{13}\n14\t$exon{14}\n15\t$exon{15}\n16\t$exon{16}\n17\t$exon{17}\n18\t$exon{18}\n19\t$exon{19}\n20\t$exon{20}\n21\t$exon{21}\n22\t$exon{22}\n23\t$exon{23}\n24\t$exon{24}\n25\t$exon{25}\n26\t$exon{26}\n27\t$exon{27}\n28\t$exon{28}\n29\t$exon{29}\n>=30\t$i\n";
close OUT;
`$config{Rscript} $Bin/simpleBar.r --infile $odir/$key.stat.exon --outfile $odir/$key.exon --x.col 1 --y.col 2 --x.lab $key.exon --y.lab count --no.grid`;
`mv $odir/$key.exon $odir/$key.exon.png`;
###########################################################################
sub max{
	my $a=shift;
	my $b=shift;
	my $tmp;
	if ($a >= $b){
		$tmp=$a;
	}else{
		$tmp=$b;
	}
	return $tmp;
}


sub USAGE {
	        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------
   Program: 	$Script
   Discription: Extract len and exon_number from filter_final.gtf & merged.gtf\n\n
   USAGE:
   		-i	input file gtf format
		-o	output dir directory
		-k input file (mRNA or lnc)
		-h	help documents
----------------------------------------------------------------------------------------
USAGE
			print $usage;
			exit ;
		}
