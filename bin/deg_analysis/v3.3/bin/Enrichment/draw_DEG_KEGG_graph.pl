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
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")};
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
&USAGE unless ($fIn and $deg and $key and $od);
$fIn=&ABSOLUTE_DIR($fIn);
mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);
&runOrDie("ln -snf $fIn/*Kegg.ko $od/");
my $mid_file_ko=(glob("$od/*Kegg.ko"))[0];

&runOrDie("ln -snf $fIn/*Kegg.pathway $od/");
my $mid_file_pathway=(glob("$od/*Kegg.pathway"))[0];
die "Err:file *.Kegg.ko or *.Kegg.pathway does not exist\n" if(!-e $mid_file_ko or !-e $mid_file_pathway);

if (-f "$mid_file_ko" && -f "$mid_file_pathway") {
	`perl $Bin/KeggGo_enrich_map_web.v1.pl -d $deg -k $key -i $od -o $od -func kegg`;
	#`perl $Bin/KeggGo_enrich_map_web.pl -d $deg -k $key -i $en_dir -o $od -map /share/nas2/genome/pipeline_cfg/db_file.cfg`;
	`rm $od/*Kegg.ko`;
	`rm $od/*Kegg.pathway`;

}

system "perl $Bin/kegg_enrichment_plot.pl -enrich_file $od/pathway/kegg_enrichment/$key.KEGG.stat -od $od/Graph -key $key";

print "[".&GetTime."] perl $Bin/draw_KEGG_histogram.pl --ipf $od/pathway/kegg_enrichment/$key.KEGG.xls --opd $od/pathway/kegg_enrichment/ --prf $key.KEGG\n";
system "perl $Bin/draw_KEGG_histogram.pl --ipf $od/pathway/kegg_enrichment/$key.KEGG.xls --opd $od/pathway/kegg_enrichment/ --prf $key.KEGG";

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
  -i     <file>  input file,All_Database_annotation.xls,forced 
  
  -deg   <file>  deg file,forced 
  
  -k     <str>   keywords of output file,forced 
  
  -od    <file>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
