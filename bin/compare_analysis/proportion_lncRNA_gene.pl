use strict;
use warnings;
use Cwd qw(abs_path);     #»ñÈ¡¹¤×÷Â·¾¶£¬¼´µ±Ç°Ä¿±êËùÔÚµÄÂ·¾¶£¨º¯ÊýÀ¨ºÅÀïµÄ¶«Î÷£©
use Getopt::Long;         #»ñÈ¡Ñ¡Ïî
use Data::Dumper;          #¿É´òÓ¡ÒýÓÃµÄ¶«Î÷¡£eg£ºprint Dumper(\%hash \@array);
use FindBin qw($Bin $Script);  #$Bin  µ÷ÓÃ½Å±¾µÄbinÄ¿Â¼µÄÂ·¾¶£¬$Script  ½Å±¾Ãû³Æ  $RealBin µ÷ÓÃ½Å±¾µÄ¾ø¶ÔÂ·¾¶  $RealScript  Óë½Å±¾Ïà¹ØµÄ½Å±¾£¨ÓÃ²»×Å£©
use File::Basename qw(basename dirname);  #basenameº¯Êý»ñÈ¡ÎÄ¼þÃû  dirnameº¯Êý»ñÈ¡Â·¾¶  fileparseº¯Êý»ñÈ¡À©Õ¹Ãû
my $BEGIN_TIME=time();    #»ñÈ¡ÏµÍ³Ê±¼ä£¬¿ÉÁ½Õß×ö²îµÃµ½ÔËÐÐÊ±¼ä£¬µ¥Î»ÎªÃë£¨s£©
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my($input_file1,$input_file2,$output_dir); 
GetOptions(
  "in1=s"    => \$input_file1,  
  "in2=s"    => \$input_file2,
  "outdir=s" => \$output_dir,
  "help|?"   => \&USAGE,
  ) or &USAGE;
 &USAGE unless (defined $input_file1 and defined $input_file2 and defined $output_dir);  #ÅÐ¶ÏÑ¡ÏîµÄºÏ·¨ÐÔ
 &log_current_time("$Script start¡­¡­");
 `mkdir $output_dir -p` if(!-d "$output_dir");
 my $tmp_path=abs_path($output_dir);
 #print "\n\n$tmp_path\n\n";
 #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 #load input file,save the result
 #print result
 #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

open READ_1,"$input_file1"||die "can't open $!\n";#Éú³ÉÈ¾É«ÌåµÄ±àºÅÎÄ¼þ,Í³¼Æ¹ýÂËÇ°µÄlncRNAºÍgeneµÄÊýÄ¿
open WRITE_1,">$tmp_path/before_filter_all_count.out"||die "can't create $!\n";
open WRITE_11,">$tmp_path/before_filter_all_proportion.out"||die "can't create $!\n";
open WRITE_2,">$tmp_path/before_filter_lncRNA.out"||die "can't create $!\n";
open WRITE_3,">$tmp_path/before_filter_mRNA.out"||die "can't create $!\n";
my (%chrom_list,@chrom_list,%lncRNA_for_chrom_before,%gene_for_chrom_before,%count_lncRNA_before,%count_gene_before,$sum_lncRNA_before,$sum_gene);#%chrom_listÈ¾É«ÌåµÄhash£¬@chrom_listÈ¾É«ÌåµÄÁÐ±í£¬%lncRNA_for_chrom_before¹ýÂËÆ÷Ç°µÄlncRNAÖ¸ÏòÈ¾É«ÌåµÄhash
#%gene_for_chrom_before»ùÒòÖ¸ÏòÈ¾É«Ìå,%count_lncRNA,%count_gene·Ö±ð±íÊ¾¼ÆÊý
while(<READ_1>){
	$_=~/(.*?)\sCufflinks.*?gene_id \"(.*?)\"\; transcript_id \"(.*?)\"/;
	if(length($1)>6){
		$lncRNA_for_chrom_before{$3}="others";
		$gene_for_chrom_before{$2}="others";
	}
	else{
		$lncRNA_for_chrom_before{$3}=$1;
		$gene_for_chrom_before{$2}=$1;
	}
	my @array=split;
	next if(length($array[0])>6);
	$chrom_list{$array[0]}=1;
}
@chrom_list=keys %chrom_list;
push(@chrom_list,"others");   #ÒìÈ¾É«ÖÊ¹éÎªothers
@chrom_list=sort{$a<=>$b}@chrom_list;
#print "@chrom_list\n";      #Èç¹ûºóÐø»¹ÐèÒªµÄ»°¿ÉÒÔ½«ÆäÊä³öµ½Ò»¸öÎÄ±¾ÖÐ×÷ÎªÒ»¸öÐòÁÐµÄÐÅÏ¢
foreach(keys %lncRNA_for_chrom_before){
	$count_lncRNA_before{$lncRNA_for_chrom_before{$_}}+=1;
}
$sum_lncRNA_before=scalar(keys %lncRNA_for_chrom_before);
foreach(keys %gene_for_chrom_before){
	$count_gene_before{$gene_for_chrom_before{$_}}+=1;
}
$sum_gene=scalar(keys %gene_for_chrom_before);
print WRITE_1"chr	num	class\n";
print WRITE_11"chr	\%	class\n";
foreach(@chrom_list){
	my ($tmp1,$tmp2);
	if(exists $count_lncRNA_before{$_}){
		$tmp1=$count_lncRNA_before{$_};
	}
	else{
		$tmp1=0;
	}
	if(exists $count_gene_before{$_}){
		$tmp2=$count_gene_before{$_};
	}
	else{
		$tmp2=0;
	}
	my $tmp3=$tmp1/$sum_lncRNA_before;
	my $tmp4=$tmp2/$sum_gene;
	print WRITE_1"$_	$tmp1	lncRNA_before\n$_	$tmp2	mRNA\n";
	print WRITE_11"$_	$tmp3	lncRNA_before\n$_	$tmp4	mRNA\n";
	print WRITE_2"$_	$tmp1\n";
	print WRITE_3"$_	$tmp2\n";
}
close READ_1;
close WRITE_1;
close WRITE_11;
close WRITE_2;
close WRITE_3;

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# draw picture by R
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
`$config{Rscript} $Bin/simpleBar_one_var.r --infile $tmp_path/before_filter_lncRNA.out --outfile $tmp_path/before_filter_lncRNA --x.col 1 --y.col 2 --x.lab Chromosome --y.lab "Number of lncRNA" --number_label`;
`$config{Rscript} $Bin/simpleBar_one_var.r --infile $tmp_path/before_filter_mRNA.out --outfile $tmp_path/before_filter_mRNA --x.col 1 --y.col 2 --x.lab Chromosome --y.lab "Number of mRNA" --number_label`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/before_filter_all_count.out --outfile $tmp_path/before_filter_all_count --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Count`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/before_filter_all_proportion.out --outfile $tmp_path/before_filter_all_proportion --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Proportion`;


open READ_2,"$input_file2"||die "can't open $!\n";
open WRITE_4,">$tmp_path/filter_all_count.out"||die "can't create $!\n";
open WRITE_41,">$tmp_path/filter_all_proportion.out"||die "can't create $!\n";
open WRITE_5,">$tmp_path/filter_and_before_filter_lncRNA.out"||die "can't create $!\n";
open WRITE_51,">$tmp_path/filter_and_before_filter_lncRNA_proportion.out"||die "can't create $!\n";
open WRITE_6,">$tmp_path/filter_lncRNA_count.out"||die "can't create $!\n";
my (%lncRNA_for_chrom_filter,%gene_for_chrom_filter,%count_lncRNA_filter,%count_gene_filter,$sum_lncRNA_filter);
while(<READ_2>){
	$_=~/(.*?)\sCufflinks.*?gene_id \"(.*?)\"\; transcript_id \"(.*?)\"/;
	if(length($1)>6){
		$lncRNA_for_chrom_filter{$3}="others";
		$gene_for_chrom_filter{$2}="others";
	}
	else{
		$lncRNA_for_chrom_filter{$3}=$1;
		$gene_for_chrom_filter{$2}=$1;
	}
	my @array=split;
	next if(length($array[0])>6);
	$chrom_list{$array[0]}=1;
}
foreach(keys %lncRNA_for_chrom_filter){
	$count_lncRNA_filter{$lncRNA_for_chrom_filter{$_}}+=1;
}
$sum_lncRNA_filter=scalar(keys %lncRNA_for_chrom_filter);
print WRITE_4"chr	num	class\n";
print WRITE_5"chr	num	class\n";
print WRITE_41"chr	\%	class\n";
print WRITE_51"chr	\%	class\n";
foreach(@chrom_list){
	my ($tmp1,$tmp3,$tmp4);
	if(exists $count_lncRNA_filter{$_}){
		$tmp1=$count_lncRNA_filter{$_}
	}
	else{
		$tmp1=0;
	}
	if(exists $count_lncRNA_before{$_}){
		$tmp3=$count_lncRNA_before{$_}
	}
	else{
		$tmp3=0;
	}
	if(exists $count_gene_before{$_}){
		$tmp4=$count_gene_before{$_}
	}
	else{
		$tmp4=0;
	}
	my $tmp5=$tmp1/$sum_lncRNA_filter;
	my $tmp6=$tmp3/$sum_lncRNA_before;
	my $tmp7=$tmp4/$sum_gene;
	print WRITE_4"$_	$tmp1	lncRNA_filter	\n$_	$tmp4	mRNA\n";
	print WRITE_41"$_	$tmp5	lncRNA_filter	\n$_	$tmp7	mRNA\n";
	print WRITE_5"$_	$tmp3	lncRNA_before	\n$_	$tmp1	lncRNA_filter\n";
	print WRITE_51"$_	$tmp6	lncRNA_before	\n$_	$tmp5	lncRNA_filter\n";
	print WRITE_6"$_	$tmp1\n";
}
close READ_2;
close WRITE_4;
close WRITE_41;
close WRITE_5;
close WRITE_51;
close WRITE_6;

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# draw picture by R
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
`$config{Rscript} $Bin/simpleBar_one_var.r --infile $tmp_path/filter_lncRNA_count.out --outfile $tmp_path/filter_lncRNA_count --x.col 1 --y.col 2 --x.lab Chromosome --y.lab "Number of lncRNA" --number_label`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/filter_all_count.out --outfile $tmp_path/filter_all_count --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Count`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/filter_all_proportion.out --outfile $tmp_path/filter_all_proportion --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Proportion`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/filter_and_before_filter_lncRNA.out --outfile $tmp_path/filter_and_before_filter_lncRNA --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Count`;
`$config{Rscript} $Bin/simpleBar_two_var.r --infile $tmp_path/filter_and_before_filter_lncRNA_proportion.out --outfile $tmp_path/filter_and_before_filter_lncRNA_proportion --x.col 1 --y.col 2 --z.col 3 --x.lab Chromosome --y.lab Proportion`;


###################################################################################################################################################################################
 &log_current_time("$Script end¡­¡­");    #µ÷ÓÃÊ±¼äº¯Êý
 my $run_time=time()-$BEGIN_TIME;
 print "the program run time is:$run_time.s\n";
 
  sub log_current_time {
     # get parameter    #»ñÈ¡²ÎÊý
     my ($info) = @_;

     # get current time with string   #»ñÈ¡µ±Ç°´®µÄÊ±¼ä
     my $curr_time = &date_time_format(localtime(time()));   #¸ñÊ½»¯»ñÈ¡Ê±¼ä±íÊ¾·½·¨

     # print info with time
     print "[$curr_time] $info\n";
}

sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

###################################################################################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: huangxy <huangxy\@biomarker.com.cn> 
      Date: 2015-08-04
      note: When the length of the chromosome id is longer than six,we classified it as \"others\"

     Usage:
            --in1      <FILE>  before filter 
            --in2      <FILE>  filter
            --outdir   <DIR>   analysis output directory
   Example:
            perl $Script --in1 merged.gtf --in2 filter_final.gtf --outdir Analysis/

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
