#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#use Algorithm::Combinatorics qw(combinations permutations);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($infile1,$infile2,$outfile);
GetOptions(
				"help|?" =>\&USAGE,
				"in1:s"=>\$infile1,
				"in2:s"=>\$infile2,
				"out:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile1 and $infile2);
&log_current_time("$Script start……");

open READ_1,"$infile1"||die "can't open the file: $!\n";
open READ_2,"$infile2"||die "can't open the file: $!\n";
open WRITE,">$outfile"||die"can't open the file :$_\n";
`sort -nr -k 2 $infile2>READ_3`;

#########先输出行名############
my (@array,@err,$final);
my $link=(<READ_1>);chomp $link;
print WRITE "$link\t";
my @links=split /\t/,<READ_2>;chomp @links;
print WRITE "$links[15]\t$links[16]\t$links[17]\t$links[18]\t$links[19]\t$links[20]\n";

open READ_3,"READ_3";
	while(<READ_3>) {
		@err=split /\t/,$_;
		if (/$array[0]/) {
			push @array,($err[15],$err[16],$err[17],$err[18],$err[19],$err[20]);
			$final=join("\t",@array);chomp $final;
#			print "$final\n";
			print WRITE "$final\n";
			last;
		}
	}

while(<READ_1>) {
	@array=split /\t/,$_;chomp(@array);
	open READ_3,"READ_3";
	while(<READ_3>) {
		@err=split /\t/,$_;
		if (/$array[0]/) {
			push @array,($err[15],$err[16],$err[17],$err[18],$err[19],$err[20]);
			$final=join("\t",@array);chomp $final;
#			print "$final\n";
			print WRITE "$final\n";
			last;
		}
	}
}

close READ_1;
close READ_2;
close READ_3;
close WRITE;

#----------------------------------------------------------------------------------
#ending of work
#----------------------------------------------------------------------------------
&log_current_time("$Script end……");    #调用时间函数
my $run_time=time()-$BEGIN_TIME;
print "$Script run time :$run_time\.s\n";

#----------------------------------------------------------------------------------
#function 
#----------------------------------------------------------------------------------
sub log_current_time {     #获取时间的函数，格式化输出
	my ($info) = @_;    #获取参数（一般为    XXX程序开始执行…………）
	my $curr_time = &date_time_format(localtime(time()));   #格式化获取时间表示方法
	print "[$curr_time] $info\n";    #输出打印
}
##########################################################################################################################
sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);  #获取格式化的时间信息，不打印（sprintf的用法）
}
########################################################################################################################## 
sub USAGE{   #选项帮助函数
	my $usage=<<"__USAGE__";#从下一行开始，知道碰到__USAGE__为止，所有的东西按统一的格式存入变量中（包括注释），#__USAGE__为结束符号，类似用eof，不会输出也不会输入，仅仅只是作为结束标志
#-----------------------------------------------------------
 Program:$Script
 Version:$version
 Contact:<luml\@biomarker.com.cn>
 Data:2015-09-08
Function:assembly sequence
   USAGE:
         --in1       <STR>   input file1  [Must]
         --in2       <STR>   input file2  [Must]
         --out       <STR>   output file  [Must]
         --help          show the docment and exit
 Example:
    perl $Script --in1 All.DEG_final.xls --in2 final_all_out.xls --out final_DEG_annotation.xls
#---------------------------------------------------------
__USAGE__
   print $usage;  #打印帮助信息
   exit;      #当调用该子程序时表示程序无法正常执行，使用exit函数强行退出该程序
}