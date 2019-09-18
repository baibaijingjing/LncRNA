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
&log_current_time("$Script start����");

open READ_1,"$infile1"||die "can't open the file: $!\n";
open READ_2,"$infile2"||die "can't open the file: $!\n";
open WRITE,">$outfile"||die"can't open the file :$_\n";
`sort -nr -k 2 $infile2>READ_3`;

#########���������############
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
&log_current_time("$Script end����");    #����ʱ�亯��
my $run_time=time()-$BEGIN_TIME;
print "$Script run time :$run_time\.s\n";

#----------------------------------------------------------------------------------
#function 
#----------------------------------------------------------------------------------
sub log_current_time {     #��ȡʱ��ĺ�������ʽ�����
	my ($info) = @_;    #��ȡ������һ��Ϊ    XXX����ʼִ�С���������
	my $curr_time = &date_time_format(localtime(time()));   #��ʽ����ȡʱ���ʾ����
	print "[$curr_time] $info\n";    #�����ӡ
}
##########################################################################################################################
sub date_time_format {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);  #��ȡ��ʽ����ʱ����Ϣ������ӡ��sprintf���÷���
}
########################################################################################################################## 
sub USAGE{   #ѡ���������
	my $usage=<<"__USAGE__";#����һ�п�ʼ��֪������__USAGE__Ϊֹ�����еĶ�����ͳһ�ĸ�ʽ��������У�����ע�ͣ���#__USAGE__Ϊ�������ţ�������eof���������Ҳ�������룬����ֻ����Ϊ������־
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
   print $usage;  #��ӡ������Ϣ
   exit;      #�����ø��ӳ���ʱ��ʾ�����޷�����ִ�У�ʹ��exit����ǿ���˳��ó���
}