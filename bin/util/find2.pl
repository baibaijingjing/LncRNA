use strict;
use warnings;
use Cwd qw(abs_path);     #获取工作路径，即当前目标所在的路径（函数括号里的东西）
use Getopt::Long;         #获取选项
use Data::Dumper;          #可打印引用的东西。eg：print Dumper(\%hash \@array);
use FindBin qw($Bin $Script);  #$Bin  调用脚本的bin目录的路径，$Script  脚本名称  $RealBin 调用脚本的绝对路径  $RealScript  与脚本相关的脚本（用不着）
use File::Basename qw(basename dirname);  #basename函数获取文件名  dirname函数获取路径  fileparse函数获取扩展名
my $BEGIN_TIME=time();    #获取系统时间，可两者做差得到运行时间，单位为秒（s）
my $version="1.0.0";


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my ($cfg,%option_hash);  #选项的hash表
GetOptions(
  "cfg=s"    => \$cfg,   #=s表示串，eg：=i表示整数
  "help|?"   => \&USAGE,
  ) or &USAGE;
&USAGE unless (defined $cfg);  #判断选项的合法性
&log_current_time("$Script start……");


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#load input file,save the result
#print result
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
open READ_1,"$cfg"||die "can't open the file: $_\n";
while(<READ_1>){
	chomp;
	next if(/#/);
	s/\s+$//;
	my @arr=split /:/,$_,2;
	$option_hash{$arr[0]}=$arr[1];
}
my $output_dir=dirname($option_hash{output});
`mkdir $output_dir -p` if(!-d "$output_dir");
my $tmp_file_in=$option_hash{input};
$tmp_file_in="$tmp_file_in.xls";
my $tmp_file_out=$option_hash{output};
$tmp_file_out="$tmp_file_out.txt"; #生成输入输出路径
open READ_2,"$tmp_file_in" ||die "can't open the file: $_\n";
open WRITE_1,">$tmp_file_out"||die "can't create the file: $_\n";
my ($num_sample,$compare); #num_sample存入物种数目，compare存入要比对的列
my $tmp=<READ_2>;  #对第一行操作，获取要对第几列筛选，获取物种的数目
$tmp=~s/\#//g;  #去除第一列可能出现的#（可以省略）。
print WRITE_1"$tmp";#输出第一行，为后续做准备
my @tmp_array=split/\t/,$tmp; #切分第一行
#print "@tmp_array\n";
for my $k(0..$#tmp_array){
	if($option_hash{type}=~/$tmp_array[$k]/){
		$compare=$k;#将要比对的列存入变量中
	}
	if($tmp_array[$k]=~/FPKM/){
		$num_sample+=1;
	}
}
#print "$compare\n";
my @gene_list=split /\;/,$option_hash{list};
#print "@gene_list\n";
if(defined $compare){
	while(<READ_2>){
		my @tmp_all=split /\t/,$_;
		for my $j(@gene_list){
			if($tmp_all[$compare]=~/$j/){
				print WRITE_1"$_";
			}
		}
	}
}
else{
	print "sorry!no result create\n";
	print WRITE_1"sorry!no result create\n";
	exit;
}
close READ_1;
close READ_2;
close WRITE_1;


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#draw barplot pictrue by R
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my $barplot=<<__BARPLOT__;
library(gplots)
barpt<-read.table(file="$option_hash{output}.txt",header=T,sep="\\t")

barpt2=as.matrix(barpt[,c(seq(6,by=2,length.out=$num_sample))])
row.names(barpt2)=barpt[,1]


tiff(filename = "$output_dir/$option_hash{picture}.tiff",height = 1200,width = 1800,res = 200)
par(mar=c(5,4,4,15))
mp <- barplot2(barpt2, beside = TRUE,#legend.text=T,#xlab=colnames(barpt),
               col = rainbow(length(barpt[,1])),#c("skyblue", scales::alpha('red',.5)),
               #ylim = c(0, 0.1),
               main = "Barplot of Select Gene", xlab = "Sample",ylab = "FPKM",
               cex.names = 1, plot.ci = F,
               plot.grid = TRUE)
legend(x=par("usr")[2], y=par("usr")[4],legend = barpt[,1],fill=rainbow(length(barpt[,1])), box.col="white", xpd=TRUE)
box()
dev.off()

png(filename = "$output_dir/$option_hash{picture}.png",height = 1200,width = 1800,res = 200)
par(mar=c(5,4,4,15))
mp <- barplot2(barpt2, beside = TRUE,#legend.text=T,#xlab=colnames(barpt),
               col = rainbow(length(barpt[,1])),#c("skyblue", scales::alpha('red',.5)),
               #ylim = c(0, 0.1),
               main = "Barplot of Select Gene", xlab = "Sample",ylab = "FPKM",
               cex.names = 1, plot.ci = F,
               plot.grid = TRUE)
legend(x=par("usr")[2], y=par("usr")[4],legend = barpt[,1],fill=rainbow(length(barpt[,1])), box.col="white", xpd=TRUE)
box()
dev.off()

pdf(file = "$output_dir/$option_hash{picture}.pdf",height = 6,width = 9)
par(mar=c(5,4,4,15))
mp <- barplot2(barpt2, beside = TRUE,#legend.text=T,#xlab=colnames(barpt),
               col = rainbow(length(barpt[,1])),#c("skyblue", scales::alpha('red',.5)),
               #ylim = c(0, 0.1),
               main = "Barplot of Select Gene", xlab = "Sample",ylab = "FPKM",
               cex.names = 1, plot.ci = F,
               plot.grid = TRUE)
legend(x=par("usr")[2], y=par("usr")[4],legend = barpt[,1],fill=rainbow(length(barpt[,1])), box.col="white", xpd=TRUE)
box()
dev.off()
__BARPLOT__
open WRITE_2,">$output_dir/barplot.r"||die "can't create the file: $_\n";
print WRITE_2$barplot;
close WRITE_2;
system "Rscript $output_dir/barplot.r";


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#function
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################################################################################################################################################################################
 &log_current_time("$Script end……");    #调用时间函数
 my $run_time=time()-$BEGIN_TIME;
 print "the program run time is:$run_time.s\n";
 
  sub log_current_time {
     # get parameter    #获取参数
     my ($info) = @_;

     # get current time with string   #获取当前串的时间
     my $curr_time = &date_time_format(localtime(time()));   #格式化获取时间表示方法

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
  function: filter data and draw three picture

     Usage:
            --cfg      <FILE>  configuration files(contain:
                                type:transcript_id(Filter condition)
                                list:ENSRNOT0000005563;ENSRNOT0000000006;(support fuzzy search and precise query)
                                output:/share/nas1/huangxy/huangxy/jfsh/123(path and file,needn't expanded-name)
                                input:/share/nas1/huangxy/huangxy/tools/out(path and file,needn't expanded-name)
                                picture:456(file name,needn't expanded-name))
   Example:
            perl $Script --cfg  cfg.txt

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
