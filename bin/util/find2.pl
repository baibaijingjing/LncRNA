use strict;
use warnings;
use Cwd qw(abs_path);     #��ȡ����·��������ǰĿ�����ڵ�·��������������Ķ�����
use Getopt::Long;         #��ȡѡ��
use Data::Dumper;          #�ɴ�ӡ���õĶ�����eg��print Dumper(\%hash \@array);
use FindBin qw($Bin $Script);  #$Bin  ���ýű���binĿ¼��·����$Script  �ű�����  $RealBin ���ýű��ľ���·��  $RealScript  ��ű���صĽű����ò��ţ�
use File::Basename qw(basename dirname);  #basename������ȡ�ļ���  dirname������ȡ·��  fileparse������ȡ��չ��
my $BEGIN_TIME=time();    #��ȡϵͳʱ�䣬����������õ�����ʱ�䣬��λΪ�루s��
my $version="1.0.0";


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my ($cfg,%option_hash);  #ѡ���hash��
GetOptions(
  "cfg=s"    => \$cfg,   #=s��ʾ����eg��=i��ʾ����
  "help|?"   => \&USAGE,
  ) or &USAGE;
&USAGE unless (defined $cfg);  #�ж�ѡ��ĺϷ���
&log_current_time("$Script start����");


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
$tmp_file_out="$tmp_file_out.txt"; #�����������·��
open READ_2,"$tmp_file_in" ||die "can't open the file: $_\n";
open WRITE_1,">$tmp_file_out"||die "can't create the file: $_\n";
my ($num_sample,$compare); #num_sample����������Ŀ��compare����Ҫ�ȶԵ���
my $tmp=<READ_2>;  #�Ե�һ�в�������ȡҪ�Եڼ���ɸѡ����ȡ���ֵ���Ŀ
$tmp=~s/\#//g;  #ȥ����һ�п��ܳ��ֵ�#������ʡ�ԣ���
print WRITE_1"$tmp";#�����һ�У�Ϊ������׼��
my @tmp_array=split/\t/,$tmp; #�зֵ�һ��
#print "@tmp_array\n";
for my $k(0..$#tmp_array){
	if($option_hash{type}=~/$tmp_array[$k]/){
		$compare=$k;#��Ҫ�ȶԵ��д��������
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
 &log_current_time("$Script end����");    #����ʱ�亯��
 my $run_time=time()-$BEGIN_TIME;
 print "the program run time is:$run_time.s\n";
 
  sub log_current_time {
     # get parameter    #��ȡ����
     my ($info) = @_;

     # get current time with string   #��ȡ��ǰ����ʱ��
     my $curr_time = &date_time_format(localtime(time()));   #��ʽ����ȡʱ���ʾ����

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
