#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($t,$c,$all,$de,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"treated:s"=>\$t,
				"control:s"=>\$c,
				"all:s"=>\$all,
				"de:s"=>\$de,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($t and $c and $all and $de and $o);
&USAGE unless ($o=~/\.png$/);

$all=&ABSOLUTE_DIR($all);
$de=&ABSOLUTE_DIR($de);

###### set env variables
my $temp = `echo \$PATH`;
print "PATH=$temp\n";
$ENV{"PATH"} = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/:" . $ENV{"PATH"};
$temp = `echo \$PATH`;
print "PATH=$temp\n";

my $str=$t.'_'.$c;
$str=~s/_/","/g;
my %DE;
open (IN,$de) or die $!;
while (<IN>) {
	next if /^\#/;
	chomp;
	my $name=(split/\s+/,$_)[0];
	my $type=(split/\s+/,$_)[-1];
	$DE{$name}=$type;
}
close IN;
my @TF;
open (IN,$all) or die $!;
while (<IN>) {
	next if /^\#/;
	my $name1=(split/\s+/,$_)[0];
	my $type=(exists $DE{$name1})?"$DE{$name1}":"normal";
	push @TF,$type;
}
close IN;
my $TF=join '","',@TF;


open (OUT,">$o.r") or die $!;
print OUT "#!$config{Rscript}\n";
print OUT <<END;
library(ggplot2)
ALL<- read.delim("$all", row.names = 1, header=TRUE,check.names =F)
Significant<-c("$TF")


v1<-strsplit("$t","_")
v2<-strsplit("$c","_")
v1<-as.vector(v1[[1]])
v2<-as.vector(v2[[1]])
	if( length(v1)== 1 ){
		lev1 <- ALL[,v1]
	}else{
		lev1 <- rowMeans( ALL[,v1] )
	}
	if( length(v2) == 1 ){
		lev2 <- ALL[,v2]
	}else{
		lev2 <- rowMeans( ALL[,v2] )
	}

p<-qplot(log10( lev2 ), log10( lev1 ), xlab="$c", ylab="$t", main="Correlation Plot", size=I(0.5), colour=Significant)
p <- p + theme_classic()
p<- p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p <- p + geom_abline(slope=1, size=0.5, colour=6, alpha=0.5)
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
print(p)
dev.off()

END

my $cmd;
#$cmd = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript $o.r";
$cmd = "$config{Rscript} $o.r";
print "$cmd\n";
$temp = `$cmd`;





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
Program Date:   2013.10.14
Usage:
  Options:
  -treated    <str>    names of treated samples,"T1_T2",forced 
  
  -control    <str>    names of control samples,"T3_T4",forced 
  
  -all        <file>   All genes expression list of Group,forced 
  
  -de         <file>   different genes expression list of Group,forced 
  
  -o          <file>   output file,*.png,forced 
 
  -h           Help

USAGE
	print $usage;
	exit;
}
