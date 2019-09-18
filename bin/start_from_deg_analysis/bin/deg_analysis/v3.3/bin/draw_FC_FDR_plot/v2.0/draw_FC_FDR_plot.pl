#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../../../../Config/lncRNA_pip.cfg")};
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($all,$de,$o,$thresh);
GetOptions(
				"help|?" =>\&USAGE,
				"all:s"=>\$all,
				"de:s"=>\$de,
				"o:s"=>\$o,
				"th:s"=>\$thresh,
				) or &USAGE;
&USAGE unless ($all and $de and $o);
&USAGE unless ($o=~/\.png$/);
$thresh||="2,0.01";
my ($fc,$p)=split/,/,$thresh;
$all=&ABSOLUTE_DIR($all);
$de=&ABSOLUTE_DIR($de);

###### set env variables
my $temp = `echo \$PATH`;
print "PATH=$temp\n";
$ENV{"PATH"} = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/:" . $ENV{"PATH"};
$temp = `echo \$PATH`;
print "PATH=$temp\n";

my %DE;
open (IN,$de) or die $!;
my %num;my $total;
while (<IN>) {
	next if /^\#/ or /^\s*$/;
	$_=~s/\s*$//;
	my ($name,$regulated)=(split/\s+/,$_)[0,-1];
	$DE{$name}=$regulated;
	$num{$regulated}++;
}
close IN;
my @TF;
chomp($total=`less -S $all|wc -l`);
$num{up}=0 if (!exists $num{up});
$num{down}=0 if (!exists $num{down});
$num{unchange}=$total-1-$num{up}-$num{down};
open (IN,$all) or die $!;
while (<IN>) {
	next if /^\#/;
	my $name1=(split/\s+/,$_)[0];
	my $type=(exists $DE{$name1})?"$DE{$name1}":"unchange";
	push @TF,"$type regulated:$num{$type}" if $type ne "unchange";
	push @TF,"$type:$num{unchange}" if $type eq "unchange";
}
close IN;
my $TF=join '","',@TF;

open (OUT,">$o.r") or die $!;
print OUT "#!$config{Rscript}\n";
print OUT <<END;
library(ggplot2)
ALL<- read.delim("$all", row.names = 1, header=TRUE,check.names =F)
Significant<-c("$TF")
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
#qplot( ALL[,c("log2FC")], -log10( ALL[,c("FDR")]), xlab="log2(FC)", ylab="-log10(FDR)", main="Volcano plot", size=I(0.5), colour=Significant)
df=data.frame(x=ALL[,c("log2FC")], y=-log10( ALL[,c("FDR")]),lab=factor(Significant,levels=unique(Significant)))
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
p=ggplot(data=df,mapping=aes(x=ALL[,c("log2FC")], y=-log10( ALL[,c("FDR")]),color=lab))
p=p+geom_point(size=1)
p=p+xlab("log2(FC)")+ylab("-log10(FDR)")+ggtitle("Volcano plot")+theme(plot.title=element_text(face="bold",size=14))
p=p+geom_hline(yintercept=-log10($p),linetype="longdash",size=0.2,colour="blue")+geom_vline(colour="blue",size=0.2,xintercept=c(log2($fc),log2(1/$fc)),linetype="longdash")
p=p+theme_classic()
p=p+scale_color_manual(name="Significant",values=c("red","green","black"),limits=c("up regulated:$num{up}","down regulated:$num{down}","unchange:$num{unchange}"))
print(p)

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
  -all        <file>   All genes expression list of Group,forced 
  
  -de         <file>   different genes expression list of Group,forced 
  
  -o          <file>   output file,*.png,forced 
  -th 		  <str>    fc and FDR thresh value ,separated by comma
  -h           Help

USAGE
	print $usage;
	exit;
}
