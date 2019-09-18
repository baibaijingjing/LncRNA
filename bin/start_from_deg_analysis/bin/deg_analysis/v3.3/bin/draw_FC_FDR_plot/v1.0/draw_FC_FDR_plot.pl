#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($all,$de,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"all:s"=>\$all,
				"de:s"=>\$de,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($all and $de and $o);
&USAGE unless ($o=~/\.png$/);

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
while (<IN>) {
	next if /^\#/;
	my $name=(split/\s+/,$_)[0];
	$DE{$name}=1;
}
close IN;
my @TF;
open (IN,$all) or die $!;
while (<IN>) {
	next if /^\#/;
	my $name1=(split/\s+/,$_)[0];
	my $type=(exists $DE{$name1})?"TRUE":"FALSE";
	push @TF,$type;
}
close IN;
my $TF=join '","',@TF;

open (OUT,">$o.r") or die $!;
print OUT "#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript\n";
print OUT <<END;
library(ggplot2)
ALL<- read.delim("$all", row.names = 1, header=TRUE)
Significant<-c("$TF")
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
qplot( ALL[,c("log2FC")], -log10( ALL[,c("FDR")]), xlab="log2(FC)", ylab="-log10(FDR)", main="Volcano plot", size=I(0.5), colour=Significant)
dev.off()

END

my $cmd;
#$cmd = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript $o.r";
$cmd = "/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript $o.r";
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
 
  -h           Help

USAGE
	print $usage;
	exit;
}
