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
my ($id,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $od);


my %DEG;
my @DE_DIR=glob "$id/*_vs_*";
my $de_limit=@DE_DIR;

# check group number
if( ($de_limit>=2) && ($de_limit<=5) ) { 
	my %SX;
	# stat
	foreach my $dir (@DE_DIR) {
		my $name=basename $dir;
		$DEG{1}{$name}=$name;
		open (IN,"$dir/$name.DEG_final.xls") or die $!;
		while (<IN>) {
			next if $.==1;
			$DEG{$.}{$name}=(split/\s/,$_)[0];
			$SX{$name}=$.;
		}
		close IN;
	}

	# output
	open(OUT,">$od/All_DEG_veen.genes") or die $!;
	foreach my $num (sort {$a<=>$b} keys %DEG) {
		my $line;
		foreach my $name (sort {$SX{$b}<=>$SX{$a}} keys %SX) {
			$line.="$DEG{$num}{$name}\t" if exists $DEG{$num}{$name};
			$line.="\t" unless exists $DEG{$num}{$name};
		}
		$line=~s/\t$/\n/;
		print OUT $line;
	}
	close OUT;

	# drawing
	`perl $Bin/veen_draw.pl -f $od/All_DEG_veen.genes -o $od/All_DEG_veen.svg -font1 10 -font 8`;
}else{
	print "\nWaring:cluster number for venn plot should be 2<= n <=5 \n";
}

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
Program Date:   2013.10.10
Usage:
  Options:
  -id   <dir>  input dir,DEG_Analysis,forced 
  
  -od   <dir>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
