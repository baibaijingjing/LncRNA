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
my ($q,$t,$od,$inter);
GetOptions(
				"help|?" =>\&USAGE,
				"q:s"=>\$q,
				"t:s"=>\$t,
				"inter:i"=>\$inter,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($q and $t and $od);
$inter= $inter|| 100000;
mkdir $od unless -d $od;
my %Genes;
$/=">";
open (IN,$t) or die $!;
<IN>;
while (<IN>) {
	chomp;
	my @A=split/ /,(split/\n/,$_)[0];
	my ($start,$end)=split/-/,$A[-1];
	$Genes{$A[1]}{$start}{$end}=(split/\n/,$_)[0];
}
close IN;

open (IN,$q) or die $!;
open (OUT,">$od/lncRNA_target.xls") or die $!;
print OUT "#lncRNA\tGenes\n";
<IN>;
while (<IN>) {
	chomp;
	my @Target;
	my @A=split/;|\s/,(split/\n/,$_)[0];
	my ($start,$end)=split/-/,$A[-1];
	foreach my $key1 (keys %{$Genes{$A[-2]}}) {
		if ($key1>=$start-$inter) {
			foreach my $key2 (keys %{$Genes{$A[-2]}{$key1}}) {
				if ($key2<=$end+$inter) {
					if ($key2<$start || $key1>$end) {
						push @Target,$Genes{$A[-2]}{$key1}{$key2};
					}
				}
			}
		}
	}
	my $cis=(split/\n/,$_)[0];
	#my $c=(split/\t/,$cis)[0];
	print OUT "$cis\t".join(";",@Target)."\n" if @Target>0;
	print OUT "$cis\n" unless @Target>0;
}
close IN;
close OUT;
$/="\n";
my @DIR=qw(Cis_target );#Trans_target);
my %T;
my %G;
#foreach my $dir (@DIR) {
	open (IN,"$od/lncRNA_target.xls") or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		next unless /\t/;
		my ($lnc_all,$target)=split/\t/,$_,2;
		my $lnc=(split /\s+/,$lnc_all)[0];
		foreach my $tar (split/;/,$target) {
			my $t = (split /\s+/,$tar)[0];
			$T{$lnc}{$t}=1;
			$G{$tar}=1;
		}
	}
#}
close IN;

open (OUT,">$od/Cis_target_result.txt") or die $!;
print OUT "#lncRNA\tGenes\n";
foreach my $key1 (keys %T) {
        print OUT "$key1\t";
        #my $k1=(split /\s+/,$key1)[0];
        my @T;
        foreach my $key2 (keys %{$T{$key1}}) {
                push @T,$key2;
                #my $k2=(split /\s+/,$key2)[0];
        }
        print OUT join(";",@T)."\n";
}
close OUT;

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
Program Date:   2013.7.13
Usage:
  Options:
  -q  <file>  lncRNA file,fasta format>trans_name;u;chr 256-368,forced 
  
  -t  <file>  gene file,fasta format>gene_name chr 2560-2700,forced 
  
  -od <dir>   output dir,forced

  -h         Help

USAGE
	print $usage;
	exit;
}
