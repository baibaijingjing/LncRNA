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
my ($fIn,$file,$o,$type);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"s:s"=>\$file,
				"o:s"=>\$o,
				"type:s"=>\$type,
				) or &USAGE;
&USAGE unless ($fIn and $file and $o);
mkdir "$o" unless (-d "$o");
my %H;
$file=&ABSOLUTE_DIR($file);
if (-f $file) {
    open(IN,$file) or die $!;
	while (<IN>) {
        chomp;
		next if (/^#/ || /^$/);
		my ($lnc,$tar)=(split /\s+/)[0,1];
		my @ta=split /;/,$tar;
		foreach my $gene (@ta){
			$H{$gene}=1;
		}
    }
    close IN;
}

if (-d $file ) {
	if ($type eq "cis") {
        my @file=(glob "$file/*_vs_*/Cis_Anno_enrichment/*lncRNA_DEG_targene.xls");
		
		foreach my $file (@file) {
			open (IN,$file) or die $!;
			while (<IN>) {
				next if /^\#/;
				my $name=(split/\t/,$_)[0];
				$H{$name}=1;
			}
			close IN;
		}
    }
    if ($type eq "trans") {
        my @file=(glob "$file/*_vs_*/Trans_Anno_enrichment/*lncRNA_DEG_targene.xls");
		foreach my $file (@file) {
			open (IN,$file) or die $!;
			while (<IN>) {
				next if /^\#/;
				my $name=(split/\t/,$_)[0];
				$H{$name}=1;
			}
			close IN;
		}
    }
    
}



open (IN,$fIn) or die $!;
open (OUT,">$o/All_lnc_target.list") or die $!;
while (<IN>) {
	print OUT $_ if /^\#/;
	my $name=(split/\t/,$_)[0];
	print OUT $_ if exists $H{$name};
}

close (IN) ;
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
Usage:
  Options:
  -i   <file>   input file,All_gene_fpkm.list,forced 
  
  -s   <files>  lncRNA target file or lnc_Diff_analysis ,forced 
  
  -o   <file>   output file all_lnc_target.list,forced
  -type	<str>	cis or trans
  
  -h         Help

USAGE
	print $usage;
	exit;
}
