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
my ($lnc,$gene,$dist,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"lnc:s"=>\$lnc,
				"gene:s"=>\$gene,
				"dist:s"=>\$dist,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($lnc and $gene and $od);
my $BedTools = "/share/nas2/genome/biosoft/bedtools/current/bin/bedtools";
mkdir $od unless -d $od;
my %targets;
$dist ||=100000;
my $lnc_temp = "$od/lnc.temp";
my $gene_temp = "$od/gene.temp";
my $result_temp = "$od/result.temp";
my $result = "$od/Cis_target_result.txt";

&GFFType($lnc,'mRNA', $lnc_temp);
&GFFType($gene,'gene', $gene_temp);

my $cmd = " $BedTools window  -w $dist -a $lnc_temp  -b $gene_temp >$result_temp" ;
&cmd_call($cmd);

############################### result Format ################################
open TEMP,"<$result_temp" or die $!;



while (<TEMP>) {
	next if (/^$/ or /^\#/) ;
	chomp;
	my @each = split /\t/,$_;
	my $lnc_q = (split /=/,(split /;/,$each[8])[0])[1];
	my $gene_t = (split /=/,$each[-1])[-1] ;
	push @{$targets{$lnc_q}}, $gene_t;
#	print OUT "$lnc_q\t$gene_t\n";

}
close TEMP;
open OUT, ">$result" or die $!;
print  OUT "#LncRNA\tGenes\n";
foreach my $l (keys %targets) {
	my $genes= join(';',@{$targets{$l}});
	print OUT "$l\t$genes\n";
}

close OUT;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}

sub GFFType{
	my $gff = shift;
	my $type = shift;
	my $out= shift;
	open GFF,"<$gff" or die $!;
	open OUT,">$out" or die $!;
	while (<GFF>) {
		next if (/^$/ or /^\#/) ;
		my $line=$_;
		my @each = split /\t/,$line;
		if ($each[2] eq $type) {
			print OUT $line;

		}
	}
	close GFF;
	close OUT;
}


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
  -lnc  <file>  lncRNA GFF file ,forced 
  
  -gene  <file>  gene GFF file ,forced 
  -dist			Distance between lnc and gene   
  -od <dir>   output dir,forced

  -h         Help

USAGE
	print $usage;
	exit;
}
