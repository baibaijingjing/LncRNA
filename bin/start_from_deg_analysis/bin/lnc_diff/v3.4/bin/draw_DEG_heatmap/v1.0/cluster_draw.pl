#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$cluster,$k,$id,$fn,$sample_cluster,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
#				"cluster:s"=>\$cluster,
#				"k:s"=>\$k,
#				"id"=>\$id,
#				"e"=>\$sample_cluster,
				"fn:i"=>\$fn,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $od);
$fIn=&ABSOLUTE_DIR($fIn);
#$cluster||=1;
$fn||=0;

$od ||="./";


mkdir $od unless -d $od;
chdir $od;
$od=`pwd`;
chomp $od;
my $od_name=(split/\//,$od)[-2];
#(my $od_name=$od)=~s/[^\/]*\///g;		视情况而定

=pod
if ($cluster==2||$cluster==3) {
	&USAGE unless (defined($k));
	&MKDIR("$od/K_means");
	chdir "$od/K_means";
	`perl $Bin/bin/K_means.pl -i $fIn -k $k -abs -od ./ `;
}
=cut
#if ($cluster==1||$cluster==3) {
	mkdir "$od/hierarchical"; 
	chdir "$od/hierarchical";
	open (IN,$fIn) or die $!;
	open (OUT,">change0") or die $!;
	while (<IN>) {
		chomp;
		my @A=split/\t/,$_;
		if($fn)
		{
			splice(@A,-$fn);
		}
		if(/^#/)
		{
			print OUT (join "\t",@A)."\n";
			next;
		}
		my $name=shift @A;
		print OUT $name;
		foreach my $value (@A) {
			$value+=1;
			print OUT "\t$value";
		}
		print OUT "\n";
	}
	close IN;
	close OUT;


system "Rscript $Bin/cluster_draw.r --infile change0 --outfile $od/hierarchical/$od_name.DEG.cluster.png";

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


sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}


sub USAGE {#
	my $usage=<<"USAGE";
---------------------------------------------------------------------------------------------------------------
ProgramName:
Version:	$version
Contact:	Wen Ping <wenp\@biomarker.com.cn> 
Program Date:   2013.12.3
Description:	this program is used to cluster data and draw heatmap
Usage:
Usage:
  Options:\n
        -i                input file              forced\n
        -fn               int                     number of columns ignored, count from the last column\n
        -od               output dir              forced\n
        -h                Help\n
---------------------------------------------------------------------------------------------------------------

USAGE
	print $usage;
	exit;
}
