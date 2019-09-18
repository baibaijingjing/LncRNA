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
my ($in,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"in:s"=>\$in,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($in and $o);
$in=&ABSOLUTE_DIR($in);
$o=&ABSOLUTE_DIR($o);

my %hash;

	

open (IN,$in) or die $!;
while (<IN>) {
        chomp;
        next if /^\#/;
        my @tem=split /\s+/,$_;
        if (exists $hash{$tem[0]}){
                if ($tem[1]=~/;/){
                        my @tar=split /;/,$tem[1];
                        foreach my $gene (@tar){
                                if (!exists $hash{$tem[0]}{$gene}){
                                        $hash{$tem[0]}{$gene}=1;
                                }
                        }
                }else{
                        $hash{$tem[0]}{$tem[1]}=1;
                }
        }else{
                if ($tem[1]=~/;/){
                        my @tar=split /;/,$tem[1];
                        foreach my $gene (@tar){
                                if (!exists $hash{$tem[0]}{$gene}){
                                        $hash{$tem[0]}{$gene}=1;
                                }
                        }
                }else{
                        $hash{$tem[0]}{$tem[1]}=1;
                }
        }		
}

open OUT,">$o" or die;
foreach my $lnc (sort keys %hash){
	foreach my $t (keys %{$hash{$lnc}}){
		print OUT "$t\n";
		}
	}

close OUT;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Usage:
  Options:
  -in <dir>  input dir,forced 
  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
