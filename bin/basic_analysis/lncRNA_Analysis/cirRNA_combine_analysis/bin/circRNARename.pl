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
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$fOut,
                "i:s"=>\$fIn,
                ) or &USAGE;
&USAGE unless ($fIn and $fOut);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------

my %hash;
open (IN,$fIn) or die $!;
open (OUT,">$fOut") or die $!;
my $head = <IN>;
print OUT "New_ID\t";
print OUT "$head";


#$/=">";
my $new_id;
while (<IN>) {
    chomp;
	my @lines=split(/\t/,$_);
	next if($lines[8]=~m/intergenic region/);
	my $gene = $lines[9];
	$gene=~s/,$//;
	$gene=~s/,/_/g;
	if (exists $hash{$gene}) {
		$hash{$gene}++;
		$new_id="$gene"."_"."$hash{$gene}";
		print OUT "$new_id\t$_\n";

	}else{
		$hash{$gene}=1;
		$new_id="$gene"."_"."$hash{$gene}";
		print OUT "$new_id\t$_\n";
	}

}

close IN;
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

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn> 
Program Date:   2012.07.02
      Modify:   
 Description:   This program is used to ......
       Usage:
        Options:
        -i <file>   input file,xxx format,forced
        -o <file>   output file,optional
        -h      help

USAGE
    print $usage;
    exit;
}
