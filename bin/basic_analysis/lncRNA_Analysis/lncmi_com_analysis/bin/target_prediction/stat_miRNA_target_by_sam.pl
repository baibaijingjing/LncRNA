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
my ($o,$i,$t);

GetOptions(
                "help|?" =>\&USAGE,
				"i:s"=>\$i,
				"t:s"=>\$t,
				"o:s"=>\$o,
                ) or &USAGE;
&USAGE unless ($o and $i and $t);

my %mir_target;
open (IN,$t) or die $!;
while (<IN>) {
    chomp;
    my ($mir,$targets)=split(/\t/,$_);
    my @target = split(/\;/,$targets);
    $mir_target{$mir}=\@target;
}
close IN;


my %TPM;
my %own_target;
my %Target;
my @samples;
my $number;;
open (IN,"$i") or die $!;
while (<IN>) {
    chomp;
    if (/\#miRNA/) {
        @samples=split(/\s+/,$_);
        shift @samples;
        pop   @samples;
        $number=@samples;
    }else{
        my @lines=split(/\s+/,$_);
        my $miRNA=shift @lines;
        my $lenth=pop   @lines;
        for (my $i=0;$i<$number ;$i++) {
            if ($lines[$i]!=0) {
                $TPM{$samples[$i]}++;
                if (exists $mir_target{$miRNA}) {
                    $own_target{$samples[$i]}++;
                    foreach my $targ (@{$mir_target{$miRNA}}) {
                        $Target{$samples[$i]}{$targ}=1;
                    }
                }
            }else{
                next;
            }
        }
    }

}
close IN;

open (OUT,">$o") or die $!;
print OUT "BMK-ID\tAll miRNA\tmiRNA with Target\tTarget gene\n";
foreach my $sample (@samples) {
    my $num=keys %{$Target{$sample}};
#     foreach my $ttaa (keys %{$Target{$sample}}) {
#         print "$ttaa\n";
#     }
    print OUT "$sample\t$TPM{$sample}\t$own_target{$sample}\t$num\n";
#die;
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

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Wangmeiyan <wangmy\@biomarker.com.cn> 
Program Date:   2015.2.9
      Modify:   
 Description:   This program is used to stat miRNA target by samples
       Usage:
        Options:
        -t  <file>   *.mir2target.list,forced
	-i  <file>   All_miRNA.expression.list,forced
	-o  <file>   output file,forced
        -h      help

USAGE
    print $usage;
    exit;
}


