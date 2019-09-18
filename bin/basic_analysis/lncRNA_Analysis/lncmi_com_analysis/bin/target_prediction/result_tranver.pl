#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Data::Dumper;
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my $version="v2.0.3";
#my %config=%{readconf("$Bin/../../cfg.cfg")};
my ($mi_tar,$m_tar,$od,$in1,$in2);
GetOptions(
    "mi_tar:s" =>\$mi_tar,
    "out:s"=>\$od,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($mi_tar and $od);
#&MKDIR($od);
#$od=&ABSOLUTE_DIR($od);
$mi_tar=&ABSOLUTE_DIR($mi_tar);

open IN,"$mi_tar" or die $!;
open OUT1, ">$od";
print OUT1 "#lncRNA_ID\ttarget_miRNA\n";
my %hash;
while (<IN>){
	chomp;
	next if (/^#/ or /^$/);
	my ($mi,$tem)=(split /\s+/,$_)[0,1];
	#print "$mi\n";
	if ($tem=~/;/){
                        my @tar=split /;/,$tem;
                        foreach my $gene (@tar){
                                if (!exists $hash{$gene}){
					my @mio=();
					$mio[0]=$mi;
                                        $hash{$gene}=\@mio;
                                }else{
					my $mir=$hash{$gene};
					push(@$mir,$mi);
				}
                        }
                
	}else{
		if (!exists $hash{$tem}){
			my @mio=();
			$mio[0]=$mi;
			$hash{$tem}=\@mio;
		}else{
			my $mir=$hash{$tem};
			push(@$mir,$mi);
		}
	}
	
}
close IN;

foreach my $l (keys %hash){
        my $mirna=$hash{$l};
        
        print OUT1 "$l\t".join(";", @$mirna)."\n";
        
}
sub MKDIR{
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}



sub ABSOLUTE_DIR{ 
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
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: 1.0
     Usage:
            -mi_tar      <FILE>  ./*.mir2target.list
            -out          <FILE>
	    -h 		help documents

   Example:
            perl $Script -mi_tar ./*.mir2target.list  -out ./lncRNA_target2mirna.list 
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
