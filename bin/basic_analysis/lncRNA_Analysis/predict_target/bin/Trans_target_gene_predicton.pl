#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Data::Dumper;
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")};
my ($od,$in);
GetOptions(
    "od:s"=>\$od,
    "in:s"=>\$in,
    #"in2:s"=>\$in2,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in and $od);
$in = &ABSOLUTE_DIR($in);

=c
my $qhost=`qhost`;
my $mem=0;
my $memus=0;
my $compute_id;
my $compute_tmp;
my $bigger;
foreach my $qhost_line (split /\n/,$qhost) {
	if ($qhost_line=~/compute/) {
		my @aline=(split /\s+/,$qhost_line);
        next if $aline[6]=~/-/;
		$aline[7]=~ s/G//g;
#		print "---------------$aline[7]----------------\n";
		if ($aline[7] >= $mem) {
			$mem         = $aline[7];
#            print "==============$mem===============\n";
			$memus       = $aline[8];
#            $compute_tmp = $aline[0];
			$compute_id  = $aline[0];
#            print "==============$compute_id===============\n";
#            print "==============$compute_tmp===============\n";
			if ($aline[7]=$mem) {
                $compute_tmp = $compute_id;
				if ($memus=~/G/) {############## 第一种情况 已用内存单位为G
					if ($aline[8]=~ /M/) {
						$aline[8]=~ s/M//;
						$memus =~ s/G//;
						$aline[8]=$aline[8]/1024;
						if ($aline[8] <= $memus) {
                            $bigger = $aline[0];
							$compute_id=$aline[0];
 #                           print "==============$compute_id===============\n";
						}else{
                            $compute_id=$compute_tmp;
                            $bigger = $compute_tmp;
#                            print "==============$compute_id===============\n";
                        }
					}
					if ($aline[8]=~ /G/) {
						$aline[8]=~ s/G//;
						$memus=~ s/G//;
						if ($aline[8] <= $memus) {
                            $bigger = $aline[0];
							$compute_id=$aline[0];
#                            print "==============$compute_id===============\n";
						}else{
                            $compute_id=$compute_tmp;
                            $bigger = $compute_tmp;
#                            print "==============$compute_id===============\n";
                        }
					}
				}
				if ($memus=~/M/) {################  第二种情况 已用内存单位为M
					if ($aline[8]=~/M/) {
						$aline[8]=~ s/M//;
						$memus=~ s/M//;
						if ($aline[8]<$memus) {
                            $bigger = $aline[0];
							$compute_id=$aline[0];
#                            print "==============$compute_id===============\n";
						}else{
                            $compute_id=$compute_tmp;
                            $bigger = $compute_tmp;
#                            print "==============$compute_id===============\n";
                        }
					}
					if ($aline[8]=~/G/) {
						$aline[8]=~ s/G//;
						$memus=~ s/M//;
						$memus=$memus/1024;
						if ($aline[8]<$memus) {
                            $bigger = $aline[0];
							$compute_id=$aline[0];
 #                           print "==============$compute_id===============\n";
						}else{
                            $compute_id=$compute_tmp;
                            $bigger = $compute_tmp;
#                            print "==============$compute_id===============\n";
                        }
					}
				}
                $compute_tmp = $bigger;
			}
            $compute_tmp = $bigger;
		}
	}
}
my $count_sample=`cd $in && ls -l |grep ".geneExpression.xls\$" |wc -l`;#$od=**/Basic_Analysis/geneExpression
&MKDIR($od);
system "ssh $compute_id";#ssh cluster &&
#system "cd $od";
print "Rscript $Bin/WGCNA_v1.2.R --indir $in --outdir $od --meanFPKM 0.5 -f 0.2 -n $count_sample";
system "Rscript $Bin/WGCNA_v1.2.R --indir $in --outdir $od --meanFPKM 0.5 -f 0.2 -n $count_sample";
=cut




################################# 选择最大可用内存计算节点做WGCNA分析 2016/01/05 ################################################	
my $count_sample=`cd $in && ls |grep ".geneExpression.xls\$" |wc -l`;
&MKDIR($od);
my $qhost=`qhost`;
#print "$qhost";
my $memus=0;
my $compute_id;

foreach my $qhost_line (split /\n/,$qhost) {
#       print "$qhost_line\n";
        if ($qhost_line=~/$config{'compute_node'}/) {
                my @aline=split /\s+/,$qhost_line;
		#print "$aline[-3]\n";
                next if($aline[-3]=~/-/);
#               print "@aline\n";
                if ($aline[-3]=~/G/) {
                        $aline[-4]=~s/G//;
                        $aline[-3]=~s/G//;
                }if($aline[-3]=~/M/){
                        $aline[-4]=~s/G//;
                        $aline[-3]=~s/M//;
                        $aline[-3]=$aline[-3]/1024;
                }
                if($aline[-4]-$aline[-3]>$memus){
                        $memus=$aline[-4]-$aline[-3];
                        $compute_id=$aline[0];
                }
        }
}
print "WGCNA work node in :-----$compute_id-----\n";
open OUT4, ">$od/WGCNA_cmd.sh";
print OUT4 "ssh $config{head_node} << EOF\n";
print OUT4 "ssh $compute_id << EOF\n";
print OUT4 "$config{Rscript} $Bin/WGCNA_v1.2.R --indir $in --outdir $od/ --meanFPKM 5 -f 0.2 -n $count_sample";###od add "/"
print OUT4 "EOF\n";
print OUT4 "EOF\n";
close OUT4;
system "sh $od/WGCNA_cmd.sh";









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
            -od          dir
	    -in         gene_exp 
	    -h 		help documents

   Example:
            perl $Script  -in gene_exp  -od  ./
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
