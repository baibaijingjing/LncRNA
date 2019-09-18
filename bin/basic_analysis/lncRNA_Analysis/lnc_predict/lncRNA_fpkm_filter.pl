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
## ------------------------------------------------------------------
## GetOptions
## ------------------------------------------------------------------
my ($cfg,$out,$indir,$lnc_list,$p);
GetOptions(
                                "help|?" =>\&USAGE,
                                #"cfg:s"=>\$cfg,
                                "p:s"=>\$p,
                                "out:s"=>\$out,
				"indir:s"=>\$indir,
				"lnc_list:s"=>\$lnc_list,
                                ) or &USAGE;
&USAGE unless ($out and $indir and $lnc_list);
#$cfg=&ABSOLUTE_DIR($cfg);
$indir=&ABSOLUTE_DIR($indir);
$lnc_list=&ABSOLUTE_DIR($lnc_list);
#mkdir $od unless (-d $od);
#$od=&ABSOLUTE_DIR($od);
$p||=0.1;
#my (%data_cfg);
#my $gene_track="$indir/genes.fpkm_tracking";
my $lnc_track="$indir/isoforms.fpkm_tracking";
# read data config
#&data_cfg_read($cfg,\%data_cfg);
#print $sample_num;
#my @sam=keys%data_cfg;
my $num;
#print "$num\n";
my (%lnc_id,%gene_id);
open IN,$lnc_list or die "$!: $lnc_list\n";
while (<IN>){
	chomp;
	#my $id = $_=~/(\S+)/g;
	#print "$_\n";
	$lnc_id{$_}=1;
}
close IN;
open  IN1,$lnc_track or die "$!: $lnc_track \n";
while (<IN1>){
	chomp;
	next if (/^\s+/ or /^#/);
	my @str = split(/\s+/, $_);
	if (!exists $lnc_id{$str[0]}){
		$gene_id{$str[4]}=1;
		
	}
	
}
close IN1;
open IN2,$lnc_track or die "$!: $lnc_track\n";
#open OUT,">$od/lnc_filter_id.list";
open OUT1,">$out" or die $!;
print OUT1 "#ID\n";
while (<IN2>){
	chomp;
	next if (/^\s+/ or /^#/);
	my @str = split /\s+/,$_;
	my $n=@str;
	$num=($n-9)/4;
	if (exists $lnc_id{$str[0]}){
		my $q=0;
		my $line="$str[0]";
		#print "$line";
		for (my $i = 1;$i <= $num;$i++){
			if ($str[($i*4)+5] < $p){
				$q++;
			}
			$line.="\t$str[($i*4)+5]";
        	}
		if ($q<$num){
			#print OUT "$line\n";
			print OUT1 "$str[0]\n";
		}
	} 
		
	
}
close IN2;
#close OUT1;

#####get cfg.cfg 





sub data_cfg_read {
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
			$data_cfg->{$sample_id}=1	
        }
    }
    close CFG;
	#my ($sample_num);
	#my @sam=keys%data_cfg;	
	#$sample_num=@sam;
#    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};
 #   print "sample_number: $data_cfg->{sample_num}\n";

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
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: niulg <niulg\@biomarker.com.cn>
      Date:

     Usage:
            --od        <DIR>   analysis output directory
	    --indir     <DIR>	cuffnorm dir
            --lnc_list   <file>  lncRNA id list
            --h                 help documents

   Example:
	perl $Script --od fpkm outdir/ --indir ./Cuffnorm/ --lnc_list lncRNA_id_list 

----------------------------------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
