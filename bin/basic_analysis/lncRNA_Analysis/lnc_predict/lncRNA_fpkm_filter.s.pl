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
my ($cfg,$out,$indir,$lnc_list,$p,$fa,$idr,$od);
GetOptions(
                                "help|?" =>\&USAGE,
                                #"cfg:s"=>\$cfg,
                                "p:s"=>\$p,
                                "out:s"=>\$out,
				"od:s"=>\$od,
				"indir:s"=>\$indir,
				"fa:s"=>\$fa,
				"idr:s"=>\$idr,
				"lnc_list:s"=>\$lnc_list,
                                ) or &USAGE;
&USAGE unless ($out and $indir and $fa and $od);
#$cfg=&ABSOLUTE_DIR($cfg);
$fa=&ABSOLUTE_DIR($fa);
$indir=&ABSOLUTE_DIR($indir);
#$lnc_list=&ABSOLUTE_DIR($lnc_list);
mkdir $od unless (-d $od);
$od=&ABSOLUTE_DIR($od);
$p||=0.1;


#####
my %FA;
&LOAD_SEQ($fa,\%FA);
######################get count and fpkm
system "perl $Bin/bin/abstract_read_count_from_dir.pl -indir $indir -out $od/All_trans_counts.list";
print "perl $Bin/bin/abstract_read_count_from_dir.pl -indir $indir -out $od/All_trans_counts.list \n";
system "perl $Bin/bin/count_to_expression.pl -i $od/All_trans_counts.list -o $od/all_trans_fpkm.list";
print "perl $Bin/bin/count_to_expression.pl -i $od/All_trans_counts.list -o $od/all_trans_fpkm.list \n";
#system "perl $Bin/data_extract_by_ids.pl -idfile $lnc_list -destfile $od/all_trans_fpkm.list -out $od/LncRNA_fpkm.list";
#######################
#filetr the fpkm 
open IN2,"$od/all_trans_fpkm.list" or die "$!: $od/all_trans_fpkm.list\n";
while (<IN2>){
	chomp;
	next if (/^\s+/ or /^#/);
	my @A=split/\t/,$_;
        my $name=shift @A;
        #my $len=pop @A;
	my $num=@A;
		my $q=0;
		for (my $i=0;$i<@A ;$i++){
			if ($A[$i] < $p){
				$q++;
			}
        	}
		if ($q=$num){
			if (exists $FA{$name}){
				delete $FA{$name};
			}	
		}	
	
}
close IN2;


###################
open OUT1,">$out" or die $!;
#print OUT1 "#ID\n";
open IN,$fa or die $!;
$/='>';
<IN>;
while (<IN>) {
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	my ($id,$seq)=split /\n+/,$_,2;
	my $seq_id=(split /\s+/,$id)[0];
	if (exists $FA{$seq_id}){
		print OUT1 ">$id\n$seq\n";
	}
}
$/="\n";
close IN ;
close OUT1;


#`rm $od/All_trans_counts.list $od/all_trans_fpkm.list`;


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


sub LOAD_SEQ{
	my ($fa,$info) = @_;
	open IN, $fa || die;
	$/='>';
	<IN>;
	while (<IN>) {
    	chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
    	my ($id,$seq)=split /\n+/,$_,2;
   		my $seq_id=(split /\s+/,$id)[0];
		$info->{$seq_id}=$seq;
	}
	$/="\n";
	close IN ;
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
            --od        <DIR>   analysis output directory ./Lnc_filter
	    --indir     <DIR>	LncExpression dir
            --fa	<file>  merged_filter_tmp.fa
	    --out	<file>	merged_filter.fa
            --h                 help documents

   Example:
	perl $Script --od fpkm outdir/ --indir ./Cuffnorm/ --lnc_list lncRNA_id_list 

----------------------------------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
