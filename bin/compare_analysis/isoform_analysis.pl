#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Data::Dumper;
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};
my ($merged,$lnc_gtf,$od);
GetOptions(
    "merged:s" =>\$merged,
    "lnc_gtf:s"=>\$lnc_gtf,
    #"mRNA:s"=>\$mRNA_gtf,
	"od:s"=>\$od,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($merged and $lnc_gtf and $od);
$merged=&ABSOLUTE_DIR($merged);
$lnc_gtf=&ABSOLUTE_DIR($lnc_gtf);
#$mRNA_gtf=&ABSOLUTE_DIR($mRNA_gtf);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
my (%mha,%lncha,%mRNAha);
&get_geneID($merged,\%mha);
&get_geneID($lnc_gtf,\%lncha);
foreach my $id (keys %mha){
	if (!exists $lncha{$id}){
		$mRNAha{$id}=$mha{$id};	
	}
}
open (OUT2,">$od/isoform_plot.txt");
open (OUT,">$od/mRNA_isoform.txt");
foreach my $key (sort keys %mRNAha){
       my @list=keys$mRNAha{$key};
        my $num=@list;
        print OUT "$key\t$num\n";
	print OUT2 "mRNA\t$num\n";
}
close OUT;
open (OUT1,">$od/lncRNA_isoform.txt");
foreach my $key (sort keys %lncha){
        my @list=keys$lncha{$key};
        my $num=@list;
        print OUT1 "$key\t$num\n";
	print OUT2 "lncRNA\t$num\n";
}
close OUT1;
#my %hash;
close OUT2;
`$config{Rscript} $Bin/iso_oneFactorBox.r --infile $od/isoform_plot.txt --outfile $od/lnc_vs_mRNA.isoform.png --value.col 2 --x.col 1 --x.lab type --y.lab isoform --title.lab lnc_vs_mRNA`;
sub get_geneID{
	my ($file,$hash)=@_;
	my %hash;
	open(IN,$file )or die "$!";
    #$\="\n"
	while (<IN>) {
		chomp;
		next if (/^\s+/ or /^#/);
        #push(@sam,$_);
    		my @tem=split("\"",$_);
		my $geneid=$tem[1];
		my $tranid=$tem[3];
		if (!exists $hash{$geneid}){
			$hash->{$geneid}{$tranid}=1;
		}
		if (exists $hash{$geneid}){
			if (!exists $hash{$geneid}{$tranid}){
				$hash->{$geneid}{$tranid}=1;
			}else{$hash->{$geneid}{$tranid}++;}
		}

	
	}
	close IN;
}
sub MKDIR{
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
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
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: Script
   Version: version
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            -merged      <GTF FILE>  ./merged.gtf
            -lnc_gtf	<GTF FILE>   ./lnc_filter_final.gtf
            -od          dir
	    -h 		help documents

   Example:
            perl Script -merged merged.gtf -lnc_gtf lnc_filter_final.gtf  -od  ./
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
    
