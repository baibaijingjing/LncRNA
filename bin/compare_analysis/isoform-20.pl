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
my (%mha,%lncha,%mRNAha,%lnc,%mrna);
&get_geneID($merged,\%mha);
&get_geneID($lnc_gtf,\%lncha);
foreach my $id (keys %mha){
	if (!exists $lncha{$id}){
		$mRNAha{$id}=$mha{$id};	
	}
}
open (OUT2,">$od/isoform_plot.txt");
print OUT2 "#Sample\tCount\tGeneNum\n";
open (OUT,">$od/mRNA_isoform.txt");
foreach my $key (sort keys %mRNAha){
       my @list=keys$mRNAha{$key};
        my $num=@list;
	if (!exists $mrna{$num}){
		$mrna{$num}=1;
	} else {
		$mrna{$num}++;
	}
        #print OUT "$key\t$num\n";
	#print OUT2 "mRNA\t$num\n";
}
foreach my $key (sort keys %mrna) {
	print OUT "$key\t$mrna{$key}\n";
	}
close OUT;
open (OUT1,">$od/lncRNA_isoform.txt");
foreach my $key (sort keys %lncha){
        my @list=keys$lncha{$key};
        my $num=@list;
	if (!exists $lnc{$num}){
                $lnc{$num}=1;
        } else {
                $lnc{$num}++;
        }
        #print OUT1 "$key\t$num\n";
	#print OUT2 "lncRNA\t$num\n";
}
foreach my $key (sort keys %lnc) {
        print OUT1 "$key\t$lnc{$key}\n";
}
close OUT1;
my(%mm,%ll);
&HIST("mRNA","$od/mRNA_isoform.txt","$od/mRNA_iso_count.txt",\%mm);
&HIST("lncRNA","$od/lncRNA_isoform.txt","$od/lncRNA_iso_count.txt",\%ll);

#foreach my $key (sort keys %mrna){
#	if (exists $lnc{$key}){
#		print OUT2 "mRNA\t$key\t$mrna{$key}\nlncRNA\t$key\t$lnc{$key}\n";
#
#	} else {
#		print OUT2 "mRNA\t$key\t$mrna{$key}\nlncRNA\t$key\t0\n";
#	}
	
#}
#close OUT1;
#my %hash;
#close OUT2;
#my @num=keys%mm;
#my $n=@num;
#print $n;
 
#foreach my $key (sort keys %mm){
#	if (exists $ll{$key}){
#		print OUT2 "mRNA\t$key\t$mm{$key}\nlncRNA\t$key\t$ll{$key}\n";
#	}
	
#}
`cat $od/mRNA_iso_count.txt $od/lncRNA_iso_count.txt >> $od/isoform_plot.txt`;
close OUT2;
#foreach my $dir (glob "$od/*iso_count.txt"){

#}
`$config{Rscript} $Bin/dodgedBar.2.r --infile $od/isoform_plot.txt --outfile $od/lnc_vs_mRNA.isoform.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"Isoform Number per Gene\" --group.lab \"Sample\" --y.lab \"log2(Number of Gene)\" --title.lab \"Isoform Density\" --legend.col 1 --no.grid`;
sub HIST{
	my ($type,$file,$out,$len)=@_;
	my %len;
	open IN ,"$file";
	$len{20}=$len{40}=$len{60}=$len{80}=$len{100}=$len{120}=$len{140}=$len{160}=$len{180}=$len{200}=$len{220}=$len{240}=$len{260}=$len{280}=$len{300}=$len{'>=300'}=0;
	while (<IN>){
        	my @tmp= (split /\s+/,$_);
	        if ($tmp[0]>0 && $tmp[0]<=20){$len{20}=$len{20}+int($tmp[1]);}     if ($tmp[0]>20 && $tmp[0]<=40){$len{40}=$len{40}+int($tmp[1]);}
        	if ($tmp[0]>40 && $tmp[0]<=60){$len{60}=$len{60}+int($tmp[1]);}   if ($tmp[0]>60 && $tmp[0]<=80){$len{80}=$len{80}+int($tmp[1]);}
	        if ($tmp[0]>80 && $tmp[0]<=100){$len{100}=$len{100}+int($tmp[1]);} if ($tmp[0]>100 && $tmp[0]<=120){$len{120}=$len{120}+int($tmp[1]);}
        	if ($tmp[0]>120 && $tmp[0]<=140){$len{140}=$len{140}+int($tmp[1]);}        if ($tmp[0]>140 && $tmp[0]<=160){$len{160}=$len{160}+int($tmp[1]);}
	        if ($tmp[0]>160 && $tmp[0]<=180){$len{180}=$len{180}+int($tmp[1]);}        if ($tmp[0]>180 && $tmp[0]<=200){$len{200}=$len{200}+int($tmp[1]);}
        	if ($tmp[0]>200 && $tmp[0]<=220){$len{220}=$len{220}+int($tmp[1]);}        if ($tmp[0]>220 && $tmp[0]<=240){$len{240}=$len{240}+int($tmp[1]);}
	        if ($tmp[0]>240 && $tmp[0]<=260){$len{260}=$len{260}+int($tmp[1]);}        if ($tmp[0]>260 && $tmp[0]<=280){$len{280}=$len{280}+int($tmp[1]);}
        	if ($tmp[0]>280 && $tmp[0]<=300){$len{300}=$len{300}+int($tmp[1]);}        if ($tmp[0]>=300){$len{'>=300'}=$len{'>=300'}+int($tmp[1]);}

	}
	close IN;
	open OUT ,">$out";
	#print OUT "#Sample\tCount\tGeneNum\n";
	print OUT "$type\t0-20\t$len{20}\n$type\t21-40\t$len{40}\n$type\t41-60\t$len{60}\n$type\t61-80\t$len{80}\n$type\t81-100\t$len{100}\n$type\t101-120\t$len{120}\n$type\t121-140\t$len{140}\n$type\t141-160\t$len{160}\n$type\t161-180\t$len{180}\n$type\t181-200\t$len{200}\n$type\t201-220\t$len{220}\n$type\t221-240\t$len{240}\n$type\t241-260\t$len{260}\n$type\t261-280\t$len{280}\n$type\t281-300\t$len{300}\n$type\t>=300\t$len{'>=300'}\n";
	#foreach my $key (sort keys %len){
	#	print OUT "$type\t$key\t$len{$key}\n";
	#}
	close OUT;
	return %len;
} 
#`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/dodgedBar.r --infile $od/isoform_plot.txt --outfile $od/all.isoform.png --x.col 2 --group.col 1 --y.col 3 --x.lab \"Isoform Number per Gene\" --group.lab \"Sample\" --y.lab \"Number of Gene\" --title.lab \"Isoform Density\" --legend.col 1`;
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
    
