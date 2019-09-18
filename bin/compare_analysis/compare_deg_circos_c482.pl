#!/usr/bin/perl -w
use strict;
#use newPerlBase;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $Title = "compare_deg_circos";
my $version="1.0.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($deg_result,$all_gene,$odir,$lncdeg_result,$lncall_gene,$chr,$outputfile);

GetOptions(
    "chr:s"     =>\$chr,
    "deg:s"     =>\$deg_result,
    "all:s"     =>\$all_gene,
    "lncdeg:s"  =>\$lncdeg_result,
    "lncall:s"  =>\$lncall_gene,
    "od:s"      =>\$odir,
    "outfile:s" =>\$outputfile,
    "help|h"    =>\&USAGE,
    ) or &USAGE;
&USAGE unless ( $deg_result and $lncdeg_result and $odir and $all_gene and $lncall_gene and $chr);

#$outputfile=$outputfile|| circos;

$odir=&ABSOLUTE_DIR($odir);
system "mkdir -p $odir" unless (-d $odir);

my $link = "$odir/$outputfile.locus_mRNA.txt";
my $lnc_link = "$odir/$outputfile.locus_lnc.txt";

my %gene_locus_hash;
my %gene_locus_hash_lnc;
my %chr_hash_num;
my %chr_hash_chanum;
my %chr_hash_cha;
my @tmp;


###### sort chr #####
open(IN1,"$chr")||die "open $chr failed!\n";
while (<IN1>)
{
        chomp;
        next if ($_=~/\#/);
        next if ($_=~/^$/);
        next if ($_=~/^\s*$/);
        my @data = split /\s+/,$_;
        if ($data[0] =~ /^\d+$/){
        	$chr_hash_num{$data[0]}=$data[1] ;
        }elsif(@tmp = $data[0] =~ /^\D+(\d+)$/){
        	$chr_hash_chanum{$tmp[0]}{$data[0]} = $data[1];
        }else{
        	$chr_hash_cha{$data[0]}=$data[1] ;
        }
        #$chr_hash{$data[0]}=$data[1] ;
}
close(IN1);
#print Dumper(\%chr_hash_num);
#print Dumper(\%chr_hash_chanum);

my $chr_sort = "$odir/chr_sort.txt";
open (CHR,">$chr_sort")||die "open or creat $chr_sort failed!\n";
if( %chr_hash_num){
	foreach my $chr_tmp (sort  {$a <=> $b} keys %chr_hash_num){
		print CHR "$chr_tmp\t$chr_hash_num{$chr_tmp}\n";
	}
}
if(%chr_hash_chanum){
	foreach my $chr_tmp (sort  {$a <=> $b} keys %chr_hash_chanum){
		foreach my $chr_key (keys %{$chr_hash_chanum{$chr_tmp}}){
			print CHR "$chr_key\t$chr_hash_chanum{$chr_tmp}{$chr_key}\n";
		}
	}
}
if( %chr_hash_cha){
	foreach my $chr_tmp (sort keys %chr_hash_cha){
		print CHR "$chr_tmp\t$chr_hash_cha{$chr_tmp}\n";
	}
}
#foreach my $chr_tmp (sort  {lc($a) cmp lc($b)} keys %chr_hash)
#{
#	print CHR "$chr_tmp\t$chr_hash{$chr_tmp}\n";
#}
close(CHR);


###### mRNA deg_result #####
open(IN2,"$all_gene")||die "open $all_gene failed!\n";
while (<IN2>)
{
        chomp;
        next if ($_=~/\#/);
        next if ($_=~/^$/);
        next if ($_=~/^\s*$/);
        my @data=split/\t/,$_;
        $gene_locus_hash{$data[0]}=$data[1] ;
}
close(IN2);


open(IN3,"$deg_result")||die "open $deg_result failed!\n";
open (OUT,">$link")||die "open or creat $link failed!\n";
while (<IN3>) {
	chomp;
	next if ( /^\#/ || /^$/ );
	my @deg = split /\t/,$_;
	if ($gene_locus_hash{$deg[0]}) {
		my $locus = $gene_locus_hash{$deg[0]};
		my @chr = split/\:/,$locus;
		my @position = split/\-/,$chr[1];
		my $start = $position[0];
		my $end = $position[1];
		my $fill_color ;
		my $stroke_color ;
		if ($deg[-1] =~ /up/){
			$fill_color = "fill_color=dred";
			$stroke_color = "stroke_color=dred";
		}elsif($deg[-1] =~ /down/){
			$fill_color = "fill_color=dgreen";
			$stroke_color = "stroke_color=dgreen";
		}
		print OUT "$chr[0]\t$start\t$end\t$fill_color,$stroke_color\n";
	}
}

close IN3;
close OUT;


###### lncRNA deg_result #####
open(IN3,"$lncall_gene")||die "open $lncall_gene failed!\n";
while (<IN3>)
{
        chomp;
        next if ($_=~/\#/);
        next if ($_=~/^$/);
        next if ($_=~/^\s*$/);
        my @data=split/\t/,$_;
        $gene_locus_hash_lnc{$data[0]}=$data[1] ;
}
close(IN3);


open(IN4,"$lncdeg_result")||die "open $lncdeg_result failed!\n";
open (OUT2,">$lnc_link")||die "open or creat $lnc_link failed!\n";
while (<IN4>) {
	chomp;
	next if ( /^\#/ || /^$/ );
	my @deg = split /\t/,$_;
	if ($gene_locus_hash_lnc{$deg[0]}) {
		my $locus = $gene_locus_hash_lnc{$deg[0]};
		my @chr = split/\:/,$locus;
		my @position = split/\-/,$chr[1];
		my $start = $position[0];
		my $end = $position[1];
		my $fill_color ;
		my $stroke_color ;
		if ($deg[-1] =~ /up/){
			$fill_color = "fill_color=vdyellow";
			$stroke_color = "stroke_color=vdyellow";
		}elsif($deg[-1] =~ /down/){
			$fill_color = "fill_color=dblue";
			$stroke_color = "stroke_color=dblue";
		}
		print OUT2 "$chr[0]\t$start\t$end\t$fill_color,$stroke_color\n";
	}
}

close IN4;
close OUT2;


print "\nCMD: perl $Bin/circos_v1.pl --chr $chr_sort --circle $link --type highlight --circle $lnc_link --type highlight --od $odir --outfile $outputfile \n\n";
`perl $Bin/circos_v1.pl --chr $chr_sort --circle $link --type highlight --circle $lnc_link --type highlight --od $odir --outfile $outputfile`;





#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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

#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}



#############################################################################################################
#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------

   Program: $Script
   Version: $version
   Contact: linhj <linhj\@biomarker.com.cn> 
      Date: 2016.04.29

     Usage:
            --chr       <FILE>   chrosome file               [ ./Ref_Genome/genome_size.txt ]
            --od        <DIR>    analysis output directory:  [ ./Compare_Analysis/DEG/T04_T05_vs_T10_T11 ]
            --deg       <file>   the mRNA deg result:        [ ./DEG_analysis/T04_T05_vs_T10_T11.DEG_final.xls ]
            --all       <file>   the mRNA expression file:   [ ./geneExpression/AllSample.genes_expression.xls ]
            --lncdeg    <file>   the lncRNA deg result:      [ ./Lnc_Diff_analysis/T04_T05_vs_T10_T11.DEG_final.xls ]
            --lncall    <file>   the lncRNA expression file: [ ./LncExpression/AllSample.isoforms_expression.xls ]
            --h                  help documents

   Example:
            perl $Script --chr genome_size.txt --od output_dir --deg T04_vs_T10.DEG_final.xls --all AllSample.genes_expression.xls \\
                 --lncdeg lnc_T04_vs_T10.DEG_final.xls --lncall AllSample.isoforms_expression.xls
                 
----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
