#!/usr/bin/perl -w
#
# Copyright (c) BMK 2011
# Writer:         He hua <heh@biomarker.com.cn>
# Program Date:   2011.
# Modifier:       He hua <heh@biomarker.com.cn>
# Last Modified:  2011.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"fa=s","g=s","intro=s","l=s","gene","cds=s","h");
if (!defined($opts{fa})||!defined($opts{g})||(!defined($opts{intro})&& !defined$opts{cds})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my $gff=$opts{g};
my $fa=$opts{fa};
my $out=$opts{intro} if(defined $opts{intro});
my $extend=$opts{l} || 0;
my $geneRNA=0;
$geneRNA=1 if(defined $opts{gene});
my %site;
my %strand;
my $GR=0;
open GFF,"$gff"||die "$!";
my $geneids;
while (<GFF>) {
	chomp;
	next if (/\#/||/^$/) ;
	my ($chr,$type,$start,$end,$strand,$ID)=(split /\s+/,$_)[0,2,3,4,6,8];
	$GR=1 if($type=~m/gene/ && $geneRNA==1);
	if ($GR==1 && $type=~m/gene/i) {
		$geneids=(split/\;|=/,$ID)[1];
	}elsif($GR==0 && $type=~m/mRNA/i){
		$geneids=(split/\;|=/,$ID)[1];
	}
	if ($type=~m/exon|CDS/i) {
		$site{$chr}{$geneids}{$start}{$end}=1;
		$strand{$chr}{$geneids}=$strand;
	}
}
close GFF;
my %Introsite;my %cds;
foreach my $chr (keys %site) {
	foreach my $geneid (keys %{$site{$chr}}) {
		my @starts=sort{$a<=>$b} keys %{$site{$chr}{$geneid}};
		my @sites;
		foreach my $start (@starts) {
			my $exon_end;
			my @end=keys %{$site{$chr}{$geneid}{$start}};
			if (@end>1) {
				$exon_end=(sort {$a<=>$b} @end)[-1];
			}else{
				$exon_end=$end[0];
			}
			push @sites,$start,$exon_end;
		}
		$cds{$chr}{$geneid}=[@sites] if(defined $opts{cds});
		if (@sites<=2) {
			next;
		}
		#print Dumper @sites;die;
		my $end=$sites[1];
		for (my $i=2;$i<@sites;$i+=2) {
			my $length=$sites[$i]-$end;
			#print "$end\t$length\n";
			if ($length<2) {
				if ($end<$sites[$i+1]) {
					$end=$sites[$i+1];
				}
				next;
			}
			$Introsite{$chr}{$geneid}{$end}=$length-1;
			#print "$end\n";
			$end=$sites[$i+1];
		#	print "$end\t",$length-1,"\n";
		}
	}
}
##foreach  (sort keys %Introsite) {
#	foreach my $geneid (sort keys %{$Introsite{"1"}}) {
#		foreach my $site(keys %{$Introsite{"1"}{$geneid}}){
#			print "1\t$geneid\t$site\t",$Introsite{"1"}{$geneid}{$site},"\n";
#		}
#	}
##}
#die;
print "Read file done!\n";
if(defined $opts{cds}){
	open CDS,">$opts{cds}" ||die "$!";
}
if(defined $opts{intro}){
	open OUT,">$out"||die "$!";
}
open FA,"$fa"||die "$!";
$/=">";
while (<FA>) {
	chomp;
	next if(/^$/);
	my ($chr,$seq)=split /\n/,$_,2;
		$chr=(split/\s+/,$chr)[0];
		$seq=~s/\s+//ig;
		my $chr_len=length$seq;
		foreach my $geneid (sort keys %{$cds{$chr}}) {
			if(defined $opts{cds}){
				my $cds;
				next if(!$cds{$chr}{$geneid});
				if(${$cds{$chr}{$geneid}}[0]-$extend<=0){
					$cds=substr($seq,0,${$cds{$chr}{$geneid}}[0]) ;
					print "$chr\t$geneid\t${$cds{$chr}{$geneid}}[0]\n";
				}
				$cds=substr($seq,${$cds{$chr}{$geneid}}[0]-$extend-1,$extend) if(${$cds{$chr}{$geneid}}[0]-$extend>0);
				for (my $i=0;$i<@{$cds{$chr}{$geneid}} ;$i+=2) {
					$cds.=substr($seq,${$cds{$chr}{$geneid}}[$i]-1,${$cds{$chr}{$geneid}}[$i+1]-${$cds{$chr}{$geneid}}[$i]+1);
				}
				if(${$cds{$chr}{$geneid}}[-1]+$extend>$chr_len){
					$cds.=substr($seq,${$cds{$chr}{$geneid}}[-1],$chr_len-${$cds{$chr}{$geneid}}[-1]) ;
					print "$chr\t$geneid\t${$cds{$chr}{$geneid}}[-1]\n";
				}
				$cds.=substr($seq,${$cds{$chr}{$geneid}}[-1],$extend) if(${$cds{$chr}{$geneid}}[-1]+$extend<=$chr_len);
				if ($strand{$chr}{$geneid} eq "-") {
					$cds=~tr/ATCG/TAGC/;
					$cds=reverse($cds);
				}
				my $gene_len=length $cds;
				print CDS ">$geneid $chr $gene_len\n$cds\n";
			}
			if (defined $opts{intro}) {
				my $seqs;
				next if(!$Introsite{$chr}{$geneid});
				foreach my $site (sort {$a<=>$b}keys %{$Introsite{$chr}{$geneid}}) {
					 $seqs.=substr($seq,$site,$Introsite{$chr}{$geneid}{$site});
				}
				if ($strand{$chr}{$geneid} eq "-") {
					if(!$seqs){
						print "$geneid\n" ;
						print Dumper %{$Introsite{$chr}{$geneid}};
						die;
					}
					$seqs=~tr/ATCG/TAGC/;

					$seqs=reverse($seqs);
				}
				next if($seqs!~/A|G|C|T/i);
				print OUT ">$geneid\n";
				print OUT $seqs,"\n";
			}
		}
}
close FA;
close OUT if(defined $opts{intro});
close CDS if(defined $opts{cds});
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub ABSOLUTE_DIR
{
        my ($in,$cur_dir)=@_;
        my $return="";

        if(-f $in)
        {
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }
        elsif(-d $in)
        {
                chdir $in;$return=`pwd`;chomp $return;
        }
        else
        {
                warn "Warning just for file and dir\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime
{ # &Runtime($BEGIN);
        my ($t1)=@_;
        my $t=time()-$t1;
        print "\nTotal elapsed time: ${t}s\n";
}
sub help{
	print << "	Usage End.";
	Description:The process get cds or Intro sequence from genome fasta file base on gff file.
		version:$ver
	Usage:

		-g     <STR>  gff file                         must be given;
		-fa    <STR>  genome fasta file                must be given;
		-gene         base on gene(defualt:mRNA)            choice;
		-cds   <STR>  cds sequence outfile                  choice;
		-intro <STR>  Intro sequence outfile                choice;
		-l     <INT>  Upstream and downstream extend<INT>bp choice;

	Usage End.
		exit;
}