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
my ($gtf,$track,$out,$isoform_len,$isoform_num,$isoform_cov);
GetOptions(
				"help|?" =>\&USAGE,
				"gtf:s"=>\$gtf,
				"track:s"=>\$track,
				"l:s"=> \$isoform_len,
				"n:s"=> \$isoform_num,
				"c:s"=> \$isoform_cov,
				"out:s"=>\$out,
				) or &USAGE;
&USAGE unless ($gtf and $out and $track);

$isoform_len||=200;
$isoform_num||=	2;
$isoform_cov||= 2;
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my %trans;
open TRACK,"<$track" ||die "$!";
while (<TRACK>) {
	chomp;
	next if (/\#/||/^$/||$.==1) ;
	my ($trans_id,$trans_len,$trans_cov,$trans_fpkm) = (split /\t/,$_)[0,7,8,9];
#	print join "\t" ,($trans_id,$trans_len,$trans_cov,"\n");
	$trans{$trans_id}{'length'} = $trans_len;
	$trans{$trans_id}{'coverage'} = $trans_cov;
	$trans{$trans_id}{'fpkm'} = $trans_fpkm;
}
close TRACK;

#print Dumper(\%trans);die;


open GTF,"$gtf"||die "$!";
while (<GTF>) {
	chomp;
	next if (/\#/||/^$/) ;
	my $line = $_;
	my ($trans_id) = $line =~/transcript_id "(\S+)";/;
	my $feature =  (split '\t',$line)[2];
	if ($feature eq 'transcript') {
		$trans{$trans_id}{'exon'} = 0;
	}
	elsif($feature eq 'exon'){
		$trans{$trans_id}{'exon'} += 1;
	}
}
close GTF;

my $trans_num = keys(%trans);
print "There are $trans_num Transcriptions in the Lib\n";
foreach my $trans_id (keys (%trans)) {
	if ($trans{$trans_id}{'length'} <= $isoform_len) {
		delete $trans{$trans_id};
	}
	elsif ($trans{$trans_id}{'coverage'} < $isoform_cov) {
		delete $trans{$trans_id};
	}
	elsif($trans{$trans_id}{'exon'} <= $isoform_num) {
		delete $trans{$trans_id};
	}
}

$trans_num = keys(%trans);
print "Step1:(Length and Exon Number and coverage):\n Ther are $trans_num Transcriptions left\n";

open GTF,"<$gtf"||die "$!";
open OUT,">$out"||die $!;
while (<GTF>) {
	chomp;
	next if (/\#/||/^$/) ;
	my $line = $_;
	my ($trans_id) = $line =~/transcript_id "(\S+)";/;
	if ( exists $trans{$trans_id} ) {
		print OUT "$line\n";
	}
}
close GTF;
close OUT ;

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
sub USAGE{
	print << "	Usage End.";
Description:The process get cds or Intro sequence from genome fasta file base on gff file.
version:$ver
Usage:
	-gtf	<STR>	gtf file                         must be given;
	-track	<STR>	track file                       must be given;
	-l	<NUM>	transcription length;
	-n	<NUM>	Exon Num;
	-c	<NUM>	transcription coverage;
	-out	<STR>	Out file;

	Usage End.
		exit;
}