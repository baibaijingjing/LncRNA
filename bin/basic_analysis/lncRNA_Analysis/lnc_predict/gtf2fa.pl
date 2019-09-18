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
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $program_name=basename($0);


my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"fa=s","gtf=s","o=s","h");
if (!defined($opts{fa})||!defined($opts{gtf})||(!defined($opts{o}))) {
	&help();
}
my $GTF_TO_FASTA = $config{gtf_to_fasta};
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my $gtf=$opts{gtf};
my $fa=$opts{fa};
my $out=$opts{o};
my $temp = "temp_file.fa";
`$GTF_TO_FASTA 	$gtf $fa $temp `;

open OUT, ">$out" or die $!;
my %fasta = FaParse($temp);
#print Dumper \%fasta;
foreach my $head (keys(%fasta)) {
#	print "$head \n";
#	my $oldid = $head;
#	my ($id) = (split /\s+/, $head)[1];
my ($locus,$id) = (split /\s+/, $head)[2,1];
	$locus=~s/\+//g;
	$locus=~s/-//g;
	print OUT  ">$id:$locus\n$fasta{$head}\n";
}
#`rm  $temp `;

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub FaParse{
	my $fa = shift;
	my %fa;
	open IN, $fa || die;
	$/='>';
	<IN>;
	while (<IN>) {
    chomp;
 	my ($id,$seq)=split /\n+/,$_,2;
#	my $seq_id=(split /\s+/,$id)[0];
	$seq=~s/\s+//g;
	#$seq = ~s/atcguUnN/ATCGTTAA/;
	$fa{$id} = $seq	;
	}
	$/='\n';
	close IN ;
	return %fa;
}

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

		-gtf     <STR>  gtf file                         must be given;
		-fa    <STR>  genome fasta file                must be given;
		-o     <STR>  output File						must be given;


	Usage End.
		exit;
}
