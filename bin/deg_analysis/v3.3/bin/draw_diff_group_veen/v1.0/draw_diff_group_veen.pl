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
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);
#my $vennDiagram	= "/share/nas2/genome/bmksoft/tool/VennDiagram/v1.1/vennDiagram.pl";
my $vennDiagram	= "$Bin/vennDiagram.pl";

my $ver="1.0";
############################################
my %opts;
my ($id,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($od and $id);

$od=abs_path($od);
mkdir $od unless (-d $od);
$id=&ABSOLUTE_DIR($id);

###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
#print Dumper(\%trans);die;
my %DEG;
my %ven;
my @DE_DIR=glob "$id/*_vs_*";
my $de_limit=@DE_DIR;
if ($de_limit==1 || $de_limit>5){
	print "WARNING: too more DEG sets, please select sets to plot.\n";
	exit;
}
if( ($de_limit>=2) && ($de_limit<=5) ) {
        my %SX;
	foreach my $dir (@DE_DIR) {
                my $name=basename $dir;
                $DEG{1}{$name}=$name;
                open (IN,"$dir/$name.DEG_final.xls") or die $!;
                while (<IN>) {
			chomp;
			next if (/^$/ || /^#/);
                        next if $.==1;
			$DEG{$.}{$name}=(split/\s+/,$_)[0];
                        $SX{$name}=$.;
                        my $ID=(split/\s+/,$_)[0];
                        push (@{$ven{$name}},$ID);
                }
                close IN;
        }
	open(OUT,">$od/All_DEG_veen.genes") or die $!;
	my $line;
        foreach my $num (sort {$a<=>$b} keys %DEG) {
                my $line;
                foreach my $name (sort {$SX{$b}<=>$SX{$a}} keys %SX) {
                        $line.="$DEG{$num}{$name}\t" if exists $DEG{$num}{$name};
                        $line.="\t" unless exists $DEG{$num}{$name};
                }
                $line=~s/\t$/\n/;
                print OUT $line;
        }
        close OUT;
}
#print Dumper(\%lnc_ids);die;
my $temp_file = "$od/temp"; 
`mkdir $temp_file ` unless (-d $temp_file);
my @cmd;
my $j = 0;
foreach my $deg (keys(%ven)) {
	print "$deg\n";
	my $len=@{$ven{$deg}};
	if ($len!=0){
		$j +=1 ;
		my $outfile = "$temp_file/$deg";
		open  OUT ,">$outfile" or die;
		foreach my $deg_id (@{$ven{$deg}}) {
			print OUT "$deg_id\n";
		}
		close OUT ;
		push (@cmd ,"-T$j",$outfile,"-ID$j",$deg);
	}
}
`perl $vennDiagram @cmd  -od $od `;

`rm -r 	$temp_file ` ;
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
	$id = (split ":",$id)[0];
    $seq=~s/\s+//g;
	$seq = ~tr/atcguUnN/ATCGTTAA/;
	$fa{$id} = $seq	;
	}
	$/='\n';
	close IN ;
	return %fa;
}
sub ABSOLUTE_DIR
{
        my ($in)=@_;
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
        #chdir $cur_dir;
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
	-id	<dir>	input dir  ;
	-od	<STR>	Out Dir;

	Usage End.
		exit;
}
