
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($in,$lnc,$out,$k);
GetOptions(
				"help|?" =>\&USAGE,
				"in:s"=>\$in,
				"lnctar:s"=>\$lnc,
				"out:s"=>\$out,
				"k:s"=>\$k,
				) or &USAGE;
&USAGE unless ($in and $lnc and $k);
$in=&ABSOLUTE_DIR($in);
#$out=&ABSOLUTE_DIR($out);
$lnc=&ABSOLUTE_DIR($lnc);

open LN,"$lnc" or die $!;
 my %lnta;
while (<LN>){
    chomp;
    next if (/^#/ || /^$/);
    my ($l,$t)=(split/\s+/,$_)[0,1];
    my @tar=split /;/,$t;
    #chop(@tar);
    $lnta{$l}=\@tar;
    
}
close LN;


my $dir=dirname($in);
#print "$dir\n";
my $outfile="$dir/$k".".deg_lnc2Target_m_Cytoscape.input.txt";
#print "$outfile\n";
open IN,$in or die $!;
open OUT,">$outfile";
while (<IN>){
    chomp;
    next if (/^#/ || /^$/);
    my $id=(split /\s+/,$_)[0];
    if (exists $lnta{$id}){
	my @tt=@{$lnta{$id}};
	#print Dumper(\@tt);
	#shift(@tt);
	my $len=@tt;
	foreach my $tid (@tt){
	    print OUT "$id\t$tid\t$len\n";
	}
    }
}
close IN;
close OUT;
###############################################################
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
   Program: $Script
   Version: $version
     Usage:
            --in	<DIR>	input file group_name.DEG_final.xls 		must bu given
	    --lnctar	<file>	lncRNA target file			must be given
	    --k		<str>	L01_vs_L02
            --h                 help documents

   Example:
            perl $Script --in  Lnc_Diff_Analysis --lnctar ./novel_lncRNA_target.xls --k L01_vs_L02

----------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
