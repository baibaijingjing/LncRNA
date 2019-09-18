use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
my ($in,$out,$fa,$pfam,$cpat,$cnci);
GetOptions(
    "i:s" =>\$in,
    "o:s"=>\$out,
	"fa:s"=>\$fa,
	"pfam:s"=>\$pfam,
	"cpat:s"=>\$cpat,
	"cnci:s"=>\$cnci,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in and $out and $fa);
$in=&ABSOLUTE_DIR($in);
#$out=&ABSOLUTE_DIR($out);
$fa=&ABSOLUTE_DIR($fa);

#&MKDIR("$out/Pfam");
my %ha;
open IN,"$in";
while (<IN>){
	chomp;
	next if (/^$/ || /^#/);
	my @tem=split /\s+/,$_;
	my $id=(split /:/,$tem[0])[0];
	if (!exists $ha{$id}){
		$ha{$id}=1;
	}
}
close IN;
print Dumper(\%ha);
my $num=keys%ha;
print "there $num candicate fa seq have pfam protein seq or domain \n";
open FA,$fa;
$/='>';
open OUT,">$out";
<FA>;
while (<FA>){
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split/\n+/,$_,2;
	print "$head\n";
	my $ID=(split/:/,$head)[0];
	print "$ID\n";

	if (defined $cpat || defined $cnci){
		if (exists $ha{$ID}){
			print OUT ">$head\n$seq";
		}
	}

	if (defined $pfam){
		if (!exists $ha{$ID}){
                        print OUT ">$head\n$seq\n";
                }
	}

}
$/="\n";
close OUT;
close FA;



sub MKDIR{
        my $dir=@_;
        mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir=`pwd`;chomp($cur_dir);
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
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: version

     Usage:
            -i     predict result of pfam or cpat
		-o 	filtet.fa
		-fa	pre_filter.fa
		-cpat	filter the cpat result
		-pfam	filter the pfam result
            -h                 help documents

   Example:
            perl $Script -i pfam.txt|cpat.txt -o filter.fa -fa pre_filter.fa -cpat|-pfam
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
