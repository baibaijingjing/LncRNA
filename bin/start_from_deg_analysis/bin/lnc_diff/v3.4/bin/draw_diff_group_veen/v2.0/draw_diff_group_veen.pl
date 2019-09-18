#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Algorithm::Combinatorics qw(combinations permutations);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $od);

system "mkdir -p $od" unless (-d $od);
$od = &ABSOLUTE_DIR($od);
$id = &ABSOLUTE_DIR($id);

print STDOUT "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n\n";
# ------------------------------------------------------------------
# traversal 2 to 5 sets venn
# ------------------------------------------------------------------
my @DEGsetName;
my @combinations;

for my $dir (glob "$id/*_vs_*") {
    my $name = basename $dir;
    push @DEGsetName, $name;
}

if (@DEGsetName>5) {
    print "WARNING: too more DEG sets, please select sets to plot.\n";
    exit;
} elsif (@DEGsetName>1) {
    for my $n (2..$#DEGsetName+1) {
        next if ($n>5);
        my $iter = combinations(\@DEGsetName,$n); 

        while (my $c = $iter->next) {
             push @combinations, join(";",@$c);
        }
    }

    while (@combinations) {
        my $combination = shift @combinations;
        my $cmd = "perl $Bin/draw_DEGsets_venn_diagram.pl ";
        my (@pre,$pre);

        for my $name (split /;/,$combination) {
            $cmd.= "--lst $id/$name/$name.DEG_final.xls --lab $name ";
            push @pre, $name;
        }

        $pre = join "_n_",@pre;
        $cmd.= "--pre $pre --od $od/$pre ";

        print "[".&GetTime(time())."] $cmd\n";
        system "$cmd";
    }
} elsif (@DEGsetName==1) {
    print "WARNING: only one DEG set, no venn plot.\n";
    exit;
} else {
    print "ERROR: no any DEG set.\n";
    exit;
}


#######################################################################################
print STDOUT "\n[".&GetTime(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
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

################################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:   $Script
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn>/<simonyoung8824\@gmail.com>
Program Date:   2014-10-29
 Description:   This program is used to draw venn diagram of 2-5 differentially expressed gene sets traversally.
       Usage:
        --id       input dir, DEG_Analysis, forced
        --od       output dir, forced

     Example:
        perl $Script --id DEG_Analysis --od DEG_Analysis/Venn

USAGE
    print $usage;
    exit;
}
