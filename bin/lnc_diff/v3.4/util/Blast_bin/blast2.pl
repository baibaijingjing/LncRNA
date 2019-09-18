#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($idir,$query,$program,$od,$view,$value,$Blast_cpu);
GetOptions(
			"help|?"=>\&USAGE,
			"idir:s"=>\$idir,
			"i:s"=>\$query,
			"p:s"=>\$program,
			"od:s"=>\$od,
            "m:s"=>\$view,
			"e:s"=>\$value,
			"cpu:s"=>\$Blast_cpu,
			)or &USAGE;
&USAGE unless ($query and $program and $od and $idir);
`mkdir $od` unless (-d $od);
`mkdir $od/Mid_dir ` unless (-d "$od/Mid_dir");
`mkdir $od/work_sh ` unless (-d "$od/work_sh");
$od = ABSOLUTE_DIR($od);

# $program ||="blastn";
$view ||="0";
$value ||="10";
$Blast_cpu||=20;
my $notename=`hostname`;chomp $notename;
$idir=ABSOLUTE_DIR($idir);
$query=ABSOLUTE_DIR($query);
my $name=basename$query;


################################ Program
#####################ÇÐ¸îÎÄ¼þ
print "Cut fa to dir\nperl /share/bioCloud/wangyj/fa_tools/Blast_bin/cut_fa_file_to_dir.pl -mRNA $query -od $od/Mid_dir\n\n";
`perl /share/bioCloud/wangyj/fa_tools/Blast_bin/cut_fa_file_to_dir.pl -mRNA $query -od $od/Mid_dir `;

#####################################
my @subfiles;
@subfiles = glob("$od/Mid_dir/*.fa");

open OUT,">$od/work_sh/$program.sh" || die "fail $od/$program.sh";
foreach my $subfile (@subfiles) {
	my $basename=basename $subfile;
	print OUT "blastall -b 20 -v 20 -p $program -e $value -F F -d $idir -i $subfile -m $view -a 2 -o $od/Mid_dir/$basename.$program.txt && \n";
}
close OUT;

&Cut_shell_qsub("$od/work_sh/$program.sh",$Blast_cpu,"6G","general.q");
#####################################
my @blastfiles=glob("$od/Mid_dir/*.txt");

open TAB1,">$od/work_sh/Covert_2Tabbest.sh" || die "$!";
open TAB2,">$od/work_sh/Covert_2Tab.sh" || die $!;

    foreach my $blastfile (sort @blastfiles) {
        my $basename = basename($blastfile);
        print TAB1 "perl /share/bioCloud/wangyj/fa_tools/Blast_bin/blast_parser.pl  -eval $value -tophit 1 -m $view -topmatch 1 $blastfile > $od/Mid_dir/$basename.tab.best \n" ;
        print TAB2 "perl /share/bioCloud/wangyj/fa_tools/Blast_bin/blast_parser.pl  -eval $value -tophit 20 -m $view -topmatch 1 $blastfile > $od/Mid_dir/$basename.tab \n" ;
    }
    close TAB1;
    close TAB2;
&Cut_shell_qsub("$od/work_sh/Covert_2Tabbest.sh",$Blast_cpu,"1G","general.q");
&Cut_shell_qsub("$od/work_sh/Covert_2Tab.sh",$Blast_cpu,"1G","general.q");

#############################################

system("cat $od/Mid_dir/*.tab.best >$od/$name.$program.best.xls");
system("cat $od/Mid_dir/*.tab >$od/$name.$program.all.xls");
system("cat $od/Mid_dir/*.txt >$od/$name.$program.raw.txt");

`rm -r $od/Mid_dir `;
`rm -r $od/work_sh `;
##############################################


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

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = system "wc -l $shell";
	if ($line<=1000) {
		if ($notename=~/login-0-4/) {
			system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell ";
		}
		else
		{
			system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell ";
		}
	}
	if ($line>=1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=1) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/login-0-4/) {
				system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file";
			}
			else
			{
				system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file ";
			}
		}
	}
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  blast_supply
Version:	$version
Contact:	lv ran <lvr\@biomarker.com.cn> 
Program Date:   2014/7/30
Usage:
Options:
-idir   input database 
-i      query file
-p      Type of program
	    blastp,blastn,blastx,tblastn,tblastx
-od      output dir
-m      Type of result  
	    0:pairwise
	    8:tabular
-e      Comparing expectations [0-+¡Þ]
-cpu    Blast cpu
-h      Help

USAGE
	print $usage;
	exit;
}
