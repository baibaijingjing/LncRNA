#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($od,$query,$program,$view,$value,$Blast_cpu,$fyes,$fno,$cfg);
GetOptions(
			"help|?"=>\&USAGE,
			"od:s"=>\$od,
			"yes:s"=>\$fyes,
			"no:s"=>\$fno,
			"fa:s"=>\$query,
			"p:s"=>\$program,
			"m:s"=>\$view,
			"cpu:s"=>\$Blast_cpu,
			"e:s"=>\$value,
			"cfg:s"=>\$cfg,
			)or &USAGE;
&USAGE unless ($query and $od);
&USAGE unless (defined $fyes or defined $fno) ;
`mkdir $od ` unless (-d $od);

$value ||="10";
$program||="blastn";
$view||="0";
$Blast_cpu||=30;


$od=&ABSOLUTE_DIR($od);
$query=&ABSOLUTE_DIR($query);
$cfg=&ABSOLUTE_DIR($cfg);
my $alignment=basename$query;

if ($fyes) {
	`mkdir $od/tmp ` unless (-d "$od/tmp");
	print "Formatdb fa \nperl $Bin/Blast_bin/formatdb.pl -i $fyes -odir $od/tmp\n\n";
	`perl $Bin/Blast_bin/formatdb.pl -i $fyes -odir $od/tmp `;
	my $name=basename$fyes;
	print "Blast and Process \nperl $Bin/Blast_bin/blast.pl -idir $od/tmp/$name  -i $query -p $program -od $od -m $view -e $value -cpu $Blast_cpu\n\n ";
	`perl $Bin/Blast_bin/blast.pl -cfg $cfg -idir $od/tmp/$name  -i $query -p $program -od $od -m $view -e $value -cpu $Blast_cpu `;
	 `rm -r $od/tmp`;
	 `perl $Bin/Blast_bin/get_coverage_from_alignment.pl -od $od -subject $fyes -align $od/$alignment.$program.all.xls -col c5,c7,c8 -key $alignment.coverage`;
}
if ($fno) {
	print "Blast and Process \nperl $Bin/Blast_bin/blast.pl -idir $fno  -i $query -p $program -od $od -m $view -e $value -cpu $Blast_cpu\n\n ";
	`perl $Bin/Blast_bin/blast.pl -cfg $cfg -idir $fno  -i $query -p $program -od $od -m $view -e $value -cpu $Blast_cpu `;
	#`perl $Bin/Blast_bin/get_coverage_from_alignment.pl -od $od -subject $fno -align $od/$alignment.$program.all.xls -col c5,c7,c8 -key $alignment.coverage`;
	
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

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:  blast_supply
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2014/12/09
Usage:
Options:
-od     output dir
-fa     query file
-p      Type of program
	    blastp,blastn,blastx,tblastn,tblastx
  -yes     <str>   need formatdb ,fa format file,-yes or -no needed 
  -no      <str>   need not formatdb ,formatdbed database, -yes or -no needed
-e      Comparing expectations [0-+¡Þ]
-m      Type of result  
	    0:pairwise
	    7:xml
-cpu    Blast cpu

-h      Help

USAGE
	print $usage;
	exit;
}
