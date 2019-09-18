#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $version    = "1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ( $detai_cfg, $analysis_dir, $step, $odir, $data_cfg, $sbs );

GetOptions(
    "cfg1:s" => \$detai_cfg,
    "idir:s" => \$analysis_dir,
    "odir:s" => \$odir,
    "step:i" => \$step,
    "sbs:s"  => \$sbs,
    "cfg2:s" => \$data_cfg,
    "help|h" => \&USAGE,
) or &USAGE;
&USAGE unless ( $detai_cfg and $analysis_dir and $odir and $data_cfg );

$detai_cfg    = abs_path($detai_cfg);
$data_cfg     = abs_path($data_cfg);
$analysis_dir = abs_path($analysis_dir);

&log_current_time("$Script start...");

# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
$step ||= 1;

# make primary dir
$odir = abs_path($odir);
mkdir $odir unless ( -d $odir );
my $web_report_dir = "$odir/Web_Report";
my $sh_dir         = "$web_report_dir/work_sh";
mkdir $web_report_dir unless ( -d $web_report_dir );
mkdir $sh_dir         unless ( -d $sh_dir );
my $HTML = "$web_report_dir/Web_Report";
mkdir $HTML unless ( -d $HTML );
my $cmd;

######################### analysis result abstract
if ( $step == 1 ) {
    $cmd = "perl $Bin/bin/Result_extract.pl -id $analysis_dir -od $web_report_dir -cfg $data_cfg  \n";
    &step_cmd_process( $cmd, "1.web_report.sh", $sh_dir );
    $step++ unless ($sbs);
}

######################### work report
if ( $step == 2 ) {
    $cmd = "cd $HTML \n";
    $cmd .= "perl $Bin/bin/LncRNA_xml.v2.pl -indir .. -config $detai_cfg -xml $web_report_dir/LncRNA_Webreport.xml \n";
    &step_cmd_process( $cmd, "2.xml.sh", $sh_dir );
    $step++ unless ($sbs);
}

if ( $step == 3 ) {
    $cmd = "cd $HTML \n";
    $cmd .= "iconv -f 'GB2312' -t 'utf-8' $web_report_dir/LncRNA_Webreport.xml -o $web_report_dir/tmp.xml \n";
#   $cmd .= "/share/nas2/genome/biosoft/Python/2.7.8/bin/python /share/nas11/yaob/research/htmlConvert/xml2HtmlConverter.py -i $web_report_dir/tmp.xml -o $HTML \n";
    $cmd .= "/share/nas2/genome/biosoft/Python/2.7.8/bin/python $Bin/htmlConvert/xml2HtmlConverter.py -i $web_report_dir/tmp.xml -o $HTML \n";
    &step_cmd_process( $cmd, "3.html.sh", $sh_dir );
}

system "rm $web_report_dir/data.txt";
system "rm $web_report_dir/LncRNA_Webreport.xml";
system "mv $web_report_dir/tmp.xml $web_report_dir/lncRNA_Bioinfor_Webreport.xml";
#system "rm $web_report_dir/LncRNA_Webreport.xml";
#system "rm -r $web_report_dir/work_sh";
#######################################################################################
my $elapse_time = ( time() - $BEGIN_TIME ) . "s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub step_cmd_process {
    my ( $cmd, $sh_name, $sh_dir ) = @_;
    my $sh_file    = "$sh_dir/$sh_name";
    my $log_file   = "$sh_file.log";
    my $flag       = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");
    &log_current_time("CMD: $cmd");

    if ( -e $sh_file ) {
        system "cat $sh_file >> $sh_file.bak";
        open( SH, ">$sh_file" ) or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }
    else {
        open( SH, ">$sh_file" ) or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("sh $sh_file > $log_file");
    if ( $flag != 0 ) {
        log_current_time("Error: command failed: $cmd");
        exit(1);
    }
    else {
        my $escaped_time = ( time() - $start_time ) . "s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
}

#############################################################################################################
sub log_current_time {

    # get parameter
    my ($info) = @_;

    # get current time with string
    my $curr_time = date_time_format( localtime( time() ) );

    # print info with time
    print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

#############################################################################################################
sub USAGE {
    my $usage = <<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Yao Bei <yaob\@biomarker.com.cn>
      Date: 2015-09-23

     Usage: (only cufflinks assembly)
		--cfg1	<FILE>	detail config(forced)
		--idir	<DIR>	analysis input directory
		--odir	<DIR>	analysis output directory
		--cfg2	<FILE>  data.cfg, (forced)
		--step		<INT>   step to start from       [1]
				1		web report
				2		xml
				3		html
		--sbs		step by step
		--h				help documents


   Example:
            perl $Script --cfg1 detail.cfg --idir Analysis/ --odir Report --cfg2 data.cfg

----------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
