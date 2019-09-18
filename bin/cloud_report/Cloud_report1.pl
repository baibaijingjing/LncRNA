#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME = time();
my $version    = "1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ( $detail_cfg, $data_cfg,$analysis_dir, $step, $data_analyzer, $customer_service_exe ,$oneStepOnly);

GetOptions(
    "cfg2:s"  => \$detail_cfg,
	"cfg:s" =>\$data_cfg,
    "dir:s"  => \$analysis_dir,
    "step:i" => \$step,
	"oneStepOnly:s"=>\$oneStepOnly,
    "PL:s"   => \$data_analyzer,
    "CSE:s"  => \$customer_service_exe,
    "help|h" => \&USAGE,
) or &USAGE;
&USAGE unless ( $detail_cfg and $analysis_dir and $data_cfg );

$detail_cfg   = abs_path($detail_cfg);
$analysis_dir = abs_path($analysis_dir);
my $odir = "$analysis_dir/Analysis_Report";
mkdir $odir unless ( -d $odir );

&log_current_time("$Script start...");

# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
$step ||= 1;

# read detail config
my %detail_cfg;
&detail_cfg_read( $detail_cfg, \%detail_cfg );

# make primary dir
my $sh_dir = "$odir/work_sh";

#my $qc_report_dir = "$odir/QC_Report";
my $web_report_dir = "$odir/Web_Report";
my $HTML = "$web_report_dir/biomarker_htmlReport";
mkdir $sh_dir unless ( -d $sh_dir );

#mkdir $qc_report_dir unless (-d $qc_report_dir);
mkdir $web_report_dir unless ( -d $web_report_dir );
my $project_id = $detail_cfg{Project_id};
my $cmd;

# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------

######################### analysis result abstract
if ( $step == 1 ) {
    $cmd = "perl $Bin/bin/Result_extract22.pl --id $analysis_dir --od $web_report_dir \n";
	$cmd .="perl $Bin/bin1/Result_extract.pl --id $analysis_dir --od $web_report_dir -cfg $data_cfg -cfg2 $detail_cfg ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Final Analysis Result and Web Report is prepared. Please check!\\n\\nResult Directory: $web_report_dir \";' |mail -s \"$project_id Final Analysis Result and Web Report\" $data_analyzer\@biomarker.com.cn " if ( $data_analyzer and $customer_service_exe );    #mail result to PL
    &step_cmd_process( $cmd, "1.1.web_report.sh", $sh_dir );
    $step++ unless($oneStepOnly);
}

######################### xml report
if ( $step == 2 ) {
    $cmd = "perl $Bin/bin1/LncRNA_xml4.pl -indir $web_report_dir -config $detail_cfg -xml $web_report_dir/configtest_raw.xml  ";
    $cmd .= "  &&  iconv -f 'GB18030' -t 'utf-8' $web_report_dir/configtest_raw.xml -o $web_report_dir/configtest_long.xml \n";
  &step_cmd_process( $cmd, "2.xml_report.sh", $sh_dir );
    $step++ unless($oneStepOnly);
}
######################### transform xml format
if ( $step == 3 ) {
    open( XML, "$web_report_dir/configtest_long.xml" ) or die $!;
    open( XML1, ">$web_report_dir/configtest_xmlconvert.xml") or die $!;
    open( XML2, ">$web_report_dir/configtest.xml") or die $!;
    while (<XML>) {
        chomp;
        #next if ( /^$/ || /^\#/ );
        $_ =~ s/$web_report_dir\///g;
        my $xml_temp1 = $_;
        print XML1 "$xml_temp1\n";
        $_ =~ s/(\w+\.KEGG\.list\.html)/$1\.cloud/; 
        #$_ =~ s/(testForDEU\.html)/$1\.cloud/;
        my $xml_temp2 = $_;
        print XML2 "$xml_temp2\n";
    }
    close XML;
    close XML1;
    close XML2;
    system "mv $web_report_dir/data.txt $web_report_dir/Template";
    system "rm $web_report_dir/configtest_long.xml";
    $step++ unless($oneStepOnly);
}

######################### xml to html 
if ( $step == 4 ) {
    mkdir $HTML unless ( -d $HTML );
    #system "$config{python} $Bin/bin_utf/htmlConvert/xml2HtmlConverter.py -i $web_report_dir/configtest_xmlconvert.xml -o $HTML ";
    $cmd="$config{python} $Bin/bin1/htmlConvert/xml2HtmlConverter.py -i $web_report_dir/configtest_xmlconvert.xml -o $HTML";
    &step_cmd_process( $cmd, "3.html_report.sh", $sh_dir );
    system "mv $web_report_dir/configtest_xmlconvert*.xml $web_report_dir/Template";
    open( Html, "$HTML/index.html" ) or die $!;
    open( Html1, ">$HTML/index_solid.html") or die $!;
    while (<Html>) {
        chomp;
        #next if ( /^$/ || /^\#/ );
=pod
        $_ =~ s/href=\"BMK_3_mRNA\/BMK_4_Alt_splice\//href1=\"BMK_3_mRNA\/BMK_4_Alt_splice\//g;
        $_ =~ s/href=\"BMK_3_mRNA\/BMK_5_SNP_Analysis\//href1=\"BMK_3_mRNA\/BMK_5_SNP_Analysis\//g;
        $_ =~ s/href=\"BMK_3_mRNA\/BMK_2_geneExpression\//href1=\"BMK_3_mRNA\/BMK_2_geneExpression\//g;
        $_ =~ s/href=\"BMK_2_LncRNA\//href1=\"BMK_2_LncRNA\//g;
        $_ =~ s/href=\"BMK_3_mRNA\/BMK_1_NewGene\//href1=\"BMK_3_mRNA\/BMK_1_NewGene\//g;
=cut
        $_ =~ s/href=\"HTML\//href1=\"HTML\//g;
        my $html_temp1 = $_;
        print Html1 "$html_temp1\n";
    }
    close Html;
    close Html1;
}

#system "rm $web_report_dir/tmp.xml";
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
sub data_cfg_read {
    my ( $cfg_file, $data_cfg ) = @_;
    my $sample_id;

    open( CFG, $cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if ( /^\s+/ or /^#/ );

        if (/^Sample/) {
            $sample_id = ( split /\s+/, $_ )[1];
        }
        if ( $_ =~ /^fq1/ || $_ =~ /^fq2/ ) {
            my $file = ( split /\s+/, $_ )[1];
            die "$file is not exist!\n" unless ( -e $file );

            $data_cfg->{rawdata}{$sample_id}{fq1} = $file if $_ =~ /^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2} = $file if $_ =~ /^fq2/;
        }
    }
    close CFG;
    &log_current_time("data config done.");
}

#############################################################################################################
sub detail_cfg_read {
    my ( $cfg_file, $detail_cfg ) = @_;

    open( CFG, $cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        s/^\s+//;s/\s+$//;s/\r$//;
        next if ( /^\s+/ or /^#/ || /^$/ );
        my ( $key, $value ) = ( split /\s+/ )[0,1];   
            $detail_cfg->{$key} = $value;
    }
    close CFG;
    &log_current_time("detail config done.");
}

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
   Contact: Simon Young <yangxh\@biomarker.com.cn>
      Date: 2014-11-13

     Usage:
            --cfg      <FILE> data config
            --cfg2       <FILE>  detail config, analysis parameters
            
            --dir       <DIR>   analysis output directory

            --step      <INT>   step to start from                  [1]
                          1     web report
                          2     xml report
                          3     en xml report

            --PL        <STR>   abbr. of Project Leader\'s name
            --CSE       <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents

   Example:
            perl $Script --cfg detail.cfg --dir Analysis/

----------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
