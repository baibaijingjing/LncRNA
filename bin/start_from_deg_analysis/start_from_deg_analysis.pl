#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $prev_analysis_report_dir, $odir, $bdir, $step, $step_by_step, $data_analyzer, $customer_service_exe);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "idir:s" =>\$prev_analysis_report_dir,
    "od:s"   =>\$odir,
    "bd:s"   => \$bdir,
    "step:s" =>\$step,
    "sbs"    =>\$step_by_step,
    "PL:s"   =>\$data_analyzer,
    "CSE:s"  =>\$customer_service_exe,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $prev_analysis_report_dir and $odir and $bdir);

system "mkdir -p $odir" unless (-d $odir);
system "mkdir -p $bdir" unless ( -d $bdir );
$odir       = abs_path($odir);
$bdir       = abs_path($bdir);
$data_cfg   = abs_path($data_cfg);
$detail_cfg = abs_path($detail_cfg);
#my $notename = `hostname`; chomp $notename;

&log_current_time("$Script start...");
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------

my (%data_cfg, %detail_cfg);

# read data config
&data_cfg_read($data_cfg,\%data_cfg);

# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);

# read steps
my %step;
$step ||= ($step_by_step) ? '4' : join ',',(4..8);
&steps_process($step,$step_by_step,\%step);

# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
# cfg backup dir
my $cfg_dir = "$odir/Config";
mkdir $cfg_dir unless (-d $cfg_dir);
system "cp $data_cfg $cfg_dir" unless ($cfg_dir eq dirname($data_cfg));
system "cp $detail_cfg $cfg_dir" unless ($cfg_dir eq dirname($detail_cfg));

# work shell backup dir
my $sh_dir = "$odir/work_sh";
mkdir $sh_dir unless (-d $sh_dir);
my $cmd;

my $index = $detail_cfg{Project_key};
my $project_id = $detail_cfg{Project_id};
my $known_unigene = $detail_cfg{Known_unigene};
my $known_anno = $detail_cfg{Known_anno};
my $novel_unigene = "$odir/Basic_Analysis/geneExpression/final_track/$index.newGene.longest_transcript.fa";

######################### Preparation
## check previous analysis result directory specified by --id
## and copy previous upstream analysis result

if (-d $prev_analysis_report_dir and (-d "$prev_analysis_report_dir/Web_Report/" or -d "$prev_analysis_report_dir/$project_id") ) {
    $prev_analysis_report_dir = abs_path($prev_analysis_report_dir);
    system "mkdir -p $odir/Analysis_Report/ " unless (-d "$odir/Analysis_Report");
#    system "cp -rf $prev_analysis_report_dir/Allgene_annoNseq $odir/Analysis_Report/ ";

    if (-d "$prev_analysis_report_dir/Web_Report/") {

        system "cp -rf $prev_analysis_report_dir/Web_Report/ $odir/Analysis_Report/ ";
    } else {
        system "cp -rf $prev_analysis_report_dir/$project_id/ $odir/Analysis_Report ";
        system "mv $odir/Analysis_Report/$project_id $odir/Analysis_Report/Web_Report ";
		print "$prev_analysis_report_dir/Web_Report/ is not exited!\n";
    }

    system "cd $odir/Analysis_Report/Web_Report/ && rm -rf ./DEG_Analysis/ ./Lnc_Diff_Analysis ./DEU_analysis ./template/configtest.xml ./template/configtest_en.xml";
	system "rm -rf $odir/Analysis_Report/Web_Report/report.html" if (-f "$odir/Analysis_Report/Web_Report/report.html");
} else {
    print STDOUT "ERROR: the directory specified by --id is illegal or may be empty. \n";
    &USAGE;
}

#my $all_unigene = "$odir/../Basic_Analysis/geneExpression/final_track/All.longest_transcript.fa";
my $all_unigene = "$odir/../Needed_Data/All.longest_transcript.fa";
if (   exists $step{4}
    || exists $step{5}
    || exists $step{6}  )
{
    open OUT, ">$sh_dir/step_4-6_merged.sh" || die;
}
######################### DEG Analysis
if ( $step{4} ) {
    $cmd = "perl $Bin/bin/deg_analysis/v3.3/DEG_Analysis.3.pl -idir $odir/Analysis_Report/Web_Report/geneExpression -cfg $detail_cfg -enrichment  $odir/../Needed_Data/Anno_Integrate/Allgene_Anno/Result/All_Database_annotation.xls -od $odir/DEG_Analysis ";
    $cmd .= "&& perl $Bin/bin/deg_analysis/v3.3/Protein_to_protein.pl -uni $all_unigene  -id $odir/DEG_Analysis --clade eukaryota -od $odir/DEG_Analysis ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id DEG Analysis is done. The result statistic is shown below:\"; system \"perl $Bin/bin/deg_analysis/v3.3/util/stat_deg.pl -id $odir/DEG_Analysis \"; print \"Result Directory: $odir/DEG_Analysis/ \";' |mail -s \"$project_id DEG Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    &step_cmd_process_1( $cmd, \%step, 4, $sh_dir );
    print OUT "$cmd\n";
}

###########################Lnc Diff Ananlysics
if ( $step{5} ) {
    $cmd = "perl $Bin/bin/lnc_diff/v3.4/DEG_Analysis.4.pl -idir $odir/../Basic_Analysis/LncExpression -idir2 $odir/DEG_Analysis/ -lnclist $odir/../Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list -lnctar $odir/../Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/novel_lncRNA_target.xls -cfg $detail_cfg -od $odir/Lnc_Diff_Analysis -enrichment $odir/../Needed_Data/Anno_Integrate/Allgene_Anno/Result/All_Database_annotation.xls ";
    $cmd .= "&& perl $Bin/bin/lnc_diff/v3.4/Protein_to_protein.lncRNA.t.pl -uni $all_unigene  -id $odir/Lnc_Diff_Analysis --clade eukaryota -od $odir/Lnc_Diff_Analysis -lnc $odir/Analysis_Report/Web_Report/LncRNA/lncRNA_target.xls ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Lnc_Diff Analysis is done. The result statistic is shown below:\"; system\"perl $Bin/bin/lnc_diff/v3.4/util/stat_deg.pl -id $odir/Lnc_Diff_Analysis \"; print \"Result Directory: $odir/Lnc_Diff_Analysis/ \";' |mail -s \"$project_id Lnc_Diff_Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    &step_cmd_process_1( $cmd, \%step, 5, $sh_dir );
    print OUT "$cmd\n";
}

######################## DEU analysis
if ( exists $detail_cfg{Sep} ) {
    if ( $step{6} ) {
        $cmd = "perl $Bin/bin/DEU_analysis/DEU_analysis.pl -cfg $detail_cfg -tophat $odir/../Needed_Data/Tophat -gtfdir $odir/../Needed_Data/Ref_Genome -od $odir/DEU_analysis";
        $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id DEU_analysis is done.\\n\\nResult Directory:$odir/DEU_analysis/ \";' |mail -s \"$project_id DEU_analysis Result \" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
        &step_cmd_process_1( $cmd, \%step, 6, $sh_dir );
        print OUT "$cmd\n";
    }
}

if (   exists $step{4}
    || exists $step{5}
	|| exists $step{6} )
{
    if ( -e "$sh_dir/step_4-6_merged.sh" ) {
        &Cut_shell_qsub( "$sh_dir/step_4-6_merged.sh", 9, "$detail_cfg{Memory}", "$detail_cfg{Queue_type}" );
        &Check_qsub_error("$sh_dir/step_4-6_merged.sh");
    }
}

######################### Analysis Reports
if ( $step{7} ) {
#	$cmd = "perl $Bin/../bioinfor_report/LncRNA_WebReport.pl --cfg1 $detail_cfg --idir $odir --odir $odir/Bioinfor_Report --cfg2 $data_cfg";
#    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Bioinfor_Report is done.\\n\\nResult Directory:$odir/Bioinfor_Report/ \";' |mail -s \"$project_id Bioinfor_Report Result \" $data_analyzer\@biomarker.com.cn \n";
    $cmd = "perl $Bin/../cloud_report/analysis_report.pl --cfg $detail_cfg --dir $odir --PL $data_analyzer --CSE $customer_service_exe ";
#    $cmd .= "&& perl $Bin/bin/cloud_report/cloud_report/lncRNA_xml_report_v1.0.pl --id $odir/Analysis_Report/Web_Report/ --cfg $detail_cfg ";
	$cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id cloud_report is done.\\n\\nResult Directory:$odir/Analysis_Report/Web_Report/ \";' |mail -s \"$project_id cloud_report Result \" $data_analyzer\@biomarker.com.cn ";
    &step_cmd_process( $cmd, \%step, 7, $sh_dir );
}

######################## Final_Report for biocloud: get gene symbol and gene name
if ( $step{8} ) {
    system "mkdir $odir/Needed_Data" unless ( -d "$odir/Needed_Data" );
    system "cp -r $odir/Config $odir/Needed_Data ";
    system "cp -r $odir/../Needed_Data/Anno_Integrate $odir/Needed_Data";
#	system "cp    $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/novel_lncRNA_target.xls $odir/Needed_Data";
    system "cp -r $odir/../Needed_Data/Tophat $odir/Needed_Data";
    system "cp -r $odir/../Needed_Data/Ref_Genome $odir/Needed_Data";
    system "cp -r $odir/../Needed_Data/All.longest_transcript.fa $odir/Needed_Data";
=pod
#    system "perl $Bin/../final_step/get_sample.pl --cfg1 $detail_cfg --cfg2 $Bin/bin/final_step/local_and_biomaRt_database.txt --od $odir";
    unless (
        -f "$odir/Analysis_Report/Web_Report/geneExpression/final_all_out.xls" )
    {
        system "cp $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls $odir/Anno_Integrate/Allgene_Anno/final_all_out.xls";
        system "cp $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls $odir/Analysis_Report/Web_Report/geneExpression/final_all_out.xls";
    }
    else {
        system "cp $odir/Analysis_Report/Web_Report/geneExpression/final_all_out.xls $odir/Anno_Integrate/Allgene_Anno";
    }
=cut
#    system "perl $Bin/../final_step/final_DEG_annotation.pl --in1 $odir/DEG_Analysis/All_DEG/All.DEG_final.xls --in2 $odir/Analysis_Report/Web_Report/geneExpression/final_all_out.xls --out $odir/Analysis_Report/Web_Report/geneExpression/final_DEG_annotation.xls";    ################### only sample & annotation
    system "cp $odir/Analysis_Report/Web_Report/geneExpression/final_DEG_annotation.xls $odir/DEG_Analysis/All_DEG";
    system "cp -r $odir/Analysis_Report/Web_Report $odir";
    $cmd = "perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Final_Report is done.\\n\\nResult Directory:$odir/Final_Report/ \";' |mail -s \"$project_id Final_Report Result \" $data_analyzer\@biomarker.com.cn ";
    &step_cmd_process( $cmd, \%step, 8, $sh_dir );


}

#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    &log_current_time("data config check:");
    my ($cfg_file, $data_cfg) = @_;
    my $sample_id;

    open (CFG, $cfg_file) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/ or /^$/);#2016/1/28
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
#           die "$file is not exist!\n" unless (-e $file);
            print "WARNING: $file is not exist!\n" unless (-e $file);

            $data_cfg->{rawdata}{$sample_id}{fq1}=$file if $_=~/^fq1/;
            $data_cfg->{rawdata}{$sample_id}{fq2}=$file if $_=~/^fq2/;
        }
    }
    close CFG;

    if (defined $data_cfg->{Qphred}) {
        print "Qphred: $data_cfg->{Qphred}\n";
    } else {
        $data_cfg->{Qphred} = 33;
        print "Qphred: $data_cfg->{Qphred} [ASCII encoding type of quality score of rawdata is unknown, and default is 33.]\n";
    }

    $data_cfg->{sample_num} = scalar keys %{$data_cfg->{rawdata}};
    print "sample_number: $data_cfg->{sample_num}\n";

    for my $s (sort keys %{$data_cfg->{rawdata}}) {
        print "${s}_fq1: $data_cfg->{rawdata}{$s}{fq1}\n${s}_fq2: $data_cfg->{rawdata}{$s}{fq2}\n";
    }
    &log_current_time("data config check done.\n");
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ( $cfg_file, $detail_cfg ) = @_;

    open( CFG, $cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if ( /^\s+/ or /^#/ or /^$/);#2016/1/28
        my ( $key, $value ) = ( split /\s+/ )[ 0, 1 ];

        if (   $key eq 'Project_name'
            or $key eq 'Customer_info'
            or $key eq 'Project_id'
            or $key eq 'Project_key' )
        {
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'Known_unigene' or $key eq 'Known_pep' ) {
            die "$key: $value is not exist!\n" unless ( -e $value );
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'Known_anno' ) {
            die "$key: $value is not illegal!\n"
              unless ( -e "$value/02.gene-annotation"
                and -e "$value/Result" );
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'Ref_seq' or $key eq 'Ref_ann' ) {
            die "$key: $value is not exist!\n" unless ( -e $value );
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'type' or $key eq 'db' ) {
            $detail_cfg->{$key} = $value;
        }

        #if ($key eq 'Com' or $key eq 'Sep'){
        #	die "$key: $value is not exist!\n" unless (-e $value);
        #	$detail_cfg->{$key} = $value;
        #}
        if ( $key =~ /^SNP_/ ) {
            $detail_cfg->{$key} = $value;
        }
        if (   $key eq 'nr'
            or $key eq 'Swissprot'
            or $key eq 'Kegg'
            or $key eq 'Pfam'
            or $key eq 'Cog'
            or $key eq 'Kog' )
        {
            die "$key: $value is not exist!\n" unless ( -e $value );
            $detail_cfg->{anno}{$key} = $value;
        }
        if ( $key eq 'fold' or 'FDR' ) {
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'Memory' or $key eq 'Queue_type' ) {
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if ( exists $detail_cfg->{$key} );
    }
    close CFG;

    &log_current_time("detail config check done.");
}


#############################################################################################################
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
#        '0' => 'Data_Assess',
#        '1' => 'Basic_Analysis',
#        '2' => 'Anno_Integrate',
#        '3' => 'DEG_Analysis',
#        '4' => 'Lnc_Diff_Analysis',
#		'5' => 'Alitsplice_Analysis',
        '4' => 'DEG_Analysis',
        '5' => 'Lnc_Diff_Analysis',
		'6' => 'DEU analysis',
		'7' => 'Analysis_Report',
        '8' => 'Final_Report for biocloud',
#		'10'=> 'Final_Report',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ($step_str =~/^\d+$/) {
            for my $s ($step_str..7) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            die "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
        }
    } else {
        print "step-by-step: OFF\n";
        for my $s (split /,/,$step_str) {
#           if ($s =~/^[0-9]$/) {
            if ($s =~/^\d+$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                die "ERROR: illegal steps specified by --step.\n";
            }
        }
    }

    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &log_current_time("step check done.\n");
}


#############################################################################################################
sub step_cmd_process {
    my ($cmd, $step_hash, $step_n, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$step_n.$step_hash->{$step_n}.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("step$step_n. $step_hash->{$step_n} start ...");
    &log_current_time("CMD: $cmd");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("$cmd > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("step$step_n. $step_hash->{$step_n} done, escaped time: $escaped_time.\n");
    }
}

##############################################################################################################
sub step_cmd_process_1 {
    my ( $cmd, $step_hash, $step_n, $sh_dir ) = @_;
    my $sh_file    = "$sh_dir/$step_n.$step_hash->{$step_n}.sh";
    my $log_file   = "$sh_file.log";
    my $flag       = 0;
    my $start_time = time();
    &log_current_time("step$step_n. $step_hash->{$step_n} start ...");
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
}
#########################################################################################################################
sub Cut_shell_qsub {    #Cut shell for qsub 1000 line one file
                        # &Cut_shell_qsub($shell,$cpu,$vf,$queue);
    my $shell = shift;
    my $cpu   = shift;
    my $vf    = shift;
    my $queue = shift;
#    my $notename;
	my $notename = `hostname`; chomp $notename;

	my $line = `less -S $shell |wc -l `;
    if ( $line <= 1000 ) {
        if ( $notename =~ /login\-0\-4/ ) {
            system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
        }
        else {
            system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
        }
    }
    if ( $line > 1000 ) {
        my @div = glob "$shell.div*";
        foreach (@div) {
            if ( -e $_ ) {
                system "rm $_";
            }
        }
        @div = ();
        my $div_index = 1;
        my $line_num  = 1;
        open IN, "$shell" || die;
        while (<IN>) {
            chomp;
            open OUT, ">>$shell.div.$div_index.sh" || die;
            if ( $line_num < 1000 ) {
                print OUT "$_\n";
                $line_num++;
            }
            else {
                print OUT "$_\n";
                $div_index++;
                $line_num = 1;
                close OUT;
            }
        }
        if ( $line_num != 1 ) {
            close OUT;
        }
        @div = glob "$shell.div*";
        foreach my $div_file (@div) {
            if ( $notename =~ /login\-0\-4/ ) {
                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
            else {
                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
        }
    }
}

sub Check_qsub_error {    #
                          # Check The qsub process if error happend
    my $sh         = shift;
    my @Check_file = glob "$sh*.qsub/*.Check";
    my @sh_file    = glob "$sh*.qsub/*.sh";

    if ( $#sh_file != $#Check_file ) {
        print "Their Some Error Happend in $sh qsub, Please Check..\n";
        die;
    }
    else {
        print "$sh qsub is Done!\n";
    }
}
#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: Simon Young <yangxh\@biomarker.com.cn> ;
            Modificated by luml <luml\@biomarker.com.cn>
      Date: 2015-11-4

     Usage:
            --cfg1      <FILE>  data config, rawdata path
            --cfg2      <FILE>  detail config, analysis parameters & refseq path
            --idir      <DIR>   input directory, previous analysis report directory
            --od        <DIR>   analysis output directory
            --bd        <DIR>   backup output directory
            --step      <INT>   step to start from or steps to run (split by comma) [6]
                          4     DEG_Analysis
                          5     Lnc Diff Analysics
                          6     DEU analysis
                          7     Analysis_Report
                          8     Final_Report
            --sbs               step-by-step analysis, start from the step specified by -step
                                or not, only run the steps specified by -step
            --PL        <STR>   abbr. of Project Leader\'s name
            --CSE       <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents
   Example:
            perl $Script --cfg1 data.cfg --cfg2 detail.cfg --idir Analysis/Analysis_Report --od Analysis/ -PL linhj -CSE linhj

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
