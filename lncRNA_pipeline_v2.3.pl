#!/usr/bin/perl -w
use newPerlBase;
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $Title="LncRNA";  
my $version="v2.3";  
my %config=%{readconf("$Bin/Config/lncRNA_pip.cfg")};
#my %config=%{selectconf("/share/nas1/niulg/pipline/v2.2.5/Config")}; 
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ( $data_cfg, $detail_cfg, $odir, $bdir, $step, $step_by_step,
    $data_analyzer, $customer_service_exe ,$test,$serial,$single_exon);

GetOptions(
    "cfg1:s" => \$data_cfg,
    "cfg2:s" => \$detail_cfg,
    "od:s"   => \$odir,
    "bd:s"   => \$bdir,
    "step:s" => \$step,
    "sbs"    => \$step_by_step,
    "PL:s"   => \$data_analyzer,
    "CSE:s"  => \$customer_service_exe,
    "test" =>\$test,
    "serial"=>\$serial,
    "single_exon"=>\$single_exon,
    "help|h" => \&USAGE,
) or &USAGE;
&USAGE
  unless ( $data_cfg
    and $detail_cfg
    and $odir
    and $bdir
    and $data_analyzer
    and $customer_service_exe );
    
##########create log file


###################
system "mkdir -p $odir" unless ( -d $odir );
system "mkdir -p $bdir" unless ( -d $bdir );
$odir       = abs_path($odir);
$bdir       = abs_path($bdir);
$data_cfg   = abs_path($data_cfg);
$detail_cfg = abs_path($detail_cfg);
my $notename = `hostname`; chomp $notename;

&log_current_time("$Script start...");
createLog($Title,$version,$$,"$odir",$test);
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------
my ( %data_cfg, %detail_cfg );
#-------------------------------------------------------------------
#check Phred 
#-------------------------------------------------------------------
system "perl $Bin/bin/Phred_Change/Phred_Change.pl -cfg $data_cfg -od $odir ";
# read data config
&data_cfg_read( $data_cfg, \%data_cfg );

# read detail config
&detail_cfg_read( $detail_cfg, \%detail_cfg );
 
# read steps
my %step;
$step ||= ($step_by_step) ? '1' : join ',', ( 1 .. 13 );
&steps_process( $step, $step_by_step, \%step );
#####DEU analysis
my $deu;
$deu ||= (exists $detail_cfg{'DEU'}) ? $detail_cfg{'DEU'} : '1';
if ($deu==0){
	print "DEU_Analysis off !";
	delete $step{'8'};
}
# ------------------------------------------------------------------
# main pipline
# ------------------------------------------------------------------
# cfg backup dir
my $cfg_dir = "$odir/Config";
mkdir $cfg_dir                   unless ( -d $cfg_dir );
system "cp $data_cfg $cfg_dir"   unless ( $cfg_dir eq dirname($data_cfg) );
system "cp $detail_cfg $cfg_dir" unless ( $cfg_dir eq dirname($detail_cfg) );

# work shell backup dir
my $sh_dir = "$odir/work_sh";
mkdir $sh_dir unless ( -d $sh_dir );
my $cmd;
my $log_step=1;
my $index      = $detail_cfg{Project_key};
my $project_id = $detail_cfg{Project_id};

######################### Data Assesss
if ( $step{1} ) {
	open SH1,">$odir/work_sh/S1_Data_Assess.sh";
    $cmd = "perl $Bin/bin/data_assess/rna_seq_data_assess.pl -config $data_cfg -detailfig $detail_cfg -outdir $odir/Data_Assess ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Data Assessment is done. The result is shown below:\\n\"; system \"cat $odir/Data_Assess/AllSample_GC_Q.stat\"; print \"\\nResult Directory: $odir/Data_Assess/\";' |mail -s \"$project_id Data Assessment Result\" $data_analyzer\@biomarker.com.cn ";
    print SH1 $cmd;
    close SH1;
    stepStart(1,"Data_Assess");
	&show_log2("This project will be analysed in 13 steps.")if($detail_cfg{Sep});
	&show_log2("This project will be analysed in 12 steps.")unless($detail_cfg{Sep});
	&show_log2("step_$log_step: Data assess start.");
    runOrDie("$odir/work_sh/S1_Data_Assess.sh");  
	
    stepTime(1);                                                                                                                                          
    #&step_cmd_process( $cmd, \%step, 1, $sh_dir );
	&show_log2("step_$log_step: Data assess finished.");
	$log_step++;
}

######################### Basic Analysis
if ( $step{2} ) {
	open SH2,">$odir/work_sh/S2_Basic_Analysis.sh";
    $cmd = "perl $Bin/bin/basic_analysis/basic_analysis.pl --cfg1 $data_cfg --cfg2 $detail_cfg --od $odir/Basic_Analysis ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Basic Analysis is done. The result statistic is shown below:\"; system \"perl $Bin/bin/util/stat_basic.pl -id $odir/Basic_Analysis \"; print \"Result Directory: $odir/Basic_Analysis/ \";' |mail -s \"$project_id Basic Analysis Result\" $data_analyzer\@biomarker.com.cn ";
    print SH2 $cmd;
    close SH2;
    stepStart(2,"Basic_Analysis");
	&show_log2("step_$log_step: Basic analysis including alignment, assembly start.");
    runOrDie("$odir/work_sh/S2_Basic_Analysis.sh");
    stepTime(2);                                         
    #&step_cmd_process( $cmd, \%step, 2, $sh_dir );
	&show_log2("step_$log_step: Basic analysis including alignment, assembly finished.");
	$log_step++;
}
###########################################
if (   exists $step{3}
    || exists $step{4}
    || exists $step{5}
    || exists $step{7}
    || exists $step{6}
    || exists $step{8} )
{
    open OUT, ">$sh_dir/step_3-8_merged.sh" || die;
}
if ( $step{3} ) {
	open SH3,">$odir/work_sh/S3_lncRNA_analysis.sh";
    $cmd = "perl $Bin/bin/basic_analysis/lncRNA_Analysis/Lncrna_analysis.v3.pl -gtf $odir/Basic_Analysis/Tophat_Cufflinks/Compare/Compare.gtf -od $odir/Basic_Analysis/Tophat_Cufflinks -cfg $detail_cfg  -LncExp $odir/Basic_Analysis/LncExpression   -gene_exp $odir/Basic_Analysis/geneExpression   ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id lncRNA Analysis is done. ' |mail -s \"$project_id Basic Analysis Result\" $data_analyzer\@biomarker.com.cn ";
    if (exists $detail_cfg{CircRNA} && $detail_cfg{CircRNA} eq "on") {
        $cmd.="\n perl $Bin/bin/basic_analysis/lncRNA_Analysis/Personality_analysis.pl -cfg1 $data_cfg -cfg2 $detail_cfg -od $odir/CircRNA_Analysis  &&  perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id CircRNA Analysis is done. ' |mail -s \"$project_id CircRNA Analysis Result\" $data_analyzer\@biomarker.com.cn "
    }
    
    print SH3 $cmd;
    close SH3;
    stepStart(3,"lncRNA_Analysis") if (defined $serial);
    qsubOrDie("$sh_dir/S3_lncRNA_analysis.sh",$detail_cfg{Queue_type},10,$detail_cfg{Memory}) if (defined $serial);
    qsubCheck("$sh_dir/S3_lncRNA_analysis.sh") if (defined $serial);
    #runOrDie("$odir/work_sh/S3.lncRNA_analysis.sh") if (defined $serial);
    stepTime(3) if (defined $serial);                                      
    print OUT "$cmd\n";
}
######################### annotation
my $all_unigene = "$odir/Basic_Analysis/geneExpression/final_track/All.longest_transcript.fa";
if ( $step{4} ) {
	open SH4,">$odir/work_sh/S4_Anno_Integrate.sh";
	mkdirOrDie("$odir/Anno_Integrate/New_Anno") unless (-d "$odir/Anno_Integrate/New_Anno");
	
    my $known_unigene = $detail_cfg{Known_unigene};
    my $known_anno    = $detail_cfg{Known_anno};
    my $novel_unigene = "$odir/Basic_Analysis/geneExpression/final_track/$index.newGene.longest_transcript.fa";
    `cat $known_unigene $novel_unigene >$all_unigene `;
	`cp $detail_cfg  $odir/Anno_Integrate/New_Anno/Gene_Func_Anno_Pipline.cfg`;
	 my $anno_cfg="$odir/Anno_Integrate/New_Anno/Gene_Func_Anno_Pipline.cfg";
     `echo mRNA  $novel_unigene  >>$anno_cfg `;

    $cmd = "perl $Bin/bin/annotation/Gene_Anno/Gene_Func_Anno_Pipline.pl  --cfg $anno_cfg --od $odir/Anno_Integrate/New_Anno  -queue  $detail_cfg{Queue_type} ";

    # load specified database
    $cmd.= "--nr " if (exists $detail_cfg{anno}{nr});
    $cmd.= "--swissprot " if (exists $detail_cfg{anno}{Swissprot});
    $cmd.= "--eggNOG ";
    $cmd.= "--kegg " if (exists $detail_cfg{anno}{Kegg});
    $cmd.= "--GO " if (exists $detail_cfg{anno}{nr});
    $cmd.= "--pfam ";
    $cmd.= "--kog ";
    $cmd.= "--cog " if (exists $detail_cfg{anno}{Cog});

    $cmd.= "&& perl $Bin/bin/annotation/annotation_integrate/anno_integrated2.pl -i ${index}_Unigene -gene $known_unigene,$novel_unigene -anno $known_anno,$odir/Anno_Integrate/New_Anno -od $odir/Anno_Integrate/Allgene_Anno";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Gene Function Annotation & Integration is done. The result statistic is shown below:\\n\"; system \"cat $odir/Anno_Integrate/Allgene_Anno/Result/Function_Annotation.stat.xls \"; print \"\\nResult Directory: $odir/Anno_Integrate/Allgene_Anno/Result/ \";' |mail -s \"$project_id Gene Function Annotation & Integration Result\" $data_analyzer\@biomarker.com.cn ";
    print SH4 $cmd;
    close SH4;
    stepStart(4,"Anno_Integrate") if (defined $serial);
    runOrDie("$odir/work_sh/S4_Anno_Integrate.sh") if (defined $serial);
    stepTime(4) if (defined $serial);
    print OUT "$cmd\n";
}
######################### Alternative Splicing Analysis
if ( $step{5} ) {
	open SH5,">$odir/work_sh/S5_Alitsplice_Analysis.sh";
	$cmd = "perl $Bin/bin/altsplice_analysis/v1.2/altsplice_analysis.pl -cfg $detail_cfg -cufflinks $odir/Basic_Analysis/Tophat_Cufflinks/Cufflinks -Ref_Genome $odir/Basic_Analysis/Tophat_Cufflinks/Ref_Genome -od $odir/Alitsplice_Analysis ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Alternative Splicing Analysis is done. \\n\\nResult Directory: $odir/Alitsplice_Analysis/ \";' |mail -s \"$project_id Alternative Splicing Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH5 $cmd;
    close SH5;
    stepStart(5,"Alitsplice_Analysis") if (defined $serial);
    runOrDie("$odir/work_sh/S5_Alitsplice_Analysis.sh") if (defined $serial);
    stepTime(5) if (defined $serial);
    print OUT "$cmd\n";
}

######################### SNP Analysis
if ( $step{6} ) {
	open SH6,">$odir/work_sh/S6_SNP_Analysis.sh";
    $cmd = "perl $Bin/bin/snp_analysis/v2.0/snp_analysis.pl -cfg $detail_cfg -tophat $odir/Basic_Analysis/Tophat_Cufflinks/Tophat -gff $odir/Basic_Analysis/geneExpression/final_track/$index.newGene_final.filtered.gff -od $odir/SNP_Analysis ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id SNP Analysis is done. The result statistic is shown below:\"; system \"perl $Bin/bin/util/stat_snp.pl -id $odir/SNP_Analysis \"; print \"Result Directory: $odir/SNP_Analysis/ \";' |mail -s \"$project_id SNP Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH6 $cmd;
    close SH6;
    stepStart(6,"SNP_Analysis") if (defined $serial);
    runOrDie("$odir/work_sh/S6_SNP_Analysis.sh") if (defined $serial);
    stepTime(5) if (defined $serial);
    print OUT "$cmd\n";
}

######################### Genic Structure Optimize
if ( $step{7} ) {
	open SH7,">$odir/work_sh/S7_Gene_Structure_Optimize.sh";
    $cmd = "perl $Bin/bin/gene_optimize/gene_structure_optimize.pl $odir/Basic_Analysis/geneExpression/final_track/$index.final.gff $detail_cfg{Ref_ann} $odir/Gene_Structure_Optimize/$index ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Gene Structure Optimization Analysis is done. \\n\\nResult Directory: $odir/Gene_Structure_Optimize/ \";' |mail -s \"$project_id Gene Structure Optimization Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH7 $cmd;
    close SH7;
    stepStart(7,"Gene_Structure_Optimize") if (defined $serial);
    runOrDie("$odir/work_sh/S7_Gene_Structure_Optimize.sh") if (defined $serial);
    stepTime(7) if (defined $serial);
    print OUT "$cmd\n";
}

######################## DEU analysis
if ($step{8}) {
	open SH8,">$odir/work_sh/S8_DEU_analysis.sh";
        $cmd = "perl $Bin/bin/DEU_analysis/DEU_analysis.pl -cfg $detail_cfg -tophat $odir/Basic_Analysis/Tophat_Cufflinks/Tophat -gtfdir $odir/Basic_Analysis/Tophat_Cufflinks/Ref_Genome -od $odir/DEU_analysis";
        $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id DEU_analysis is done.\\n\\nResult Directory:$odir/DEU_analysis/ \";' |mail -s \"$project_id DEU_analysis Result \" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH8 $cmd;
	close SH8;
	stepStart(8,"DEU_analysis") if (defined $serial);
	runOrDie("$odir/work_sh/S8_DEU_analysis.sh") if (defined $serial);
    stepTime(8) if (defined $serial);
    print OUT "$cmd\n";
    
}
if (!defined $serial) {
    if (   exists $step{4}
        || exists $step{5}
        || exists $step{6}
        || exists $step{7}
        || exists $step{8}
        || exists $step{3} )
    {
        if ( -e "$sh_dir/step_3-8_merged.sh" ) {
        stepStart(3,"LncRNA Analysis");
        stepStart(4,"Annotation Analysis");
        stepStart(5,"Alternative Splicing Analysis");
        stepStart(6,"SNP Analysis");
        stepStart(7,"Genic Structure Optimize");
        stepStart(8,"DEU analysis");
            &show_log2("step_3-8: step 3 to step 8 run in the same time start.");
            #&Cut_shell_qsub( "$sh_dir/step_3-8_merged.sh", 10, "$detail_cfg{Memory}", "$detail_cfg{Queue_type}" );
            #&Check_qsub_error("$sh_dir/step_3-8_merged.sh");
            qsubOrDie("$sh_dir/step_3-8_merged.sh",$detail_cfg{Queue_type},10,$detail_cfg{Memory});
            qsubCheck("$sh_dir/step_3-8_merged.sh");
            &show_log2("step_3-8: step 4 to step 9 run in the same time finished.");
        stepTime(3);
        stepTime(4);
        stepTime(5);
        stepTime(6);
        stepTime(7);
        stepTime(8);
        
        }
    }
}


######################### DEG Analysis
my ($cmd4,$cmd5);
if ( $step{9} ) {
	open SH9,">$odir/work_sh/S9_DEG_Analysis.sh";
    $cmd4 = "perl $Bin/bin/deg_analysis/v3.3/DEG_Analysis.pl -idir $odir/Basic_Analysis/geneExpression -cfg $detail_cfg -enrichment  $odir/Anno_Integrate/Allgene_Anno/Result/All_Database_annotation.xls -od $odir/DEG_Analysis ";
    $cmd4 .= "&& perl $Bin/bin/deg_analysis/v3.3/Protein_to_protein.pl -cfg $detail_cfg -uni $all_unigene  -id $odir/DEG_Analysis --clade eukaryota -od $odir/DEG_Analysis ";
    $cmd4 .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id DEG Analysis is done. The result statistic is shown below:\"; system \"perl $Bin/bin/util/stat_deg.pl -id $odir/DEG_Analysis \"; print \"Result Directory: $odir/DEG_Analysis/ \";' |mail -s \"$project_id DEG Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH9 $cmd4;
    close SH9;
    stepStart(9,"DEG_Analysis") if (defined $serial);
    runOrDie("$odir/work_sh/S9_DEG_Analysis.sh") if (defined $serial);
    stepTime(9) if (defined $serial);
    
    #print OUT "$cmd\n";
}

###########################Lnc Diff Ananlysics


if ( $step{10} ) {
	open SH10,">$odir/work_sh/S10_Lnc_Diff_Analysis.sh";
    $cmd5 = "perl $Bin/bin/lnc_diff/v3.4/DEG_Analysis.pl -idir $odir/Basic_Analysis/LncExpression -idir2 $odir/DEG_Analysis/ -lnclist $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list -lnctar $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict -cfg $detail_cfg -od $odir/Lnc_Diff_Analysis -enrichment $odir/Anno_Integrate/Allgene_Anno/Result/All_Database_annotation.xls ";
    $cmd5 .= "&& perl $Bin/bin/lnc_diff/v3.4/Protein_to_protein.lncRNA.t.pl -cfg $detail_cfg -uni $all_unigene  -id $odir/Lnc_Diff_Analysis --clade eukaryota -od $odir/Lnc_Diff_Analysis -lnc $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/novel_lncRNA_target.xls ";
    $cmd5 .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Lnc_Diff Analysis is done. The result statistic is shown below:\"; system\"perl $Bin/bin/util/stat_deg.pl -id $odir/Lnc_Diff_Analysis \"; print \"Result Directory: $odir/Lnc_Diff_Analysis/ \";' |mail -s \"$project_id Lnc_Diff_Analysis Result\" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH10 $cmd5;
    close SH10;
    stepStart(10,"Lnc_Diff_Analysis") if (defined $serial);
    runOrDie("$odir/work_sh/S10_Lnc_Diff_Analysis.sh") if (defined $serial);
    stepTime(10) if (defined $serial);
}
if (exists $step{9} || exists $step{10}) {
    open(DEG,">$odir/work_sh/step_9-10_DEG_merged.sh") or die $!;
    if (exists $step{9} && exists $step{10}){print DEG "$cmd4 && $cmd5  ";stepStart(9,"DEG Analysis");stepStart(10,"Lnc Diff Ananlysics");}
    if (exists $step{9} && !exists $step{10}){print DEG "$cmd4 ";stepStart(9,"DEG Analysis");}
    if (!exists $step{9} && exists $step{10}){print DEG "$cmd5 ";stepStart(10,"Lnc Diff Ananlysics");}
    #&Cut_shell_qsub( "$sh_dir/step_9-10_DEG_merged.sh", 10, "$detail_cfg{Memory}", "$detail_cfg{Queue_type}" );
    qsubOrDie("$sh_dir/step_9-10_DEG_merged.sh",$detail_cfg{Queue_type},10,$detail_cfg{Memory});
    qsubCheck("$sh_dir/step_9-10_DEG_merged.sh");
    #&Check_qsub_error("$sh_dir/step_9-10_DEG_merged.sh");
    
}
if ($deu==1) {$log_step=11;}else{$log_step=10;}

########################## Compare analysis of mRNA and lncRNA
if ( $step{11} ) {
	open SH11,">$odir/work_sh/S11_Compare_Analysis.sh";
    $cmd = " \n perl $Bin/bin/compare_analysis/compare_analysis.pl -m $odir/Basic_Analysis/Tophat_Cufflinks/Compare/Compare.gtf -l $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gtf -mfa $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/gene.fa -lncfa $odir/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_final.fa -e $odir/Lnc_Diff_Analysis/all_trans_fpkm.list -cfg $detail_cfg -od $odir/Compare_analysis ";
    $cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Compare_analysis is done.\\n\\nResult Directory:$odir/Compare_analysis/ \";' |mail -s \"$project_id Compare_analysis Result \" $data_analyzer\@biomarker.com.cn ";    #mail result to data analyzer
    print SH11 $cmd;
	close SH11;
	stepStart(11,"Compare_Analysis");
	runOrDie("$odir/work_sh/S11_Compare_Analysis.sh");
	&show_log2("step_$log_step: Compare analysis of mRNA and lncRNA start.");
	#&cmd_process("$odir/work_sh/S11_Compare_Analysis.sh");
	stepTime(11);
	&show_log2("step_$log_step: Compare analysis of mRNA and lncRNA finished.");
	$log_step++;
}

######################### Analysis Reports
if ( $step{12} ) {
	open SH12,">$odir/work_sh/S12_Analysis_Report.sh";

    $cmd = "perl $Bin/bin/cloud_report/Cloud_report.pl --cfg2 $detail_cfg --cfg $data_cfg --dir $odir --PL $data_analyzer --CSE $customer_service_exe ";
	$cmd .= "&& perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id cloud_report is done.\\n\\nResult Directory:$odir/Analysis_Report/Web_Report/ \";' |mail -s \"$project_id cloud_report Result \" $data_analyzer\@biomarker.com.cn ";
    print SH12 $cmd;
	close SH12;
	stepStart(12,"Analysis_Report");
	runOrDie("$odir/work_sh/S12_Analysis_Report.sh");
	&show_log2("step_$log_step: Bioinfor_report and cloud_report start.");
	stepTime(12);
	&show_log2("step_$log_step: Bioinfor_report and cloud_report finished.");
	$log_step++;
}

######################### Final_Report for biocloud: get gene symbol and gene name
if ( $step{13} ) {
	&show_log2("step_$log_step: Get gene symbol and gene name for biocloud start.");
	stepStart(13,"Final_Report_for_biocloud");
	
	open SH13,">$odir/work_sh/S13_Final_Report_for_biocloud.sh";
    $cmd="perl $Bin/bin/package/copy_Needed_Data.pl -in  $odir -od $odir/Needed_Data -cfg $detail_cfg  &&  ";
    $cmd .="perl $Bin/bin/package/package.pl -in $odir  &&  ";
    $cmd .= "perl -e 'print \"Dear $data_analyzer,\\n\\t$project_id Final_Report is done.\\n\\nResult Directory:$odir/Final_Report/ \";' |mail -s \"$project_id Final_Report Result \" $data_analyzer\@biomarker.com.cn ";
    print SH13 $cmd;
    close SH13;
    runOrDie("$odir/work_sh/S13_Final_Report_for_biocloud.sh");
    stepTime(13);
	&show_log2("step_$log_step:  For biocloud finished.");
	&show_log2("Analysis of this project is completed!");
}

totalTime();

#######################################################################################
my $elapse_time = ( time() - $BEGIN_TIME ) . "s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
sub data_cfg_read {
    &log_current_time("data config check:");
    my ( $cfg_file, $data_cfg ) = @_;
    my $sample_id;

    open( CFG, $cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
        next if ( /^\s+/ or /^#/ or /^$/);
        if (/^Qphred/) {
            $data_cfg->{Qphred} = ( split /\s+/, $_ )[1];
        }
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

    if ( defined $data_cfg->{Qphred} ) {
        print "Qphred: $data_cfg->{Qphred}\n";
    }
    else {
        $data_cfg->{Qphred} = 33;
        print "Qphred: $data_cfg->{Qphred} [ASCII encoding type of quality score of rawdata is unknown, and default is 33.]\n";
    }

    $data_cfg->{sample_num} = scalar keys %{ $data_cfg->{rawdata} };
    print "sample_number: $data_cfg->{sample_num}\n";

    for my $s ( sort keys %{ $data_cfg->{rawdata} } ) {
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
	s/^\s+//;s/\s+$//;s/\r$//;
        next if ( /^\s+/ or /^#/ or /^$/);
	s/\s$//;
        my ($key,$value) = split /\s+/,$_,2;

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

        if ($key eq 'DEU' or $key eq 'CircRNA'){
        	#die "$key: $value is not exist!\n" unless (-e $value);
        	$detail_cfg->{$key} = $value;
        }
        if ( $key =~ /^SNP_/ ) {
            $detail_cfg->{$key} = $value;
        }
        if (   $key eq 'nr'
            or $key eq 'Swissprot'
            or $key eq 'Kegg'
            or $key eq 'Pfam'
            or $key eq 'Cog'
            or $key eq 'Kog'
		or $key eq 'eggNOG' )
        {
            die "$key: $value is not exist!\n" unless ( -e $value );
            $detail_cfg->{anno}{$key} = $value;
        }
        if ( $key eq 'fold' or 'FDR' ) {
            $detail_cfg->{$key} = $value;
        }
        if ( $key eq 'Memory' or $key eq 'Queue_type' or $key eq 'CPU' ) {
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
    my ( $step_str, $step_by_step, $step_hash ) = @_;
    my %_step_ = (
        '1'  => 'Data_Assess',
        '2'  => 'Basic_Analysis',
        '3' =>  'lncRNA_Analysis',
        '4'  => 'Anno_Integrate',
        '5'  => 'Alitsplice_Analysis',
        '6'  => 'SNP_Analysis',
        '7'  => 'Gene_Structure_Optimize',
        '8' => 'DEU_analysis',
        '9'  => 'DEG_Analysis',
        '10'  => 'Lnc_Diff_Analysis',
        '11' => 'Compare_Analysis',
        '12' => 'Analysis_Report',
		'13' => 'Final_Report_for_biocloud',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ( $step_str =~ /^\d+$/ ) {
            for my $s ( $step_str .. 13 ) {
                $step_hash->{$s} = $_step_{$s};
            }
        }
        else {
            die "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
        }
    }
    else {
        print "step-by-step: OFF\n";
        for my $s ( split /,/, $step_str ) {
            if ( $s =~ /^\d+$/ ) {
                $step_hash->{$s} = $_step_{$s};
            }
            else {
                die "ERROR: illegal steps specified by --step.\n";
            }
        }
    }

    print "steps_to_run: "
      . (
        join ", ",
        ( map { sprintf "$_.$step_hash->{$_}" } sort keys %{$step_hash} )
      ) . ".\n";
    &log_current_time("step check done.\n");
}

#############################################################################################################
sub step_cmd_process {
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

    $flag = system("$cmd > $log_file");
    if ( $flag != 0 ) {
        log_current_time("Error: command failed: $cmd");
        exit(1);
    }
    else {
        my $escaped_time = ( time() - $start_time ) . "s";
        &log_current_time("step$step_n. $step_hash->{$step_n} done, escaped time: $escaped_time.\n");
    }
}

sub cmd_process {
    my ( $sh_file ) = @_;
    my $log_file   = "$sh_file.log";
    my $flag       = 0;
    my $start_time = time();
	my $sh_name=basename($sh_file);
	my $step=(split /_/,$sh_name)[0];
    if ( -e $sh_file ) {
        system "cat $sh_file >> $sh_file.bak";
        
    }
    $flag = system("sh $sh_file > $log_file");
    if ( $flag != 0 ) {
        log_current_time("Error: command failed: $sh_file");
        exit(1);
    }
    else {
        my $escaped_time = ( time() - $start_time ) . "s";
        &log_current_time("step_$step  done, escaped time: $escaped_time.\n");
    }
}
sub Cut_shell_qsub {    #Cut shell for qsub 1000 line one file
                        # &Cut_shell_qsub($shell,$cpu,$vf,$queue);
    my $shell = shift;
    my $cpu   = shift;
    my $vf    = shift;
    my $queue = shift;

    my $line = `less -S $shell |wc -l `;
    chomp $line;
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
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$odir/../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################
sub USAGE {    #
    my $usage = <<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: niulg <niulg\@biomarker.com.cn>
      Date: 2016.10.10

     Usage:
            --cfg1      <FILE>  data config, rawdata & refseq path
            --cfg2      <FILE>  detail config, analysis parameters
            --od        <DIR>   analysis output directory
            --bd        <DIR>   backup output directory

            --step      <INT>   step to start from or steps to run (split by comma) [1]
                          1     Data_Assess
                          2     Basic_Analysis
                          3     lncRNA_Analysis
                          4     Anno_Integrate
                          5     Alitsplice_Analysis
                          6     SNP_Analysis
                          7     Gene_Structure_Optimize
                          8     DEU_analysis
                           9     DEG_Analysis
                         10     Lnc_Diff_Analysis
                          11    Compare_Analysis
                          12    Analysis_Report
                          13    Final_step
            --sbs               step-by-step analysis, start from the step specified by -step
                                or not, only run the steps specified by -step
            --serial            analysis the step serial ;default Parallel ;

            --PL        <STR>   abbr. of Project Leader\'s name
            --CSE       <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents

   Example:
            perl $Script --cfg1 data.cfg --cfg2 detail.cfg --od Analysis/ --bd Backup/ -PL niulg -CSE niulg

----------------------------------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
