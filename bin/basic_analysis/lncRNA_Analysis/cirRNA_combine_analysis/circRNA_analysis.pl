#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $odir,$gtf,$samFile,$step);

GetOptions(
    "data_cfg:s" =>\$data_cfg,
    "detail_cfg:s" =>\$detail_cfg,
    "gtf:s" =>\$gtf,
    "sam:s" =>\$samFile ,
    "od:s"   =>\$odir,
	"step:s" =>\$step,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir and $gtf);

system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$data_cfg=abs_path($data_cfg);
$detail_cfg=abs_path($detail_cfg);
$gtf=abs_path($gtf);
$samFile =abs_path($samFile);
&log_current_time("$Script start...");
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------
my (%data_cfg, %detail_cfg);

# read data config
&data_cfg_read($data_cfg,\%data_cfg);

# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);

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
#my $index = $detail_cfg{Project_key};
#my $project_id = $detail_cfg{Project_id};
my $Rscript = $config{Rscript};
my $GFFREAD_BIN     = $config{GFFREAD_BIN};
$detail_cfg{Queue_type}= "general.q" unless (defined $detail_cfg{Queue_type});
$detail_cfg{Memory}= "general.q" unless (defined $detail_cfg{Memory});
$step = $step  || 1;
my $in;
my @samples;
for my $sample (sort keys %{$data_cfg{rawdata}}) 
{
	push @samples,$sample;
}
######################### circRNA identify
if ($step==1) {
	$samFile =~s/[\r\n]$//;
	$cmd = "perl $Bin/bin/circRNA_identify.v2.pl -samFile  $samFile -ref $detail_cfg{Ref_seq} -gff $gtf -od $odir/circRNA_identify && perl $Bin/bin/exact_circRNA_gene.pl $odir/circRNA_identify/All_CircRNA.xls $odir/circRNA_identify/circRNA_gene.list";
	&step_cmd_process($cmd,"identify",$sh_dir);
	$step = 2;
}
$in = glob("$odir/*/All_CircRNA.*");
######################### circRNA expression analysis
if ($step==2) {
	
	my $sample = join(",",@samples);
	my $readCount = glob("$odir/../../totalRead.stat.xls");
	$cmd = "perl $Bin/bin/circRNAexpression.pl -i $in -sample $sample -dir $odir/expression && perl  $Bin/bin/count_normalization.pl -i $odir/expression/All_gene_counts.list -dir $odir/circRNA_identify/Data_preprocess -o $odir/expression/norm_count.xls && perl $Bin/bin/circRNA_vennFile.pl -i $odir/expression/All_gene_counts.list -o $odir/expression ";
	&step_cmd_process($cmd,"expression",$sh_dir);	
	my @venn = glob("$odir/expression/*.venn.txt");
	my $venn = join (",",@venn);
	my $name;
	foreach(@venn)
	{
		my $tmp = basename($_);
		$tmp =~s/.venn.txt//;
		$name .="$tmp,";
	}
	$name =~s/,$//;
	$cmd="$Rscript $Bin/bin/R/Vennerable_v1.0.R ", $cmd.="-d $venn -n $name -o $odir/expression -p sample",`$cmd && touch $odir/expression/venn.ok` if(@venn>5 && @venn<=9 && !-e "$odir/expression/venn.ok");
	$cmd = "$Rscript $Bin/bin/R/vennDiagram_v1.0.R ", $cmd.="-d $venn -n $name -O $odir/expression -o sample",`$cmd && touch $odir/expression/venn.ok` if(@venn<=5 && !-e "$odir/expression/venn.ok");
	$step = 3;
}
######################### circRNA sequence
if ($step ==3) {
    $cmd = "perl $Bin/bin/circRNAseq.pl -ref $detail_cfg{Ref_seq} -circ $in -o $odir/circRNA.fa";
    &step_cmd_process($cmd,"circRNA_sequence",$sh_dir);
	$step = 4;
}
######################### circRNA statistics
if ($step ==4) {
        mkdir("$odir/statistics") if(!-d "$odir/statistics");
    $cmd = "$Rscript $Bin/bin/R/pie.r $in $odir/statistics/circRNA.type 1 && $Rscript $Bin/bin/R/simpleBar.r --infile $in --outfile $odir/statistics/chromosome_distrbution.png --x.col 2 --y.col 5 --x.lab \"chromosome\" --y.lab \"chromosome distrbution\" --skip 1 && $Rscript $Bin/bin/R/statistics.r $odir/expression/All_gene_counts.list $odir/statistics/statistic.xls && perl $Bin/bin/circRNAlength.pl -i $odir/circRNA_identify/All_CircRNA.xls -o $odir/statistics/forlength.list";
    &step_cmd_process($cmd,"circRNA_statistics",$sh_dir);
        if(exists $detail_cfg{circBase})
        {
                `perl $Bin/bin/newCircRNA.pl -ref $detail_cfg{circBase} -c $odir/circRNA.fa -od $odir/statistics/ -queue $detail_cfg{Queue_type}`;
                `$Rscript $Bin/bin/R/pie.r $odir/statistics/forpie.list $odir/statistics/circ_pred 0`;
        }
        $step = 5;
}
######################### compare circRNA with mRNA
if ($step ==5) {
	my $dir = "$odir/mRNA_vs_circRNA/";
	mkdir($dir) if(!-d "$dir");
	if(!exists $detail_cfg{EmRNA})
	{
		$cmd = "perl $Bin/bin/mRNAvscircRNAlocation.pl -i $in -gff $detail_cfg{Ref_ann} -o $dir/mRNA_vs_circRNAlocation.xls ";
	}
	else
	{
		my $mRNA = $detail_cfg{EmRNA};
		$cmd = "perl $Bin/bin/mRNAVScircRNACompare.pl -i $odir/expression/All_gene_counts_detail.xls -lnc $mRNA -o $dir/mRNA_vs_circRNA.xls ";
		$cmd .= "&& $Rscript $Bin/bin/R/rpkm_density_plot_func.r $dir/mRNA_vs_circRNA.xls $dir/mRNA_vs_circRNA ";
		$cmd .= "&& perl $Bin/bin/mRNAvscircRNAlocation.pl -i $in -gff $detail_cfg{Ref_ann} -o $dir/mRNA_vs_circRNAlocation.xls ";
	}
    &step_cmd_process($cmd,"circRNA_vs_mRNA",$sh_dir);
	$step = 6;
}
######################### circRNA alternative splicing analysis
if ($step ==6) {
	mkdir("$odir/overlap_alitisplice/") if(!-d "$odir/overlap_alitisplice/");
    $cmd = "perl $Bin/bin/overlap_alitisplice.pl -i $in -o $odir/overlap_alitisplice/overlap_alitisplice.xls";
    &step_cmd_process($cmd,"alternative_splicing",$sh_dir);
	$step = 7;
}
######################### circRNA Rename with gene_id
if ($step ==7) {
	mkdir("$odir/new_name/") if(!-d "$odir/new_name/");

    $cmd = "perl $Bin/bin/circRNARename.pl -i $in -o $odir/new_name/circRNA_newname.xls";
    &step_cmd_process($cmd,"circRNA_Rename",$sh_dir);
	$step = 8;
}
###############################circRNA conservation analysis
if ($step ==8) {
	if(exists $detail_cfg{phastCons}){
	mkdir("$odir/conservation/") if(!-d "$odir/conservation/");
    $cmd = "perl $Bin/bin/conservation.pl -bw $detail_cfg{phastCons} -circ $odir/circRNA_identify/All_CircRNA.xls -od $odir/conservation";
    &step_cmd_process($cmd,"conservation",$sh_dir);
	}
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
        next if (/^\s+/ or /^#/);
        if (/^Qphred/) {
            $data_cfg->{Qphred} = (split /\s+/,$_)[1];
        }
        if (/^Sample/) {
            $sample_id=(split/\s+/,$_)[1];
        }
        if ($_=~/^fq1/ || $_=~/^fq2/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);

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

sub creat_new_config {
	my $config = shift;
	my $key = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $para = shift;
	my $conf = shift;

	$conf .= "Index\t$key\n";
	$conf .= "FQ1\t$fq1\n";
	$conf .= "FQ2\t$fq2\n\n";

	foreach my $para_meter (sort keys %{$para}) {
		$para_meter =~ /^para_(.*)/;
		$conf .= "$1\t$para->{$para_meter}\n";
	}
	open OUT,">$config" || die $!;
	print OUT "$conf";
	close OUT;
}

#############################################################################################################
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
	 s/^\s+//;s/\s+$//;s/\r$//;
        next if ( /^\s+/ or /^#/ or /^$/);
        my ($key, $value) = (split /\s+/)[0,1];
		$detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_anno') {
            die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
        }

        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;

    &log_current_time("detail config check done.");
}

#############################################################################################################
sub step_cmd_process {
    my ($cmd, $step,$sh_dir) = @_;
	
    my $sh_file = "$sh_dir/$step.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
	if(-e "$sh_file.finish")
	{
		return;
	}
	else
	{
		&log_current_time("$step step start ...");
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
			&log_current_time("$step step done, escaped time: $escaped_time.\n");
			`touch $sh_file.finish`;
		}
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
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            --data_cfg      <FILE>  data config, rawdata path
            --detail_cfg      <FILE>  detail config, analysis parameters & refseq path
            --gtf                  <FILE>  genome.gtf file
            --sam               All.sam
            --od        <DIR>   analysis output directory
            --h                 help documents

   Example:
            perl $Script --data_cfg data.cfg --detail_cfg detail.cfg --od CircRNA_identify/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
