#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.4.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg, $detail_cfg, $odir);

GetOptions(
    "data_cfg:s" =>\$data_cfg,
    "detail_cfg:s" =>\$detail_cfg,
    "od:s"   =>\$odir,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir);

system "mkdir -p $odir" unless (-d $odir);
$odir=abs_path($odir);
$data_cfg=abs_path($data_cfg);
$detail_cfg=abs_path($detail_cfg);
my $thread ||= 10;
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
my $bwa_dir = "/share/nas2/genome/biosoft/bwa/current/bwa";
my $CIRI_dir = "/share/nas34/songmm/circRNA_process/CIRI_v1.2.pl";
#my $index = $detail_cfg{Project_key};
#my $project_id = $detail_cfg{Project_id};
my $step = 1;
######################### Data Preprocessing
if ($step==1) {
    $cmd = "perl $Bin/bin/sample_Preprocess_v2.pl -data_cfg $data_cfg -o $odir/Data_preprocess ";
    &step_cmd_process($cmd,"data_preprocess",$sh_dir);
	my $config = "$odir/Config/All_sample.config";
	my $conf;
	&creat_new_config($config,"All_Combination","$odir/Data_preprocess/All_left.fq","$odir/Data_preprocess/All_right.fq",\%Assembly_Para,$conf);
	$step = 2;
}

######################### Generate Sam file
if ($step==2) {
	$cmd = "$bwa_dir index –a bwtsw $detail_cfg{Ref_seq} && $bwa_dir mem -T 19 -t $thread $detail_cfg{Ref_seq} $odir/Data_preprocess/All_left.fq $odir/Data_preprocess/All_right.fq > $odir/All.sam";
    &step_cmd_process($cmd,"blast",$sh_dir);
	$step = 3;
}
######################### CircRNA Identify
if ($step ==3) {
	my $GFFREAD_BIN     = "/share/nas2/genome/biosoft/cufflinks/2.2.1/gffread";
	$gff = basename($detail_cfg{Ref_ann});
    $gtf = basename($detail_cfg{Ref_ann});
    $gtf =~s/\.gff3?$//i;
    $gtf.= ".gtf";
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gtf";
    $cmd = "perl $CIRI_dir -P -I $odir/All.sam -O $odir/All_CircRNA.outfile -F $detail_cfg{Ref_seq} -A $gtf -RE $detail_cfg{ciri_RE}";
    &step_cmd_process($cmd,"circRNA_identify",3,$sh_dir);
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
        next if (/^\s+/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];

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
   Contact: Simon Young <yangxh\@biomarker.com.cn> 
      Date: 2014-11-13

     Usage:
            --data_cfg      <FILE>  data config, rawdata path
            --detail_cfg      <FILE>  detail config, analysis parameters & refseq path
            --od        <DIR>   analysis output directory
            --h                 help documents

   Example:
            perl $Script --data_cfg data.cfg --detail_cfg detail.cfg --od CircRNA_identify/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
