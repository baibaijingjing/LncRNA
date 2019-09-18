#!/usr/bin/perl -w
#
my $ver="1.0";
my $title = "fusionmap";
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg,$detail_cfg,$odir,$script_cfg,$test);

GetOptions(
    "cfg1:s" =>\$data_cfg,
    "cfg2:s" =>\$detail_cfg,
    "script_cfg:s"=>\$script_cfg,
    "od:s"   =>\$odir,
    "test"=>\$test,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($data_cfg and $detail_cfg and $odir and $script_cfg);
#$script_cfg=&ABSOLUTE_DIR($script_cfg);
my %script_config=%{readconf("$script_cfg")};
system "mkdir -p $odir" unless (-d $odir);

$odir=&ABSOLUTE_DIR($odir);
$data_cfg=&ABSOLUTE_DIR($data_cfg);
$detail_cfg=&ABSOLUTE_DIR($detail_cfg);

createLog($title,$ver,$$,"$odir/log/",1);
timeLog ();
my %para;
my %sample;
open (IN,"cat $data_cfg $detail_cfg|") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);

	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/Sample/) {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;


my $bowtie_path =$script_config{bowtie};
my $fusionmap_path =$script_config{fusionmap};
my $mono_path =$script_config{mono};
my $samtools_path =$script_config{samtools};
my $queue = $script_config{que};
my $vf = $script_config{vf};
my $cpu = $script_config{cpu};
my $step;




stepStart(1,"bowtie");
# ------------------------------------------------------------------
# read configs and steps
# ------------------------------------------------------------------
my (%data_cfg, %detail_cfg);
# read data config
&data_cfg_read($data_cfg,\%data_cfg);
# read detail config
&detail_cfg_read($detail_cfg,\%detail_cfg);
my $work_sh = "$odir/work_sh";
mkdir $work_sh;
my $shfile = "$work_sh/bowite.sh" ;
open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
mkdirOrDie("$odir/bowtie_result");
foreach my $sam (sort keys %sample) {
  print SH "$bowtie_path --all -p 15 --chunkmbs 512 -S /share/nas1/backup_techserver/gaom/database/database/human/hg19/genome -1 $sample{$sam}{FQ1} -2 $sample{$sam}{FQ2} > $odir/bowtie_result/$sam\.sam && ";
  print SH "$samtools_path view -S -b $odir/bowtie_result/$sam\.sam > $odir/bowtie_result/$sam\.bam\n";
  print "$sam\t$sample{$sam}\t$sample{$sam}{FQ1}\n";
}
close SH;
qsubOrDie("$odir/work_sh/bowite.sh",$queue,$cpu,$vf);
stepTime(1,"bowite");

stepStart(2,"fuisonmap_Analysis");
$shfile = "$work_sh/fusionmap.sh" ;
open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
mkdirOrDie("$odir/fusionmap_result");
`cp -r $Bin/Fusion $odir/fusionmap_result/`;
`cp -r $Bin/ReferenceLibrary $odir/fusionmap_result/`;
foreach my $id (sort keys %sample){
  mkdirOrDie("$odir/fusionmap_result/$id");
  print SH "perl $Bin/fusion_config/FusionMapconfigCreater.pl --bamfile $odir/bowtie_result/$id\.bam --outputpath $odir/fusionmap_result/$id --outputname $id --datatype PE && ";
  print SH "$mono_path $fusionmap_path --pereport $odir/fusionmap_result Human.hg19 RefGene $odir/fusionmap_result/$id/FusionPE.config > $odir/fusionmap_result/$id/log.txt\n";
}
close SH;
qsubOrDie("$odir/work_sh/fusionmap.sh",$queue,$cpu,$vf);
stepTime(2,"fusionmap_Analysis");

stepStart(3,"circos_plot");
$shfile = "$work_sh/circos.sh";
open (SH,">$shfile") || die "Can't creat $shfile, $!\n" ;
mkdirOrDie("$odir/circos_result/png");
mkdrrOrDie("$odir/circos_result/reprot");
foreach my $id (sort keys %sample){
 mkdirOrDie("$odir/circos_result/$id");
 print SH "perl $Bin/circos/link.pl --cfg $Bin/circos/path.cfg -od $odir/circos_result/$id --fr $odir/fusionmap_result/$id/$id\.PairedEndFusionReport.txt ";
 print SH "&& cp $odir/circos_result/$id/circos\.png $odir/circos_result/png/$id\_circos\.png ";
 print SH "&& cp $odir/fusionmap_result/$id/$id\.PairedEndFusionReport.txt $odir/circos_result/reprot/$id\_FusionReport.txt\n";
 # perl circos/link.pl --cfg circos/path.cfg -od ../SRR1657561_fumap/circos/T01/ --fr ../SRR1657561_fumap/fusionmap_result/T01/T01.PairedEndFusionReport.txt
}
close SH;
qsubOrDie("$odir/work_sh/circos.sh",$queue,$cpu,$vf);
stepTime(3,"circos_plot");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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
		if ($_=~/^Basecall_stat/) {
            my $file=(split/\s+/,$_)[1];
            die "$file is not exist!\n" unless (-e $file);
            $data_cfg->{Basecall_stat}= $file;
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
            #die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' ) {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
	#	die "Kegg database is wrong ,should not be kobas data\n " if ($key eq 'Kegg' and $value!~/kobas/);
        }
        if ($key eq 'Queue_type') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'gff') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'gtf') {
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;
    die "Must choose Queue_type !\n" unless (exists $detail_cfg->{Queue_type});
    &log_current_time("detail config check done.");
}

#############################################################################################################
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
        '1' => 'Data_Assess',
        '2' => 'Basic_Analysis',
        '3' => 'Anno_Integrate',
        '4' => 'DEG_Analysis',
        '5' => 'Alitsplice_Analysis',
        '6' => 'SNP_Analysis',
        '7' => 'Gene_Structure_Optimize',
        '8' => 'Analysis_Report',
        '9' => 'More',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ($step_str =~/^[1-9]$/) {
            for my $s ($step_str..9) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            print "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
            die;
        }
    } else {
        if ($step eq join ',',(1..9)) {
            print "step-by-step: ON\n";
        } else {
            print "step-by-step: OFF\n";
        }

        for my $s (split /,/,$step_str) {
            if ($s =~/^[1-9]$/) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
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
   Program: $title
   Version: $ver

     Usage:
            --cfg1      <FILE>  data config, rawdata path
            --cfg2      <FILE>  detail config, analysis parameters & refseq path
            --od        <DIR>   analysis output directory
            --script_cfg <FILE> script.cfg,contain software path and queue ;

            --h                 help documents

   Example:
            perl $Script --cfg1 ref_trans.data.cfg --cfg2 ref_trans.detail.cfg --od Analysis/ --script_cfg script.config

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
