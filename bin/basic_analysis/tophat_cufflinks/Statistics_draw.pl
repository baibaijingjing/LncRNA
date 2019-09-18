#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME    = time();
my @Original_ARGV = @ARGV;
my $version="2.4.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ( $cfg1, $cfg2, $od, $step, $log, $oneStepOnly,$test,$in );
GetOptions(
    "help|?"        => \&USAGE,
    "cfg1:s"        => \$cfg1,
    "cfg2:s"        => \$cfg2,
	"in:s" =>\$in,
    "od:s"          => \$od,
    "s:s"           => \$step,
    "oneStepOnly:s" => \$oneStepOnly,
    "test"=>\$test,
) or &USAGE;
&USAGE unless ( $cfg1 and $cfg2 and $od );
################################

my $notename = `hostname`;
chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od = &ABSOLUTE_DIR($od);
#$in = &ABSOLUTE_DIR($in);
&MKDIR("$od/work_sh");

$step = $step || 1;

#
# log file
#
my $startTime = GetTime();
my $user      = `whoami`;
chomp $user;
my $workDir = `pwd`;
chomp $workDir;
my $task = "perl $Bin/$Script " . join( "\t", @Original_ARGV );

open( $log, ">", "$od/Statistics_draw." . time() . ".log" ) or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $cfg1\n";
print $log "detail config file:  $cfg2\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# load config file
#==================================================================

my %total_read;
my %para;
my %sample;

open( IN, "cat $cfg1 $cfg2|" ) || die "$!\n";
while (<IN>) {
    chomp;
    s/\r$//;
    s/^\s+//;
    s/\s+$//;
    next if ( /^\#/ || /^$/ );

    my @tmp = split /\s+/, $_;
    if ( $tmp[0] =~ m/Sample/ ) {
        my $fq1 = <IN>;
        chomp $fq1;
        my $fq2 = <IN>;
        chomp $fq2;
        my @fq_1 = split /\s+/, $fq1;
        $sample{ $tmp[1] }{FQ1} = $fq_1[1];
        my @fq_2 = split /\s+/, $fq2;
        $sample{ $tmp[1] }{FQ2} = $fq_2[1];
    }
    $para{ $tmp[0] } = $tmp[1];
}
close IN;

if (-f "$od/totalRead.stat.xls" ) {
    open IN, "$od/totalRead.stat.xls" || die;
    while (<IN>) {
        chomp;
        next if (/^$/);
        my @tmp = split /\s+/, $_;
        $total_read{ $tmp[0] } = $tmp[1];
    }
    close IN;
}

#==================================================================
# pipeline
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome =basename( $para{Ref_seq} );
my $genome_size;
my $gtf;
my $gff;
my @chromosome;
my $idx_prefix;


###############################
#
# step 6: Statistic bam files
#
# ######
if ( $step == 1 ) {
    print STDOUT "=== Statistic bam files  ===\n";
    print STDOUT "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
    print STDOUT "shell file: $od/work_sh/genome_bam2depth.sh\n";
    print STDOUT "shell file: $od/work_sh/genome_Checkgraph.sh\n";
    print $log "=== Statistic bam files  ===\n";
    print $log "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
    print $log "shell file: $od/work_sh/genome_bam2depth.sh\n";
    print $log "shell file: $od/work_sh/genome_Checkgraph.sh\n";
    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "30";
    #
    # # write shell
    #
    open OUT1, ">$od/work_sh/Tophat_bam_stat.sh"   || die;
    open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
    open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
    open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
    &MKDIR("$od/Map_Stat");
    my @str_type_stat;
    foreach my $sam ( sort keys %sample ) {
        push @str_type_stat, "$od/Map_Stat/$sam.type.stat";
        print OUT1 "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Tophat/$sam/accepted_hits.bam -totalRead $total_read{$sam} -od $od/Map_Stat\n";    #insert
        print OUT2 "$config{samtools} depth $od/Tophat/$sam/accepted_hits.bam >$od/Tophat/$sam/$sam.sort.bam.depth &&\n";
        print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $para{Ref_ann} -i $od/Tophat/$sam/$sam.sort.bam.depth -od $od/Map_Stat -index $sam  &&  ";    #randcheck
        print OUT3 "perl $Bin/bin/draw_total_random.pl -id $od/Map_Stat -od $od/Map_Stat  \n";
        print OUT4 "perl $Bin/bin/plotReadDensity.pl -vf $para{Memory} -q $para{Queue_type} -cpu $para{CPU} -a $od/Ref_Genome/$genome.fai -f bam -i $od/Tophat/$sam/accepted_hits.bam -o $od/Map_Stat/ -k $sam &&\n";    #2015/08/10,modify by niulg,map_stat
        print OUT4 "$config{Rscript} $Bin/bin/pie.R infile=$od/Map_Stat/$sam.type.stat outfile=$od/Map_Stat/$sam.type.png legend.col=1  value.col=2 skip=1 sep=t &&\n";    #2015/08/10,modify bu niulg,type
    }
    close OUT1;
    close OUT2;
    close OUT3;
    close OUT4;

    # #qsub
    #`sh $od/work_sh/Tophat_bam_stat.sh`;
   
    &Cut_shell_qsub( "$od/work_sh/Tophat_bam_stat.sh","$para{CPU}", "$para{Memory}", "$para{Queue_type}" );
    &Cut_shell_qsub( "$od/work_sh/genome_bam2depth.sh","$para{CPU}", "$para{Memory}", "$para{Queue_type}");
    &Cut_shell_qsub( "$od/work_sh/genome_Checkgraph.sh","$para{CPU}", "$para{Memory}", "$para{Queue_type}" );
    &Cut_shell_qsub( "$od/work_sh/plot_ReadDensity.sh","$para{CPU}", "$para{Memory}", "$para{Queue_type}" );

    #&Check_qsub_error("$od/work_sh/Tophat_bam_stat.sh");

    my $str_type_stat = join " ", @str_type_stat;    #2015/08/10,modify bu niulg
`perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png` ;                                              #2015/08/10,modify bu niulg
    print STDOUT "\n";
    print $log "\n";
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
print $log "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
#######################################################################################
close($log);

#==================================================================
# subs
#==================================================================
sub fasta {
    my $fa = shift;
    my %Fasta;
    open IN, $fa || die;
    $/ = '>';
    <IN>;
    while (<IN>) {
        chomp;
        my ( $id, $seq ) = split /\n+/, $_, 2;
        my $seq_id = ( split /\s+/, $id )[0];

        #$seq=~s/\s+//g;
        $Fasta{$seq_id} = $seq;
    }
    close IN;
    return %Fasta;

}

sub LOAD_PARA {
    my $para_file = shift;
    my $para      = shift;

    my $error_status = 0;
    open IN, $para_file || die "fail open: $para_file";
    while (<IN>) {
        chomp;
        s/^\s+//;
        s/\s+$//;
        s/\r$//;
        next if ( /^$/ or /^\#/ );
        my ( $para_key, $para_value ) = split( /\s+/, $_ );
        $para->{$para_key} = $para_value;
        if ( !defined $para_value ) {
            warn "Non-exist: $para_key\n";
            $error_status = 1;
        }
    }
    close IN;
    die "\nExit due to error Undefined parameter\n" if ($error_status);
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
        if ( $notename =~ /cluster/ ) {    ####2015-09-25
                                           #if ($notename=~/login\-0\-4/) {
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
            if ( $notename =~ /cluster/ ) {    ####2015-09-25
                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
            else {
                system	"sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
        }
    }
}

sub Check_qsub_error {                         #
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

sub GetTMR {    #
     #Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
    my $fStat = shift;
    open( IN, "<", $fStat ) or die $!;
    while (<IN>) {
        if (/^Mapped Reads\s(\d+)/) {
            close(IN);
            return $1;
        }
    }
    close(IN);
    die "Error Reads Stat file.\n";
}

sub MKDIR {    # &MKDIR($out_dir);
    my ($dir) = @_;

    #	rmdir($dir) if(-d $dir);
    mkdir($dir) if ( !-d $dir );
}

sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir = `pwd`;
    chomp($cur_dir);
    my ($in) = @_;
    my $return = "";

    if ( -f $in ) {
        my $dir  = dirname($in);
        my $file = basename($in);
        chdir $dir;
        $dir = `pwd`;
        chomp $dir;
        $return = "$dir/$file";
    }
    elsif ( -d $in ) {
        chdir $in;
        $return = `pwd`;
        chomp $return;
    }
    else {
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}

sub GetTime {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

sub USAGE {    #
    my $usage = <<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination£ºDesigned for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:				bam file  Map Stat
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	
                            
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
    print $usage;
    exit;
}
