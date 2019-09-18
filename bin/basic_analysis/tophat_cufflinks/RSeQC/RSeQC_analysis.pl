use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use newPerlBase;
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $version="1.0";
my ($in,$out,$genome,$step,$oneStepOnly,$q,$cpu,$m);
GetOptions(
    "i:s" =>\$in,
    "o:s"=>\$out,
    "genome:s"=>\$genome,
    "step:s"=>\$step,
    "oneStepOnly:s"=>\$oneStepOnly,
    "queue:s"=>\$q,
    "cpu:s"=>\$cpu,
    "m:s"=>\$m,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($in and $out);
#&MKDIR("$out");
mkdir $out unless -d $out;
mkdir "$out/work_sh" unless -d "$out/work_sh";
#&MKDIR("$out/work_sh");
my $work_sh="$out/work_sh";
$out=&ABSOLUTE_DIR($out);
$in=&ABSOLUTE_DIR($in);

$genome=&ABSOLUTE_DIR($genome);
my $gff_name=basename($genome);
$gff_name =~ s/\.gff$//ig;
$step= $step ||1;
$q= $q|| "general.q";
$cpu=$cpu || "20";
$m=$m || "20G";
if ($step==1) {
	#`export PATH=\$PATH:$Bin/bin/`;
	$ENV{"PATH"} = "$Bin/bin:" . $ENV{"PATH"};
    system "$config{sortBed} -i $genome | $Bin/bin/gff2bed > $out/$gff_name.bed";
    $step++ unless (defined $oneStepOnly);
}
if ($step==2) {
        open(SH,">$work_sh/bam_qc.sh") or die $!;
        my @bam=glob("$in/*");
        foreach my $dir (@bam){
                print "$dir\n";
                my $sam=basename($dir);
                my $bamfile="$dir/accepted_hits.bam";
                my $bin_dir=dirname($config{RSeQC});
                print SH "$config{RSeQC} -i $bamfile -r $out/$gff_name.bed > $out/$sam.RSeQC.txt \n" if (-f $bamfile);
                #print SH "cd $bin_dir && infer_experiment.py -i $bamfile -r $out/$gff_name.bed > $out/$sam.RSeQC.txt \n" if (-f $bamfile);
        }
        close SH;
        system "sh $config{qsub_sh} --queue $q --maxproc $cpu  --resource vf=$m --reqsub --independent $work_sh/bam_qc.sh";
        &Check_qsub_error("$work_sh/bam_qc.sh");
        $step++ unless (defined $oneStepOnly);
}
my %Stat;
open(OUT,">$out/fq_map.xls") or die $!;
print OUT "#sample_id\tread1_vs_gene\tread2_vs_gene\tundefined\n";
if ($step==3) {
    my @stat=glob("$out/*RSeQC.txt");
    
    foreach my $stat_file (@stat){
        print "$stat_file\n";
        my $file=basename($stat_file);
        $file=~/(.*?).RSeQC.txt/;
        my $name=$1;
        open(IN,$stat_file) or die $!;
        while (<IN>) {
            chomp;
            next  if (/^$/ || /^#/ || /This/);
            my @tem=split /\s+/,$_;
            if (/failed/) {
                $Stat{$name}{un}=$tem[-1]; 
            }
            if (/1\+\+/) {
                $Stat{$name}{left_fq}=$tem[-1];
            }
            if (/2\+\+/) {
                $Stat{$name}{right_fq}=$tem[-1];
            }
        }
        close IN;
        #print Dumper(\%Stat);
        
    }
}
foreach my $key (sort keys %Stat){
                print OUT "$key\t$Stat{$key}{left_fq}\t$Stat{$key}{right_fq}\t$Stat{$key}{un}\n" if (exists $Stat{$key}{un});
                print OUT "$key\t$Stat{$key}{left_fq}\t$Stat{$key}{right_fq}\t0.00\n" if (!exists $Stat{$key}{un});
                if ($Stat{$key}{right_fq} < $Stat{$key}{left_fq}) {
                    `echo "the $key rawdata lib_type is not fr_firststrand !!!!!!!" >> $out/warn_libtype.file`;
                }
                
        }
close OUT;
sub MKDIR{
        my $dir=@_;
        mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir=`pwd`;chomp($cur_dir);
    my ($in)=@_;
    my $return="";

    if(-f $in)
    {
        my $dir=dirname($in);
        my $file=basename($in);
        chdir $dir;$dir=`pwd`;chomp $dir;
        $return="$dir/$file";
    }
    elsif(-d $in)
    {
        chdir $in;$return=`pwd`;chomp $return;
    }
    else
    {
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}
sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version

     Usage:
            -i      <dir>       Tophat/
            -genome     <file>          *.gff3
            -o       <dir>           output dir 
            -queue		qsub queue ; general.q
            -cpu 			qsub maxproc ;30
            -m				qsub vf ;15G;
            -h                 help documents

   Example:
            perl Script -i Tophat/ -genome genome.gff3 -o outdir 
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
