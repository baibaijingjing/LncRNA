#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use newPerlBase;
my %config=%{readconf("$Bin/../db_file.cfg")};


#########根据你集群的情况进行修改
my $python=$config{python1};
my $parse_script="$Bin/bin/eggNOGblast.py";
my $eggNOGdb;
#my $eggNOGdb="/share/nas28/jintt/database/eggNOGdb/";
##########
my ($id,$cpu,$Blast_e,$queues,$help);
my $odir;

GetOptions(
    "help|?"=>sub { usage() },
    "id:s"=>\$id,
    "Database:s"=>\$eggNOGdb,
    "od:s"=>\$odir,
    "cpu:s"=>\$cpu,
    "Blast_e:s"=>\$Blast_e,
    "queue:s"=>\$queues,
           );
sub usage{
    print qq{
Optional:
--id      the input of nucleotide of code sequence(fasta)
--Database  the eggnog.db file,force
--od      the output directory
--cpu             cpu number,default 50
--queue                 the queue is used for qsub jobs 

-h      print the help
};
exit;
}
if (!$id||!$eggNOGdb||!$odir ||$help) {
   &usage;
}
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";


my $Q_name=(glob "$id/*.fa")[0];
$Q_name=basename $Q_name;
if ($Q_name=~/(.+)\.\d+\.fa$/) {
	$Q_name=$1;
}
else {print "Your file name is Wrong!\n";die;}


&MKDIR($odir);
$odir=&ABSOLUTE_DIR($odir);
$id=&ABSOLUTE_DIR($id);
$eggNOGdb=&ABSOLUTE_DIR($eggNOGdb);
$cpu||=50;
$queues||="middle.q";

&MKDIR("$odir/Result");
my $Result_dir="$odir/Result";
&MKDIR("$odir/eggNOG_Dir");
&MKDIR("$odir/work_sh");
&MKDIR("$odir/02.gene-annotation");
my $Tab_dir="$odir/02.gene-annotation";
&MKDIR("$odir/work_sh");
&MKDIR("$odir/work_sh/eggNOG_sh");
my $sh_dir="$odir/work_sh/eggNOG_sh";

my $blast_shell_file = "$odir/work_sh/eggNOG_sh/eggNOG1.sh";
my @subfiles;
@subfiles = glob("$id/*.fa");

#first split the  input_fasta file into small file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	my $name=basename $subfile;
	print OUT "$python $parse_script -f $subfile -e $eggNOGdb -c $config{blast}/blastx -o $odir/eggNOG_Dir &&\n";
}
close OUT;
&qsubOrDie("$blast_shell_file",$queues,$cpu,"6G");
#&Cut_shell_qsub("$blast_shell_file",$cpu,"6G",$queues);
&Check_qsub_error("$blast_shell_file");

my %Class=(
	"J" => [1,"Translation, ribosomal structure and biogenesis"],
	"A" => [2,"RNA processing and modification"],
	"K" => [3,"Transcription"],
	"L" => [4,"Replication, recombination and repair"],
	"B" => [5,"Chromatin structure and dynamics"],
	"D" => [6,"Cell cycle control, cell division, chromosome partitioning"],
	"Y" => [7,"Nuclear structure"],
	"V" => [8,"Defense mechanisms"],
	"T" => [9,"Signal transduction mechanisms"],
	"M" => [10,"Cell wall/membrane/envelope biogenesis"],
	"N" => [11,"Cell motility"],
	"Z" => [12,"Cytoskeleton"],
	"W" => [13,"Extracellular structures"],
	"U" => [14,"Intracellular trafficking, secretion, and vesicular transport"],
	"O" => [15,"Posttranslational modification, protein turnover, chaperones"],
	"C" => [16,"Energy production and conversion"],
	"G" => [17,"Carbohydrate transport and metabolism"],
	"E" => [18,"Amino acid transport and metabolism"],
	"F" => [19,"Nucleotide transport and metabolism"],
	"H" => [20,"Coenzyme transport and metabolism"],
	"I" => [21,"Lipid transport and metabolism"],
	"P" => [22,"Inorganic ion transport and metabolism"],
	"Q" => [23,"Secondary metabolites biosynthesis, transport and catabolism"],
	"R" => [24,"General function prediction only"],
	"S" => [25,"Function unknown"],
);

system "echo 'Query\tMatch\tNOG\thsp_score\tFunctional Category\tDescription'>$odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab";
system "echo 'cat $odir/eggNOG_Dir/*functions.txt >>$odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab'>$odir/work_sh/eggNOG_sh/eggNOG2.sh";
system "rm $odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab" if (-d "$odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab");
system "sh $odir/work_sh/eggNOG_sh/eggNOG2.sh";
system "cp $odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab $odir/02.gene-annotation/$Q_name.eggNOG.blast.tab";

open(OTU,"$odir/eggNOG_Dir/$Q_name.eggNOG.blast.tab");
my $t=0;
my %hash1;
my %hash2;
my %hash3;
while (<OTU>)
{
    $t++;
    chomp;
    if ($t !=1)
    {
      my @array=split;
      if (exists$hash1{$array[0]} and $hash1{$array[0]}<$array[3])
      {
        $hash2{$array[0]}=$_;
        $hash3{$array[0]}=$array[4];
      }
      if (!$hash1{$array[0]})   
      {
        $hash1{$array[0]}=$array[3];
        $hash2{$array[0]}=$_;
        $hash3{$array[0]}=$array[4];
      }  
    }
}
close(OTU);

open (OUT,">$odir/02.gene-annotation/$Q_name.eggNOG_class")||die $!;
print OUT "#Query\tMatch\tNOG\thsp_score\tFunctional Category\tDescription\tFunction class defination\n";
foreach my $key (keys %hash1)
{
    print OUT $hash2{$key},"\t$Class{$hash3{$key}}[1]\n" if  (exists $Class{$hash3{$key}});
	print OUT $hash2{$key},"\t--\n" if (!exists $Class{$hash3{$key}});
}
close (OUT);

system "cp $odir/02.gene-annotation/$Q_name.eggNOG_class $odir/Result/$Q_name.eggNOG_class.txt";
system "perl $Bin/bin/eggNOGClassDrawer.pl -i $odir/Result/$Q_name.eggNOG_class.txt -o $odir/Result/$Q_name.eggNOG.cluster"; 
################################################################################################################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
			#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";

	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
				#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
                system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
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
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR { #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if (-f $in) {
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	} elsif(-d $in) {
		chdir $in;$return=`pwd`;chomp $return;
	} else {
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}
