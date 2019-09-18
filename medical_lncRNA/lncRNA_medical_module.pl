use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use Switch;

use newPerlBase;
my ($conf, $outpath, $inpath, $step);
GetOptions(
        "h|?"           =>\&USAGE,
	"conf:s"	=>\$conf,
	"i:s"     	=>\$inpath,
	"o:s"		=>\$outpath,	
	"step:s"	=>\$step,	
)or &USAGE;
&USAGE unless ($conf and $inpath);

############################################################################
#
#				Main pipeline
############################################################################

$outpath ||=$inpath;
`mkdir -p $outpath`	if(!-e $outpath);

$outpath=abs_path($outpath);
$inpath=abs_path($inpath);
$conf=abs_path($conf);

my %config=&readConfig($conf);
open(SAMPLE,$config{sample})||die $!;
my @samples=();
while(<SAMPLE>){
	next if($_!~/^Sample/);
	chomp;
	push @samples,(split(/\s+/,$_))[1];
}
close(SAMPLE);
print @samples,"\n\n";
my $sh="$outpath/medical_work_sh";
`mkdir $sh`	if(!-e $sh);
$step ||=join(",",1..4);
my @steps=split(/,/,$step);

my @known_gff_species=("GRCh37","GRCh38","GRCm38","mm9","Rnor_6.0");
my $db=$config{db};
foreach my $step (sort {$a<=>$b} @steps){
	if($step==1){
		my $n=grep /^$db$/, @known_gff_species;
		next	if($n == 0);
		print "Begin distinguish the new lncRNA and known lncRNA!\n";
		open (SH,">$sh/known_lncRNA_identify.sh")||die $!;
		my $command="perl $Bin/Known_lncRNA/distinguish_lncRNA_knowornew.pl -i $inpath/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gff -d $Bin/Known_lncRNA/lncRNA_gff/$db.lncRNA.gff -o $outpath/Known_lncRNA";
		print SH "$command\n";
		close(SH);
		runOrDie($command);
		print "End distinguish!\n";	
	}

	if($step==2){
		print "Begin Annotation for DEGs!\n";
		open (SH,">$sh/anno.sh")||die $!;
		my $tmp_dir="$inpath/DEG_Analysis";
		opendir(DEG,$tmp_dir)||die $!;
		my @vs=grep{/vs/ && -d "$tmp_dir/$_"} readdir(DEG);
		closedir(DEG);
		foreach my $v (@vs){
			my $command="perl $Bin/anno/medical_anno.pl -deg $tmp_dir/$v/$v.DEG_final.xls -db $db -o $outpath/DEG_Analysis/$v";
			print SH "$command\n";
			runOrDie($command);
	
		}
		close(SH);
		print "End:Annotation for DEGs!\n";
	}
	if($step==3){
# /share/nas1/longn/project/lncRNA/BMK151121-A382/Analysis/Basic_Analysis/Tophat_Cufflinks/Tophat/
		print "Begin gene fusion analysis!\n";
		open(SH,">$sh/gene_fusion.sh")||die $!;
		my $command="perl $Bin/gene_fusion/fusionmap3.pl --indir $inpath/Basic_Analysis/Tophat_Cufflinks/Tophat --type $config{Geno_type}  --od  $outpath/Gene_Fusion --script_cfg $Bin/gene_fusion/script.cfg   --cfg2 $conf"; ##$conf should contain Pfam and Ref_seq　parameter
		print SH "$command\n";
		close(SH);
		print "$command\n";
		runOrDie ($command);
		print "End gene fusion analysis!\n";
	}
	if($step==4){
		print "Begin produce new web report based on the regular analysis and medical module!\n";
		#########################################################################################
		####################      produce new map png contain whole chromsome

		my $mappath="$inpath/Basic_Analysis/Tophat_Cufflinks/Map_Stat";
		my $newpath="$outpath/Basic_Analysis/Tophat_Cufflinks/Map_Stat";
		`mkdir -p $newpath`	if(-d $newpath);
		foreach my $s(@samples){
			&sort_by_chr("$mappath/$s.plus.txt.stat","$newpath/$s.plus.txt.stat.tmp");
			&sort_by_chr("$mappath/$s.minus.txt.stat","$newpath/$s.minus.txt.stat.tmp");
			my $replot="Rscript $Bin/tools/plot_coverage_distribution.r --input1 $newpath/$s.plus.txt.stat.tmp --input2 $newpath/$s.minus.txt.stat.tmp --output $newpath/$s.map.png --bg F -w 10000 -C $Bin/bin/chr.list -x \"Chromsome position\" -y \"Median of read deinsity(log(2))\" -t \"Genomewide distribution of read coverage\" --col \"steelblue1,palegreen\"";
			runOrDie($replot);
		}
		#########################################################################################
		###################     Begin to produce new webreport(medical)
		
		open(SH,">$sh/web_report.sh")||die $!;
		my $command="perl $Bin/web_report/web_report.pl -regularod $inpath -medicalod $outpath";
		print SH "$command";
		close(SH);
		runOrDie ($command);

		print "End web report!\n";		
	}

	die "Error: step out of range: $_\n"	if($step<1 || $step>4);
}
totalTime();


############################################################################
#
#				sub function
############################################################################

sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}

sub sort_by_chr{
	my ($input,$output)=@_;
	my $dir=dirname $output;
	`mkdir -p $dir`	if(!-d $dir);
	`grep '^[1-9]' $input|sort -n -k 1 -k 2 >$output`;
	`grep '^X' $input|sort -n -k 2 >>$output`;
	`grep '^Y' $input|sort -n -k 2 >>$output`;
}

sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-conf	<file>	input config file, forced
	-i	<path>	regular analysis output path, forced
	-o	<path>	medical path, default -i(regular analysis output path)
	-step	<int>	analysis step, default 1,2,3,4
			1. Known lncRNA identification
			2. DEG annotation (TF, COSMIC)
			3. Fusion gene
			4. final report

	-h	Help

Example:
	
USAGE
	print $usage;
	exit;
}


