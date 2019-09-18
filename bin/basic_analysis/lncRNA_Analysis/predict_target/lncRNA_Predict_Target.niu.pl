#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($lnc,$lncGFF,$geneGFF,$od,$gene,$dist,$step,$q,$m,$oneStepOnly,$gene_exp);
####
my ($genome);

GetOptions(
				"help|?" =>\&USAGE,
				"lnc:s"=>\$lnc,
				"lncGFF:s"=>\$lncGFF,
				"geneGFF:s"=>\$geneGFF,
				"gene:s"=>\$gene,
				"gene_exp:s"=>\$gene_exp,
				"dist:s"=>\$dist,
				"step:s"=>\$step,
				"od:s"=>\$od,
				"q:s"=>\$q,
				"m:s"=>\$m,
				"genome:s"=>\$genome,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($od and $lnc and $lncGFF and $gene and $genome and $gene_exp);

$dist ||= 100000;
$step ||=1;
$q ||= "general.q";
$m ||= "15G";
my $LNCTAR="$Bin/bin/Predict_LncTar.pl";
my %targets;
mkdir $od unless -d $od;
$lnc=&ABSOLUTE_DIR($lnc);
$lncGFF=&ABSOLUTE_DIR($lncGFF);
$gene=&ABSOLUTE_DIR($gene);
$geneGFF=&ABSOLUTE_DIR($geneGFF);
$genome=&ABSOLUTE_DIR($genome);
$gene_exp=&ABSOLUTE_DIR($gene_exp);

mkdir "$od/work_sh" unless -d "$od/work_sh";
open OUT,">$od/work_sh/lncRNA_tar_predict.sh";
########################Cis Target Predict ###############################
my $cmd;
if ($step==1) {
	$cmd="perl $Bin/bin/get_gene_by_gff.pl -genome $genome -gff $geneGFF -o $od/gene.fa &&";
	$cmd .="perl $Bin/bin/cis_target_find.pl -inter $dist -q $lnc -t $od/gene.fa -od $od/Cis_target \n";
	#my $cmd="perl $Bin/bin/cis_target_find.s.pl -lnc $lncGFF  -gene $geneGFF  -d $dist -od $od/Cis_target\n";
#	print  OUT "$cmd\n";
#	system $cmd;
	$step++ unless ($oneStepOnly);
}

################################# Trans Predict  by LncTar based on sequence similar############################
if ($step==2) {
	$cmd .= "perl $LNCTAR -lncSeq $lnc -gene $gene -odir $od/LncTar -q $q -m $m \n";
	print OUT "$cmd\n";
#	system $cmd1;

	$step++ unless ($oneStepOnly);
}
if ($step==3){
	my $cmd2="perl $Bin/bin/Trans_target_gene_predicton.pl -in $gene_exp -od $od/Trans_target ";
	print OUT "$cmd2\n";
}
close OUT;
################################
&Cut_shell_qsub( "$od/work_sh/lncRNA_tar_predict.sh", 10, $m, $q);
&Check_qsub_error("$od/work_sh/lncRNA_tar_predict.sh");
################################


################################# Target intergation ###################################################

my $cis_target = "$od/Cis_target/Cis_target_result.txt";
my $lncTar_target = "$od/LncTar/LncTar_basepair.target";

	open CIS ,"<$cis_target" or die $!;
	while (<CIS>) {
	
		next if (/^$/ or /^\#/) ;
		chomp;
		my ($lnc_id,$target) = (split /\t/,$_);
		my @targets	= split(/;/, $target);
		push @{$targets{$lnc_id}} ,@targets;
	}	
	close CIS;
	open LNCTAR,"<$lncTar_target" or die $!;
	while(<LNCTAR>){
		next if (/^$/ or /^\#/) ;
		chomp;
		my ($lnc_id,$target) = (split /\t/,$_);
		my @targets	= split(/;/, $target);
		push @{$targets{$lnc_id}} ,@targets;
	}
	close LNCTAR;

	open OUT,">$od/novel_lncRNA_target.xls";
	print OUT "#lncRNA_ID\tGene_ID\n";
	foreach my $l (keys %targets) {
		my @targets = @{$targets{$l}};
		my %saw; 
		@saw{ @targets } = ( ); 
		my @uniq_target = sort keys %saw;
		my $genes=join(";",@uniq_target);
		print OUT "$l\t$genes\n"
	}
	close OUT;







#print "cat $od/lncRNA_position.target $od/LncTar/lncRNA_basepair.target >$od/novel_lncRNA_target.xls\n";
#`cat $od/lncRNA_position.target $od/LncTar/lncRNA_basepair.target >$od/novel_lncRNA_target.xls`;

#$/=">";
#open (IN,$t) or die $!;
#open (OUT,">$od/lnc_novel_target_gene.fa") or die $!;
#<IN>;
#while (<IN>) {
#	chomp;
#	if (exists $G{(split/\n/,$_)[0]}) {
#		print OUT ">$_";
#	}
#}
#close IN;
#close OUT;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

################################################################################################################
sub Cut_shell_qsub{
        my $sh = shift;
	my $cpu= shift;
	my $men =shift;
        my $qeue = shift;
        `sh $config{qsub_sh} --queue $qeue --reqsub -maxproc $cpu --independent $sh `;
}
###########################################################################################################
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
#############################################################################################
sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.7.13
Usage:
  Options:
  -lnc  <file>  lncRNA file,fasta format,forced;	 must be given; 
  -lncGFF       lncRNA GFF file
  -gene  <file>  gene file,fasta format,forced 	 must be given;
  -geneGFF       gene GFF file;	 must be given;
-genome		ref_fa;		must be given;
-gene_exp   	Basic_Analysis/geneExpression ; must be given;
  -dist			distance between lncRNA and genes
  -step			 step to run 
  -od <dir>   output dir,forced;	 must be given;
 
  -h         Help

  Ex:
  
	perl $Script -lnc lnc_filter_final.fa -genome genome.fa   -gene Known.longest_transcript.fa -geneGFF gene.gff3 -od Lnc_target_predict -gene_exp ./geneExpression
  
USAGE
	print $usage;
	exit;
}
