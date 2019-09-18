#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($lnc,$lncGFF,$geneGFF,$od,$gene,$dist,$step,$q,$m,$oneStepOnly,$gene_ex);
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
`mkdir $od/work_sh unless -d $od/work_sh`;
open OUT,">$od/work_sh/lncRNA_tar_predict.sh";
########################Cis Target Predict ###############################
my $cmd;
if ($step==1) {
	$cmd="perl $Bin/bin/get_gene_by_gff.pl -genome $genome -gff $geneGFF -o $od/gene.fa &&";
	$cmd .="perl $Bin/bin/cis_target_find.pl -inter $dist -q $lnc -t $od/gene.fa -od $od/Cis_target\n";
	#my $cmd="perl $Bin/bin/cis_target_find.s.pl -lnc $lncGFF  -gene $geneGFF  -d $dist -od $od/Cis_target\n";
#	print  OUT "$cmd\n";
#	system $cmd;
	$step++ unless ($oneStepOnly);
}

################################# Trans Predict  by LncTar based on sequence similar############################
if ($step==2) {
	$cmd .= "perl $LNCTAR -lncSeq $lnc -gene $gene -odir $od/LncTar\n";
	print OUT "$cmd\n";
#	system $cmd1;

	$step++ unless ($oneStepOnly);
}
#if ($step==3){
#	my $cmd2="perl $Bin/Trans_target_gene_predicton.pl -in $gene_exp -od $od/Trans_target ";
#	print OUT "$cmd2\n";
#}
close OUT;
################################
&Cut_shell_qsub( "$od/work_sh/lncRNA_tar_predict.sh", 10, $m, $q);
&Check_qsub_error("$od/work_sh/lncRNA_tar_predict.sh");
################################

if ($step==3){
    my $count_sample=`cd $od/../../geneExpression && ls |grep ".geneExpression.xls" |wc -l`;#$od=**/Basic_Analysis/geneExpression
    `mkdir $od/Trans_target -p` if(!-d "$od/Trans_target");
    ################################# 选择最大内存计算节点为 WGCNA 2015/11/30 ################################################	
    my $qhost=`qhost`;
    my $mem=0;
    my $memus=0;
    my $compute_id;
    my $compute_tmp;
    my $bigger;
    foreach my $qhost_line (split /\n/,$qhost) {
        if ($qhost_line=~/compute/) {
            my @aline=(split /\s+/,$qhost_line);
            next if $aline[6]=~/-/;
            $aline[7]=~ s/G//g;
    #		print "---------------$aline[7]----------------\n";
            if ($aline[7] >= $mem) {
                $mem         = $aline[7];
    #            print "==============$mem===============\n";
                $memus       = $aline[8];
    #            $compute_tmp = $aline[0];
                $compute_id  = $aline[0];
    #            print "==============$compute_id===============\n";
    #            print "==============$compute_tmp===============\n";
                if ($aline[7]=$mem) {
                    $compute_tmp = $compute_id;
                    if ($memus=~/G/) {############## 第一种情况 已用内存单位为G
                        if ($aline[8]=~ /M/) {
                            $aline[8]=~ s/M//;
                            $memus =~ s/G//;
                            $aline[8]=$aline[8]/1024;
                            if ($aline[8] <= $memus) {
                                $bigger = $aline[0];
                                $compute_id=$aline[0];
     #                           print "==============$compute_id===============\n";
                            }else{
                                $compute_id=$compute_tmp;
                                $bigger = $compute_tmp;
    #                            print "==============$compute_id===============\n";
                            }
                        }
                        if ($aline[8]=~ /G/) {
                            $aline[8]=~ s/G//;
                            $memus=~ s/G//;
                            if ($aline[8] <= $memus) {
                                $bigger = $aline[0];
                                $compute_id=$aline[0];
    #                            print "==============$compute_id===============\n";
                            }else{
                                $compute_id=$compute_tmp;
                                $bigger = $compute_tmp;
    #                            print "==============$compute_id===============\n";
                            }
                        }
                    }
                    if ($memus=~/M/) {################  第二种情况 已用内存单位为M
                        if ($aline[8]=~/M/) {
                            $aline[8]=~ s/M//;
                            $memus=~ s/M//;
                            if ($aline[8]<$memus) {
                                $bigger = $aline[0];
                                $compute_id=$aline[0];
    #                            print "==============$compute_id===============\n";
                            }else{
                                $compute_id=$compute_tmp;
                                $bigger = $compute_tmp;
    #                            print "==============$compute_id===============\n";
                            }
                        }
                        if ($aline[8]=~/G/) {
                            $aline[8]=~ s/G//;
                            $memus=~ s/M//;
                            $memus=$memus/1024;
                            if ($aline[8]<$memus) {
                                $bigger = $aline[0];
                                $compute_id=$aline[0];
     #                           print "==============$compute_id===============\n";
                            }else{
                                $compute_id=$compute_tmp;
                                $bigger = $compute_tmp;
    #                            print "==============$compute_id===============\n";
                            }
                        }
                    }
                    $compute_tmp = $bigger;
                }
                $compute_tmp = $bigger;
            }
        }
    }

    system "ssh $compute_id";
    system "Rscript $Bin/bin/WGCNA_v1.2.R --indir $od/../../geneExpression --outdir $od/Trans_target/ --meanFPKM 0.5 -n $count_sample -f 0.2";
}


################################# Target intergation ###################################################

my $cis_target = "$od/Cis_target/Cis_target_result.txt";
my $lncTar_target = "$od/LncTar/LncTar_basepair.target";

if ($step==3) {
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
	foreach my $l (keys %targets) {
		my @targets = @{$targets{$l}};
		my %saw; 
		@saw{ @targets } = ( ); 
		my @uniq_target = sort keys %saw;
		my $genes=join(";",@uniq_target);
		print OUT "$l\t$genes\n"
	}
	close OUT;
}







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
        my $qeue = shift;
        my $cpu = shift;
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
  -lnc  <file>  lncRNA file,fasta format,forced 
  -lncGFF       lncRNA GFF file
  -gene  <file>  gene file,fasta format,forced 
  -geneGFF       gene GFF file
  -dist			distance between lncRNA and genes
  -step			 step to run 
  -od <dir>   output dir,forced
 
  -h         Help

  Ex:
  
	perl $Script -lnc lnc_filter_final.fa -lncGFF filter_final.gff   -gene gene.fa  Known.longest_transcript.fa -geneGFF gene.gff3 -od Lnc_target_predict
  
USAGE
	print $usage;
	exit;
}
