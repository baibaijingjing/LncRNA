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
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ( $gtf, $od, $step, $key, $type, $db, $log,$LncExp,$oneStepOnly,$cfg,$gene_exp);
GetOptions(
				"help|?" =>\&USAGE,
				#"fa:s"  =>\$fa,
				"gtf:s"  =>\$gtf,
				"od:s"   =>\$od,
				"cfg:s"    =>\$cfg,
				"type:s"=>\$type,
				"db:s"=>\$db,
				"LncExp:s"=>\$LncExp,
				"gene_exp:s"=>\$gene_exp,
				"step:s"=>\$step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($gtf and $cfg and $od and $LncExp and $gene_exp) ;
################################
&show_log2("step_3:lncRNA analysis start.");
#########################
my $notename = `hostname`; chomp $notename;
$cfg = &ABSOLUTE_DIR($cfg);
#$fa = &ABSOLUTE_DIR($fa);
$gtf = &ABSOLUTE_DIR($gtf);
$gene_exp=&ABSOLUTE_DIR($gene_exp);
$LncExp=&ABSOLUTE_DIR($LncExp);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");


$step = $step || 1;
my %para;
open (IN,"$cfg") || die "$!\n";
while (<IN>) {
        chomp;
        s/\r$//;s/^\s+//;s/\s+$//;
        next if (/^\#/ || /^$/);

        my @tmp=split /\s+/,$_;
        $para{$tmp[0]}=$tmp[1];
}

if (exists $para{'type'}){
	$type=$para{'type'};
}else{
	$type=$type || "ve";
}
if (exists $para{'db'}){
	$db=$para{'db'};
}else{
	$db=$db || "MAM";
}
#$type = $type || "MAM";
#$db = $db || "ve";
#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/lncRNA_identfy.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

#print $log "data config file:  $fa\n";
print $log "detail config file:  $gtf\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";


#==================================================================
&MKDIR("$od/Lnc_filter");
&MKDIR("$od/Lnc_filter/code_filter");
my $cmd;
if ($step==1) {
	print STDOUT "=== candidate lncRNA filter   ===\n";
	print STDOUT "candidate lncRNA filter shell file: $od/work_sh/lncRNA_candidate_filter.sh\n";
	print $log "=== candidate lncRNA filter ===\n";
	print $log "candidate lncRNA filter shell file: $od/work_sh/lncRNA_candidate_filter.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/lncRNA_candidate_filter.sh" || die;
	#my $cmd;
	
	$cmd="perl $Bin/lnc_predict/stat.class_code.pl -i $gtf -o $od/Lnc_filter \n";
	$cmd.= "perl $Bin/lnc_predict/lnc_merger_filter.pl  -gtf  $gtf  -out $od/Lnc_filter/merged_filter.gtf  ";
	if (exists $para{'exon'}){
		$cmd.= " -n $para{'exon'}  ";
	}
	if (exists $para{'No_sense'}){
		$cmd.="-no_sense yes ";
	}
	$cmd.= "   &&   perl $Bin/lnc_predict/gtf2fa.pl  -gtf $od/Lnc_filter/merged_filter.gtf  -fa $para{'Ref_seq'}   -o $od/Lnc_filter/merged_filter_tmp.fa  &&  ";
	if (exists $para{'fpkm'}){
		$cmd.= "perl $Bin/lnc_predict/lncRNA_fpkm_filter.s.pl --p $para{'fpkm'} --od $od/Lnc_filter --indir $LncExp  --fa $od/Lnc_filter/merged_filter_tmp.fa --out $od/Lnc_filter/merged_filter.fa  ";
	}else{
		$cmd.= "perl $Bin/lnc_predict/lncRNA_fpkm_filter.s.pl --od $od/Lnc_filter --indir $LncExp  --fa $od/Lnc_filter/merged_filter_tmp.fa --out $od/Lnc_filter/merged_filter.fa";
	}
	
	print OUT "$cmd";
	close OUT;
	&Cut_shell_qsub("$od/work_sh/lncRNA_candidate_filter.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/lncRNA_candidate_filter.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}


if ($step==2) {
	print STDOUT "===  lncRNA prediction  ===\n";
	print STDOUT "lncRNA prediction shell file: $od/work_sh/lncRNA prediction1.sh\n";
	print $log "=== lncRNA prediction ===\n";
	print $log "lncRNA prediction shell file: $od/work_sh/lncRNA prediction1.sh\n";
	open OUT,">$od/work_sh/lncRNA_prediction1.sh" || die;
	#my $cmd;
	my $lnc_fa="$od/Lnc_filter/merged_filter.fa";
	if (exists $para{'DATA'}){
                $cmd= "perl $Bin/lnc_predict/Lnc_predict.v6.pl -cfg $cfg -fa $lnc_fa -type $type -db $db -od $od/Lnc_filter/code_filter -key lnc_code_filter -pfam $para{'DATA'} && ";
        }else {
                $cmd= "perl $Bin/lnc_predict/Lnc_predict.v6.pl -cfg $cfg -fa $lnc_fa -type $type -db $db -od $od/Lnc_filter/code_filter -key lnc_code_filter && ";
        }
	#  $cmd= "perl $Bin/lnc_predict/Lnc_predict.v5.pl -fa $lnc_fa -type $type -db $db -od $od/Lnc_filter/code_filter -key lnc_code_filter && ";
	print OUT "$cmd ";
        close OUT;
        &Cut_shell_qsub("$od/work_sh/lncRNA_prediction1.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
        &Check_qsub_error("$od/work_sh/lncRNA_prediction1.sh");

        $step++ unless ($oneStepOnly) ;
        print STDOUT "\n";
        print $log "\n";
}
if ($step==3) {
        print STDOUT "===  lncRNA prediction and statistics  ===\n";
        print STDOUT "lncRNA prediction shell file: $od/work_sh/lncRNA prediction2.sh\n";
        print $log "=== lncRNA prediction ===\n";
        open OUT,">$od/work_sh/lncRNA_prediction2.sh" || die;
        #my $cmd;
        my $lnc_fa="$od/Lnc_filter/merged_filter.fa";
	$cmd= "grep -v 'NA' $od/Lnc_filter/code_filter/list.txt |cut -f 1 >$od/Lnc_filter/lnc_filter_id.list  && ";
	$cmd.= "perl $Bin/lnc_predict/abstract_gtf_seq_by_transid.pl -i $od/Lnc_filter/lnc_filter_id.list -gtf $od/Lnc_filter/merged_filter.gtf -o $od/Lnc_filter/filter_final.gtf && ";
	$cmd.= "perl $Bin/lnc_predict/gtf_to_gff.pl -i $od/Lnc_filter/filter_final.gtf -o $od/Lnc_filter/filter_final.gff && ";
	$cmd.= "perl $Bin/lnc_predict/lnc_gtf2fa.pl -fa $para{'Ref_seq'} -gtf $od/Lnc_filter/filter_final.gtf -o $od/Lnc_filter/lnc_filter_final.fa -gff $od/Lnc_filter/filter_final.gff &&  " ;
	$cmd.="perl $Bin/lnc_predict/lnc.stat.class_code.pl -i $od/Lnc_filter/filter_final.gtf -o $od/Lnc_filter  && ";
	$cmd.="perl $Bin/lnc_predict/lncRNA_draw_circos.pl -genome $od/Ref_Genome/genome_size.txt -gff $od/Lnc_filter/filter_final.gff -gtf  $od/Lnc_filter/filter_final.gtf -od $od/Lnc_filter/circos  &&  ";
	print OUT "$cmd ";
	close OUT;
	&Cut_shell_qsub("$od/work_sh/lncRNA_prediction2.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/lncRNA_prediction2.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

if ($step==4) {
	print STDOUT "=== lncRNA targetgene prediction   ===\n";
	print STDOUT "lncRNA targetgene prediction shell file: $od/work_sh/Lnc_target_predict.sh\n";
	print $log "=== lncRNA_targetgene_predict ===\n";
	print $log "lncRNA_targetgene_predict shell file: $od/work_sh/Lnc_target_predict.sh\n";

#
# write shell
#
	open OUT,">$od/work_sh/Lnc_target_predict.sh" || die;
	&MKDIR("$od/Lnc_target_predict");
	#my $cmd;
	if (exists $para{'Cis_dist'}){
		$cmd= "perl $Bin/predict_target/lncRNA_Predict_Target.niu.pl -gene_exp $gene_exp -dist $para{'Cis_dist'} -genome $para{Ref_seq} -lnc $od/Lnc_filter/lnc_filter_final.fa -lncGFF  $od/Lnc_filter/filter_final.gff  -gene  $para{Known_unigene}  -geneGFF $para{'Ref_ann'}  -od $od/Lnc_target_predict -q $para{'Queue_type'} -m $para{'Memory'} -gene_exp $gene_exp\n";
	}else{
		$cmd= "perl $Bin/predict_target/lncRNA_Predict_Target.niu.pl -gene_exp $gene_exp -genome $para{Ref_seq} -lnc $od/Lnc_filter/lnc_filter_final.fa -lncGFF  $od/Lnc_filter/filter_final.gff  -gene  $para{Known_unigene}  -geneGFF $para{'Ref_ann'}  -od $od/Lnc_target_predict -q $para{'Queue_type'} -m $para{'Memory'} -gene_exp  $gene_exp\n";
	}
	#$cmd.= "perl $Bin/predict_target/lncRNA_Predict_Target.v1.pl -q $od/lnc_filter_final.fa -t $od/Lnc_target_predict/gene.fa -od $od/Lnc_target_predict -mRNA $para{Known_unigene}  &&";
	
	#########################
	print OUT "$cmd";
	close OUT;
	&Cut_shell_qsub("$od/work_sh/Lnc_target_predict.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
		&Check_qsub_error("$od/work_sh/Lnc_target_predict.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}
if ($step==5) {
		print STDOUT "=== lncRNA and miRNA combine analysis   ===\n";
		print STDOUT "lncRNA and miRNA combine analysis shell file: $od/work_sh/lncRNA_com_miRNA.sh\n";
		print $log "===  lncRNA and miRNA combine analysis ===\n";
		print $log " lncRNA and miRNA combine analysis shell file: $od/work_sh/lncRNA_com_miRNA.sh\n";
		
		#
		# write shell
		#
		open OUT,">$od/work_sh/lncRNA_com_miRNA.sh" || die;
		&MKDIR("$od/LncRNA_com_miRNA");
		#my $cmd;
		$cmd="perl $Bin/lncmi_com_analysis/lncRNA_com_miRNA.pl -cfg $cfg -od $od/LncRNA_com_miRNA -lnc_rna $od/Lnc_filter/lnc_filter_final.fa -lnc_tar $od/Lnc_target_predict/novel_lncRNA_target.xls \n";
		#$cmd.= "perl $Bin/predict_target/lncRNA_Predict_Target.v1.pl -q $od/lnc_filter_final.fa -t $od/Lnc_target_predict/gene.fa -od $od/Lnc_target_predict -mRNA $para{Known_unigene}  &&";
		
		#########################
		print OUT "$cmd";
		close OUT;
		&Cut_shell_qsub("$od/work_sh/lncRNA_com_miRNA.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
		&Check_qsub_error("$od/work_sh/lncRNA_com_miRNA.sh");
		print STDOUT "\n";
		print $log "\n";
}
	


	

&show_log2("step_3: lncRNA analysis finished.");
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	chomp $line;
	if ($line<=1000) {
		if ($notename=~/login\-0\-4/) {
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
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
			open OUT,">>$shell.div.$div_index.sh" || die;
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
		if ($line_num!=1) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/login\-0\-4/) {
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
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

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
#	rmdir($dir) if(-d $dir);
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");

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

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
}
#################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$od/../../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
####################


sub USAGE {#
	my $usage=<<"USAGE";

Description:
	This program is a Procedure deal with lncRNA prediction 


Usage:
	-cfg	detail.cfg ,	must be given;
	-gtf	the compare gtf file ,	must be given;
	-od 	output dir ,	./Lnc_filter must be given;
	-LncExp		LncExpression dir,	must be given;
	-gene_exp 	geneExpression
	-step		
		1:candidate lncRNA filter;
		2:lncRNA prediction ;
		3:lncRNA  statistics;
		4:lncRNA targetgene prediction;
		5:lncRNA combine miRNA analysis
	-oneStepOnly		 
	-type   DataBase For CNCI predict:"ve" for vertebrate species, "pl" for plat species;
        -db     DataBase For CPC predict:"PLN,BCT,INV,MAM,PHG,PRI,ROD,SYN,UNA,VRL,VRT,ENV";
	
USAGE
	print $usage;
	exit;
}
