#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
#use Algorithm::Combinatorics qw(combinations permutations);
my $BEGIN_TIME=time();
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$cfg,$enrich,$od,$lnclist,$lnctar,$id2);
GetOptions(
				"help|?" =>\&USAGE,
				"idir:s"=>\$id,
				"idir2:s"=>\$id2,
				"cfg:s"=>\$cfg,
				"enrichment:s"=>\$enrich,
				"lnclist:s"=>\$lnclist,
				"lnctar:s"=>\$lnctar,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $cfg and $od );

mkdir $od unless -d $od;
mkdir "$od/work_sh" unless -d "$od/work_sh";
$id=&ABSOLUTE_DIR($id);
$id2=&ABSOLUTE_DIR($id2);
$cfg=&ABSOLUTE_DIR($cfg);
$od=&ABSOLUTE_DIR($od);
$lnclist=&ABSOLUTE_DIR($lnclist);
$lnctar=&ABSOLUTE_DIR($lnctar);
my $cis_tar="$lnctar/Cis_target/Cis_target_result.txt";
my $trans_tar="$lnctar/LncTar/LncTar_basepair.target";
my $Rscript=$config{Rscript};
my $notename=`hostname`;chomp $notename;
open (SH1,">$od/work_sh/select.sh") or die $!;
open (SH2,">$od/work_sh/draw.sh") or die $!;
if (defined $enrich) {
	$enrich=dirname($enrich);
	open (SH3,">$od/work_sh/enrich.sh") or die $!;
}

################# bioinfor_pipeline.log ###################
&show_log2("step_10: Difference expression genes analysis of lncRNA start.");
###########################################################

##############abstract_read_count_from_dir
print "abstract_read_count_from_dir\nperl $Bin/bin/abstract_read_count_from_dir/v1.0/abstract_read_count_from_dir.pl -indir $id -out $od/All_gene_counts.list\n\n";
system "perl $Bin/bin/abstract_read_count_from_dir/v1.0/abstract_read_count_from_dir.pl -indir $id -out $od/All_trans_counts.list";
#`perl $Bin/data_extract_by_ids.pl -idfile $id/../Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list -destfile $od/All_trans_counts.list -out $od/LncRNA_counts.list`;
system"perl $Bin/data_extract_by_ids.pl -idfile $lnclist -destfile $od/All_trans_counts.list -out $od/LncRNA_counts.list";
#`rm $od/All_gene_counts.list`;
print "Calculate Expression\nperl $Bin/bin/count_to_expression/v1.0/count_to_expression.pl -i $od/All_trans_counts.list -o $od/all_trans_fpkm.list && perl $Bin/data_extract_by_ids.pl -idfile $id/../Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list -destfile $od/all_trans_fpkm.list -out $od/All_trans_fpkm.list \n\n";
system "perl $Bin/bin/count_to_expression/v1.0/count_to_expression.pl -i $od/All_trans_counts.list -o $od/all_trans_fpkm.list ";
system "perl $Bin/data_extract_by_ids.pl -idfile $lnclist -destfile $od/all_trans_fpkm.list -out $od/LncRNA_fpkm.list";
#############相关性分析及作图
print "cor_Analysis\nperl $Bin/bin/fpkm_cor_plot/v1.0/fpkm_cor_plot.pl -i $od/LncRNA_counts.list -od $od/density\n\n";
system "perl $Bin/bin/fpkm_cor_plot/v1.0/fpkm_cor_plot.pl -i $od/LncRNA_counts.list -od $od/density";

#############密度分布图
print "fpkm_density\nperl $Bin/bin/fpkm_density_plot/v1.0/fpkm_density_plot.pl -i $od/LncRNA_counts.list -od $od/density\n\n";
system "perl $Bin/bin/fpkm_density_plot/v1.0/fpkm_density_plot.pl -i $od/LncRNA_counts.list -od $od/density";

############盒型图
print "fpkm_box\nperl $Bin/bin/fpkm_box_plot/v1.0/fpkm_box_plot.pl -i $od/LncRNA_counts.list -od $od/density\n\n";
system "perl $Bin/bin/fpkm_box_plot/v1.0/fpkm_box_plot.pl -i $od/LncRNA_counts.list -od $od/density";




#######################################################主体

############计算所有基因FDR和FC

print "Calculate All Gene FDR and FC\nperl $Bin/bin/gene_counts_to_FDR_FC/v1.0/gene_counts_to_FDR_FC.pl  -i $od/LncRNA_counts.list -cfg $cfg -od $od/\n\n";
#`perl $Bin/bin/gene_counts_to_FDR_FC/v1.0/gene_counts_to_FDR_FC.pl  -i $od/LncRNA_counts.list -cfg $cfg -od $od/`;
system "perl $Bin/bin/gene_counts_to_FDR_FC/v1.0/gene_counts_to_FDR_FC.pl  -i $od/All_trans_counts.list -cfg $cfg -od $od/";


############筛选差异基因及作图
my(%detail_cfg);
&detail_cfg_read($cfg,\%detail_cfg);
my ($FC,$FDR,$genome,$index,%kmean);
if (exists $detail_cfg{'fold'}) {
	$FC=$detail_cfg{'fold'};
	if ($FC<=0) {
		print "Fold Change Value Should Greater than zero!";die;
	}
}
if (!exists $detail_cfg{'fold'}){
	$FC=2;
}

if (exists $detail_cfg{'FDR'}) {
	$FDR=$detail_cfg{'FDR'};
}else{
	$FDR=0.01;
}
if (exists $detail_cfg{'Project_key'}){
	$index=$detail_cfg{'Project_key'};
	
}

my $de_files;
my $i=1;
mkdir "$od/All_DEG" unless -d "$od/All_DEG";
open (DAT,">$od/All_DEG/group.dat") or die $!;
foreach my $dir (glob "$od/*_vs_*") {
	next unless -d $dir;
	mkdir "$dir/DEG_Cluster" unless -d "$dir/DEG_Cluster";
	mkdir "$dir/DEG_Cluster/hierarchical" unless -d "$dir/DEG_Cluster/hierarchical";
	my $group_name=basename $dir;
	my $cmd="perl $Bin/bin/filter_by_names/v1.0/filter_by_names.l.pl -i $dir/$group_name.tmp.xls -lnc $lnclist -o $dir/$group_name.final.xls && ";	
	$cmd.="perl $Bin/bin/filter_by_names/v1.0/filter_by_names.l.pl -i $dir/$group_name.DEG_tmp.xls -lnc $lnclist -o $dir/$group_name.DEG_final.xls";

	$de_files.="$dir/$group_name.DEG_final.xls ";
	
	$cmd.="  &&  perl $Bin/bin/Enrichment/get_deg_lnc2m_regulat.pl --in $dir/$group_name.DEG_final.xls --lnctar $lnctar -k $group_name ";
	print SH1 "$cmd  &&  ";
	my ($control,$treated)=split/_vs_/,$group_name;
	print DAT "$control\tgroup$i\n";
        print DAT "$treated\tgroup$i\n";
        $i++;
#	#########表达量散点图
#	print SH2 "perl $Bin/bin/draw_cor_plot/v1.0/draw_cor_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -o $dir/$group_name.cor.png \n";
	#########表达量彩图
	print SH2 "perl $Bin/bin/draw_DEG_cor_plot/v1.0/draw_DE_cor_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.DEG_cor.png \n";
	#########MA图
	print SH2 "perl $Bin/bin/draw_FC_count_plot/v2.0/draw_FC_count_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.FC_count.png \n";
	#########火山图
	print SH2 "perl $Bin/bin/draw_FC_FDR_plot/v2.0/draw_FC_FDR_plot.pl -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -th $FC,$FDR -o $dir/$group_name.FC_FDR.png \n";
	#########聚类作图
	#print SH2 "perl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -fn 3 -i $dir/$group_name.DEG_final.xls -od $dir/DEG_Cluster \n";
	print SH2 "$Rscript $Bin/bin/draw_cluster_heatmap/heatmap2.R --infile $dir/$group_name.DEG_final.xls --outfile $dir/DEG_Cluster/hierarchical/$group_name.DEG.cluster --rowname F --height 2000 --width 1500 --color \"blue,white,red\" \n";#modified by niulg
	#print SH2 "perl $Bin/bin/draw_DEG_heatmap/v1.0/cluster_draw.pl -fn 3 -i $dir/$group_name.DEG_final.xls -od $dir/DEG_Cluster \n";
    ############得到差异lncRNA对应的靶基因列表
	#print SH3 "perl $Bin/bin/Enrichment/get_diff_lnctarget_gene.pl -i $dir/$group_name.DEG_final.xls -t $od/../Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/novel_lncRNA_target.xls -o $dir/Anno_enrichment/$group_name.target.list \n\ ";
	#########富集分析
	if (defined $enrich) {
		$enrich=&ABSOLUTE_DIR($enrich);
		###########得到差异lncRNA对应的靶基因列表
		print SH1 "perl $Bin/bin/Enrichment/get_lnc_target_gene.pl -l $dir/$group_name.DEG_final.xls -tg $trans_tar -out $dir/Trans_Anno_enrichment -k $group_name -m $id2/$group_name/$group_name.final.xls  &&";
		print SH1 "perl $Bin/bin/Enrichment/get_lnc_target_gene.pl -l $dir/$group_name.DEG_final.xls -tg $cis_tar -out $dir/Cis_Anno_enrichment -k $group_name -m $id2/$group_name/$group_name.final.xls  \n";
		
		print SH3 "perl $Bin/bin/Enrichment/select_DEG_Anno.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -o $dir/Trans_Anno_enrichment/$group_name.annotation.xls\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_GO_graph.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Trans_Anno_enrichment\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_top_GO_graph.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Trans_Anno_enrichment/\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_KEGG_graph.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Trans_Anno_enrichment\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_COG_graph.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Trans_Anno_enrichment\n\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_eggNOG_graph.pl -i $enrich -deg $dir/Trans_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Trans_Anno_enrichment\n\n";
		
		print SH3 "perl $Bin/bin/Enrichment/select_DEG_Anno.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -o $dir/Cis_Anno_enrichment/$group_name.annotation.xls\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_GO_graph.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Cis_Anno_enrichment\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_top_GO_graph.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Cis_Anno_enrichment/\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_KEGG_graph.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Cis_Anno_enrichment\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_COG_graph.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Cis_Anno_enrichment\n\n";
		print SH3 "perl $Bin/bin/Enrichment/draw_DEG_eggNOG_graph.pl -i $enrich -deg $dir/Cis_Anno_enrichment/$group_name.lncRNA_DEG_targene.xls -k $group_name -od $dir/Cis_Anno_enrichment\n\n";
			
	}
}
print "Selecting DEG by qsub:\n$od/work_sh/select.sh\n\n";
&Shell_qsub ("$od/work_sh/select.sh",$detail_cfg{'Queue_type'},$detail_cfg{'CPU'});
close SH1;
print "Drawing Graphs by qsub:\n$od/work_sh/draw.sh\n\n";
&Shell_qsub ("$od/work_sh/draw.sh",$detail_cfg{'Queue_type'},$detail_cfg{'CPU'});
close SH2;
if (defined $enrich) {
	print "Enriching Annotation by qsub:\n$od/work_sh/enrich.sh\n\n";
	&Shell_qsub ("$od/work_sh/enrich.sh",$detail_cfg{'Queue_type'},$detail_cfg{'CPU'});
	close SH3;
}



########Kmeans analysis 
if( exists $detail_cfg{Kmeans} ) {
        open (SH,">$od/work_sh/kmeans.sh") or die $!;
        my $kmean_groups = $detail_cfg{Kmeans};
        my $len = @$kmean_groups;
        for(my $i=0; $i<$len; $i++) {
                mkdir "$od/kmean_$i" unless (-d "$od/kmean_$i");
                #mkdir "$od/kmean" unless (-d "$od/kmean");
                my $cmd;
                my @tem=split /,/,$$kmean_groups[$i];
                my $num=@tem;
                if ($num==2){
                        $cmd="perl $Bin/bin/kmeans_analysis/fpkm_prepare.pl -i $od/LncRNA_fpkm.list  -o $od/kmean_$i/kmean_group_fpkm.list -order \"$$kmean_groups[$i]\" && $Rscript $Bin/bin/kmeans_analysis/kmeans_log.R  --infile $od/kmean_$i/kmean_group_fpkm.list --outdir $od/kmean_$i";
                }else {
			$cmd="perl $Bin/bin/kmeans_analysis/fpkm_prepare.pl -i $od/LncRNA_fpkm.list  -o $od/kmean_$i/kmean_group_fpkm.list -order \"$$kmean_groups[$i]\" && $Rscript $Bin/bin/kmeans_analysis/kmeans_scale.R  --infile $od/kmean_$i/kmean_group_fpkm.list --outdir $od/kmean_$i";
		}
		print SH "$cmd\n";
	}
	print "Kmeans Analysis by qsub:\n$od/work_sh/kmeans.sh\n\n";
        &Shell_qsub ("$od/work_sh/kmeans.sh",$detail_cfg{'Queue_type'},$detail_cfg{'CPU'});
        close SH;
}

############整合所有差异基因
mkdir "$od/All_DEG" unless -d "$od/All_DEG";
mkdir "$od/All_DEG/DEG_Cluster/" unless -d "$od/All_DEG/DEG_Cluster/";
mkdir "$od/All_DEG/DEG_Cluster/hierarchical" unless -d "$od/All_DEG/DEG_Cluster/hierarchical";
print "Get All DEG\nperl $Bin/bin/filter_by_names/v1.0/filter_by_names.pl -i $od/LncRNA_fpkm.list -s $de_files -o $od/All_DEG/All.DEG_final.xls \n\n";
system "perl $Bin/bin/filter_by_names/v1.0/filter_by_names.pl -i $od/LncRNA_fpkm.list -s $de_files -o $od/All_DEG/All.DEG_final.xls ";
############所有差异基因注释聚类图
print "perl $Bin/bin/draw_anno_cluster/get_deg_file.pl -i $id2/All_gene_fpkm.list -type cis -s $od -o $od/Anno_Cluster \n";
system "perl $Bin/bin/draw_anno_cluster/get_deg_file.pl -i $id2/All_gene_fpkm.list -type cis -s $od -o $od/Cis_Anno_Cluster/ ";
print "perl $Bin/bin/draw_anno_cluster/get_deg_file.pl -i $id2/All_gene_fpkm.list -type trans -s $od -o $od/Trans_Anno_Cluster/ \n";
system "perl $Bin/bin/draw_anno_cluster/get_deg_file.pl -i $id2/All_gene_fpkm.list -type trans -s $od -o $od/Trans_Anno_Cluster/ ";
##################
print  "perl $Bin/bin/draw_anno_cluster/draw_anno_cluster.pl -i $enrich -deg $od/Cis_Anno_Cluster/All_lnc_target.list -od $od/Cis_Anno_Cluster/ \n";
system "perl $Bin/bin/draw_anno_cluster/draw_anno_cluster.pl -i $enrich -deg $od/Cis_Anno_Cluster/All_lnc_target.list -od $od/Cis_Anno_Cluster/ ";
print  "perl $Bin/bin/draw_anno_cluster/draw_anno_cluster.pl -i $enrich -deg $od/Trans_Anno_Cluster/All_lnc_target.list -od $od/Trans_Anno_Cluster/ \n";
system "perl $Bin/bin/draw_anno_cluster/draw_anno_cluster.pl -i $enrich -deg $od/Trans_Anno_Cluster/All_lnc_target.list -od $od/Trans_Anno_Cluster/ ";
############所有差异基因聚类作图
#print "Cluster_draw\nperl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -i $od/All_DEG/All.DEG_final.xls -od $od/All_DEG/DEG_Cluster \n\n";
print "Cluster_draw\n$Rscript $Bin/bin/draw_cluster_heatmap/heatmap.R --infile $od/All_DEG/All.DEG_final.xls --outfile $od/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster  --width 1500 -height 2000 --rowname F --color \"blue,white,red\"\n\n";
system "$Rscript $Bin/bin/draw_cluster_heatmap/heatmap.R --infile $od/All_DEG/All.DEG_final.xls --outfile $od/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster --height 2000 --width 1500 --rowname F --color \"blue,white,red\" ";
#`perl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -i $od/All_DEG/All.DEG_final.xls -od $od/All_DEG/DEG_Cluster `;
############veen图
print "draw_diff_group_veen\nperl $Bin/bin/draw_diff_group_veen/v1.0/draw_diff_group_veen.pl -id $od -od $od \n\n";
system "perl $Bin/bin/draw_diff_group_veen/v1.0/draw_diff_group_veen.pl -id $od -od $od ";
print "draw_diff_group_veen\nperl $Bin/bin/draw_diff_group_veen/v2.0/draw_diff_group_veen.pl -id $od -od $od/Venn \n\n";
system "perl $Bin/bin/draw_diff_group_veen/v2.0/draw_diff_group_veen.pl -id $od -od $od/Venn ";



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
sub detail_cfg_read {
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/ or /^$/);
        my ($key, $value) = (split /\s+/)[0,1];

        if ($key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'fold') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'FDR') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
	if ( $key eq 'Memory' or $key eq 'Queue_type' or $key eq 'CPU' ) {
            $detail_cfg->{$key} = $value;
        }
	if ($key eq 'Kmeans'){
		#my $ke=(split/\s+/,$_)[0];
		my @str = $_=~/(\S+)/g;
		my $len = @str;
                #if( $len != 2 ){ print "$in: $line: the number of world != 2\n"; exit; }
		if (exists $detail_cfg->{$str[0]}){
			my $tmp = $detail_cfg->{$str[0]};
                        push(@$tmp, $str[1]);
			
		}else{
			my @tmp = ();
                        $tmp[0] = $str[1];
                        $detail_cfg->{$str[0]} = \@tmp;
		}

		
	}
    }
    close CFG;
    #&log_current_time("detail config done.");
}
##########################################################################################################################
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

sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu);
	my $sh = shift;
	my $qeue = shift;
	my $cpu = shift;

	if ($notename=~/cluster/)
	{
		`sh $config{qsub_sh} --queue $qeue --reqsub -maxproc $cpu --independent $sh `;
	}
	else
	{
		`sh $config{qsub_sh} --queue $qeue --reqsub --maxproc $cpu --independent $sh  `;
	}
}


sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log1("cmd")
sub show_log1()
{
    open LOG, ">$od/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$od/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.18
Usage:
  Options:
  -idir          isoExpress files dir                       		must be given;
  -idir2	mRNA DEG_analysis dir					choice;
  -lnclist	lncRNA list file ./Lnc_filter/lnc_filter_id.list		choice;
  -lnctar	lncRNA target_gene list ./Lnc_target_predict/novel_lncRNA_target.xls	choice;
  -cfg           soft parameter to DE mRNA Analysis         			must be given;

  -enrichment    All_Database_annotation.xls                 choice

  -od            Out dir                                     		must be given;

  -h                help document

USAGE
	print $usage;
	exit;
}
