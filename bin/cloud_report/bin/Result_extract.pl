#!/usr/bin/perl -w
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $ver="2.0";

############################################
my ($id,$od,$data_cfg,$detail);
GetOptions(
	"id:s"=>\$id,
	"od:s"=>\$od,
	"cfg:s"=>\$data_cfg,
	"cfg2:s"=>\$detail,
	"h|?"=>\&help,
)or &help;
&help unless ($id and $od and $data_cfg and $detail);

###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";

#############read cfg

my %ID;

open(IN,$data_cfg) or die $!;
while(<IN>){
	chomp;
	next if (/^#/ || /^$/);
	if (/^Sample/){
		#my ($bmk_id,$sampleid)=(split /\s+/,$_,3)[1,2];
        my @IDline = split("\t",$_); 
        my $IDnum = @IDline;
        if ($IDnum == 2) {
		my $bmk_id=(split /\s+/,$_)[1];
		$ID{$bmk_id} = $bmk_id;
        }
		if($IDnum == 3){
			my ($bmk_id,$sampleid)=(split /\s+/,$_,3)[1,2];
		    $ID{$bmk_id} = $sampleid;
        }
	}
}
close IN;
my %detail_cfg;
open( CFG, $detail ) or die "$!: $detail\n";
while (<CFG>) {
	chomp;
	s/^\s+//;s/\s+$//;s/\r$//;
	next if ( /^\s+/ or /^#/ or /^$/);
	s/\s$//;
	my ($key,$value) = split /\s+/,$_,2;
	$detail_cfg{$key} = $value;
}
close CFG;
my $index;
if (exists $detail_cfg{"Project_key"}){
  $index=$detail_cfg{"Project_key"};
}
###############
$id = &ABSOLUTE_DIR($id);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
###########BMK_1
system "cp $Bin/PDF/Web_Report.pdf $od/Readme.pdf";
my $QC="$od/BMK_1_rawData";
&MKDIR($QC);
my $QC_data="$QC/BMK_1_Data_Assess";
my $QC_map="$QC/BMK_2_Mapped_Statistics";
my $QC_lib="$QC/BMK_3_Library_Assessment";
&MKDIR("$QC_data");
&MKDIR("$QC_map");
&MKDIR("$QC_lib");
############BMK_3
my @dirge;
my $gene="$od/BMK_3_mRNA";
&MKDIR($gene);
my $newgene="$od/BMK_3_mRNA/BMK_1_NewGene";push(@dirge,$newgene);
my $newgene_ann="$od/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno";push(@dirge,$newgene_ann);
my $newgene_ann_a="$od/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_1_Annotation";push(@dirge,$newgene_ann_a);
my $newgene_ann_st="$od/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_2_Statistic";push(@dirge,$newgene_ann_st);
my $newgene_ann_kegg="$od/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno/BMK_3_KEGG_map";push(@dirge,$newgene_ann_kegg);

my $gene_exp="$od/BMK_3_mRNA/BMK_2_geneExpression";push(@dirge,$gene_exp);
my $gene_exp_cor="$od/BMK_3_mRNA/BMK_2_geneExpression/BMK_1_Sample_Correlation";push(@dirge,$gene_exp_cor);
my $gene_exp_stat="$od/BMK_3_mRNA/BMK_2_geneExpression/BMK_2_Expression_Statistics";push(@dirge,$gene_exp_stat);

my $gene_deg="$od/BMK_3_mRNA/BMK_6_DEG_Analysis";push(@dirge,$gene_deg);
my $gene_deg_all="$od/BMK_3_mRNA/BMK_6_DEG_Analysis/BMK_1_All_DEG";push(@dirge,$gene_deg_all);
my $gene_deg_clust="$od/BMK_3_mRNA/BMK_6_DEG_Analysis/BMK_2_Anno_Cluster";push(@dirge,$gene_deg_clust);
my $gene_deg_ppi="$od/BMK_3_mRNA/BMK_6_DEG_Analysis/BMK_3_DEG_PPI";push(@dirge,$gene_deg_ppi);

my $as="$gene/BMK_4_Alt_splice";push(@dirge,$as);
my $snp="$gene/BMK_5_SNP_Analysis";push(@dirge,$snp);
my $snp_type="$gene/BMK_5_SNP_Analysis/BMK_3_SNP_type";push(@dirge,$snp_type);
my $snp_ann="$gene/BMK_5_SNP_Analysis/BMK_2_SNP_anno";push(@dirge,$snp_ann);
my $snp_indel="$gene/BMK_5_SNP_Analysis/BMK_1_InDel_anno";push(@dirge,$snp_indel);

my $structure="$gene/BMK_3_Gene_Structure_Optimize";push(@dirge,$structure);
my $deu;
if (-d "$id/DEG_Analysis" && -d "$id/DEU_analysis"){
  $deu="$gene/BMK_7_DEU";push(@dirge,$deu);
}

foreach my $d (@dirge){
  &T($d);
  &MKDIR($d);
  
}

################BMK_2
my $LncRNA="$od/BMK_2_LncRNA";
&MKDIR($LncRNA);
my @dirlnc;
my $lnc_assem="$od/BMK_2_LncRNA/BMK_1_Assembly_Result";push (@dirlnc,$lnc_assem);

my $lnc_pre="$od/BMK_2_LncRNA/BMK_2_LncRNA_Prediction";push (@dirlnc,$lnc_pre);
my $lnc_pre_soft="$od/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result";push (@dirlnc,$lnc_pre_soft);

my $lnc_tar="$od/BMK_2_LncRNA/BMK_4_LncRNA_Target";push (@dirlnc,$lnc_tar);
my $lnc_tar_trans="$od/BMK_2_LncRNA/BMK_4_LncRNA_Target/WGCNA_Result";push (@dirlnc,$lnc_tar_trans);

my $lnc_exp="$od/BMK_2_LncRNA/BMK_3_LncRNA_Expression";push (@dirlnc,$lnc_exp);
my $lnc_exp_cor="$od/BMK_2_LncRNA/BMK_3_LncRNA_Expression/BMK_1_Sample_Correlation";push (@dirlnc,$lnc_exp_cor);
my $lnc_exp_stat="$od/BMK_2_LncRNA/BMK_3_LncRNA_Expression/BMK_2_Expression_Statistics";push (@dirlnc,$lnc_exp_stat);
my ($lnc_deg_all, $lnc_deg,$lnc_deg_ppi,$lnc_deg_cluster,$lnc_deg_cis,$lnc_deg_trans);
if (-d "$id/Lnc_Diff_Analysis"){
	 $lnc_deg="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis";push (@dirlnc,$lnc_deg);
	$lnc_deg_all="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis/BMK_1_All_DEG";push (@dirlnc,$lnc_deg_all);
	$lnc_deg_ppi="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis/BMK_2_DEG_PPI";push (@dirlnc,$lnc_deg_ppi);
	$lnc_deg_cluster="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis/BMK_3_Anno_Cluster";push (@dirlnc,$lnc_deg_cluster);
	$lnc_deg_cis="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis/BMK_3_Anno_Cluster/BMK_1_Cis_Cluster";push (@dirlnc,$lnc_deg_cis);
	$lnc_deg_trans="$od/BMK_2_LncRNA/BMK_6_DEG_Analysis/BMK_3_Anno_Cluster/BMK_2_Trans_Cluster";push (@dirlnc,$lnc_deg_trans);
}

my $lnc_vs="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA";push (@dirlnc,$lnc_vs);
my $lnc_vs_len="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_1_Length_Compare";push (@dirlnc,$lnc_vs_len);
my $lnc_vs_exon="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_2_Exon_Compare";push (@dirlnc,$lnc_vs_exon);
my $lnc_vs_orf="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_3_ORF_Compare";push (@dirlnc,$lnc_vs_orf);
my $lnc_vs_exp="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_4_Expression_Compare";push (@dirlnc,$lnc_vs_exp);
my $lnc_vs_iso="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_5_Isoform_Compare";push (@dirlnc,$lnc_vs_iso);
my $lnc_vs_deg="$od/BMK_2_LncRNA/BMK_5_LncRNA_vs_mRNA/BMK_6_DEG_Compare";push (@dirlnc,$lnc_vs_deg);

foreach my $d (@dirlnc){&T($d);&MKDIR($d);}


############### Extract 
if (-d "$id/Data_Assess") {
	my $allsample=glob("$id/Data_Assess/AllSample_GC_Q.stat");
	open (IN,$allsample) or die $!;
	open (OUT,">$QC_data/AllSample_GC_Q.stat") or die $!;
	print OUT "#SampleID\tBMK-ID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ30(%)\n";
	while(<IN>){
		chomp;
		next if (/^#/||/^\s+$/);
		my ($bmkid)=(split /\t+/,$_)[0];
		print OUT "$ID{$bmkid}\t$_\n";
	}
	close IN;
	close OUT;
	system "cp -r $id/Data_Assess/PNG $QC_data/";
	#system "cp $Bin/readme/readme/rawdata.README.txt $rawdata_od";
}

if (-d "$id/Basic_Analysis/geneExpression") {
	system "cp $id/Basic_Analysis/geneExpression/*map.png $QC_map";
	system "cp $id/Basic_Analysis/geneExpression/*type.png $QC_map";
	system "cp $id/Basic_Analysis/geneExpression/*mappedStat.xls $QC_map";
	system "cp $id/Basic_Analysis/geneExpression/*insertSize.r.png $QC_lib";
	system "cp $id/Basic_Analysis/geneExpression/*randcheck.png $QC_lib";
	system "cp $id/Basic_Analysis/geneExpression/*Saturation.png $QC_lib";
	my %map_stat;
	my @map = glob "$id/Basic_Analysis/geneExpression/*.mappedStat.xls";
	open (OUT,">$QC_map/Total_mapped.stat") or die $!;
	print OUT "BMK-ID\tTotal Reads\tMapped Reads\tUniq Mapped Reads\tMultiple Mapped Reads\tReads Map to '+'\tReads Map to '-'\n";
	foreach my $file(@map){
		open (IN,$file) or die $!;
		basename($file)=~/^(\w*)\./;
		my $name=$1;
		while(<IN>){
			chomp;
			my @div=(split /\t/,$_);
			if (/^Total Reads/){
				$map_stat{$name}{'Total Reads'}=$div[1];
			}
			if (/^mapped Reads/){
				$map_stat{$name}{'mapped Reads'}="$div[1]($div[2])";
			}
			if (/^Uniq Map/){
				$map_stat{$name}{'Uniq Map'}="$div[1]($div[2])";
			}
			if (/^Multiple Map/){
				$map_stat{$name}{'Multiple Map'}="$div[1]($div[2])";
			}
			if (/^Only Map Plus Strand/){
				$map_stat{$name}{'Left Map'}="$div[1]($div[2])";
			}
			if (/^Only Map Minus Strand/){
				$map_stat{$name}{'Right Map'}="$div[1]($div[2])";
			}
		}
		close IN;
	}
	foreach my $k (sort keys  %map_stat){
		print OUT "$k\t$map_stat{$k}{'Total Reads'}\t$map_stat{$k}{'mapped Reads'}\t$map_stat{$k}{'Uniq Map'}\t$map_stat{$k}{'Multiple Map'}\t$map_stat{$k}{'Left Map'}\t$map_stat{$k}{'Right Map'}\n";
	}
	close OUT;
	my $fa=glob("$id/Basic_Analysis/geneExpression/final_track/*.newGene.longest_transcript.fa");
	$fa=~/final_track\/(.*).newGene.longest_transcript.fa$/;
	 if (!defined $index){$index=$1;}
	system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene.longest_transcript.fa $newgene/$index.newGene.longest_transcript.fa";
	system "sed -i 's/\t/ /' $newgene/$index.newGene.longest_transcript.fa";###20160115
	system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene_final.filtered.gff $newgene/$index.newGene.gff";
	#system "cp $Bin/readme/readme/geneExpression.README.txt $od/geneExpression";
}
if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lib_type_stat/") {
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lib_type_stat/fq_map.xls $QC_lib/Lib_assess.txt";
}

if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/") {
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lncRNA_filter_pie.stat.png $lnc_pre/LncRNA_classification.png";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/circos/circos.png $lnc_pre/LncRNA_circos.png" if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/circos/circos.png");
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gff $lnc_pre/LncRNA.gff";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gtf $lnc_pre/LncRNA.gtf";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list $lnc_pre/LncRNA.id";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_final.fa $lnc_pre/LncRNA.fa";
	
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPC/lnc_code_filter.result.txt $lnc_pre_soft/CPC.txt";
	
	my $pfamfile="$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/Pfam/Pfam_result.txt" ;
	open IN,"$pfamfile";
	open OUT,">$lnc_pre_soft/Pfam.txt";
	print OUT "trans_id\thmm_acc\thmm_name\thmm_start\thmm_end\thmm_length\tbit_score\tE-value\n";
	while (<IN>){
		chomp;
		next if /^#/ || /^$/;
		my @tmp=split /\s+/,$_;
		my $id=(split /\:/,$tmp[0])[0];
		print OUT "$id\t$tmp[5]\t$tmp[6]\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\n";
	}
	close OUT;
	close IN;
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CNCI/CNCI.index $lnc_pre_soft/CNCI.txt";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt $lnc_pre_soft/CPAT.txt"  if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt");
	
	#############################
	
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/venn.png $lnc_pre_soft";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/list.txt $lnc_pre_soft/Software_veen.txt";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Cis_target/Cis_target_result.txt  $lnc_tar/LncRNA_Cis_target.txt";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target $lnc_tar/LncRNA_LncTar_target.txt" if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target");

	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.txt  $lnc_tar_trans/" if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/");  ####by linhj
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.png  $lnc_tar_trans/" if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/");  ####by linhj
	my @assembly=glob"$id/Basic_Analysis/Tophat_Cufflinks/Cufflinks/*/transcripts.gtf";
	foreach my $ass (@assembly){
		$ass=~/\/Cufflinks\/(.*)\/transcripts.gtf/;
		my $name=$1;
		system "cp $ass $lnc_assem/$name.Cufflinks.transcripts.gtf";
		system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $lnc_assem/$name.Cufflinks.transcripts.gtf";
	}
	if(-d "$id/Basic_Analysis/Tophat_Cufflinks/Scripture/"){
		my @scripture=glob"$id/Basic_Analysis/Tophat_Cufflinks/Scripture/*.gtf";
		foreach my $scr(@scripture){
			if ($scr=~/\/Scripture\/(.*).gtf/){
				my $sampl=$1;
				system "cp -r $scr $lnc_assem/$sampl.Scripture.transcripts.gtf";
				system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $lnc_assem/$sampl.Scripture.transcripts.gtf";
			}
		}
	}
}

if (-d "$id/Alitsplice_Analysis") {
	system "cp -r $id/Alitsplice_Analysis/*.png $as";
	my @aslist=glob("$id/Alitsplice_Analysis/*.fpkm");
	foreach my $fpkm (@aslist){
	  $fpkm=~/Alitsplice_Analysis\/(.*).fpkm$/;
	  my $assam=$1;
	  system	"cut -f 1-8 $fpkm > $as/$assam.AS.list";
	}
	#system "cp $Bin/readme/readme/Alt_splice.README.txt $od/Alt_splice";
}

if (-d "$id/DEG_Analysis") {
	my @DEG_dir=glob "$id/DEG_Analysis/*";
	my %Anno;
	my %Stat;
	system "cp  $id/DEG_Analysis/All_DEG_veen.* $gene_deg_all" if (-f "$id/DEG_Analysis/All_DEG_veen.png");
	system "cp $id/DEG_Analysis/All_gene_counts.list $gene_exp" if (-f "$id/DEG_Analysis/All_gene_counts.list");
	system "cp $id/DEG_Analysis/All_gene_fpkm.list $gene_exp" if (-f "$id/DEG_Analysis/All_gene_fpkm.list");
	system "perl $Bin/anno_stat.pl -i $id/DEG_Analysis -anno all -od $gene_deg_all";
	system "cp -r $id/DEG_Analysis/Anno_Cluster/*_Cluster/GO* $gene_deg_clust " if (-d "$id/DEG_Analysis/Anno_Cluster/");
	system "cp -r $id/DEG_Analysis/Anno_Cluster/*_Cluster/KEGG* $gene_deg_clust " if (-d "$id/DEG_Analysis/Anno_Cluster/");
	my $n=4;
	foreach my $dir (@DEG_dir) {
		if (-d $dir) {
			if ($dir=~/\/All_DEG$/) {
			  system "cp $dir/All.DEG_final.xls $gene_deg_all/All.DEG.Expression.xls";
			  system "cp $dir/final_DEG_annotation.xls $gene_deg_all/All.DEG_gene_annotation.xls";
			  
			  system "cp -r $dir/DEG_Cluster/hierarchical/* $gene_deg_all/" ;
			}
			elsif ($dir=~ /\/kmean_0$/) {
				system "cp -r $dir/k-means.png $gene_deg_all/";
			}
			elsif ($dir=~/\/density$/) {
			  `awk '{for(i=1;i<=NF/2+1;i++) {if(NR==1){if(i==1) printf(\"\%s\\t\%s\\t\",\"sample\",\$i);if(i>1 && i<NF/2+1) printf(\"\%s\\t\",\$i) } else {printf(\"\%s\\t\",\$i)}} printf \"\\n\"}' $dir/free_com.cor|sed 's/decor.//g' > $gene_exp_cor/Sample_correlation.list`;
				system "cp -r $dir/cor_plot/*.png $gene_exp_cor/";
				system "cp  $dir/sample_cluster.png $gene_exp_cor/" if (-f "$dir/sample_cluster.png");
				system "cp -r $dir/*fpkm*png $gene_exp_stat";
            }
			elsif ($dir=~/_vs_/){
				my $nam=basename($dir);
				my $vs_dir="$gene_deg/BMK_".$n."_$nam";
				
				&MKDIR("$vs_dir");
				&MKDIR("$vs_dir/BMK_1_Statistics_Visualization");
				&MKDIR("$vs_dir/BMK_2_DEG_Annotation");
				&MKDIR("$vs_dir/BMK_3_GO_Enrichment");
				&MKDIR("$vs_dir/BMK_4_Pathway_Enrichment");
				system "cp $Bin/PDF/BMK_3_mRNA/BMK_vs.pdf $vs_dir";
				#&MKDIR("$vs_dir/BMK_4_DEG_Cluster");
				
				system "cp -r $dir/*.png $vs_dir/BMK_1_Statistics_Visualization" ;
				system "cp -r $dir/*.DEG_final.xls $vs_dir/BMK_1_Statistics_Visualization" unless $dir=~/.*_vs_.*_vs_.*/;
				system "cp -r $dir/*.DEG_gene.fa $vs_dir/BMK_1_Statistics_Visualization" unless $dir=~/.*_vs_.*_vs_.*/;
				system "cp -r $dir/DEG_Cluster/hierarchical/* $vs_dir/BMK_1_Statistics_Visualization"; 
				system "cp -r $dir/Anno_enrichment/*_vs_*.annotation.xls $vs_dir/BMK_2_DEG_Annotation/";
				system "cp -r $dir/Anno_enrichment/Cog_Anno/* $vs_dir/BMK_2_DEG_Annotation/";
				system "cp -r $dir/Anno_enrichment/eggNOG_Anno/* $vs_dir/BMK_2_DEG_Annotation/";
				system "cp -r $dir/Anno_enrichment/pathway/* $vs_dir/BMK_4_Pathway_Enrichment/";
				system "cp -r $dir/Anno_enrichment/Graph/*KEGG* $vs_dir/BMK_4_Pathway_Enrichment/kegg_enrichment/";
				system "rm $vs_dir/BMK_4_Pathway_Enrichment/kegg_enrichment/*.stat";
				system "rm $vs_dir/BMK_4_Pathway_Enrichment/kegg_enrichment/*.svg";
				system "cp -r $dir/Anno_enrichment/Graph/*topGO* $vs_dir/BMK_3_GO_Enrichment/";
				system "cp -r $dir/Anno_enrichment/go_enrichment/*.GO_enrichment.stat.xls $vs_dir/BMK_3_GO_Enrichment/$nam.GO_stat.xls";
				system "cp -r $dir/Anno_enrichment/go_enrichment/*png $vs_dir/BMK_3_GO_Enrichment/";
				
				system "cp -r $dir/$nam.DEG.CytoscapeInput.txt $gene_deg_ppi";

				&rename_recursion("$vs_dir/BMK_1_Statistics_Visualization","FC_count","MA");&rename_recursion("$vs_dir/BMK_1_Statistics_Visualization","FC_FDR","Volcano");
				#system "rm -r $DEG/$nam/Kog_Anno" if (-d "$DEG/$nam/Kog_Anno");
				$n++;
			}
		}
	}
}
&rename_recursion($gene,"Kegg","KEGG");&rename_recursion($gene,"Cog","COG");
if (-d "$id/Anno_Integrate/New_Anno/Result/") {
	
    system "cp -r $id/Anno_Integrate/New_Anno/Result/Function_Annotation.stat.xls $newgene_ann/";
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/All_Database_annotation.xls $newgene_ann/";
	#$newgene_ann_a;$newgene_ann_st;$newgene_ann_kegg;
	 my $pre=(glob "$id/Anno_Integrate/New_Anno/Result/*longest_transcript.fa.nr.anno.txt")[0];
	 $pre=basename($pre);
	 $pre=~/(.*).nr.anno.txt/;
	 my $new=$1;
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/*png $newgene_ann_st/";
	 #system "cp -r $id/Anno_Integrate/New_Anno/Result/*.fa.GO.pdf $newgene_ann_st/";
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/*Cog_class.txt $newgene_ann_st/";
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/*GO_enrichment.stat.xls $newgene_ann_st/$new.GO.stat";
	
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/*.nr.lib.stat $newgene_ann_st/";
	 #system "cp -r $id/Anno_Integrate/New_Anno/Result/*GO_enrichment.stat.xls $newgene_ann_st/";
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/*.txt $newgene_ann_a/";
	 system "cp $id/Anno_Integrate/New_Anno/Result/*Kegg.* $newgene_ann_a/";
	 system "cp -r $id/Anno_Integrate/New_Anno/Result/Kegg_map/* $newgene_ann_kegg/ ";
	 #system "cp -r $id/Anno_Integrate/New_Anno/Result/*GO_tree.stat.xls $newgene_ann_st/";
	 system "rm $newgene_ann_a/*Kegg.globe" if (-f "$newgene_ann_a/$new.Kegg.globe");
	 #system "rm $newgene_ann_a/*GO.list.txt " if (-f "$newgene_ann_a/$new.GO.list.txt");
	 
	&rename_recursion($newgene_ann_a,"Kegg","KEGG");&rename_recursion($newgene_ann_st,"Kegg","KEGG");
	&rename_recursion($newgene_ann_a,"Kog","KOG");&rename_recursion($newgene_ann_st,"Kog","KOG");
	&rename_recursion($newgene_ann_a,"Cog","COG");&rename_recursion($newgene_ann_st,"Cog","COG");
	&rename_recursion($newgene_ann_a,"nr","NR");&rename_recursion($newgene_ann_st,"nr","NR");
	&rename_recursion($newgene_ann_a,"Swissprot","SwissProt");&rename_recursion($newgene_ann_st,"Swissprot","SwissProt");

	my $function_anno = "$newgene_ann/Function_Annotation.stat.xls";
	my $tmp = "$newgene_ann/tmp.txt";
	open (IN,$function_anno) or die $!;
	open (OUT,">$tmp") or die $!;
	print OUT "Annotated databases\tNew Gene Number\n";
	while (<IN>) {
		chomp;
		next if (/^#/);
		my @div=(split /\t/,$_);
		print OUT "$div[0]\t$div[1]\n";
	}
	close IN;
	close OUT;
	system "rm $function_anno";
	system "mv $tmp $function_anno";
	
	#system "cp $Bin/readme/readme/NewGene.README.txt $od/NewGene";
}
if (-d "$id/Lnc_Diff_Analysis") {
	my @Lnc_Diff_dir=glob "$id/Lnc_Diff_Analysis/*";
	my %Anno;
	my %Stat;
        system "cp $id/Lnc_Diff_Analysis/*.png $lnc_deg_all/" if (-f "$id/Lnc_Diff_Analysis/All_DEG_venn.png");
        system "cp $id/Lnc_Diff_Analysis/*.genes $lnc_deg_all/" if (-f "$id/Lnc_Diff_Analysis/All_DEG_venn.genes");
        
        system "cp $id/Lnc_Diff_Analysis/LncRNA_counts.list $lnc_exp/" if (-f "$id/Lnc_Diff_Analysis/LncRNA_counts.list");
        system "cp $id/Lnc_Diff_Analysis/LncRNA_fpkm.list $lnc_exp/" if (-f "$id/Lnc_Diff_Analysis/LncRNA_fpkm.list");
	system "perl $Bin/anno_stat.pl -i $id/Lnc_Diff_Analysis -anno cis -od $lnc_deg_all/";
	system "perl $Bin/anno_stat.pl -i $id/Lnc_Diff_Analysis -anno trans -od $lnc_deg_all/";
	system "cp -r $id/Lnc_Diff_Analysis/Cis_Anno_Cluster/GO_Cluster/GO*  $lnc_deg_cis/ " if (-d "$id/Lnc_Diff_Analysis/Cis_Anno_Cluster/");
	system "cp -r $id/Lnc_Diff_Analysis/Cis_Anno_Cluster/KEGG_Cluster/KEGG*  $lnc_deg_cis/ " if (-d "$id/Lnc_Diff_Analysis/Cis_Anno_Cluster/");
	system "cp -r $id/Lnc_Diff_Analysis/Trans_Anno_Cluster/GO_Cluster/GO*  $lnc_deg_trans/" if (-d "$id/Lnc_Diff_Analysis/Trans_Anno_Cluster/");
	system "cp -r $id/Lnc_Diff_Analysis/Trans_Anno_Cluster/KEGG_Cluster/KEGG*  $lnc_deg_trans/" if (-d "$id/Lnc_Diff_Analysis/Trans_Anno_Cluster/");
	 

	my $ln=4;
	foreach my $dir (@Lnc_Diff_dir) {
            
		if (-d $dir ){
			if ($dir=~/\/All_DEG$/) {
				system "cp $dir/All.DEG_final.xls $lnc_deg_all/All.DEG.Expression.xls";
			 #system "cp $dir/final_DEG_annotation.xls $gene_deg_all/All.DEG_gene_annotation.xls";
			  
			  #system "cp  $dir/DEG_Cluster/hierarchical/* $lnc_deg_all/" if (-f  "$dir/DEG_Cluster/hierarchical/*");
			  system "cp  $dir/DEG_Cluster/hierarchical/* $lnc_deg_all/";
			}
			
			elsif($dir=~/\/kmean_0$/){
				system"cp -r $dir/k-means.png $lnc_deg_all" if (-f "$dir/k-means.png");
				#system "rm -r $LncRNA/DEG/kmean_0";
			}
			elsif ($dir=~/\/density$/) {
				`awk '{for(i=1;i<=NF/2+1;i++) {if(NR==1){if(i==1) printf(\"\%s\\t\%s\\t\",\"sample\",\$i);if(i>1 && i<NF/2+1) printf(\"\%s\\t\",\$i) } else {printf(\"\%s\\t\",\$i)}} printf \"\\n\"}' $dir/free_com.cor|sed 's/decor.//g' > $lnc_exp_cor/correlation.txt`;
				system "cp -r $dir/cor_plot/*.png $lnc_exp_cor/";
				system "cp -r $dir/sample_cluster.png $lnc_exp_cor/";
				system "cp -r $dir/*fpkm*png $lnc_exp_stat";
			}
			elsif ($dir=~/_vs_/){
				my $nam=basename($dir);
				my $vs_dir="$lnc_deg/BMK_".$ln."_$nam";
				
				&MKDIR("$vs_dir");
				&MKDIR("$vs_dir/BMK_1_Statistics_Visualization");
				&MKDIR("$vs_dir/BMK_2_Cis_Anno_enrichment");
				&MKDIR("$vs_dir/BMK_2_Cis_Anno_enrichment/BMK_1_DEG_Annotation");
				&MKDIR("$vs_dir/BMK_2_Cis_Anno_enrichment/BMK_2_GO_Enrichment");
				&MKDIR("$vs_dir/BMK_2_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment");
				&MKDIR("$vs_dir/BMK_3_Trans_Anno_enrichment");
				&MKDIR("$vs_dir/BMK_3_Trans_Anno_enrichment/BMK_1_DEG_Annotation");
				&MKDIR("$vs_dir/BMK_3_Trans_Anno_enrichment/BMK_2_GO_Enrichment");
				&MKDIR("$vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment");
				system "cp $Bin/PDF/BMK_2_LncRNA/BMK_vs.pdf $vs_dir";
				
				system "cp -r $dir/DEG_Cluster/hierarchical/* $vs_dir/BMK_1_Statistics_Visualization"; 
				system "cp  $dir/$nam.DEG_final.xls $vs_dir/BMK_1_Statistics_Visualization/$nam.DEG.xls";
                system "cp  $dir/$nam.deg_lnc2Target_m_Cytoscape.input.txt $vs_dir/BMK_1_Statistics_Visualization/$nam.DEG_lncRNA2Target_Cytoscape.input.txt" if (-f  "$dir/$nam.deg_lnc2Target_m_Cytoscape.input.txt");
                system "cp  $dir/$nam.DEG.CytoscapeInput.txt $lnc_deg_ppi" if (-f  "$dir/$nam.DEG.CytoscapeInput.txt");
				system "cp -r $dir/*.png $vs_dir/BMK_1_Statistics_Visualization" unless $dir=~/.*_vs_.*_vs_.*/;
				######Cis result
				system "cp -r $dir/Cis_Anno_enrichment/*_vs_*.annotation.xls $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_1_DEG_Annotation/$nam.target_gene_annotation.xls" if (-f "$dir/Cis_Anno_enrichment/$nam.annotation.xls");
				system "cp -r $dir/Cis_Anno_enrichment/*_vs_*.lncRNA_DEG_targene.xls $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_1_DEG_Annotation/" if (-f "$dir/Cis_Anno_enrichment/$nam.lncRNA_DEG_targene.xls");
				system "cp -r $dir/Cis_Anno_enrichment/Cog_Anno/* $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_1_DEG_Annotation/";
				system "cp -r $dir/Cis_Anno_enrichment/eggNOG_Anno/* $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_1_DEG_Annotation/";
				system "cp -r $dir/Cis_Anno_enrichment/pathway/* $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/";
				system "cp -r $dir/Cis_Anno_enrichment/Graph/*KEGG* $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/";
				system "rm $vs_dir/BMK_3_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.stat" if (-f "$vs_dir/BMK_3_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.stat");
				system "rm $vs_dir/BMK_3_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.svg" if (-f "$vs_dir/BMK_3_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.svg");
				system "cp -r $dir/Cis_Anno_enrichment/Graph/*vs*topGO* $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_2_GO_Enrichment/";
				
				system "cp -r $dir/Cis_Anno_enrichment/go_enrichment/*.GO_enrichment.stat.xls $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_2_GO_Enrichment/$nam.GO_stat.xls";
				system "cp -r $dir/Cis_Anno_enrichment/go_enrichment/*png $vs_dir/BMK_2_Cis_Anno_enrichment/BMK_2_GO_Enrichment/";
				########Trans result
				system "cp -r $dir/Trans_Anno_enrichment/*_vs_*.annotation.xls $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_1_DEG_Annotation/$nam.target_gene_annotation.xls";
				system "cp -r $dir/Trans_Anno_enrichment/*_vs_*.lncRNA_DEG_targene.xls $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_1_DEG_Annotation/";
				system "cp -r $dir/Trans_Anno_enrichment/Cog_Anno/* $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_1_DEG_Annotation/";
				system "cp -r $dir/Trans_Anno_enrichment/eggNOG_Anno/* $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_1_DEG_Annotation/";
				system "cp -r $dir/Trans_Anno_enrichment/pathway/* $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/";
				system "cp -r $dir/Trans_Anno_enrichment/Graph/*KEGG* $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/";
				system "rm $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.stat" if (-f "$vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.stat");
				system "rm $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.svg" if (-f "$vs_dir/BMK_3_Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/$nam.KEGG.svg");
				system "cp -r $dir/Trans_Anno_enrichment/Graph/*topGO* $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_2_GO_Enrichment";
				system "cp -r $dir/Trans_Anno_enrichment/go_enrichment/*.GO_enrichment.stat.xls $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_2_GO_Enrichment/$nam.GO_stat.xls";
				system "cp -r $dir/Trans_Anno_enrichment/go_enrichment/*png $vs_dir/BMK_3_Trans_Anno_enrichment/BMK_2_GO_Enrichment/";
				$ln++;
				&rename_recursion("$vs_dir/BMK_1_Statistics_Visualization","FC_count","MA");&rename_recursion("$vs_dir/BMK_1_Statistics_Visualization","FC_FDR","Volcano");
			}
		}
	}
}
&rename_recursion($LncRNA,"Kegg","KEGG");&rename_recursion($LncRNA,"Cog","COG");
if(-d "$id/SNP_Analysis") {
	system "cp -r $id/SNP_Analysis/stat/AllSample.* $snp";
	system "cp -r $id/SNP_Analysis/stat/*type* $snp_type";
	system "cp -r $id/SNP_Analysis/stat/final.*.anno.gatk.all.list $snp";
	system "cp -r $id/SNP_Analysis/stat/indel_anno/*.png  $snp_indel";
	system "cp -r $id/SNP_Analysis/stat/indel_anno/*.stat  $snp_indel";
	system "cp -r $id/SNP_Analysis/stat/snp_anno/*.png  $snp_ann";
	system "cp -r $id/SNP_Analysis/stat/snp_anno/*.stat  $snp_ann";
	&rename_recursion($snp,"snp","SNP");&rename_recursion($snp,"indel","InDel");
	&rename_recursion($snp_ann,"snp","SNP");&rename_recursion($snp_ann,"indel","InDel");
	&rename_recursion($snp_type,"snp","SNP");&rename_recursion($snp_type,"indel","InDel");
	&rename_recursion($snp_indel,"snp","SNP");&rename_recursion($snp_indel,"Indel","InDel");
}

if (-d "$id/Gene_Structure_Optimize") {
	system "cp -r $id/Gene_Structure_Optimize/*xls $structure/";
	#system "cp $Bin/readme/readme/Gene_Structure_Optimize.README.txt $od/Gene_Structure_Optimize";
}

if (-d "$id/Compare_analysis"){

	system "cp -r $id/Compare_analysis/lncRNA_mRNA_DEG/*png $lnc_vs_deg";
    #system "cp -r $id/Compare_analysis/chr_num_stat $LncRNA/Compare/";
	system "cp -r $id/Compare_analysis/*RNA.len.png $lnc_vs_len";
	system "cp -r $id/Compare_analysis/*exon.png $lnc_vs_exon";
	system "cp -r $id/Compare_analysis/gene.orf.png $lnc_vs_orf/mRNA.ORF.png";
	system "cp -r $id/Compare_analysis/lnc_filter_final.orf.png $lnc_vs_orf/lncRNA.ORF.png";
	system "cp -r $id/Compare_analysis/*FPKM.xls $lnc_vs_exp";
	system "cp -r $id/Compare_analysis/*fpkm.png $lnc_vs_exp";
	system "cp -r $id/Compare_analysis/*isoform.png $lnc_vs_iso";
	#system "cp -r $id/Compare_analysis/*RNA_isoform.txt $lnc_vs_iso";	
	&rename_recursion("$lnc_vs_deg","FC_count","MA");&rename_recursion("$lnc_vs_deg","FC_FDR","Volcano");
}

if (-d "$id/DEU_analysis") {
	#&MKDIR("$GENE/DEU");
	my @DEU_dir = glob"$id/DEU_analysis/*";
	foreach my $dir (@DEU_dir) {
		if ($dir =~ /_vs_/) {
			my $name = basename($dir);
			&MKDIR("$deu/$name");
			system "cp -r $dir/DEU_Result_Final.xls $deu/$name/$name.DEU_Result_Final.xls";
			system "cp -r $dir/DEU_Result_All.xls $deu/$name/$name.DEU_Result_All.xls";
			system "cp -r $dir/DEXSeqReport $deu/$name/" if (-d "$dir/DEXSeqReport");
			#system "cp -r $dir/DEXSeqReport/testForDEU.html $deu/$name/DEXSeqReport/testForDEU.html.cloud";
		}
	}
}


my $per=1;
system "rm -r $od/BMK_4_Personality " if (-d "$od/BMK_4_Personality");
&MKDIR("$od/BMK_4_Personality");
my $per_dir="$od/BMK_4_Personality";
if ( -d "$id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/Known_lncRNA"){
	
	my $line=`less $id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/Known_lncRNA/known_lncRNA.result.txt |wc -l`;
	chomp $line;
	&MKDIR("$per_dir/BMK_$per"."_Known_LncRNA") if ($line > 1);
	my $known="$per_dir/BMK_$per"."_Known_LncRNA";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/Known_lncRNA/known_lncRNA.result.txt $known/" if ($line > 1);
	$per++ if ($line > 1 );
	}
if (-d "$id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/miRNA_Target2LncRNA/"){
	&MKDIR("$per_dir/BMK_$per"."_miRNA_Target2LncRNA");
	my $s2l="$per_dir/BMK_$per"."_miRNA_Target2LncRNA";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/miRNA_Target2LncRNA/lncRNA_target2mirna.list $s2l";
	$per++;
	#system "cp $id/Personality/miRNA_Target2LncRNA/*.mir2target.stat $od/Personality/miRNA_Target2LncRNA";
}
if (-d "$id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/precursor"){
	&MKDIR("$per_dir/BMK_$per"."_Precursor_Result");
	my $precursor="$per_dir/BMK_$per"."_Precursor_Result";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/precursor/lncRNA_precursor.txt $precursor";
	$per++;
}
if (-d "$id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/mRNA_lncRNA_miRNA"){
	&MKDIR("$per_dir/BMK_$per"."_mRNA_lncRNA_miRNA_Regulation");
	my $thr_regulat="$per_dir/BMK_$per"."_mRNA_lncRNA_miRNA_Regulation";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/LncRNA_com_miRNA/mRNA_lncRNA_miRNA/* $thr_regulat";
	$per++;
}
if (-d "$id/CircRNA_Analysis") {
	&MKDIR("$per_dir/BMK_$per"."_CircRNA_Analysis");
	my $circRNA="$per_dir/BMK_$per"."_CircRNA_Analysis";
    system "cut -f 1-10 $id/CircRNA_Analysis/circRNA_identify/All_CircRNA.xls > $circRNA/All_CircRNA.xls" if (-f "$id/CircRNA_Analysis/circRNA_identify/All_CircRNA.xls");
	system "cp $id/CircRNA_Analysis/expression/sample_venn.png $circRNA" if (-f "$id/CircRNA_Analysis/expression/sample_venn.png");
	system "cp $id/CircRNA_Analysis/statistics/circRNA.type.pie.png $circRNA" if (-f "$id/CircRNA_Analysis/statistics/circRNA.type.pie.png");
	system "cp $id/CircRNA_Analysis/statistics/chromosome_distrbution.png $circRNA" if (-f "$id/CircRNA_Analysis/statistics/chromosome_distrbution.png");
	system "cp $id/CircRNA_Analysis/expression/All_gene_counts.list $circRNA" if (-f "$id/CircRNA_Analysis/expression/All_gene_counts.list");
	$per++;
}
if (-d "$id/Personality"){	
	if (-d "$id/Personality/Conservation"){
		&MKDIR("$per_dir/BMK_$per"."_Conservation");
		my $conservation="$per_dir/BMK_$per"."_Conservation";
		system "cp $id/Personality/Conservation/*.png $conservation";
		$per++;
	}
	
	if (-d "$id/Personality/TF"){
		&MKDIR("$per_dir/BMK_$per"."_TF");
		my $tf="$per_dir/BMK_$per"."_TF";
		system "cp $id/Personality/TF/final.xls $tf/TF.txt";
		$per++;
	}
	if (-d "$id/Personality/Tissue_specific"){
		&MKDIR("$per_dir/BMK_$per"."_Tissue_specific");
		my $ts="$per_dir/BMK_$per"."_Tissue_specific";
		system "cp $id/Personality/Tissue_specific/*.png $ts";
		$per++;
	}
	if(-d "$id/Personality/PCA"){
		&MKDIR("$per_dir/BMK_$per"."_PCA");
		my $pca="$per_dir/BMK_$per"."_PCA";
		system "cp $id/Personality/PCA/*.png $pca";
		$per++;
	}
	
}
#system "cp $Bin/../README $od/";
################copy the readme.pdf to $od
system "perl $Bin/Readme_copy.pl -in $od ";
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $Script Time :[$Time_End]\n\n";
###############Subs
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
sub rename_recursion
{
        my($dir,$key,$replace)=@_;
        my %options;
        #目录路径
        my @cases;
        if (-d $dir) {#判断目录是否存在
                my @files;
                my $dh;
                push(@files, $dir);
				push @cases,$dir;
                while (@files) {
                        if (-d $files[0])
                        {#若是目录执行以下操作
                                opendir $dh, $files[0] or die $!;#打开目录句柄,若失败打印错误信息
                                @_ = grep { /^[^\.]/ } readdir $dh;#过滤掉以"."和".."的文件,即UNIX下的隐藏文件
                                foreach (@_) {
                                        if ($_ !~ /\./) {
                                                s/\&/\\&/ if($_=~/\&/);
                                                push @files,"$files[0]/".$_."/";
                                                push @cases,"$files[0]/".$_."/"; ##
                                        }
                                }
                                closedir $dh;
                        }
                        shift @files;
                }
        }

`rename $key $replace $_/*` foreach @cases;#打印文件列表
#print "$_\n" foreach @cases;
}
sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}
sub T{
	my $id=shift;
	$id = '"'.$id.'"';
	return $id;
}
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract LncRNA Analysis Results for Html Process;
	version:$ver
	Usage:
		-id  <STR>   input dir, analysis output directory   force
		-od  <STR>   result output dir                      force
		-cfg		data.cfg			force
		-cfg2		detail.cfg
		-h           help
	Usage End.
		exit;
}
