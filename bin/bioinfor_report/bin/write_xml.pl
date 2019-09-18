#!/usr/bin/perl -w
use autodie qw(:all);
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw (abs_path getcwd);
use File::Spec;
use File::Basename qw (basename dirname);
my $BEGIN_TIME = time ();
my $version = "1.0.0";
use RTF::Writer;
use Encode qw (decode);

###############################################################

# ----------------------------------------------------------------
# GetOptions
# -----------------------------------------------------------------
my ($dir,$config,$xml);
GetOptions (
	"h|?" => \&USAGE,
	"indir:s" => \$dir,
	"config:s" => \$config,
	"xml:s"=>\$xml,
) or &USAGE;
&USAGE unless ($dir and $config and $xml);

my $type="type1";
my $dir_abs = abs_path ($dir);
my $dir_template ="$Bin/Template";
my $tem="$dir/tem";
mkdir $tem unless (-d $tem);
#system "cp -r $template $dir_abs";

#my $dir_template = "$dir/Template";

my %data;
my %table_info;
my %pic_info;
my %data_among_words;
my %Ref_info;
my %config;
my %reference;

open (IN,"$Bin/LncRNA_reference.txt") or die $!;
while (<IN>) {
	chomp;
	next if (/^#/);
	my ($id,$link)=(split /\t/,$_);
	$link=~s/\s+$//;
	$reference{$id}=$link;
}
close IN;

system "perl $Bin/data.pl -indir $dir_abs -od $dir_abs";



open (IN,"$dir/data.txt") or die $!;
while (<IN>) {
	chomp;
	my ($k1,$v1) = split /\t+/,$_;
	$data{$k1}=$v1;
}
close IN;
my @com;
my @sep;
open (IN,"$config") or die $!;
while (<IN>) {
	chomp;
	s/\s+$//;s/\r$//;
	next if (/^#/ || /^\s+$/);
	my ($k,$v) = (split /\s+/,$_)[0,1];
	if ($k =~/^Com/){
		push @com,$v;
	}
	if($k =~ /^Sep/){
		push @sep,$v;
	}
	if ($k =~ /^Project_name/){
		$config{$k}=$v;
	}
	if ($k =~ /^min/){
		$config{$k}=$v;
	}
	if($k =~ /^test_time/){
		$config{$k}=$v;
	}
	if ($k =~ /^info_time/){
		$config{$k}=$v;
	}
	if ($k =~ /^start_time/){
		$config{$k}=$v;
	}
	if ($k =~ /^finish_data/){
		$config{$k}=$v;
	}
	if ($k=~/^ref_gene_name/){
		$config{$k}=$v;
	}
	if ($k=~/^ref_gene_addr/){
		$config{$k}=$v;
	}
}
close IN;

my $com_num=scalar (@com);
my $sep_num=scalar(@sep);
my $total=$com_num+$sep_num;
print "$com_num\n";
print "$sep_num\n";
my %templat;



my $ftemplate;
my @ftemplates;
if ($total>1 && $total<6){
	$ftemplate="$Bin/1.1.Basic_tem.txt";
}
if ($total<2 || $total>5){
	$ftemplate="$Bin/1.1.Basic_tem_no_deg_veen.txt";
}
push @ftemplates,$ftemplate;
print "$ftemplate\n";
if ($sep_num>0){
	#system "cat $ftemplate $Bin/2.DEU_tem.txt >>$dir/template";
	$ftemplate="$Bin/2.DEU_tem.txt";
	push @ftemplates,$ftemplate;
}
if(-d "$dir/Personality/Conservation"){
	#system "cat $Bin/4.conservation.txt >>$dir/template";
	$ftemplate="$Bin/4.conservation.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/Known_LncRNA"){
	#system "cat $Bin/6.Known_lncRNA.txt >>$dir/template";
	$ftemplate="$Bin/6.Known_lncRNA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/Tissue_specific"){
	system "cat $Bin/5.tissue_specify.txt >>$dir/template";
}
if (-d "$dir/Personality/precursor"){
	#system "cat $Bin/7.precusor.txt >>$dir/template";
	$ftemplate="$Bin/7.precusor.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/miRNA_Target2LncRNA"){
	#system "cat $Bin/8.miRNA_target2lncRNA.txt >>$dir/template";
	$ftemplate="$Bin/8.miRNA_target2lncRNA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/PCA"){
	$ftemplate="$Bin/10.PCA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/WGCNA"){
	$ftemplate="$Bin/11.WGCNA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/TF"){
	$ftemplate="$Bin/9.TF.txt";
	push @ftemplates,$ftemplate;
}

#system "cat $Bin/3.attach_tem.txt >>$dir/template";
$ftemplate="$Bin/3.attach_tem.txt";
push @ftemplates,$ftemplate;

`cat @ftemplates >$dir/template`;
$ftemplate="$dir/template";
print "$ftemplate\n";



my $finish_data=&time_format($config{finish_data});
my $sample_num = $data{Sample_num};
my $base_num = $data{base_num};
my $min_data = $data{min_data};
my $min_Q30 = $data{min_Q30};
my $min_map_ratio = $data{min_map_ratio};
my $max_map_ratio = $data{max_map_ratio};
my $new_gene_num = $data{new_gene_num};
my $new_gene_ann_num = $data{new_gene_ann_num};
my $deg_num = $data{deg_num};
my $lnc_num = $data{lnc_num};
my $diff_lnc = $data{diff_lnc};
my $optimized_gene_num = $data{optimized_gene_num};
my $prefix = $config{Project_name};
my $min = $config{min};
my $test_time = &time_format($config{test_time});
my $info_time = &time_format($config{info_time});
my $start_time = &time_format($config{start_time});
my $ref_gene_addr=$config{ref_gene_addr};
my $ref_gene_name=$config{ref_gene_name};
&get_data_among_words();




chomp (my $user = `whoami`);
my $user_addr = $user."\@biomarker.com.cn";
my $report_version = &T ("1.2");
my $report_name = &substitute($prefix);
$report_name=&T($prefix);
my $report_code = &T("XXX");
my $report_user = &T($user);
my $report_user_addr = &T($user_addr);
my $report_time = &T(&GetTime);






# -------------------------------------------------------------------------
# output xml
# -------------------------------------------------------------------------
open (OUT,">$xml") or die $!;

print OUT "<?xml version=\"1.0\" encoding=\"utf-8\"?> \n";
print OUT "<report> \n";
print OUT "<report_version value=$report_version	\/> \n";
print OUT "<report_name value=$report_name	\/> \n";
print OUT "<report_code value=$report_code	\/> \n";
print OUT "<report_user value=$report_user	\/> \n";
print OUT "<report_user_addr value=$report_user_addr	\/> \n";
print OUT "<report_time value=$report_time	\/> \n";
print OUT "<report_abstract value='' \/> \n";



open (TEM, $ftemplate) or die $!;
$/="》";
my $i=0;
while (<TEM>){
	chomp;
	next if (/^$/ || /^\#/);
	my (undef,$format,$context)=split /\n/,$_,3;
	$i++;
#	print "$i:\t $context\n";
	$format =~ s/\s+$//;
	$context =~ s/\s+$//;
	&xml_write($format,$type,$context);
}
print OUT "<\/report> \n";

close TEM;
close OUT;
$/="\n"	;    ### 
	
sub xml_write {
	my ($format,$type,$text) = @_;
	if ($format eq "表格") {### 表格
		print OUT &table_write($format,$type,$text);
	}
	elsif ($format =~ /图片/) {### 图片
		print OUT &picture_write($format,$type,$text);
	}
	#elsif ($format =~ /链接/) { ## 链接
	#	print OUT &Link($format,$type,$text);
	#}
	elsif ($format =~ /公式/) {### 公式
		print OUT &picture_write($format,$type,$text);
	}
	elsif ($format =~ /级标题/) {### 标题
		print OUT &Head($format,$type,$text);
	}
	elsif ($format =~ /正文/) {### 正文
		print OUT &Text($format,$type,$text);
	}
	elsif ($format eq '注释') {
		print OUT &Anno($format,$type,$text);
	}
	elsif ($format eq '参考文献') {
		print OUT &Reference ($format,$type,$text);
	}
	elsif ($format eq '附录标题') {
		print OUT &Attach($format,$type,$text);
	}
}



sub table_write {
	my ($name,$type,$desc) = @_;
	my @tables = glob (&get_table($desc));
	my $content;
	my $nam=$desc;
	$type = "type1";
	$type = &T($type);
	$desc = &T($desc);
	chomp $desc;
	#print "$desc\n";
	my $table_num = @tables;
	if (1<$table_num) {
		if ($desc=~/DEU分析结果示意表/){
			$content="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				print "$table\n";
				$table=~/\/(\w+_vs_\w+)\//;
				my $table_name=$1;
				#print "$table_name\n";
				$table_name=$table_name.".DEU";
				$table_name=&T($table_name);
				my $path=&T($table);
				my $action=&T("type1");
				$content=$content."<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .='</file_list>'."\n";
		}
		elsif($desc=~/DEU分析结果/){
			$content="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				$table=~/\/(\w+_vs_\w+)\//;
				my $table_name=$1;
				$table_name=$table_name.".DEU";
				$table_name=&T($table_name);
				my $path=&T($table);
				my $action=&T("type1");
				$content .= "<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .='</file_list>'."\n" ;
		}
		else{
			$content = "<file_list name=$desc type=$type >\n";
			foreach my $table (@tables) {
				my $table_name = basename($table);
				&table_html($table,"$dir/tem/$table_name.html","$nam","","","../Web_Report/src","") ;
				my $path = &T("$dir/tem/$table_name.html") if ($table!~/html$/);
				my $action = &T("xls");
				$table_name=&T($table_name);
				$content = $content."<file	name=$table_name	type=$type	path=$path	action=$action	/> \n";
			}
			$content = $content.'</file_list>'."\n";
		}
	}
	elsif($table_num==1) {
		my $path = $tables[0];
		if (-f "$path"){
			my $table_name = $desc;
			my $line=`wc -l $path`; chomp($line);
			my $line_num = (split /\s+/,$line)[0];
			$path=&T($path);
			if (30>$line_num) {
				$content = "<table	name=$table_name	type=$type	path=$path	\/>\n";
			}
			else {
				my $action = &T("xls");
				$content = "<file	name=$table_name	type=$type	path=$path	action=$action	/> \n";
			}
		}
	}
	return $content;
}




sub get_table {
	chomp (my $id = shift);
	%table_info=(
		"碱基质量值与碱基识别出错的概率的对应关系表" => "$dir_template/error.txt",
		"样品测序数据评估统计表" => "$dir/QC/Clean_Data/AllSample_GC_Q.stat",
		"样品测序数据与所选参考基因组的序列比对结果统计表" => "$dir/QC/Map_assess/Total_mapped.stat",
		"SNP位点信息" => "$dir/mRNA/Structure/SNP/snp_anno/final.snp.anno.gatk.*.list",
		"InDel位点信息" => "$dir/mRNA/Structure/SNP/indel_anno/final.indel.anno.gatk.*.list",
		"SNP位点统计表" => "$dir/mRNA/Structure/SNP/AllSample.snp.stat",
		"可变剪切事件结构和表达量统计表" => "$dir/mRNA/Structure/Alt_splice/*.fpkm",
		"基因结构优化结果" => "$dir/mRNA/Gene_Structure_Optimize/*.geneStructure.optimize.xls",
		"新基因的GFF文件" => "$dir/mRNA/NewGene/*.newGene_final.filtered.gff",
		"新基因序列FASTA文件" => "$dir/mRNA/NewGene/*.newGene.longest_transcript.fa",
		"新基因功能注释结果统计" => "$dir/mRNA/NewGene/Annotation/Function_Annotation.stat.xls",
		"基因表达量结果文件" => "$dir/mRNA/Expression/*.geneExpression.xls",
		"差异表达分析结果" => "$dir/mRNA/DEG/*_vs_*/*_vs_*.DEG_final.xls",
		"差异表达基因数目统计表" => "$dir/mRNA/DEG/DEG.stat",
		"注释的差异表达基因数量统计表" => "$dir/mRNA/DEG/DEG.anno.stat",
		"差异表达基因topGO富集结果" => "$dir/mRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.xls",
		"差异表达基因的KEGG富集结果" => "$dir/mRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.stat",
		#"差异表达基因的KEGG通路注释" => "$dir/mRNA/DEG/*_vs_*/pathway/kegg_map/*.html",
		"DEU分析结果示意表" => "$dir/mRNA/DEU/*_vs_*/DEU_Result_Final.xls",
		"DEU分析结果"=>"$dir/mRNA/DEU/*_vs_*/DEXSeqReport/testForDEU.html",
		"Cufflinks拼接结果" => "$dir/LncRNA/Assembly/*.Cufflinks.transcripts.gtf",
		"Scripture拼接结果" => "$dir/LncRNA/Assembly/*.Scripture.transcripts.gtf",
		"CPC分析结果统计" => "$dir/LncRNA/Identify/CPC.txt",
		"CNCI分析结果统计" => "$dir/LncRNA/Identify/CNCI.txt",
		"pfam分析结果统计" => "$dir/LncRNA/Identify/Pfam.txt",
		"CPAT分析结果统计" => "$dir/LncRNA/Identify/cpat.txt",
		"基于位置关系的LncRNA靶基因预测结果" => "$dir/LncRNA/Target/lncRNA_position.target",
		"基于互补序列的LncRNA靶基因预测结果"=>"$dir/LncRNA/Target/lncRNA_basepair.target",
		"lncRNA表达量结果" => "$dir/LncRNA/DEG/LncRNA_fpkm.list",
		"差异表达lncRNA结果" => "$dir/LncRNA/DEG/*_vs_*/*_vs_*.DEG_final.xls",
		"差异表达lncRNA数目统计表" => "$dir/LncRNA/DEG/Diff_Lnc.stat",
		"注释的差异表达lncRNA靶基因数量统计表" => "$dir/LncRNA/DEG/Diff_lnc.anno.stat",
		"差异表达lncRNA靶基因topGO富集结果" => "$dir/LncRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.xls",
		"差异表达lncRNA靶基因KEGG富集结果" => "$dir/LncRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.stat",
		#"差异表达lncRNA靶基因的KEGG通路注释"=>"$dir/LncRNA/DEG/*_vs_*/pathway/kegg_map/*.html",
		"附表1 软件列表" => "$dir_template/software_used_list.txt",
		"附表2 数据库列表" => "$dir_template/database_used_list.txt",
		"附表3 核酸编码表" => "$dir_template/Nucleic_acids_encoding_table.txt",
		"LncRNA的gff文件"=>"$dir/LncRNA/Identify/LncRNA.gff",
		"已知lncRNA鉴定统计表"=>"$dir/Personality/Known_LncRNA/known_lncRNA.result.txt",
		"miRNA前体的lncRNA统计"=>"$dir/Personality/precursor/lncRNA_precursor.txt",
		"miRNA靶向lncRNA分析表"=>"$dir/Personality/miRNA_Target2LncRNA/*.mir2target.list",
		"差异基因转录因子分析表"=>"$dir/Personality/TF/TF.txt",
	);
	return $table_info{$id};
}




sub picture_write {
	my ($name,$type,$desc) = @_;
	#$desc =~ s/\s+$//;
	my @picts = glob(&get_picture($desc));
	my $pict_num = @picts;
	my $content;
	$type="type1";
	$type = &T($type);
	if (1<$pict_num) {
		$name = &T($name);
		$desc = &T($desc);
		$content = "<pic_list	name=$desc	type=$type > \n";
		foreach my $pict(@picts) {
			my $pict_name = &T (basename($pict));
			my $path = &T($pict);
			$content = $content."<pic	name=$pict_name	path=$path \/> \n";
		}
		$content = $content.'</pic_list>'."\n";
	}
	if($pict_num==1) {
		my $path=$picts[0];
		if (-f "$path"){
			if ($desc=~/长链非编码测序实验流程图/){
				my $pic_name=&T($desc);
				$path = &T($path);
				$type="img-width-normal";
				$type = &T($type);
				$content = "<pic name=$pic_name type=$type path=$path \/>\n";
		}
		else{
			my $pict_name = &T($desc);
			$path = &T($path);
			$content="<pic name=$pict_name type=$type path=$path \/> \n";
			}	
		}
	}
	return $content;
}



sub get_picture {
	my $name=shift;
	%pic_info = (
		"长链非编码测序实验流程图" => "$dir_template/P01_RNA-Seq_experimental_workflow.png",
		"生物信息分析流程图（mRNA部分）"=> "$dir_template/P02_RNA-Seq_analysis_workflow.png",
		"长链非编码RNA生物信息分析流程图" => "$dir_template/P02_LncRNA_analysis_workflow.png",
		"FASTQ格式文件示意图" => "$dir_template/P03_FASTQ_format.png",
		"碱基质量值分布图" => "$dir_abs/QC/Clean_Data/PNG/*.quality.png",
		"ATGC含量分布图" => "$dir_abs/QC/Clean_Data/PNG/*.acgtn.png",
		"公式1 质量值计算公式" => "$dir_template/F01_Qscore_formula.png",
		"公式2 FPKM计算公式" => "$dir_template/F02_FPKM_formula.png",
		"TopHat2分析流程" => "$dir_template/TopHat2_workflow.png",
		"Mapped Reads在参考基因组上的位置及覆盖深度分布图" => "$dir_abs/QC/Map_assess/L??.png",
		"基因组不同区域Reads分布图" => "$dir_abs/QC/Map_assess/*.type.png",
		"IGV浏览器界面" => "$dir_template/P07_IGV_interface.png",
		"Mapped Reads在mRNA上的位置分布图" => "$dir_abs/QC/Map_assess/Total.randcheck.png",
		"插入片段长度模拟分布图" => "$dir_abs/QC/Map_assess/*.insertSize.r.png",
		"长链非编码测序数据饱和度模拟图" => "$dir_abs/QC/Map_assess/*.Saturation.png",
		#"长链非编码测序数据饱和度模拟图"=>"$dir_abs/QC/Map_assess/Total.gene_tag.png",
		"SNP突变类型分布图" => "$dir_abs/mRNA/Structure/SNP/All.snp.type.png",
		"SNP密度分布图" => "$dir_abs/mRNA/Structure/SNP/AllSample.SNP_density.png",
		"SNP注释分类图" => "$dir_abs/mRNA/Structure/SNP/snp_anno/all.anno.stat.png",
		"InDel注释分类图" => "$dir_abs/mRNA/Structure/SNP/indel_anno/all.anno.stat.png",
		"可变剪切类型统计图" => "$dir_abs/mRNA/Structure/Alt_splice/*.png",
		"各样品FPKM密度分布对比图" => "$dir_abs/mRNA/Expression/all.fpkm_density.png",
		"各样品FPKM箱线图" => "$dir_abs/mRNA/Expression/all.fpkm_box.png",
		"样品间相关性图" => "$dir_abs/mRNA/Expression/sample_cluster.png",
		"差异表达火山图" => "$dir_abs/mRNA/DEG/*_vs_*/*_vs_*.FC_FDR.png",
		"差异表达MA图" => "$dir_abs/mRNA/DEG/*_vs_*/*_vs_*.FC_count.png",
		"差异表达基因维恩图" => "$dir_abs/mRNA/DEG/All_DEG_veen.png",
		"差异表达基因聚类图" => "$dir_abs/mRNA/DEG/All_DEG/all_sample_DEG_cluster.png",
		"差异表达基因聚类折线图" => "$dir_abs/mRNA/DEG/k-means.png",
		"差异表达基因GO注释分类统计图" => "$dir_abs/mRNA/DEG/*_vs_*/go_enrichment/*_vs_*.GO.png",
		"差异表达基因topGO富集有向无环图" => "$dir_abs/mRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.png",
		"差异表达基因COG注释分类统计图" => "$dir_abs/mRNA/DEG/*_vs_*/Cog_Anno/*_vs_*.Cog.classfy.png",
		"差异表达基因KEGG分类图" => "$dir_abs/mRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.png",
		"差异表达基因的KEGG通路注释图" =>"$dir_abs/mRNA/DEG/*_vs_*/pathway/kegg_map/*.png",
		"差异表达基因KEGG通路富集散点图" => "$dir_abs/mRNA/DEG/*_vs_*/Graph/*_vs_*.KEGG.Phase.png",
		"差异表达基因蛋白质互作网络图" => "$dir_template/P26_pp_network.png",
		"预测长链非编码RNA统计图" => "$dir_abs/LncRNA/Identify/lnc.class_code.stat.png",
		"预测方法维恩图" => "$dir_abs/LncRNA/Identify/venn.png",
		"DEU分析结果示意图"=>"$dir_template/DEU.png",
		"差异表达lncRNA聚类图" => "$dir_abs/LncRNA/DEG/All_DEG/all_sample_DEG_cluster.png",
		"差异表达lncRNA靶基因聚类折线图"=>"$dir_abs/LncRNA/DEG/k-means.png",
		"差异表达lncRNA靶基因GO注释分类统计图" => "$dir_abs/LncRNA/DEG/*_vs_*/go_enrichment/*_vs_*.GO.png",
		"差异表达lnc靶基因topGO有向无环图" => "$dir_abs/LncRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.png",
		"差异表达lncRNA靶基因COG注释分类统计图" => "$dir_abs/LncRNA/DEG/*_vs_*/Cog_Anno/*_vs_*.Cog.classfy.png",
		"差异表达lncRNA靶基因KEGG分类图" => "$dir_abs/LncRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.png",
		"差异表达lncRNA靶基因KEGG通路富集散点图" => "$dir_abs/LncRNA/DEG/*_vs_*/Graph/*_vs_*.KEGG.Phase.png",
		"差异表达lncRNA靶基因的KEGG通路注释图"=> "$dir_abs/LncRNA/DEG/*_vs_*/pathway/kegg_map/*.png",
		"差异表达lncRNA靶基因蛋白互作网络图" => "$dir_template/P26_pp_network_lncRNA.png",
		"mRNA长度统计图" => "$dir_abs/LncRNA/Compare/Length/mRNA.len.png",
		"lncRNA长度统计图" => "$dir_abs/LncRNA/Compare/Length/lncRNA.len.png",
		"mRNA的exon个数统计图" => "$dir_abs/LncRNA/Compare/Exon/mRNA.exon.png",
		"lncRNA的exon个数统计图" => "$dir_abs/LncRNA/Compare/Exon/lncRNA.exon.png",
		"mRNA的orf长度统计图" => "$dir_abs/LncRNA/Compare/ORF/gene.orf.png",
		"lncRNA的orf长度统计图" => "$dir_abs/LncRNA/Compare/ORF/lnc_filter_final.orf.png",
		"lncRNA和mRNA表达量比较图" => "$dir_abs/LncRNA/Compare/Expression/lnc_vs_mRNA.fpkm.png",
		"lncRNA和mRNA可变剪切异构体比较图"=>"$dir_abs/LncRNA/Compare/Isoform/lnc_vs_mRNA.isoform.png",
		"code区和lncRNA的序列保守性分析图"=>"$dir_abs/Personality/Conservation/phastcons.cdf.png",
		"lncRNA和mRNA位点保守型示意图"=>"$dir_template/Conservation_position.png",
		"Jensen–Shannon divergence公式"=>"$dir_template/JS.png",
		"lncRNA组织特异性表达图"=>"$dir_abs/Personality/Tissue_specific/*.png",
		"差异基因主成分分析图"=>"$dir_abs/Personality/PCA/*.png",
		"mRNA与lncRNA热图"=>"$dir_abs/Personality/WGCNA/*.png",
		
	);
	return $pic_info{$name};
}

=c
sub Link{
	my ($format,$type,$text)=@_;
	my @KEGG=glob(&get_picture($text));
	$type=&T($type);
	$text=&T($text);
	my $content="<pic_list name=$text  type=$type >\n";
	foreach my $k(@KEGG){
		my $name=basename($k);
		$k=~/(^\S+)\./; ####   ####  ####
		my $pat=$1;
		$pat="$pat.html";
		$pat=&T($pat);
		$k=&T($k);
		$name=&T($name);
		$content .= "<pic name=$name path=$k\/>\n";
	}
	$content .= "</pic_list>\n";
	return $content;
}
=cut

sub Head{
	my ($name,$type,$desc) = @_;
	my $class;
	if ($name eq '一级标题'){
		$class='h1';
	}
	elsif ($name eq '二级标题'){
		$class = 'h2';
	}
	elsif ($name eq '三级标题') {
		$class = 'h3';
	}
	elsif ($name eq '四级标题') {
		$class = 'h4';
	}
	
	$type = &T($type);
	$desc = &T($desc);
	my $content = "<$class name=$desc type=$type desc=$desc \/> \n";
	return $content;
}




sub Text{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "<p type=$type desc=$desc \/> \n";
	return $content;
}



sub Anno{
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "<p type=$type desc=$desc \/> \n";
	return $content;
}









sub tab2lines {
	my $text =shift;
	my $line;
	my @lines =split(/\t/,$text);
	foreach my $each (@lines) {
		$each = $each.'<br>';
		$line= $line.$each
	}
	return $line;
}




sub Reference{
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	my @lines = split /\n/,$desc;
	foreach my $line (@lines) {
		my @each = split /\t/,$line;
		$Ref_info{$each[0]} = $each[1];
	}
	my $content = "<ref_list name=\"参考文献\" type=$type > \n";
	foreach my $id (sort {$a<=>$b} keys %Ref_info) {
		my $ref_name = $Ref_info{$id};
		my $ref_lin = $reference{$id};
		$id = &T($id);
		$ref_name = &T($ref_name);
		$ref_lin=&T($ref_lin);
		$content = $content."<ref id=$id name=$ref_name link=$ref_lin \/> \n";
	}
	$content = $content."</ref_list>\n";
	return $content;
	
}




sub Attach {
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	my $content = "<Attach name=\"附录标题\" type=$type  \/> \n";
	return $content;
}






sub Trans{
	my $text = shift;
	my (@data_among_word) = $text =~ /(\$[a-z,A-Z,_,0-9]+)/g;
	if (@data_among_word !=0) {
		for (my $j=0;$j<@data_among_word;$j++) {
			my $data_id = $data_among_word[$j];
			$data_among_word[$j] =~ s/\$/\\\$/;
			$data_id =~ s/\$//;
			$text =~ s/$data_among_word[$j]/$data_among_words{$data_among_word[$j]}/;
		}
	}
	my (@ref_among_word) = $text =~ /\[(\d*)\]/g;
	if (@ref_among_word != 0) {
		for (my $j=0;$j<@ref_among_word;$j++) {
			my $ref=$ref_among_word[$j];
			my $ref_link = &href($ref);
			$ref= quotemeta($ref);
			$text =~ s/\[$ref\]/$ref_link/;
		}
	}
	return $text;
}



sub href {
	my $text =shift;
	my $name = '&lt;a href=&quot;#ref'.$text.'&quot;&gt;['.$text.']&lt;/a&gt;';
	return $name;
}



sub get_data_among_words {
	
	%data_among_words = (
		"\\\$sample_num" => format_figure($sample_num),
		"\\\$finish_date" => $finish_data,
		"\\\$base_num" => format_figure($base_num)."Gb",
		"\\\$min_data" => format_figure($min_data)."Gb",
		"\\\$min_Q30" => $min_Q30,
		"\\\$min_map_ratio" => $min_map_ratio,
		"\\\$max_map_ratio" => $max_map_ratio,
		"\\\$new_gene_num" => format_figure($new_gene_num),
		"\\\$new_gene_ann_num" => format_figure($new_gene_ann_num),
		"\\\$deg_num" => format_figure($deg_num),
		"\\\$lnc_num" => format_figure($lnc_num),
		"\\\$diff_lnc" => format_figure($diff_lnc),
		"\\\$optimized_gene_num" => format_figure($optimized_gene_num),
		"\\\$min" => format_figure($min),
		"\\\$test_time" => format_figure($test_time),
		"\\\$info_time" => format_figure($info_time),
		"\\\$start_time"=>$start_time,
		"\\\$ref_gene_addr"=>$ref_gene_addr,
		"\\\$ref_gene_name"=> $ref_gene_name,
	)
	
	
}


sub table_html{
	my($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	my @inputs=();
	$/="\n";
	open (IN,$input) || die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @tmp=split /\t+/;
		push @inputs,\@tmp;
	}
	$/="》";
	my $titles=basename $input;
	open HH,">$outHtml.tmp" or die $!;
	print HH <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtmll/DTD/xhtmll-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/htm15shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width,initial-scale=1" name="viewport"></meta>
	<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
	<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
	<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/nav.js" type="text/javascript"></script>
	<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
	<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
	<script src="$srcPath/js/ready.js" type="text/javascript"></script>
	<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
	</head>
	<body>
	<div class="container shadow"><header><img src="$srcPath/images/logo.jpg" class="pull-right" />
	<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
	</header>
	</div>
HTML
	if($text){
		print HH "<div class=\"table-responsive\"><p></div>\n";
	}
	print HH "<div class=\"table-responsive\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info\">\n";
	for (my $i=0;$i<=$#{$inputs[0]};$i++){
		print HH "<th>$inputs[0][$i]</th>\n";
	}
	print HH "</tr></thead>\n<tbody>\n";
	for (my $k=1;$k<=$#inputs;$k++){
		print HH "<tr>";
		for (my $i=0;$i<=$#{$inputs[$k]};$i++){
			if($linkColNum){
				my $j=$i+1;
				if($linkColNum=~/,$j,/){
					print HH "<td><a href=\"$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
				}
				else{
					print HH "<td>$inputs[$k][$i]</td>";
				}
			}
			else{
				print HH "<td>$inputs[$k][$i]</td>";
			}
		}
		print HH "</tr>\n";
	}
print HH <<XGL;
	</tbody>
	</table>
	</body>
	</html>
XGL
	close HH;
	`iconv -f "GB2312" -t "UTF-8" $outHtml.tmp -o $outHtml`;
	`rm -r $outHtml.tmp`;
}

sub format_figure {
	my $figure = shift;
	$figure =~ s/\r//;
	$figure =~ s/%//;
	if (!defined $figure) {
		die;
	}
	if ($figure =~ /\./) {
		if ($figure ==100) {
			$figure=100;
		}else {
			$figure = sprintf ("%.2f",$figure);
		}
	}else {
		$figure = Integer_Three_Digit($figure);
	}
	return $figure;
}



sub Integer_Three_Digit {
	my $integer = shift;
	$integer =~ s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $integer;
}


	
sub T{
	my $id=shift;
	$id = '"'.$id.'"';
	return $id;
}



sub substitute {
	my $text = shift;
	$text =~ s/&/&amp;/g;
	$text =~ s/'/&apos;/g;
	$text =~ s/"/&quot;/g;
	$text =~ s/</&lt;/g;
	$text =~ s/>/&gt;/g;
	
	return $text;
}



sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d年%02d月%02d日", $year+1900, $mon+1, $day);
}

sub time_format {
	my $time = shift;
	my ($year,$mon,$day) = (split /\//,$time);
	return (sprintf("%4d年%02d月%02d日", $year, $mon, $day));
}


sub USAGE {
	my $usage =<< "USAGE";
Program: $Script
Version: $version
Usage:(only cufflinks assembly)
	Options:
	-indir		<dir>
	-config	<config.txt>
	-xml	<xml report>
	-h			help
USAGE
	print $usage;
	exit;
}
