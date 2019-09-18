#!/usr/bin/perl -w
use autodie qw(:all);
use strict;
use warnings;
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
`rm -r $dir/Template` if (-d "$dir/Template");
system "cp -r $Bin/Template $dir";

my $type="type1";
my $desc_anno="";
$desc_anno=&T($desc_anno);
my $table_type="full";
$table_type=&T($table_type);

my $dir_abs = abs_path ($dir);
`rm -r $dir/HTML` if (-d "$dir/HTML");
mkdir "$dir/HTML" unless -d "$dir/HTML";
my $HTML="$dir_abs/HTML";
my $dir_template = "$dir_abs/Template";

my %data;
my %table_info;
my %pic_info;
my %data_among_words;
my %Ref_info;
my %config;
my %file_all_info; ##2016-03-03
my %reference;
my %second_links;

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

open (IN,"$config") or die $!;
while (<IN>) {
	chomp;
	s/\s+$//;s/\r$//;
	next if (/^#/ || /^$/);
	my ($k,$v) = (split /\s+/,$_)[0,1];
	if ($k =~ /^Project_name/){
		$config{$k}=$v;
	}
    if ($k =~ /^Project_id/){
		$config{$k}=$v;
	}
	if($k=~/^Sep/){
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
    if ($k=~/^Ref_seq/){
		$config{$k}=$v;
	}
}
close IN;

my $finish_data;
if ( !defined $config{finish_data}) {
    $finish_data = &GetDate;
}else{
	$finish_data=&time_format($config{finish_data});
}

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
my $min;
if ( !defined $config{min}) {
    $min = "XXXX";
}else{
    $min = &time_format($config{min});
}
my $test_time;
if ( !defined $config{test_time}) {
    $test_time = "XXXX/XX/XX";
}else{
    $test_time = &time_format($config{test_time});
}
my $info_time;
if ( !defined $config{info_time}) {
    $info_time = "XXXX/XX/XX";
}else{
    $info_time = &time_format($config{info_time});
}
my $start_time;
if ( !defined $config{start_time}) {
    $start_time = "XXXX/XX/XX";
}else{
    $start_time = &time_format($config{start_time});
}
my $ref_gene_addr;
if ( !defined $config{ref_gene_addr}) {
    my $tem1 = $config{Ref_seq};
#    print "==========$config{Ref_seq}===================\n";
    my ($Ref_seq1,$Ref_seq2) = (split /\//, $tem1)[6,7];#   split /\n/, $_;
    $ref_gene_addr = "$Ref_seq1"."."."$Ref_seq2";
}else{
    $ref_gene_addr = $config{ref_gene_addr};
}

my $ref_gene_name;
if ( !defined $config{ref_gene_name}) {
    my $Ref_seq = (split /\//,$config{Ref_seq})[6];
    $ref_gene_name = $Ref_seq;
}else{
    $ref_gene_name = $config{ref_gene_name};
}
&get_data_among_words();
#&get_picture();


my $ftemplate="$Bin/1.Basic_tem.v4.txt";
my @ftemplates;
push @ftemplates,$ftemplate;
if (-d "$dir/BMK_3_mRNA/BMK_7_DEU/"){
	$ftemplate="$Bin/2.DEU_tem.txt";
	push @ftemplates,$ftemplate;
}
my @persona=glob("$dir/BMK_4_Personality/*");
foreach my $per_dir (@persona){
  my $name=basename($per_dir);
  if($name=~/Conservation/){
	$ftemplate="$Bin/4.conservation.txt";
	push @ftemplates,$ftemplate;
	}
	if ($name=~/Known_LncRNA/){
		$ftemplate="$Bin/6.Known_lncRNA.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/Tissue_specific/){
		$ftemplate="$Bin/5.tissue_specify.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/Precursor_Result/){
		$ftemplate="$Bin/7.precusor.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/miRNA_Target2LncRNA/){
		$ftemplate="$Bin/8.miRNA_target2lncRNA.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/PCA/){
		$ftemplate="$Bin/10.PCA.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/TF/){
		$ftemplate="$Bin/9.TF.txt";
		push @ftemplates,$ftemplate;
	}
	if ($name=~/CircRNA_Analysis/){
		$ftemplate="$Bin/13.circRNA.txt";
		push @ftemplates,$ftemplate;
	}
}
if ($sample_num >= 4) {
		$ftemplate="$Bin/11.WGCNA.txt";
		push @ftemplates,$ftemplate;
}
$ftemplate="$Bin/3.attach_tem.txt";
push @ftemplates,$ftemplate;

`cat @ftemplates >$dir/Template/template`;
$ftemplate="$dir/Template/template";
print "$ftemplate\n";



chomp (my $user = `whoami`);
my $user_addr = $user."\@biomarker.com.cn";
my $report_version = &T ("v1.2");
my $report_name = &substitute($prefix);
$report_name=&T($prefix);
$report_name = &Trans($report_name);
my $report_code = $config{Project_id};#&T("XXX");
$report_code=&T($report_code);
$report_code = &Trans($report_code);
my $report_user = &T($user);
my $report_user_addr = &T($user_addr);
my $report_time = &GetDate;#&time_format(&GetTime);
my $alltime =join(";",$test_time,$info_time,$start_time,$report_time);
$alltime=&T($alltime);
$alltime=&Trans($alltime);
$report_time = &T($report_time);




# -------------------------------------------------------------------------
# output xml
# -------------------------------------------------------------------------
open (OUT,">$xml") or die $!;

print OUT "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
print OUT "<report>\n";
print OUT "\t<report_version value=$report_version \/>\n";
print OUT "\t<report_name value=$report_name \/>\n";
print OUT "\t<report_code value=$report_code \/>\n";
print OUT "\t<report_user value=$report_user \/>\n";
print OUT "\t<report_user_addr value=$report_user_addr \/>\n";
print OUT "\t<report_time value=$alltime \/>\n";
#print OUT "\t<report_abstract value='' \/>\n";



open (TEM, $ftemplate) or die $!;
$/="》";
my $i=0;
while (<TEM>){
	chomp;
	next if (/^$/ || /^\#/);
	my (undef ,$format,@context)=split /\n/,$_;;
	$i++;
	my $num=@context;
	print "$.\t$context[0]\n";
	if ($num > 1){
	  #print Dumper(\@context);
	  if ($format =~ /图片/ ){print OUT &picture_write($format,$type,$context[0],$context[1]);}
	  elsif ($format =~ /表格/){print OUT &table_write($format,$type,$context[0],$context[1]); }
	  elsif($format =~ /正文/){print OUT &Text_mul($format,$type,\@context);}
	  elsif($format =~ /参考文献/){print OUT &Reference($format,$type,\@context);}
	}else{
	  &xml_write($format,$type,$context[0]);
	}
	$format =~ s/\s+$//;
}
print OUT "<\/report>\n";

close TEM;
close OUT;
$/="\n";    ### 


my @DEU_file;
my $file;
my $cmd;
if(-d "$dir/BMK_3_mRNA/BMK_7_DEU/"){
    @DEU_file = glob("$dir/BMK_3_mRNA/BMK_*_DEU/*_vs_*/DEXSeqReport/testForDEU.html");
    for $file (@DEU_file) {
        $cmd = "sed 's/href=\"/href=\"project_DEXSeqReport_Path\\//' $file > $file.cloud ";
        system "$cmd";
    }
}



sub xml_write {
	my ($format,$type,$text) = @_;
	if ($format eq "表格") {### 表格
		print OUT &table_write($format,$type,$text);
    }
	elsif($format=~/文件集合/) {## 公式
        print OUT &file_all_write($format,$text);
	}
	elsif ($format =~ /图片/) {### 图片
		print OUT &picture_write($format,$type,$text);
	}
	elsif ($format =~ /分析结果概述/) { ## 分析结果概述
		print OUT &abstract($format,$type,$text);
	}
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
	my ($name,$type,$desc,$anno) = @_;
	my @tables = glob (&get_table($desc));
	my $content;
	$type = "xls"; ###linhj
	$type = &T($type);
	
	$anno=&T($anno);
	#chomp $desc;
	my $table_num = @tables;
	my $txt=$tables[0];
	my $txt_line=`less -S $txt|wc -l`; chomp $txt_line;
	next if ($table_num ==0);
	$txt=$tables[1] unless ($txt_line !=1);
	if ($desc!~/DEU分析结果$/ ){
	  $desc = &T($desc);
=pod
	  if ($txt_line < 20 && $table_num==1){
		$txt=&T($txt);
		$content="\t<table name=$desc type=$table_type path=$txt desc=$anno \/>\n";
	  }
=cut
		if ($txt_line < 20 &&  $txt_line >1){
			
			if ($desc=~/碱基质量值与碱基识别出错的概率的对应关系表|附表/){
				$txt=&T($txt);
				$content="\t<table name=$desc type=$table_type path=$txt desc=$anno \/>\n";
			}else{
			  #$txt=&T($txt);
			  print "$desc\n";
				my $table_name=basename($txt);
				my $txt_tmp=shift;
				if(!-f "$dir/Template/$table_name"){
					`head -n 6 $txt >$dir/Template/$table_name`;
					$txt_tmp="$dir/Template/$table_name";
					}else{
						`head -n 6 $txt >$dir/Template/$table_name.tmp`;
						$txt_tmp="$dir/Template/$table_name.tmp";
					}
					$txt_tmp=&T($txt_tmp);
				  $content = "\t<table name=$desc type=$table_type path=$txt_tmp desc=$anno \/>\n";
				 $content .="\t<file_list name=$desc type=$type desc=\"\" >\n";
					foreach my $path (@tables){
					my $title=&Title_write($path);
					if (!-f "$dir/HTML/$title.html"){
					  &table_html($path,"$dir/HTML/$title.html",$desc,"","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
					  &table_html_cloud($path,"$dir/HTML/$title.html.cloud",$desc,"","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
					  $path="$dir/HTML/$title.html" if ($path!~/html$/);
					}else{
					  &table_html_cloud($path,"$dir/HTML/$title.tmp.html.cloud",$desc,"","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
					  &table_html($path,"$dir/HTML/$title.tmp.html",$desc,"","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
					  $path="$dir/HTML/$title.tmp.html" if ($path!~/html$/);
					}
					my $html_title=basename($path);
					$content .="\t\t".'<file name="'.$html_title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
				  }
				  $content .="\t".'</file_list>'."\n";
		  
			}
		}
		elsif ($txt_line>=20){
		  #$txt=&T($txt);
		  my $table_name=&Title_write($txt);
		  my $txt_tmp=shift;
		  if(!-f "$dir/Template/$table_name"){
			`head -n 6 $txt >$dir/Template/$table_name`;
			$txt_tmp="$dir/Template/$table_name";
		  }else{
			`head -n 6 $txt >$dir/Template/$table_name.tmp`;
			$txt_tmp="$dir/Template/$table_name.tmp";
		  }
		  #$txt_tmp=&T($txt_tmp);
		  if ($desc=~/Scripture拼接结果/){
			my $tit="Scripture拼接"; $tit=&T($tit);
			my $tit_class="h3";
			my $scripture="Scripture 拼接基于统计学分段模型区分表达位点和实验背景噪音，比较适用于长转录本的拼接。Scripture 拼接结果如下表：";
			$scripture=&substitute($scripture); $scripture=&T($scripture);
			$content ="\t<$tit_class name=$tit type=$type desc=$tit \/>\n";
			$content .="\t<p type=\"type1\" desc=$scripture \/>\n";
			$txt_tmp=&T($txt_tmp);
			$content .= "\t<table name=$desc type=$table_type path=$txt_tmp desc=$anno \/>\n";
			$content .="\t<file_list name=$desc type=$type desc=\"\" >\n";
			foreach my $path (@tables){
				my $title=&Title_write($path);
				#print "niulg 0823 modify $title\n";
				  &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$desc","","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
				  &table_html($path,"$dir/HTML/$title.html","$desc","","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
				  $path="$dir/HTML/$title.html" if ($path!~/html$/);
				my $html_title=basename($path);
				$content .="\t\t".'<file name="'.$html_title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
			}
		  }else{
			$txt_tmp=&T($txt_tmp);
			$content = "\t<table name=$desc type=$table_type path=$txt_tmp desc=$anno \/>\n";
			$content .="\t<file_list name=$desc type=$type desc=\"\" >\n";
			foreach my $path (@tables){
				my $title=&Title_write($path);
				if (!-f "$dir/HTML/$title.html"){
				  &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$desc","","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
				  &table_html($path,"$dir/HTML/$title.html","$desc","","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
				  $path="$dir/HTML/$title.html" if ($path!~/html$/);
				}else{
				  &table_html_cloud($path,"$dir/HTML/$title.tmp.html.cloud","$desc","","","project_Template_Path/../src",$anno,"logo_image_path") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
					&table_html($path,"$dir/HTML/$title.tmp.html","$desc","","","../src/",$anno) if ($path!~/html$/ and $path!~/KEGG\.list/);
					$path="$dir/HTML/$title.tmp.html" if ($path!~/html$/);
				}
				my $html_title=basename($path);
				$content .="\t\t".'<file name="'.$html_title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
			}
		  }
		  $content .="\t".'</file_list>'."\n";
		  
		}
	}  
  
	elsif($desc=~/DEU分析结果$/){
	  $desc = &T($desc);
			$content="\t<file_list name=$desc type=$type desc=$desc_anno >\n";
			foreach my $table(@tables){
				$table=~/\/(\w+_vs_\w+)\//;
				my $table_name=$1;
				$table_name=$table_name.".DEU.html";
				$table_name=&T($table_name);
				my $path=&T($table);
				my $action=&T("type1");
				$content .= "\t\t<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .= "\t".'</file_list>'."\n" ;
		}
	return $content;
  }




sub get_table {
	chomp (my $id = shift);
	%table_info=(
		"碱基质量值与碱基识别出错的概率的对应关系表" => "$dir_template/error.txt",
		"样品测序数据评估统计表" => "$dir/BMK_1_rawData/BMK_1_Data_Assess/AllSample_GC_Q.stat",
		"样品测序数据与所选参考基因组的序列比对结果统计表" => "$dir/BMK_*_rawData/BMK_*_Mapped_Statistics/Total_mapped.stat",
		"SNP位点信息" => "$dir/BMK_*_mRNA/BMK_*_SNP_Analysis/final.SNP.anno.gatk.all.list",
		"InDel位点信息" => "$dir/BMK_*_mRNA/BMK_*_SNP_Analysis/final.InDel.anno.gatk.all.list",
		"SNP位点统计表" => "$dir/BMK_*_mRNA/BMK_*_SNP_Analysis/AllSample.SNP.stat",
		"可变剪切事件结构和表达量统计表" => "$dir/BMK_*_mRNA/BMK_*_Alt_splice/*.AS.list",
		"基因结构优化结果" => "$dir/BMK_*_mRNA/BMK_*_Gene_Structure_Optimize/*.geneStructure.optimize.xls",
		"新基因的GFF文件" => "$dir/BMK_3_mRNA/BMK_1_NewGene/*.newGene.gff",
		"新基因序列FASTA文件" => "$dir/BMK_3_mRNA/BMK_*_NewGene/*.newGene.longest_transcript.fa",
		"新基因功能注释结果统计" => "$dir/BMK_3_mRNA/BMK_1_NewGene/BMK_1_NewGene_Anno/Function_Annotation.stat.xls",
		"基因表达量结果文件" => "$dir/BMK_3_mRNA/BMK_2_geneExpression/All_gene_fpkm.list",
		"差异表达分析结果" => "$dir/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*vs*/BMK_1_Statistics_Visualization/*_vs_*.DEG_final.xls",
		"差异表达基因数目统计表" => "$dir/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_All_DEG/DEG.stat",
		"注释的差异表达基因数量统计表" => "$dir/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_All_DEG/DEG.anno.stat",
		"差异表达基因topGO富集结果" => "$dir/BMK_3_mRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_*_GO_Enrichment/*_vs_*.topGO_*.xls",
		"差异表达基因的KEGG富集结果" => "$dir/ BMK_3_mRNA/BMK_*_DEG_Analysis/*_vs*/BMK_4_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.stat",
		#"差异表达基因的KEGG通路注释" => "$dir/DEG_Analysis/*_vs_*/pathway/kegg_map/*.html",
		"DEU分析结果表" => "$dir/BMK_3_mRNA/BMK_*_DEU/*_vs_*/*DEU_Result_Final.xls",
		"DEU分析结果"=>"$dir/BMK_3_mRNA/BMK_*_DEU/*_vs_*/DEXSeqReport/testForDEU.html",
		"Cufflinks拼接结果" => "$dir/BMK_2_LncRNA/BMK_1_Assembly_Result/*.Cufflinks.transcripts.gtf",
		"Scripture拼接结果" => "$dir/BMK_2_LncRNA/BMK_1_Assembly_Result/*.Scripture.transcripts.gtf",
		"CPC分析结果统计" => "$dir/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CPC.txt",
		"CNCI分析结果统计" => "$dir/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CNCI.txt",
		"pfam分析结果统计" => "$dir/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/Pfam.txt",
		"CPAT分析结果统计" => "$dir/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/CPAT.txt",
		"基于位置关系的LncRNA靶基因预测结果" => "$dir/BMK_2_LncRNA/BMK_4_LncRNA_Target/LncRNA_Cis_target.txt",
		"基于互补序列的LncRNA靶基因预测结果"=>"$dir/BMK_2_LncRNA/BMK_4_LncRNA_Target/LncRNA_LncTar_target.txt", #if (-f "$dir/LncRNA/lncRNA_basepair.target"),
		"lncRNA表达量结果" => "$dir/BMK_2_LncRNA/BMK_3_LncRNA_Expression/LncRNA_fpkm.list",
		"差异表达lncRNA结果" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.DEG.xls",
		"差异表达lncRNA数目统计表" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/DEG_lncRNA.stat",
		"注释的差异表达lncRNA顺式靶基因数量统计表" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/DEG_lncRNA.cis_anno.stat",
		"注释的差异表达lncRNA反式靶基因数量统计表" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/DEG_lncRNA.trans_anno.stat",
		"差异表达lncRNA顺式靶基因topGO富集结果" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.topGO_*.xls",
		"差异表达lncRNA反式靶基因topGO富集结果" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_*_Trans_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.topGO_*.xls",
		"差异表达lncRNA顺式靶基因KEGG富集结果" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_*_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.stat",
		"差异表达lncRNA反式靶基因KEGG富集结果" => "$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/*_vs_*/BMK_*_Trans_Anno_enrichment/BMK_*_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.stat",
		
		#"差异表达lncRNA靶基因的KEGG通路注释"=>"$dir/Lnc_Diff_Analysis/*_vs_*/pathway/kegg_map/*.html",
		"附表1 软件列表" => "$dir_template/software_used_list.txt",
		"附表2 数据库列表" => "$dir_template/database_used_list.txt",
		"附表3 核酸编码表" => "$dir_template/Nucleic_acids_encoding_table.txt",
		"LncRNA的gff文件"=>"$dir/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/LncRNA.gtf",
		"已知lncRNA鉴定统计表"=>"$dir/BMK_4_Personality/*_Known_lncRNA/known_lncRNA.result.txt",
		"miRNA前体的lncRNA统计"=>"$dir/BMK_4_Personality/*_Precursor_Result/lncRNA_precursor.txt",
		"miRNA靶向lncRNA分析表"=>"$dir/BMK_4_Personality/*_miRNA_Target2LncRNA/lncRNA_target2mirna.list",
		"差异基因转录因子分析表"=>"$dir/BMK_4_Personality/BMK_*_TF/TF.txt",
		#"miRNA前体的lncRNA统计示意表"=>"$dir_template/precursor.txt",
		"circRNA鉴定表"=>"$dir/BMK_4_Personality/*_CircRNA_Analysis/All_CircRNA.xls",
		"circRNA表达量统计表" =>"$dir/BMK_4_Personality/*_CircRNA_Analysis/All_gene_counts.list",
		"数据链特异性建库评估信息表"=>"$dir/BMK_1_rawData/BMK_*_Library_Assessment/Lib_assess.txt",
	);
	return $table_info{$id};
    #print Dumper(\%table_info);
}
sub get_file {#
    chomp (my $id = shift);
            ## 文件 差异表达基因的KEGG富集结果（网页）
            # width
            if ($id=~/差异表达基因的KEGG富集网页结果/) {
                my @kegg_enrich=glob("$dir/BMK_*_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Pathway_Enrichment/kegg_enrichment/*.KEGG.list");
               $file_all_info{"差异表达基因的KEGG富集网页结果"} =[@kegg_enrich];
               $/ = "\n";
               foreach my $degs (@kegg_enrich) {
                    my $key_name=basename$degs;
                    my $key_dir=dirname$degs;
                    open IN,$degs ||die $!;
                    while (<IN>) {
                        next if ($_=~/^\#|^\s*$/) ;
                        my ($keywords,$koid)=split(/\t+/,$_);
                        my $mid_path="$key_dir/../kegg_map/$koid.html";
                        $mid_path=~s/$dir\///;#print "bbbbbbbbbbbb\t$mid_path\n";
                        $second_links{KEGG}{$key_name}{$keywords}="../$mid_path";
                    }
                }
                $/ = "》";
                return $file_all_info{$id};
                return $second_links{KEGG};
            }
}
sub get_file_lnc {#
    chomp (my $id = shift);
            if ($id=~/差异表达lncRNA顺式靶基因KEGG富集网页结果/) {
               my @kegg_enrich=glob("$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_*_Pathway_Enrichment/kegg_enrichment/*.KEGG.list");
               $file_all_info{"差异表达lncRNA顺式靶基因KEGG富集网页结果"} =[@kegg_enrich];
               $/ = "\n";
               foreach my $degs (@kegg_enrich) {
                    my $key_name=basename$degs;
                    $key_name="lnc_cis.$key_name";
                    my $key_dir=dirname$degs;
                    open IN,$degs ||die $!;
                    while (<IN>) {
                        next if ($_=~/^\#|^\s*$/) ;
                        my ($keywords,$koid)=split(/\t+/,$_);
                        my $mid_path="$key_dir/../kegg_map/$koid.html";
                        $mid_path=~s/$dir\///;#print "aaaaaaaaaa\t$mid_path\n";
                        $second_links{lnc_cisKEGG}{$key_name}{$keywords}="../$mid_path";
                    }
                }
                $/ = "》";
                return $file_all_info{$id};
                return $second_links{lnc_cisKEGG};
			}
			if ($id=~/差异表达lncRNA反式靶基因KEGG富集网页结果/) {
               my @kegg_enrich=glob("$dir/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Trans_Anno_enrichment/BMK_*_Pathway_Enrichment/kegg_enrichment/*.KEGG.list");
               $file_all_info{"差异表达lncRNA反式靶基因KEGG富集网页结果"} =[@kegg_enrich];
               $/ = "\n";
               foreach my $degs (@kegg_enrich) {
                    my $key_name=basename$degs;
                    $key_name="lnc_trans.$key_name";
                    my $key_dir=dirname$degs;
                    open IN,$degs ||die $!;
                    while (<IN>) {
                        next if ($_=~/^\#|^\s*$/) ;
                        my ($keywords,$koid)=split(/\t+/,$_);
                        my $mid_path="$key_dir/../kegg_map/$koid.html";
                        $mid_path=~s/$dir\///;#print "aaaaaaaaaa\t$mid_path\n";
                        $second_links{lnc_transKEGG}{$key_name}{$keywords}="../$mid_path";
                    }
                }
                $/ = "》";
                return $file_all_info{$id};
                return $second_links{lnc_transKEGG};
            }
            print Dumper(\%second_links);
}




sub picture_write {
	my ($name,$type,$desc,$anno) = @_;
	#$desc =~ s/\s+$//;
	my @picts = glob(&get_picture($desc));
	my $pict_num = @picts;
	my $content;
	$type="type1";
	$type = &T($type);
	$anno=&T($anno);
	if (1<$pict_num) {
		if ($desc=~/原始数据组成/){
			my $tit="测序质量控制"; $tit=&T($tit);
			my $tit_class="h3"; 
			my $ann="在进行数据分析之前，首先需要确保这些Reads有足够高的质量，以保证后续分析的准确。百迈客对数据进行严格的质量控制，进行如下过滤方式：";
#			$ann .="\n"."在进行数据分析之前，首先需要确保这些Reads有足够高的质量，以保证后续分析的准确。百迈客对数据进行严格的质量控制，进行如下过滤方式：";
			$ann .="<br>"."(1) 去除含有接头的Reads；";
			$ann .="<br>"."(2) 去除低质量的Reads（包括去除N的比例大于10%的Reads；去除质量值Q≤10的碱基数占整条Read的50%以上的Reads）。";
			$ann .="<br>"."经过上述一系列的质量控制之后得到的高质量的Clean Data，以FASTQ格式提供。";
			$ann=&substitute($ann); $ann=&T($ann);
			my $ann1="注：Adapter related：过滤掉的含有接头Reads数占总Raw Reads数的比例。Low quality：过滤掉的低质量Reads数占总Raw Reads数的比例。Clean Reads：经过以上过滤得到的Clean Reads 数占总Raw Reads 数的比例。";
			$ann1=&substitute($ann1); $ann1=&T($ann1);
			$content = "\t<$tit_class name=$tit type=$type desc=$tit \/>\n";
			$content .="\t<p type=$type desc=$ann \/>\n";
			$desc=&T($desc);
			$content .="\t<pic_list name=$desc type=$type desc=$ann1 >\n";
			foreach my $pict(@picts){
				my $pict_name=&T(basename($pict));
				my $path=&T($pict);
				$content .="\t\t<pic name=$pict_name desc=$desc_anno path=$path \/>\n";
			}
			$content .= "\t".'</pic_list>'."\n";
		}
		else{
			$desc = &T($desc);
			$content = "\t<pic_list name=$desc type=$type desc=$anno >\n";
			foreach my $pict(@picts) {
				my $pict_name = &T (basename($pict));
				my $path = &T($pict);
				$content = $content."\t\t<pic name=$pict_name desc=$desc_anno path=$path \/>\n";
			}
		$content = $content."\t".'</pic_list>'."\n";
		}
	}
	if($pict_num==1) {
		my $path=$picts[0];
        my $desc_anno1=&T(" ");
		if (-f "$path"){
			if ($desc=~/长链非编码测序实验流程图/){
				my $pic_name=&T($desc);
				$path = &T($path);
				$type="type1";
				$type = &T($type);
				$content = "\t<pic name=$pic_name type=$type desc=$desc_anno1 path=$path \/>\n";
			}
			elsif ($desc=~/差异表达基因维恩图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="将每组差异基因进行维恩图绘制，见下图。图中展示了各比较组特有的差异基因的个数，以及比较组间的共有的差异基因个数。"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $tit="差异表达基因维恩图";	$tit=&T($tit);
				my $tit_class="h3";	
#                my $desc_anno1=&T(" ");
				$content="\t<$tit_class name=$tit type=$type desc=$tit \/>\n";
				$content .="\t<p type=$type desc=$ann \/>\n";
				$content .="\t<pic name=$picture_name type=$type desc=$desc_anno1 path=$path \/>\n";
			}
			elsif ($desc=~/表达基因聚类折线图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="注：x轴表示实验条件，y轴表示标准化的FPKM。黑色线条表示这个cluster中的所有基因在不同实验条件下相对表达量的平均值的折线图。"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="对所有基因的表达量做K-means聚类分析，对基因的表达量水平值log2(FPKM+1) 进行聚类分析,得到基因在不同实验条件下的表达模式，表达模式相同或相近的基因聚集成类。由于同类的基因可能具有相似的功能，或是共同参与同一代谢过程或细胞通路。因此，通过对同一个cluster的基因进行功能注释及富集分析，可预测未知基因的功能或已知基因的未知功能。";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content ="\t<p type=$type desc=$ann1 \/>\n";
				$content .="\t<pic name=$picture_name type=$type desc=$ann path=$path \/>\n";
				#$content .="\t<p type=$type desc=$ann \/>\n";
				
			}
			elsif ($desc=~/差异表达lncRNA聚类图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="对筛选出的差异表达lncRNA做层次聚类分析，将具有相同或相似表达行为的lncRNA进行聚类，差异表达lncRNA聚类结果如下图：" ; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1 ="注：图中不同的列代表不同的样品，不同的行代表不同的lncRNA。颜色代表了lncRNA在样品中的表达量水平log2（FPKM+1）。"; 
				###0818
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				my $tit="差异表达lncRNA聚类分析"; $tit=&T($tit);
				my $tit_class='h3'; 
				$content ="\t<$tit_class name=$tit type=$type desc=$tit \/>\n";
				$content .="\t<p type=$type desc=$ann \/>\n";
				$content .="\t<pic name=$picture_name type=$type desc=$ann1 path=$path \/>\n";
				#$content .="\t<p type=$type desc=$ann1 \/>\n";
			}
			elsif ($desc=~/lncRNA聚类折线图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="注：x轴表示实验条件，y轴表示标准化的FPKM。黑色线条表示这个cluster中的所有基因在不同实验条件下相对表达量的平均值的折线图。";
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="对筛选出的所有lncRNA的表达量做层次聚类分析，得到lncRNA在不同实验条件下的表达模式，表达模式相同或相近的lncRNA聚集成类。见上图：";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content = "\t<pic name=$picture_name type=$type desc=$ann path=$path \/>\n";
				#$content .="\t<p type=$type desc=$ann \/>\n";
				$content .="\t<p type=$type desc=$ann1 \/>\n";
			}
			elsif ($desc=~/样品间相关性图/){
				if (exists $config{'Sep'}){
					my $picture_name=&T($desc);
					my $tit="重复相关性评估"; $tit=&T($tit);
					my $tit_class="h3";
					my $ann="研究表明，基因的表达在不同的个体间存在生物学可变性（Biological Variability），不同的基因之间表达的可变程度存在差异，而长链非编码测序技术、qPCR以及生物芯片等技术都不能消除这种可变性。为了寻找真正感兴趣的差异表达基因，需要考虑和处理因生物学可变性造成的表达差异。目前最常用且最有效的方法是在实验设计中设立生物学重复（Biological Replicates）。重复条件限制越严格，重复样品数目越多，寻找到的差异表达基因越可靠。对于设立生物学重复的项目，评估生物学重复的相关性对于分析转录组测序数据非常重要。生物学重复的相关性不仅可以检验生物学实验操作的可重复性；还可以评估差异表达基因的可靠性和辅助异常样品的筛查。";
					$ann .="将皮尔逊相关系数r（Pearson’s Correlation Coefficient）作为生物学重复相关性的评估指标。r2越接近1，说明两个重复样品相关性越强。见下图";
					my $ann1="注：横坐标表示样品名称，纵坐标表示对应的样品名称。颜色代表r2值大小。";
					###0818
					$ann1=&substitute($ann1); $ann1=&T($ann1);
					$ann=&substitute($ann); $ann=&T($ann);
					$path=&T($path);
					$content="\t<$tit_class name=$tit type=$type desc=$tit \/>\n";
					$content .="\t<p type=$type desc=$ann \/>\n";
					$content .="\t<pic name=$picture_name type=$type desc=$ann1 path=$path \/>\n";
					#$content .="\t<p type=$type desc=$ann1 \/>\n";
				}
			}
			else{
				my $pict_name = &T($desc);
				$path = &T($path);
				$content="\t<pic name=$pict_name type=$type desc=$anno path=$path \/>\n";
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
		"碱基错误率分布图" => "$dir_abs/BMK_1_rawData/BMK_1_Data_Assess/PNG/*.quality.png",
		"ATGC含量分布图" => "$dir_abs/BMK_1_rawData/BMK_1_Data_Assess/PNG/*.acgtn.png",
		"公式1 质量值计算公式" => "$dir_template/F01_Qscore_formula.png",
		"公式2 FPKM计算公式" => "$dir_template/F02_FPKM_formula.png",
		"TopHat2分析流程" => "$dir_template/TopHat2_workflow.png",
		"Mapped Reads在参考基因组上的位置及覆盖深度分布图" => "$dir_abs/BMK_1_rawData/BMK_2_Mapped_Statistics/*.map.png",
		"基因组不同区域Reads分布图" => "$dir_abs/BMK_1_rawData/BMK_2_Mapped_Statistics/*.type.png",
		"IGV浏览器界面" => "$dir_template/P07_IGV_interface.png",
		"Mapped Reads在mRNA上的位置分布图" => "$dir_abs/BMK_1_rawData/BMK_3_Library_Assessment/Total.randcheck.png",
		"插入片段长度模拟分布图" => "$dir_abs/BMK_1_rawData/BMK_3_Library_Assessment/*.insertSize.r.png",
		"长链非编码测序数据饱和度模拟图" => "$dir_abs/BMK_1_rawData/BMK_3_Library_Assessment/*.Saturation.png",
		#"长链非编码测序数据饱和度模拟图"=>"$dir_abs/geneExpression/Total.gene_tag.png",
		"SNP突变类型分布图" => "$dir_abs/BMK_3_mRNA/BMK_5_SNP_Analysis/BMK_3_SNP_type/All.SNP.type.png",
		"SNP密度分布图" => "$dir_abs/BMK_3_mRNA/BMK_5_SNP_Analysis/AllSample.SNP_density.png",
		"SNP注释分类图" => "$dir_abs/BMK_3_mRNA/BMK_5_SNP_Analysis/BMK_2_SNP_anno/all.anno.stat.png",
		"InDel注释分类图" => "$dir_abs/BMK_3_mRNA/BMK_5_SNP_Analysis/BMK_1_InDel_anno/all.anno.stat.png",
		"可变剪切类型统计图" => "$dir_abs/BMK_3_mRNA/BMK_4_Alt_splice/*.png",
		"各样品FPKM密度分布对比图" => "$dir_abs/BMK_3_mRNA/BMK_2_geneExpression/BMK_2_Expression_Statistics/all.fpkm_density.png",
		"各样品FPKM箱线图" => "$dir_abs/BMK_3_mRNA/BMK_2_geneExpression/BMK_2_Expression_Statistics/all.fpkm_box.png",
		"样品间相关性图" => "$dir_abs/BMK_3_mRNA/BMK_2_geneExpression/BMK_1_Sample_Correlation/sample_cluster.png",
		"差异表达火山图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.Volcano.png",
		"差异表达MA图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*_vs_*.MA.png",
		"差异表达基因维恩图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_All_DEG/All_DEG_veen.png",
		"差异表达基因聚类图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/all_sample_DEG_cluster.png",
		"表达基因聚类折线图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/k-means.png",
		"差异表达基因GO注释分类统计图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_3_GO_Enrichment/*_vs_*.GO.png",
		"差异表达基因topGO富集有向无环图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_3_GO_Enrichment/*_vs_*.topGO_*.png",
		"差异表达基因COG注释分类统计图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK*_vs_*/BMK_2_DEG_Annotation/*_vs_*.COG.classfy.png",
		"差异表达基因KEGG分类图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_4_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.png",
		"差异表达基因的KEGG通路注释图" =>"$dir_template/P20_pathway_example.png",
		"差异表达基因KEGG通路富集散点图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_4_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.Phase.png",
		"差异表达基因蛋白质互作网络图" => "$dir_template/P26_pp_network.png",
		"预测长链非编码RNA统计图" => "$dir_abs/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/LncRNA_classification.png",
		"预测方法维恩图" => "$dir_abs/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/BMK_1_Software_Result/venn.png",
		"lncRNA染色体分布环状图" =>"$dir_abs/BMK_2_LncRNA/BMK_2_LncRNA_Prediction/LncRNA_circos.png",
		"DEU分析结果示意图"=>"$dir_template/DEU.png",
		"差异表达lncRNA聚类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/all_sample_DEG_cluster.png",
		"lncRNA聚类折线图"=>"$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_1_All_DEG/k-means.png",
		"差异表达lncRNA顺式靶基因GO注释分类统计图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.GO.png",
		"差异表达lncRNA顺式靶基因topGO有向无环图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.topGO_*.png",
		"差异表达lncRNA顺式靶基因COG注释分类统计图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_1_DEG_Annotation/*_vs_*.COG.classfy.png",
		"差异表达lncRNA顺式靶基因KEGG分类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.png",
		"差异表达lncRNA顺式靶基因KEGG通路富集散点图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK*_vs_*/BMK_*_Cis_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.Phase.png",
		"差异表达lncRNA顺式靶基因的KEGG通路注释图"=> "$dir_template/P19_pathway_example_cis.png",
		"差异表达lncRNA顺式靶基因GO富集聚类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/BMK_*_Cis_Cluster/GO_cluster.png",
		"差异表达lncRNA顺式靶基因KEGG富集聚类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/BMK_*_Cis_Cluster/KEGG_cluster.png",
		"差异表达lncRNA反式靶基因GO注释分类统计图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*Trans_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.GO.png",
		"差异表达lncRNA反式靶基因topGO有向无环图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*Trans_Anno_enrichment/BMK_2_GO_Enrichment/*_vs_*.topGO_*.png",
		"差异表达lncRNA反式靶基因COG注释分类统计图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_vs_*/BMK_*Trans_Anno_enrichment/BMK_1_DEG_Annotation/*_vs_*.COG.classfy.png",
		"差异表达lncRNA反式靶基因KEGG分类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK*_vs_*/BMK_*Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichmentt/*_vs_*.KEGG.png",
		"差异表达lncRNA反式靶基因KEGG通路富集散点图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK*_vs_*/BMK_*Trans_Anno_enrichment/BMK_3_Pathway_Enrichment/kegg_enrichment/*_vs_*.KEGG.Phase.png",
		"差异表达lncRNA反式靶基因的KEGG通路注释图"=> "$dir_template/P18_pathway_example_trans.png",
		"差异表达lncRNA反式靶基因GO富集聚类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/BMK_*_Trans_Cluster/GO_cluster.png",
		"差异表达lncRNA反式靶基因KEGG富集聚类图" => "$dir_abs/BMK_2_LncRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/BMK_*_Trans_Cluster/KEGG_cluster.png",
		"差异表达lncRNA靶基因蛋白互作网络图" => "$dir_template/P26_pp_network_lncRNA.png",
		"mRNA长度统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Length_Compare/mRNA.len.png",
		"lncRNA长度统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Length_Compare/lncRNA.len.png",
		"mRNA的exon个数统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Exon_Compare/mRNA.exon.png",
		"lncRNA的exon个数统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Exon_Compare/lncRNA.exon.png",
		"mRNA的orf长度统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_ORF_Compare/mRNA.ORF.png",
		"lncRNA的orf长度统计图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_ORF_Compare/lncRNA.ORF.png",
		"lncRNA和mRNA表达量比较图" => "$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Expression_Compare/lnc_vs_mRNA.fpkm.png",
		"lncRNA和mRNA可变剪切异构体比较图"=>"$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_Isoform_Compare/lnc_vs_mRNA.isoform.png",
		"差异表达lncRNA和mRNA火山交互式图"=>"$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_DEG_Compare/*.Volcano.png",
		"差异表达lncRNA和mRNAMA交互式图"=>"$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_DEG_Compare/*.MA.png",
		"差异表达lncRNA和mRNA染色体分布环状图"=>"$dir_abs/BMK_*_LncRNA/BMK_*_LncRNA_vs_mRNA/BMK_*_DEG_Compare/*.circos.png",
		"code区和lncRNA的序列保守性分析图"=>"$dir_abs/BMK_*_Personality/BMK_*_Conservation/phastcons.cdf.png",
		"lncRNA和mRNA位点保守型示意图"=>"$dir_template/Conservation_position.png",
		"Jensen_Shannon divergence公式"=>"$dir_template/JS.png",
		"lncRNA特异性表达图"=>"$dir_abs/BMK_*_Personality/BMK_*_Tissue_specific/*.png",
		"差异基因主成分分析图"=>"$dir_abs/Personality/PCA/*.png",
		"mRNA与lncRNA共表达网络分析示意图"=>"$dir_abs/BMK_2_LncRNA/BMK_*_LncRNA_Target/WGCNA_Result/Fig_8_1_networkHeatmap.png",
		"原始数据组成"=>"$dir_abs/BMK_1_rawData/BMK_*_Data_Assess/PNG/Raw_data_ratio_*.png",
		"差异表达基因GO富集聚类图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/GO_cluster.png",
		"差异表达基因KEGG富集聚类图" => "$dir_abs/BMK_3_mRNA/BMK_*_DEG_Analysis/BMK_*_Anno_Cluster/KEGG_cluster.png",
		"CircRNA来源分布图"=>"$dir_abs/BMK_*_Personality/*CircRNA_Analysis/circRNA.type.pie.png",
		"CircRNA染色体分布图"=>"$dir_abs/BMK_*_Personality/*CircRNA_Analysis/chromosome_distrbution.png",
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
sub Title_write{
  my ($path)=@_;
  my $title=basename($path);
  if ($path=~/Cis_Anno_enrichment/){
	$title="LncRNA_Cis.".$title;
	#print "niulg $title\n";
  }elsif($path=~/Trans_Anno_enrichment/){
	$title="LncRNA_Trans.".$title;
	#print "niulg $title\n";
  }
  return $title;
}
#####2016-03-03
sub file_all_write{# 
	my ($format,$context) = @_;
#    my @context_file = &get_file($context);
#    my $context_file=@context_file;
#    print "aaa  $context\n";
#    print("@context_file\n");
    #push @context2,$context;
#    my $line= scalar @context_file;
#	my ($name)= split /\s+/,$context;
#	my @file_info = @{$file_all_info{$name}};
    if ($context=~/差异表达基因的KEGG富集网页结果/){
        my @context_file = &get_file($context);
        my $context_file=@context_file;
        my $line= scalar @context_file;
        my ($name)= split /\s+/,$context;
        my @file_info = @{$file_all_info{$name}};
        my $content ="\t".'<file_list name="'.$context.'" type="'.'xls'.'" desc="">'."\n";
        foreach my $path (@file_info) {
            my $title=basename$path;
            &table_html($path,"$dir/HTML/$title.html","$name","","","../src/",$context) if ($path!~/html$/ and $path!~/KEGG\.list/);
            &table_html($path,"$dir/HTML/$title.html","$name",",1,",\%{$second_links{KEGG}},"../src/",$context) if ($path=~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name","","","project_Template_Path/../src",$context,"project_Template_Path/../src/images/logo.jpg") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name",",1,",\%{$second_links{KEGG}},"project_Template_Path/../src",$context,"logo_image_path") if ($path=~/KEGG\.list/);

            $path="$dir/HTML/$title.html" if ($path!~/html$/);
            $content .="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
        }
        $content .="\t".'</file_list>'."\n";
        return $content;
    }
    elsif ($context=~/差异表达lncRNA顺式靶基因KEGG富集网页结果/){
        my @context_file = &get_file_lnc($context);
        my $context_file=@context_file;
        my $line= scalar @context_file;
        my ($name)= split /\s+/,$context;
        my @file_info = @{$file_all_info{$name}};
        my $content1 ="\t".'<file_list name="'.$context.'" type="'.'xls'.'" desc="">'."\n";
        foreach my $path (@file_info) {
            my $title=basename$path;
            $title="lnc_cis.$title";
            &table_html($path,"$dir/HTML/$title.html","$name","","","../src/",$context) if ($path!~/html$/ and $path!~/KEGG\.list/);
            &table_html($path,"$dir/HTML/$title.html","$name",",1,",\%{$second_links{lnc_cisKEGG}},"../src/",$context) if ($path=~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name","","","project_Template_Path/../src",$context,"project_Template_Path/../src/images/logo.jpg") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name",",1,",\%{$second_links{lnc_cisKEGG}},"project_Template_Path/../src",$context,"logo_image_path") if ($path=~/KEGG\.list/);

            $path="$dir/HTML/$title.html" if ($path!~/html$/);
            $content1 .="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
        }
        $content1 .="\t".'</file_list>'."\n";
        return $content1;
    }
    elsif ($context=~/差异表达lncRNA反式靶基因KEGG富集网页结果/){
        my @context_file = &get_file_lnc($context);
        my $context_file=@context_file;
        my $line= scalar @context_file;
        my ($name)= split /\s+/,$context;
        my @file_info = @{$file_all_info{$name}};
        my $content1 ="\t".'<file_list name="'.$context.'" type="'.'xls'.'" desc="">'."\n";
        foreach my $path (@file_info) {
            my $title=basename$path;
            $title="lnc_trans.$title";
            &table_html($path,"$dir/HTML/$title.html","$name","","","../src/",$context) if ($path!~/html$/ and $path!~/KEGG\.list/);
            &table_html($path,"$dir/HTML/$title.html","$name",",1,",\%{$second_links{lnc_transKEGG}},"../src/",$context) if ($path=~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name","","","project_Template_Path/../src",$context,"project_Template_Path/../src/images/logo.jpg") if ($path!~/html(\.cloud)?$/ and $path!~/KEGG\.list/);
            &table_html_cloud($path,"$dir/HTML/$title.html.cloud","$name",",1,",\%{$second_links{lnc_transKEGG}},"project_Template_Path/../src",$context,"logo_image_path") if ($path=~/KEGG\.list/);

            $path="$dir/HTML/$title.html" if ($path!~/html$/);
            $content1 .="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
        }
        $content1 .="\t".'</file_list>'."\n";
        return $content1;
    }

}

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
	my $content = "\t<$class name=$desc type=$type desc=$desc \/>\n";
	return $content;
}

sub table_html{
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	my @inputs=();
	$/="\n";
	open (IN,$input)||die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @tmp=split/\t+/;
		push@inputs,\@tmp;
	}
	$/="》";
    my $titles;
    if($input=~/Cis_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_cis.$titles";
    }
    elsif($input=~/Trans_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_trans.$titles";
    }else{
	    $titles=basename$input;
    }
    $title=~s/"//g;
	open HH,">$outHtml" or die "$!";
	print HH <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
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
	<body style=\"overflow:hidden\"> 
	<div class="container shadow"  id="box"><header><img src="$srcPath/images/logo.jpg" class="pull-right" />
	<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
	</header>
	</div>
HTML
        if($text){
		 	$text=~s/"//g;
                print HH "<div class=\"table-responsive\" id=\"textbox\"><p>$text</p></div>\n";
        }
        print HH "<div class=\"table-responsive\" id=\"box2\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info\">\n";
        if($text){
                print HH <<HTML;
<script type="text/javascript">
    var textbox=\$("#textbox").height();
    var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8-textbox;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
		 var height2 = height*0.8-textbox;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML
        }
        else{
        print HH <<HTML;
<script type="text/javascript">
        var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML

        }

       for (my $i=0;$i<=$#{$inputs[0]};$i++){
			print HH "<th>$inputs[0][$i]</th>\n";	
		}
        print HH "</tr></thead>\n<tbody>\n";
        for (my $k=1;$k<=$#inputs ;$k++) {
                print HH "<tr>";
                
				for (my $i=0;$i<=$#{$inputs[$k]};$i++){
                        if($linkColNum){
                                my $j=$i+1;
                                if($linkColNum=~/,$j,/){
                                        #print "out:$outHtml\n$k  $i $$input[$k][$i]\n"  if(!exists $$linkHash{$$input[$k][$i]});exit;
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
</div>
</body>
</html>
XGL
	close HH;

}

sub table_html_cloud{
    my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text,$logo_path)=@_;
    my @inputs=();
    $/="\n";
    open (IN,$input)||die $!;
    while(<IN>){
        chomp;
        next if /^\s*$/;
        my @tmp=split/\t+/;
        push@inputs,\@tmp;
    }
    $/="》";
    #my $titles = basename($input);
    my $titles;
    if($input=~/Cis_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_cis.$titles";
    }
    elsif($input=~/Trans_Anno_enrichment/){
        $titles=basename$input;
        $titles = "lnc_trans.$titles";
    }else{
	    $titles=basename$input;
    }
    open HH,">$outHtml" or die "$!";
	$title=~s/"//g;
    print HH <<HTML;
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
    <meta charset="UTF-8"></meta>
    <!--[if lt IE 9]>
    <script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
    <![endif]-->
    <meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
    <meta content="width=device-width, initial-scale=1" name="viewport"></meta>
    <!--added css-->
    <link rel="stylesheet" href="$srcPath/css/amend.css"/>
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
    <div class="container shadow"><header><img src="$logo_path" class="pull-right" />
    <div role="main" ><header><h2 id="title" class="text-center amend-title">$title</h2>
    </header>
    </div>
HTML
        if($text){
		 	$text=~s/"//g;
                print HH "<div class=\"table-responsive\" id=\"textbox\"><p>$text</p></div>\n";
        }
        print HH "<div class=\"table-responsive\" id=\"box2\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info amend\">\n";
		      if($text){
                print HH <<HTML;
<script type="text/javascript">
    var textbox=\$("#textbox").height();
    var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8-textbox;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
		 var height2 = height*0.8-textbox;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML
        }
        else{
        print HH <<HTML;
<script type="text/javascript">
        var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML

        }
        for (my $i=0;$i<=$#{$inputs[0]};$i++){
            print HH "<th>$inputs[0][$i]</th>\n";    
        }
        print HH "</tr></thead>\n<tbody>\n";
        for (my $k=1;$k<=$#inputs ;$k++) {
            print HH "<tr>";
            for (my $i=0;$i<=$#{$inputs[$k]};$i++){
                if($linkColNum){
                    my $j=$i+1;
                    if($linkColNum=~/,$j,/){
#                        print HH "<td><a href=\"$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
                        print HH "<td><a href=\"project_Template_Path/$$linkHash{$titles}{$inputs[$k][$i]}\" target=\"_blank\">$inputs[$k][$i]</a></td>";
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
</div>
</body>
</html>
XGL
    close HH;

}

sub fa_html{
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	open FA,"$input" ||die $!;
	my %inputs=();
	$/=">";
	while(<FA>){
		chomp;
		next if ($_=~/^\s*$/);
		my ($id,$seq)=split(/\n+/,$_,2);
		$inputs{$id}=$seq;
	}
	close FA;
	$/="》";
	open H,">$outHtml";
	print H <<HTML;
	<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
	<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
	<meta charset="UTF-8"></meta>
	<!--[if lt IE 9]>
	<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
	<![endif]-->
	<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
	<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
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
	<div style="word-wrap:break-word;word-break:break-all">
HTML
		if($text){
			print H "<p>$text</p>\n";
		}
		for my $key(keys %inputs){
			print H ">$key<br/>\n$inputs{$key}<br/>\n";	
		}
print H <<XGL;
	</div>
	</body>
	</html>
XGL
	close H;

}


sub Text{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "\t<p type=$type desc=$desc \/>\n";
	return $content;
}
sub Text_mul{
	my ($name,$type,$desc) = @_;
	$type = &T($type);
	my $content;
	foreach my $de (@$desc){
	  #print "niuling d####$de\n";	
		$de = &substitute($de);
		$de = &T($de);
		$de = &Trans($de);
		$content .= "\t<p type=$type desc=$de \/>\n";
	}
	
	
	return $content;
}


sub Anno{
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	$desc = &substitute($desc);
	$desc = &T($desc);
	$desc = &Trans($desc);
	my $content = "\t<p type=$type desc=$desc \/>\n";
	return $content;
}



sub abstract{
    my ($name,$type,$desc) = @_;
	$type = &T($type);
#	$desc = &substitute($desc);
#	$desc = &T($desc);
#	$desc = &Trans($desc);
    my $p_abstract = join("","&lt;p class=&quot; p-abstract&quot; &gt;",$desc,"&lt;/p&gt;");
#    $p_abstract = &substitute($p_abstract);
    $p_abstract = &Trans($p_abstract);
    $p_abstract = &T($p_abstract);
    my $content = "\t<report_abstract value=$p_abstract \/>\n";
#    "\t<p type=$type desc=$desc \/>\n";
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
	#my @lines = split /\n/,$desc;
	my @lines = @$desc;
	foreach my $line (@lines) {
		my @each = split /\t/,$line;
		$Ref_info{$each[0]} = $each[1];
	}
	my $content = "\t<ref_list name=\"参考文献\" type=$type desc=$desc_anno >\n";
	foreach my $id (sort {$a<=>$b} keys %Ref_info) {
		my $ref_name = $Ref_info{$id};
		my $ref_lin = $reference{$id};
		$id = &T($id);
		$ref_name = &T($ref_name);
		$ref_lin=&T($ref_lin);
		$content = $content."\t\t<ref id=$id name=$ref_name link=$ref_lin \/>\n";
	}
	$content = $content."\t</ref_list>\n";
	return $content;
	
}




sub Attach {
	my ($name,$type,$desc) =@_;
	$type = &T($type);
	my $content = "\t<Attach name=\"附录标题\" type=$type desc=$desc_anno \/>\n";
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

sub GetDate {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime( time() );
    return sprintf( "%4d\/%02d\/%02d", $year + 1900, $mon + 1, $day );
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
	Usage:
	Options:
	-indir		<dir>
	-config	<config.txt>
	-xml	<xml report>
	-h			help
USAGE
	print $usage;
	exit;
}
