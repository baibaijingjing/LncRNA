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
my $desc_anno="";
$desc_anno=&T($desc_anno);
my $table_type="full";
$table_type=&T($table_type);

my $dir_abs = abs_path ($dir);
my $dir_template ="$Bin/Template";
#system "cp -r $template $dir_abs";

my $template="$dir/template";
mkdir $template unless (-d $template);

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

system "perl $Bin/data.tmp.pl -indir $dir_abs -od $dir_abs";



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

#my $com_num=scalar (@com);
#my $sep_num=scalar(@sep);
#my $total=$com_num+$sep_num;
#print "$com_num\n";
#print "$sep_num\n";

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
	$min =  &time_format($config{min});
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
    my ($Ref_seq1,$Ref_seq2) = (split/\//,$config{Ref_seq})[6,7];
    $ref_gene_addr = $Ref_seq1.$Ref_seq2;
}else{
	$ref_gene_addr=$config{ref_gene_addr};
}

my $ref_gene_name;
if ( !defined $config{ref_gene_name}) {
    my $Ref_seq = (split/\//,$config{Ref_seq})[6];
    $ref_gene_name = $Ref_seq;
}else{
	$ref_gene_name=$config{ref_gene_name};
}
&get_data_among_words();


my $ftemplate="$Bin/1.Basic_tem.txt";
my @ftemplates;
push @ftemplates,$ftemplate;
if (-d "$dir/mRNA/DEU"){
	$ftemplate="$Bin/2.DEU_tem.txt";
	push @ftemplates,$ftemplate;
}
if(-d "$dir/Personality/Conservation"){
	$ftemplate="$Bin/4.conservation.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/Known_LncRNA"){
	$ftemplate="$Bin/6.Known_lncRNA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/Tissue_specific"){
	$ftemplate="$Bin/5.tissue_specify.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/precursor"){
	$ftemplate="$Bin/7.precusor.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/miRNA_Target2LncRNA"){
	$ftemplate="$Bin/8.miRNA_target2lncRNA.txt";
	push @ftemplates,$ftemplate;
}
if (-d "$dir/Personality/PCA"){
	$ftemplate="$Bin/10.PCA.txt";
	push @ftemplates,$ftemplate;
}
#if (-d "$dir/Personality/WGCNA"){
#	$ftemplate="$Bin/11.WGCNA.txt";
#	push @ftemplates,$ftemplate;
if ($sample_num >= 4) {
	$ftemplate="$Bin/11.WGCNA.txt";
	push @ftemplates,$ftemplate;
}

if (-d "$dir/Personality/TF"){
	$ftemplate="$Bin/9.TF.txt";
	push @ftemplates,$ftemplate;
}


$ftemplate="$Bin/3.attach_tem.txt";
push @ftemplates,$ftemplate;

`cat @ftemplates >$dir/Template`;
$ftemplate="$dir/Template";
print "$ftemplate\n";




chomp (my $user = `whoami`);
my $user_addr = $user."\@biomarker.com.cn";
my $report_version = &T ("1.2");
my $report_name = &substitute($prefix);
$report_name=&T($prefix);
my $report_code = $config{Project_id};#&T("XXX");
#my $report_code = &T("XXX");
$report_code = &T($report_code);
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
		print "$text\n";
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
	$type = "type1";
	$type = &T($type);
	chomp $desc;
	#print "$desc\n";
	my $table_num = @tables;
	if (1<$table_num) {
		my $txt=$tables[0];
		$desc=&T($desc);
		my $txt_line=`less -S $txt|wc -l`; chomp $txt_line;
		if ($desc=~/DEU分析结果表/){
			my $txt_name="DEU分析结果示意表"; $txt_name=&T($txt_name);
			if($txt_line>30){
					`head -n 11 $txt >$dir/Template/DEU.txt`;
					my $txt_path="$dir/Template/DEU.txt"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			}else{
					my $path=&T($txt);
					$content="<table name=$desc type=$type path=$path \/>\n";
				}
			$content="\t<file_list name=$desc type=$type desc=$desc_anno >\n";
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
		elsif ($desc=~/SNP位点信息/){
			my $txt_name="SNP/InDel位点信息部分结果"; 
			$txt_name=&T($txt_name);
			my $txt_tmp=shift;
			if ($txt_line<30){
				$txt_tmp=$txt;
			}else {
				`head -n 11 $txt >$dir/template/snp.anno`;
				$txt_tmp="$dir/template/snp.anno";
			}
			$txt_tmp=&T($txt_tmp);
			$content = "<table name=$txt_name type=$type path=$txt_tmp \/>\n";
			my $anno="注：Chr：SNP/InDel位点所在染色体编号；Pos：SNP/InDel位点在染色体上的位置；Gene_id：SNP/InDel位点所在的基因或原来未注释的基因区（表中用Intergenic表示）；Ref：所选参考基因组中的SNP/InDel等位；Alt：测序样品中识别到的其他的SNP/InDel等位；L**：样品L**中SNP/InDel位点的分型；Depth：样品L**中SNP/InDel位点的测序深度；AlleDp：样品L**中SNP/InDel位点的各等位测序深度；Effect：SNP/InDel所在区域或类型；Codon_change：编码改变方式，未改变用点表示。核酸编码表见附表3，Effect具体说明详见：http://snpeff.sourceforge.net/SnpEff_manual.html。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/> \n";
			my $anno1="各样品的SNP/InDel位点详细信息见下表：";
			$anno1=&substitute($anno1); $anno1=&T($anno1);
			$content .="<p type=$type desc=$anno1 \/>\n";
			$content .="<file_list name=$desc type=$type >\n";
			foreach my $table (@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif ($desc=~/可变剪切事件结构和表达量统计表/){
			my $txt_name="可变剪切事件结构和表达量部分结果";
			$txt_name=&T($txt_name);
			my $txt_tmp=shift;
			if ($txt_line<30){
				$txt_tmp=$txt;
			}else {
				`head -n 11 $txt >$dir/template/AS.fpkm`;
				$txt_tmp="$dir/template/AS.fpkm";
			}
			$txt_tmp=&T($txt_tmp);
			$content = "<table name=$txt_name type=$type path=$txt_tmp \/>\n";
			my $anno="注：(1) event_id: AS事件编号；(2) event_type: AS事件类型； (3) gene_id: cufflink组装结果中的基因编号；(4) chrom: 染色体编号；(5) event_start: AS事件起始位置；(6) event_end: AS事件结束位置；(7) event_pattern: AS事件特征； (8) strand: 基因正负链信息；(9) fpkm: 此AS类型所在转录本的表达量。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/> \n";
			$content .="<file_list name=$desc type=$type >\n";
			foreach my $table (@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif ($desc=~/基因表达量结果文件/){
			my $txt_name="基因表达量部分结果"; $txt_name=&T($txt_name);
			my $anno="注： GeneID：基因ID号；Length：基因长度；FPKM：基因表达量；Locus：基因位置；Strand：正负链；Count：比对到该基因上的片段数；Normalized：校正后的片段数。";
			$anno=&substitute($anno); $anno=&T($anno);
			my $path_tmp=shift;
			if ($txt_line>30){
				`head -n 11 $txt >$dir/template/gene.fpkm`;
				$path_tmp="$dir/template/gene.fpkm";
			}
			else {
				$path_tmp=$txt;	
			}			
			$path_tmp=&T($path_tmp);
			$content="<table name=$txt_name type=$type path=$path_tmp \/>\n";
			$content .="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .="</file_list> \n";
		}
		elsif ($desc=~/差异表达分析结果/){
			my $txt_name="差异表达分析部分结果"; $txt_name=&T($txt_name);
			my $txt_tmp=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/deg.final`;
				$txt_tmp="$dir/template/deg.final";
			}else{
				$txt_tmp=$txt;
			}
			$txt_tmp=&T($txt_tmp);
			$content="<table name=$txt_name type=$type path=$txt_tmp \/>\n";
			my $anno="注：ID：基因编号；FDR：错误发现率；log2FC：表达量差异倍数的对数值；regulated：上调基因（up）还是下调基因（down）。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif($desc=~/差异表达基因topGO富集结果/){
			my $txt_name="差异表达基因topGO富集部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/topGO.result`;
				$txt_path="$dir/template/topGO.result";
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：GO.ID：GO term的ID；Term：GO功能；Annotated：所有基因注释到该功能的基因数；Significant：DEG注释到该功能的基因数；Expected：注释到该功能DEG数目的期望值；KS：富集Term的显著性统计，KS值越小，表明富集越显著。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .="</file_list>\n";
		}
		elsif($desc=~/差异表达基因的KEGG富集结果/){
			my $txt_name="差异表达基因的KEGG富集部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/KEGG.result`;
				$txt_path="$dir/template/KEGG.result";
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：#Kegg_pathway：KEGG通路名称；ko_id：KEGG通路ID；Cluter_frequency：通路中差异表达基因数占总的差异表达基因数的频率；Genome_frequency：通路中总的基因数占所有基因数的频率；P-value：富集显著性，P-value值越小，富集越显著；Corrected_P-value：修正后的富集显著性。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .="</file_list>\n";
		}
		elsif($desc=~/Cufflinks拼接结果/){
			my $txt_name="Cufflinks拼接部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>10){
				`head -n 6 $txt >$dir/template/Cufflinks.result`;
				$txt_path="$dir/template/Cufflinks.result";
			}
			else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：第一列：Chr：染色体序号；第二列：Source：来源；第三列：Type：转录本或者其外显子；第四列：Start Site：起始坐标；第五列： End Site：终止坐标；第六列：Score：对应位置得分；第七列：Strand：链的信息；第八列：frame：frame的信息；第九列：Description：描述信息。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type>\n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action />\n";
			}
			$content .="</file_list>\n";
		}
		elsif ($desc=~/Scripture拼接结果/){
			my $tit="Scripture拼接"; $tit=&T($tit);
			my $tit_class="h3";
			my $ann="Scripture 拼接基于统计学分段模型区分表达位点和实验背景噪音，比较适用于长转录本的拼接。Scripture 拼接结果如下表：";
			$ann=&substitute($ann); $ann=&T($ann);
			my $ann1="注：第一列：Chr：染色体序号；第二列：Source：来源；第三列：Type：转录本或者其外显子；第四列：Start Site：起始坐标；第五列： End Site：终止坐标；第六列：Score：对应位置得分；第七列：Strand：链的信息；第八列：frame：frame的信息；第九列：Description：描述信息。";
			$ann1=&substitute($ann1); $ann1=&T($ann1);
			$content ="<$tit_class name=$tit type=$type desc=$tit \/> \n";
			$content .="<p type=$type desc=$ann \/> \n";
			my $name1="Scripture拼接部分结果"; $name1=&T($name1);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 6 $txt >$dir/template/Scripture.result`;
				$txt_path="$dir/template/Scripture.result";
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content .="<table name=$name1 type=$type path=$txt_path \/> \n";
			$content .="<p type=$type desc=$ann1 \/> \n";
			$content .="<file_list name=$desc type=$type> \n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif ($desc=~/差异表达lncRNA结果/){
			my $txt_name="差异表达lncRNA部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc.deg`;
				$txt_path="$dir/template/lnc.deg"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：ID：lncRNA编号；FDR：错误发现率；log2FC：表达量差异倍数以2为底的对数值；regulated：上调lncRNA（up）还是下调lncRNA（down）。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content.="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type> \n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif($desc=~/差异表达lncRNA靶基因topGO富集结果/){
			my $txt_name="差异表达lncRNA靶基因topGO富集部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc_topGO.result`;
				$txt_path="$dir/template/lnc_topGO.result"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：GO.ID：GO term的ID；Term：GO功能；Annotated：所有注释到该功能的lncRNA靶基因数；Significant：DEG注释到该功能的lncRNA靶基因数；Expected：注释到该功能DEG数目的期望值；KS：富集Term的显著性统计，KS值越小，表明富集越显著。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content.="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type> \n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		elsif($desc=~/差异表达lncRNA靶基因KEGG富集结果/){
			my $txt_name="差异表达lncRNA靶基因KEGG富集部分结果"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc_KEGG.result`;
				$txt_path="$dir/template/lnc_KEGG.result"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="注：#Kegg_pathway：KEGG通路名称；ko_id：KEGG通路ID；Cluter_frequency：通路中差异表达基因数占总的差异表达基因数的频率；Genome_frequency：通路中总的基因数占所有基因数的频率；P-value：富集显著性，P-value值越小，富集越显著；Corrected_P-value：修正后的富集显著性。";
			$anno=&substitute($anno); $anno=&T($anno);
			$content.="<p type=$type desc=$anno \/>\n";
			$content .="<file_list name=$desc type=$type> \n";
			foreach my $table(@tables){
				my $table_name=&T(basename($table));
				my $path=&T($table);
				my $action=&T("xls");
				$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
			}
			$content .="</file_list> \n";
		}
		
		
		
		else{
			$content = "<file_list name=$desc type=$type >\n";
			foreach my $table (@tables) {
				my $table_name = &T(basename($table));
				my $path = &T($table);
				my $action = &T("xls");
				$content = $content."<file	name=$table_name	type=$type	path=$path	action=$action	/> \n";
			}
			$content = $content.'</file_list>'."\n";
		}
	}
	if ($table_num==1) {
		my $path = $tables[0];
		my $txt_line=`less -S $path|wc -l`; chomp $txt_line;
		my $table_name=&T($desc);
		if (-f "$path"){
			if ($desc=~/基因结构优化结果/){
				my $tmp_name="基因结构优化部分结果"; $tmp_name=&T($tmp_name);
				my $anno="注：GeneID：基因ID；Locus：基因座，格式为“染色体编号:起点坐标-终点坐标”；Strand：正负链；Site：优化的位置，3'或5'UTR；OriginalRegion：原来注释的第一个或最后一个外显子的起止坐标；OptimizedRegion：延伸之后的第一个或最后一个外显子的起止坐标。";
				$anno=&substitute($anno); $anno=&T($anno);
				if ($txt_line>30){
					`head -n 11 $path >$dir/template/gene.structure.optimize`;
					my $path_tmp="$dir/template/gene.structure.optimize";
					$path_tmp=&T($path_tmp);
					$path=&T($path);
					$content ="<table name=$tmp_name type=$type path=$path_tmp \/> \n";
					$content .="<p type=$type desc=$anno \/> \n";
					$content .="<file name=$table_name type=$type path=$path \/>\n";
				}else {
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif ($desc=~/新基因的GFF文件/){
				my $tmp_name="新基因的GFF文件部分结果"; $tmp_name=&T($tmp_name);
				my $anno="注：#Seq_ID：染色体号；Source：注释信息的来源，Cufflinks软件；Type：注释特征（Feature）类型；Start/End：特征序列的起止位置；Score：得分，数字，注释信息可能性的说明，“.”表示缺失值；Strand：特征序列所在的正负链；Phase：仅对注释类型为CDS有效，表示起始编码的位置，有效值为0、1、2，“.”表示缺失值；Attributes：以多个键值对组成的注释信息描述。";
				$anno=&substitute($anno); $anno=&T($anno);
				if ($txt_line>30){
					`head -n 11 $path >$dir/template/newgene.gff`;
					my $path_tmp="$dir/template/newgene.gff";
					$path_tmp=&T($path_tmp);
					$path=&T($path);
					$content ="<table name=$tmp_name type=$type path=$path_tmp \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$content .="<file name=$table_name type=$type path=$path \/>\n";
				}
				else {
					$path=&T($path);
					$content ="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/新基因序列FASTA文件/){
				my $tmp_name="新基因的FA文件部分结果"; $tmp_name=&T($tmp_name);
				my $anno="注：FASTA格式每一个序列单元以“>”开头，直到出现下一个“>”之前为止。“>”开头的行为序列ID行，后面紧接着基因ID；下面一行或多行为该基因的碱基序列。";
				$anno=&substitute($anno); $anno=&T($anno);
				if ($txt_line>30){
					`head -n 11 $path >$dir/template/newgene.fa`;
					my $path_tmp="$dir/template/newgene.fa";
					$path_tmp=&T($path_tmp);
					$path=&T($path);
					$content ="<table name=$tmp_name type=$type path=$path_tmp \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$content .="<file name=$table_name type=$type path=$path \/>\n";
				}
				else {
					$path=&T($path);
					$content ="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/CPC分析结果统计/){
				my $txt_name="CPC分析部分结果"; $txt_name=&T($txt_name);
				my $anno="注：#transcrits_id：转录本ID；length：ORF长度；feature：转录本类型；score:转录本得分，当score＜0时，为Noncoding。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/cpc.result`;
					my $txt_path="$dir/template/cpc.result";
					$txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content .="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/CNCI分析结果统计/){
				my $txt_name="CNCI分析部分结果"; $txt_name=&T($txt_name);
				my $anno="注：第1列：Transcript_ID：转录本ID；第2列：index：转录本类型；第3列：score:转录本得分，当score＜0时，为Noncoding；第4列：start：ORF起始位置；第5列：end：ORF终止位置；第6列：length：转录本长度。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/cnci.result`;
					my $txt_path="$dir/template/cnci.result"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content.="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content .="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content.="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/CPAT分析结果统计/){
				my $txt_name="CPAT分析部分结果"; $txt_name=&T($txt_name);
				my $anno="注：#ID：转录本ID；mRNA_size：转录本长度；ORF_size：ORF长度；Fickett_score：Fickett得分；Hexamer_score：Hexamer得分；coding_prob：编码能力得分；coding_noncoding：编码or非编码，no代表非编码。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/cpat.result`;
					my $txt_path="$dir/template/cpat.result"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content.="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/pfam分析结果统计/){
				my $txt_name="pfam分析部分结果"; $txt_name=&T($txt_name);
				my $anno="注：第1列：trans_id：转录本ID；第2列：hmm_acc：比对到pfam结构域ID；第3列： hmm_name：pfam结构域名称；第4列：hmm_start：比对到结构域的起始位置；第5列：hmm_end：比对到结构域的终止位置；第6列：hmm_length：pfam结构域的长度；第7列： bit_score比对打分值；第八列：E-value：比对的E值，pfam结构域筛选的条件E-value＜0.001。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/pfam.result`;
					my $txt_path="$dir/template/pfam.result"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/基于位置关系的LncRNA靶基因预测结果/){
				my $txt_name="基于位置关系的LncRNA靶基因预测部分结果"; $txt_name=&T($txt_name);
				my $anno="注：#lncRNA：lncRNA的id号；Genes：lncRNA对应的靶基因id号。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/lnc2target.position`;
					my $txt_path="$dir/template/lnc2target.position"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content.="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/基于互补序列的LncRNA靶基因预测结果/){
				my $txt_name="基于互补序列的LncRNA靶基因预测部分结果"; $txt_name=&T($txt_name);
				my $anno="注：#LncRNA_ID：lncRNA的id号；Target_mRNA_ID：lncRNA对应的靶基因ID。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/lnc2target.basepair`;
					my $txt_path="$dir/template/lnc2target.basepair"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			elsif($desc=~/lncRNA表达量结果/){
				my $txt_name="lncRNA表达量部分结果"; $txt_name=&T($txt_name);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/lnc.fpkm`;
					my $txt_path="$dir/template/lnc.fpkm"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
				}
			}
			elsif ($desc=~/DEU分析结果表/){
				my $txt_name="DEU分析结果示意表"; $txt_name=&T($txt_name);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/DEU.txt`;
					my $txt_path="$dir/template/DEU.txt"; $txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
				}
				
			}
			elsif($desc=~/DEU分析结果/){
				my $path=&T($path);
				my $action=&T("xls");
				$content .= "<file name=$table_name type=$type path=$path action=$action />\n";
			}
			elsif ($desc=~/差异表达分析结果/){
				my $txt_name="差异表达分析部分结果"; $txt_name=&T($txt_name);
				my $txt_path=shift;
				my $anno="注：ID：基因编号；FDR：错误发现率；log2FC：表达量差异倍数的对数值；regulated：上调基因（up）还是下调基因（down）。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/deg.final`;
					$txt_path="$dir/template/deg.final"; 
				
					$txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					$content.="<file name=$table_name type=$type path=$path \/>\n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
				
			}
			elsif ($desc=~/差异表达lncRNA结果/){
				my $txt_name="差异表达lncRNA部分结果"; $txt_name=&T($txt_name);
				my $txt_path=shift;
				my $anno="注：ID：lncRNA编号；FDR：错误发现率；log2FC：表达量差异倍数以2为底的对数值；regulated：上调lncRNA（up）还是下调lncRNA（down）。";
				$anno=&substitute($anno); $anno=&T($anno);
				if($txt_line>30){
					`head -n 11 $path >$dir/template/lnc.deg`;
					$txt_path="$dir/template/lnc.deg"; 
					
					$txt_path=&T($txt_path);
					$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
					$content.="<p type=$type desc=$anno \/>\n";
					$path=&T($path);
					my $action=&T("xls");
					$content .="<file name=$table_name type=$type path=$path action=$action /> \n";
				}else{
					$path=&T($path);
					$content="<table name=$table_name type=$type path=$path \/>\n";
					$content .="<p type=$type desc=$anno \/>\n";
				}
			}
			else {
				$path=&T($path);
				$content="<table name=$table_name type=$type path=$path \/>\n";
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
		"Cufflinks拼接结果" => "$dir/LncRNA/Assembly/*.Cufflinks.transcripts.gtf",
		"Scripture拼接结果" => "$dir/LncRNA/Assembly/*.Scripture.transcripts.gtf",
		"CPC分析结果统计" => "$dir/LncRNA/Identify/CPC.txt",
		"CNCI分析结果统计" => "$dir/LncRNA/Identify/CNCI.txt",
		"CPAT分析结果统计" => "$dir/LncRNA/Identify/cpat.txt",
		"pfam分析结果统计" => "$dir/LncRNA/Identify/Pfam.txt",
		"基于位置关系的LncRNA靶基因预测结果" => "$dir/LncRNA/Target/lncRNA_position.target",
		"基于互补序列的LncRNA靶基因预测结果"=>"$dir/LncRNA/Target/lncRNA_basepair.target",
		"lncRNA表达量结果" => "$dir/LncRNA/DEG/LncRNA_fpkm.list",
		"差异表达lncRNA结果" => "$dir/LncRNA/DEG/*_vs_*/*_vs_*.DEG_final.xls",
		"差异表达lncRNA数目统计表" => "$dir/LncRNA/DEG/Diff_Lnc.stat",
		"注释的差异表达lncRNA靶基因数量统计表" => "$dir/LncRNA/DEG/Diff_lnc.anno.stat",
		"差异表达lncRNA靶基因topGO富集结果" => "$dir/LncRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.xls",
		"差异表达lncRNA靶基因KEGG富集结果" => "$dir/LncRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.stat",
		"DEU分析结果表" => "$dir/mRNA/DEU/*_vs_*/DEU_Result_Final.xls",
		"DEU分析结果"=>"$dir/mRNA/DEU/*_vs_*/DEXSeqReport/testForDEU.html",
		#"差异表达lncRNA靶基因的KEGG通路注释"=>"$dir/LncRNA/DEG/*_vs_*/pathway/kegg_map/*.html",
		"附表1 软件列表" => "$dir_template/software_used_list.txt",
		"附表2 数据库列表" => "$dir_template/database_used_list.txt",
		"附表3 核酸编码表" => "$dir_template/Nucleic_acids_encoding_table.txt",
		"LncRNA的gff文件"=>"$dir/LncRNA/Identify/LncRNA.gff",
		"已知lncRNA鉴定统计表"=>"$dir/Personality/Known_LncRNA/known_lncRNA.result.txt",
		"miRNA前体的lncRNA统计"=>"$dir/Personality/precursor/lncRNA_precursor.txt",
		"miRNA靶向lncRNA分析表"=>"$dir/Personality/miRNA_Target2LncRNA/*.mir2target.list",
		"差异基因转录因子分析表"=>"$dir/Personality/TF/TF.txt",
		#"SNP/InDel位点信息示意表"=>"$dir_template/snp_indel.anno.gatk.list",
		#"可变剪切事件结构和表达量示意表"=>"$dir_template/ASprofile.fpkm",
		#"基因结构优化部分结果"=>"$dir_template/gene.structure.optimize.txt",
		#"新基因的GFF文件示意表"=>"$dir_template/newgene.gff",
		#"新基因序列FASTA文件示意表"=>"$dir_template/newgene.fa",
		#"基因表达量结果示意表"=>"$dir_template/geneExpression.txt",
		"miRNA前体的lncRNA统计示意表"=>"$dir_template/precursor.txt",
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
		if ($desc=~/原始数据组成/){
			my $tit="测序质量控制"; $tit=&T($tit);
			my $tit_class="h3"; 
			my $ann="在进行数据分析之前，首先需要确保这些Reads有足够高的质量，以保证后续分析的准确。百迈客对数据进行严格的质量控制，进行如下过滤方式：";
			$ann .="\n"."在进行数据分析之前，首先需要确保这些Reads有足够高的质量，以保证后续分析的准确。百迈客对数据进行严格的质量控制，进行如下过滤方式：";
			$ann .="\n"."(1) 去除含有接头的Reads；";
			$ann .="\n"."(2) 去除低质量的Reads（包括去除N的比例大于10%的Reads；去除质量值Q≤10的碱基数占整条Read的50%以上的Reads）。";
			$ann .="\n"."经过上述一系列的质量控制之后得到的高质量的Clean Data，以FASTQ格式提供。";
			$ann=&substitute($ann); $ann=&T($ann);
			my $ann1="注：Adapter related：过滤掉的含有接头Reads数占总Raw Reads数的比例。Low quality：过滤掉的低质量Reads数占总Raw Reads数的比例。Clean Reads：经过以上过滤得到的Clean Reads 数占总Raw Reads 数的比例。";
			$ann1=&substitute($ann1); $ann1=&T($ann1);
			$content = "<$tit_class name=$tit type=$type desc=$tit \/> \n";
			$content .="<p type=$type desc=$ann \/> \n";
			$desc=&T($desc);
			$content .="<pic_list name=$desc type=$type> \n";
			foreach my $pict(@picts){
				my $pict_name=&T(basename($pict));
				my $path=&T($pict);
				$content .="<pic name=$pict_name path=$path \/> \n";
			}
			$content .="</pic_list> \n";
		}
		else{
			$desc = &T($desc);
			$content = "<pic_list	name=$desc	type=$type > \n";
			foreach my $pict(@picts) {
				my $pict_name = &T (basename($pict));
				my $path = &T($pict);
				$content = $content."<pic	name=$pict_name	path=$path \/> \n";
			}
		$content = $content.'</pic_list>'."\n";
		}
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
			elsif ($desc=~/差异表达基因维恩图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="将每组差异基因进行维恩图绘制，见下图。图中展示了各比较组特有的差异基因的个数，以及比较组间的共有的差异基因个数。"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $tit="差异表达基因维恩图";	$tit=&T($tit);
				my $tit_class="h3";	
				$content="<$tit_class name=$tit type=$type desc=$tit \/>\n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/>\n";
			}
			elsif ($desc=~/差异表达基因聚类折线图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="注：x轴表示实验条件，y轴表示标准化的FPKM。黑色线条表示这个cluster中的所有基因在不同实验条件下相对表达量的平均值的折线图。"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="对筛选出的所有差异表达基因的表达量做K-means聚类分析，对差异基因的表达量水平值log2(FPKM+1) 进行聚类分析,得到差异基因在不同实验条件下的表达模式，表达模式相同或相近的基因聚集成类。由于同类的基因可能具有相似的功能，或是共同参与同一代谢过程或细胞通路。因此，通过对同一个cluster的基因进行功能注释及富集分析，可预测未知基因的功能或已知基因的未知功能。";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content ="<p type=$type desc=$ann1 \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				
			}
			elsif ($desc=~/差异表达lncRNA聚类图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="对筛选出的差异表达lncRNA做层次聚类分析，将具有相同或相似表达行为的lncRNA进行聚类，差异表达lncRNA聚类结果如下图：" ; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1 ="注：图中不同的列代表不同的样品，不同的行代表不同的lncRNA。颜色代表了lncRNA在样品中的表达量水平log2（FPKM+1）。"; 
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				my $tit="差异表达lncRNA聚类分析"; $tit=&T($tit);
				my $tit_class='h3'; 
				$content ="<$tit_class name=$tit type=$type desc=$tit \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/> \n";
				$content .= "<p type=$type desc=$ann1 \/> \n";
			}
			elsif ($desc=~/差异表达lncRNA靶基因聚类折线图/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="注：x轴表示实验条件，y轴表示标准化的FPKM。黑色线条表示这个cluster中的所有基因在不同实验条件下相对表达量的平均值的折线图。";
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="对筛选出的所有lncRNA的表达量做层次聚类分析，得到lncRNA在不同实验条件下的表达模式，表达模式相同或相近的lncRNA聚集成类。";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content = "<pic name=$picture_name type=$type path=$path \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<p type=$type desc=$ann1 \/> \n";
			}
			elsif ($desc=~/样品间相关性图/){
				if (exists $config{'Sep'}){
					my $picture_name=&T($desc);
					my $tit="重复相关性评估"; $tit=&T($tit);
					my $tit_class="h3";
					my $ann="研究表明，基因的表达在不同的个体间存在生物学可变性（Biological Variability），不同的基因之间表达的可变程度存在差异，而长链非编码测序技术、qPCR以及生物芯片等技术都不能消除这种可变性。为了寻找真正感兴趣的差异表达基因，需要考虑和处理因生物学可变性造成的表达差异。目前最常用且最有效的方法是在实验设计中设立生物学重复（Biological Replicates）。重复条件限制越严格，重复样品数目越多，寻找到的差异表达基因越可靠。对于设立生物学重复的项目，评估生物学重复的相关性对于分析转录组测序数据非常重要。生物学重复的相关性不仅可以检验生物学实验操作的可重复性；还可以评估差异表达基因的可靠性和辅助异常样品的筛查。";
					$ann .="\n"."将皮尔逊相关系数r（Pearson’s Correlation Coefficient）作为生物学重复相关性的评估指标。r2越接近1，说明两个重复样品相关性越强。";
					my $ann1="注：横坐标表示样品名称，纵坐标表示对应的样品名称。颜色代表r2值大小。";
					$ann1=&substitute($ann1); $ann1=&T($ann1);
					$ann=&substitute($ann); $ann=&T($ann);
					$path=&T($path);
					$content="<$tit_class name=$tit type=$type desc=$tit \/> \n";
					$content .="<p type=$type desc=$ann \/> \n";
					$content .="<pic name=$picture_name type=$type path=$path \/> \n";
					$content .="<p type=$type desc=$ann1 \/> \n";
				}
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
		"碱基错误率分布图" => "$dir_abs/QC/Clean_Data/PNG/*.quality.png",
		"ATGC含量分布图" => "$dir_abs/QC/Clean_Data/PNG/*.acgtn.png",
		"公式1 质量值计算公式" => "$dir_template/F01_Qscore_formula.png",
		"公式2 FPKM计算公式" => "$dir_template/F02_FPKM_formula.png",
		"TopHat2分析流程" => "$dir_template/TopHat2_workflow.png",
		"Mapped Reads在参考基因组上的位置及覆盖深度分布图" => "$dir_abs/QC/Map_assess/*.map.png",
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
		"差异表达基因聚类图" => "$dir_abs/mRNA/DEG/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster.png",
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
		"差异表达lncRNA聚类图" => "$dir_abs/LncRNA/DEG/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster.png",
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
		"Jensen_Shannon divergence公式"=>"$dir_template/JS.png",
		"lncRNA特异性表达图"=>"$dir_abs/Personality/Tissue_specific/*.png",
		"差异基因主成分分析图"=>"$dir_abs/Personality/PCA/*.png",
		"mRNA与lncRNA共表达网络分析示意图"=>"$dir_abs/LncRNA/Target/Trans_target/Fig_8_1_networkHeatmap.png",
		"原始数据组成"=>"$dir_abs/QC/Clean_Data/PNG/Raw_data_ratio_*.png",
		
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
