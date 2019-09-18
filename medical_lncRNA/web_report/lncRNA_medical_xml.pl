#!/usr/bin/perl
use warnings;

use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);
use XML::Writer;
use IO::File;
use Encode;

my ($old,$output);
GetOptions(
        "h|?"           =>\&USAGE,
	"oldxml:s"	=>\$old,
	"o:s"		=>\$output,
	#"conf:s"	=>\$conf,
)or &USAGE;
&USAGE unless ($old);

#my %config=&readConfig($conf);
my $inpath=dirname $old;
$output ||= "$inpath/configtest_medical.xml";

my @images=();
my @tables=();
my ($h1,$h2,$h3,$h1id,$h2id,$h3id,$degid);

open(XML,$old)||die $!;
open(OUT,">$output")||die $!;
#my @tmp=();
#my @known_gff_species=("GRCh37","GRCh38","GRCm38","mm9","Rnor_6.0");
#open(TEST,">$inpath/test.out")||die $!;
#my $num=grep /^$config{db}$/,@know_gff_species;


$h1id=1;
$h2id=1;
$h3id=1;
$h1="";
$h2="";
$h3="";
my $known_flag=0;
while(my $write=<XML>){
	chomp($write);

	$h1=(split(/\"/,$write))[1]	if($write=~/\<h1/);
	$h2=(split(/\"/,$write))[1]	if($write=~/\<h2/);
	$h3=(split(/\"/,$write))[1]	if($write=~/\<h3/);

	
	if($h2 =~/lncRNA靶基因预测/){
		next if($known_flag ne 0);
		$known_flag=1;
		next if(!-e "$inpath/BMK_medical/BMK_1_Known_lncRNA/lncRNA_mRNA_id.list");
		print OUT &pattern("h3","type1","区分已知和新的lncRNA","区分已知和新的lncRNA"),"\n";
		print OUT &pattern("p","type1","整合当前数据库的已知lncRNA信息（人，GENCODE、NONCODE、lncipedia；小鼠，GENCODE、NONCODE；大鼠，NONCODE）来区分出的lncRNA中的已知lncRNA和新的lncRNA。预测转录本与已知数据库有相交的转录本，且两个转录本中存在至少一个相同的外显子或者预测的外显子在已知外显子里面。当预测转录本对应多个已知转录本，给每个已知转录本打分，预测转录本被认为是最高打分的已知转录本。初始分数为０，打分标准为：",""),"\n";
		print OUT &pattern("p","type1","1. 预测转录本在已知转录本中或完全相同，+1分",""),"\n";
		print OUT &pattern("p","type1","2. 当满足条件1且预测转录本与已知转录本有相同数目的外显子，+1分",""),"\n";
		print OUT &pattern("p","type1","3. 每当预测转录本中的一个外显子与已知转录本的外显子相交，+0.5分",""),"\n";
		print OUT &pattern("p","type1","4. 每当预测转录本中的一个外显子与已知转录本的一个外显子里面或相同，+0.5分",""),"\n";

		print OUT &pattern("file","xls","注：第一列为预测出的lncRNA转录本编号，第二列为数据库中已知的编号。--则说明不是已知的lncRNA。","已知的lncRNA列表","BMK_medical/BMK_1_Known_lncRNA/lncRNA_mRNA_id.list","xls"),"\n";
	}
	if($write =~/差异表达基因蛋白质互作网络图/){
		print OUT "$write";
		$h3="转录因子注释";
		print OUT &pattern("h3","type1","转录因子注释","转录因子注释"),"\n";
		print OUT &pattern("p","type1","转录水平调控是基因表达调控的重要环节，其中转录因子（Transcriptionfactor，TF）通过结合基因上游特异性核苷酸序列实现对该基因的转录调控，许多生物学功能都通过调节特异性转录因子实现对下游功能基因的激活或抑制，因此有必要注释不同分组间差异表达基因的转录组因子,使用动物转录因子数据库AnimalTFDB对差异表达基因中的转录因子进行鉴定。",""),"\n";
#		print OUT &pattern("p","type1","Transcription_factor分析结果",""),"\n";
		print OUT &pattern("file_list","xls","Transcription_factor分析结果","Transcription_factor分析结果"),"\n";
		@tables=glob("$inpath/BMK_medical/BMK_*_DEG_anno/*DEG.TF.xls");		
		foreach my $i(@tables){
			my $base=basename $i;
			my $relatepath=join("/",("BMK_medical",(split(/\//,$i))[-2],$base));
			print OUT &pattern("file","xls",$base,$base,$relatepath,"xls"),"\n";
		}	
		print OUT '</file_list>',"\n";
		print OUT &pattern("p","type1","注:ID，Ensembl基因ID；Gene_Symbol，基因通用名称；log2FC，样品间差异倍数变化；Regulated，差异上下调状态；Family，转录因子所属的家族。",""),"\n";


		@tables=glob("$inpath/BMK_medical/BMK_*_DEG_anno/*DEG.cosmic.xls");	
		my $cancer_num=@tables;
		if($cancer_num>0){
#			my $newid=join(".",($d1,$d2,($d3+2)));
			$h3="癌症基因功能注释";
			print OUT &pattern("h3","type1","癌症基因功能注释","癌症基因功能注释"),"\n";
			print OUT &pattern("p","type1","原癌基因（Proto-oncogene）是参与细胞生长、分裂和分化的正常基因，当其发生突变后（如基因序列被改变）就会变成致癌基因（Oncogene）。通常在肿瘤或恶性细胞系中某些特异性癌基因会上调表达，通过了解癌基因在不同实验组的表达情况有助于深入认识疾病的发病机理。Cosmic(https://cancer.sanger.ac.uk/cosmic)是英国Sanger实验室开发并维护的癌基因及相关注释数据库，有较高的权威性和可信度，通过Cosmic数据库，可对差异表达基因中的癌基因部分进行注释。",""),"\n";
#			print OUT &pattern("p","type1","Oncogenes分析结果",""),"\n";

			print OUT &pattern("file_list","xls","Oncogenes分析结果","Oncogenes分析结果"),"\n";
			foreach my $i(@tables){
				my $base=basename $i;
				my $relatepath=join("/",("BMK_medical",(split(/\//,$i))[-2],$base));
				print OUT &pattern("file","xls",$base,$base,$relatepath,"xls"),"\n";
			}	
			print OUT '</file_list>',"\n";
			print OUT &pattern("p","type1","注：ID，Ensembl基因ID；Gene_Symbol，基因通用名称；log2FC，样品间差异倍数变化；Regulated，差异上下调状态；Description，基因描述；Tumor Type(Somatic)，体细胞癌症类型；Tumor Type(Germline)，生殖细胞系癌症类型。",""),"\n";
		}
		next;
	}

	if($h2=~/基因表达量分析/){
#		my $currenth3=(split(/\s+/,$h3))[0];
#		@tmp=split(/\./,$currenth3);
#		
#		my $currenth2=(split(/\s+/,$h2))[0];
#		print "当前的h2标题是$h2\nid是$currenth2\n";
#		@tmp=split(/\./,$currenth2);
#		my $newid=join(".",($tmp[0],($tmp[1]+1)));
		$h2="融合基因分析";
		print OUT &pattern("h2","type1","融合基因分析","融合基因分析"),"\n";
		print OUT &pattern("p","type1","融合基因是指将两个或多个基因的编码区首尾相连，置于同一套调控序列(包括启动子、增强子、核糖体结合序列、终止子等)控制之下，构成的嵌合基因。融合基因的表达产物为融合蛋白。我们使用Fusionmap在转录组中研究基因融合事件。Fusionmap首先通过比对到基因组和转录本中双末端(pairend)关系的序列寻找候选的基因融合，然后采用通过与nt等数据库比较，过滤掉假阳性结果。",""),"\n";

		$pattern=&pattern("pic_list","type1","注：红色的线代表同一染色体上发生的融合事件，绿色的线代表不同染色体上发生的融合事件。","检测到的基因融合事件");
		print OUT "$pattern\n";
		@images=glob("$inpath/BMK_medical/BMK_*_Gene_fusion/png/*.png");
		foreach my $i(@images){
			my $base=basename $i;
			my $relatepath=join("/",("BMK_medical",(split(/\//,$i))[-3],"png",$base));
			print OUT &pattern("pic","type1",$base,$base,$relatepath),"\n";
		}	
		print OUT '</pic_list>',"\n";
		print OUT &pattern("h3","type1","融合基因统计","融合基因统计"),"\n";
		print OUT &pattern("file_list","xls","融合基因事件统计","融合基因事件统计",""),"\n";
		@tables=glob("$inpath/BMK_medical/BMK_*_Gene_fusion/report/*.xls");		
		my $fusiontitle=(split(/\//,$tables[0]))[-3];
		foreach my $i(@tables){
			my $base=basename $i;
			print OUT &pattern("file","xls",$base,$base,"BMK_medical/$fusiontitle/report/$base","xls"),"\n";
		}		
		print OUT '</file_list>',"\n";
		print OUT &pattern("p","type1","注：Fusion ID：基因融合事件编号；Unique Mapping Position：映射到唯一位点的基因融合数目；Count：基因融合位点总数；Gene：发生融合的基因；Strand：基因融合发生在正义链“+”还是反义链“-”；Chromosome：融合基因位于的染色体；Start：基因融合起始位点；End：基因融合终止位点；Filter：InFamilyList代表基因对属于同一家族；InParalogueList代表基因对来自同一旁系同源组；SameEnsembleGene代表基因对有一个super ensemble基因交集；InBlackList代表基因对包含线粒体和核糖体基因或者假基因。",""),"\n";

		my $pfam="$inpath/BMK_medical/$fusiontitle/Fusion_gene_Pfam.anno";
		if(-e $pfam){
			print OUT &pattern("h3","type1","融合基因结构域分析","融合基因结构域分析"),"\n";
			print OUT &pattern("p","type1","融合基因可能会导致具有新的或不同功能的基因产物生成，这些异常活性的蛋白质会发挥癌基因还可能与一个新的启动子发生融合，从而激活下游癌基因的表达。后者常见于淋巴肿瘤当中。我们对融合基因进行结构域及相关功能预测，有利于进一步探索其与肿瘤发生、发展或者转移过程中的作用机制。Pfam蛋白结构域搜索结果展示如下：",""),"\n";
#			print OUT &pattern("p","type1","融合基因的Pfam注释结果",""),"\n";
			print OUT &pattern("file","xls","注：Fusion Gene_ID，fusionmap产生的融合基因ID；Pfam acc，比对到pfam结构域的ID；Pfam name，pfam结构域名称；Pfam description，pfam结构域详细描述；Bit score，比对打分分值；E-value，比对evalue。","融合基因的Pfam注释全部结果","BMK_medical/$fusiontitle/Fusion_gene_Pfam.anno","xls"),"\n";
		}

	}
	print OUT "$write";
}


close(XML);
close(OUT);



########################################################################################################
#
#						Sub function
########################################################################################################



sub pattern{
	my($title,$type,$desc,$name,$path,$action)=@_;
	my $pattern="\<$title type=\"$type\" desc=\"$desc\" name=\"$name\"";
	$pattern .= " path =\"$path\" "		if($path);
	$pattern .= " action=\"$action\" "	if($action);

	if($title =~/list/){
		$pattern .="\>";
	}else{
		$pattern .="\/\>";
	}
	return $pattern;
}



sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}




sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	Options:
	-oldxml	<file>	input the xml file produced from the regular analysis, forced
	-o	<str>	output xml file, default the same path with old xml called configtest_medical.xml

	-h	Help

Example:

USAGE
	print $usage;
	exit;
}


