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
$/="��";
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
	if ($format eq "���") {### ���
		print OUT &table_write($format,$type,$text);
		print "$text\n";
	}
	elsif ($format =~ /ͼƬ/) {### ͼƬ
		print OUT &picture_write($format,$type,$text);
	}
	#elsif ($format =~ /����/) { ## ����
	#	print OUT &Link($format,$type,$text);
	#}
	elsif ($format =~ /��ʽ/) {### ��ʽ
		print OUT &picture_write($format,$type,$text);
	}
	elsif ($format =~ /������/) {### ����
		print OUT &Head($format,$type,$text);
	}
	elsif ($format =~ /����/) {### ����
		print OUT &Text($format,$type,$text);
	}
	elsif ($format eq 'ע��') {
		print OUT &Anno($format,$type,$text);
	}
	elsif ($format eq '�ο�����') {
		print OUT &Reference ($format,$type,$text);
	}
	elsif ($format eq '��¼����') {
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
		if ($desc=~/DEU���������/){
			my $txt_name="DEU�������ʾ���"; $txt_name=&T($txt_name);
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
		elsif($desc=~/DEU�������/){
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
		elsif ($desc=~/SNPλ����Ϣ/){
			my $txt_name="SNP/InDelλ����Ϣ���ֽ��"; 
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
			my $anno="ע��Chr��SNP/InDelλ������Ⱦɫ���ţ�Pos��SNP/InDelλ����Ⱦɫ���ϵ�λ�ã�Gene_id��SNP/InDelλ�����ڵĻ����ԭ��δע�͵Ļ�������������Intergenic��ʾ����Ref����ѡ�ο��������е�SNP/InDel��λ��Alt��������Ʒ��ʶ�𵽵�������SNP/InDel��λ��L**����ƷL**��SNP/InDelλ��ķ��ͣ�Depth����ƷL**��SNP/InDelλ��Ĳ�����ȣ�AlleDp����ƷL**��SNP/InDelλ��ĸ���λ������ȣ�Effect��SNP/InDel������������ͣ�Codon_change������ı䷽ʽ��δ�ı��õ��ʾ���������������3��Effect����˵�������http://snpeff.sourceforge.net/SnpEff_manual.html��";
			$anno=&substitute($anno); $anno=&T($anno);
			$content .="<p type=$type desc=$anno \/> \n";
			my $anno1="����Ʒ��SNP/InDelλ����ϸ��Ϣ���±�";
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
		elsif ($desc=~/�ɱ�����¼��ṹ�ͱ����ͳ�Ʊ�/){
			my $txt_name="�ɱ�����¼��ṹ�ͱ�������ֽ��";
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
			my $anno="ע��(1) event_id: AS�¼���ţ�(2) event_type: AS�¼����ͣ� (3) gene_id: cufflink��װ����еĻ����ţ�(4) chrom: Ⱦɫ���ţ�(5) event_start: AS�¼���ʼλ�ã�(6) event_end: AS�¼�����λ�ã�(7) event_pattern: AS�¼������� (8) strand: ������������Ϣ��(9) fpkm: ��AS��������ת¼���ı������";
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
		elsif ($desc=~/������������ļ�/){
			my $txt_name="�����������ֽ��"; $txt_name=&T($txt_name);
			my $anno="ע�� GeneID������ID�ţ�Length�����򳤶ȣ�FPKM������������Locus������λ�ã�Strand����������Count���ȶԵ��û����ϵ�Ƭ������Normalized��У�����Ƭ������";
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
		elsif ($desc=~/������������/){
			my $txt_name="������������ֽ��"; $txt_name=&T($txt_name);
			my $txt_tmp=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/deg.final`;
				$txt_tmp="$dir/template/deg.final";
			}else{
				$txt_tmp=$txt;
			}
			$txt_tmp=&T($txt_tmp);
			$content="<table name=$txt_name type=$type path=$txt_tmp \/>\n";
			my $anno="ע��ID�������ţ�FDR���������ʣ�log2FC����������챶���Ķ���ֵ��regulated���ϵ�����up�������µ�����down����";
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
		elsif($desc=~/���������topGO�������/){
			my $txt_name="���������topGO�������ֽ��"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/topGO.result`;
				$txt_path="$dir/template/topGO.result";
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="ע��GO.ID��GO term��ID��Term��GO���ܣ�Annotated�����л���ע�͵��ù��ܵĻ�������Significant��DEGע�͵��ù��ܵĻ�������Expected��ע�͵��ù���DEG��Ŀ������ֵ��KS������Term��������ͳ�ƣ�KSֵԽС����������Խ������";
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
		elsif($desc=~/����������KEGG�������/){
			my $txt_name="����������KEGG�������ֽ��"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/KEGG.result`;
				$txt_path="$dir/template/KEGG.result";
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="ע��#Kegg_pathway��KEGGͨ·���ƣ�ko_id��KEGGͨ·ID��Cluter_frequency��ͨ·�в����������ռ�ܵĲ������������Ƶ�ʣ�Genome_frequency��ͨ·���ܵĻ�����ռ���л�������Ƶ�ʣ�P-value�����������ԣ�P-valueֵԽС������Խ������Corrected_P-value��������ĸ��������ԡ�";
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
		elsif($desc=~/Cufflinksƴ�ӽ��/){
			my $txt_name="Cufflinksƴ�Ӳ��ֽ��"; $txt_name=&T($txt_name);
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
			my $anno="ע����һ�У�Chr��Ⱦɫ����ţ��ڶ��У�Source����Դ�������У�Type��ת¼�������������ӣ������У�Start Site����ʼ���ꣻ�����У� End Site����ֹ���ꣻ�����У�Score����Ӧλ�õ÷֣������У�Strand��������Ϣ���ڰ��У�frame��frame����Ϣ���ھ��У�Description��������Ϣ��";
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
		elsif ($desc=~/Scriptureƴ�ӽ��/){
			my $tit="Scriptureƴ��"; $tit=&T($tit);
			my $tit_class="h3";
			my $ann="Scripture ƴ�ӻ���ͳ��ѧ�ֶ�ģ�����ֱ��λ���ʵ�鱳���������Ƚ������ڳ�ת¼����ƴ�ӡ�Scripture ƴ�ӽ�����±�";
			$ann=&substitute($ann); $ann=&T($ann);
			my $ann1="ע����һ�У�Chr��Ⱦɫ����ţ��ڶ��У�Source����Դ�������У�Type��ת¼�������������ӣ������У�Start Site����ʼ���ꣻ�����У� End Site����ֹ���ꣻ�����У�Score����Ӧλ�õ÷֣������У�Strand��������Ϣ���ڰ��У�frame��frame����Ϣ���ھ��У�Description��������Ϣ��";
			$ann1=&substitute($ann1); $ann1=&T($ann1);
			$content ="<$tit_class name=$tit type=$type desc=$tit \/> \n";
			$content .="<p type=$type desc=$ann \/> \n";
			my $name1="Scriptureƴ�Ӳ��ֽ��"; $name1=&T($name1);
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
		elsif ($desc=~/������lncRNA���/){
			my $txt_name="������lncRNA���ֽ��"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc.deg`;
				$txt_path="$dir/template/lnc.deg"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="ע��ID��lncRNA��ţ�FDR���������ʣ�log2FC����������챶����2Ϊ�׵Ķ���ֵ��regulated���ϵ�lncRNA��up�������µ�lncRNA��down����";
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
		elsif($desc=~/������lncRNA�л���topGO�������/){
			my $txt_name="������lncRNA�л���topGO�������ֽ��"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc_topGO.result`;
				$txt_path="$dir/template/lnc_topGO.result"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="ע��GO.ID��GO term��ID��Term��GO���ܣ�Annotated������ע�͵��ù��ܵ�lncRNA�л�������Significant��DEGע�͵��ù��ܵ�lncRNA�л�������Expected��ע�͵��ù���DEG��Ŀ������ֵ��KS������Term��������ͳ�ƣ�KSֵԽС����������Խ������";
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
		elsif($desc=~/������lncRNA�л���KEGG�������/){
			my $txt_name="������lncRNA�л���KEGG�������ֽ��"; $txt_name=&T($txt_name);
			my $txt_path=shift;
			if($txt_line>30){
				`head -n 11 $txt >$dir/template/lnc_KEGG.result`;
				$txt_path="$dir/template/lnc_KEGG.result"; 
			}else{
				$txt_path=$txt;
			}
			$txt_path=&T($txt_path);
			$content="<table name=$txt_name type=$type path=$txt_path \/>\n";
			my $anno="ע��#Kegg_pathway��KEGGͨ·���ƣ�ko_id��KEGGͨ·ID��Cluter_frequency��ͨ·�в����������ռ�ܵĲ������������Ƶ�ʣ�Genome_frequency��ͨ·���ܵĻ�����ռ���л�������Ƶ�ʣ�P-value�����������ԣ�P-valueֵԽС������Խ������Corrected_P-value��������ĸ��������ԡ�";
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
			if ($desc=~/����ṹ�Ż����/){
				my $tmp_name="����ṹ�Ż����ֽ��"; $tmp_name=&T($tmp_name);
				my $anno="ע��GeneID������ID��Locus������������ʽΪ��Ⱦɫ����:�������-�յ����ꡱ��Strand����������Site���Ż���λ�ã�3'��5'UTR��OriginalRegion��ԭ��ע�͵ĵ�һ�������һ�������ӵ���ֹ���ꣻOptimizedRegion������֮��ĵ�һ�������һ�������ӵ���ֹ���ꡣ";
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
			elsif ($desc=~/�»����GFF�ļ�/){
				my $tmp_name="�»����GFF�ļ����ֽ��"; $tmp_name=&T($tmp_name);
				my $anno="ע��#Seq_ID��Ⱦɫ��ţ�Source��ע����Ϣ����Դ��Cufflinks�����Type��ע��������Feature�����ͣ�Start/End���������е���ֹλ�ã�Score���÷֣����֣�ע����Ϣ�����Ե�˵������.����ʾȱʧֵ��Strand�������������ڵ���������Phase������ע������ΪCDS��Ч����ʾ��ʼ�����λ�ã���ЧֵΪ0��1��2����.����ʾȱʧֵ��Attributes���Զ����ֵ����ɵ�ע����Ϣ������";
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
			elsif($desc=~/�»�������FASTA�ļ�/){
				my $tmp_name="�»����FA�ļ����ֽ��"; $tmp_name=&T($tmp_name);
				my $anno="ע��FASTA��ʽÿһ�����е�Ԫ�ԡ�>����ͷ��ֱ��������һ����>��֮ǰΪֹ����>����ͷ����Ϊ����ID�У���������Ż���ID������һ�л����Ϊ�û���ļ�����С�";
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
			elsif($desc=~/CPC�������ͳ��/){
				my $txt_name="CPC�������ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע��#transcrits_id��ת¼��ID��length��ORF���ȣ�feature��ת¼�����ͣ�score:ת¼���÷֣���score��0ʱ��ΪNoncoding��";
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
			elsif($desc=~/CNCI�������ͳ��/){
				my $txt_name="CNCI�������ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע����1�У�Transcript_ID��ת¼��ID����2�У�index��ת¼�����ͣ���3�У�score:ת¼���÷֣���score��0ʱ��ΪNoncoding����4�У�start��ORF��ʼλ�ã���5�У�end��ORF��ֹλ�ã���6�У�length��ת¼�����ȡ�";
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
			elsif($desc=~/CPAT�������ͳ��/){
				my $txt_name="CPAT�������ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע��#ID��ת¼��ID��mRNA_size��ת¼�����ȣ�ORF_size��ORF���ȣ�Fickett_score��Fickett�÷֣�Hexamer_score��Hexamer�÷֣�coding_prob�����������÷֣�coding_noncoding������or�Ǳ��룬no����Ǳ��롣";
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
			elsif($desc=~/pfam�������ͳ��/){
				my $txt_name="pfam�������ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע����1�У�trans_id��ת¼��ID����2�У�hmm_acc���ȶԵ�pfam�ṹ��ID����3�У� hmm_name��pfam�ṹ�����ƣ���4�У�hmm_start���ȶԵ��ṹ�����ʼλ�ã���5�У�hmm_end���ȶԵ��ṹ�����ֹλ�ã���6�У�hmm_length��pfam�ṹ��ĳ��ȣ���7�У� bit_score�ȶԴ��ֵ���ڰ��У�E-value���ȶԵ�Eֵ��pfam�ṹ��ɸѡ������E-value��0.001��";
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
			elsif($desc=~/����λ�ù�ϵ��LncRNA�л���Ԥ����/){
				my $txt_name="����λ�ù�ϵ��LncRNA�л���Ԥ�ⲿ�ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע��#lncRNA��lncRNA��id�ţ�Genes��lncRNA��Ӧ�İл���id�š�";
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
			elsif($desc=~/���ڻ������е�LncRNA�л���Ԥ����/){
				my $txt_name="���ڻ������е�LncRNA�л���Ԥ�ⲿ�ֽ��"; $txt_name=&T($txt_name);
				my $anno="ע��#LncRNA_ID��lncRNA��id�ţ�Target_mRNA_ID��lncRNA��Ӧ�İл���ID��";
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
			elsif($desc=~/lncRNA��������/){
				my $txt_name="lncRNA��������ֽ��"; $txt_name=&T($txt_name);
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
			elsif ($desc=~/DEU���������/){
				my $txt_name="DEU�������ʾ���"; $txt_name=&T($txt_name);
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
			elsif($desc=~/DEU�������/){
				my $path=&T($path);
				my $action=&T("xls");
				$content .= "<file name=$table_name type=$type path=$path action=$action />\n";
			}
			elsif ($desc=~/������������/){
				my $txt_name="������������ֽ��"; $txt_name=&T($txt_name);
				my $txt_path=shift;
				my $anno="ע��ID�������ţ�FDR���������ʣ�log2FC����������챶���Ķ���ֵ��regulated���ϵ�����up�������µ�����down����";
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
			elsif ($desc=~/������lncRNA���/){
				my $txt_name="������lncRNA���ֽ��"; $txt_name=&T($txt_name);
				my $txt_path=shift;
				my $anno="ע��ID��lncRNA��ţ�FDR���������ʣ�log2FC����������챶����2Ϊ�׵Ķ���ֵ��regulated���ϵ�lncRNA��up�������µ�lncRNA��down����";
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
		"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��" => "$dir_template/error.txt",
		"��Ʒ������������ͳ�Ʊ�" => "$dir/QC/Clean_Data/AllSample_GC_Q.stat",
		"��Ʒ������������ѡ�ο�����������бȶԽ��ͳ�Ʊ�" => "$dir/QC/Map_assess/Total_mapped.stat",
		"SNPλ����Ϣ" => "$dir/mRNA/Structure/SNP/snp_anno/final.snp.anno.gatk.*.list",
		"InDelλ����Ϣ" => "$dir/mRNA/Structure/SNP/indel_anno/final.indel.anno.gatk.*.list",
		"SNPλ��ͳ�Ʊ�" => "$dir/mRNA/Structure/SNP/AllSample.snp.stat",
		"�ɱ�����¼��ṹ�ͱ����ͳ�Ʊ�" => "$dir/mRNA/Structure/Alt_splice/*.fpkm",
		"����ṹ�Ż����" => "$dir/mRNA/Gene_Structure_Optimize/*.geneStructure.optimize.xls",
		"�»����GFF�ļ�" => "$dir/mRNA/NewGene/*.newGene_final.filtered.gff",
		"�»�������FASTA�ļ�" => "$dir/mRNA/NewGene/*.newGene.longest_transcript.fa",
		"�»�����ע�ͽ��ͳ��" => "$dir/mRNA/NewGene/Annotation/Function_Annotation.stat.xls",
		"������������ļ�" => "$dir/mRNA/Expression/*.geneExpression.xls",
		"������������" => "$dir/mRNA/DEG/*_vs_*/*_vs_*.DEG_final.xls",
		"�����������Ŀͳ�Ʊ�" => "$dir/mRNA/DEG/DEG.stat",
		"ע�͵Ĳ������������ͳ�Ʊ�" => "$dir/mRNA/DEG/DEG.anno.stat",
		"���������topGO�������" => "$dir/mRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.xls",
		"����������KEGG�������" => "$dir/mRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.stat",
		#"����������KEGGͨ·ע��" => "$dir/mRNA/DEG/*_vs_*/pathway/kegg_map/*.html",
		"Cufflinksƴ�ӽ��" => "$dir/LncRNA/Assembly/*.Cufflinks.transcripts.gtf",
		"Scriptureƴ�ӽ��" => "$dir/LncRNA/Assembly/*.Scripture.transcripts.gtf",
		"CPC�������ͳ��" => "$dir/LncRNA/Identify/CPC.txt",
		"CNCI�������ͳ��" => "$dir/LncRNA/Identify/CNCI.txt",
		"CPAT�������ͳ��" => "$dir/LncRNA/Identify/cpat.txt",
		"pfam�������ͳ��" => "$dir/LncRNA/Identify/Pfam.txt",
		"����λ�ù�ϵ��LncRNA�л���Ԥ����" => "$dir/LncRNA/Target/lncRNA_position.target",
		"���ڻ������е�LncRNA�л���Ԥ����"=>"$dir/LncRNA/Target/lncRNA_basepair.target",
		"lncRNA��������" => "$dir/LncRNA/DEG/LncRNA_fpkm.list",
		"������lncRNA���" => "$dir/LncRNA/DEG/*_vs_*/*_vs_*.DEG_final.xls",
		"������lncRNA��Ŀͳ�Ʊ�" => "$dir/LncRNA/DEG/Diff_Lnc.stat",
		"ע�͵Ĳ�����lncRNA�л�������ͳ�Ʊ�" => "$dir/LncRNA/DEG/Diff_lnc.anno.stat",
		"������lncRNA�л���topGO�������" => "$dir/LncRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.xls",
		"������lncRNA�л���KEGG�������" => "$dir/LncRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.stat",
		"DEU���������" => "$dir/mRNA/DEU/*_vs_*/DEU_Result_Final.xls",
		"DEU�������"=>"$dir/mRNA/DEU/*_vs_*/DEXSeqReport/testForDEU.html",
		#"������lncRNA�л����KEGGͨ·ע��"=>"$dir/LncRNA/DEG/*_vs_*/pathway/kegg_map/*.html",
		"����1 ����б�" => "$dir_template/software_used_list.txt",
		"����2 ���ݿ��б�" => "$dir_template/database_used_list.txt",
		"����3 ��������" => "$dir_template/Nucleic_acids_encoding_table.txt",
		"LncRNA��gff�ļ�"=>"$dir/LncRNA/Identify/LncRNA.gff",
		"��֪lncRNA����ͳ�Ʊ�"=>"$dir/Personality/Known_LncRNA/known_lncRNA.result.txt",
		"miRNAǰ���lncRNAͳ��"=>"$dir/Personality/precursor/lncRNA_precursor.txt",
		"miRNA����lncRNA������"=>"$dir/Personality/miRNA_Target2LncRNA/*.mir2target.list",
		"�������ת¼���ӷ�����"=>"$dir/Personality/TF/TF.txt",
		#"SNP/InDelλ����Ϣʾ���"=>"$dir_template/snp_indel.anno.gatk.list",
		#"�ɱ�����¼��ṹ�ͱ����ʾ���"=>"$dir_template/ASprofile.fpkm",
		#"����ṹ�Ż����ֽ��"=>"$dir_template/gene.structure.optimize.txt",
		#"�»����GFF�ļ�ʾ���"=>"$dir_template/newgene.gff",
		#"�»�������FASTA�ļ�ʾ���"=>"$dir_template/newgene.fa",
		#"�����������ʾ���"=>"$dir_template/geneExpression.txt",
		"miRNAǰ���lncRNAͳ��ʾ���"=>"$dir_template/precursor.txt",
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
		if ($desc=~/ԭʼ�������/){
			my $tit="������������"; $tit=&T($tit);
			my $tit_class="h3"; 
			my $ann="�ڽ������ݷ���֮ǰ��������Ҫȷ����ЩReads���㹻�ߵ��������Ա�֤����������׼ȷ�������Ͷ����ݽ����ϸ���������ƣ��������¹��˷�ʽ��";
			$ann .="\n"."�ڽ������ݷ���֮ǰ��������Ҫȷ����ЩReads���㹻�ߵ��������Ա�֤����������׼ȷ�������Ͷ����ݽ����ϸ���������ƣ��������¹��˷�ʽ��";
			$ann .="\n"."(1) ȥ�����н�ͷ��Reads��";
			$ann .="\n"."(2) ȥ����������Reads������ȥ��N�ı�������10%��Reads��ȥ������ֵQ��10�ļ����ռ����Read��50%���ϵ�Reads����";
			$ann .="\n"."��������һϵ�е���������֮��õ��ĸ�������Clean Data����FASTQ��ʽ�ṩ��";
			$ann=&substitute($ann); $ann=&T($ann);
			my $ann1="ע��Adapter related�����˵��ĺ��н�ͷReads��ռ��Raw Reads���ı�����Low quality�����˵��ĵ�����Reads��ռ��Raw Reads���ı�����Clean Reads���������Ϲ��˵õ���Clean Reads ��ռ��Raw Reads ���ı�����";
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
			if ($desc=~/�����Ǳ������ʵ������ͼ/){
				my $pic_name=&T($desc);
				$path = &T($path);
				$type="img-width-normal";
				$type = &T($type);
				$content = "<pic name=$pic_name type=$type path=$path \/>\n";
			}
			elsif ($desc=~/���������ά��ͼ/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="��ÿ�����������ά��ͼ���ƣ�����ͼ��ͼ��չʾ�˸��Ƚ������еĲ������ĸ������Լ��Ƚ����Ĺ��еĲ�����������"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $tit="���������ά��ͼ";	$tit=&T($tit);
				my $tit_class="h3";	
				$content="<$tit_class name=$tit type=$type desc=$tit \/>\n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/>\n";
			}
			elsif ($desc=~/����������������ͼ/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="ע��x���ʾʵ��������y���ʾ��׼����FPKM����ɫ������ʾ���cluster�е����л����ڲ�ͬʵ����������Ա������ƽ��ֵ������ͼ��"; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="��ɸѡ�������в��������ı������K-means����������Բ������ı����ˮƽֵlog2(FPKM+1) ���о������,�õ���������ڲ�ͬʵ�������µı��ģʽ�����ģʽ��ͬ������Ļ���ۼ����ࡣ����ͬ��Ļ�����ܾ������ƵĹ��ܣ����ǹ�ͬ����ͬһ��л���̻�ϸ��ͨ·����ˣ�ͨ����ͬһ��cluster�Ļ�����й���ע�ͼ�������������Ԥ��δ֪����Ĺ��ܻ���֪�����δ֪���ܡ�";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content ="<p type=$type desc=$ann1 \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				
			}
			elsif ($desc=~/������lncRNA����ͼ/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="��ɸѡ���Ĳ�����lncRNA����ξ����������������ͬ�����Ʊ����Ϊ��lncRNA���о��࣬������lncRNA����������ͼ��" ; 
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1 ="ע��ͼ�в�ͬ���д���ͬ����Ʒ����ͬ���д���ͬ��lncRNA����ɫ������lncRNA����Ʒ�еı����ˮƽlog2��FPKM+1����"; 
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				my $tit="������lncRNA�������"; $tit=&T($tit);
				my $tit_class='h3'; 
				$content ="<$tit_class name=$tit type=$type desc=$tit \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<pic name=$picture_name type=$type path=$path \/> \n";
				$content .= "<p type=$type desc=$ann1 \/> \n";
			}
			elsif ($desc=~/������lncRNA�л����������ͼ/){
				my $picture_name=&T($desc);
				$path=&T($path);
				my $ann="ע��x���ʾʵ��������y���ʾ��׼����FPKM����ɫ������ʾ���cluster�е����л����ڲ�ͬʵ����������Ա������ƽ��ֵ������ͼ��";
				$ann=&substitute($ann); $ann=&T($ann);
				my $ann1="��ɸѡ��������lncRNA�ı��������ξ���������õ�lncRNA�ڲ�ͬʵ�������µı��ģʽ�����ģʽ��ͬ�������lncRNA�ۼ����ࡣ";
				$ann1=&substitute($ann1); $ann1=&T($ann1);
				$content = "<pic name=$picture_name type=$type path=$path \/> \n";
				$content .="<p type=$type desc=$ann \/> \n";
				$content .="<p type=$type desc=$ann1 \/> \n";
			}
			elsif ($desc=~/��Ʒ�������ͼ/){
				if (exists $config{'Sep'}){
					my $picture_name=&T($desc);
					my $tit="�ظ����������"; $tit=&T($tit);
					my $tit_class="h3";
					my $ann="�о�����������ı���ڲ�ͬ�ĸ�����������ѧ�ɱ��ԣ�Biological Variability������ͬ�Ļ���֮����Ŀɱ�̶ȴ��ڲ��죬�������Ǳ����������qPCR�Լ�����оƬ�ȼ����������������ֿɱ��ԡ�Ϊ��Ѱ����������Ȥ�Ĳ����������Ҫ���Ǻʹ���������ѧ�ɱ�����ɵı����졣Ŀǰ���������Ч�ķ�������ʵ���������������ѧ�ظ���Biological Replicates�����ظ���������Խ�ϸ��ظ���Ʒ��ĿԽ�࣬Ѱ�ҵ��Ĳ��������Խ�ɿ���������������ѧ�ظ�����Ŀ����������ѧ�ظ�������Զ��ڷ���ת¼��������ݷǳ���Ҫ������ѧ�ظ�������Բ������Լ�������ѧʵ������Ŀ��ظ��ԣ��������������������Ŀɿ��Ժ͸����쳣��Ʒ��ɸ�顣";
					$ann .="\n"."��Ƥ��ѷ���ϵ��r��Pearson��s Correlation Coefficient����Ϊ����ѧ�ظ�����Ե�����ָ�ꡣr2Խ�ӽ�1��˵�������ظ���Ʒ�����Խǿ��";
					my $ann1="ע���������ʾ��Ʒ���ƣ��������ʾ��Ӧ����Ʒ���ơ���ɫ����r2ֵ��С��";
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
		"�����Ǳ������ʵ������ͼ" => "$dir_template/P01_RNA-Seq_experimental_workflow.png",
		"������Ϣ��������ͼ��mRNA���֣�"=> "$dir_template/P02_RNA-Seq_analysis_workflow.png",
		"�����Ǳ���RNA������Ϣ��������ͼ" => "$dir_template/P02_LncRNA_analysis_workflow.png",
		"FASTQ��ʽ�ļ�ʾ��ͼ" => "$dir_template/P03_FASTQ_format.png",
		"��������ʷֲ�ͼ" => "$dir_abs/QC/Clean_Data/PNG/*.quality.png",
		"ATGC�����ֲ�ͼ" => "$dir_abs/QC/Clean_Data/PNG/*.acgtn.png",
		"��ʽ1 ����ֵ���㹫ʽ" => "$dir_template/F01_Qscore_formula.png",
		"��ʽ2 FPKM���㹫ʽ" => "$dir_template/F02_FPKM_formula.png",
		"TopHat2��������" => "$dir_template/TopHat2_workflow.png",
		"Mapped Reads�ڲο��������ϵ�λ�ü�������ȷֲ�ͼ" => "$dir_abs/QC/Map_assess/*.map.png",
		"�����鲻ͬ����Reads�ֲ�ͼ" => "$dir_abs/QC/Map_assess/*.type.png",
		"IGV���������" => "$dir_template/P07_IGV_interface.png",
		"Mapped Reads��mRNA�ϵ�λ�÷ֲ�ͼ" => "$dir_abs/QC/Map_assess/Total.randcheck.png",
		"����Ƭ�γ���ģ��ֲ�ͼ" => "$dir_abs/QC/Map_assess/*.insertSize.r.png",
		"�����Ǳ���������ݱ��Ͷ�ģ��ͼ" => "$dir_abs/QC/Map_assess/*.Saturation.png",
		#"�����Ǳ���������ݱ��Ͷ�ģ��ͼ"=>"$dir_abs/QC/Map_assess/Total.gene_tag.png",
		"SNPͻ�����ͷֲ�ͼ" => "$dir_abs/mRNA/Structure/SNP/All.snp.type.png",
		"SNP�ܶȷֲ�ͼ" => "$dir_abs/mRNA/Structure/SNP/AllSample.SNP_density.png",
		"SNPע�ͷ���ͼ" => "$dir_abs/mRNA/Structure/SNP/snp_anno/all.anno.stat.png",
		"InDelע�ͷ���ͼ" => "$dir_abs/mRNA/Structure/SNP/indel_anno/all.anno.stat.png",
		"�ɱ��������ͳ��ͼ" => "$dir_abs/mRNA/Structure/Alt_splice/*.png",
		"����ƷFPKM�ܶȷֲ��Ա�ͼ" => "$dir_abs/mRNA/Expression/all.fpkm_density.png",
		"����ƷFPKM����ͼ" => "$dir_abs/mRNA/Expression/all.fpkm_box.png",
		"��Ʒ�������ͼ" => "$dir_abs/mRNA/Expression/sample_cluster.png",
		"�������ɽͼ" => "$dir_abs/mRNA/DEG/*_vs_*/*_vs_*.FC_FDR.png",
		"������MAͼ" => "$dir_abs/mRNA/DEG/*_vs_*/*_vs_*.FC_count.png",
		"���������ά��ͼ" => "$dir_abs/mRNA/DEG/All_DEG_veen.png",
		"������������ͼ" => "$dir_abs/mRNA/DEG/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster.png",
		"����������������ͼ" => "$dir_abs/mRNA/DEG/k-means.png",
		"���������GOע�ͷ���ͳ��ͼ" => "$dir_abs/mRNA/DEG/*_vs_*/go_enrichment/*_vs_*.GO.png",
		"���������topGO���������޻�ͼ" => "$dir_abs/mRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.png",
		"���������COGע�ͷ���ͳ��ͼ" => "$dir_abs/mRNA/DEG/*_vs_*/Cog_Anno/*_vs_*.Cog.classfy.png",
		"���������KEGG����ͼ" => "$dir_abs/mRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.png",
		"����������KEGGͨ·ע��ͼ" =>"$dir_abs/mRNA/DEG/*_vs_*/pathway/kegg_map/*.png",
		"���������KEGGͨ·����ɢ��ͼ" => "$dir_abs/mRNA/DEG/*_vs_*/Graph/*_vs_*.KEGG.Phase.png",
		"��������򵰰��ʻ�������ͼ" => "$dir_template/P26_pp_network.png",
		"Ԥ�ⳤ���Ǳ���RNAͳ��ͼ" => "$dir_abs/LncRNA/Identify/lnc.class_code.stat.png",
		"Ԥ�ⷽ��ά��ͼ" => "$dir_abs/LncRNA/Identify/venn.png",
		"DEU�������ʾ��ͼ"=>"$dir_template/DEU.png",
		"������lncRNA����ͼ" => "$dir_abs/LncRNA/DEG/All_DEG/DEG_Cluster/hierarchical/all_sample_DEG_cluster.png",
		"������lncRNA�л����������ͼ"=>"$dir_abs/LncRNA/DEG/k-means.png",
		"������lncRNA�л���GOע�ͷ���ͳ��ͼ" => "$dir_abs/LncRNA/DEG/*_vs_*/go_enrichment/*_vs_*.GO.png",
		"������lnc�л���topGO�����޻�ͼ" => "$dir_abs/LncRNA/DEG/*_vs_*/Graph/*_vs_*.topGO_*.png",
		"������lncRNA�л���COGע�ͷ���ͳ��ͼ" => "$dir_abs/LncRNA/DEG/*_vs_*/Cog_Anno/*_vs_*.Cog.classfy.png",
		"������lncRNA�л���KEGG����ͼ" => "$dir_abs/LncRNA/DEG/*_vs_*/pathway/kegg_enrichment/*_vs_*.KEGG.png",
		"������lncRNA�л���KEGGͨ·����ɢ��ͼ" => "$dir_abs/LncRNA/DEG/*_vs_*/Graph/*_vs_*.KEGG.Phase.png",
		"������lncRNA�л����KEGGͨ·ע��ͼ"=> "$dir_abs/LncRNA/DEG/*_vs_*/pathway/kegg_map/*.png",
		"������lncRNA�л��򵰰׻�������ͼ" => "$dir_template/P26_pp_network_lncRNA.png",
		"mRNA����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/Length/mRNA.len.png",
		"lncRNA����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/Length/lncRNA.len.png",
		"mRNA��exon����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/Exon/mRNA.exon.png",
		"lncRNA��exon����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/Exon/lncRNA.exon.png",
		"mRNA��orf����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/ORF/gene.orf.png",
		"lncRNA��orf����ͳ��ͼ" => "$dir_abs/LncRNA/Compare/ORF/lnc_filter_final.orf.png",
		"lncRNA��mRNA������Ƚ�ͼ" => "$dir_abs/LncRNA/Compare/Expression/lnc_vs_mRNA.fpkm.png",
		"lncRNA��mRNA�ɱ�����칹��Ƚ�ͼ"=>"$dir_abs/LncRNA/Compare/Isoform/lnc_vs_mRNA.isoform.png",
		"code����lncRNA�����б����Է���ͼ"=>"$dir_abs/Personality/Conservation/phastcons.cdf.png",
		"lncRNA��mRNAλ�㱣����ʾ��ͼ"=>"$dir_template/Conservation_position.png",
		"Jensen_Shannon divergence��ʽ"=>"$dir_template/JS.png",
		"lncRNA�����Ա��ͼ"=>"$dir_abs/Personality/Tissue_specific/*.png",
		"����������ɷַ���ͼ"=>"$dir_abs/Personality/PCA/*.png",
		"mRNA��lncRNA������������ʾ��ͼ"=>"$dir_abs/LncRNA/Target/Trans_target/Fig_8_1_networkHeatmap.png",
		"ԭʼ�������"=>"$dir_abs/QC/Clean_Data/PNG/Raw_data_ratio_*.png",
		
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
	if ($name eq 'һ������'){
		$class='h1';
	}
	elsif ($name eq '��������'){
		$class = 'h2';
	}
	elsif ($name eq '��������') {
		$class = 'h3';
	}
	elsif ($name eq '�ļ�����') {
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
	my $content = "<ref_list name=\"�ο�����\" type=$type > \n";
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
	my $content = "<Attach name=\"��¼����\" type=$type  \/> \n";
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
	return sprintf("%4d��%02d��%02d��", $year+1900, $mon+1, $day);
}

sub time_format {
	my $time = shift;
	my ($year,$mon,$day) = (split /\//,$time);
	return (sprintf("%4d��%02d��%02d��", $year, $mon, $day));
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
