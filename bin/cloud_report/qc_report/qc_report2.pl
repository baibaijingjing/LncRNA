#!usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

my ($od,$cfg,$indir,$data_type,$cut,$sep);
GetOptions(
	"h|?"=>\&USAGE,
	"od:s"=>\$od,
	"cfg:s"=>\$cfg,
	"indir:s"=>\$indir,
	"types:s"=>\$data_type,
	"sep:s"=>\$sep,
	"cut:s"=>\$cut,
)or &USAGE;
&USAGE unless ($od and $cfg and $indir);

$indir=abs_path($indir);
$od=abs_path($od);
mkdir $od unless (-d $od);
$cut ||=15;
$sep ||=2;
$data_type ||= 'pe';

my %CFG;
my %data;
&read_cfg();

my $Project_type=$CFG{'Project_type'};
my $Project_ID=$CFG{'Project_ID'};
my $species_type=$CFG{'species_type'};
my $RefSeq=$CFG{'RefSeq?'};
my $RequiredData=$CFG{'RequiredData'};
my $RequiredQuality=$CFG{'Q30'};
my $FinishedDate=&GetDate;
my $Deadline=$CFG{'Deadline'};

my %BaseSum;
my %Q30;
my %GC_CONTENT;
my @GC_Con;
my $QC_PATH=glob("$indir/Data_Assess/AllSample_GC_Q.stat");
open(IN,$QC_PATH) or die $!;
while(<IN>){
	chomp;
	next if (/^#/ || /^\s*$/);
	my($sample_id,$basesum,$quality,$GC_content)=(split /\s+/,$_)[0,2,5,3];
	$basesum=sprintf "%0.2f",$basesum/1000000000;
	$BaseSum{$sample_id}=$basesum;
	$Q30{$sample_id}=$quality."%";
	$GC_CONTENT{$sample_id}=(sprintf "%.2f",$GC_content)."%";
	push @GC_Con,$GC_content;
}
close IN;

my ($min_GC_content,$max_GC_content)= (sort {$a<=>$b} @GC_Con)[0,-1];
my $GC_div=$max_GC_content-$min_GC_content;
$GC_div=&percent($GC_div/100);

my @base_assess=glob("$indir/Data_Assess/*.acgtn");
foreach my $base_asses (@base_assess){
	my ($sample_id)=basename($base_asses)=~/^(\S+)\.acgtn$/;
	my ($base_iso,$base_wave)=&base_assess($base_asses,$data_type,$sep,$cut);
	$base_iso=&percent($base_iso/100);
	$base_wave=&percent($base_wave/100);
	push @{$data{$sample_id}},($base_iso,$base_wave);
}


my %insert_around_per;
my %insert_mode;
&read_insertSize();

my %randCheck;
&read_randCheck();

my $library_qc=(-d "/share/nas1/basecall/Config_dir/")?"/share/nas1/basecall/Config_dir":"/share/nas7/no_ref/Analysis/Current_Data/heying/Config_Dir";
my %library_info;
&read_library() if (-d $library_qc);

my(%map_ratio,$map_ratio_deviation);
&read_mapstat();

my $lnc_pre_ratio;
#my $none=0;
my $lnc_final=glob("$indir/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_final.fa");
my $insect=`less -S $lnc_final |wc -l`; chomp $insect;
$insect=$insect/2;


my $lnc_pre=glob("$indir/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/merged_filter.fa");
my $union=`less $lnc_pre|wc -l`; chomp($union);
$union=$union/2;
$lnc_pre_ratio=(sprintf("%.2f",$insect/$union*100))."%";

my $lnc_fpkm=glob("$indir/Lnc_Diff_Analysis/LncRNA_fpkm.list");
my $total_lnc=`less $lnc_fpkm|wc -l`;
chomp ($total_lnc);
$total_lnc=$total_lnc-1;
my $non_expressed=0;
open(IN,$lnc_fpkm) or die $!;
while(<IN>){
	chomp;
	next if (/^#/||/^\s*$/);
	my @fpkms=(split /\s+/);
	#for (@fpkms[1..$#fpkms]){$non_expressed++ if $_==0;}
	my @non;
	for (my $i=0;$i<$#fpkms;$i++){$non[$i]=0;}
	#print "@non\n"; 
	#print "@fpkms[1..$#fpkms]\n";
	if(join(" ", @fpkms[1..$#fpkms]) eq join(" ", @non)){
		$non_expressed++;
	}
}
close IN;
print "$non_expressed\n";
my $fpkm_rato=(sprintf("%.2f",(($total_lnc-$non_expressed)/$total_lnc*100)))."%";

my $deg_num_deviation;
&read_deg();

my ($lnc_deg_deviation,$lnc_deg_ratio);
&read_lnc_deg();

my ($known_gene_num,$total_gene_num,$new_gene_num);
&read_gff();

my %as_ratio;
&read_altersplice();


my %sensitivity;
&read_geneexpression();


open(OUT,">$od/$Project_ID.QC.stat") or die $!;
print OUT join("\t",'Project_type','Project_ID','Sample_ID','species_type','RefSeq?','RequiredData','ObtainedData','RequiredQuality','UseData','Qubit','LibraryQuality','BaseInfo','BaseWave','NT%(fungus)','NT%(irrelative)','rRNA%','NEB?','InsertAround%','InserMode','GC%','GC%_Range','MaxDepthPer','Mapped%','Mapped%_range','LncPredic_Insection/Union','Lnc_expressed/Total_transcript','DEG_num','Lnc_DEG_num','Lnc_DEG_num/Lnc_total','AS/Gene(%)','Gene_expressed/Gene_total','NewGene','Correlation_min','Correlation_max','Deadline','FinishedDate')."\n";

foreach my $sample_id(sort keys %map_ratio){
	my $obtainedData=$BaseSum{$sample_id};
	my ($id,$qubit,$libraryquality)=(exists $library_info{$sample_id})?@{$library_info{$sample_id}}:('-','-','-');
	my ($base_info,$base_wave)=@{$data{$sample_id}};
	my $insertAround=$insert_around_per{$sample_id};
	my $inserMode=$insert_mode{$sample_id};
	my $gc=$GC_CONTENT{$sample_id};
	my $rand_check=$randCheck{$sample_id};
	my $map_Ratio=$map_ratio{$sample_id};
	my $ASRatio=$as_ratio{$sample_id};
	my $gene_ratio=$sensitivity{$sample_id};
	print OUT join("\t",$Project_type,$Project_ID,$sample_id,$species_type,$RefSeq,$RequiredData,$obtainedData,$RequiredQuality,$obtainedData,$qubit,$libraryquality,$base_info,$base_wave,'-','-','-','ÊÇ',$insertAround,$inserMode,$gc,$GC_div,$rand_check,$map_Ratio,$map_ratio_deviation,$lnc_pre_ratio,$fpkm_rato,$deg_num_deviation,$lnc_deg_deviation,$lnc_deg_ratio,$ASRatio,$gene_ratio,$new_gene_num,'-','-',$Deadline,$FinishedDate)."\n";
}
close OUT;

sub read_geneexpression{
	my @geneExpression_files=glob("$indir/Basic_Analysis/geneExpression/*.geneExpression.xls");
	foreach my $geneExpression_file(@geneExpression_files){
		my ($sample_id)=basename($geneExpression_file)=~/^(\S+)\.geneExpression\.xls/;
		my $detected_num;
		open(IN,$geneExpression_file) or die $!;
		while(<IN>){
			chomp;
			next if (/^#/||/^\s*$/);
			my @data=split /\t+/;
			$detected_num++ if $data[2]>0;
		}
		close IN;
		my $sensitivity=(sprintf "%.2f",($detected_num/$total_gene_num*100))."%";
		$sensitivity{$sample_id}=$sensitivity;
	}
}


sub read_altersplice{
	my @as_files=glob("$indir/Alitsplice_Analysis/*.fpkm");
	foreach my $as_file(@as_files){
		my ($sample_id)=basename($as_file)=~/^(\S+)\.fpkm$/;
		my $total_as_num=`less $as_file|cut -f 3|sort|uniq|wc -l`;
		chomp $total_as_num;
		my $total_as_ratio=(sprintf "%.2f",(($total_as_num-1)/$known_gene_num)*100)."%";
		$as_ratio{$sample_id}=$total_as_ratio;
	}
}

sub read_gff{
	my $gff=glob("$indir/Basic_Analysis/Tophat_Cufflinks/Ref_Genome/*.gff3");
	$known_gene_num=`less $gff|awk '\$3=="gene"'|wc -l`;
	chomp $known_gene_num;
	my $newgene_gff=glob("$indir/Basic_Analysis/geneExpression/final_track/*.newGene_final.filtered.gff");
	$new_gene_num=`less $newgene_gff|awk '\$3=="gene"'|wc -l`;
	chomp $new_gene_num;
	$total_gene_num=$known_gene_num+$new_gene_num;
	return($known_gene_num,$total_gene_num,$new_gene_num);
}

sub read_lnc_deg{
	my @lnc_deg_files=glob("$indir/Lnc_Diff_Analysis/*vs*/*.DEG_final.xls");
	my @deg_lnc_nums;
	foreach my $lnc_deg_file(@lnc_deg_files){
		my $deg_lnc_num=`less $lnc_deg_file|wc -l`;
		chomp $deg_lnc_num;
		$deg_lnc_num--;
		push @deg_lnc_nums,$deg_lnc_num;
	}
	my ($max_lnc_degnum,$min_lnc_degnum)=(sort {$a<=>$b}@deg_lnc_nums)[-1,0];
	$lnc_deg_deviation=join("-",$min_lnc_degnum,$max_lnc_degnum);
	$lnc_deg_ratio=(sprintf("%.2f",($max_lnc_degnum/$insect*100)))."%";
	return($lnc_deg_deviation,$lnc_deg_ratio);
}

sub read_deg{
	my @deg_files=glob("$indir/DEG_Analysis/*vs*/*.DEG_final.xls");
	my @deg_nums;
	foreach my $deg_file(@deg_files){
		my $deg_num=`less $deg_file|wc -l`;
		chomp $deg_num;
		$deg_num = $deg_num-1;
		push @deg_nums,$deg_num;
	}
	my ($max_deg_num,$min_deg_num)=(sort {$a<=>$b}@deg_nums)[-1,0];
	$deg_num_deviation=join("-",$min_deg_num,$max_deg_num);
}

sub read_mapstat{
	my @mapstat_files=glob("$indir/Basic_Analysis/Tophat_Cufflinks/Map_Stat/*.mappedStat.xls");
	my @map_ratio;
	foreach my $mapstat_file(@mapstat_files){
		my ($sample_id)=basename($mapstat_file)=~/^(\S+)\.mappedStat\.xls/;
		open(IN,$mapstat_file) or die $!;
		my @lines=<IN>;
		close IN;
		my $mapratio=(split /\s+/,$lines[1])[-1];
		chomp $mapratio;
		$map_ratio{$sample_id}=$mapratio;
		chop $mapratio;
		push @map_ratio,$mapratio;
	}
	my ($min_map_ratio,$max_map_ratio)=(sort {$a<=>$b}@map_ratio)[0,-1];
	$map_ratio_deviation=&format_figure($max_map_ratio-$min_map_ratio);
	$map_ratio_deviation=$map_ratio_deviation."%";
}

sub read_randCheck{
	my @randCheck_files=glob("$indir/Basic_Analysis/Tophat_Cufflinks/Map_Stat/*.randcheck_per.list");
	foreach my $randCheck_file (@randCheck_files){
		my ($sample_id)=basename($randCheck_file)=~/^(\S+)\.randcheck_per\.list/;
		$randCheck{$sample_id}=0;
		open(IN,$randCheck_file) or die "Open $randCheck_file failed!\n";
		while(<IN>){
			chomp;
			if(/^[\d\.]+:([\.\d]+)$/){
				$randCheck{$sample_id}=($randCheck{$sample_id}>$1)?$randCheck{$sample_id}:$1;
			}
		}
		close IN;
		$randCheck{$sample_id}=&format_figure($randCheck{$sample_id}).'%';
	}
}

sub read_insertSize{
	my @insertSize_files = glob ("$indir/Basic_Analysis/Tophat_Cufflinks/Map_Stat/*.insertSize.list");
	foreach my $insertSize_file(@insertSize_files){
		open (IN,$insertSize_file) || die "Open $insertSize_file failed!\n";
		my ($sample_id)=basename($insertSize_file)=~/^(\S+)\.insertSize\.list/;
		my %insert_info;
		while(<IN>){
			chomp;
			next unless /^(\d+):(\d+)$/;
			$insert_info{$1}=$2;
		}
		close IN;
		
		my @len=sort {$insert_info{$a}<=>$insert_info{$b}} (keys %insert_info);
		$insert_mode{$sample_id}=$len[-1];
		my ($count,$sum);
		for (my $i=$len[-1]-75;$i<=$len[-1]+75;$i++){
			$insert_info{$i} ||=0;
			$count += $insert_info{$i};
		}
		foreach my $len (keys %insert_info){
			$sum+=$insert_info{$len};
		}
		$insert_around_per{$sample_id}=&percent($count/$sum);
	}
}

sub base_assess{
	my($fIn,$data_type,$sep,$cut)=@_;
	my (%GC,$base_iso,$GC_sep);
	open (IN,$fIn) or die $!;
	while(<IN>){
		chomp;
		next if (/^#/);
		my @in=split /\t/,$_;
		$GC{$in[0]}{A}=$in[1];
		$GC{$in[0]}{C}=$in[2];
		$GC{$in[0]}{G}=$in[3];
		$GC{$in[0]}{T}=$in[4];
	}
	close IN;
	my $sum_AT=0;
	my $sum_GC=0;
	my $num=0;
	if($data_type eq 'se'){
		foreach my $key (keys %GC){
			if($key>$cut){
				$sum_AT += abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC += abs ($GC{$key}{G}-$GC{$key}{C});
				if (abs ($GC{$key}{A}-$GC{$key}{T})>$sep or abs ($GC{$key}{G}-$GC{$key}{C})>$sep) {
					$num++;
				}
			}
		}
		$base_iso=&max($sum_AT,$sum_GC)/(101-$cut);
		$GC_sep=$num/(101-$cut);
		return ($base_iso,$GC_sep);
	}else{
		foreach my $key (keys %GC){
			if($key>$cut and $key<=101){
				$sum_AT +=abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC +=abs ($GC{$key}{G}-$GC{$key}{C});
				if ((abs ($GC{$key}{A}-$GC{$key}{T}))>=$sep or (abs ($GC{$key}{G}-$GC{$key}{C}))>=$sep) {
					$num++;
				}
			}elsif ($key>101+$cut){
				$sum_AT += abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC += abs ($GC{$key}{G}-$GC{$key}{C});
				if (abs ($GC{$key}{A}-$GC{$key}{T})>$sep or abs ($GC{$key}{G}-$GC{$key}{C})>$sep){
					$num++;
				}
			}
		}
		my $aa=$sum_AT > $sum_GC ? $sum_AT : $sum_GC;
		$base_iso=$aa/(202-$cut*2);
		$GC_sep=$num/(202-$cut*2);
		return ($base_iso,$GC_sep);
	}
}

sub read_library{
	opendir (LIBRARY,$library_qc) or die "Opendir $library_qc failed!\n";
	my $flag;
	while(my $file=readdir LIBRARY){
		next unless $file=~/.*\.txt$/;
		$file="$library_qc/$file";
		open(FILE,$file) or die "Open file $file failed!\n";
		while(<FILE>){
			chomp;
			s/\r//;
			next unless (/$Project_ID/&&/$Project_type/);
			my @data=split /\s+/,$_;
			my ($sample_id)=$data[3]=~/^[a-zA-Z0-9]+-([a-zA-Z0-9]+)-.*/;
			my $qubit;
			if ($data[14] ne '-'){
				$qubit=sprintf "%.2f",$data[14];
			}
			else{
				$qubit=$data[14];
			}
			$library_info{$sample_id}=[$data[3],$qubit,$data[19]];
			$flag++;
		}
		close FILE;
	}
	closedir LIBRARY;
	if (scalar (keys %library_info) ==0){
		print "This project library information is lost in library_qc dir!\n";
	}
}

sub read_cfg{
	open(IN,$cfg) or die $!;
	while(<IN>){
		chomp;
		next if (/^#/ || /^\s*$/);
		my ($name,$anno)=split /\s+/;
		$CFG{$name}=$anno;
	}
	close IN;
}

sub percent{
	my $num=shift;
	if($num=~s/\%//){
		$num=sprintf "%.2f",$num;
		$num .='%';
	}else{
		$num=sprintf "%.2f",$num*100;
		$num .='%';
	}
}

####ÕûÊý¸ñÊ½£º3Î»Ò»¸ö¶ººÅ
####Ð¡Êý¸ñÊ½£ºÐ¡ÊýµãºóÁ½Î»
sub format_figure{
	my $figure=shift;
	if(!defined $figure){
		die;
	}
	if($figure=~/\./){
		if($figure==100){
			$figure=100;
		}else{
			$figure=sprintf("%.2f",$figure);
		}
	}else{
		$figure=&Interger_Three_Digit($figure);
	}
	return $figure;
}

#####ÕûÊý¸ñÊ½£º3Î»Ò»¸ö¶ººÅ
sub Interger_Three_Digit{
	my $interger=shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

sub GetDate{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d",$year+1900,$mon+1,$day);
}

sub USAGE{
	my $usage=<<"USAGE";
#######################################################
Usage:
	Options:
	-h		Help
	-od		output directory
	-cfg		configure file
	-indir		Analysis
	-types		data type pe or se (default pe)
	-sep		the separation of AT or GC more than the number (default 2)
	-cut		the number of removing bases per read(default 15)
#######################################################
USAGE
	print $usage;
	exit;
}