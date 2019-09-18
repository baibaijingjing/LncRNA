#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use GD;
use GD::Polyline;
use Text::NSP::Measures::2D::Fisher::right;
use Encode;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";   
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")};
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"k:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);
$deg=&ABSOLUTE_DIR($deg);
$fIn=&ABSOLUTE_DIR($fIn);

my $group="$deg/group.dat" if (-f "$deg/group.dat");
my ($anno,$list_file,$pvalue_file,$detail_file,$pre);
unless ($key eq "kegg"){
    $anno=(glob "$fIn/*GO_tree.stat*")[0];
    $list_file="$od/GO_term.list";
    $pvalue_file="$od/DEG_GO_term_pvalue.list";
    $detail_file="$od/DEG_GO_term_detail.list";
    $pre="GO_cluster";
}
unless ($key eq "go"){
    $anno = (glob "$fIn/*Kegg.pathway*")[0];
    $list_file="$od/KEGG_pathway.list";
    $pvalue_file="$od/DEG_KEGG_term_pvalue.list";
    $detail_file="$od/DEG_KEGG_term_detail.list";
    $pre="KEGG_cluster";
}


my (%gene,%pathmy ,%deglist,%kegg,%sort,%path);
###get diff gene id list and hash


open IN,"$deg" or die $!;
while (<IN>){
    chomp;
    next if (/^#/ || /^$/);
    my $gene=(split /\s+/,$_)[0];
    $deglist{$gene}=1;
}
close IN;
my @degid=keys %deglist;
my $deg_num=@degid;

=pod
注意cut的输出可能有换行符
my @degid=`cut -f 1 $deg`;
shift @degid;
#print "shuzu :".Dumper(\@degid);
my $deg_num=@degid;
my %deglist = map{$_ => 1} @degid;
=cut
#print Dumper(\%deglist);
print "total diff gene id nummber :$deg_num\n";
#######################get per sample gene id list

open DE,"$deg" or die $!;
my $head_line=<DE>;
my @Samples=split/\t/,$head_line;
shift @Samples;
my $sam=@Samples;
while (<DE>){
    chomp;
    next if (/^#/ || /^$/);
    my @A=split /\s+/,$_;
    my $id=shift @A;

    for (my $i=0;$i<@A ;$i++) {
        if ($A[$i]>0){
            if (!exists $gene{$Samples[$i]}){
                my @tem=();
                push (@tem,$id);
                $gene{$Samples[$i]}=\@tem;
            }
            if (exists $gene{$Samples[$i]}){
                my $tem=$gene{$Samples[$i]};
                push (@$tem,$id);
                #$gene{$Samples[$i]}=\@tem;
            }
        }
        
    }
}


#print Dumper(\%gene);
#########get pathway hash 
open (KEGG,"$anno") or die $!;


while (<KEGG>) {
	next if /^\#/;
	chomp;
	next if (/Biological Process\s/ || /Cellular Component\s/ || /Molecular Function\s/);
	my @t = split /\t/,$_;
        my $name;
        my @p;
        unless ($key eq "go"){
            $name=$t[0];
            @p = split /;/,($t[3]);
        }
        unless ($key eq "kegg"){
            my @tem=split /\s+/,$t[0];
            pop @tem;
            $name=join(" ",@tem[2..$#tem]);
            @p = split /;/,($t[2]);
        }
        #print "way name:$name\n";
        #print Dumper(\@p);
        $path{$name}=\@p;
	my $nu = 0;#pathways total nummbere
        my (@p_new,@q_new);
        for my $i (0..$#p) {
                $kegg{all}{$p[$i]}=1;
                if (exists $deglist{$p[$i]}) {
                        push @p_new,$p[$i];
                        
                        $nu += 1;
                        $kegg{match}{$p[$i]}=1;
                }
        }
	
	
}
$kegg{num}{match}=keys %{$kegg{match}};
$kegg{num}{all}=keys %{$kegg{all}};
close KEGG;
print "total pathway gene id list : $kegg{'num'}{'all'} \n";
print "total diff gene id list (2):$kegg{num}{match}";

##########out put result

my $N=$kegg{num}{all};
my $M=$kegg{num}{match};

foreach my $way (sort keys %path){
    my %phash;
    my @way_diff=();
    foreach my $sam (sort keys %gene){
	
        my @ways=@{$path{$way}};
        my @way1=@ways;
        my $num_way=@ways;### one pathway contain gene bumber====$n
        my @ge=@{$gene{$sam}};
        my %ways = map{$_ => 1} @ways;
        my %ge = map{$_ => 1} @ge;
        #print "$sam".Dumper(\%ways);
        #print "$sam".Dumper(\%ge);
        my @share = grep {$ge{$_}} @ways;
        my $num=@share;### one pathway diff gene ====$m
        #push @way_diff,$share[$_] for(0..$#share);
        #print "$sam".Dumper(\@share)."\n";
        push @way_diff,@share;
        $kegg{$way}{$sam}{'diff'}=$num;
        $kegg{$way}{$sam}{'total'}=$num_way;
	
    }
    #print Dumper(\@way_diff);
    @way_diff = grep { ++$phash{$_} < 2 } @way_diff;
    my $tm_num=@way_diff;
    $sort{$way}=$tm_num;
}
#print Dumper(\ %sort);

open OUT,">$list_file";
open OUT1,">$pvalue_file" or die $!;
open OUT2,">$detail_file" or die $!;
my $sam_head=join("\t",@Samples);
chomp $sam_head;
print OUT "#path_way\t$sam_head\n";
print OUT1 "#path_way\t$sam_head\n";
print OUT2 "#path_way\t$sam_head\tway_total_diff\tway_gene\tall_diff_gene\tall_gene\n";
foreach my $final (sort{$sort{$b}<=>$sort{$a}} keys %sort){
    next if ($final eq "num" || $final eq "all" || $final eq "match");
    my @len=split /\s+/,$final;
    my $len=@len;
    my ($line,$line1,$line2);
    if ($len > 7){
        my $title=join(" ",@len[0..4])."...";
        $line= "$title($sort{$final})";
        $line1="$title($sort{$final})";
    }else{
        $line= "$final($sort{$final})";
        $line1="$final($sort{$final})";
    }
    
    #my $line1="$final($sort{$final})";
    $line2="$final($sort{$final})";
    my $q=0;
    foreach my $sm (sort keys %gene){
        $q++ if ($kegg{$final}{$sm}{'diff'}==0);
	$line.= "\t$kegg{$final}{$sm}{'diff'}";
	my $p_value = &hyper($kegg{$final}{$sm}{'diff'},$M,$kegg{$final}{$sm}{'total'},$N);
	$line1.= "\t$$p_value";
        $line2.="\t$kegg{$final}{$sm}{'diff'}";
    }
    $line.="\n";
   $line1.= "\n";
   $line2.="\t$sort{$final}\t$kegg{$final}{$Samples[0]}{'total'}\t$M\t$N\n";
   print OUT $line unless ($q==$sam);
   print OUT1 $line1;
   print OUT2  $line2;
}
close OUT;
close OUT1;
close OUT2;
#system "$config{Rscript} $Bin/anno_cluster_heatmap.r --infile $list_file --outfile $od/$pre --rowname T --height 2000 --width 2000 --color \"blue,white,red\" --size 4 --groupfile $group";
system "$config{Rscript} $Bin/anno_cluster_heatmap.r --infile $list_file --outfile $od/$pre --rowname T --height 2000 --width 2000 --color \"blue,white,red\" --size 4 --zero 0 ";
print "$config{Rscript} $Bin/anno_cluster_heatmap.r --infile $list_file --outfile $od/$pre --rowname T --height 2000 --width 2000 --color \"blue,white,red\" --size 4 ";
#$kegg{num}{match}=keys %{$kegg{match}};
#$kegg{num}{all}=keys %{$kegg{all}};


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
sub hyper{ #calculate  Fisher's exact test
	my ($n11,$n1p,$np1,$npp)=@_;
	my ($out,$errorCode,$right_value);
	$right_value = calculateStatistic(
		n11=>$n11,
		n1p=>$n1p,
		np1=>$np1,
		npp=>$npp
	);
	if( ( $errorCode = getErrorCode() ) ) {
		$out = $errorCode." - ".getErrorMessage();
	}
	else {
		$out = $right_value;
	}
	return \$out;
}
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
Contact:	niulg <niulg\@biomarker.com.cn> 
Usage:
  Options:
  -i     <file>  input file,All_Database_annotation.xls,forced 
  
  -deg   <file>  deg file,forced 
  
  -k     <str>   keywords of output file,forced 
  
  -od    <file>  output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
 
