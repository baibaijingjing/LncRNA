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
my $BEGIN_TIME=time();

#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($flinks,$unigene,$ftest,$od,$protein,$blast,$id,$clade_or_taxid,$cfg);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$od,
				"clade:s"=>\$clade_or_taxid,
				"uni:s"=>\$unigene,
				"id:s"=>\$id,
				"cfg:s"=>\$cfg,
			) or &USAGE;
&USAGE unless ( $unigene and $od and $id);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
#my $STRING_db = "/share/nas2/nas/nas1/database/STRING/v10"; 
my $STRING_db = $config{STRING};
my ($protein_links, $protein_actions, $protein_seq);

if ($clade_or_taxid eq 'all') {
    $protein_links = "$STRING_db/protein.links.v10.txt";
    $protein_actions = "$STRING_db/protein.actions.v10.txt";
    $protein_seq = "$STRING_db/protein.sequences.v10.fa";
} elsif ($clade_or_taxid =~/^\d+$/) {
    $protein_links = "$STRING_db/split/protein_links/eukaryota.$clade_or_taxid.protein.links.txt";
    $protein_actions = "$STRING_db/split/protein_actions/eukaryota.$clade_or_taxid.protein.actions.txt";
    $protein_seq = "$STRING_db/split/protein_seq/eukaryota.$clade_or_taxid.protein.sequences.fa";
} else {
    $protein_links = "$STRING_db/split/protein_links/$clade_or_taxid.protein.links.txt";
    $protein_actions = "$STRING_db/split/protein_actions/$clade_or_taxid.protein.actions.txt";
    $protein_seq = "$STRING_db/split/protein_seq/$clade_or_taxid.protein.sequences.fa";
}

unless (-f $protein_actions) {
    print "ERROR: Illegal arguments of --clade is given! \n";
    print "ERROR: the taxid $clade_or_taxid is illegal, maybe it isn't an eukaryota, or not exist actually.\n\n" if ($clade_or_taxid =~/^\d+$/);
    &USAGE;
}

#$flinks="/share/nas15/database/STRING/v9.1/protein.links.v9.1.txt";
#$protein="/share/nas15/database/STRING/v9.1/protein.sequences.v9.1.fa";

system "mkdir -p $od" unless (-d $od);
$id=&ABSOLUTE_DIR($id);
$od=&ABSOLUTE_DIR($od);
$cfg=&ABSOLUTE_DIR($cfg);
my $name=basename$unigene;

if (-d $id ) {
	` cat $id/*_vs_*/*DEG*final.xls |grep -v "#" |cut -f 1|sort |uniq >$od/used_gene.list `;
}
elsif (-f $id) {
	`cat $id |grep -v "#" |cut -f 1|sort |uniq >$od/used_gene.list `;
} 

my $cmd="perl $Bin/util/extract_list.pl -od $od -data $unigene -list $od/used_gene.list -type extract ";
print `$cmd ` unless (-e "$od/used_gene.list.fa") ;


$cmd="perl $Bin/util/blast_process.pl -cfg $cfg -od $od -fa $od/used_gene.list.fa  -e 1e-5 -p blastx -no $protein_seq ";
print `$cmd ` unless (-e "$od/used_gene.list.fa.blastx.all.xls") ;



$blast="$od/used_gene.list.fa.blastx.all.xls";
#$blast="/share/bioCloud/wangyj/research/DGE_pipeline/v1.2/Maize.Unigene.fa.blastx.best1.xls";


open (OUT,">$od/used_blastx_result.xls") or die $!;
open (IN,$blast) or die $!;
<IN>;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=(split /\t/,$_)[0,1,4,5,8,11];
	$lines[4]=~s/.*\((\S+)\)/$1/;
	if (($lines[4] >75) ) {
		if ( ($lines[5] > 0.7*$lines[1]) || ($lines[5] >0.7*$lines[3]) ) {
			print OUT $_,"\n";
		}	
	}
}
close IN;
close OUT;

$blast="$od/used_blastx_result.xls";
my $length=`less $od/used_blastx_result.xls|wc -l`;
exit if ($length < 2);
my $species=`cut -f 5 $blast|cut -d '.' -f 1 |sort |uniq -c |sort -rn |head -1 `;
	$species =~s/^\s+//g;
	(undef,$species)=split(/\s+/,$species);


my %blast_link;
my %temp;
open (IN,$blast) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @lines=(split /\t/,$_)[0,1,4,5,8,11,12];
	next unless ($lines[2]=~/^$species/) ;
	push @{$blast_link{$lines[2]}},$lines[0],$lines[6];
=c
	if (exists $blast_link{$lines[2]} ) {
		next if $temp{$lines[2]}<$lines[-1];
		$blast_link{$lines[2]}=$lines[0];
	}
	else {
		$temp{$lines[2]}=$lines[-1];
		$blast_link{$lines[2]}=$lines[0];
	}
=cut
}
close IN;

open (IN,$protein_links) or die $!;
#my %LINK;
my %new_LINK;

<IN>;
while (<IN>) {
	chomp;
	next if (/^$/);
	my ($one,$two,$three)=split /\s+/,$_,3;
	next unless ($one=~/$species/ and $two=~/$species/) ;
	$new_LINK{$one}{$two}=$three;
}
close IN;

###########################################
if (-d $id){
foreach my $dir (glob "$id/*_vs_*") {
	my $group_name=basename $dir;
	my %deg=();
	my @genes=glob("$dir/$group_name.DEG*final.xls");
	my $len=@genes;
	next if ($len==0);
	open (LIST,$genes[0])||die $!;
	while(<LIST>){
		chomp;
		next if (/^\s*$|^#/);
		my ($gene)=split /\s+/,$_;
		$deg{$gene}=1;	
	}
	close(LIST);
my %test;
	open (OUT2,">$dir/$group_name.DEG.CytoscapeInput.txt")||open (OUT2,">$od/$group_name.DEG.CytoscapeInput.txt");
	foreach my $former(keys %new_LINK){
		foreach my $latter (sort {$new_LINK{$former}{$b}<=>$new_LINK{$former}{$a}} keys %{$new_LINK{$former}}){
			if (exists $blast_link{$former} and exists $blast_link{$latter}){
				my @form=@{$blast_link{$former}};
				my @latt=@{$blast_link{$latter}};
				foreach my $temp1 (@form) {
					foreach my $temp2 (@latt) {
						if (exists $deg{$temp1} and exists $deg{$temp2}){
							next if (exists $test{$temp1}{$temp2});
							$test{$temp1}{$temp2}=1;
							next if (exists $test{$temp2}{$temp1});
							print OUT2 "$temp1\t$former\t$temp2\t$latter\t$new_LINK{$former}{$latter}\n";

						}
					next;next;
					}
				next;next;
					}
				}
			}
	
		}
	#close(OUT);
	close(OUT2);
 }
}
elsif(-f $id){
	my %deg=();
	open (LIST,$id)||die $!;
	while(<LIST>){
		chomp;
		next if (/^\s*$|^#/);
		my ($gene)=split /\s+/,$_;
		$deg{$gene}=1;	
	}
	close(LIST);
	
	my %test;	
	my $nam=basename$id;
	open (OUT,">$od/$nam.CytoscapeInput.txt")||die $!;

	foreach my $former(keys %new_LINK){
		foreach my $latter (sort {$new_LINK{$former}{$b}<=>$new_LINK{$former}{$a}} keys %{$new_LINK{$former}}){
			if (exists $blast_link{$former} and exists $blast_link{$latter}){

				foreach my $temp1 (@{$blast_link{$former}}) {
					foreach my $temp2 (@{$blast_link{$latter}}) {
						if (exists $deg{$temp1} and exists $deg{$temp2}){
							next if (exists $test{$temp1}{$temp2});
							$test{$temp1}{$temp2}=1;
							next if (exists $test{$temp2}{$temp1});
							#print OUT "$temp1\t$former\t$temp2\t$latter\t$new_LINK{$former}{$latter}\n";
							print OUT2 "$temp1\t$former\t$temp2\t$latter\t$new_LINK{$former}{$latter}\n";

						}
					}
				}
				}
			}
	
		}
	close(OUT);

}
########################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

################# bioinfor_pipeline.log ###################
&show_log2("step_4: Difference expression genes analysis of mRNA finished.");
###########################################################
#close LOG;
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

#############################################################################################################

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

sub USAGE {
	my $usage = <<"USAGE";
 ProgramName:
     Version:	$version
     Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:	2015.09.06
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-uni <file>		fa file,Unigene ,Known and New gene, forced
		-id <dir>		DEG_Analysis directory,like ./DEG_Analysis/,
		                or gene list file,forced
	
		--clade <STR>   organism category or certain eukaryota\'s taxid to transfer interaction
                        'all' | 'archaea' | 'bacteria' | 'eukaryota' | TAXID,           ['eukaryota']
                        eg. 'eukaryota' or 9606 for Homo sapiens
                        [tip1: option 'all' is time-consuming! ]
                        [tip2: your can search TAXID of a species in file $STRING_db/species.taxid.txt!]
		-od <dir>		output directory,forced
		-h		help

USAGE
	print $usage;
	exit;
}
