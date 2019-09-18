#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};

##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %opts;
GetOptions(\%opts,"fa=s","od=s","h");
if (!defined($opts{fa}) || !defined($opts{od}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));

######################
my $notename=`hostname`;chomp $notename;
my $seq_fa=&ABSOLUTE_DIR ($opts{fa});
my $name=basename ($seq_fa);
$name=~m/(.*)\.fa/;
my $key=$1;
`mkdir $opts{od} ` if(!-d $opts{od});
my $odir=&ABSOLUTE_DIR ($opts{od});

chdir "$odir";
`getorf -find 1 $seq_fa $key.orf `;

#my $orf=glob "$odir/*.orf";
my $orf="$odir/$key.orf";
open IN ,"$orf" || die;
$/='>';
<IN>;
#my %len;
my %L;
while(<IN>){
	my ($head,$seq)=split/\n+/,$_,2;
	my $a=(split /\s+/,$head)[1];
	my $b=(split /\s+/,$head)[3];
	my $start=(split /\[/,$a)[1];
	my $end=(split /\]/,$b)[0];
	my $length=(abs($start-$end)+1)/3;
	my $id=(split /\s+/,$head)[0];
	my @idd=(split /\_/,$id);
        my $ID=join ("\_",@idd[0..$#idd-1]);
	#my $ID=(split /\_/,$id)[0];
	if (exists $L{$ID}){
		$L{$ID}=&max($L{$ID},$length);
	}else {
		$L{$ID}=$length;
	}
	#my $leng=$end-$start;
	#if($leng>0){ $length=($leng+1)/3;}
	#if (exists $len{$length}){
	#	$len{$length}++;
	#}else { $len{$length}=1	}
}
close IN;
open OUT1,">$odir/$key.orf.len.xls";
for my $l(sort keys %L){
	print OUT1 "$l\t$L{$l}\n";
}
close OUT1;
my %count;
$count{100}=$count{200}=$count{300}=$count{400}=$count{500}=$count{600}=$count{700}=$count{800}=$count{900}=$count{1000}=$count{1100}=$count{1200}=$count{1300}=$count{1400}=$count{1500}=$count{1600}=$count{1700}=$count{1800}=$count{1900}=$count{2000}=$count{2100}=0;
open OUT ,">$odir/$key.orf.txt";
for my $a (sort values %L){
	if ($a>0 && $a<100){ $count{100}++;}	if ($a>=100 && $a<200) { $count{200}++;}
	if ($a>=200 && $a<300) { $count{300}++;}	if ($a>=300 && $a<400) { $count{400}++;}
	if ($a>=400 && $a<500) { $count{500}++;}	if ($a>=500 && $a<600) { $count{600}++;}
	if ($a>=600 && $a<700) { $count{700}++;}	if ($a>=700 && $a<800) { $count{800}++;}
	if ($a>=800 && $a<900) { $count{900}++;}	if ($a>=900 && $a<1000) { $count{1000}++;}
	if ($a>=1000 && $a<1100) { $count{1100}++;}	if ($a>=1100 && $a<1200) { $count{1200}++;}
	if ($a>=1200 && $a<1300) { $count{1300}++;}	if ($a>=1300 && $a<1400) { $count{1400}++;}
	if ($a>=1400 && $a<1500) { $count{1500}++;}	if ($a>=1500 && $a<1600) { $count{1600}++;}
	if ($a>=1600 && $a<1700) { $count{1700}++;}	if ($a>=1700 && $a<1800) { $count{1800}++;}
	if ($a>=1800 && $a<1900) { $count{1900}++;}	if ($a>=1900 && $a<2000) { $count{2000}++;}
	if ($a>=2000){ $count{2100}++;}
}
print OUT "100\t$count{100}\n200\t$count{200}\n300\t$count{300}\n400\t$count{400}\n500\t$count{500}\n600\t$count{600}\n700\t$count{700}\n800\t$count{800}\n900\t$count{900}\n1000\t$count{1000}\n1100\t$count{1100}\n1200\t$count{1200}\n1300\t$count{1300}\n1400\t$count{1400}\n1500\t$count{1500}\n1600\t$count{1600}\n1700\t$count{1700}\n1800\t$count{1800}\n1900\t$count{1900}\n2000\t$count{2000}\n>=2000\t$count{2100}\n";
close OUT;
`$config{Rscript} $Bin/simpleBar.r --infile $odir/$key.orf.txt --outfile $odir/$key.orf.png --x.col 1 --y.col 2 --x.lab $key.orf.length --y.lab count --axis.size 10`;
#print %len;

#open IN1,"$seq_fa" || die;
#$/='>';
#<IN1>;
#my %seq_length;
#my %seq_line;
#while (<IN1>) {
#	chomp;
#	next if (/^$/ || /^\#/);
#	my ($head,$seq)=split/\n+/,$_,2;
#	my $id=(split/\s+/,$head)[0];
#	$seq=~s/\s+//g;
#	$seq_line{$id}=$seq;
#	my $length=length($seq);
#	$seq_length{$id}=$length;
#}

#close IN1;

#open IN2,"$orf" || die;
#<IN2>;
#my %cds_pos;
#my %pep;
#my %cds_strand;
#while (<IN2>) {
#	chomp;
#	next if (/^$/ || /^\#/);
#	my ($head,$seq)=split /\n+/,$_,2;
#	$seq=~s/\s+//g;
#	my $length=length($seq);
#	my @tmp=split/\s+/,$head;
#	$tmp[0]=~/(.*)_\d+/;
#	my $id=$1;
#	$tmp[1]=~/\[(\d+)/;
#	my $start_pos=$1;
#	$tmp[3]=~/(\d+)\]/;
#	my $end_pos=$1;
#	my $strand="+";
#	if ($head=~/REVERSE\sSENSE/) {
#		$strand="-";
#	}
#	if (!defined $pep{$id}) {
#		$pep{$id}=$seq;
#		$cds_pos{$id}{$start_pos}=$end_pos;
#		$cds_strand{$id}=$strand;
#	}
#	elsif (defined $pep{$id} && $length>length ($pep{$id})) {
#		$pep{$id}=$seq;
#		foreach my $s (keys %{$cds_pos{$id}}) {
#			delete $cds_pos{$id};
#		}
#		$cds_pos{$id}{$start_pos}=$end_pos;
#		$cds_strand{$id}=$strand;
#	}
#}
#$/="\n";
#close IN2;

#open OUT1,">$odir/$key.pep.fa" || die;
#open OUT2,">$odir/$key.cds.fa" || die;
#open OUT3,">$odir/$key.cds_pep.stat.xls" || die;
#foreach my $id (keys %pep) {
#	print OUT1 ">$id\n$pep{$id}\n";
#	my $cds;
#	foreach my $start (keys %{$cds_pos{$id}}) {
#		if ($cds_strand{$id} eq "+") {
#			print OUT3 ">$id\tlength=$seq_length{$id}\tstrand=\'$cds_strand{$id}\'\tstart=$start\tend=$cds_pos{$id}{$start}\n";
#			$cds=substr($seq_line{$id},$start-1,$cds_pos{$id}{$start}-$start+1);
#		}
##=		else {
#			my $cds_seq=substr($seq_line{$id},$cds_pos{$id}{$start}-1,$start-$cds_pos{$id}{$start}+1);
#			print OUT3 ">$id\tlength=$seq_length{$id}\tstrand=\'$cds_strand{$id}\'\tstart=$cds_pos{$id}{$start}\tend=$start\n";
#			$cds=reverse $cds_seq;
#			$cds=~tr/ATCG/TAGC/;
#		}
#	}
#	print OUT2 ">$id\n$cds\n";
#	print OUT3 "$cds\n$pep{$id}\n";
#}

#close OUT1;
#close OUT2;
#close OUT3;
sub max{
	my $a=shift;
	my $b=shift;
	my $tmp;
	if ($a>=$b){
		$tmp=$a;
	}else{
		$tmp=$b
	}
	return $tmp;
}

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);

###########subs
sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

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

sub help
{
	print <<"	Usage End.";
	Description:
		Function : use Getorf predict transcript cDNA and the Pep sequence;
		Version  : $ver.
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-fa
		   Unigene seq fa;
		-od
		   pep seq Out dir ;
		-h
		    Help document
		Attention The Program Cat't Run Backplat;
	Usage End.

	exit;
}
