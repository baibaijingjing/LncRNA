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
	my $a1=(split /\_/,$id)[0];
	my $a2=(split /\_/,$id)[1];
	my $ID="$a1\_$a2";
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
$count{50}=$count{100}=$count{150}=$count{200}=$count{250}=$count{300}=$count{350}=$count{400}=$count{450}=$count{500}=$count{550}=$count{600}=$count{650}=$count{700}=$count{750}=$count{800}=$count{850}=$count{900}=$count{950}=$count{1000}=0;
open OUT ,">$odir/$key.orf.txt";
for my $a (sort values %L){
        if ($a>0 && $a<50) { $count{50}++;}	if ($a>=50 && $a<100){ $count{100}++;}    
	if ($a>=100 && $a<150) { $count{150}++;}	if ($a>=150 && $a<200) { $count{200}++;}
	if ($a>=200 && $a<250) { $count{250}++;}        if ($a>=250 && $a<300) { $count{300}++;}
	if ($a>=300 && $a<350) { $count{350}++;}        if ($a>=350 && $a<400) { $count{400}++;}
        if ($a>=400 && $a<450) { $count{450}++;}        if ($a>=450 && $a<500) { $count{500}++;}
        if ($a>=500 && $a<550) { $count{550}++;}        if ($a>=550 && $a<600) { $count{600}++;}
        if ($a>=600 && $a<650) { $count{650}++;}        if ($a>=650 && $a<700) { $count{700}++;}
        if ($a>=700 && $a<750) { $count{750}++;}        if ($a>=750 && $a<800) { $count{800}++;}
        if ($a>=800 && $a<850) { $count{850}++;}     if ($a>=850 && $a<900) { $count{900}++;}
        if ($a>=900 && $a<950) { $count{950}++;}     if ($a>=950) { $count{1000}++;}
}
print OUT "50\t$count{50}\n100\t$count{100}\n150\t$count{150}\n200\t$count{200}\n250\t$count{250}\n300\t$count{300}\n350\t$count{350}\n400\t$count{400}\n450\t$count{450}\n500\t$count{500}\n550\t$count{550}\n600\t$count{600}\n650\t$count{650}\n700\t$count{700}\n750\t$count{750}\n800\t$count{800}\n850\t$count{850}\n900\t$count{900}\n950\t$count{950}\n>=1000\t$count{1000}\n";
close OUT;
`$config{Rscript} $Bin/simpleBar.r --infile $odir/$key.orf.txt --outfile $odir/$key.orf.png --x.col 1 --y.col 2 --x.lab $key.orf.length --y.lab count --axis.size 10 --no.grid`;
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
