my $ver="1.0.0";
use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
my $programe_dir=basename($0);
my $path=dirname($0);

our %opts;
GetOptions(\%opts,"i=s","out=s","od=s","enrichment=s","anno=s","h" );


if( !defined($opts{i}) || !defined($opts{od}) || !defined($opts{anno})|| defined($opts{h})){
	&help();
	exit;
}

my $dir=&ABSOLUTE_DIR($opts{i});
my $outdir=&ABSOLUTE_DIR($opts{od});
my $anno=$opts{anno};
my (@DIR,@DEG,$out,$out1);
if ($anno eq "trans"){
	@DEG=glob("$dir/*vs*/*.DEG_final.xls");
	@DIR=glob("$dir/*vs*/Trans_Anno_enrichment/*.annotation.xls") ;
	$out="$outdir/DEG_lncRNA.trans_anno.stat";
	$out1="$outdir/DEG_lncRNA.stat";
	
}
if ($anno eq "cis"){
	@DEG=glob("$dir/*vs*/*.DEG_final.xls");
	@DIR=glob("$dir/*vs*/Cis_Anno_enrichment/*.annotation.xls") ;
	$out="$outdir/DEG_lncRNA.cis_anno.stat";
	$out1="$outdir/DEG_lncRNA.stat";
}
if ($anno eq "all"){
	@DEG=glob("$dir/*vs*/*.DEG_final.xls");
	@DIR=glob("$dir/*vs*/Anno_enrichment/*.annotation.xls") ;
	$out="$outdir/DEG.anno.stat";
	$out1="$outdir/DEG.stat";
}
my %Site;
my %Anno;
#my @DIR=glob("$dir/Cis_Anno_enrichment/*.annotation.xls") ;
foreach my $file (@DIR){
	my $nam=basename($file);
	$nam=~s/.annotation.xls//;
	open (IN,$file) or die $!;
	while (<IN>) {
		chomp;
		if (/^\#/) {
			my @Anno=split/\s+/,$_;
			for (my $s=0;$s<@Anno ;$s++) {
				if ($Anno[$s] eq 'COG_class') {
					$Site{'COG'}=$s;
				}
				if ($Anno[$s] eq 'KOG_class') {
					$Site{'KOG'}=$s;
				}
				if ($Anno[$s] eq 'Swissprot_annotation') {
					$Site{'Swiss-Prot'} = $s;
				}
				elsif ($Anno[$s]=~/^([^_]+)_annotation/) {
					$Site{$1}=$s;
				}
			}
			$Anno{$nam}{'Annotated'} = 0;
		}
		else{
			my @Info=split /\t+/,$_;
			foreach my $key (keys %Site) {
				$Anno{$nam}{$key}||=0;
				$Anno{$nam}{$key}++ unless ($Info[$Site{$key}] eq '--');
			}
			$Anno{$nam}{'Annotated'} ++;
		}
	}
	close IN;
}
#print Dumper(\%Anno);
my %Stat;
foreach my $deg (@DEG){
	my $nam=basename($deg);
	#print "$nam\n";
	$nam=~s/.DEG_final.xls//;
	#print "$nam\n";
	open (IN,$deg) or die $!;
	while (<IN>) {
		chomp;
		next if /^\#/;
		my $type=(split/\s+/,$_)[-1];
		$Stat{$nam}{up}++ if $type eq 'up';
		$Stat{$nam}{down}++ if $type eq 'down';
		$Stat{$nam}{total}++;
	}
close IN;
}



open (OUT,">$out") or die $!;
my $limit_anno=0;
foreach my $key (sort keys %Anno) {
    if ($Anno{$key}{'Annotated'}!=0){
	  if ($limit_anno==0) {
		  print OUT "#DEG Set";
		  foreach my $key1 (sort keys %{$Anno{$key}}) {
			  print OUT "\t$key1";
		  }
		  print OUT "\n";
		  $limit_anno++;
	  }
	  print OUT "$key";
	  foreach my $key1 (sort keys %{$Anno{$key}}) {
		  print OUT "\t$Anno{$key}{$key1}";
	  }
	  print OUT "\n";
  }else{
	  if ($limit_anno==0) {
		  print OUT "#DEG Set\tAnnotated\tCOG\tGO\tKEGG\tSwiss-prot\tTrEMBL\tnr\tnt\n";
		  $limit_anno++;
	  }
		  print OUT "$key\t0\t0\t0\t0\t0\t0\t0\t0\n";
	  }
}
close OUT;

open (OUT,">$out1") or die $!;
	print OUT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		if ($Stat{$key}{total}!=0){
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}else{
		print OUT "$key\t0\t0\t0\n";
	}
	}
close OUT;



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

sub help{
print <<"Usage End.";
Description: Differencial Express Analysis pipeline;
Version: $ver

Usage:
-i                lnc_diff dir                             must be given;
-anno              cis or trans              must be given;
-od               LncRNA/DEG                                         must be given;
-h                help document
Usage End.
exit;
}
