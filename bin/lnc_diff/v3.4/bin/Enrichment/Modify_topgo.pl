#!/usr/bin/perl -w
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir,$links);
GetOptions(
				"help|?" =>\&USAGE,
				"links:s"=>\$links,
				"indir:s"=>\$indir,
				) or &USAGE;
&USAGE unless ($indir);


############################

open IN,$links or die $!;

my %link;
while(<IN>){
	chomp;
	s/\r//;
	next if (/^$/);
	if(!/^\#/){
		my @line=split /\t/,$_,2;
		$link{$line[0]}=$line[1];
	}
}
close IN;
my @file;
@file=glob "$indir/*topGO*.xls";
my $length=@file;
if($length==0){
	print STDERR "Error:topGO files is not found!\n";exit;
}
foreach my $val(@file){
	next if $val =~/_gene.xls/;
    my $base=basename($val);
    my $new=(split("_",$base))[-1];
    my $dir=dirname($val);
    `mv $val $dir/$new`;
    open IN,"$dir/$new" or die $!;
    open OUT,">$dir/$base" or die $!;
    print OUT "GO.ID\tTerm\tAnnotated\tSignificant\tExpected\tKS\n";
    <IN>;
    while(<IN>){
	chomp;
	s/\r//;
	next if (/^$/);
	my @tmp=split /\t/,$_;
	print OUT "$tmp[0]\t$link{$tmp[0]}\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\n" if(defined $link{$tmp[0]});
	print OUT "$_\n" if(!defined $link{$tmp[0]});
      }
      close OUT;
      close IN;
	  system "rm $dir/$new";
}





#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath
{		#��ȡָ��Ŀ¼���ļ��ľ���·��
        my ($type,$input) = @_;

        my $return;
	$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:cui ye<cuiy\@biomarker.com.cn> 

Usage:
  Options:
  -indir       target file 
  -links        links file
  -h         Help

USAGE
	print $usage;
	exit;
}



