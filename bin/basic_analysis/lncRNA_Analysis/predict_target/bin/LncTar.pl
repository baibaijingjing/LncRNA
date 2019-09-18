#!/usr/bin/perl -w 
use strict;
use warnings;
use Getopt::Std;
use List::Util qw[min max];
use FindBin qw($Bin $Script);
use lib ("$Bin");
use Genecare::PerlOligos;
my %opt=();
my $lncprimerfile='';
my $mRNAprimerfile='';
my $inputfile='';
my $delta=0;
my $kind=0;
my @res = ();

my $pdm_result = "";
my $whetherpic='';
my ($fprimer,$rprimer);
my @primers;
my @temp;
my @para;
my @pdm_result;
my $list_len = 0;
my $lncPrimerCount=0;

my $usage_1="Usage: perl LncTar.pl [-p kind] [-l lncfile] [-m mrnafile] [-d delta] [-s WhetherPic ] [-o outfile] \n";
my $usage_2="Usage: perl LncTar.pl [-p kind] [-f filename] [-d delta] [-s WhetherPic]  [-o outfile] \n";

if(!(@ARGV==12 || @ARGV==10))
{
   print "param error ,please check out and run again\n";
   print "$usage_1";
   print "$usage_2";
   exit;
}


if(@ARGV==10)
{ 
   if($ARGV[0] eq '-p')
    { $para[1]=$ARGV[1];}
   if($ARGV[2] eq '-p')
    { $para[1]=$ARGV[3];}
   if($ARGV[4] eq '-p')
    { $para[1]=$ARGV[5];}
   if($ARGV[6] eq '-p')
    { $para[1]=$ARGV[7];}
   if($ARGV[8] eq '-p')
    { $para[1]=$ARGV[9];}   
    if($ARGV[0] eq '-f')
    { $para[3]=$ARGV[1];}
    if($ARGV[2] eq '-f')
    { $para[3]=$ARGV[3];}
    if($ARGV[4] eq '-f')
    { $para[3]=$ARGV[5];}
    if($ARGV[6] eq '-f')
    { $para[3]=$ARGV[7];}
    if($ARGV[8] eq '-f')
    { $para[3]=$ARGV[9];}   
    if($ARGV[0] eq '-s')
    { $para[5]=$ARGV[1];}
    if($ARGV[2] eq '-s')
    { $para[5]=$ARGV[3];}
    if($ARGV[4] eq '-s')
    { $para[5]=$ARGV[5];}
    if($ARGV[6] eq '-s')
    { $para[5]=$ARGV[7];}  
    if($ARGV[8] eq '-s')
    { $para[5]=$ARGV[9];}      
    if($ARGV[0] eq '-d')
    { $para[7]=$ARGV[1];}
    if($ARGV[2] eq '-d')
    { $para[7]=$ARGV[3];}
    if($ARGV[4] eq '-d')
    { $para[7]=$ARGV[5];}
    if($ARGV[6] eq '-d')
    { $para[7]=$ARGV[7];}
    if($ARGV[8] eq '-d')
    { $para[7]=$ARGV[9];}   
    if($ARGV[0] eq '-o')
    { $para[9]=$ARGV[1];}
    if($ARGV[2] eq '-o')
    { $para[9]=$ARGV[3];}
    if($ARGV[4] eq '-o')
    { $para[9]=$ARGV[5];}
    if($ARGV[6] eq '-o')
    { $para[9]=$ARGV[7];}
    if($ARGV[8] eq '-o')
    { $para[9]=$ARGV[9];}
}
if(@ARGV==12)
{ 
   if($ARGV[0] eq '-p')
    { $para[1]=$ARGV[1];}
   if($ARGV[2] eq '-p')
    { $para[1]=$ARGV[3];}
   if($ARGV[4] eq '-p')
    { $para[1]=$ARGV[5];}
   if($ARGV[6] eq '-p')
    { $para[1]=$ARGV[7];}
   if($ARGV[8] eq '-p')
    { $para[1]=$ARGV[9];}
   if($ARGV[10] eq '-p')
    { $para[1]=$ARGV[11];}  
    if($ARGV[0] eq '-l')
    { $para[3]=$ARGV[1];}
    if($ARGV[2] eq '-l')
    { $para[3]=$ARGV[3];}
    if($ARGV[4] eq '-l')
    { $para[3]=$ARGV[5];}
    if($ARGV[6] eq '-l')
    { $para[3]=$ARGV[7];}
    if($ARGV[8] eq '-l')
    { $para[3]=$ARGV[9];}
   if($ARGV[10] eq '-l')
    { $para[3]=$ARGV[11];} 
    if($ARGV[0] eq '-m')
    { $para[5]=$ARGV[1];}
    if($ARGV[2] eq '-m')
    { $para[5]=$ARGV[3];}
    if($ARGV[4] eq '-m')
    { $para[5]=$ARGV[5];}
    if($ARGV[6] eq '-m')
    { $para[5]=$ARGV[7];}
    if($ARGV[8] eq '-m')
    { $para[5]=$ARGV[9];}
   if($ARGV[10] eq '-m')
    { $para[5]=$ARGV[11];}  
    if($ARGV[0] eq '-d')
    { $para[7]=$ARGV[1];}
    if($ARGV[2] eq '-d')
    { $para[7]=$ARGV[3];}
    if($ARGV[4] eq '-d')
    { $para[7]=$ARGV[5];}
    if($ARGV[6] eq '-d')
    { $para[7]=$ARGV[7];}
    if($ARGV[8] eq '-d')
    { $para[7]=$ARGV[9];}
   if($ARGV[10] eq '-d')
    { $para[7]=$ARGV[11];}
     if($ARGV[0] eq '-s')
    { $para[9]=$ARGV[1];}
    if($ARGV[2] eq '-s')
    { $para[9]=$ARGV[3];}
    if($ARGV[4] eq '-s')
    { $para[9]=$ARGV[5];}
    if($ARGV[6] eq '-s')
    { $para[9]=$ARGV[7];}
    if($ARGV[8] eq '-s')
    { $para[9]=$ARGV[9];}
   if($ARGV[10] eq '-s')
    { $para[9]=$ARGV[11];}
     if($ARGV[0] eq '-o')
    { $para[11]=$ARGV[1];}
    if($ARGV[2] eq '-o')
    { $para[11]=$ARGV[3];}
    if($ARGV[4] eq '-o')
    { $para[11]=$ARGV[5];}
    if($ARGV[6] eq '-o')
    { $para[11]=$ARGV[7];}
    if($ARGV[8] eq '-o')
    { $para[11]=$ARGV[9];}
   if($ARGV[10] eq '-o')
    { $para[11]=$ARGV[11];}
}
if(!($para[1]==1 ||$para[1]==2))
{
	print "the first param should be 1 or 2!\n";
	exit;
}
if($para[1]==1)    
{
	if(@ARGV!=12)
	{
		 print "param error.\n";
		 print "$usage_1";
		 exit;
  }
        $lncprimerfile=$para[3];
        $mRNAprimerfile=$para[5];
        $delta=$para[7];
        $whetherpic=$para[9];
        $pdm_result=$para[11];
   if((!$lncprimerfile)||(!$mRNAprimerfile)||(!$delta) ||(!$whetherpic)||(!$pdm_result))
   { print "param error,please check carefully.\n";
		 print "$usage_1";
		 exit;
  }
        
}
if($para[1]==2)  
{
	if(@ARGV!=10)
	{
		 print "param error.\n";
		 print "$usage_2";
		 exit;
        }
        $lncprimerfile=$para[3];               
        $whetherpic=$para[5];
        $delta=$para[7];
        $pdm_result=$para[9];
        if((!$lncprimerfile)||(!$delta) ||(!$whetherpic)||(!$pdm_result))
   { print "param error,please check carefully.\n";
		 print "$usage_2";
		 exit;
  }
        
}
open (LOGF, ">$pdm_result")||die"can not open pdm_result!\n"; 
close LOGF;     #clear the file of pdm_result

if($para[1]==1)
{
	      push (@res,"Query\t\t Length_Query\t\t Target\t\t Length_Target\t\t dG\t\t ndG\t\t Start_Position_Query\t\t End_Position_Query\t\t Start_Position_Target\t\t End_Position_Target\n");
	      open(PDM,">>$pdm_result") or die "Can't open file $pdm_result : No such file or directory!\n";
			  print PDM @res;	
			  close (PDM); 
open (PRIMER_FILE, "<$lncprimerfile") or die "Can't open file $lncprimerfile : No such file or directory!\n";
while (<PRIMER_FILE>) {	   
	$list_len=0;
        my $line = $_;
	my @temp = split (/[\t ]/, $line);
	 if($temp[1] eq '')
  {
  	print "the format of file (-l) is error,Please check.\n ";
  	exit;
  }
        $temp[1]=~s///g; 
        $temp[1]=~s/ //g; 
        $temp[1]=~s/u/t/g;
		    $temp[1]=~s/U/T/g;
        $temp[1]=~s/\n//g; 
                       
	push (@primers, [@temp]); 
	$list_len++;
	open (MRNA_FILE,"<$mRNAprimerfile") or die "Can't open file $mRNAprimerfile : No such file or directory!\n";
	while(<MRNA_FILE>)
	{
		my $line1 = $_;
		my @temp1 = split (/[\t ]/, $line1);
			if($temp1[1] eq '')
  {
  	print "the format of file (-m) is error,Please check.\n ";
  	exit;
  }  
		$temp1[1]=~s///g;
		$temp1[1]=~s/u/t/g;
		$temp1[1]=~s/U/T/g;
		$temp1[1]=~s/ //g; 
		$temp1[1]=~s/\n//g;    
	     
		push (@primers, [@temp1]); 
		$list_len++;
	}

	my $index=$lncPrimerCount/50;
	my $nowCount=0;
	while($index>0)
	{
		$index--;
		$nowCount++;
	}
	if($index<0)
	{
		$nowCount--;
	}
	for (my $i = 0; $i <=0; $i++) {   
			for (my $j =$i+1; $j <= $list_len-1; $j++) {
				cal_dimer($primers[$i][0], $primers[$i][1], $primers[$j][0], $primers[$j][1],$delta,$whetherpic,$pdm_result);
				}	
		}
	$lncPrimerCount++;
	@temp=();
	@primers=();
	close MRNA_FILE;	
  }
close PRIMER_FILE;
}  




if($para[1]==2)
{
        push (@res,"Query\t\t Length_Query\t\t Target\t\t Length_Target\t\t dG\t\t ndG\t\t Start_Position_Query\t\t End_Position_Query\t\t Start_Position_Target\t\t End_Position_Target\n");
	      open(PDM,">>$pdm_result") or die "Can't open file: $!";
			  print PDM @res;	
			  close (PDM);
open(PRIMER_FILE, "<$lncprimerfile") or die "Can't open file $lncprimerfile : No such file or directory!\n";;
while (<PRIMER_FILE>) {	
	$list_len=0;
        my $line = $_;
	      my @temp = split (/[\t ]/, $line);
	      	if($temp[1] eq '' ||$temp[1] eq '' )
  {
  	print "the format of file (-f) is error,Please check.\n ";
  	exit;
  }  
        $temp[1]=~s///g; 
        $temp[1]=~s/ //g; 
        $temp[1]=~s/u/t/g;
		    $temp[1]=~s/U/T/g;
        $temp[3]=~s///g;
        $temp[3]=~s/ //g;
        $temp[3]=~s/u/t/g;
		    $temp[3]=~s/U/T/g;
        $temp[1]=~s/\n//g;
        $temp[3]=~s/\n//g;                       
	      
	cal_dimer_2($temp[0],$temp[1],$temp[2],$temp[3],$delta,$whetherpic,$pdm_result);		
	$lncPrimerCount++;
	@temp=();
	@primers=();

  }
close PRIMER_FILE;
}
close LOGF;
