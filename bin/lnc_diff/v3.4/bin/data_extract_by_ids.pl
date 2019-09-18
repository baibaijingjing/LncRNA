#!/usr/bin/perl -w
# 
# Copyright (c) BMK 2013
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2014/3/28.
# Modifier:       lium <lium@biomarker.com.cn>
# Last Modified:  2014/8/27.
# v1.1 vs v1.0:   1) add filter option to control data filter and extract
#
# ------------------------------------------------------------------
# add package dir
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# setting version number
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# load package
# ------------------------------------------------------------------
use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require bmkPerlBase;
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};


unshift(@INC, "$config{bmkperl}");

# ------------------------------------------------------------------
# new bmkPerlBase
# ------------------------------------------------------------------
our $bmk = bmkPerlBase->new();


# ------------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------------
# Time start
my $BEGIN=time();
$bmk->timeLog(-info=>"Start Time");

# pipeline
main_pipeline();

# Time end
$bmk->timeLog(-info=>"End Time");
$bmk->totalRunTime(-start_time=>"$BEGIN");



# ------------------------------------------------------------------
# the define of main_pipeline
# it contains main step in this program
# ------------------------------------------------------------------
sub main_pipeline{
	# init variable
	my %opts;
	my %RunArgv = ();

	# get option and check
	get_options(\%opts);

	# check input file
	para_check(\%opts, \%RunArgv);
	$bmk->timeLog(-info=>"para_check is done");

	# check exp and create awk exp
	get_id(\%opts, \%RunArgv);
	$bmk->timeLog(-info=>"get_id is done");

	# calling awk
	data_extract(\%opts, \%RunArgv);
	$bmk->timeLog(-info=>"data_extract is done");
}


# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
sub get_options {
	my ($opts) = @_;
	GetOptions($opts, "idfile=s", "destfile=s", "out=s", "idcol=i", "destcol=i", "filter", "h" ) or USAGE();

	# check
	USAGE() if(!defined($opts->{"idfile"}) || !defined($opts->{"destfile"}) || !defined($opts->{"out"}) 
		|| defined($opts->{"h"}));

	# set default value
	$opts->{"idcol"} = 1 if( !defined($opts->{"idcol"}) );
	$opts->{"destcol"} = 1 if( !defined($opts->{"destcol"}) );
	# check
	if( $opts->{"idcol"} < 1 ){
		$bmk->logAndDie(-info=>"idcol must >= 1");
	}
	if( $opts->{"destcol"} < 1 ){
		$bmk->logAndDie(-info=>"destcol must >= 1");
	}
}


# ------------------------------------------------------------------
# the regular expression create and check 
# ------------------------------------------------------------------
sub get_id {
	# get parameter
	my ($opts, $run_argv) = @_;

	# open file
	open(IN, "$run_argv->{idfile}") or $bmk->logAndDie(-info=>"open idfile ($run_argv->{idfile}) is failed");

	# iter
	my %id = ();
	my $idcol = $opts->{idcol} - 1;
	my $line = "";
	while($line=<IN>){
		# get id
		chomp $line;
		my @str = split(/\t/,$line);
		# check and store
		if(exists $id{$str[$idcol]}){
			$bmk->logAndDie(-info=>"repeat id ($str[$idcol])");
		}else{
			$id{$str[$idcol]} = 1;
		}
	}

	# close
	close(IN);

	# store
	$run_argv->{"id"} = \%id;
}



# ------------------------------------------------------------------
# check parameter 
# ------------------------------------------------------------------
sub para_check {
	# get parameter
	my ($opts, $run_argv) = @_;

	# check input file
	if(!-f $opts->{"idfile"}){
		$bmk->logAndDie(-info=>"Error: idfile \"$opts->{idfile}\" not exist");
	}
	if(!-f $opts->{"destfile"}){
		$bmk->logAndDie(-info=>"Error: destfile \"$opts->{destfile}\" not exist");
	}

	# get absolute file
	# in
	$run_argv->{"idfile"} = $bmk->absolutePath(-in=>$opts->{"idfile"});
	$run_argv->{"destfile"} = $bmk->absolutePath(-in=>$opts->{"destfile"});
	# out
	$bmk->runOrDie(-cmd=>"touch $opts->{out}");
	$run_argv->{"out"} = $bmk->absolutePath(-in=>$opts->{"out"});
	# filter
	$run_argv->{"filter"} = 0;
	if( defined($opts->{"filter"}) ) {
		$run_argv->{"filter"} = 1;
	}
}



# ------------------------------------------------------------------
# calling awk
# ------------------------------------------------------------------
sub data_extract {
	# get parameter
	my ($opts, $run_argv) = @_;

	# open file
	open(IN, "$run_argv->{destfile}") 
		or $bmk->logAndDie(-info=>"open destfile ($run_argv->{destfile}) is failed");
	open(OUT, ">$run_argv->{out}") 
		or $bmk->logAndDie(-info=>"open or create output file ($run_argv->{out}) is failed");

	# iter
	my $id = $run_argv->{id};
	my $line = "";
	my $destcol = $opts->{destcol} - 1;
	# iter each record
	# extract data
	if( $run_argv->{"filter"} == 0 ) {
		while($line=<IN>){
			chomp $line;
			my @str = split(/\t/,$line);
			# check and output
			if(exists $id->{$str[$destcol]}){
				print OUT "$line\n";
			}
		}
	# filter data
	} else {
		while($line=<IN>){
			chomp $line;
			my @str = split(/\t/,$line);
			# check and output
			if(!exists $id->{$str[$destcol]}){
				print OUT "$line\n";
			}
		}
	}

	# close
	close(OUT);
	close(IN);
}




# ------------------------------------------------------------------
# print usage information and exit
# ------------------------------------------------------------------
sub USAGE {
	my $usage=<<"USAGE";
Version: 	Version$version
Writer: 	lium <lium\@biomarker.com.cn>
Description: 	extracting specified columns from file, 
Usage:
  Example:1) perl data_extract_by_ids.pl -idfile id.txt -destfile dest.txt -out out.txt -idcol 2 -destcol 3
  		2) perl data_extract_by_ids.pl -idfile id.txt -destfile dest.txt -out out.txt -idcol 2 -filter
  Options:
  forced parameters:
  -idfile 	<str> 	the file name of id
  -destfile 	<str> 	the file name of dest file
  -out 		<str> 	the file name of output
  optional parameters:
  -idcol 	<int> 	the colomn of id in idfile [default 1]
  -destcol 	<int> 	the colomn of id in destfile [default 1]
  -filter 	<none> 	filtering id record from dest file
  -h 		<none> 	Help

USAGE
	print $usage;
	exit(1);
}










