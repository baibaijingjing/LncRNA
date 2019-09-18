#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use File::Spec;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($in);

GetOptions(
				"help|?" =>\&USAGE,
				"in:s"=>\$in,
				) or &USAGE;
&USAGE unless ($in);
$in=&ABSOLUTE_DIR($in);



#$in = '/media/niulg/A834496E34494114/lncRNA/8月份计划/change_name/';
my @dir=glob("$in/BMK*");
foreach my $in (@dir){
    my $bmk=basename($in);
    my @cases;
    &READDIR($in,\@cases);
    foreach my $file (@cases){
        my $name=basename($file);
        my $readme="$Bin/PDF/$bmk/$name.pdf";
        if (-f $readme) {
            system "cp $readme $file/Readme.pdf";
        }   
    }
}



sub READDIR {#判断目录是否存在
    my($in,$cases)=@_;
    my @files;
    my $dh;
    push(@files, $in);
    while (@files) {
        if (-d $files[0]) {#若是目录执行以下操作
            opendir $dh, $files[0] or die $!;#打开目录句柄,若失败打印错误信息
            @_ = grep { /^[^\.]/ } readdir $dh;#过滤掉以"."和".."的文件,即UNIX下的隐藏文件
            foreach (@_) {
                my $dir=File::Spec->catfile ($files[0], $_);
                if (-d $dir) {
                    push(@files, File::Spec->catfile ($files[0], $_));#连接目录名和文件名形成一个完整的文件路径:
                    push(@$cases, File::Spec->catfile ($files[0], $_));
                }
            }
            closedir $dh;
        }
        shift @files;
    }
    
}


#################################################
sub LOAD_PARA {
    my $para_file = shift;
    my $para      = shift;

    my $error_status = 0;
    open IN, $para_file || die "fail open: $para_file";
    while (<IN>) {
        chomp;
        s/^\s+//;
        s/\s+$//;
        s/\r$//;
        next if ( /^$/ or /^\#/ );
        my ( $para_key, $para_value ) = split( /\s+/, $_ );
        $para->{$para_key} = $para_value;
        }
    close IN;
    die "\nExit due to error Undefined parameter\n" if ($error_status);
}
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

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version

Usage:copy the readme.pdf to BMK_*_
  Options:
  -in       Web_Report
  
  -h         Help

USAGE
	print $usage;
	exit;
}
