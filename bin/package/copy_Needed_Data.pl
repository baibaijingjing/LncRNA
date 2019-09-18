#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";


#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$od,$cfg);

GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$in,
                "od:s"=>\$od,
                "cfg:s"=>\$cfg,
                ) or &USAGE;
&USAGE unless ($in);
# ------------------------------------------------------------------
#
# ---indir is program folder---------------------------------------------------------------
$in = ABSOLUTE_DIR("$in");
mkdir $od unless (-d $od);
$od=ABSOLUTE_DIR("$od");
$cfg=ABSOLUTE_DIR("$cfg");
mkdir "$od/Tophat" unless (-d "$od/Tophat");

system "cp -r $in/Config $od/ ";
system "cp -r $in/Anno_Integrate $od/";
my @Sam=(glob "$in/Basic_Analysis/Tophat_Cufflinks/Tophat/*");
foreach my $dir (@Sam){
    my $sam=basename($dir);
    &MKDIR("$od/Tophat/$sam");
    system "cp $in/Basic_Analysis/Tophat_Cufflinks/Tophat/$sam/*.bam $od/Tophat/$sam/";
    system "cp $in/Basic_Analysis/Tophat_Cufflinks/Tophat/$sam/*.bam.bai $od/Tophat/$sam/";
    system "cp $in/Basic_Analysis/Tophat_Cufflinks/Tophat/$sam/*.txt $od/Tophat/$sam/";
}

system "cp -r $in/Basic_Analysis/Tophat_Cufflinks/Ref_Genome $od/";
system "cp -r $in/Basic_Analysis/geneExpression/final_track/All.longest_transcript.fa $od/";
system "perl $Bin/bin/final_step/get_sample.pl --cfg1 $cfg --cfg2 $Bin/bin/final_step/local_and_biomaRt_database.txt --od $in";
unless ( -f "$in/Analysis_Report/Web_Report/geneExpression/final_all_out.xls" ) {
    system "cp $in/Anno_Integrate/Allgene_Anno/Result/tmp.xls $in/Anno_Integrate/Allgene_Anno/final_all_out.xls";
    system "cp $in/Anno_Integrate/Allgene_Anno/Result/tmp.xls $in/Analysis_Report/Web_Report/geneExpression/final_all_out.xls";
}
else {
    system "cp $in/Analysis_Report/Web_Report/geneExpression/final_all_out.xls $in/Anno_Integrate/Allgene_Anno";
}
system "perl $Bin/bin/final_step/final_DEG_annotation.pl --in1 $in/DEG_Analysis/All_DEG/All.DEG_final.xls --in2 $in/Analysis_Report/Web_Report/geneExpression/final_all_out.xls --out $in/Analysis_Report/Web_Report/geneExpression/final_DEG_annotation.xls";    ################### only sample & annotation
system "cp $in/Analysis_Report/Web_Report/geneExpression/final_DEG_annotation.xls $in/DEG_Analysis/All_DEG";
system "cp -r $in/Analysis_Report/Web_Report $in";




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub MKDIR {    # &MKDIR($out_dir);
    my ($dir) = @_;

    #	rmdir($dir) if(-d $dir);
    mkdir($dir) if ( !-d $dir );
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

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
 Description:   This program is used to Copy Data  to Need_Data......
       Usage:
        Options:
        -in <dir>   Input directory,Analysis, forced.
        -od <dir>   Input directory,Analysis/Needed_Data, forced.
        -cfg <file> detail.cfg,forced.
        -h          help

USAGE
    print $usage;
    exit;
}
