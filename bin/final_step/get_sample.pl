use strict;
use warnings;
use Cwd qw(abs_path); #获取工作路径，即当前目标所在的路径（函数括号里的东西）
use Getopt::Long; #获取选项
use Data::Dumper; #可打印引用的东西。eg：print Dumper(\%hash \@array);
use FindBin qw($Bin $Script); #$Bin  调用脚本的bin目录的路径，$Script  脚本名称  $RealBin 调用脚本的绝对路径  $RealScript  与脚本相关的脚本（用不着）
use File::Basename qw(basename dirname); #basename函数获取文件名  dirname函数获取路径  fileparse函数获取扩展名
my $BEGIN_TIME = time(); #获取系统时间，可两者做差得到运行时间，单位为秒（s）
my $version = "1.0.0";
use newPerlBase;
my $Title="LncRNA";
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my ( $cfg1, $cfg2, %option_for_flag, $species, $odir, %species_for_database )
  ;    #选项的hash表
GetOptions(
    "cfg1=s" => \$cfg1,    #=s表示串，eg：=i表示整数
    "cfg2=s" => \$cfg2,    #=s表示串，eg：=i表示整数
    "od:s"   => \$odir,
    "help|?" => \&USAGE,
) or &USAGE;
&USAGE unless ( defined $cfg1 and defined $cfg2 );    #判断选项的合法性
&log_current_time("$Script start……");
open READ_1, "$cfg1" || die "can't open the file:$_\n";
while (<READ_1>) {
    chomp;
    next if (/#/);
    my @arr = split /\s+/, $_, 2;
    $option_for_flag{ $arr[0] } = $arr[1];
}
if ( $option_for_flag{Ref_seq} =~ /\/.+?\/.+?\/.+?\/.+?\/.+?\/(.+?)\// ) {
    $species = $1;
}
close READ_1;
open READ_2, "$cfg2" || die "can't open the file:$_\n";
while (<READ_2>) {
    chomp;
    my @arr = split /\s+/, $_;
    $species_for_database{ $arr[0] } = [ @arr[ 1 .. 3 ] ];
}
###############2015/11/18 by luml
unless ( ${ $species_for_database{$species} }[0] ) {
    print "$species is not exist in $cfg2!";
    system "perl $Bin/two_file_Integration.pl --in2 $odir/Basic_Analysis/geneExpression/AllSample.genes_expression.xls --in2_col 0 --in1 $odir/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls --in1_col 0 --out $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls";    ######## intermediate document to final_all_out.xls

    #	next;
}
else {
    my $database = ${ $species_for_database{$species} }[0];
    my $dataset  = ${ $species_for_database{$species} }[1];
    my $host     = ${ $species_for_database{$species} }[2]
      ;    ##biomaRt change 2015-11-10 by linhj

#print "$database  $dataset\n";
#system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset"; ##记得要改输入输出文件的路径
    system "$config{Rscript} $Bin/id_conversion.r --infile $odir/Basic_Analysis/geneExpression/AllSample.genes_expression.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset --host $host";    ##biomaRt change 2015-11-10 by linhj
    system "perl $Bin/two_file_Integration.pl -in1 $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out.xls -in1_col 0 -in2 $odir/Basic_Analysis/geneExpression/AllSample.genes_expression.xls -in2_col 0 -out $odir/Basic_Analysis/geneExpression/tmp2.xls";    ######## intermediate document to final_all_out.xls
    system "perl $Bin/two_file_Integration.pl -in1 $odir/Anno_Integrate/Allgene_Anno/Result/Integrated_Function.annotation.xls -in1_col 0 -in2 $odir/Basic_Analysis/geneExpression/tmp2.xls -in2_col 0 -out $odir/Analysis_Report/Web_Report/geneExpression/final_all_out.xls";    #######################include gene samble & name,annotation
}

=pod
unless (${$species_for_database{$species}}[0]) {### 2015-11-17
	print "$species is not exist in $cfg2!";
	next;
}else{
    my $database=${$species_for_database{$species}}[0];
    my $dataset=${$species_for_database{$species}}[1];
    my $host=${$species_for_database{$species}}[2];##biomaRt change 2015-11-10 by linhj
    #print "$database  $dataset\n";
    #system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset"; ##记得要改输入输出文件的路径
    system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset --host $host"; ##biomaRt change 2015-11-10 by linhj
    #system "rm $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls";
    }
=cut

###################################################################################################################################################################################
&log_current_time("$Script end……");    #调用时间函数
my $run_time = time() - $BEGIN_TIME;
print "the program run time is:$run_time.s\n";

sub log_current_time {

    # get parameter    #获取参数
    my ($info) = @_;

    # get current time with string   #获取当前串的时间
    my $curr_time = &date_time_format( localtime( time() ) )
      ;                                    #格式化获取时间表示方法

    # print info with time
    print "[$curr_time] $info\n";
}

sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

###################################################################################################################################################################################
sub USAGE {
    my $usage = <<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: huangxy <huangxy\@biomarker.com.cn>
      Date: 2015-08-13

     Usage:
            --cfg1      <FILE>   Must(detail.cfg)
            --cfg2      <FILE>   Must(local_and_biomaRt_database.txt)
   Example:
            perl $Script --cfg1 cfg.txt --cfg2 cfg2.txt

----------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
