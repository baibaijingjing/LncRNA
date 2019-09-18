use strict;
use warnings;
use Cwd qw(abs_path); #��ȡ����·��������ǰĿ�����ڵ�·��������������Ķ�����
use Getopt::Long; #��ȡѡ��
use Data::Dumper; #�ɴ�ӡ���õĶ�����eg��print Dumper(\%hash \@array);
use FindBin qw($Bin $Script); #$Bin  ���ýű���binĿ¼��·����$Script  �ű�����  $RealBin ���ýű��ľ���·��  $RealScript  ��ű���صĽű����ò��ţ�
use File::Basename qw(basename dirname); #basename������ȡ�ļ���  dirname������ȡ·��  fileparse������ȡ��չ��
my $BEGIN_TIME = time(); #��ȡϵͳʱ�䣬����������õ�����ʱ�䣬��λΪ�루s��
my $version = "1.0.0";
use newPerlBase;
my $Title="LncRNA";
my %config=%{readconf("$Bin/../../Config/lncRNA_pip.cfg")};

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GetOptions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my ( $cfg1, $cfg2, %option_for_flag, $species, $odir, %species_for_database )
  ;    #ѡ���hash��
GetOptions(
    "cfg1=s" => \$cfg1,    #=s��ʾ����eg��=i��ʾ����
    "cfg2=s" => \$cfg2,    #=s��ʾ����eg��=i��ʾ����
    "od:s"   => \$odir,
    "help|?" => \&USAGE,
) or &USAGE;
&USAGE unless ( defined $cfg1 and defined $cfg2 );    #�ж�ѡ��ĺϷ���
&log_current_time("$Script start����");
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
#system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset"; ##�ǵ�Ҫ����������ļ���·��
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
    #system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset"; ##�ǵ�Ҫ����������ļ���·��
    system "Rscript $Bin/id_conversion.r --infile $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls --outfile $odir/Anno_Integrate/Allgene_Anno/Result/trans_id_out --id_name ensembl_gene_id --id.col 1 --database $database --dataset $dataset --host $host"; ##biomaRt change 2015-11-10 by linhj
    #system "rm $odir/Anno_Integrate/Allgene_Anno/Result/tmp.xls";
    }
=cut

###################################################################################################################################################################################
&log_current_time("$Script end����");    #����ʱ�亯��
my $run_time = time() - $BEGIN_TIME;
print "the program run time is:$run_time.s\n";

sub log_current_time {

    # get parameter    #��ȡ����
    my ($info) = @_;

    # get current time with string   #��ȡ��ǰ����ʱ��
    my $curr_time = &date_time_format( localtime( time() ) )
      ;                                    #��ʽ����ȡʱ���ʾ����

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
