use strict;
use warnings;
use Cwd qw(abs_path)
  ; #��ȡ����·��������ǰĿ�����ڵ�·��������������Ķ�����
use Getopt::Long; #��ȡѡ��
use Data::Dumper; #�ɴ�ӡ���õĶ�����eg��print Dumper(\%hash \@array);
use FindBin qw($Bin $Script)
  ; #$Bin  ���ýű���binĿ¼��·����$Script  �ű�����  $RealBin ���ýű��ľ���·��  $RealScript  ��ű���صĽű����ò��ţ�
use File::Basename qw(basename dirname)
  ; #basename������ȡ�ļ���  dirname������ȡ·��  fileparse������ȡ��չ��
my $BEGIN_TIME = time()
  ; #��ȡϵͳʱ�䣬����������õ�����ʱ�䣬��λΪ�루s��
my $version = "1.0.0";

#---------------------------------------------------------------------------
#GetOptions
#---------------------------------------------------------------------------
my ( $input_file1, $input_file2, $output_file, $input1_ncol, $input2_ncol );
GetOptions(
    "in1=s" => \$input_file1,    #=s��ʾ����eg��=i��ʾ����
    "in1_col=i" =>
      \$input1_ncol,    #�������input_file1ӵ����ͬ������������
    "in2=s" => \$input_file2,
    "in2_col=i" =>
      \$input2_ncol,    #�������input_file2ӵ����ͬ������������
    "out=s"  => \$output_file,
    "help|?" => \&USAGE,
) or &USAGE;
&USAGE
  unless ( defined $input_file1
    and defined $input1_ncol
    and defined $input_file2
    and defined $input2_ncol
    and defined $output_file );    #�ж�ѡ��ĺϷ���
&log_current_time("$Script start����");    #����ʱ�亯��

#----------------------------------------------------------------------------
#load input file and print the result
#----------------------------------------------------------------------------
open READ_1, "$input_file1" || die "can't open the file: $!\n";
open READ_2, "$input_file2" || die "can't open the file: $!\n";

#$output_file ||="./".basename($input_file).".out";
open WRITE, ">$output_file" || die "can't open the file :$_\n";

##�Ȼ��ܱ��⣬��file1�н���ͬ��һ������ȥ��,��ӡ���Ϻ�Ĺ�ͬ������
##���Ϻ������˳��Ϊfile2-file1
my @array1 = split /\t/, <READ_1>;
chomp(@array1);
my $file1_ncol = scalar(@array1);
my @array2 = split /\t/, <READ_2>;
chomp(@array2);
my $file2_ncol = scalar(@array2);
splice( @array1, $input1_ncol, 1 );
my $tmp_str1 = join( "\t", ( @array2, @array1 ) );
print WRITE"$tmp_str1\n";

##��ȡ�ļ������ڴ�����ͬ�������ظ�ֵ������hash��hash
#input_file1������Ϊkey1,����ͬ��Ϊkey2����������ϢΪvalue
my %geneid1;
my $col_1 = 1;
while (<READ_1>) {
    chomp;
    my @array = split /\t/, $_;
    my $tmp_id = splice( @array, $input1_ncol, 1 );
    my $tmp_file1 = join( "\t", @array );
    $geneid1{$col_1}{$tmp_id} = $tmp_file1;
    $col_1++;
}

##��file2Ϊ˳�򣬶���ͬ���н���ƥ�����������ϵ��ļ�

while (<READ_2>) {
    chomp;
    my @array = split /\t/, $_;
    my $tmp_id = $array[$input2_ncol];
    my $tmp;
    my $flag = 0;    #���ڼ�¼�Ƿ����ظ�ֵ
    foreach ( keys %geneid1 ) {
        if ( exists $geneid1{$_}{$tmp_id} ) {
            $tmp = $geneid1{$_}{$tmp_id}
              ;      #���ڴ洢��Ҫ�ӵ�����ļ���file1����Ϣ
            delete $geneid1{$_}{$tmp_id}
              ;      #����ӡ֮���file1�е�hash����ɾ��
            my $tmp_str = join( "\t", ( @array, $tmp ) );
            print WRITE"$tmp_str\n";
            $flag++;
            $tmp = '';
        }
    }
    if ( $flag == 0 ) {
        $tmp = "--\t";
        my $i;
        for ( $i = 2 ; $i < $file1_ncol ; $i++ ) {
            $tmp .= "--\t";
        }
        my $tmp_str = join( "\t", ( @array, $tmp ) );
        print WRITE"$tmp_str\n";
        $tmp = '';
    }
}
##��file1��δ��file2ƥ������ݴ�ӡ����
my ( $col, @tmp, $id );
foreach $col ( sort { $geneid1{$a} <=> $geneid1{$b} } keys %geneid1 ) {
    foreach $id ( sort keys %{ $geneid1{$col} } ) {
        my $i;
        splice(@tmp);
        for ( $i = 0 ; $i < $file2_ncol ; $i++ ) {
            $tmp[$i] .= "--";
        }
        $tmp[$input2_ncol] = $id;
        my $tmp_str = join( "\t", ( @tmp, $geneid1{$col}{$id} ) );
        print WRITE "$tmp_str\n";
    }
}
close READ_2;
close READ_1;
close WRITE;

#----------------------------------------------------------------------------------
#ending of work
#----------------------------------------------------------------------------------
&log_current_time("$Script end����");    #����ʱ�亯��
my $run_time = time() - $BEGIN_TIME;
print "$Script run time :$run_time\.s\n";

#----------------------------------------------------------------------------------
#function
#----------------------------------------------------------------------------------
sub log_current_time {    #��ȡʱ��ĺ�������ʽ�����
    my ($info) =
      @_;    #��ȡ������һ��Ϊ    XXX����ʼִ�С���������
    my $curr_time = &date_time_format( localtime( time() ) )
      ;                              #��ʽ����ȡʱ���ʾ����
    print "[$curr_time] $info\n";    #�����ӡ
}
##########################################################################################################################
sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );    #��ȡ��ʽ����ʱ����Ϣ������ӡ��sprintf���÷���
}
##########################################################################################################################
sub USAGE {    #ѡ���������
    my $usage =
      <<"__USAGE__"; #����һ�п�ʼ��֪������__USAGE__Ϊֹ�����еĶ�����ͳһ�ĸ�ʽ��������У�����ע�ͣ���#__USAGE__Ϊ�������ţ�������eof���������Ҳ�������룬����ֻ����Ϊ������־
#-----------------------------------------------------------
 Program:$Script
 Version:$version
 Contact:BMK
 Data:2015-08-12
Function:assembly sequence
   USAGE:
         --in1       <STR>   input file1    [Must]
         --in1_col   <STR>  the same col number in input file1  [Must]
	       --in2       <STR>   input file2  [Must]
	       --in2_col   <STR>  the same col number in input file2 [Must]
         --out       <STR>   output file ,the head order is the file2-file1   [Must]
         --help          show the docment and exit
 Example:
    perl $Script --in1 test1.xls --in1_col 1 --in2  test2.xls --in2_col --out out.xls
#---------------------------------------------------------
__USAGE__
    print $usage;
    exit;
}
