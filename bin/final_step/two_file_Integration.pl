use strict;
use warnings;
use Cwd qw(abs_path)
  ; #获取工作路径，即当前目标所在的路径（函数括号里的东西）
use Getopt::Long; #获取选项
use Data::Dumper; #可打印引用的东西。eg：print Dumper(\%hash \@array);
use FindBin qw($Bin $Script)
  ; #$Bin  调用脚本的bin目录的路径，$Script  脚本名称  $RealBin 调用脚本的绝对路径  $RealScript  与脚本相关的脚本（用不着）
use File::Basename qw(basename dirname)
  ; #basename函数获取文件名  dirname函数获取路径  fileparse函数获取扩展名
my $BEGIN_TIME = time()
  ; #获取系统时间，可两者做差得到运行时间，单位为秒（s）
my $version = "1.0.0";

#---------------------------------------------------------------------------
#GetOptions
#---------------------------------------------------------------------------
my ( $input_file1, $input_file2, $output_file, $input1_ncol, $input2_ncol );
GetOptions(
    "in1=s" => \$input_file1,    #=s表示串，eg：=i表示整数
    "in1_col=i" =>
      \$input1_ncol,    #输入的是input_file1拥有相同的列名的列数
    "in2=s" => \$input_file2,
    "in2_col=i" =>
      \$input2_ncol,    #输入的是input_file2拥有相同的列名的列数
    "out=s"  => \$output_file,
    "help|?" => \&USAGE,
) or &USAGE;
&USAGE
  unless ( defined $input_file1
    and defined $input1_ncol
    and defined $input_file2
    and defined $input2_ncol
    and defined $output_file );    #判断选项的合法性
&log_current_time("$Script start……");    #调用时间函数

#----------------------------------------------------------------------------
#load input file and print the result
#----------------------------------------------------------------------------
open READ_1, "$input_file1" || die "can't open the file: $!\n";
open READ_2, "$input_file2" || die "can't open the file: $!\n";

#$output_file ||="./".basename($input_file).".out";
open WRITE, ">$output_file" || die "can't open the file :$_\n";

##先汇总标题，在file1中将相同的一列名称去掉,打印整合后的共同的列名
##整合后的列名顺序为file2-file1
my @array1 = split /\t/, <READ_1>;
chomp(@array1);
my $file1_ncol = scalar(@array1);
my @array2 = split /\t/, <READ_2>;
chomp(@array2);
my $file2_ncol = scalar(@array2);
splice( @array1, $input1_ncol, 1 );
my $tmp_str1 = join( "\t", ( @array2, @array1 ) );
print WRITE"$tmp_str1\n";

##读取文件，由于存在相同列中有重复值，存入hash的hash
#input_file1以行名为key1,以相同列为key2，以其他信息为value
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

##以file2为顺序，对相同的列进行匹配输出最后到整合的文件

while (<READ_2>) {
    chomp;
    my @array = split /\t/, $_;
    my $tmp_id = $array[$input2_ncol];
    my $tmp;
    my $flag = 0;    #用于记录是否有重复值
    foreach ( keys %geneid1 ) {
        if ( exists $geneid1{$_}{$tmp_id} ) {
            $tmp = $geneid1{$_}{$tmp_id}
              ;      #用于存储需要加到输出文件的file1的信息
            delete $geneid1{$_}{$tmp_id}
              ;      #将打印之后的file1中的hash内容删除
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
##将file1中未与file2匹配的内容打印出来
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
&log_current_time("$Script end……");    #调用时间函数
my $run_time = time() - $BEGIN_TIME;
print "$Script run time :$run_time\.s\n";

#----------------------------------------------------------------------------------
#function
#----------------------------------------------------------------------------------
sub log_current_time {    #获取时间的函数，格式化输出
    my ($info) =
      @_;    #获取参数（一般为    XXX程序开始执行…………）
    my $curr_time = &date_time_format( localtime( time() ) )
      ;                              #格式化获取时间表示方法
    print "[$curr_time] $info\n";    #输出打印
}
##########################################################################################################################
sub date_time_format {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );    #获取格式化的时间信息，不打印（sprintf的用法）
}
##########################################################################################################################
sub USAGE {    #选项帮助函数
    my $usage =
      <<"__USAGE__"; #从下一行开始，知道碰到__USAGE__为止，所有的东西按统一的格式存入变量中（包括注释），#__USAGE__为结束符号，类似用eof，不会输出也不会输入，仅仅只是作为结束标志
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
