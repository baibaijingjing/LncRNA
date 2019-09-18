use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename;
use Cwd qw(abs_path getcwd);

my ($regularod,$medicalod);
GetOptions(
        "h|?"           =>\&USAGE,

	"regularod:s"	=>\$regularod,
	"medicalod:s"	=>\$medicalod,
)or &USAGE;
&USAGE unless ($regularod);

$medicalod ||= $regularod;
$regularod=abs_path($regularod);
$medicalod=abs_path($medicalod);



if($regularod ne $medicalod){
	`mkdir -p $medicalod/Analysis_Report/Web_Report`	if(!-d "$medicalod/Analysis_Report/Web_Report");
#	`mkdir -p $resultcloud`	if(-d $resultcloud);
	`cp -r $regularod/Analysis_Report/Web_Report/BMK* $medicalod/Analysis_Report/Web_Report`;
	`cp -r $regularod/Analysis_Report/Web_Report/HTML $medicalod/Analysis_Report/Web_Report`;
	`cp $regularod/Analysis_Report/Web_Report/configtest.xml $medicalod/Analysis_Report/Web_Report`;
        `cp -r $regularod/Analysis_Report/Web_Report/Template $medicalod/Analysis_Report/Web_Report`;
}
my $resultpath="$medicalod/Analysis_Report/Web_Report/BMK_medical";  #####只考虑返回老师的结果
#my $resultcloud="$medicalod/Web_Report";

`mkdir -p $resultpath`	if(!-d $resultpath);
my $id=1;
my $folder;
#####################################################################################
###################替换read覆盖深度图
`cp $medicalod/Basic_Analysis/Tophat_Cufflinks/Map_Stat/*.map.png $medicalod/Analysis_Report/Web_Report/BMK_1_rawData/BMK_2_Mapped_Statistics`;


#####################################################################################
######################## known and new lncRNA distinguish

if(-d "$medicalod/Known_lncRNA"){
	$folder=join("_",("BMK",$id,"Known_lncRNA"));
	`mkdir $resultpath/$folder`	if(!-d "$resultpath/$folder");
	$id ++;	
	`cp $medicalod/Known_lncRNA/* $resultpath/$folder`;
	`cp $Bin/bin/Readme_known.txt $resultpath/$folder/Readme.txt`;
}

#####################################################################################
########################extract cosmic/tf result table

opendir(DIR,"$medicalod/DEG_Analysis")||die $!;
my @vs=grep {/vs/ && -d "$medicalod/DEG_Analysis/$_"} readdir(DIR);
closedir(DIR);
$folder=join("_",("BMK",$id,"DEG_anno"));
`mkdir $resultpath/$folder`	if(!-d "$resultpath/$folder");
$id ++;


foreach my $v(@vs){	
	`cp $medicalod/DEG_Analysis/$v/$v.DEG.cosmic.txt $resultpath/$folder/$v.DEG.cosmic.xls`;
	`cp $medicalod/DEG_Analysis/$v/$v.DEG.TF.txt $resultpath/$folder/$v.DEG.TF.xls`;
}
`cp $Bin/bin/Readme_anno.txt $resultpath/$folder/Readme.txt`;
#####################################################################################
##########################extract fusion gene result

$folder=join("_",("BMK",$id,"Gene_fusion"));
`mkdir  $resultpath/$folder`	if(!-d "$resultpath/$folder");
$id ++;
`cp -r $medicalod/Gene_Fusion/result/* $resultpath/$folder`;
`cp $Bin/bin/Readme_fusion.txt $resultpath/$folder/Readme.txt`;
#####################################################################################
########################sub function
sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage: 
        Options:
	-regularod	<path>	regular analysis outpath, forced
        -medicalod      <path>	medical analysis outpath, default regularod

        -h      Help

Example:

USAGE
        print $usage;
        exit;
}

