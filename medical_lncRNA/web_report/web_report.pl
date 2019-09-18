use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename;

use newPerlBase;

my ($regularod,$medicalod);
GetOptions(
        "h|?"           =>\&USAGE,

	"regularod:s"	=>\$regularod,
	"medicalod:s"	=>\$medicalod,
	"conf:s"	=>\$conf,
)or &USAGE;
&USAGE unless ($regularod);

#my %config=&readConfig($conf);
$medicalod ||=$regularod;
$medicalod =abs_path($medicalod);
$regularod =abs_path($regularod);

####################### first step :extract result

runOrDie("perl $Bin/extract_result.pl -regularod $regularod -medicalod $medicalod");


####################### build new xml


runOrDie("perl $Bin/lncRNA_medical_xml.pl -oldxml $medicalod/Analysis_Report/Web_Report/configtest.xml -o $medicalod/Analysis_Report/Web_Report/configtest_medical.xml");

####################### convent to html

runOrDie("python $Bin/xmltohtml/xml2HtmlConverter.py -i $medicalod/Analysis_Report/Web_Report/configtest_medical.xml -o  $medicalod/Analysis_Report/Web_Report -n index_medical");


#####################################################################################
########################sub function
=head
my $formatDateTime = sub{
          my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
          my $format_time = sprintf("%4d-%02d-%02d %02d:%02d:%02d",
                    $year+1900, $mon+1, $day, $hour, $min, $sec);
          return $format_time;
};

sub timeLog {
          my $detail = shift;
          # get current time with string
          my $curr_time = &$formatDateTime(localtime(time()));
          # print info with time
          writeLog("[$curr_time] $detail");
}
sub writeLog{
          my $detail=shift;
          if(defined $userlog and defined $syslog){
                    open OUT1,">>$userlog" or die("$userlog PATH ERROR!\n");
                    open OUT2,">>$syslog" or die("$syslog PATH ERROR!\n");
                    print OUT1 "$detail \n";
                    print OUT2 "$detail \n";
                    close OUT1;
                    close OUT2;
          }elsif(defined $userlog and !defined $syslog){
                    open OUT1,">>$userlog" or die("$userlog PATH ERROR!\n");
                    print OUT1 "$detail \n";
                    close OUT1;
          }else{
                    print "$detail \n";
          }
}
sub runOrDie
{
          my ($cmd) = @_ ;
          if($cmd!~/\s/){
                    $cmd=abs_path($cmd);
                    $cmd="sh $cmd"; # for shell file
          }
          &timeLog($cmd);                                                                                                                                                                                            
          my $begintime=time();
          my $flag = system($cmd);
          if ($flag != 0){
                    &timeLog("Error: command fail: $cmd");
                    exit(1);
          }
          infoTime($begintime,"command done!");
          return ;
}
=cut

sub USAGE{
        my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage: 
        Options:
	-regularod	<path>	regular analysis outpath, forced
        -medicalod      <path>	medical analysis outpath, default regularod
	-db		<str>	the ref genome, optional (GRCh37,GRCh38,GRCm38,mm9,Rnor_6.0)
				default GRCh38
				this paramter whether identify the known lncRNA

        -h      Help

Example:

USAGE
        print $usage;
        exit;
}

