use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../../../Config/lncRNA_pip.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
my ($ref,$circRNAFa,$od);
my $BEGIN_TIME=time();
my $bigWigSummary =$config{bigWigSummary};
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($bw, $gtf,$type,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"bw:s"  =>\$bw,
				"circ:s"  =>\$circ,
				"od:s"   =>\$od,
				) or &USAGE;
&USAGE unless ($bw and $circ and $od) ;
mkdir $od if(!-e $od);
&conservation($circ,"$od/conservation.mean.result","mean");
&conservation($circ,"$od/conservation.max.result","max");
&conservation($circ,"$od/conservation.std.result","std");
&conservation($circ,"$od/conservation.min.result","min");
system("paste $od/conservation.mean.result $od/conservation.max.result $od/conservation.std.result $od/conservation.min.result > $od/All_conservation.result");
system("less $od/All_conservation.result|cut -f6,7,14,21,28 > $od/score.txt");
`$config{Rscript} $Bin/R/fpkm_density_plot_func.r $od/score.txt $od`;
# print "################################################\n";
# print "Draw boxplot of circRNA scorephastcons\n";
# if (-e "$od/scores.txt") {
# 	`rm $od/scores.txt`;
# }
# my $cmd1=qx(echo "Type    Score" >>$od/scores.txt);
# my $cmd3 = qx(cut -f 2,7 $od/conservation.result >>$od/scores.txt);
# #my $cmd4= qx(grep -v \'intergenic region\' $od/scores.txt > $od/scores.txt1);
# `sed -i 's/intergenic region/intergenic/g' $od/scores.txt`;
# my $cmd5 = "$config{Rscript} $Bin/R/phastcons_cdf.r $od/scores.txt $od ";
# $cmd5 .= ">$od/plot_box.log 2>&1";
# `$cmd5`;

#####################################################
#####sub

sub cmd_call {
	my ($shell,$cpu,$vf,$queue)=@_;
	system "sh $config{qsub_sh} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent $shell";
}

sub conservation{
	my ($circ,$out,$type) = @_;
	open(IN, "$circ") or die("Could not open input circular result file!\n"); 
	open(OUTPUT, ">$out") or die("Could not open output file!\n");
	print OUTPUT "chr\tfeature\tstart\tend\tlength\tcircRNA_ID\t${type}_score\n";
	while(<IN>)
	{
		next if(/#/);
		my @line = split/\t/,$_;
		($chr,$feature,$start,$end) = ("chr".$line[1],$line[8],$line[2],$line[3]);
		my $len= $end -$start;
		$trans = $line[0];
		system("$bigWigSummary -type=$type $bw $chr $start $end 1 >$od/$$.$type.tmp 2>&1");
		open(TEMP, "$od/$$.$type.tmp") or die("Cannot open result temp file!\n");
		my $score = '';
		while(<TEMP>) {
			$score .= $_;
			chomp($score);
		}
		close(TEMP);
		if($score =~ /^no/) { 
			$score = 0;
		}
		print OUTPUT "$chr\t$feature\t$start\t$end\t$len\t$trans\t$score\n";
		system("rm $od/$$.$type.tmp");
	}
	close(IN);
	close(OUTPUT);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: Get PhastCons from BigWig file downloaded from UCSC 
Version: $version
Contact: Zhang Qiuxue <zhangqx\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Get PhastCons from BigWig file downloaded from UCSC

Usage:
	-bw				BigWig file 
	-circ				circular result file 
	-od				output file 
 
USAGE
	print $usage;
	exit;
}
