#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3"; 
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($lnc_fa,$cfg1,$cfg2,$od, $step, $log,$oneStepOnly,$step_by_step,$gtf);
GetOptions(
				"help|?" =>\&USAGE,
				#"fa:s"  =>\$fa,
				"lnc_fa:s"  =>\$lnc_fa,
				"od:s"   =>\$od,
				"cfg1:s"    =>\$cfg1,
				"cfg2:s"	=>\$cfg2,
				"step:s"=>\$step,
				"sbs:s"=>\$step_by_step,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ( $cfg1 and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;
$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
$lnc_fa = &ABSOLUTE_DIR($lnc_fa) if (defined $lnc_fa);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

my %step;
$step ||= ($step_by_step) ? '1' : join ',', ( 1 .. 3 );
&steps_process( $step, $step_by_step, \%step );

delete $step{3};
my %para;
open (IN,"$cfg2") || die "$!\n";
while (<IN>) {
        chomp;
        s/\r$//;s/^\s+//;s/\s+$//;
        next if (/^\#/ || /^$/);

        my @tmp=split /\s+/,$_;
        $para{$tmp[0]}=$tmp[1];
}
print Dumper(\%step);
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/personality_analysis.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

#print $log "data config file:  $fa\n";

print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";


#==================================================================
my $cmd;



if (exists $step{1}) {
	print STDOUT "===  circRNA prediction  ===\n";
	print STDOUT "circRNA prediction shell file: $od/work_sh/circRNA_prediction.sh\n";
	print $log "=== circRNA prediction ===\n";
	print $log "circRNA predictionshell file: $od/work_sh/circRNA_prediction.sh\n";
	open OUT,">$od/work_sh/circRNA_prediction.sh" || die;
	#my $cmd;
	$cmd="perl $Bin/cirRNA_combine_analysis/cirRNA_identify.pl -cfg1 $cfg1 -cfg2 $cfg2 -od $od";
	print OUT "$cmd ";
        close OUT;
        &Cut_shell_qsub("$od/work_sh/circRNA_prediction.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
        &Check_qsub_error("$od/work_sh/circRNA_prediction.sh");

        $step++ unless ($oneStepOnly);
        print STDOUT "\n";
        print $log "\n";
}
if (exists $step{2}) {
        print STDOUT "===  cirRNA target miRNA analysis   ===\n";
        print STDOUT "cirRNA target miRNA analysis $od/work_sh/cirRNA_target.sh\n";
        print $log "=== cirRNA target miRNA analysis===\n";
        open OUT,">$od/work_sh/cirRNA_target.sh" || die;
        #my $cmd;
        &MKDIR("$od/Target_miRNA");
	if (exists $para{'miRNA_fa'} and -e $para{'miRNA_fa'}) {
            print OUT "perl $Bin/lncmi_com_analysis/bin/target_prediction/mir2target.cut.pl -miR $para{'miRNA_fa'} -lnc $od/circRNA.fa  -type $para{'SPECIES_TYPE'} -cfg $cfg2 -key $para{'Project_key'} -od $od/Target_miRNA  \n";
    }
    if (!exists $para{'miRNA_fa'}) {
        print OUT "perl $Bin/lncmi_com_analysis/bin/target_prediction/select_fa.pl -fa $config{miRBase} -i $para{'Key'} -o $od/Target_miRNA/$para{'Key'}_sRNA.fa &&  ";
        print OUT "perl $Bin/lncmi_com_analysis/bin/target_prediction/mir2target.cut.pl -miR $od/Target_miRNA/$para{'Key'}_sRNA.fa -lnc $od/circRNA.fa -type $para{'SPECIES_TYPE'} -cfg $cfg2 -key $para{'Project_key'} -od $od/Target_miRNA  \n";
    }
	close OUT;
	&Cut_shell_qsub("$od/work_sh/cirRNA_target.sh",$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	&Check_qsub_error("$od/work_sh/cirRNA_target.sh");
	print STDOUT "\n";
	print $log "\n";
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log  "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close($log);
#==================================================================
# subs 
#==================================================================
sub steps_process {
    print "\n";
    &log_current_time("step check:");
    my ( $step_str, $step_by_step, $step_hash ) = @_;
    my %_step_ = (
        '1'  => 'CircRNA_Analysis',
        '2'  => 'CircRNA_target_miRNA',
    );

    if ($step_by_step) {
        print "step-by-step: ON\n";
        if ( $step_str =~ /^\d+$/ ) {
            for my $s ( $step_str .. 12 ) {
                $step_hash->{$s} = $_step_{$s};
            }
        }
        else {
            die "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
        }
    }
    else {
        print "step-by-step: OFF\n";
        for my $s ( split /,/, $step_str ) {
            if ( $s =~ /^\d+$/ ) {
                $step_hash->{$s} = $_step_{$s};
            }
            else {
                die "ERROR: illegal steps specified by --step.\n";
            }
        }
    }

    print "steps_to_run: "
      . (
        join ", ",
        ( map { sprintf "$_.$step_hash->{$_}" } sort keys %{$step_hash} )
      ) . ".\n";
    &log_current_time("step check done.\n");
}
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	chomp $line;
	if ($line<=1000) {
		if ($notename=~/login\-0\-4/) {
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=1) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/login\-0\-4/) {
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
}

sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
#	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub step_cmd_process {
    my ($cmd, $sh_name, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$sh_name";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("$sh_name start...");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("sh $sh_file > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("$sh_name done, escaped time: $escaped_time.");
    }
}
sub log_current_time {

    # get parameter
    my ($info) = @_;

    # get current time with string
    my $curr_time = date_time_format( localtime( time() ) );

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

sub USAGE {#
	my $usage=<<"USAGE";

Description:
	This program is a Procedure deal with lncRNA prediction 


Usage:
	-cfg1	data.cfg ,	must be given;
	-cfg2	detail.cfg ,	must be given;
	-od 	output dir ,	./Lnc_filter must be given;
	-step		
		1:cirRNA idendify;
		2:cirRNA target miRNA prediction ;
	-sbs
	
USAGE
	print $usage;
	exit;
}
