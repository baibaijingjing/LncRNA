#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
#print "$Bin/../../../Config/lncRNA_pip.cfg" if (-e "$Bin/../../../Config/lncRNA_pip.cfg");
my %config=%{readconf("$Bin/../../../Config/lncRNA_pip.cfg")};

my $USAGE = "

#*******************************************************************************************
Options:
-cufflinks   <dir>         Cufflinks OUT dir                               [required]
-od          <dir>         AS output dir                                   [default: (./)]
-Ref_Genome  <dir>         genome dir,*.hdrs and *.gtf should include      [required]
-step        <int>         The step of program                             [default 1]
              1             Extract AS
              2             Stat AS Result


#******************************************************************************************

";

my $extract_as="$Bin/software/ASprofile.b-1.0.4/extract-as";
my $summary_as="$Bin/software/ASprofile.b-1.0.4/summarize_as.pl";
my $as_fpkm="$Bin/software/ASprofile.b-1.0.4/extract-as-fpkm";
my $collect_fpkm="$Bin/software/ASprofile.b-1.0.4/collect_fpkm.pl";

my ($cufflinks, $fout, $step,$Ref_Genome,$cfg);

GetOptions (  "cufflinks:s"  =>\$cufflinks, "Ref_Genome:s"   =>\$Ref_Genome,  "od:s"   =>\$fout,    "step:n" =>\$step, "cfg:s" =>\$cfg,) || die $USAGE;

die $USAGE unless ($cufflinks and $Ref_Genome);

$fout ||= "./";
mkdir $fout unless (-d $fout);
$fout       = abs_path($fout);
$Ref_Genome = abs_path($Ref_Genome);
$cufflinks  = abs_path($cufflinks);
$cfg = abs_path($cfg);
#my $genome =(glob("$Ref_Genome/*.genome.fa"))[0];
#`cd $Ref_Genome && grep '>' $genome > $genome.hdrs`;
my $hdrs=(glob("$Ref_Genome/*.hdrs"))[0];
my $gtf =(glob("$Ref_Genome/*.gtf"))[0];

################# bioinfor_pipeline.log ###################
&show_log2("step_5: Alternative splicing analysis start.");
###########################################################

$step ||= 1;
my %para=%{readconf($cfg)};
my $note_name = `hostname`;  chomp $note_name;
if ($note_name =~ /cluster/)
{
	print "This jobs cannot be down on cluster host!\n";
	exit;
}

my @cufflinks=glob("$cufflinks/*");


my $work_sh = "$fout/work_sh";       mkdir $work_sh unless (-d $work_sh);
my $Step_1_sh = "$work_sh/Step_1_Extract_AS.sh";
my $Step_2_sh = "$work_sh/Step_2_Extract_fpkm.sh";
my $Step_3_sh = "$work_sh/Step_3_fpkm_compare.sh";

### Step 1 : Extract Alt gene list from cuffcmp.tracking file

if ($step==1)
{
	open (OUT1, ">", $Step_1_sh) || die "Open $Step_1_sh failed!\n";
	open (OUT2, ">", $Step_2_sh) || die "Open $Step_2_sh failed!\n";
	foreach my $cufflink (sort @cufflinks) 
	{
		my $sam=basename$cufflink;
		next if ($sam!~/\w+/);
		mkdir "$fout/$sam" unless (-d "$fout/$sam");
		print OUT1 "cd $fout/$sam && perl $Bin/bin/changename.pl $cufflink/cuffcmp.transcripts.gtf.tmap $cufflink/transcripts.gtf $fout/$sam && ";
		print OUT1 "$extract_as transcripts.gtf $hdrs  >$sam.tmap.as && ";
		print OUT1 "perl $summary_as  transcripts.gtf $sam.tmap.as  -p $sam && \n";
		print OUT2 "cd $fout/$sam && $as_fpkm  transcripts.gtf $hdrs $sam.as.nr  -W 9 >$sam.fpkm \n";
	}
	close OUT1;
	close OUT2;
	&Cut_shell_qsub($Step_1_sh,$para{'CPU'},$para{'Memory'},$para{'Queue_type'});
	#&Check_qsub_error($Step_1_sh);
	print "\n\nStep 1 :";
	print "\n\tExtract_AS is Done!\n\n";
	print "sh $Step_2_sh ";
	`sh $Step_2_sh `;
	print "\n\tExtract_fpkm is Done!\n\n";
	$step++;
}

### Step 2 : AS stat and draw
my %stat;
my %sample;
if ($step==2) 
{

	foreach my $cufflink (sort @cufflinks) 
	{
		my $sam=basename$cufflink;
		next if ($sam!~/\w+/);
		$sample{$sam}=1;
		my $stats="$fout/$sam/$sam.fpkm";
		open (IN,"$stats")||die $!;
		while (<IN>) {
			chomp;
			next if ($_=~/event|_OFF/);
			my $type=(split(/\t+/,$_))[1];
			$type=~s/_ON|_OFF//g;
			$stat{$type}{$sam}++;
		}
		`mv $fout/$sam/$sam.fpkm $fout `;

	}
	
#	open (OUT, ">", $Step_3_sh) || die "Open $Step_3_sh failed!\n";
#	print OUT "cd $fout && perl $collect_fpkm ";
#	@samples=sort keys %sample;
#	print OUT join(",",@samples),"  -s W9.fpkm > All_sample.W9.diff-exp";
#	close OUT;
#	`sh $Step_3_sh `;
#	print "\n\nStep 3 :";
#	print "\n\t Sample Fpkm Compare is Done!\n\n";
	
	open (OUT, ">$fout/As_Result.stat") || die "Open $fout/As_Result.stat  failed!\n";
	print OUT "AS\t",join("\t",sort {$a cmp $b} keys %sample),"\n";
	foreach my $type (sort {$b cmp $a} keys %stat) {
		print OUT "$type";
		foreach  my $sam(sort {$a cmp $b} keys %sample) {
			print OUT  "\t$stat{$type}{$sam}";
		}
		print OUT "\n";
	}
	close OUT;
	my $cmd="$config{Rscript} $Bin/bin/plot_stat_AS_events.r -i $fout/As_Result.stat -o $fout/As_event_stat -b FALSE";
	print $cmd,"\n";
	system $cmd;

	print "\n\nStep 4 :";
	print "\n\t Stat AS Result is Done!\n\n";

}
##############################
&show_log2("step_5: Alternative splicing analysis finished.");
#close LOG;
##############################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
			system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
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
				system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
}
=c
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
=cut
#############################################################################################################

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log1("cmd")
sub show_log1()
{
    open LOG, ">$fout/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$fout/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################