#!/usr/local/bin/perl -w

my $ver="1.0.2";
use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
######################.........................................

my %opts;
GetOptions(\%opts,"miR=s","lnc=s","cfg=s","od=s","key=s","type=s","step=s","h" );

if(!defined($opts{miR}) || !defined($opts{lnc}) || !defined($opts{cfg}) || !defined($opts{key}) || !defined($opts{type}) || !defined($opts{od}) || defined($opts{h})){
	&help();
	exit;
}

#==================================================================
# bin
#==================================================================
my $RNAHybrid_BIN = "/share/nas2/genome/biosoft/RNAhybrid/current/bin/RNAhybrid";
my $miranda_BIN = "/share/nas2/genome/biosoft/miranda/bin/miranda";

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################
my $mirna=&ABSOLUTE_DIR($opts{miR});
my $gene=&ABSOLUTE_DIR($opts{lnc});
my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $type=$opts{type};
my $key=$opts{key};
&MKDIR($opts{od});
my $odir=&ABSOLUTE_DIR($opts{od});
my $work_sh="$odir/work_sh";
my $miR_div="$odir/miR_div";
&MKDIR($work_sh);
&MKDIR($miR_div);
`ln -s $mirna $odir`;
my $step=$opts{step} || 0;

print "start from step $step\n";

###########
my %para;
&LOAD_PARA($cfg,\%para);

## if the user forgot setting software options, choose RNAhybrid as default analysis tool
$para{RNAhybrid} = 1 unless (defined $para{RNAhybrid}) ;

## if the user forget setting parameter for miRanda, using the following default values
$para{miRanda_SCORE} = 50.0 unless (defined $para{miRanda_SCORE}) ;
$para{miRanda_EN} = -20 unless (defined $para{miRanda_EN}) ;
$para{miRanda_scale} = 4.0 unless (defined $para{miRanda_scale}) ;
$para{miRanda_go} = -2 unless (defined $para{miRanda_go}) ;
$para{miRanda_ge} = -8 unless (defined $para{miRanda_ge}) ;


my %seq;
&Cut_mir($mirna,$miR_div,$key,\%seq);

my $gene_utr="$odir/mRNA_utr.fa";
if ($type==0) {
	&Gene_utr($gene,$gene_utr);
}

############Target Predict
if ($step==0) {
	print "Target prediction start from the beginning...\n";
	if ($type==0) {
		print "Choose animal Predict Use RNAhybrid or miranda or both\n\n";
		my @files=glob "$miR_div/*.fa";
		open SH,">$work_sh/$key.TargetFinder.sh" || die "$!";
		if (defined $para{RNAhybrid} && $para{RNAhybrid} == 1) {
			foreach my $file (@files) {
				my $name=basename($file);
				print SH "cd $miR_div && ";
				print SH "$RNAHybrid_BIN -d 1.9,0.28 -b 1 -e $para{RNAhybrid_EN} -t $gene_utr -q $file >$name.RNAhybrid.txt && ";
				print SH "perl $Bin/bin/RNAhybrid_result.pl $miR_div/$name.RNAhybrid.txt $miR_div $name &&";
				print SH "perl $Bin/bin/testsvm.pl $name.RNAhybrid.aln.txt $name.RNAhybrid.aln_svm.txt \n";
			}
		}
		
		if (defined $para{miRanda} && $para{miRanda} == 1) {
			foreach my $file (@files) {
				my $__miranda_out__ = "$file.miRanda.txt";
				my $name = basename $file;
				my $miranda_dir = dirname $file;
				
				my $cmd = "$miranda_BIN $file $gene -sc $para{miRanda_SCORE} -en $para{miRanda_EN} -scale $para{miRanda_scale} -go $para{miRanda_go} -ge $para{miRanda_ge} -out $__miranda_out__ -quiet ";
				$cmd .= " -strict " if (defined $para{miRanda_strict} && $para{miRanda_strict}) ;
				$cmd .= "&& ";
				$cmd .= "perl $Bin/bin/miRanda_result.pl $__miranda_out__ $miranda_dir $name.miRanda.aln ";

				print SH $cmd, "\n";
			}
		}

		close SH;
		print "Begin to Target prediction qsub ...\n\n";
		#if(`hostname` =~/cluster/){
			#`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc 20 --independent $work_sh/$key.TargetFinder.sh `;
			#&cmd_call("sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $work_sh/$key.TargetFinder.sh --reqsub -maxproc 20 --independent");
		#	&cmd_call("sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc 20 --independent $work_sh/$key.TargetFinder.sh");
		#	$step=1;
		#}else{
			#`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc 20 --independent $work_sh/$key.TargetFinder.sh `;
			#&cmd_call("sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $work_sh/$key.TargetFinder.sh --reqsub -maxproc 20 --independent");
		
		&Cut_shell_qsub("$work_sh/$key.TargetFinder.sh","$para{Cpu}","$para{Queue_type}");
		$step=1;
		#}
		print "Target prediction proccess done!\n\n";
	}
	if ($type==1) {
		print "Choose plant Predict Use TargetFinder\n\n";
		open SH,">$work_sh/$key.TargetFinder.sh" || die "$!";
		foreach my $id (keys %seq) {
			if ($para{TargetFinder_rev}==0) {
				print SH "perl $Bin/bin/targetfinder.pl -s $seq{$id} -d $gene -q $id -c $para{TargetFinder_score} >$miR_div/$id.TargetFinder.txt && ";
				print SH "perl $Bin/bin/targetFinder_result.pl $miR_div/$id.TargetFinder.txt $miR_div $id \n";
			}
			else {
				print SH "perl $Bin/bin/targetfinder.pl -s $seq{$id} -d $gene -q $id -c $para{TargetFinder_score} -r >$miR_div/$id.TargetFinder.txt &&";
				print SH "perl $Bin/bin/targetFinder_result.pl $miR_div/$id.TargetFinder.txt $miR_div $id \n";
			}
		}
		close SH;
		print "Begin to TargetFinder qsub ...\n\n";
		&Cut_shell_qsub("$work_sh/$key.TargetFinder.sh","$para{Cpu}","$para{Queue_type}");
		$step=1;
		print "TargetFinder proccess done!\n\n";
=cut
		if(`hostname` =~/cluster/){
			`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc 20 --independent $work_sh/$key.TargetFinder.sh `;
			$step=1;
		}else{
			`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc 20 --independent $work_sh/$key.TargetFinder.sh `;
			$step=1;
		}
=cut
		
	}
}

############ cat target prediction result

if ($step==1 && $type==0) {
	if (defined $para{RNAhybrid} && $para{RNAhybrid} == 1) {
		`cat $miR_div/*.RNAhybrid.aln_svm.txt >$odir/$key.RNAhybrid.aln_svm.txt`;
	}
	if (defined $para{miRanda} && $para{miRanda} == 1) {
		`cat $miR_div/*.miRanda.aln.txt > $odir/$key.miRanda.aln.txt`;
	}
	
	$step=2;
}

if ($step==1 && $type==1) {
	`cat $miR_div/*.TargetFinder.aln.txt >$odir/$key.TargetFinder.aln.txt`;
	$step=2;
}

`rm -r $miR_div/`;
############## target stat
my %mir_tar;
my %target_gene;
my $mir_num;
my $target_gene_num;
if ($step==2 && $type==0) {
	if (defined $para{RNAhybrid} && $para{RNAhybrid} == 1) {
		&_animal_target_stat_RNAhybrid("$odir/$key.RNAhybrid.aln_svm.txt",\%mir_tar,\%target_gene);
	}

	if (defined $para{miRanda} && $para{miRanda} == 1) {
		&_animal_target_stat_miRanda("$odir/$key.miRanda.aln.txt",\%mir_tar,\%target_gene);
	}

	my $target_list_file = "$odir/$key.mir2target.list";
	my $output_mode = ($para{RNAhybrid} && $para{miRanda}) ? 'both' : ($para{RNAhybrid} ? 'rnahybrid' : 'miranda');
	&_animal_output_taget_gene_list(\%mir_tar, \%target_gene, $target_list_file, $output_mode, \$mir_num, \$target_gene_num);
	
#	open OUT,">$odir/$key.mir2target.list\n" || die $!;
#	foreach my $mir (sort keys %mir_tar) {
#		print OUT "$mir\t";
#		print OUT join (";",keys %{$mir_tar{$mir}}),"\n";
#	}
#	close OUT;
}

if ($step==2 && $type==1) {
	&_plant_target_stat("$odir/$key.TargetFinder.aln.txt",\%mir_tar,\%target_gene);
	open OUT,">$odir/$key.mir2target.list\n" || die $!;
	print OUT "miRNA\tTarget\n";
	foreach my $mir (sort keys %mir_tar) {
		print OUT "$mir\t";
		print OUT join(";",keys %{$mir_tar{$mir}}),"\n";
	}
	close OUT;

	$mir_num=keys %mir_tar;
	$target_gene_num=keys %target_gene;
}

## mir num and target num 
open OUT,">$odir/$key.mir2target.stat\n" || die $!;
print OUT "miR_num\t$mir_num\n";
print OUT "target_gene_num\t$target_gene_num\n";
close OUT;

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);

###########subs
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	my $shell = shift;
	my $cpu = shift;
	# my $vf = shift;
	 my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		#print "$shell\n";
		&cmd_call( "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc $cpu --queue $queue --independent $shell");
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
		 my $count=0;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_ \n";
				$line_num++;
			}
			else{
				$line_num=1;
				$div_index++;
				print OUT "$_ \n";
				close OUT;
			} 
			
		}
		if ($line_num!=1) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			#my $dir=dirname($shell);
			
				# &cmd_call( "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --reqsub -maxproc 20 --independent");
				&cmd_call( "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub -maxproc $cpu --queue $queue --independent $div_file");
				#print "$div_file \n";

		}
	}
}




sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}

sub LOAD_PARA
{
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub Cut_mir 
{ # &cut_mir($mir,$cut_dir,$prefix,\%seq);
	my ($fa,$cut_dir,$prefix,$seq_info)=@_;

	open IN,"$fa" || die $!;
	$/='>';
	my $num=1;
	while (<IN>) {
		chomp;
		s/\s+$//;
		next if (/^$/) ;
		open OUT,">$cut_dir/$prefix.$num.fa" || die;
		print OUT ">$_\n";
		close OUT;
		my ($id,$seq)=split/\n+/,$_,2;
		$id=~s/\s+.*$//;
		$seq_info->{$id}=$seq;
		$num++;
	}
	close IN;
	$/="\n";
}

sub Gene_utr 
{ # for animal use RNAhybrid : the gene length must <=2000
  # &Gene_utr($gene,$gene_utr);
  my ($gene_fa,$utr_fa)=@_;

  open IN,"$gene_fa" || die $!;
  open OUT,">$utr_fa" || die $!;
  $/='>';
  while (<IN>) {
      chomp;
      next if (/^$/);
      my ($id,$seq)=split /\n+/,$_,2;
      $id=~s/\s.*$//;
      $seq=~s/\s+//g;
      my $length=length $seq;
      my $new_seq;
      if ($length>2000) {
        $new_seq=substr ($seq,0,2000);
        print OUT ">$id\n$new_seq\n";
      }
      else {
        print OUT ">$id\n$seq\n";
      }
  }
  close IN;
  close OUT;
  $/="\n";
}

sub _animal_target_stat_RNAhybrid
{
	my ($aln,$stat,$gene)=@_;

	open IN,"$aln" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		if (/^>/) {
			$_=~s/^>//;
			my @tmp=split/\s+/,$_;
			$stat->{$tmp[0]}->{'RNAhybrid'}{$tmp[2]}=1;
			$gene->{'RNAhybrid'}{$tmp[2]}=1;
		}
	}
	close IN;
}


sub _animal_target_stat_miRanda {#
	my ($aln,$stat,$gene)=@_;

	open IN,"$aln" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		if (/^>>/) {
			$_=~s/^>>//g;
			my @tmp=split/\s+/,$_;
			$stat->{$tmp[0]}->{'miRanda'}{$tmp[1]}=1;
			$gene->{'miRanda'}{$tmp[1]}=1;
		}
	}
	close IN;
}

sub _animal_output_taget_gene_list {#
	my ($h_mir, $h_gene, $f_target, $m, $n_mir, $n_tar_gene) = @_;
	
	open (OUT,">$f_target") or die $!;
	if ($m eq 'both') {
		my %tmp = ();
		foreach my $mir (sort keys %{$h_mir}) {
			if (defined $h_mir->{$mir}{'RNAhybrid'} && defined $h_mir->{$mir}{'miRanda'}) {
				my $target_gene_intersection = intersection([keys %{$h_mir->{$mir}{'RNAhybrid'}}],[keys %{$h_mir->{$mir}{'miRanda'}}]);
				if (scalar @$target_gene_intersection > 0) {
					print OUT "$mir\t";
					print OUT join(";", @$target_gene_intersection),"\n";

					$$n_mir++;
					map {$tmp{$_}=1} @$target_gene_intersection;
				}
			}
		}

		$$n_tar_gene = scalar keys %tmp; 
	}elsif($m eq 'rnahybrid'){
		foreach my $mir (sort keys %{$h_mir}) {
			print OUT "$mir\t";
			print OUT join(";", keys %{$h_mir->{$mir}{'RNAhybrid'}}),"\n";
			$$n_mir++;
		}
		$$n_tar_gene = scalar keys %{$h_gene->{'RNAhybrid'}};
	}elsif($m eq 'miranda'){
		foreach my $mir (sort keys %{$h_mir}) {
			print OUT "$mir\t";
			print OUT join(";", keys %{$h_mir->{$mir}{'miRanda'}}),"\n";
			$$n_mir++;
		}
		$$n_tar_gene = scalar keys %{$h_gene->{'miRanda'}};
	}else{
		die "unknow mode\n";
	}
	close (OUT) ;
}

sub intersection {#
   my ($A,$B)=@_;
   my %uniqA=map {$_,1} @{$A};
   my %uniqB=map {$_,1} @{$B};
   my %merge=();
   my %overlap=();
   foreach  (keys %uniqA,keys %uniqB) {
           $merge{$_}++ && $overlap{$_}++;
   }
   my @result = keys %overlap;
   return \@result;
}

sub _plant_target_stat
{
	my ($aln,$stat,$gene)=@_;

	open IN,"$aln" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		if (/^>/) {
			$_=~s/^>//;
			my @tmp=split/\s+/,$_;
			$stat->{$tmp[0]}->{$tmp[1]}=1;
			$gene->{$tmp[1]}=1;
		}
	}
	close IN;
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
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

sub help
{
	print <<"	Usage End.";
	Description: Predict miRNA Target proccess (for plant and animals);
		Version: $ver
		Writer by Mengf <mengf\@biomarker.com.cn>
	Usage:
		-miR              miRNA fa                                        must be given;

		-gene             gene sequence                                   must be given;

		-type             species type (0 animal;1 plant)                 must be given;

		-cfg              soft parameter to select Target                 must be given;

		-key              Sample name                                     must be given;

		-od               Out dir                                         must be given;

		-step             process step (0 or 1 default 0)                 choice;
		                  0  from RNAhybrid or TargetFinder
		                  1  from RNAhybrid or TargetFinder Result tackle
		                  2  stat target predicton

		-h                help document
	Usage End.
	exit;
}

