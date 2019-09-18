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
my ($cfg1, $cfg2, $od, $step, $log, $oneStepOnly,$test);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg1:s"  =>\$cfg1,
				"cfg2:s"  =>\$cfg2,
				"od:s"   =>\$od,
				"s:s"    =>\$step,
				"test"	 =>\$test,
				"oneStepOnly:s"=>\$oneStepOnly,
				) or &USAGE;
&USAGE unless ($cfg1 and $cfg2 and $od) ;
################################

my $notename = `hostname`; chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

$step = $step || 3;

#
# log file 
#
my $startTime = GetTime();
my $user = `whoami`;  chomp $user;
my $workDir = `pwd`;  chomp $workDir;
my $task = "perl $Bin/$Script ".join("\t", @Original_ARGV);

open ($log,">", "$od/Tophat_cufflinks.".time().".log") or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $cfg1\n";
print $log "detail config file:  $cfg2\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins 
#==================================================================


#my $GFFREAD_BIN ="/share/nas2/genome/biosoft/cufflinks/2.2.1/gffread";    # 2014-12-17 ~
my $GFFREAD_BIN=$config{GFFREAD_BIN};
#my $TOPHAT_BIN = "/share/nas2/genome/biosoft/tophat/2.0.13_fix/tophat2";    # 2015-10-09 ~
my $TOPHAT_BIN=$config{TOPHAT_BIN};
#my $GTF2FA_BIN = "/share/nas2/genome/biosoft/tophat/2.0.7/gtf_to_fasta";
my $GTF2FA_BIN=$config{GTF2FA_BIN};
## cufflinks suite
#my $CUFFLINKS_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cufflinks";    # 2014-12-17 ~
my $CUFFLINKS_BIN=$config{CUFFLINKS_BIN};
#my $CUFFCOMPARE_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffcompare";    # 2014-12-17 ~
my $CUFFCOMPARE_BIN=$config{CUFFCOMPARE_BIN};
#my $CUFFMERGE_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffmerge";      # 2014-12-17 ~
my $CUFFMERGE_BIN=$config{CUFFMERGE_BIN};
#my $CUFFDIFF_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffdiff";       # 2014-12-17 ~
my $CUFFDIFF_BIN=$config{CUFFDIFF_BIN};

my $SCRIPTURE_BIN = "$Bin/software/scripture-beta2.jar";

#my $CUFFNORM_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffnorm";       # 2014-12-17 ~
my $CUFFNORM_BIN=$config{CUFFNORM_BIN};
#my $CUFFQUANT_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffquant";
my $CUFFQUANT_BIN=$config{CUFFQUANT_BIN};

#==================================================================
# load config file 
#==================================================================

my %total_read;
my %para;
my %sample;

open (IN,"cat $cfg1 $cfg2|") || die "$!\n";
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	
	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/Sample/) {
		my $fq1=<IN>;  chomp $fq1;
		my $fq2=<IN>;  chomp $fq2;
		my @fq_1=split /\s+/,$fq1;
		$sample{$tmp[1]}{FQ1}=$fq_1[1];
		if (!-f "$od/totalRead.stat.xls") {
			my $total_line=`less -S $sample{$tmp[1]}{FQ1} | wc -l `;chomp $total_line;
			$total_read{$tmp[1]}=$total_line/2;
		}
		my @fq_2=split /\s+/,$fq2;
		$sample{$tmp[1]}{FQ2}=$fq_2[1];
	}
	$para{$tmp[0]}=$tmp[1];
}
close IN;

if (!-f "$od/totalRead.stat.xls") {
	open OUT,">$od/totalRead.stat.xls" || die;
	foreach my $sam (sort keys %total_read) {
		print OUT "$sam\t$total_read{$sam}\n";
	}
	close OUT;
}
else {
	open IN,"$od/totalRead.stat.xls" || die;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @tmp=split/\s+/,$_;
		$total_read{$tmp[0]}=$tmp[1];
	}
}

#==================================================================
# pipeline 
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome;
my $genome_size;
my $gtf;
my $gff;
my @chromosome;
my $idx_prefix;

if ($step!=1) {
	&MKDIR("$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome=basename ($para{Ref_seq});
	$idx_prefix = basename($para{Ref_seq});
	$idx_prefix =~s/.fa$//;
	my $genome_path = dirname($para{Ref_seq});

	if (!-f "$od/Ref_Genome/$genome".".hdrs"){
	system "grep '>' $para{Ref_seq} > $od/Ref_Genome/$genome.hdrs";#2015/08/19,modified by niulg
	}
	$genome_size = "$od/Ref_Genome/genome_size.txt";
	my $gff_name=basename ($para{Ref_ann});
        $gff = $gff_name;
        $gff_name=~s/\.gff3?$//i;
        #system "ln -s $para{Ref_ann} ./";
        #system "$GFFREAD_BIN $gff -T -o $gff_name.gtf";
        $gtf = "$gff_name.gtf";
        chdir "../";

        #$step++ unless ($oneStepOnly) ;
        print STDOUT " check and build bowtie index  Done\n";
        print $log "check and build bowtie index  Done\n";
}
if ($step==1) {
	print STDOUT "=== check and build bowtie index ===\n";
	print $log  "=== check and build bowtie index ===\n";
	&MKDIR("$od/Ref_Genome");
	chdir "$od/Ref_Genome";
	$genome=basename ($para{Ref_seq});
	$idx_prefix = basename($para{Ref_seq});
    $idx_prefix =~s/.fa$//;
	if (!-f "$od/Ref_Genome/$genome".".hdrs"){
        system "grep '>' $para{Ref_seq} > $od/Ref_Genome/$genome.hdrs";#2015/08/19,modified by niulg
        }
	if (!-f "$idx_prefix.1.bt2") {
        my $genome_abs_path = $para{Ref_seq};
		my $genome_path = dirname($para{Ref_seq});

        if (-e "$genome_path/$idx_prefix.1.bt2") {
            system "ln -s $genome_abs_path $genome_path/$idx_prefix.*.bt2 ./";
		} elsif (-e "$genome_path/$idx_prefix.1.bt2l") {   # for large genomes more than about 4 billion nucleotides in length
            system "ln -s $genome_abs_path $genome_path/$idx_prefix.*.bt2l ./";
        } else {
            system "ln -s $genome_abs_path ./";
            system "$config{bowtie2build} $genome $idx_prefix";
        }
	}

	$genome_size = "$od/Ref_Genome/genome_size.txt";
	open OUT, ">$genome_size";
	my %genomes = fasta($genome);
	foreach my $chro (keys %genomes){
		next if length($chro)>11;##2015-11-12 for PSI genome by linhj

		my $genome_len = length($genomes{$chro});
		if ($genome_len>1000000){             ######2015.10.15 
			push @chromosome ,$chro;
			if (!-e "$od/Ref_Genome/$chro.fa" ){
				open OUT1, ">$od/Ref_Genome/$chro.fa";
				print OUT1 ">$chro\n$genomes{$chro}\n";
				close OUT1;
			}
            print OUT "$chro\t$genome_len\n";
		#push @chromosome ,$chro;
		#if (!-e "$od/Ref_Genome/$chro.fa" ){
			#open OUT1, ">$od/Ref_Genome/$chro.fa";
			#print OUT1 ">$chro\n$genomes{$chro}\n";
			#close OUT1;
		};
		#my $genome_len = length($genomes{$chro});
		#if ($genome_len>10000||length($chro)<6){
		#	print OUT "$chro\t$genome_len\n";
		#	}else{
		#		next;
		#	}
		};
		#print @chromosome ;
		#die;
	close OUT;
	
	################################## gff2gtf
	my $gff_name=basename ($para{Ref_ann});
	$gff = $para{Ref_ann};#$gff_name;
	$gff_name=~s/\.gff3?$//i;
	system "ln -s $para{Ref_ann} ./" unless (-f $gff);
	system "$GFFREAD_BIN $gff -T -o $gff_name.gtf";
	$gtf = "$gff_name.gtf";
	chdir "../";

	$step++ unless ($oneStepOnly) ;
	print STDOUT " check and build bowtie index  Done\n";
	print $log "check and build bowtie index  Done\n";

}

#################################### 
#
# step 2: Align the RNA-seq read to genome using Tophat2 + bowtie2
#
#########

if ($step==2) {
	print STDOUT "=== Align the RNA-seq read to genome using Tophat2 + bowtie2  ===\n";
	print STDOUT "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	print $log "=== Align the RNA-seq read to genome using Tophat2 + bowtie2  ===\n";
	print $log "Tophat mapping shell file: $od/work_sh/Tophat.sh\n";
	
	#
	# write shell
	#
	open OUT,">$od/work_sh/Tophat.sh" || die;
	&MKDIR("$od/Tophat");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Tophat/$sam");
		print OUT "cd $od/Tophat/$sam && ";
        print OUT "$TOPHAT_BIN --output-dir ./ --read-mismatches $para{Mismatch} --read-edit-dist $para{Mismatch} --max-intron-length 5000000 --library-type $para{Lib_type} --num-threads 8 --GTF $od/Ref_Genome/$gtf --mate-inner-dist $para{Insert_size} $od/Ref_Genome/$idx_prefix $sample{$sam}{FQ1} $sample{$sam}{FQ2} && ";
		print OUT "samtools index accepted_hits.bam accepted_hits.bam.bai &&\n";
	}
	close OUT;
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type}||="middle.q";
	$para{CPU} ||= "30";
	&Cut_shell_qsub("$od/work_sh/Tophat.sh",30, "50G", "$para{Queue_type}" );
	&Check_qsub_error("$od/work_sh/Tophat.sh");
	

	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

#####################################
#
# step 3: Scripture Assembly Analysis
#
########

if ($step==3) {
	print STDOUT "shell file: $od/work_sh/Cufflinks.sh\n";
	print $log "=== Cufflinks Assembly Analysis  ===\n";
	print $log "shell file: $od/work_sh/Cufflinks.sh\n";
	open OUT,">$od/work_sh/Cufflinks.sh" || die;
	&MKDIR ("$od/Cufflinks");
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cufflinks/$sam");
                print OUT "cd $od/Cufflinks/$sam && ";
                print OUT "$CUFFLINKS_BIN -o ./ -p 4 -g $od/Ref_Genome/$gtf --library-type $para{Lib_type} -u -L $sam $od/Tophat/$sam/accepted_hits.bam \n";
		
	}		
	close OUT;
	$para{Memory}||="15G";
	$para{Queue_type}||="general.q";
	$para{CPU} ||= "30";
	&Cut_shell_qsub("$od/work_sh/Cufflinks.sh",50,$para{'Memory'},"$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/Cufflinks.sh");
	
	$step++ unless ($oneStepOnly) ;
	print STDOUT "Cufflinks Finished \n";
	print $log "Cufflinks Finished\n";
}	
	


#################################### 
#
# step 4 & 5: Cuffcompare & Cuffmerge Analysis
#
########
##未对转录本进行二选其一过滤，样本多的时候不宜使用
if ($step==4) {
	print STDOUT "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";
	print STDOUT "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	print $log "=== Cuffcompare: compare gtf of each sample with known gene model ===\n";	
	print $log "cuffcompare shell file: $od/work_sh/Cuffcompare.sh\n";
	
	#
	# write shell 
	#
	open OUT1,">$od/work_sh/Cuffcompare.sh" || die;
	&MKDIR("$od/Cuffcompare");
	open OUT2,">$od/assembly_GTF_list.txt" || die;
	open OUT3,">$od/work_sh/Compare.sh" || die;
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffcompare/$sam");
		print OUT1 "cd $od/Cuffcompare/$sam && ";
		&MKDIR("$od/Compare");
		print OUT1 "$CUFFCOMPARE_BIN -r $od/Ref_Genome/$gtf $od/Cufflinks/$sam/transcripts.gtf \n";
		print OUT2 "$od/Cufflinks/$sam/transcripts.gtf\n";
	}
	print OUT3 "cd $od/Compare/ && ";
	print OUT3 "$CUFFCOMPARE_BIN -i $od/assembly_GTF_list.txt -o first -r $od/Ref_Genome/$gtf  && ";
	#print OUT3 "perl $Bin/bin/filter_transcript.2.pl -in1 $od/Compare/first.tracking -in2 $od/Compare/first.combined.gtf -out $od/Compare/middle.gtf &&";
	#print OUT3 "$CUFFCOMPARE_BIN -r $od/Ref_Genome/$gtf -o $para{Project_key} -s $od/Ref_Genome/$genome middle.gtf && ";
	print OUT3 "cp first.combined.gtf Compare.gtf";
	close OUT1;
	close OUT2;
	close OUT3;
	
	#
	# qsub
	#
	$para{Memory}||="15G";
	$para{Queue_type}||="general.q";
	$para{CPU} ||= "30";
	&Cut_shell_qsub("$od/work_sh/Cuffcompare.sh",$para{CPU},"30G",$para{Queue_type});
	&Check_qsub_error("$od/work_sh/Cuffcompare.sh");
	&Cut_shell_qsub("$od/work_sh/Compare.sh","$para{CPU}",$para{Memory},"$para{Queue_type}");
	&Check_qsub_error("$od/work_sh/Compare.sh");
	$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}


##################################### 
#
# step 5: Cuffquant and Cuffnorm Analysis
#
#########

if ($step==5) {
	print STDOUT "=== Cuffquant and Cuffnorm analysis(only use gene FPKM result) ===\n";
	print STDOUT "Cuffquant shell file: $od/work_sh/Cuffquant.sh\n";
	print STDOUT "Cuffnorm shell file: $od/work_sh/Cuffnorm.sh\n";
	#print STDOUT "bam2sam shell file : $od/work_sh/bam2sam.sh\n";
	print $log "=== Cuffquant and Cuffnormanalysis(only use gene FPKM result) ===\n";
	print $log "Cuffquant shell file: $od/work_sh/Cuffquant.sh\n";
	print STDOUT "Cuffnorm shell file: $od/work_sh/Cuffnorm.sh\n";
	#print $log "bam2sam shell file : $od/work_sh/bam2sam.sh\n";
	
	#
	# write shell
	# 
	open OUT1,">$od/work_sh/Cuffquant.sh" || die;
	open OUT2,">$od/work_sh/Cuffnorm.sh" || die;
	#open OUT2,">$od/work_sh/bam2sam.sh" || die;
	&MKDIR("$od/Cuffquant");
	&MKDIR("$od/Cuffnorm");
	print OUT2 "$CUFFNORM_BIN --library-type $para{Lib_type} -o $od/Cuffnorm -p 6 -total-hits-norm -no-update-check -output-format cuffdiff $od/Compare/Compare.gtf ";
        print OUT2 "-L ";
        my $lable;
        my $sam_str;
	foreach my $sam (sort keys %sample) {
		&MKDIR("$od/Cuffquant/$sam");
                print OUT1 "cd $od/Cuffquant/$sam && ";
		print OUT1 "$CUFFQUANT_BIN -o ./ -p 6 -b $od/Ref_Genome/$genome --library-type $para{Lib_type}  $od/Compare/Compare.gtf $od/Tophat/$sam/accepted_hits.bam \n";
		$lable.="$sam".",";
                $sam_str.="$od/Cuffquant/$sam/abundances.cxb ";
	}
	$lable=~s/,$//;
        $sam_str=~s/\s+$//;
        print OUT2 "$lable $sam_str\n";
	close OUT1;
	close OUT2;
	
	#
	# qsub 
	#
	$para{Memory}||="15G";
	$para{Queue_type}||="general.q";
	$para{CPU} ||= "30";
	&Cut_shell_qsub("$od/work_sh/Cuffquant.sh",$para{CPU},"30G",$para{Queue_type});
	&Check_qsub_error("$od/work_sh/Cuffquant.sh");
	&Cut_shell_qsub("$od/work_sh/Cuffnorm.sh",$para{CPU},"30G",$para{Queue_type});
	&Check_qsub_error("$od/work_sh/Cuffnorm.sh");
	
	#$step++ unless ($oneStepOnly) ;
	print STDOUT "\n";
	print $log "\n";
}

###############################
#
# step 6: Statistic bam files
#
# ######
if ($step==6) {
    print STDOUT "=== Statistic bam files  ===\n";
    print STDOUT "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
    print STDOUT "shell file: $od/work_sh/genome_bam2depth.sh\n";
    print STDOUT "shell file: $od/work_sh/genome_Checkgraph.sh\n";
    print $log "=== Statistic bam files  ===\n";
    print $log "shell file: $od/work_sh/Tophat_bam_stat.sh\n";
    print $log "shell file: $od/work_sh/genome_bam2depth.sh\n";
    print $log "shell file: $od/work_sh/genome_Checkgraph.sh\n";
	$para{Memory}||="15G";
	$para{Queue_type}||="general.q";
	$para{CPU} ||= "30";
	#
	# # write shell 
	#        
    open OUT1, ">$od/work_sh/Tophat_bam_stat.sh"   || die;
    open OUT2, ">$od/work_sh/genome_bam2depth.sh"  || die;
    open OUT3, ">$od/work_sh/genome_Checkgraph.sh" || die;
    open OUT4, ">$od/work_sh/plot_ReadDensity.sh"  || die;
	&MKDIR("$od/Map_Stat");
	my @str_type_stat;
	foreach my $sam (sort keys %sample) {
		push @str_type_stat, "$od/Map_Stat/$sam.type.stat";
		print OUT1 "perl $Bin/bin/bam2map_stat.2.pl -i $sam -bam $od/Tophat/$sam/accepted_hits.bam -totalRead $total_read{$sam} -od $od/Map_Stat\n";#insert 
		print OUT2 "$config{samtools} depth $od/Tophat/$sam/accepted_hits.bam >$od/Tophat/$sam/$sam.sort.bam.depth &&\n";
		print OUT3 "perl $Bin/bin/get_percent_of_exon_intro_inter_by_gff_v1.4.pl -gff $para{Ref_ann} -i $od/Tophat/$sam/$sam.sort.bam.depth -od $od/Map_Stat -index $sam  &&  ";#randcheck
		print OUT3 "perl $Bin/bin/draw_total_random.pl -id $od/Map_Stat -od $od/Map_Stat  \n";
		print OUT4 "perl $Bin/bin/plotReadDensity.2.pl -vf $para{Memory} -q $para{Queue_type} -cpu $para{CPU} -a $od/Ref_Genome/$genome.fai -f bam -i $od/Tophat/$sam/accepted_hits.bam -o $od/Map_Stat/ -k $sam &&\n";#2015/08/10,modify by niulg,map_stat
		print OUT4 "$config{Rscript} $Bin/bin/pie.R infile=$od/Map_Stat/$sam.type.stat outfile=$od/Map_Stat/$sam.type.png legend.col=1  value.col=2 skip=1 sep=t &&\n";#2015/08/10,modify bu niulg,type
	}
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;
	# #qsub 
	#`sh $od/work_sh/Tophat_bam_stat.sh`;
	
	&Cut_shell_qsub("$od/work_sh/Tophat_bam_stat.sh",  "$para{CPU}","$para{Memory}","$para{Queue_type}");
	&Cut_shell_qsub("$od/work_sh/genome_bam2depth.sh",  "$para{CPU}","$para{Memory}","$para{Queue_type}");
	&Cut_shell_qsub("$od/work_sh/genome_Checkgraph.sh", "$para{CPU}","$para{Memory}","$para{Queue_type}");
	&Cut_shell_qsub("$od/work_sh/plot_ReadDensity.sh",  "$para{CPU}","$para{Memory}","$para{Queue_type}");
	#&Check_qsub_error("$od/work_sh/Tophat_bam_stat.sh");

	my $str_type_stat = join " ", @str_type_stat;#2015/08/10,modify bu niulg
	`perl $Bin/bin/map_total_type.pl -i $str_type_stat -o $od/Map_Stat/Total_Type.png`;#2015/08/10,modify bu niulg
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
sub fasta{
	my $fa = shift;
	my %Fasta;
	open IN, $fa || die;
	$/='>';
	<IN>;
	while (<IN>) {
    	chomp;
   		my ($id,$seq)=split /\n+/,$_,2;
    	my $seq_id=(split /\s+/,$id)[0];
    	#$seq=~s/\s+//g;
		$Fasta{$seq_id} = $seq;
	}
	close IN ;
	return %Fasta;

}
sub LOAD_PARA {
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

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	chomp $line;
	if ($line<=1000) {
        if ($notename=~/cluster/) { ####2015-09-25
		#if ($notename=~/login\-0\-4/) {
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
            		if ($notename=~/cluster/) { ####2015-09-25
			#if ($notename=~/login\-0\-4/) {
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

sub USAGE {#
	my $usage=<<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination：Designed for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            #1  Check the index of Genome & build index for alignment
                            #2  run Tophat analysis
                            3  Cufflinks Assembly Analysis
                            4  run Cuffcompare analysis Alt splice
                            5  run Cuffquant and Cuffnorm
                            #6  Map Stat
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
	print $usage;
	exit;
}
