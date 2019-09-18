use 5.012;
use Getopt::Long;
my ($fq1, $fq2, $out, $gtf, $coverage, $coverage2, $rand_mode, $rand_mode2, $read_length, $seq_err, $insert_length, $ref_dir, $help, $if_chr);
Getopt::Long::GetOptions(
	'1=s'	=>	\$fq1,
	'2=s'	=>	\$fq2,
	'O=s'	=>	\$out,
	'G=s'	=>	\$gtf,
	'C=i'	=>	\$coverage,
	'LC=i'	=>	\$coverage2,
	'R=i'	=>	\$rand_mode,
	'LR=i'	=>	\$rand_mode2,
	'L=i'	=>	\$read_length,
	'E=i'	=>	\$seq_err,
	'I=i'	=>	\$insert_length,
	'D=s'	=>	\$ref_dir,
	'CHR1=i'	=>	\$if_chr,
	'H!'	=>	\$help
);

my $if_die;
#my ($fq1, $fq2, $gtf) = (">>./simulate_80_10X_80bp_200_1.fq",">>./simulate_80_10X_80bp_200_2.fq",">>./gtf_80_10X_80bp_200.out");
#my $coverage = 10;
#my $coverage2 = 100;
#my $rand_mode = 1;
#my $rand_mode2 = 2;
#my $read_length = 80;
my $cRNA_size = 0;
#my $if_PE = 2;
#my $seq_err = 1;
#my $insert_length = 200;
#my $ref_dir = "/panfs/home/zhao/gaoyuan/bwaphage/hg19/";

if(defined($help)){
	print "This is CIRI-simulator, a simulation tool for circRNAs. Welcome!\n\n";
	print "Written by Yuan Gao. Any questions please mail to gaoyuan06\@mails.ucas.ac.cn.\n\n";
	print "Arguments (all required):\n";
	print "\t-1\t\toutput simulated PE reads file 1 name\n";
	print "\t-2\t\toutput simulated PE reads file 2 name\n";
	print "\t-O\t\toutput simulated reads list name\n";
	print "\t-G\t\tinput gtf formatted annotation file name\n";
	print "\t-C\t\tset coverage or max coverage (when choosing -R 2) for circRNAs\n";
	print "\t-LC\t\tset coverage or max coverage (when choosing -LR 2) for linear transcripts\n";
	print "\t-R\t\tset random mode for circRNAs: 1 for constant coverage; 2 for random coverage\n";
	print "\t-LR\t\tset random mode for linear transcripts: 1 for constant coverage; 2 for random coverage\n";
	print "\t-L\t\tread length of simulated reads (e.g. 100)\n";
	print "\t-E\t\tpercentage of sequencing error (e.g. 2)\n";
	print "\t-I\t\tinsertion length (should be larger than read length) (e.g. 350)\n";
	print "\t-D\t\tdirectory of reference sequence(s) (please make sure all references referred in gtf file are included in the directory)\n";
	print "\t-CHR1\t\tif only choose chr1 to simulate sequencing reads: 1 for yes; 0 for no\n";
	print "\t-H\t\tshow help information\n";
	$if_die = 1;
}elsif(!defined($fq1) or !defined($fq2) or !defined($out) or !defined($gtf) or!defined($coverage) or !defined($coverage2) or !defined($rand_mode) or !defined($rand_mode2) or !defined($read_length) or !defined($seq_err) or !defined($insert_length) or !defined($ref_dir)or !defined($if_chr)){
	$if_die = 1;
	print "Please input complete arguments.\n";
}elsif($insert_length <= $read_length){
	$if_die = 1;
	print "Insertion length should be larger than read length.\n";
}
die if $if_die == 1;

$fq1 = ">>".$fq1;
$fq2 = ">>".$fq2;
#$out = ">>./".$out;
my %chr_gene_trsc_exon;
my $pre_gene = '';
my @gene_anno;
my @chr;
my $seqID;
my $sim_total;
$ref_dir = $ref_dir."/" unless rindex($ref_dir, "/") == length($ref_dir) - 1;
open GTF, "<", $gtf or die "cannot open gtf file: $!";
open OUT, ">>", $out or die;
while(<GTF>){
	chomp;
	next if /^#/;
	my @line = split /\t/;
	last if ($if_chr == 1 and $line[0] ne "chr1");
	my @atr = split '; ', $line[8];
	if($pre_gene ne $atr[0] and $pre_gene ne ''){
		&split_transcript(@gene_anno);
		@gene_anno = ();
	}
	push @gene_anno, $_;
	$pre_gene = $atr[0];
}

#for my $chr(@chr){
#	for my $gene(keys %{$chr_$gene_trsc_exon{$chr}}){
#		for my $trsc(keys %{$chr_$gene_trsc_exon{$chr}{$gene}}){
#			for $exon(@{$chr_$gene_trsc_exon{$chr}{$gene}{$trsc}}){
#				print OUT "$chr\t$gene\t$trsc\t@$exon\n";
#			}
#		}
#	}
#}
for my $chromo(@chr){
	open CHR, "<", $ref_dir."$chromo.fa" or open CHR, "<", $ref_dir."$chromo.fasta" or die "cannot open the chr fasta file $chromo: $!";
	my $uni_seq = 0;
	my $chr_seq;
	while(<CHR>){
		chomp;
		if(/^>/ and $uni_seq == 0){
			$uni_seq = 1;
		}elsif(/^>/){
			die "There are more than one sequence in $chromo file. Please check!";
		}else{
			$chr_seq .= $_;
		}
	}
	for my $gene(keys %{$chr_gene_trsc_exon{$chromo}}){
		for my $trsc(keys %{$chr_gene_trsc_exon{$chromo}{$gene}}){
			my ($rand_exon1, $rand_exon2);
			my %rand_num1;
			my %rand_num2;
			my ($bingo_num1, $bingo_num2);
			my $cRNA_seq;
			my $trsc_seq;
			my $if_cRNA = int(rand(1)+1);	#the probability of cRNA generated from this transcrpt
			for my $i(0 .. $#{$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}){
				my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$i];
				push @$exon, substr ($chr_seq, $$exon[0]-1, $$exon[1]-$$exon[0]+1);
				push @$exon, substr ($chr_seq, $$exon[0]-3, 2);
				push @$exon, substr ($chr_seq, $$exon[1], 2);
				#print "$chromo\t$gene\t$trsc\t@$exon\n";
				$trsc_seq .= $$exon[-3];
				if ( $if_cRNA == 1 and ( ($$exon[2] eq '-' and $$exon[-2] =~ /AC/i) or ($$exon[2] eq '+' and $$exon[-2] =~ /AG/i) ) ){
					$rand_exon1 ++;
					$rand_num1{$rand_exon1} = $i;
				}
				if ( $if_cRNA == 1 and ( ($$exon[2] eq '-' and $$exon[-1] =~ /CT/i) or ($$exon[2] eq '+' and $$exon[-1] =~ /GT/i) ) ){
					$rand_exon2 ++;
					$rand_num2{$rand_exon2} = $i;
				}
			}
			&simulate_reads2( $rand_mode2, $trsc_seq, $coverage2 ) if length($trsc_seq) > $insert_length;
			if( $if_cRNA == 1 and ( $rand_exon1>=1 and $rand_exon2>=1 ) ){
				$bingo_num1 = int(rand($rand_exon1)+1);
				$bingo_num2 = int(rand($rand_exon2)+1);
				if($rand_num1{$bingo_num1} <= $rand_num2{$bingo_num2}){
					for my $j( $rand_num1{$bingo_num1} .. $rand_num2{$bingo_num2} ){
						my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
						$cRNA_seq .= $$exon[-3];
						#print ">$chromo\t$gene\t$trsc\t$j\t$$exon[0]\t$$exon[1]\t$$exon[2]\t$$exon[-2]\t$$exon[-1]\n";
					}
					#print "?", $cRNA_seq, "\n";
					if ($cRNA_size == 1){
						if ( length($cRNA_seq) <= $read_length ){	#simulate the situation where insert length shorter than read length later
							print OUT "$chromo\t$gene\t$trsc\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num1{$bingo_num1}][0], "\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num2{$bingo_num2}][1], "\n";
							&simulate_reads( $rand_mode, $cRNA_seq, $coverage );
							$sim_total++;
						}
					}elsif($cRNA_size == 2){
						if ( length($cRNA_seq) <= $insert_length/2 ){	#simulate the situation where insert length shorter than read length later
							print OUT "$chromo\t$gene\t$trsc\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num1{$bingo_num1}][0], "\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num2{$bingo_num2}][1], "\n";
							&simulate_reads( $rand_mode, $cRNA_seq, $coverage );
							$sim_total++;
						}
					}else{
						print OUT "$chromo\t$gene\t$trsc\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num1{$bingo_num1}][0], "\t", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$rand_num2{$bingo_num2}][1], "\n";
						&simulate_reads( $rand_mode, $cRNA_seq, $coverage );
						$sim_total++;
					}
				}
			}
		}
	}
}
print OUT "!!total: $sim_total\n";
sub simulate_reads2{
	open FQ1, $fq1 or die;
	open FQ2, $fq2 or die;
	my $mode = shift @_;
	my $trsc_coverage;
	my $seq_length = length($_[0]);
	if ($mode == 1){
		$trsc_coverage = $_[1];
	}else{
		$trsc_coverage = rand($_[1]+1);
	}
	my ($read_num, undef) = sort{$b <=> $a}(int( $seq_length * $trsc_coverage / $read_length / 2 ),1);
	my $err_num = int( $seq_length * $trsc_coverage * $seq_err / 100 );
	my %err_read;
	for (1 .. $err_num){
		my $err_loci = int( rand( ($read_num)*2+1 ) );
		$err_read{$err_loci} ++;
	}
	for my $x( 1 .. $read_num ){
		my $strand = int( rand(2) );
		my ($seq1, $seq2);
		my $ins_len_rand = $insert_length;		#insert length can be simulated later
		my $start_loci = int( rand($seq_length - $ins_len_rand) );
		my $start_loci2 = $start_loci + $ins_len_rand - $read_length;
		$seqID ++;
		if ($strand == 0){
			$seq1 = substr( $_[0], $start_loci, $read_length );
			$seq2 = &comp_rev( substr( $_[0], $start_loci2, $read_length ) );
		}else{
			$seq1 = &comp_rev( substr( $_[0], $start_loci, $read_length ) );
			$seq2 = substr( $_[0], $start_loci2, $read_length );
		}
		my @errs1;
		for (1 .. $err_read{$x * 2 - 1}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs1;
			push @errs1, $err_loci;
			my $ori_base = substr( $seq1, $err_loci, 1 );
			substr( $seq1, $err_loci, 1 ) = &simulate_seq_error($ori_base);
		}
		my @errs2;
		for (1 .. $err_read{$x * 2}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs2;
			push @errs2, $err_loci;
			my $ori_base = substr( $seq2, $err_loci, 1 );
			substr( $seq2, $err_loci, 1 ) = &simulate_seq_error($ori_base);
		}
		print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
		print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
		print FQ1 "$seq1\n";
		print FQ2 "$seq2\n";
		print FQ1 "+\n";
		print FQ2 "+\n";
		print FQ1 ("!" x $read_length) . "\n";
		print FQ2 ("!" x $read_length) . "\n";
	}
}
sub simulate_reads{
	open FQ1, $fq1 or die;
	open FQ2, $fq2 or die;
	open OUT, ">>", $out or die;
	my $mode = shift @_;
	my $seq_length = length($_[0]);
	my $cRNA_coverage;
	my $seq4substr;
	if ($mode == 1){
		$cRNA_coverage = $_[1];
	}else{
		$cRNA_coverage = rand($_[1]+1);
	}
	#my $err_num = int($seq_length*$seq_err/100);
	#my @errs;
	#for (1 .. $err_num){
	#	my $err_loci = int( rand($seq_length) );
	#	redo if $err_loci ~~ @errs;
	#	push @errs, $err_loci;
	#	my $ori_base = substr( $_[0], $err_loci, 1 );
	#	substr( $_[0], $err_loci, 1) = &simulate_seq_error($ori_base);
	#}
	if($cRNA_size == 1 or $cRNA_size == 0){
		$seq4substr = $_[0] x 10;	#different according to $seq_length	my $seq4substr = $_[0] . substr( $_[0], 0, $read_length-1 );
	}else{
		$seq4substr = $_[0] x 10;
	}
	my $read_num = int( $seq_length * $cRNA_coverage / $read_length / 2 );
	my $err_num = int( $seq_length * $cRNA_coverage * $seq_err / 100 );
	my %err_read;
	for (1 .. $err_num){
		my $err_loci = int( rand( ($read_num)*2+1 ) );
		$err_read{$err_loci} ++;
	}
	for my $x( 1 .. $read_num ){
		my $start_loci = int( rand($seq_length) );
		my $strand = int( rand(2) );
		my ($seq1, $seq2);
		my $ins_len_rand = $insert_length;		#insert length can be simulated later
		my $start_loci2 = ( $start_loci + $ins_len_rand - $read_length ) % $seq_length;
		$seqID ++;
		if ($strand == 0){
			$seq1 = substr( $seq4substr, $start_loci, $read_length );
			$seq2 = &comp_rev( substr( $seq4substr, $start_loci2, $read_length ) );
		}else{
			$seq1 = &comp_rev( substr( $seq4substr, $start_loci, $read_length ) );
			$seq2 = substr( $seq4substr, $start_loci2, $read_length );
		}
		my @errs1;
		for (1 .. $err_read{$x * 2 - 1}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs1;
			push @errs1, $err_loci;
			my $ori_base = substr( $seq1, $err_loci, 1 );
			substr( $seq1, $err_loci, 1) = &simulate_seq_error($ori_base);
		}
		my @errs2;
		for (1 .. $err_read{$x * 2}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs2;
			push @errs2, $err_loci;
			my $ori_base = substr( $seq2, $err_loci, 1 );
			substr( $seq2, $err_loci, 1) = &simulate_seq_error($ori_base);
		}
		print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
		print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
		print FQ1 "$seq1\n";
		print FQ2 "$seq2\n";
		print FQ1 "+\n";
		print FQ2 "+\n";
		print FQ1 ("!" x $read_length) . "\n";
		print FQ2 ("!" x $read_length) . "\n";
		print OUT ">\t$x\t$seqID\n";
		print OUT "**\t1\n" if $seq_length - $start_loci <= 69;
		print OUT "**\t2\n" if $seq_length - $start_loci2 <= 69;
	}
}
sub simulate_seq_error{
	my $ori_base = $_[0];
	my @base = ('A', 'T', 'C', 'G');
	my $err_base_index;
	for my $i( 0 .. $#base ){
		if ($base[$i] =~ /$ori_base/i){
			while(1){
				$err_base_index = int(rand(4));
				last unless $err_base_index  == $i;
			}
			last;
		}
	}
	$base[$err_base_index];
}
sub comp_rev{
	my $seq = reverse($_[0]);
	$seq =~ s/[Aa]/X/g;
	$seq =~ s/[Tt]/A/g;
	$seq =~ s/X/T/g;
	$seq =~ s/[Cc]/Y/g;
	$seq =~ s/[Gg]/C/g;
	$seq =~ s/Y/G/g;
	$seq;
}
sub split_transcript{
	#my $gene;
	for (@_){
		my @line = split /\t/;
		if($line[2] eq 'exon'){
			my @atr = split ('; ', $line[8], 3);
			if ($atr[1] =~ /transcript_id \"(\w+\.\w+)\"/){
				push @{$chr_gene_trsc_exon{$line[0]}{$atr[0]}{$1}}, [ $line[3], $line[4], $line[6] ];
				#$gene = $atr[0];
			}else{
				print "error: no transcript_id found for $atr[0]!\n";
			}
		}
	}
	my @line2 = split (/\t/, $_[0], 2);
	push @chr, $line2[0] unless $line2[0] ~~ @chr;
}

