#!/usr/bin/perl -w
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $ver="2.0";

############################################
=c
my %opts;
GetOptions(\%opts,"id=s","od=s","allsamplecfg=s","h");
if (!defined($opts{id})||!defined($opts{od})||!defined($opts{allsamplecfg})||defined($opts{h})) {
	&help();
}
=cut
my ($id,$od,$data_cfg);
GetOptions(
	"id:s"=>\$id,
	"od:s"=>\$od,
	"cfg:s"=>\$data_cfg,
	"h|?"=>\&help,
)or &help;
&help unless ($id and $od and $data_cfg);

###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";


###############
$id = &ABSOLUTE_DIR($id);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);

my $QC="$od/QC";
&MKDIR($QC);
&MKDIR("$QC/Clean_Data");
&MKDIR("$QC/Map_assess");
my $GENE="$od/mRNA";
&MKDIR($GENE);
&MKDIR("$GENE/NewGene");
&MKDIR("$GENE/NewGene/Annotation");
&MKDIR("$GENE/Structure");
my $DEG="$GENE/DEG";
&MKDIR($DEG);
&MKDIR("$GENE/Expression");
my $LncRNA="$od/LncRNA";
&MKDIR($LncRNA);
&MKDIR("$LncRNA/Identify");
&MKDIR("$LncRNA/Target");
&MKDIR("$LncRNA/Target/Trans_target/"); ###by linhj
&MKDIR("$LncRNA/Expression");
&MKDIR("$LncRNA/Assembly");
&MKDIR("$LncRNA/DEG");
&MKDIR("$LncRNA/Compare");

my %ID;
open(IN,$data_cfg) or die $!;
while(<IN>){
	chomp;
	next if (/^#/||/^\s+$/);
	if (/^Sample/){
		my ($bmk_id,$sampleid)=(split /\s+/,$_,3)[1,2];
        my @IDline = split("\t",$_); 
        my $IDnum = @IDline;
        if ($IDnum == 2) {
            $ID{$bmk_id} = $bmk_id;
        }else{
		    $ID{$bmk_id} = $sampleid;
        }
	}
}
close IN;

############### Extract 
if (-d "$id/Data_Assess") {
	#system "cp $id/Data_Assess/AllSample_GC_Q.stat $QC/Clean_Data";
	my $allsample=glob "$id/Data_Assess/AllSample_GC_Q.stat";
	open (IN,$allsample) or die $!;
	open (OUT,">$QC/Clean_Data/AllSample_GC_Q.stat") or die $!;
	print OUT "#SampleID\tBMK-ID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ30(%)\n";
	while(<IN>){
		chomp;
		next if (/^#/||/^\s+$/);
		my ($bmkid)=(split /\t+/,$_)[0];
		print OUT "$ID{$bmkid}\t$_\n";
	}
	close IN;
	close OUT;
	#`awk '{print \$1}' $allsamplecfg |paste - $allsample >$QC/Clean_Data/AllSample_GC_Q.stat`;
	#`iconv -f 'GB2312' -t 'utf-8' $QC/Clean_Data/AllSample_GC_Q.stat -o $QC/Clean_Data/AllSample_GC_Q.stat.custom`;
	system "cp -r $id/Data_Assess/PNG $QC/Clean_Data";
	#system "cp $Bin/readme/readme/rawdata.README.txt $rawdata_od";
}

if (-d "$id/Basic_Analysis/geneExpression") {
	system "cp $id/Basic_Analysis/geneExpression/*.png $QC/Map_assess";
	system "cp -r $id/Basic_Analysis/geneExpression/*.xls $GENE/Expression";
	system "rm $GENE/Expression/*isoExpression.xls";
	my %map_stat;
	my @map = glob "$id/Basic_Analysis/geneExpression/*.mappedStat.xls";
	open (OUT,">$QC/Map_assess/Total_mapped.stat") or die $!;
	print OUT "BMK-ID\tTotal Reads\tMapped Reads\tUniq Mapped Reads\tMultiple Mapped Reads\tReads Map to '+'\tReads Map to '-'\n";
	foreach my $file(@map){
		open (IN,$file) or die $!;
		basename($file)=~/^(.*)\.mappedStat.xls/;
		my $name=$1;
		while(<IN>){
			chomp;
			my @div=(split /\t/,$_);
			if (/^Total Reads/){
				$map_stat{$name}{'Total Reads'}=$div[1];
			}
			if (/^mapped Reads/){
				$map_stat{$name}{'mapped Reads'}="$div[1]($div[2])";
			}
			if (/^Uniq Map/){
				$map_stat{$name}{'Uniq Map'}="$div[1]($div[2])";
			}
			if (/^Multiple Map/){
				$map_stat{$name}{'Multiple Map'}="$div[1]($div[2])";
			}
			if (/^Left Map/){
				$map_stat{$name}{'Left Map'}="$div[1]($div[2])";
			}
			if (/^Right Map/){
				$map_stat{$name}{'Right Map'}="$div[1]($div[2])";
			}
		}
		close IN;
	}
	foreach my $k (sort keys  %map_stat){
		print OUT "$k\t$map_stat{$k}{'Total Reads'}\t$map_stat{$k}{'mapped Reads'}\t$map_stat{$k}{'Uniq Map'}\t$map_stat{$k}{'Multiple Map'}\t$map_stat{$k}{'Left Map'}\t$map_stat{$k}{'Right Map'}\n";
	}
	close OUT;
	system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene.longest_transcript.fa $GENE/NewGene";
	system "sed -i 's/\t/ /'  $GENE/NewGene/*.newGene.longest_transcript.fa";###20160115
	system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene_final.filtered.gff $GENE/NewGene";
	#system "cp $Bin/readme/readme/geneExpression.README.txt $od/geneExpression";
}

if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/") {
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/*.png $LncRNA/Identify";
	system"rm $LncRNA/Identify/class_code.stat.png" if (-f "$LncRNA/Identify/class_code.stat.png");
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPC/lnc_code_filter.result.txt $LncRNA/Identify/CPC.txt";
	
	my $pfamfile="$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/Pfam/Pfam_result.txt" ;
	open IN,"$pfamfile";
	open OUT,">$LncRNA/Identify/Pfam.txt";
	print OUT "trans_id\thmm_acc\thmm_name\thmm_start\thmm_end\thmm_length\tbit_score\tE-value\n";
	while (<IN>){
		chomp;
		next if /^#/ || /^$/;
		my @tmp=split /\s+/,$_;
		my $id=(split /\:/,$tmp[0])[0];
		print OUT "$id\t$tmp[5]\t$tmp[6]\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\n";
	}
	close OUT;
	close IN;
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CNCI/CNCI.index $LncRNA/Identify/CNCI.txt";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt $LncRNA/Identify/cpat.txt"  if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt");
	
	#############################
	
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/venn.png $LncRNA/Identify";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/list.txt $LncRNA/Identify/CPC_CNCI_Pfam_CPAT.txt";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list $LncRNA/Identify/LncRNA.id";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_final.fa $LncRNA/Identify/LncRNA.fa";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gtf $LncRNA/Identify/LncRNA.gtf";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gff $LncRNA/Identify/LncRNA.gff";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Cis_target/Cis_target_result.txt  $LncRNA/Target/lncRNA_position.target";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.txt  $LncRNA/Target/Trans_target/" if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/");  ####by linhj
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.png  $LncRNA/Target/Trans_target/" if (-d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/");  ####by linhj
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target $LncRNA/Target/lncRNA_basepair.target" if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target");
	my @assembly=glob"$id/Basic_Analysis/Tophat_Cufflinks/Cufflinks/*/transcripts.gtf";
	foreach my $ass (@assembly){
		$ass=~/\/(\w\d+)\//;
		my $name=$1;
		system "cp $ass $LncRNA/Assembly/$name.Cufflinks.transcripts.gtf";
		system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $LncRNA/Assembly/$name.Cufflinks.transcripts.gtf";
	}
	my @scripture=glob"$id/Basic_Analysis/Tophat_Cufflinks/Scripture/*.gtf";
	foreach my $scr(@scripture){
		if ($scr=~/\/(\w\d+).gtf/){
			my $sampl=$1;
			system "cp -r $scr $LncRNA/Assembly/$sampl.Scripture.transcripts.gtf";
			system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $LncRNA/Assembly/$sampl.Scripture.transcripts.gtf";
		}
	}
}

if (-d "$id/Alitsplice_Analysis") {
	&MKDIR("$GENE/Structure/Alt_splice");
	system "cp -r $id/Alitsplice_Analysis/*.png $GENE/Structure/Alt_splice";
	system "cp -r $id/Alitsplice_Analysis/*.fpkm $GENE/Structure/Alt_splice";
	#system "cp $Bin/readme/readme/Alt_splice.README.txt $od/Alt_splice";
}

if (-d "$id/DEG_Analysis") {
	my @DEG_dir=glob "$id/DEG_Analysis/*";
	my %Anno;
	my %Stat;
	system "cp  $id/DEG_Analysis/*.png $DEG" if (-f "$id/DEG_Analysis/All_DEG_veen.png");
	system "cp $id/DEG_Analysis/All_gene_counts.list $DEG" if (-f "$id/DEG_Analysis/All_gene_counts.list");
	system "cp $id/DEG_Analysis/All_gene_fpkm.list $DEG" if (-f "$id/DEG_Analysis/All_gene_fpkm.list");
	foreach my $dir (@DEG_dir) {
		if (-d $dir) {
			if ($dir=~/\/All_DEG$/) {
				my $name=basename($dir);
				&MKDIR("$DEG/$name");
				system "cp -r $dir/* $DEG/$name/";
				system "rm $DEG/All_DEG/group.dat" if (-f "$DEG/All_DEG/group.dat");
				system "rm $DEG/All_DEG/Rplots.pdf" if (-f "$DEG/All_DEG/Rplots.pdf");
			}
			elsif ($dir=~ /\/kmean_0$/) {
				system "cp -r $dir/*.png $DEG";
			}
			elsif ($dir=~/\/density$/) {
				system "cp -r $dir/*.png $GENE/Expression/";
				system "cp -r $dir/*.cor $GENE/Expression/";
				system "cp -r $dir/cor_plot $GENE/Expression/";
            }
			elsif ($dir=~/_vs_/){
				my $nam=basename($dir);
				&MKDIR("$DEG/$nam");
				system "cp -r $dir/Anno_enrichment/* $DEG/$nam/";
                system "cp -r $dir/DEG_Cluster $DEG/$nam"; ###by linhj
				system "cp -r $dir/*.DEG_final.xls $DEG/$nam";
				system "cp -r $dir/*.DEG_gene.fa $DEG/$nam";
				system "cp -r $dir/*.DEG.CytoscapeInput.txt $DEG/$nam";
				system "cp -r $dir/*.png $DEG/$nam" unless $dir=~/.*_vs_.*_vs_.*/;
				system "rm $DEG/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat" if (-f "$DEG/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat");
				system "rm $DEG/$nam/pathway/kegg_enrichment/*.Kegg.ko" if (-f "$DEG/$nam/pathway/kegg_enrichment/*.Kegg.ko");
				system "rm -r $DEG/$nam/Kog_Anno" if (-d "$DEG/$nam/Kog_Anno");
				#system "rm $DEG/$nam/go_enrichment/*.svg.list" if (-f "$DEG/$nam/go_enrichment/*.svg.list");
                system "rm $DEG/$nam/Cog_Anno/*.plot.log" if (-f "$DEG/$nam/Cog_Anno/*.plot.log");
                system "rm $DEG/$nam/Graph/topGO.*";
				system "rm -r $DEG/$nam/*_File_*" if (-f "$DEG/$nam/*_File_*");

				if (-d "$dir/Graphic") {
					system "cp $dir/Graphic/*.png $DEG/$nam";
					system "cp $dir/Graphic/*.pdf $DEG/$nam";
				}
				my %Site;
				my $file=(glob "$dir/Anno_enrichment/*.annotation.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					if (/^\#/) {
						my @Anno=split/\s+/,$_;
						for (my $s=0;$s<@Anno ;$s++) {
							if ($Anno[$s] eq 'COG_class') {
								$Site{'COG'}=$s;
							}
							if ($Anno[$s] eq 'KOG_class') {
								$Site{'KOG'}=$s;
							}
                            if ($Anno[$s] eq 'Swissprot_annotation') {
								$Site{'Swiss-Prot'} = $s;
                            }
							elsif ($Anno[$s]=~/^([^_]+)_annotation/) {
								$Site{$1}=$s;
							}
						}
                        $Anno{$nam}{'Annotated'} = 0;
					}
					else{
						my @Info=split /\t+/,$_;
						foreach my $key (keys %Site) {
							$Anno{$nam}{$key}||=0;
							$Anno{$nam}{$key}++ unless ($Info[$Site{$key}] eq '--');
						}
                        $Anno{$nam}{'Annotated'} ++;
					}
				}
				close IN;
				$file=(glob "$dir/*.DEG_final.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					next if /^\#/;
					my $type=(split/\s+/,$_)[-1];
					$Stat{$nam}{up}++ if $type eq 'up';
					$Stat{$nam}{down}++ if $type eq 'down';
					$Stat{$nam}{total}++;
				}
				close IN;
			}
		}
	}
	open (OUT,">$DEG/DEG.anno.stat") or die $!;
	my $limit_anno=0;
	foreach my $key (sort keys %Anno) {#!
		if ($Anno{$key}{'Annotated'}!=0){
		if ($limit_anno==0) {
			print OUT "#DEG Set";
			foreach my $key1 (sort keys %{$Anno{$key}}) {
				print OUT "\t$key1";
			}
			print OUT "\n";
			$limit_anno++;
		}
		print OUT "$key";
		foreach my $key1 (sort keys %{$Anno{$key}}) {
			print OUT "\t$Anno{$key}{$key1}";
		}
		print OUT "\n";
	}else{
		if ($limit_anno==0) {
			print "#DEG Set\tAnnotated\tCOG\tGO\tKEGG\tSwiss-prot\tTrEMBL\tnr\tnt\n";
			$limit_anno++;
			}
			print "$key\t0\t0\t0\t0\t0\t0\t0\t0\n";
		}
	}
	close OUT;
	open (OUT,">$DEG/DEG.stat") or die $!;
	print OUT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		if ($Stat{$key}{total}!=0){
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}else{
		print OUT "$key\t0\t0\t0\n";
	}
	}
	close OUT;
	#system "cp $Bin/readme/readme/DEG_Analysis.README.txt $DEG_od";
}

if (-d "$id/Anno_Integrate/New_Anno/Result/") {
    system "cp -r $id/Anno_Integrate/New_Anno/Result/* $GENE/NewGene/Annotation";
    system "rm -r $GENE/NewGene/Annotation/Blast2go" if (-d "$GENE/NewGene/Annotation/Blast2go");
    system "rm $GENE/NewGene/Annotation/*.annot" if (-d "$GENE/NewGene/Annotation/*.annot");
	system "rm $GENE/NewGene/Annotation/*.svg";
    system "rm $GENE/NewGene/Annotation/*.svg.list" if (-f "$GENE/NewGene/Annotation/*.svg.list");
	my $function_anno = "$GENE/NewGene/Annotation/Function_Annotation.stat.xls";
	my $tmp = "$GENE/NewGene/Annotation/tmp.txt";
	open (IN,$function_anno) or die $!;
	open (OUT,">$tmp") or die $!;
	print OUT "Annotated databases\tNew Gene Number\n";
	while (<IN>) {
		chomp;
		next if (/^#/);
		my @div=(split /\t/,$_);
		print OUT "$div[0]\t$div[1]\n";
	}
	close IN;
	close OUT;
	system "rm $function_anno";
	system "mv $tmp $function_anno";
	#system "cp $Bin/readme/readme/NewGene.README.txt $od/NewGene";
}
if (-d "$id/Lnc_Diff_Analysis") {
	my @Lnc_Diff_dir=glob "$id/Lnc_Diff_Analysis/*";
	my %Anno;
	my %Stat;
        system "cp $id/Lnc_Diff_Analysis/*.png $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/All_DEG_venn.png");
        system "cp $id/Lnc_Diff_Analysis/*.genes $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/All_DEG_venn.genes");
        #system "cp $id/Lnc_Diff_Analysis/All_trans_counts.list $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/All_trans_counts.list");
        #system "cp $id/Lnc_Diff_Analysis/all_trans_fpkm.list $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/all_trans_fpkm.list");
        system "cp $id/Lnc_Diff_Analysis/LncRNA_counts.list $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/LncRNA_counts.list");
        system "cp $id/Lnc_Diff_Analysis/LncRNA_fpkm.list $LncRNA/DEG" if (-f "$id/Lnc_Diff_Analysis/LncRNA_fpkm.list");
	foreach my $dir (@Lnc_Diff_dir) {
            
		if (-d $dir ){
			if ($dir=~/\/All_DEG$/) {
				#my $name=basename($dir);
				#&MKDIR("$LncRNA/DEG/$name");
				system "cp -r $dir $LncRNA/DEG/";
				system"rm $LncRNA/DEG/All_DEG/group.dat" if (-f "$LncRNA/DEG/All_DEG/group.dat");
			}
			elsif($dir=~/\/kmean_0$/){
				system"cp -r $dir/k-means.png $LncRNA/DEG/" if (-f "$dir/k-means.png");
				#system "rm -r $LncRNA/DEG/kmean_0";
			}
			#elsif ($dir=~/\/Venn$/) {
			#	system "rm -r $LncRNA/DEG/Venn";
			#}
			elsif ($dir=~/\/density$/) {
				system "cp -r $dir/*.png $LncRNA/Expression/" if (-f "$dir/*.png");
				system "cp -r $dir/*.cor $LncRNA/Expression/" if (-f "$dir/*.cor");
				system "cp -r $dir/cor_plot $LncRNA/Expression/" if (-f "$dir/cor_plot");
				#system "rm -r $LncRNA/DEG/density/";
			}
			elsif ($dir=~/_vs_/){
				my $nam=basename($dir);
				&MKDIR("$LncRNA/DEG/$nam");
				system "cp -r $dir/Anno_enrichment/* $LncRNA/DEG/$nam/";
                system "cp -r $dir/DEG_Cluster $LncRNA/DEG/$nam"; ###by linhj
				system "cp -r $dir/*.DEG_final.xls $LncRNA/DEG/$nam/";
                                system "cp -r $dir/*deg_lnc2Target_m_Cytoscape.input.txt $LncRNA/DEG/$nam/" if (-f  "$dir/*deg_lnc2Target_m_Cytoscape.input.txt");
                                system "cp -r $dir/*DEG.CytoscapeInput.txt $LncRNA/DEG/$nam/" if (-f  "$dir/*DEG.CytoscapeInput.txt");
				system "cp -r $dir/*.png $LncRNA/DEG/$nam" unless $dir=~/.*_vs_.*_vs_.*/;
				system "rm -r $LncRNA/DEG/$nam/pathway/kegg_enrichment/*.Kegg.ko" unless (-e "$LncRNA/DEG/$nam/pathway/kegg_enrichment/*.Kegg.ko");
				system "rm -r $LncRNA/DEG/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat" unless (-e "$LncRNA/DEG/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat");
				system "rm $LncRNA/DEG/$nam/go_enrichment/*.svg.list" unless (-e "$LncRNA/DEG/$nam/go_enrichment/*.svg.list");
				system "rm $LncRNA/DEG/$nam/Cog_Anno/*.plot.log" unless (-e "$LncRNA/DEG/$nam/Cog_Anno/*.plot.log");
				system "rm $LncRNA/DEG/$nam/Graph/topGO.*";
				system "rm -r $LncRNA/DEG/$nam/*_File*" if (-f "$LncRNA/DEG/$nam/*_File*");

				my %Site;
				my $file=(glob "$dir/Anno_enrichment/*.annotation.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					if (/^\#/) {
						my @Anno=split/\s+/,$_;
						for (my $s=0;$s<@Anno ;$s++) {
							if ($Anno[$s] eq 'COG_class') {
								$Site{'COG'}=$s;
							}
							if ($Anno[$s] eq 'KOG_class') {
								$Site{'KOG'}=$s;
							}
							if ($Anno[$s] eq 'Swissprot_annotation') {
								$Site{'Swiss-Prot'} = $s;
							}
							elsif ($Anno[$s]=~/^([^_]+)_annotation/) {
								$Site{$1}=$s;
							}
						}
						$Anno{$nam}{'Annotated'} = 0;
					}
					else{
						my @Info=split /\t+/,$_;
						foreach my $key (keys %Site) {
							$Anno{$nam}{$key}||=0;
							$Anno{$nam}{$key}++ unless ($Info[$Site{$key}] eq '--');
						}
						$Anno{$nam}{'Annotated'} ++;
					}
				}
				close IN;
				$file=(glob "$dir/*.DEG*final.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					next if /^\#/;
					my $type=(split/\s+/,$_)[-1];
					$Stat{$nam}{up}++ if $type eq 'up';
					$Stat{$nam}{down}++ if $type eq 'down';
					$Stat{$nam}{total}++;
				}
				close IN;
			}
		}
	}
	open (OUT,">$LncRNA/DEG/Diff_lnc.anno.stat") or die $!;
	my $limit_anno=0;
	foreach my $key (sort keys %Anno) {
		if ($Anno{$key}{'Annotated'}!=0){
		if ($limit_anno==0) {
			print OUT "#DEG Set";
			foreach my $key1 (sort keys %{$Anno{$key}}) {
				print OUT "\t$key1";
			}
			print OUT "\n";
			$limit_anno++;
		}
		print OUT "$key";
		foreach my $key1 (sort keys %{$Anno{$key}}) {
			print OUT "\t$Anno{$key}{$key1}";
		}
		print OUT "\n";
	}else{
		if ($limit_anno==0) {
			print OUT "#DEG Set\tAnnotated\tCOG\tGO\tKEGG\tSwiss-prot\tTrEMBL\tnr\tnt\n";
			$limit_anno++;
		}
			print OUT "$key\t0\t0\t0\t0\t0\t0\t0\t0\n";
		}
	}
	close OUT;
	open (OUT,">$LncRNA/DEG/Diff_Lnc.stat") or die $!;
	print OUT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		if ($Stat{$key}{total}!=0){
			$Stat{$key}{up}||=0;
			$Stat{$key}{down}||=0;
			print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
		}else{
			print OUT "$key\t0\t0\t0\n";
		}
	}
	close OUT;
	#system "cp $Bin/readme/readme/DEG_Analysis.README.txt $Lnc_Diff_od";
}

if(-d "$id/SNP_Analysis") {
	system "cp -r $id/SNP_Analysis/stat $GENE/Structure/";
	system "mv -f $GENE/Structure/stat/ $GENE/Structure/SNP/";
	#system "cp $Bin/readme/readme/SNP.README.txt $od/SNP_Analysis";
}

if (-d "$id/Gene_Structure_Optimize") {
	system "cp -r $id/Gene_Structure_Optimize $GENE/";
	#system "cp $Bin/readme/readme/Gene_Structure_Optimize.README.txt $od/Gene_Structure_Optimize";
}

if (-d "$id/Compare_analysis"){
	&MKDIR("$LncRNA/Compare/Length");
	&MKDIR("$LncRNA/Compare/Exon");
	&MKDIR("$LncRNA/Compare/ORF");
	&MKDIR("$LncRNA/Compare/Expression/");
	&MKDIR("$LncRNA/Compare/Isoform");
    system "cp -r $id/Compare_analysis/chr_num_stat $LncRNA/Compare/";
	system "cp -r $id/Compare_analysis/*RNA.len.png $LncRNA/Compare/Length";
	system "cp -r $id/Compare_analysis/*exon.png $LncRNA/Compare/Exon";
	system "cp -r $id/Compare_analysis/*orf.png $LncRNA/Compare/ORF";
	system "cp -r $id/Compare_analysis/*FPKM.xls $LncRNA/Compare/Expression/";
	system "cp -r $id/Compare_analysis/*fpkm.png $LncRNA/Compare/Expression/";
	system "cp -r $id/Compare_analysis/*isoform.png $LncRNA/Compare/Isoform";
	system "cp -r $id/Compare_analysis/*RNA_isoform.txt $LncRNA/Compare/Isoform";	
}

if (-d "$id/DEU_analysis") {
	&MKDIR("$GENE/DEU");
	my @DEU_dir = glob"$id/DEU_analysis/*";
	foreach my $dir (@DEU_dir) {
		if ($dir =~ /_vs_/) {
			my $name = basename($dir);
			#&MKDIR("$GENE/DEU/$name/DEXSeqReport");
			#&MKDIR("$GENE/DEU/$name/DEXSeqReport/files");
			system "cp -r $dir $GENE/DEU";
			system "rm $GENE/DEU/$name/DEU_Result_All.xls" if (-f "$GENE/DEU/$name/DEU_Result_All.xls");
			#system "rm -r $GENE/DEU/$name/DEXSeqReport/files/*.svg" if (-f "$GENE/DEU/$name/DEXSeqReport/files/*.svg");
			#system"cp -r $dir/DEU_Result_Final.xls $GENE/DEU/$name";
			#system"cp -r $dir/DEXSeqReport/testForDEU.html $GENE/DEU/$name/DEXSeqReport";
			#system"cp -r $dir/DEXSeqReport/files/*.html $GENE/DEU/$name/DEXSeqReport/files";
		}
	}
}

if (-d "$id/Personality"){
	&MKDIR("$od/Personality");
	if (-d "$id/Personality/Conservation"){
		mkdir "$od/Personality/Conservation" unless (-d "$od/Personality/Conservation");
		system "cp $id/Personality/Conservation/*.png $od/Personality/Conservation";
	}
	if (-d "$id/Personality/Known_LncRNA"){
		&MKDIR("$od/Personality/Known_LncRNA");
		system "cp $id/Personality/Known_LncRNA/known_lncRNA.result.txt $od/Personality/Known_LncRNA";
	}
	if (-d "$id/Personality/miRNA_Target2LncRNA"){
		&MKDIR("$od/Personality/miRNA_Target2LncRNA");
		system "cp $id/Personality/miRNA_Target2LncRNA/*.mir2target.list $od/Personality/miRNA_Target2LncRNA";
		system "cp $id/Personality/miRNA_Target2LncRNA/*.mir2target.stat $od/Personality/miRNA_Target2LncRNA";
	}
	if (-d "$id/Personality/precursor"){
		&MKDIR("$od/Personality/precursor");
		system "cp $id/Personality/precursor/lncRNA_precursor.txt $od/Personality/precursor";
	}
	if (-d "$id/Personality/TF"){
		&MKDIR("$od/Personality/TF");
		system "cp $id/Personality/TF/final.xls $od/Personality/TF/TF.txt";
	}
	if (-d "$id/Personality/Tissue_specific"){
		&MKDIR("$od/Personality/Tissue_specific");
		system "cp $id/Personality/Tissue_specific/*.png $od/Personality/Tissue_specific";
	}
	if(-d "$id/Personality/PCA"){
		&MKDIR("$od/Personality/PCA");
		system "cp $id/Personality/PCA/*.png $od/Personality/PCA";
	}
	if(-d "$id/Personality/WGCNA"){
		&MKDIR("$od/Personality/WGCNA");
		system "cp $id/Personality/WGCNA/result/wgcna/*.png $od/Personality/WGCNA";
	}
}
#system "cp $Bin/../README $od/";

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $Script Time :[$Time_End]\n\n";
###############Subs
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

sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract LncRNA Analysis Results for Html Process;
	version:$ver
	Usage:
		-id  <STR>   input dir, analysis output directory   force
		-od  <STR>   result output dir                      force
		-cfg		data.cfg			force
		-h           help
	Usage End.
		exit;
}
