#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name = basename($0);

my $ver = "2.0";
############################################
my %opts;
GetOptions( \%opts, "id=s", "od=s", "h" );
if ( !defined( $opts{id} ) || !defined( $opts{od} ) || defined( $opts{h} ) ) {
    &help();
}
###############Time
my $BEGIN = time();
my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\n[$Time_Start] $Script start ... \n\n";

###############
#my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $id = &ABSOLUTE_DIR( $opts{id} );
&MKDIR( $opts{od} );
my $od     = &ABSOLUTE_DIR( $opts{od} );
my $DEG_od = "$od/DEG_Analysis";
&MKDIR($DEG_od);
my $Lnc_Diff_od = "$od/Lnc_Diff_Analysis";
&MKDIR($Lnc_Diff_od);
my $rawdata_od = "$od/rawdata";
&MKDIR($rawdata_od);
&MKDIR("$od/LncRNA");
&MKDIR("$od/geneExpression");
&MKDIR("$od/LncRNA/LncExpression");
&MKDIR("$od/NewGene");
&MKDIR("$od/Allgene");
&MKDIR("$od/NewGene/NewGene_Anno");
&MKDIR("$od/Allgene/Allgene_Anno");
&MKDIR("$od/SNP_Analysis");
&MKDIR("$od/Compare_Analysis");
&MKDIR("$od/Compare_Analysis/Length");
&MKDIR("$od/Compare_Analysis/Exon");
&MKDIR("$od/Compare_Analysis/Orf_len");
&MKDIR("$od/Compare_Analysis/Expression/");
&MKDIR("$od/Compare_Analysis/Isoform/");
&MKDIR("$od/LncRNA/Trans_target");
&MKDIR("$od/Compare_Analysis/lncRNA_mRNA_DEG");

############### Extract Assembly dir
if ( -d "$id/Data_Assess" ) {
    system "cp $id/Data_Assess/AllSample_GC_Q.stat $rawdata_od";
    system "cp -r $id/Data_Assess/PNG $rawdata_od";
    system "cp $Bin/readme/readme/rawdata.README.txt $rawdata_od";
}

if ( -d "$id/Basic_Analysis/geneExpression" ) {
    system "cp $id/Basic_Analysis/geneExpression/*.png $od/geneExpression";
    system "cp $id/Basic_Analysis/geneExpression/*.xls $od/geneExpression";
	my %map_stat;
	my @map = glob "$id/Basic_Analysis/geneExpression/*.mappedStat.xls";
	open (OUT,">$od/geneExpression/Total_mapped.stat") or die $!;
	print OUT "BMK-ID\tTotal Reads\tMapped Reads\tUniq Mapped Reads\tMultiple Mapped Reads\tReads Map to '+'\tReads Map to '-'\n";
	foreach my $file(@map){
		open (IN,$file) or die $!;
		basename($file)=~/^(\w*)\./;
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
			if (/^Only Map Plus Strand/){
				$map_stat{$name}{'Left Map'}="$div[1]($div[2])";
			}
			if (/^Only Map Minus Strand/){
				$map_stat{$name}{'Right Map'}="$div[1]($div[2])";
			}
		}
		close IN;
	}
	foreach my $k (sort keys  %map_stat){
		print OUT "$k\t$map_stat{$k}{'Total Reads'}\t$map_stat{$k}{'mapped Reads'}\t$map_stat{$k}{'Uniq Map'}\t$map_stat{$k}{'Multiple Map'}\t$map_stat{$k}{'Left Map'}\t$map_stat{$k}{'Right Map'}\n";
	}
	close OUT;
    system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene.longest_transcript.fa $od/NewGene";
    system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene_final.filtered.gff $od/NewGene";
    system "cp $Bin/readme/readme/geneExpression.README.txt $od/geneExpression";
}
if ( -d "$id/Basic_Analysis/LncExpression" ) {
    system "cp $id/Basic_Analysis/LncExpression/*xpression.xls $od/LncRNA/LncExpression";
}

if ( -d "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/" ) {
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/*.png $od/LncRNA/"; #    任务1，饼图的绘制改动地方，修改文件的指向路径--------------------------------------------------------------------------------------------------------------------------------------------
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/*.pdf $od/LncRNA/";   #
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPC/lnc_code_filter.result.txt $od/LncRNA/CPC.xls";

#system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/Pfam/Pfam_result.txt  $od/LncRNA/Pfam.xls";
    my $pfamfile = "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/Pfam/Pfam_result.txt"; #unless (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/Pfam/Pfam_result.txt");
    open IN,  "$pfamfile";
    open OUT, ">$od/LncRNA/Pfam.xls";
    print OUT "#trans_id\thmm_acc\thmm_name\thmm_start\thmm_end\thmm_length\tbit_score\tE-value\n";
    while (<IN>) {
        chomp;
        next if /^#/ || /^$/;
		my @tmp=split /\s+/,$_;
        my $id = ( split /\:/, $tmp[0] )[0];
        print OUT "$id\t$tmp[5]\t$tmp[6]\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t$tmp[12]\n";
    }
    close OUT;
    close IN;
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CNCI/CNCI.index $od/LncRNA/CNCI.xls";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt $od/LncRNA/cpat.txt"  if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/CPAT/cpat.txt");
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/venn.png $od/LncRNA";
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/code_filter/list.txt $od/LncRNA/CPC_CNCI_Pfam_CPAT.txt";#code_filter.result.txt";
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_id.list $od/LncRNA/LncRNA.id";
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/lnc_filter_final.fa $od/LncRNA/LncRNA.fa";
    system "cp -r $id/Compare_analysis/chr_num_stat/ $od/LncRNA"; #  任务2lncRNAcount画图，改动的地方，新的文件指向新的路径----------------------------------------------------------------------------------------------------------------------------------------------------------------
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_filter/filter_final.gtf $od/LncRNA/LncRNA.gtf";
    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/novel_lncRNA_target.xls  $od/LncRNA/lncRNA_target.xls";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Cis_target/Cis_target_result.txt  $od/LncRNA/lncRNA_position.target";
	system "cp $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target $od/LncRNA/lncRNA_basepair.target" if (-f "$id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/LncTar/LncTar_basepair.target");

#	system "cp -r $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target  $od/LncRNA/";#      任务1WGCNA结果，改动的地方文件输出的指向路径-------------------------------------------------------------------------------------------------------------------------------------------------
    system "cp -r $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/Fig* $od/LncRNA/Trans_target";
    system "cp -r $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.txt $od/LncRNA/Trans_target";
    system "cp -r $id/Basic_Analysis/Tophat_Cufflinks/Lnc_target_predict/Trans_target/*.xls $od/LncRNA/Trans_target";

#    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Cufflinks/*1/transcripts.gtf $od/LncRNA/Cufflinks.partial.transcripts.gtf";    #修改 2015/10/13 by luml
#    system "cp $id/Basic_Analysis/Tophat_Cufflinks/Scripture/*1.gtf $od/LncRNA/Scripture.partial.transcripts.gtf";
	my @assembly=glob"$id/Basic_Analysis/Tophat_Cufflinks/Cufflinks/*/transcripts.gtf";
	system "cp $assembly[0] $od/LncRNA/Cufflinks.partial.transcripts.gtf";
	foreach my $ass (@assembly){
		$ass=~/\/Cufflinks\/(.*)\/transcripts.gtf/;
		my $name=$1;
		system "cp $ass $od/LncRNA/$name.Cufflinks.transcripts.gtf";
		system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $od/LncRNA/$name.Cufflinks.transcripts.gtf";
	}
	my @scripture=glob"$id/Basic_Analysis/Tophat_Cufflinks/Scripture/*.gtf";
	system "cp $scripture[0] $od/LncRNA/Scripture.partial.transcripts.gtf";
	foreach my $scr(@scripture){
		if ($scr=~/\/Scripture\/(.*).gtf/){
			my $sampl=$1;
			system "cp -r $scr $od/LncRNA/$sampl.Scripture.transcripts.gtf";
			system "sed -i '1i Chr\tSource\tType\tStart Site\tEnd Site\tScore\tStrand\tframe\tDescription' $od/LncRNA/$sampl.Scripture.transcripts.gtf";
		}
	}
}

if ( -d "$id/Alitsplice_Analysis" ) {
=pod
    &MKDIR("$od/Alt_splice");
    my @ALT_dir = glob "$id/Alitsplice_Analysis/*";
    foreach my $dir (@ALT_dir) {
        next if ( $dir =~ /work_sh/ );
        if ( -d $dir ) {
            $dir =~ m/.*\/(\S+)/;
            my $nam = $1;
            &MKDIR("$od/Alt_splice/$nam");

            #			system "mkdir $od/Alt_splice/$nam/gff";
            #system "cp -r $dir/AS_Result/gff $od/Alt_splice/$nam/";
            system "cp $dir/AS_Result/*.SpliceGrapher.gff $od/Alt_splice/$nam/";
            system "mkdir $od/Alt_splice/$nam/pdf"
              unless ( -d "$od/Alt_splice/$nam/pdf" );
            system "cp -r $dir/AS_Result/pdf/*/* $od/Alt_splice/$nam/pdf";

#----------------------------------------- by Simon Young 2014-01-18 ---------------------------------------------------
# it takes too long time.
#            system "gs -q -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=$od/Alt_splice/$nam/$nam.SpliceGrapher.pdf -dBATCH $dir/AS_Result/pdf/*/*.pdf"; # combine all pdf files to one per sample
            system "cp -r $dir/AS_Result/*.xls $od/Alt_splice/$nam/"
              ;    # NOTE: here is no content
        }
    }
    system "cp $Bin/readme/readme/Alt_splice.README.txt $od/Alt_splice";
=cut
	&MKDIR("$od/Alt_splice");
	system "cp -r $id/Alitsplice_Analysis/*.png $od/Alt_splice";
	system "cp -r $id/Alitsplice_Analysis/*.fpkm $od/Alt_splice";
}
if ( -d "$id/DEG_Analysis" ) {
    my @DEG_dir = glob "$id/DEG_Analysis/*";
    my %Anno;
    my %Stat;
	system "cp -r $id/DEG_Analysis/*.png $DEG_od" if (-e "$id/DEG_Analysis/*.png");
        system "cp -r $id/DEG_Analysis/kmean_0/k-means.png $DEG_od" if (-e "$id/DEG_Analysis/kmean_0/k-means.png");
	system "cp $id/DEG_Analysis/All_gene_counts.list $DEG_od" if (-f "$id/DEG_Analysis/All_gene_counts.list");
	system "cp $id/DEG_Analysis/All_gene_fpkm.list $DEG_od" if (-f "$id/DEG_Analysis/All_gene_fpkm.list");
    system "perl $Bin/anno_stat.pl -i $id/DEG_Analysis -anno all -od $DEG_od";
	system "cp -r $id/DEG_Analysis/Anno_Cluster/ $DEG_od " if (-d "$id/DEG_Analysis/Anno_Cluster/");
	system "rm $DEG_od/Anno_Cluster/*_File_*";
	system "rm $DEG_od/Anno_Cluster/*pdf ";
	system "rm -r $DEG_od/Anno_Cluster/work_sh" if (-d "$DEG_od/Anno_Cluster/work_sh");
	foreach my $dir (@DEG_dir) {
        if ( -d $dir ) {
            $dir =~ m/.*\/(\S+)/;
            my $nam = $1;
            &MKDIR("$DEG_od/$nam");
            if ( $dir =~ /\/All_DEG$/ ) {
                system "cp -r $dir/* $DEG_od/$nam/";
            }
            elsif ( $dir =~ /\/work_sh$/ ) {
                `rm -r $DEG_od/$nam`;
            }
			elsif ($dir=~ /\/kmean_0$/) {
				system "cp -r $dir/ $DEG_od";
			}
            elsif ( $dir =~ /\/density$/ ) {
                system "cp -r $dir/*.png $od/geneExpression/";
                system "cp -r $dir/*.cor $od/geneExpression/";
                system "cp -r $dir/cor_plot $od/geneExpression/";
                #system "rm -r $DEG_od/density/";
            }
            elsif ( $dir =~ /_vs_/ ) {
				my $nam=basename($dir);
				&MKDIR("$DEG_od/$nam");
                system "cp -r $dir/Anno_enrichment/* $DEG_od/$nam/";
                system "cp -r $dir/DEG_Cluster $DEG_od/$nam";
                system "cp -r $dir/*.DEG_final.xls $DEG_od/$nam";
				system "cp -r $dir/*.DEG.CytoscapeInput.txt $DEG_od/$nam";
                system "cp -r $dir/*.png $DEG_od/$nam" unless $dir =~ /.*_vs_.*_vs_.*/;
				system "rm $DEG_od/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat" if (-f "$DEG_od/$nam/pathway/kegg_enrichment/*.KEGG.tree.stat");
                system "rm $DEG_od/$nam/DEG_Cluster/hierarchical/change0" if ( -f "$DEG_od/$nam/DEG_Cluster/hierarchical/change0" );
                system "rm $DEG_od/$nam/Graph/topGO.*";

                if ( -d "$dir/Graphic" ) {
                    system "cp $dir/Graphic/*.png $DEG_od/$nam";
                    system "cp $dir/Graphic/*.pdf $DEG_od/$nam";
                }
            }
        }
    }
    
  
	my $tmp=0;
	open(OUT,">$DEG_od/DEG.anno.stat_add");
	while(<IN>){
		if($_=~/#DEG Set/){
			last;
		}
		else{
			if($tmp==0){
				print OUT "#DEG Set\tAnnotated\tCOG\tGO\tKEGG\tSwiss-prot\tTrEMBL\tnr\tnt\n";
			}
			print OUT "$_";
			$tmp++;
			
		}
	}
	close OUT;
	close IN;
	unless (-z "$DEG_od/DEG.anno.stat_add") { ###### 删除中间文件DEG.anno.stat_add
		`cp "$DEG_od/DEG.anno.stat_add" "$DEG_od/DEG.anno.stat"`;
		`rm "$DEG_od/DEG.anno.stat_add"`;
	}
    
    system "cp $Bin/readme/readme/DEG_Analysis.README.txt $DEG_od";
}

if ( -d "$id/Basic_Analysis/geneExpression/final_track/" ) {
    system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene_final.filtered.gff $od/NewGene";
    system "cp $id/Basic_Analysis/geneExpression/final_track/*.newGene.longest_transcript.fa $od/NewGene";
}

if ( -d "$id/Anno_Integrate/New_Anno/Result/" ) {
    system "cp -r $id/Anno_Integrate/New_Anno/Result/* $od/NewGene/NewGene_Anno";
	system "cp -r $id/Anno_Integrate/Allgene_Anno/02.gene-annotation/*_Unigene.fa $od/Allgene/Allgene_Anno";
    system "rm -r $od/NewGene/NewGene_Anno/Blast2go" if ( -d "$od/NewGene/NewGene_Anno/Blast2go" );
    system "rm $od/NewGene/NewGene_Anno/*.annot"; system "rm $od/NewGene/NewGene_Anno/*.svg.list";
	my $function_anno = "$od/NewGene/NewGene_Anno/Function_Annotation.stat.xls";
	my $tmp = "$od/NewGene/NewGene_Anno/tmp";
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
#	system "rm $function_anno";
	system "mv $tmp $function_anno.2";
    system "cp $Bin/readme/readme/NewGene.README.txt $od/NewGene";
}
if ( -d "$id/Lnc_Diff_Analysis" ) {
    my @Lnc_Diff_dir = glob "$id/Lnc_Diff_Analysis/*";
    my %Anno;
    my %Stat;
	system "perl $Bin/anno_stat.pl -i $id/Lnc_Diff_Analysis -anno cis -od $Lnc_Diff_od";
	system "perl $Bin/anno_stat.pl -i $id/Lnc_Diff_Analysis -anno trans -od $Lnc_Diff_od";
	print "perl $Bin/anno_stat.pl -i $id/Lnc_Diff_Analysis -anno cis -od $Lnc_Diff_od";
	system "cp -r $id/Lnc_Diff_Analysis/Cis_Anno_Cluster/ $Lnc_Diff_od " if (-d "$id/Lnc_Diff_Analysis/Cis_Anno_Cluster/");
	system "cp -r $id/Lnc_Diff_Analysis/Trans_Anno_Cluster/ $Lnc_Diff_od " if (-d "$id/Lnc_Diff_Analysis/Trans_Anno_Cluster/");
	system "rm $Lnc_Diff_od/*Anno_Cluster/*_File_*";
	system "rm $Lnc_Diff_od/*Anno_Cluster/*pdf ";
	system "rm -r $Lnc_Diff_od/Cis_Anno_Cluster/work_sh" if (-d "$Lnc_Diff_od/Cis_Anno_Cluster/work_sh");
    system "rm -r $Lnc_Diff_od/Trans_Anno_Cluster/work_sh" if (-d "$Lnc_Diff_od/Trans_Anno_Cluster/work_sh");
	system "cp -r $id/Lnc_Diff_Analysis/kmean_0/k-means.png $Lnc_Diff_od " if (-e "$id/Lnc_Diff_Analysis/kmean_0/k-means.png");
    foreach my $dir (@Lnc_Diff_dir) {
        if ( -f $dir ) {
            next if $dir =~ /\.svg$/;
            next if $dir =~ /\.log$/;
            system "cp $dir $Lnc_Diff_od";
            system "mv $Lnc_Diff_od/All_gene_expression.list $Lnc_Diff_od/All_trans_fpkm.list" if ( -f "$Lnc_Diff_od/All_gene_expression.list" );
        }
        if ( -d $dir ) {
            $dir =~ m/.*\/(\S+)/;
            my $nam = $1;
            &MKDIR("$Lnc_Diff_od/$nam");
            if ( $dir =~ /\/All_DEG$/ ) {
                system "cp -r $dir/ $Lnc_Diff_od/";
            }
			elsif($dir=~/\/kmean_0$/){
				system "cp -r $dir/ $Lnc_Diff_od/";
				#system "rm -r $Lnc_Diff_od/kmean_0";
			}
            elsif ( $dir =~ /\/work_sh$/ ) {
                `rm -r $Lnc_Diff_od/$nam`;
            }
            elsif ( $dir =~ /\/density$/ ) {
                system "cp -r $dir/*.png $od/LncRNA/LncExpression/";
                system "cp -r $dir/*.cor $od/LncRNA/LncExpression/";
                system "cp -r $dir/cor_plot $od/LncRNA/LncExpression/";
                system "rm -r $Lnc_Diff_od/density/";
            }
            elsif ( $dir =~ /_vs_/ ) {
                
                system "cp -r $dir/DEG_Cluster $Lnc_Diff_od/$nam";
                system "cp -r $dir/*.DEG_final.xls $Lnc_Diff_od/$nam/";
                system "cp -r $dir/*.png $Lnc_Diff_od/$nam" unless $dir =~ /.*_vs_.*_vs_.*/;
                system "cp -r $dir/Cis_Anno_enrichment/ $Lnc_Diff_od/$nam/";
                system "cp -r $dir/Trans_Anno_enrichment/ $Lnc_Diff_od/$nam/";
				  system "rm $Lnc_Diff_od/$nam/DEG_Cluster/hierarchical/change0"
                  if (
                    -f "$Lnc_Diff_od/$nam/DEG_Cluster/hierarchical/change0" );
                system "rm $Lnc_Diff_od/$nam/Graph/topGO.*";

                if ( -d "$dir/Graphic" ) {
                    system "cp $dir/Graphic/*.png $Lnc_Diff_od/$nam";
                    system "cp $dir/Graphic/*.pdf $Lnc_Diff_od/$nam";
                }
            }
        }
    }
    system "cp $Bin/readme/readme/DEG_Analysis.README.txt $Lnc_Diff_od";
}

if ( -d "$id/SNP_Analysis" ) {
    my $SNP_dir = "$id/SNP_Analysis";
    system "cp -r $SNP_dir/stat $od/SNP_Analysis";
	system "mv -f $od/SNP_Analysis/stat $od/SNP_Analysis/SNP";
}

if ( -d "$id/Gene_Structure_Optimize" ) {
    system "cp -r $id/Gene_Structure_Optimize $od";
    system "cp $Bin/readme/readme/Gene_Structure_Optimize.README.txt $od/Gene_Structure_Optimize";
}

if ( -d "$id/Compare_analysis" ) {
    system "cp -r $id/Compare_analysis/*RNA.len* $od/Compare_Analysis/Length" ; #改动的地方不需要改---------------------------------------------------------------------------------------------------------------------------------------
    system "cp -r $id/Compare_analysis/chr_num_stat $od/Compare_Analysis/";
	system "cp -r $id/Compare_analysis/lncRNA_mRNA_DEG/*png $od/Compare_Analysis/lncRNA_mRNA_DEG";
    system "cp -r $id/Compare_analysis/*exon* $od/Compare_Analysis/Exon";
    system "cp $id/Compare_analysis/*orf* $od/Compare_Analysis/Orf_len";
    system "cp $id/Compare_analysis/*FPKM* $od/Compare_Analysis/Expression/";
    system "cp $id/Compare_analysis/*fpkm* $od/Compare_Analysis/Expression/";
	system "cp -r $id/Compare_analysis/*isoform.png $od/Compare_Analysis/Isoform";
	system "cp -r $id/Compare_analysis/*RNA_isoform.txt $od/Compare_Analysis/Isoform";	
}

if (-d "$id/DEU_analysis") {
	&MKDIR("$od/DEU_analysis");
	my @DEU_dir = glob"$id/DEU_analysis/*";
	foreach my $dir (@DEU_dir) {
		if ($dir =~ /_vs_/) {
			my $name = basename($dir);
			my $sample=$name;
			$sample=~/\/(\w+_vs_\w+)\//;
			&MKDIR("$od/DEU_analysis/$name/DEXSeqReport");
			&MKDIR("$od/DEU_analysis/$name/DEXSeqReport/files");
			system "cp -r $dir $od/DEU_analysis";
			system "rm $od/DEU_analysis/$name/DEU_Result_All.xls" if (-f "$od/DEU_analysis/$name/DEU_Result_All.xls");
			system "rm -r $od/DEU_analysis/$name/DEXSeqReport/files/*.svg" if (-f "$od/DEU_analysis/$name/DEXSeqReport/files/*.svg");
			system "mv $od/DEU_analysis/$name/DEU_Result_Final.xls $od/DEU_analysis/$name/$sample.DEU_Result_Final.xls";
			#system "cp -r $dir/DEXSeqReport/testForDEU.html $od/DEU_analysis/$name/DEXSeqReport";
			#system "cp -r $dir/DEXSeqReport/files/*.html $od/DEU_analysis/$name/DEXSeqReport/files";
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
system "cp $Bin/readme/readme.txt $od/";

################Time
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir = `pwd`;
    chomp($cur_dir);
    my ($in) = @_;
    my $return = "";

    if ( -f $in ) {
        my $dir  = dirname($in);
        my $file = basename($in);
        chdir $dir;
        $dir = `pwd`;
        chomp $dir;
        $return = "$dir/$file";
    }
    elsif ( -d $in ) {
        chdir $in;
        $return = `pwd`;
        chomp $return;
    }
    else {
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}

sub para_load {
    my ( $file, $para ) = @_;
    open IN, $file || die "$!";
    while (<IN>) {
        chomp;
        s/\r+//g;
        next if ( /^$/ || /^\#/ );
        my ( $key, $info ) = ( split /\s+/, $_ )[ 0, 1 ];
        if ( !$key ) { print "$_\n"; die; }
        $para->{$key} = $info;
    }
    close IN;
}

sub sub_format_datetime {    #Time calculation subroutine
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
    $wday = $yday = $isdst = 0;
    sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

sub Runtime {                # &Runtime($BEGIN);
    my ($t1) = @_;
    my $t = time() - $t1;
    print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR {                  # &MKDIR($out_dir);
    my ($dir) = @_;
    rmdir($dir) if ( -d $dir );
    mkdir($dir) if ( !-d $dir );
}

sub help {
    print << "	Usage End.";
Description: Extract lncRNA Analysis Reaults for Html Process;
    version:$ver
      Usage:
        --id  <STR>   input dir, analysis output directory   force
        --od  <STR>   result output dir                      force
        --h           help
	Usage End.
    exit;
}
