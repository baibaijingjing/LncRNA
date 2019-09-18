#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use newPerlBase;
my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 

use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
my $version="1.4.0";
my $BEGIN_TIME=time();

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($cfg,$odir);
my ($nr,$nt,$Swissprot,$TrEMBL,$GO,$Cog,$Kog,$Pfam,$Kegg,$all);
my ($Query,$Blastp,$Blast_cpu,$Blast_e,$blast_cut,$hmmscan_cpu);
my ($sh_dir,$Result_dir,$Tab_dir,$div_dir);
my $step;

GetOptions(
		"cfg:s"=>\$cfg,
		"qnuc:s"=>\$Query,
		"all"=>\$all,
		"nr"=>\$nr,
		"nt"=>\$nt,
		"swissprot"=>\$Swissprot,
		"trembl"=>\$TrEMBL,
		"cog"=>\$Cog,
		"kog"=>\$Kog,
		"pfam"=>\$Pfam,
		"kegg"=>\$Kegg,
		"GO"=>\$GO,
		"od:s"=>\$odir,
        "step=i" =>\$step,
		"help"=>\&USAGE
	) or &USAGE;
&USAGE unless ($cfg and $odir and $Query) ;

# ------------------------------------------------------------------
# load arguments form config file & init parameters.
# ------------------------------------------------------------------
&MKDIR($odir);
$odir=&ABSOLUTE_DIR($odir);
$cfg=&ABSOLUTE_DIR($cfg);

# load arguments form config file
my %CFG;
&LOAD_PARA($cfg,\%CFG);

#$Query=$CFG{mRNA};
$Blast_cpu=$CFG{blast_cpu};
$hmmscan_cpu=$CFG{hmmscan_cpu};
$Blast_e=$CFG{blast_e};
$blast_cut=$CFG{blast_cut};

# init parameters
my $Q_name=basename $Query;
$step ||= 1;

# deal with illegal parameters
if (( !defined $nr)&& ( !defined $nt) && ( !defined $Cog) && ( !defined $Kog) && ( !defined $Pfam) && ( !defined $Swissprot) && ( !defined $TrEMBL) && ( !defined $Kegg) && (! defined $GO) && (!defined $all) ) {
	print BOLD RED "Error: You must choose database to align.\n";
	exit;
}
if (-f $Query) {
    print "[".&date_time_format(time())."] make sure your query FASTA file: $Query\n";
} else {
    print BOLD RED "Error: Can't find your query file to annotaion.\n";
    exit;
}
unless ($step>=1 and $step<=4) {
	print BOLD RED "Error: start step must belong to {1,2,3,4}.\n";
	exit;
}

#######################################################################################
print "\n[".&date_time_format($BEGIN_TIME)."] $Script start ...\n";
#######################################################################################
# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
# make basic directories
&MKDIR("$odir/mid");
my $mid_dir="$odir/mid";
&MKDIR("$odir/Result");
$Result_dir="$odir/Result";
&MKDIR("$odir/work_sh");
$sh_dir="$odir/work_sh";
#system("cp -f /home/xugl/.kobasrc ~/");
#system("export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/R/3.1.1/lib64/R/lib/:\$LD_LIBRARY_PATH");
# ------------------------------------------------------------------
# cut query FASTA file
# ------------------------------------------------------------------
if ($step == 1) {
    print "[".&date_time_format(time())."] cut query FASTA file:\n";
    print "[".&date_time_format(time())."] perl $Bin/bin/cut_fa_file_to_dir.pl -mRNA $Query -od $mid_dir -cut $blast_cut\n\n";
    `perl $Bin/bin/cut_fa_file_to_dir.pl -mRNA $Query -od $mid_dir -cut $blast_cut`;
    $step = 2;
}

# ------------------------------------------------------------------
# align with databases and alignment result format.
# ------------------------------------------------------------------
my $database_num=0;
if ($step == 2) {
    print "[".&date_time_format(time())."] now align with databases:\n";
    my $SH="$odir/work_sh/Align_Database.sh";
    open (SH,">$SH") or die $!;
    if (defined $Cog||defined $all) {
		$database_num++ if (defined $Cog);
        print "       COG: $CFG{Cog} \n";
    #	print "Anno Cog Database\nperl $Bin/bin/Gene_Cog_Anno.pl -id $mid_dir -Database $CFG{Cog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_Cog_Anno.pl -id $mid_dir -Database $CFG{Cog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2}\n";
    }
    if (defined $Kog||defined $all) {
		$database_num++ if (defined $Kog);
        print "       KOG: $CFG{Kog} \n";
    #	print "Anno Kog Database\nperl $Bin/bin/Gene_Kog_Anno.pl -id $mid_dir -Database $CFG{Kog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_Kog_Anno.pl -id $mid_dir -Database $CFG{Kog} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2} \n";
    }
    if (defined $Pfam||defined $all) {
		$database_num++ if (defined $Pfam);
        print "      Pfam: $CFG{Pfam} \n";
    #	print "Anno Kog Database\nperl $Bin/bin/Gene_Pfam_Anno.pl -id $mid_dir -Database $CFG{Pfam} --od $odir --cpu $hmmscan_cpu \n";
        print SH "perl $Bin/bin/Gene_Pfam_Anno.pl --idir $mid_dir --pfam $CFG{Pfam} --odir $odir --cpu $hmmscan_cpu -queue  $CFG{Queue_type2}\n";
    }
    if (defined $Swissprot||defined $all) {
		$database_num++ if (defined $Swissprot);
        print "Swiss-Prot: $CFG{Swissprot} \n";
    #	print "Anno Swissprot Database\nperl $Bin/bin/Gene_Swissprot_Anno.pl -id $mid_dir -Database $CFG{Swissprot} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_Swissprot_Anno.pl -id $mid_dir -Database $CFG{Swissprot} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2}\n";
    }
    if (defined $Kegg||defined $all) {
		$database_num++ if (defined $Kegg);
        print "      KEGG: $CFG{Kegg} \n";
    #	print "Anno Kegg Database\nperl $Bin/bin/Gene_KEGG_Anno.pl -id $mid_dir -Database $CFG{Kegg} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_KEGG_Anno.pl -id $mid_dir -Database $CFG{Kegg} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2} \n";
    }
    if (defined $TrEMBL||defined $all) {
		$database_num++ if (defined $TrEMBL);
        print "    TrEMBL: $CFG{TrEMBL} \n";
    #	print "Now Anno TrEMBL Database\nperl $Bin/bin/Gene_TrEMBL_Anno.pl -id $mid_dir -Database $CFG{TrEMBL} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_TrEMBL_Anno.pl -id $mid_dir -Database $CFG{TrEMBL} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2} \n";
    }
    if (defined $nr||defined $all) {
        print "        nr: $CFG{nr} \n";
    #	print "Now Anno Nr Database\nperl $Bin/bin/Gene_Nr_Anno.pl -id $mid_dir -Database $CFG{nr} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
     print SH "perl $Bin/bin/Gene_Nr_Anno.pl -id $mid_dir -Database $CFG{nr} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2}\n";
		unless (defined $nt || defined $all ){
			print SH "\n" x ($database_num-1);
		}

    }
	if (defined $nt||defined $all) {
        print "        nt: $CFG{nt} \n";
    #	print "Now Anno Nt Database\nperl $Bin/bin/Gene_Nt_Anno.pl -id $mid_dir -Database $CFG{nt} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e \n";
        print SH "perl $Bin/bin/Gene_Nt_Anno.pl -id $mid_dir -Database $CFG{nt} -od $odir -Blast_cpu $Blast_cpu -Blast_e $Blast_e -queue  $CFG{Queue_type2} \n";
		unless (defined $nr || defined $all ) {
			print SH "\n" x ($database_num-1);
		}
		if (defined $nr and !defined $all) {
			print SH "\n" x $database_num;
		}
		if (defined $all) {
			print SH "\n" x 6;
		}
    }

    close SH;
    print "\n[".&date_time_format(time())."] start to qsub $SH.\n";
    &qsubOrDie($SH,$CFG{Queue_type2},8,"10G");
    $step = 3;
}

# ------------------------------------------------------------------
# Blast2GO based on nr database blast result
# ------------------------------------------------------------------
if ($step == 3 && ($GO||$all)) {
    print "[".&date_time_format(time())."] start Blast2GO based on nr blast result:\n";
    print "[".&date_time_format(time())."] perl $Bin/bin/Blast2GO.pl -id $odir/Nr_Dir/ -od $Result_dir -sh_dir $sh_dir -k $Q_name  -queue  $CFG{Queue_type2}  \n\n";
    `perl $Bin/bin/Blast2GO.pl -id $odir/Nr_Dir/ -od $Result_dir -sh_dir $sh_dir -k $Q_name -queue  $CFG{Queue_type2} `;
    $step = 4;
}

# ------------------------------------------------------------------
# annotation result tackle, plot and stat
# ------------------------------------------------------------------
if ($step == 4) {
    ## result integrate
    if((defined $Cog || defined $Kog || defined $nr || defined $nt || defined $Swissprot || defined $TrEMBL || defined $Kegg) || defined $all){
        print "[".&date_time_format(time())."] Anno Integrate:\n";
        print "[".&date_time_format(time())."] perl $Bin/bin/anno_integrate.pl -gene $Query -id $Result_dir -od $Result_dir -key $Q_name \n\n";
        `perl $Bin/bin/anno_integrate.pl -gene $Query -id $Result_dir -od $Result_dir -key $Q_name `;
    }

    ## GO histogram plot
    if ($GO||$all) {
        print "[".&date_time_format(time())."] Draw GO Graph:\n";
        print "[".&date_time_format(time())."] perl $Bin/bin/draw_GO_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir \n\n";
        `perl $Bin/bin/draw_GO_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir `;
    }

    ## COG and KOG histogram plot
    if ($Cog||$Kog||$all) {
        print "[".&date_time_format(time())."] Draw COG & KOG Graph:\n";
        print "[".&date_time_format(time())."] perl $Bin/bin/draw_COG_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir \n\n" if ($Cog||$all);
        `perl $Bin/bin/draw_COG_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir ` if ($Cog||$all);
        print "[".&date_time_format(time())."] perl $Bin/bin/draw_KOG_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir \n\n" if ($Kog||$all);
        `perl $Bin/bin/draw_KOG_graph.pl -i $Result_dir/All_Database_annotation.xls -k $Q_name -od $Result_dir ` if ($Kog||$all);
    }

    ## nr pie plot
    if ($nr||$all) {
        `perl $Bin/bin/util/nr_pie_stat.pl -i $Result_dir/$Q_name.nr.anno.txt -o $Result_dir/$Q_name.nr.lib.stat -m 10 `;
    }
}

#######################################################################################
print "\n[".&date_time_format(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub LOAD_PARA {
	my $para_file= shift;
	my $para= shift;

    # switch node in which store databases to search.
    # add node to configs file when needed.
    my %db_node;
    my $db_idx = 0;
    open(NODE,"$Bin/../../../config/db_file.cfg") or die;
    while (<NODE>) {
        next unless (/^db_node/);
        my $node = (split /\s+/,$_)[1];
        $db_node{$db_idx++} = $node;
    }
    close NODE;
    my $random = int(rand(scalar keys %db_node));

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
        $para_value =~s/\/share|\/lustre/$db_node{$random}/ unless ($para_key eq 'mRNA');
        #$para_value =~s/\/lustre/$db_node{$random}/ unless ($para_key eq 'mRNA');
        $para->{$para_key} = $para_value;

		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
        if ($para_key eq "nr" or $para_key eq "Cog" or $para_key eq "Kog" or $para_key eq "Pfam" or $para_key eq "Swissprot" or $para_key eq "Kegg" or $para_key eq "nt" or $para_key eq "TrEMBL") {
            unless (-f $para_value) {
                warn "Error: Can't find $para_key database: $para_value.\n";
                $error_status = 1;
            }
			die "Kegg database is wrong ,should include kobas data\n " if ($para_key eq 'Kegg' and $para_value!~/kobas/);
        }
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub LOAD_SEQ {
	my ($fa,$info) = @_;

	open IN,"$fa" || die $!;
	$/='>';
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/);
		my ($head,$seq)=split/\n+/,$_,2;
		my $id=(split/\s+/,$head)[0];
		$info->{$id}=$seq;
	}
	$/="\n";
	close IN;
}



sub cut_str {
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
}

sub parse_config { # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+//;s/\s+$//;s/\r$//;
		next if(/$/ or /\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}

sub Runtime { # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub date_time_format {#Time calculation subroutine
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	`mkdir -p $dir ` if(!-d $dir);
}

sub ABSOLUTE_DIR { #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if (-f $in) {
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	} elsif(-d $in) {
		chdir $in;$return=`pwd`;chomp $return;
	} else {
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}

sub USAGE {
    my $usage =<<"_USAGE_";
#-------------------------------------------------------------------------------------------------
    Program: $Script
    Version: $version
     Writer: Meng Fei <mengf\@biomarker.com.cn>
       Data: 2012-00-00
   Modifier: Simon Young <simonyoung8824\@gmail.com>
       Data: 2014-09-28
Description: alignment with fuctional database and annotate genes. v1.1.0-based,
             1) add KOG and Pfam database annotation, additionally predict CDS;
             2) deal with disturbed work_sh output and excessive qsub;
             3) fix several Blast2GO empty annotation bug;
             4) fix excel overflow bug by spliting files with more than 65,000 lines.

      Usage:
             --all             search against all database

               --nr            search against nr database
               --swissprot     search against Swiss-Prot database
               --cog           search against COG database
               --kog           search against KOG database for eukaryote
               --kegg          search against KEGG database
               --pfam          search against Pfam database
               --GO            run Blast2GO start GO Annotation (--nr is required)
               --nt            search against nt database
               --trembl        search against TrEMBL database

             --step            program start step, (1|2|3|4)
               1               start from the beginning, default
               2               start from alignment with databases and alignment result format
               3               start from Blast2GO
               4               start from result tackle, plot and stat

             --cfg             database config
			 --qnuc            query (fa file) 
             --od              output dir
             --help            show help document

	Example:
	  perl $Script --all --qnuc new_gene.fa --cfg Gene_Func_Anno_Pipline.cfg --od ./

#-------------------------------------------------------------------------------------------------
_USAGE_
    print $usage;
    exit;
}
