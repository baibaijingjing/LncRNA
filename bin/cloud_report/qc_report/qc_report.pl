#!usr/bin/perl -w
use strict;

use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

my ( $cfg, $od, $data_type, $sep, $cut,$in );

GetOptions(
    "cfg:s"   => \$cfg,
    "od:s"    => \$od,
    "types:s" => \$data_type,
    "a:i"     => \$sep,
    "b:i"     => \$cut,
    "in:s"=>\$in,
) || die &USAGE;

die &USAGE unless ( $cfg && $od  && $in);

$cfg = abs_path($cfg);
system "mkdir -p $od" unless ( -d $od );
$cut       ||= 15;
$sep       ||= 2;
$data_type ||= 'pe';

my %para;
read_para();

my $library_qc = "/share/nas1/basecall/Config_dir/";

my $gff                = $para{'Ref_ann'};
my $project_id         = $para{'Project_id'};

my $data_assess     ="$in/Data_Assess";
my $deg_dir            = "$in/DEG_Analysis";
my $map_stat_dir       = "$in/Basic_Analysis/Tophat_Cufflinks/Map_Stat/";
my $altSplice_dir      = "$in/Alitsplice_Analysis/";
my $geneExpression_dir = "$in/Basic_Analysis/geneExpression/";
my $new_gene_gff       = "$in/Basic_Analysis/geneExpression/final_track/$para{'Project_key'}.newGene_final.filtered.gff";

my $lnc_deg_dir="$in/Lnc_Diff_Analysis";
my $lnc_pre_dir="$in/Basic_Analysis/Tophat_Cufflinks/Lnc_filter"; 

my $new_gene_num;
read_newGene();
my $known_gene_num;
my $total_gene_num;
read_gff();

my $GC_div;

my %library_info;
read_library();
my %data;
data_assess();
my ($deg_num_deviation,$all_deg_num);
read_deg();
my ($lnc_deg_num_deviation,$lnc_deg_all_num);
read_lnc_deg();
my ($lncRNA_ratio,$lnc_exp_ratio,$lnc_deg_ratio);
read_lnc_predict();
my ( %map_ratio, $map_ratio_deviation );
read_mapStat();
my %insert_around_per;
my %insert_mode;
read_insertSize();
my %randCheck;
read_randCheck();
my %as_ratio;
read_alterSplice();
my %sensitivity;
read_geneExpression();

open( OUT, ">$od/$project_id.QC.stat.xls" ) or die;
print OUT join( "\t",
    '#ProjectType', 'ProjectID',        'SampleID', 'SpeciesType',  'RefSeq?', 'RequiredData',
    'ObtainedData', 'RequiredQuality',  'ObtainQuality','UseData',    'Qubit',        'LibraryQuality',   'BaseIso',    'BaseWave',
    'NT%(contaminate)', 'NT%(irrelative)',    'rRNA%',        'NEBorNOT',
    'InsertAround%',   'GCcontent',  'MaxDepthPer',  'Mapped%',      'Mapped%Range',
    'N50Length',    '<300nt%',          'Tissue',    'assemble_div',
    'LncRNA_ratio','DEG_lncRNA','DEG_lncRNA_ratio','DEGnum', 'TotalGene', 'AlterSpliceGene%',
    'Saturation',   'NovelGeneNum',      'IntraSq.Rmin', 'InterSq.Rmax',     'Deadline',    'FinishDate' )
  . "\n";

foreach my $sample_id ( sort keys %map_ratio ) {
    my ( $id, $qubit, $quality ) =( exists $library_info{$sample_id} ) ? @{ $library_info{$sample_id} } : ( '?', '?', '?' );
    my (
         $obtained_base, $GC_content,
        $Q_20,          $Q_30,          $base_iso,
        $base_wave,     $rRNA_per  ) = @{ $data{$sample_id} };
    my $GC_div_              = $GC_div;
    my $map_ratio            = $map_ratio{$sample_id};
    my $map_ratio_deviation_ = $map_ratio_deviation;
    my $insert_around_per    = $insert_around_per{$sample_id};
    my $insert_mode          = $insert_mode{$sample_id};
    my $randCheck            = $randCheck{$sample_id};
    my $sensitivity          = $sensitivity{$sample_id};
    my $deg_num_deviation_   = &format_figure($deg_num_deviation);
    my $as_ratio             = $as_ratio{$sample_id};
    my $new_gene_num_        = &format_figure($new_gene_num);
    my $FinishDate           = &GetDate;

    print OUT join( "\t",
        'Lnc',               $project_id,  $sample_id,        '?',                   'Ref',        "10G",
        $obtained_base,        '85.00%',     $Q_30,$obtained_base,        '--',                '--',     $base_iso,        $base_wave,
        '?',          '?',        $rRNA_per,             '?',
        $insert_around_per, $GC_content,     $randCheck,   $map_ratio,        $map_ratio_deviation_,
        '-',          '-',        '-',     '-',
        $lncRNA_ratio,$lnc_deg_all_num,$lnc_deg_ratio, $all_deg_num,$total_gene_num,        $as_ratio,
        $sensitivity, $new_gene_num_,        '?',                   '?',          '?',   $FinishDate )
      . "\n";
}

close OUT;

if ( -d "$data_assess/PNG" ) {
    system "cp -r $data_assess/PNG $od";
    system "cp $geneExpression_dir/Total.*.png $geneExpression_dir/*.insertSize.r.png $od/PNG";
}

sub data_assess {
    my @data_statistics = glob ("$data_assess/*.stat");
    my @GC_content;
    foreach my $data_statistics (@data_statistics) {
        
        my $stat_name=basename($data_statistics);
        next if ($stat_name=~/AllSample/);
        $stat_name=~/(\S+).stat/;
        my $sample_id=$1;
        open( DATA, $data_statistics ) || die "Open $data_statistics failed!\n";
        <DATA>;
        <DATA>;
        <DATA>;
        while (<DATA>) {
            chomp;
            my @data          = split /\t+/;
            #my $sample_id     = $data[0];
            my $obtained_read = $data[1];
            my $obtained_base = &format_figure( $data[2] );
            my $GC_content    = $data[3];
    
            my $Q_20 = $data[5] . '%';
            my $Q_30 = $data[7] . '%';
            push @{ $data{$sample_id} },
              (
                $obtained_read,
                sprintf( "%.2f", $GC_content ) . '%',
                $Q_20, $Q_30
              );
    
            push @GC_content, $GC_content;
        }
        close DATA;
    }
    

    my ( $min_GC_content, $max_GC_content ) =( sort { $a <=> $b } @GC_content )[ 0, -1 ];
    $GC_div = $max_GC_content - $min_GC_content;
    $GC_div = &percent( $GC_div / 100 );

    my @base_assess = glob("$data_assess/*.acgtn");
    foreach my $base_assess (@base_assess) {
        my ($sample_id) = basename($base_assess) =~ /^(\S+)\.acgtn$/;
        my ( $base_iso, $base_wave ) =
          base_assess( $base_assess, $data_type, $sep, $cut );
        $base_iso  = percent( $base_iso / 100 );
        $base_wave = percent( $base_wave / 100 );
        push @{ $data{$sample_id} }, ( $base_iso, $base_wave );
    }
    if (-f "$data_assess/AllSample.data.stat") {
        open( ASS, "$data_assess/AllSample.data.stat" ) or die $!;
        while (<ASS>) {
            chomp;
            next if ( /^\s+$/ || $. == 1 || /^#/ );
            my @col = ( split /\t/ );
            #die "Note your file: $data_assess/AllSample.data.stat \n" if ( @col != 11 );
            my $sample_id = $col[0];
            push (@{ $data{$sample_id} },     &format_figure( $col[3] ) . '%' );
        }
        close ASS;
    }
    
    

}

sub read_library {
    opendir( LIBRARY, $library_qc ) || die "Opendir $library_qc failed!\n";
    my $flag;
    while ( my $file = readdir LIBRARY ) {
        next unless $file =~ /.*\.txt$/;
        $file = "$library_qc/$file";
        open( FILE, $file ) || die "Open file $file failed!\n";
        while (<FILE>) {
            chomp;
            s/\r//;
            next unless ( /$project_id/ && /Lnc-Str/ );
            my @data = split /\t+/, $_;

#			my ($sample_id) = $data[3] =~/.*-(.*)-.*/;      # error when H25-T07-2-I
#----------------------------- by Simon Young 2015-01-08 -----------------------------------
            my ($sample_id) = $data[3] =~ /^[a-zA-Z0-9]+-([a-zA-Z0-9]+)-.*/;
            my $qubit = sprintf "%.2f", $data[14];
            $library_info{$sample_id} = [ $data[3], $qubit, $data[19] ];

            $flag++;
        }
        close FILE;

    #		last if $flag;          # error when samples' info in more than one file.
    }
    closedir LIBRARY;

    if ( scalar( keys %library_info ) == 0 ) {
        print "This project library infomation is lost in library_qc dir!\n";
        #die;
    }
}

sub read_insertSize {
    my @insertSize_files = glob("$map_stat_dir/*.insertSize.list");
    foreach my $insertSize_file (@insertSize_files) {
        open( IN, $insertSize_file ) || die "Open $insertSize_file failed!\n";
        my ($sample_id) =
          basename($insertSize_file) =~ /^(\S+)\.insertSize\.list/;
        my %insert_info;
        while (<IN>) {
            chomp;
            next unless /^(\d+):(\d+)$/;
            $insert_info{$1} = $2;
        }
        close IN;

        my @len =
          sort { $insert_info{$a} <=> $insert_info{$b} } ( keys %insert_info );
        $insert_mode{$sample_id} = $len[-1];
        my ( $count, $sum );
        for ( my $i = $len[-1] - 75 ; $i <= $len[-1] + 75 ; $i++ ) {
            $insert_info{$i} ||= 0;
            $count += $insert_info{$i};
        }
        foreach my $len ( keys %insert_info ) {
            $sum += $insert_info{$len};
        }
        $insert_around_per{$sample_id} = &percent( $count / $sum );
    }
}

sub read_deg {
    my @deg_files = glob("$deg_dir/*vs*/*.DEG_final.xls");
    my @deg_num;
    foreach my $deg_files (@deg_files) {
        my $deg_num = `wc -l $deg_files`;
        chomp $deg_num;
        $deg_num = ( split /\s+/, $deg_num )[0];
        $deg_num--;
        push @deg_num, $deg_num;
    }

    my ( $deg_num_max, $deg_num_min ) =
      ( sort { $a <=> $b } @deg_num )[ -1, 0 ];
    $deg_num_deviation = $deg_num_max - $deg_num_min;
    $all_deg_num=`less $deg_dir/All_DEG/All.DEG_final.xls |wc -l`;
    chomp $all_deg_num;
}
sub read_lnc_predict{
    my $candi_lnc =`grep '>' $lnc_pre_dir/merged_filter.fa|wc -l`;
    chomp $candi_lnc;
    my $lncRNA=`grep '>' $lnc_pre_dir/lnc_filter_final.fa|wc -l`;
    chomp $lncRNA;
    $lncRNA_ratio=( sprintf "%.2f", ( $lncRNA/$candi_lnc)*100)."%";
    $lnc_deg_ratio=( sprintf "%.2f", ( $lnc_deg_all_num/$lncRNA)*100)."%";
    
}
sub read_lnc_deg{
    my @lnc_deg_files=glob ("$lnc_deg_dir/*vs*/*.DEG_final.xls");
    my @lnc_num;
    foreach my $deg_files (@lnc_deg_files) {
        my $deg_num = `wc -l $deg_files`;
        chomp $deg_num;
        $deg_num = ( split /\s+/, $deg_num )[0];
        $deg_num--;
        push @lnc_num, $deg_num;
    }

    my ( $deg_num_max, $deg_num_min ) =      ( sort { $a <=> $b } @lnc_num )[ -1, 0 ];
    $lnc_deg_num_deviation = $deg_num_max - $deg_num_min;
    $lnc_deg_all_num=`less $lnc_deg_dir/All_DEG/All.DEG_final.xls|wc -l`;
    chomp $lnc_deg_all_num;
    $lnc_deg_all_num=$lnc_deg_all_num-1;
    
}

sub read_geneExpression {
    my @geneExpression_files =
      glob("$geneExpression_dir/*\.geneExpression.xls");
    foreach my $geneExpression_file (@geneExpression_files) {
        my ($sample_id) =
          basename($geneExpression_file) =~ /^(\S+)\.geneExpression\.xls/;
        my $detected_num;
        open( IN, $geneExpression_file )
          || die "Open $geneExpression_file failed!\n";
        while (<IN>) {
            chomp;
            next if ( /^#/ || /^\s*$/ );
            my @data = split /\t+/;
            $detected_num++ if $data[2] > 0;
        }
        close IN;

        my $sensitivity =
          ( sprintf "%.2f", ( $detected_num / $total_gene_num * 100 ) ) . "%";
        $sensitivity{$sample_id} = $sensitivity;
    }
}

sub read_mapStat {
    my @mapStat_files = glob("$map_stat_dir/*.mappedStat.xls");
    my @map_ratio;
    foreach my $mapStat_file (@mapStat_files) {
        my ($sample_id) = basename($mapStat_file) =~ /^(\S+)\.mappedStat\.xls/;
        open( IN, $mapStat_file ) || die "Open $mapStat_file failed!\n";
        my @lines = <IN>;
        close IN;

        my $map_ratio = ( split /\s+/, $lines[1] )[-1];
        chomp $map_ratio;
        $map_ratio{$sample_id} = $map_ratio;
        chop $map_ratio;
        push @map_ratio, $map_ratio;
    }
    my ( $map_ratio_max, $map_ratio_min ) = ( sort { $a <=> $b } @map_ratio )[ -1, 0 ];
    $map_ratio_deviation = &format_figure( $map_ratio_max - $map_ratio_min );
    $map_ratio_deviation .= "%";
}

sub read_randCheck {
    my @randCheck_files = glob("$map_stat_dir/*randcheck_per.list");
    foreach my $randCheck_file (@randCheck_files) {
        my ($sample_id) =
          basename($randCheck_file) =~ /^(\S+)\.randcheck_per\.list/;
        $randCheck{$sample_id} = 0;
        open( IN, $randCheck_file ) || die "Open $randCheck_file failed!\n";
        while (<IN>) {
            chomp;
            if (/^[\d\.]+:([\.\d]+)$/) {

                #				$randCheck{$sample_id} += ($1-1)**2;  #DepthVariance
                $randCheck{$sample_id} = ( $randCheck{$sample_id} > $1 ) ? $randCheck{$sample_id}: $1;    #MaxDepthPer
            }
        }
        close IN;
        $randCheck{$sample_id} = &format_figure( $randCheck{$sample_id} ) . '%';
    }
}

sub read_alterSplice {
    my @mapStat_files = glob("$altSplice_dir/*fpkm");
    foreach my $mapStat_file (@mapStat_files) {
        my ($sample_id) =
          basename($mapStat_file) =~ /^(\S+)\.fpkm$/;
        
        my $total_as_number =`less $mapStat_file|cut -f 3|sort|uniq|wc -l`;
        chomp $total_as_number;
        my $total_as_ratio =( sprintf "%.2f", ( $total_as_number / $known_gene_num ) * 100 ). "%";
        $as_ratio{$sample_id} = $total_as_ratio;
    }
}

sub read_newGene {
    $new_gene_num = `less $new_gene_gff|awk '\$3=="gene"'|wc -l`;
    chomp $new_gene_num;
}

sub read_gff {
    $known_gene_num = `less $gff|awk '\$3=="gene"'|wc -l`;
    chomp $known_gene_num;
    $total_gene_num = $known_gene_num + $new_gene_num;   # total = known + novel
}

sub read_para {
    open( IN, $cfg ) || die "Open $cfg failed!\n";
    while (<IN>) {
        chomp;
        next if ( /^#/ || /^\s*$/ );
        my ( $mark, $path ) = split /\s+/;
        $para{$mark} = $path;
    }
    close IN;
}

sub base_assess {
    my ( $fIn, $data_type, $sep, $cut ) = @_;
    my ( %GC, $base_iso, $GC_sep );
    open IN, $fIn or die $!;
    <IN>;
    while (<IN>) {
        my @in = split /\t/, $_;
        $GC{ $in[0] }{A} = $in[1];
        $GC{ $in[0] }{C} = $in[2];
        $GC{ $in[0] }{G} = $in[3];
        $GC{ $in[0] }{T} = $in[4];
    }
    close IN;
    my $sum_AT = 0;
    my $sum_GC;
    my $num = 0;
    if ( $data_type eq 'se' ) {
        foreach my $key ( keys %GC ) {
            if ( $key > $cut ) {
                $sum_AT += abs( $GC{$key}{A} - $GC{$key}{T} );
                $sum_GC += abs( $GC{$key}{G} - $GC{$key}{C} );
                if (   abs( $GC{$key}{A} - $GC{$key}{T} ) > $sep
                    or abs( $GC{$key}{G} - $GC{$key}{C} ) > $sep )
                {
                    $num++;
                }
            }
        }
        $base_iso = &max( $sum_AT, $sum_GC ) / ( 101 - $cut );
        $GC_sep = $num / ( 101 - $cut );
        return ( $base_iso, $GC_sep );
    }
    else {
        foreach my $key ( keys %GC ) {
            if ( $key > $cut and $key <= 101 ) {
                $sum_AT += abs( $GC{$key}{A} - $GC{$key}{T} );
                $sum_GC += abs( $GC{$key}{G} - $GC{$key}{C} );
                if (   ( abs( $GC{$key}{A} - $GC{$key}{T} ) ) >= $sep
                    or ( abs( $GC{$key}{G} - $GC{$key}{C} ) ) >= $sep )
                {
                    $num++;
                }
            }
            elsif ( $key > 101 + $cut ) {
                $sum_AT += abs( $GC{$key}{A} - $GC{$key}{T} );
                $sum_GC += abs( $GC{$key}{G} - $GC{$key}{C} );
                if (   abs( $GC{$key}{A} - $GC{$key}{T} ) > $sep
                    or abs( $GC{$key}{G} - $GC{$key}{C} ) > $sep )
                {
                    $num++;
                }
            }
        }
        my $aa = $sum_AT > $sum_GC ? $sum_AT : $sum_GC;
        $base_iso = $aa /  ( 202 - $cut * 2 );
        $GC_sep   = $num / ( 202 - $cut * 2 );
        return ( $base_iso, $GC_sep );
    }
}

sub percent {
    my $num = shift;
    if ( $num =~ s/\%// ) {
        $num = sprintf "%.2f", $num;
        $num .= '%';
    }
    else {
        $num = sprintf "%.2f", $num * 100;
        $num .= '%';
    }
}

################################################################################################################
#ÕûÊý¸ñÊ½£º3Î»Ò»¸ö¶ººÅ
sub Integer_Three_Digit {    #
    my $interger = shift;
    $interger =~ s/(?<=\d)(?=(\d\d\d)+$)/,/g;
    return $interger;
}

################################################################################################################
#ÕûÊý¸ñÊ½£º3Î»Ò»¸ö¶ººÅ
#Ð¡Êý¸ñÊ½£ºÐ¡ÊýµãºóÁ½Î»
sub format_figure {    #
    my $figure = shift;
    if ( !defined $figure ) {
        die;
    }
    if ( $figure =~ /\./ ) {
        if ( $figure == 100 ) {
            $figure = 100;
        }
        else {
            $figure = sprintf( "%.2f", $figure );
        }
    }
    else {
        $figure = Integer_Three_Digit($figure);
    }
    return $figure;
}

################################################################################################################
sub GetDate {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf( "%4d-%02d-%02d", $year + 1900, $mon + 1, $day );
}

################################################################################################################
sub USAGE {    #
    my $usage = <<"USAGE";
#---------------------------------------------------------------------------------------------------
Usage:
	Options:
	-cfg   <file>   config file,detail.cfg ,must be given;
	-od     <dir>    output dir,must be given;
    -in     <dir>      input dir ,Analysis,must be given;
	-type  <str>    data type pe or se(default pe);
	-a     <number> input number,the separation of AT or GC more than the number(default 2)
	-b     <number> input number,the number of removing bases per read(default 15)
	-c     <number> the output batches (default 1)

	-h         Help
#---------------------------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}
