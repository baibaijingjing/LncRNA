use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my $version="v2.0.3";
my %config=%{readconf("$Bin/../../../../Config/lncRNA_pip.cfg")};
my $ver="1.0.0";
my ($gtf,$gff,$genome,$od,$win,$cut);
my ($units,$show_label,$sense_color,$linc_color,$anti_sense_color,$intronic_color);
GetOptions(
	"genome:s"=>\$genome,
	"od:s"=>\$od,
    "gtf:s"=>\$gtf,
    "gff:s"=>\$gff,
    "cut:s"=>\$cut,
	"win:s"=>\$win,
    "units:s"=>\$units,
    "show_lable:s"=>\$show_label,
    "sense_color:s"=>\$sense_color,
    "linc_color:s"=>\$linc_color,
    "anti_sense_color:s"=>\$anti_sense_color,
    "intronic_color:s"=>\$intronic_color,
	"h|?"=>\&help,
)or &help;
&help unless ($genome and $gtf and $gff and $od);
#open(OUT,">/run/media/bmk/A834496E34494114/perl/test/circos/cir1.txt") or die $!;
$win =$win || 1000000;
$cut = $cut || 10000000;
$units =$units || 1000000;
$show_label = $show_label || "yes";
$sense_color = $sense_color || "vdgreen";
$linc_color =$linc_color || "vdred";
$anti_sense_color =$anti_sense_color || "vdgrey";
$intronic_color =$intronic_color || "vdblue";
$genome=&ABSOLUTE_DIR($genome);
$gtf=&ABSOLUTE_DIR($gtf);
$gff=&ABSOLUTE_DIR($gff);
&MKDIR($od);
#&MKDIR("$od/circos");
$od=&ABSOLUTE_DIR($od);
my %type = (
    "j" => "sense_lncRNA",
    "e" => "sense_lncRNA",
    "o" => "sense_lncRNA",
    "u" => "lincRNA",
    "i" => "intronic_RNA",
    "x" => "anti_sense_RNA",
    
    
);
my %Fasta;
open(IN,$genome) or die $!;
open(OUT,">$od/genome_chrom.txt") or die $!;
if ($genome =~/.txt$/) {
    while (<IN>) {
        chomp;
        my ($chr,$length)=split /\s+/,$_,2;
        if ($length > $cut ) {
            print OUT "$chr\t$length\n";
            $Fasta{$chr}=$length;
        }
    }
}
if ($genome =~/.fa$/) {
    $/='>';
    <IN>;
    while (<IN>) {
        chomp;
        my ($head,$seq)=split /\n/,$_,2;
        my $chr=(split /\s+/,$head)[0];
        $seq=~s/\n//g;
        my @len=split //,$seq;
        my $length=@len;
        if ($length > $cut ) {
            print OUT "$chr\t$length\n";
            $Fasta{$chr}=$length;
        }
    }
    $/="\n";
}



close IN;
close OUT;
my $chr=&format_chr("$od/genome_chrom.txt");

open(IN,$gff) or die $!;
my %locu;
while (<IN>) {
    chomp;
    my @tmp=split /\t/,$_;
    next if ($tmp[2]!~/mRNA/);
    #print "$tmp[8]\n";
    $tmp[8]=~/^ID=([^\=]+);Parent/;
    my $start=&min($tmp[3],$tmp[4]);
    my $end=&max($tmp[3],$tmp[4]);
    my $id=$1;
    #print "$id\n";
    $locu{$id}{star}=$start;
    $locu{$id}{end}=$end;
    $locu{$id}{chr}=$tmp[0];
    #print OUT "$tmp[0]\t$start\t$end\t$";
}
open(IN1,$gtf) or die $!;
my %lnc;
my %sta;#tong ji hash
my %total;
while (<IN1>) {
    chomp;
    my @tem=split /\t+/,$_;
    $tem[8]=~/gene_id\s\"([^\"]+)\";\stranscript_id\s\"([^\"]+)\";\sexon_number\s\"(\d+)\";/;
    my $id=$2;
    #print "$2\n";
    $tem[8]=~/class_code\s\"([^\"+])\";/;
    my $code=$1;
    $lnc{$type{$code}}{$id}=$locu{$id};
    my $st=$locu{$id}{star};
    my $ed=$locu{$id}{end};
    $sta{$type{$code}}{$tem[0]}{$st}=$ed;
    $total{$tem[0]}{$st}=$ed;
    
}
#print Dumper(\%sta);
#print Dumper(\%lnc);
my @ke=sort keys %lnc;
#print Dumper(\@ke);
foreach my $ty (@ke){
    open(OUT,">$od/$ty"."_locu.txt") or die $!;
    open(OUT1,">$od/$ty"."_density.txt") or die $!;
    foreach my $l (sort keys %{$lnc{$ty}}){
        print OUT "$l\t$lnc{$ty}{$l}{chr}\t$lnc{$ty}{$l}{star}\t$lnc{$ty}{$l}{end}\n";
    }
    foreach my $id (sort keys %{$sta{$ty}}){
        #print "$id\n";
        my %hash1=%{$sta{$ty}};
        #print Dumper(\%hash1);
        if (defined $hash1{$id}) {
                my $num=0;
                my $i=$win-1;
                my @shuzu=sort {$a<=>$b} (keys %{$hash1{$id}});
                my $h=0;
                my $len=$Fasta{$id};
                #print "$len";
                while($i-$win+1<$len) {
                    foreach my $key1 ($h ... $#shuzu){
                        if ($shuzu[$key1] < $i) {
                            $num++;
                        }
                        else
                        {
                            print OUT1 "$id\t",$i-$win+1,"\t","$i\t$num\n";
                            $h=$key1;
                            $num=0;
                            last;
                        }
                    }
                print OUT1 "$id\t",$i-$win+1,"\t","$i\t$num\n" if($i>$len);
                $i=$i+$win;
                }
                $h=0;
                
            }
        }
    
    close OUT;
}
################################
open(OUT2,">$od/lncRNA_denstith.txt") or die $!;
foreach my $chro (sort keys %total){
    if (defined $Fasta{$chro}) {
                my $num=0;
                my $i=$win-1;
                my @shuzu=sort {$a<=>$b} (keys %{$total{$chro}});
                my $h=0;
                my $len=$Fasta{$chro};
                #print "$len";
                while($i-$win+1<$len) {
                    foreach my $key1 ($h ... $#shuzu){
                        if ($shuzu[$key1] < $i) {
                            $num++;
                        }
                        else
                        {
                            print OUT2 "$chro\t",$i-$win+1,"\t","$i\t$num\n";
                            $h=$key1;
                            $num=0;
                            last;
                        }
                    }
                print OUT2 "$chro\t",$i-$win+1,"\t","$i\t$num\n" if($i>$len);
                $i=$i+$win;
                }
                $h=0;
            }
}
close OUT2;
#################################################################
open(CONF,">$od/circos.conf") or die $!;

print CONF <<END;
chromosomes_units=$units
<ideogram>
    fill=yes
    label_font=default
    label_parallel=yes
    label_radius=dims(ideogram,radius)-60p
    label_size=20
    radius=0.40r
    show_label=$show_label
    <spacing>
        default=0.005r
    </spacing>
    stroke_color=dgrey
    stroke_thickness=2p
    thickness=40p
</ideogram>
karyotype=$od/chr.info
<plots>
    <axes axis>
        color=grey
        spacing=0.2r
        thickness=1
    </axes>
    <backgrounds background>
        color=vvlgrey
    </backgrounds>
    <plot>
        color=$anti_sense_color
        extend_bin=no
        file=$od/anti_sense_RNA_density.txt
        max=1122
        min=0
        r0=0.4r
        r1=0.47r+20p
        stroke_type=bin
        thickness=2
        type=histogram
    </plot>
    <plot>
        color=$linc_color
        extend_bin=no
        file=$od/lincRNA_density.txt
        max=1520
        min=0
        r0=0.66r
        r1=0.73r+20p
        stroke_type=bin
        thickness=2
        type=histogram
    </plot>
    <plot>
        color=$sense_color
        extend_bin=no
        file=$od/sense_lncRNA_density.txt
        max=1122
        min=0
        r0=0.79r
        r1=0.86r+20p
        stroke_type=bin
        thickness=2
        type=histogram
    </plot>
    <plot>
        color=$intronic_color
        extend_bin=no
        file=$od/intronic_RNA_density.txt
        max=804
        min=0
        r0=0.53r
        r1=0.6r+20p
        stroke_type=bin
        thickness=2
        type=histogram
    </plot>
    type=histogram
</plots>
show_tick_labels=yes
show_ticks=yes
spacing=10u
<ticks>
    color=black
    format=%d
    multiplier=1e-5
    radius=1r
    thickness=2p
    <tick>
        size=10p
        spacing=10u
    </tick>
    <tick>
        color=black
        format=%d
        label_offset=10p
        label_size=20p
        show_label=yes
        size=15p
        spacing=50u
        thickness=4p
    </tick>
</ticks>
<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
#<<include etc/colors_fonts_patterns.conf>>
#<<include colors.ucsc.conf>>
#<<include colors.hsv.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<image>
<<include etc/image.conf>>
</image>
<<include etc/housekeeping.conf>>
max_links* =70000
max_points_per_track* =70000
END
`$config{circos} -conf $od/circos.conf -outputdir $od -outputfile circos`;


sub format_chr{
        my$file=shift @_;
        my$i=0;
        my$col_len;
        open INC,"$file"  or die "$!";
        open OCHR,">$od/chr.info" or die "$!";
        while (my $line = <INC> ){
                chomp $line;
                my @tmp=split(/\s+/,$line);
                $col_len=@tmp if $.==1;
                last if ($col_len==7);
                next if($line=~/^#/ or $line=~/^$/ or  $line=~/^\s*$/);
                $i++;
                $i=1 if ($i>24);
                print OCHR "chr\t-\t$tmp[0]\t$tmp[0]\t0\t$tmp[1]\tchr$i\n";

        }
        close INC;
        close OCHR;
        if($col_len==7){
                return "$file";
                `rm $od/chr.info`;
        }else{
                return "$od/chr.info";
        }

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


sub max{#&max(lists or arry);
	#ÇóÁÐ±íÖÐµÄ×î´óÖµ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#ÇóÁÐ±íÖÐµÄ×îÐ¡Öµ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

sub help
{
	print <<"	Usage End.";
	Description:sefla
		Function : use draw the circos picture;
            Contact : niulg\@biomarker.com.cn
		Version  : $ver
		Usage    :
		-genome
		   genome fa or Ref_Genome/genome_size.txt ;must be given;
		-od
		   out put dir ;must be given;
        -gff
            lncRNA gff file;must be given;
        -gtf
            lncRNA grf file;must be given;
        -win
            stastic unit of chrome ;default 1000000;choise 
        -cut
            filter the chr length;default 100000000  ;choise,
        ----------------------------------------------------------------------------------
        for cicos png chosie
        -units      chr lable
        -show_lable     show chr nummber :default yes
        -sense_color        default vdgreen
        -linc_color         default vdred
        -anti_sense_color       default vdgrey
        -intronic_color                 default vdblue
        
	----------------------------------------------------------------------------------
		-h
		    Help document
        for example:
        perl $Script -genome genome.fa -gff lncRNA.gff -gtf lncRNA.gtf -od ./circos -win 1000000 -cut 100000000
	Usage End.

	exit;
}
