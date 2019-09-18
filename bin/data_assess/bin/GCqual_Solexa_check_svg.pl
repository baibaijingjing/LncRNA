#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

if (@ARGV<3) {
	print "Usage: perl $0 qufile gcfile outdir\n";
	exit;
}
###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
############################################

my @units;
my $i = 0;
my $j = 0;
my %info;
my $tmp;
my $qu = $ARGV[0];
my $gc = $ARGV[1];
my $out= $ARGV[2];
mkdir $out unless (-d "$out");
$out =~s/\/$//;
print "$qu\n$gc\n";


my $Space=10;
my $frame_width=400;
my $frame_height=300;
my $X_len=$frame_width-$Space;
my $Y_len=$frame_height-$Space;
my $left_side=20;
my $right_side=20;
my $up_side=20;
my $down_side=20;

my $TitleSize=16;
my $YtextSize=14;
my $YnumSize=12;
my $YspaceLen=3;
my $XtextSize=14;
my $XnumSize=12;
my $XspaceLen=3;

my $IntextSize=10;

my $height_paper=$frame_height+$up_side+$TitleSize+$XtextSize+$XnumSize+$down_side+2*$XspaceLen;
my $width_paper=$frame_width+$left_side+$YtextSize+$YnumSize+$right_side+$YspaceLen;

#### base quality file:
open IN,$qu or die $!;
my %quality;
my @title=split /\s+/,<IN>;
my $cycle_num=0;
while (<IN>) {
	chomp;
	next if (/^$/);
	$cycle_num++;
	my @units=split /\s+/,$_;
	for (my $i=1;$i<@units ;$i++) {
		$quality{$cycle_num}{$title[$i]}=$units[$i];
	}
}
close IN;

my @tmp=keys %quality;
my $len=600/@tmp;

#### svg output:
my $graph_qu = basename($qu).".svg";
print "$out/$graph_qu\n";
open OUT,">$out/$graph_qu" or die $!;
print OUT &svg_paper($width_paper,$height_paper);
print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width/2,$up_side+$TitleSize/2,$TitleSize,"black","Quality Distribution");  ## Title
print OUT &svg_mid_txt($left_side+$YtextSize/2-$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height/2,$YtextSize,"black","Quality Score",3);  ## Y text
print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width/2,$up_side+$TitleSize+$XspaceLen+$frame_height+$XspaceLen+$XnumSize+$XtextSize+$XspaceLen,$XtextSize,"black","Cycle Number");   ## X text

## 矩形
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen,"black",0.5);   ## up line
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## left line
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$Y_len/8-1.5,"black",0.5);   ## left line (up)
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$Y_len/8+1.5,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## left line (down)

print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## right line
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,$left_side+$YtextSize+$YnumSize+$XspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## down line

## Y轴刻度
my $Ystep=10;
my $Y_base=$Y_len/$title[-1];
my $YstepNum=$title[-1]/$Ystep;
for (my $i=0;$i<=$YstepNum ;$i++) {
	print OUT &svg_mid_txt($left_side+$YtextSize,$up_side+$TitleSize+$XspaceLen+$YnumSize*0.4+($Y_len-$i*$Ystep*$Y_base),$YnumSize,"black",$i*$Ystep);
	print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+($Y_len-$i*$Ystep*$Y_base),$left_side+$YtextSize+$YnumSize+$YspaceLen+$YspaceLen,$up_side+$TitleSize+$XspaceLen+($Y_len-$i*$Ystep*$Y_base),"black",0.5);
}

## Y轴不定坐标线
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen-3,$up_side+$TitleSize+$XspaceLen+$Y_len/8-3,$left_side+$YtextSize+$YnumSize+$YspaceLen+3,$up_side+$TitleSize+$XspaceLen+$Y_len/8,"black",0.5);
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen-3,$up_side+$TitleSize+$XspaceLen+$Y_len/8,$left_side+$YtextSize+$YnumSize+$YspaceLen+3,$up_side+$TitleSize+$XspaceLen+$Y_len/8+3,"black",0.5);

## X轴刻度
my $Xstep=50;
my $XstepNum=$cycle_num/$Xstep;
my $X_base=$X_len/$cycle_num;
for (my $i=0;$i<=$XstepNum ;$i++) {
	print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Space+$Y_len+$XnumSize*4/3,$XnumSize,"black",$i*$Xstep);
	print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Y_len+$Space,$left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Y_len+$Space-$XspaceLen,"black",0.5);
}

## reads中间虚线
print OUT &svg_dashed($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$X_len/2,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$X_len/2,$up_side+$TitleSize+$XspaceLen+$frame_height,"gray","1",0.5);

## 碱基质量分布
foreach my $cycle (sort {$a<=>$b} keys %quality) {
	foreach my $qual (sort {$a<=>$b} keys %{$quality{$cycle}}) {
		my $opacity=sqrt($quality{$cycle}{$qual}/50);
#		print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+($cycle-1)*$X_base,$up_side+$TitleSize+$XspaceLen+($Y_len-$qual*$Y_base),$X_base-1,$Y_base,"blue",$opacity);
		print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+($cycle-1)*$X_base,$up_side+$TitleSize+$XspaceLen+($Y_len-$qual*$Y_base),$X_base*0.9,$Y_base,"blue",$opacity);
	}
}

print OUT &svg_end();
close OUT;




##########################################################
#Drawing The Reads GC distribution
################Define main variables here:
my @names;
%info = ();
my %color = (
	"A(%)"=>"green",
	"T(%)"=>"red",
	"C(%)"=>"blue",
	"G(%)"=>"black",
	"N(%)"=>"gray"
);
################
open(IN,"$gc")||die"Can't open $gc\n";
@names = split(/\s+/,<IN>);
my $reads_len;
while (<IN>) {
	@units = split;
	$reads_len++;
	$info{$names[1]}{$units[0]} = $units[1];
	$info{$names[2]}{$units[0]} = $units[2];
	$info{$names[3]}{$units[0]} = $units[3];
	$info{$names[4]}{$units[0]} = $units[4];
	$info{$names[5]}{$units[0]} = $units[5];
}
close IN;

my $graph_gc = basename($gc).".svg";
print "$out/$graph_gc\n";

my @YMark=qw(0 10 20 30 40 50);

open(O,">$out/$graph_gc")||die"$!";
print O <<_FLAG_;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="800px" height="640px" version="1.1"	xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 640">
  <desc>Solexa Data GC Distribution</desc>
  <desc>lpan\@bmk.com</desc>
  <g font-family="ArialNarrow-Bold" font-size="20"  font-weight="bold" stroke-width="2" fill="none">
	<g id="top" transform="translate(100,0)">
		<rect width="800" height="80" stroke-fill="black"/>
		<text x="300" y="60" fill="black" font-size="28" style="text-anchor:middle">Base Distribution</text>
	</g>
	<g id="content" transform="translate(100,80)">
		<rect x="0" y="0" width="600" height="480" stroke-width="1" stroke="black"/>
			<g id="grid">
				<text x="-40" y="10" fill="black">$YMark[5]</text>
				<text x="-40" y="106" fill="black">$YMark[4]</text>
				<text x="-40" y="202" fill="black">$YMark[3]</text>
				<text x="-40" y="298" fill="black">$YMark[2]</text>
				<text x="-40" y="394" fill="black">$YMark[1]</text>
				<text x="-40" y="490" fill="black">$YMark[0]</text>

				<line x1="0" y1="96" x2="10" y2="96" stroke="black" stroke-width="1"/>
				<line x1="0" y1="192" x2="10" y2="192" stroke="black" stroke-width="1"/>
				<line x1="0" y1="288" x2="10" y2="288" stroke="black" stroke-width="1"/>
				<line x1="0" y1="384" x2="10" y2="384" stroke="black" stroke-width="1"/>

			</g>
			<g id="right" transform="translate(480,0)">
				<rect x="20" y="10" width="90" height="170" stroke="black" stroke-width="1"/>
				<g transform="translate(30,45)" font-family="ArialNarrow-Bold" font-weight="bold" font-size="18" fill="none">
					<text x="0" y="0" fill="green">A</text>
					<line x1="25" y1="-15" x2="60" y2="-15" stroke="green" stroke-width="1" />
					<text x="0" y="30" fill="red">T</text>
					<line x1="25" y1="15" x2="60" y2="15" stroke="red" stroke-width="1"/>
					<text x="0" y="60" fill="blue">C</text>
					<line x1="25" y1="45" x2="60" y2="45" stroke="blue" stroke-width="1"/>
					<text x="0" y="90" fill="black">G</text>
					<line x1="25" y1="75" x2="60" y2="75" stroke="black" stroke-width="1"/>
					<text x="0" y="120" fill="gray">N</text>
					<line x1="25" y1="105" x2="60" y2="105" stroke="gray" stroke-width="1"/>
				</g>
			</g>
			<g id="story">
_FLAG_

my $dis;
if ($reads_len%6==0) {
	$dis=$reads_len/6;
}
elsif ($reads_len%7==0) {
	$dis=$reads_len/7;
}
elsif ($reads_len%5==0) {
	$dis=$reads_len/5;
}

if ($reads_len%7!=0 && $reads_len%5!=0 && $reads_len%6!=0) {
	if ($reads_len>=180) {
		$dis=50;
	}
	elsif ($reads_len>=150 && $reads_len<180) {
		$dis=30;
	}
	elsif($reads_len>90 && $reads_len<150){
		$dis=20;
	}
	elsif($reads_len<=90){
		$dis=15;
	}
}
my $step=int(@tmp/$dis);
my @XMark=qw(0 1 2 3 4 5 6 7 8 9 10);
for ($i=0; $i<=$step; $i++) {
	$XMark[$i] = $i*$dis;
}
for (my $n=0;$n<=$step;$n++) {
	my $xtxt=-15+$n*600*$dis/@tmp;
	print O "<text x=\"$xtxt\" y=\"508\" fill=\"black\">$XMark[$n]</text>\n";
	my $xline;
	if ($n!=$step) {
		$xline=$n*600*$dis/@tmp+600*$dis/@tmp;
		print O "<line x1=\"$xline\" y1=\"470\" x2=\"$xline\" y2=\"480\" stroke=\"black\" stroke-width=\"0.5\"/>\n";
	}
}

foreach my $key1 (sort keys %info) {
	print O "<polyline points=\"";
	foreach my $key2 (sort {$a<=>$b} keys %{$info{$key1}}) {
		my $var1 = $key2 * $len;
		my $var2 = (50 - $info{$key1}{$key2}) * 9.6;
		print O "$var1,$var2\t";
	}
	print O "\" stroke-width=\"0.5\" stroke=\"$color{$key1}\"\/>\n";
}

print O <<_FLAG_;
		</g>
	</g>
	<g id="bottom" transform="translate(100,560)">
		<text x="300" y="60" fill="black" font-size="24" style="text-anchor:middle">Cycle Number</text>
	</g>
	<g id="left" transform="translate(0,0)">
		<text x="-120" y="0" style="text-anchor:middle" fill="black" transform="rotate(-90 120,80) " font-size="24">Percentage</text>
	</g>
  </g>
</svg>
_FLAG_
close O;

#`perl /share/nas1/chenx/bin/tools/Svg/svg2xxx_release/svg2xxx ./ ./PNG`;

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $program_name Time :[$Time_End]\n\n";
##################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub svg_paper (){#&svg_paper(width,height,[color])
	#my $svg_drawer = getlogin()."@".(`hostname`);
	my $svg_drawer = "chenx"."@"."biomarker\.com\.cn";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub max 
{
	my ($x1,$x2)=@_;
	my $max;
	if ($x1 > $x2) {
		$max=$x1;
	}
	else {
		$max=$x2;
	}
	return $max;
}

sub min 
{
	my ($x1,$x2)=@_;
	my $min;
	if ($x1 < $x2) {
		$min=$x1;
	}
	else {
		$min=$x2;
	}
	return $min;
}

sub color_gradient  #datanow,data1,data2,color1,color2
{
	my @svg_x=@_;
	my $out_color;
	if ($svg_x[0] >=$svg_x[2]) {
		$out_color=$svg_x[4];
	}
	elsif ($svg_x[0] <=$svg_x[1]) {
		$out_color=$svg_x[3];
	}
	else {
		my $tmp_red1=&hex2ten(substr($svg_x[3],1,2));
		my $tmp_gre1=&hex2ten(substr($svg_x[3],3,2));
		my $tmp_blu1=&hex2ten(substr($svg_x[3],5,2));
		my $tmp_red2=&hex2ten(substr($svg_x[4],1,2));
		my $tmp_gre2=&hex2ten(substr($svg_x[4],3,2));
		my $tmp_blu2=&hex2ten(substr($svg_x[4],5,2));
		my $new_red=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_red2-$tmp_red1)+$tmp_red1);
		my $new_gre=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_gre2-$tmp_gre1)+$tmp_gre1);
		my $new_blu=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_blu2-$tmp_blu1)+$tmp_blu1);
		$new_red=&ten2hex($new_red);$new_red="0$new_red" if(length($new_red)==1);
		$new_gre=&ten2hex($new_gre);$new_gre="0$new_gre" if(length($new_gre)==1);
		$new_blu=&ten2hex($new_blu);$new_blu="0$new_blu" if(length($new_blu)==1);
		$out_color="#$new_red$new_gre$new_blu";
	}
	return $out_color;
}

sub ten2hex  #ten
{
	my $tmp_ten=$_[0];   
	my $hex_value=uc(sprintf("%lx",$tmp_ten));
	return $hex_value;
}

sub hex2ten  #hex
{
	my $tmp_hex=$_[0];
	my $ten_value=0;
	my $tmp_i=0;
	my $tmp_j=0;
	my @tmp_x=split(//,$tmp_hex);
	my %hash_hex=(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,'A',10,'B',11,'C',12,'D',13,'E',14,'F',15);
	for ($tmp_i=@tmp_x-1;$tmp_i>=0 ;$tmp_i--) {
		$ten_value+=$hash_hex{$tmp_x[$tmp_i]}*(16**$tmp_j);
		$tmp_j++;
	}
	return $ten_value;
}

sub svg_polygon  #colorfill,colorstroke,coloropacity,point1,point2,...
{
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_circle  #&svg_circle(x,y,r,color,[info])
{
	my @svg_x=@_;
	my $line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" />\n";
	if (defined $svg_x[4]) {
		$line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" onclick=\"alert('$svg_x[4]')\" onmousemove=\"window.status='$svg_x[4]'\" />\n";
	}
	return $line;
}

sub svg_txt  #&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_mid_txt #&svg_mid_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" text-anchor=\"middle\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_dashed  #&svg_line(x1,y1,x2,y2,color,"10 5",[width])
{
	my @svg_x=@_;
	my $line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4];stroke-width:$svg_x[6]\"/>\n";
	}
	return $line;
}
sub svg_line  #&svg_line(x1,y1,x2,y2,color,[width])
{
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_rect  #&svg_rest(x,y,width,height,color,[opacity])
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
	return $line;
}

sub svg_end  #end
{
	return "</svg>\n";
}
