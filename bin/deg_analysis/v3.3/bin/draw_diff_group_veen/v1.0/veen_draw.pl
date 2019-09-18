#!/usr/bin/perl -w
#
# Copyright (c) BMK 2013
# Writer:         He hua <heh@biomarker.com.cn>
# Program Date:   2013-6-4.
# Modifier:       He hua <heh@biomarker.com.cn>
# Last Modified:  2013-6-4.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"i=s","f=s","font1=s","font=s","t=s","sample=s","o=s","h");
if ((!defined($opts{i}) && !defined $opts{f})||!defined($opts{o})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my $font1=$opts{font1}||36;
my $font=$opts{font}||24;
$font1.="px";$font.="px";
my $infile;
if (defined $opts{f}) {
	open IN,$opts{f}||die;
	my $head=<IN>;chomp$head;
	my @sample=split/\s+/,$head;
	if (@sample>5) {
		print "\tThe Sample more than 5,this process only for 2 to 5.\n";
		exit;
	}
	my %info;
	while (<IN>) {
		chomp;
		next if(/^$/);
		my @info=split/\s+/,$_;
		for (my $i=0;$i<@info ;$i++) {
			$info{$info[$i]}{$sample[$i]}++;
		}
	}
	close IN;
	my %com;
    my %element;

	foreach my $key (sort keys %info) {
		my $comm=join ",",sort keys %{$info{$key}};
		$com{$comm}++;
        $element{$comm}{$key} = 1;
	}
	
	$infile=$opts{o};
	$infile=~s/svg$/stat/;
	open OUT,">$infile"||die;
	foreach my $key (sort keys %com) {
        my $elements = join ";",(sort keys %{$element{$key}});
        print $elements."\n";
		print OUT "$key\t$com{$key}\t$elements\n";
	}
	close OUT;
}else{
	$infile=$opts{i};
}

my @ABCDE=("A","B","C","D","E");
my $count=0;my (%tmp,%sample,%sample1);
my %value=("A"=>0,"B"=>0,"C"=>0,"D"=>0,"E"=>0,"AB"=>0,"AC"=>0,"AD"=>0,"AE"=>0,"BC"=>0,"BD"=>0,"BE"=>0,"CD"=>0,"CE"=>0,"DE"=>0,"ABC"=>0,"ABD"=>0,"ABE"=>0,"ACD"=>0,"ACE"=>0,"ADE"=>0,"BCD"=>0,"BCE"=>0,"BDE"=>0,"CDE"=>0,"ABCD"=>0,"ABCE"=>0,"ABDE"=>0,"ACDE"=>0,"BCDE"=>0,"ABCDE"=>0);
open IN,$infile||die;
while (<IN>) {
	chomp;
	next if(/^$/||/^\#/);
	my ($sample,$value)=split/:|\s+/,$_;
	my @sample=split/\,/,$sample;
	for (my $i=0;$i<@sample ;$i++) {
		next if(exists $sample1{$sample[$i]});
		$count++;
		$sample{$ABCDE[$count-1]}=$sample[$i];
		$sample1{$sample[$i]}=$ABCDE[$count-1];
	}
	$tmp{$sample}=$value;
}
close IN;

my %total_value;
foreach my $sample (sort keys %tmp) {
	my @sample=split/\,/,$sample;
	my @NABCDE;
	for (my $i=0;$i<@sample;$i++) {
		push @NABCDE,$sample1{$sample[$i]};
		$total_value{$sample1{$sample[$i]}}+=$tmp{$sample} if(exists $sample1{$sample[$i]});
	}
#	print "$sample\n";
	my $NABCDE=join "",sort @NABCDE;
	$value{$NABCDE}=$tmp{$sample};
}

my $svgout;my $svgout1;
if ($count==2) {
	$font1=14;$font=9;
	$svgout=&veen2($font1,$font,$sample{"A"},$sample{"B"},$total_value{"A"},$total_value{"B"},$value{"A"},$value{"B"},$value{"AB"});
}elsif($count==3){
	$svgout=&veen3($font1,$font,$sample{"A"},$sample{"B"},$sample{"C"},$total_value{"A"},$total_value{"B"},$total_value{"C"},$value{"A"},$value{"B"},$value{"C"},$value{"AB"},$value{"AC"},$value{"BC"},$value{"ABC"});
}elsif($count==4){
	$svgout1=&veen4($font1,$font,$sample{"A"},$sample{"B"},$sample{"C"},$sample{"D"},$total_value{"A"},$total_value{"B"},$total_value{"C"},$total_value{"D"},$value{"A"},$value{"B"},$value{"C"},$value{"D"},$value{"AB"},$value{"AC"},$value{"AD"},$value{"BC"},$value{"BD"},$value{"CD"},$value{"ABC"},$value{"ABD"},$value{"ACD"},$value{"BCD"},$value{"ABCD"});
	$svgout=&veen4_1($font1,$font,$sample{"A"},$sample{"B"},$sample{"C"},$sample{"D"},$total_value{"A"},$total_value{"B"},$total_value{"C"},$total_value{"D"},$value{"A"},$value{"B"},$value{"C"},$value{"D"},$value{"AB"},$value{"AC"},$value{"AD"},$value{"BC"},$value{"BD"},$value{"CD"},$value{"ABC"},$value{"ABD"},$value{"ACD"},$value{"BCD"},$value{"ABCD"});
#                                                                                            $A,         $B,         $C,          $D,         $AB,          $AC,        $AD,         $BC,        $BD,         $CD,         $ABC,         $ABD,         $ACD,         $BCD,          $ABCD
}else{
	$svgout=&veen5($font1,$font,$sample{"A"},$sample{"B"},$sample{"C"},$sample{"D"},$sample{"E"},$total_value{"A"},$total_value{"B"},$total_value{"C"},$total_value{"D"},$total_value{"E"},$value{"A"},$value{"B"},$value{"C"},$value{"D"},$value{"E"},$value{"AB"},$value{"AC"},$value{"AD"},$value{"AE"},$value{"BC"},$value{"BD"},$value{"BE"},$value{"CD"},$value{"CE"},$value{"DE"},$value{"ABC"},$value{"ABD"},$value{"ABE"},$value{"ACD"},$value{"ACE"},$value{"ADE"},$value{"BCD"},$value{"BCE"},$value{"BDE"},$value{"CDE"},$value{"ABCD"},$value{"ABCE"},$value{"ABDE"},$value{"ACDE"},$value{"BCDE"},$value{"ABCDE"});
#																										$A,			$B,			$C,			$D,			$E,			$AB,			$AC,		$AD,		$AE,		$BC,		$BD,			$BE,		$CD,		$CE,		$DE,		$ABC,			$ABD,			$ABE,		$ACD,			$ACE,		$ADE,			$BCD,		$BCE,			$BDE,		$CDE,			$ABCD,			$ABCE,		$ABDE,			$ACDE,			$BCDE,		$ABCDE
}

open OUT,">$opts{o}"||die;
print OUT $svgout;
close OUT;
$opts{o}=&path($opts{o});
my $dir=dirname$opts{o};
my $outfile=basename($opts{o});
chdir $dir;

my $type=$opts{t}|| "png";
`perl /share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx $outfile -t $type `;
if ($count==4) {
	$outfile=~s/.svg$/_1.svg/;
	open OUT,">$outfile"||die;
	print OUT $svgout1;
	close OUT;
	`/share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx $outfile -t $type `;
}
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub path{
	my ($in)=@_;
	my $return;
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
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
sub help{
	print << "	Usage End.";
	Description:Draw Venn\'s construction for 2 to 5 sets base on data.
		version:$ver
	Usage:
		-i            input file            
		-f            input list file       
		-font1        Sample text font size(36)
		-font         digit font size(24)
		-t            type [jpg|jpeg|jpe|png|tiff|pdf|ps](png)
		-o            out file
	eg: 
	   perl $program_name -i input.file -o out.svg
	   perl $program_name -f input.list.file -o out.svg

	input.file format:               input.list.file format:
	Pe:123                             Sample1    Sample2
	Ab,Pe:643                           gene1     gene2
	Ab:592                              gene2     gene3
	                                    gene3

	Usage End.
		exit;
}
##################sub svg###############
sub veen2{
my ($fontpx1,$fontpx,$TA,$TB,$TAn,$TBn,$A,$B,$AB)=@_;
my $SVG=<<"SVGList";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->
<svg xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.0" width="250.21001" height="197.71001" id="svg2">
  <defs
	 id="defs4">
	<linearGradient
	   id="linearGradient4553">
	  <stop
		 style="stop-color:#00ffdc;stop-opacity:0.79215693"
		 offset="0"
		 id="stop4555" />
	  <stop
		 style="stop-color:#13c1ec;stop-opacity:0.79215693"
		 offset="1"
		 id="stop4557" />
	</linearGradient>
	<linearGradient
	   id="linearGradient3661">
	  <stop
		 style="stop-color:#ff8100;stop-opacity:0.7904762"
		 offset="0"
		 id="stop3663" />
	  <stop
		 style="stop-color:#ecd913;stop-opacity:0.7904762"
		 offset="1"
		 id="stop3665" />
	</linearGradient>
	<linearGradient
	   id="linearGradient2762">
	  <stop
		 style="stop-color:#ff8100;stop-opacity:0.7904762"
		 offset="0"
		 id="stop2764" />
	  <stop
		 style="stop-color:#ecd913;stop-opacity:0.7904762"
		 offset="1"
		 id="stop2766" />
	</linearGradient>
	<linearGradient
	   x1="84.285713"
	   y1="108.07647"
	   x2="51.428574"
	   y2="55.933605"
	   id="linearGradient2768"
	   xlink:href="#linearGradient2762"
	   gradientUnits="userSpaceOnUse" />
	<linearGradient
	   x1="84.285713"
	   y1="108.07647"
	   x2="51.428574"
	   y2="55.933605"
	   id="linearGradient3659"
	   xlink:href="#linearGradient4553"
	   gradientUnits="userSpaceOnUse" />
	<linearGradient
	   x1="84.285713"
	   y1="108.07647"
	   x2="51.428574"
	   y2="55.933605"
	   id="linearGradient4579"
	   xlink:href="#linearGradient2762"
	   gradientUnits="userSpaceOnUse" />
	<linearGradient
	   x1="84.285713"
	   y1="108.07647"
	   x2="51.428574"
	   y2="55.933605"
	   id="linearGradient4581"
	   xlink:href="#linearGradient4553"
	   gradientUnits="userSpaceOnUse" />
  </defs>
  <g
	 transform="translate(-19.28572,-20.21933)"
	 id="layer1">
	<g
	   transform="translate(3.675958,2.042199)"
	   id="g4567">
	  <path
		 d="M 170 95.576469 A 75.35714 75.35714 0 1 1  19.285721,95.576469 A 75.35714 75.35714 0 1 1  170 95.576469 z"
		 style="opacity:1;fill:url(#linearGradient4579);fill-opacity:1;stroke:#343434;stroke-width:0;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1"
		 id="path2760" />
	  <path
		 d="M 170 95.576469 A 75.35714 75.35714 0 1 1  19.285721,95.576469 A 75.35714 75.35714 0 1 1  170 95.576469 z"
		 transform="translate(92.49999,0)"
		 style="opacity:1;fill:url(#linearGradient4581);fill-opacity:1;stroke:#343434;stroke-width:0;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1"
		 id="path3657" />
	  <text
		 x="60.714287"
		 y="103.12558"
		 style="font-size:$fontpx;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4559"
		 xml:space="preserve"><tspan
		   x="60.714287"
		   y="103.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4561">$A</tspan></text>
	  <text
		 x="213.65625"
		 y="103.12558"
		 style="font-size:$fontpx;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="213.65625"
		   y="103.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$B</tspan></text>
	  <text
		 x="140.65625"
		 y="103.12558"
		 style="font-size:$fontpx;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="140.65625"
		   y="103.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$AB</tspan></text>
	  <text
		 x="100.65625"
		 y="190.12558"
		 style="font-size:$fontpx1;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="100.65625"
		   y="190.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$TA</tspan></text>
	  <text
		 x="100.65625"
		 y="210.12558"
		 style="font-size:$fontpx;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="100.65625"
		   y="210.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$TAn</tspan></text>
	  <text
		 x="190.65625"
		 y="190.12558"
		 style="font-size:$fontpx1;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="190.65625"
		   y="190.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$TB</tspan></text>
	  <text
		 x="190.65625"
		 y="210.12558"
		 style="font-size:$fontpx;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
		 id="text4563"
		 xml:space="preserve"><tspan
		   x="190.65625"
		   y="210.12558"
		   style="font-weight:bold;text-align:center;text-anchor:middle"
		   id="tspan4565">$TBn</tspan></text>
	</g>
  </g>
</svg>
SVGList
	return $SVG;
}
sub veen3{
	my ($fontpx1,$fontpx,$TA,$TB,$TC,$TAn,$TBn,$TCn,$A,$B,$C,$AB,$AC,$BC,$ABC)=@_;
	my $SVG=<<"SVGList";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->
<svg xmlns:svg="http://www.w3.org/2000/svg" xmlns="http://www.w3.org/2000/svg" version="1.0" width="600" height="600" id="svg2837">
  <defs
	 id="defs2839" />
  <g
	 transform="translate(76.892021,-667.37832)"
	 id="layer1">
	<path
	   d="M 300.01761,904.03602 C 298.25491,913.0464 298.56029,922.56147 296.33758,931.59089 C 289.51657,966.24729 270.10015,997.7607 242.28249,1019.5389 C 237.18789,1023.9771 230.99795,1027.3024 225.8435,1031.6253 C 229.82357,1035.3331 235.59033,1037.3064 240.57566,1039.6788 C 259.3003,1047.2478 279.34269,1052.533 299.70511,1051.4021 C 314.09126,1051.7628 328.47817,1050.1429 342.22317,1045.8073 C 353.40077,1042.4419 364.61104,1037.9755 374.43799,1031.8777 C 374.75968,980.36255 346.05583,929.52871 301.20064,903.90209 C 300.81686,903.75164 300.36286,903.71354 300.01761,904.03602 L 300.01761,904.03602 z"
	   id="path3690"
	   style="fill:#00ffff;fill-opacity:1;fill-rule:nonzero;stroke:none" />
	<path
	   d="M 218.32118,777.09406 C 180.6431,801.77205 154.71884,843.7888 151.12957,888.89855 C 151.14617,892.93111 149.51169,897.37709 150.73189,901.24584 C 154.14527,902.28867 157.29452,898.61466 160.48636,897.78602 C 181.00431,887.92631 203.74253,882.54335 226.64701,883.48703 C 250.94378,882.91255 274.6615,889.76595 295.976,900.90353 C 297.54618,902.47489 299.7691,901.34406 299.19171,899.12531 C 299.21601,849.04114 270.67066,800.78611 228.20804,774.81184 C 224.53679,770.69805 221.74042,775.40097 218.32118,777.09406 L 218.32118,777.09406 z"
	   id="path3686"
	   style="fill:#ff00ff;fill-opacity:1;fill-rule:nonzero;stroke:none" />
	<path
	   d="M 148.9685,903.54495 C 106.88819,928.32456 78.183152,974.54978 75.772032,1023.5209 C 76.029539,1026.9345 73.662998,1032.5602 78.522093,1033.6898 C 96.097968,1043.1378 115.55096,1049.335 135.44171,1051.2235 C 153.43203,1051.6596 171.79605,1052.1073 189.14197,1046.6311 C 201.42843,1043.1567 213.60804,1038.4483 224.50421,1031.9824 C 221.81108,1028.5521 216.77757,1026.8654 213.32118,1023.9244 C 177.32478,998.60924 153.93112,956.26523 150.99629,912.29959 C 150.27983,909.51442 151.447,905.71415 149.43725,903.50031 L 148.9685,903.54491 L 148.9685,903.54495 z"
	   id="path3688"
	   style="fill:#ffff00;fill-opacity:1;fill-rule:nonzero;stroke:none" />
	<path
	   d="M 150,752.875 C 67.433431,752.875 0.5,819.80843 0.5,902.375 C 0.5,957.80941 30.682968,1006.1848 75.5,1032 C 75.635483,976.72214 105.78388,928.47579 150.5,902.71875 C 150.49974,902.60347 150.5,902.49034 150.5,902.375 C 150.5,847.14452 180.45541,798.88546 225,773 C 202.95285,760.1881 177.33609,752.875 150,752.875 z"
	   id="path3670"
	   style="fill:#ff0000;fill-opacity:1;stroke:#000000;stroke-width:1;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
	<path
	   d="M 225,773 C 269.54459,798.88546 299.5,847.14452 299.5,902.375 C 299.5,902.49034 299.50026,902.60347 299.5,902.71875 C 344.21612,928.47579 374.36452,976.72214 374.5,1032 C 419.31703,1006.1848 449.5,957.80941 449.5,902.375 C 449.5,819.80843 382.56657,752.875 300,752.875 C 272.66391,752.875 247.04715,760.1881 225,773 z"
	   id="path3668"
	   style="fill:#0000ff;fill-opacity:1;stroke:#000000;stroke-width:1;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
	<path
	   d="M 374.5,1032 C 352.56448,1044.6351 327.13216,1051.875 300,1051.875 C 272.66391,1051.875 247.04715,1044.5306 225,1031.7188 C 202.95285,1044.5306 177.33609,1051.875 150,1051.875 C 122.86784,1051.875 97.435522,1044.6351 75.5,1032 C 75.499694,1032.1247 75.5,1032.2502 75.5,1032.375 C 75.5,1114.9416 142.43343,1181.875 225,1181.875 C 307.56657,1181.875 374.5,1114.9416 374.5,1032.375 C 374.5,1032.2502 374.50031,1032.1247 374.5,1032 z"
	   id="path3666"
	   style="fill:#00ff00;fill-opacity:1;stroke:#000000;stroke-width:1;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
	<path
	   d="M 225,1031.7188 C 269.45156,1005.8874 299.37523,957.80466 299.5,902.71875 C 277.56448,890.08362 252.13216,882.875 225,882.875 C 197.86784,882.875 172.43552,890.08362 150.5,902.71875 C 150.62477,957.80466 180.54844,1005.8874 225,1031.7188 z"
	   id="path2847"
	   style="fill:#ffffff;fill-opacity:1;stroke:#000000;stroke-width:1;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
	<path
	   d="M 424.46112,127.06245 C 424.50609,477.3992 127.72168,761.42607 -238.39599,761.42607 C -604.51367,761.42607 -901.29808,477.3992 -901.25311,127.06245 C -901.29808,-223.2743 -604.51367,-507.30117 -238.396,-507.30117 C 127.72168,-507.30117 424.50609,-223.2743 424.46112,127.06245 L 424.46112,127.06245 z"
	   transform="matrix(-0.2557112,0,0,-0.2672297,89.039464,936.31704)"
	   id="path3623"
	   style="opacity:0.5;fill:none;stroke:none" />
	<path
	   d="M 424.46112,127.06245 C 424.50609,477.3992 127.72168,761.42607 -238.39599,761.42607 C -604.51367,761.42607 -901.29808,477.3992 -901.25311,127.06245 C -901.29808,-223.2743 -604.51367,-507.30117 -238.396,-507.30117 C 127.72168,-507.30117 424.50609,-223.2743 424.46112,127.06245 L 424.46112,127.06245 z"
	   transform="matrix(0,-0.2557112,0.2672297,0,266.04513,841.40164)"
	   id="path3627"
	   style="opacity:0.5;fill:none;stroke:none" />
	<path
	   d="M 424.46112,127.06245 C 424.50609,477.3992 127.72168,761.42607 -238.39599,761.42607 C -604.51367,761.42607 -901.29808,477.3992 -901.25311,127.06245 C -901.29808,-223.2743 -604.51367,-507.30117 -238.396,-507.30117 C 127.72168,-507.30117 424.50609,-223.2743 424.46112,127.06245 L 424.46112,127.06245 z"
	   transform="matrix(0.2370913,9.579111e-2,-0.100106,0.2477711,294.24133,1023.7161)"
	   id="path3631"
	   style="opacity:0.5;fill:none;stroke:none" />

	<text
	   x="-6.714287"
	   y="885.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx1;font-style:normal;font-weight:normal;text-align:center;text-anchor:end;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="-6.714287"
		 y="885.07648"
		 id="tspan4233"
		 style="font-size:$fontpx1;text-align:center;text-anchor:end;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TA</tspan></text>
	<text
	   x="455.714287"
	   y="885.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx1;font-style:normal;font-weight:normal;text-align:center;text-anchor:left;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="455.714287"
		 y="885.07648"
		 id="tspan4233"
		 style="font-size:$fontpx1;text-align:center;text-anchor:left;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TB</tspan></text>
	<text
	   x="223.714287"
	   y="1215.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx1;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="223.714287"
		 y="1215.07648"
		 id="tspan4233"
		 style="font-size:$fontpx1;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TC</tspan></text>

	<text
	   x="-6.714287"
	   y="908.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:end;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="-6.714287"
		 y="908.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:end;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TAn</tspan></text>
	<text
	   x="455.714287"
	   y="908.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:left;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="455.714287"
		 y="908.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:left;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TBn</tspan></text>
	<text
	   x="223.714287"
	   y="1238.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="223.714287"
		 y="1238.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$TCn</tspan></text>

	<text
	   x="223.714287"
	   y="845.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="223.714287"
		 y="845.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$AB</tspan></text>
	<text
	   x="223.714287"
	   y="945.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="223.714287"
		 y="945.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$ABC</tspan></text>
	<text
	   x="125.714287"
	   y="1015.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="125.714287"
		 y="1015.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$AC</tspan></text>
	<text
	   x="315.714287"
	   y="1015.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="315.714287"
		 y="1015.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#000000;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$BC</tspan></text>

	<text
	   x="79.714287"
	   y="885.07648"
	   id="text4231"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="79.714287"
		 y="885.07648"
		 id="tspan4233"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$A</tspan></text>
	<text
	   x="371.71429"
	   y="885.07648"
	   id="text4235"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="371.71429"
		 y="885.07648"
		 id="tspan4237"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$B</tspan></text>
	<text
	   x="223.71429"
	   y="1109.0764"
	   id="text4239"
	   xml:space="preserve"
	   style="font-size:$fontpx;font-style:normal;font-weight:normal;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;stroke:none;font-family:Sans"><tspan
		 x="223.71429"
		 y="1109.0764"
		 id="tspan4241"
		 style="font-size:$fontpx;text-align:center;text-anchor:middle;fill:#ffffff;fill-opacity:1;font-family:Sans;-inkscape-font-specification:Sans">$C</tspan></text>
  </g>
</svg>
SVGList
	return $SVG;
}

sub veen4{
	my ($fontpx1,$fontpx,$TA,$TB,$TC,$TD,$TAn,$TBn,$TCn,$TDn,$A,$B,$C,$D,$AB,$AC,$AD,$BC,$BD,$CD,$ABC,$ABD,$ACD,$BCD,$ABCD)=@_;
	my $SVG=<<"SVGList";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->
<svg xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.w3.org/2000/svg" height="481.65" width="518.79" version="1.1" xmlns:cc="http://creativecommons.org/ns#" xmlns:dc="http://purl.org/dc/elements/1.1/">
 <metadata>
  <rdf:RDF>
   <cc:Work rdf:about="">
	<dc:format>image/svg+xml</dc:format>
	<dc:type rdf:resource="http://purl.org/dc/dcmitype/StillImage"/>
	<dc:title/>
   </cc:Work>
  </rdf:RDF>
 </metadata>
 <g transform="translate(-116.89 -293.82)">
  <g transform="translate(2.7143 2.2857)" stroke="#000" stroke-linecap="round" stroke-width="1.7158">
   <path opacity="0.641" d="m441.43 500.93c0 82.054-66.518 148.57-148.57 148.57-82.054 0-148.57-66.518-148.57-148.57 0-82.054 66.518-148.57 148.57-148.57 82.05 0 148.57 66.52 148.57 148.57z" transform="matrix(.96163 0 0 .96163 5.619 -18)" fill="#ff5c5c"/>
   <path opacity="0.641" d="m441.43 500.93c0 82.054-66.518 148.57-148.57 148.57-82.054 0-148.57-66.518-148.57-148.57 0-82.054 66.518-148.57 148.57-148.57 82.05 0 148.57 66.52 148.57 148.57z" transform="matrix(.96163 0 0 .96163 177.52 -18)" fill="#739aff"/>
   <path opacity="0.641" d="m441.43 500.93c0 82.054-66.518 148.57-148.57 148.57-82.054 0-148.57-66.518-148.57-148.57 0-82.054 66.518-148.57 148.57-148.57 82.05 0 148.57 66.52 148.57 148.57z" transform="matrix(.96163 0 0 .96163 5.619 119)" fill="#ffff6a"/>
   <path opacity="0.641" d="m441.43 500.93c0 82.054-66.518 148.57-148.57 148.57-82.054 0-148.57-66.518-148.57-148.57 0-82.054 66.518-148.57 148.57-148.57 82.05 0 148.57 66.52 148.57 148.57z" transform="matrix(.96163 0 0 .96163 177.52 119)" fill="#73ff73"/>
  </g>
 </g>

  <text x="35" y="45" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TA</text>
  <text x="480" y="45" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TB</text>
  <text x="480" y="420" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TC</text>
  <text x="35" y="420" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TD</text>
  <text x="35" y="70" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TAn</text>
  <text x="480" y="70" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TBn</text>
  <text x="480" y="450" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TCn</text>
  <text x="35" y="450" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TDn</text>

  <text x="125" y="95" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$A</text>
  <text x="375" y="95" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$B</text>
  <text x="258" y="135" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$AB</text>
  <text x="375" y="385" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$C</text>
  <text x="125" y="385" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$D</text>
  <text x="375" y="245" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$BC</text>
  <text x="258" y="345" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$CD</text>
  <text x="125" y="245" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$AD</text>

  <text x="290" y="200" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABC</text>
  <text x="290" y="288" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$BCD</text>
  <text x="228" y="287" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ACD</text>
  <text x="228" y="200" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABD</text>
  <text x="258" y="244" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABCD</text>

</svg>
SVGList
	return $SVG;
}

sub veen4_1{
	my ($fontpx1,$fontpx,$TA,$TB,$TC,$TD,$TAn,$TBn,$TCn,$TDn,$A,$B,$C,$D,$AB,$AC,$AD,$BC,$BD,$CD,$ABC,$ABD,$ACD,$BCD,$ABCD)=@_;
	my $SVG=<<"SVGList";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="640" height="420">
<defs>
	<ellipse id="elli" rx="187.33333" ry="115.33333" />
</defs>

<g style="fill-opacity:0.2; stroke:black; stroke-width:2.4"
   transform="translate(-145,-775)">

	<g transform="translate(529 1024)">
		<use transform="rotate(-40)" fill="#6fff05" xlink:href="#elli" />
	</g>
	<g transform="translate(457 938)">
		<use transform="rotate(-40)" fill="#ff6405" xlink:href="#elli" />
	</g>
	<g transform="translate(460 938)">
		<use transform="rotate(40)" fill="#0525ff" xlink:href="#elli" />
	</g>
	<g transform="translate(388 1024)">
		<use transform="rotate(40)" fill="#1e1e1e" xlink:href="#elli" />
	</g>
  </g>
  
  <text x="35" y="160" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TA</text>
  <text x="125" y="30" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TB</text>
  <text x="505" y="30" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TC</text>
  <text x="578" y="160" dy="0.7ex" font-size="$fontpx1" text-anchor="middle">$TD</text>
  <text x="35" y="189" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TAn</text>
  <text x="125" y="59" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TBn</text>
  <text x="505" y="59" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TCn</text>
  <text x="578" y="189" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$TDn</text>

  <text x="115" y="157" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$A</text>
  <text x="204" y="60" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$B</text>
  <text x="415" y="60" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$C</text>
  <text x="494" y="157" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$D</text>

  <text x="175" y="117" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$AB</text>
  <text x="310" y="90" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$BC</text>
  <text x="453" y="117" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$CD</text>
  <text x="185" y="248" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$AC</text>
  <text x="310" y="357" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$AD</text>
  <text x="430" y="248" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$BD</text>
   <text x="230" y="175" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABC</text>
  <text x="384" y="175" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$BCD</text>
  <text x="255" y="295" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ACD</text>
  <text x="370" y="295" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABD</text>
  <text x="310" y="245" dy="0.7ex" font-size="$fontpx" text-anchor="middle">$ABCD</text>

</svg>
SVGList
	return $SVG;
}

sub veen5{
	my ($fontpx1,$fontpx,$TA,$TB,$TC,$TD,$TE,$TAn,$TBn,$TCn,$TDn,$TEn,$A,$B,$C,$D,$E,$AB,$AC,$AD,$AE,$BC,$BD,$BE,$CD,$CE,$DE,$ABC,$ABD,$ABE,$ACD,$ACE,$ADE,$BCD,$BCE,$BDE,$CDE,$ABCD,$ABCE,$ABDE,$ACDE,$BCDE,$ABCDE)=@_;
	my $SVG=<<"SVGList";
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="940" height="840" viewBox="-390 -450 800 800">
 <title>Radially-symmetrical Five-set Venn Diagram</title>
 <desc>Devised by Branko Gruenbaum and rendered by CMG Lee.</desc>
 <defs>
  <ellipse id="ellipse" cx="36" cy="-56" rx="160" ry="320" />
  <g id="ellipses">
   <use xlink:href="#ellipse" fill="#0000ff" />
   <use xlink:href="#ellipse" fill="#0099ff" transform="rotate(72)" />
   <use xlink:href="#ellipse" fill="#00cc00" transform="rotate(144)" />
   <use xlink:href="#ellipse" fill="#cc9900" transform="rotate(216)" />
   <use xlink:href="#ellipse" fill="#ff0000" transform="rotate(288)" />
  </g>
 </defs>
 <use xlink:href="#ellipses" fill-opacity="0.3" />
 <use xlink:href="#ellipses" fill-opacity="0" stroke="#000" stroke-width="2" />
 <g text-anchor="middle" font-family="sans-serif" font-size="16">

  <text x="30"   y="-425" dy="0.7ex" text-anchor="middle" font-size="$fontpx1">$TA</text>
  <text x="410"  y="-92"  dy="0.7ex" font-size="$fontpx1">$TB</text>
  <text x="260"  y="280"  dy="0.7ex" font-size="$fontpx1">$TC</text>
  <text x="-340" y="240"  dy="0.7ex" text-anchor="middle" font-size="$fontpx1">$TD</text>
  <text x="-395" y="-175" dy="0.7ex" text-anchor="middle" font-size="$fontpx1">$TE</text>
  <text x="30"   y="-393" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$TAn</text>
  <text x="410"  y="-60"  dy="0.7ex" font-size="$fontpx">$TBn</text>
  <text x="260"  y="312"  dy="0.7ex" font-size="$fontpx">$TCn</text>
  <text x="-340" y="272"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$TDn</text>
  <text x="-395" y="-140" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$TEn</text>

  <text x="30"   y="-300" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$A</text>
  <text x="300"  y="-60"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$B</text>
  <text x="160"  y="280"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$C</text>
  <text x="-220" y="220"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$D</text>
  <text x="-280" y="-130" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$E</text>
  <text x="180"  y="-130" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$AB</text>
  <text x="40"   y="230"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$AC</text>
  <text x="100"  y="-200" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$AD</text>
  <text x="-80"  y="-215" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$AE</text>
  <text x="190"  y="125"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BC</text>
  <text x="-190" y="120"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BD</text>
  <text x="230"  y="40"   dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BE</text>
  <text x="-60"  y="220"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$CD</text>
  <text x="-170" y="-150" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$CE</text>
  <text x="-222" y="0"    dy="0.7ex" text-anchor="middle" font-size="$fontpx">$DE</text>
  <text x="90"   y="150"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABC</text>
  <text x="148"  y="-153" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABD</text>
  <text x="170"  y="-20"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABE</text>
  <text x="-33"  y="208"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ACD</text>
  <text x="-93"  y="-193" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ACE</text>
  <text x="20"   y="-180" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ADE</text>
  <text x="-120" y="120"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BCD</text>
  <text x="190"  y="100"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BCE</text>
  <text x="-211" y="32"   dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BDE</text>
  <text x="-150" y="-80"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$CDE</text>
  <text x="-30"  y="160"  dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABCD</text>
  <text x="140"  y="80"   dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABCE</text>
  <text x="120"  y="-100" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABDE</text>
  <text x="-60"  y="-140" dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ACDE</text>
  <text x="-160" y="20"   dy="0.7ex" text-anchor="middle" font-size="$fontpx">$BCDE</text>
  <text x="0"    y="0"    dy="0.7ex" text-anchor="middle" font-size="$fontpx">$ABCDE</text>
 </g>
</svg>
SVGList
	return $SVG;
}
