#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Data::Dumper;
use File::Basename qw(basename dirname);
my ($mi_lnc,$lnc_m,$mi_m,$od,$in1,$in2,$mim);
GetOptions(
    "miRNA2lncRNA:s" =>\$mi_lnc,
    "lncRNA2mRNAtar:s"=>\$lnc_m,
    "miRNA2mRNA:s"=>\$mi_m,
    "od:s"=>\$od,
    "in1:s"=>\$in1,
    "in2:s"=>\$in2,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($mi_lnc and $lnc_m  and $in1 and $in2 and $od);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
$in1=&ABSOLUTE_DIR($in1);
$in2=&ABSOLUTE_DIR($in2);
$mi_lnc=&ABSOLUTE_DIR($mi_lnc);
$lnc_m=&ABSOLUTE_DIR($lnc_m);

my (%m_s,%m_l,%l_s);
&get_hash_tar($mi_lnc,\%l_s,"$od/miRNA_lncRNA.regulation_cytoscape.input.txt","$od/miRNA_lncRNA.regulation_cytoscape.Attribute.list","sl");
&get_hash_tar($lnc_m,\%m_l,"$od/lncRNA_mRNA.regulation_cytoscape.input.txt","$od/lncRNA_mRNA.regulation_cytoscape.Attribute.list","lm");
if (defined $mi_m and $mi_m ne ""){
	$mi_m=&ABSOLUTE_DIR($mi_m);
	&get_hash_tar($mi_m,\%m_s,"$od/miRNA_mRNA.regulation_cytoscape.input.txt","$od/miRNA_mRNA.regulation_cytoscape.Attribute.list","sm");
	#print Dumper(\%m_s);
	open SHARE,">$od/miRNA2lncRNA_share_target.list" or die $!;
	foreach my $gene (keys %m_s){
		if (exists $m_l{$gene}){
			my $ll=$m_l{$gene};
			my $ss=$m_s{$gene};

			print SHARE "$gene\t".join(",",@$ll).";".join(",",@$ss)."\n";
			#print SHARE"$gene\t$lcl;$scs";
			}
	}
}
my (%lncrna,%mrna, %srna);##Attribute file
open OUT2, ">$od/miRNA_lncRNA_mRNA.Attribute.list";
my @d=glob("$in1/*_vs_*");
foreach my $dir (@d){
	my $name=basename($dir);
	my $degl="$dir/$name.DEG_final.xls";
	print "$degl\n";
	my $n=$name;
	$n=~s/L/S/g;
	my $degmi="$in2/$n/$n.DEG_final.xls";
	print "$degmi\n";
	open DEG,$degl;
	my %id;
	while (<DEG>){
		chomp;
		my $i=(split /\s+/,$_)[0];
		$id{$i}=1;		
	}
	close DEG;
	open MI,$degmi;
	my %d_mi;
	while (<MI>){
		chomp;
		my $j=(split /\s+/,$_)[0];
		$d_mi{$j}=1;
	}
	close MI;
	open OUT1,">$od/$name.miRNA_lncRNA_mRNA.Interaction.list";
	#open OUT2, ">$od/$name.miRNA_lncRNA_mRNA.Attribute.list";
	foreach my $lg (keys %m_l){
		if (!exists $mrna{$lg}){
			print OUT2 "$lg\tm\n";
			$mrna{$lg}=1;
		}
		my @lrna=@{$m_l{$lg}};
		foreach my $lnc (@lrna){
			if (!exists $lncrna{$lnc}){
			print OUT2 "$lnc\tl\n";
			$lncrna{$lnc}=1;
			}
			if (exists $id{$lnc}){
				print OUT1 "$lnc\t$lg\tml\n";
			}
				
		}
	}
	
	if (defined $mi_m and $mi_m ne ""){
		foreach my $sg (keys %m_s){
			if (!exists $mrna{$sg}){

				print OUT2 "$sg\tm\n";
				$mrna{$sg}=1;
			}
			my @mirna=@{$m_s{$sg}};
			foreach my $s (@mirna){
				if (!exists $srna{$s}) {
				
					print OUT2 "$s\ts\n"; 
					$srna{$s}=1;
				}
				if (exists $d_mi{$s}){
					print OUT1 "$s\t$sg\tsm\n";
				}
			}
		}
	}	
	foreach my $key (keys %l_s){
		if (!exists $srna{$key}) {
			print OUT2 "$key\ts\n";
			$srna{$key}=1;
		}
		if (exists $id{$key}){
			my @star=@{$l_s{$key}};
			foreach my $val (@star){
				if (!exists $lncrna{$val}) {

	                                print OUT2 "$val\tl\n";
					$lncrna{$val}=1;
				}
				if (exists $d_mi{$val}){
					print OUT1 "$key\t$val\tls\n";
					}
				}
			}
	}
	close OUT1;
	
} 

close OUT2;




sub get_hash_tar {
	my ($tar,$hash,$out_inter,$out_att,$type)=@_;
	open IN,"$tar" or die $!;
	open OUT1, ">$out_inter";
	open OUT2, ">$out_att";
	my %at;#shu xing 
	while (<IN>){
		chomp;
		next if (/^#/ or /^$/);
		my ($t,$tem)=(split /\s+/,$_)[0,1];
		print OUT2 "$t\tt\n";
		if ($tem=~/;/){
			my @tar=split /;/,$tem;
			foreach my $gene (@tar){
				print OUT1 "$t\t$gene\t$type\n";
				if (!exists $at{$gene}){
					print OUT2 "$gene\tg\n";
					$at{$gene}=1;
					}
				if (!exists $hash->{$gene}){
					my @mio=();
					$mio[0]=$t;
					$hash->{$gene}=\@mio;
				
				}else{
					my $mir=$hash->{$gene};
					push(@$mir,$t);
				}
			}
	                
		}else{
			if (!exists $hash->{$tem}){
				my @mio=();
				$mio[0]=$t;
				$hash->{$tem}=\@mio;
			}else{
				my $mir=$hash->{$tem};
				push(@$mir,$t);
			}
		}
		
	}
	close IN;
	close OUT1;
	close OUT2;
}




sub MKDIR{
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}



sub ABSOLUTE_DIR{ 
        my $cur_dir=`pwd`;chomp($cur_dir);
        my ($in)=@_;
        my $return="";
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
sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: 1.0
     Usage:
            -mi_lnc      <FILE>  ./lncRNA_target2mirna.list
            -lnc_m       <FILE>   ./novel_lncRNA_target.xls
            -mi_m        <FILE>   *.mir2target.list
	          
            -od          dir
	         -in1          Lnc_Diff_Analysis dir
	       -in2	miRNA  Diff_Analysis/
	        -h 		help documents

   Example:
            perl $Script -mi_m ./*.mir2target.list -lnc_m ./novel_lncRNA_target.xls  -mi_lnc /lncRNA_target2mirna.list -in1 Lnc_Diff_Analysis -in2 Diff_Analysis/ -od  ./
----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
