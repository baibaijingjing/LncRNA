#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use experimental qw(smartmatch);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Times=localtime();
my $Year=$Times[5]+1990;
my $Month=$Times[4]+1;
my $Day=$Times[3];
#######################################################################################

sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   ziy <ziy\@biomarker.com.cn> 
Program Date:   2016.11.2
      Modify:   
 Description:   This program is used to ......
       Usage:   perl $Script -i filter_final.gff  -d database_lncRNA.gff -o ./
        Options:
        -i <file>   input gff file,forced

	-d <file>   input the species lncRNA database gff file,forced

        -o <dir>   output dir,optional, default[./]

        -h      help

USAGE
    print $usage;
    exit;
}

#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($infile,$outdir,$db_gff);

GetOptions(
                "help|?" =>\&USAGE,
                "o:s"=>\$outdir,
                "i:s"=>\$infile,
		"d:s"=>\$db_gff,
                ) or &USAGE;
&USAGE unless ($infile and $db_gff);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
$outdir ||='./';
`mkdir -p $outdir` unless (-d $outdir);
$infile=abs_path($infile);
$outdir=abs_path($outdir);
$db_gff=abs_path($db_gff);

my (%db_mrna,%db_exon);

#############################################
open (IN,$db_gff) or die $!;
while (<IN>) { 
    chomp;
    next if (/^\s*$|^\#/);
    my @lines=split /\t/,$_;
    die "Error:$db_gff please input useful lncRNA database! must be gff format\n" if ($lines[3]!~/^\d+$/ || $lines[4]!~/^\d+$/ || $lines[6]!~/^[\.\-\+]$/);
    if ($lines[2] eq 'mRNA' and $lines[-1]=~/ID=([^;\s]+).*Parent=([^;\s]+)/){
	push (@{$db_mrna{$lines[0]}{$lines[6]}},[$lines[3],$lines[4],$2,$1,$lines[1]])
    }elsif($lines[2] eq 'exon' and $lines[-1]=~/Parent=([^;\s]+)/){
	push (@{$db_exon{$1}},[$lines[3],$lines[4]]);
    }
}
close IN;

if (!%db_mrna or !%db_exon){
	die "Error:$db_gff please input useful lncRNA database! must be gff format\n";
}

##############################################
my (%ref_mrna,%ref_exon,%ref_db);

open (IN,$infile) or die $!;
while (<IN>){
	chomp;
	next if (/^\s*$|^\#/);
	my @lines=split /\t/,$_;
	die "Error:$infile have some problems,please check whether to gff format or other problems\n" if ($lines[3]!~/^\d+$/ || $lines[4]!~/^\d+$/ || $lines[6]!~/^[\.\-\+]$/);
	if ($lines[2] eq 'mRNA' and $lines[-1]=~/ID=([^;\s]+).*Parent=([^;\s]+)/){
		push  (@{$ref_mrna{$lines[0]}},[$lines[3],$lines[4],$lines[6],$1,$2]);
	}elsif($lines[2] eq 'exon' and $lines[-1]=~/Parent=([^;\s]+)/){
		push (@{$ref_exon{$1}},[$lines[3],$lines[4],$_]);
	}
}
close IN;

if (!%ref_mrna or !%ref_exon){
	die "Error:$infile have some problems,please check whether to gff format or other problems\n";	
}

foreach my $chr(sort keys %ref_mrna){
    	foreach my $m(sort{$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$ref_mrna{$chr}}){
	 	if (exists $db_mrna{$chr}{$m->[2]}){
			foreach my $dm(sort{$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$db_mrna{$chr}{$m->[2]}}){
			    	my $flag=0;
				my $score=0;
				last if ($dm->[0] > $m->[1]);
				if(($m->[0] >= $dm->[0] and $m->[0] <= $dm->[1]) || ($m->[0] < $dm->[0] and $m->[1] > $dm->[0])){
					if($m->[0] >= $dm->[0] and $m->[1] <= $dm->[1]){
						$score+=1;
						$score+=1 if (@{$db_exon{$dm->[3]}} == @{$ref_exon{$m->[3]}});
					}
					foreach my $e(sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$ref_exon{$m->[3]}}){
						foreach my $de(sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$db_exon{$dm->[3]}}){
							if (($e->[0] >= $de->[0] and $e->[0] <= $de->[1]) || ($e->[0] < $de->[0] and $e->[1] > $de->[0])){
								$score+=0.5;
								if($e->[0] >= $de->[0] and $e->[1] <= $de->[1]){
									$score+=0.5;
									$flag=1;
								}
						    	}
						}
					}
				}
				if ($flag==1){
					push (@{$ref_db{$m->[3]}},[$dm->[3],$score,$dm->[-1],$dm->[2]]);
				}
			}
		}
	}

}

########################################################
my (%uniq_mrna,%uniq_gene);
open (OUT1,">$outdir/lncRNA_known.gtf") or die $!;
open (OUT2,">$outdir/lncRNA_new.gtf") or die $!;
open (OUT3,">$outdir/lncRNA_mRNA_id.list") or die $!;
open (OUT4,">$outdir/lncRNA_gene_id.list") or die $!;

foreach my $chr(sort keys %ref_mrna){
        foreach my $m(sort{$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$ref_mrna{$chr}}){
		foreach my $e(sort{$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$ref_exon{$m->[3]}}){
			my @tmp=split /\t/,$e->[2];
			if (exists $ref_db{$m->[3]}){
				my (@db_from,@gene_alias,@transcript_alias,@transcript_score);
				my $flag=0;
				foreach my $t(sort{$b->[1] <=> $a->[1] || $b->[2] cmp $a->[2] || $a->[0] cmp $b->[0] || $a->[3] cmp $b->[3]} @{$ref_db{$m->[3]}}){
				    	push @db_from,$t->[2] if (!(@db_from~~/^$t->[2]$/));
					push @gene_alias,$t->[3] if (!(@gene_alias~~/^$t->[3]$/));
					if (!(@transcript_alias~~/^$t->[0]$/)){
						push @transcript_alias,$t->[0];
						push @transcript_score,$t->[1];
				        }
				}
				print OUT1 "$chr\t",join(",",@db_from),"\t",join("\t",@tmp[2..7]);
				print OUT1 "\tgene_id \"$m->[-1]\"; transcript_id \"$m->[3]\";";
				for(my $i=0;$i<@gene_alias;$i++){
				    	print OUT1 ' gene_alias_',$i+1," \"$gene_alias[$i]\";";
				}
				for(my $i=0;$i<@transcript_alias;$i++){
					print OUT1 ' transcript_alias_',$i+1," \"$transcript_alias[$i]\";";
				}
				for(my $i=0;$i<@transcript_alias;$i++){
					print OUT1 ' transcript_score_',$i+1," \"$transcript_score[$i]\";";
				}
				print OUT1 "\n";
				if (!exists $uniq_mrna{$m->[3]}{$transcript_alias[0]}){
					$uniq_mrna{$m->[3]}{$transcript_alias[0]}=1;
					print OUT3 "$m->[3]\t$transcript_alias[0]\n";
				}
				if (!exists $uniq_gene{$m->[-1]}{$gene_alias[0]}){
					$uniq_gene{$m->[-1]}{$gene_alias[0]}=1;
					print OUT4 "$m->[-1]\t$gene_alias[0]\n";
				}
					
			}else{
			    	print OUT2 join("\t",@tmp[0..7]),"\tgene_id \"$m->[-1]\"; transcript_id \"$m->[3]\";\n";
				if (!exists $uniq_mrna{$m->[3]}{'--'}){
					print OUT3 "$m->[3]\t--\n";
					$uniq_mrna{$m->[3]}{'--'}=1;
				}
				if (!exists $uniq_gene{$m->[-1]}{'--'}){
					print OUT4 "$m->[-1]\t--\n";
					$uniq_gene{$m->[-1]}{'--'}=1;
			        } 
			}
		}
	}
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
#######################################################################################
print STDOUT "\nperl $Script Done. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my $Time = &sub_format_datetime(localtime($time));
        print "$Time:\t$txt\n" ;
        return ($time) ;
}
#############################################################
#&run_or_die($cmd);
sub run_or_die()
{
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
        return ;
}
######################################################################
## qsub
sub qsub()
{
        my ($shfile, $queue, $ass_maxproc) = @_ ;
        $queue ||= 'general.q' ;
        $ass_maxproc ||= 20 ;
        if (`hostname` =~ /cluster/){
                my $cmd = "/share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $ass_maxproc --reqsub $shfile --queue $queue --independent" ;
                &run_or_die($cmd);
        }
        else{
                my $cmd = "/share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $ass_maxproc --reqsub $shfile --queue $queue --independent" ;
                &run_or_die($cmd);
        }

        return ;
}
#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
        sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
#######################################################################################
#                             .       .
#                            / `.   .' \
#                    .---.  <    > <    >  .---.
#                    |    \  \ - ~ ~ - /  /    |
#                     ~-..-~             ~-..-~
#                 \~~~\.'                    `./~~~/
#       .-~~^-.    \__/                        \__/
#     .'  O    \     /               /       \  \
#    (_____,    `._.'               |         }  \/~~~/
#     `----.          /       }     |        /    \__/
#           `-.      |       /      |       /      `. ,~~|
#               ~-.__|      /_ - ~ ^|      /- _      `..-'   f: f:
#                    |     /        |     /     ~-.     `-. _||_||_
#                    |_____|        |_____|         ~ - . _ _ _ _ _>
#

