use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my $version="v1.0.2";
my %config=%{readconf("$Bin/../../../../../../Config/lncRNA_pip.cfg")};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$k,$out);
GetOptions(
    "h|?" =>\&USAGE,
    "out|o:s"=>\$out,
    "fa|i:s"=>\$fa,
    "key|k:s"=>\$k,
    ) or &USAGE;
&USAGE unless ($fa and $k and $out );
$fa=&ABSOLUTE_DIR($fa);
&MKDIR("$out/");
$out=&ABSOLUTE_DIR($out);
#&MKDIR("$out/Known_LncRNA");
my $outdir="$out";
#$gff=&ABSOLUTE_DIR($gff);
#mkdir $out unless -d $out;
my $BLASTN=$config{blastn};
my $DBDIR=$config{NONCODE};
my %hash;
&get_dafa_ha($fa,\%hash);
open OUT,">$outdir/known_lncRNA.result.txt";
print OUT "query_id\tdb_id\tidentify\tcoverage\tdatabase\tquery_length\tdb_length\tmatch_length\tmis-match\tgap\tquery_star\tquery_end\tdb_star\tdb_end\te-value\tscore\n";

if ($k){
    if ($k=~/,/){
        my @key=split /,/,$k;
        foreach my $s(@key){
	    print "$s\n";
            my $db=$DBDIR.$s."_lncRNA.fa";
            `$BLASTN -evalue 1e-5 -db $db -query $fa -num_threads 2 -out $outdir/$s.known.blast.txt -outfmt 6`;
	    my%ha;
	    &get_dafa_ha($db,\%ha);
	    open IN1,"$outdir/$s.known.blast.txt";
	    $/="\n";
	    while (<IN1>){
		chomp;
		
		my @b=split /\s+/,$_;
		if ($hash{$b[0]}){
		    my $lnc=$hash{$b[0]}{'len'};
		    my $dbl=$ha{$b[1]}{'len'};
		    #my $lengt=&min($lnc,$dbl);
			my $lengt=$lnc;
		    my $length=$lengt-$b[5];
		    my $cov=sprintf "%0.2f",(int($b[3])/$length);
		    
		    if ($cov > 0.8){
			if ($b[2]>80){
			    my $s1=join"\t",@b[0..2];
			    my $s2=join"\t",@b[4..$#b];
			    print OUT "$s1\t$cov\t$s\t$lnc\t$dbl\t$b[3]\t$s2\n";
			}
			
		    }
		   
		}
		
	    }
	    
        }
    
    }else{
        my $db=$DBDIR.$k."_lncRNA.fa";
        `$BLASTN -evalue 1e-5 -db $db -query $fa -num_threads 2 -out $outdir/$k.known.blast.txt -outfmt 6`;
	my%ha;
            &get_dafa_ha($db,\%ha);
            open IN1,"$outdir/$k.known.blast.txt";
            $/="\n";
            while (<IN1>){
                chomp;

                my @b=split /\s+/,$_;
                if ($hash{$b[0]}){
                    my $lnc=$hash{$b[0]}{'len'};
                    my $dbl=$ha{$b[1]}{'len'};
                    my $lengt=$lnc;
                    my $length=$lengt-$b[5];
                    my $cov=sprintf "%0.2f",(int($b[3])/$length);
                    if ($cov > 0.8){
                        if ($b[2]>80){
                            my $s1=join"\t",@b[0..2];
                            my $s2=join"\t",@b[4..$#b];
                            print OUT "$s1\t$cov\t$k\t$lnc\t$dbl\t$b[3]\t$s2\n";
                        }
                        
                    }
                   
                }
                
            }
    }
}else{
    warn "the known lncRNA skip!"
}

sub get_lncfa_ha {
    my ($file,$hash)=@_;
    open IN, "$file";
    $/=">";
    #my %hash;
    <IN>;
    while (<IN>){
        chomp;
        my ($head,$fa)=(split /\n/,$_)[0,1];
        my $id=(split /\s+/,$head)[0];
        my $lous=(split /\s+/,$head)[1].":".(split /\s+/,$head)[2];
        my @tem=split //,$fa;
        my $len=@tem;
        #print "$id\t$len\n";
        $hash->{$id}->{'lous'}=$lous;
        $hash->{$id}->{'len'}=$len;
    }
}


sub get_dafa_ha {
    my ($file,$ha)=@_;
    open IN2, "$file";
    $/=">";
    #my %ha;
    <IN2>;
    while (<IN2>){
        chomp;
        my ($head,$fa);
	$head=(split /\n/,$_)[0];
	$fa=(split /\n/,$_)[1];
        my $id=(split /\s+/,$head)[0];
        #my $lous=(split /\s+/,$head)[1].":".(split /\s+/,$head)[2];
        my @tem=split //,$fa;
        my $len=@tem;
        #print "$id\t$len\n";
        #$ha{$id}{'lous'}=$lous;
        $ha->{$id}->{'len'}=$len;
    }
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
sub min{
    my $min=shift;
    my $temp;
    while (@_) {
        $temp=shift;
        $min=$min<$temp?$min:$temp;
    }
    return $min;
}

sub MKDIR{ # &MKDIR($out_dir);
    my ($dir)=@_;
#rmdir($dir) if(-d $dir);
    mkdir($dir) if(!-d $dir);
}
sub USAGE {
    my $usage="USAGE;
#-------------------------------------------------
Contact: lgniu
   Data: 2015-07-12
Fuction: the script is used to ...
  Usage:
    -fa  <STR>    lncRNA fasta file.
    -key  <STR>    noncode_fa_database  name eg.noncode,lncrnadb 
    -out    <STR>   output
    -h  help.
Example:
    perl Script -fa lnc_filter_final.fa -key noncode ,lncrnadb  -out output

#-------------------------------------------------
USAGE";
    print $usage;
    exit;
};
