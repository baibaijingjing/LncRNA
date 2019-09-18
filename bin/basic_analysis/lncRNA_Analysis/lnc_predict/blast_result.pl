my $odir="/share/nas31/niulg/test/blast/blastn";
my $key="osa";
my (%hp_lnc,%lnc_hp);
&blast_result_do("$odir/precursor.txt",\%hp_lnc,\%lnc_hp);
open OUT1,">$odir/$key.precursor.list\n" || die $!;
#open OUT1,">$odir/$key.target2mir.list\n" || die $!;
foreach my $tcon (sort keys %lnc_hp) {
	print OUT1 "$tcon\t";
        my $tcon_mir_num=keys%{lnc_hp{$tcon}};
 	print OUT1 "$tcon_mir_num\t";
	foreach my $hp (keys %{lnc_hp{$tcon}}){
		my $num=${lnc_hp{$tcon}}{$hp};
		print OUT1 "$hp\t$num;";
	}
        print OUT1 "\n";
	
}


close OUT1;

$mir_num=keys %mir_tar;
$target_gene_num=keys %target_gene;	



sub blast_result_do{

        my ($txt,$hp_lnc,$lnc_hp)=@_;

        open IN,"$txt" || die $!;
        while (<IN>) {
                chomp;
                next if (/^$/);
                $_=~s/^>//;
                my @tmp=split/\s+/,$_;
		if (exists $hp_lnc{$tmp[0]}){
			if (exists $hp_lnc{$tmp[0]}{$tmp[1]}){
				$hp_lnc->{$tmp[0]}->{$tmp[1]}++;	
			}else{$hp_lnc->{$tmp[0]}->{$tmp[1]}=1;}
			 
		}else{
			$hp_lnc->{$tmp[0]}->{$tmp[1]}=1
		}
                





		if (exists $lnc_hp{$tmp[1]}){
                         if (exists $lnc_hp{$tmp[1]}{$tmp[0]}){
                                 $lnc_hp->{$tmp[1]}->{$tmp[0]}++;
                         }else{$lnc_hp->{$tmp[1]}->{$tmp[0]}=1;}
                      
                 }else{
                         $lnc_hp->{$tmp[1]}->{$tmp[0]}=1
                 }
                 
               
        }
        close IN;
}

