#!/usr/bin/perl -w
package predict_one_duplex;

use strict;
use warnings;
use Exporter;
use Cwd;
use FindBin qw($Bin $Script);
use vars qw(@ISA @EXPORT @EXPORT_OK);
@ISA       = qw(Exporter);
@EXPORT    = qw(predict_one_duplex);
@EXPORT_OK = qw();

#$target, $miRNA
sub predict_one_duplex{
		my $target = shift;
		my $miRNA = shift;
		my $aln = shift;
		my $current_time = time;
		$current_time.=$aln;
		# -----------------------
		# ------- remove "-" from the 3' end of the miRNA
		# -----------------------
		if(substr($miRNA, 0, 1) eq "-")			# 3' end of miRNA duplex is '-'
		{
			my $have_dash = 1;
			
			# count how many '-' are their in the 3' end of miRNA
			my $dash_count = 0;
			while (substr($miRNA, $dash_count, 1) eq '-')
			{
				$dash_count++;
			}
			
			# remove "-" from the left of the duplex		
			$target = substr($target, $dash_count);
			$miRNA = substr($miRNA, $dash_count);
		}
		
		# -----------------------
		# -------  remove "-" from the 5' end of the miRNA
		# -----------------------
		if(substr($miRNA, length($miRNA)-1, 1) eq "-")			# 5' end of miRNA duplex is '-'
		{
			my $have_dash = 1;
			
			# count how many '-' are their in the 5' end of miRNA
			my $dash_count = 0;
			while (substr($miRNA, length($miRNA)-1-$dash_count, 1) eq '-')
			{
				$dash_count++;
			}
			
			# remove "-" from the right of the duplex		
			$target = substr($target, 0, length($target)-$dash_count);
			$miRNA = substr($miRNA, 0, length($miRNA)-$dash_count);
		}
		# -----------------------
		# -------  encode duplex
		# -----------------------
		my @duplex_code;
		for (my $j=0; $j<length($miRNA); $j++)
		{
			# encode the duplex_code first
			# AU | UA:	0
			# GC | CG:	1
			# GU | UG:	2
			# gap:		3
			# mismatch:	4
			my $tmp = substr($miRNA, length($miRNA)-1-$j, 1).substr($target, length($target)-1-$j, 1);
			if ($tmp eq "AU" || $tmp eq "UA")
			{
				$duplex_code[$j] = 0;
			}
			elsif($tmp eq "GC" || $tmp eq "CG")
			{	
				$duplex_code[$j] = 1;
			}
			elsif($tmp eq "GU" || $tmp eq "UG")
			{	
				$duplex_code[$j] = 2;
			}
			elsif(index($tmp, "-")!=-1)		# search for the first occurrence of "-"
			{	
				$duplex_code[$j] = 3;
			}
			else							# mismatch
			{
				$duplex_code[$j] = 4;
			}
		}
		
		# -----------------------
		# -------  compute length of seed region. seed region is defined as the first 9 nucleotides counting from 5' of miRNA in the duplex
		# -----------------------
		my $seed_length = 0;
		my $nt_count = 0;
		{
			my $j;
			for ($j=0; $j<length($miRNA); $j++)
			{
				if (substr($miRNA, length($miRNA)-1-$j, 1) ne "-")
				{
					$nt_count++;
				}
				if ($nt_count == 9)		# got 9 nucleotides, end of seed region reached
				{
					last;
				}
			}
			$seed_length = $j+1;		# $seed_length is the actual length, so is the array index + 1
		}
		
		# -----------------------
		# -------  compute pair compositions
		# -----------------------
		
		#-----1. entire ------#
		my $num_AU_entire = 0;
		my $num_GC_entire = 0;
		my $num_GU_entire = 0;
		my $num_gap_entire = 0;
		my $num_mismatch_entire = 0;
		
		for (my $j=0; $j<=$#duplex_code; $j++)
		{
			if ($duplex_code[$j] == 0)
			{
				$num_AU_entire++;
			}
			elsif($duplex_code[$j] == 1)
			{
				$num_GC_entire++;
			}
			elsif($duplex_code[$j] == 2)
			{
				$num_GU_entire++;
			}
			elsif($duplex_code[$j] == 3)
			{
				$num_gap_entire++;
			}
			else
			{
				$num_mismatch_entire++;
			}
		}
		
		#-----2. seed ------#
		my $num_AU_seed = 0;
		my $num_GC_seed = 0;
		my $num_GU_seed = 0;
		my $num_gap_seed = 0;
		my $num_mismatch_seed = 0;
		for (my $j=0; $j<$seed_length; $j++)
		{
			if ($duplex_code[$j] == 0)
			{
				$num_AU_seed++;
			}
			elsif($duplex_code[$j] == 1)
			{
				$num_GC_seed++;
			}
			elsif($duplex_code[$j] == 2)
			{
				$num_GU_seed++;
			}
			elsif($duplex_code[$j] == 3)
			{
				$num_gap_seed++;
			}
			else
			{
				$num_mismatch_seed++;
			}
		}
		
		#-----3. non_seed ------#
		my $num_AU_non_seed = 0;
		my $num_GC_non_seed = 0;
		my $num_GU_non_seed = 0;
		my $num_gap_non_seed = 0;
		my $num_mismatch_non_seed = 0;
		for (my $j=$seed_length; $j<=$#duplex_code; $j++)
		{
			if ($duplex_code[$j] == 0)
			{
				$num_AU_non_seed++;
			}
			elsif($duplex_code[$j] == 1)
			{
				$num_GC_non_seed++;
			}
			elsif($duplex_code[$j] == 2)
			{
				$num_GU_non_seed++;
			}
			elsif($duplex_code[$j] == 3)
			{
				$num_gap_non_seed++;
			}
			else
			{
				$num_mismatch_non_seed++;
			}
		}
		
		# -----------------------
		# -------  append pair compositions in the feature vector
		# -----------------------
			
		# append pair compositions in the feature vector
		my @feature_vector;
		my $index = 0;
		my $duplex_length = $#duplex_code+1;
		my $non_seed_length = $duplex_length - $seed_length;
		$feature_vector[$index++] = $num_AU_entire/$duplex_length;
		$feature_vector[$index++] = $num_GC_entire/$duplex_length;
		$feature_vector[$index++] = $num_GU_entire/$duplex_length;
		$feature_vector[$index++] = $num_gap_entire/$duplex_length;
		$feature_vector[$index++] = $num_mismatch_entire/$duplex_length;
		
		$feature_vector[$index++] = $num_AU_seed/$seed_length;
		$feature_vector[$index++] = $num_GC_seed/$seed_length;
		$feature_vector[$index++] = $num_GU_seed/$seed_length;
		$feature_vector[$index++] = $num_gap_seed/$seed_length;
		$feature_vector[$index++] = $num_mismatch_seed/$seed_length;
		
		$feature_vector[$index++] = $num_AU_non_seed/$non_seed_length;
		$feature_vector[$index++] = $num_GC_non_seed/$non_seed_length;
		$feature_vector[$index++] = $num_GU_non_seed/$non_seed_length;
		$feature_vector[$index++] = $num_gap_non_seed/$non_seed_length;
		$feature_vector[$index++] = $num_mismatch_non_seed/$non_seed_length;
		# -----------------------
		# -------  compute di-pair compositions
		# -----------------------
		
		#-----1. entire ------#
		# initialize the di_pair_composition 5*5 2D array with all elements 0
		my @di_pair_composition;
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$di_pair_composition[$j*5+$k] = 0;	# width 5
			}
		}
		
		# scan the duplex_code and compute the composition of each di-pair
		for (my $j=0; $j<=$#duplex_code-1; $j++)		# note $j can only reach the second last duplex_code
		{
			$di_pair_composition[$duplex_code[$j]*5 + $duplex_code[$j+1]]++;
		}
		
		# append to feature vector
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$feature_vector[$index++] = $di_pair_composition[$j*5 + $k]/($duplex_length-1);	# denorminator is duplex length - 1
			}
		}
		
		#-----2. seed ------#
		# initialize the di_pair_composition 2D array with all elements 0
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$di_pair_composition[$j*5+$k] = 0;
			}
		}
		
		# scan the seed region of the duplex_code and compute the composition of each di-pair
		for (my $j=0; $j<$seed_length-1; $j++)				# note that $seed_length is the actual length. $j can only reach the 2nd last
		{
			$di_pair_composition[$duplex_code[$j]*5+$duplex_code[$j+1]]++;
		}
		
		# append to feature vector
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$feature_vector[$index++] = $di_pair_composition[$j*5+$k]/($seed_length-1);	# denorminator is seed length - 1
			}
		}
		
		
		#-----3. non_seed ------#
		# initialize the di_pair_composition 2D array with all elements 0
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$di_pair_composition[$j*5+$k] = 0;
			}
		}
		
		# scan the non_seed region of duplex_code and compute the composition of each di-pair
		for (my $j=$seed_length; $j<=$#duplex_code-1; $j++)			# note $seed_length is the actual length, which is the starting index of the non_seed region
		{
			$di_pair_composition[$duplex_code[$j]*5+$duplex_code[$j+1]]++;
		}
		
		# append to feature vector
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				$feature_vector[$index++] = $di_pair_composition[$j*5+$k]/($duplex_length-$seed_length-1);	# denorminator is non-seed length - 1
			}
		}
			
		
		# -----------------------
		# -------  compute tri-pair compositions
		# -----------------------
		
		#-----1. entire ------#
		# initialize the tri_pair_composition 3D array with all elements 0
		my @tri_pair_composition;
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$tri_pair_composition[$j*25+$k*5+$p] = 0;
				}
			}
		}
		
		# scan the duplex_code and compute the composition of each tri-pair
		for (my $j=0; $j<=$#duplex_code-2; $j++)			# note $j can only reach the third last duplex_code
		{
			$tri_pair_composition[$duplex_code[$j]*25+$duplex_code[$j+1]*5+$duplex_code[$j+2]]++;
		}
		
		# append to feature vector
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$feature_vector[$index++] = $tri_pair_composition[$j*25+$k*5+$p]/($duplex_length-2);	# denorminator is duplex length - 2
				}
			}
		}
		
		#-----2. seed ------#
		# initialize the tri_pair_composition 3D array with all elements 0
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$tri_pair_composition[$j*25+$k*5+$p] = 0;
				}
			}
		}
		
		# scan the seed region of the duplex_code and compute the composition of each tri-pair
		for (my $j=0; $j<$seed_length-2; $j++)		# note that $seed_length is the actual length, which is 1 more than the index, and $j can only reach the 3nd last
		{
			$tri_pair_composition[$duplex_code[$j]*25+$duplex_code[$j+1]*5+$duplex_code[$j+2]]++;
		}
		
		# append to feature vecto
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$feature_vector[$index++] = $tri_pair_composition[$j*25+$k*5+$p]/($seed_length-2);	# denorminator is seed length - 2
				}
			}
		}
		
		
		#-----3. non_seed ------#
		# initialize the tri_pair_composition 3D array with all elements 0
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$tri_pair_composition[$j*25+$k*5+$p] = 0;
				}
			}
		}
		
		# scan the non_seed region of duplex_code and compute the composition of each tri-pair
		for (my $j=$seed_length; $j<=$#duplex_code-2; $j++)			# note $seed_length is the actual length, which is the starting index of the non_seed region
		{
			$tri_pair_composition[$duplex_code[$j]*25+$duplex_code[$j+1]*5+$duplex_code[$j+2]]++;
		}
		
		# append to feature vector
		for (my $j=0; $j<=4; $j++)
		{
			for (my $k=0; $k<=4; $k++)
			{
				for (my $p=0; $p<=4; $p++)
				{
					$feature_vector[$index++] = $tri_pair_composition[$j*25+$k*5+$p]/($duplex_length-$seed_length-2);	# denorminator is non-seed length - 2
				}
			}
		}
		
		# -----------------------
		# -------  get the informative features, based on $feature_f_score, 1 indicates f_score > 0.1, 0 otherwise
		# -----------------------
		my @feature_f_score = (1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1);
		my $k=0;
		my @feature_to_svm;
		for (my $j=0; $j<=$#feature_f_score; $j++)
		{
			if ($feature_f_score[$j] == 1)
			{
				$feature_to_svm[$k++] = $feature_vector[$j];
			}
		}
				
		
		# -----------------------
		# -------  construct the input file to svm
		# -----------------------
		my $seq_file_name = $current_time.'.seq';
		my $predict_file_name = $current_time.'.predict';
		
		open OUT, ">$seq_file_name" or die $!;
		
		print OUT "0 ";				# label 
		
		for (my $i=1; $i<=$#feature_to_svm+1; $i++)
		{
			print OUT $i;
			print OUT ":";
			print OUT $feature_to_svm[$i-1];	# feature
			print OUT " ";
		}
		print OUT "\n";		# this is important as SVM would otherwise complains
		
		close OUT;
		
		# -----------------------
		# -------  run svm
		# -----------------------
		my $shell_command = "$Bin/svm_classify ".$seq_file_name." $Bin/model_f_score_over_0.1.txt ".$predict_file_name.">/dev/null";
		my $shell_output = system($shell_command);
		
		
		# -----------------------
		# -------  read svm predict file, and report result
		# -----------------------
		open IN, "$predict_file_name" or die $!;
		my $result;
		my $prediction_score = <IN>;
                chomp $prediction_score;
		close IN;
		
		# -----------------------
		# -------  delete seq file and predict files
		# -----------------------
		$shell_command = 'rm '.$seq_file_name;
		system($shell_command);
		$shell_command = 'rm '.$predict_file_name;
		system($shell_command);
		return $prediction_score;
}

1;
