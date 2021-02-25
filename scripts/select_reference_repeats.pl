#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use Cwd 'abs_path';
use Term::ANSIColor;

use FindBin qw($Bin);

#---- Logic: Process the combined_crispr_file to identify CRISPR repeats with
#------- 1. There should be at least 5 unique CRISPR arrays associated to the repeat with clear boundaries.
#------- 2. CRISPR quality score >3
#------- 3. Not >3 mutations in the array, or < 25% of the repeats has mutation(s)



#my 
#our $refseq_or_metagenomic="REFSEQ";

#my $combined_crispr_file="$Bin/../REF_CRISPR_FILES/combined_archaea_and_bacterial_crisprs_from_metagenomes_v1.txt";
#our $refseq_or_metagenomic="METAGENOMIC";


#-

#--- define some parameters
our $min_no_of_seq_in_a_cluster=2;
our $keep_seq_id=0;
our $preserve_cluster_index=0;

our $remove_identical_sequences=1;

our $identify_cluster_direction_from_the_sequence_header=1; #--- special case to handle CRISPRDetect predicted repeats
our $identify_cluster_direction_using_majority_rule=0;
		



print "\nThis script requires specific files and folders. Check code before running.\n"; exit;



#--- start processing the combined array file and get the repeats

my @arr_crispr_sources;
push(@arr_crispr_sources,"REFSEQ");
push(@arr_crispr_sources,"METAGENOMIC");
push(@arr_crispr_sources,"REFSEQ_AND_METAGENOMIC");

foreach my $refseq_or_metagenomic(@arr_crispr_sources)
	{	
		print "Processing: $refseq_or_metagenomic ..\n";
		
		#--- define the source
		my $input_crispr_file;
		
		if($refseq_or_metagenomic eq "REFSEQ")
			{
				$input_crispr_file="$Bin/../REF_CRISPR_FILES/combined_archaea_and_bacterial_crisprs_from_refseq_86.txt";
			}
		elsif($refseq_or_metagenomic eq "METAGENOMIC")
			{
				$input_crispr_file="$Bin/../REF_CRISPR_FILES/combined_archaea_and_bacterial_crisprs_from_metagenomes_v1.txt";
			}
		elsif($refseq_or_metagenomic eq "REFSEQ_AND_METAGENOMIC")
			{
				$input_crispr_file="$Bin/../REF_CRISPR_FILES/combined_archaea_and_bacterial_crisprs_from_REFSEQ_AND_METAGENOMIC_v0.txt";
			}	
		#-----
		
		
		
		
		
		
		#--- specify file names to be created --	
		my $selected_repeat_sequences="$Bin/../REF_BLAST_DB_FILES/selected_repeats.$refseq_or_metagenomic.fa";
		my $blast_db_of_good_repeats="$Bin/../REF_BLAST_DB_FILES/REF_REPEATS.$refseq_or_metagenomic.db";
		
		
		#---Step1:  process the CRISPRDetect output file and create sequence file with good repeats --	
		&identify_good_repeats($refseq_or_metagenomic,$input_crispr_file,$selected_repeat_sequences);
		
			
		#--- Step2: run cd-hit to cluster the repeats --
		system("cd-hit-est -i $selected_repeat_sequences -o $selected_repeat_sequences.NR -n 3 -c 0.95 -aL 0.90 -aS 0.90 -M 0 -T 0 -d 0");
			
		#---- Step3: process the cluster file
		my $cluster_file="$selected_repeat_sequences.NR.clstr";
		our $cluster_prefix="REPEAT_$refseq_or_metagenomic\_";
		
		my %hash_of_no_of_entries_and_cluster_id;		
		my %hash_of_cluster_id_and_ordered_entries;			
		&process_cluster_file_and_load_hashes($cluster_prefix,$selected_repeat_sequences,$cluster_file,\%hash_of_no_of_entries_and_cluster_id,\%hash_of_cluster_id_and_ordered_entries);
		
		
		#--- sort the cluster_ids based on no_of_entries
		my @sorted_cluster_ids;
		foreach my $nof_entries(sort{$b<=>$a} keys %hash_of_no_of_entries_and_cluster_id)
			{
				foreach my $cluster_id(sort{length($hash_of_no_of_entries_and_cluster_id{$nof_entries}{$b})<=>length($hash_of_no_of_entries_and_cluster_id{$nof_entries}{$a})} keys %{$hash_of_no_of_entries_and_cluster_id{$nof_entries}})
					{
						push(@sorted_cluster_ids,$cluster_id);
					}
			}	
					
		
		
		#---- Step4: now create the db file		
		open(CR,">$blast_db_of_good_repeats");
		
		my %hash_of_already_written_sequences;
		my $new_c_id=0;
		foreach my $c_id(@sorted_cluster_ids)
			{		
				$new_c_id++;
				
				#-- get the cluster_size
				my $cluster_size=keys %{$hash_of_cluster_id_and_ordered_entries{$c_id}};
				foreach my $seq_index(sort{$a<=>$b} keys %{$hash_of_cluster_id_and_ordered_entries{$c_id}})
					{
						my @arr_t1=split('\|',$hash_of_cluster_id_and_ordered_entries{$c_id}{$seq_index});
						my $seq=$arr_t1[1];
						
						if($remove_identical_sequences >0 and $hash_of_already_written_sequences{$seq})
							{
								next;
							}
						
										
						my $seq_id;
						
						if($preserve_cluster_index !=0)
							{
								$seq_id=$cluster_prefix."CLUSTER_".$new_c_id."-SIZE_".$cluster_size."-S_".$seq_index;
							}
						else{
								$seq_id=$cluster_prefix."CLUSTER_".$c_id."-SIZE_".$cluster_size."-S_".$seq_index;
							}	
						
						if($keep_seq_id!=0)
							{
								$seq_id="$seq_id|$arr_t1[0]";
							}
							
						
						
						print CR ">$seq_id\n$seq\n";
						
						$hash_of_already_written_sequences{$seq}=1;
					}
				#print "\n";	
			}		
		close(CR);	
		
		#--- now create the blast DB
		system("makeblastdb -in $blast_db_of_good_repeats -dbtype nucl -parse_seqids");
	}

#--- cleanup some files --
system("rm $Bin/../REF_BLAST_DB_FILES/*.NR*");


############# subs 


sub identify_good_repeats()
	{
		my($refseq_or_metagenomic,$input_crispr_file,$selected_repeat_sequences)=@_;
		
		#--- process the CRISPRs
		open(RD,"$input_crispr_file");
		my @arr_crispr_file=<RD>;
		close(RD);
				
				
		#--- now load all repeats in a hash with associated parameters
		my %hash_of_repeat_index_and_repeat_seq_with_parameters;
		 
		
		open(WR,">$selected_repeat_sequences");	
		my $repeat_index=0;	
		for(my $i=0;$i<=$#arr_crispr_file;$i++)
			{
				my $line=$arr_crispr_file[$i]; chomp $line; 						
						
						
				if($line=~ /^>/)     #------ first level--------------------------------------------------------
					{
								
						my $array_start_line_index=$i;
						my @arr_current_array;
						
						my $is_questionable_array=0;
						my $k=$i;
						while($arr_crispr_file[$k]!~/\/\//)
							{
								#--- check to skip "# Questionable array : YES"
								if($line=~/# Questionable array : YES/)
									{
										$is_questionable_array++;
									}
								#----	
								push(@arr_current_array,$arr_crispr_file[$k]);
								$k++;
							}
						$i=$k;
						
						
						
						#--- process the current array and get parameters --
						if($is_questionable_array==0)
						{
							my($model_repeat,$array_direction,$directional_confidence,$score,$nof_repeats,$nof_perfect_repeats,$nof_mutated_repeats,$total_insertions,$total_mutations)=&process_array_and_load_hashes(\@arr_current_array);
							
							#print "$model_repeat\t$array_direction\t$directional_confidence\t$score\t$nof_repeats\t$total_insertions\t$total_mutations\t$nof_mutated_repeats\n";
							
							my $bad_array=0;
							
							if($score < 3)
								{
									$bad_array++;
								}
							if($nof_repeats < 3)
								{
									$bad_array++;
								}
							if($nof_perfect_repeats<3)
								{
									$bad_array++;
								}	
							if(length($model_repeat) < 23)
								{
									$bad_array++;
								}	
							if($total_insertions >3 and ($nof_mutated_repeats/$nof_repeats) >0.25 )
								{
									$bad_array++;
								}
							if($total_mutations > 3 and ($nof_mutated_repeats/$nof_repeats) >0.25)
								{
									$bad_array++;
								}
							
							
							#---- if all conditions are met, write the repeat sequence
							if($bad_array==0)
								{
									$repeat_index++;
									
									print WR ">REPEAT_$refseq_or_metagenomic\_$repeat_index-$array_direction-$directional_confidence-$score-$nof_repeats-$total_insertions-$total_mutations-$nof_mutated_repeats\n$model_repeat\n";
								}
						}		
					}	
								
			}
		close(WR);
		
		return 1;
	}


sub process_array_and_load_hashes()
	{
		my($current_array)=@_;
		
		
		#my $array_start; my $array_stop;my $array_seq="";my $crispr_index; my $left_flank=""; my $right_flank="";
		my $model_repeat="NA";
		my $array_direction="NA";
		my $directional_confidence="NA";
		my $score=0;
		my $nof_repeats=0;			
		my $total_insertions=0;
		my $total_mutations=0;		
		my $nof_mutated_repeats=0;
		my $nof_perfect_repeats=0;
		
		for(my $i=0;$i<=$#{$current_array};$i++)
			{
				my $line=$$current_array[$i]; chomp $line; 						
					
				#---- now get the direction ------------------------------------------------------------------------------
				if($line=~/Array_Orientation: (\S+)/)
					{
						$array_direction=$1; chomp $array_direction; $array_direction=~s/\r+//g;
						my @arr_d1=split('',$array_direction); #Forward|Reverse|Unconfirmed
						
						$array_direction=$arr_d1[0];
					}	
				#--- get the representative repeat ---
				if($line=~/# Primary repeat :     (\S+)/)
					{
						$model_repeat=$1; chomp $model_repeat; $model_repeat=~s/\r+//g;
					}	
				#--- get the array score ---
				if($line=~/# Questionable array : NO/)
					{
						my @arr_s1=split('\s+',$line);
						$score=$arr_s1[$#arr_s1];
					}
				if($line=~/Final direction:/)		##       Final direction:         R [0,1.15   Confidence: HIGH] 
					{
						my @arr_s1=split('\s+',$line);
						$directional_confidence=$arr_s1[$#arr_s1]; $directional_confidence=~s/\]//g;
						
					}
				#-- now check the individual repeats							
				if($line=~ /^>/)     
					{		
													
						#---now process all the dotted repeat, and add to sequence_strand on the go
						my $j=$i+4;
						while($$current_array[$j]!~/====/)
							{
								my $current_line=$$current_array[$j]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;
										
								my @arr_t1=split('\t',$current_line);
										
										my $r_start=$arr_t1[0];$r_start=~s/\s+//g;
										
										my $r_length=$arr_t1[1];$r_length=~s/\s+//g;
										my $s_length=$arr_t1[3];$s_length=~s/\s+//g;
										my $r_seq=$arr_t1[4]; $r_seq=~s/\s+//g; 
										my $s_seq=$arr_t1[5]; $s_seq=~s/\s+//g;
										
										my $comment=$arr_t1[6];		$s_length=~s/^\s+//;							
										
										
										
										
										my $tmp_r_seq=$r_seq;
										#--- get total mutations (includes deletions) ---
										my $nof_mutations=0;
										$tmp_r_seq=~s/\.//g;
										if($tmp_r_seq and $tmp_r_seq=~/\S+/)
											{
												$nof_mutations=length($tmp_r_seq);
											}
										
										
												
										#---- dont forget to check insertion bases in comment ---------------------------										
										my $no_of_insertions=0;
										if($arr_t1[6] and $arr_t1[6]!~/^Del/)
											{
												$comment=$arr_t1[6]; chomp $comment; $comment=~s/^\s+//;
															#if($comment=~/^Del/){next;}												
												my @tmp_arr1=split(' ',$comment);
												my $insertion_bases=$tmp_arr1[0];
												my $insertion_positions=$tmp_arr1[1];
												
												$insertion_bases=~s/,//g;
												$no_of_insertions=length($insertion_bases);	
												
											}
										
										#---------------------------------------------------------
										$total_mutations= $total_mutations + $nof_mutations;
										$total_insertions=$total_insertions+$no_of_insertions;
										
										#--- get the nof_repeats
										$nof_repeats++;
										
										#--- get the nof_mutated_repeats ----
										if($nof_mutations >0 or $no_of_insertions >0)
											{
												$nof_mutated_repeats++;
											}
										else{
												$nof_perfect_repeats++;
											}	
										
										$j++;
										
										if($j>=$#{$current_array}){last;}
									}
								$i=$j;
					}

							
			}
			
			
		return ($model_repeat,$array_direction,$directional_confidence,$score,$nof_repeats,$nof_perfect_repeats,$nof_mutated_repeats,$total_insertions,$total_mutations);
	}


sub process_cluster_file_and_load_hashes()
	{
		my($cluster_prefix,$selected_repeat_sequences,$cluster_file,$hash_of_no_of_entries_and_cluster_id,$hash_of_cluster_id_and_ordered_entries)=@_;
		
		#---- load seq IDs and sequences in a hash -------------------
		my %hash_of_seq_id_and_seq;
		my @arr_input_seq_file=`cat $selected_repeat_sequences >&1`;
				
		for(my $i=0;$i<$#arr_input_seq_file;$i++)
			{
				if($arr_input_seq_file[$i]=~/^>/)
					{
						my $id=$arr_input_seq_file[$i]; 	chomp $id;	$id=~s/\r//g; $id=~s/>//;
						my $seq=$arr_input_seq_file[$i+1];	chomp $seq;	$seq=~s/\r//g;
								
						$hash_of_seq_id_and_seq{$id}=$seq;
					}
			}	


		#---- now process the cluster file
		my @arr_clstr_file=`cat $selected_repeat_sequences.NR.clstr >&1`;
		
		
		my $cluster_index=0;
		#open(SINGLE,">$output_seq_file.single");		
		for(my $i=0;$i<$#arr_clstr_file;$i++)
			{
				my $c_line=$arr_clstr_file[$i];
				
				if($c_line=~/>Cluster/)
					{							
							
							
							# if identify_cluster_direction_from_the_sequence_header is 1; identify the optimal direction for the entire cluster
							my $optimal_cluster_direction="plus";
							if($identify_cluster_direction_from_the_sequence_header==1 or $identify_cluster_direction_using_majority_rule==1)
								{
									#--- load the lines belongs to this cluster and get the direction
									my $c=$i+1;
									my $no_of_seqs_in_cluster=0;
									
									
															
									my %hash_cluster_lines;
									my %hash_of_current_seqIDs;
									while($arr_clstr_file[$c] !~ />Cluster/)
										{
											my $current_line=$arr_clstr_file[$c]; chomp $current_line;$current_line=~s/\r//g;
											
											
											#--- get the rep_seq for this cluster --------------
											if($current_line=~/>/)
												{									
													
													my @arr_t1=split('>',$current_line);
													my @arr_t2=split('\.\.\.',$arr_t1[$#arr_t1]);											
													my $seq_id=$arr_t2[0];
													
													$hash_of_current_seqIDs{$seq_id}=1;
													
													$hash_cluster_lines{$current_line}=1;											
													$no_of_seqs_in_cluster++;
												}
											$c++;
											
											if($c > $#arr_clstr_file){last;}									
										}
									
									
									if($no_of_seqs_in_cluster < $min_no_of_seq_in_a_cluster) #-- as the clusters with just 2 sequences (one with HIGH confidence) will work
										{
											#--- write the leftover repeat in the single file
											#foreach my $seq_id(keys %hash_temp_seq_id_and_seq)
											#	{
											#		my $seq=$hash_of_seq_id_and_seq{$seq_id};	
											#		print SINGLE ">$seq_id\n$seq\n";
											#	}
											
											next;
										}
									
									#------
									if($no_of_seqs_in_cluster >=2 && $identify_cluster_direction_from_the_sequence_header==1) #-- as the clusters with just 2 sequences (one with HIGH confidence) will work
										{
											($optimal_cluster_direction)=&get_optimal_cluster_direction_from_the_sequence_header($optimal_cluster_direction,\%hash_cluster_lines);
											
											#if($optimal_cluster_direction eq "minus")
											#	{
													#print "All sequences from this cluster will be changed to $optimal_cluster_direction\n";
													#foreach my $l(keys %hash_cluster_lines)
													#	{
													#		print "$l\n";
													#	}
											#	}
										}
									elsif($no_of_seqs_in_cluster >=3 && $identify_cluster_direction_using_majority_rule==1)
										{
											($optimal_cluster_direction)=&get_optimal_cluster_direction_using_majority_rule($optimal_cluster_direction,\%hash_cluster_lines);
											
											#if($optimal_cluster_direction eq "minus")
											#	{
											#		#print "All sequences from this cluster [ $c_line ] will be changed to $optimal_cluster_direction\n";
											#	}
										}
										
										
								}
							
							#---- get no. of sequences in the cluster ---
							my $c=$i+1;
							my $no_of_seqs_in_cluster=0;
							
							my $rep_seq_id="";
							my $rep_seq="";
							
							my %hash_temp_seq_id_and_seq;
							while($arr_clstr_file[$c]!~/>Cluster/)
								{
									my $current_line=$arr_clstr_file[$c]; chomp $current_line;$current_line=~s/\r//g;
									
									
									#--- get the rep_seq for this cluster --------------
									if($current_line=~/>/)
										{									
											
											my @arr_t1=split('>',$current_line);
											my @arr_t2=split('\.\.\.',$arr_t1[$#arr_t1]);
											
											my $seq_id=$arr_t2[0];
											my $seq=$hash_of_seq_id_and_seq{$seq_id};									
											
											
											
											#---find the strand ---
											my $strand="plus";
											if($current_line=~/ -\//)
												{
													$strand='minus';
												}	
											
											#---- compare the current strand against the optimal direction
											if($strand ne $optimal_cluster_direction)
												{
													$seq=reverse($seq);$seq=~tr/ACGT/TGCA/;
												}
												
											if($current_line=~/ \*/)
												{
													$rep_seq_id=$seq_id;
													$rep_seq=$seq;
												}
											else{
													$hash_temp_seq_id_and_seq{$seq_id}=$seq;
												}		
																						
											#print "\t$seq\t[Direction: $strand]\n";									
										}
									#---------------------------------------------------
									if($current_line=~/>/)
										{
											$no_of_seqs_in_cluster++;
										}
									$c++;
									
									if($c > $#arr_clstr_file){last;}	
									
								}
							
											

							#----------------------------------------------------------------	
							
							
							#my $cluster_id="CLUSTER_".$cluster_index;
							$cluster_index++;
							my $cluster_id=$cluster_index;
							
							$hash_of_no_of_entries_and_cluster_id->{$no_of_seqs_in_cluster}->{$cluster_id}=length($rep_seq);
							
							#--- first load the rep_seq_id and seq -----
							$hash_of_cluster_id_and_ordered_entries->{$cluster_id}->{'0'}="$rep_seq_id|$rep_seq";
							
							#---- now load the sequences that belong to this cluster in a hash
							my $seq_index=1;
							foreach my $seq_id(sort{length($hash_temp_seq_id_and_seq{$b})<=>length($hash_temp_seq_id_and_seq{$a})} keys %hash_temp_seq_id_and_seq)
								{
									my $t_seq=$hash_temp_seq_id_and_seq{$seq_id};
									$hash_of_cluster_id_and_ordered_entries->{$cluster_id}->{$seq_index}="$seq_id|$t_seq";
									
									#print ">$seq_id\n$t_seq\n";
									$seq_index++;
								}

					}
				
			}	
		#close(SINGLE);




	
		
		return 1;
	}




sub get_optimal_cluster_direction_from_the_sequence_header()
	{
		my($optimal_cluster_direction,$hash_cluster_lines)=@_;
		
		#--- load the sensitivity values --
		my %hash_of_sensitivity_values;
		$hash_of_sensitivity_values{'HIGH'}=0.75;
		$hash_of_sensitivity_values{'MEDIUM'}=0.50;
		$hash_of_sensitivity_values{'LOW'}=0.25;
		
		
		#----
		my %hash_of_directions;
		
		
		my $line_index=1;
		foreach my $current_line(keys %{$hash_cluster_lines})
			{
				my @arr_t1=split('>',$current_line);
				
				my @arr_t2=split('\.\.\.',$arr_t1[$#arr_t1]);
									
				my $seq_id=$arr_t2[0];
									
				#----					
									
				
				#---find the strand ---
				my $strand="plus";
				if($current_line=~/ \*/)
					{
						$strand="plus";
					}
				if($current_line=~/ -\//)
					{
						$strand='minus';
				
					}
				#----	
				my $prediction_sensitivity="NA";
				if($current_line=~m/-(HIGH|MEDIUM|LOW|NA)-/) # >REPEAT_90600-R-HIGH-5.86-5-0-3-2
					{						
						$prediction_sensitivity=$1;
					}
				
				#--- get the score 
				my $sensitivity_score_for_this_direction=0;
				if($hash_of_sensitivity_values{$prediction_sensitivity})
					{
						$sensitivity_score_for_this_direction=$hash_of_sensitivity_values{$prediction_sensitivity};
					}
				
					
				if($hash_of_directions{$strand})
					{	
						$hash_of_directions{$strand}=$hash_of_directions{$strand}+$sensitivity_score_for_this_direction;	
					}
				else{
						$hash_of_directions{$strand}=$sensitivity_score_for_this_direction;
					}
						
				$line_index++;	
			}
		
		#-- now check which direction is higher
		if($hash_of_directions{'minus'} and $hash_of_directions{'minus'} >$hash_of_directions{'plus'})
			{
				$optimal_cluster_direction="minus";
			}
		else{
				$optimal_cluster_direction="plus";
			}	
		
		return ($optimal_cluster_direction);
	}





