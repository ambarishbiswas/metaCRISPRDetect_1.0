#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use Cwd 'abs_path';
use Term::ANSIColor;

use FindBin qw($Bin);




my $input_cluster_file="";

my $input_seq_file="";
my $output_seq_file="";

my $min_no_of_seq_in_a_cluster=10;
my $keep_seq_id=0;
my $preserve_cluster_index=0;
my $cluster_prefix="";
my $remove_identical_sequences=0;

my $identify_cluster_direction_from_the_sequence_header=0; #--- special case to handle CRISPRDetect predicted repeats
my $identify_cluster_direction_using_majority_rule=0;

#----------------------------------------------------
my $no_of_threads=5;

#use lib '/PROJECTS/CRISPRHost/lib1';
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager($no_of_threads);  

#-----------------------------------------------------

		

for(my $i=0;$i<=$#ARGV;$i++)
	{
						
				#--------------------------- input options ------------------------------------------------------------
				if($ARGV[$i]=~/-c$/)
					{
						$input_cluster_file=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-i$/)
					{
						$input_seq_file=$ARGV[$i+1];					
					}					
				elsif($ARGV[$i]=~/-o$/)
					{
						$output_seq_file=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-n$/)
					{
						$min_no_of_seq_in_a_cluster=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-k$/)
					{
						$keep_seq_id=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-p$/)
					{
						$preserve_cluster_index=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-cp$/)
					{
						$cluster_prefix=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-rm$/)
					{
						$remove_identical_sequences=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-identify_cluster_direction_from_the_sequence_header$/)
					{
						$identify_cluster_direction_from_the_sequence_header=$ARGV[$i+1];					
					}
				elsif($ARGV[$i]=~/-identify_cluster_direction_using_majority_rule$/)
					{
						$identify_cluster_direction_using_majority_rule=$ARGV[$i+1];					
					}							
	}

if($input_cluster_file!~/\S/ or $output_seq_file!~/\S/)
	{
		print "Error: please provide the input and out put file in this manner: -c input_cluster_file -i input_sequences_file -o output_file -n min_no_of_seq_in_a_cluster -p 1 -k 1 \n\n"; exit;
	}

			


		
#---- load seq IDs and sequences in a hash -------------------
my %hash_of_seq_id_and_seq;

my @arr_input_seq_file=`less $input_seq_file >&1`;
		
		
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


my @arr_clstr_file=`less $input_cluster_file >&1`;


my %hash_of_no_of_entries_and_cluster_id;		
my %hash_of_cluster_id_and_ordered_entries;
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
					
					$hash_of_no_of_entries_and_cluster_id{$no_of_seqs_in_cluster}{$cluster_id}=length($rep_seq);
					
					#--- first load the rep_seq_id and seq -----
					$hash_of_cluster_id_and_ordered_entries{$cluster_id}{'0'}="$rep_seq_id|$rep_seq";
					
					#---- now load the sequences that belong to this cluster in a hash
					my $seq_index=1;
					foreach my $seq_id(sort{length($hash_temp_seq_id_and_seq{$b})<=>length($hash_temp_seq_id_and_seq{$a})} keys %hash_temp_seq_id_and_seq)
						{
							my $t_seq=$hash_temp_seq_id_and_seq{$seq_id};
							$hash_of_cluster_id_and_ordered_entries{$cluster_id}{$seq_index}="$seq_id|$t_seq";
							
							#print ">$seq_id\n$t_seq\n";
							$seq_index++;
						}

			}
		
	}	
#close(SINGLE);




#---- now create the output file		
open(CR,">$output_seq_file");

#--- sort the cluster_ids based on no_of_entries
my @sorted_cluster_ids;
foreach my $nof_entries(sort{$a<=>$b} keys %hash_of_no_of_entries_and_cluster_id)
	{
		foreach my $cluster_id(sort{length($hash_of_no_of_entries_and_cluster_id{$nof_entries}{$b})<=>length($hash_of_no_of_entries_and_cluster_id{$nof_entries}{$a})} keys %{$hash_of_no_of_entries_and_cluster_id{$nof_entries}})
			{
				push(@sorted_cluster_ids,$cluster_id);
			}
	}


#foreach my $c_id(sort{$a<=>$b} keys %hash_of_cluster_id_and_ordered_entries)
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


################# subs

sub get_optimal_cluster_direction_using_majority_rule()
	{
		my($optimal_cluster_direction,$hash_cluster_lines)=@_;
		
		my %hash_of_directions;
		my $line_index=1;
		foreach my $current_line(keys %{$hash_cluster_lines})
			{
				my @arr_t1=split('>',$current_line);
				
				my @arr_t2=split('\.\.\.',$arr_t1[$#arr_t1]);
									
				my $seq_id=$arr_t2[0];
									
									
									
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
					
				if($hash_of_directions{$strand})
					{	
						$hash_of_directions{$strand}=$hash_of_directions{$strand}+1;	
					}
				else{
						$hash_of_directions{$strand}=1;
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



exit;
	
