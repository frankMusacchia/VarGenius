
#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2017>  <Francesco Musacchia>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
package LIB::picard_management;
## picard_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of the Picard tools for the pipeline

BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_PICARD_AddOrReplaceReadGroups run_PICARD_MarkDuplicates_gen run_PICARD_MergeSamFiles_gen
									run_PICARD_MarkDuplicates_chrom run_PICARD_MergeSamFiles_chrom run_PICARD_BedToIntervalList
									run_PICARD_BuildBamIndex run_PICARD_SortSam run_PICARD_AddOrReplaceReadGroups
									run_PICARD_CreateSequenceDictionary run_PICARD_ReorderSam
									run_PICARD_AddOrReplaceReadGroups_gen run_PICARD_SortVcf
									run_PICARD_GatherVcfs);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit  execute_threads build_input_name_from_executed);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present extract_name);


#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);

#########################################################PICARD



=head2 run_PICARD_MergeSamFiles_chrom

 Title   : run_PICARD_MergeSamFiles_chrom
 Usage   : run_PICARD_MergeSamFiles_chrom(   );

 Function: Uses MergeSAmFiles from PICARD to merge together BAM files
						It works both if the files have been sub-divided per-chromosomes
						and if it is not.
						
						Takes the distinct read files for the sample given in input
 Returns : nothing
=cut
sub run_PICARD_MergeSamFiles_chrom{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $in_fold_suff = shift;
	my $out_fold_suff = shift;#Used to know the folder where to take the input bam files
	my $step = shift;
	my $step_needed = shift;
	my $steps_array = shift;
	my $steps_array2 = shift;
	my $perchrom = shift;#Execution per chromosome (0 or 1)
	
	my $tempFolder = $cfg_hash->{'scratch_f'};		
		
	my $main_name = $cfg_hash->{'db_readf_name'};		

	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);		
						    
	#Here we get the read file ids associated with this sample. They will be used for the input file names
	my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
						$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
							
	#Be sure the previous step has been performed checking in the samples table
	my $step_needed_done = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$step_needed,$cfg_hash->{'db_sample_id'},$sample_id);
						
	if ( $step_needed_done == 1){
		#my $step = 'mergesam_step';
		print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);
		
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
		
		#Program execution
		my $java_call = $cfg_hash->{'java_path'};
		if (defined $cfg_hash->{'javamem_MSF'}){
			$java_call .= " -".$cfg_hash->{'javamem_MSF'}." ";
		}
		$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
		
		my $picard_params = "";
		if ( defined $cfg_hash->{'picard_USE_THREADING'} ){
			$picard_params .= " USE_THREADING=".$cfg_hash->{'picard_USE_THREADING'};
		}

					    		
		#If it was performed per chromosome, read file have already been merged into the PrintReads
		#program. So we just need to get the BAM for each chromosome
		my $numExec = 1;
		if ( $perchrom == 1 ){
				#Getting all chromosomes
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));

				#In folder
				my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
				my $params;
				getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
								$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

				#Prepare output name for files for data recalibration
				my $recal_bam = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'print_reads_step'},
								$params->{$cfg_hash->{'db_sample_name'}},$steps_array2)."_".$cfg_hash->{'print_reads_step'}.".".$cfg_hash->{'bam_ext'};
		
				#Get a single bam file from separated
				my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
									$sample_name,$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};
																		
				#Here get the samples names to be merged
				my $samples_to_merge  = "";
								
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
					#Build i-th bam output name
					my $ith_in_bam = $outFolder."/".extract_name($recal_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					$samples_to_merge .= " I=$ith_in_bam ";				
				}

				my $command = "$java_call TMP_DIR=$tempFolder $samples_to_merge OUTPUT=$outFile $picard_params";
				print_and_log( "Executing command: $command\n",$log_file);
				try_exec_command($command) or die "Unable to execute command: $command\n";
		}
		#If it hasn't been performed per-chromosome: execute only one merge for all the samples
		else{
			#Here get the samples names to be merged
			my $samples_to_merge  = "";
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			my $params;
			foreach my $readf_id (sort keys %{$distinct_readf_ids}) {
				#print_and_log("Sample id is: $readf_id..\n",$log_file);#DEBUGCODE
				getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
									$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
				 #Build the input file name using as stop point the base recalibration step
				 my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};
				$samples_to_merge .= " I=$input_file ";
			}
			
			#Get a single bam file from separated
			my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
								$sample_name,$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};
			
			my $command = "$java_call TMP_DIR=$tempFolder $samples_to_merge OUTPUT=$outFile $picard_params";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";		
		}	
	}else{
			print_and_log( "WARNING: Cannot run $prog_used because $step_needed has not been run..\n",$log_file);
	}												
}

=head2 run_PICARD_MergeSamFiles_gen

 Title   : run_PICARD_MergeSamFiles_gen
 Usage   : run_PICARD_MergeSamFiles_gen(   );

 Function: run_PICARD_MergeSamFiles_gen.
						When running this program the sample Id will be always the sample
						name not separated per lanes
 Returns : nothing
=cut
sub run_PICARD_MergeSamFiles_gen{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $samples_to_merge = shift;
	my $outFile = shift;
	my $log_file = shift;
		
	my $tempFolder = $cfg_hash->{'scratch_f'};
		#print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);#DEBUGCODE
	
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_MSF'}){
		$java_call .= " -".$cfg_hash->{'javamem_MSF'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";
	if ( defined $cfg_hash->{'picard_USE_THREADING'} ){
		$picard_params .= " USE_THREADING=".$cfg_hash->{'picard_USE_THREADING'};
	}
		
	my $command = "$java_call TMP_DIR=$tempFolder $samples_to_merge OUTPUT=$outFile $picard_params";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";												
	
}


=head2 run_PICARD_MarkDuplicates_chrom 

 Title   : run_PICARD_MarkDuplicates_chrom
 Usage   : run_PICARD_MarkDuplicates_chrom(   );

 Function: run_PICARD_MarkDuplicates_chrom
			Is a tool from PICARD which permits to tag reads that are duplicate
			but coming from the same fragment. 
			Duplicates are identified as read pairs having identical 5' 
			positions (coordinate and strand) for both reads in a mate pair
			If REMOVE_DUPLICATES = true is used then the duplicates are not only
			flagged but also removed.
			If METRICS_FILE = filepath then a file with metrics about the duplicates
			is also given
			
 Returns : nothing
=cut
sub run_PICARD_MarkDuplicates_chrom{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $in_fold_suff = shift;
	my $out_fold_suff = shift;#Used to know the folder where to take the input bam files
	my $step = shift;
	my $prev_step = shift;
	my $steps_array_in = shift;
	my $steps_array_out = shift;
	my $perchrom = shift;#Execution per chromosome (0 or 1)
		
	my $main_name = $cfg_hash->{'db_readf_name'};	
	my $tempFolder = $cfg_hash->{'scratch_f'};
	#inFolder for the bam files
	my $inFolder = $cfg_hash->{$group_id.'_'.$in_fold_suff.'_out_f'};
	#Out folder for the bam files cleaned
	my $outFolder = $cfg_hash->{$group_id.'_'.$out_fold_suff.'_out_f'};
		
	print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);
	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);		
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	#Getting the previous step status
	my $prev_step_done = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$prev_step,
						$cfg_hash->{'db_sample_id'},$sample_id);
	if ( $prev_step_done == 1){
				    																			
		#Program execution
		my $java_call = $cfg_hash->{'java_path'};
		if (defined $cfg_hash->{'javamem_MD'}){
			$java_call .= " -".$cfg_hash->{'javamem_MD'}." ";
		}
		$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
		
		my $picard_params = "";
		if ( defined $cfg_hash->{'picard_VALIDATION_STRINGENCY'} ){
			$picard_params .= " VALIDATION_STRINGENCY=".$cfg_hash->{'picard_VALIDATION_STRINGENCY'};
		}
		if ( defined $cfg_hash->{'picard_CREATE_INDEX'} ){
			$picard_params .= " CREATE_INDEX=".$cfg_hash->{'picard_CREATE_INDEX'};
		}
		if ( defined $cfg_hash->{'picard_REMOVE_DUPLICATES'} ){
			$picard_params .= " REMOVE_DUPLICATES=".$cfg_hash->{'picard_REMOVE_DUPLICATES'};
		}
		if ( defined $cfg_hash->{'picard_METRICS_FILE'} ){
			$picard_params .= " METRICS_FILE=".$cfg_hash->{'picard_METRICS_FILE'};
		}
		if ( defined $cfg_hash->{'picard_ASSUME_SORTED'} ){
			$picard_params .= " ASSUME_SORTED=".$cfg_hash->{'picard_ASSUME_SORTED'};
		}	
		if ( defined $cfg_hash->{'picard_MAX_FILE_HANDLES_FOR_READ_ENDS_MAP'} ){
			$picard_params .= " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=".$cfg_hash->{'picard_MAX_FILE_HANDLES_FOR_READ_ENDS_MAP'};
		}

		#If it was performed per chromosome
		my $numExec = 1;
		if ( $perchrom == 1 ){
				##Using also Chromosomex X and Y
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
				#Get a single bam file from separated
				my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
									$sample_name,$steps_array_in)."_".$exec.".".$cfg_hash->{'bam_ext'};
				#Get a single bam file from separated
				my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
									$sample_name,$steps_array_out)."_".$step."_".$exec.".".$cfg_hash->{'bam_ext'};											
				#Checks if the SAM output exists and is non-zero
				if ( file_not_present($input_file)  > 0){ log_and_exit( "ERROR: Cannot find $input_file. Please check mrdup_groups task\n", $log_file)}
			
				#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
				#my $command = $java_call." TMP_DIR=".$cfg_hash->{$group_id.'_'.$out_fold_suff.'_out_f'}." $picard_params INPUT=$input_file OUTPUT=$outFile";		 
				my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
						
				print_and_log( "Executing command: $command\n",$log_file);
				#try_exec_command($command) or die "Unable to execute command: $command\n";
				push(@commands,$command);
			}
			execute_threads(\@commands);
		}else{
			#Get a single bam file from separated
			my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
								$sample_name,$steps_array_in).".".$cfg_hash->{'bam_ext'};
			#Get a single bam file from separated
			my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
								$sample_name,$steps_array_out)."_".$step.".".$cfg_hash->{'bam_ext'};											
			#Checks if the SAM output exists and is non-zero
			if ( file_not_present($input_file) > 0){ log_and_exit( "ERROR: Cannot find $input_file. Please check mrdup_groups task\n", $log_file)}
		
			#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
			#my $command = $java_call." TMP_DIR=".$cfg_hash->{$group_id.'_'.$out_fold_suff.'_out_f'}." $picard_params INPUT=$input_file OUTPUT=$outFile";		 
			my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";		
		}
	}else{
			print_and_log( "WARNING: Cannot run $prog_used because $prev_step has not been run..\n",$log_file);
	}											
	
}


=head2 run_PICARD_AddOrReplaceReadGroups 

 Title   : run_PICARD_AddOrReplaceReadGroups
 Usage   : run_PICARD_AddOrReplaceReadGroups(   );

 Function: run_PICARD_AddOrReplaceReadGroups
			Add Read Group to the alignmed bam
			
 Returns : nothing
=cut
sub run_PICARD_AddOrReplaceReadGroups_gen{

	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;	
	my $analysis_id = shift;	
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	
	my $main_name = $cfg_hash->{'db_readf_name'};	

	#print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);#DEBUGCODE

	#get the sample name						    	
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    				    																			
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_MD'}){
		$java_call .= " -".$cfg_hash->{'javamem_MD'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";
	if ( defined $cfg_hash->{'picard_SORT_ORDER'} ){
		$picard_params .= " SORT_ORDER=".$cfg_hash->{'picard_SORT_ORDER'};
	}
	if ( defined $cfg_hash->{'picard_RGID'} and  $cfg_hash->{'picard_RGID'} eq 'YES' ){
		$picard_params .= " RGID=".$sample_name;
	}
	if ( defined $cfg_hash->{'picard_RGLB'} and  $cfg_hash->{'picard_RGLB'} eq 'YES' ){
		$picard_params .= " RGLB=".$targetbed;
	}
	if ( defined $cfg_hash->{'picard_RGSM'} and  $cfg_hash->{'picard_RGSM'} eq 'YES' ){
		$picard_params .= " RGSM=".$sample_name;
	}		
	if ( defined $cfg_hash->{'picard_RGPU'} and  $cfg_hash->{'picard_RGPU'} eq 'YES' ){
		$picard_params .= " RGPU=".$analysis_id;
	}
	if ( defined $cfg_hash->{'picard_RGPL'} ){
		$picard_params .= " RGPL=".$cfg_hash->{'picard_RGPL'};
  }

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_PICARD_AddOrReplaceReadGroups 

 Title   : run_PICARD_AddOrReplaceReadGroups
 Usage   : run_PICARD_AddOrReplaceReadGroups(   );

 Function: run_PICARD_AddOrReplaceReadGroups
			Add Read Group to the alignmed bam
			
 Returns : nothing
=cut
sub run_PICARD_AddOrReplaceReadGroups{

	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;	
	my $analysis_id = shift;	
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	
	my $main_name = $cfg_hash->{'db_readf_name'};	

	#print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);#DEBUGCODE

	#Getting the sample id...
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_readf_id'},$readf_id);
	#..to get the sample name						    	
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	#Get flowcell and lane from the database for the read file
	my $flowcell = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_flowcell'},
						    $cfg_hash->{'db_readf_id'},$readf_id);
	my $lane = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_lane'},
						    $cfg_hash->{'db_readf_id'},$readf_id);
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    				    																			
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_MD'}){
		$java_call .= " -".$cfg_hash->{'javamem_MD'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";
	if ( defined $cfg_hash->{'picard_SORT_ORDER'} ){
		$picard_params .= " SORT_ORDER=".$cfg_hash->{'picard_SORT_ORDER'};
	}
	if ( defined $cfg_hash->{'picard_RGID'} and  $cfg_hash->{'picard_RGID'} eq 'YES' ){
		$picard_params .= " RGID=".$flowcell.".".$lane;
	}
	if ( defined $cfg_hash->{'picard_RGLB'} and  $cfg_hash->{'picard_RGLB'} eq 'YES' ){
		$picard_params .= " RGLB=".$targetbed;
	}
	if ( defined $cfg_hash->{'picard_RGSM'} and  $cfg_hash->{'picard_RGSM'} eq 'YES' ){
		$picard_params .= " RGSM=".$sample_name;
	}		
	if ( defined $cfg_hash->{'picard_RGPU'} and  $cfg_hash->{'picard_RGPU'} eq 'YES' ){
		$picard_params .= " RGPU=".$analysis_id;
	}
	if ( defined $cfg_hash->{'picard_RGPL'} ){
		$picard_params .= " RGPL=".$cfg_hash->{'picard_RGPL'};
  }

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 run_PICARD_MarkDuplicates_gen 

 Title   : run_PICARD_MarkDuplicates_gen
 Usage   : run_PICARD_MarkDuplicates_gen(   );

 Function: run_PICARD_MarkDuplicates_gen
			Is a tool from PICARD which permits to tag reads that are duplicate
			but coming from the same fragment. 
			Duplicates are identified as read pairs having identical 5' 
			positions (coordinate and strand) for both reads in a mate pair
			If REMOVE_DUPLICATES = true is used then the duplicates are not only
			flagged but also removed.
			If METRICS_FILE = filepath then a file with metrics about the duplicates
			is also given
			
 Returns : nothing
=cut
sub run_PICARD_MarkDuplicates_gen{

	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	
	my $main_name = $cfg_hash->{'db_readf_name'};	

	#print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);#DEBUGCODE

				    																			
		#Program execution
		my $java_call = $cfg_hash->{'java_path'};
		if (defined $cfg_hash->{'javamem_MD'}){
			$java_call .= " -".$cfg_hash->{'javamem_MD'}." ";
		}
		$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
		
		my $picard_params = "";
		if ( defined $cfg_hash->{'picard_VALIDATION_STRINGENCY'} ){
			$picard_params .= " VALIDATION_STRINGENCY=".$cfg_hash->{'picard_VALIDATION_STRINGENCY'};
		}
		if ( defined $cfg_hash->{'picard_CREATE_INDEX'} ){
			$picard_params .= " CREATE_INDEX=".$cfg_hash->{'picard_CREATE_INDEX'};
		}
		if ( defined $cfg_hash->{'picard_REMOVE_DUPLICATES'} ){
			$picard_params .= " REMOVE_DUPLICATES=".$cfg_hash->{'picard_REMOVE_DUPLICATES'};
		}
		if ( defined $cfg_hash->{'picard_METRICS_FILE'} ){
			$picard_params .= " METRICS_FILE=$input_file\_".$cfg_hash->{'picard_METRICS_FILE'};
		}
		if ( defined $cfg_hash->{'picard_ASSUME_SORTED'} ){
			$picard_params .= " ASSUME_SORTED=".$cfg_hash->{'picard_ASSUME_SORTED'};
		}	
		if ( defined $cfg_hash->{'picard_MAX_FILE_HANDLES_FOR_READ_ENDS_MAP'} ){
			$picard_params .= " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=".$cfg_hash->{'picard_MAX_FILE_HANDLES_FOR_READ_ENDS_MAP'};
		}

		#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
		my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
		print_and_log( "Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 run_PICARD_BedToIntervalList 

 Title   : run_PICARD_BedToIntervalList
 Usage   : run_PICARD_BedToIntervalList(   );

 Function: run_PICARD_BedToIntervalList
			Is a tool from PICARD which permits create the interval list file 
			from the BED file. This file is useful in other tools like GATK
			
 Returns : nothing
=cut
sub run_PICARD_BedToIntervalList{

	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $input_file = shift;
	my $outFile = shift;

	print_and_log( "Starting $prog_used using picard tools..\n ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_MD'}){
		$java_call .= " -".$cfg_hash->{'javamem_MD'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." $picard_params I=$input_file O=$outFile SD=".extract_name($cfg_hash->{'hum_ref'},'noext').".".$cfg_hash->{'dict_ext'};		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

#=head2 run_PICARD_SortVcf 

 #Title   : run_PICARD_SortVcf
 #Usage   : run_PICARD_SortVcf(   );

 #Function: run_PICARD_SortSam
			#Is a tool from PICARD which permits to sort a VCF file using
			#an external sequence dictionary. If no external dictionary is supplied,
			 #the VCF file headers of multiple inputs must have the same sequence dictionaries.
			
 #Returns : nothing
#=cut
#sub run_PICARD_SortVcf{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $input_files = shift;
	#my $outFile = shift;
	#my $log_file = shift;
	#my $dict = shift;
	
	#my $tempFolder = $cfg_hash->{'scratch_f'};	

	##print_and_log( "Starting $prog_used using PICARD tools.. ",$log_file);
																		
	##Program execution
	#my $java_call = $cfg_hash->{'java_path'};
	#if (defined $cfg_hash->{'javamem_SS'}){
		#$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	#}
	#$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	

	#my $picard_params = "";
	
	#my @input_files = split(",",$input_files);
	#foreach my $input_file (@input_files){
		#$picard_params .= " INPUT=$input_file";
	#} 
		
	#if ( defined $dict ){
		#$picard_params .= " SEQUENCE_DICTIONARY=$dict";
	#}

	##print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	#my $command = $java_call." TMP_DIR=$tempFolder $picard_params OUTPUT=$outFile";		 
	#print_and_log( "Executing command: $command\n",$log_file);
	#try_exec_command($command) or die "Unable to execute command: $command\n";
#}

=head2 run_PICARD_ReorderSam 

 Title   : run_PICARD_ReorderSam
 Usage   : run_PICARD_ReorderSam(   );

 Function: run_PICARD_ReorderSam
			Is a tool from PICARD which permits to reorder a bam or a sam file 
			to match a given reference 
			
 Returns : nothing
=cut
sub run_PICARD_ReorderSam{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $reference = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	#print_and_log( "Starting $prog_used using PICARD tools.. ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_SS'}){
		$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";
	if ( defined $reference ){
		$picard_params .= " REFERENCE=$reference ";
	}

	#if ( defined $cfg_hash->{'picard_CREATE_INDEX'} ){
		#$picard_params .= " CREATE_INDEX=".$cfg_hash->{'picard_CREATE_INDEX'};
	#}
		
	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile CREATE_INDEX=TRUE";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_PICARD_SortSam 

 Title   : run_PICARD_SortSam
 Usage   : run_PICARD_SortSam(   );

 Function: run_PICARD_SortSam
			Is a tool from PICARD which permits to sort a bam or a sam file using
			a specific sort order that is required!
			Possible values: {unsorted, queryname, coordinate, duplicate, unknown}
			
 Returns : nothing
=cut
sub run_PICARD_SortSam{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	#print_and_log( "Starting $prog_used using PICARD tools.. ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_SS'}){
		$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";
	if ( defined $cfg_hash->{'picard_SORT_ORDER'} ){
		$picard_params .= " SORT_ORDER=".$cfg_hash->{'picard_SORT_ORDER'};
	}

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file OUTPUT=$outFile";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_PICARD_SortVcf 

 Title   : run_PICARD_SortVcf
 Usage   : run_PICARD_SortVcf(   );

 Function: run_PICARD_SortVcf
			Is a tool from PICARD which permits to sort a vcf file using
			the refernce dictionary file
			
 Returns : nothing
=cut
sub run_PICARD_SortVcf{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	#print_and_log( "Starting $prog_used using PICARD tools.. ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_SS'}){
		$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $dictionary = extract_name($cfg_hash->{'hum_ref'},"noext").".".$cfg_hash->{'dict_ext'};

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder SEQUENCE_DICTIONARY=$dictionary INPUT=$input_file OUTPUT=$outFile";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}



=head2 run_PICARD_GatherVcfs 

 Title   : run_PICARD_GatherVcfs
 Usage   : run_PICARD_GatherVcfs(   );

 Function: run_PICARD_GatherVcfs
			Merge VCFs input file in one VFC output file 
			
 Returns : returns the command if needed to run in a job outside
=cut
sub run_PICARD_GatherVcfs{

	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_vcfs = shift;
	my $outFile = shift;
	my $log_file = shift;
	my $get_command = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};		

	print_and_log( "Starting $prog_used using PICARD tools..\n ",$log_file);#DEBUGCODE
						  								
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_CV'}){
		$java_call .= " -".$cfg_hash->{'javamem_CV'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	#build input VCFs string
	my @input_vcfs = split(",",$input_vcfs);
	my $input_str = "";
	foreach my $input_vcf (@input_vcfs){
			$input_str .= " I=$input_vcf ";
	}
	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $input_str OUTPUT=$outFile";	
	if ($get_command ne 'GET_CMD'){	 
	 print_and_log( "Executing command: $command\n",$log_file);
	 try_exec_command($command) or die "Unable to execute command: $command\n";
	}
	return $command
}

=head2 run_PICARD_BuildBamIndex 

 Title   : run_PICARD_BuildBamIndex
 Usage   : run_PICARD_BuildBamIndex(   );

 Function: run_PICARD_BuildBamIndex
			Is a tool from PICARD which permits to index a bam  file.
			This tool creates an index file for the input BAM that allows fast
			look-up of data in a BAM file, lke an index on a database. 
			Note that this tool cannot be run on SAM files, 
			and that the input BAM file must be sorted in coordinate order.
			
			If you do not want to give the ouput name use NONE
 Returns : prints the BAI file as $input_file.bam.bai if $out_file is not given
=cut
sub run_PICARD_BuildBamIndex{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	print_and_log( "Starting $prog_used using PICARD tools..\n ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_SS'}){
		$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params INPUT=$input_file";
	#Add output if is given
	if ( $outFile ne 'NONE'){
		$command .= " OUTPUT=$outFile ";
	}
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_PICARD_CreateSequenceDictionary

 Title   : run_PICARD_CreateSequenceDictionary
 Usage   : run_PICARD_CreateSequenceDictionary(   );

 Function: run_PICARD_CreateSequenceDictionary
			Is a tool from PICARD which permits to generate a dictionary.
			
			If you do not want to give the ouput name use NONE
 Returns : prints the BAI file as $input_file.bam.bai if $out_file is not given
=cut
sub run_PICARD_CreateSequenceDictionary{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	print_and_log( "Starting $prog_used using PICARD tools..\n ",$log_file);
																		
	#Program execution
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_SS'}){
		$java_call .= " -".$cfg_hash->{'javamem_SS'}." ";
	}
	$java_call .= " -jar ".$cfg_hash->{'picard_path'}." $prog_used ";
	
	my $picard_params = "";

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $java_call." TMP_DIR=$tempFolder $picard_params REFERENCE=$input_file";
	#Add output if is given
	if ( $outFile ne 'NONE'){
		$command .= " OUTPUT=$outFile ";
	}
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}


1;
