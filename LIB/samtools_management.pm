
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
    
package LIB::samtools_management;
## samtools_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of the samtools for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(  run_SAMTOOLS_Index run_SAMTOOLS_Sort run_SAMTOOLS_View
									run_SAMTOOLS_View_gen run_SAMTOOLS_Sort_gen run_SAMTOOLS_Index_gen
									run_SAMTOOLS_flagstat run_SAMTOOLS_mpileup);
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
use LIB::files_management qw( save_hash load_hash file_not_present  extract_name);

use LIB::files_manipulation qw();

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);


##################



=head2 run_SAMTOOLS_View

 Title   : run_SAMTOOLS_View
 Usage   : run_SAMTOOLS_View(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs samtools sort -b -S file.sam > file.bam
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_View{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	#my $final_step = $cfg_hash->{'sort_idx_step'};
	my $step_needed = $cfg_hash->{'align_step'};
	
	#Get infos for subject sample from the database
	#print_and_log("Getting information for $sample_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
				$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
				$cfg_hash-> {'db_readf_id'},$sample_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	
	#Check from the DB if the Alignment was performed
	if ( $params->{$step_needed} == 1 ){
		
			my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
			#Set the input file name
			my $input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},$steps_array).".".$cfg_hash->{'sam_ext'};
			
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1).".".$cfg_hash->{'bam_ext'};
			my $param_str = " -b -S ";	
			
			#Execute the command
			my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'samview_prog'}." $param_str $input_file > $outFile";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
	}else {
		print_and_log( "ERROR: The Alignment for $sample_id has not been performed. Cannot execute $prog_used..\n",$log_file);
	}
}


=head2 run_SAMTOOLS_Sort 

 Title   : run_SAMTOOLS_Sort
 Usage   : run_SAMTOOLS_Sort(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs samtools sort
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_Sort{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	#Here I call this 'final' step because when I construct the name the .bam 
	#output file will be the last one. After that, index will produce a .bai associated to 
	#ŧhis output
	#my $final_step = $cfg_hash->{'sort_idx_step'};
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $step_needed = $cfg_hash->{'align_step'};#samview_step

	#Get parameters for subject sample from the database
	print_and_log("Getting information for $sample_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
				$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
				$cfg_hash-> {'db_readf_id'},$sample_id);
	
	#Check from the DB if the Alignment was performed
	if ( $params->{$step_needed} == 1 ){
		
			my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
			#Set the input file name
			my $input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$step.".".$cfg_hash->{'bam_ext'};
			my $param_str = " ";	
			
			if (defined $cfg_hash->{'align_cpus'}){
				$param_str .= " -@ ".$cfg_hash->{'align_cpus'}." ";
			}
			if (defined $cfg_hash->{'align_threads_mem'}){
				$param_str .= " -m ".$cfg_hash->{'align_threads_mem'}." ";
			}

			#Execute the command
			my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'sort_prog'}." $param_str $input_file > $outFile";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
	}else {
		print_and_log( "ERROR: The Alignment for $sample_id has not been performed. Cannot execute $prog_used..\n",$log_file);
	}
}


=head2 run_SAMTOOLS_Index

 Title   : run_SAMTOOLS_Index
 Usage   : run_SAMTOOLS_Index(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs samtools index
					
	INPUT: needs bam output from the alignment against the reference
		Here the input file is the same as from the sorting since we consider
		Sort and indexing as the same operation				
 Returns : nothing
=cut
sub run_SAMTOOLS_Index{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	#my $final_step = $cfg_hash->{'sort_idx_step'};
	my $step_needed = $cfg_hash->{'align_step'};#samsort_step
	
	#Get infos for subject sample from the database
	print_and_log("Getting information for $sample_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$sample_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	#Check from the DB if the Alignment was performed
	if ( $params->{$step_needed} == 1 ){
		
			my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
			#Set the input file name. Here the input file is the same as from the sorting since we consider
			#Sort and indexing as the same operation
			my $input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},$steps_array)."_".$step.".".$cfg_hash->{'bam_ext'};
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#The output will be the indexed file .bai
	
			#Execute the command
			my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'index_prog'}." $input_file ";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
	}else {
		print_and_log( "ERROR: The Alignment for $sample_id has not been performed. Cannot execute $prog_used..\n",$log_file);
	}
}

=head2 run_SAMTOOLS_mpileup

 Title   : run_SAMTOOLS_flagstat
 Usage   : run_SAMTOOLS_flagstat(   );

 Function: Executes SAMTOOLS-mpileup for variant calling
								
					INPUT: needs a list of bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_mpileup{
	my $cfg_hash = shift;
	my $bam_inputs = shift;
	my $outFile = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'mpileup_prog'}." -go $outFile -f ".$cfg_hash->{'hum_ref'}." $bam_inputs  > $outFile";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

=head2 run_SAMTOOLS_flagstat

 Title   : run_SAMTOOLS_flagstat
 Usage   : run_SAMTOOLS_flagstat(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs samtools flagstat
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_flagstat{
	my $cfg_hash = shift;
	my $bam_input = shift;
	my $outFile = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'flagstat_prog'}."  $bam_input  > $outFile";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}



=head2 run_SAMTOOLS_View_gen

 Title   : run_SAMTOOLS_View_gen
 Usage   : run_SAMTOOLS_View_gen(   );

 Function: Converts SAM input in a BAM output
					
					INPUT: needs SAM output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_View_gen{
	my $cfg_hash = shift;
	my $sam_input = shift;
	my $bam_output = shift;
	my $log_file = shift;

	my $param_str = " -b -S ";	
	
	#Execute the command
	my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'samview_prog'}." $param_str $sam_input > $bam_output";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
							
}

=head2 run_SAMTOOLS_Sort_gen

 Title   : run_SAMTOOLS_Sort_gen
 Usage   : run_SAMTOOLS_Sort_gen(   );

 Function: Sorts a bam file in input
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_Sort_gen{
	my $cfg_hash = shift;
	my $bam_input = shift;
	my $bam_output = shift;
	my $log_file = shift;

	my $param_str = " ";	
	
	if (defined $cfg_hash->{'samsort_threads'}){
		$param_str .= " -@ ".$cfg_hash->{'samsort_threads'}." ";
	}
	if (defined $cfg_hash->{'samsort_mem'}){
		$param_str .= " -m ".$cfg_hash->{'samsort_mem'}." ";
	}

	#Execute the command
	my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'sort_prog'}." $param_str $bam_input > $bam_output";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
							
}

=head2 run_SAMTOOLS_Index_gen

 Title   : run_SAMTOOLS_Index_gen
 Usage   : run_SAMTOOLS_Index_gen(   );

 Function: Indexes the bam output from sort
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_SAMTOOLS_Index_gen{
	my $cfg_hash = shift;
	my $bam_input = shift;
	my $log_file = shift;

	#Execute the command
	my $command = $cfg_hash->{'samtools_path'}." ".$cfg_hash->{'index_prog'}." $bam_input ";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
							
}

1;
