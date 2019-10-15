
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
    
package LIB::bwakit_management;
## bwakit_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of the bwa-kit tools for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_BWAKIT_altmapping run_BWAKIT_samblaster run_BWAKIT_bwa run_BWAKIT_seqtk_mergepe 
									run_BWAKIT_seqtk_seq run_STAR);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 
use File::Copy;#To manage files

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit  execute_threads build_input_name_from_executed);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present extract_name);

use LIB::files_manipulation qw();

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);



##############


=head2 run_BWAKIT_seqtk

 Title   : run_BWAKIT_seqtk
 Usage   : run_BWAKIT_seqtk(   );

 Function: Runs bwa_kit seqtk
 Returns : nothing
=cut
sub run_BWAKIT_seqtk_seq{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	my $prev_task = "qc";	
	#Get infos for subject sample
	#print_and_log("Getting information for $readf_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash->{'db_readf_id'},$readf_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	
	my $read1 = "";
	my $read2 = "";
	my $outTemp1 = "";
	my $outTemp2 = "";
		
	#print_and_log( "group: $group_id and task: $task..\n",$log_file);
	my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
	my $inFolder = $cfg_hash->{$group_id.'_'.$prev_task.'_out_f'};
			
	#Check from the DB if the Trimming was performed
	if ( $params->{$cfg_hash->{'trim_step'}} == 1 ){
			
		##Build the names of the two trimmed files
		#The names that is giving VarGenius	
		my $trimmed1 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
		my $trimmed2 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
	
		#Check if they exist
		if (file_not_present($trimmed1)  > 0){ die "Cannot proceed with $prog_used! Check: $trimmed1.\n";}
		if (file_not_present($trimmed2)  > 0){ die "Cannot proceed with $prog_used! Check: $trimmed2.\n";}
		#Assign the names
		$read1 = $trimmed1;
		$read2 = $trimmed2;
		$outTemp1 = $trimmed1.".temp";
		$outTemp2 = $trimmed2.".temp";
	}else {
		print_and_log( "(The trimming for $readf_id has not been performed. I will use the raw.)",$log_file);
		my $fq1 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1".$cfg_hash->{'trim_out_ext'};
		my $fq2 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2".$cfg_hash->{'trim_out_ext'};
		#Assign the names
		$read1 = $fq1;
		$read2 = $fq2;
		$outTemp1 = $fq1.".temp";
		$outTemp2 = $fq2.".temp";
	}
	#Now that you have the names, convert the scores
	if ( defined $cfg_hash->{'align_type'} and $cfg_hash->{'align_type'} eq 'PE'){
		my $command = $cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'seqtk_prog'}." seq -Q64 -V $read1 | gzip > $outTemp1";
		print_and_log( "Executing command: $command - ",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
		$command = $cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'seqtk_prog'}." ".$cfg_hash->{'seqtk_seq_prog'}." -Q64 -V $read2 | gzip > $outTemp2";
		print_and_log( "Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
		move($outTemp1,$read1);
		move($outTemp2,$read2);
	}

}

=head2 run_BWAKIT_seqtk_mergepe

 Title   : run_BWAKIT_seqtk_mergepe
 Usage   : run_BWAKIT_seqtk_mergepe(   );

 Function: Runs bwa_kit mergepe that merges the paired end in a 
						unique file that will be used in BWA
						
 Returns : nothing
=cut
sub run_BWAKIT_seqtk_mergepe{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	my $prev_task = "qc";
	
	#Get infos for subject sample
	#print_and_log("Getting information for $readf_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
	my $inFolder = $cfg_hash->{$group_id.'_'.$prev_task.'_out_f'};
		
	my $read1 = "";
	my $read2 = "";
	
		#Check from the DB if the Trimming was performed
	if ( $params->{$cfg_hash->{'trim_step'}} == 1 ){
		#The names that is giving VarGenius	
		$read1 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
		$read2 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
	}else {
		print_and_log( "The trimming for $readf_id has not been performed. I will use the raw.\n",$log_file);
		$read1 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1".$cfg_hash->{'trim_out_ext'};
		$read2 = $inFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2".$cfg_hash->{'trim_out_ext'};	
	}			
	#Check if they exist
	if (file_not_present($read1)  > 0){ die "Cannot proceed with $prog_used! Check: $read1.\n";}
	if (file_not_present($read2)  > 0){ die "Cannot proceed with $prog_used! Check: $read2.\n";}
	
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergepe_step'},
		$params->{$main_name},$steps_array)."_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'};
		
	if ( defined $cfg_hash->{'align_type'} and $cfg_hash->{'align_type'} eq 'PE'){
			my $command = $cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'seqtk_prog'}." ".$cfg_hash->{'seqtk_mergepe_prog'}." $read1 $read2 > $outFile";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
	}

}


=head2 run_STAR 

 Title   : run_STAR
 Usage   : run_STAR(   );

 Function: Runs STAR 
 Returns : nothing
=cut
sub run_STAR{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	my $param_in = shift;
	my $runMode = shift;
		
	#Get infos for subject sample
	#print_and_log("Getting information for $readf_id from ".$cfg_hash->{'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
				$cfg_hash->{'db_readf_id'},$readf_id);
	#Getting the group name
	my $group_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	#Getting the sample id...
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_readf_id'},$readf_id);
	#..to get the sample name						    	
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $prev_task = "qc";
	
	#Initialize the parameter string for the input files
	my $param_str = "";
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
					    
	my $outFile = "";
	
	#If we are not producing indexing file for STAR, execute this code
	#to get the actual file names (trimmed or not)				    
	if ( $runMode ne 'genomeGenerate'){
		if ( defined $cfg_hash->{'align_type'} and $cfg_hash->{'align_type'} eq 'PE'){
			my $trimmed1 = "";
			my $trimmed2 = "";
			
			#Sample sheet variables (These I will take from the DB)
			my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
			my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
						
			#Build output name 
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $qcFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
			
						
			#Gives error if merge paired ends is erroneously used
			if ( $params->{$cfg_hash->{'mergepe_step'}} == 1){
				log_and_exit( "ERROR: ".$cfg_hash->{'mergepe_step'}." step cannot be used in STAR \n",$log_file);	
			}else{
				#print_and_log( "The program to make interleaved reads has not been run...\n",$log_file);#DEBUGCODE
				#Check from the DB if the Trimming was performed
				if ( $params->{$cfg_hash->{'trim_step'}} == 1 ){
					#Build the names of the two trimmed files
					$trimmed1 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
					$trimmed2 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
				}else{
					print_and_log( "No trimming, aligning the raw fastq...\n",$log_file);	
					$trimmed1 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1".$cfg_hash->{'trim_out_ext'};
					$trimmed2 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2".$cfg_hash->{'trim_out_ext'};
				}
				#Set the name of the output file
				$outFile = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'align_step'},
						$params->{$main_name},$steps_array)."_".$cfg_hash->{'align_step'};
			}
			
			#Check if they exist
			if ($trimmed1  ne '' and file_not_present($trimmed1) > 0 ){ die "Cannot proceed with $prog_used! Check: $trimmed1.\n";}
			if ($trimmed2  ne '' and file_not_present($trimmed2) > 0 ){ die "Cannot proceed with $prog_used! Check: $trimmed2.\n";}

			#Give read files input		
			$param_str .= " --readFilesIn $trimmed1 $trimmed2";		
			
			#if read files are in .gz format use the readFilesCommand parameter
			if ($trimmed1 =~ /\.gz$/ and $trimmed2 =~ /\.gz$/){
					$param_str .= " --readFilesCommand zcat ";
			}
			
			#Set the outname prefix
			$param_str .= " --outFileNamePrefix $outFile ";
		}			
	}else{
		#Set the run mode
		$param_str .= " --runMode genomeGenerate ";
			#Add reference genome
			if (defined $cfg_hash->{'hum_ref'}){
					$param_str .= " --genomeFastaFiles ".$cfg_hash->{'hum_ref'}." ";
			}		
	}



	#Threads
	if (defined $cfg_hash->{'align_threads'}){
		$param_str .= " --runThreadN ".$cfg_hash->{'align_threads'}." ";
	}

	my $command = $cfg_hash->{'star_path'}." $param_str $param_in ";
	print_and_log( "Executing command: $command ",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	if ( $runMode ne 'genomeGenerate'){
		#Change the name of the output file 
		print_and_log( "Moving $outFile".$cfg_hash->{'STAR_alignout'}."  to $outFile.".$cfg_hash->{'sam_ext'}." \n",$log_file);#DEBUGCODE		
		move($outFile.$cfg_hash->{'STAR_alignout'}, $outFile.".".$cfg_hash->{'sam_ext'}) or die "ERROR: Cannot move $outFile".$cfg_hash->{'STAR_alignout'}."  to $outFile.".$cfg_hash->{'sam_ext'}."\n";
	}
	
	return $outFile;
}



=head2 run_BWAKIT_bwa 

 Title   : run_BWAKIT_bwa
 Usage   : run_BWAKIT_bwa(   );

 Function: Runs BWA 
 
 @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
where the ’\t’ stands for the tab character
The @RG identifies that this is the “read group” tag. Following the tabs we have the following keys:
ID – globally unique string identifying this run. In general we use an 
	identifier linked to the lane where the data was run (for illumina data) 
	or instrument where lane is not applicable.
SM – the name associated with the DNA sample in this file. In general this 
	will be a sample identifier such as NA12878, but it can be any string. 
	This is the most important field for the GATK because all analysis is 
	done by sample, and this is what the tools will use to define what data
	 belongs to which sample.
PL – the platform used. The list of supported platforms can be found in
 the GATK’s online manual as it is ever growing. For example you can use 
 “illumina”, “pacbio” or “iontorrent” here.
LB – an identifier of the library from which this DNA was sequenced. This 
	is important for future reference and quality control. In case errors
	 are associated with the DNA prep, this will be the field linking the 
	 data to the laboratory process.
PU – the platform unit identifier for this run. Any generic identifier 
	that will allow any investigation to go back to the very machine, and 
	time where this data was run. We typically use the “flowcell-barcode.lane” 
	unique identifier for Illumina runs.
The read group information is key for downstream GATK functionality. 
The GATK will not work without a read group tag. Make sure to enter as 
much metadata as you know about your data in the read group fields provided.
 For more information about all the possible fields in the @RG tag, take a
  look at the SAM specification (http://samtools.sourceforge.net/SAM1.pdf).

 Returns : nothing
=cut
sub run_BWAKIT_bwa{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	#Get infos for subject sample
	#print_and_log("Getting information for $readf_id from ".$cfg_hash->{'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
				$cfg_hash->{'db_readf_id'},$readf_id);
	#Getting the group name
	my $group_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	#Getting the sample id...
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_readf_id'},$readf_id);
	#..to get the sample name						    	
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $prev_task = "qc";
	
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
						    
						    
		if ( defined $cfg_hash->{'align_type'} and $cfg_hash->{'align_type'} eq 'PE'){
			my $trimmed1 = "";
			my $trimmed2 = "";
			
			#Sample sheet variables (These I will take from the DB)
			my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
			my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
						
			#Build output name 
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $qcFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
			my $outFile = "";
						
			#If SEQKT was run then there is a unique file for Paired ends in the alignment folder
			if ( $params->{$cfg_hash->{'mergepe_step'}} == 1){
				#$trimmed1 = $outFolder."/".$params->{$main_name}."_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'};
				$trimmed1 = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'align_step'},
						$params->{$main_name},$steps_array).".".$cfg_hash->{'fastq_ext'};
			
				$outFile = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'align_step'},
						$params->{$main_name},$steps_array)."_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'sam_ext'};
						
				print_and_log( "(Mergepe step executed using only $trimmed1 for the alignment ",$log_file);
			}else{
				#print_and_log( "The program to make interleaved reads has not been run...\n",$log_file);#DEBUGCODE
				#Check from the DB if the Trimming was performed
				if ( $params->{$cfg_hash->{'trim_step'}} == 1 ){
					#Build the names of the two trimmed files
					$trimmed1 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
					$trimmed2 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
					print_and_log( "Trimming was executed using  $trimmed1 and $trimmed2 for the alignment.. ",$log_file);						
				}else{
					
					$trimmed1 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1".$cfg_hash->{'trim_out_ext'};
					$trimmed2 = $qcFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2".$cfg_hash->{'trim_out_ext'};
					print_and_log( "(No trimming, aligning the raw fastq) $trimmed1 and $trimmed2 ",$log_file);	
				}

				$outFile = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'align_step'},
						$params->{$main_name},$steps_array)."_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'sam_ext'};
			}
			
			#Check if they exist
			if ($trimmed1  ne '' and file_not_present($trimmed1) > 0 ){ die "Cannot proceed with $prog_used! Check: $trimmed1.\n";}
			if ($trimmed2  ne '' and file_not_present($trimmed2) > 0 ){ die "Cannot proceed with $prog_used! Check: $trimmed2.\n";}
						
			my $param_str = " ";
			#Put Information about group id if exist
			if (defined $cfg_hash->{'vcf_rg_id'} or defined $cfg_hash->{'vcf_rg_lb'} or defined $cfg_hash->{'vcf_rg_sm'}
						or defined $cfg_hash->{'vcf_rg_pu'} or defined $cfg_hash->{'vcf_rg_pl'}){
						$param_str .= ' -R"@RG\t';		
			
				if (defined $cfg_hash->{'vcf_rg_id'}){
					$param_str .= $cfg_hash->{'vcf_rg_id'}.":".$flowcell.".".$lane.'\t';
				}
				if (defined $cfg_hash->{'vcf_rg_lb'}){
					$param_str .= $cfg_hash->{'vcf_rg_lb'}.":".$targetbed.'\t';
				}
				if (defined $cfg_hash->{'vcf_rg_sm'}){
					$param_str .= $cfg_hash->{'vcf_rg_sm'}.":".$sample_name.'\t';
				}				
				if (defined $cfg_hash->{'vcf_rg_pu'}){
					$param_str .= $cfg_hash->{'vcf_rg_pu'}.":".$analysis_id.'\t';
				}			
				if (defined $cfg_hash->{'vcf_rg_pl'}){
					$param_str .= $cfg_hash->{'vcf_rg_pl'}.":".$cfg_hash->{'vcf_rg_pl_val'}.'\t';
				}
				chop($param_str);#Remove \
				chop($param_str);#Remove t
				$param_str .= '"'; 
			}
			#Add reference genome
			if (defined $cfg_hash->{'hum_ref'}){
					$param_str .= " ".$cfg_hash->{'hum_ref'}." ";
			}			
			
			#Paired end reads interleaved or not
			if ( $params->{$cfg_hash->{'mergepe_step'}} == 1 ){
				$param_str .= " -p ";
			}
			if ( $trimmed1 ne ''){
					$param_str .= " $trimmed1 ";
			}	
			if ( $trimmed2 ne ''){
				if ( $params->{$cfg_hash->{'mergepe_step'}} == 0 ){
						$param_str .= " $trimmed2 ";
				}#else{
						#print_and_log( "WARNING: Something went wrong with parameters. You are using both the reads for BWA in BWA-KIT".
							#" as well as keeping seqtk'} == 1 in the configuration file...\n",$log_file);
				#}
			}
			#Threads
			if (defined $cfg_hash->{'align_threads'}){
				$param_str .= " -t ".$cfg_hash->{'align_threads'}." ";
			}
	
			my $command = $cfg_hash->{'bwa_path'}." mem $param_str > $outFile";
			print_and_log( "Executing command: $command \n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}
}


=head2 run_BWAKIT_samblaster

 Title   : run_BWAKIT_samblaster
 Usage   : run_BWAKIT_samblaster(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bwa_kit samblaster
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BWAKIT_samblaster{
	my $cfg_hash = shift;
	my $prog_used = shift;
	#my $readf_id = shift;
	my $group_id = shift;
	my $input_file = shift;
	my $outFile = shift;
	#my $task = shift;		
	#my $steps_array = shift;
	#my $step = shift;	
	my $log_file = shift;

	
	##Get infos for subject sample from the database
	#print_and_log("Getting information for $readf_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);
	#my $params;
	#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash->{'db_dsn'},
						#$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
	
	#my $main_name = $cfg_hash->{'db_readf_name'};
	#Check from the DB if the Alignment was performed
	if ( -e  $input_file){
		
			#my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
			##Set the input file name
			#my $input_file = "";#Input file
			#if ( $params->{$cfg_hash->{'mergepe_step'}} == 1){
				#$input_file = $outFolder."/".$params->{$main_name}."_".$cfg_hash->{'mergepe_step'}.
										#"_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'sam_ext'};
			#}else{
				#$input_file = $outFolder."/".$params->{$main_name}."_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'sam_ext'};
			#}
			
			#$input_file = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
							#$params->{$main_name},$steps_array).".".$cfg_hash->{'sam_ext'};
			
			#if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			##Prepare output name
			#my $outFile = $outFolder."/".extract_name($input_file,1)."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'sam_ext'};
			 
			my $param_str = " ";	
			
			#If active will remove duplicates
			if (defined $cfg_hash->{'sambl_rem_dups'}){
				$param_str .=  " --".$cfg_hash->{'sambl_rem_dups'}." ";
			}
			
			#Execute the command
			my $command = $cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'samblaster_prog'}." --input $input_file  $param_str > $outFile";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
	}else {
		print_and_log( "ERROR: Cannot find $input_file. Cannot execute $prog_used..\n",$log_file);
	}
}


=head2 run_BWAKIT_altmapping

 Title   : run_BWAKIT_altmapping
 Usage   : run_BWAKIT_altmapping(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs both k8 and bwa-postalt.js
					
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BWAKIT_altmapping{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $readf_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;		
	my $steps_array = shift;
	my $step = shift;
	
	
	#Get infos for subject sample from the database
	#print_and_log("Getting information for $readf_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);#
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	#Check from the DB if the Alignment was performed
	if ( $params->{$cfg_hash->{'align_step'}} == 1 ){
			print "Settings for gr: $group_id and task $task\n";
			my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
			#Set the input file name
			my $input_file = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'map_alt_step'},
							$params->{$main_name},$steps_array).".".$cfg_hash->{'sam_ext'};
			
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$cfg_hash->{'map_alt_step'}.".".$cfg_hash->{'postalt_prefix'};
			my $param_str = " ";	
			$param_str .=  " -p $outFile ";

			
			#Execute the command
			my $command = $cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'k8_prog'}." ".
										$cfg_hash->{'bwakit_path'}."/".$cfg_hash->{'map_alt_prog'}." $param_str $input_file ".$cfg_hash->{'hum_ref'}." > $outFile";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
						
	}else {
		print_and_log( "ERROR: The Alignment for $readf_id has not been performed. Cannot execute $prog_used..\n",$log_file);
	}
}


1;
