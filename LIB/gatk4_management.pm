
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
    
package LIB::gatk4_management;
## gatk_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of the gatk tools for the pipeline

BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_GATK4_RealignerTargetCreator run_GATK4_IndelRealigner 
            run_GATK4_BaseRecalibrator_perlane run_GATK4_HaplotypeCaller_parallel
            run_GATK4_CatVariants run_GATK4_GenotypeGVCFs run_GATK4_VariantFiltration
            run_GATK4_VariantRecalibration run_GATK4_ApplyRecalibration 
            run_GATK4_PhaseByTransmission run_GATK4_DepthOfCoverage
            run_GATK4_BaseRecalibrator_persample run_GATK4_PrintReads_persample
            run_GATK4_SelectVariants run_GATK4_CombineVariants
            run_GATK4_VariantAnnotator run_GATK4_CalculateGenotypePosteriors 
            run_GATK4_SplitNCigarReads run_GATK4_GenotypeGVCFs_parallel
            run_GATK4_ApplyBQSR run_GATK_MergeVcfs run_GATK4_GenomicsDBImport_parallel
            run_GATK4_ApplyVQSR);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit  execute_threads build_input_name_from_executed
				execute_jobs);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present delete_file 
								extract_name split_bedfile);
																
use LIB::files_manipulation qw();

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);

##########################################################GATK TOOLS COMING FROM PICARD


=head2 run_GATK_MergeVcfs 

 Title   : run_GATK_MergeVcfs
 Usage   : run_GATK_MergeVcfs(   );

 Function: run_GATK_MergeVcfs
			This is a tool from PICARD and integrated into GATK which permits to merge VCF 
			files. This function takes in input 
			input_files:  a string with VCFs paths separated by comma
			outFile: the VCF output path
			
 Returns : nothing
=cut
sub run_GATK_MergeVcfs{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_files = shift;
	my $outFile = shift;
	my $log_file = shift;

	my $tempFolder = $cfg_hash->{'scratch_f'};	

	#print_and_log( "Starting $prog_used using GATK-PICARD tools.. ",$log_file);
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_SS'} or defined $cfg_hash->{'scratch_f'} ){
		#Open string
		$jatk_call .= " --java-options '";
		if (defined $cfg_hash->{'javamem_SS'}){
			$jatk_call .= " -".$cfg_hash->{'javamem_SS'}." ";
		}	
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}	
							
		#Close string
		$jatk_call .= "'";
	}
	#######END JAVA PARAMETERS
	
	my @input_files = split(",",$input_files);
	my $input_str = "";
	foreach my $input_file (@input_files){
			$input_str .= " -I $input_file ";
	}
	
	my $picard_params = "";
	#if ( defined $cfg_hash->{'picard_SORT_ORDER'} ){
		#$picard_params .= " SORT_ORDER=".$cfg_hash->{'picard_SORT_ORDER'};
	#}

	#print_and_log( "$java_call - picard params: $picard_params\n",$log_file);
	my $command = $jatk_call."  $picard_params $input_str --OUTPUT $outFile";		 
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

##############################################################


=head2 run_GATK4_SplitNCigarReads

 Title   : run_GATK4_SplitNCigarReads
 Usage   : run_GATK4_SplitNCigarReads(   );

 Function: run_GATK4_SplitNCigarReads
			Is a tool from GATK which splits reads into exon segments 
			(getting rid of Ns but maintaining grouping information) and 
			hard-clip any sequences overhanging into the intronic regions.
 Returns : nothing
=cut
sub run_GATK4_SplitNCigarReads{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_bam =shift;
	my $out_bam = shift;
	my $log_file = shift;
	
	#my $previous_task = 'align';
	#my $step_needed = $cfg_hash->{'sort_idx_step'};
	
	#Check from the DB if the Sorting and indexing was performed
	#if ( $params->{$step_needed} == 1 ){
																	
		print_and_log( "Starting SplitNCigarReads for $input_bam using GATK tools..\n ",$log_file);
		
		my $param_str = " ";	
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_SNCR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_SNCR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_SNCR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS


		if ( defined $cfg_hash->{'SNCR_read_filter'} ){
			$param_str .= " --read_filter ".$cfg_hash->{'SNCR_read_filter'}." ";
		}
		if ( defined $cfg_hash->{'SNCR_reassign_mapping_quality_from'} and correct_type($cfg_hash->{'SNCR_reassign_mapping_quality_from'},"positiveint") ){
					$param_str .= " --reassign_mapping_quality_from ".$cfg_hash->{'SNCR_reassign_mapping_quality_from'};
		}
			
		if ( defined $cfg_hash->{'SNCR_reassign_mapping_quality_to'} and correct_type($cfg_hash->{'SNCR_reassign_mapping_quality_to'},"positiveint") ){
					$param_str .= " --reassign_mapping_quality_to ".$cfg_hash->{'SNCR_reassign_mapping_quality_to'};
		}
		if ( defined $cfg_hash->{'SNCR_unsafe'} ){
			$param_str .= " --unsafe ".$cfg_hash->{'SNCR_unsafe'}." ";
		}		
		

		
		if (file_not_present($input_bam) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_bam.\n";}
		my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bam ";			 
		print_and_log("Executing command: $command \n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";

	#}else{print_and_log( "ERROR: Step $step_needed has not been run for sample: $sample_id.f\n ",$log_file);}
	
}


=head2 run_GATK4_RealignerTargetCreator

 Title   : run_GATK4_RealignerTargetCreator
 Usage   : run_GATK4_RealignerTargetCreator(   );

 Function: run_GATK4_RealignerTargetCreator
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			This function takes in input the bam file from the alignment and produces a file 
			to be used for the indel realignment later. Hence in this step 
 Returns : nothing
=cut
sub run_GATK4_RealignerTargetCreator{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;
  my $steps_array = shift;      
	my $step = shift;
	
	my $previous_task = 'align';
	my $step_needed = $cfg_hash->{'sort_idx_step'};
	
	#Get infos for subject sample
	print_and_log("Getting information for $sample_id from ".$cfg_hash->{'db_name'}."...\n",$log_file);
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$sample_id);
	
	my $main_name = $cfg_hash->{'db_readf_name'};
	#Check from the DB if the Sorting and indexing was performed
	if ( $params->{$step_needed} == 1 ){
		#InFolder is the alignment one at this stage
		my $inFolder = $cfg_hash->{$group_id.'_'.$previous_task.'_out_f'};
		my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$group_id);
						    			
		#Set the input file name
		my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};
				
		print_and_log( "Starting RealignerTargetCreator for $input_bam using GATK tools..\n ",$log_file);
		
		#Prepare output name
		my $targ_list = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array)."_".$step.".".$cfg_hash->{'targ_int_ext'};
		my $param_str = " ";	
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_RTC'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_RTC'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_RTC'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

		#if ( defined $cfg_hash->{'nt_RTC'} ){
			#$param_str .= " -nt ".$cfg_hash->{'nt_RTC'};
		#}
		if ( $cfg_hash->{'filter_mismatching_base_and_quals_RTC'} eq 'YES'){
			$param_str .= " --filter_mismatching_base_and_quals ";
		}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}  
                
		#Adds known indel snp to the list of known
		if ( defined $cfg_hash->{'known_RTC'} ){
			my @known_indel_snp  = split($cfg_hash->{'word_sep'},$cfg_hash->{'known_RTC'});
			foreach my $known (@known_indel_snp){
				$param_str .= " -known ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
			}
		}

			
		my $numExec = 1;
		if ( $perchrom == 1 ){
			##Using also Chromosomex X and Y that Romans did not use!!!
			$numExec = $cfg_hash->{'chromosomes_num'};
			my @chroms = 1..$numExec;
			push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
						
			#An array will contain all the commands to be executed
			my @commands = ();
			
			#Executing one GATK command for each target chromosome using bed files
			foreach my $exec (@chroms){
			#for (my $exec = 1; $exec <= $numExec; $exec++){
					my $l_string .= " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
					my $ith_targ = $outFolder."/".extract_name($targ_list,1)."_$exec.".$cfg_hash->{'targ_int_ext'};
					my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_targ ";			 
					print_and_log("Executing command: $command\n",$log_file);
					push(@commands,$command);
			}
			#Run a different thread for each command
			execute_threads(\@commands);
		}else{
			# When perchrom is not used the HaplotypeCaller uses the target regions as intervals	
			#Getting the target bed file using the folder for targets and the name contained in the database
			my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
										$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
										$cfg_hash->{'db_analysis_id'},$group_id);
			if (file_not_present($target_bed) > 0 ){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
			$param_str .= " --intervals $target_bed ";			
		
			if (file_not_present($input_bam) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_bam.\n";}
			my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $targ_list ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}
	}else{print_and_log( "ERROR: Step sort_idx_step has not been run for  sample: $sample_id.f\n ",$log_file);}
	
}


=head2 run_GATK4_IndelRealigner

 Title   : run_GATK4_IndelRealigner
 Usage   : run_GATK4_IndelRealigner(   );

 Function: run_GATK4_IndelRealigner
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			
 Returns : nothing
=cut
sub run_GATK4_IndelRealigner{
    my $cfg_hash = shift;
    my $prog_used = shift;
    my $sample_id = shift;
    my $group_id = shift;
    my $log_file = shift;
    my $task = shift;
    my $steps_array = shift;

    my $previous_task = 'align';#Used to know the folder where to take the input bam file
    my $step = $cfg_hash->{'indel_real_step'};
    my $prev_step = $cfg_hash->{'real_targ_step'};
    
    #Get infos for subject sample
    print_and_log("Getting information for $sample_id from ".$cfg_hash-> {'db_name'}."...\n",$log_file);
    my $params;
    getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
                                 $cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
                                 $cfg_hash-> {'db_readf_id'},$sample_id);
    
    my $main_name = $cfg_hash->{'db_readf_name'};
    #Check from the DB if the Sorting and indexing was performed
    if ( $params->{$prev_step} == 1 ){
			#inFolder for the bam file
			my $inFolder = $cfg_hash->{$group_id.'_'.$previous_task.'_out_f'};
			#Out folder
			my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};

			#Verifying that there has been an execution per chromosome
			my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
					    ,$cfg_hash->{'db_analysis_id'},$group_id);
					                    
			#Build the input bam file name using as stop point the realigner target 
			#creator step because at that time is constructed only a list of intervals
			my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
																			$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};

			#Add here the target intervals coming from RealignerTargetCreator
			my $targ_list = $outFolder."/".build_input_name_from_executed($params,$step,
																			$params->{$main_name},$steps_array).".".$cfg_hash->{'targ_int_ext'};

			#Prepare output name for bam files
			my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
																			$params->{$main_name},$steps_array)."_".$step.".".$cfg_hash->{'bam_ext'};
			
			#print_and_log( "Starting $prog_used for $input_bam using GATK tools..\n ",$log_file);#DEBUGCODE
			
			#Set the Java call
			my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_IR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_IR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_IR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

			#Parameters string construction
			my $param_str = " ";
			#Number of threads to use CANNOT USE IN INDELREALIGNER
			#if ( defined $cfg_hash->{'refine_threads'} ){
							#$param_str .= " -nt ".$cfg_hash->{'refine_threads'};
			#}
			if ( $cfg_hash->{'filter_bases_not_stored'} eq 'YES'){
							$param_str .= " --filter_bases_not_stored ";
			}
			if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
				$param_str .= " --fix_misencoded_quality_scores ";
			}
			if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
				$param_str .= " --allow_potentially_misencoded_quality_scores ";
			}

			#Adds known indel snp to the list of known
			if ( defined $cfg_hash->{'known_IR'} ){
				my @known_indel_snp  = split($cfg_hash->{'word_sep'},$cfg_hash->{'known_IR'});
				foreach my $known (@known_indel_snp){
					$param_str .= " -known ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
				}
			}

			my $numExec = 1;
			if ( $perchrom == 1){
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
					#Build i-th target file name
					my $ith_targ = $outFolder."/".extract_name($targ_list,1)."_".$exec.".".$cfg_hash->{'targ_int_ext'};
					#Filenotpresent > 1 means that the size could be zero!
					if ( file_not_present($ith_targ) > 1 ){ die "Cannot proceed with $prog_used! Check: $ith_targ.\n";}
					my $targ_interv = " --targetIntervals $ith_targ ";
					#Build ith chromosome target name
					my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
					#Build i-th bam output name
					my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string $targ_interv -o $ith_out_bam ";			 
					print_and_log("Executing a thread for: $command\n",$log_file);
					#try_exec_command($command) or die "Unable to execute command: $command\n";
					push(@commands,$command);
				}
				execute_threads(\@commands);
			}else{
					if ( file_not_present($input_bam) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_bam.\n";}
					if ( file_not_present($targ_list) > 0 ){ die "Cannot proceed with $prog_used! Check: $targ_list.\n";}
					#add target list to parameters
					$param_str .= " --targetIntervals $targ_list ";
					#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bam ";			 
					print_and_log("Executing command: $command\n",$log_file);
					try_exec_command($command) or die "Unable to execute command: $command\n";
			}
    }else{print_and_log( "ERROR: Step $step has not been run for sample: $sample_id.f\n ",$log_file);}
    
}


=head2 run_GATK4_ApplyBQSR

 Title   : run_GATK4_ApplyBQSR
 Usage   : run_GATK4_ApplyBQSR(   );

 Function: run_GATK4_ApplyBQSR
			Replaces the PrintReads for GATK4 to execute after BaseRecalibrato
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			avoid the next step of MergeSamFiles
 
 Returns : nothing
=cut
sub run_GATK4_ApplyBQSR{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $out_bam = shift;
	my $prev_task = shift;#Used to know the folder where to take the input bam file
	my $task = shift;
	my $steps_array = shift;      
	my $steps_array2 = shift;      
	my $log_file = shift;
	
	#my $previous_task = 'refine';
	my $step = $cfg_hash->{'base_recal_step'};
	my $prev_step = $cfg_hash->{'base_recal_step'};

	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
		
	my $db_sample_name = $cfg_hash->{'db_sample_name'};
	my $db_readf_name = $cfg_hash->{'db_readf_name'};
	
	#Check from the DB if the Sorting and indexing was performed
	if ( $params->{$prev_step} == 1 ){
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
				
		#The base quality score recal file generated previously
		my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2)."_$step.".$cfg_hash->{'recal_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		#print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";


		####### JAVA PARAMETERS	
		if (defined $cfg_hash->{'javamem_PR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
			$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_PR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_PR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
		#######END JAVA PARAMETERS
	

		#Parameters string construction
		my $param_str = " ";
		#Not used in GATK4
		#baq
		#if ( defined $cfg_hash->{'baq_PR'} ){
			#$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		#}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}		
		my $numExec = 1;
		
		#After the update of GATK that the IndelRealignment should not be done the
		#per-chromosome mode is deactivated here because starting from the alignment for
		#sure we do not have that separation. Hence check if indel_real was executed
		#Getting the indel_real_step status
		my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'}
														,$cfg_hash->{'db_sample_id'},$sample_id);	
		if ( $perchrom == 1  and $indel_real_step == 1 ){
			$numExec = $cfg_hash->{'chromosomes_num'};
			my @chroms = 1..$numExec;
			push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
						
			#An array will contain all the commands to be executed
			my @commands = ();
			
			#Executing one GATK command for each target chromosome using bed files
			foreach my $exec (@chroms){
				my $in_bam_set = "";
				
				##Construct a command with all reads files associated with this sample for each chromosome
				#I defined params here because I need it in the construction of output file name
				#there we need just one among all
				my $params;
				foreach my $readf_id (
					sort keys %{$distinct_readf_ids}
					)
					{
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used. The step before here is the step used because
						#the bam files did not change after base recalibration
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
						#Build i-th input name
						my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
														
						if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						$in_bam_set .= " -I $ith_in_bam ";									
					}
					#Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					#be better because all the bases are considered together
					my  $bqsr = " --bqsr-recal-file ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					#Build i-th bam output name
					my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					my $command = $jatk_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -O $ith_out_bam ";			 
					print_and_log("Executing a thread for: $command\n",$log_file);
					#try_exec_command($command) or die "Unable to execute command: $command\n";
					push(@commands,$command);
			}				
			execute_threads(\@commands);
		}else{
			my $in_bam_set = "";
			
			##Construct the input set using all the samples and all the chromosomes
			#Get the reads files associated with this sample
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			my $params;
			foreach my $readf_id (
				sort keys %{$distinct_readf_ids}
				)
				{
					#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
					getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
									$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
					#Build the input file name. The step will not be there in the array, hence all the suffixes 
					#(from steps performed) will be used
					my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
									$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
		
				if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
				$in_bam_set .= " -I $input_bam ";									
			}
			#Check input files
			if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." --bqsr-recal-file $bqsr_file $param_str -O $out_bam ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}

	}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
}



=head2 run_GATK4_PrintReads_persample

 Title   : run_GATK4_PrintReads_persample
 Usage   : run_GATK4_PrintReads_persample(   );

 Function: run_GATK4_PrintReads_persample
			Replaces the PrintReads for GATK4 to execute after BaseRecalibrato
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			avoid the next step of MergeSamFiles
 
 Returns : nothing
=cut
sub run_GATK4_ApplyBQSROLD{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $prev_task = shift;#Used to know the folder where to take the input bam file
	my $task = shift;
	my $steps_array = shift;      
	my $steps_array2 = shift;      
	
	
	#my $previous_task = 'refine';
	my $step = $cfg_hash->{'print_reads_step'};
	my $prev_step = $cfg_hash->{'base_recal_step'};
	
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	
	my $db_sample_name = $cfg_hash->{'db_sample_name'};
	my $db_readf_name = $cfg_hash->{'db_readf_name'};
	
	#Check from the DB if the Sorting and indexing was performed
	if ( $params->{$prev_step} == 1 ){
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
				
		#Build the base quality score recal file name using as stop point the current step
		my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};
				
		#Prepare output name for files for data recalibration
		my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		#print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";


		####### JAVA PARAMETERS	
		if (defined $cfg_hash->{'javamem_PR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
			$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_PR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_PR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
		#######END JAVA PARAMETERS
	

		#Parameters string construction
		my $param_str = " ";
		#baq
		if ( defined $cfg_hash->{'baq_PR'} ){
			$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		}
		##Number of threads to use CANNOT USE IN INDELREALIGNER
		#if ( defined $cfg_hash->{'nct_PR'} ){
			#$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		#}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}		
		my $numExec = 1;
		#After the update of GATK that the IndelRealignment should not be done the
		#per-chromosome mode is deactivated here because starting from the alignment for
		#sure we do not have that separation. Hence check if indel_real was executed
		#Getting the indel_real_step status
		my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'}
																			,$cfg_hash->{'db_sample_id'},$sample_id);	
		if ( $perchrom == 1  and $indel_real_step == 1 ){
			$numExec = $cfg_hash->{'chromosomes_num'};
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
				my $in_bam_set = "";
				
				##Construct the input set using all the samples for each chromosome
				#Get the reads files associated with this sample
				#I defined params here because I need it in the construction of output file name
				#there we need just one among all
				my $params;
				foreach my $readf_id (
					sort keys %{$distinct_readf_ids}
					)
					{
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used. The step before here is the step used because
						#the bam files did not change after base recalibration
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
						#Build i-th input name
						my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
														
						if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						$in_bam_set .= " -I $ith_in_bam ";									
					}
					#Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					#be better because all the bases are considered together
					my  $bqsr = " --bqsr-recal-file ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					#Build i-th bam output name
					my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					my $command = $jatk_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -O $ith_out_bam ";			 
					print_and_log("Executing a thread for: $command\n",$log_file);
					#try_exec_command($command) or die "Unable to execute command: $command\n";
					push(@commands,$command);
			}				
			execute_threads(\@commands);
		}else{
				my $in_bam_set = "";
				
				##Construct the input set using all the samples and all the chromosomes
				#Get the reads files associated with this sample
				#I defined params here because I need it in the construction of output file name
				#there we need just one among all
				my $params;
				foreach my $readf_id (
					sort keys %{$distinct_readf_ids}
					)
					{
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
					if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
					$in_bam_set .= " -I $input_bam ";									
				}
			#Check input files
			if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." --bqsr-recal-file $bqsr_file $param_str -O $out_bam ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}

	}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
}

=head2 run_GATK4_BaseRecalibrator_persample

 Title   : run_GATK4_BaseRecalibrator_persample
 Usage   : run_GATK4_BaseRecalibrator_persample(   );

 Function: run_GATK4_BaseRecalibrator_persample
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			
 Returns : nothing
=cut
sub run_GATK4_BaseRecalibrator_persample{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $prev_task = shift;
	my $task = shift;
	my $prev_step = shift;
  my $steps_array = shift;
  my $steps_array2 = shift;
  
  
	#my $previous_task = 'align';#Used to know the folder where to take the input bam file
	my $curr_step = $cfg_hash->{'base_recal_step'};
	#my $prev_step = $cfg_hash->{'indel_real_step'};
	
	#Get infos for subject sample
	print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	
	my $db_sample_name = $cfg_hash->{'db_sample_name'};
	my $db_readf_name = $cfg_hash->{'db_readf_name'};
	#Check from the DB if the Indel realignment was performed
	#if ( $params->{$prev_step} == 1 ){
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$group_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
		
		#Prepare output name for files for data recalibration
		my $out_bqsr = $outFolder."/".build_input_name_from_executed($params,$curr_step,
						$params->{$db_sample_name},$steps_array2)."_".$curr_step.".".$cfg_hash->{'recal_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$group_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$group_id,$cfg_hash->{'db_sample_id'},$sample_id);
												    		
		#print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}."  using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

		####### JAVA PARAMETERS	
		if (defined $cfg_hash->{'javamem_BR'} or defined $cfg_hash->{'scratch_f'}){
			#Open string
			$jatk_call .= " --java-options '";
				#Add temporary folder
				if ( defined $cfg_hash->{'scratch_f'} ){
					$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
				}	
				if (defined $cfg_hash->{'javamem_BR'}){
					$jatk_call .= " -".$cfg_hash->{'javamem_BR'}." ";
				}			
				#Close string
				$jatk_call .= "'";
			}
		#######END JAVA PARAMETERS

		
		
		
		#Parameters string construction
		my $param_str = " ";
		##Number of threads to use CANNOT USE IN INDELREALIGNER
		#if ( defined $cfg_hash->{'refine_threads'} ){
			#$param_str .= " -nct ".$cfg_hash->{'refine_threads'}." ";
		#}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}
		#if ( $cfg_hash->{'filter_mismatching_base_and_quals_BR'} eq 'YES'){
			#$param_str .= " --filter_mismatching_base_and_quals ";
		#}		
		
		#print_and_log("Known BR  ".$cfg_hash->{'known_BR'}."...\n",$log_file);#DEBUGCODE
		#Adds known indel snp to the list of known
		if ( defined $cfg_hash->{'known_BR'} ){
			my @known_indel_snp  = split($cfg_hash->{'word_sep'},$cfg_hash->{'known_BR'});
			foreach my $known (@known_indel_snp){
				#$param_str .= " -knownSites:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
				$param_str .= " --known-sites ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
			}
		}
		#print_and_log("param_str $param_str...\n",$log_file);#DEBUGCODE
		##Adds covariates to use
		#if ( defined $cfg_hash->{'covariates_BR'} ){
			#print "Covariates: ".$cfg_hash->{'covariates_BR'}."\n";
			#my @covariates  = split(",",$cfg_hash->{'covariates_BR'});
			#foreach my $covariate (@covariates){
				#$param_str .= " -cov $covariate ";
			#}
		#}		
		#if ( defined $cfg_hash->{'BR_plots'} ){
			#$param_str .= " -plots ".$cfg_hash->{'BR_plots'}." ";
		#}

		#Getting the target bed file using the folder for targets and the name contained in the database
		my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
							,$cfg_hash->{'db_analysis_id'},$group_id);
									
		my $numExec = 1;
		#After the update of GATK that the IndelRealignment should not be done the
		#per-chromosome mode is deactivated here because starting from the alignment for
		#sure we do not have that separation. Hence check if indel_real was the previous step
		if ( $perchrom == 1  and $prev_step eq $cfg_hash->{'indel_real_step'} ){
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				
			#Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
			#be better because all the bases are considered together
			#Build ith chromosome target name
			
			if (file_not_present($target_bed) > 0 ){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
			#my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'all_chr'}.".".$cfg_hash->{'bed_ext'}." ";#ROMAN FILE
			my $l_string = " --intervals $target_bed ";#USE THIS 
			
			#Build i-th bam output name
			my $bqsr_file = $outFolder."/".extract_name($out_bqsr,1).".".$cfg_hash->{'recal_ext'};
			
			my $in_bam_set = "";
			
			##Construct the input set using all the samples and all the chromosomes
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			my $params;
			foreach my $readf_id (
				sort keys %{$distinct_readf_ids}
				)
				{
					#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
					getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readf_id);
					#Build the input file name. The step will not be there in the array, hence all the suffixes 
					#(from steps performed) will be used
					my $input_bam = $inFolder."/".build_input_name_from_executed($params,$curr_step,
									$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
					#Executing one GATK command for each target chromosome using bed files
					foreach my $exec (@chroms){
					#for (my $exec = 1; $exec <= $numExec; $exec++){
							#Build i-th input name
							my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
							if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $ith_in_bam.\n",$log_file);}
							$in_bam_set .= " -I $ith_in_bam ";
					}		
				}

			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $in_bam_set $param_str $l_string -O $bqsr_file ";			 
			print_and_log("Executing a thread for: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
			push(@commands,$command);
			#Execute the thread
			execute_threads(\@commands);
		}else{
			# When perchrom is not used, it uses the target regions as intervals	

			if (file_not_present($target_bed) > 0 ){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
			$param_str .= " --intervals $target_bed ";	
			
			my $in_bam_set = "";
			##Construct the input set using all the samples and all the chromosomes
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			my $params;
			foreach my $readf_id (
				sort keys %{$distinct_readf_ids}
				)
				{
					#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
					getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
									$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
					#Build the input file name. The step will not be there in the array, hence all the suffixes 
					#(from steps performed) will be used
					my $input_bam = $inFolder."/".build_input_name_from_executed($params,$curr_step,
									$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};

				if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
				$in_bam_set .= " -I $input_bam ";				
			}	
			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." $param_str -O $out_bqsr ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}
	#}else{log_and_exit( "ERROR: Step $curr_step has not been run for sample: $sample_id.f\n ",$log_file);}
	
}



#=head2 run_GATK4_PrintReads_persample

 #Title   : run_GATK4_PrintReads_persample
 #Usage   : run_GATK4_PrintReads_persample(   );

 #Function: run_GATK4_PrintReads_persample
			#Replaces the PrintReads for GATK4 to execute after BaseRecalibrato
			#Is a tool from GATK which permits to realign the reads towards well assessed indels
			#Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			#and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			#recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			#avoid the next step of MergeSamFiles
 
 #Returns : nothing
#=cut
#sub run_GATK4_ApplyBQSR{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $analysisid = shift;
	#my $out_bam = shift;
	#my $log_file = shift;
	#my $prev_task = shift;#Used to know the folder where to take the input bam file
	#my $task = shift;
	#my $steps_array = shift;      
	#my $steps_array2 = shift;      
	
	
	##my $previous_task = 'refine';
	#my $step = $cfg_hash->{'print_reads_step'};
	#my $prev_step = $cfg_hash->{'base_recal_step'};

	##Get infos for subject sample
	##print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	#my $params;
	#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						#$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
		
	#my $db_sample_name = $cfg_hash->{'db_sample_name'};
	#my $db_readf_name = $cfg_hash->{'db_readf_name'};
	
	##Check from the DB if the Sorting and indexing was performed
	#if ( $params->{$prev_step} == 1 ){
		##inFolder for the bam file
		#my $inFolder = $cfg_hash->{$analysisid.'_'.$prev_task.'_out_f'};
		##Out folder
		#my $outFolder = $cfg_hash->{$analysisid.'_'.$task.'_out_f'};
				
		##Build the base quality score recal file name using as stop point the current step
		#my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						#$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};

		##Verifying that there has been an execution per chromosome
		#my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    #,$cfg_hash->{'db_analysis_id'},$analysisid);

		##Here we get the read file ids associated with this sample. They will be used for the input file names
		#my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						#,$cfg_hash->{'db_analysis_id'},$analysisid,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		##print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		##Set the Java call
		#my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";


		######## JAVA PARAMETERS	
		#if (defined $cfg_hash->{'javamem_PR'} or defined $cfg_hash->{'scratch_f'}){
		##Open string
			#$jatk_call .= " --java-options '";
			##Add temporary folder
			#if ( defined $cfg_hash->{'scratch_f'} ){
				#$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			#}	
			#if (defined $cfg_hash->{'javamem_PR'}){
				#$jatk_call .= " -".$cfg_hash->{'javamem_PR'}." ";
			#}			
			##Close string
			#$jatk_call .= "'";
		#}
		########END JAVA PARAMETERS
	

		##Parameters string construction
		#my $param_str = " ";
		##baq
		#if ( defined $cfg_hash->{'baq_PR'} ){
			#$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		#}
		###Number of threads to use CANNOT USE IN INDELREALIGNER
		##if ( defined $cfg_hash->{'nct_PR'} ){
			##$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		##}
		#if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --fix_misencoded_quality_scores ";
		#}
		#if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --allow_potentially_misencoded_quality_scores ";
		#}		
		#my $numExec = 1;
		##After the update of GATK that the IndelRealignment should not be done the
		##per-chromosome mode is deactivated here because starting from the alignment for
		##sure we do not have that separation. Hence check if indel_real was executed
		##Getting the indel_real_step status
		#my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'}
																			#,$cfg_hash->{'db_sample_id'},$sample_id);	
		#if ( $perchrom == 1  and $indel_real_step == 1 ){
			#$numExec = $cfg_hash->{'chromosomes_num'};
				###Using also Chromosomex X and Y that Romans did not use!!!
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				##An array will contain all the commands to be executed
				#my @commands = ();
				
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
				##for (my $exec = 1; $exec <= $numExec; $exec++){
				#my $in_bam_set = "";
				
				###Construct the input set using all the samples for each chromosome
				##Get the reads files associated with this sample
				##I defined params here because I need it in the construction of output file name
				##there we need just one among all
				#my $params;
				#foreach my $readf_id (
					#sort keys %{$distinct_readf_ids}
					#)
					#{
						##print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										#$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						##Build the input file name. The step will not be there in the array, hence all the suffixes 
						##(from steps performed) will be used. The step before here is the step used because
						##the bam files did not change after base recalibration
						#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										#$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
						##Build i-th input name
						#my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
														
						#if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						#$in_bam_set .= " -I $ith_in_bam ";									
					#}
					##Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					##be better because all the bases are considered together
					#my  $bqsr = " --bqsr-recal-file ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					##Build i-th bam output name
					#my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					#my $command = $jatk_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -O $ith_out_bam ";			 
					#print_and_log("Executing a thread for: $command\n",$log_file);
					##try_exec_command($command) or die "Unable to execute command: $command\n";
					#push(@commands,$command);
			#}				
			#execute_threads(\@commands);
		#}else{
				#my $in_bam_set = "";
				
				###Construct the input set using all the samples and all the chromosomes
				##Get the reads files associated with this sample
				##I defined params here because I need it in the construction of output file name
				##there we need just one among all
				#my $params;
				#foreach my $readf_id (
					#sort keys %{$distinct_readf_ids}
					#)
					#{
						##print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										#$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						##Build the input file name. The step will not be there in the array, hence all the suffixes 
						##(from steps performed) will be used
						#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
										#$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
					#if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
					#$in_bam_set .= " -I $input_bam ";									
				#}
			##Check input files
			#if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			#my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." --bqsr-recal-file $bqsr_file $param_str -O $out_bam ";			 
			#print_and_log("Executing command: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
		#}

	#}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
#}



#=head2 run_GATK4_PrintReads_persample

 #Title   : run_GATK4_PrintReads_persample
 #Usage   : run_GATK4_PrintReads_persample(   );

 #Function: run_GATK4_PrintReads_persample
			#Replaces the PrintReads for GATK4 to execute after BaseRecalibrato
			#Is a tool from GATK which permits to realign the reads towards well assessed indels
			#Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			#and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			#recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			#avoid the next step of MergeSamFiles
 
 #Returns : nothing
#=cut
#sub run_GATK4_ApplyBQSROLD{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $analysisid = shift;
	#my $log_file = shift;
	#my $prev_task = shift;#Used to know the folder where to take the input bam file
	#my $task = shift;
	#my $steps_array = shift;      
	#my $steps_array2 = shift;      
	
	
	##my $previous_task = 'refine';
	#my $step = $cfg_hash->{'print_reads_step'};
	#my $prev_step = $cfg_hash->{'base_recal_step'};
	
	##Get infos for subject sample
	##print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	#my $params;
	#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						#$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	
	#my $db_sample_name = $cfg_hash->{'db_sample_name'};
	#my $db_readf_name = $cfg_hash->{'db_readf_name'};
	
	##Check from the DB if the Sorting and indexing was performed
	#if ( $params->{$prev_step} == 1 ){
		##inFolder for the bam file
		#my $inFolder = $cfg_hash->{$analysisid.'_'.$prev_task.'_out_f'};
		##Out folder
		#my $outFolder = $cfg_hash->{$analysisid.'_'.$task.'_out_f'};
				
		##Build the base quality score recal file name using as stop point the current step
		#my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						#$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};
				
		##Prepare output name for files for data recalibration
		#my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
						#$params->{$db_sample_name},$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};

		##Verifying that there has been an execution per chromosome
		#my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    #,$cfg_hash->{'db_analysis_id'},$analysisid);

		##Here we get the read file ids associated with this sample. They will be used for the input file names
		#my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						#,$cfg_hash->{'db_analysis_id'},$analysisid,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		##print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		##Set the Java call
		#my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";


		######## JAVA PARAMETERS	
		#if (defined $cfg_hash->{'javamem_PR'} or defined $cfg_hash->{'scratch_f'}){
		##Open string
			#$jatk_call .= " --java-options '";
			##Add temporary folder
			#if ( defined $cfg_hash->{'scratch_f'} ){
				#$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			#}	
			#if (defined $cfg_hash->{'javamem_PR'}){
				#$jatk_call .= " -".$cfg_hash->{'javamem_PR'}." ";
			#}			
			##Close string
			#$jatk_call .= "'";
		#}
		########END JAVA PARAMETERS
	

		##Parameters string construction
		#my $param_str = " ";
		##baq
		#if ( defined $cfg_hash->{'baq_PR'} ){
			#$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		#}
		###Number of threads to use CANNOT USE IN INDELREALIGNER
		##if ( defined $cfg_hash->{'nct_PR'} ){
			##$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		##}
		#if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --fix_misencoded_quality_scores ";
		#}
		#if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --allow_potentially_misencoded_quality_scores ";
		#}		
		#my $numExec = 1;
		##After the update of GATK that the IndelRealignment should not be done the
		##per-chromosome mode is deactivated here because starting from the alignment for
		##sure we do not have that separation. Hence check if indel_real was executed
		##Getting the indel_real_step status
		#my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'}
																			#,$cfg_hash->{'db_sample_id'},$sample_id);	
		#if ( $perchrom == 1  and $indel_real_step == 1 ){
			#$numExec = $cfg_hash->{'chromosomes_num'};
				###Using also Chromosomex X and Y that Romans did not use!!!
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				##An array will contain all the commands to be executed
				#my @commands = ();
				
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
				##for (my $exec = 1; $exec <= $numExec; $exec++){
				#my $in_bam_set = "";
				
				###Construct the input set using all the samples for each chromosome
				##Get the reads files associated with this sample
				##I defined params here because I need it in the construction of output file name
				##there we need just one among all
				#my $params;
				#foreach my $readf_id (
					#sort keys %{$distinct_readf_ids}
					#)
					#{
						##print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										#$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						##Build the input file name. The step will not be there in the array, hence all the suffixes 
						##(from steps performed) will be used. The step before here is the step used because
						##the bam files did not change after base recalibration
						#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										#$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
						##Build i-th input name
						#my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
														
						#if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						#$in_bam_set .= " -I $ith_in_bam ";									
					#}
					##Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					##be better because all the bases are considered together
					#my  $bqsr = " --bqsr-recal-file ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					##Build i-th bam output name
					#my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					#my $command = $jatk_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -O $ith_out_bam ";			 
					#print_and_log("Executing a thread for: $command\n",$log_file);
					##try_exec_command($command) or die "Unable to execute command: $command\n";
					#push(@commands,$command);
			#}				
			#execute_threads(\@commands);
		#}else{
				#my $in_bam_set = "";
				
				###Construct the input set using all the samples and all the chromosomes
				##Get the reads files associated with this sample
				##I defined params here because I need it in the construction of output file name
				##there we need just one among all
				#my $params;
				#foreach my $readf_id (
					#sort keys %{$distinct_readf_ids}
					#)
					#{
						##print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										#$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						##Build the input file name. The step will not be there in the array, hence all the suffixes 
						##(from steps performed) will be used
						#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
										#$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
					#if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
					#$in_bam_set .= " -I $input_bam ";									
				#}
			##Check input files
			#if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			#my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." --bqsr-recal-file $bqsr_file $param_str -O $out_bam ";			 
			#print_and_log("Executing command: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
		#}

	#}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
#}


=head2 run_GATK4_PrintReads_persample

 Title   : run_GATK4_PrintReads_persample
 Usage   : run_GATK4_PrintReads_persample(   );

 Function: run_GATK4_PrintReads_persample
			
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			avoid the next step of MergeSamFiles
 
 Returns : nothing
=cut
sub run_GATK4_PrintReads_persample{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $prev_task = shift;#Used to know the folder where to take the input bam file
	my $task = shift;
	my $steps_array = shift;      
	my $steps_array2 = shift;      
	
	
	#my $previous_task = 'refine';
	my $step = $cfg_hash->{'print_reads_step'};
	my $prev_step = $cfg_hash->{'base_recal_step'};
	
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	
	my $db_sample_name = $cfg_hash->{'db_sample_name'};
	my $db_readf_name = $cfg_hash->{'db_readf_name'};
	
	#Check from the DB if the Sorting and indexing was performed
	if ( $params->{$prev_step} == 1 ){
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$group_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
				
		#Build the base quality score recal file name using as stop point the current step
		my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};
				
		#Prepare output name for files for data recalibration
		my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$group_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$group_id,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		#print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";


####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_PR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_PR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_PR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	

		#Parameters string construction
		my $param_str = " ";
		#baq
		if ( defined $cfg_hash->{'baq_PR'} ){
			$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		}
		##Number of threads to use CANNOT USE IN INDELREALIGNER
		#if ( defined $cfg_hash->{'nct_PR'} ){
			#$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		#}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}		
		my $numExec = 1;
		#After the update of GATK that the IndelRealignment should not be done the
		#per-chromosome mode is deactivated here because starting from the alignment for
		#sure we do not have that separation. Hence check if indel_real was executed
		#Getting the indel_real_step status
		my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'}
																			,$cfg_hash->{'db_sample_id'},$sample_id);	
		if ( $perchrom == 1  and $indel_real_step == 1 ){
			$numExec = $cfg_hash->{'chromosomes_num'};
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
				my $in_bam_set = "";
				
				##Construct the input set using all the samples for each chromosome
				#Get the reads files associated with this sample
				#I defined params here because I need it in the construction of output file name
				#there we need just one among all
				my $params;
				foreach my $readf_id (
					sort keys %{$distinct_readf_ids}
					)
					{
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used. The step before here is the step used because
						#the bam files did not change after base recalibration
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
						#Build i-th input name
						my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
														
						if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						$in_bam_set .= " -I $ith_in_bam ";									
					}
					#Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					#be better because all the bases are considered together
					my  $bqsr = " -BQSR ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					#Build i-th bam output name
					my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
					my $command = $jatk_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -o $ith_out_bam ";			 
					print_and_log("Executing a thread for: $command\n",$log_file);
					#try_exec_command($command) or die "Unable to execute command: $command\n";
					push(@commands,$command);
			}				
			execute_threads(\@commands);
		}else{
				my $in_bam_set = "";
				
				##Construct the input set using all the samples and all the chromosomes
				#Get the reads files associated with this sample
				#I defined params here because I need it in the construction of output file name
				#there we need just one among all
				my $params;
				foreach my $readf_id (
					sort keys %{$distinct_readf_ids}
					)
					{
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
					if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
					$in_bam_set .= " -I $input_bam ";									
				}
			#Check input files
			if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." -BQSR $bqsr_file $param_str -o $out_bam ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}

	}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
}


#=head2 run_GATK4_BaseRecalibrator_perlane

 #Title   : run_GATK4_BaseRecalibrator_perlane
 #Usage   : run_GATK4_BaseRecalibrator_perlane(   );

 #Function: run_GATK4_BaseRecalibrator_perlane
			#Is a tool from GATK which permits to realign the reads towards well assessed indels
			
 #Returns : nothing
#=cut
#sub run_GATK4_BaseRecalibrator_perlane{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $group_id = shift;
	#my $log_file = shift;
	#my $task = shift;
  #my $steps_array = shift;
        
	##my $previous_task = 'align';#Used to know the folder where to take the input bam file
	#my $curr_step = $cfg_hash->{'base_recal_step'};
	#my $prev_step = $cfg_hash->{'indel_real_step'};
	
	##Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);
	#my $params;
	#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						#$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$sample_id);
	
	#my $main_name = $cfg_hash->{'db_readf_name'};
	##Check from the DB if the Sorting and indexing was performed
	#if ( $params->{$prev_step} == 1 ){
		##inFolder for the bam file
		#my $inFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
		##Out folder
		#my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
		##Build the input file name using as stop point the indel realigner step
		#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$curr_step,
						#$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};
		
		##Prepare output name for files for data recalibration
		#my $out_bqsr = $outFolder."/".build_input_name_from_executed($params,$curr_step,
						#$params->{$main_name},$steps_array)."_".$curr_step.".".$cfg_hash->{'recal_ext'};

		##Verifying that there has been an execution per chromosome
		#my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    #,$cfg_hash->{'db_analysis_id'},$group_id);
						    		
		#print_and_log( "Starting $prog_used for $input_bam using GATK tools..\n ",$log_file);
		
		##Set the Java call
		#my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
		##Set the java memory
		#if (defined $cfg_hash->{'javamem_BR'}){
			#$jatk_call .= " -".$cfg_hash->{'javamem_BR'}." ";
		#}
		#
		#$jatk_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

		##Parameters string construction
		#my $param_str = " ";
		##Number of threads to use CANNOT USE IN INDELREALIGNER
		#if ( defined $cfg_hash->{'nct_BR'} ){
			#$param_str .= " -nct ".$cfg_hash->{'nct_BR'}." ";
		#}
		#if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --fix_misencoded_quality_scores ";
		#}
		#if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			#$param_str .= " --allow_potentially_misencoded_quality_scores ";
		#}
		##Adds known indel snp to the list of known
		#if ( defined $cfg_hash->{'known_BR'} ){
			#my @known_indel_snp  = split($cfg_hash->{'word_sep'},$cfg_hash->{'known_BR'});
			#foreach my $known (@known_indel_snp){
				##$param_str .= " -knownSites:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
				#$param_str .= " -knownSites ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
			#}
		#}
		##Adds covariates to use
		#if ( defined $cfg_hash->{'covariates_BR'} ){
			#print "Covariates: ".$cfg_hash->{'covariates_BR'}."\n";
			#my @covariates  = split(",",$cfg_hash->{'covariates_BR'});
			#foreach my $covariate (@covariates){
				#$param_str .= " -cov $covariate ";
			#}
		#}		
		
		#my $numExec = 1;
		#if ( $perchrom == 1){
				###Using also Chromosomex X and Y that Romans did not use!!!
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				##An array will contain all the commands to be executed
				#my @commands = ();
				

			##Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
			##be better because all the bases are considered together
			#if ( $cfg_hash->{'BQSR_uniq'} eq 'YES'){
				##Build ith chromosome target name
				#my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'all_chr'}.".".$cfg_hash->{'bed_ext'}." ";
				##Build i-th bam output name
				#my $bqsr_file = $outFolder."/".extract_name($out_bqsr,1).".".$cfg_hash->{'recal_ext'};
				
				#my $in_bam_set = "";
				##QUI DOVREI CONTROLLARE SE I FILE BED SONO POI EFFETTIVAMENTE QUELLI...#TODO
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
				##for (my $exec = 1; $exec <= $numExec; $exec++){
						##Build i-th input name
						#my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
						#if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $ith_in_bam.\n",$log_file);}
						#$in_bam_set .= " -I $ith_in_bam ";
				#}				
				##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
				#my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $in_bam_set $param_str $l_string -o $bqsr_file ";			 
				#print_and_log("[gr:$group_id]Executing a thread for: $command\n",$log_file);
				##try_exec_command($command) or die "Unable to execute command: $command\n";
				#push(@commands,$command);
			#}
			##Analysis with a BQSR for each BAM per chromosome. This kind of analysis should
			##be faster but the bases are not considered together			
			#else{
			 
				##QUI DOVREI CONTROLLARE SE I FILE BED SONO POI EFFETTIVAMENTE QUELLI...#TODO
				##Executing one GATK command for each target chromosome using bed files
				#for (my $exec = 1; $exec <= $numExec; $exec++){
						##Build i-th input name
						#my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
						#if (file_not_present($ith_in_bam)  > 0){ log_and_exit( "Cannot proceed with $prog_used! Check: $ith_in_bam.\n",$log_file);}
						##Build ith chromosome target name
						#my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
						##Build i-th bam output name
						#my $ith_bqsr = $outFolder."/".extract_name($out_bqsr,1)."_$exec.".$cfg_hash->{'recal_ext'};
						##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
						#my $command = $jatk_call." -I $ith_in_bam -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_bqsr ";			 
						#print_and_log("Executing a thread for: $command\n",$log_file);
						##try_exec_command($command) or die "Unable to execute command: $command\n";
						#push(@commands,$command);
				#}				
			#}

			#execute_threads(\@commands);
		#}else{
			#if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
			##print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			#my $command = $jatk_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bqsr ";			 
			#print_and_log("[gr:$group_id] Executing command: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
		#}
	#}else{log_and_exit( "ERROR: Step $curr_step has not been run for sample: $sample_id.f\n ",$log_file);}
	
#}




=head2 run_GATK4_HaplotypeCaller_parallel

 Title   : run_GATK4_HaplotypeCaller_parallel
 Usage   : run_GATK4_HaplotypeCaller_parallel(   );

 Function: run_GATK4_HaplotypeCaller_parallel
			Is a tool from GATK which permits to call variants. Haplotype caller works
			per-sample hence if until this moment we have worked with files separated
			by lanes a step of Merge is necessary.
			A check for this is made when the distinct samples are fetched if they exist.
			In that case all the steps for each file should be completed.
			
			I have written only a "parallel" version of GATK4 HaplotypeCaller
			
 Returns : nothing
=cut
sub run_GATK4_HaplotypeCaller_parallel{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $in_fold_suff = shift;
	my $out_fold_suff = shift;#Used to know the folder where to take the input bam files
	my $step = shift;
	my $step_needed  = shift;
	my $steps_array = shift;
	my $steps_array2 = shift;
	my $log_file = shift;

	#my $previous_task = '';
	my $steps_arr1_in_use;
	my $paramsIn;
	my $paramsOut;
	my $main_name = $cfg_hash->{'db_sample_name'};

	#Case IR only or no refinement step: there is a single file and only IndelRealignment or no refinement step was executed: 
	#construct the name using readfile level fields
	if ( $step_needed eq $cfg_hash->{'indel_real_step'} or $step_needed eq $cfg_hash->{'sort_idx_step'} ){
		$steps_arr1_in_use = $steps_array;
		my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
								$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
		foreach my $dist_sample ( keys %{$distinct_samples}){
			getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$dist_sample);
		}
	}	
	#CASE: IR-BR or IR-BR-M or IR-M: there are multiple read files or the last step is BR
	#construct the name using sample level fields
	elsif ($step_needed eq $cfg_hash->{'mrdup_groups_step'} or $step_needed eq $cfg_hash->{'print_reads_step'}){
		$steps_arr1_in_use = $steps_array2;
		getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	}
	
	#Get the sample parameters from the the db for the output file
	getSampleConfiguration_locked(\$paramsOut,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);		
	
	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
					    $cfg_hash->{'db_analysis_id'},$analysis_id);
																	
	#InFolder is the alignment one at this stage
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};

	#Getting the pedigree file
	my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);	


	#Getting the flag indicating if the joint genotype has to be run
	my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $out_ext;
	#If it has been run, take the output form the genotype, else from the varcalling
	if ($dogenot){
		$out_ext = $cfg_hash->{'gvcf_ext'};
	}else{
		$out_ext = $cfg_hash->{'vcf_ext'};
	}
																				
	#Set the input and output file names
	my $input_file = $inFolder."/".build_input_name_from_executed($paramsIn,$step,$sample_name,$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'};
	my $out_file = $outFolder."/".build_input_name_from_executed($paramsOut,$step,$sample_name,$steps_array2)."_".$step.".".$out_ext;
	
	#print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Prepare output name
	my $param_str = " ";	
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_HAPC'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_HAPC'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_HAPC'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	##nct: number of processors to use
	#if ( defined $cfg_hash->{'varcall_threads'} ){
					#$param_str .= " -nct ".$cfg_hash->{'varcall_threads'};
	#}
	if ( defined $cfg_hash->{'emitRefConfidence'} and $dogenot == 1){
					$param_str .= " --emit-ref-confidence ".$cfg_hash->{'emitRefConfidence'};
	}
	if ( defined $cfg_hash->{'variant_index_type'} ){
					$param_str .= " --variant_index_type ".$cfg_hash->{'variant_index_type'};
	}
	if ( defined $cfg_hash->{'variant_index_parameter'} ){
					$param_str .= " --variant_index_parameter ".$cfg_hash->{'variant_index_parameter'};
	}
	if ( defined $cfg_hash->{'HC_stand_call_conf'} and correct_type($cfg_hash->{'HC_stand_call_conf'},"positiveint") ){
					$param_str .= " -stand-call-conf ".$cfg_hash->{'HC_stand_call_conf'};
	}
	if ( defined $cfg_hash->{'HC_stand_emit_conf'} and correct_type($cfg_hash->{'HC_stand_emit_conf'},"positiveint") ){
					$param_str .= " -stand_emit_conf ".$cfg_hash->{'HC_stand_emit_conf'};
	}
	if ( defined $cfg_hash->{'HC_min_base_quality_score'} and correct_type($cfg_hash->{'HC_min_base_quality_score'},"positiveint") ){
					$param_str .= " --min_base_quality_score ".$cfg_hash->{'HC_min_base_quality_score'};
	}
	if ( defined $cfg_hash->{'HC_dontUseSoftClippedBases'} and  $cfg_hash->{'HC_dontUseSoftClippedBases'} eq 'YES' ){
					$param_str .= " --dontUseSoftClippedBases ";
	}
	if ( defined $cfg_hash->{'hapc_dbsnp'} ){
					$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	}	
	

	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'}
									,$cfg_hash->{'db_analysis_id'},$analysis_id);			

	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);	
	#Get target bed
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}		
	}else{$target_bed = "-";}

	if ( defined $cfg_hash->{'hapc_dbsnp'} ){
		$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	}		
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'hapc_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'hapc_annotation'});
		foreach my $annotation (@annotations){
			$param_str .= " -A $annotation ";
		}
	}
	
	my $numExec = 1;
	if ( $perchrom == 1 ){
		##Using also Chromosomex X and Y
		$numExec = $cfg_hash->{'chromosomes_num'};
		my @chroms = 1..$numExec;
		push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
					
		#An array will contain all the commands to be executed
		my @commands = ();
		my $out_files_str = "";
		
		#Executing one GATK command for each target chromosome using bed files
		foreach my $exec (@chroms){
			#Build ith chromosome target name
			my $l_string = "";
			#WES
			if (  $target_bed ne "-" ) {
				$l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
			}else{
				$l_string = " -L ".$cfg_hash->{'chr_suff'}."$exec ";
			}
			
			my $ith_out = $outFolder."/".extract_name($out_file,1)."_$exec.".$out_ext;
			$out_files_str .= $ith_out.",";	
			
			my $ith_input = $inFolder."/".extract_name($input_file,1)."_".$exec.".".$cfg_hash->{'bam_ext'};
			if (file_not_present($ith_input) > 0 ){ die "Cannot proceed with $prog_used! Check: $ith_input.\n";}
			#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
			my $command = $jatk_call." -I $ith_input -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_out ";			 
			print_and_log("Executing command: $command \n",$log_file);
			push(@commands,$command);
		}
		execute_threads(\@commands);
		return "no_depend",$out_files_str,$out_file
	}else{
		if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
		#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
		my $command = $jatk_call." -I $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -O $out_file ";			 
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
	}
}



=head2 run_GATK4_CatVariants

 Title   : run_GATK4_CatVariants
 Usage   : run_GATK4_CatVariants(   );

 Function: run_GATK4_CatVariants
			Is a tool from GATK which permits to merge together gvcf files from
			Haplotype caller so that we will have only one VCF file per sample
			
 Returns : nothing
=cut
sub run_GATK4_CatVariants{
	my $cfg_hash = shift;
	my $prog_used = "GatherVcfs";
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;
  my $step = shift;
	my $steps_array = shift;
	
	my $previous_task = 'varcall';
	my $step_needed = $cfg_hash->{'varcall_step'};

	#Getting the group name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'}
						,$cfg_hash->{'db_sample_id'},$sample_id);		
																		
	my $params;
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
				    $cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_sample_table'},
																					$cfg_hash-> {'db_sample_id'},$sample_id);

	my $prev_step = $cfg_hash->{'varcall_step'};
	#Getting the group name
	my $prev_step_done = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$prev_step 
						,$cfg_hash->{'db_sample_id'},$sample_id);
	if ( $prev_step_done == 1){
	
		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$group_id);

		#Getting the flag indicating if the joint genotype has to be run
		my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
						$cfg_hash->{'db_analysis_id'},$group_id);
		my $out_ext;
		#If it has been run, take the output form the genotype, else from the varcalling
		if ($dogenot){
			$out_ext = $cfg_hash->{'gvcf_ext'};
		}else{
			$out_ext = $cfg_hash->{'vcf_ext'};
		}
											
		#InFolder is the variant call at this stage
		my $inFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
		my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
			
		#Set the unique output
		my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
						$sample_name,$steps_array)."_".$step.".".$out_ext;

		#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Prepare output name
		my $param_str = " ";	
		
		#Set the Java call
		my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_CV'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_CV'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_CV'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS


		#$jatk_call .= " -cp ".$cfg_hash->{'gatk_path'}." ".$cfg_hash->{'gatk_tools_path'}.".$prog_used ";

		if ( $cfg_hash->{'assumeSorted'} eq 'YES'){
			$param_str .= " -assumeSorted ";
		}

		#Constructing input string with files to be merged 
		my $vars_to_cat = "";
		my @f_to_remove = ();
		my $numExec = 1;
		if ( $perchrom == 1){
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));

				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
					#Set the input file name to be used as component for the single files to be merged
					my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,$sample_name,$steps_array)."_$exec.".$out_ext;
					if (file_not_present($input_file) > 0 ){ log_and_exit("Cannot proceed with $prog_used! Check: $input_file.\n",$log_file);}
					$vars_to_cat .=  " -I $input_file ";
					push(@f_to_remove,$input_file);
			}
		}
		#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
		my $command = $jatk_call." $vars_to_cat --REFERENCE ".$cfg_hash->{'hum_ref'}." $param_str --OUTPUT $outFile ";			 
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";	
		
		if ( $cfg_hash->{'remove_temp'} eq 'YES'){
			foreach my $f (@f_to_remove){
				delete_file($f);
			}		
		}
	}else{log_and_exit("ERROR: step $prev_step_done has not been run for $sample_id. Cannot run $step process..\n",$log_file);}
}



=head2 run_GATK4_GenomicsDBImport_parallel

Title   : run_GATK4_GenomicsDBImport_parallel
Usage   : run_GATK4_GenomicsDBImport_parallel;

Function:  run_GATK4_GenomicsDBImport

		   Import single-sample GVCFs into one GVCF before joint genotyping.
		   It is able to divide the job in multiple jobs (returning the job ids)
		   and also to execute a single command.
			
Returns : The genomics db fo GenotypeGVCFs

=cut			
sub run_GATK4_GenomicsDBImport_parallel{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $sample_id = shift;
	my $task = shift;
	my $step = shift;
	my $steps_array2 = shift;
	my $log_file = shift;
	
	
	print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file); ##DEBUGCODE
	#defining $previous_task
	my $previous_task = 'varcall';

	print_and_log("Previous task is $previous_task using GATK tools..\n ",$log_file);##DEBUGCODE

	#Getting the analysis name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);	

																													
	#Here we get the samples ids associated with that analysis id
	print_and_log( "Here we get the samples ids associated with that analysis id..\n ",$log_file);##DEBUGCODE	
	my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
					    $cfg_hash->{'db_analysis_id'},$analysis_id);

	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);	
	#Get target bed
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}		
	}else{$target_bed = "-";}
						    							
	#input vcfs folder 
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
    #Set the output folder for chr for final vcf 
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	#Prepare parameters string
	my $param_str = " --overwrite-existing-genomicsdb-workspace ";		
	#Constructing input string with files to be merged 
	my $input_file_str = "";

	my $params;	
	if ( (scalar keys %{$distinct_samples}) > 1){
		#Check from the DB if the previous step has been performed for all the samples
		foreach my $dist_sample ( keys %{$distinct_samples}){
				getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$dist_sample);
				$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
									$params->{$cfg_hash->{'db_sample_name'}},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
				#print_and_log("Sample $input_file_str  will be used...\n",$log_file);#DEBUGCODE
		}
	}else{
		#Get infos for the single sample associated to this group from the 'samples' table
		print_and_log("Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
		$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$cfg_hash->{'db_sample_name'}},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	}

	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_CV'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}	
		if (defined $cfg_hash->{'javamem_CV'}){
			$jatk_call .= " -".$cfg_hash->{'javamem_CV'}." ";
		}			
		#Close string
		$jatk_call .= "'";
	}
	#######END JAVA PARAMETERS
	
	##Generation of the commands to be launched in parallel using the chromosomes
	#Using also Chromosomex X and Y
	my $numExec = $cfg_hash->{'chromosomes_num'};
	my @chroms = 1..$numExec;
	push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	#An array will contain all the commands to be executed
	my @commands = ();
    				
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	
	#Execution of parallel jobs divided per chromosome
	if ( $perchrom == 1 ){
		#Check if some chromosomes is so huge to use splitted versions for the jobs them in two different parts
		my @chroms_to_split = ();
		if ( $cfg_hash->{'chroms_to_split'} ne 'none' ){
			my @chroms_to_split = split(",",$cfg_hash->{'chroms_to_split'});
			print_and_log( " Chromosomes to split: ".@chroms_to_split." ..\n ",$log_file);
		}
		#Executing one GATK command for each target chromosome using bed files
		foreach my $exec (@chroms){
			#Check if this chromosome has been splitted to separate job
			print_and_log( " Check if this chromosome has been splitted $exec ..\n ",$log_file);
			#if ( grep $_ == $exec , @chroms_to_split ){
			if ( $exec == 1 ){
				print_and_log( " This chromosome $exec will be split..\n ",$log_file);
				my @chr_parts = (1,2);
				foreach my $chr_part (@chr_parts){
					my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec\_$chr_part";
					#folder created
					print_and_log( " The genomics db folder out for chr$exec part $chr_part is $genomics_db_f ..\n ",$log_file);
					#targetbed  for both parts chr1
					my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec\_$chr_part.".$cfg_hash->{'bed_ext'}." ";			 
					my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $input_file_str --genomicsdb-workspace-path $genomics_db_f $l_string $param_str";
					push(@commands,$command);	
				}
			}else{		
				#folder for genomics db
				my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec";
				print_and_log( " The genomics db folder out is $genomics_db_f ..\n ",$log_file);
				#target
				my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." "; 
				my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $input_file_str --genomicsdb-workspace-path $genomics_db_f  $l_string $param_str";
				print_and_log( " The command is $command ..\n ",$log_file);
				#print_and_log("Preparing command: $command\n",$log_file);
				push(@commands,$command);
			}
		}	
		#print_and_log("END",$log_file);
		#Run a different job for each command
		my $jobs_to_wait_dbi = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$task,$task,$log_file);
		#Return the job ids for each job and the list of output file and the final output name
		return $jobs_to_wait_dbi;		
	}else{
		#folder for genomics db
		my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'};		
		print_and_log( " The genomics db folder out is $genomics_db_f ..\n ",$log_file);
		#target
		my $l_string = " --intervals $target_bed "; 
		my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $input_file_str --genomicsdb-workspace-path $genomics_db_f  $l_string $param_str";
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
		#Return the command (Just to return something...)
		return $command;
	}


}

=head2 run_GATK4_GenomicsDBImport

Title   : run_GATK4_GenomicsDBImport
Usage   : run_GATK4_GenomicsDBImport;

Function:  run_GATK4_GenomicsDBImport

		   Import single-sample GVCFs into one GVCF before joint genotyping
			
Returns : Return one GVCF before joint genotyping   (????)

=cut			
sub run_GATK4_GenomicsDBImport{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $sample_id = shift;
	my $task = shift;
	my $step = shift;
	my $steps_array2 = shift;
	my $log_file = shift;
	
	
	print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file); ##DEBUGCODE
	#defining $previous_task
	my $previous_task = 'varcall';

	print_and_log("Previous task is $previous_task using GATK tools..\n ",$log_file);##DEBUGCODE

	#Getting the analysis name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);	

																													
	#Here we get the samples ids associated with that analysis id
	print_and_log( "Here we get the samples ids associated with that analysis id..\n ",$log_file);##DEBUGCODE	
	my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);
							
	#input vcfs folder 
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
    #Set the output folder for chr for final vcf 
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	#Prepare parameters string
	my $param_str = " --overwrite-existing-genomicsdb-workspace ";		
	#Constructing input string with files to be merged 
	my $input_file_str = "";

	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'}
									,$cfg_hash->{'db_analysis_id'},$analysis_id);			

	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);	
	#Get target bed
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}		
	}else{$target_bed = "-";}
	
	my $params;	
	if ( (scalar keys %{$distinct_samples}) > 1){
		#Check from the DB if the previous step has been performed for all the samples
		foreach my $dist_sample ( keys %{$distinct_samples}){
				getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$dist_sample);
				$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
									$params->{$cfg_hash->{'db_sample_name'}},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
				#print_and_log("Sample $input_file_str  will be used...\n",$log_file);#DEBUGCODE
		}
	}else{
		#Get infos for the single sample associated to this group from the 'samples' table
		print_and_log("Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
		$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$cfg_hash->{'db_sample_name'}},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	}

	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_CV'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}	
		if (defined $cfg_hash->{'javamem_CV'}){
			$jatk_call .= " -".$cfg_hash->{'javamem_CV'}." ";
		}			
		#Close string
		$jatk_call .= "'";
	}
	#######END JAVA PARAMETERS
	
	##Generation of the commands to be launched in parallel using the chromosomes
	#Using also Chromosomex X and Y
	my $numExec = $cfg_hash->{'chromosomes_num'};
	my @chroms = 1..$numExec;
	push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	#An array will contain all the commands to be executed
	my @commands = ();
    				
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	
	#Check if some chromosomes is so huge to use splitted versions for the jobs them in two different parts
	my @chroms_to_split = ();
	if ( $cfg_hash->{'chroms_to_split'} ne 'none' ){
		my @chroms_to_split = split(",",$cfg_hash->{'chroms_to_split'});
		print_and_log( " Cromosomes to split: ".@chroms_to_split." ..\n ",$log_file);
	}
	#Executing one GATK command for each target chromosome using bed files
	foreach my $exec (@chroms){
		#Check if this chromosome has been splitted to separate job
		print_and_log( " Check if this chromosome has been splitted $exec ..\n ",$log_file);
		#if ( grep $_ == $exec , @chroms_to_split ){
		if ( $exec == 1 ){
			print_and_log( " This chromosome $exec will be split..\n ",$log_file);
			my @chr_parts = (1,2);
			foreach my $chr_part (@chr_parts){
				my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec\_$chr_part";
				#folder created
				print_and_log( " The genomics db folder out for chr$exec part $chr_part is $genomics_db_f ..\n ",$log_file);
				#targetbed  for both parts chr1
				my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec\_$chr_part.".$cfg_hash->{'bed_ext'}." ";			 
				my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $input_file_str --genomicsdb-workspace-path $genomics_db_f $l_string $param_str";
				push(@commands,$command);	
			}
		}else{		
			#folder for final vcf
			my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec";
			print_and_log( " The genomics db folder out is $genomics_db_f ..\n ",$log_file);
			#target
			my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." "; 
			my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $input_file_str --genomicsdb-workspace-path $genomics_db_f  $l_string $param_str";
			print_and_log( " The command is $command ..\n ",$log_file);
			#print_and_log("Preparing command: $command\n",$log_file);
			push(@commands,$command);
		}
	}
	#print_and_log("END",$log_file);
	#Run a different job for each command
	my $jobs_to_wait_dbi = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$task,$task,$log_file);
	#Return the job ids for each job and the list of output file and the final output name
	return $jobs_to_wait_dbi;
}


=head2 run_GATK4_GenotypeGVCFs

 Title   : run_GATK4_GenotypeGVCFs
 Usage   : run_GATK4_GenotypeGVCFs(   );

 Function: run_GATK4_GenotypeGVCFs
			Given the gvcf files from all the samples in a groups will make a Joint
			Genotype using GATK tools
			
 Returns : nothing
=cut
sub run_GATK4_GenotypeGVCFs{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;
  my $step = shift;
  my $step_needed = shift;
  my $steps_array2 = shift;
  my $steps_array3 = shift;

        
	my $previous_task = 'varcall';
	
	#Getting the sample name
	my $sample_name = $cfg_hash->{'db_sample_name'}; 
	#Getting the group name
	my $group_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'}
						,$cfg_hash->{'db_analysis_id'},$group_id);	
																																			
	#Here we get the samples names associated with that group id
	my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							$cfg_hash->{'db_analysis_id'},$group_id);
	my $params;

	#input folder setting
	my $inFolder = $cfg_hash->{$group_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
	
	
	#Prepare parameters string
	my $param_str = " ";	
		
	#If there are more samples associated to this group
	#Constructing input string with files to be merged 
	my $input_file_str = "";
	
	#If GenotypeGVCFs is executed for a single sample we reduce the threshold for confidence call
	#(https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller)
	#Set the phred-scaled Qscore threshold to separate high confidence from low confidence calls
	if (defined $cfg_hash->{'ggvcf_single_sam_stand_call_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_call_conf'},"positiveint") ) {
		$param_str .= " -stand-call-conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
	}
	##Set the java memory
	#if (defined $cfg_hash->{'ggvcf_single_sam_stand_emit_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_emit_conf'},"positiveint") ){
		#$param_str .= " -stand_emit_conf ".$cfg_hash->{'ggvcf_single_sam_stand_emit_conf'}." ";
	#}	
	
	if ( (scalar keys %{$distinct_samples}) > 1){
		#Check from the DB if the previous step has been performed for all the samples
		foreach my $dist_sample ( keys %{$distinct_samples}){
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$dist_sample);
			$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
		}
	}else{

		
		#Get infos for the single sample associated to this group from the 'samples' table
		#print_and_log("Getting information for $group_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
							$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash-> {'db_analysis_id'},$group_id);	
		$input_file_str = " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	}
	
	#Set the unique output name
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					$group_name,$steps_array3)."_".$step.".".$cfg_hash->{'gvcf_ext'};

	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_GGVCF'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_GGVCF'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	#if ( defined $cfg_hash->{'ggvcf_dbsnp'} ){
	#	$param_str .= " -dbsnp ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'ggvcf_dbsnp'}};
	#}		
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'ggvcf_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'ggvcf_annotation'});
		foreach my $annotation (@annotations){
			$param_str .= " -A $annotation ";
		}
	}

	#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
	my $command = $jatk_call." $input_file_str --reference ".$cfg_hash->{'hum_ref'}." $param_str --output $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

#=head2 run_GATK4_GenotypeGVCFs_parallel

 #Title   : run_GATK4_GenotypeGVCFs_parallel
 #Usage   : run_GATK4_GenotypeGVCFs_parallel(   );

 #Function: run_GATK4_GenotypeGVCFs_parallel
			#Given the gvcf files from all the samples in a groups will make a Joint
			#Genotype using GATK tools
			
 #Returns : Return the job ids for each job and the list of output file and the final output name
#=cut
#sub run_GATK4_GenotypeGVCFs_parallel{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $analysis_id = shift;
	#my $log_file = shift;
	#my $task = shift;
  #my $step = shift;
  #my $step_needed = shift;
  #my $steps_array2 = shift;
  #my $steps_array3 = shift;

        
	#my $previous_task = 'varcall';
	
	##Getting the sample name
	#my $sample_name = $cfg_hash->{'db_sample_name'}; 
	##Getting the group name
	#my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						#$cfg_hash->{'db_analysis_id'},$analysis_id);	
																																			
	##Here we get the samples names associated with that group id
	#my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							#$cfg_hash->{'db_analysis_id'},$analysis_id);
	#my $params;

	##input folder setting
	#my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	#my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	
	
	##Prepare parameters string
	#my $param_str = " ";	
		
	##If there are more samples associated to this group
	##Constructing input string with files to be merged 
	#my $input_file_str = "";
	
	##If GenotypeGVCFs is executed for a single sample we reduce the threshold for confidence call
	##(https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller)
	##Set the phred-scaled Qscore threshold to separate high confidence from low confidence calls
	#if (defined $cfg_hash->{'ggvcf_single_sam_stand_call_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_call_conf'},"positiveint") ) {
		#$param_str .= " -stand-call-conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
	#}
	###Set the java memory
	##if (defined $cfg_hash->{'ggvcf_single_sam_stand_emit_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_emit_conf'},"positiveint") ){
		##$param_str .= " -stand_emit_conf ".$cfg_hash->{'ggvcf_single_sam_stand_emit_conf'}." ";
	##}	
	
	#if ( (scalar keys %{$distinct_samples}) > 1){
		##Check from the DB if the previous step has been performed for all the samples
		#foreach my $dist_sample ( keys %{$distinct_samples}){
			#getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$dist_sample);
			#$input_file_str .= " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								#$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
		#}
	#}else{

		
		##Get infos for the single sample associated to this group from the 'samples' table
		##print_and_log("Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		#getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
							#$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash-> {'db_analysis_id'},$analysis_id);	
		#$input_file_str = " --variant ".$inFolder."/".build_input_name_from_executed($params,$step,
								#$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	#}
	


	##print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);

	
	##Set the Java call
	#my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
	######## JAVA PARAMETERS	
	#if (defined $cfg_hash->{'javamem_GGVCF'} or defined $cfg_hash->{'scratch_f'}){
		##Open string
		#$jatk_call .= " --java-options '";
			##Add temporary folder
			#if ( defined $cfg_hash->{'scratch_f'} ){
				#$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			#}	
			#if (defined $cfg_hash->{'javamem_GGVCF'}){
				#$jatk_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
			#}			
			##Close string
			#$jatk_call .= "'";
		#}
	########END JAVA PARAMETERS

	##if ( defined $cfg_hash->{'ggvcf_dbsnp'} ){
	##	$param_str .= " -dbsnp ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'ggvcf_dbsnp'}};
	##}		
	##Adds known indel snp to the list of known
	#if ( defined $cfg_hash->{'ggvcf_annotation'} ){
		#my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'ggvcf_annotation'});
		#foreach my $annotation (@annotations){
			#$param_str .= " -A $annotation ";
		#}
	#}


	###Using also Chromosomex X and Y
	#my $numExec = $cfg_hash->{'chromosomes_num'};
	#my @chroms = 1..$numExec;
	#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	##An array will contain all the commands to be executed
	#my @commands = ();

	##Set the output name
	#my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					#$analysis_name,$steps_array3)."_".$step.".".$cfg_hash->{'gvcf_ext'};
					
	##Path to the script to execute multiple jobs using an array of commands
	#my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	#my $out_files_str = "";
	##Executing one GATK command for each target chromosome using bed files
	#foreach my $exec (@chroms){
			#my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
			
			##Set the output name
			#my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
					#$analysis_name,$steps_array3)."_".$step."_$exec.".$cfg_hash->{'gvcf_ext'};
			#$out_files_str .= $nth_out.",";					
			#my $command = $jatk_call." $input_file_str --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $nth_out ";			 
			##print_and_log("Preparing command: $command\n",$log_file);
			#push(@commands,$command);
	#}
	##Run a different job for each command
	#my $jobs_to_wait = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$task,$log_file);
	
	##Return the job ids for each job and the list of output file and the final output name
	#chop($out_files_str);
	#return $jobs_to_wait,$out_files_str,$outFile;
#}

=head2 run_GATK4_GenotypeGVCFs_parallel

 Title   : run_GATK4_GenotypeGVCFs_parallel
 Usage   : run_GATK4_GenotypeGVCFs_parallel(   );

 Function: run_GATK4_GenotypeGVCFs_parallel
 
		   Given the gvcf file  it will make a Joint Genotype using GATK tools		
		   It executes several separate jobs in parallel using execute_jobs
		   
		   If GET_VCF is used as last parameter, it does not execute the jobs
		   but only returns the paths of the VCFs. Useful when you already have
		   ran the jobs and only want to gather them
		   
Input :  GenomicsDB workspace created by GenomicsDBImport	
		
Output : A final VCF in which all samples have been jointly genotyped

Returns : Return the job ids for each job and the list of output file and the final output name
=cut
sub run_GATK4_GenotypeGVCFs_parallel{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $sample_id = shift;
	my $task = shift;
	my $step = shift;
	my $step_needed = shift;
	my $steps_array2 = shift;
	my $steps_array3 = shift;
	my $jobs_to_wait_dbi = shift;
	my $log_file = shift;
	my $get_vcf = shift; 
	
    
    #print  "$log_file";
	my $previous_task = 'varcall';
	
	#Getting the sample name
	my $sample_name = $cfg_hash->{'db_sample_name'}; 
	#Getting the analysis name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);	
	#Verifying if there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
					    $cfg_hash->{'db_analysis_id'},$analysis_id);
					    																																			
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
			$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_analyses_table'},$cfg_hash-> {'db_analysis_id'},$analysis_id);	

	#input and output folder setting
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};  
	my $param_str = " ";	
		
	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);	
	#Get target bed
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}		
	}else{$target_bed = "-";}
		
	#If GenotypeGVCFs is executed for a single sample we reduce the threshold for confidence call
	#(https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller)
	#Set the phred-scaled Qscore threshold to separate high confidence from low confidence calls
	if (defined $cfg_hash->{'ggvcf_single_sam_stand_call_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_call_conf'},"positiveint") ) {
	     $param_str .= " -stand_call_conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
	}
	
	print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_GGVCF'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_GGVCF'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'ggvcf_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'ggvcf_annotation'});
		foreach my $annotation (@annotations){
			$param_str .= " -A $annotation ";
		}
	}

	#NB: In GATK4 calculation of Mapping Quality has been updated
	if ( $cfg_hash->{'ggvcf_allow_old_rms'} ne 'NO'){
		$param_str .= " --allow-old-rms-mapping-quality-annotation-data ";
	}
	
	
 	#Set the output name
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					$analysis_name,$steps_array3)."_".$step.".".$cfg_hash->{'vcf_ext'};
					   
	##Using also Chromosomex X and Y
	my $numExec = $cfg_hash->{'chromosomes_num'};
	my @chroms = 1..$numExec;
	push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	#An array will contain all the commands to be executed
	my @commands = ();
    my @commands2 = ();

					
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	my $out_files_str = "";
	
	#Check if some chromosomes is so huge to use splitted versions for the jobs them in two different parts
	my @chroms_to_split = ();
	if ($cfg_hash->{'chroms_to_split'} ne 'none' ){
		my @chroms_to_split = split(",",$cfg_hash->{'chroms_to_split'});
	}
	
	#Parallel execution
	if ( $perchrom == 1 ){
		#Executing one GATK command for chr1 target chromosome using bed file
		foreach my $exec (@chroms){
			#Check if this chromosome has been splitted to separate job
			#if ( grep {/\b$exec\b/} @chroms_to_split){
			if ( $exec == 1 ){
				my @chr_parts = (1,2);
				foreach my $chr_part (@chr_parts){
					my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec\_$chr_part";
					my $genomics_db_path = "-V gendb://".$genomics_db_f;
					#folder created
					print_and_log( " The genomics db folder out for chr$exec part $chr_part is $genomics_db_f ..\n ",$log_file);
					#targetbed  for both parts chr1
					my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec\_$chr_part.".$cfg_hash->{'bed_ext'}." ";			 
					
					#Set the output name (job)
					my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
							$analysis_name,$steps_array3)."_".$step."_$exec\_$chr_part.".$cfg_hash->{'gvcf_ext'};		
							
					$out_files_str .= $nth_out.",";
					my $command = $jatk_call." $genomics_db_path --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $nth_out ";
					print_and_log("Preparing command: $command\n",$log_file);
					push(@commands,$command);	
				}
			}else{		
				my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
				my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec";
				#set input file/folder from which the program takes the merged vcfs
				my $genomics_db_path = "-V gendb://".$genomics_db_f;
				#Set the output name of the final vcf for each chr
				my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
						$analysis_name,$steps_array3)."_".$step."_$exec.".$cfg_hash->{'gvcf_ext'};
				$out_files_str .= $nth_out.",";					
				my $command = $jatk_call." $genomics_db_path --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $nth_out ";			 
				print_and_log("Preparing command: $command\n",$log_file);
				push(@commands,$command);
			} 
		}
		#print_and_log("END",$log_file);
		my $jobs_to_wait = "no_depend";
		if ($get_vcf ne "GET_VCF_LIST"){
			#Run a different job for each command
			$jobs_to_wait = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,$jobs_to_wait_dbi,$task,$task,$log_file);
		 }
		 
		chop($out_files_str);
		#Return the job ids for each job and the list of output file and the final output name
		return $jobs_to_wait, $out_files_str, $outFile;
	}else{
		my $l_string = " --intervals $target_bed "; 
		my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'};
		#set input file/folder from which the program takes the merged vcfs
		my $genomics_db_path = "-V gendb://$genomics_db_f ";
		#Set the output name of the final vcf for each chr
		my $out = $outFolder."/".build_input_name_from_executed($params,$step,
				$analysis_name,$steps_array3).".".$cfg_hash->{'gvcf_ext'};
		my $command = $jatk_call." $genomics_db_path --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $out ";			 
		print_and_log("Preparing command: $command\n",$log_file);
	}

}


=head2 run_GATK4_GenotypeGVCFs_parallel

 Title   : run_GATK4_GenotypeGVCFs_parallel
 Usage   : run_GATK4_GenotypeGVCFs_parallel(   );

 Function: run_GATK4_GenotypeGVCFs_parallel
 
		   Given the gvcf file  it will make a Joint Genotype using GATK tools		
		   It executes several separate jobs in parallel using execute_jobs
		   
		   If GET_VCF is used as last parameter, it does not execute the jobs
		   but only returns the paths of the VCFs. Useful when you already have
		   ran the jobs and only want to gather them
		   
Input :  GenomicsDB workspace created by GenomicsDBImport	
		
Output : A final VCF in which all samples have been jointly genotyped

Returns : Return the job ids for each job and the list of output file and the final output name
=cut
sub run_GATK4_GenotypeGVCFs_parallelOLD{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $sample_id = shift;
	my $task = shift;
	my $step = shift;
	my $step_needed = shift;
	my $steps_array2 = shift;
	my $steps_array3 = shift;
	my $jobs_to_wait_dbi = shift;
	my $log_file = shift;
	my $get_vcf = shift; 
	
    
    #print  "$log_file";
	my $previous_task = 'varcall';
	
	#Getting the sample name
	my $sample_name = $cfg_hash->{'db_sample_name'}; 
	#Getting the analysis name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);	
																																			
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
			$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_analyses_table'},$cfg_hash-> {'db_analysis_id'},$analysis_id);	

	#input and output folder setting
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};  
	my $param_str = " ";	
		
	
	#If GenotypeGVCFs is executed for a single sample we reduce the threshold for confidence call
	#(https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller)
	#Set the phred-scaled Qscore threshold to separate high confidence from low confidence calls
	if (defined $cfg_hash->{'ggvcf_single_sam_stand_call_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_call_conf'},"positiveint") ) {
	     $param_str .= " -stand_call_conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
	}
	
	print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_GGVCF'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_GGVCF'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'ggvcf_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'ggvcf_annotation'});
		foreach my $annotation (@annotations){
			$param_str .= " -A $annotation ";
		}
	}

	#NB: In GATK4 calculation of Mapping Quality has been updated
	if ( $cfg_hash->{'ggvcf_allow_old_rms'} ne 'NO'){
		$param_str .= " --allow-old-rms-mapping-quality-annotation-data ";
	}
	
	
 	#Set the output name
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					$analysis_name,$steps_array3)."_".$step.".".$cfg_hash->{'vcf_ext'};
					   
	##Using also Chromosomex X and Y
	my $numExec = $cfg_hash->{'chromosomes_num'};
	my @chroms = 1..$numExec;
	push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	#An array will contain all the commands to be executed
	my @commands = ();
    my @commands2 = ();

					
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	my $out_files_str = "";
	
	#Check if some chromosomes is so huge to use splitted versions for the jobs them in two different parts
	my @chroms_to_split = ();
	if ($cfg_hash->{'chroms_to_split'} ne 'none' ){
		my @chroms_to_split = split(",",$cfg_hash->{'chroms_to_split'});
	}
		
	#Executing one GATK command for chr1 target chromosome using bed file
	foreach my $exec (@chroms){
		#Check if this chromosome has been splitted to separate job
		#if ( grep {/\b$exec\b/} @chroms_to_split){
		if ( $exec == 1 ){
			my @chr_parts = (1,2);
			foreach my $chr_part (@chr_parts){
				my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec\_$chr_part";
				my $genomics_db_path = "-V gendb://".$genomics_db_f;
				#folder created
				print_and_log( " The genomics db folder out for chr$exec part $chr_part is $genomics_db_f ..\n ",$log_file);
				#targetbed  for both parts chr1
				my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec\_$chr_part.".$cfg_hash->{'bed_ext'}." ";			 
				
				#Set the output name (job)
				my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
						$analysis_name,$steps_array3)."_".$step."_$exec\_$chr_part.".$cfg_hash->{'gvcf_ext'};		
						
				$out_files_str .= $nth_out.",";
				my $command = $jatk_call." $genomics_db_path --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $nth_out ";
				print_and_log("Preparing command: $command\n",$log_file);
				push(@commands,$command);	
			}
		}else{		
			my $l_string = " --intervals ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
			my $genomics_db_f = $outFolder."/".$cfg_hash->{'gatk_genomicsdb_name'}."_chr$exec";
			#set input file/folder from which the program takes the merged vcfs
			my $genomics_db_path = "-V gendb://".$genomics_db_f;
			#Set the output name of the final vcf for each chr
			my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
					$analysis_name,$steps_array3)."_".$step."_$exec.".$cfg_hash->{'gvcf_ext'};
			$out_files_str .= $nth_out.",";					
			my $command = $jatk_call." $genomics_db_path --reference ".$cfg_hash->{'hum_ref'}." $param_str $l_string --output $nth_out ";			 
			print_and_log("Preparing command: $command\n",$log_file);
			push(@commands,$command);
		} 
	}
	#print_and_log("END",$log_file);
	my $jobs_to_wait = "no_depend";
	if ($get_vcf ne "GET_VCF_LIST"){
		#Run a different job for each command
		$jobs_to_wait = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,$jobs_to_wait_dbi,$task,$task,$log_file);
     }
     
	chop($out_files_str);
	#Return the job ids for each job and the list of output file and the final output name
	return $jobs_to_wait, $out_files_str, $outFile;
}


=head2 run_GATK4_VariantFiltration

 Title   : run_GATK4_VariantFiltration
 Usage   : run_GATK4_VariantFiltration(   );

 Function: run_GATK4_VariantFiltration
			Given the gvcf files for the single group it will filter base on 
			parameters given in input.
			The group name here is in input since we do not need sample informations
			
 Returns : nothing 
=cut
sub run_GATK4_VariantFiltration{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $param_str = shift;	
	my $filters = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_VF'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_VF'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_VF'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Get the filter set from input to filter the INFO field
	my @filter_set = split(",",$filters);
	foreach my $filter ( @filter_set){
		if ( defined $cfg_hash->{'filterExpression'.$filter} ){
			$param_str .= " --filter-expression ".$cfg_hash->{'filterExpression'.$filter}." ";
		}	
		if ( defined $cfg_hash->{'filterName'.$filter} ){
			$param_str .= " --filter-name ".$cfg_hash->{'filterName'.$filter}." ";
		}		
		#Genotype filters variables have _gX
		if ( defined $cfg_hash->{'filterExpression_'.$filter} ){
			$param_str .= " --genotype-filter-expression ".$cfg_hash->{'filterExpression_'.$filter}." ";
		}	
		if ( defined $cfg_hash->{'filterName_'.$filter} ){
			$param_str .= " --genotype-filter-name ".$cfg_hash->{'filterName_'.$filter}." ";
		}				
	}

	my $command = $jatk_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str --output $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK4_SelectVariants

 Title   : run_GATK4_SelectVariants
 Usage   : run_GATK4_SelectVariants(   );

 Function: run_GATK4_SelectVariants
			Subsets a VCF in order to facilitate analyses
			
			It is used here to select INDELS and SNP to two different files
			but it has many purposes
			
 Returns : nothing 
=cut
sub run_GATK4_SelectVariants{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $param_str = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_SV'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_SV'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_SV'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	my $command = $jatk_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str --output $outFile ";			 
	print_and_log("[gr:$analysis_id] Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

=head2 run_GATK4_CalculateGenotypePosteriors

 Title   : run_GATK4_CalculateGenotypePosteriors
 Usage   : run_GATK4_CalculateGenotypePosteriors(   );

 Function: run_GATK4_CalculateGenotypePosteriors
					This function derives the posteriors of genotype calls in our callset, 
					recalibratedVariants.vcf, which just came out of the VQSR filtering step;
					 it contains among other samples a trio of individuals (mother, father and child) 
					 whose family structure is described in the pedigree file trio.ped (which you need to supply).
					 To do this, we are using the most comprehensive set of high confidence SNPs available to us,
					 a set of sites from Phase 3 of the 1000 Genomes project (available in our resource bundle),
					 which we pass via the --supporting argument.
			
 Returns : nothing 
=cut
sub run_GATK4_CalculateGenotypePosteriors{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $ped_file = shift;
	my $outFile = shift;
	my $param_str = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);##DEBUGCODE
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_CGP'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_CGP'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_CGP'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	#Execute command
	my $command = $jatk_call." --variant $input_file --reference ".$cfg_hash->{'hum_ref'}." $param_str --pedigree $ped_file --output $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK4_VariantAnnotator

 Title   : run_GATK4_VariantAnnotator
 Usage   : run_GATK4_VariantAnnotator(   );

 Function: run_GATK4_VariantAnnotator
					Executes VariantAnnotator
			
 Returns : nothing 
=cut
sub run_GATK4_VariantAnnotator{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;	
	my $ped_file = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_VA'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_VA'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_VA'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Set the filters by searching in the user file
	my $param_str = " ";	

	#Adds annotations to use
	if ( defined $cfg_hash->{'GRW_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'GRW_annotation'});
		foreach my $annotation (@annotations){
					$param_str .= " --annotation $annotation ";
		}
	}
		
	my $command = $jatk_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str --pedigree $ped_file --output $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}





=head2 run_GATK4_CombineVariants

 Title   : run_GATK4_CombineVariants
 Usage   : run_GATK4_CombineVariants(   );

 Function: run_GATK4_CombineVariants
						combines multiple variant records present at the same site in the different 
						input sources into a single variant record in the output. If sample names overlap
						, then they are "uniquified" by default, which means a suffix is appended to make them unique
			
 Returns : nothing 
=cut
sub run_GATK4_CombineVariants{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_vcfs = shift;
	my $outFile = shift;
	my $log_file = shift;
  my $get_command  =shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_CV'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_CV'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_CV'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Set the filters by searching in the user file
	my $param_str = " ";	
	
	my @vcf_set = split(",",$input_vcfs);
	
	foreach my $vcf ( @vcf_set ){
		$param_str .= " --variant $vcf ";		
	}

	if (defined $cfg_hash->{'cv_genotypeMergeOptions'}){
		$param_str .= " -genotypeMergeOptions ".$cfg_hash->{'cv_genotypeMergeOptions'}." ";
	}
	
	my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." $param_str --output $outFile ";			 
	
	if ($get_command ne 'YES'){
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";	
  }
	return $command;
}

=head2 run_GATK4_VariantRecalibration

 Title   : run_GATK4_VariantRecalibration
 Usage   : run_GATK4_VariantRecalibration(   );

 Function: run_GATK4_VariantRecalibration
			Is a tool from GATK which permits to recalibrate variants. It can work
			per-sample 
			
 Returns : nothing
=cut
sub run_GATK4_VariantRecalibration{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $mode = shift;
	my $input_file = shift;
	my $recal_file = shift;
	my $tranches_file = shift;
	my $rscript_file = shift;
	my $log_file = shift;

	
	#print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Prepare output name
	my $param_str = " ";	
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_VR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_VR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_VR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	if ( defined $cfg_hash->{'vr_std'} ){
					$param_str .= " --stdThreshold ".$cfg_hash->{'vr_std'};
	}
	#Max Gaussians is the number of different cluster to made with variants 
	if ( defined $cfg_hash->{"vr_max_gaussians_$mode"} ){
					$param_str .= " --max-gaussians ".$cfg_hash->{"vr_max_gaussians_$mode"};
	}	
	my $numres = $cfg_hash->{"vr_numres_$mode"};

	#Get resources strings
	my $res_string = "";
	for (my $r = 1; $r <= $numres; $r++){
		$res_string .= "-resource:".$cfg_hash->{"vr_resources_$mode".$r}." ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{"vr_resources_$mode".$r.'_file'}}." ";
	}

	#Adds known indel snp to the list of known
	my $ann_string = "";
	if ( defined $cfg_hash->{'vr_'.$mode.'_annotation'} ){
					my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'vr_'.$mode.'_annotation'});
					foreach my $annotation (@annotations){
									$ann_string .= " -an $annotation ";
					}
	}


	my $command = $jatk_call." --reference ".$cfg_hash->{'hum_ref'}." --variant $input_file $res_string $ann_string -mode $mode --output $recal_file --tranches-file $tranches_file ".
								" --rscript-file $rscript_file  $param_str";			 
	print_and_log("Executing command: $command\n",$log_file);
	
	if (file_not_present($input_file) > 0 ){ warn "Could not proceed with $prog_used! Check: $input_file.\n";}
	#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_GATK4_ApplyVQSR

 Title   : run_GATK4_ApplyVQSR
 Usage   : run_GATK4_ApplyVQSR(   );

 Function: run_GATK4_ApplyVQSR
			Is a tool from GATK which permits to recalibrate variants.
			
 Returns : nothing
=cut
sub run_GATK4_ApplyVQSR{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $mode = shift;
	my $input_file = shift;
	my $recal_file = shift;
	my $tranches_file = shift;
	my $rscript_file = shift;
	my $out_vcf = shift;
	my $log_file = shift;
			
	#print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Prepare output name
	my $param_str = " ";	
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
	
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_AR'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_AR'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_AR'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	if ( defined $cfg_hash->{'ar_ts_filter_level_'.$mode} ){
		$param_str .= " --truth-sensitivity-filter-level ".$cfg_hash->{'ar_ts_filter_level_'.$mode};
	}

	#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
	my $command = $jatk_call." --variant $input_file --reference ".$cfg_hash->{'hum_ref'}." -mode $mode --recal-file $recal_file --tranches-file $tranches_file ".
								"$param_str  --output $out_vcf ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 run_GATK4_PhaseByTransmission

 Title   : run_GATK4_PhaseByTransmission
 Usage   : run_GATK4_PhaseByTransmission(   );

 Function: run_GATK4_PhaseByTransmission
			Given the filtered gvcf files and parent/child information for the single group it will make the analysis
			of phasing
			
 Returns : nothing
=cut
sub run_GATK4_PhaseByTransmission{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $group_name = shift;
	my $group_id = shift;
	my $outFolder = shift;
	my $input_file = shift;
  my $out_vcf = shift;
  my $ped_f_path = shift;
	my $log_file = shift;
	
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";

####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_PBT'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_PBT'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_PBT'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS
	
	#Set the filters by searching in the user file
	my $param_str = " ";	

	#Cluster Window size
	if ( $cfg_hash->{'FatherAlleleFirst'} eq 'YES'){
		$param_str .= " --FatherAlleleFirst ";
	}
	if (defined $cfg_hash->{'pedigreeValidationType'}){
		$param_str .= " --pedigreeValidationType ".$cfg_hash->{'pedigreeValidationType'}." ";
	}
	if (defined $cfg_hash->{'MendelianViolationsFile'}){
		$param_str .= " --MendelianViolationsFile ".$outFolder."/".extract_name($out_vcf,1).".".$cfg_hash->{'MendelianViolationsFile'}." ";
	}
		
	#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
	my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." --variant $input_file -ped $ped_f_path $param_str -o $out_vcf ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or log_and_exit( "Unable to execute command: $command\n",$log_file);	
}


=head2 run_GATK4_DepthOfCoverage

 Title   : run_GATK4_DepthOfCoverage
 Usage   : run_GATK4_DepthOfCoverage(   );

 Function: run_GATK4_DepthOfCoverage
						Given a BAM file of aligned reads and a file with target regions
						it uses DepthOfCoverage to generate statistics 
			
 Returns : nothing
=cut
sub run_GATK4_DepthOfCoverage{
	my $cfg_hash = shift;
	my $bam_input = shift;#This can also be a file with paths
	my $out_suffix = shift;
	my $group_id = shift;
	my $log_file = shift;
	  
  #Program name
  my $prog_used = $cfg_hash->{'depth_coverage_prog'};
	#Set the Java call
	my $jatk_call = $cfg_hash->{'gatk_path'}." $prog_used ";
	
	
	####### JAVA PARAMETERS	
	if (defined $cfg_hash->{'javamem_DOC'} or defined $cfg_hash->{'scratch_f'}){
		#Open string
		$jatk_call .= " --java-options '";
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$jatk_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}	
			if (defined $cfg_hash->{'javamem_DOC'}){
				$jatk_call .= " -".$cfg_hash->{'javamem_DOC'}." ";
			}			
			#Close string
			$jatk_call .= "'";
		}
	#######END JAVA PARAMETERS

	
	#Set the filters by searching in the user file
	my $param_str = " ";	

	if ( $cfg_hash->{'filter_mismatching_base_and_quals_DOC'} eq 'YES'){
		$param_str .= " --filter_mismatching_base_and_quals ";
	}	
			
	if (defined $cfg_hash->{'DOC_genelist'}){
		my $gene_list_file = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'DOC_genelist'};
		die "ERROR: Cannot proceed with $prog_used because cannot find $gene_list_file\n" unless (-e $gene_list_file);
		$param_str .= " -geneList ".$cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'DOC_genelist'}." ";
	}
	#Put the coverage thresholds if present
	if (defined $cfg_hash->{'DOC_coverage_thresholds'}){
		my @trhs = split(",",$cfg_hash->{'DOC_coverage_thresholds'});
		foreach my $thr (@trhs){
				$param_str .= " -ct $thr ";
		}
	}
	
	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = $cfg_hash->{'target_reg_f'}."/".
						get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						$cfg_hash->{'db_analysis_id'},$group_id);
	if (file_not_present($target_bed)  > 0){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
	
	#Construct the target interval list file from the target file
	my $target_int_list = $target_bed.".".$cfg_hash->{'target_inter_list_ext'};
	
	#print_and_log( "$jatk_call - gatk params: $picard_params\n",$log_file);
	my $command = $jatk_call." -R ".$cfg_hash->{'hum_ref'}." -I $bam_input -L $target_int_list $param_str -o $out_suffix ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or log_and_exit( "Unable to execute command: $command\n",$log_file);	
}

1;
