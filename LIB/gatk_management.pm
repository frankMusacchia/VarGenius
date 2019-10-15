
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
    
package LIB::gatk_management;
## gatk_management.pm
#Permits the management of the gatk tools for the pipeline

BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_GATK_RealignerTargetCreator run_GATK_IndelRealigner 
            run_GATK_BaseRecalibrator_perlane run_GATK_HaplotypeCaller
            run_GATK_CatVariants run_GATK_GenotypeGVCFs run_GATK_VariantFiltration
            run_GATK_VariantRecalibration run_GATK_ApplyRecalibration 
            run_GATK_PhaseByTransmission run_GATK_DepthOfCoverage
            run_GATK_BaseRecalibrator_persample run_GATK_PrintReads_persample
            run_GATK_SelectVariants run_GATK_CombineVariants
            run_GATK_VariantAnnotator run_GATK_CalculateGenotypePosteriors 
            run_GATK_SplitNCigarReads run_GATK_GenotypeGVCFs_parallel
            run_GATK_GCContentByInterval run_GATK_DepthOfCoverage_gen
            run_GATK_HaplotypeCaller_parallel);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit  execute_threads build_input_name_from_executed
				execute_jobs JOB_get_jobid_form_jobname alter_job update_dependencies_hash);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present 
																delete_file extract_name);
																
#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);

##########################################################GATK

=head2 run_GATK_SplitNCigarReads

 Title   : run_GATK_SplitNCigarReads
 Usage   : run_GATK_SplitNCigarReads(   );

 Function: run_GATK_SplitNCigarReads
			Is a tool from GATK which splits reads into exon segments 
			(getting rid of Ns but maintaining grouping information) and 
			hard-clip any sequences overhanging into the intronic regions.
 Returns : nothing
=cut
sub run_GATK_SplitNCigarReads{
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
		my $java_call = $cfg_hash->{'java_path'};
		#Set the java memory
		if (defined $cfg_hash->{'javamem_SNCR'}){
			$java_call .= " -".$cfg_hash->{'javamem_SNCR'}." ";
		}
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}
		#Set the program to use
		$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";


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
		my $command = $java_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bam ";			 
		print_and_log("Executing command: $command \n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";

	#}else{print_and_log( "ERROR: Step $step_needed has not been run for sample: $sample_id.f\n ",$log_file);}
	
}


=head2 run_GATK_RealignerTargetCreator

 Title   : run_GATK_RealignerTargetCreator
 Usage   : run_GATK_RealignerTargetCreator(   );

 Function: run_GATK_RealignerTargetCreator
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			This function takes in input the bam file from the alignment and produces a file 
			to be used for the indel realignment later. Hence in this step 
 Returns : nothing
=cut
sub run_GATK_RealignerTargetCreator{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
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
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);
						    			
		#Set the input file name
		my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array).".".$cfg_hash->{'bam_ext'};
				
		print_and_log( "Starting RealignerTargetCreator for $input_bam using GATK tools..\n ",$log_file);
		
		#Prepare output name
		my $targ_list = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array)."_".$step.".".$cfg_hash->{'targ_int_ext'};
		my $param_str = " ";	
		
		#Set the Java call
		my $java_call = $cfg_hash->{'java_path'};
		#Set the java memory
		if (defined $cfg_hash->{'javamem_RTC'}){
			$java_call .= " -".$cfg_hash->{'javamem_RTC'}." ";
		}
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}			
		#Set the program to use
		$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

		if ( defined $cfg_hash->{'nt_RTC'} ){
			$param_str .= " -nt ".$cfg_hash->{'nt_RTC'};
		}
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
			##Using also Chromosomex X and Y
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
					my $command = $java_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_targ ";			 
					print_and_log("Executing command: $command\n",$log_file);
					push(@commands,$command);
			}
			#Run a different thread for each command
			execute_threads(\@commands);
		}else{
			# When perchrom is not used the HaplotypeCaller uses the target regions as intervals	
			#Getting the target bed file using the folder for targets and the name contained in the database
			my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);
			if ($target_bed ne '-'){
				$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
			}else{$target_bed = "-";}
												
			if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}
			if ( $target_bed ne '-'){ 
				$param_str .= " --intervals $target_bed ";	
			}			
		
			if (file_not_present($input_bam) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_bam.\n";}
			my $command = $java_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $targ_list ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}
	}else{print_and_log( "ERROR: Step sort_idx_step has not been run for  sample: $sample_id.f\n ",$log_file);}
	
}


=head2 run_GATK_IndelRealigner

 Title   : run_GATK_IndelRealigner
 Usage   : run_GATK_IndelRealigner(   );

 Function: run_GATK_IndelRealigner
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			
 Returns : nothing
=cut
sub run_GATK_IndelRealigner{
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
			my $java_call = $cfg_hash->{'java_path'};
			#Set the java memory
			if (defined $cfg_hash->{'javamem_IR'}){
							$java_call .= " -".$cfg_hash->{'javamem_IR'}." ";
			}
			#Add temporary folder
			if ( defined $cfg_hash->{'scratch_f'} ){
				$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
			}			
			#Set the program to use
			$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

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
				##Using also Chromosomex X and Y
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
					#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					my $command = $java_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string $targ_interv -o $ith_out_bam ";			 
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
					#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					my $command = $java_call." -I $input_bam -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bam ";			 
					print_and_log("Executing command: $command\n",$log_file);
					try_exec_command($command) or die "Unable to execute command: $command\n";
			}
    }else{print_and_log( "ERROR: Step $step has not been run for sample: $sample_id.f\n ",$log_file);}
    
}

=head2 run_GATK_BaseRecalibrator_persample

 Title   : run_GATK_BaseRecalibrator_persample
 Usage   : run_GATK_BaseRecalibrator_persample(   );

 Function: run_GATK_BaseRecalibrator_persample
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			
 Returns : nothing
=cut
sub run_GATK_BaseRecalibrator_persample{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
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
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			
		
		#Prepare output name for files for data recalibration
		my $out_bqsr = $outFolder."/".build_input_name_from_executed($params,$curr_step,
						$params->{$db_sample_name},$steps_array2)."_".$curr_step.".".$cfg_hash->{'recal_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
												    		
		#print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}."  using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $java_call = $cfg_hash->{'java_path'};
		#Set the java memory
		if (defined $cfg_hash->{'javamem_BR'}){
			$java_call .= " -".$cfg_hash->{'javamem_BR'}." ";
		}
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}		
		#Set the program to use
		$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

		#Parameters string construction
		my $param_str = " ";
		#Number of threads to use CANNOT USE IN INDELREALIGNER
		if ( defined $cfg_hash->{'refine_threads'} ){
			$param_str .= " -nct ".$cfg_hash->{'refine_threads'}." ";
		}
		if ( $cfg_hash->{'fix_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --fix_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'allow_potentially_misencoded_quality_scores'} eq 'YES'){
			$param_str .= " --allow_potentially_misencoded_quality_scores ";
		}
		if ( $cfg_hash->{'filter_mismatching_base_and_quals_BR'} eq 'YES'){
			$param_str .= " --filter_mismatching_base_and_quals ";
		}		
		#Adds known indel snp to the list of known
		if ( defined $cfg_hash->{'known_BR'} ){
			my @known_indel_snp  = split($cfg_hash->{'word_sep'},$cfg_hash->{'known_BR'});
			foreach my $known (@known_indel_snp){
				#$param_str .= " -knownSites:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
				$param_str .= " -knownSites ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$known}." ";
			}
		}
		#Adds covariates to use
		if ( defined $cfg_hash->{'covariates_BR'} ){
			print "Covariates: ".$cfg_hash->{'covariates_BR'}."\n";
			my @covariates  = split(",",$cfg_hash->{'covariates_BR'});
			foreach my $covariate (@covariates){
				$param_str .= " -cov $covariate ";
			}
		}		

		#if ( defined $cfg_hash->{'BR_plots'} ){
			#$param_str .= " -plots ".$cfg_hash->{'BR_plots'}." ";
		#}

		#Getting the target bed file using the folder for targets and the name contained in the database
		my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
							,$cfg_hash->{'db_analysis_id'},$analysis_id);
		if ($target_bed ne '-'){
			$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		}else{$target_bed = "-";}
		
									
		my $numExec = 1;
		
		#OBSOLETE CODE FROM GATK3.8
		#After the update of GATK that the IndelRealignment should not be done the
		#per-chromosome mode is deactivated here because starting from the alignment for
		#sure we do not have that separation. Hence check if indel_real was the previous step
		
		#if ( $perchrom == 1  and $prev_step eq $cfg_hash->{'indel_real_step'} )
		#{
				###Using also Chromosome X and Y 
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				##An array will contain all the commands to be executed
				#my @commands = ();
				
			##Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
			##be better because all the bases are considered together
			##Build ith chromosome target name
			
			#if (file_not_present($target_bed) > 0 ){ print_and_log("WARNING no target bed file was uses using $prog_used!\n",$log_file);}
			##my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'all_chr'}.".".$cfg_hash->{'bed_ext'}." ";#ROMAN FILE
			#my $l_string = "";#USE THIS 
			
			#if ( $target_bed ne '-'){ 
				#$l_string = " --intervals $target_bed ";	
			#}
			
			##Build i-th bam output name
			#my $bqsr_file = $outFolder."/".extract_name($out_bqsr,1).".".$cfg_hash->{'recal_ext'};
			
			#my $in_bam_set = "";
			
			###Construct the input set using all the samples and all the chromosomes
			##I defined params here because I need it in the construction of output file name
			##there we need just one among all
			#my $params;
			#foreach my $readf_id ( sort keys %{$distinct_readf_ids} ) {
				##print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
				#getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readf_id);
				##Build the input file name. The step will not be there in the array, hence all the suffixes
				##(from steps performed) will be used
				#my $input_bam = $inFolder."/".build_input_name_from_executed($params,$curr_step,
								#$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
					##Build i-th input name
					#my $ith_in_bam = $inFolder."/".extract_name($input_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#if (file_not_present($ith_in_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $ith_in_bam.\n",$log_file);}
					#$in_bam_set .= " -I $ith_in_bam ";
				#}		
			#}

			##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			#my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." $in_bam_set $param_str $l_string -o $bqsr_file ";			 
			#print_and_log("Executing a thread for: $command\n",$log_file);
			##try_exec_command($command) or die "Unable to execute command: $command\n";
			#push(@commands,$command);
			##Execute the thread
			#execute_threads(\@commands);
		#}	
		#else{
			# When perchrom is not used, it uses the target regions as intervals	
			if (file_not_present($target_bed) > 0 ){ print_and_log(" no target bed file was used with $prog_used.\n",$log_file);}
			if ( $target_bed ne '-'){ 
				$param_str .= " --intervals $target_bed ";	
			}
			
			my $in_bam_set = "";
			##Construct the input set using all the samples and all the chromosomes
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			foreach my $readf_id (
				sort keys %{$distinct_readf_ids}
				)
				{
					#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
					getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
					#Build the input file name. The step will not be there in the array, hence all the suffixes 
					#(from steps performed) will be used
					my $input_bam = $inFolder."/".build_input_name_from_executed($params,$curr_step,
									$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};

				if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
				$in_bam_set .= " -I $input_bam ";				
			}	
			#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			my $command = $java_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_bqsr ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		#}
	#}else{log_and_exit( "ERROR: Step $curr_step has not been run for sample: $sample_id.f\n ",$log_file);}
	
}



=head2 run_GATK_PrintReads_persample

 Title   : run_GATK_PrintReads_persample
 Usage   : run_GATK_PrintReads_persample(   );

 Function: run_GATK_PrintReads_persample
			
			Is a tool from GATK which permits to realign the reads towards well assessed indels
			Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			avoid the next step of MergeSamFiles
			
			In summer 2019 I extended it to be executed per-chromosome so that this process can be parallelized
			on multiple nodes for WGS analysis.
 
 Returns : nothing
=cut
sub run_GATK_PrintReads_persample{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $out_bam = shift;
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

	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
	}else{$target_bed = "-";}
		
	#Check from the DB if the Sorting and indexing was performed
	if ( $params->{$prev_step} == 1 ){
		#inFolder for the bam file
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$prev_task.'_out_f'};
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};##
				
		#Build the base quality score recal file name using as stop point the current step
		my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};

		#Verifying that there has been an execution per chromosome
		my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);

		#Here we get the read file ids associated with this sample. They will be used for the input file names
		my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		#print_and_log( "Starting $prog_used for sample  [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		#Set the Java call
		my $java_call = $cfg_hash->{'java_path'};
		#Set the java memory
		if (defined $cfg_hash->{'javamem_PR'}){
			$java_call .= " -".$cfg_hash->{'javamem_PR'}." ";
		}
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}			
		#Set the program to use
		$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

		#Parameters string construction
		my $param_str = " ";
		#baq
		if ( defined $cfg_hash->{'baq_PR'} ){
			$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		}
		#Number of threads to use CANNOT USE IN INDELREALIGNER
		if ( defined $cfg_hash->{'nct_PR'} ){
			$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		}
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
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'},
																			$cfg_hash->{'db_sample_id'},$sample_id);	
		#if ( $perchrom == 1  and $indel_real_step == 1 ){
		print_and_log( "perchrom = $perchrom for sample  [$sample_id]..\n ",$log_file);#DEBUGCODE
		if ( $perchrom == 1  ){
		
			$numExec = $cfg_hash->{'chromosomes_num'};
				##Using also Chromosomex X and Y
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
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
						
										
					my $in_bam_set = "";
					
					##Construct the input set using all the samples for each chromosome
					#Get the reads files associated with this sample
					#I defined params here because I need it in the construction of output file name
					#there we need just one among all
					my $params;
					foreach my $readf_id (sort keys %{$distinct_readf_ids}){
						#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
						getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
						#Build the input file name. The step will not be there in the array, hence all the suffixes 
						#(from steps performed) will be used. The step before here is the step used because
						#the bam files did not change after base recalibration
						my $input_bam = $inFolder."/".build_input_name_from_executed($params,$prev_step,
										$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
			
						if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
						$in_bam_set .= " -I $input_bam ";							
					}
					#Analysis with one single BQSR for all BAM per chromosome. This kind of analysis should
					#be better because all the bases are considered together
					my $bqsr = " -BQSR ".$outFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					#Build i-th bam output name
					my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					my $command = "$java_call $in_bam_set -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr $l_string -o $ith_out_bam ";			 
					print_and_log("Executing a job for: $command\n",$log_file);
					push(@commands,$command);
			}
			#Path to the script to execute multiple jobs using an array of commands
			my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
						
			my $pr_jobs = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$task,"PR",$log_file);

			#Now alter the dependencies of the variant calling job: it has to wait the print reads job
			my $new_depend = "-W depend=afterok:$pr_jobs ";
			my $job_name = $cfg_hash->{'varcall_task_s'}."_s".$sample_id;
			my $VC_jobid = JOB_get_jobid_form_jobname($job_name,$cfg_hash->{'qsub_username'});
			print_and_log("Altering dependencies of $job_name ($VC_jobid):\n ",$log_file);
			alter_job($VC_jobid,$new_depend,$log_file);

			#Update the hash with dependencies
			my $deps_f = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'jc_dep_hash_f'};
			update_dependencies_hash($deps_f,$job_name,$VC_jobid,$pr_jobs,$log_file);
						
		}
		##FORSE QUESTA PARTE NON MI SERVE. NEL CASO DI MULTIPLE READS IL MERGE VIENE FATTO DURANTE 
		#PRINTREADS. SERVIREBBE SOLO SE PRINTREADS NON VIENE ESEGUITO 
		else{
			my $in_bam_set = "";
			
			##Construct the input set using all the samples and all the chromosomes
			#Get the reads files associated with this sample
			#I defined params here because I need it in the construction of output file name
			#there we need just one among all
			my $params;
			foreach my $readf_id ( sort keys %{$distinct_readf_ids}	) {
					#print_and_log("Readf id is: $readf_id..\n",$log_file);#DEBUGCODE
					getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readf_id);
					#Build the input file name. The step will not be there in the array, hence all the suffixes 
					#(from steps performed) will be used
					my $input_bam = $inFolder."/".build_input_name_from_executed($params,$step,
									$params->{$db_readf_name},$steps_array).".".$cfg_hash->{'bam_ext'};
		
				if (file_not_present($input_bam) > 0 ){ log_and_exit( "Cannot proceed with $prog_used! Check: $input_bam.\n",$log_file);}
				$in_bam_set .= " -I $input_bam ";									
			}
			#Check input files
			if (file_not_present($bqsr_file) ){ die "Cannot proceed with $prog_used! Check: $bqsr_file.\n";}
			#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			my $command = $java_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." -BQSR $bqsr_file $param_str -o $out_bam ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
		}

	}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
}




#=head2 run_GATK_PrintReads_persample

 #Title   : run_GATK_PrintReads_persample
 #Usage   : run_GATK_PrintReads_persample(   );

 #Function: run_GATK_PrintReads_persample
			
			#Is a tool from GATK which permits to realign the reads towards well assessed indels
			#Print recalibrated reads per-sample takes in input all the reads files (e.g. separated per-lane)
			#and the .grp file from BaseRecalibrator and generates a single .bam per sample with the reads
			#recalibrated. This step merges the .bam files per lane in a single .bam and hence permits to
			#avoid the next step of MergeSamFiles
 
 #Returns : nothing
#=cut
#sub run_GATK_PrintReads_persampleOLD{
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $group_id = shift;
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
		#my $inFolder = $cfg_hash->{$group_id.'_'.$prev_task.'_out_f'};
		##Out folder
		#my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
				
		##Build the base quality score recal file name using as stop point the current step
		#my $bqsr_file = $outFolder."/".build_input_name_from_executed($params,$step,
						#$params->{$db_sample_name},$steps_array2).".".$cfg_hash->{'recal_ext'};
				
		##Prepare output name for files for data recalibration
		#my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
						#$params->{$db_sample_name},$steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};

		##Verifying that there has been an execution per chromosome
		#my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    #,$cfg_hash->{'db_analysis_id'},$group_id);

		##Here we get the read file ids associated with this sample. They will be used for the input file names
		#my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
						#,$cfg_hash->{'db_analysis_id'},$group_id,$cfg_hash->{'db_sample_id'},$sample_id);						    		
		
		##print_and_log( "Starting $prog_used for sample ".$params->{$db_sample_name}." [$sample_id] using GATK tools..\n ",$log_file);#DEBUGCODE
		
		##Set the Java call
		#my $java_call = $cfg_hash->{'java_path'};
		##Set the java memory
		#if (defined $cfg_hash->{'javamem_PR'}){
			#$java_call .= " -".$cfg_hash->{'javamem_PR'}." ";
		#}
		##Add temporary folder
		#if ( defined $cfg_hash->{'scratch_f'} ){
			#$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		#}			
		##Set the program to use
		#$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

		##Parameters string construction
		#my $param_str = " ";
		##baq
		#if ( defined $cfg_hash->{'baq_PR'} ){
			#$param_str .= " -baq ".$cfg_hash->{'baq_PR'}." ";
		#}
		##Number of threads to use CANNOT USE IN INDELREALIGNER
		#if ( defined $cfg_hash->{'nct_PR'} ){
			#$param_str .= " -nct ".$cfg_hash->{'nct_PR'}." ";
		#}
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
				###Using also Chromosomex X and Y
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
					#my  $bqsr = " -BQSR ".$inFolder."/".extract_name($bqsr_file,1).".".$cfg_hash->{'recal_ext'};
								
					##Build i-th bam output name
					#my $ith_out_bam = $outFolder."/".extract_name($out_bam,1)."_$exec.".$cfg_hash->{'bam_ext'};
					##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					#my $command = $java_call." $in_bam_set  -R ".$cfg_hash->{'hum_ref'}." $param_str $bqsr -o $ith_out_bam ";			 
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
			##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			#my $command = $java_call." $in_bam_set -R ".$cfg_hash->{'hum_ref'}." -BQSR $bqsr_file $param_str -o $out_bam ";			 
			#print_and_log("Executing command: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
		#}

	#}else{log_and_exit( "ERROR: Step $prev_step has not been run for sample: $sample_id.\n ",$log_file);}
	
#}





=head2 run_GATK_HaplotypeCaller

 Title   : run_GATK_HaplotypeCaller
 Usage   : run_GATK_HaplotypeCaller(   );

 Function: run_GATK_HaplotypeCaller
			Is a tool from GATK which permits to call variants. Haplotype caller works
			per-sample hence if until this moment we have worked with files separated
			by lanes a step of Merge is necessary.
			A check for this is made when the distinct samples are fetched if they exist.
			In that case all the steps for each file should be completed.
			
			PARALLEL version of HC permits to return jobs to launch per-chromosome
			
 Returns : nothing
=cut
sub run_GATK_HaplotypeCaller_parallel{
  my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $input_file = shift;
	my $out_file = shift;
	my $inFolder = shift;
	my $outFolder = shift;
	my $log_file = shift;

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

	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
								$cfg_hash->{'db_sample_id'},$sample_id);		
	
	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);


	#Getting the pedigree file
	my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
								$cfg_hash->{'db_analysis_id'},$analysis_id);	
											
	#print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Prepare output name
	my $param_str = " ";	
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_HAPC'}){
		$java_call .= " -".$cfg_hash->{'javamem_HAPC'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	#nct: number of processors to use
	if ( defined $cfg_hash->{'varcall_threads'} ){
					$param_str .= " -nct ".$cfg_hash->{'varcall_threads'};
	}
	if ( defined $cfg_hash->{'emitRefConfidence'} and $dogenot == 1){
					$param_str .= " --emitRefConfidence ".$cfg_hash->{'emitRefConfidence'};
	}
	if ( defined $cfg_hash->{'variant_index_type'} ){
					$param_str .= " --variant_index_type ".$cfg_hash->{'variant_index_type'};
	}
	if ( defined $cfg_hash->{'variant_index_parameter'} ){
					$param_str .= " --variant_index_parameter ".$cfg_hash->{'variant_index_parameter'};
	}
	if ( defined $cfg_hash->{'HC_stand_call_conf'} and correct_type($cfg_hash->{'HC_stand_call_conf'},"positiveint") ){
					$param_str .= " -stand_call_conf ".$cfg_hash->{'HC_stand_call_conf'};
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
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);			
	
	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);

	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
		if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}		
	}else{$target_bed = "-";}
											
	# When perchrom is not used the HaplotypeCaller uses the target regions as intervals	
	if ( $perchrom == 0 ){
				
		if ( $target_bed ne '-'){ 
			$param_str .= " --intervals $target_bed ";	
		}
	}

	if ( defined $cfg_hash->{'hapc_dbsnp'} ){
		$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	}		
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'hapc_annotation'} ){
					my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'hapc_annotation'});
					foreach my $annotation (@annotations){
							#If annotation is inbreedingCoeff we need the ped file
							if ( $annotation eq 'InbreedingCoeff'){
								if ($ped_file ne 'none'){
									my $ped_f_path = $cfg_hash->{$analysis_id.'_data_fold'}."/$ped_file";
									if (  -e  $ped_f_path ){
										$param_str .= " -A $annotation ";
										$param_str .= " -ped $ped_f_path ";
									}else{
											print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									}									
								}else{
											print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									}	
						}else{
								$param_str .= " -A $annotation ";
							}
					}
	}

	#Preparation of jobs for parallel execution per-chromosome
	my $numExec = 1;
	if ( $perchrom == 1 ){
				##Using also Chromosomex X and Y
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));

									
				#An array will contain all the commands to be executed
				my @commands = ();

				#Path to the script to execute multiple jobs using an array of commands
				my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
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
					#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					my $command = $java_call." -I $ith_input -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_out ";			 
					print_and_log("Preparing command: $command \n",$log_file);
					push(@commands,$command);
			}
			#Run a different job for each command
			my $task = $cfg_hash->{'varcall_task_s'};
			my $jobs_to_wait = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$cfg_hash->{'varcall_task'},$cfg_hash->{'varcall_task_s'},$log_file);
			
			#Return the job ids for each job and the list of output file and the final output name
			chop($out_files_str);
			return $jobs_to_wait,$out_files_str,$out_file;		
	}
	#All the target is used
	else{
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			my $command = $java_call." -I $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_file ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
	}
}

=head2 run_GATK_HaplotypeCaller

 Title   : run_GATK_HaplotypeCaller
 Usage   : run_GATK_HaplotypeCaller(   );

 Function: run_GATK_HaplotypeCaller
			Is a tool from GATK which permits to call variants. Haplotype caller works
			per-sample hence if until this moment we have worked with files separated
			by lanes a step of Merge is necessary.
			A check for this is made when the distinct samples are fetched if they exist.
			In that case all the steps for each file should be completed.
			
 Returns : nothing
=cut
sub run_GATK_HaplotypeCaller{
  my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
	my $input_file = shift;
	my $out_file = shift;
	my $inFolder = shift;
	my $outFolder = shift;
	my $log_file = shift;

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

	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
								$cfg_hash->{'db_sample_id'},$sample_id);		
	
	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);


	#Getting the pedigree file
	my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
								$cfg_hash->{'db_analysis_id'},$analysis_id);	
											
	#print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Prepare output name
	my $param_str = " ";	
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_HAPC'}){
					$java_call .= " -".$cfg_hash->{'javamem_HAPC'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	#nct: number of processors to use
	if ( defined $cfg_hash->{'varcall_threads'} ){
					$param_str .= " -nct ".$cfg_hash->{'varcall_threads'};
	}
	if ( defined $cfg_hash->{'emitRefConfidence'} and $dogenot == 1){
					$param_str .= " --emitRefConfidence ".$cfg_hash->{'emitRefConfidence'};
	}
	if ( defined $cfg_hash->{'variant_index_type'} ){
					$param_str .= " --variant_index_type ".$cfg_hash->{'variant_index_type'};
	}
	if ( defined $cfg_hash->{'variant_index_parameter'} ){
					$param_str .= " --variant_index_parameter ".$cfg_hash->{'variant_index_parameter'};
	}
	if ( defined $cfg_hash->{'HC_stand_call_conf'} and correct_type($cfg_hash->{'HC_stand_call_conf'},"positiveint") ){
					$param_str .= " -stand_call_conf ".$cfg_hash->{'HC_stand_call_conf'};
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
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);			
	
	# When perchrom is not used the HaplotypeCaller uses the target regions as intervals	
	if ( $perchrom == 0 ){
			my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
								,$cfg_hash->{'db_analysis_id'},$analysis_id);
			if ($target_bed ne '-'){
				$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
			}else{$target_bed = "-";}
			
			if (file_not_present($target_bed) > 0 ){ print_and_log("no target bed file was used with $prog_used!\n",$log_file);}
				
		if ( $target_bed ne '-'){ 
			$param_str .= " --intervals $target_bed ";	
		}
	}

	if ( defined $cfg_hash->{'hapc_dbsnp'} ){
		$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	}		
	#Adds known indel snp to the list of known
	if ( defined $cfg_hash->{'hapc_annotation'} ){
					my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'hapc_annotation'});
					foreach my $annotation (@annotations){
							#If annotation is inbreedingCoeff we need the ped file
							if ( $annotation eq 'InbreedingCoeff'){
								if ($ped_file ne 'none'){
									my $ped_f_path = $cfg_hash->{$analysis_id.'_data_fold'}."/$ped_file";
									if (  -e  $ped_f_path ){
										$param_str .= " -A $annotation ";
										$param_str .= " -ped $ped_f_path ";
									}else{
											print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									}									
								}else{
											print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									}	
						}else{
								$param_str .= " -A $annotation ";
							}
					}
	}


	#my $numExec = 1;
	#if ( $perchrom == 1 ){
				###Using also Chromosomex X and Y
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));

									
				##An array will contain all the commands to be executed
				#my @commands = ();
				
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
					##Build ith chromosome target name
					#my $l_string = "";
					##WES
					#if (  $target_bed ne "-" ) {
						#$l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
					#}else{
						#$l_string = " -L ".$cfg_hash->{'chr_suff'}."$exec ";
					#}

					##my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
					#my $ith_out = $outFolder."/".extract_name($out_file,1)."_$exec.".$out_ext;
					#my $ith_input = $inFolder."/".extract_name($input_file,1)."_".$exec.".".$cfg_hash->{'bam_ext'};
					#if (file_not_present($ith_input) > 0 ){ die "Cannot proceed with $prog_used! Check: $ith_input.\n";}
					##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					#my $command = $java_call." -I $ith_input -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_out ";			 
					#print_and_log("Executing command: $command \n",$log_file);
					#push(@commands,$command);
			#}
			#execute_threads(\@commands);
	#}else{
			if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			my $command = $java_call." -I $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_file ";			 
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";
	#}
}



#=head2 run_GATK_HaplotypeCaller

 #Title   : run_GATK_HaplotypeCaller
 #Usage   : run_GATK_HaplotypeCaller(   );

 #Function: run_GATK_HaplotypeCaller
			#Is a tool from GATK which permits to call variants. Haplotype caller works
			#per-sample hence if until this moment we have worked with files separated
			#by lanes a step of Merge is necessary.
			#A check for this is made when the distinct samples are fetched if they exist.
			#In that case all the steps for each file should be completed.
			
 #Returns : nothing
#=cut
#sub run_GATK_HaplotypeCallerOLD{
  #my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $sample_id = shift;
	#my $group_id = shift;
	#my $log_file = shift;
	#my $in_fold_suff = shift;
	#my $out_fold_suff = shift;#Used to know the folder where to take the input bam files
  #my $step = shift;
  #my $step_needed  = shift;
  #my $steps_array = shift;
  #my $steps_array2 = shift;
       
	##my $previous_task = '';
	#my $steps_arr1_in_use;
	#my $paramsIn;
	#my $paramsOut;
	#my $main_name = $cfg_hash->{'db_sample_name'};

	##Case IR only or no refinement step: there is a single file and only IndelRealignment or no refinement step was executed: 
	##construct the name using readfile level fields
	#if ( $step_needed eq $cfg_hash->{'indel_real_step'} or $step_needed eq $cfg_hash->{'sort_idx_step'} ){
		#$steps_arr1_in_use = $steps_array;
		#my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
								#$cfg_hash->{'db_analysis_id'},$group_id,$cfg_hash->{'db_sample_id'},$sample_id);
		#foreach my $dist_sample ( keys %{$distinct_samples}){
			#getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							#$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$dist_sample);
		#}
	#}	
	##CASE: IR-BR or IR-BR-M or IR-M: there are multiple read files or the last step is BR
	##construct the name using sample level fields
	#elsif ($step_needed eq $cfg_hash->{'mrdup_groups_step'} or $step_needed eq $cfg_hash->{'print_reads_step'}){
		#$steps_arr1_in_use = $steps_array2;
		#getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	#}
	
	##Get the sample parameters from the the db for the output file
	#getSampleConfiguration_locked(\$paramsOut,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	##Getting the sample name
	#my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    #$cfg_hash->{'db_sample_id'},$sample_id);		
	
	##Verifying that there has been an execution per chromosome
	#my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
					    #$cfg_hash->{'db_analysis_id'},$group_id);
																	
	##InFolder is the alignment one at this stage
	#my $inFolder = $cfg_hash->{$group_id.'_'.$in_fold_suff.'_out_f'};
	#my $outFolder = $cfg_hash->{$group_id.'_'.$out_fold_suff.'_out_f'};

	##Getting the pedigree file
	#my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
						    #$cfg_hash->{'db_analysis_id'},$group_id);	


	##Getting the flag indicating if the joint genotype has to be run
	#my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
						#$cfg_hash->{'db_analysis_id'},$group_id);
	#my $out_ext;
	##If it has been run, take the output form the genotype, else from the varcalling
	#if ($dogenot){
		#$out_ext = $cfg_hash->{'gvcf_ext'};
	#}else{
		#$out_ext = $cfg_hash->{'vcf_ext'};
	#}
																				
	##Set the input and output file names
	#my $input_file = $inFolder."/".build_input_name_from_executed($paramsIn,$step,$sample_name,$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'};
	#my $out_file = $outFolder."/".build_input_name_from_executed($paramsOut,$step,$sample_name,$steps_array2)."_".$step.".".$out_ext;
	
	##print_and_log( "Starting $prog_used for $input_file using GATK tools..\n ",$log_file);#DEBUGCODE
	
	##Prepare output name
	#my $param_str = " ";	
	
	##Set the Java call
	#my $java_call = $cfg_hash->{'java_path'};
	##Set the java memory
	#if (defined $cfg_hash->{'javamem_HAPC'}){
					#$java_call .= " -".$cfg_hash->{'javamem_HAPC'}." ";
	#}
	##Add temporary folder
	#if ( defined $cfg_hash->{'scratch_f'} ){
		#$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	#}	
	##Set the program to use
	#$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	##nct: number of processors to use
	#if ( defined $cfg_hash->{'varcall_threads'} ){
					#$param_str .= " -nct ".$cfg_hash->{'varcall_threads'};
	#}
	#if ( defined $cfg_hash->{'emitRefConfidence'} and $dogenot == 1){
					#$param_str .= " --emitRefConfidence ".$cfg_hash->{'emitRefConfidence'};
	#}
	#if ( defined $cfg_hash->{'variant_index_type'} ){
					#$param_str .= " --variant_index_type ".$cfg_hash->{'variant_index_type'};
	#}
	#if ( defined $cfg_hash->{'variant_index_parameter'} ){
					#$param_str .= " --variant_index_parameter ".$cfg_hash->{'variant_index_parameter'};
	#}
	#if ( defined $cfg_hash->{'HC_stand_call_conf'} and correct_type($cfg_hash->{'HC_stand_call_conf'},"positiveint") ){
					#$param_str .= " -stand_call_conf ".$cfg_hash->{'HC_stand_call_conf'};
	#}
	#if ( defined $cfg_hash->{'HC_stand_emit_conf'} and correct_type($cfg_hash->{'HC_stand_emit_conf'},"positiveint") ){
					#$param_str .= " -stand_emit_conf ".$cfg_hash->{'HC_stand_emit_conf'};
	#}
	#if ( defined $cfg_hash->{'HC_min_base_quality_score'} and correct_type($cfg_hash->{'HC_min_base_quality_score'},"positiveint") ){
					#$param_str .= " --min_base_quality_score ".$cfg_hash->{'HC_min_base_quality_score'};
	#}
	#if ( defined $cfg_hash->{'HC_dontUseSoftClippedBases'} and  $cfg_hash->{'HC_dontUseSoftClippedBases'} eq 'YES' ){
					#$param_str .= " --dontUseSoftClippedBases ";
	#}
	#if ( defined $cfg_hash->{'hapc_dbsnp'} ){
					#$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	#}	
	

	#my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'}
									#,$cfg_hash->{'db_analysis_id'},$group_id);			
	## When perchrom is not used the HaplotypeCaller uses the target regions as intervals	
	#if ( $perchrom == 0 ){
		##Getting the target bed file using the folder for targets and the name contained in the database
		#my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
						#,$cfg_hash->{'db_analysis_id'},$group_id);
			#if (file_not_present($target_bed) > 0 ){ print_and_log("WARNING no target bed file was uses using $prog_used!\n",$log_file);}
				
		#if ( $target_bed ne '-'){ 
			#$param_str .= " --intervals $target_bed ";	
		#}
	#}

	#if ( defined $cfg_hash->{'hapc_dbsnp'} ){
					#$param_str .= " -D:dbsnp,VCF ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'hapc_dbsnp'}}
	#}		
	##Adds known indel snp to the list of known
	#if ( defined $cfg_hash->{'hapc_annotation'} ){
					#my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'hapc_annotation'});
					#foreach my $annotation (@annotations){
							##If annotation is inbreedingCoeff we need the ped file
							#if ( $annotation eq 'InbreedingCoeff'){
								#if ($ped_file ne 'none'){
									#my $ped_f_path = $cfg_hash->{$group_id.'_data_fold'}."/$ped_file";
									#if (  -e  $ped_f_path ){
										#$param_str .= " -A $annotation ";
										#$param_str .= " -ped $ped_f_path ";
									#}else{
											#print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									#}									
								#}else{
											#print_and_log("Pedigree file does not exist. Cannot use annotation: $annotation\n",$log_file);
									#}	
						#}else{
								#$param_str .= " -A $annotation ";
							#}
					#}
	#}


	#my $numExec = 1;
	#if ( $perchrom == 1 ){
				###Using also Chromosomex X and Y
				#$numExec = $cfg_hash->{'chromosomes_num'};
				#my @chroms = 1..$numExec;
				#push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				##An array will contain all the commands to be executed
				#my @commands = ();
				
				##Executing one GATK command for each target chromosome using bed files
				#foreach my $exec (@chroms){
				##for (my $exec = 1; $exec <= $numExec; $exec++){
					#my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
					#my $ith_out = $outFolder."/".extract_name($out_file,1)."_$exec.".$out_ext;
					#my $ith_input = $inFolder."/".extract_name($input_file,1)."_".$exec.".".$cfg_hash->{'bam_ext'};
					#if (file_not_present($ith_input) > 0 ){ die "Cannot proceed with $prog_used! Check: $ith_input.\n";}
					##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
					#my $command = $java_call." -I $ith_input -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $ith_out ";			 
					#print_and_log("Executing command: $command \n",$log_file);
					#push(@commands,$command);
			#}
			#execute_threads(\@commands);
	#}else{
			#if (file_not_present($input_file) > 0 ){ die "Cannot proceed with $prog_used! Check: $input_file.\n";}
			##print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
			#my $command = $java_call." -I $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -o $out_file ";			 
			#print_and_log("Executing command: $command\n",$log_file);
			#try_exec_command($command) or die "Unable to execute command: $command\n";
	#}
#}

=head2 run_GATK_CatVariants

 Title   : run_GATK_CatVariants
 Usage   : run_GATK_CatVariants(   );

 Function: run_GATK_CatVariants
			Is a tool from GATK which permits to merge together gvcf files from
			Haplotype caller so that we will have only one VCF file per sample
			
 Returns : nothing
=cut
sub run_GATK_CatVariants{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $task = shift;
	my $step = shift;
	my $steps_array = shift;
	my $get_command  =shift;	

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
		my $java_call = $cfg_hash->{'java_path'};
		#Set the java memory
		if (defined $cfg_hash->{'javamem_CV'}){
			$java_call .= " -".$cfg_hash->{'javamem_CV'}." ";
		}
		#Add temporary folder
		if ( defined $cfg_hash->{'scratch_f'} ){
			$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
		}		
		#Set the program to use
		$java_call .= " -cp ".$cfg_hash->{'gatk_path'}." ".$cfg_hash->{'gatk_tools_path'}.".$prog_used ";

		if ( $cfg_hash->{'assumeSorted'} eq 'YES'){
			$param_str .= " -assumeSorted ";
		}

		#Constructing input string with files to be merged 
		my $vars_to_cat = "";
		my @f_to_remove = ();
		my $numExec = 1;
		if ( $perchrom == 1){
				##Using also Chromosomex X and Y
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));

				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
					#Set the input file name to be used as component for the single files to be merged
					my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,$sample_name,$steps_array)."_$exec.".$out_ext;
					$vars_to_cat .=  " -V $input_file ";
					push(@f_to_remove,$input_file);
			}
		}
		if ( $cfg_hash->{'remove_temp'} eq 'YES'){
			foreach my $f (@f_to_remove){
				delete_file($f);
			}		
		}
		
		#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
		my $command = $java_call." $vars_to_cat -R ".$cfg_hash->{'hum_ref'}." $param_str --outputFile $outFile ";			 

		if ($get_command ne 'GET'){
			print_and_log("Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";	
		}
		return $command;			

	}else{log_and_exit("ERROR: step $prev_step_done has not been run for $sample_id. Cannot run $step process..\n",$log_file);}
}


=head2 run_GATK_GenotypeGVCFs

 Title   : run_GATK_GenotypeGVCFs
 Usage   : run_GATK_GenotypeGVCFs(   );

 Function: run_GATK_GenotypeGVCFs
			Given the gvcf files from all the samples in a groups will make a Joint
			Genotype using GATK tools
			
 Returns : nothing
=cut
sub run_GATK_GenotypeGVCFs{
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
		$param_str .= " -stand_call_conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
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
			$input_file_str .= " -V ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
		}
	}else{

		
		#Get infos for the single sample associated to this group from the 'samples' table
		#print_and_log("Getting information for $group_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
							$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash-> {'db_analysis_id'},$group_id);	
		$input_file_str = " -V ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	}
	
	#Set the unique output name
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					$group_name,$steps_array3)."_".$step.".".$cfg_hash->{'gvcf_ext'};

	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);

	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_GGVCF'}){
		$java_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

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

	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	my $command = $java_call." $input_file_str -R ".$cfg_hash->{'hum_ref'}." $param_str -o $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

=head2 run_GATK_GenotypeGVCFs_parallel

 Title   : run_GATK_GenotypeGVCFs_parallel
 Usage   : run_GATK_GenotypeGVCFs_parallel(   );

 Function: run_GATK_GenotypeGVCFs_parallel
			Given the gvcf files from all the samples in a groups will make a Joint
			Genotype using GATK tools
			
 Returns : Return the job ids for each job and the list of output file and the final output name
=cut
sub run_GATK_GenotypeGVCFs_parallel{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $analysis_id = shift;
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
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);	
																																			
	#Here we get the samples names associated with that group id
	my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $params;

	#input folder setting
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	
	
	#Prepare parameters string
	my $param_str = " ";	
		
	#If there are more samples associated to this group
	#Constructing input string with files to be merged 
	my $input_file_str = "";
	
	#If GenotypeGVCFs is executed for a single sample we reduce the threshold for confidence call
	#(https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller)
	#Set the phred-scaled Qscore threshold to separate high confidence from low confidence calls
	if (defined $cfg_hash->{'ggvcf_single_sam_stand_call_conf'} and correct_type($cfg_hash->{'ggvcf_single_sam_stand_call_conf'},"positiveint") ) {
		$param_str .= " -stand_call_conf ".$cfg_hash->{'ggvcf_single_sam_stand_call_conf'}." ";
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
			$input_file_str .= " -V ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
		}
	}else{

		
		#Get infos for the single sample associated to this group from the 'samples' table
		#print_and_log("Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
							$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash-> {'db_analysis_id'},$analysis_id);	
		$input_file_str = " -V ".$inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$sample_name},$steps_array2).".".$cfg_hash->{'gvcf_ext'}." ";
	}
	


	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);

	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_GGVCF'}){
		$java_call .= " -".$cfg_hash->{'javamem_GGVCF'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

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


	##Using also Chromosomex X and Y
	my $numExec = $cfg_hash->{'chromosomes_num'};
	my @chroms = 1..$numExec;
	push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
				
	#An array will contain all the commands to be executed
	my @commands = ();

	#Set the output name
	my $outFile = $outFolder."/".build_input_name_from_executed($params,$step,
					$analysis_name,$steps_array3)."_".$step.".".$cfg_hash->{'gvcf_ext'};
					
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	my $out_files_str = "";
	#Executing one GATK command for each target chromosome using bed files
	foreach my $exec (@chroms){
			my $l_string = " -L ".$cfg_hash->{'chrom_reg_f'}."/".$cfg_hash->{'chr_suff'}."$exec.".$cfg_hash->{'bed_ext'}." ";
			
			#Set the output name
			my $nth_out = $outFolder."/".build_input_name_from_executed($params,$step,
					$analysis_name,$steps_array3)."_".$step."_$exec.".$cfg_hash->{'gvcf_ext'};
			$out_files_str .= $nth_out.",";					
			my $command = $java_call." $input_file_str -R ".$cfg_hash->{'hum_ref'}." $param_str $l_string -o $nth_out ";			 
			#print_and_log("Preparing command: $command\n",$log_file);
			push(@commands,$command);
	}
	#Run a different job for each command
	my $jobs_to_wait = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,"no_depend",$task,$task,$log_file);
	
	#Return the job ids for each job and the list of output file and the final output name
	chop($out_files_str);
	return $jobs_to_wait,$out_files_str,$outFile;
}

=head2 run_GATK_VariantFiltration

 Title   : run_GATK_VariantFiltration
 Usage   : run_GATK_VariantFiltration(   );

 Function: run_GATK_VariantFiltration
			Given the gvcf files for the single analysis it will filter based on 
			parameters given in input.
			The analysis name here is in input since we do not need sample informations
			
 Returns : nothing 
=cut
sub run_GATK_VariantFiltration{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $param_str = shift;	
	my $filters = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_VF'}){
		$java_call .= " -".$cfg_hash->{'javamem_VF'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}		
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
	#Get the filter set from input to filter the INFO field
	my @filter_set = split(",",$filters);
	foreach my $filter ( @filter_set){
		if ( defined $cfg_hash->{'filterExpression'.$filter} ){
			$param_str .= " --filterExpression ".$cfg_hash->{'filterExpression'.$filter}." ";
		}	
		if ( defined $cfg_hash->{'filterName'.$filter} ){
			$param_str .= " --filterName ".$cfg_hash->{'filterName'.$filter}." ";
		}		
		#Genotype filters variables have _gX
		if ( defined $cfg_hash->{'filterExpression_'.$filter} ){
			$param_str .= " --genotypeFilterExpression ".$cfg_hash->{'filterExpression_'.$filter}." ";
		}	
		if ( defined $cfg_hash->{'filterName_'.$filter} ){
			$param_str .= " --genotypeFilterName ".$cfg_hash->{'filterName_'.$filter}." ";
		}				
	}

	my $command = $java_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -o $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK_SelectVariants

 Title   : run_GATK_SelectVariants
 Usage   : run_GATK_SelectVariants(   );

 Function: run_GATK_SelectVariants
			Subsets a VCF in order to facilitate analyses
			
			It is used here to select INDELS and SNP to two different files
			but it has many purposes
			
 Returns : nothing 
=cut
sub run_GATK_SelectVariants{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $param_str = shift;
	my $input_file = shift;
	my $outFile = shift;
	my $log_file = shift;

	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_SV'}){
		$java_call .= " -".$cfg_hash->{'javamem_SV'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}		
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	my $command = $java_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -o $outFile ";			 
	print_and_log("[gr:$analysis_id] Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

=head2 run_GATK_GCContentByInterval

 Title   : run_GATK_GCContentByInterval
 Usage   : run_GATK_GCContentByInterval(   );

 Function: run_GATK_GCContentByInterval
					Executes GATK GCContentByInterval
			
 Returns : nothing 
=cut
sub run_GATK_GCContentByInterval{
	my $cfg_hash = shift;
	my $targetbed = shift;
	my $outFile = shift;	
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_VA'}){
		$java_call .= " -".$cfg_hash->{'javamem_VA'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}		
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T ".$cfg_hash->{'gatk_gccontent_prog'}." ";
	
	
	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}."  --intervals $targetbed -o $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK_CalculateGenotypePosteriors

 Title   : run_GATK_CalculateGenotypePosteriors
 Usage   : run_GATK_CalculateGenotypePosteriors(   );

 Function: run_GATK_CalculateGenotypePosteriors
					This function derives the posteriors of genotype calls in our callset, 
					recalibratedVariants.vcf, which just came out of the VQSR filtering step;
					 it contains among other samples a trio of individuals (mother, father and child) 
					 whose family structure is described in the pedigree file trio.ped (which you need to supply).
					 To do this, we are using the most comprehensive set of high confidence SNPs available to us,
					 a set of sites from Phase 3 of the 1000 Genomes project (available in our resource bundle),
					 which we pass via the --supporting argument.
			
 Returns : nothing 
=cut
sub run_GATK_CalculateGenotypePosteriors{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $group_name = shift;
	my $group_id = shift;
	my $input_file = shift;
	my $ped_file = shift;
	my $outFile = shift;
	my $param_str = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);##DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_CGP'}){
		$java_call .= " -".$cfg_hash->{'javamem_CGP'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
	my $command = $java_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -ped $ped_file -o $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK_VariantAnnotator

 Title   : run_GATK_VariantAnnotator
 Usage   : run_GATK_VariantAnnotator(   );

 Function: run_GATK_VariantAnnotator
					Executes VariantAnnotator
			
 Returns : nothing 
=cut
sub run_GATK_VariantAnnotator{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $group_name = shift;
	my $group_id = shift;
	my $input_file = shift;
	my $outFile = shift;	
	my $ped_file = shift;
	my $log_file = shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_VA'}){
		$java_call .= " -".$cfg_hash->{'javamem_VA'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}		
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
	#Set the filters by searching in the user file
	my $param_str = " ";	

	#Adds annotations to use
	if ( defined $cfg_hash->{'GRW_annotation'} ){
		my @annotations  = split($cfg_hash->{'word_sep'},$cfg_hash->{'GRW_annotation'});
		foreach my $annotation (@annotations){
					$param_str .= " -A $annotation ";
		}
	}
		
	my $command = $java_call." --variant $input_file -R ".$cfg_hash->{'hum_ref'}." $param_str -ped $ped_file -o $outFile ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 run_GATK_CombineVariants

 Title   : run_GATK_CombineVariants
 Usage   : run_GATK_CombineVariants(   );

 Function: run_GATK_CombineVariants
						combines multiple variant records present at the same site in the different 
						input sources into a single variant record in the output. If sample names overlap
						, then they are "uniquified" by default, which means a suffix is appended to make them unique
						
						If get_command is given and is "YES", it returns a string that is the command to run. This
						is useful when you want to run a job with that command only.
						
			
 Returns : nothing 
=cut
sub run_GATK_CombineVariants{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $input_vcfs = shift;
	my $outFile = shift;
	my $log_file = shift;
  my $get_command  =shift;
	
	#print_and_log( "Starting $prog_used using GATK tools..\n ",$log_file);#DEBUGCODE
	
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_CV'}){
		$java_call .= " -".$cfg_hash->{'javamem_CV'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
	#Set the filters by searching in the user file
	my $param_str = " ";	
	
	my @vcf_set = split(",",$input_vcfs);
	
	foreach my $vcf ( @vcf_set ){
		$param_str .= " --variant $vcf ";		
	}

	if (defined $cfg_hash->{'cv_genotypeMergeOptions'}){
		$param_str .= " -genotypeMergeOptions ".$cfg_hash->{'cv_genotypeMergeOptions'}." ";
	}
	
	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." $param_str -o $outFile ";			 
	
	if ($get_command ne 'YES'){
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";	
  }
	return $command;
}

=head2 run_GATK_VariantRecalibration

 Title   : run_GATK_VariantRecalibration
 Usage   : run_GATK_VariantRecalibration(   );

 Function: run_GATK_VariantRecalibration
			Is a tool from GATK which permits to recalibrate variants. It can work
			per-sample 
			
 Returns : nothing
=cut
sub run_GATK_VariantRecalibration{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
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
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_VR'}){
					$java_call .= " -".$cfg_hash->{'javamem_VR'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}		
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	if ( defined $cfg_hash->{'vr_std'} ){
					$param_str .= " --stdThreshold ".$cfg_hash->{'vr_std'};
	}
	#Max Gaussians is the number of different cluster to made with variants 
	if ( defined $cfg_hash->{"vr_max_gaussians_$mode"} ){
					$param_str .= " --maxGaussians ".$cfg_hash->{"vr_max_gaussians_$mode"};
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


	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." -input $input_file $res_string $ann_string -mode $mode -recalFile $recal_file -tranchesFile $tranches_file ".
								" -rscriptFile $rscript_file  $param_str";			 
	print_and_log("Executing command: $command\n",$log_file);
	
	if (file_not_present($input_file) > 0 ){ warn "Could not proceed with $prog_used! Check: $input_file.\n";}
	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_GATK_ApplyRecalibration

 Title   : run_GATK_ApplyRecalibration
 Usage   : run_GATK_ApplyRecalibration(   );

 Function: run_GATK_ApplyRecalibration
			Is a tool from GATK which permits to recalibrate variants. It can work
			per-sample 
			
 Returns : nothing
=cut
sub run_GATK_ApplyRecalibration{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $sample_id = shift;
	my $group_id = shift;
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
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_AR'}){
					$java_call .= " -".$cfg_hash->{'javamem_AR'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";

	#if ( defined $cfg_hash->{'refine_threads'} ){
					#$param_str .= " -nt ".$cfg_hash->{'refine_threads'};
	#}
	if ( defined $cfg_hash->{'ar_ts_filter_level_'.$mode} ){
		$param_str .= " --ts_filter_level ".$cfg_hash->{'ar_ts_filter_level_'.$mode};
	}
	#Number of threads to use CANNOT BE USED IN ApplyRecalibration
	#if ( defined $cfg_hash->{'nct_PR'} ){
		#$param_str .= " -nct ".$cfg_hash->{'nct_PR'};
	#}
	
	#my $retval = file_not_present($input_file);
	#if ( $retval  > 0 ){ log_and_exit( "Cannot proceed with $prog_used (ERROR: $retval)! Check: $input_file.\n",$log_file);}
		
	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	my $command = $java_call." --input $input_file -R ".$cfg_hash->{'hum_ref'}." -mode $mode -recalFile $recal_file -tranchesFile $tranches_file ".
								"$param_str  -o $out_vcf ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 run_GATK_PhaseByTransmission

 Title   : run_GATK_PhaseByTransmission
 Usage   : run_GATK_PhaseByTransmission(   );

 Function: run_GATK_PhaseByTransmission
			Given the filtered gvcf files and parent/child information for the single group it will make the analysis
			of phasing
			
 Returns : nothing
=cut
sub run_GATK_PhaseByTransmission{
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
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_PBT'}){
		$java_call .= " -".$cfg_hash->{'javamem_PBT'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
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
		
	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." --variant $input_file -ped $ped_f_path $param_str -o $out_vcf ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or log_and_exit( "Unable to execute command: $command\n",$log_file);	
}


=head2 run_GATK_DepthOfCoverage

 Title   : run_GATK_DepthOfCoverage
 Usage   : run_GATK_DepthOfCoverage(   );

 Function: run_GATK_DepthOfCoverage
						Given a BAM file of aligned reads and a file with target regions
						it uses DepthOfCoverage to generate statistics 
			
 Returns : nothing
=cut
sub run_GATK_DepthOfCoverage{
	my $cfg_hash = shift;
	my $bam_input = shift;#This can also be a file with paths
	my $out_suffix = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	  
  #Program name
  my $prog_used = $cfg_hash->{'depth_coverage_prog'};
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_DOC'}){
		$java_call .= " -".$cfg_hash->{'javamem_DOC'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	
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
	my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	if ($target_bed ne '-'){
		$target_bed = $cfg_hash->{'target_reg_f'}."/".$target_bed;
	}else{$target_bed = "-";}
	if (file_not_present($target_bed)  > 0){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
	
	#Construct the target interval list file from the target file
	my $target_int_list = $target_bed.".".$cfg_hash->{'target_inter_list_ext'};
	
	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." -I $bam_input -L $target_int_list $param_str -o $out_suffix ";			 
	print_and_log("Executing command: $command\n",$log_file);
	try_exec_command($command) or log_and_exit( "Unable to execute command: $command\n",$log_file);	
}

=head2 run_GATK_DepthOfCoverage_gen

 Title   : run_GATK_DepthOfCoverage_gen
 Usage   : run_GATK_DepthOfCoverage_gen(   );

 Function: run_GATK_DepthOfCoverage_gen
						Given a BAM file of aligned reads and a file with target regions
						it uses DepthOfCoverage to generate statistics 
			
 Returns : nothing
=cut
sub run_GATK_DepthOfCoverage_gen{
	my $cfg_hash = shift;
	my $bam_input = shift;#This can also be a file with paths
	my $param_str = shift;
	my $targetbed = shift;
	my $out_suffix = shift;
	my $analysisid = shift;
	my $log_file = shift;
	my $get_command = shift;
	
	
  #Program name
  my $prog_used = $cfg_hash->{'depth_coverage_prog'};
	#Set the Java call
	my $java_call = $cfg_hash->{'java_path'};
	#Set the java memory
	if (defined $cfg_hash->{'javamem_DOC'}){
		$java_call .= " -".$cfg_hash->{'javamem_DOC'}." ";
	}
	#Add temporary folder
	if ( defined $cfg_hash->{'scratch_f'} ){
		$java_call .= " -Djava.io.tmpdir=".$cfg_hash->{'scratch_f'}." ";
	}	
	#Set the program to use
	$java_call .= " -jar ".$cfg_hash->{'gatk_path'}." -T $prog_used";
	

	if (file_not_present($targetbed)  > 0){ die "Cannot proceed with $prog_used! Check: $targetbed.\n";}
	

	#print_and_log( "$java_call - gatk params: $picard_params\n",$log_file);
	my $command = $java_call." -R ".$cfg_hash->{'hum_ref'}." -I $bam_input -L $targetbed $param_str -o $out_suffix ";	
	
	if ( !(defined $get_command) or ($get_command ne 'GET') ){
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or log_and_exit( "Unable to execute command: $command\n",$log_file);			
	}		 
	return $command;
}

1;
