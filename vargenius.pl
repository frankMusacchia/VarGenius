#!/usr/bin/perl

#vargenius.pl path
use lib '/pico/work/TELET_TIGEM/VarGeniusBeta/';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

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
    
    
=head1 NAME
VarGenius - Variant Discovery and Annotation Tool
   

=head1 SYNOPSIS

VarGenius is a tool for the detection and annotation of variants 

Please read carefully the user Readme and User Guide before to start executing it.
 
=cut


=head1 OPTIONS

=over 8

=item B<--help>
Shows the brief help informations

=item B<--version>
Shows the version

=back

=cut


=head1 EXAMPLES

  Write the config.txt file and give it in input to this script. In the config file you have to specify the locations of the data
  that you want to manage.

  The following is an example of call to this script:

 vargenius.pl  -c config.txt

=cut


=head1 DESCRIPTION
  
  Here you have the descriptions of the elements that you have to write on the file config.txt
           
   
=cut


=head1 AUTHORS

 Francesco Musacchia 
 --
 http://www.tigem.it

 $Id: VarGenius, 2017/06/01 11:00:12 Exp $
=cut

# NEEDED LIBRARIES
use strict;
use warnings;
#use FindBin;#To search libraries
#use lib $FindBin::Bin;#To search libraries
use Getopt::Long;
use Cwd;#To change work directory
use Data::Dumper;#To print the hashes
use Pod::Usage;#Used to write an usage
use File::Copy;#To manage files
use IO::Handle;#To immediately print with autoflush 
use Time::HiRes qw( time ); #To compute the running time of jobs

#Using a library to manage files
use LIB::files_management qw(  delete_file save_hash load_hash file_not_present
								file_name_wostrange_chars print_file extract_name
								append_file_2_file separate_bed_auto_and_sex_chr
								get_only_first_on_col_from_BED split_bedfile);

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(hash2ConfigFile configFile2Hash try_exec_command 
				try_exec_job initialize_folders print_and_log log_and_exit
				separate_input_ids remove_sample_sheet show_analyses_from_sample_sheet
				get_ped_file remove_analyses clean_folders save_analyses checkConfigVariables
				separate_bed_perchr show_analysesids_from_analysesnames show_samples_from_group
				check_config_parameters check_genedb_links get_flowcell_and_lane
				get_ped_file_for_enlarged_analysis store_and_compress_analyses
				exec_store_and_compress_analyses check_datasets_links
				show_samplesids_from_samplesnames get_extended_regions
				JOB_get_status_field JOB_is_active_queue);

#Using a library to manage configuration files and their global variables
use LIB::db_management qw(createDBAndTables getSampleConfiguration copyTable2DB
				exists_ss_locked insert_into_table_locked get_id_if_exists_from_db
				select_all_from_db db_present select_distinct_samples createTables
				 do_query_select_all update_table update_table_2
				insert_only_if_not_exists_locked get_count_from_selected);

#Using a library for picard functions management
use LIB::picard_management qw(run_PICARD_BedToIntervalList);

#Using a library for bedtools functions management
use LIB::bedtools_management qw(run_BEDTOOLS_intersectBed run_BEDTOOLS_sortBed);

#Using a library for standard utilities								
use LIB::std_lib qw(print_hash print_array correct_type);

#Using a library for managing CNVs
use LIB::cnv_management qw( calculate_gc_content );

								
#GLOBAL VARIABLES COMMON TO MANY PROGRAMS
my $configHash;#This hash substain all the life of the program
my $program_name = "vargenius.pl";
my $readme = "README.me";#Readme filename
my $run = 1;#Determines if the program should execute tasks or not (used to print the help)
my $program_config = "program_config.txt";#The path to the configuration file for the program
my $variables_file = "variables.txt";#The file with all the variables that must be present in the config files
my $config_folder = "CONFIGURATION";
my $license_file = "COPYRIGHT";
my $working_folder = '';#Folder where to work
my $session_folder = '';#Folder where all the analyses will be written
my $program_folder = '';#Folder of the program
my $foldersFile = "folders.txt";
my $program_version = "0.1";
my $libFolder = "LIB/";#folder with all the libraries
my $logFolder = "LOG";
my $dataFolder = "DATA";
my $main_data_folder = "";
my $main_log_folder = "";
my $log_file = "";
my $timesFile = "";
my $user_config = "";

#LOG nice printings
my $niceClosure = "\n=================================================\n";
my $configClosure = "\n*************************************************\n";
my $vargenius_title = ' _  _  _______  ______  ______ _______ __   _ _____ _     _ _______'."\n".
											' \  /  |_____| |_____/ |  ____ |______ | \  |   |   |     | |______'."\n".
											'  \/   |     | |    \_ |_____| |______ |  \_| __|__ |_____| ______|'."\n";

#Global variables
my $inputFolder = "";
my $paramFile = "";
my $sample_sheet = "";
my $create_database = "";
my $create_tables = "";
my $qualityCheck = "";
my $trimming = "";
my $assemQual = "";
my $alignment = "";
my $refining = "";
my $mergesam = "";
my $varcall = "";
my $genot = "";
my $varfilt = "";
my $phasing = "";
my $stats = "";
my $anstats = "";
my $annotate = "";
my $finalout = "";
my $fastfout = "";
my $vcfimport = "";
my $start_all = "";
my $genedb = "";
my $dldatasets = "";
my $geneann = "";
my $reannote = "";
my $remove_sample_sheet = "";
my $remove_analyses = "";
my $store_analyses = "";
my $vg_store_analyses = "";#Used to re-execute the command as a job
my $save_analyses = "";
my $show_analyses = "";
my $show_samples = "";
my $get_analyses_ids = "";
my $get_samples_ids = "";
my $clean_analyses = "";
my $clean_folder = "";
my $overwrite_conf = "";
my $update_freqs = "";
my $cnv_detection = "";#Flag to activate the CNV pipeline
my $cnv_out = "";#Flag to activate the CNV output generation

my $profile = "";
my $license = "";# Shows Open Source GNU GPL license
	
my $pipeline = "pipeline";

#This variable is useful to understand if the user runs an analysis or something
#other (like: --genedb )
my $do_pipeline_tasks = 1;

my $outFolder;
my $cpus = 0;
my $run_samples = "";#samples identifiers
my $run_analyses = "";#analyses identifiers
my $force_run = "";#Use this to force runs
my $run_name = "";#Name of the run with which you are working
my $exec_utils = 0;

#Configuration hash
my $cfg_hash;
my $config_file = "";
#Sample information hash
my $ss_hash;
my $sampleID_field = "sampleID";
my @ss_samples = ();
#Hash that contains the analyses and the samples to execute
my $exec_hash;
#Dependencies hash
my $deps_f;

my $qsub_command = "qsub";
my $programs_runner = "programs_runner.pl";
#my $jobs_cleaner = "jobs_cleaner.pl";

###############DB TABLES
my $sample_conf_table = "sample_conf";

#Tasks names
my @pipeline_tasks = ('qc','trim','align','refine','MS','varcall','genot','varfilt','phasing','stats',
											'anstats','finalout','fastfout');
#Tasks to exec
my @tasks_to_exec = ();

#####################################################MAIN#############################################

######################################################################################################
#LOG FILE CREATION
#Here we take the current time to get a log name different at each computation  
my $time = scalar(localtime);
$time =~ s/ /_/g;
$time =~ s/:/-/g;


#$log_file = $program_name."_exec_$time.log";
#$timesFile = $program_name."_times_$time.log";

STDOUT->autoflush(1);#This makes STDOUT "hot" in the sense that everything will be print immediately

parse_command_line_args();

#VarGenius permits both the usage of some utilities and the run of an analysis pipeline
#we set here the log file name depending on what is doing

if ( $exec_utils > 0){
	$log_file = "utils_$time.log";
}else{
	$log_file = join("_",@tasks_to_exec)."_$time.log";
}

if ( -e $user_config and $user_config ne ''){
	#print_and_log("\nStore configuration files parameters..\n",$log_file);#DEBUGCODE
	#Sets useful directories and stores the content of the configuration files
	set_global_var_and_confighash($user_config);

 #Print license
 if ( $license ){
			print_file($program_folder."/".$license_file);
			log_and_exit("Exiting..\n",$log_file);
	}  
  	
	#Utilities execution
	if ( $exec_utils > 0){
			
			#Other tasks to execute like utilities
			print_and_log("No tasks to execute! Checking if you want to make utils...\n\n",$log_file);
			
			#Remove completely a run
			if ( $remove_sample_sheet ne '' ){
				remove_sample_sheet($cfg_hash,$log_file,$remove_sample_sheet);
			}
			#Remove a single analysis
			if ( $remove_analyses ne '' ){
				remove_analyses($cfg_hash,$log_file,$remove_analyses);
			}
			#Given a run name (sample sheet file) gives all the analysis ids
			if ( $show_analyses ne '' ){
				show_analyses_from_sample_sheet($cfg_hash,$show_analyses,$log_file);
			}
			#Given a analysis name or a analysis id returns all the sample ids
			if ( $show_samples ne '' ){
				show_samples_from_group($cfg_hash,$show_samples,$log_file);
			}			
			#Given a set of analysis names gives all the analysis ids
			if ( $get_analyses_ids ne '' ){
				show_analysesids_from_analysesnames($cfg_hash,$get_analyses_ids,$log_file);
			}	
			#Given a set of sample names gives all the samples ids
			if ( $get_samples_ids ne '' ){
				show_samplesids_from_samplesnames($cfg_hash,$get_samples_ids,$log_file);
			}						
			#Given analysis ids stores and compresses the output files into the storage area
			if ( $vg_store_analyses ne '' ){
				exec_store_and_compress_analyses($cfg_hash,$store_analyses,$force_run,$log_file);
			}
			#Given analysis ids stores and compresses the output files into the storage area
			if ( $store_analyses ne '' ){
				store_and_compress_analyses($cfg_hash,$working_folder,$store_analyses,$force_run,$log_file);
			}			
			#Given analysis ids stores the output files into the storage area and makes a symbolic link
			if ( $save_analyses ne '' ){
				save_analyses($cfg_hash,$working_folder,$store_analyses,$force_run,$log_file);
			}					
			#Clean all the folders refine,varcall,genotype completely and
			#the useless bam and sam in alignment
			if ( $clean_folder ne '' ){
				my $sep = ",";
				my @analyses = split($sep,separate_input_ids($clean_folder,$sep));
				clean_folders($cfg_hash,\@analyses,$working_folder,$log_file);
			}

	}else { 
			
			preliminary_configuration();
			#print_and_log("The selected pipeline is $pipeline...\n",$log_file);
			if (scalar ( @tasks_to_exec) > 0 ) {
					if ( $pipeline eq 'pipeline'){
						print_and_log("Pipeline execution...\n",$log_file);
						pipeline_execution();
					}
			}
	}
}else{
	log_and_exit("Cannot open configuration file. Please give a correct parameter for user_config! Exiting...\n",$log_file);
}

#We also remove all the other remaining log file in the folder created by user's interruptions
system("rm -f $working_folder/*.log"); 


#####################################################################################################
																																																		#
#####################################################################################################

=head2 make_analysis_folders
 Title   : make_analysis_folders
 Usage   : make_analysis_folders(   );

 Function:  Creates all the folders needed for the analyses. Needs the config
						has been generated.
 
 Returns : nothing

=cut
sub make_analysis_folders{
	$cfg_hash->{'all_analysis_ids'} = "";
	
	#For each analysis generates a folder with all the analyses
	foreach my $analysis_id (keys %$exec_hash ){
		
		#Fill a variable in cfg_hash with all the analysis ids
		$cfg_hash->{'all_analysis_ids'} .= "$analysis_id,";
		#print "Creating folders for analysis: $analysis_id\n";#DEBUGCODE
		#Obtain the analysis name from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
		#Create the analysis folder
		$cfg_hash->{$analysis_id."_f"} = $working_folder."/".$analysis_name;
		#Check if directory exists, otherwise it creates it
		unless(-d $cfg_hash->{$analysis_id."_f"} ){
			print_and_log("Creating folder ".$cfg_hash->{$analysis_id."_f"}."...\n",$log_file);
			mkdir $cfg_hash->{$analysis_id."_f"} or die "ERROR: can't create folder ".$cfg_hash->{$analysis_id."_f"}.". Check permissions. \n";
		 }
		 
		#data folder will go in the analysis folder	
		my $dataFolder = $cfg_hash->{$analysis_id."_f"}."/".$cfg_hash->{'data_fold'};
		$cfg_hash->{$analysis_id."_data_fold"} = $dataFolder;
		#Check if directory exists, otherwise it creates it
		unless(-d $dataFolder){
			print_and_log("Creating folder $dataFolder...\n",$log_file);
			mkdir $dataFolder or die "ERROR: can't create folder $dataFolder. Check permissions. \n";
		}				 
		
		#Set a file name to be used to list the links to output files
		$cfg_hash->{$analysis_id."_outlist_file"} = $dataFolder."/".$analysis_name."_".$cfg_hash->{'outlist_file'};
		#Set a file name to be used to enrich the output page with statistics obtained during the computation
		$cfg_hash->{$analysis_id."_outstats_file"} = $dataFolder."/".$analysis_name."_".$cfg_hash->{'outstats_file'};
		
		#print_and_log("Setting paths for outlist and outstats files for analysis $analysis_id:\n",$log_file);#DEBUGCODE
		#print_and_log("Outlist: ".$cfg_hash->{$analysis_id."_outlist_file"}."\n",$log_file);#DEBUGCODE
		#print_and_log("Outstats: ".$cfg_hash->{$analysis_id."_outstats_file"}."\n",$log_file);#DEBUGCODE
		
		
		#For each task generates a folder if it does not exist
		foreach my $task (@pipeline_tasks){
			#print "Creating the folder for $analysis_id and task: $task\n";#DEBUGCODE
			my $outFolder = $cfg_hash->{$analysis_id."_f"}."/".$cfg_hash->{$task.'_out_f'};
			$cfg_hash->{$analysis_id.'_'.$task.'_out_f'} = $outFolder;
			#Check if directory exists, otherwise it creates it
			unless(-d $outFolder){
				print_and_log("Creating folder $outFolder...\n",$log_file);
				mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
			 }
			#Log folder will go in the output folder	
			my $logFolder = $outFolder."/".$cfg_hash->{'log_fold'};
			$cfg_hash->{$analysis_id.'_'.$task.'_log_fold'} = $logFolder;
			#Check if directory exists, otherwise it creates it
			unless(-d $logFolder){
				print_and_log("Creating folder $logFolder...\n",$log_file);
				mkdir $logFolder or die "ERROR: can't create folder $logFolder. Check permissions. \n";
			}
			
			#Only into the final output folder, make also the folder for images
			if ( $task eq 'finalout'){
				#Log folder will go in the output folder	
				my $imgFolder = $outFolder."/".$cfg_hash->{'img_fold'};
				$cfg_hash->{$analysis_id.'_'.$task.'_img_fold'} = $imgFolder;
				#Check if directory exists, otherwise it creates it
				unless(-d $imgFolder){
					print_and_log("Creating folder $imgFolder...\n",$log_file);
					mkdir $imgFolder or die "ERROR: can't create folder $imgFolder. Check permissions. \n";
				}
				
				#CNV analysis folder
				my $cnvFolder = $outFolder."/".$cfg_hash->{'cnv_out_f'};
				$cfg_hash->{$analysis_id.'_cnv_out_f'} = $cnvFolder;
				#Check if directory exists, otherwise it creates it
				unless(-d $cnvFolder){
					print_and_log("Creating folder $cnvFolder...\n",$log_file);
					mkdir $cnvFolder or die "ERROR: can't create folder $cnvFolder. Check permissions. \n";
				}	

				#CNV analysis folder: results for SEX chromosomes
				my $sexFolder = $cnvFolder."/".$cfg_hash->{'cnv_sex_f'};
				#Check if directory exists, otherwise it creates it
				unless(-d $sexFolder){
					print_and_log("Creating folder $sexFolder...\n",$log_file);
					mkdir $sexFolder or die "ERROR: can't create folder $sexFolder. Check permissions. \n";
				}
				#CNV analysis folder: results for AUTOsomes
				my $autoFolder = $cnvFolder."/".$cfg_hash->{'cnv_auto_f'};
				#Check if directory exists, otherwise it creates it
				unless(-d $autoFolder){
					print_and_log("Creating folder $autoFolder...\n",$log_file);
					mkdir $autoFolder or die "ERROR: can't create folder $autoFolder. Check permissions. \n";
				}				
								
				#CNV analysis Log folder will go in the CNV folder	
				$logFolder = $cnvFolder."/".$cfg_hash->{'log_fold'};
				$cfg_hash->{$analysis_id.'_cnv_log_fold'} = $logFolder;
				#Check if directory exists, otherwise it creates it
				unless(-d $logFolder){
					print_and_log("Creating folder $logFolder...\n",$log_file);
					mkdir $logFolder or die "ERROR: can't create folder $logFolder. Check permissions. \n";
				}							
			}		

		}
	}
	chop($cfg_hash->{'all_analysis_ids'});
	
}


=head2 get_modified_config_hash
 Title   : get_modified_config_hash
 Usage   : get_modified_config_hash(   );

 Function:  
 		#Depending by the parameters something in the configuration may be changed
		#The cfg_hash will be overwritten
		
 Returns : returns the hash updated depending by the different situations
						Sequencing type, research group particular interest, generic about configuration
						
						at the end of the process if requested, overwrites the variables that
						user wanted to change by using the command line

=cut
sub get_modified_config_hash{
	  my $cfg_hash = shift;
	  my $analysis_id = shift;
	  my ($cfg_hash_updated) = shift;
	  #my $overwrite_conf = shift;
	  
	  #Copy the main in the updated hash so that it has all the initial values
	  foreach my $key (keys %$cfg_hash){
				$$cfg_hash_updated->{$key} = $cfg_hash->{$key};
		}
		
		if ( $analysis_id eq 'dummy'){
			#my $dummy_cfg_f = $config_folder."/".$cfg_hash->{'dummy_profile'};
			#configFile2Hash($dummy_cfg_f,\$$cfg_hash_updated);
			#Add the dummy profile with the specific time and resources requests for the jobcleaner
			if ( $profile ne '' ){
			$profile .= ",".$cfg_hash->{'dummy_profile'};
			}else{
			$profile = $cfg_hash->{'dummy_profile'};	
			}
			#Get the configuration for the jobcleaner from input
			my @profiles = split(",",$profile); 
			foreach my $profile (@profiles){
				my $new_config = $config_folder."/".$profile;
				if ( -e $new_config ){
					print_and_log("Changing the configuration hash using the specific profile at $new_config..\n",$log_file);
					configFile2Hash($new_config,\$$cfg_hash_updated);					
				}else{
					print_and_log("WARNING: Cannot open non existent profile $new_config\n",$log_file);					
				}					
			}			
		}else{
			
			#Research group specific requests
			my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
			if ( correct_type($userid,"positiveint")){
				my @pieces = split(/\./, $cfg_hash->{'profile_file_suff'});
				my $param_cfg_f = $config_folder."/".$pieces[0]."_$userid.".$pieces[1];
				#The specific configuration file should exist or nothing will change
				if ( -e $param_cfg_f){
					#print_and_log("Changing the configuration hash for User id: $userid using $param_cfg_f..\n",$log_file);#DEBUGCODE
					configFile2Hash($param_cfg_f,\$$cfg_hash_updated);					
				}else{
					print_and_log("WARNING: Cannot find a profile for user id: $userid using $param_cfg_f.\n",$log_file);					
				}
			}
			
			#Other specific profiles that the user wants to load using the input variable --profile
			if ( $profile ne '' ){
				#They could be separated by comma
				my @profiles = split(",",$profile); 
				foreach my $profile (@profiles){
					my $new_config = $config_folder."/".$profile;
					if ( -e $new_config ){
						print_and_log("Changing the configuration hash using the specific profile at $new_config..\n",$log_file);
						configFile2Hash($new_config,\$$cfg_hash_updated);					
					}else{
						print_and_log("WARNING: Cannot open non existent profile $new_config\n",$log_file);					
					}					
				}
			}
			
			####CONFIGURATION ADJUSTMENT
			#Depending by the parameters something in the configuration may be changed
			#The cfg_hash will be overwritten
			#Here we change configuration depending by the run per chromosome or per sample
			my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
			if ($perchrom == 1){
				my $mem_cfg_f = $config_folder."/".$cfg_hash->{'perchrom_profile'};
				configFile2Hash($mem_cfg_f,\$$cfg_hash_updated);
			}else{
				my $mem_cfg_f = $config_folder."/".$cfg_hash->{'persample_profile'};
				configFile2Hash($mem_cfg_f,\$$cfg_hash_updated);
			}
			
			###Sequencing Types
			#Right now we have only exomes and amplicons panels sequencing
			my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);		
			#Overwrite the config hash with parameters specific for Haloplex Technology
			if ($sequencingtype eq 'targeted'){
				print_and_log("Changing the configuration hash for sequencing type: $sequencingtype..\n",$log_file);
				my $param_cfg_f = $config_folder."/".$cfg_hash->{'targeted_profile'};
				configFile2Hash($param_cfg_f,\$$cfg_hash_updated);
			}
			elsif ($sequencingtype eq 'rna'){
				print_and_log("Changing the configuration hash for sequencing type: $sequencingtype..\n",$log_file);
				my $param_cfg_f = $config_folder."/".$cfg_hash->{'rna_profile'};
				configFile2Hash($param_cfg_f,\$$cfg_hash_updated);
			}elsif ($sequencingtype eq 'wgs'){
				print_and_log("Changing the configuration hash for sequencing type: $sequencingtype..\n",$log_file);
				my $param_cfg_f = $config_folder."/".$cfg_hash->{'wgs_profile'};
				configFile2Hash($param_cfg_f,\$$cfg_hash_updated);
			}			
			#print "Hash obtained after configuration modification: \n";
			#print Dumper\$$cfg_hash_updated;
		}
		
		#Overwrite the configuration using input command
		if ($overwrite_conf ne ''){
			print_and_log("Changing the configuration hash using inputed variables..\n",$log_file);
				my $sep = $cfg_hash->{'mult_ann_sep'};
				my @vars_to_ow = split(/\Q$sep\E/,$overwrite_conf);
				
				foreach my $var_to_ow (@vars_to_ow){
					$var_to_ow =~ s/\s//;
					my @parts = split("=",$var_to_ow);
					if ( defined $$cfg_hash_updated->{$parts[0]} ){
						if ( defined $parts[1] ){
							$$cfg_hash_updated->{$parts[0]} = $parts[1];
							print_and_log($parts[0]." = ".$parts[1]."...\n",$log_file);
						}else{
								delete($$cfg_hash_updated->{$parts[0]});
						}
					}else{
						print_and_log("Parameter ".$parts[0]." does not exist...\n",$log_file);
					}
				}
		}

	print_and_log("Finished modifying configuration..\n",$log_file);#DEBUGCODE
	#Set additional variables of the config hash now with the new params
	set_other_global_var_into_hash($$cfg_hash_updated);
	#After having modified the configuration check it
	check_config_parameters($$cfg_hash_updated,$log_file);
	#############
	
}


=head2 pipeline_execution
 Title   : pipeline_execution
 Usage   : pipeline_execution(   );

 Function:  Executes the pipeline
 
					N.B. THE CONFIG HASH HERE ONLY CHANGES TO cfg_hash_gr BECAUSE
					FOR EACH analysis THERE WILL BE A DIFFERENT ONE.
					PLEASE USE ALWAYS THAT AFTER THE LOOP ON analysis IDS
					
 Returns : nothing

=cut
sub pipeline_execution{
	#An hash that considers for each job its dependencies before than it can be run
	my $dependence_hash;
	my $dep_clean_hash;
	
	#Path to the programs runner
	my $program_call = $program_folder."/".$cfg_hash->{'programs_runner_script'};
	#Path to the job cleaner script
	my $jobclean_prog = $program_folder."/".$cfg_hash->{'jobclean_prog'};

	
	#Put variables in a configuration file. It will be used in each job
	#save_hash(\$cfg_hash,$config_file);

	
	
	#The following are tasks that must be executed with jobs but they are not analysis
	#tasks. Thus, if they are executed the pipeline tasks cannot be
	###N.B. FOR THESE JOBS THAT ARE RUN INDEPENDENTLY THE CONFIG FILE AT THE END OF THE 
	#       SUBROUTINE MUST NOT BE REMOVED!!!
	#####################DOWNLOAD THE DATASETS###################
	if ($dldatasets){
		
		#Creation of support tables used for the output generation 
		my $task = "dldatasets";
						
		print_and_log($niceClosure,$log_file);
		print_and_log("(".scalar(localtime)."): $task step \n",$log_file);
		
		#Checks the links containing information for gene annotation
		check_datasets_links($cfg_hash,$log_file);
		
		#Get dependencies
		my $dependency = "no_depend";
		print_and_log("Launching a job for $task \n",$log_file);
		#Call the program using qsub and giving it parameters for the job and
		my $job_name = $task;
		my $job_err = $main_log_folder."/".$task.'.err';
		my $job_log = $main_log_folder."/".$task.'.log';
		my $env_vars = "CONFIG_FILE=$config_file,TASK=$task,SAMPLE_ID=-1,ANALYSIS_ID=-1,LOG_FILE=$log_file,PIPE=$pipeline ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file);
		print_and_log("Qsub id: $job_id\n",$log_file);
		##Put this job id in the hash for dependencies
		$dependence_hash->{$job_name} = $job_id;		
		print_and_log($niceClosure,$log_file);
	}	
	##################################################
	
	#####################GENE TABLES###################
	if ($genedb){
		
		#Creation of support tables used for the output generation 
		my $task = "genedb";
						
		print_and_log($niceClosure,$log_file);
		print_and_log("(".scalar(localtime)."): $task step \n",$log_file);
		
		#Checks the links containing information for gene annotation
		check_genedb_links($cfg_hash,$log_file);
		
		#Get dependencies
		my $dependency = "no_depend";
		print_and_log("Launching a job for $task \n",$log_file);
		#Call the program using qsub and giving it parameters for the job and
		my $job_name = $task;
		my $job_err = $main_log_folder."/".$task.'.err';
		my $job_log = $main_log_folder."/".$task.'.log';
		my $env_vars = "CONFIG_FILE=$config_file,TASK=$task,SAMPLE_ID=-1,ANALYSIS_ID=-1,LOG_FILE=$log_file,PIPE=$pipeline ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file);
		print_and_log("Qsub id: $job_id\n",$log_file);
		##Put this job id in the hash for dependencies
		$dependence_hash->{$job_name} = $job_id;		
		print_and_log($niceClosure,$log_file);
	}	
	##################################################

	#####################GENE ANNOTATION###################
	if ($geneann){
		
		#Creation of support tables used for the output generation 
		my $task = "geneann";
						
		print_and_log($niceClosure,$log_file);
		print_and_log("(".scalar(localtime)."): $task step \n",$log_file);
		
		#Get dependencies
		my $dependency = "no_depend";
		print_and_log("Launching a job for $task \n",$log_file);
		#Call the program using qsub and giving it parameters for the job and
		my $job_name = $task;
		my $job_log = $main_log_folder."/".$task.'.log';
		my $env_vars = "CONFIG_FILE=$config_file,TASK=$task,SAMPLE_ID=-1,ANALYSIS_ID=-1,LOG_FILE=$log_file,PIPE=$pipeline ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_log,$log_file);
		print_and_log("Qsub id: $job_id\n",$log_file);
		##Put this job id in the hash for dependencies
		$dependence_hash->{$job_name} = $job_id;		
		print_and_log($niceClosure,$log_file);
	}	
	##################################################
		
	#####################REANNOTATION###################
	if ($reannote){
		#Creation of support tables used for the output generation 
		my $task = "reannote";
				
		print_and_log($niceClosure,$log_file);
		print_and_log("(".scalar(localtime)."): $task step \n",$log_file);

		#Get dependencies
		my $dependency = "no_depend";
		print_and_log("Launching a job for $task \n",$log_file);
		#Call the program using qsub and giving it parameters for the job and
		my $job_name = $task;
		my $job_log = $main_log_folder."/".$task.'.log';
		my $env_vars = "CONFIG_FILE=$config_file,TASK=$task,SAMPLE_ID=-1,ANALYSIS_ID=-1,LOG_FILE=$log_file,PIPE=$pipeline ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_log,$log_file);
		print_and_log("Qsub id: $job_id\n",$log_file);
		#Put this job id in the hash for dependencies
		$dependence_hash->{$job_name} = $job_id;		
		print_and_log($niceClosure,$log_file);
	}	
	##################################################
	

	#####################RECOMPUTE FREQUENCIES###################
	if ($update_freqs){
		my $task = "update_freqs";


		print_and_log($niceClosure,$log_file);
		print_and_log("(".scalar(localtime)."): $task step \n",$log_file);

		#Get dependencies
		my $dependency = "no_depend";
		print_and_log("Launching a job for $task \n",$log_file);
		#Call the program using qsub and giving it parameters for the job and
		my $job_name = $task;
		my $job_log = $main_log_folder."/".$task.'.log';
		$cfg_hash->{'qsub_walltime'} = "100:00:00";
		my $env_vars = "CONFIG_FILE=$config_file,TASK=$task,SAMPLE_ID=-1,ANALYSIS_ID=-1,LOG_FILE=$log_file,PIPE=$pipeline ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_log,$log_file);
		print_and_log("Qsub id: $job_id\n",$log_file);
		#Put this job id in the hash for dependencies
		$dependence_hash->{$job_name} = $job_id;		
		print_and_log($niceClosure,$log_file);
	}	
	##################################################

	#N.B. THE CONFIG HASH HERE CHANGES TO cfg_hash_gr AND THE LOG FILE
	#CHANGE TO log_file_gr BECAUSE
	#FOR EACH GROUP THERE WILL BE A DIFFERENT ONE.
	#PLEASE USE ALWAYS THAT AFTER THE LOOP ON GROUP IDS	
					
	#An array to know which folder are to clean
	my @tasks_for_clean = ();
	
	#########################
	##########Start the Analyses
	#########################
	my $tot_analyses = scalar ( keys %{$exec_hash});
	my $num_analysis = 0;
	#Execute the tasks only if no other job is executed: update_freqs, gene_db, etc..
	if ( $do_pipeline_tasks ){
		#print_and_log("\n\nStarting the execution of the pipeline...\n",$log_file);#DEBUGCODE
		#For each analysis it executes the pipeline separately
		foreach my $analysis_id ( keys %{$exec_hash}){
			my $analysis_indication2 = "[an:$analysis_id ](".scalar(localtime)."): ";
			#print_and_log("Working with analysis: $analysis_id..\n",$log_file);
			
			##CHECKS AND CONFIGURATION PER-ANALYSIS
			#Obtain an hash specific for the analysis we are executing <---N.B.!!	
			#print_and_log("[an:$analysis_id ](".scalar(localtime).") Changing the configuration hash \n",$log_file);#DEBUGCODE
			my $cfg_hash_gr;
			get_modified_config_hash($cfg_hash,$analysis_id,\$cfg_hash_gr);
			
			##Parameters that are modified based on the current set of analyses
			$num_analysis++;
			if ( $num_analysis < $tot_analyses){
				  #If there is more than one analysis in progress,
					#the update of frequencies must be performed only at the last analysis
					$cfg_hash_gr->{'update_frequencies'} = "NO";
			}
			
			#Generate a log file to be used per group that will go in the data folder
			my $log_file_gr = $cfg_hash->{$analysis_id."_data_fold"}."/$analysis_id\_$time.log";
			#print_and_log("Copy $log_file to $log_file_gr ..\n",$log_file);#DEBUGCODE
			#Copy the content of the main LOG file to the group log file
			append_file_2_file($log_file,$log_file_gr);
			#Here I remove the main log file because there are the group ones
			#and jobclean wil use its own
			delete_file($log_file);
			
		 ###PEDIGREE FILE
			#################DB
			#Get the count of distinct samples where the genotype has been predicted (must be in use for frequency calculation and for this type of sequencing)
			my $query_num_samples =  "SELECT ".$cfg_hash->{'db_sample_id'}." FROM ".$cfg_hash->{'db_sample_table'}.
													 " WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id";
			#Connect to database
			my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
			# PERL DBI CONNECT
			my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
													 
			my $num_samples = get_count_from_selected($dbh,$query_num_samples," ",$log_file);	
			#Disconnect db
			$dbh->disconnect();
			#################
			#Obtain the analysis name from the database given the analysis id
			my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
																$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
																$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
			#Set the path of the file   
			my $out_ped = $cfg_hash_gr->{$analysis_id."_data_fold"}."/$analysis_name.".$cfg_hash->{'ped_ext'};	
			if ( $cfg_hash_gr->{'get_ped_file'} eq 'YES' ){
				#Verifying if your are doing a joint analysis of many samples		
				if ( $num_samples > $cfg_hash_gr->{'min_samples_4_joint'}){
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Obtaining a PED file for $num_samples samples : $out_ped \n",$log_file);
					get_ped_file_for_enlarged_analysis($cfg_hash_gr,$analysis_id,$out_ped,$log_file);					
					
					#Change the time needed for the final steps to the maximum
					$cfg_hash_gr->{'genot_walltime'} = $cfg_hash_gr->{'qsub_walltime_max'};
					$cfg_hash_gr->{'varfilt_walltime'} = $cfg_hash_gr->{'qsub_walltime_max'};
					$cfg_hash_gr->{'finalout_walltime'} = $cfg_hash_gr->{'qsub_walltime_max'};
					
					#Here we set perchrom=1 because if the analysis is a JOINT one, the user is constrict to run it per chromosome
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Since you are using $num_samples  samples, the database will be updated with: ".$cfg_hash->{'db_perchrom'}."=1  \n",$log_file_gr);
					my $fields = $cfg_hash_gr->{'db_perchrom'};
					my $values = "1";
					update_table_2($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
					$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$fields,$values,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);							
				}
				#I create here the pedfile, if there is a trio (M,F,P) then also the phasing can be done
					#unless the user decides to avoid it
				else{
					my $trios = get_ped_file($cfg_hash_gr,$analysis_id,$out_ped,$log_file_gr);
					#print_and_log("Created the Pedigree file $out_ped..\n",$log_file_gr);#DEBUGCODE
					if( $cfg_hash_gr->{'phasing'} eq 'YES' and $trios == 2){
						print_and_log("[an:$analysis_id ](".scalar(localtime).") The analysis is of a TRIO. Phasing will be performed \n",$log_file_gr);
						$phasing = 1;
					}else{
						print_and_log("[an:$analysis_id ](".scalar(localtime).") The analysis is not of a TRIO. Phasing will not be performed \n",$log_file_gr);
						$phasing = 0;	
					}									
				}				
			}


			####Management of operation with a target file.
			#Avoid to do this if the sequencing type is WGS
			my $sequencingtype = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_sequencingtype'},
									$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
									
			if ( $sequencingtype ne 'wgs'){
								
				####TARGET FILE
				#Check if the target file is present in the targets_table
				my $target_bed = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_targetbed'},
									$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
				my $target_bed_path = $cfg_hash_gr->{'target_reg_f'}."/".$target_bed;
				#Check if it exists								
				if (file_not_present($target_bed_path) > 0 ){ print_and_log("[an:$analysis_id ](".scalar(localtime).") WARNING: $target_bed_path does not exist. Please put the target file into ".$cfg_hash_gr->{'target_reg_f'}." ...\n",$log_file_gr);}		
				
				##TARGET EXTENSION
				#If the user wants to use the extended target, in this step we call a function to create a new target
				#the database is also updated
				my $target_extended = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_targetextended'},
						$cfg_hash_gr->{'db_analysis_id'},$analysis_id); 
				if ( $cfg_hash_gr->{'target_extension'} eq 'YES'  ){
						if ($target_extended == 0){
							my $target_bed_ext_path = extract_name($target_bed_path,"noext")."_plus".$cfg_hash_gr->{'target_extens_dim'}.".".$cfg_hash_gr->{'bed_ext'};
							if (!(-e $target_bed_ext_path) or (-z $target_bed_ext_path) ) {
								print_and_log("[an:$analysis_id ](".scalar(localtime).") Generating the extended target into: $target_bed_ext_path  \n",$log_file_gr);
								print_and_log("[an:$analysis_id ](".scalar(localtime).") Extending ".$cfg_hash_gr->{'target_extens_dim'}."nt on both 3' an 5' \n",$log_file_gr);
								get_extended_regions($target_bed_path,$target_bed_ext_path,$cfg_hash_gr->{'target_extens_dim'});						
							}else{
								print_and_log("[an:$analysis_id ](".scalar(localtime).") Extended target file $target_bed_ext_path have been already produced\n",$log_file_gr);		
							}

							
							#Update table
							$target_bed = extract_name($target_bed_ext_path,0);
							$target_bed_path = $cfg_hash_gr->{'target_reg_f'}."/".$target_bed;
							print_and_log("[an:$analysis_id ](".scalar(localtime).") The database will be updated with: $target_bed  \n",$log_file_gr);
							my $fields = $cfg_hash_gr->{'db_targetbed'}.",".$cfg_hash_gr->{'db_targetextended'};
							my $values = "'$target_bed',1";
							update_table_2($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$fields,$values,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
							
						}else{
							#The target extended should already exist
							if (!(-e $target_bed_path) or (-z $target_bed_path) ) {
								print_and_log("[an:$analysis_id ](".scalar(localtime).") The extended target $target_bed_path does not exist. Please correct... \n",$log_file_gr);
								print_and_log("[an:$analysis_id ](".scalar(localtime).") Update the db with: 1. row target file 2. targetextended=0 \n",$log_file_gr);
							}
							print_and_log("[an:$analysis_id ](".scalar(localtime).") $target_bed has been already extended.\n",$log_file_gr);	
						}
				}
				
				#Getting the target bed file using the folder for targets and the name contained in the database
				my $target_id = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},								
												$cfg_hash_gr->{'db_user'},$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$cfg_hash_gr->{'db_targetid'},
												$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");
				#If it does not exist, insert
				if ($target_id < 0 ){
					 $target_id = insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_target_files_table'},$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");
				}								
						

				#PER-CHROMOSOME: Any analysis could be executed per-chromosome, hence we generate the bed files for each cromosome
				#Get the perchrom flag
				my $perchrom = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_perchrom'},
									$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
				#Get the flag for status
				my $db_target_perchrdiv = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
								$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$cfg_hash_gr->{'db_target_perchrdiv'},
								$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");
				my $dest_folder = $main_data_folder."/".$cfg_hash_gr->{'target_perchr_fold'}."/$target_bed\_".$cfg_hash_gr->{'chrom_reg_f'};
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Target $target_bed per chromomsome folder: in $dest_folder (flag $db_target_perchrdiv)...\n",$log_file_gr);
				if ( $db_target_perchrdiv == 0 ){
					#If folder does not exists, creates it
					unless(-d $dest_folder){
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Creating folder $dest_folder...\n",$log_file_gr);
						mkdir $dest_folder or die "[an:$analysis_id ](".scalar(localtime).") ERROR: can't create folder $dest_folder. Check permissions. \n";
					}	
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Separating  $target_bed_path per chromomsome in $dest_folder...\n",$log_file_gr);
					
					separate_bed_perchr($target_bed_path,$dest_folder);
					
					#Additionally some chromosome can be further separated in two parts
					if ( length($cfg_hash_gr->{'chroms_to_split'}) > 0 ){
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Separating in two bed chr:".$cfg_hash_gr->{'chroms_to_split'}."...\n",$log_file_gr);

						my @chroms_to_split = split(",",$cfg_hash_gr->{'chroms_to_split'});
						foreach my $chrom_to_split (@chroms_to_split) {
							my $bed_in = $dest_folder."/".$cfg_hash->{'chr_suff'}."$chrom_to_split.".$cfg_hash->{'bed_ext'};
							my $bed_out1 = $dest_folder."/".$cfg_hash->{'chr_suff'}."$chrom_to_split\_1.".$cfg_hash->{'bed_ext'};
							my $bed_out2 = $dest_folder."/".$cfg_hash->{'chr_suff'}."$chrom_to_split\_2.".$cfg_hash->{'bed_ext'};
							#QUESTA FUNZIONE ANDRA FATTA IN MODO DA FARE UNO SPLIT IN N PARTI INVECE CHE SOLO 2
							split_bedfile ($bed_in,$bed_out1,$bed_out2);
						}
					}
					#Update the chromosome regions folder for this analysis
					$cfg_hash_gr->{'chrom_reg_f'} = $dest_folder;
					
					#Update table
					my $fields = $cfg_hash_gr->{'db_target_perchrdiv'}.";".$cfg_hash_gr->{'db_target_perchrdivdir'};
					my $values = "1;'$dest_folder'";
					update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
					$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$fields,$values,$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");
				}else{
					##Getting the target_bed
					$dest_folder = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$cfg_hash_gr->{'db_target_perchrdivdir'},
						$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");	
				}		
				#Set the chromosome regions folder into the config hash
				$cfg_hash_gr->{'chrom_reg_f'} = $dest_folder;
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Chromosomal regions will be pickedfrom: ".$cfg_hash_gr->{'chrom_reg_f'}."..\n",$log_file_gr);	    
			
				#INTERVAL LIST FILE :Construct the target interval list file from the target file (This is needed ALWAYS in downstream analyses)
				my $target_int_list = $target_bed_path.".".$cfg_hash_gr->{'target_inter_list_ext'};
				if (file_not_present($target_int_list) > 0 ){
					my $prog_used = $cfg_hash_gr->{'bedtointerval_prog'};
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Executing $prog_used to generate the ".$cfg_hash_gr->{'target_inter_list_ext'}." file for $target_bed : ",$log_file_gr);
					run_PICARD_BedToIntervalList($cfg_hash_gr,$prog_used,$log_file_gr,$target_bed_path,$target_int_list);
					#Update table
					update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
					$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$cfg_hash_gr->{'db_target_bedtointerval'},"1",$cfg_hash_gr->{'targetname'},"'".$target_bed."'");
				}	

				##TARGET BED FOR CNV ANALYSIS
				if ( $cnv_detection ){
					print_and_log("[an:$analysis_id ](".scalar(localtime).") CNV detection pipeline active...\n",$log_file_gr);
					
					#SEX AND AUTO CHROMOSOMES TARGET BED FILE
					#For all programs for CNV based on depth
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Generating a target file only for autosomes and one for sex chromosomes...\n",$log_file_gr);
					#Give a name to the target with only sex chromosomes
					$cfg_hash_gr->{'target_bed_sex'} = extract_name($target_bed_path,"noext")."_".$cfg_hash_gr->{'cnv_sex_f'}.".".$cfg_hash_gr->{'bed_ext'};
					#and to that with only autosomes
					$cfg_hash_gr->{'target_bed_auto'} = extract_name($target_bed_path,"noext")."_".$cfg_hash_gr->{'cnv_auto_f'}.".".$cfg_hash_gr->{'bed_ext'}; 
					
					#If the target has not been extended, clean it
					my $target_bed_path_temp = "";
					#if ( $target_extended == 0 ) {
						#Leave only one description on the last column of the target
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Removing descriptions from $target_bed_path into $target_bed_path.temp...\n",$log_file_gr);			
						if ( get_only_first_on_col_from_BED($target_bed_path,$target_bed_path.".temp",4) > 0 ){
							#It happens that the BED target file has descriptions empty with '-'. Which can cause problem during
							#the process. Hence, change '-' to NONE
							my $command = "sed 's/-/NONE/' $target_bed_path.temp > $target_bed_path.temp2";
							print_and_log("[an:$analysis_id ](".scalar(localtime).") Executing $command...\n",$log_file_gr);				
							(system( $command )) == 0 or die "Unable to execute: $command!\n";
							$target_bed_path_temp = "$target_bed_path.temp2";					
						}else{
							$target_bed_path_temp =  $target_bed_path;
						}
					#}else{
						#$target_bed_path_temp =  $target_bed_path;
					#}			
					#
					if ( (!(-e $cfg_hash_gr->{'target_bed_sex'}) or (-z $cfg_hash_gr->{'target_bed_sex'})) or (!(-e $cfg_hash_gr->{'target_bed_auto'}) or (-z $cfg_hash_gr->{'target_bed_auto'})) ){
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Separating SEX and AUTO chromosomes in two target files...\n",$log_file_gr);
						separate_bed_auto_and_sex_chr($target_bed_path_temp,$cfg_hash_gr->{'target_bed_sex'},$cfg_hash_gr->{'target_bed_auto'});
					}else{
						print_and_log("[an:$analysis_id ](".scalar(localtime).") SEX and AUTO chromosomes target files for $target_bed have been already produced\n",$log_file_gr);		
					}
					
					if ($cfg_hash_gr->{'remove_temp'} eq 'YES' ){
						delete_file("$target_bed_path.temp");
						delete_file("$target_bed_path.temp2");					
					}

					
					#Get flag to see if the gc_content file was already produced for this target
					my $db_target_gccontentf = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_target_files_table'},$cfg_hash_gr->{'db_target_gccontentf'},
									$cfg_hash_gr->{'db_targetname'},"'".$target_bed."'");
					if ( $db_target_gccontentf == 0 ){
						#XHMM
						#Calculate GC content must be performed only once per-target
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Generating a gc content reference file for $target_bed (this operation is ".
						" performed only once per-each target)\n",$log_file_gr);				
						calculate_gc_content($cfg_hash_gr,$target_bed,$log_file);					
					}else{
						print_and_log("[an:$analysis_id ](".scalar(localtime).") GC content reference file for $target_bed has been already produced\n",$log_file_gr);		
					}
				}
				
				##TARGET BED for EXON COVERAGE
				#Prepare Target Exons:
				my $ucsc_exons = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_exons_bed'};
				my $target_exons = extract_name($target_bed_path,"noext")."_".$cfg_hash->{'ucsc_exons_bed'};
				
				if ( !(-e $target_exons) or (-z $target_exons) ){
					print_and_log( "[an:$analysis_id ](".scalar(localtime).") Getting the exons on the target...",$log_file);#DEBUGCODE

					#1. Intersect the exons with the target so that we have the exons on-target
					#to see how much is the overlap and also when overlap is missing	
					#Parameters used are:
					 #  -a the exons
					  # -b the target file
						# -wa which allows to write the original entry in A . 
					my $inters_params = " -wa ";
					print_and_log( "...intersecting UCSC exons with target (intersectBed)...",$log_file);#DEBUGCODE
					run_BEDTOOLS_intersectBed($cfg_hash,$ucsc_exons,$target_bed_path,$target_exons.".temp",$inters_params,$log_file);#######	
					
					#Get only the chrom and interval, sort and get uniq
					print_and_log( "Get only the chrom and interval, sort and get unique lines...",$log_file);#DEBUGCODE
					my $command =	"  cut -f1-3 $target_exons.temp | sort -k1,1 -k2,2n | uniq > $target_exons.temp2";
					print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
					try_exec_command($command) or die "Unable to execute command: $command\n";		
					
					#Sort using bedtools and a reference chromosome order
					my $chr_order_f = "chr_order.txt";
					my $sortparams = " -faidx $main_data_folder/$chr_order_f ";
					run_BEDTOOLS_sortBed($cfg_hash,"$target_exons.temp2",$sortparams,$target_exons,$log_file);
						
					if ($cfg_hash_gr->{'remove_temp'} eq 'YES' ){
						delete_file("$target_exons.temp");
						delete_file("$target_exons.temp2");
					}				
				}else{
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Target BED using exons $ucsc_exons has been already produced\n",$log_file_gr);		
				}	
			}else{
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Analysis on WGS data. Target bed and accessories will not be used/produced\n",$log_file_gr);		
			}	
			
			
			#Store the config hash in the data folder for the analysis
			my $cfg_hash_f = $cfg_hash->{$analysis_id."_data_fold"}."/cfg_$analysis_id.hash";
			print_and_log("[an:$analysis_id ](".scalar(localtime).") Storing the configuration hash in $cfg_hash_f..\n",$log_file_gr);
			$cfg_hash_gr->{'cfg_hash_f'} = $cfg_hash_f;
			save_hash(\$cfg_hash_gr,$cfg_hash_f);

      ################################
      #PIPELINE EXECUTION FOR EACH SAMPLE
      ###################################
      
      #For each sample in each analysis
			while (my ($sample_id, $value) = each %{ $exec_hash->{$analysis_id} } ) {
				#print_and_log($niceClosure,$log_file_gr);
				#print_and_log("\nStarting the pipeline for analysis: $analysis_id and sample $sample_id...\n",$log_file_gr);#DEBUGCODE
				#print_and_log("Files to analyze are indexed with: ".$exec_hash->{$analysis_id}->{$sample_id}." \n",$log_file_gr);#DEBUGCODE
				
				#Get the samples to be analyzed
				my @readf_ids = ();
				#Do it with an array
				@readf_ids = split($cfg_hash_gr->{'parameters_sep'},$exec_hash->{$analysis_id}->{$sample_id});	
				
				#Quality Check includes the quality check and the trimming
				if ($qualityCheck){
					my $sample_info_hash;
					my $task = "qc";
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					print_and_log($niceClosure,$log_file_gr);
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id ](".scalar(localtime)."): ";
					#print_and_log("[an:$analysis_id ](".scalar(localtime)."): $task step \n",$log_file_gr);#DEBUGCODE

					#Get dependencies
					my $dependency = "no_depend";
					#Loop through sample ids so that it will run a job for each sample
					foreach my $readf_id (@readf_ids){
						$analysis_indication2 ="[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task for read file $readf_id : ",$log_file_gr);

						#Getting the read file name
						my $readf_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_readf_table'},$cfg_hash_gr->{'db_readf_name'},
						$cfg_hash_gr->{'db_readf_id'},$readf_id);	
						
						#Call the program using qsub and giving it parameters for the job and
						my $job_name = $task."_r".$readf_id;
						#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_".$cfg_hash_gr->{$task.'_log'};
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'err_ext'};						
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$readf_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."):  - step:$task - dependencies:$dependency - jobid:$job_id]\n",$log_file_gr);
						#Put this job id in the hash for dependencies
						$dependence_hash->{$job_name} = $job_id;		
						
					}
					print_and_log($niceClosure,$log_file_gr);
				}

				#Trimming includes the step of trimming of sequences
				if ($trimming){
					my $sample_info_hash;
					my $task = "trim";
				
					#Managing the dependencies
					my $task_to_wait = "qc";
										
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					print_and_log($niceClosure,$log_file_gr);
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id ](".scalar(localtime)."): ";
					#print_and_log("[an:$analysis_id ](".scalar(localtime)."): $task step \n",$log_file_gr);#DEBUGCODE

					#Loop through sample ids so that it will run a job for each sample
					foreach my $readf_id (@readf_ids){
						$analysis_indication2 ="[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task for read file $readf_id : ",$log_file_gr);
						#Get dependencies
						my $dependency = "no_depend";
						if (defined $dependence_hash->{$task_to_wait."_r".$readf_id}){
							$dependency = $dependence_hash->{$task_to_wait."_r".$readf_id};	
						}
					
						#Getting the read file name
						my $readf_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_readf_table'},$cfg_hash_gr->{'db_readf_name'},
						$cfg_hash_gr->{'db_readf_id'},$readf_id);	
						
						#Call the program using qsub and giving it parameters for the job and
						my $job_name = $task."_r".$readf_id;
						#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_".$cfg_hash_gr->{$task.'_log'};
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'err_ext'};						
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$readf_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id ](".scalar(localtime).") - step:$task - dependencies:$dependency - jobid:$job_id]\n",$log_file_gr);
						#Put this job id in the hash for dependencies
						$dependence_hash->{$job_name} = $job_id;		
						
					}
					print_and_log($niceClosure,$log_file_gr);
				}
				
								
				#Alignment against reference
				if ($alignment){

					my $task = "align";
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";				
					print_and_log($niceClosure,$log_file_gr);
					#print_and_log("[an:$analysis_id ](".scalar(localtime).") : $task step \n",$log_file_gr);#DEBUGCODE
				
					#Managing the dependencies
					my $task_to_wait = "trim";
					#Loop through sample ids so that it will run a job for each sample
					foreach my $readf_id (@readf_ids){
						$analysis_indication2 ="[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
						
						#Get dependencies
						my $dependency = "no_depend";
						if (defined $dependence_hash->{$task_to_wait."_r".$readf_id}){
							$dependency = $dependence_hash->{$task_to_wait."_r".$readf_id};	
						}
						#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
						
						#Getting the read file name
						my $readf_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_readf_table'},$cfg_hash_gr->{'db_readf_name'},
						$cfg_hash_gr->{'db_readf_id'},$readf_id);	
													
						#Call the program using qsub and giving it parameters for the job and
						my $job_name = $task."_r".$readf_id;
						#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_".$cfg_hash_gr->{$task.'_log'};
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'err_ext'};						
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$readf_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id] (".scalar(localtime)."): [sample:$sample_id - readfid:$readf_id - step:$task - dependencies:$dependency - jobid:$job_id]\n",$log_file_gr);
						
						#Put this job id in the hash for dependencies
						$dependence_hash->{$job_name} = $job_id;
						#The STATS step will use the overall sample aligned files. Hence here we set the 
						#waiting for the sample
						if ( ! defined $dependence_hash->{$task."_s".$sample_id}){
							$dependence_hash->{$task."_s".$sample_id} = $job_id.":";		 	
						}else{
							$dependence_hash->{$task."_s".$sample_id} .= $job_id.":";
						}
						#The GR_STATS step will use the overall group aligned files. Hence here we set the 
						#waiting for the group
						if ( ! defined $dependence_hash->{$task."_a".$analysis_id}){
							$dependence_hash->{$task."_a".$analysis_id} = $job_id.":";		 	
						}else{
							$dependence_hash->{$task."_a".$analysis_id} .= $job_id.":";
						}
						#Fill the hash for cleaning with this jobname and dependencies
						if ($dependency ne "no_depend"){
								$dep_clean_hash->{$job_name}->{'id'} = $job_id;
								$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
						}
					}
				}

						
				#refining alignments with Base recalibration
				if ($refining){
					my $main_task = "refine";
					#print "Values: RTC".$cfg_hash_gr->{'realign_target_creator'}." IR".$cfg_hash_gr->{'indel_realigner'}." BR ".$cfg_hash_gr->{'base_recalibrator'}." PR ".$cfg_hash_gr->{'print_reads'}."\n";
					if ($cfg_hash_gr->{'realign_target_creator'} eq 'YES'  or $cfg_hash_gr->{'indel_realigner'} eq 'YES'){
						
						#INDEL REALIGNMENT PER-LANE
						my $task = "IR";
						#Put in the config hash a variable with all the tasks to be executed
						if ( not grep {/\b$task\b/} @tasks_for_clean){
							push(@tasks_for_clean,$task);
						}
						my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";				
						print_and_log($niceClosure,$log_file_gr);
						#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step \n",$log_file_gr);#DEBUGCODE
						
						#Managing the dependencies
						my $task_to_wait = "align";
						#Set an id for the next job to wait
						my $next_job_to_wait = $task."_s".$sample_id;
						#Loop through sample ids so that it will run a job for each sample
						foreach my $readf_id (@readf_ids){
							$analysis_indication2 ="[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
							

							#Get dependencies
							my $dependency = "no_depend";
							if (defined $dependence_hash->{$task_to_wait."_r".$readf_id}){
								$dependency = $dependence_hash->{$task_to_wait."_r".$readf_id};	
							}
							#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
							print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
							
							#Getting the read file name
							my $readf_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_readf_table'},$cfg_hash_gr->{'db_readf_name'},
							$cfg_hash_gr->{'db_readf_id'},$readf_id);	
														
							#Call the program using qsub and giving it parameters for the job and
							my $job_name = $task."_r".$readf_id;
							#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$readf_name."_".$cfg_hash_gr->{$main_task.'_log'};
							my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'log_ext'};
							my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$readf_name."_$task.".$cfg_hash_gr->{'err_ext'};								
							my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$readf_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
							my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$main_task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
							#print_and_log("Qsub id: $job_id\n",$log_file_gr);
							print_and_log("[an:$analysis_id ](".scalar(localtime).") - step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

							#Put this job id in the hash for dependencies
							if ( ! defined $dependence_hash->{$next_job_to_wait}){
								$dependence_hash->{$next_job_to_wait} = $job_id.":";		 	
							}else{
								$dependence_hash->{$next_job_to_wait} .= $job_id.":";
							}
							#Fill the hash for cleaning with this jobname and dependencies
							if ($dependency ne "no_depend"){
										$dep_clean_hash->{$job_name}->{'id'} = $job_id;
										$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
								}
							}
							chop($dependence_hash->{$next_job_to_wait});
					}
					
					if ($cfg_hash_gr->{'base_recalibrator'} eq 'YES'  or $cfg_hash_gr->{'print_reads'} eq 'YES'){
						##############################
						#BASE RECALIBRATION PER SAMPLE
						##############################
						my $task = "BR";
						#Put in the config hash a variable with all the tasks to be executed
						if ( not grep {/\b$task\b/} @tasks_for_clean){
							push(@tasks_for_clean,$task);
						}
						my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";				
						print_and_log($niceClosure,$log_file_gr);
						#print_and_log("[an:$analysis_id ](".scalar(localtime).")  $task step \n",$log_file_gr);#DEBUGCODE
						
						#Managing the dependencies
						my $task_to_wait;
						#If the indel realigner has been run, wait for that..
						if ($cfg_hash_gr->{'realign_target_creator'} eq 'YES'  or $cfg_hash_gr->{'indel_realigner'} eq 'YES'){
							$task_to_wait  = "IR";
						}#Otherwise wait for the alignment
						else{
							$task_to_wait  = "align";
						}
						
						#The name here will be reduced since it must be shorter than 15 chars
						my $next_job_to_wait = $task."_s".$sample_id;

						#Get dependencies
						my $dependency = "no_depend";
						if (defined $dependence_hash->{$task_to_wait."_s".$sample_id}){
							$dependency = $dependence_hash->{$task_to_wait."_s".$sample_id};	
						}
						#print_and_log("[an:$analysis_id ](".scalar(localtime).")  This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
						print_and_log("[an:$analysis_id ](".scalar(localtime).")  Launching a job for $task : ",$log_file_gr);
						
						#Getting the sample name
						my $sample_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_sample_name'},
							$cfg_hash_gr->{'db_sample_id'},$sample_id);	
													
						#Call the program using qsub and giving it parameters for the job and
						my $job_name = $task."_s".$sample_id;
						#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$sample_name."_".$cfg_hash_gr->{$main_task.'_log'};
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$main_task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'err_ext'};							
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$sample_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$main_task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id ](".scalar(localtime).") - step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

						#Put this job id in the hash for dependencies for samples
						if ( ! defined $dependence_hash->{$next_job_to_wait}){
							$dependence_hash->{$next_job_to_wait} = $job_id.":";		 	
						}else{
							$dependence_hash->{$next_job_to_wait} .= $job_id.":";
						}
						#2019
						#The CNV step will use the files after base recalibration. Hence here we set the 
						#waiting for the analysis
						if ( ! defined $dependence_hash->{$task."_a".$analysis_id}){
							$dependence_hash->{$task."_a".$analysis_id} = $job_id.":";		 	
						}else{
							$dependence_hash->{$task."_a".$analysis_id} .= $job_id.":";
						}
						#Fill the hash for cleaning with this jobname and dependencies
						if ($dependency ne "no_depend"){
								$dep_clean_hash->{$job_name}->{'id'} = $job_id;
								$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
						}
						chop($dependence_hash->{$next_job_to_wait});
					}
				}
				
				#Merge read files in unique files for each sample. This task will be executed only if the refinement wasn't executed
				#Here I have commented the check on the input parameter and included the check
				#only using the parameter from the config file so that if kept, the pipeline will always 
				#run this step when the variable multiple or perchrom is 1 . In this way the user does not have to remember
				#to give the input parameter
				if ( $varcall and $cfg_hash_gr->{'print_reads'} eq 'NO' and
					($cfg_hash_gr->{'mrdup_groups'} eq 'YES' or $cfg_hash_gr->{'index_after_merge'} eq 'YES') ){
					my $task = $cfg_hash_gr->{'merge_task_s'};
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";	
					print_and_log($niceClosure,$log_file_gr);
					#print_and_log("[an:$analysis_id ](".scalar(localtime).")  $task step \n",$log_file_gr);#DEBUGCODE
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Checking if is useful to apply mergesam to your samples..\n",$log_file_gr);
					my $multiple = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_multiple'},
							$cfg_hash_gr->{'db_sample_id'},$sample_id);
					#PER-CHROMOSOME: If this analysis is executed per-chromosome, then we must create the bed files for each cromosome
					#Verifying that there has been an execution per chromosome
					my $perchrom = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
								$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_perchrom'},
								$cfg_hash_gr->{'db_analysis_id'},$analysis_id);					
					#Get sequencing type
					my $sequencingtype = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
										$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_sequencingtype'},
										$cfg_hash_gr->{'db_analysis_id'},$analysis_id);		

					#Verifying if the refinement was executed
					my $printreads = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'print_reads_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    					
					#Verify if there are multiple read files and refinement step was not executed launch the merge
					if ( $multiple == 1 and $printreads == 0){
						$mergesam = 1;#To activate the wait in variant calling step
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Sample $sample_id is subdivided in multiple read files. Merging results..\n",$log_file_gr);
		
						#Managing the dependencies
						my $task_to_wait;
						#If no refine step must be executed before, then the task to wait is the alignment
						if ($cfg_hash_gr->{'base_recalibrator'} eq 'NO' and $cfg_hash_gr->{'print_reads'} eq 'NO' and $sequencingtype eq 'targeted' ){
							$task_to_wait = "align";
						}elsif ($cfg_hash_gr->{'base_recalibrator'} eq 'YES'  or $cfg_hash_gr->{'print_reads'} eq 'YES'){
							$task_to_wait = "BR";
						}else{
							$task_to_wait = "align";	
						}
						#Get dependencies 					
						#The name here will be reduced since it must be shorter than 15 chars
						my $job_to_wait = $task_to_wait."_s".$sample_id;
						my $dependency = "no_depend";
						if (defined $dependence_hash->{$job_to_wait}){
							$dependency = $dependence_hash->{$job_to_wait};	
						}
						#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
						#Getting the sample name
						my $sample_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_sample_name'},
							$cfg_hash_gr->{'db_sample_id'},$sample_id);	
													
						#Call the program using qsub and giving it parameters for the job and
						#my $job_name = substr($task."_".$sample_id,0,$cfg_hash_gr->{'qsub_name_len'});
						my $job_name = $task."_s".$sample_id;
						#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_".$cfg_hash_gr->{$task.'_log'};
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'err_ext'};						
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$sample_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id ](".scalar(localtime).")  - step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

						#Put this job id in the hash for dependencies
						if ( ! defined $dependence_hash->{$job_name}){
							$dependence_hash->{$job_name} = $job_id.":";		 	
						}else{
							$dependence_hash->{$job_name} .= $job_id.":";
						}
						chop($dependence_hash->{$job_name});	#code
						#Fill the hash for cleaning with this jobname and dependencies
						if ($dependency ne "no_depend"){
								$dep_clean_hash->{$job_name}->{'id'} = $job_id;
								$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
						}
					}else{ $mergesam = "";}
				}#MERGEGROUPS END
					
				#VARIANT CALL STEP	
				if ( $varcall ){
					my $task = "varcall";
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					
					###For each analysis we check if there are multiple samples. It it is, a Joint Genotyping will be run###
					#Get the samples ids involved for the group. If they are more than one, do JG
					my $query = "SELECT ".$cfg_hash_gr->{'db_sample_id'}." FROM  ".$cfg_hash_gr->{'db_sample_table'}." WHERE ".$cfg_hash_gr->{'db_analysis_id'}."=$analysis_id;";	
					print_and_log( "Executing: $query\n",$log_file_gr);
					my $group_sam = do_query_select_all($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},$cfg_hash_gr->{'db_pass'},$query,$cfg_hash_gr->{'db_sample_id'});
					
					#GenotypeGVCF is executed for multiple samples and if the user chooses it for a single sample
					if (scalar(keys %{$group_sam}) > 1 or  $cfg_hash_gr->{'do_genotypeGVCF_4_single_sample'} eq 'YES'){
					#if (scalar(keys %{$group_sam}) > 0){
						print_and_log("Setting dogenot=1... (Joint Genotype will be executed after variant calling (GATK BP)) ...\n",$log_file_gr);
						#Update table
						update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_dogenot'},"1",$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
						#Activating the variable that runs the job
						$genot = 1;
					}else{
						print_and_log("Setting dogenot=0...(There is a single sample and GenotypeGVCF will not be run (GATK BP)) ...\n",$log_file_gr);
						#Update table
						update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_dogenot'},"0",$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
						#Activating the variable that runs the job
						$genot = 0;					
					}
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";						
					print_and_log($niceClosure,$log_file_gr);
					#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step \n",$log_file_gr);#DEBUGCODE		
					
					#Get sequencing type
					my $sequencingtype = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
										$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_sequencingtype'},
										$cfg_hash_gr->{'db_analysis_id'},$analysis_id);	
										
					#Managing the dependencies
					my $task_to_wait = "";			
					if ( $mergesam ) {
						$task_to_wait = $cfg_hash_gr->{'merge_task_s'};
					}else{
						#If no refine step must be executed before, then the task to wait is the alignment
						if (( 
							#$cfg_hash_gr->{'realign_target_creator'} eq 'NO' and $cfg_hash_gr->{'indel_realigner'} eq 'NO' and
								$cfg_hash_gr->{'base_recalibrator'} eq 'NO' and $cfg_hash_gr->{'print_reads'} eq 'NO'
								and $sequencingtype eq 'targeted' )){
							$task_to_wait = "align";
						}	
						#If base recalibration was executed than wait it 
						elsif ($cfg_hash_gr->{'base_recalibrator'} eq 'YES'  or $cfg_hash_gr->{'print_reads'} eq 'YES'){
							$task_to_wait = $cfg_hash_gr->{'baserecal_task_s'};
						}
						##Otherwise wait for the alignment
						else{
							$task_to_wait = "align";
						}
						
					}

					#Get dependencies
					my $dependency = "no_depend";
					my $job_to_wait = $task_to_wait."_s".$sample_id;
					if (defined $dependence_hash->{$job_to_wait}){
						$dependency = $dependence_hash->{$job_to_wait};	
					}
					#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : \n",$log_file_gr);

					#Getting the sample name
					my $sample_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_sample_name'},
							$cfg_hash_gr->{'db_sample_id'},$sample_id);	
										
					#Call the program using qsub and giving it parameters for the job and
					my $job_name = $task."_s".$sample_id;
					#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_".$cfg_hash_gr->{$task.'_log'};
					my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'log_ext'};
					my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'err_ext'};						
					my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$sample_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
					my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
					#print_and_log("Qsub id: $job_id\n",$log_file_gr);
					print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);
					#Put this job id in the hash for dependencies
					#The next step will use the overall analysis. Hence here we set the 
					#waiting for the analysis
					if ( ! defined $dependence_hash->{$task."_a".$analysis_id}){
						$dependence_hash->{$task."_a".$analysis_id} = $job_id.":";		 	
					}else{
						$dependence_hash->{$task."_a".$analysis_id} .= $job_id.":";
					}
					#Fill the hash for cleaning with this jobname and dependencies
					if ($dependency ne "no_depend"){
							$dep_clean_hash->{$job_name}->{'id'} = $job_id;
							$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
					}
				}#VARIANT CALL END	
				
				#Statistics on the aligned bam files
				if ( $stats ){
					my $task = "stats";
					#Put in the config hash a variable with all the tasks to be executed
					if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
					my $analysis_indication2 = "[an:$analysis_id sam: $sample_id](".scalar(localtime)."): ";				
					print_and_log($niceClosure,$log_file_gr);
					#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step \n",$log_file_gr);#DEBUGCODE				
				
					#Managing the dependencies
					#my $task_to_wait = "align";			
					#Since January 2019 the statistics are performed on the recalibrated BAM, hence
					#the task to wait is "BR"
					my $task_to_wait = "BR";
					
					#Get dependencies
					my $dependency = "no_depend";
					my $job_to_wait = $task_to_wait."_s".$sample_id;
					if (defined $dependence_hash->{$job_to_wait}){
						$dependency = $dependence_hash->{$job_to_wait};	
					}
					#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : \n",$log_file_gr);
				
					#Getting the sample name
					my $sample_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_sample_name'},
							$cfg_hash_gr->{'db_sample_id'},$sample_id);	
												
					#Call the program using qsub and giving it parameters for the job and
					my $job_name = $task."_s".$sample_id;
					#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_".$cfg_hash_gr->{$task.'_log'};
					my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'log_ext'};
					my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$sample_name."_$task.".$cfg_hash_gr->{'err_ext'};						
					my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$sample_id,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
					my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
					#print_and_log("Qsub id: $job_id\n",$log_file_gr);
					print_and_log("[an:$analysis_id ](".scalar(localtime).")  step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

					#Put this job id in the hash for dependencies
					#The next step will use the overall analysis. Hence here we set the 
					#waiting for the analysis
					if ( ! defined $dependence_hash->{$task."_a".$analysis_id}){
						$dependence_hash->{$task."_a".$analysis_id} = $job_id.":";		 	
					}else{
						$dependence_hash->{$task."_a".$analysis_id} .= $job_id.":";
					}
					#Fill the hash for cleaning with this jobname and dependencies
					if ($dependency ne "no_depend"){
							$dep_clean_hash->{$job_name}->{'id'} = $job_id;
							$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
					}
				}#STAT END	
			}#Finished processing samples here
	    ################################
      #END PIPELINE EXECUTIONG FOR EACH SAMPLE
      ###################################
			
			#####################################
			######CONSENSUS VARIANT CALLING START
			######################################
			#In this block of instructions we check if the user wants to execute the consensus approach
			#if YES VarGenius runs two jobs for Freebayes and Samtools-mpileup pipelines
			#and changes the dependencies of the fout task that will wait for these two jobs
			#N.B. I decided to put it here because it can start after the statistics step has been run and created 
			#the merged bam files
			#If the consensus approach was used in a previous analysis, it will be always used independently from the
			#configuration file
			#Getting the consensus flag
			my $consensus = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
											$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_consensus'},
											$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
			if ( $varcall and ($cfg_hash_gr->{'consensus_vc_approach'} ne 'NO' or $consensus == 1) ){
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";	
				
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Starting Consensus Variant Calling approach: the BAM file aligned".
				"against the genome will be used for variant calling with other VC tools.\n",$log_file);	
				
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step \n",$log_file_gr);#DEBUGCODE				

				#Managing the dependencies
				my $task_to_wait;
				#Get the samples ids involved for the analysis
				
				#If there is at least 1 sample that has multiple readfile then this job has to wait the
				#stats step for this group
				my $query = "SELECT ".$cfg_hash_gr->{'db_sample_id'}." FROM  ".$cfg_hash_gr->{'db_sample_table'}." WHERE ".$cfg_hash_gr->{'db_analysis_id'}."=$analysis_id;";	
				#print_and_log( "[an:$analysis_id ](".scalar(localtime).") Executing: $query\n",$log_file_gr);#DEBUGCODE
				my $group_sam = do_query_select_all($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},$cfg_hash_gr->{'db_pass'},$query,$cfg_hash_gr->{'db_sample_id'});
				
				my $multiple = 0;
				foreach my $sample_id (keys %{$group_sam}){
				#Verifying that there has been an execution of multiple samples
				 if (get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_multiple'},
									$cfg_hash_gr->{'db_sample_id'},$sample_id) == 1){
					$multiple = 1;
				 }		
				}
				#If there are multiple read files than it must wait for the merge doing in stats
				if ( $multiple == 1) {
					$task_to_wait = "stats";
				}else{#Otherwise go on after the alignment
					$task_to_wait = "align";
				}
								
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);	
				#Get dependencies
				my $dependency = "no_depend";
				my $job_to_wait = $task_to_wait."_a".$analysis_id;
				if (defined $dependence_hash->{$job_to_wait}){
					$dependency = $dependence_hash->{$job_to_wait};	
				}
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE

				my $task_4_fold = $cfg_hash->{'varfilt_task'};	
				my $varfil_jobname = $task_4_fold."_a".$analysis_id;
				
				my @vc_algorithms = split(",",$cfg_hash->{'vc_algorithms'});
				foreach my $vc_algorithm (@vc_algorithms){
					#Run a job for Samtools pipeline
					if ( $vc_algorithm eq 'samtools' ){
						#Call a job for Samtools-mpileup using qsub and giving it parameters for the job 
						my $task = $cfg_hash->{'VCSAM_task'};	
						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
						
						my $job_name = $task."_a".$analysis_id;
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task_4_fold.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task_4_fold.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};						
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);
						
						#Add the jobid for this task to the varfilt list because
						#every task that has to wait the variant filtering to start
						#must wait also the consensus variant calling
						if ( ! defined $dependence_hash->{$varfil_jobname}){
							$dependence_hash->{$varfil_jobname} = $job_id.":";		 	
						}else{
							$dependence_hash->{$varfil_jobname} .= $job_id.":";
						}
						#Fill the hash for cleaning with this jobname and dependencies
						if ($dependency ne "no_depend"){
								$dep_clean_hash->{$job_name}->{'id'} = $job_id;
								$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
						}					
					}
					
					#Run a job for  Freebayes pipeline
					if ( $vc_algorithm eq 'freebayes' ){
						#Call the program using qsub and giving it parameters for the job
						my $task = $cfg_hash->{'VCFB_task'};	

						print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);				
						#Call a job for FreeBayes using qsub and giving it parameters for the job
						my $job_name = $task."_a".$analysis_id;
						my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task_4_fold.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
						my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task_4_fold.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};							
						my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
						my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
						#print_and_log("Qsub id: $job_id\n",$log_file_gr);
						print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

						#Add the jobid for this task to the varfilt list because
						#every task that has to wait the variant filtering to start
						#must wait also the consensus variant calling

						if ( ! defined $dependence_hash->{$varfil_jobname}){
							$dependence_hash->{$varfil_jobname} = $job_id.":";		 	
						}else{
							$dependence_hash->{$varfil_jobname} .= $job_id.":";
						}
						#Fill the hash for cleaning with this jobname and dependencies
						if ($dependency ne "no_depend"){
								$dep_clean_hash->{$job_name}->{'id'} = $job_id;
								$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
						}					
					}							
				}
			
								

				
				#Since this is an analysis with a consensus method update the database with this info
				my $fields = $cfg_hash_gr->{'db_consensus'};
				my $values = "1";
				#Update
				update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$fields,$values,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);			
			}######CONSENSUS VARIANT CALLING END
			
			  
			##########################   		
			#STATISTICS PER-ANALYSIS #
			##########################
			if ( $anstats ){
				my $task = "anstats";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
						push(@tasks_for_clean,$task);
					}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";		
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step \n",$log_file_gr);#DEBUGCODE				

				#Managing the dependencies
				my $task_to_wait;
				#Get the samples ids involved for the analysis
				
				#If there is at least 1 sample that has multiple readfile then this job has to wait the
				#stats step for this group
				my $query = "SELECT ".$cfg_hash_gr->{'db_sample_id'}." FROM  ".$cfg_hash_gr->{'db_sample_table'}." WHERE ".$cfg_hash_gr->{'db_analysis_id'}."=$analysis_id;";	
				#print_and_log( "[an:$analysis_id ](".scalar(localtime).") Executing: $query\n",$log_file_gr);#DEBUGCODE
				my $group_sam = do_query_select_all($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},$cfg_hash_gr->{'db_pass'},$query,$cfg_hash_gr->{'db_sample_id'});
				
				my $multiple = 0;
				foreach my $sample_id (keys %{$group_sam}){
				#Verifying that there has been an execution of multiple samples
				 if (get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
									$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_sample_table'},$cfg_hash_gr->{'db_multiple'}
									,$cfg_hash_gr->{'db_sample_id'},$sample_id) == 1){
					$multiple = 1;
				 }		
				}
				##If there are multiple read files than it must wait for the merge doing in stats
				#if ( $multiple == 1) {
					#$task_to_wait = "stats";
				#}else{#Otherwise go on after the refinement
					#$task_to_wait = "refine";
				#}
				#The statistics on the analysis are performed after the BaseRecalibrator task
				$task_to_wait = "BR";	
									
				#Getting the group name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);	
				#Get dependencies
				my $dependency = "no_depend";
				my $job_to_wait = $task_to_wait."_a".$analysis_id;
				if (defined $dependence_hash->{$job_to_wait}){
					$dependency = $dependence_hash->{$job_to_wait};	
				}
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
				
				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_".$cfg_hash_gr->{$task.'_log'};
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};				
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

				#Put this job id in the hash for dependencies
				#The next step will use the overall group. Hence here we set the 
				#waiting for the groups
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}
			}#STAT END
			
			##############
			#GENOTYPING	
			###############
			
			#The genotyping step must be used only when there are multiple samples in the same analysis
			if ( $genot ){	

				my $task = "genot";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";	
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step\n",$log_file_gr);#DEBUGCODE

				#Managing the dependencies
				my $task_to_wait = "varcall";
				#Getting the group name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);

				#Get dependencies
				my $dependency = "no_depend";
				if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
					$dependency = $dependence_hash->{$task_to_wait."_a".$analysis_id};
					chop($dependency);
				}
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);

				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_".$cfg_hash_gr->{$task.'_log'};
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};
								
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

				#Put this job id in the hash for dependencies
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				chop($dependence_hash->{$job_name});
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}
			}#GENOTYPING END
			
			###################
			#VARIANT FILTRATION 
			##################
			if ($varfilt){
				my $task = "varfilt";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";				
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step\n",$log_file_gr);#DEBUGCODE

				#Managing the dependencies
				#If the genotype step has been run, the task to wait is that, otherwise is the varcall
				my $task_to_wait;
				if ($genot){
						$task_to_wait ="genot";
				}else{
						$task_to_wait ="varcall"
				}
				
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
																				$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
																				$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
						
				#Get dependencies
				my $dependency = "no_depend";
				if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
					$dependency = $dependence_hash->{$task_to_wait."_a".$analysis_id};	
				}
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
				

				#For this step I pass the group name to the script using the SAMPLE_ID environment variable
				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				#my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_".$cfg_hash_gr->{$task.'_log'};
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};
								
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);
				#Put this job id in the hash for dependencies
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				chop($dependence_hash->{$job_name});
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}
			}#VAR_FILTERING END
			
			#######################
			##PHASING
			####################### 
			#This step will be done automatically  if the user asks for variant filtering
			#and if there are three samples with mother, father and proband
			if ( $varfilt and $phasing ){

				my $task = "phasing";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";				
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step\n",$log_file_gr);#DEBUGCODE	
				#print_and_log("Phasing for analysis $analysis_id..\n",$log_file_gr);

				#Managing the dependencies
				my $task_to_wait ="varfilt";
				
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
																				$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
																				$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
						
				#Get dependencies
				my $dependency = "no_depend";
				if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
					$dependency = $dependence_hash->{$task_to_wait."_a".$analysis_id};	
				}
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
					
				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};				
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id] (".scalar(localtime)."): [step:$task - dependencies:$dependency - jobid:$job_id]\n",$log_file_gr);
				#Put this job id in the hash for dependencies
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				chop($dependence_hash->{$job_name});
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}
			}#PHASING END

			####################
			#FAST FINAL OUTPUT
			####################	
			#Fast Final output step
			if ($fastfout){
				my $maintask = "finalout";
				my $task = "fastfout";
				
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";				
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).")  $task step\n",$log_file_gr);#DEBUGCODE	
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") Generating a VCF and annotating it ..\n",$log_file_gr);

				#Managing the dependencies
				my $task_to_wait;
				
				#If the phasing is currently being running, wait that
				if ( $phasing ){
					$task_to_wait	= 'phasing';
				}else{
					$task_to_wait	= 'varfilt';	
				}				
					
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
																				$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
																				$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
				#Get dependencies
				my $dependency = "no_depend";
				if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
					$dependency = $dependence_hash->{$task_to_wait."_a".$analysis_id};	
				}
					
				#Only for the final out step it must be executed concurrently because
				#the contemporary access to the db of many processes may cause strange result
				#Hence here I get the dbbusy field from the database which says if the database is busy
				my $dbbusy = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$cfg_hash_gr->{'db_busy'},
						$cfg_hash_gr->{'db_info_row'},"1");
				#If the database is busy then another field contains the identifier of the 
				#job to wait. This job will wait only for that, the next job will wait this waiter and so on..
				if ($dbbusy == 1 and $cfg_hash_gr->{'block_db'} eq 'YES' ) {
					print_and_log("Dbbusy: $dbbusy\n",$log_file_gr);
					#Get the id of the job that is using the db
					my $jobin = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$cfg_hash_gr->{'db_jobin'},
						$cfg_hash_gr->{'db_info_row'},"1");
					#The process is alone and has no dependencies
					if ($dependency eq 'no_depend'){
						$dependency = $jobin;
					}else{#or it is waiting for something other 
						$dependency .= ":".$jobin;
					}				
				}					
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);
					
				#Call the program using qsub and giving it parameters for the job and
				#The name must be shorter than 15 chars
				my $job_name = $task."_a".$analysis_id;
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$maintask.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$maintask.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};				
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);
				
				if ($job_id ne ''){
					#Now that the job has been launched set the new jobID to wait to this one and set the groupid
					#The flag will be set once that the job will start (in programs_runner.pl)
					#Get the sample sheet id
					my $ss_id = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_sample_sheet_id'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Setting analysis $analysis_id (jobid:$job_id) to be the last waiter of the db for the final output ...\n",$log_file_gr);
					
					#BLOCKDB
					#If the output is printed in the mode blocking the db than set it as busy and put the jobid
					if ( $cfg_hash_gr->{'block_db'} eq 'YES' ) {
						#UPDATE fields = values
						my $fields = $cfg_hash_gr->{'db_jobin'}.";".$cfg_hash_gr->{'db_analysis_id'}.";".$cfg_hash_gr->{'db_sample_sheet_id'};
						my $values = $job_id.";".$analysis_id.";".$ss_id;
						#If the db was not busy now it is
						if ($dbbusy == 0) {
							$fields .= ";".$cfg_hash_gr->{'db_busy'};
							$values .= ";1";
						}
						#WHERE 
						my $fields2 = $cfg_hash_gr->{'db_info_row'};
						my $values2 = "1";
						#Update
						update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
								$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$fields,$values,$fields2,$values2);					
					}

								
					#Put this job id in the hash for dependencies
					if ( ! defined $dependence_hash->{$job_name}){
						$dependence_hash->{$job_name} = $job_id.":";		 	
					}else{
						$dependence_hash->{$job_name} .= $job_id.":";
					}
					chop($dependence_hash->{$job_name});
					#Fill the hash for cleaning with this jobname and dependencies
					if ($dependency ne "no_depend"){
							$dep_clean_hash->{$job_name}->{'id'} = $job_id;
							$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
					}
				}else{
					log_and_exit("[an:$analysis_id ](".scalar(localtime).") ERROR: Cannot get a correct jobid. Check if there are problems with job execution.\n",$log_file_gr);
				}

			}#Fast Final Output END
			
			####################
			#CNV 
			####################
			#CNV are detected using ExomeDepth and XHMM.
			#They both need the BAM files at the end of the process of BaseRecalibration
			if ( $cnv_detection){
				my $task = "cnv";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";	
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step\n",$log_file_gr);#DEBUGCODE

				#Managing the dependencies
				#This job has to wait the BR task
				my @tasks_to_wait = ("BR");
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);

				#Get dependencies
				my $dependency = "";
				
				foreach my $task_to_wait (@tasks_to_wait){
					if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
						$dependency .= $dependence_hash->{$task_to_wait."_a".$analysis_id}.":";
						
					}					
				}
				#If empty, set nodepend
				if ( $dependency eq "" ){$dependency = "no_depend";}else{chop($dependency);}
			
								
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);

				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$task.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

				#Put this job id in the hash for dependencies
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				chop($dependence_hash->{$job_name});
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}

				#update the database to say that CNV analysis was performed
				my $fields = $cfg_hash_gr->{'db_analysis_cnv'};
				my $values = "1";
				#Update
				update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$fields,$values,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);	
			}
			
			#CNV OUTPUT : generate a tabular human readable output from the output
			#of the programs for CNV detection
			if ( $cnv_out){
				my $maintask = "cnv";
				my $task = "cnvo";
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";	
				print_and_log($niceClosure,$log_file_gr);
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") $task step\n",$log_file_gr);#DEBUGCODE

				#Managing the dependencies
				#This job has to wait the cnv detection job
				my @tasks_to_wait = ("cnv");
				#Getting the analysis name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);

				#Get dependencies
				my $dependency = "";
				
				foreach my $task_to_wait (@tasks_to_wait){
					if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
						$dependency .= $dependence_hash->{$task_to_wait."_a".$analysis_id}.":";	
					}					
				}
				#If empty, set nodepend
				if ( $dependency eq "" ){$dependency = "no_depend";}else{chop($dependency);}
			
								
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task : ",$log_file_gr);

				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$maintask.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$maintask.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);

				#Put this job id in the hash for dependencies
				if ( ! defined $dependence_hash->{$job_name}){
					$dependence_hash->{$job_name} = $job_id.":";		 	
				}else{
					$dependence_hash->{$job_name} .= $job_id.":";
				}
				chop($dependence_hash->{$job_name});
				#Fill the hash for cleaning with this jobname and dependencies
				if ($dependency ne "no_depend"){
						$dep_clean_hash->{$job_name}->{'id'} = $job_id;
						$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
				}

				#Commented because I do not know if to use an additional flag for the cnv output
				##update the database to say that CNV analysis was performed
				#my $fields = $cfg_hash_gr->{'db_analysis_cnv'};
				#my $values = "1";
				##Update
				#update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						#$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$fields,$values,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);			
				
			}
			####################
			#CNV 
			####################
			
			
			
			###################
			# VCF IMPORT
			##################
			#VCF import step
			if ($vcfimport){
				my $fold_suff = 'finalout';
				my $task = $cfg_hash_gr->{'vcfimport_task_s'};
				#Put in the config hash a variable with all the tasks to be executed
				if ( not grep {/\b$task\b/} @tasks_for_clean){
					push(@tasks_for_clean,$task);
				}
				my $analysis_indication2 = "[an:$analysis_id](".scalar(localtime)."): ";				
				print_and_log($niceClosure,$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).")  $task step\n",$log_file_gr);	#DEBUGCODE
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") Importing VCF in DB for analysis $analysis_id..\n",$log_file_gr);#DEBUGCODE

				
				#If variants must be inserted, check here if the variants from this analysis were already inserted
				#and alert the user if true. She/he cannot run this command
				if ($cfg_hash_gr->{'vcfimport'} eq 'YES'){
					my $query = "SELECT ".$cfg_hash_gr->{'db_analysis_id'}." FROM  ".$cfg_hash_gr->{'db_genotype_sample_table'}." WHERE ".$cfg_hash_gr->{'db_analysis_id'}." = $analysis_id;";	
					#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
					my $anal_ids_fetch = do_query_select_all($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},$cfg_hash_gr->{'db_pass'},$query,$cfg_hash_gr->{'db_analysis_id'});
					#Now check if they are present
					if ( scalar(keys %{$anal_ids_fetch}) > 0 ){
						print_and_log("[an:$analysis_id ](".scalar(localtime).") WARNING: Variants from analysis  $analysis_id will not be imported. I found ".scalar(keys %{$anal_ids_fetch}).
						" genotypes into the ".$cfg_hash_gr->{'db_genotype_sample_table'}." table. Please check this information. Exiting...\n",$log_file_gr);
					}					
				}
				
				#Managing the dependencies
				my $task_to_wait = 'fastfout';
					
				#Getting the group name
				my $analysis_name = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
																				$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_analysis_name'}
																				,$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
							
				#Get dependencies
				my $dependency = "no_depend";
				if (defined $dependence_hash->{$task_to_wait."_a".$analysis_id}){
					$dependency = $dependence_hash->{$task_to_wait."_a".$analysis_id};	
				}
		
				#Only this step to import the VCF must be executed concurrently if there are several 
				#analyses running because
				#the contemporary access to the db of many processes may cause strange result
				#Hence here I get the dbbusy field from the database which says if the database is busy
				#and add that dependency
				my $dbbusy = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$cfg_hash_gr->{'db_busy'},
						$cfg_hash_gr->{'db_info_row'},"1");
				#If the database is busy then another field contains the identifier of the 
				#job to wait. This job will wait only for that, the next job will wait this waiter and so on..
				#EXPERIMENTAL: VCFIMPORT WILL ALWAYS WAIT IF DBBUSY=1
				#if ($dbbusy == 1) {
				if ($dbbusy == 1 and $cfg_hash_gr->{'always_block_db_vcfi'} eq 'YES' ) {
					
					#print_and_log("Dbbusy: $dbbusy\n",$log_file_gr);#DEBUGCODE
					#Get the id of the job that is using the db
					my $jobin = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
						 $cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$cfg_hash_gr->{'db_jobin'},
						 $cfg_hash_gr->{'db_info_row'},"1");

					#Check if the job is in an active queue
					my $active_q = JOB_is_active_queue($jobin);
					if ( $active_q ){
						#The process is alone and has no dependencies
						if ($dependency eq 'no_depend'){
							$dependency = $jobin;
						}else{
							#If the fastfout is the job 
							#or it is waiting for something other 
							$dependency .= ":".$jobin;
						}		
					}else{
						$dbbusy = 0;
					}
					##Check if the job has finished
					#my $status_field = 9;
					#my $jobstat = JOB_get_status_field($jobin,$status_field,$log_file);
					#if ( $jobstat =~ /[RQ]/){
						##The process is alone and has no dependencies
						#if ($dependency eq 'no_depend'){
							#$dependency = $jobin;
						#}else{
							##If the fastfout is the job 
							##or it is waiting for something other 
							#$dependency .= ":".$jobin;
						#}							
					#}else{
						#$dbbusy = 0;
					#}
			
				}					
				#print_and_log("[an:$analysis_id ](".scalar(localtime).") This job will wait for these jobs: $dependency\n",$log_file_gr);#DEBUGCODE
				print_and_log("[an:$analysis_id ](".scalar(localtime).") Launching a job for $task :",$log_file_gr);
					
				#Call the program using qsub and giving it parameters for the job and
				my $job_name = $task."_a".$analysis_id;
				my $job_log = $cfg_hash_gr->{$analysis_id.'_'.$fold_suff.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'log_ext'};
				my $job_err = $cfg_hash_gr->{$analysis_id.'_'.$fold_suff.'_log_fold'}."/".$analysis_name."_$task.".$cfg_hash_gr->{'err_ext'};
				my $env_vars = "CONFIG_FILE=$cfg_hash_f,TASK=$task,SAMPLE_ID=$analysis_name,ANALYSIS_ID=$analysis_id,LOG_FILE=$log_file_gr,PIPE=$pipeline ";
				my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$program_call,$job_name,$dependency,$job_log,$job_err,$log_file_gr);						
				#print_and_log("Qsub id: $job_id\n",$log_file_gr);
				print_and_log("[an:$analysis_id ](".scalar(localtime).") step:$task - dependencies:$dependency - jobid:$job_id\n",$log_file_gr);
				
				if ($job_id ne ''){
					#Now that the job has been launched set the new jobID to wait to this one and set the analysisid
					#The flag will be set once that the job will start (in programs_runner.pl)
					#Get the sample sheet id
					my $ss_id = get_id_if_exists_from_db($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
							$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_analyses_table'},$cfg_hash_gr->{'db_sample_sheet_id'},
							$cfg_hash_gr->{'db_analysis_id'},$analysis_id);
					print_and_log("[an:$analysis_id ](".scalar(localtime).") Setting analysis $analysis_id (jobid:$job_id) to be the last waiter of the db for the final output ...\n",$log_file_gr);
					
					#The VCF import can be always blocked with 'always_block_db_vcfi' parameter.
					if ( $cfg_hash_gr->{'always_block_db_vcfi'} eq 'YES' ) {
						#UPDATE fields = values
						my $fields = $cfg_hash_gr->{'db_jobin'}.";".$cfg_hash_gr->{'db_analysis_id'}.";".$cfg_hash_gr->{'db_sample_sheet_id'};
						my $values = $job_id.";".$analysis_id.";".$ss_id;

						#If the db was not busy now it is
						if ($dbbusy == 0) {
							$fields .= ";".$cfg_hash_gr->{'db_busy'};
							$values .= ";1";
						}
						#WHERE 
						my $fields2 = $cfg_hash_gr->{'db_info_row'};
						my $values2 = "1";
						#Update
						update_table($cfg_hash_gr->{'db_name'},$cfg_hash_gr->{'db_dsn'},$cfg_hash_gr->{'db_user'},
								$cfg_hash_gr->{'db_pass'},$cfg_hash_gr->{'db_info_table'},$fields,$values,$fields2,$values2);					
					}
					
								
					#Put this job id in the hash for dependencies
					if ( ! defined $dependence_hash->{$job_name}){
						$dependence_hash->{$job_name} = $job_id.":";		 	
					}else{
						$dependence_hash->{$job_name} .= $job_id.":";
					}
					chop($dependence_hash->{$job_name});
					#Fill the hash for cleaning with this jobname and dependencies
					if ($dependency ne "no_depend"){
							$dep_clean_hash->{$job_name}->{'id'} = $job_id;
							$dep_clean_hash->{$job_name}->{'deps'} = $dependency;
					}
				}else{
					log_and_exit("[an:$analysis_id ](".scalar(localtime).") ERROR: Cannot get a correct jobid. Check if there are problems with job execution.\n",$log_file_gr);
				}
			}#VCF import END
		
		}#Finished processing groups here		

		
		
		########################
		## JOB CLEANER TASK
		########################
		#This task starts only whether jobs dependencies are present
		if ( defined $dep_clean_hash or $clean_analyses){
			$cfg_hash->{'tasks_to_exec'} = join(",",@tasks_for_clean); 
			#Obtain an hash specific for the cleaning
			my $cfg_hash_gr;
			get_modified_config_hash($cfg_hash,"dummy",\$cfg_hash_gr);
			
			#Define the time needed for the job cleaner as the sum of the WallTime
			#requested for the other jobs
			my $hours_need = (split(":",$cfg_hash->{'qsub_walltime'}))[0];
			my $tasks_req;
			if ( scalar(@tasks_for_clean) >0 ){
					$tasks_req = scalar(@tasks_for_clean);
			}else{$tasks_req = 1;}
			#$cfg_hash_gr->{'qsub_walltime'} = "24:00:00"; #($hours_need * $tasks_req).":00:00";
			print_and_log("Walltime needed for JobCleaner: ".$cfg_hash_gr->{'qsub_walltime'}."\n",$log_file);
			
			#Dependencies
			if ( defined $dep_clean_hash ){
				print_and_log("Saving the hash with dependencies..\n",$log_file);
				#Saving the hash with the dependencies
				save_hash(\$dep_clean_hash,$deps_f);
			}
			#Saving the general hash
			my $cl_hash_f = $main_data_folder."/".$cfg_hash->{'jc_hash_f'};
			save_hash(\$cfg_hash_gr,$cl_hash_f);		

			#Running a job for jobs cleaning
			my $task = "jobclean";
			print_and_log($niceClosure,$log_file);
			print_and_log("(".scalar(localtime)."): Running $task \n",$log_file);	
			my $job_log = "$log_file.$task";
			my $env_vars = "CONFIG_FILE=$cl_hash_f,DEP_HASH_F=$deps_f,TASK=$task,LOG_FILE=$log_file,JC_TIME=".
						$cfg_hash->{'jobclean_time'}.",QS_USER=".$cfg_hash->{'qsub_username'};
			my $job_id = try_exec_job( $cfg_hash_gr,$env_vars,$task,$jobclean_prog,$task,'no_depend',$job_log,$job_log,$job_log);	
		}		
	}		

	
	#Remove the config file saved
	if(!$genedb and !$reannote and !$update_freqs and !$geneann and !$geneann and !$dldatasets){
			delete_file($config_file);
	}
	
}



###################################CHECK SUBROUTINES



=head2 check_sample_sheet

 Title   : check_sample_sheet
 Usage   : check_sample_sheet( 	sample_sheet => sample sheet path,
																ss_head_file => file with the header
																sep => separator used
							);

 Function:  checks the sample sheet before to start the pipeline.
						It uses the file 'ss_head_file' containing the fields which
						have to be present in the header of the sample sheet.
						This fields are the same as those in the database.
 
 Returns : nothing

=cut
sub check_sample_sheet {
	my $sample_sheet = shift;
	my $sep = shift;
	my $ss_head_file = shift;
	my $ss_head_sep = shift;
	
	#To check the analysisname field
	my $analysis_name_fld = $cfg_hash->{'db_analysis_name'};
	my $sample_gender_fld = $cfg_hash->{'db_sample_gender'};
	my $targetbed_fld = $cfg_hash->{'db_targetbed'};
	my $kinship_fld = $cfg_hash->{'db_sample_kinship'};
	my $fqdir_fld = $cfg_hash->{'db_readf_fqdir'};
	my $fq1_fld = $cfg_hash->{'db_readf_fq1'};
	my $fq2_fld = $cfg_hash->{'db_readf_fq2'};
	
	print_and_log( "Checking sample sheet $sample_sheet...",$log_file);
	
	if ( -e $sample_sheet){
		#Here I check if the same sample sheet was used previously
		if (exists_ss_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_ss_table'},$cfg_hash->{'db_ss_name'},extract_name($sample_sheet,0)) == -1){
			#Here we check if the sample sheet header correspond to the one we have in input
			#opening the ss header from the configuration folder
			open (ss_head,$ss_head_file) or die "ERROR: Cannot open $ss_head_file. The program will exit..\n";
			open (s_sheet,$sample_sheet) or die "ERROR: Cannot open $sample_sheet. The program will exit..\n";
			#Take the first line
			my $line_ss = <s_sheet>;
			my $line_ssh = <ss_head>;
			#Remove the newline
			chomp($line_ss);
			chomp($line_ssh);
			#Put fields in an array
			my @head_ss = split($sep,$line_ss);
			my @head_ssh = split($ss_head_sep,$line_ssh);
			#Check if number of fields are the same
			my $lines_ss = scalar(@head_ss);
			my $lines_ssh = scalar(@head_ssh);
			if( $lines_ss == $lines_ssh){
				my $num_line = 0;
				#Checks field by field if they are the same
				foreach my $field_ss (@head_ss ){
					#The order must be the same as that in the header file
					if( $field_ss ne $head_ssh[$num_line] ){
						log_and_exit("$field_ss is not a correct field for the sample sheet".
								" or not in the correct position. Check $ss_head_file sorting..",$log_file);					
					}
					$num_line++;
				}
				#Search all the lines in search of errors
				$num_line = 0;
				while (my $row = <s_sheet>) {
							my @fields = split($sep,$row);
							my $num_col  = 0;
							
							my $fqdir = "";
							foreach my $field (@fields){
								#Check the analysisname
								if (  $head_ssh[$num_col] =~ /$analysis_name_fld/){
									#print_and_log("Checking Analysis name $field!\n",$log_file);#DEBUGCODE
									#The analysisname must be there...
									log_and_exit("No ".$head_ssh[$num_line]."... no party! Please add it in the sample sheet (line $num_line). Exiting..\n",$log_file)
									unless ($field ne '' and $field ne '-');
									
									#Check if the analysisname is present into the db, unless force_run is used (e.g. new samples for a Joint analysis)
									unless ( $force_run ){
										#Getting the analysis name
										my $analysisname = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																					$cfg_hash->{'db_analysis_name'},"'".$field."'");			
										log_and_exit("ERROR: ( line ".($num_line+1).") $analysis_name_fld $field exists already into the database! Please change this name as $program_name cannot ".
										" manage analysis with the same name. Exiting..\n",$log_file)
										unless ($analysisname eq "-1");	
									}
													
								}
								#The gender must be one among M and F
								if (  $head_ssh[$num_col] =~ /$sample_gender_fld/){
									#print_and_log("Checking Gender $field!\n",$log_file);#DEBUGCODE
									if ( ($field ne 'M') and ( $field ne 'F') and ( $field ne '-') ){
										log_and_exit("ERROR: ( line ".($num_line+1).") The value $field is invalid as $sample_gender_fld at line $num_line. Please use only M or F) . Exiting..\n",$log_file)
									}
								}	
								#Check the kinship that must be one among kinship_type
								if (  $head_ssh[$num_col] =~ /$kinship_fld/){
									#print_and_log("Checking kinship $field!\n",$log_file);#DEBUGCODE
									my $kinship_types = join("",split(",",$cfg_hash->{'kinship_types'}));
									if ( $field !~ /[$kinship_types]/){
											log_and_exit("ERROR: ( line ".($num_line+1).")The value $field is invalid as $kinship_fld at line $num_line. Please use one among ".$cfg_hash->{'kinship_types'}.") . Exiting..\n",$log_file)
									}
								}
								#The target bed file must be present
								if (  $head_ssh[$num_col] =~ /$targetbed_fld/){
									#print_and_log( "Checking the target file name: $field ...\n",$log_file);#DEBUGCODE
									#Check length and chars of the target file name
									file_name_wostrange_chars($field,$cfg_hash->{'bed_ext'});
									#if ( ! (-e $cfg_hash->{'target_reg_f'}."/".$field) ){
									#	log_and_exit("ERROR: ( line ".($num_line+1).") The target bed file $field does not exist in ".$cfg_hash->{'target_reg_f'}.". Please put it there!!\n",$log_file)
									#}
								}	
								#The fastq files must be present
								if (  $head_ssh[$num_col] =~ /$fqdir_fld/){
									$fqdir = $field;
								}
								if (  $head_ssh[$num_col] =~ /$fq1_fld/ or $head_ssh[$num_col] =~ /$fq2_fld/){
									#print_and_log( "Checking the fastq file name: $field ...\n",$log_file);#DEBUGCODE
									unless (-e $fqdir."/".$field ){
										log_and_exit("ERROR: ( line ".($num_line+1).") The fastq file $field does not exist in $fqdir. Please put it there!!\n",$log_file)
									}								
								}																						
								$num_col++;
							}
		
					$num_line++;
				}
			}else{
					log_and_exit("Sample sheet $sample_sheet has $lines_ss fields ".
							"while they should be $lines_ssh. Please check the error...\n",$log_file);
			}

			
			print_and_log("...DONE!\n",$log_file);
		}else{
			print_and_log("The sample sheet $sample_sheet is already present in the database\n",$log_file);
		}		
	}
}


=head2 sample_sheet_2_db

 Title   : sample_sheet_2_db
 Usage   : sample_sheet_2_db( 	sample_sheet => sample sheet path,
								ss_hash => hash that will be used to store all the information about the samples
								sep => separator used
							);

 Function:  read the sample sheet and assigns an id (sequential number)
						to each sample. Then fills the database in the corresponding table
 
 Returns : nothing

=cut
sub sample_sheet_2_db{
	my $sample_sheet = shift;
	my $sep = shift;

	#The field where there is the sample id
	my $id_field = -1;
	my $sample_id = 0;
	
	#Used to remember analyses already inserted in db
	my $analyses_inserted;
	my $samples_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	print_and_log("Reading the sample sheet $sample_sheet and storing in database ".$cfg_hash->{'db_name'}."...\n ",$log_file);
	
	#Here I check if the same sample sheet was used previously
	if (exists_ss_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
		$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_name'},extract_name($sample_sheet,0)) == -1){
		my $ss_name = "'".extract_name($sample_sheet,0)."'";
		###############################
		#INSERT SAMPLE SHEET NAME			#
		###############################
		#Insert the sample sheet id into the database and gets the id generated
		my $ss_id = insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_name'},$ss_name);
		#The last id is positive if the sample sheet name is succesfully inserted
		if ($ss_id >=0 ){
			#Store the analysis id
			$cfg_hash->{'ss_id'} = $ss_id;
			open (s_sheet,"<$sample_sheet") or log_and_exit("ERROR: Cannot open $sample_sheet. The program will exit...\n",$log_file);
							
			#Go line by line and get the field values for each sample
			my $row_num = 0;
			while (my $line = <s_sheet>){
				#print "Analyzing line: ".$line;
				chomp($line);
				
				#######READ HEADER
				#Get the field names from the header
				if ($row_num == 0){
					#Put fields in a way to be used
					my @head = split($sep,$line);
					my $n_field = 0;
					#each field will be associated to the index
					foreach my $field (@head){
							$head_hash->{$n_field} = $field;
							$n_field++;
					}
				}else{#The row is composed of data
					my $analysis_id = -1;
					my @fields = split($sep,$line);
					my $n_field = 0;
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}
					#print "\n";
					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					#Putting the analysisname and analysis id in the analysis table
					my $analysisname = $ss_hash->{'analysisname'};
					if ( ! defined $analyses_inserted->{$analysisname} ){
							my $fields = $cfg_hash->{'db_analysis_name'}.",".$cfg_hash->{'db_sample_sheet_id'};
							my $values = "'$analysisname',$ss_id";
							#The second values are used to evaluate if the analysisname is present
							my $fields2 = $cfg_hash->{'db_analysis_name'};
							my $values2 = "'$analysisname'";
														
							#Check if the research group identifier is present
							if ( $ss_hash->{$cfg_hash->{'db_analysis_resgroupid'}} ne '-'){
									$fields .= ",".$cfg_hash->{'db_analysis_resgroupid'};
									$values .= ",".$ss_hash->{$cfg_hash->{'db_analysis_resgroupid'}};
							}
							#Check if the user identifier is present
							if ( $ss_hash->{$cfg_hash->{'db_analysis_userid'}} ne '-'){
								$fields .= ",".$cfg_hash->{'db_analysis_userid'};
								$values .= ",".$ss_hash->{$cfg_hash->{'db_analysis_userid'}};
							}

							#Check if the user wants to run this analysis per-chromosome and set the db field for that
							if ( $ss_hash->{$cfg_hash->{'db_perchrom'}} == 1){
								$fields .= ",".$cfg_hash->{'db_perchrom'};
								$values .= ",1";
							}
							#Check if the user wants to run this analysis per-chromosome and set the db field for that
							if ( $ss_hash->{$cfg_hash->{'db_targetbed'}} ne '-' ){
								$fields .= ",".$cfg_hash->{'db_targetbed'};
								$values .= ",'".$ss_hash->{'targetbed'}."'";
							}
							if ( $ss_hash->{$cfg_hash->{'db_sequencingtype'}} ne '-'){
								$fields .= ",".$cfg_hash->{'db_sequencingtype'};
								$values .= ",'".$ss_hash->{'seqtype'}."'";
							}
							#Field for the step of conversion of the quality scores
							if ( $ss_hash->{$cfg_hash->{'db_convscores'}} ne '-'){
								$fields .= ",".$cfg_hash->{'db_convscores'};
								$values .= ",".$ss_hash->{'convscores'}."";
							}
							#Field for the step of conversion of the quality scores
							if ( $ss_hash->{$cfg_hash->{'db_analysis_infreq'}} ne '-'){
								$fields .= ",".$cfg_hash->{'db_analysis_infreq'};
								$values .= ",".$ss_hash->{$cfg_hash->{'db_analysis_infreq'}}."";
							}					
							#print "DB: Inserting $fields = $values in table ".$cfg_hash->{'db_analyses_table'}."\n";#DEBUGCODE
							#Inserting target file into the target file table
							insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_target_files_table'},$cfg_hash->{'db_targetid'},
													$cfg_hash->{'db_targetname'},"'".$ss_hash->{'targetbed'}."'",$cfg_hash->{'db_targetname'},"'".$ss_hash->{'targetbed'}."'");
													
							my ($new_analysis_id,$exist_analysis_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
													$fields,$values,$fields2,$values2);
							#Save in the hash
							$analyses_inserted->{$analysisname}->{'anal_id'} = $ss_id;
							#If the groupid exists insert that, otherwise the new one
							if ( $new_analysis_id >= 0){
								$analyses_inserted->{$analysisname}->{'analysis_id'} = $new_analysis_id;
								$analysis_id = $new_analysis_id;
							}elsif ( $exist_analysis_id >= 0) {
								$analyses_inserted->{$analysisname}->{'analysis_id'} = $exist_analysis_id;
								$analysis_id = $exist_analysis_id;
							}else{
								log_and_exit("ERROR: during the import of the analysis_id for $analysisname. Cannot understand what is happen with the db...\n",$log_file);
							}
					}else{
							$analysis_id = $analyses_inserted->{$analysisname}->{'analysis_id'};
					}
					###############################
					#INSERT SAMPLE INFO						#
					###############################					
					#Putting the sample name and analysis id in the samples table
					my $sampleName = $ss_hash->{'samplename'};
					if ( ! defined $samples_inserted->{$sampleName} ){
							####Update the samples_info table
							my $fields = $cfg_hash->{'db_sample_name'}.",";
							my $values = "'$sampleName',";
							#Fields to use to check if the person do not exist yet
							my $fields2 = $cfg_hash->{'db_sample_name'};
							my $values2 = "'$sampleName'";
							
							#Get the family id from the sample name
							#If there is '_'
							my $familyid = $analysisname;					
							if ( $sampleName =~ '_' ){$familyid = (split("_",$sampleName))[1];}
							
							if (defined $familyid ){
								$fields .= $cfg_hash->{'db_familyid'}.",";
								$values .= "'$familyid',";
						  }
							#Add Date of birth if exists
							if ( $ss_hash->{$cfg_hash->{'db_sample_dob'}} ne '-' ){
								$fields .= $cfg_hash->{'db_sample_dob'}.",";
								$values .= "'".$ss_hash->{$cfg_hash->{'db_sample_dob'}}."',";
							}
							#Add place of birth if exists
							if ( $ss_hash->{$cfg_hash->{'db_sample_pob'}} ne '-' ){
									$fields .= $cfg_hash->{'db_sample_pob'}.",";
									$values .= "'".$ss_hash->{$cfg_hash->{'db_sample_pob'}}."',";
							}
							#Add Kinship if exists
							if ( $ss_hash->{$cfg_hash->{'db_sample_kinship'}} ne '-' ){
								$fields .= "".$cfg_hash->{'db_sample_kinship'}.",";
								$values .= "'".$ss_hash->{$cfg_hash->{'db_sample_kinship'}}."',";
							}
							#Check if the gender identifier is present
							if ( $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne '-' ){
									$fields .= "".$cfg_hash->{'db_sample_gender'}.",";
									$values .= "'".$ss_hash->{$cfg_hash->{'db_sample_gender'}}."',";
							}
							#Check if the sample affected flag is present
							if ( $ss_hash->{$cfg_hash->{'db_sample_affected'}} ne '-' ){
									$fields .= $cfg_hash->{'db_sample_affected'}.",";
									$values .= $ss_hash->{$cfg_hash->{'db_sample_affected'}}.",";
									#my $affected = $ss_hash->{$cfg_hash->{'db_sample_affected'}};
									#if ($affected == 1){ 
										#$values .= "true,";
									#} else{
										#$values .= "false,";	
									#}
							}
							##Check if the heritability is present
							#if ( $ss_hash->{$cfg_hash->{'db_sample_heritability'}} ne '-' ){
									#$fields .= ",".$cfg_hash->{'db_sample_heritability'};
									#$values .= ",".$ss_hash->{$cfg_hash->{'db_sample_heritability'}}."";
							#}

							
							chop($fields);
							chop($values);
							#print "DB: Inserting $fields = $values in table ".$cfg_hash->{'db_analyses_table'}."\n";#DEBUGCODE
							my $personid = -1;
							
							my ($new_personid,$exist_personid) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_info_table'},$cfg_hash->{'db_samples_personid'},
													$fields,$values,$fields2,$values2);
					
							#If the personid exists insert that, otherwise the new one
							if ( $new_personid >= 0){
								$personid = $new_personid;
							}elsif ( $exist_personid >= 0) {
								$personid = $exist_personid;
							}else{
								log_and_exit("ERROR: during the import of the personid for $sampleName. Cannot understand what is happen with the db...\n",$log_file);
							}		
												
							###Now insert data into the samples table													
							$fields = $cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_samples_personid'}.",".$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_sheet_id'};
							$values = "'$sampleName',$personid,$analysis_id,$ss_id";
							#If the sample has multiple files it adds the flag (1) to the multiple field
							if ( $ss_hash->{'samplename'} ne $ss_hash->{'readfname'}){
								$fields .= ",".$cfg_hash->{'db_multiple'};
								$values .= ",1";
							}
							#Add Kinship if exists
							if ( $ss_hash->{'kinship'} ne '-' ){
								$fields .= ",".$cfg_hash->{'db_sample_kinship'};
								$values .= ",'".$ss_hash->{'kinship'}."'";
							}
							#Check if the gender identifier is present
							if ( $ss_hash->{$cfg_hash->{'db_sample_gender'}} ne '-' ){
									$fields .= ",".$cfg_hash->{'db_sample_gender'};
									$values .= ",'".$ss_hash->{$cfg_hash->{'db_sample_gender'}}."'";
							}
							#print "DB: Inserting $fields = $values in table ".$cfg_hash->{'db_analyses_table'}."\n";#DEBUGCODE
							$sample_id = insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$fields,$values);
													
							#Save in the hash
							$samples_inserted->{$sampleName}->{'ss_id'} = $ss_id;
							$samples_inserted->{$sampleName}->{'sample_id'} = $sample_id;
					}else{
							$sample_id = $samples_inserted->{$sampleName}->{'sample_id'};
					}
					###############################
					#INSERT READFILE INFO					#
					###############################
					#To insert read file info, we need that both the group and sample infos have been inserted
					#print "The analysis_id  is: $analysis_id\n";#DEBUGCODE	
					if ( $analysis_id >= 0 and $sample_id >= 0 ){
						#Inserting each sample in the database
						my $fields = $cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_readf_fqdir'}.",".$cfg_hash->{'db_readf_fq1'}.",".
												 $cfg_hash->{'db_readf_fq2'}.",".$cfg_hash->{'db_readf_name'}.",".
												 $cfg_hash->{'db_sample_sheet_id'}.",".$cfg_hash->{'db_analysis_id'};
						my $values = "$sample_id,'".$ss_hash->{'fqdir'}."','".$ss_hash->{'fq1'}."','".
													$ss_hash->{'fq2'}."','".$ss_hash->{'readfname'}."',$ss_id,$analysis_id";

						my $readf_id = insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$fields,$values);						
					
						#Here I get information about the flowcell for each fastq
						my ($instrument,$flowcell,$lane) = get_flowcell_and_lane($cfg_hash,$ss_hash->{'fqdir'}."/".$ss_hash->{'fq1'},$log_file);
						#..and put them into the db
						#print_and_log("Updating ".$cfg_hash->{'db_readf_table'}." putting flowcell: $flowcell and lane $lane for $readf_id ..\n",$log_file);#DEBUGCODE
						update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
												$cfg_hash->{'db_readf_instrument'},"'".$instrument."'",$cfg_hash->{'db_readf_id'},$readf_id);
						update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
												$cfg_hash->{'db_readf_flowcell'},"'".$flowcell."'",$cfg_hash->{'db_readf_id'},$readf_id);
						update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
												$cfg_hash->{'db_readf_lane'},"'".$lane."'",$cfg_hash->{'db_readf_id'},$readf_id);						
					}
				}
				$row_num++;
			}
			close(s_sheet);
		}else{
			log_and_exit("There was some error inserting $sample_sheet in the database",$log_file);
		}
	}else{
		log_and_exit("The sample sheet $sample_sheet is already present in the database",$log_file);
	}
}


=head2 set_global_var_and_confighash

 Title   : set_global_var_and_confighash
 Usage   : set_global_var_and_confighash(   );

 Function:  Just sets important variables and imports the parameters from the configuration 	
					file into the hash
 
 Returns : nothing

=cut
sub set_global_var_and_confighash{
	my $config_user = shift;
	
	print_and_log($vargenius_title,$log_file);
	##Set the folders name using inizialize_folders
	print_and_log($configClosure,$log_file);	
	($working_folder,$program_folder) = initialize_folders($foldersFile);
	$config_folder = $program_folder."/".$config_folder;
	$program_config = $config_folder."/".$program_config;
	$variables_file = $config_folder."/".$variables_file;

	##Set important variables
	$config_file = $working_folder."/conf_$time.hash";
	
	#Check if config files are ok (1 and 2 indicates the line into the variables file
	print_and_log("Checking $program_config ",$log_file);	
	#print_and_log("using line 1 of $variables_file...\n",$log_file);	#DEBUGCODE
	checkConfigVariables($program_config,$variables_file,1);
	print_and_log("Checking $config_user ",$log_file);	 
	#print_and_log("using line 2 of $variables_file...\n",$log_file);	#DEBUGCODE
	checkConfigVariables($config_user,$variables_file,2);
	
	#loads the parameters in the configuration hash as from the configuration files
	print_and_log("Loading the parameters for the run as from $program_config and $config_user...\n",$log_file);	
	configFile2Hash($program_config,\$cfg_hash);
	configFile2Hash($config_user,\$cfg_hash);

	#Same some variable in the config hash
	$cfg_hash->{'work_fold'} = $working_folder;
	$cfg_hash->{'program_folder'} = $program_folder;
	$cfg_hash->{'program_name'} = $program_name;
	$cfg_hash->{'logFolder'} = $logFolder;
	$cfg_hash->{'lib_folder'} = $program_folder."/".$cfg_hash->{'lib_folder'};
	

	
	#print Dumper\$cfg_hash;#DEBUGCODE
}

=head2 preliminary_configuration

 Title   : preliminary_configuration
 Usage   : preliminary_configuration(   );

 Function:  preliminary_configuration of the pipeline:
							- set the folders name using inizialize_folders
							- loads the parameters in the hash 'cfg_hash' as from the configuration files
							- creates or fills or check if the database exists
							- uses the samplesheet to generate a configuration for the execution
							- if run_analyses and run_samples are used the configuration is modified
								using 'exec_hash'
							- creates the folders for the analyses as given by the configuration
							- sets the reference genome
							- saves the configuration hash 'cfg_hash' in a common file 'log_file'
							
 Returns : nothing

=cut
sub preliminary_configuration{

	print_and_log("Creating the database if needed..\n",$log_file);	
	##creates or fills or check if the database exists
	if ( $create_database ){
		#Create the database if asked
		if ( (db_present($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'}) == 0) or ($force_run) ){
			print_and_log("Creating the database ".$cfg_hash->{'db_name'}.". The old one will be removed..\n",$log_file);
			#Create database and tables
			createDBAndTables($cfg_hash->{'db_name'},$config_folder."/".$cfg_hash->{'db_schema'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_host'});
			#Fill the INFO table
			configure_database($cfg_hash);
		}else{
			log_and_exit("The database ".$cfg_hash->{'db_name'}." exists. Please use the parameter --force_run to remove it..\n",$log_file);
			}
			exit;
	}#Elsif the user wants to create only the tables in an existing database
	elsif ( $create_tables ){
		
		if ( db_present($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'}) == 1 ){
				#Fills tables
				print_and_log("Filling the database ".$cfg_hash->{'db_name'}." using the sql file..\n",$log_file);
				createTables($cfg_hash->{'db_name'},$config_folder."/".$cfg_hash->{'db_schema'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_host'});
				#Fill the INFO table
				configure_database($cfg_hash);
		}else{
				log_and_exit("The database ".$cfg_hash->{'db_name'}." does not exist. Please use --create_database to make a newer one..\n",$log_file);
		}
		exit;
	}#Otherwise the database must exist
	else{
		#Check if exists
		if ( db_present($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'}) == 0){
			log_and_exit("The database ".$cfg_hash->{'db_name'}." does not exist. Please use --create_database to make a newer one..\n",$log_file);
		}
	}
	
	##If there are no utils to execute, uses the samplesheet to generate 
	#a configuration for the execution of the analyses
	if ( $do_pipeline_tasks == 1){
		if ($sample_sheet ne ''){
			#print_and_log("Checking and reading the sample sheet..\n",$log_file);		
			check_sample_sheet($sample_sheet,$cfg_hash->{'ss_separator'},$config_folder."/".$cfg_hash->{'ss_head_file'},$cfg_hash->{'ss_head_sep'});
			#Fills sample, group and analsysis table
			sample_sheet_2_db($sample_sheet,$cfg_hash->{'ss_separator'});
			$cfg_hash->{'sample_sheet'} = $sample_sheet;
			#If the user does not choose groups to run.. take them all
			if ($run_analyses eq ''){
				get_samples_exec_config("all",$cfg_hash->{'ss_id'});
				#print_and_log("Samples to execute:\n",$log_file);#DEBUGCODE
				#print_hash($exec_hash)#DEBUGCODE
			}
			#log_and_exit("The sample sheet has now been stored. Restart without the parameter. Exiting...\n",$log_file);
		}else{
			#if run_analyses and run_samples are used the configuration is modified
			#print_and_log("Checking samples to analyze..\n",$log_file);#DEBUGCODE	
			if ($run_analyses ne ''){
				get_samples_exec_config("input");
				#print_and_log("Samples to execute:\n".print_hash($exec_hash),$log_file);#DEBUGCODE
			}elsif ($run_samples ne ''){
					log_and_exit("Cannot start analyzing samples without defining an analysis. Please use --run_analyses as well. Exiting...\n",$log_file);
					#CHECK THE SAMPLES GIVEN. ALSO THE FACT THAT THEY MUST BE SEPARATED BY COMMA #TOCHECK
			}elsif ($sample_sheet eq '') {
				log_and_exit("Cannot start an analysis without the sample sheet or the analysis id. Exiting...\n",$log_file);
			}		
		}
	}

	print_and_log("Organizing folders..\n",$log_file);	
	#Creates the needed folders for all the analyses
	make_analysis_folders();

	#Set the human reference for the alignment program
	#$cfg_hash->{'hum_ref_align'} =  $cfg_hash->{'ref_align_f'}."/".$cfg_hash->{'human_genome'};

	#Set main log and data folders
	$main_log_folder = $working_folder."/".$logFolder;
	$main_data_folder = $working_folder."/".$dataFolder;
	$cfg_hash->{'main_data_folder'} = $main_data_folder;

	#Set the folders for the data that the programs use
	$cfg_hash->{'references_fold'} = $main_data_folder."/".$cfg_hash->{'references_fold'};
	
	#Variables for GATK databases
	if (defined $cfg_hash->{'known_mills_link'} ){$cfg_hash->{'known_mills'} = extract_name($cfg_hash->{'known_mills_link'},'gz');}
	if (defined $cfg_hash->{'known_db_snp_link'} ){$cfg_hash->{'known_db_snp'} = extract_name($cfg_hash->{'known_db_snp_link'},'gz');}
	if (defined $cfg_hash->{'known_hapmap_link'} ){	$cfg_hash->{'known_hapmap'} = extract_name($cfg_hash->{'known_hapmap_link'},'gz');}
	if (defined $cfg_hash->{'known_1000g_link'} ){$cfg_hash->{'known_1000g'} = extract_name($cfg_hash->{'known_1000g_link'},'gz');}
	if (defined $cfg_hash->{'known_omni_link'} ){$cfg_hash->{'known_omni'} = extract_name($cfg_hash->{'known_omni_link'},'gz');}

	if (defined $cfg_hash->{'known_sites_mills_link'} ){$cfg_hash->{'known_sites_mills'} = extract_name($cfg_hash->{'known_sites_mills_link'},'gz');}
	if (defined $cfg_hash->{'known_sites_db_snp_link'} ){$cfg_hash->{'known_sites_db_snp'} = extract_name($cfg_hash->{'known_sites_db_snp_link'},'gz');}
	if (defined $cfg_hash->{'known_sites_hapmap_link'} ){$cfg_hash->{'known_sites_hapmap'} = extract_name($cfg_hash->{'known_sites_hapmap_link'},'gz');}
	if (defined $cfg_hash->{'known_sites_1000g_link'} ){$cfg_hash->{'known_sites_1000g'} = extract_name($cfg_hash->{'known_sites_1000g_link'},'gz');}
	if (defined $cfg_hash->{'known_sites_omni_link'} ){$cfg_hash->{'known_sites_omni'} = extract_name($cfg_hash->{'known_sites_omni_link'},'gz');	}
	if (defined $cfg_hash->{'known_sites_1000gP3_link'} ){$cfg_hash->{'known_sites_1000gP3'} = extract_name($cfg_hash->{'known_sites_1000gP3_link'},'gz');}
			
	
	#Check if directory exists, otherwise it creates it
	unless(-d $main_log_folder){
		print_and_log("Creating folder $main_log_folder...\n",$log_file);
		mkdir $main_log_folder or die "ERROR: can't create folder $main_log_folder. Check permissions. \n";
	 }
	#Create the main data folder as a subfolder of the working dir
	unless(-d $main_data_folder){
		print_and_log("Creating folder $main_data_folder...\n",$log_file);
		mkdir $main_data_folder or die "ERROR: can't create folder $main_data_folder. Check permissions. \n";
	 }
	#Create the target per-chromosome folder as a subfolder of the main data folder
	my $target_perchr_fold = $main_data_folder."/".$cfg_hash->{'target_perchr_fold'};
	unless(-d $target_perchr_fold){
		print_and_log("Creating folder $target_perchr_fold...\n",$log_file);
		mkdir $target_perchr_fold or die "ERROR: can't create folder $target_perchr_fold. Check permissions. \n";
	 }
	#Create the gene db folder as a subfolder of the main data folder
	my $gene_db_fold = $main_data_folder."/".$cfg_hash->{'gene_db_fold'};
	unless(-d $gene_db_fold){
		print_and_log("Creating folder $gene_db_fold...\n",$log_file);
		mkdir $gene_db_fold or die "ERROR: can't create folder $gene_db_fold. Check permissions. \n";
	 }
	#Create the gene panels folder as a subfolder of the main data folder
	my $gene_panels_fold = $main_data_folder."/".$cfg_hash->{'gene_panels_fold'};
	unless(-d $gene_panels_fold){
		print_and_log("Creating folder $gene_panels_fold...\n",$log_file);
		mkdir $gene_panels_fold or die "ERROR: can't create folder $gene_panels_fold. Check permissions. \n";
	 }

	#Create the references folder as a subfolder of the main data folder
	my $references_fold = $cfg_hash->{'references_fold'};
	unless(-d $references_fold){
		print_and_log("Creating folder $references_fold...\n",$log_file);
		mkdir $references_fold or die "ERROR: can't create folder $references_fold. Check permissions. \n";
	 }
	 
	#REFERENCE DATABASES FOLDERS. 
	#The users can define their own folders into the user_config file. If defined there
	#default folders won't be used.
	if ( ! (defined $cfg_hash->{'ref_genome_f'}) ){
		$cfg_hash->{'ref_genome_f'} = $cfg_hash->{'references_fold'}."/".$cfg_hash->{'ref_genome_f_def'};
	}
	#Set the human reference for all the programs
	$cfg_hash->{'hum_ref'} =  $cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'};	
	if ( ! (defined $cfg_hash->{'target_reg_f'}) ){
		$cfg_hash->{'target_reg_f'} = $cfg_hash->{'references_fold'}."/".$cfg_hash->{'target_reg_f_def'};
	}
	if ( ! (defined $cfg_hash->{'annovar_db_f'}) ){
		$cfg_hash->{'annovar_db_f'} = $cfg_hash->{'references_fold'}."/".$cfg_hash->{'annovar_db_f_def'};
	}
	if ( ! (defined $cfg_hash->{'gatk_ref_f'}) ){
		$cfg_hash->{'gatk_ref_f'} = $cfg_hash->{'references_fold'}."/".$cfg_hash->{'gatk_ref_f_def'};
	}		

	#Set scratch and storage folders
	if ( ! (defined $cfg_hash->{'scratch_f'}) ){
		$cfg_hash->{'scratch_f'} = $main_data_folder."/".$cfg_hash->{'scratch_f_def'};
	}
	if ( ! (defined $cfg_hash->{'storage_f'}) ){
		$cfg_hash->{'storage_f'} = $main_data_folder."/".$cfg_hash->{'storage_f_def'};
	}	

	 
	#Create the targets folder as a subfolder of the references folder
	my $targets_fold = $cfg_hash->{'target_reg_f'};
	unless(-d $targets_fold){
		print_and_log("Creating folder $targets_fold...\n",$log_file);
		mkdir $targets_fold or die "ERROR: can't create folder $targets_fold. Check permissions. \n";
	 }
	#Create the genome folder as a subfolder of the references folder
	my $genome_fold = $cfg_hash->{'ref_genome_f'};
	unless(-d $genome_fold){
		print_and_log("Creating folder $genome_fold...\n",$log_file);
		mkdir $genome_fold or die "ERROR: can't create folder $genome_fold. Check permissions. \n";
	 }
	#Create the targets folder as a subfolder of the references folder
	my $gatk_dbs_fold = $cfg_hash->{'gatk_ref_f'};
	unless(-d $gatk_dbs_fold){
		print_and_log("Creating folder $gatk_dbs_fold...\n",$log_file);
		mkdir $gatk_dbs_fold or die "ERROR: can't create folder $gatk_dbs_fold. Check permissions. \n";
	 }
	#Create the folders for annovar databases
	my $annovar_db_f = $cfg_hash->{'annovar_db_f'};
	unless(-d $annovar_db_f){
		print_and_log("Creating folder $annovar_db_f...\n",$log_file);
		mkdir $annovar_db_f or die "ERROR: can't create folder $annovar_db_f. Check permissions. \n";
	 }

	#Scratch and storage folder
	my $scratch_f = $cfg_hash->{'scratch_f'};
	unless(-d $scratch_f){
		print_and_log("Creating folder $scratch_f...\n",$log_file);
		mkdir $scratch_f or die "ERROR: can't create folder $scratch_f. Check permissions. \n";
	 }
	#Scratch and storage folder
	my $storage_f = $cfg_hash->{'storage_f'};
	unless(-d $storage_f){
		print_and_log("Creating folder $storage_f...\n",$log_file);
		mkdir $storage_f or die "ERROR: can't create folder $storage_f. Check permissions. \n";
	 }

	#Folder or WebServer for results. If not defined is the working folder
	if ( ! defined $cfg_hash->{'html_host'} ){
			$cfg_hash->{'html_host'} = $working_folder;
	}

	 	 	 	 	 	 	 	 	 
	#Set the dependence hash file
	$deps_f = $main_data_folder."/".$cfg_hash->{'jc_dep_hash_f'};
	
	#Put in the config hash a variable with all the tasks to be executed
	$cfg_hash->{'tasks_to_exec'} = "";
	
	#Shift the log file in the LOG folder
	move($log_file,$main_log_folder."/".$log_file);
	$log_file = $main_log_folder."/".$log_file;

	##Store the config hash in log_file
	#print_and_log("Storing the configuration hash in $config_file..\n",$log_file);#DEBUGCODE
	save_hash(\$cfg_hash,$config_file);
	print_and_log($configClosure,$log_file);
}

#This function sets global variables into the hash
#I wrote it because some global are not re-set into the get modified config hash
sub set_other_global_var_into_hash {
	my $cfg_hash_updated = shift;
	$cfg_hash_updated->{'hum_ref'} =  $cfg_hash_updated->{'ref_genome_f'}."/".$cfg_hash_updated->{'hum_gen_align'};	
}


=head2 get_samples_exec_config

 Title   : get_samples_exec_config
 Usage   : get_samples_exec_config(   );

 Function:  Extracts the samples to be executed in correspondence to each
						group. The input should be 
							groups: 1,2,4-8 (comma separates groups)
							samples: 1,4,5-10;2,3,5-7;3 (semicolon separates samples associated to groups)
						ORDER IS IMPORTANT!!
 
 Returns : nothing

=cut
sub get_samples_exec_config{
			my $mode = shift;
			my $ss_id = shift;
			
			my $sep = $cfg_hash->{'parameters_sep'};
			
			#Mode getting all the samples for the analysis given by the samplesheet	
			if ( $mode eq 'all' ){
				#I inizialize a parameter with the samples that must be analyzed
				#This parameter is used during the VQSR to define the samples that must be given in input
				#into the output file. If all samples are run I write: ALL
				$cfg_hash->{'run_samples'} = "ALL";				
				#Here we get the group ids associated with the analysis id
				my $groups_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'}
																,$cfg_hash->{'db_sample_sheet_id'},$ss_id);
				
				foreach my $analysis_id ( keys %{$groups_associated}){
					#Here we get the samples ids associated with the group id
					my $samples_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'}
																,$cfg_hash->{'db_analysis_id'},$analysis_id);						
					foreach my $sample_id ( keys %{$samples_associated}){
						my $readf_id_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
																,$cfg_hash->{'db_sample_id'},$sample_id);						
							foreach my $readf_id ( keys %{$readf_id_associated}){
								if ( defined $sample_id ){
									if (! defined $exec_hash->{$analysis_id}->{$sample_id} ){
										$exec_hash->{$analysis_id}->{$sample_id} = $readf_id.",";
									}else{
										$exec_hash->{$analysis_id}->{$sample_id} .= $readf_id.",";
									}
								#print "Analysis = $analysis_id sample = $sample_name readf_ids =".$exec_hash->{$analysis_id}->{$sample_name}." \n";#DEBUGCODE
								}
						}
					}
					#print "Removing the last comma from the exec_hash...\n";#DEBUGCODE
					while (my ($key, $value) = each %{ $exec_hash->{$analysis_id} } ) {
						#print "$key = $value \n";
						chop($exec_hash->{$analysis_id}->{$key});
					}
				}
			}else{
				print_and_log("Getting samples from input..\n",$log_file);	
				#Mode getting all the samples for the analysis given by the user in input	
					#Groups that the user wants to run
					my @exec_analyses = split($sep,separate_input_ids($run_analyses,$sep));

					my $exists = -1;

					print_and_log("Analyses to run: ",$log_file);#DEBUGCODE
					print_array(\@exec_analyses);
					print_and_log("\n",$log_file);#DEBUGCODE
					#Separate the samples using the ;
					my @given_samples = ();
					if ($run_samples ne ''){
						#print "Analyzed samples: $run_samples\n";#DEBUGCODE
						my $sam_sep =  $cfg_hash->{'param_groups_sep'};
						@given_samples = split($sam_sep,$run_samples);
						#print "Splitted: ";#DEBUGCODE
						#print_array(\@given_samples);#DEBUGCODE
					}
					my $group_num = 0;							
					foreach my $analysis_id (@exec_analyses){
							#Here we check if the sample sheet has been loaded and get the analysis id
							if ($group_num == 0){
								#print_and_log("Getting the analysis id...\n",$log_file);#DEBUGCODE
								$ss_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sample_sheet_id'},
																				$cfg_hash->{'db_analysis_id'},$analysis_id);
							}
							#Samples present in the db associated to a given analysis
							my @samples_present = ();
							#print_and_log("Checking the presence of analysis $analysis_id...\n",$log_file);#DEBUGCODE
							#Check if the group really exists
							$exists = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'}
																				,$cfg_hash->{'db_analysis_id'},$analysis_id);
							#Only if the group id exists goes on inserting the analysis. Otherwise this analysis is not taken													
							if ($exists ne -1){
								#Here we get the samples ids associated with that group
								my $samples_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'}
																				,$cfg_hash->{'db_analysis_id'},$analysis_id);
							#	print_and_log("Samples present for analysis $analysis_id\n",$log_file);#DEBUGCODE
								foreach my $key ( keys %{$samples_associated}){
										#print_and_log(" $key ",$log_file);
										push(@samples_present,$key);
								}
								#Finally we get the samples to analyze for each group if they are present
								if ($run_samples ne ''){
									#I inizialize a parameter with the samples that must be analyzed
									#This parameter is used during the VQSR to define the samples that must be given in input
									#into the output file
									$cfg_hash->{'run_samples'} = $run_samples;
									
									#Get the samples associated to the group
									#print_and_log("Separating: ".$given_samples[$group_num]."\n",$log_file);#DEBUGCODE
									my $group_samples = separate_input_ids($given_samples[$group_num],$sep);
									#print "Obtained: $group_samples\n ";#DEBUGCODE
									my @group_samples = split($sep,$group_samples);
									#print_and_log("Searching for samples: \n",$log_file);#DEBUGCODE
									#print_array(\@group_samples);#DEBUGCODE
									foreach my $sample_id (@group_samples ){
											if ( grep {/\b$sample_id\b/} @samples_present){
												my $readf_id_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}
																						,$cfg_hash->{'db_sample_id'},$sample_id);						
													foreach my $readf_id ( keys %{$readf_id_associated}){
														if ( !(defined $exec_hash->{$analysis_id}->{$sample_id}) ){
															$exec_hash->{$analysis_id}->{$sample_id} = $readf_id.",";
														}else{
															$exec_hash->{$analysis_id}->{$sample_id} .= $readf_id.",";
														}
														#print "analysis = $analysis_id sample = $sample_id readf_ids =".$exec_hash->{$analysis_id}->{$sample_id}." \n";#DEBUGCODE
													}
											}else{
													print_and_log("Sample $sample_id is not present for analysis $analysis_id. It will not be run..\n",$log_file);
											}
									}
								}else{
									#I inizialize a parameter with the samples that must be analyzed
									#This parameter is used during the VQSR to define the samples that must be given in input
									#into the output file. If all samples are run I write: ALL
									$cfg_hash->{'run_samples'} = "ALL";									
									#Otherwise all the samples are considered
									foreach my $sample_id ( @samples_present ){
										my $readf_id_associated = select_all_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
																				$cfg_hash->{'db_sample_id'},$sample_id);		
										foreach my $readf_id ( keys %{$readf_id_associated}){
											if (! defined $exec_hash->{$analysis_id}->{$sample_id} ){
												$exec_hash->{$analysis_id}->{$sample_id} = $readf_id.",";
											}else{
												$exec_hash->{$analysis_id}->{$sample_id} .= $readf_id.",";
											}
										 #print "analysis = $analysis_id sample = $sample_id readf_ids =".$exec_hash->{$analysis_id}->{$sample_id}." \n";#DEBUGCODE
										}
									}
								}
							#Remove the last comma from all
							#print "Removing the last comma from the exec_hash...\n";#DEBUGCODE
							while (my ($key, $value) = each %{ $exec_hash->{$analysis_id} } ) {
								#print "$key = $value \n";
								chop($exec_hash->{$analysis_id}->{$key});
							}
							
							}else{
								print_and_log("WARNING: the analysis_id $analysis_id does not exist in the database. It will not be analyzed..\n",$log_file);
							}
						#Next analysis index
						$group_num++;
					}
			}		
}

=head2 configure_database

 Title   : configure_database
 Usage   : configure_database(   );

 Function:  db configuration
 
 Returns : nothing

=cut
sub configure_database{

	##################################
	#INFO TABLE
	####################################
	#Set the information into the info table
	#Insert into table info a line useful as a semaphore
	my $fields = $cfg_hash->{'db_busy'}.",".$cfg_hash->{'db_info_row'}.",".$cfg_hash->{'db_analysis_id'}.",".
	$cfg_hash->{'db_sample_sheet_id'}.",".$cfg_hash->{'db_dbversion'}.",".$cfg_hash->{'db_progver'};
	my $values = "1,1,-1,-1,'".$cfg_hash->{'VarGeniusVer'}."','".$cfg_hash->{'dbVersion'}."'";		
	#Inserts
	insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$fields,$values);
}



=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
						If you want to add an utility you need to add:
						
						- new_u|--new_utility, among the description in $howToUse
						- new_u|new_utility=s"=> \$new_utility, among the GetOptions(..)
						- push the element into global array @utils_to_exec
							if ($new_utility){push(@utils_to_exec,"new_utility");}
						- increase the utilities to exec counter in 
							if ( $new_utility ne ''){$new_utility++;}
						- go to "#Global variables" and add a new variable with the same name at the start of this script
							my $new_utility = '';
						- go to "#Utilities execution" and add the "IF"  checking if this command is active
							and containing the function to call
							
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

	
	my $howToUse = "Use with: perl $program_name \n\n".
	"-c|--user_config: The configuration file with all the parameters needed for the execution. It should have been copied ".
									"into the working folder during the installation.\n".
  "-ss|--sample_sheet: You need to make the sample sheet indicating the information about your samples and analyses. Please ".
									" use the script get_sample_sheet.pl to this aim.\n".
  "-db_c|--create_database: Creates a database named as from the configuration file db_name variable. You need to run this ".
  " command the first time that you use $program_name. This command creates also the tables.\n".
  "-db_t|--create_tables: Fills the tables of the database named as from the configuration file db_name variable. You need to run this ".
								" command the first time that you use $program_name only if the database has been already created and you did not run ".
								" --create_database before. \n".
  "-gdb|--genedb: Fills the tables needed to store the information about the genes.\n".
  "-dlds|--dldatasets: Downloads the genome in fasta format;\n".
  "-gan|--geneann: Generates a table with annotation for all genes of the database.\n".
  "-ra|--run_analyses: Use this command to select the analysis identifiers that you want to start.\n".
  "-rs|--run_samples: Use this command to select the samples for which you want to run the analysis. Needs you also use --run_analyses \n".
  "-qc|--quality_check: Execute the quality check  steps given in the user configuration file.\n".
  "-trim|--trimming: Execute the trimming steps given in the user configuration file.\n".
  "-aln|--alignment: Executes the alignment steps given in the user configuration file\n".
  "-ref|--refining: Executes the refinement steps given in the user configuration file\n".
 # "-m|--mergesam: \n".
  "-vc|--varcall:  Executes the variant calling steps given in the user configuration file.\n".
  "-gen|--genot: \n".
  "-vf|--varfilt: Executes the variant filtering step.\n".
  #"-ph|--phasing: \n".
  "-st|--stats: Executes statistics about the coverage. Executes samtools flagstats to find what is the amount of aligned".
		" sequences and finds the non-covered regions of the target file. It also merges the bam files if they are separated per-lane".
		" for further analyses.\n".
  "-an_st|--analysis_stats: Executes GATK DepthOfCoverage using all the samples bam files to get the genes and interval coverage.\n".
  #"-ann|--annotate: DEPRECATED\n".#
  "-vcfi|--vcfimport: Import variants from VCF into the database. Needs the --varfilt step completed.\n".#
  #"-out|--finalout: So far it just imports the variants into the database. This step is needed to compute the internal frequency of variants.\n".
  "-fout|--fastfout: starts from the last VCF obtained. Takes all the variants, annotates them, builds an output, do some rearrangements on the output".
								"Creates plots and an HTML page with all the output.\n".
	"-start|--start_all: identical to run -qc -aln -ref -vc -vf -fout -st -an_st.\n".							
  #"-p|--pipeline\n".#
  "-oc|--overwrite_conf: Use this parameter to change only for this time any of the configuration file variable. Usage: ".
				" --overwrite_conf 'parameter1=value1]---[parameter2=value2]---[parameter3=value3' ".
  "-prof|--profile: use this parameter to change the configuration parameters using an additional file. Please put the ".
		"profile file into the CONFIGURATION folder of VarGenius and use only the file name.".
  "-reann|--reannote: used to reannote all the variants from the database using Annovar in a complete file.\n".
  "-upfreqs|--update_freqs: used anytime the user wants to update the frequencies. Needs less than one hour.\n".
  "-remss|--remove_sample_sheet: remove completely from the database groups from a sample sheet given the sample_sheet id or its name,\n".
  "-reman|--remove_analyses: Removes all the database data related with the analyses ids given in input\n".
  "-store|--store_analyses: Makes a tar of the output of the analyses given in input and stores it into the folder storage_f that you set into the configuration file.\n ".
  "-save|--save_analyses: Saves the output of the analyses given in input into the folder storage_f that you set into the configuration file and makes symbolic links.\n ".
  "-sh_a|--show_analyses: given the analysis id (positive integer) or the name of the sample sheet returns the analyses contained\n".
  "-sh_s|--show_samples: given the analysis id (positive integer) or the analysisname returns the samples contained\n".
  "-aid|--get_analyses_ids: given an analysis name or a set of analyses names separated by comma, returns the analyses ids\n".
  "-sid|--get_samples_ids: given a sample name or a set of sample names separated by comma, returns the sample ids\n".
  "-clean|--clean_analyses: DEPRECATED\n".
  "-clf|--clean_folder: removes SAM and BAM useless files from the alignment_out folder, all the files from refine_out and varcall_out for the given groups.\n".
  "-fr|--force_run: \n".
  "-lic|--licence: shows GNU GPL Open Source warranty and conditions.\n\n";
  
  #  Parse options
  GetOptions(
           "help" =>        \$HELP,
           "version" =>      \$VERSION,
           "license" =>      \$license,
           "c|user_config=s" => \$user_config,   #It's mandatory
           "ss|sample_sheet=s" => \$sample_sheet,   #It's mandatory
           "db_c|create_database" => \$create_database,
           "db_t|create_tables" => \$create_tables, 
			 "gdb|genedb" => \$genedb,
			 "dlds|dldatasets" => \$dldatasets,
			 "gan|geneann" => \$geneann, 
           "ra|run_analyses=s" => \$run_analyses,
           "rs|run_samples=s" => \$run_samples,
           "qc|quality_check" => \$qualityCheck,
           "trim|trimming" => \$trimming,
           "aln|alignment" => \$alignment,
           "ref|refining" => \$refining,
          # "m|mergesam" => \$mergesam,
           "vc|varcall" => \$varcall,
           "gen|genot" => \$genot,
           "vf|varfilt" => \$varfilt,
           #"ph|phasing" => \$phasing,
           "st|stats" => \$stats,
           "an_st|analysis_stats" => \$anstats,
           "ann|annotate" => \$annotate,
		 "fout|fastfout" => \$fastfout,
		 "vcfi|vcfimport" => \$vcfimport,
		 "cnv|cnv_detection" => \$cnv_detection,
		 "cnvo|cnv_out" => \$cnv_out,
		 "start|start_all" => \$start_all,
		 "p|pipeline=s" => \$pipeline, 
		 "oc|overwrite_conf=s" => \$overwrite_conf, 
		 "prof|profile=s" => \$profile, 
		 "reann|reannote"  => \$reannote, 
		 "freqs|update_freqs\n"  => \$update_freqs, 
		 "remss|remove_sample_sheet=s" => \$remove_sample_sheet,
		 "reman|remove_analyses=s" => \$remove_analyses,
		 "store|store_analyses=s" => \$store_analyses,
		  "vgstore=s" => \$vg_store_analyses,
		 "save|save_analyses=s" => \$save_analyses,
		 "sh_a|show_analyses=s" => \$show_analyses,
		 "sh_s|show_samples=s" => \$show_samples,
		 "aid|get_analyses_ids=s" => \$get_analyses_ids,
		 "sid|get_samples_ids=s" => \$get_samples_ids,
	     "clean|clean_analyses"=> \$clean_analyses,
		 "clf|clean_folder=s"=> \$clean_folder,
           "fr|force_run" => \$force_run
            );

  #Print a little help
  if ( $HELP ){
    print $howToUse;
    pod2usage(1);
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "Version: $program_name $program_version \n";
    exit;
  }
  
	#-start executes a complete standard pipeline with all tasks
	#equivalent to start: qc, trim, aln, ref, vc, vf, fout, st, an_st, vcfi
	if ($start_all){
			$qualityCheck = 1;
			$trimming = 1;
			$alignment = 1;
			$refining = 1;
			$varcall = 1;
			$varfilt = 1;
			$stats = 1;
			$anstats = 1;
			$fastfout = 1;
			$vcfimport = 1;
	}
	#To check that all the tasks chosen are consecutive, I match with a regex the
	#string composed of all the tasks with the string formed by the tasks
	#wanted
	my $all_tasks_str = "quality_check_trimming_alignment_refining_varcall_varfilt_fastfout_vcfimport";
	my $formed_str = "";
	my $doanalysis = 1;
	#Get the tasks that will be executed
	if ($qualityCheck){push(@tasks_to_exec,"quality_check");$formed_str .= "quality_check"."_";}
	if ($trimming){push(@tasks_to_exec,"trimming");$formed_str .= "trimming"."_";}
	if ($alignment){push(@tasks_to_exec,"alignment");$formed_str .= "alignment"."_";}
	if ($refining){push(@tasks_to_exec,"refining");$formed_str .= "refining"."_";}
	if ($mergesam){push(@tasks_to_exec,"MS");}
	if ($varcall){push(@tasks_to_exec,"varcall");$formed_str .= "varcall"."_";}
	if ($genot){push(@tasks_to_exec,"genot");}
	if ($varfilt){push(@tasks_to_exec,"varfilt");$formed_str .= "varfilt"."_";}
	#if ($phasing){push(@tasks_to_exec,"phasing");}
	if ($stats){push(@tasks_to_exec,"stats");$doanalysis = 0;}
	if ($anstats){push(@tasks_to_exec,"anstats");$doanalysis = 0;}
	if ($annotate){push(@tasks_to_exec,"annotate");}
	if ($fastfout){push(@tasks_to_exec,"fastfout");$formed_str .= "fastfout"."_";}
	#if ($finalout){push(@tasks_to_exec,"finalout");$formed_str .= "finalout"."_";}
	if ($vcfimport){push(@tasks_to_exec,"vcfimport");$formed_str .= "vcfimport"."_";}
	if ($cnv_detection) {push(@tasks_to_exec,"cnv_detection");$doanalysis = 0;}
	if ($cnv_out) {push(@tasks_to_exec,"cnv_out");$doanalysis = 0;}

	if ($create_database){
		$doanalysis = 0;
		#If the reannotation is run, only the that is run..
		$do_pipeline_tasks = 0;		
	}
	if ($create_tables){
		$doanalysis = 0;
		#If the reannotation is run, only the that is run..
		$do_pipeline_tasks = 0;		
	}	
	if ($genedb){
		push(@tasks_to_exec,"genedb");
		$doanalysis = 0;		#If the genedb is run, only the that is run..
		$do_pipeline_tasks = 0;
	}
	if ($dldatasets){
		push(@tasks_to_exec,"dldatasets");
		$doanalysis = 0;		#If the dldatasets is run, only the that is run..
		$do_pipeline_tasks = 0;
	}	
	if ($geneann){
		push(@tasks_to_exec,"geneann");
		$doanalysis = 0;		#If the geneann is run, only the that is run..
		$do_pipeline_tasks = 0;
	}	
	if ($reannote){
		push(@tasks_to_exec,"reannote");
		$doanalysis = 0;
		#If the reannotation is run, only the that is run..
		$do_pipeline_tasks = 0;		
	}
	if ($update_freqs){
		push(@tasks_to_exec,"update_freqs");
		$doanalysis = 0;
		#If the frequencies update is run, only the that is run..
		$do_pipeline_tasks = 0;	
		}
	if ($clean_analyses){
		push(@tasks_to_exec,"clean");
		$doanalysis = 0;
		#If the user wants to make cleaning does just that..
		$do_pipeline_tasks = 0;	
		}
	

	#Utilities that will be executed
	my @utils_to_exec = ();
	if ($create_database){push(@utils_to_exec,"create_database");}
	if ($create_tables){push(@utils_to_exec,"create_tables");}
	if ($show_analyses){push(@utils_to_exec,"show_analyses");}
	if ($show_samples){push(@utils_to_exec,"show_samples");}
	if ($get_analyses_ids){push(@utils_to_exec,"get_analyses_ids");}
	if ($remove_sample_sheet){push(@utils_to_exec,"remove_sample_sheet");}
	if ($remove_analyses){push(@utils_to_exec,"remove_analyses");}
	if ($store_analyses){push(@utils_to_exec,"store_analyses");}
	if ($vg_store_analyses){push(@utils_to_exec,"vg_store_analyses");}
	if ($save_analyses){push(@utils_to_exec,"save_analyses");}
	if ($clean_folder){push(@utils_to_exec,"clean_folder");}
	
	#Checks tasks to execute	
	if (scalar(@tasks_to_exec) ){
		#Check if taskVCFBs to execute are all consecutives and there are more than one to execute
		chop($formed_str);
		if ( $all_tasks_str !~ $formed_str and $doanalysis == 1  and scalar(@tasks_to_exec) > 1){
		  print "This command is not permitted in $program_name. Please choose all consecutive steps. Exiting... \n";
		  exit;	
		}		
		print "Tasks that will be run: \n";
		foreach my $ex_task (@tasks_to_exec){
			print "\t $ex_task\n";
		}
	}
	if ( scalar(@tasks_to_exec) == 0  and  scalar(@utils_to_exec) == 0){
		print "You did not select a task to execute. No analysis will be executed...\n";
		print $howToUse;
	}
	
	#Execution of utilities
	if ( $remove_sample_sheet ne ''){$exec_utils++;}
	if ( $remove_analyses ne ''){$exec_utils++;}
	if ( $store_analyses ne ''){$exec_utils++;}
	if ( $vg_store_analyses ne ''){$exec_utils++;}
	if ( $save_analyses ne ''){$exec_utils++;}
	if ( $show_analyses ne ''){$exec_utils++;}
	if ( $show_samples ne ''){$exec_utils++;}
	if ( $clean_folder ne ''){$exec_utils++;}
	if ( $get_analyses_ids ne ''){$exec_utils++;}
	if ( $get_samples_ids ne ''){$exec_utils++;}
}


