#!/usr/bin/perl

#VarGenius utilities
#This script is a template to access the output of VarGenius programmatically
#It has many different functions depending on which data should be used and which operation must be completed.
#
#perl VarGeni -f BAM_COPY_FROM_ARCHIVE 

use strict;
use warnings;
use Getopt::Long;
use File::Copy;#To manage files

#vargenius.pl path
use lib '/pico/work/TELET_TIGEM/VarGeniusBeta/';

#Using a library for database management
use LIB::db_management qw(get_id_if_exists_from_db do_query_select_all get_id_if_exists_from_db_woconn
												fetch_all_rows get_count_from_selected do_fetch_row_array_woconn
												do_query_select_all_woconn fetch_all_rows_woconn
												do_fetch_row_array do_query_select_all_woconn getSampleConfiguration_locked);
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash initialize_folders print_and_log separate_input_ids get_output_name_from_executed);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present
				download_file dl_and_extract delete_file extract_columns_from_file
				append_str_2_file check_presence file_num_rows 
				invert_cols_position get_col_index insert_col_in_file_table delete_columns
				extract_col_from_file file_name_wostrange_chars append_hash_to_file
				delete_rows_containing file_list_to_array extract_colnum_from_file_linux
				shift_column compress_folder extract_name append_str_2_file_if_path_notexist
				join_files_with_cat separate_elements_on_col join_vcf_files_with_cat
				vcf_insert_string_into_field vcf_filter_identical_calls
				overwrite_str_in_file_if_exists tsv_2_xls list_to_array
				count_lines_file compid_intervals_2_bed insert_line_into_file);

#Global variables
my $foldersFile = "folders.txt";
my $program_name = "vargenius_utilities.pl";
my $config_folder = "CONFIGURATION";
my $dataFolder = "DATA";
my $program_version = "0.1";
my $storage_path = 'http://UDProject:udpNGSd4t4@130.186.13.104/udp/';
#The steps array is used both to construct the input and output file names
#and to have a reference for the type of pipeline to use

#Input parameters
my $func = "";#function to execute			
my $parameters = "";								
my $user_config = "";
my $input_f = "";
my $suppl_file = "";
my $out_file = "";
my $log_file = "";

#Steps array to build output names
my @steps_array;
my @steps_array2;
my @steps_array3;
my @stats_steps_array;


						
#Get the working folder and vargenius	
#print "Inizialyzing folders ...\n";
my ($working_folder,$program_folder) = initialize_folders($foldersFile);
my $program_config = $program_folder."/".$config_folder."/program_config.txt";
my $main_data_folder = $working_folder."/".$dataFolder;

parse_command_line_args();
my $cfg_hash;

#print "Reading configuration $program_config...\n";
configFile2Hash($program_config,\$cfg_hash);
#print "Reading configuration $user_config...\n";
configFile2Hash($user_config,\$cfg_hash);

#Fill the steps arrays using the $pipeline mode given in input
#print "steps array ...\n";
fill_steps_arrays();
#Get parameters from input
#print "parse command line ...\n";




# for the analysis given, copy the BAM files from the Archive folder
#to the working folder
if ( $func eq 'BAM_COPY_FROM_ARCHIVE'){
	#print "Starting to copy files for IGV session for $parameters ...\n";
	bam_copy_from_archive($parameters);
}


#Copies the files related to the IGV session from the storage folder to the 
#working folder
sub bam_copy_from_archive{
	my $analysis_ids = shift;
	
	die "ERROR: analysisid is needed to get files \n" unless ($parameters ne '');
	die "ERROR: storage folder is needed \n" unless ($suppl_file ne '');

	#Get an array of group ids to store
	my $sep = $cfg_hash->{'parameters_sep'};
	my @ans_to_store = split($sep,separate_input_ids($analysis_ids,$sep));

	#Execute for each analysisid to store	 	
	foreach my $analysis_id (@ans_to_store){
		
		#Obtain the analysis name from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
															$cfg_hash->{'db_analysis_id'},$analysis_id);
		#print "Starting to copy files for IGV session for $analysis_name ...\n";														
		my $vcfFolder = $working_folder."/".$analysis_name."/".$cfg_hash->{'varfilt_out_f'}."/";
		my $bamFolder = $working_folder."/".$analysis_name."/".$cfg_hash->{'refine_out_f'}."/";
		
		#Get the samples ids involved for the analysis
		my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
		#print "Executing: $query\n";#DEBUGCODE
		my $res_group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
		my @group_sam = ();
		
		my $storage_f = $suppl_file;
		
		foreach my $sample_id (keys %{$res_group_sam}){
			
			my $params;
			#######SAMPLE LEVEL
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
					
			my $recal_bam = $storage_f."/".$analysis_name."/".get_output_name_from_executed($params,$cfg_hash->{'print_reads_step'},
										$params->{$cfg_hash->{'db_sample_name'}},\@steps_array2).".".$cfg_hash->{'bam_ext'};
			my $recal_bai = extract_name($recal_bam,"noext").".".$cfg_hash->{'bai_ext'};
			print "Sample ".$params->{$cfg_hash->{'db_sample_name'}}.":\n copying $recal_bam to $bamFolder...\n";
			copy($recal_bam,$bamFolder) or print "WARNING: Cannot copy the BAM file ($recal_bam) to $bamFolder!";
			
			print "Sample ".$params->{$cfg_hash->{'db_sample_name'}}.":\n copying $recal_bai to $bamFolder...\n";
			copy($recal_bai,$bamFolder) or print "WARNING: Cannot copy the BAM file ($recal_bai) to $bamFolder!";
			#Given the analysis id copies the BAM files and the VCF from STORAGE into the WORK folder

		}					
			
		my $params;
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
		$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
			 
		my $vcf_file = $storage_f."/".$analysis_name."/".get_output_name_from_executed($params,$cfg_hash->{'phasing_step'},
									$params->{$cfg_hash->{'db_analysis_name'}},\@steps_array3).".".$cfg_hash->{'vcf_ext'};	
		my $vcfi_file = extract_name($vcf_file,"noext").".".$cfg_hash->{'vcfi_ext'};									
		print "copying $vcf_file to $vcfFolder...\n";
		copy($vcf_file,$vcfFolder) or print "WARNING: Cannot copy the VCF file ($vcf_file) to $vcfFolder!";
		print "copying $vcfi_file to $vcfFolder...\n";
		copy($vcfi_file,$vcfFolder) or print "WARNING: Cannot copy the VCF file ($vcfi_file) to $vcfFolder!";
	}	
}

=head2 fill_steps_arrays

 Title   : fill_steps_arrays
 Usage   : fill_steps_arrays(   );

 Function: Fills the steps array using the type of pipeline that must be 
					executed
					
 Returns : 
=cut
sub fill_steps_arrays{
	#Get the parameters from input
	#my $cfg_hash = load_hash($config_file);
	
	my $comma = ",";
	my $pipe_ind = 1;
	#Split the variables from the configuration file and fill the arrays
	@steps_array = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st1'});
	@steps_array2 = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st2'});
	@steps_array3 = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st3'});
	@stats_steps_array = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_stats'});
	#print "LEGGO:" .$cfg_hash->{'pipe'.$pipe_ind.'_st2'} ;
	#print "LEGGO:" .$cfg_hash->{'pipe'.$pipe_ind.'_st3'};
	#print "LEGGO: $cfg_hash->{'pipe'.$pipe_ind.'_st3'}";
	#print "LEGGO: $cfg_hash->{'pipe'.$pipe_ind.'_stats'}";
	
}

=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

	my $howToUse = "Use with: \nperl $program_name \n\n".
	"-c|--user_config: The configuration file with all the parameters needed for the execution. It should have been copied ".
									"into the working folder during the installation.\n".
	"\t-f|--func: functions: \n".
					"\t\tBAM_COPY_FROM_ARCHIVE (Copies the needed files for IGV session from the Storage folder to the WORK folder\n".
					"\t\t\t To use this function you need to give in input:\n".
	"\t\t\t -p: a list of analysisid\n".
          
	"\t-s|--suppl_file: supplementary file to use in the functions.\n".
	"\t-o|--out_file: output file to use in the functions .\n".
	"\t-p|--parameters: a generic argument to use parameters into the functions .\n".
	"\t-log|--log_file: log file to be used .\n".
	"\t-i|--input_file: input file to use\n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "version" => \$VERSION,
           "c|user_config=s" => \$user_config,   #It's mandatory
           "f|func=s" => \$func,   #It's mandator
           "i|input_file=s" => \$input_f,
           "s|suppl_file=s" => \$suppl_file,
           "p|parameters=s" => \$parameters,
           "log|log_file=s" => \$log_file,
           "o|out_file=s" => \$out_file
            );

  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "version: $program_version \n";
    exit;
  }
  
  #Func and user config are mandatory
  if ($user_config eq ''){
		print "configuration file in input is needed. Please use --user_config. Exiting... \n";
    exit;
	}

  #Func and user config are mandatory
  if ($func eq ''){
		print "function to use is needed. Please use --func. Exiting... \n";
    exit;
 }
}
