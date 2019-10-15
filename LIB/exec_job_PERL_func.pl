#!/usr/bin/perl

use lib '/pico/work/TELET_TIGEM/VarGeniusBeta/';

#Libraries
use strict;
use warnings;

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit );

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash delete_file );

#Using a library for database management
use LIB::db_management qw( get_id_if_exists_from_db );

####################
#
#This script is used to launch a PERL subroutine
# 
#N.B. You can insert whatever!

#Get the file with the parameters for the module
my $config_file = $ENV{CONFIG_FILE};
my $func = $ENV{FUNC};
#Parameters must be separated by semicolon (;)
my $params = $ENV{PARAMS};
#Gets the log file
my $log_file = $ENV{LOG_FILE};

my $parameters;

#CNV Pipeline
#This function  is useful for the CNV pipeline when XHMM jobs are executed and before to launch
# a job for the XHMM pipeline, this one is launched to check if the DOC jobs are completed and well done.
# It constructs the string of sample_interval_summary files to be merged together for XHMM
# into an output file given as parameter that will be used into the SHELL script
if ($func eq "CHECK_CNV_INPUT"){
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	#separate the parameters
	my @params = split(";",$params);				
	#Put parameters into the hash
	foreach my $param (@params){
		#split variable name and value
		my @vals = split("=",$param);
		$parameters->{$vals[0]} = $vals[1];
	}
	my @samples_ids = split(/-/,$parameters->{"samples"});
	print_and_log("Creating a job to search for available sample_interval_summary files..\n",$log_file);

	open (OUT,">".$parameters->{"out_f"}) or die "Cannot open $config_file\n";
	
	
	#search available sample_interval_summary files and construct the input string for both SEX and AUTO
	my $samp_int_sums_paths_sex = "";
	my $samp_int_sums_paths_auto = "";
	foreach my $sample_id (@samples_ids){
		#Get the analysis id and name for the given sample
		my $curr_anal_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
		my $curr_anal_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$curr_anal_id);	
		my $curr_sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
		
		#The input folder is given by the working folder  (vargenius_analyses) and the analysis name

		my $statsFolder = $cfg_hash->{'work_fold'}."/".$curr_anal_name."/".$cfg_hash->{'stats_out_f'};
				
		my $int_sum_f_auto = "$statsFolder/$curr_sample_name.XHMM.AUTO.DOC.".$cfg_hash->{'DOC_sam_interval_sum'};
		print_and_log("Searching for $int_sum_f_auto..\n",$log_file);
		if ( -e $int_sum_f_auto and !(-z $int_sum_f_auto) ){
			#print DOC_O_AUTO $int_sum_f_auto."\n";
			$samp_int_sums_paths_auto .= " --GATKdepths $int_sum_f_auto ";
			print_and_log("$int_sum_f_auto will be used..\n",$log_file);
			
			#Cleaning operation (only the interval_summary file is needed in XHMM, hence all other file can be removed)
			if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
				print "Removing DoC unused output..\n";	
				my @xhmm_remove_DOC_outputs = split(",",$cfg_hash->{'xhmm_remove_DOC_outputs'});
				foreach my $xhmm_remove (@xhmm_remove_DOC_outputs){
					delete_file("$statsFolder/$curr_sample_name.XHMM.AUTO.DOC.".$cfg_hash->{$xhmm_remove});	
				}
			}		
		}
		my $int_sum_f_sex = "$statsFolder/$curr_sample_name.XHMM.SEX.DOC.".$cfg_hash->{'DOC_sam_interval_sum'};
		print_and_log("Searching for $int_sum_f_sex..\n",$log_file);
		if ( -e $int_sum_f_sex and !(-z $int_sum_f_sex) ){
			#print DOC_O_SEX  $int_sum_f_sex."\n";
			$samp_int_sums_paths_sex .= " --GATKdepths $int_sum_f_sex ";
			print_and_log("$int_sum_f_sex will be used..\n",$log_file);
			
			#Cleaning operation (only the interval_summary file is needed in XHMM, hence all other file can be removed)
			if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
				print_and_log("Removing  DoC unused output..\n",$log_file);	
				my @xhmm_remove_DOC_outputs = split(",",$cfg_hash->{'xhmm_remove_DOC_outputs'});
				foreach my $xhmm_remove (@xhmm_remove_DOC_outputs){
					delete_file("$statsFolder/$curr_sample_name.XHMM.SEX.DOC.".$cfg_hash->{$xhmm_remove});	
				}
			}
		}		
	}
	#Print the two lines of paths into the output file
	print OUT $samp_int_sums_paths_auto."\n";
	print OUT $samp_int_sums_paths_sex."\n";
	close (OUT);
}
