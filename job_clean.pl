#!/usr/bin/perl
#To work at CINECA
##!/opt/software/perl/bin/perl
#To work at TIGEM

#Works at CINECA
use lib '/pico/work/TELET_TIGEM/VarGeniusBeta'; #PICO CINECA
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5/x86_64-linux-thread-multi';
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5';
#Parallel::ForkManager
#use lib qw(/pico/home/userexternal/fmusacch/bin/perl_lib/Parallel-ForkManager-1.17/lib/);
use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

#TIGEM
#use lib '/home/users/musacchia/VarGenius';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED


# Package VarGenius - Variant Discovery and Annotation Tool
# Author: Francesco Musacchia (2019)
#
#This software is a pipeline to discover SNPs and Variants through clinical
#sequencing data and to update our internal databases

#Libraries
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(try_exec_command print_and_log log_and_exit kill_job
																JOB_get_info_from_jobname);
															
#Using a library to manage files
use LIB::files_management qw(   delete_file save_hash load_hash);

#Environmental variables
#Get the file with the saved hash of dependencies
my $deps_f = $ENV{DEP_HASH_F};
#Get the separator to use for data
#my $sep = $ENV{SEP};
#Get the task name to run
my $task = $ENV{TASK};
#Gets the log file
my $log_file = $ENV{LOG_FILE};
#Gets the log file
my $sleep_time = $ENV{JC_TIME};
#Gets the username to access the jobs
my $qsub_username = $ENV{QS_USER};
#Get the file with the parameters for the module
my $config_file = $ENV{CONFIG_FILE};


my $sep = ":";
my $status_field = 9;
my $id_field = 0;


#my $qstat_out = $ARGV[0];
#my $sleep_time = 3;
#Variables for computing execution times
my $globalStart = time;
my $partTime  = time;
my $partDuration = undef;

#Get the config hash
my $cfg_hash = load_hash($config_file);

#Continue checking until there are jobs in the hash of 
#dependencies
if ( -e $deps_f ){
	#Load the hash of dependencies
	my $deps = load_hash($deps_f);
	print Dumper\$deps;
	while ( scalar(keys %{$deps}) > 0 ){
		print "(".scalar(localtime)."): RUNNING $task..\n";

		#Printing qstat result in a file
		my $qstat_out = "qstat_out.txt";
		my $command = "qstat -n1 -t -u".$qsub_username." -w > $qstat_out";
		try_exec_command($command) or die "Unable to execute command: $command\n";

		my $start =0 ;
		my $stat_hash;
		#print "Looping...\n";
		open (QSTAT_OUT,"<$qstat_out") or die "ERROR: Cannot open $qstat_out. The program will exit..\n";
		while (my $line = <QSTAT_OUT>){
				if( $start == 1 ){
					#print $line."\n";
					my @fields = split(/\s+/,$line);
					my @name_parts = split(/\./,$fields[$id_field]);
					#print " id1: ".$fields[0]."\n";
					$stat_hash->{$name_parts[0]} = $fields[$status_field];
					#print "1: ".$fields[0]." - 2: ".$fields[1]." - 3:".$fields[2]." " ;#DEBUGCODE
				}
				if ($line =~ /---------------/){
					$start = 1;
				}
		}
		close(QSTAT_OUT);
		print Dumper\$stat_hash;


		# If some job in the hold status is waiting for nothing, remove it and all 
		# jobs depending by it
		foreach my $job_name ( keys %{$deps}){
			#print "Searching for $job_name..\n";
			#If job has finished, it can be removed from the list of deps
			if ( ! (defined $stat_hash->{$deps->{$job_name}->{'id'}} ) ){
				#print_and_log ( "Deleting: ".$deps->{$job_name}."\n",$log_file);				
				delete($deps->{$job_name});	
				#Saving the hash with the ne dependencies
				save_hash(\$deps,$deps_f);				
			}
			#If the job is still in hold status then check all dependencies. 
			#If they have finished (or gone out) the subject job will be killed
			elsif ( $stat_hash->{$deps->{$job_name}->{'id'}} eq 'H'){
				print_and_log( "$job_name is in holding status...\n",$log_file);
				print_and_log( "Its dependencies are: ".$deps->{$job_name}->{'deps'}."\n",$log_file);
				my @dep_jobs = split($sep,$deps->{$job_name}->{'deps'});
				my $num_dep_jobs = scalar(@dep_jobs);
				my $exited = 0;
				
				foreach my $dep_job ( @dep_jobs ){
					if ( ! (defined $stat_hash->{$dep_job}) ){
						$exited++;
					}
				}
				print_and_log( "num_dep_jobs: $num_dep_jobs and exited: $exited\n",$log_file);
				if ( $exited == $num_dep_jobs){
						kill_job($deps->{$job_name}->{'id'});
						#Get the information associated with the jobname to remove			
						my ($task,$type,$name) = JOB_get_info_from_jobname($cfg_hash,$job_name);
						print_and_log ( "ERROR: Killed job $job_name because all depending jobs (".$deps->{$job_name}->{'deps'}.") exited.\n ",$log_file);
						print_and_log ( "Task: $task - Type: $type - Name: $name\n ",$log_file);
						
						#Delete also from hashes
						#print "Deleting ".$deps->{$job_name}->{'id'}."\n ";
						delete($stat_hash->{$deps->{$job_name}->{'id'}});	
						delete($deps->{$job_name});
						#Saving the hash with the ne dependencies
						save_hash(\$deps,$deps_f);
						
				}#else{print_and_log(  "Cannot kill job $job_name because it is waiting some working job.\n ",$log_file);}
			}

		}
		sleep($sleep_time);

		#Reload  the hash with the dependencies because some job may have
		#altered some other
		$deps = load_hash($deps_f);
	}
}



###################################################
#Before	to leave, remove useless files
###################################################

#Before exiting clean all the useless files if the user wants
if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
	
	#Get all the group ids 
	my @analysis_ids = split(",",$cfg_hash->{'all_analysis_ids'});
	
	#Get all the tasks that have been executed for this run
	my @tasks_to_exec = split(",",$cfg_hash->{'tasks_to_exec'});
	
	foreach my $analysis_id ( @analysis_ids ){
		print_and_log("Removing temporary files for analysis id: $analysis_id.\n",$log_file);
		print_and_log("Tasks executed:".$cfg_hash->{'tasks_to_exec'}."\n",$log_file);
		#Remove useless files after alignment
		my $task = $cfg_hash->{'align_step'};
		if( grep {/\b$task\b/} @tasks_to_exec ){
			
			#Get the folder where all the files have been inserted
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			print_and_log("Removing all SAM files for task $task from $outFolder ".
				      $outFolder."/*.".$cfg_hash->{'sam_ext'}.
				      "...\n",$log_file);
			#Remove all SAM files generated by BWA
			delete_file($outFolder."/*.".$cfg_hash->{'sam_ext'});
			#Remove all bam files not sorted and not indexed
			print_and_log("Removing all bam files not sorted and not indexed for task $task from $outFolder:".
				      $outFolder."/*_".$cfg_hash->{'align_step'}."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'}.
				      "...\n",$log_file);
			delete_file($outFolder."/*_".$cfg_hash->{'align_step'}."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'});
			#Remove all fastq useless files
			print_and_log("Removing all fastq useless files for task $task from $outFolder:".
				      $outFolder."/*_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'}.
				      "...\n",$log_file);
			delete_file($outFolder."/*_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'});
		}
		
		#Remove useless files after refine
		$task = $cfg_hash->{'refine_task'};
		if( grep {/\b$task\b/} @tasks_to_exec ){
			#print_and_log("Files for task $task.\n",$log_file);
			#Get the folder where all the files have been inserted
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			#Remove the target intervals from RealignerTargetCreator
			print_and_log("Removing files for task $task: ".$outFolder."/*.".$cfg_hash->{'targ_int_ext'}."....\n",$log_file);
			try_remove($outFolder."/*.".$cfg_hash->{'targ_int_ext'});
			
			#If the user ran also the IndelRealigner then the bam files produced are now useless
			#Remove also the .grp recalibration tables output of BaseRecalibrator
			print_and_log("Removing files for task $task: ".$outFolder."/*_".$cfg_hash->{'indel_real_step'}.".".$cfg_hash->{'bam_ext'}."....\n",$log_file);
			try_remove($outFolder."/*_".$cfg_hash->{'indel_real_step'}.".".$cfg_hash->{'bam_ext'});
			#Remove the recalibration table ouput of BaseRecalibrator
			print_and_log("Removing files for task $task: ".$outFolder."/*_".$cfg_hash->{'base_recal_step'}.".".$cfg_hash->{'recal_ext'}."....\n",$log_file);
			try_remove($outFolder."/*_".$cfg_hash->{'base_recal_step'}.".".$cfg_hash->{'recal_ext'});		
		}

	}
}




################################################Support subroutines
sub try_remove {
		my $path_to_rem = shift;
		
		my $del_ret = delete_file($path_to_rem);
		
		if ( $del_ret == 1){
			print "$path_to_rem succesfully deleted!\n";
		}elsif ( $del_ret == -1){
			print "ERROR: $path_to_rem was not deleted!\n";
			}elsif ($del_ret == -2){
				print "ERROR: cannot find $path_to_rem..\n";
				}
}



