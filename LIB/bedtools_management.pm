
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
    
package LIB::bedtools_management;
## bedtools_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of the bedtools for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(  run_BEDTOOLS_coverageBed run_BEDTOOLS_mergeBed run_BEDTOOLS_intersectBed
									run_BEDTOOLS_genomeCoverageBed run_BEDTOOLS_subtractBed
									run_BEDTOOLS_coverageBed_gen run_BEDTOOLS_sortBed);

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

use LIB::files_manipulation qw();

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db);


##################



=head2 run_BEDTOOLS_coverageBed

 Title   : run_BEDTOOLS_coverageBed
 Usage   : run_BEDTOOLS_coverageBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : the constructed output name
=cut
sub run_BEDTOOLS_coverageBed{
	my $cfg_hash = shift;
	my $bed_input = shift;
	my $outFile = shift;
	my $log_file = shift;
	my $target_bed = shift;
	
	my $param_str = " -a ".$target_bed;#$cfg_hash->{$group_id.'_target_bed'};	####TO BE SET!
	
	#Set the mode
	if (defined $cfg_hash->{'covbed_mode'}){
		$param_str .= " -".$cfg_hash->{'covbed_mode'}." ";
		#Set the output extension: .cov_hist for hist and cov_d for d
		$outFile .= ".".$cfg_hash->{'cov_'.$cfg_hash->{'covbed_mode'}.'_ext'};
	}else{
		log_and_exit("ERROR: Parameter covbed_mode in user_config must be set [hist|d]. Please check it..\n",$log_file);
	}

	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'coverbed_prog'}." $param_str -b $bed_input > $outFile";;
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
	
	return $outFile;
}

=head2 run_BEDTOOLS_mergeBed

 Title   : run_BEDTOOLS_mergeBed
 Usage   : run_BEDTOOLS_mergeBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_mergeBed{
	my $cfg_hash = shift;
	my $bed = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'mergebed_prog'}." -i $bed  $params > $outFile";
	print_and_log( " [ Executing command: $command ] ",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

=head2 run_BEDTOOLS_coverageBed_gen

 Title   : run_BEDTOOLS_coverageBed_gen
 Usage   : run_BEDTOOLS_coverageBed_gen(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_coverageBed_gen{
	my $cfg_hash = shift;
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = "";
	#if ( $bed1 =~ /\.bam$/ and $bed2 =~ /\.bed$/ ){
		##$params .= " -bed ";
		#$files_str = " -abam $bed1 -b $bed2 ";
	#}elsif ($bed2 =~ /\.bam$/ and $bed1 =~ /\.bed$/){
		#$files_str = " -abam $bed2 -b $bed1 ";
	#}else{
		$files_str = " -a $bed1 -b $bed2 ";
	#}
		
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'coverbed_prog'}." $files_str $params > $outFile";
	if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

=head2 run_BEDTOOLS_intersectBed

 Title   : run_BEDTOOLS_intersectBed
 Usage   : run_BEDTOOLS_intersectBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_intersectBed{
	my $cfg_hash = shift;
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = "";
	if ( $bed1 =~ /\.bam$/ and $bed2 =~ /\.bed$/ ){
		#$params .= " -bed ";
		$files_str = " -abam $bed1 -b $bed2 ";
	}elsif ($bed2 =~ /\.bam$/ and $bed1 =~ /\.bed$/){
		$files_str = " -abam $bed2 -b $bed1 ";
	}else{
		$files_str = " -a $bed1 -b $bed2 ";
	}
		
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'intersectbed_prog'}." $files_str $params > $outFile";
	if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

=head2 run_BEDTOOLS_subtractBed

 Title   : run_BEDTOOLS_subtractBed
 Usage   : run_BEDTOOLS_subtractBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools subtract
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_subtractBed{
	my $cfg_hash = shift;
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = " -a $bed1 -b $bed2 ";
		
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'subtractbed_prog'}." $files_str $params > $outFile";
	if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}


=head2 run_BEDTOOLS_sortBed

 Title   : run_BEDTOOLS_sortBed
 Usage   : run_BEDTOOLS_sortBed(   );

 Function: Runs sortBed taking in input the bam file
								
			INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_sortBed{
	my $cfg_hash = shift;
	my $file = shift;
	my $params = shift;
	my $outFile = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'sortbed_prog'}." -i $file $params  > $outFile";
	print_and_log( "[ Executing command: $command ]",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}


=head2 run_BEDTOOLS_genomeCoverageBed

 Title   : run_BEDTOOLS_genomeCoverageBed
 Usage   : run_BEDTOOLS_genomeCoverageBed(   );

 Function: Runs genomeCoverageBed taking in input the bam file
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_genomeCoverageBed{
	my $cfg_hash = shift;
	my $bam = shift;
	my $outFile = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = $cfg_hash->{'bedtools_path'}."/".$cfg_hash->{'genomecov_prog'}." -ibam $bam -bga  > $outFile";
	print_and_log( "[ Executing command: $command ]",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

1;
