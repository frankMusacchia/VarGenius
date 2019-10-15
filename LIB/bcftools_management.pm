
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
    
package LIB::bcftools_management;
#Permits the management of the bcftools for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(  run_BCFTOOLS_call run_BCFTOOLS_filter);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

use LIB::programs_management qw(try_exec_command	print_and_log log_and_exit  );

use LIB::db_management qw(get_id_if_exists_from_db);

=head2 run_BCFTOOLS_call

 Title   : run_BCFTOOLS_call
 Usage   : run_BCFTOOLS_call(   );

 Function: Executes BCFTOOLS call for variant calling
								
					INPUT: needs a BCF file
					
 Returns : nothing
=cut
sub run_BCFTOOLS_call{
	my $cfg_hash = shift;
	my $bcf_input = shift;
	my $outFile = shift;
	my $param_str = shift;
	my $target_bed = shift;
	my $log_file = shift;
	

	if (file_not_present($target_bed) > 0 ){ print_and_log("WARNING no target bed file was usesìd using in freebayes!\n",$log_file);}
			
	if ( $target_bed ne '-'){ 
		$param_str .= " --targets-file $target_bed ";	
	}
		
	
	#Execute the command
	my $command = $cfg_hash->{'bcftools_path'}." ".$cfg_hash->{'bcftools_call_prog'}." -o $outFile $param_str $bcf_input ";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

=head2 run_BCFTOOLS_filter

 Title   : run_BCFTOOLS_filter
 Usage   : run_BCFTOOLS_filter(   );

 Function: Executes BCFTOOLS filter for variant calling
								filters is a comma separated numbers list
								
								The filters will be concatenated to $parameters, hence leave
								a white space when you do not have any input paramter
					INPUT: needs a BCF file
					
 Returns : nothing
=cut
sub run_BCFTOOLS_filter{
	my $cfg_hash = shift;
	my $bcf_input = shift;
	my $outFile = shift;
	my $param_str = shift;
	my $filters = shift;
	my $log_file = shift;

	#my @fields = split(":",$filters);
	#my @filter_set = split(",",$fields[1]);
	#foreach my $filter ( @filter_set){
		#if ( defined $cfg_hash->{$fields[0].'_filterName'.$filter} ){
			#$parameters .= " -s".$cfg_hash->{$fields[0].'_filterName'.$filter}." ";
		#}	
		#if ( defined $cfg_hash->{$fields[0].'_filterExpression'.$filter} ){
			#$parameters .= " ".$cfg_hash->{$fields[0].'_filterExpression'.$filter}." ";
		#}	
	#}	

	my @filter_set = split(",",$filters);
	foreach my $filter ( @filter_set){
		if ( defined $cfg_hash->{'filterName'.$filter} ){
			$param_str .= " -s".$cfg_hash->{'filterName'.$filter}." ";
		}	
		if ( defined $cfg_hash->{'filterExpression'.$filter} ){
			$param_str .= " ".$cfg_hash->{'filterExpression'.$filter}." ";
		}	
	}	
	
	#Execute the command
	my $command = $cfg_hash->{'bcftools_path'}." ".$cfg_hash->{'bcftools_filter_prog'}." -o $outFile $param_str $bcf_input ";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

1;
