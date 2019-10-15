
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
    
package LIB::vcftools_management;
#Permits the management of the vcftools for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_VCFTOOLS_sort);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

use LIB::programs_management qw(try_exec_command	print_and_log log_and_exit  );

use LIB::db_management qw(get_id_if_exists_from_db);

=head2 run_VCFTOOLS_sort

 Title   : run_VCFTOOLS_sort
 Usage   : run_VCFTOOLS_sort(   );

 Function: Executes VCFTOOLS sort to sort a VCF file
								
					INPUT: needs a VCF file
					
 Returns : nothing
=cut
sub run_VCFTOOLS_sort{
	my $cfg_hash = shift;
	my $vcf_input = shift;
	my $outFile = shift;
	my $param_str = shift;
	my $log_file = shift;
	

	#Execute the command
	my $command = $cfg_hash->{'vcftools_path'}."/".$cfg_hash->{'vcftools_sort_prog'}." $param_str $vcf_input > $outFile ";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

1;
