#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2018>  <Francesco Musacchia>

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
    
package LIB::fastq_quality_management;
## files_management.pm
#Author: Francesco Musacchia
#Permits the management of quality check and trimming programs
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(run_trimmomatic  run_fastqc run_trimgalore);
}
use strict;
use warnings;
use File::Copy;#To manage files

use LIB::files_management qw( load_hash delete_file file_not_present extract_name);

use LIB::db_management qw(getSampleConfiguration_locked);

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(try_exec_command print_and_log log_and_exit);
																											
=head2 run_trimmomatic

 Title   : run_trimmomatic
 Usage   : run_trimmomatic(   );

 Function:  Trimmomatic is ran reading the target file and for each replicate the forward and reverse reads
  are used in the same command. Both for Single Pair and Paired Ends can be run but parameters always have to be chosen
  by the user in the parameters file. Actually there is the problem of how to give in input the parameters before the file inputed
 
 Returns : nothing

=cut
sub run_trimmomatic{
	my $config_file = shift;
	my $task  = shift;	
	my $sample_id = shift;
	my $prog_used = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $trimmed1 = shift; 
	my $trimmed2  = shift;
	my $log_file = shift;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $params;
	
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash->{'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'}, $cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}, $sample_id);
			
	my $program_path = $cfg_hash->{$prog_used.'_path'};
	my $threads = $cfg_hash->{$task.'_threads'};
	##Sample sheet variables (These I will take from the DB)
	#my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
	#my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
	
	##Output folder
	#my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
	##The ouput names that is giving VarGenius	
	#my $trimmed1 = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
	#my $trimmed2 = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};

	my $align_type = $cfg_hash->{'align_type'};
		
	my $fu = $trimmed1."_1_unpaired.fq.gz";
	my $ru = $trimmed2."_2_unpaired.fq.gz";

	#Java call
	my $java_call = $cfg_hash->{'java_path'};
	if (defined $cfg_hash->{'javamem_TRM'}){
		$java_call .= " -".$cfg_hash->{'javamem_TRM'}." ";
	}
	$java_call .= " -jar $program_path ";
		
	my $parameters = "";
	
	if (defined $cfg_hash->{'TRM_ILLUMINACLIP'}){
		$parameters .= " ILLUMINACLIP:".$cfg_hash->{'TR_adapters'}.$cfg_hash->{'TRM_ILLUMINACLIP'};
	}
	if (defined $cfg_hash->{'TRM_SLIDINGWINDOW'}){
		$parameters .= " SLIDINGWINDOW:".$cfg_hash->{'TRM_SLIDINGWINDOW'};
	}
	if (defined $cfg_hash->{'TRM_MINLEN'}){
		$parameters .= " MINLEN:".$cfg_hash->{'TRM_MINLEN'};
	}	
	my $command = "$java_call PE -phred33 -threads $threads $fq1 $fq2 $trimmed1 $fu $trimmed2 $ru $parameters";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#Delete unpaired trimmed sequences
	delete_file($fu);
	delete_file($ru);
	
	#SINGLE END.. WHEN WILL BE USED		
  #$command = "java -jar $program_path SE -phred33 -threads 20 $readsF1 $outSE $params";
} 


=head2 run_fastqc

 Title   : run_fastqc
 Usage   : run_fastqc(   );

 Function: Runs FastQC
 Returns : nothing
=cut
sub run_fastqc{
	my $config_file = shift;
	my $task  = shift;
	my $group_id = shift;
	my $sample_id = shift;
	my $prog_used = shift;
	my $fq1 = shift;
	my $fq2 = shift;	
	my $log_file = shift;
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $params;
	
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$sample_id);

	my $program_path = $cfg_hash->{$prog_used.'_path'};

	#print "The outfolder will be formed by $group_id and $task..\n";
	my $outFolder = $cfg_hash->{$group_id.'_'.$task.'_out_f'};
	my $align_type = $cfg_hash->{'align_type'};

	#Set Parameters
	my $parameters = "";
	if ( defined $cfg_hash->{$task.'_threads'}){
		$parameters .= " -t ".$cfg_hash->{$task.'_threads'}." ";
	}
	
	print "The log file is $log_file\n";
	#print Dumper\$params;
	#print_and_log( "Running $prog_used for job: $sample_id...\n",$log_file);#DEBUGCODE
	if ( defined $align_type and $align_type eq 'PE'){
			#print_and_log( "This alignment is a paired-end one...\n",$log_file);#DEBUGCODE
			#Aligning read 1 to reference
			my $command = "$program_path $parameters -o $outFolder $fq1 $fq2 ";
			print_and_log( "Executing command: $command\n",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";

	}else{
			log_and_exit( "ERROR: Alignment type undefined. Cannot execute $prog_used..\n",$log_file);
	}

}

=head2 run_trimgalore

 Title   : run_trimgalore
 Usage   : run_trimgalore(   );

 Function: Runs TrimGalore
 Returns : nothing
=cut
sub run_trimgalore{
	my $config_file = shift;
	my $task  = shift;
	my $analysis_id = shift;
	my $sample_id = shift;
	my $prog_used = shift;
	my $log_file = shift;
	
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $params;
	
	#Get infos for subject sample
	#print_and_log("Getting information for $sample_id from database ".$cfg_hash->{'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'}, $cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'}, $sample_id);
			
	my $program_path = $cfg_hash->{$prog_used.'_path'};
	my $threads = $cfg_hash->{$task.'_threads'};
	#Sample sheet variables (These I will take from the DB)
	my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
	my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};

	#print_and_log(" fq1: $fq1  - fq2: $fq2...\n",$log_file);#DEBUGCODE
	
	#Output file names
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	#print_and_log("outFolder $outFolder \n",$log_file);#DEBUGCODE
	
	#the names that will give TrimGalore
	my $trimmed1 = $outFolder."/".extract_name($params->{'fq1'},1);
	#print_and_log(" trimmed1: $trimmed1...\n",$log_file);#DEBUGCODE
	$trimmed1 .= $cfg_hash->{'trim_out_suff'}."1".$cfg_hash->{'trim_out_ext'};
	
	my $trimmed2 = $outFolder."/".extract_name( $params->{'fq2'},1);
	$trimmed2 .= $cfg_hash->{'trim_out_suff'}."2".$cfg_hash->{'trim_out_ext'};
	#print_and_log(" trimmed2: $trimmed2  ...\n",$log_file);#DEBUGCODE

	
	#The new name that is giving VarGenius	
	my $trimmed1_new = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
	my $trimmed2_new = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};

	my $align_type = $cfg_hash->{'align_type'};
	
	#print Dumper\$params;
	#print_and_log( "Running $prog_used for sample: $sample_id...\n",$log_file);#DEBUGCODE
	if ( defined $align_type and $align_type eq 'PE'){
		#print_and_log( "This alignment is a paired-end one...\n",$log_file);#DEBUGCODE
			my $params2 = " ";
			if (defined $cfg_hash->{'trimg_min_qual'}){
				$params2 .= " -q ".$cfg_hash->{'trimg_min_qual'}." ";
			}
			if (defined $cfg_hash->{'trimg_qual_enc'}){
				$params2 .= " --".$cfg_hash->{'trimg_qual_enc'}." ";
			}
			if (defined $cfg_hash->{'clip_R1'}){
				$params2 .= " --clip_R1 ".$cfg_hash->{'clip_R1'}." ";
			}
			if (defined $cfg_hash->{'clip_R2'}){
				$params2 .= " --clip_R2 ".$cfg_hash->{'clip_R2'}." ";
			}
			if ( $cfg_hash->{'qc_after_trim'} eq 'YES' ){
				$params2 .= " --fastqc --fastqc_args '--outdir $outFolder' ";
			}
			if ( defined  $cfg_hash->{'cutadapt_path'} ){
				$params2 .= " --path_to_cutadapt ".$cfg_hash->{'cutadapt_path'}." ";
			}
			$params2 .= " -o $outFolder --paired $fq1 $fq2 ";

			#Trimming the reads 
			my $command = "$program_path $params2 ";
			print_and_log( "Executing command: $command ",$log_file);
			try_exec_command($command) or die "Unable to execute command: $command\n";

	}else{
			log_and_exit( "ERROR: Alignment type undefined. Cannot execute $prog_used..\n",$log_file);
	}

	#Change the names of the output files
	#Check if they exist
	if (file_not_present($trimmed1) > 0 ){ die "ERROR $trimmed1 has not been created.\n";}
	if (file_not_present($trimmed2) > 0 ){ die "ERROR $trimmed2 has not been created.\n";}
	move($trimmed1,$trimmed1_new);
	move($trimmed2,$trimmed2_new);
}

1;
