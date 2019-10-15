
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
    
package LIB::annotation_management;
## annotation_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of tools for the annotation
#In particular we use: Annovar
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( run_ANNOVAR_convert2annovar run_ANNOVAR_download_db  
						run_ANNOVAR_annotate_variation run_ANNOVAR_table_annovar
						run_ANNOVAR_download_WG_fasta run_ANNOVAR_retrieveseq
						all_variants_2_VCF checkIncludedAnnovarDBs);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
																print_and_log log_and_exit  execute_threads
																build_input_name_from_executed);

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present extract_name);
use LIB::files_manipulation qw();

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
													select_distinct_samples get_id_if_exists_from_db fetch_all_rows);

##################



=head2 run_ANNOVAR_convert2annovar

 Title   : run_ANNOVAR_convert2annovar
 Usage   : run_ANNOVAR_convert2annovar(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs convert2annovar
					
					This tool converts a VCF file in a convenient shape for Annovar
					
 Returns : nothing
=cut
sub run_ANNOVAR_convert2annovar{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $input_file = shift;	
	my $outFile = shift;		
	
	my $param_str = "";
	if ( defined $cfg_hash->{'ann_include_all_var'} ){
		$param_str .= "  ".$cfg_hash->{'ann_include_all_var'}." ";
	}
	if ( $cfg_hash->{'ann_includeinfo'} eq 'YES'){
		$param_str .= "  -includeinfo ";
	}
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$cfg_hash->{'conv2annov_prog'}." $param_str $input_file  > $outFile";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 checkIncludedAnnovarDBs

 Title   : checkIncludedAnnovarDBs
 Usage   : checkIncludedAnnovarDBs( - configFile -> file with the user configuration
                              - variablesFile -> the path to a file with all variables written
                              - lineToCheck -> line in the variables file to be used
          );

 Function: this subroutine checks if an Annovar database is among those "included" in VarGenius. If this is true,
 uncompress and copies the database to the annovar directory.

 Returns : nothing

=cut
sub checkIncludedAnnovarDBs {
	my $cfg_hash = shift;
	my $db_name = shift;
	my $ext = shift;
	my $log_file = shift;
	
	my @included_dbs = split(",",$cfg_hash->{'annov_included'});
	#If this database is included into the program bundle, then uncompress and move into the annovar folder
	if ( grep {/\b$db_name\b/} @included_dbs){
		print_and_log("$db_name is included in the installation. Putting it into ".$cfg_hash->{'annovar_db_f'}."\n",$log_file);
		my_extract_any_file($cfg_hash->{'main_data_folder'}."/".$db_name.".".$ext.".".$cfg_hash->{'gz_ext'},$cfg_hash->{'annovar_db_f'});
	}
}


=head2 run_ANNOVAR_download_db

 Title   : run_ANNOVAR_download_db
 Usage   : run_ANNOVAR_download_db(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs run_ANNOVAR_annotate_variation
					
					This tool executes annotate_variation tool from Annovar and can 
					download a database
					
					/pico/home/userexternal/fmusacch/bin/annovar/annotate_variation.pl
					 -downdb  -webfrom annovar -buildver hg19 refGene  
					 /pico/work/TELET_UDP/home/shared/references/annovar/humandb

					
 Returns : nothing
=cut
sub run_ANNOVAR_download_db{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $db_2_dl = shift;
	my $ann_type = shift;
	my $source = shift;
	
	my $param_str = "";
	
	#Source selected (for region base annotation the annovar source want be used)
	if ( (defined $cfg_hash->{'ann_db_source'}) and  ($ann_type ne 'r') and ($source ne 'UCSC') ){
		$param_str .= " -webfrom ".$cfg_hash->{'ann_db_source'};
	}
	#Build version
	if ( defined $cfg_hash->{'ann_build_ver'} ){
		$param_str .= " -buildver ".$cfg_hash->{'ann_build_ver'};
	}else{ print_and_log( "ERROR: ann_build_ver is not defined...\n",$log_file);}
	
	#Database name
	if ( defined $db_2_dl ){
		$param_str .= " $db_2_dl ";
	}else{ print_and_log( "ERROR: db_2_dl is not defined...\n",$log_file);}

	my $db_folder = $cfg_hash->{'annovar_db_f'};
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$cfg_hash->{'annot_var_prog'}." -downdb $param_str $db_folder";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_ANNOVAR_retrieveseq

 Title   : run_ANNOVAR_retrieveseq
 Usage   : run_ANNOVAR_retrieveseq(   );

 Function: Downloads whole genome fasta files. It is needed to have GenCode ids
					
					/pico/home/userexternal/fmusacch/bin/annovar/retrieve_seq_from_fasta.pl 
					-format genericGene  -seqdir /pico/work/TELET_UDP/home/shared/references/annovar/humandb/hg19_seq
					-outfile /pico/work/TELET_UDP/home/shared/references/annovar/humandb/hg19_wgEncodeGencodeBasicV19Mrna.fa 
					/pico/work/TELET_UDP/home/shared/references/annovar/humandb/hg19_wgEncodeGencodeBasicV19.txt
 Returns : nothing
=cut
sub run_ANNOVAR_retrieveseq{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $db_2_dl = shift;
	
	my $param_str = "";
	
	my $db_folder = $cfg_hash->{'annovar_db_f'};
		
	#Build version
	if ( defined $cfg_hash->{'annovar_seqdir'} ){
		$param_str .= " -seqdir $db_folder/".$cfg_hash->{'ann_build_ver'}."_".$cfg_hash->{'annovar_seqdir'};
	}else{ print_and_log( "ERROR: annovar_seqdir is not defined...\n",$log_file);}

	#Build version
	my $build_ver = "";
	if ( defined $cfg_hash->{'ann_build_ver'} ){
		$build_ver = $cfg_hash->{'ann_build_ver'};
	}else{ print_and_log( "ERROR: ann_build_ver is not defined...\n",$log_file);}
	
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$cfg_hash->{'retrieveseq_prog'}." -format genericGene $param_str -outfile $db_folder/$build_ver"."_".$db_2_dl."Mrna.fa  $db_folder/$build_ver"."_".$db_2_dl.".txt"; 
								
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 run_ANNOVAR_download_WG_fasta

 Title   : run_ANNOVAR_download_WG_fasta
 Usage   : run_ANNOVAR_download_WG_fasta(   );

 Function: Downloads whole genome fasta files. It is needed to have GenCode ids
					
 Returns : nothing
=cut
sub run_ANNOVAR_download_WG_fasta{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $db_2_dl = shift;
	
	my $param_str = "";
	
	#Build version
	if ( defined $cfg_hash->{'ann_build_ver'} ){
		$param_str .= " -buildver ".$cfg_hash->{'ann_build_ver'};
	}else{ print_and_log( "ERROR: ann_build_ver is not defined...\n",$log_file);}
	
	#Database name
	if ( defined $db_2_dl ){
		$param_str .= " $db_2_dl ";
	}

	my $db_folder = $cfg_hash->{'annovar_db_f'}."/".$cfg_hash->{'ann_build_ver'}."_".$db_2_dl;
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$cfg_hash->{'annot_var_prog'}." -downdb $param_str $db_folder";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 run_ANNOVAR_annotate_variation

 Title   : run_ANNOVAR_annotate_variation
 Usage   : run_ANNOVAR_annotate_variation(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs run_ANNOVAR_annotate_variation
					
					This tool executes annotate_variation tool from Annovar and can 
					annotate a Annovar format file
					
 Returns : nothing
=cut
sub run_ANNOVAR_annotate_variation{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $in_file = shift;	
	my $out_file = shift;
	my $db_type = shift;
	my $ann_type = shift;
	
	my $param_str = "";
	
	#Annotation type
	if ( defined $ann_type ){
		$param_str .= " -$ann_type ";
	}
	#Database name
	if ( defined $db_type ){
		$param_str .= " -dbtype $db_type ";
	}
	#Build version
	if ( defined $cfg_hash->{'ann_build_ver'} ){
		$param_str .= " -buildver ".$cfg_hash->{'ann_build_ver'};
	}else{ print_and_log( "ERROR: ann_build_ver is not defined...\n",$log_file);}
		

	my $db_folder = $cfg_hash->{'annovar_db_f'};
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$cfg_hash->{'annot_var_prog'}."  $param_str -out $out_file $in_file $db_folder";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}



=head2 run_ANNOVAR_table_annovar

 Title   : run_ANNOVAR_table_annovar
 Usage   : run_ANNOVAR_table_annovar(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs run_ANNOVAR_annotate_variation
					
					This tool executes annotate_variation tool from Annovar and can 
					annotate a Annovar format file
					
					The other_info parameter permits to attach to the annotation the columns
					inserted in the input file.
					
 Returns : nothing
=cut
sub run_ANNOVAR_table_annovar{
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $log_file = shift;
	my $in_file = shift;	
	my $out_file_suff = shift;
	
	my $param_str = "";

	#CSV output
	if ( $cfg_hash->{'annov_csvout'} eq 'YES'){
		$param_str .= " -csvout ";
	}
	#CSV output
	if (  $cfg_hash->{'ann_otherinfo'} eq 'YES'){
		$param_str .= " -otherinfo ";
	}
	#Nastring
	if (  $cfg_hash->{'annov_nastring'} eq 'YES'){
		$param_str .= " -nastring . ";
	}
	#Remove
	if ( $cfg_hash->{'annov_remove'} eq 'YES'){
		$param_str .= " -remove ";
	}		
	#Build version
	if ( defined $cfg_hash->{'ann_build_ver'} ){
		$param_str .= " -buildver ".$cfg_hash->{'ann_build_ver'};
	}else{ print_and_log( "ERROR: ann_build_ver is not defined...\n",$log_file);}
	
	#annov_input
	if ( $cfg_hash->{'annov_vcf_input'} eq 'YES' ){
		$param_str .= " --vcf_input ";
	}		

	#######
	my $modes_string = " -operation ";
	my $dbs_string = " -protocol ";
	#Arguments to be used for each operation
	my $op_arguments = "";
	#annov_input
	if ( $cfg_hash->{'annov_hgvs'} eq 'YES'){
		$op_arguments = " -argument  '";
	}			
	#Get the annotation types wanted
	my @ann_types = split($cfg_hash->{'word_sep'},$cfg_hash->{'annov_ann_types'});
	#For each annotation type get its databases
	foreach my $ann_type (@ann_types){
		#Check if the list of annotations is present
		if ( defined $cfg_hash->{'annov_'.$ann_type.'_dbs'}) {
			#Get the list of needed databases
			my @databases = split($cfg_hash->{'word_sep'},$cfg_hash->{'annov_'.$ann_type.'_dbs'});

			#Get annotation with all the needed databases
			foreach my $db (@databases){
				#Obtain the database name if the database is already installed
				my $db_present = "";
				$db_present = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_annovar_dbs_table'},$cfg_hash->{'db_annovar_dbs_name'}
								,$cfg_hash->{'db_annovar_dbs_name'},"'".$db."'",$cfg_hash->{'db_annovar_dbs_type'},
								"'".$ann_type."'");
	
				#Length smaller than 2 means that the return is -1, db does not exist
				if ( (length $db_present) > 2 ){
					#Parameters for operation g only
					if ( $ann_type eq 'g' ){
						#splicing_threshold
						if ( defined $cfg_hash->{'annov_splicing_threshold'} ){
							$op_arguments .= " -splicing_threshold ".$cfg_hash->{'annov_splicing_threshold'}." ";
						}	
					}
					#Condition for 1000g database because they are more than 1
					if ( $db =~ /^1000g/ ){
						my @modes = split($cfg_hash->{'word_sep'},$cfg_hash->{'annov_1000g_popul'});
						foreach my $mode (@modes){
							$modes_string .= "$ann_type,";
							$dbs_string .= "$db\_$mode,";
							#annov_hgvs
							if ( $cfg_hash->{'annov_hgvs'} eq 'YES' ){
								$op_arguments .= "-hgvs,";
							}	
						}
					}#EXTERNAL DATABASES (dbtype generic)
					elsif( $db =~ /\_EXT$/ ){
						my $db_name = $db;
						$db_name =~ s/\_EXT//;
						$modes_string .= "$ann_type,";
						$dbs_string .= "$db_name,"; 
						#parameters for annovar_annotation.pl
						$op_arguments .= "-hgvs,";
						$param_str .=	" -genericdbfile ".$cfg_hash->{'annovar_db_f'}."/".
							$cfg_hash->{'ann_build_ver'}."_".$db_name.".".$cfg_hash->{'txt_ext'};
					}
					#BED files
					#Alternatively, if you want to annotate multiple columns int tab-delimited format, just run the same operation twice
					#-protocol bed,bed -operation r,r -bedfile file,file -arg '-colWanted 4','-colWanted 5'
					#The name of columns will be: bed,bed2,bed3... and so on..
					elsif( $db =~ /\_BED$/ ){
						my $db_name = $db;
						$db_name =~ s/\_BED//;
						$modes_string .= "$ann_type,";
						$dbs_string .= "bed,"; 
						#parameters for annovar_annotation.pl
						$op_arguments .= "-colsWanted 4 -hgvs,";
						$param_str .=	" -bed ".$cfg_hash->{'ann_build_ver'}."_".$db_name.".".$cfg_hash->{'bed_ext'};
					}
					#UCSC
					else{
						#If the database has been taken from UCSC we need to remove the suffix _UCSC
						if ($db =~ /\_UCSC$/){
							$db =~ s/\_UCSC//;
						}
						$modes_string .= "$ann_type,";
						$dbs_string .= "$db,"; 
						#annov_hgvs
						if ( $cfg_hash->{'annov_hgvs'} eq 'YES' ){
							$op_arguments .= "-hgvs,";
						}	
					}
				}else{
					print_and_log( "Database $db for Annovar does not exists... \n",$log_file);
				}											
			} 
		}
	}	
	chop($modes_string);
	chop($dbs_string);
	#annov_hgvs
	if ( $op_arguments ne '' ){
		chop($op_arguments);
		$op_arguments .= "'";
	}		
	#####

	my $db_folder = $cfg_hash->{'annovar_db_f'};
	#Execute the command
	my $command = $cfg_hash->{'annovar_path'}."/".$prog_used."  $param_str $modes_string $dbs_string $op_arguments -out $out_file_suff $in_file $db_folder";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}




=head2 all_variants_2_VCF

 Title   : all_variants_2_VCF
 Usage   : all_variants_2_VCF(   );

 Function: Fetches all the variants from the database and puts them in 
 a VCF file
					
 Returns : nothing
 
=cut
sub all_variants_2_VCF{
	my $cfg_hash = shift;
	my $vcf_to_ann = shift;
	my $log_file = shift;

	my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	my $id_str = $cfg_hash->{'db_vcf_id_str'};
	my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	my $info_str = $cfg_hash->{'db_vcf_info_str'};
	
	##############################################################
	#Generates a VCF file with all the variants from the database#
	##############################################################
	my $tab = "\t";
	#my $vcf_to_ann = "all_variants.vcf";
	#Query to select all the variants
	my $query = "SELECT * FROM  ".$cfg_hash->{'db_variants_table'}.";";							
	print  "Executing: $query\n";
	#fetch_all_rows
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'},$query);
	#Writes into a VCF
	open (ALL_VAR,">$vcf_to_ann") or die "ERROR: Cannot open $vcf_to_ann. The program will exit..\n";
	#print ALL_VAR $cfg_hash->{'db_var_chrom'}.$tab.$cfg_hash->{'db_var_pos'}.$tab.$cfg_hash->{'db_var_ident'}.$tab.
	#$cfg_hash->{'db_var_ref'}.$tab.$cfg_hash->{'db_var_alt'}.$tab.".\n";
	
	print ALL_VAR join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";
	
	#For each result print all the needed fields for a VCF
	foreach my $row (@$res) {
		my @arr = @$row;
		 #print join($tab, @$row), "\t.\n";
		 
		print ALL_VAR $arr[1].$tab.$arr[2].$tab.$arr[3].$tab.$arr[4].$tab.$arr[5].$tab.".".$tab.".".$tab.$arr[0]."\n";
	}
	close(ALL_VAR);


}


1;
