
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
    
package LIB::output_management;
## output_management.pm
#Author: Francesco Musacchia  (2018)
#Permits the management of the output from VarGenius
#
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( import_vcf_2_db import_ann_out_2_db get_non_covered_regions
			 export_results get_non_cov_reg_with_bedtools
			 splicing_var_dist_2_final_out rearrange_rawtabular_out build_output_from_annotation
			 generate_vcf_to_annotate output_from_annotation_without_db get_gene_annotation
			 print_variants_stats generate_igvsession_xml get_allgene_info separate_joint_outputs
			 send_email get_gene_coverage_above get_all_samples_coverage_table
			 get_all_samples_reads_number get_statistics_tables get_all_genes_annotation
			 separate_analysis_from_joint tag_panel_genes import_annovarout_2_db
			 import_vgout_2_db import_vcf_2_db_fast import_annovarout_2_db_fast
			 print_alignment_statistics_table evaluate_quality_of_sequencing
			 get_gene_4_interval get_non_cov_exons run_BEDTOOLS_coverageBed_gen);
}

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 
use Vcf;#Use vcf_tools library
use File::Copy;#To manage files

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw( print_and_log log_and_exit try_exec_command
																sort_samples build_input_name_from_executed
																separate_input_ids);

#Using a library to manage files
use LIB::files_management qw(save_hash load_hash file_not_present 
														delete_file  extract_col_from_file
														extract_colnum_from_file_linux extract_columns_from_file
														delete_columns delete_rows_containing get_col_index
														extract_name twocolfile_to_hash file_list_to_array
														append_str_2_file_if_path_notexist
														list_to_array);


#Using a library for database management
use LIB::db_management qw(insert_if_not_exists_locked insert_into_table_locked
							get_id_if_exists_from_db do_query_select_all do_query
													do_fetch_row_array select_element_from_db insert_only_if_not_exists_locked
													fetch_all_rows get_dbelement_ifexists_woconn insert_only_if_not_exists_locked_woconn
													insert_into_table_locked_woconn do_fetch_row_array_woconn
													getSampleConfiguration_locked insert_into_table_woconn
													delete_from_db get_id_if_exists_from_db_woconn load_csv_2_db
													update_table update_table_2) ;

#Using a library for bedtools functions management
use LIB::bedtools_management qw( 	run_BEDTOOLS_genomeCoverageBed run_BEDTOOLS_mergeBed 
									run_BEDTOOLS_intersectBed run_BEDTOOLS_coverageBed_gen );

#Using a library for standard utilities								
use LIB::std_lib qw( correct_type);






=head2 import_vcf_2_db_fast

 Title   : import_vcf_2_db_fast
 Usage   : import_vcf_2_db_fast( -database => 'name of the database,
                               );

 Function: Imports the VCF file content to the given table of variants in the database
			using the COPY command. Hence it first writes a file with all the rows
			and then loads it with a query.
			
					#Split different alternative alleles using ALT field
					#Put the variants information in the variants table if the comp_id is not present
					#Generates a VCF for variants that have to be annotated
					#Then inserts information from INFO field Separating using the comma, if there are diff alleles
					#Finally generates the Table for genotypes for different samples
 Returns : 

=cut			
sub import_vcf_2_db_fast {
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $vcf_file = shift;
	my $vcf_to_ann = shift;
	my $vcf_all_var = shift;
	
  #print_and_log( "Opening VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  #Loads the VCF
	my $vcf = Vcf->new(file=>$vcf_file);
	#Reads the header
	$vcf->parse_header();
	
	#Field names and indexes
	my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	my $chr_ind = $cfg_hash->{'vcf_chr_ind'};
	my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	my $pos_ind = $cfg_hash->{'vcf_pos_ind'};
	my $id_str = $cfg_hash->{'db_vcf_id_str'};
	my $id_ind = $cfg_hash->{'vcf_id_ind'};
	my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	my $ref_ind = $cfg_hash->{'vcf_ref_ind'};
	my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	my $alt_ind = $cfg_hash->{'vcf_alt_ind'};
	my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	my $qual_ind = $cfg_hash->{'vcf_qual_ind'};
	my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	my $filter_ind = $cfg_hash->{'vcf_filter_ind'};
	my $info_str = $cfg_hash->{'db_vcf_info_str'};
	my $info_ind = $cfg_hash->{'vcf_info_ind'};
	my $genot_ind = $cfg_hash->{'vcf_genot_ind'};
	my $prob_gen_ind = $cfg_hash->{'vcf_prob_gen_ind'};
	my $comp_id_str = $cfg_hash->{'db_vcf_comp_id_str'};
	#Separators
	my $sep = ",";
	my $tab = "\t";
	my $vcf_info_sep = $cfg_hash->{'vcf_info_sep'};
	my $vcf_gen_sep = $cfg_hash->{'vcf_gen_sep'};
	
	my $var_inserted = 0;
	#VCF file of variants to be annotated
	open (VCF_ANN,">".$vcf_to_ann) or die "ERROR: Cannot open $vcf_to_ann. The program will exit..\n";
	print VCF_ANN  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";

	#VCF file with all variants
	open (VCF_ALL,">".$vcf_all_var) or die "ERROR: Cannot open $vcf_all_var. The program will exit..\n";
	print VCF_ALL  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";

	##Open also a file for the STATISTICS table to load into the db
	#my $new_statistics_table = $vcf_all_var.".statistics.tab";
	#open (VCF_STAT,">".$new_statistics_table) or die "ERROR: Cannot open $new_statistics_table. The program will exit..\n";

	##Open also a file for the GENOTYPE_SAMPLES table to load into the db
	#my $new_genotsam_table = $vcf_all_var.".statistics.tab";
	#open (VCF_GEN,">".$new_genotsam_table) or die "ERROR: Cannot open $new_genotsam_table. The program will exit..\n";
		
	#Make a query to get all the variants compid so that you can match them from an hash 
	my $compids_hash;
	#to verify if it must be inserted
	my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_var_compid'}." FROM  ".$cfg_hash->{'db_variants_table'}.";";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'},$query);
	foreach my $row (@$res) {
		#Get the row values into an array
		my @arr = @$row;
		my $fld_num = 0;
		my $compid = "";
		my $varid = -1;
		#Each element of the array is now a field of the table
		foreach  my $elem (@arr) {
			if (defined $elem and $fld_num == 0) {
				$varid = $elem;
			}
			if (defined $elem and $fld_num == 1) {
				$compid = $elem;
			}	
			$fld_num++;		
		}
		#Put the elements into an hash
		$compids_hash->{$compid} = $varid;
		#print_and_log( "Putting  $varid into hash for: $compid\n",$log_file);#DEBUGCODE
	}
	
	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
		
	#my $group_id =1;#GET THE GROUP ID....##########
	print_and_log( "Reading VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  # This will split the fields and print a list of CHR:POS
  while (my $x = $vcf->next_data_array()){
		
	 
#	if ($var_inserted < 100) {#DEBUGCODE
		   #print_and_log( "$$x[$chr_ind]:$$x[$pos_ind]:$$x[$id_ind]:$$x[$ref_ind]:$$x[$alt_ind]\n",$log_file);#DEBUGCODE

		#An array to collect the ids corresponding to each alternate allele variant
		my @alt_ids = ();   
    ###################################################
    ####Adding variants and alternative variants ######
    ###################################################  
		#Construct the string with first information from the VCF
		my $fields = $chr_str.$sep.$pos_str.$sep.$id_str.$sep.$ref_str.$sep.$alt_str.$sep.$comp_id_str;
		#Construct the string with values
		my $start_values = "'".$$x[$chr_ind]."'".$sep."'".$$x[$pos_ind]."'".$sep."'".$$x[$id_ind]."'".$sep."'".$$x[$ref_ind]."'";
	      
		     
		#Split different alternative alleles using ALT field
		my @alts = split('\,', $$x[$alt_ind]);
		my $alt_id = 0;
		foreach my $alt (@alts){
	    #Construct the comp_id
	    my $comp_id = $$x[$chr_ind]."_".$$x[$pos_ind]."_".$$x[$ref_ind]."_".$alt;
	    my $values = $start_values.$sep."'".$alt."'".$sep."'".$comp_id."'";
	    #print_and_log( "Loading fields: $fields and  values: $values...\n",$log_file);#DEBUGCODE
	    
	    #Put the variants information in the variants table if the comp_id is not present
	    my $alt_id = -1;
	    my $existing_id = -1;
	    if ( !(defined $compids_hash->{$comp_id}) ){
				#print_and_log( "Putting  a new variant: $fields values:$values\n",$log_file);#DEBUGCODE
				$alt_id = insert_into_table_woconn($dbh,$cfg_hash->{'db_variants_table'},$fields,$values); 
				#print_and_log( "Id: $alt_id\n",$log_file);#DEBUGCODE
			}else{
				 $existing_id = $compids_hash->{$comp_id};
			}
		#OLD SLOWER WAY 
	    #my ($alt_id,$existing_id) = insert_only_if_not_exists_locked_woconn($dbh,$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'}
		#				    ,$fields,$values,$cfg_hash->{'db_vcf_comp_id_str'},"'".$comp_id."'"); 	
	   
	    
	    #The variant was not there it has ben added and altid is the id
	    if ( $alt_id >= 0){
				#Print all variants
				print VCF_ALL join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$alt_id)."\n";
		    #Print the variant in the VCF if it must be annotated
		    #We put the var identifier in the INFO field that can be useful to insert annotation results
		    #my $var_id_fld = $cfg_hash->{'db_var_id'}."=".$alt_id;
		    print VCF_ANN join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$alt_id)."\n";
		    #Save variant ids in order
		    push(@alt_ids,$alt_id);
		    #print_and_log( "Variant was not there id:$alt_id...\n",$log_file);#DEBUGCODE
	    }elsif ($existing_id >= 0){
				#Print all variants
	      print VCF_ALL join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$existing_id)."\n";					  
		    #The variant was already there. Put the existing id
		    #print_and_log( "Variant was already there id:$existing_id...\n",$log_file);#DEBUGCODE
		    push(@alt_ids,$existing_id);
	    }else{
		    log_and_exit("ERROR: during the import of $comp_id. Cannot understand what is happen with the db...\n",$log_file);
	    }
	    #$alt_id++;#USED FOR TRIAL
		}
			    
			    
    ###################################################
    ####Inserting information from INFO field##########
    ###################################################
    
    #hash with alleles information
		my $alt_hash;
		#Split the INFO field
		my @info_fields = split($vcf_info_sep,$$x[$info_ind]);
    #For each field of the INFO column...
    foreach my $info_field (@info_fields){
		    my @info_f_parts = split("=",$info_field);
		    #Separate using the comma, if there are diff alleles
		    my @alt_infos = split(",",$info_f_parts[1]);
		    #If there is the alternative info, get it
		    for (my $alt_info_ind = 0; $alt_info_ind < scalar(@alt_infos); $alt_info_ind++){
			    #Store the  alternative info for each of them
			    #alt_hash->{1}->{DP} = 6
			    $alt_hash->{$alt_info_ind}->{$info_f_parts[0]} = $alt_infos[$alt_info_ind];
		    }
		}	
		
		#Taking filter and quality from the main line
		my $info_fields_base = $cfg_hash->{'db_qual'}.",".$cfg_hash->{'db_filter'}.",".$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'}.",";#.join(",",@info_fields);
		my $info_values_base = $$x[$qual_ind].",'".$$x[$filter_ind]."',$group_id";
		
		#Go through the hash with information for each allele
		foreach my $alt_info_ind ( keys %{$alt_hash}){
	    my $info_fields = $info_fields_base;
	    my $info_values = $info_values_base.",".$alt_ids[$alt_info_ind].",";
	    #Get all the INFO values for the allele
	    #The values are taken as numbers here
	    while (my ($key, $value) = each %{ $alt_hash->{$alt_info_ind} } ) {
			    #print "$key = $value \n";
			    $info_fields .= "$key,";
			    #if the value is not  a number enclose in quotes
			    if ( correct_type($value,"real") ){	
						$info_values .= "$value,";
					}else{$info_values .= "'".$value."',";}
			    
	    }
	    #Remove the last comma
	    chop($info_fields);
	    chop($info_values);

	    #print_and_log( "Loading fields: $info_fields and  values: $info_values...\n",$log_file);#DEBUGCODE
		my $alt_id = insert_into_table_woconn($dbh,$cfg_hash->{'db_var_statistics_table'},$info_fields,$info_values);		
		}
		
		###############################################################
		###########Table for genotypes for different samples###########
		###########from FORMAT field##################################
		###############################################################
		
		my $format_fields_base = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_var_ids'}.",";
		my $format_values_base = "";
		
		#Get the list of samples using the Vcf library
		my (@samples) = $vcf->get_samples();
		#Go through the genotype information per sample
    foreach my $sample (@samples){
	    my $format_fields_str = $format_fields_base;
	    
	    my $sample_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
						$cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_analysis_id'},"'".$sample."'".",".$group_id);
	    #Joins the variant ids associated to different alleles if present
	    my $var_ids = join(",",@alt_ids);
	    my $format_values_str = $group_id.",".$sample_id.",'".$var_ids."',";
	    #Get the array with genotyping names
	    my @format_fields = split($vcf_gen_sep,$$x[$genot_ind]);
	    #Get the values corresponding to this sample
	    my $sample_col = $vcf->get_column($x, $sample);
	    #For each format field
	    foreach my $format_field (@format_fields){
		    #Get the index of the genotype information
		    my $form_idx = $vcf->get_tag_index($$x[$genot_ind],$format_field,':');
		    #If it exists
		    if ($form_idx != -1) {
				    #Add the field name
				    $format_fields_str .= $format_field.",";
				    #Get the value using the index
				    my $form_val = $vcf->get_field($sample_col,$form_idx);
				    #Concatenate the value
				    $format_values_str .= "'".$form_val."',";
		    }
		    ##If it is the GT information than use that to associate the allele variant
		    #if ( $format_field eq 'GT'){
				    ##Get the value using the index
				    #my $form_val = $vcf->get_field($sample_col,$form_idx);
				    #my @gen_flds = split(/[\|\/]/,$form_val);
		    #}
	    }
	    #Remove the last comma
	    chop($format_fields_str);
	    chop($format_values_str);
	    #print_and_log( "Loading fields: $format_fields_str and  values: $format_values_str...\n",$log_file);
	    #fill the table with genotypes for each sample
	    my $var_sample_id = insert_into_table_woconn($dbh,$cfg_hash->{'db_genotype_sample_table'},$format_fields_str,$format_values_str);
	 }
#	}#DEBUGCODE 100vars
	
	$var_inserted++;
  }

 	$dbh->disconnect();
 
  close(VCF_ANN);
  close(VCF_ALL);
}


=head2 import_vcf_2_db

 Title   : import_vcf_2_db
 Usage   : import_vcf_2_db( -database => 'name of the database,
                               );

 Function: Imports the VCF file content to the given table of variants in the database
					#Split different alternative alleles using ALT field
					#Put the variants information in the variants table if the comp_id is not present
					#Generates a VCF for variants that have to be annotated
					#Then inserts information from INFO field Separating using the comma, if there are diff alleles
					#Finally generates the Table for genotypes for different samples
 Returns : 

=cut			
sub import_vcf_2_db {
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $group_id = shift;
	my $log_file = shift;
	my $vcf_file = shift;
	my $vcf_to_ann = shift;
	my $vcf_all_var = shift;
	
  #print_and_log( "Opening VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  #Loads the VCF
	my $vcf = Vcf->new(file=>$vcf_file);
	#Reads the header
	$vcf->parse_header();
	
	#Field names and indexes
	my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	my $chr_ind = $cfg_hash->{'vcf_chr_ind'};
	my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	my $pos_ind = $cfg_hash->{'vcf_pos_ind'};
	my $id_str = $cfg_hash->{'db_vcf_id_str'};
	my $id_ind = $cfg_hash->{'vcf_id_ind'};
	my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	my $ref_ind = $cfg_hash->{'vcf_ref_ind'};
	my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	my $alt_ind = $cfg_hash->{'vcf_alt_ind'};
	my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	my $qual_ind = $cfg_hash->{'vcf_qual_ind'};
	my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	my $filter_ind = $cfg_hash->{'vcf_filter_ind'};
	my $info_str = $cfg_hash->{'db_vcf_info_str'};
	my $info_ind = $cfg_hash->{'vcf_info_ind'};
	my $genot_ind = $cfg_hash->{'vcf_genot_ind'};
	my $prob_gen_ind = $cfg_hash->{'vcf_prob_gen_ind'};
	my $comp_id_str = $cfg_hash->{'db_vcf_comp_id_str'};
	#Separators
	my $sep = ",";
	my $tab = "\t";
	my $vcf_info_sep = $cfg_hash->{'vcf_info_sep'};
	my $vcf_gen_sep = $cfg_hash->{'vcf_gen_sep'};
	
	my $var_inserted = 0;
	#VCF file of variants to be annotated
	open (VCF_ANN,">".$vcf_to_ann) or die "ERROR: Cannot open $vcf_to_ann. The program will exit..\n";
	print VCF_ANN  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";

	#VCF file with all variants
	open (VCF_ALL,">".$vcf_all_var) or die "ERROR: Cannot open $vcf_all_var. The program will exit..\n";
	print VCF_ALL  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";

	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
		
	#my $group_id =1;#GET THE GROUP ID....##########
	print_and_log( "Reading VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  # This will split the fields and print a list of CHR:POS
  while (my $x = $vcf->next_data_array()){
#	if ($var_inserted < 100) {#DEBUGCODE
		   #print_and_log( "$$x[$chr_ind]:$$x[$pos_ind]:$$x[$id_ind]:$$x[$ref_ind]:$$x[$alt_ind]\n",$log_file);#DEBUGCODE

		#An array to collect the ids corresponding to each alternate allele variant
		my @alt_ids = ();   
    ###################################################
    ####Adding variants and alternative variants ######
    ###################################################  
		#Construct the string with first information from the VCF
		my $fields = $chr_str.$sep.$pos_str.$sep.$id_str.$sep.$ref_str.$sep.$alt_str.$sep.$comp_id_str;
		#Construct the string with values
		my $start_values = "'".$$x[$chr_ind]."'".$sep."'".$$x[$pos_ind]."'".$sep."'".$$x[$id_ind]."'".$sep."'".$$x[$ref_ind]."'";
	      
		     
		#Split different alternative alleles using ALT field
		my @alts = split('\,', $$x[$alt_ind]);
		my $alt_id = 0;
		foreach my $alt (@alts){
	    #Construct the comp_id
	    my $comp_id = $$x[$chr_ind]."_".$$x[$pos_ind]."_".$$x[$ref_ind]."_".$alt;
	    my $values = $start_values.$sep."'".$alt."'".$sep."'".$comp_id."'";
	    #print_and_log( "Loading fields: $fields and  values: $values...\n",$log_file);#DEBUGCODE
	    #Put the variants information in the variants table if the comp_id is not present
	    my ($alt_id,$existing_id) = insert_only_if_not_exists_locked_woconn($dbh,$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'}
						    ,$fields,$values,$cfg_hash->{'db_vcf_comp_id_str'},"'".$comp_id."'"); 	
	    
	    #The variant was not there it has ben added and altid is the id
	    if ( $alt_id >= 0){
				#Print all variants
				print VCF_ALL join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$alt_id)."\n";
		    #Print the variant in the VCF if it must be annotated
		    #We put the var identifier in the INFO field that can be useful to insert annotation results
		    #my $var_id_fld = $cfg_hash->{'db_var_id'}."=".$alt_id;
		    print VCF_ANN join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$alt_id)."\n";
		    #Save variant ids in order
		    push(@alt_ids,$alt_id);
		    #print_and_log( "Variant was not there id:$alt_id...\n",$log_file);#DEBUGCODE
	    }elsif ($existing_id >= 0){
				#Print all variants
	      print VCF_ALL join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alt,$$x[$qual_ind],".",$existing_id)."\n";					  
		    #The variant was already there. Put the existing id
		    #print_and_log( "Variant was already there id:$existing_id...\n",$log_file);#DEBUGCODE
		    push(@alt_ids,$existing_id);
	    }else{
		    log_and_exit("ERROR: during the import of $comp_id. Cannot understand what is happen with the db...\n",$log_file);
	    }
	    #$alt_id++;#USED FOR TRIAL
		}
			    
			    
    ###################################################
    ####Inserting information from INFO field##########
    ###################################################
    
    #hash with alleles information
		my $alt_hash;
		#Split the INFO field
		my @info_fields = split($vcf_info_sep,$$x[$info_ind]);
    #For each field of the INFO column...
    foreach my $info_field (@info_fields){
		    my @info_f_parts = split("=",$info_field);
		    #Separate using the comma, if there are diff alleles
		    my @alt_infos = split(",",$info_f_parts[1]);
		    #If there is the alternative info, get it
		    for (my $alt_info_ind = 0; $alt_info_ind < scalar(@alt_infos); $alt_info_ind++){
			    #Store the  alternative info for each of them
			    #alt_hash->{1}->{DP} = 6
			    $alt_hash->{$alt_info_ind}->{$info_f_parts[0]} = $alt_infos[$alt_info_ind];
		    }
		}	
		
		#Taking filter and quality from the main line
		my $info_fields_base = $cfg_hash->{'db_qual'}.",".$cfg_hash->{'db_filter'}.",".$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'}.",";#.join(",",@info_fields);
		my $info_values_base = $$x[$qual_ind].",'".$$x[$filter_ind]."',$group_id";
		
		#Go through the hash with information for each allele
		foreach my $alt_info_ind ( keys %{$alt_hash}){
	    my $info_fields = $info_fields_base;
	    my $info_values = $info_values_base.",".$alt_ids[$alt_info_ind].",";
	    #Get all the INFO values for the allele
	    #The values are taken as numbers here
	    while (my ($key, $value) = each %{ $alt_hash->{$alt_info_ind} } ) {
			    #print "$key = $value \n";
			    $info_fields .= "$key,";
			    #if the value is not  a number enclose in quotes
			    if ( correct_type($value,"real") ){	
						$info_values .= "$value,";
					}else{$info_values .= "'".$value."',";}
			    
	    }
	    #Remove the last comma
	    chop($info_fields);
	    chop($info_values);
	    #print_and_log( "Loading fields: $info_fields and  values: $info_values...\n",$log_file);#DEBUGCODE
      my $alt_id = insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_var_statistics_table'},$info_fields,$info_values);		
		}
		
		###############################################################
		###########Table for genotypes for different samples###########
		###########from FORMAT field##################################
		###############################################################
		
		my $format_fields_base = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_var_ids'}.",";
		my $format_values_base = "";
		
		#Get the list of samples using the Vcf library
		my (@samples) = $vcf->get_samples();
		#Go through the genotype information per sample
    foreach my $sample (@samples){
	    my $format_fields_str = $format_fields_base;
	    
	    my $sample_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
						$cfg_hash->{'db_sample_name'}.",".$cfg_hash->{'db_analysis_id'},"'".$sample."'".",".$group_id);
	    #Joins the variant ids associated to different alleles if present
	    my $var_ids = join(",",@alt_ids);
	    my $format_values_str = $group_id.",".$sample_id.",'".$var_ids."',";
	    #Get the array with genotyping names
	    my @format_fields = split($vcf_gen_sep,$$x[$genot_ind]);
	    #Get the values corresponding to this sample
	    my $sample_col = $vcf->get_column($x, $sample);
	    #For each format field
	    foreach my $format_field (@format_fields){
		    #Get the index of the genotype information
		    my $form_idx = $vcf->get_tag_index($$x[$genot_ind],$format_field,':');
		    #If it exists
		    if ($form_idx != -1) {
				    #Add the field name
				    $format_fields_str .= $format_field.",";
				    #Get the value using the index
				    my $form_val = $vcf->get_field($sample_col,$form_idx);
				    #Concatenate the value
				    $format_values_str .= "'".$form_val."',";
		    }
		    ##If it is the GT information than use that to associate the allele variant
		    #if ( $format_field eq 'GT'){
				    ##Get the value using the index
				    #my $form_val = $vcf->get_field($sample_col,$form_idx);
				    #my @gen_flds = split(/[\|\/]/,$form_val);
		    #}
	    }
	    #Remove the last comma
	    chop($format_fields_str);
	    chop($format_values_str);
	    #print_and_log( "Loading fields: $format_fields_str and  values: $format_values_str...\n",$log_file);
	    #fill the table with genotypes for each sample
	    my $var_sample_id = insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genotype_sample_table'},$format_fields_str,$format_values_str);
	 }
#	}#DEBUGCODE 100vars
	
	$var_inserted++;
  }

 	$dbh->disconnect();
 
  close(VCF_ANN);
  close(VCF_ALL);
}




=head2 generate_vcf_to_annotate

 Title   : generate_vcf_to_annotate
 Usage   : generate_vcf_to_annotate( -database => 'name of the database,
                               );

 Function:
					#construct two files (allvar) and allvarinfo where
					#The first contains the information about the variant: chrom, pos, ref alt,qual
					#and the second contains the INFO fields needed and the FORMAT fields for all the samples		
		
					#Split different alternative alleles using ALT field
					#Put the variants information in the variants table if the comp_id is not present
					#Generates a VCF for variants that have to be annotated
					#Then inserts information from INFO field Separating using the comma, if there are diff alleles
					#Finally generates the Table for genotypes for different samples
					
					#In 2019 I added also the possibility to get the variant called by Freebayes. This tools
					#uses the AO and RO INFO fields of the VCF to annotate the reads count on alternate and reference
					#In this script when I meet the RO field than I search the AO too and merge them into a single AD
					#field separeting with comma.
					
			
 Returns : 

=cut			
sub generate_vcf_to_annotate {
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $vcf_file = shift;
	my $vcf_all_var = shift;
	
  #print_and_log( "Opening VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  #Loads the VCF
	my $vcf = Vcf->new(file=>$vcf_file);
	#Reads the header
	$vcf->parse_header();
	
	#Field names and indexes
	my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	my $chr_ind = $cfg_hash->{'vcf_chr_ind'};
	my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	my $pos_ind = $cfg_hash->{'vcf_pos_ind'};
	my $id_str = $cfg_hash->{'db_vcf_id_str'};
	my $id_ind = $cfg_hash->{'vcf_id_ind'};
	my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	my $ref_ind = $cfg_hash->{'vcf_ref_ind'};
	my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	my $alt_ind = $cfg_hash->{'vcf_alt_ind'};
	my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	my $qual_ind = $cfg_hash->{'vcf_qual_ind'};
	my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	my $filter_ind = $cfg_hash->{'vcf_filter_ind'};
	my $info_str = $cfg_hash->{'db_vcf_info_str'};
	my $info_ind = $cfg_hash->{'vcf_info_ind'};
	my $genot_ind = $cfg_hash->{'vcf_genot_ind'};
	my $prob_gen_ind = $cfg_hash->{'vcf_prob_gen_ind'};
	my $comp_id_str = $cfg_hash->{'db_vcf_comp_id_str'};
	#Separators
	my $sep = ",";
	my $tab = "\t";
	my $vcf_info_sep = $cfg_hash->{'vcf_info_sep'};
	my $vcf_gen_sep = $cfg_hash->{'vcf_gen_sep'};
	my $semicolon = ";";
	my $comma = ",";
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	my @info_fields_needed = split(",",$cfg_hash->{'info_fields_needed'});
	
	my $var_inserted = 0;


	
	########################################SORTED SAMPLES
	#Get the list of samples using the Vcf library (order is the one of the VCF)
	#my (@samples) = $vcf->get_samples();
	

	#The genotype info will be print in the order given in the user config file (samples_order)
	#Getting the samples list associated to the groupid
	my $samples_h;#Samples hash to be reordered
	my @sort_samples = ();#The array with samples id reordered
	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
	my @samples = ();
	#################VQSR Considerations
  #If VQSR was executed	then we get the samples from those selected by the user.
  #I check the genotype refinement step into the db
 	my $genotref = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_genotref'},
									$cfg_hash->{'db_analysis_id'},$analysis_id); 
	#print_and_log( "genotref = $genotref and runsamples: ".$cfg_hash->{'run_samples'}." \n",$log_file);									
  if ( ($genotref == 1) and ($cfg_hash->{'run_samples'} ne 'ALL') ){
		my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},",");
		@samples = split(",",$samples_ids);
	}else{
		foreach my $sample_id (keys %{$res_samples}){
			push(@samples,$sample_id);
		}
	}
	#################VQSR Considerations	
		
	my $kinship_ok = 1;
	#Get the kinship, to make the resorting
	foreach my $sample_name (@samples){
		
		
		#Obtain the kinship from the database given the sample id
		my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
														$cfg_hash->{'db_sample_name'},"'".$sample_name."'");	
		#Obtain the kinship from the database given the sample id
		my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
													#	$cfg_hash->{'db_sample_id'},$sample_id);	
													$cfg_hash->{'db_sample_name'},"'".$sample_name."'");
		#Build an hash to reorder
		#If at least one sample does not specifies kinship, set a flag
		if (defined $kinship ){
			if ( $kinship ne $cfg_hash->{'vargenius_empty_val'}){	
				$samples_h->{$sample_name}->{'id'} = $sample_name;
				$samples_h->{$sample_name}->{'k'} = $kinship;
			}else{
				$kinship_ok = 0;
			}
	  }else{
			$kinship_ok = 0;
		}
	}
	#Sort the samples as in samples_order parameter, if the kinship is present
	if ( $kinship_ok == 1){
		sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	}else{
		#Get the VCF order
		@sort_samples = @samples;
	}
	#print_and_log("Samples header: " ,$log_file);#DEBUGCODE
	#foreach my $sample (@sort_samples){	
		#print_and_log(" $sample " ,$log_file);#DEBUGCODE
	#}

	#Check if the target is extended 
	my $target_extended = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetextended'},
					$cfg_hash->{'db_analysis_id'},$analysis_id); 
	#Getting the consensus flag
	my $consensus = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_consensus'},
											$cfg_hash->{'db_analysis_id'},$analysis_id);
																
	#A file for the INFO and FORMAT field
	my $other_info_f = $vcf_all_var.".info";
	open (OTH_INFO,">".$other_info_f) or die "ERROR: Cannot open $other_info_f. The program will exit..\n";
	#VCF file with all variants
	open (VCF_ALL,">".$vcf_all_var) or die "ERROR: Cannot open $vcf_all_var. The program will exit..\n";
	print VCF_ALL  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";
	
	print_and_log( "Reading VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  # This will split the fields and print a list of CHR:POS
  while (my $x = $vcf->next_data_array()){
		#hash with alleles information
		my $alt_hash;

		#Try to understand from the FORMAT field if the result is obtained with freebayes
		my @format_fields_needed = split(",",$cfg_hash->{'genotype_fields'});
		if ( $$x[$genot_ind] =~ /:RO/ ){
			@format_fields_needed = split(",",$cfg_hash->{'FB_genotype_fields'});
		}
		
		#FILTERS ON THE VCF 
		#if ($var_inserted < 100) {#DEBUGCODE
		#When using the extended target and/or the consensus approach, hundreds of variants come out
		#hence I insert here a limit on the PASS filter
		if ( $cfg_hash->{'filter_vcf'} ){
			if ( $$x[$filter_ind] ne $cfg_hash->{'VCF_minqfilter'} or $$x[$qual_ind] <  $cfg_hash->{'VCF_minqual'}){
				next;
			}
		}
		#print_and_log( $$x[$chr_ind].":".$$x[$pos_ind].":".$$x[$id_ind].":".$$x[$ref_ind].":".$$x[$alt_ind]."\n",$log_file);#DEBUGCODE
		
		#An array to collect the ids corresponding to each alternate allele variant
		my @alt_ids = ();   
    ###################################################
    ####Adding variants and alternative variants ######
    ###################################################  

		#Split different alternative alleles using ALT field
		my @alts = split('\,', $$x[$alt_ind]);
		for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
		 #Insert all the needed info for the variant into the hash for the alternative
		  $alt_hash->{$alt_all}->{'VAR'} = join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alts[$alt_all],$$x[$qual_ind],$$x[$filter_ind]);
			#print_and_log( "VAR $alt_all: ".$alt_hash->{$alt_all}->{'VAR'}."\n",$log_file);#DEBUGCODE	 
		}
			    
   
		###################################################
		####Inserting information from INFO field##########
		###################################################
		##Split the INFO field
		#my @info_fields = split($vcf_info_sep,$$x[$info_ind]);
		##For each field of the INFO column...
		#foreach my $info_field (@info_fields){
				#my @info_f_parts = split("=",$info_field);
				#my $info_name = $info_f_parts[0];
				##If this info field is needed
				#if ( grep {/\b$info_name\b/} @info_fields_needed){
				
				##Separate using the comma, if there are diff alleles
				#my @alt_infos = split(",",$info_f_parts[1]);
				##If there is the alternative info, get it
				#for (my $alt_all = 0; $alt_all < scalar(@alt_infos); $alt_all++){
					##Store the  alternative info for each of them
					#if (defined $alt_hash->{$alt_ind}->{'INFO'} ){
							#$alt_hash->{$alt_all}->{'INFO'} = $alt_infos[$alt_all].$semicolon;
					#}else{
							#$alt_hash->{$alt_all}->{'INFO'} .= $alt_infos[$alt_all].$semicolon;
					#}
				#}
			#}
		#}	
		
		###################################################
		####Inserting information from INFO field##########
		###################################################		
		#Split the INFO field
		my @info_fields = split($vcf_info_sep,$$x[$info_ind]);
		my $info_fields_hash;
		#For each field of the INFO column...
		foreach my $info_field (@info_fields){
			my @info_f_parts = split("=",$info_field);
			$info_fields_hash->{$info_f_parts[0]} = $info_f_parts[1];
		}		
		#For each INFO field needed, if it is there write it..
		for ( my $i = 0; $i < scalar(@info_fields_needed); $i++){
			#If this info field is present
			if ( defined $info_fields_hash->{$info_fields_needed[$i]} ){
				
				#Separate for diff alleles
				my @alt_infos = split(",",$info_fields_hash->{$info_fields_needed[$i]});
				#If there is the alternative info, get it
				for (my $alt_all = 0; $alt_all < scalar(@alt_infos); $alt_all++){
					#Store the  alternative info for each of them (separated by comma)
					if (! defined $alt_hash->{$alt_all}->{'INFO'} ){
							$alt_hash->{$alt_all}->{'INFO'} = $alt_infos[$alt_all].$semicolon;
					}else{
							$alt_hash->{$alt_all}->{'INFO'} .= $alt_infos[$alt_all].$semicolon;
					}
				}
				#print_and_log( "FOR ".$info_fields_needed[$i].". Setting: ".$alt_hash->{0}->{'INFO'}."\n",$log_file);#DEBUGCODE
			}else{
				#Insert an empty value "-" when info is absent
				if (! defined $alt_hash->{0}->{'INFO'} ){
						$alt_hash->{0}->{'INFO'} = $vargenius_empty_val.$semicolon;
				}else{
						$alt_hash->{0}->{'INFO'} .= $vargenius_empty_val.$semicolon;
				}
				#print_and_log( "FOR ".$info_fields_needed[$i].". Setting: ".$alt_hash->{0}->{'INFO'}."\n",$log_file);#DEBUGCODE

			}
		}
			
		#Remove last comma
		for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			if (defined $alt_hash->{$alt_all}->{'INFO'}){
				#Remove last semicolon
				chop($alt_hash->{$alt_all}->{'INFO'});
				#print_and_log( "INFO $alt_all: ".$alt_hash->{$alt_all}->{'INFO'}."\n",$log_file);#DEBUGCODE
			}	
		}		
		
		###############################################################
		###########Table for genotypes for different samples###########
		###########from FORMAT field##################################
		###############################################################


		#For each alternative allele
		for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			my $mult_flag = 0;
			#There are multiple alleles?
			if ( scalar(@alts) > 1){
				$mult_flag = 1;
			}
			#Go through the genotype information per sample sorted. 
			#For each sample will concatenate its genotype information into the same hash location
			foreach my $sample (@sort_samples){
				my $genotype_info = "";
				
				#Get the array with FORMAT info from the VCF object
				#Get the column corresponding to this sample from the VCF object
				my $sample_col = $vcf->get_column($x, $sample);
				#For each format field (in this way I can also give the order)
				foreach my $format_field (@format_fields_needed){
					#If this format field is needed
					my $lc_f_f = uc $format_field;
						#Get the index of the specific genotype information
						my $form_idx = $vcf->get_tag_index($$x[$genot_ind],$lc_f_f,':');
						my $form_val = $cfg_hash->{'vargenius_empty_val'};
						#If it exists
						if ($form_idx != -1) {
							#Get the value for the specified sample using the index
							$form_val = $vcf->get_field($sample_col,$form_idx);
						}		
						#If we have GT and there are more alleles, we calculate the genotype using PL
						if ( $lc_f_f eq 'GT' and $mult_flag == 1){
								#Get the PL index
								my $pl_idx = $vcf->get_tag_index($$x[$genot_ind],'PL',':');
								if ($pl_idx != -1) {
									#Concatenate the value using the index
									$form_val .= $mult_ann_sep.$vcf->get_field($sample_col,$pl_idx);
								}
						}
						#Freebayes VCF:
						#If the field is RO (reference reads observed then get also the AO (alternative reads observed) concatenating them
						if ( $lc_f_f eq 'RO' ){
								#Get the AO index
								my $ao_idx = $vcf->get_tag_index($$x[$genot_ind],'AO',':');
								if ($ao_idx != -1) {
									#Concatenate the value of Alt allele counts using the index
									$form_val .= ",".$vcf->get_field($sample_col,$ao_idx);
								}
								#change the name of the genotype field (so that I do not change
								#the function get_genotype_info_vcf
								$lc_f_f = 'AD';
						}						
						#Given the format field (AD,GT,GQ and DP) converts or adds information:
						#For AD prints REF_reads, ALT_reads, tot_reads, perc_ALT
						#For GT, if the value is for a multiple allele, gets the genotype using the (just concatenated) PL value
						#For DP and GQ, will return the string containing the form_val
						$genotype_info = get_genotype_info_vcf($cfg_hash,$lc_f_f,$form_val,$alt_all,$mult_flag,$semicolon,$log_file);
						#print_and_log("result of get nenotype $genotype_info \n",$log_file);#DEBUGCODE
												
						if ( defined $alt_hash->{$alt_all}->{'FORMAT'}){
							$alt_hash->{$alt_all}->{'FORMAT'} .= $genotype_info.$semicolon;
						}else{
							$alt_hash->{$alt_all}->{'FORMAT'} =	$genotype_info.$semicolon;
						}					
				}
				#Remove last semicolon and insert the VarGenius separator ]---[
				if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
					chop($alt_hash->{$alt_all}->{'FORMAT'});
					$alt_hash->{$alt_all}->{'FORMAT'} .= $mult_ann_sep;
				}
			}	
			for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
				if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
					#print_and_log( "FORMAT $alt_all: ".$alt_hash->{$alt_all}->{'FORMAT'}."\n",$log_file);#DEBUGCODE
				}
			}
			

			
	 }
	#	}#DEBUGCODE 100vars


		#Now construct two files (allvar) and allvarinfo where
		#The first contains the information about the variant: chrom, pos, ref alt,qual
		#and the second contains the INFO fields needed and the FORMAT fields for all the samples
		for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			#Final string to write in the file for each alt allele
			my $final_string = "";
			#String with info to write in the otherinfo file
			my $info_string = "";
			
			#concatenate the variant information and the index of the variant
			if (defined $alt_hash->{$alt_all}->{'VAR'}){
				$final_string .=  $alt_hash->{$alt_all}->{'VAR'}.$tab.$var_inserted;
			}else{
					die "ERROR: No information for the VAR. That is weird!\n";
			}
			
			#Get the information needed from the INFO field
			#first write the index of the variant
			$info_string .= $var_inserted.$tab;
			if (defined $alt_hash->{$alt_all}->{'INFO'}){
				$info_string .= $alt_hash->{$alt_all}->{'INFO'}.$mult_ann_sep;
			}else{
				#print_and_log("WARNING: No information for the INFO. \n",$log_file);#DEBUGCODE
				$info_string .= "-".$mult_ann_sep;	
			}
			if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
				$info_string .= $alt_hash->{$alt_all}->{'FORMAT'};
			}else{
					die "ERROR: No information for the FORMAT. That is weird!\n";
			}
			#Remove the separator at the end of the string
			$info_string = substr($info_string,0,-length($mult_ann_sep));
			#Print a different variant for each allele		
			print VCF_ALL $final_string."\n";
			print OTH_INFO $info_string."\n"; 
			$var_inserted++;
		}
	
  }
  close(OTH_INFO);
  close(VCF_ALL);
}



#=head2 generate_vcf_to_annotate

 #Title   : generate_vcf_to_annotate
 #Usage   : generate_vcf_to_annotate( -database => 'name of the database,
                               #);

 #Function:
					##construct two files (allvar) and allvarinfo where
					##The first contains the information about the variant: chrom, pos, ref alt,qual
					##and the second contains the INFO fields needed and the FORMAT fields for all the samples		
		
					##Split different alternative alleles using ALT field
					##Put the variants information in the variants table if the comp_id is not present
					##Generates a VCF for variants that have to be annotated
					##Then inserts information from INFO field Separating using the comma, if there are diff alleles
					##Finally generates the Table for genotypes for different samples
					
			
 #Returns : 

#=cut			
#sub generate_vcf_to_annotate {
	#my $cfg_hash = shift;
	#my $prog_used = shift;
	#my $group_id = shift;
	#my $log_file = shift;
	#my $vcf_file = shift;
	#my $vcf_all_var = shift;
	
  ##print_and_log( "Opening VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  ##Loads the VCF
	#my $vcf = Vcf->new(file=>$vcf_file);
	##Reads the header
	#$vcf->parse_header();
	
	##Field names and indexes
	#my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	#my $chr_ind = $cfg_hash->{'vcf_chr_ind'};
	#my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	#my $pos_ind = $cfg_hash->{'vcf_pos_ind'};
	#my $id_str = $cfg_hash->{'db_vcf_id_str'};
	#my $id_ind = $cfg_hash->{'vcf_id_ind'};
	#my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	#my $ref_ind = $cfg_hash->{'vcf_ref_ind'};
	#my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	#my $alt_ind = $cfg_hash->{'vcf_alt_ind'};
	#my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	#my $qual_ind = $cfg_hash->{'vcf_qual_ind'};
	#my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	#my $filter_ind = $cfg_hash->{'vcf_filter_ind'};
	#my $info_str = $cfg_hash->{'db_vcf_info_str'};
	#my $info_ind = $cfg_hash->{'vcf_info_ind'};
	#my $genot_ind = $cfg_hash->{'vcf_genot_ind'};
	#my $prob_gen_ind = $cfg_hash->{'vcf_prob_gen_ind'};
	#my $comp_id_str = $cfg_hash->{'db_vcf_comp_id_str'};
	##Separators
	#my $sep = ",";
	#my $tab = "\t";
	#my $vcf_info_sep = $cfg_hash->{'vcf_info_sep'};
	#my $vcf_gen_sep = $cfg_hash->{'vcf_gen_sep'};
	#my $semicolon = ";";
	#my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};


	#my @format_fields_needed = split(",",$cfg_hash->{'genotype_fields'});
	#my @info_fields_needed = split(",",$cfg_hash->{'info_fields_needed'});
	
	#my $var_inserted = 0;


	
	#########################################SORTED SAMPLES
	##Get the list of samples using the Vcf library (order is the one of the VCF)
	##my (@samples) = $vcf->get_samples();
	

	##The genotype info will be print in the order given in the user config file (samples_order)
	##Getting the samples list associated to the groupid
	#my $samples_h;#Samples hash to be reordered
	#my @sort_samples = ();#The array with samples id reordered
	##Get the samples ids involved for the group
	#my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$group_id;";	
	##print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	#my $res_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
	#my @samples = ();
	##################VQSR Considerations
  ##If VQSR was executed	then we get the samples from those selected by the user.
  ##I check the genotype refinement step into the db
 	#my $genotref = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_genotref'},
									#$cfg_hash->{'db_analysis_id'},$group_id); 
	##print_and_log( "genotref = $genotref and runsamples: ".$cfg_hash->{'run_samples'}." \n",$log_file);									
  #if ( ($genotref == 1) and ($cfg_hash->{'run_samples'} ne 'ALL') ){
		#my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},",");
		#@samples = split(",",$samples_ids);
	#}else{
		#foreach my $sample_id (keys %{$res_samples}){
			#push(@samples,$sample_id);
		#}
	#}
	##################VQSR Considerations	
		
	#my $kinship_ok = 1;
	##Get the kinship, to make the resorting
	#foreach my $sample_name (@samples){
		
		
		##Obtain the kinship from the database given the sample id
		#my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
														#$cfg_hash->{'db_sample_name'},"'".$sample_name."'");	
		##Obtain the kinship from the database given the sample id
		#my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
													##	$cfg_hash->{'db_sample_id'},$sample_id);	
													#$cfg_hash->{'db_sample_name'},"'".$sample_name."'");
		##Build an hash to reorder
		##If at least one sample does not specifies kinship, set a flag
		#if (defined $kinship ){
			#if ( $kinship ne $cfg_hash->{'vargenius_empty_val'}){	
				#$samples_h->{$sample_name}->{'id'} = $sample_name;
				#$samples_h->{$sample_name}->{'k'} = $kinship;
			#}else{
				#$kinship_ok = 0;
			#}
	  #}else{
			#$kinship_ok = 0;
		#}
	#}
	##Sort the samples as in samples_order parameter, if the kinship is present
	#if ( $kinship_ok == 1){
		#sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	#}else{
		##Get the VCF order
		#@sort_samples = @samples;
	#}
	##print_and_log("Samples header: " ,$log_file);#DEBUGCODE
	##foreach my $sample (@sort_samples){	
		##print_and_log(" $sample " ,$log_file);#DEBUGCODE
	##}
	
	##A file for the INFO and FORMAT field
	#my $other_info_f = $vcf_all_var.".info";
	#open (OTH_INFO,">".$other_info_f) or die "ERROR: Cannot open $other_info_f. The program will exit..\n";
	##VCF file with all variants
	#open (VCF_ALL,">".$vcf_all_var) or die "ERROR: Cannot open $vcf_all_var. The program will exit..\n";
	#print VCF_ALL  join($tab,$chr_str,$pos_str,$id_str,$ref_str,$alt_str,$qual_str,$filter_str,$info_str)."\n";
	
	#print_and_log( "Reading VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  ## This will split the fields and print a list of CHR:POS
  #while (my $x = $vcf->next_data_array()){
		##hash with alleles information
		#my $alt_hash;

##	if ($var_inserted < 100) {#DEBUGCODE
		##print_and_log( $$x[$chr_ind].":".$$x[$pos_ind].":".$$x[$id_ind].":".$$x[$ref_ind].":".$$x[$alt_ind]."\n",$log_file);#DEBUGCODE

		##An array to collect the ids corresponding to each alternate allele variant
		#my @alt_ids = ();   
    ####################################################
    #####Adding variants and alternative variants ######
    ####################################################  

		##Split different alternative alleles using ALT field
		#my @alts = split('\,', $$x[$alt_ind]);
		#for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
		 ##Insert all the needed info for the variant into the hash for the alternative
		  #$alt_hash->{$alt_all}->{'VAR'} = join($tab,$$x[$chr_ind],$$x[$pos_ind],$$x[$id_ind],$$x[$ref_ind],$alts[$alt_all],$$x[$qual_ind],$$x[$filter_ind]);
			##print_and_log( "VAR $alt_all: ".$alt_hash->{$alt_all}->{'VAR'}."\n",$log_file);#DEBUGCODE	 
		#}
			    
   
    ####################################################
    #####Inserting information from INFO field##########
    ####################################################
		##Split the INFO field
		#my @info_fields = split($vcf_info_sep,$$x[$info_ind]);
    ##For each field of the INFO column...
    #foreach my $info_field (@info_fields){
			##If this info field is needed
			#if ( grep {/\b$info_field\b/} @info_fields_needed){
		    #my @info_f_parts = split("=",$info_field);
		    ##Separate using the comma, if there are diff alleles
		    #my @alt_infos = split(",",$info_f_parts[1]);
		    ##If there is the alternative info, get it
		    #for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			    ##Store the  alternative info for each of them
			    #if (defined $alt_hash->{$alt_ind}->{'INFO'} ){
						#$alt_hash->{$alt_all}->{'INFO'} = $alt_infos[$alt_all].$semicolon;
					#}else{
						#$alt_hash->{$alt_all}->{'INFO'} .= $alt_infos[$alt_all].$semicolon;
					#}
		    #}
			#}
		#}	
		##for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			##print_and_log( "INFO $alt_all: ".$alt_hash->{$alt_all}->{'INFO'}."\n",$log_file);#DEBUGCODE
		##}		
		
		################################################################
		############Table for genotypes for different samples###########
		############from FORMAT field##################################
		################################################################


		##For each alternative allele
		#for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			#my $mult_flag = 0;
			##There are multiple alleles?
			#if ( scalar(@alts) > 1){
				#$mult_flag = 1;
			#}
			##Go through the genotype information per sample sorted. 
			##For each sample will concatenate its genotype information into the same hash location
			#foreach my $sample (@sort_samples){
				#my $genotype_info = "";
				
				##Get the array with FORMAT info from the VCF object
				##Get the column corresponding to this sample from the VCF object
				#my $sample_col = $vcf->get_column($x, $sample);
				##For each format field (in this way I can also give the order)
				#foreach my $format_field (@format_fields_needed){
					##If this format field is needed
					#my $lc_f_f = uc $format_field;
						##Get the index of the specific genotype information
						#my $form_idx = $vcf->get_tag_index($$x[$genot_ind],$lc_f_f,':');
						#my $form_val = $cfg_hash->{'vargenius_empty_val'};
						##If it exists
						#if ($form_idx != -1) {
							##Get the value for the specified sample using the index
							#$form_val = $vcf->get_field($sample_col,$form_idx);
						#}		
						##If we have GT and there are more alleles, we calculate the genotype using PL
						#if ( $lc_f_f eq 'GT' and $mult_flag == 1){
								##Get the PL index
								#my $pl_idx = $vcf->get_tag_index($$x[$genot_ind],'PL',':');
								#if ($pl_idx != -1) {
									##Concatenate the value using the index
									#$form_val .= $mult_ann_sep.$vcf->get_field($sample_col,$pl_idx);
								#}
						#}
						##Given the format field (AD,GT,GQ and DP) converts or adds information:
						##For AD prints REF_reads, ALT_reads, tot_reads, perc_ALT
						##For GT, if the value is for a multiple allele, gets the genotype using the (just concatenated) PL value
						##For DP and GQ, will return the string containing the form_val
						#$genotype_info = get_genotype_info_vcf($cfg_hash,$lc_f_f,$form_val,$alt_all,$mult_flag,$semicolon,$log_file);
						##print_and_log("result of get nenotype $genotype_info \n",$log_file);#DEBUGCODE
												
						#if ( defined $alt_hash->{$alt_all}->{'FORMAT'}){
							#$alt_hash->{$alt_all}->{'FORMAT'} .= $genotype_info.$semicolon;
						#}else{
							#$alt_hash->{$alt_all}->{'FORMAT'} =	$genotype_info.$semicolon;
						#}					
				#}
				##Remove last semicolon and insert the VarGenius separator ]---[
				#if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
					#chop($alt_hash->{$alt_all}->{'FORMAT'});
					#$alt_hash->{$alt_all}->{'FORMAT'} .= $mult_ann_sep;
				#}
			#}	
			#for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
				#if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
					##print_and_log( "FORMAT $alt_all: ".$alt_hash->{$alt_all}->{'FORMAT'}."\n",$log_file);#DEBUGCODE
				#}
			#}
			

			
	 #}
	##	}#DEBUGCODE 100vars


		##Now construct two files (allvar) and allvarinfo where
		##The first contains the information about the variant: chrom, pos, ref alt,qual
		##and the second contains the INFO fields needed and the FORMAT fields for all the samples
		#for (my $alt_all = 0; $alt_all < scalar(@alts); $alt_all++){
			##Final string to write in the file for each alt allele
			#my $final_string = "";
			##String with info to write in the otherinfo file
			#my $info_string = "";
			
			##concatenate the variant information and the index of the variant
			#if (defined $alt_hash->{$alt_all}->{'VAR'}){
				#$final_string .=  $alt_hash->{$alt_all}->{'VAR'}.$tab.$var_inserted;
			#}else{
					#die "ERROR: No information for the VAR. That is weird!\n";
			#}
			
			##Get the information needed from the INFO field
			##first write the index of the variant
			#$info_string .= $var_inserted.$tab;
			#if (defined $alt_hash->{$alt_all}->{'INFO'}){
				#$info_string .= $alt_hash->{$alt_all}->{'INFO'}.$mult_ann_sep;
			#}else{
				##print_and_log("WARNING: No information for the INFO. \n",$log_file);#DEBUGCODE
				#$info_string .= "-".$mult_ann_sep;	
			#}
			#if (defined $alt_hash->{$alt_all}->{'FORMAT'}){
				#$info_string .= $alt_hash->{$alt_all}->{'FORMAT'};
			#}else{
					#die "ERROR: No information for the FORMAT. That is weird!\n";
			#}
			##Remove the separator at the end of the string
			#$info_string = substr($info_string,0,-length($mult_ann_sep));
			##Print a different variant for each allele		
			#print VCF_ALL $final_string."\n";
			#print OTH_INFO $info_string."\n"; 
			#$var_inserted++;
		#}
	
  #}
  #close(OTH_INFO);
  #close(VCF_ALL);
#}






=head2 output_from_annotation_without_db

 Title   : output_from_annotation_without_db
 Usage   : output_from_annotation_without_db( -database => 'name of the database,
                               );

 Function:  Uses the output from Annovar to build a first output for VarGenius.
						It separates the information for the transcripts and also adds the genotype information.
 						-gets the needed fields for variants,statistics and genotype tables
 						
						
						Then a final output in the form of a tab separated table is constructed:
						- an Header is built using both data from the configuration file
						- gene information are in the table GENES and will be fetched associated to genes symbols
						- genotype information is fetched related to variants from the GENOTYPE_SAMPLE table hence
								   fields names for this information is collected too
						
							
 Returns : 

=cut
sub output_from_annotation_without_db {
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $vargenius_out = shift;
	my $compl_annovar_out = shift;
	my $other_info_f = shift;
	my $vcf_file = shift;
	
	my $semicolon = ";";
	my $tab = "\t";
	my $comma = ",";
	
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	#Annotations fields are obtained from the program_config and are the fields from Annovar output that
	#will be used.
	my @annotations_fields = split($comma,$cfg_hash->{'annotations_fields'});
	my $ann_fields_used = 0;#Used to check if we are using all. Otherwise there is an error
	
	#The total number of information carried from Annovar in the Otherinfo column
	my $other_info_col = $cfg_hash->{'other_info_col'};
	my $tot_info_annovar = $cfg_hash->{'tot_info_annovar'};
	my $start_info_annovar = $cfg_hash->{'start_info_ann_wo_db'};
	my $annov_empty_val = $cfg_hash->{'annov_empty_val'};

	
	#Inizialization of Variants, info and format fields to pick from the output of annovar
	my @variants_fields_needed = ();
	if (defined $cfg_hash->{'variant_fields_from_annovar'}){ @variants_fields_needed = split($comma,$cfg_hash->{'variant_fields_from_annovar'});}	
	my @info_fields_needed = ();
	if (defined $cfg_hash->{'info_fields_from_annovar'}){ @info_fields_needed = split($comma,$cfg_hash->{'info_fields_from_annovar'});}
	my @format_fields_needed = ();
	if (defined $cfg_hash->{'format_fields_from_annovar'}){ @format_fields_needed =split($comma,$cfg_hash->{'format_fields_from_annovar'});}
	
	#Index in the list of the refgene_symbol, useful later when we have to access gene info
	my $gene_field_ind = -1;
	#print_and_log( "Printing the VarGenius outputs...\n",$log_file); #DEBUGCODE
	
	#Positions of the genotype information
	my $gt_pos = $cfg_hash->{'gt_pos'};
	my $zyg_pos = $cfg_hash->{'zyg_pos'};
	my $ref_r_pos = $cfg_hash->{'ref_r_pos'};
	my $alt_r_pos = $cfg_hash->{'alt_r_pos'};
	my $tot_r_pos = $cfg_hash->{'tot_r_pos'};
	my $percalt_r_pos = $cfg_hash->{'percalt_r_pos'};
	my $dp_pos = $cfg_hash->{'dp_pos'};
	my $gq_pos = $cfg_hash->{'gq_pos'};
	
	
	####SETTINGS FOR BED FILES IN ANNOVAR
	#In Annovar, if you use a BED file you must write the name with _BED into ##Annovar DATABASES in program_config.txt
	#then Annovar when annotates uses the column name "bed" for the first element in the annotation, then "bed2,bed3,..,bedN".
	#We change the header name accordingly here bulding an hash containing for each BED Annovar header  the corresponding database name
	#First get the XX_BED from the three parameters annov_g_dbs,annov_f_dbs, annov_r_dbs (THE SAME ORDER IS IMPORTANT)
	my @ann_all_dbs = ();
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_g_dbs'}));
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_f_dbs'}));
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_r_dbs'}));
	#The smartest way to fetch just those with _BED is the following
	my @ann_bed_dbs = grep(/_BED/, @ann_all_dbs);
	#Now remove the _BED and put it into an hash
	my $ann_bed_names;
	my $bed_counter = 1;
	foreach my $ann_bed_db (@ann_bed_dbs){
			$ann_bed_db =~ s/_BED//g;
			if ($bed_counter > 1){
				$ann_bed_names->{$cfg_hash->{'bed_ext'}.$bed_counter} = $ann_bed_db;
			}else{
				$ann_bed_names->{$cfg_hash->{'bed_ext'}} = $ann_bed_db;	
			}
			$bed_counter++;
	}
	
	##### FIRST PART HEADER#########
	###Here I first save the header for the variants, statistics and genotype info then
	# the header for the annotation will be added during the read of the annovar output

	#HEADER GENOTYPE INFOS
	#The genotype info will be print in the order given in the user config file (samples_order)

	#Getting the samples list associated to the analysisid
  my $samples_h;#Samples hash to be reordered
	my @sort_samples = ();#The array with samples id reordered
	#Get the samples ids involved for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res_group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	my @group_sam = ();
	#################VQSR Considerations
	#If VQSR was executed	then we get the samples from those selected by the user.
	#I check the genotype refinement step into the db
 	my $genotref = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_genotref'},
									$cfg_hash->{'db_analysis_id'},$analysis_id); 
	#print_and_log( "genotref = $genotref and runsamples: ".$cfg_hash->{'run_samples'}." \n",$log_file);									
	if ( ($genotref == 1) and ($cfg_hash->{'run_samples'} ne 'ALL') ){
		my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},",");
		@group_sam = split(",",$samples_ids);
	}else{
		foreach my $sample_id (keys %{$res_group_sam}){
			push(@group_sam,$sample_id);
		}
	}
	#################VQSR Considerations	
	
	#Print the header for information about the genotype for the sample
	my $sample_gen_fields1 = "";
	my $sample_gen_fields2 = "";
	my $kinship_ok = 1;
	#Get the kinship, to make the resorting
	foreach my $sample_id (@group_sam){
		#Obtain the kinship from the database given the sample id
		my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
														$cfg_hash->{'db_sample_id'},$sample_id);	
		#Obtain the sample name from the database given the sample id
		my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
														$cfg_hash->{'db_sample_id'},$sample_id);	
		#Build an hash to reorder
		#If at least one sample does not specifies kinship, set a flag
		if (defined $kinship ){
			if ( $kinship ne $cfg_hash->{'vargenius_empty_val'}){	
				$samples_h->{$sample_id}->{'id'} = $sample_name;
				$samples_h->{$sample_id}->{'k'} = $kinship;
			}else{
				$kinship_ok = 0;
			}
	  }else{
			$kinship_ok = 0;
		} 
	}
	
	#Sort the samples as in samples_order parameter, if the kinship is present
	if ( $kinship_ok == 1){
		#Sort the samples as in samples_order parameter
		sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	}#Otherwise Get the samples order as in the VCF
	else{
		my $vcf = Vcf->new(file=>$vcf_file);
		$vcf->parse_header(); 
		my (@samples_vcf) = $vcf->get_samples();
		#Get the VCF order
		@sort_samples = @samples_vcf;
	}

	
	#Initialize two strings for two sets of information for genotype
	my $gt_add_1 = "";
	my $gt_add_2 = "";
	#print_and_log("Output from annotation samples header: " ,$log_file);#DEBUGCODE
	#Write header for genotypes. Go through sample names
	foreach my $sample_name (@sort_samples){
		print_and_log("Getting genotype header for the sample : $sample_name..\n" ,$log_file);#DEBUGCODE
		#print_and_log(" $sample_name " ,$log_file);#DEBUGCODE
		
		$sample_gen_fields1 .= "$sample_name\_".$cfg_hash->{'ref_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'alt_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'tot_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'percalt_r_str'}.$tab;
		
		$gt_add_1 .= "$sample_name\_".$cfg_hash->{'gt_str'}.
								 "$tab$sample_name\_".$cfg_hash->{'zyg_str'}.$tab;
		$gt_add_2 .= "$sample_name\_".$cfg_hash->{'dp_str'}.
								 "$tab$sample_name\_".$cfg_hash->{'gq_str'}.$tab;				
						
	}
	chop($sample_gen_fields1);		
	$sample_gen_fields2 .= $gt_add_1.$gt_add_2;
	chop($sample_gen_fields2);		
	#END HEADER GENOTYPE INFOS
	
	
	###############################
	# INFORMATION FROM THE DB INTO AN HASH
	####################################
	#Getting the type of sequencing. The frequencies table will be updated for that only
	my $seqtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
									
	my $db_hash;
	#Store into memory with an hash useful information from the database to be used for the annotation
	my $compid_cols = $cfg_hash->{'vcf_chr_ind'}.",".$cfg_hash->{'vcf_pos_ind'}.",".$cfg_hash->{'vcf_ref_ind'}.",".$cfg_hash->{'vcf_alt_ind'};
	#Get only the columns with the information with the VCF
	copy ($vcf_file, $vcf_file.".temp" ) or die "ERROR: Cannot copy $vcf_file into $vcf_file.onlyinfo \n";
	delete_rows_containing($vcf_file.".temp","^#",$vcf_file.".onlyinfo");
	#Get only the variants chr,pos,alt,ref
	extract_columns_from_file($vcf_file.".onlyinfo",$compid_cols,$vcf_file.".compid");
	#Put into an array
	my @compids = list_to_array($vcf_file.".compid",'NO_NEW_LINE');
	#Generate a string to be used into a query
	my $query_str_compids = "";
	foreach my $compid (@compids){
		$compid =~ s/\t/\_/g;
		$query_str_compids .= "'$compid',";
	}
	chop($query_str_compids);
	
	##ALLELIC FREQUENCIES
	#Get all the variant frequencies needed using a query for which all the compid are used
	print_and_log( "Get all the variant ids\n",$log_file);
	$query = "SELECT ".$cfg_hash->{'db_var_compid'}.",".$cfg_hash->{"db_$seqtype\_allele_freq"}.",".$cfg_hash->{"db_$seqtype\_freq_factors"}." FROM ".$cfg_hash->{'db_variants_table'}.
					" WHERE ".$cfg_hash->{'db_var_compid'}." IN ($query_str_compids);";							
	#print_and_log("Executing: $query\n",$log_file);
	print_and_log("Executing a query to fetch all the variants frequencies for this analysis from $vcf_file\n",$log_file);

	#fetch_all_rows
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
										$cfg_hash->{'db_pass'},$query);
	foreach my $var_field (@$res) {
		my @arr = @$var_field;
		#Get the database id of the variant
		my $compid = $arr[0];
		my $allele_freq = $arr[1];
		my $freq_factors = $arr[2];
		
		$db_hash->{$compid}->{'allele_freq'} = $allele_freq;
		$db_hash->{$compid}->{'freq_factors'} = $freq_factors;
	}	

	##########################################
	# END INFORMATION FROM THE DB INTO AN HASH
	####################################
		
	###########################################
	## PRINT THE FIRST PART OF THE HEADER   ###
	###########################################
	#Print the header using the variants and statistics fields and the genotype information
	#(the second part of the header consists in the annotation of Annovar and is written during the
	#annovar output reading later)
	my $head_var_stat_gen = "";
	if ((scalar @variants_fields_needed) > 0){ $head_var_stat_gen .= join($tab,@variants_fields_needed).$tab;}
	if ((scalar @info_fields_needed) > 0){ $head_var_stat_gen .= join($tab,@info_fields_needed).$tab;}
	if ($sample_gen_fields1 ne ''){ $head_var_stat_gen .= $sample_gen_fields1.$tab; }
	if ($sample_gen_fields2 ne ''){ $head_var_stat_gen .= $sample_gen_fields2.$tab; }
	
	#Build an hash with all the INFO and FORMAT fields for each variant
	my $oth_info_hash;
	open (OTH_INFO,"<".$other_info_f) or die "ERROR: Cannot open $other_info_f. The program will exit..\n";
	while (my $line = <OTH_INFO>){
		chop($line);
		my @fields = split($tab,$line);
		#Put in the hash at the variant index location the corresponding info prepared in generate_vcf
		$oth_info_hash->{$fields[0]} = $fields[1];
	}
	close(OTH_INFO);
	######END FIRST PART HEADER SAVING#####
	
	open (FILT_OUT,">".$vargenius_out) or die "ERROR: Cannot open $vargenius_out. The program will exit..\n";
		
	#####################################
	###########PRINT OUTPUT##############
	#####################################	

	#Open output from the annotation
	open (ANN_OUT,"<".$compl_annovar_out) or die "ERROR: Cannot open $compl_annovar_out. The program will exit..\n";
	#RefSeq fields: contains the transcripts information.Convert it in a readable form
	my $tr_ann_field_col = -1;
	######################
	#Genomic annotation (RefSeq, Gencode, ..) inizialization
	#Get the field of Annovar annotation  describing the AAChange (the field that must be splitted)
	my $aa_change_gen_field = $cfg_hash->{'aa_change_gen_field'};
	$aa_change_gen_field =~ s/[\-\.\+]/\_/g;
	$aa_change_gen_field = lc($aa_change_gen_field);
	#Genomic annotation fields to use (RefSeq, Gencode, etc)
	my $genomic_ann_fields = $cfg_hash->{'genomic_ann_fields'};
	$genomic_ann_fields =~ s/[\-\.\+]/\_/g;
	$genomic_ann_fields = lc($genomic_ann_fields);	
	my @genomic_ann_fields = split(",",$genomic_ann_fields);
	my @compl_gen_ann_fields = ();#Concatenate AAChange with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_gen_ann_fields,$aa_change_gen_field."_".$genomic_ann_field);
	}
	my @compl_gen_ann_cols = ();#Used for the indexes of these fields
	##Get the field names that will be used for splitting the AAChange field
	#my $tr_ann_field = $cfg_hash->{'tr_ann_field'};
	#$tr_ann_field =~ s/[\-\.\+]/\_/g;
	#$tr_ann_field = lc($tr_ann_field);
	#################

	#The index of the column Otherinfo from Annovar (with varid)
	my $other_info_col_ind = -1;
	#Information required to filter base on the frequencies in ExAC
	my $exacall_field = $cfg_hash->{'exacall_fld'};
	my $exacall_ind = -1;
	my $exacall_min_freq = $cfg_hash->{'exacall_min_freq'};
	
	#Numerates the line just to know what is the header of the file
	my $num_line = 0;
	my $tot_cols = 0;
	#The array of field names
	my @head_flds = ();
	#the indexes of the fields that are to be used
	my @flds_touse = ();
	
				
	#Go through each line of a TXT tab separated table
	while (my $line = <ANN_OUT>){ # and ($num_line < 1000){#DEBUGCODE
		#Remove last \n
		chop($line);
		#Initialize the newline
		my $newline = "";

		#Used to filter the variants
		my $take_var = 1;
		
		###########################
		##########HEADER###########
		###########################
		#If we are in the header
		if ($num_line == 0){
			#Put the first part of the header
			$newline .= $head_var_stat_gen;
			#Split the header to get the fields names
			@head_flds = split($tab,$line);
			$tot_cols = scalar(@head_flds);
			#Number of column that will be inserted in the flds_touse
			my $num_col = 0;
			#Keep in count if a field has been already inserted
			my @insd_flds = ();
			#Get each field in the fields variable
			foreach my $head_fld (@head_flds){
				
				#print_and_log( "Headfield:$head_fld \n",$log_file);
				#Other info is the last field with supplementary information
				#and contains the variant id
				if ($head_fld eq $other_info_col) {
					#print_and_log( "This is otherinfo:$head_fld \n",$log_file);
					$other_info_col_ind = $num_col;
					#push(@flds_touse,$num_col);
				}else{
						#Before to add the value I search for those strings
						#that cannot be inserted as they are written in the output from Annovar
						#
						#	print_and_log( "Field $head_fld\n",$log_file);
						#Convert all - and . into _
						$head_fld =~ s/[\-\.\+]/\_/g;
						#These are keywords in SQL, change as before
						if ( ($head_fld eq 'Start') or ($head_fld eq 'End') or ($head_fld eq 'Ref') ){
							$head_fld = "_".$head_fld;
						}
						#Convert upper to lowercase
						$head_fld = lc( $head_fld);
						#If the field starts with a number (e.g. 1000g), add 'db'
						if ($head_fld =~ /^\d/){ $head_fld= "db".$head_fld;}					
						
						#Here I filter the annotation using the frequency of Exac. 
						#This was useful for the Undiagnose Diseases. Being undiagnosed, the variants
						#are interesting only if they are not frequent in Exac
						if ( $head_fld eq $exacall_field){
							$exacall_ind = $num_col; 
							print "Exacall $head_fld ind :$exacall_ind\n"
						}
						
						#If this field is needed then it will be inserted
						if ( grep {/\b$head_fld\b/} @annotations_fields){
							#SECOND PART OF THE HEADER IS WRITTEN HERE
							$ann_fields_used++;#To check if we use all
							#check if this field was already inserted
							if ( ! grep {/\b$head_fld\b/} @insd_flds) {
								#Columns for which the name must be substituted
								#If the column is one among bed..bedN
								if ( $head_fld =~ /^bed\d*$/ ){
									$newline .= $ann_bed_names->{$head_fld}.$tab;

								}#NB:Add elsif if there are others..
								#Columns for which something must be added or they must be kept as they are
								else{
									#Keep the index of the field that keeps the RefSeq information
									#and adds the fields names for the transcript infos
									if ( grep {/\b$head_fld\b/} @compl_gen_ann_fields ){
										#Get the column index
										push(@compl_gen_ann_cols,$num_col);
										#Get the genomic annotation name removing the aachange field
										$head_fld =~ s/^$aa_change_gen_field\_// ;
										#Get the header by summing the head field obtained with all the annotation fields
										my $tr_ann_fields = $cfg_hash->{'tr_ann_fields'};
										my @tr_ann_fields = split (",",$tr_ann_fields);
										foreach my $tr_ann_field (@tr_ann_fields){
											$newline .= $head_fld."_".$tr_ann_field.$tab;
										}							
									}									
																	
									#This will be applied for any other field!
									$newline .= $head_fld.$tab;
								}
								push(@insd_flds,$head_fld);
								#Save the indexes to use
								push(@flds_touse,$num_col);	
							}else{
								print_and_log("WARNING: File $compl_annovar_out contains a repeated field: $head_fld.".
												" The column number $num_col will not be used..\n",$log_file);
							}
						}
				}
				
				$num_col++;
			}

			##Check if all the annotations field are used, otherwise stop with error
			#unless ( scalar( @annotations_fields) == $ann_fields_used){
				#log_and_exit("ERROR: File $ann_fields_used annotation fields used while they are: ".scalar( @annotations_fields).".".
				#" Check variable annotations_fields into the configuration file. Exiting...\n",$log_file);				
			#}
			#NEW INFORMATION HEADER ( N.B. the last tab will be removed and a newline inserted at the end!)
			#Add the internal frequency field 
			$newline .= $cfg_hash->{'intern_var_freq_fld'}.$tab.$cfg_hash->{'var_frequency_factors_fld'}.$tab;
			
		}##########HEADER###########
		#If we are in the info header. Second line of Annovar output
		elsif ($num_line == 1){
			#DONOTHING
		}
		
		###########################
		###########DATA############
		###########################
		else{
			#Split the fields of the Annovar annotation ouput
			my @ann_fields = split($tab,$line);
			my $tot_fields = $tot_cols + $tot_info_annovar;
			if ( scalar(@ann_fields) ne $tot_fields ){
				print_and_log("WARNING: Line $num_line contains ".scalar(@ann_fields)." fields while they should be $tot_fields..\n",$log_file);				
			}
			#Select only those variants which have exac frequency lower than a threshold
			#or that do not have Exac frequency
			if ( $ann_fields[$exacall_ind] ne '-' and  $ann_fields[$exacall_ind] ne '.' ){
				if ($ann_fields[$exacall_ind] > $exacall_min_freq ){
					$take_var = 0;
				}
			}
			#Takes or not the variant depending by the ExAC frequency
			if ( $take_var == 1 ){
				
				###########
				#CHR, ID, POS, REF, ALT, QUAL
				##########
				#To get the compid
				my $compid = "";
				my $varinfo_count = 1;
				#Here we get and print the info about the variant from the OtherInfo column of the Annovar output
				for (my $info_ind = ($tot_cols + $start_info_annovar)-1; $info_ind < (($tot_cols + $tot_info_annovar) -1); $info_ind++){
					if ($ann_fields[$info_ind] ne '-' and $ann_fields[$info_ind] ne '' ){
						$newline .= $ann_fields[$info_ind].$tab;
						
					}else{
						$newline .= "-".$tab;
					}
					#The fields (chrom_pos_ref_alt) are concatenated to set the compid The fourth is ID
					if ( $varinfo_count <= 5 and $varinfo_count != 3){
						$compid .= $ann_fields[$info_ind]."_";
					}
					$varinfo_count++;
				}
				#Remove last _ from compid
				chop($compid);
				
				######################
				##INFO and FORMAT fields using the all_var.info file
				#####################		
				#From the last field of the other info column of the Annovar annotation get the index of the all_var.info file
				#where it should fetch the genotype information
				my $info_format_field_ind = $ann_fields[($tot_cols + $tot_info_annovar)-1];
								
				#Split the fields of all_var.info  using the program separator
				my @format_fields = split(/\Q$mult_ann_sep\E/,$oth_info_hash->{$info_format_field_ind});
				
				################INFO FIELD
				#The first information from the file is the INFO field, get it and remove and print into the file
				my $info_field = shift @format_fields;
				#If the INFO field is not empty, add to the string
				#NB. different INFO fields are separated with semicolon!
				if ( $info_field ne $vargenius_empty_val){
					my @sel_info_fields = split($semicolon,$info_field);
					foreach my $sel_info_field (@sel_info_fields){
						$newline .=  $sel_info_field.$tab;	
					}
				}else{
					#Add empty val for each field needed
					foreach my $info_field_needed (@info_fields_needed){
						$newline .=  $vargenius_empty_val.$tab;	
					}
					
				}
				
				###############FORMAT FIELD
				#total number of FORMAT fields per sample
				my $tot_f_per_s = scalar(@format_fields) /scalar (@sort_samples);
				#Add ref_reads, alt_reads, perc_alt_reads, tot_reads
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$ref_r_pos].$tab.$sample_format_flds[$alt_r_pos].$tab.
						$sample_format_flds[$tot_r_pos].$tab.$sample_format_flds[$percalt_r_pos].$tab;
				}
				#Add GT and Zygosity information for each sample
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$gt_pos].$tab.$sample_format_flds[$zyg_pos].$tab;
				}
				#Add DP and GQ information for each sample
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$dp_pos].$tab.$sample_format_flds[$gq_pos].$tab;
				}				
										
				#Annotation
				#Get only those column from @ann_fields that have been selected and inserted in 
				#@annotations
				foreach my $fld_touse (@flds_touse){
					my $elem = $ann_fields[$fld_touse];
				
					
					##################################
					#PRINT THE TRANSCRIPT Genomic INFORMATION#
					##################################
					#NB:If you put additional evaluation use ELSIF!!
					if ( grep {/\b$fld_touse\b/} @compl_gen_ann_cols ){
						my $subj_field = $ann_fields[$fld_touse];
						#Transcript annotation field must not be empty or unknown
						#[E.g.: PLEKHN1:NM_001160184:exon13:c.1480A>C:p.R494R , PLEKHN1:NM_032129:exon14:c.1585A>C:p.R529R ]
						if ( ($ann_fields[$fld_touse] ne $annov_empty_val) and (lc($ann_fields[$fld_touse]) ne 'unknown') ){
							#Change quotes with nothing
							$ann_fields[$fld_touse] =~ s/[\"\']//g;
							#Separates the transcript info fields
							my @tr_fields = split($comma,$ann_fields[$fld_touse]);
							#Separates the fields of the first info 
							my @first_tr_info = split(":",shift @tr_fields);
							#after the shift the array has one element less
							# we put the values of the first element in different fields
							#if an element is not there will be substituted with $vargenius_empty_val
							#my $tr_ann_fields = $cfg_hash->{'tr_ann_fields'};#REMOVE
							my @tr_ann_fields = split($comma,$cfg_hash->{'tr_ann_fields'});
							for (my $fld = 0; $fld < scalar(@tr_ann_fields); $fld++){
								if ( $fld < scalar(@first_tr_info) ) {
									$newline .= $first_tr_info[$fld].$tab;
								}else{
									$newline .= $vargenius_empty_val.$tab;
								}
							}
							#and put the remaining transcript info in a single value if it is there 
							my $rem_tr_info = $vargenius_empty_val; 
							if (scalar(@tr_fields) > 0 ){
								$rem_tr_info = join($comma,@tr_fields);
							}
							$newline .= $rem_tr_info.$tab;
						}else{
							#the value is not there, hence
							#substitute the same number of fields with vargenius_empty_val (-)
							my @tr_fields = split($comma,$cfg_hash->{'tr_ann_fields'});
							foreach my $tr_f (@tr_fields){
								$newline .= $vargenius_empty_val.$tab;
							}
							#The dot for the base transcript info
							$newline .= $vargenius_empty_val.$tab;
						}
					}
					#NB:If you put additional evaluation use ELSIF!!				
					
					##########################
					######Any other field    #
					##########################
					else{
						#IF the value is there
						if ( $ann_fields[$fld_touse] ne $annov_empty_val) {
							#Remove single and double quotes
							$ann_fields[$fld_touse] =~ tr/"'//d;
							$newline .=  $ann_fields[$fld_touse].$tab;
						}else{
							#the value is not there, hence
							#substitute the field with vargenius_empty_val (-)
							$newline .=  $vargenius_empty_val.$tab;
						}
					}						
				}		
				#NB: If you put additional transcript information use ESLIF like here
				
				#NB: ADDITIONAL INFORMATION WILL BE PUT AT THE END OF THE FILE AND COULD BE REARRANGED LATER
				
				#################
				#ALLELE FREQUENCY
				#################
				#Getting and adding the variant internal frequency and the factors which generated the frequency
				#print_and_log("Getting the variant allele frequency and factors which generated it using the composite id: $compid\n ",$log_file);#DEBUGCODE
				if (defined $db_hash->{$compid}->{'allele_freq'} and defined $db_hash->{$compid}->{'freq_factors'} ) {
					$newline .=  $db_hash->{$compid}->{'allele_freq'}.$tab.$db_hash->{$compid}->{'freq_factors'}.$tab;
				}else{
					$newline .=  "-1".$tab.$cfg_hash->{'vargenius_empty_val'}.$tab;
				}

			}
		}
		#Put the new line
		if ($newline ne ""){
			#Remove the last tab
			chop($newline);	
			print FILT_OUT $newline."\n";
		}
		$num_line++;
	}
 	
	#Close files
	close(FILT_OUT);
	close(ANN_OUT);
}





=head2 output_from_annotation_without_db

 Title   : output_from_annotation_without_db
 Usage   : output_from_annotation_without_db( -database => 'name of the database,
                               );

 Function:  Uses the output from Annovar to build a first output for VarGenius.
						It separates the information for the transcripts and also adds the genotype information.
 						-gets the needed fields for variants,statistics and genotype tables
 						
						
						Then a final output in the form of a tab separated table is constructed:
						- an Header is built using both data from the configuration file
						- gene information are in the table GENES and will be fetched associated to genes symbols
						- genotype information is fetched related to variants from the GENOTYPE_SAMPLE table hence
								   fields names for this information is collected too
						
							
 Returns : 

=cut
sub output_from_annotation_without_dbOLD {
	my $cfg_hash = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $vargenius_out = shift;
	my $compl_annovar_out = shift;
	my $other_info_f = shift;
	my $vcf_file = shift;
	
	my $semicolon = ";";
	my $tab = "\t";
	my $comma = ",";
	
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	#Annotations fields are obtained from the program_config and are the fields from Annovar output that
	#will be used.
	my @annotations_fields = split($comma,$cfg_hash->{'annotations_fields'});
	my $ann_fields_used = 0;#Used to check if we are using all. Otherwise there is an error
	
	#The total number of information carried from Annovar in the Otherinfo column
	my $other_info_col = $cfg_hash->{'other_info_col'};
	my $tot_info_annovar = $cfg_hash->{'tot_info_annovar'};
	my $start_info_annovar = $cfg_hash->{'start_info_ann_wo_db'};
	my $annov_empty_val = $cfg_hash->{'annov_empty_val'};

	
	#Inizialization of Variants, info and format fields to pick from the output of annovar
	my @variants_fields_needed = ();
	if (defined $cfg_hash->{'variant_fields_from_annovar'}){ @variants_fields_needed = split($comma,$cfg_hash->{'variant_fields_from_annovar'});}	
	my @info_fields_needed = ();
	if (defined $cfg_hash->{'info_fields_from_annovar'}){ @info_fields_needed = split($comma,$cfg_hash->{'info_fields_from_annovar'});}
	my @format_fields_needed = ();
	if (defined $cfg_hash->{'format_fields_from_annovar'}){ @format_fields_needed =split($comma,$cfg_hash->{'format_fields_from_annovar'});}
	
	#Index in the list of the refgene_symbol, useful later when we have to access gene info
	my $gene_field_ind = -1;
	#print_and_log( "Printing the VarGenius outputs...\n",$log_file); #DEBUGCODE
	
	#Positions of the genotype information
	my $gt_pos = $cfg_hash->{'gt_pos'};
	my $zyg_pos = $cfg_hash->{'zyg_pos'};
	my $ref_r_pos = $cfg_hash->{'ref_r_pos'};
	my $alt_r_pos = $cfg_hash->{'alt_r_pos'};
	my $tot_r_pos = $cfg_hash->{'tot_r_pos'};
	my $percalt_r_pos = $cfg_hash->{'percalt_r_pos'};
	my $dp_pos = $cfg_hash->{'dp_pos'};
	my $gq_pos = $cfg_hash->{'gq_pos'};
	
	
	####SETTINGS FOR BED FILES IN ANNOVAR
	#In Annovar, if you use a BED file you must write the name with _BED into ##Annovar DATABASES in program_config.txt
	#then Annovar when annotates uses the column name "bed" for the first element in the annotation, then "bed2,bed3,..,bedN".
	#We change the header name accordingly here bulding an hash containing for each BED Annovar header  the corresponding database name
	#First get the XX_BED from the three parameters annov_g_dbs,annov_f_dbs, annov_r_dbs (THE SAME ORDER IS IMPORTANT)
	my @ann_all_dbs = ();
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_g_dbs'}));
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_f_dbs'}));
	push(@ann_all_dbs,split(",",$cfg_hash->{'annov_r_dbs'}));
	#The smartest way to fetch just those with _BED is the following
	my @ann_bed_dbs = grep(/_BED/, @ann_all_dbs);
	#Now remove the _BED and put it into an hash
	my $ann_bed_names;
	my $bed_counter = 1;
	foreach my $ann_bed_db (@ann_bed_dbs){
			$ann_bed_db =~ s/_BED//g;
			if ($bed_counter > 1){
				$ann_bed_names->{$cfg_hash->{'bed_ext'}.$bed_counter} = $ann_bed_db;
			}else{
				$ann_bed_names->{$cfg_hash->{'bed_ext'}} = $ann_bed_db;	
			}
			$bed_counter++;
	}
	
	##### FIRST PART HEADER#########
	###Here I first save the header for the variants, statistics and genotype info then
	# the header for the annotation will be added during the read of the annovar output

	#HEADER GENOTYPE INFOS
	#The genotype info will be print in the order given in the user config file (samples_order)

	#Getting the samples list associated to the analysisid
  my $samples_h;#Samples hash to be reordered
	my @sort_samples = ();#The array with samples id reordered
	#Get the samples ids involved for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res_group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	my @group_sam = ();
	#################VQSR Considerations
	#If VQSR was executed	then we get the samples from those selected by the user.
	#I check the genotype refinement step into the db
 	my $genotref = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_genotref'},
									$cfg_hash->{'db_analysis_id'},$analysis_id); 
	#print_and_log( "genotref = $genotref and runsamples: ".$cfg_hash->{'run_samples'}." \n",$log_file);									
	if ( ($genotref == 1) and ($cfg_hash->{'run_samples'} ne 'ALL') ){
		my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},",");
		@group_sam = split(",",$samples_ids);
	}else{
		foreach my $sample_id (keys %{$res_group_sam}){
			push(@group_sam,$sample_id);
		}
	}
	#################VQSR Considerations	
	
	#Print the header for information about the genotype for the sample
	my $sample_gen_fields1 = "";
	my $sample_gen_fields2 = "";
	my $kinship_ok = 1;
	#Get the kinship, to make the resorting
	foreach my $sample_id (@group_sam){
		#Obtain the kinship from the database given the sample id
		my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
														$cfg_hash->{'db_sample_id'},$sample_id);	
		#Obtain the sample name from the database given the sample id
		my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
														$cfg_hash->{'db_sample_id'},$sample_id);	
		#Build an hash to reorder
		#If at least one sample does not specifies kinship, set a flag
		if (defined $kinship ){
			if ( $kinship ne $cfg_hash->{'vargenius_empty_val'}){	
				$samples_h->{$sample_id}->{'id'} = $sample_name;
				$samples_h->{$sample_id}->{'k'} = $kinship;
			}else{
				$kinship_ok = 0;
			}
	  }else{
			$kinship_ok = 0;
		} 
	}
	
	#Sort the samples as in samples_order parameter, if the kinship is present
	if ( $kinship_ok == 1){
		#Sort the samples as in samples_order parameter
		sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	}#Otherwise Get the samples order as in the VCF
	else{
		my $vcf = Vcf->new(file=>$vcf_file);
		$vcf->parse_header(); 
		my (@samples_vcf) = $vcf->get_samples();
		#Get the VCF order
		@sort_samples = @samples_vcf;
	}

	
	#Initialize two strings for two sets of information for genotype
	my $gt_add_1 = "";
	my $gt_add_2 = "";
	#print_and_log("Output from annotation samples header: " ,$log_file);#DEBUGCODE
	#Write header for genotypes. Go through sample names
	foreach my $sample_name (@sort_samples){
		print_and_log("Getting genotype header for the sample : $sample_name..\n" ,$log_file);#DEBUGCODE
		#print_and_log(" $sample_name " ,$log_file);#DEBUGCODE
		
		$sample_gen_fields1 .= "$sample_name\_".$cfg_hash->{'ref_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'alt_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'tot_r_str'}.
														"$tab$sample_name\_".$cfg_hash->{'percalt_r_str'}.$tab;
		
		$gt_add_1 .= "$sample_name\_".$cfg_hash->{'gt_str'}.
								 "$tab$sample_name\_".$cfg_hash->{'zyg_str'}.$tab;
		$gt_add_2 .= "$sample_name\_".$cfg_hash->{'dp_str'}.
								 "$tab$sample_name\_".$cfg_hash->{'gq_str'}.$tab;				
						
	}
	chop($sample_gen_fields1);		
	$sample_gen_fields2 .= $gt_add_1.$gt_add_2;
	chop($sample_gen_fields2);		
	#END HEADER GENOTYPE INFOS
	
	###########################################
	## PRINT THE FIRST PART OF THE HEADER   ###
	###########################################
	#Print the header using the variants and statistics fields and the genotype information
	#(the second part of the header consists in the annotation of Annovar and is written during the
	#annovar output reading later)
	my $head_var_stat_gen = "";
	if ((scalar @variants_fields_needed) > 0){ $head_var_stat_gen .= join($tab,@variants_fields_needed).$tab;}
	if ((scalar @info_fields_needed) > 0){ $head_var_stat_gen .= join($tab,@info_fields_needed).$tab;}
	if ($sample_gen_fields1 ne ''){ $head_var_stat_gen .= $sample_gen_fields1.$tab; }
	if ($sample_gen_fields2 ne ''){ $head_var_stat_gen .= $sample_gen_fields2.$tab; }
	
	#Build an hash with all the INFO and FORMAT fields for each variant
	my $oth_info_hash;
	open (OTH_INFO,"<".$other_info_f) or die "ERROR: Cannot open $other_info_f. The program will exit..\n";
	while (my $line = <OTH_INFO>){
		chop($line);
		my @fields = split($tab,$line);
		#Put in the hash at the variant index location the corresponding info prepared in generate_vcf
		$oth_info_hash->{$fields[0]} = $fields[1];
	}
	close(OTH_INFO);
	######END FIRST PART HEADER SAVING#####
	
	open (FILT_OUT,">".$vargenius_out) or die "ERROR: Cannot open $vargenius_out. The program will exit..\n";
		
	#####################################
	###########PRINT OUTPUT##############
	#####################################	

	#Open output from the annotation
	open (ANN_OUT,"<".$compl_annovar_out) or die "ERROR: Cannot open $compl_annovar_out. The program will exit..\n";
	#RefSeq fields: contains the transcripts information.Convert it in a readable form
	my $tr_ann_field_col = -1;
	######################
	#Genomic annotation (RefSeq, Gencode, ..) inizialization
	#Get the field of Annovar annotation  describing the AAChange (the field that must be splitted)
	my $aa_change_gen_field = $cfg_hash->{'aa_change_gen_field'};
	$aa_change_gen_field =~ s/[\-\.\+]/\_/g;
	$aa_change_gen_field = lc($aa_change_gen_field);
	#Genomic annotation fields to use (RefSeq, Gencode, etc)
	my $genomic_ann_fields = $cfg_hash->{'genomic_ann_fields'};
	$genomic_ann_fields =~ s/[\-\.\+]/\_/g;
	$genomic_ann_fields = lc($genomic_ann_fields);	
	my @genomic_ann_fields = split(",",$genomic_ann_fields);
	my @compl_gen_ann_fields = ();#Concatenate AAChange with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_gen_ann_fields,$aa_change_gen_field."_".$genomic_ann_field);
	}
	my @compl_gen_ann_cols = ();#Used for the indexes of these fields
	##Get the field names that will be used for splitting the AAChange field
	#my $tr_ann_field = $cfg_hash->{'tr_ann_field'};
	#$tr_ann_field =~ s/[\-\.\+]/\_/g;
	#$tr_ann_field = lc($tr_ann_field);
	#################

	#The index of the column Otherinfo from Annovar (with varid)
	my $other_info_col_ind = -1;
	#Information required to filter base on the frequencies in ExAC
	my $exacall_field = $cfg_hash->{'exacall_fld'};
	my $exacall_ind = -1;
	my $exacall_min_freq = $cfg_hash->{'exacall_min_freq'};
	
	#Numerates the line just to know what is the header of the file
	my $num_line = 0;
	my $tot_cols = 0;
	#The array of field names
	my @head_flds = ();
	#the indexes of the fields that are to be used
	my @flds_touse = ();
	
	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	#################
	
	#Getting the type of sequencing. The frequencies table will be updated for that only
	my $seqtype = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);					
				
	#Go through each line of a TXT tab separated table
	while (my $line = <ANN_OUT>){ # and ($num_line < 1000){#DEBUGCODE
		#Remove last \n
		chop($line);
		#Initialize the newline
		my $newline = "";

		#Used to filter the variants
		my $take_var = 1;
		
		###########################
		##########HEADER###########
		###########################
		#If we are in the header
		if ($num_line == 0){
			#Put the first part of the header
			$newline .= $head_var_stat_gen;
			#Split the header to get the fields names
			@head_flds = split($tab,$line);
			$tot_cols = scalar(@head_flds);
			#Number of column that will be inserted in the flds_touse
			my $num_col = 0;
			#Keep in count if a field has been already inserted
			my @insd_flds = ();
			#Get each field in the fields variable
			foreach my $head_fld (@head_flds){
				
				#print_and_log( "Headfield:$head_fld \n",$log_file);
				#Other info is the last field with supplementary information
				#and contains the variant id
				if ($head_fld eq $other_info_col) {
					#print_and_log( "This is otherinfo:$head_fld \n",$log_file);
					$other_info_col_ind = $num_col;
					#push(@flds_touse,$num_col);
				}else{
						#Before to add the value I search for those strings
						#that cannot be inserted as they are written in the output from Annovar
						#
						#	print_and_log( "Field $head_fld\n",$log_file);
						#Convert all - and . into _
						$head_fld =~ s/[\-\.\+]/\_/g;
						#These are keywords in SQL, change as before
						if ( ($head_fld eq 'Start') or ($head_fld eq 'End') or ($head_fld eq 'Ref') ){
							$head_fld = "_".$head_fld;
						}
						#Convert upper to lowercase
						$head_fld = lc( $head_fld);
						#If the field starts with a number (e.g. 1000g), add 'db'
						if ($head_fld =~ /^\d/){ $head_fld= "db".$head_fld;}					
						
						#Here I filter the annotation using the frequency of Exac. 
						#This was useful for the Undiagnose Diseases. Being undiagnosed, the variants
						#are interesting only if they are not frequent in Exac
						if ( $head_fld eq $exacall_field){
							$exacall_ind = $num_col; 
							print "Exacall $head_fld ind :$exacall_ind\n"
						}
						
						#If this field is needed then it will be inserted
						if ( grep {/\b$head_fld\b/} @annotations_fields){
							#SECOND PART OF THE HEADER IS WRITTEN HERE
							$ann_fields_used++;#To check if we use all
							#check if this field was already inserted
							if ( ! grep {/\b$head_fld\b/} @insd_flds) {
								#Columns for which the name must be substituted
								#If the column is one among bed..bedN
								if ( $head_fld =~ /^bed\d*$/ ){
									$newline .= $ann_bed_names->{$head_fld}.$tab;

								}#NB:Add elsif if there are others..
								#Columns for which something must be added or they must be kept as they are
								else{
									#Keep the index of the field that keeps the RefSeq information
									#and adds the fields names for the transcript infos
									if ( grep {/\b$head_fld\b/} @compl_gen_ann_fields ){
										#Get the column index
										push(@compl_gen_ann_cols,$num_col);
										#Get the genomic annotation name removing the aachange field
										$head_fld =~ s/^$aa_change_gen_field\_// ;
										#Get the header by summing the head field obtained with all the annotation fields
										my $tr_ann_fields = $cfg_hash->{'tr_ann_fields'};
										my @tr_ann_fields = split (",",$tr_ann_fields);
										foreach my $tr_ann_field (@tr_ann_fields){
											$newline .= $head_fld."_".$tr_ann_field.$tab;
										}							
									}									
																	
									#This will be applied for any other field!
									$newline .= $head_fld.$tab;
								}
								push(@insd_flds,$head_fld);
								#Save the indexes to use
								push(@flds_touse,$num_col);	
							}else{
								print_and_log("WARNING: File $compl_annovar_out contains a repeated field: $head_fld.".
												" The column number $num_col will not be used..\n",$log_file);
							}
						}
				}
				
				$num_col++;
			}

			##Check if all the annotations field are used, otherwise stop with error
			#unless ( scalar( @annotations_fields) == $ann_fields_used){
				#log_and_exit("ERROR: File $ann_fields_used annotation fields used while they are: ".scalar( @annotations_fields).".".
				#" Check variable annotations_fields into the configuration file. Exiting...\n",$log_file);				
			#}
			#NEW INFORMATION HEADER ( N.B. the last tab will be removed and a newline inserted at the end!)
			#Add the internal frequency field 
			$newline .= $cfg_hash->{'intern_var_freq_fld'}.$tab.$cfg_hash->{'var_frequency_factors_fld'}.$tab;
			
		}##########HEADER###########
		#If we are in the info header. Second line of Annovar output
		elsif ($num_line == 1){
			#DONOTHING
		}
		
		###########################
		###########DATA############
		###########################
		else{
			#Split the fields of the Annovar annotation ouput
			my @ann_fields = split($tab,$line);
			my $tot_fields = $tot_cols + $tot_info_annovar;
			if ( scalar(@ann_fields) ne $tot_fields ){
				print_and_log("WARNING: Line $num_line contains ".scalar(@ann_fields)." fields while they should be $tot_fields..\n",$log_file);				
			}
			#Select only those variants which have exac frequency lower than a threshold
			#or that do not have Exac frequency
			if ( $ann_fields[$exacall_ind] ne '-' and  $ann_fields[$exacall_ind] ne '.' ){
				if ($ann_fields[$exacall_ind] > $exacall_min_freq ){
					$take_var = 0;
				}
			}
			#Takes or not the variant depending by the ExAC frequency
			if ( $take_var == 1 ){
				
				###########
				#CHR, ID, POS, REF, ALT, QUAL
				##########
				#To get the compid
				my $compid = "";
				my $varinfo_count = 1;
				#Here we get and print the info about the variant from the OtherInfo column of the Annovar output
				for (my $info_ind = ($tot_cols + $start_info_annovar)-1; $info_ind < (($tot_cols + $tot_info_annovar) -1); $info_ind++){
					if ($ann_fields[$info_ind] ne '-' and $ann_fields[$info_ind] ne '' ){
						$newline .= $ann_fields[$info_ind].$tab;
						
					}else{
						$newline .= "-".$tab;
					}
					#The fields (chrom_pos_ref_alt) are concatenated to set the compid The fourth is ID
					if ( $varinfo_count <= 5 and $varinfo_count != 3){
						$compid .= $ann_fields[$info_ind]."_";
					}
					$varinfo_count++;
				}
				#Remove last _ from compid
				chop($compid);
				
				######################
				##INFO and FORMAT fields using the all_var.info file
				#####################		
				#From the last field of the other info column of the Annovar annotation get the index of the all_var.info file
				#where it should fetch the genotype information
				my $info_format_field_ind = $ann_fields[($tot_cols + $tot_info_annovar)-1];
								
				#Split the fields of all_var.info  using the program separator
				my @format_fields = split(/\Q$mult_ann_sep\E/,$oth_info_hash->{$info_format_field_ind});
				
				################INFO FIELD
				#The first information from the file is the INFO field, get it and remove and print into the file
				my $info_field = shift @format_fields;
				#If the INFO field is not empty, add to the string
				#NB. different INFO fields are separated with semicolon!
				if ( $info_field ne $vargenius_empty_val){
					my @sel_info_fields = split($semicolon,$info_field);
					foreach my $sel_info_field (@sel_info_fields){
						$newline .=  $sel_info_field.$tab;	
					}
				}else{
					#Add empty val for each field needed
					foreach my $info_field_needed (@info_fields_needed){
						$newline .=  $vargenius_empty_val.$tab;	
					}
					
				}
				
				###############FORMAT FIELD
				#total number of FORMAT fields per sample
				my $tot_f_per_s = scalar(@format_fields) /scalar (@sort_samples);
				#Add ref_reads, alt_reads, perc_alt_reads, tot_reads
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$ref_r_pos].$tab.$sample_format_flds[$alt_r_pos].$tab.
						$sample_format_flds[$tot_r_pos].$tab.$sample_format_flds[$percalt_r_pos].$tab;
				}
				#Add GT and Zygosity information for each sample
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$gt_pos].$tab.$sample_format_flds[$zyg_pos].$tab;
				}
				#Add DP and GQ information for each sample
				for ( my $sample_ind = 0; $sample_ind < scalar (@sort_samples); $sample_ind++){
					#Get the FORMAT for the sample
					my $sample_format = $format_fields[$sample_ind];
					#Get all the FORMAT fields
					my @sample_format_flds =  split($semicolon,$sample_format);
					$newline .=  $sample_format_flds[$dp_pos].$tab.$sample_format_flds[$gq_pos].$tab;
				}				
										
				#Annotation
				#Get only those column from @ann_fields that have been selected and inserted in 
				#@annotations
				foreach my $fld_touse (@flds_touse){
					my $elem = $ann_fields[$fld_touse];
				
					
					##################################
					#PRINT THE TRANSCRIPT Genomic INFORMATION#
					##################################
					#NB:If you put additional evaluation use ELSIF!!
					if ( grep {/\b$fld_touse\b/} @compl_gen_ann_cols ){
						my $subj_field = $ann_fields[$fld_touse];
						#Transcript annotation field must not be empty or unknown
						#[E.g.: PLEKHN1:NM_001160184:exon13:c.1480A>C:p.R494R , PLEKHN1:NM_032129:exon14:c.1585A>C:p.R529R ]
						if ( ($ann_fields[$fld_touse] ne $annov_empty_val) and (lc($ann_fields[$fld_touse]) ne 'unknown') ){
							#Change quotes with nothing
							$ann_fields[$fld_touse] =~ s/[\"\']//g;
							#Separates the transcript info fields
							my @tr_fields = split($comma,$ann_fields[$fld_touse]);
							#Separates the fields of the first info 
							my @first_tr_info = split(":",shift @tr_fields);
							#after the shift the array has one element less
							# we put the values of the first element in different fields
							#if an element is not there will be substituted with $vargenius_empty_val
							#my $tr_ann_fields = $cfg_hash->{'tr_ann_fields'};#REMOVE
							my @tr_ann_fields = split($comma,$cfg_hash->{'tr_ann_fields'});
							for (my $fld = 0; $fld < scalar(@tr_ann_fields); $fld++){
								if ( $fld < scalar(@first_tr_info) ) {
									$newline .= $first_tr_info[$fld].$tab;
								}else{
									$newline .= $vargenius_empty_val.$tab;
								}
							}
							#and put the remaining transcript info in a single value if it is there 
							my $rem_tr_info = $vargenius_empty_val; 
							if (scalar(@tr_fields) > 0 ){
								$rem_tr_info = join($comma,@tr_fields);
							}
							$newline .= $rem_tr_info.$tab;
						}else{
							#the value is not there, hence
							#substitute the same number of fields with vargenius_empty_val (-)
							my @tr_fields = split($comma,$cfg_hash->{'tr_ann_fields'});
							foreach my $tr_f (@tr_fields){
								$newline .= $vargenius_empty_val.$tab;
							}
							#The dot for the base transcript info
							$newline .= $vargenius_empty_val.$tab;
						}
					}
					#NB:If you put additional evaluation use ELSIF!!				
					
					##########################
					######Any other field    #
					##########################
					else{
						#IF the value is there
						if ( $ann_fields[$fld_touse] ne $annov_empty_val) {
							#Remove single and double quotes
							$ann_fields[$fld_touse] =~ tr/"'//d;
							$newline .=  $ann_fields[$fld_touse].$tab;
						}else{
							#the value is not there, hence
							#substitute the field with vargenius_empty_val (-)
							$newline .=  $vargenius_empty_val.$tab;
						}
					}						
				}		
				#NB: If you put additional transcript information use ESLIF like here
				
				#NB: ADDITIONAL INFORMATION WILL BE PUT AT THE END OF THE FILE AND COULD BE REARRANGED LATER
				
				#################
				#ALLELE FREQUENCY
				#################
				#Getting and adding the variant internal frequency and the factors which generated the frequency
				#print_and_log("Getting the variant allele frequency and factors which generated it using the composite id: $compid\n ",$log_file);#DEBUGCODE
	
				my $var_query = "SELECT ".$cfg_hash->{"db_$seqtype\_allele_freq"}.",".$cfg_hash->{"db_$seqtype\_freq_factors"}.
												" FROM ".$cfg_hash->{'db_variants_table'}." WHERE ".$cfg_hash->{'db_var_compid'}." = '$compid';";	
				#print_and_log( "Executing: $var_query\n",$log_file);		
				my @res = ();
				do_fetch_row_array_woconn($dbh,$var_query,\@res);
				if ( scalar(@res) > 0 ){
					if ( (defined $res[0]) and (defined $res[1]) ){
						if ($res[0] ne '' and $res[1] ne ''){
							#print_and_log("allele_freq: for $compid are ".$res[0]." ".$res[1]." \n ",$log_file);	
							$newline .=  $res[0].$tab.$res[1].$tab;
						}else{
							$newline .=  $cfg_hash->{'vargenius_empty_val'}.$tab.$cfg_hash->{'vargenius_empty_val'}.$tab;
						}	#code
					}else{
						$newline .=  $cfg_hash->{'vargenius_empty_val'}.$tab.$cfg_hash->{'vargenius_empty_val'}.$tab;
					}
				}else{
					$newline .=  $cfg_hash->{'vargenius_empty_val'}.$tab.$cfg_hash->{'vargenius_empty_val'}.$tab;
				}
			}
		}
		#Put the new line
		if ($newline ne ""){
			#Remove the last tab
			chop($newline);	
			print FILT_OUT $newline."\n";
		}
		$num_line++;
	}
	#Disconnect db
  $dbh->disconnect(); 	
	#Close files
	close(FILT_OUT);
	close(ANN_OUT);
}





=head2 separate_analysis_from_joint

 Title  : separate_analysis_from_joint
 Usage  : separate_analysis_from_joint( - filePath => 'the path to the file'
															- func_column => the column with the type of variant
															- detail_col => the column with the detail of the variant);
 
 Function: This script takes in input the output for multiple joined samples and
					creates  a specific output for a set of samples. 
					It basically goes through the header of the output and:
					- gets the starting column numbers of: reads info, genotype info and quality info
					- at the quality info puts the sample name into an array @samplenames
					- prints the variants info and annotation as it is and only the header for those samples
					  where grep is succesfull for $cfg_hash->{'run_samples'}
					  
					  
					Then goes through the data:
					- prints variants info and annotation as it is and...
					- loops through all sample names in @samplenames and when encounters
						one of those from the selected uses the column numbers for read, genotype and info qual
						to get its specific information from the row
						
 Returns : an output containing genotype information for the specified samples

=cut		
sub separate_analysis_from_joint {
	my $cfg_hash = shift;
	my $needed_samples = shift;
	my $analysis_id = shift;
	my $complete_out = shift;	
	my $out_path = shift;
	my $log_file = shift;
	
	my $sep =",";
	my $samples_ids = separate_input_ids($needed_samples,",");
	my @samples_ids = split($sep,$samples_ids);
	#Array with sample names
	my @sample_names = ();
	#Array with sample names to keep
	my @sample_names_tokeep = ();	
	foreach my $sel_sample_id (@samples_ids){
			my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
							$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'},"$analysis_id,$sel_sample_id"); 
		push(@sample_names_tokeep,$sample_name);
		#print_and_log( " $sample_name...",$log_file);#DEBUGCODE
	}
	#print_and_log( " \n",$log_file);#DEBUGCODE
	#print_and_log( "Opening $complete_out and getting the output for ".$cfg_hash->{'run_samples'}."..\n",$log_file);#DEBUGCODE
	open (C_OUT, "<$complete_out") or die ("Cannot open $complete_out\n");
	open (I_OUT, ">$out_path") 	or die ("Cannot open $out_path\n");
							
	#Get the string that indicates when the three blocks of information  start
	my $r_start = "_".$cfg_hash->{'ref_r_str'};
	my $gt_start = "_".$cfg_hash->{'gt_str'};
	my $qual_start = "_".$cfg_hash->{'dp_str'};
	
	my $tab = "\t";
	#Number of infos for reads, genotype and quality
	my $reads_infos = 4;
	my $gt_infos = 2;
	my $qual_infos = 2;
	#Index where the three sections start
	my $rinfo_start = 0;
	my $gtinfo_start = 0;
	my $qualinfo_start = 0;
			
	#Numerates the line just to know what is the header of the file
	my $num_line = 0;			
	
	#Go through the file and create a new output for each different sample
	while ( my $row = <C_OUT> ) {
		#Remove newline
		chop($row);
		##########HEADER###########
		#If we are in the header
		if ($num_line == 0){
			#print "Reading header..\n";
			#Split the header to get the fields names
			my @head_flds = split($tab,$row);
			my $num_col = 0;
			my $num_file = 1;
			
			#Flags to get the index of the first element of each block
			my $rfound = 0;
			my $gtfound = 0;
			my $qfound = 0;
			
			#Get each field in the fields variable
			foreach my $head_fld (@head_flds){
				#Find the index where the reads info starts				
				if ($head_fld =~ /$r_start$/) {
					#print "Matched $head_fld with $r_start\n";
					if ($rfound == 0){
						$rinfo_start = $num_col;
						$rfound = 1;
					}
				}
				#Find the index where the genotype info starts
				if ($head_fld =~ /$gt_start$/) {
					#print "Matched $head_fld with $gt_start\n";
					if ($gtfound == 0){
						$gtinfo_start = $num_col;
						$gtfound = 1;
					}
				}
				#Find the index where the quality info starts
				if ($head_fld =~ /$qual_start$/) {
					#print "Matched $head_fld with $qual_start\n";
					if ($qfound == 0){
						$qualinfo_start = $num_col;
						$qfound = 1;
					}
					#When also the last information is found then save the sample name...
					my $sample_name = $head_fld;
					$sample_name =~ s/$qual_start//;
					#Fill an array with all sample names contained
					push(@sample_names,$sample_name);

					#Increment number of samples
					$num_file++;
					#print "sample_name: $sample_name r: $rinfo_start; gt:$gtinfo_start; q:$qualinfo_start\n";
				}
				#increment column 
				$num_col++;
			}
				#Print headers
				#Get header until the sample info starts
				my $first_part  = "";
				for (my $i = 0; $i < $rinfo_start; $i++){
						$first_part .= $head_flds[$i].$tab;
				}
				chop($first_part);	
				
				#Get header after the sample info end
				my $last_part  = "";
				my $last_part_st = $qualinfo_start + ($qual_infos * scalar(@sample_names));
				for (my $i = $last_part_st; $i < scalar(@head_flds); $i++){
						$last_part .= $head_flds[$i].$tab;
				}	
				chop($last_part);	
				#print header using the first and last part
				#and building the sample genotype header using the same method used
				#for the print of the output (output_from_annotation_without_db	)
				my $sample_gen_fields = "";
				my $gt_add_1 = "";
				my $gt_add_2 = "";
				
				foreach my $sample_name (@sample_names_tokeep){
					$sample_gen_fields .= "$sample_name\_".$cfg_hash->{'ref_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'alt_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'tot_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'percalt_r_str'}.
													 $tab;
					
					$gt_add_1 .= "$sample_name\_".$cfg_hash->{'gt_str'}.
											 "$tab$sample_name\_".$cfg_hash->{'zyg_str'}.
												$tab;
					$gt_add_2 .= "$sample_name\_".$cfg_hash->{'dp_str'}.
												"$tab$sample_name\_".$cfg_hash->{'gq_str'}.
												$tab;
					
					#print_and_log( "Creating the header for $sample_name ...\n",$log_file);
				}

				print I_OUT $first_part.$tab.$sample_gen_fields.$gt_add_1.$gt_add_2.$last_part."\n";
		##################
		#PRINT OF THE DATA
		##################
		}else{
			#Split the row to get the fields
			my @flds = split($tab,$row);
			#print "Reading $row..\n";
			#################################
			#Print the first part Variant info
			my $first_part  = "";
			for (my $i = 0; $i < $rinfo_start; $i++){
					$first_part .= $flds[$i].$tab;
			}

			print I_OUT $first_part;

			#print "$first_part\n";
			my $num_sample = 0;
			
			###############################
			#Print the second part: Genotype
			#print_and_log( "Sample names are: ".scalar(@sample_names)."...\n",$log_file);#DEBUGCODE
			my $second_part  = "";
			my $sample_gen_fields = "";
			my $gt_add_1 = "";
			my $gt_add_2 = "";			
			foreach my $sample_name (@sample_names){
				#print_and_log( "Going through info for sample $sample_name...\n",$log_file);#DEBUGCODE
				if (grep {/\b$sample_name\b/} @sample_names_tokeep){

					#print_and_log( "Printing info for sample $sample_name...\n",$log_file);#DEBUGCODE
					#Reads info
					my $rinfo_s =  $rinfo_start + ($reads_infos * $num_sample);
					#print_and_log( "Reads info starting at $rinfo_s. End: ".($rinfo_s + $reads_infos - 1)."\n",$log_file);#DEBUGCODE
					for (my $i = $rinfo_s; $i < ($rinfo_s + $reads_infos); $i++){
						$sample_gen_fields .= $flds[$i].$tab;
					}
					#Genotype info
					my $ginfo_s = $gtinfo_start + ($gt_infos * $num_sample);
					#print_and_log( "Genotype info starting at $ginfo_s. End: ".($ginfo_s + $gt_infos - 1)."\n",$log_file);#DEBUGCODE
					for (my $i = $ginfo_s; $i < ($ginfo_s + $gt_infos); $i++){
						$gt_add_1 .= $flds[$i].$tab;
					}
					#Quality info
					my $qinfo_s = $qualinfo_start + ($qual_infos * $num_sample);
					#print_and_log( "Quality info starting at $qinfo_s. End: ".($qinfo_s + $qual_infos - 1)."\n",$log_file);#DEBUGCODE
					for (my $i = $qinfo_s; $i < ($qinfo_s + $qual_infos); $i++){
						$gt_add_2 .= $flds[$i].$tab;
					}
				}
				#print_and_log( "2nd part:$second_part\n",$log_file);#DEBUGCODE
				$num_sample++;
			}		
			$second_part = 	$sample_gen_fields.$gt_add_1.$gt_add_2;
			print I_OUT $second_part;
						
			########################
			#Print the last part: Annotation
			my $last_part  = "";
			my $last_part_st = $qualinfo_start + ($qual_infos * scalar(@sample_names));
			for (my $i = $last_part_st; $i < scalar(@flds); $i++){
					$last_part .= $flds[$i].$tab;
			}
			chop($last_part);
			$last_part .= "\n";
			print I_OUT $last_part;
			#print "$last_part\n";
		}
		$num_line++;
	}
	
	close(C_OUT);
}


				

				
=head2 tag_panel_genes

 Title  : tag_panel_genes
 Usage  : tag_panel_genes ();
 
 Function: Given the user id, searches the data folder for the userid_data folder 
					and starts an R script that adds a flag to each row if the gene is one
					of the panel list in input or one of the candidates.
					Then shifts the flag near the genes information and converts the file in
					XLS format. Finally adds the path to the outlist genes
  
 Returns : the final output

=cut		
sub tag_panel_genes {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $vargenius_out = shift;
	my $panel = shift;
	my $user_id = shift;
	my $outFolder = shift;
	my $log_file = shift;
	
	my $panel_path = $cfg_hash->{'main_data_folder'}."/$user_id\_data/$panel";
	my $cand_path = $cfg_hash->{'main_data_folder'}."/$user_id\_data/$panel\_cand";

	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
															$cfg_hash->{'db_analysis_id'},$analysis_id);
																
	#print_and_log("Tagging panel genes..\n",$log_file);
	my $RLogPath = $outFolder."/".$cfg_hash->{'log_fold'}."/$analysis_name\_tag_panel_genes".$cfg_hash->{'R_log_file'};
	my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_utils'};
	my $out_path = extract_name($vargenius_out,"noext")."_$panel.".$cfg_hash->{'txt_ext'};
	my $empty_val = $cfg_hash->{'vargenius_empty_val'};
	my $gene_field =  $cfg_hash->{'gene_field'};
	my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args TAG_PANEL_GENES $vargenius_out $panel $panel_path $cand_path $empty_val $gene_field $out_path'  $R_script $RLogPath";
	
	print_and_log("The command is: $command\n",$log_file);
	try_exec_command($command) or R_die($command,$R_script);			
			
	return $out_path;
}


=head2 separate_joint_outputs

 Title  : separate_joint_outputs
 Usage  : separate_joint_outputs( - filePath => 'the path to the file'
															- func_column => the column with the type of variant
															- detail_col => the column with the detail of the variant);
 
 Function: This script takes in input the output for multiple joined samples and
					creates  a different output for each of them. This is a very simple script
					that could create many different files occupying more space but the outputs
					will be more readable
  
 Returns : a distinct output for each sample

=cut		
sub separate_joint_outputs {
	my $cfg_hash = shift;
	my $complete_out = shift;	
	my ($out_paths) = shift;
	my $log_file = shift;
	
	print "Opening $complete_out\n";
	
	open (C_OUT, "<$complete_out") or die ("Cannot open $complete_out\n");
	
	#Get the string that indicates when the three blocks of information  start
	my $r_start = "_".$cfg_hash->{'ref_r_str'};
	my $gt_start = "_".$cfg_hash->{'gt_str'};
	my $qual_start = "_".$cfg_hash->{'dp_str'};
	
	my $tab = "\t";
	#Number of infos for reads, genotype and quality
	my $reads_infos = 4;
	my $gt_infos = 2;
	my $qual_infos = 2;
	#Index where the three sections start
	my $rinfo_start = 0;
	my $gtinfo_start = 0;
	my $qualinfo_start = 0;
		

	#Array with sample names
	my @sample_names = ();
	
	#Numerates the line just to know what is the header of the file
	my $num_line = 0;			
	
	#Go through the file and create a new output for each different sample
	while ( my $row = <C_OUT> ) {
		#Remove newline
		chop($row);
		##########HEADER###########
		#If we are in the header
		if ($num_line == 0){
			print "Reading header..\n";
			#Split the header to get the fields names
			my @head_flds = split($tab,$row);
			my $num_col = 0;
			my $num_file = 1;
			
			#Flags to get the index of the first element of each block
			my $rfound = 0;
			my $gtfound = 0;
			my $qfound = 0;
			
			#Get each field in the fields variable
			foreach my $head_fld (@head_flds){
				#Find the index where the reads info starts				
				if ($head_fld =~ /$r_start$/) {
					#print "Mathced $head_fld with $r_start\n";
					if ($rfound == 0){
						$rinfo_start = $num_col;
						$rfound = 1;
					}
				}
				#Find the index where the genotype info starts
				if ($head_fld =~ /$gt_start$/) {
					#print "MAthced $head_fld with $gt_start\n";
					if ($gtfound == 0){
						$gtinfo_start = $num_col;
						$gtfound = 1;
					}
				}
				#Find the index where the quality info starts
				if ($head_fld =~ /$qual_start$/) {
					#print "MAthced $head_fld with $qual_start\n";
					if ($qfound == 0){
						$qualinfo_start = $num_col;
						$qfound = 1;
					}
					#When also the last information is found then save the sample name...
					my $sample_name = $head_fld;
					$sample_name =~ s/$qual_start//;
					push(@sample_names,$sample_name);
					#...and the paths of the new files
					push(@$out_paths,extract_name($complete_out,'noext')."_".$sample_name.".".$cfg_hash->{'txt_ext'});
					#Increment number of samples
					$num_file++;
					#print "sample_name: $sample_name r: $rinfo_start; gt:$gtinfo_start; q:$qualinfo_start\n";
				}
				#Print headers
				#Get header until the sample info starts
				my $first_part  = "";
				for (my $i = 0; $i < $rinfo_start; $i++){
						$first_part .= $head_flds[$i].$tab;
				}
				chop($first_part);	
				#Get header after the sample info end
				my $last_part  = "";
				my $last_part_st = $qualinfo_start + ($qual_infos * scalar(@sample_names));
				for (my $i = $last_part_st; $i < scalar(@head_flds); $i++){
						$last_part .= $head_flds[$i].$tab;
				}	
				chop($last_part);	
				#print header using the first and last part
				#and building the sample genotype header using the same method used
				#for the print of the output (output_from_annotation_without_db	)
				foreach my $sample_name (@sample_names){
					my $sample_gen_fields .= "$sample_name\_".$cfg_hash->{'ref_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'alt_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'tot_r_str'}.
													"$tab$sample_name\_".$cfg_hash->{'percalt_r_str'};
					
					my $gt_add_1 .= "$sample_name\_".$cfg_hash->{'gt_str'}.
											 "$tab$sample_name\_".$cfg_hash->{'zyg_str'};
					my $gt_add_2 .= "$sample_name\_".$cfg_hash->{'dp_str'}.
							 "$tab$sample_name\_".$cfg_hash->{'gq_str'};
							 		
					open (I_OUT, ">".extract_name($complete_out,'noext')."_".$sample_name.".".$cfg_hash->{'txt_ext'}) 
						or die ("Cannot open ".extract_name($complete_out,'noext')."_".$sample_name.".".$cfg_hash->{'txt_ext'}."\n");
					print I_OUT $first_part.$tab.$sample_gen_fields.$tab.$gt_add_1.$tab.$gt_add_2.$tab.$last_part."\n";
					close(I_OUT);
				}
				 
				#increment column 
				$num_col++;
			}
		##################
		#PRINT OF THE DATA
		##################
		}else{
			#Split the row to get the fields
			my @flds = split($tab,$row);
			#print "Reading $row..\n";
			#################################
			#Print the first part Variant info
			my $first_part  = "";
			for (my $i = 0; $i < $rinfo_start; $i++){
					$first_part .= $flds[$i].$tab;
			}
			foreach my $out_path (@$out_paths){
							open (I_OUT, ">>$out_path") or die ("Cannot open $out_path\n");
							print I_OUT $first_part;
							close(I_OUT);
			}
			#print "$first_part\n";
			my $num_sample = 0;
			
			###############################
			#Print the second part Genotype
			foreach my $sample_name (@sample_names){
				open (S_OUT, ">>".extract_name($complete_out,'noext')."_".$sample_name.".".$cfg_hash->{'txt_ext'})
				 or die ("Cannot open ".extract_name($complete_out,'noext')."_".$sample_name.".".$cfg_hash->{'txt_ext'}."\n");
				my $second_part  = "";
				print_and_log( "Printing info for sample $sample_name...\n",$log_file);
				#Reads info
				my $rinfo_s =  $rinfo_start + ($reads_infos * $num_sample);
				#print_and_log( "Reads info starting at $rinfo_s. End: ".($rinfo_s + $reads_infos - 1)."\n",$log_file);#DEBUGCODE
				for (my $i = $rinfo_s; $i < ($rinfo_s + $reads_infos); $i++){
					$second_part .= $flds[$i].$tab;
				}
				#Genotype info
				my $ginfo_s = $gtinfo_start + ($gt_infos * $num_sample);
				#print_and_log( "Genotype info starting at $ginfo_s. End: ".($ginfo_s + $gt_infos - 1)."\n",$log_file);#DEBUGCODE
				for (my $i = $ginfo_s; $i < ($ginfo_s + $gt_infos); $i++){
					$second_part .= $flds[$i].$tab;
				}
				#Quality info
				my $qinfo_s = $qualinfo_start + ($qual_infos * $num_sample);
				#print_and_log( "Quality info starting at $qinfo_s. End: ".($qinfo_s + $qual_infos - 1)."\n",$log_file);#DEBUGCODE
				for (my $i = $qinfo_s; $i < ($qinfo_s + $qual_infos); $i++){
					$second_part .= $flds[$i].$tab;
				}
				print S_OUT $second_part;
				close(S_OUT);
				#print "$second_part\n";
				$num_sample++;
			}			
			
			########################
			#Print the last part
			my $last_part  = "";
			my $last_part_st = $qualinfo_start + ($qual_infos * scalar(@sample_names));
			for (my $i = $last_part_st; $i < scalar(@flds); $i++){
					$last_part .= $flds[$i].$tab;
			}
			chop($last_part);
			$last_part .= "\n";
			foreach my $out_path (@$out_paths){
							open (I_OUT, ">>$out_path") or die ("Cannot open $out_path\n");
							print I_OUT $last_part;
							close(I_OUT);
			}
			#print "$last_part\n";
		}
		$num_line++;
	}
	
	close(C_OUT);
}




=head2 get_genotype_info

 Title   : get_genotype_info
 Usage   : get_genotype_info( -database => 'name of the database,
                               );

 Function: Gets information from the database related with the genotype for a given
						sample and a given variant
 
 					For multiallelic variant  (i.e. there are multiple variant 
 					ids in varids) it calls the function get_gt_from_pl
					to infer the biallelic genotype
					
					Puts a tag to describe the zygosity
					If alleles are different the sample is heterozygous	
					If the alleles are the same the sample is homozygous
 Returns : a string with 
						if a value is missing will add $vargenius_empty_val

=cut
sub get_genotype_info {
	my $cfg_hash = shift;
	my $log_file = shift;
	my $varid = shift;
	my $sampleid = shift;
	
	my $gt_ind  = $cfg_hash->{'gt_ind'};
	my $ad_ind  = $cfg_hash->{'ad_ind'};
	my $dp_ind  = $cfg_hash->{'dp_ind'};
	my $gq_ind  = $cfg_hash->{'gq_ind'};
	my $sep = "\t";
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	my $out_str = "";#Final output string
	my $out_str1 = "";
	my $out_str2 = "";
	my $out_str3 = "";
	#Gives the position of the alternative variant in the varids field
	#of genotype_sample table where the reference is always the first. Hence starts from 1
	my $var_pos = 1;
	my $genot_table = $cfg_hash->{'db_genotype_sample_table'};
	#Get the selected genotype fields from the configuration file
	my $sel_f_str = $cfg_hash->{'db_var_ids'}.",".$cfg_hash->{'genotype_fields'};
	
	if ( $sampleid ne $vargenius_empty_val and $sampleid ne '.'){
			#print_and_log("Varid: $varid sampleid: $sampleid -> ",$log_file);
		#select * from genotype_sample where varids like any (values('238,%'),('238'),('%,238,%'),('%,238'));
		my $var_query = "SELECT $sel_f_str FROM $genot_table WHERE ".$cfg_hash->{'db_var_ids'}." LIKE ".
			"ANY ( values('$varid'),('$varid,%'),('%,$varid'),('%,$varid,%') ) AND sampleid='$sampleid' ;";	
		#print_and_log( "Executing: $var_query\n",$log_file);		
		my @res = ();
		do_fetch_row_array($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$var_query,\@res);
		if (scalar (@res) > 0){
			#Check what is the position of this variant among the varids
			#and set the var_pos
			my @varids = split(",",$res[0]);
			foreach my $id ( @varids){
					if ( $id != $varid){
							$var_pos++;
					}else{
						last;
					}
			}
			my $fld_num = 1;
			foreach my $elem (@res){
				###################AD##########################
				#If the field is AD insert the first element, the element corrisponding
				#to the alternative allele and the last element
				# number of reference reads,number of alterate reads, total reads number, percentage alterate reads
				if ( $fld_num == $ad_ind   ){	
					#print_and_log("AD is: $elem - varpos = $var_pos\t",$log_file);
					if (defined $elem  and ($elem ne $vargenius_empty_val) and ($elem ne '.')){	
						my @all_dpts = split(",",$elem);
						#Compute the percentage of alterate reads respect to the sum
						#of the number of altered reads and non altered
						my $perc_alterate = 0;
						my $readsum = $all_dpts[0] + $all_dpts[$var_pos];
						if (  ($readsum ne 0) ){#($all_dpts[$var_pos] ne 0) and
							$perc_alterate = $all_dpts[$var_pos]/$readsum;
						}else{
							$perc_alterate = $vargenius_empty_val;
						}
						$out_str1 .=  $all_dpts[0].$sep.$all_dpts[$var_pos].$sep.$readsum.$sep.$perc_alterate.$sep;
					}else{
						$out_str1 .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
					}
				}
				###################GT##########################
				#If the field is GT than convert any number greater than zero to 1
				#because that number indicates the allele position but we have only one position
				#in the database
				if ( $fld_num == $gt_ind   ){
					#If this is a multiallelic variant then call the function get_gt_from_pl
					#to infer the biallelic genotype
					#I see if there are multiple variant ids in varids
					if( $res[0] =~ /,/ ){
							#print_and_log("Var:$varid Varids: ".$res[0]." applying get_gtfrompl..\n ",$log_file);
							$elem = get_gt_from_pl($cfg_hash,$sampleid,$res[0],$var_pos,$log_file);
					}
					#print_and_log("GT is: $elem\t",$log_file);
					if (defined $elem  and ($elem ne $vargenius_empty_val) and ($elem ne '.')){
						#Whatever is related with an additional alternative change to 1
						$elem =~ s/[2-9]/1/g;
						$out_str2 .=  $elem.$sep;
						#Now determine the Zygosity. Separate the GT using the slash or |
						my @gens = split("[\|,\/]",$elem);
						#If the alleles are the same the sample is homozygous
						if ( $gens[0] eq $gens[1]){
							#0/0
							if ( $gens[0] eq 0){
								$out_str2 .=  "HOMREF".$sep;
							}#Situation with ./.
							elsif ( $gens[0] eq '.' ) {
								$out_str2 .=  $vargenius_empty_val.$sep;
							}#1/1
							else{
									$out_str2 .=  "HOM".$sep;
							}
						#If alleles are different the sample is heterozygous	
						}else{#0/1 or 1/0 or 0|1 or 1|0
							$out_str2 .=  "HET".$sep;
						}
					}else{
							$out_str2 .=  $vargenius_empty_val.$sep;
					}
				}

				###################DP and GQ##########################
				#If the field is DP or GQ just insert the value and the sep
				if ( $fld_num == $gq_ind  or  $fld_num == $dp_ind  ){	
					#print_and_log("DP or GQ is: $elem\t",$log_file);
					if (defined $elem  and ($elem ne $vargenius_empty_val) and ($elem ne '.')){	
						$out_str3 .=  $elem.$sep;
					}else{
						$out_str3 .=  $vargenius_empty_val.$sep;
					}
				}
				$fld_num++;
				#print_and_log("\n",$log_file);
			}
		}else{
			print_and_log("WARNING: No result for the query $var_query...\t",$log_file);
			$out_str1 .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
			$out_str2 .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
			$out_str3 .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
		}#No result for the query
	}else{
					print_and_log("Sample id is not detectable...\t",$log_file);	
			}	
	#$out_str = $out_str1.$out_str2.$out_str3;
	return $out_str1,$out_str2,$out_str3;
}

=head2 get_gt_from_pl_vcf

 Title   : get_gt_from_pl_vcf
 Usage   : get_gt_from_pl_vcf( -database => 'name of the database,
                               );

 Function: Infers a biallelic genotype for a given alternate variant
					using the PL field of the VCF.
					Even if the genotype has been phased (indicated with |), if we make this computation
					it will loose the phasing (only / is used).
 
 Returns : a string with a genotype (0/1

=cut
sub get_gt_from_pl_vcf {
	my $cfg_hash = shift;
	my $gt = shift;
	my $pl = shift;	
	my $alt_ind = shift;
	my $log_file = shift;

	#PL threshold to infer the genotype
	my $pl_thr = $cfg_hash->{'pl_thr'};

	#print_and_log("Varids: $varids sampid= $sample_id varpos = $var_pos. GT= $gt, PL= $pl..\n ",$log_file);			
	
	#THe new genotype
	my $new_gt = "";
			
	if ( (defined $pl) and ($pl ne '') and ($gt ne '') and ($pl ne '.')){
		#Get the separator of the two elements of the genotype
		my $sep = "";
		if ( $gt =~ /\|/){
			$sep = "|";
		}else{$sep= "/";}
		#Get the elements
		my @g_els = split($sep,$gt);
		#Get the phred scaled likelihoods
		my @pl_vals = split(",",$pl);
		
		#Get the location where to start to take the pl values
		#(subtract 1 because it starts from 1)
		$alt_ind = ($alt_ind-1) * 3;
		#Get the PL values associated with this var_pos
		my $pl0 = $pl_vals[$alt_ind];
		my $pl1 = $pl_vals[$alt_ind+1];
		my $pl2 = $pl_vals[$alt_ind+2];
		
		#Get the minimum among the three that will specify the new genotype
		my $min = 100000000000;

		if ($pl0 < $min) {$min = $pl0;$new_gt = "0/0";}
		if ($pl1 < $min) {$min = $pl1;$new_gt = "0/1";}
		if ($pl2 < $min) {$min = $pl2;$new_gt = "1/1";}
		
		#Return a genotype only if the minimum reached is lower than a threshold
		if ($min > $pl_thr){
			$new_gt = "./.";
		}		
	}else{
			$new_gt = $gt;
	}								
	#print_and_log("Returning with new GT= $new_gt..\n ",$log_file);
	return $new_gt;
}


=head2 get_gt_from_pl

 Title   : get_gt_from_pl
 Usage   : get_gt_from_pl( -database => 'name of the database,
                               );

 Function: Infers a biallelic genotype for a given alternate variant
					using the PL field of the VCF.
					Even if the genotype has been phased (indicated with |), if we make this computation
					it will loose the phasing (only / is used).
 
 Returns : a string with a genotype (0/1

=cut
sub get_gt_from_pl {
	my $cfg_hash = shift;
	my $sample_id = shift;
	my $varids = shift;
	my $var_pos = shift;
	my $log_file = shift;

	#PL threshold to infer the genotype
	my $pl_thr = $cfg_hash->{'pl_thr'};
	
	my @fields = ($cfg_hash->{'db_sample_id'}, $cfg_hash->{'db_var_ids'});
	my @values = ($sample_id,"'".$varids."'");
	#print_and_log("Query: fields: ".$fields[0].",".$fields[1]." values= ".$values[0].",".$values[1]." ...\n ",$log_file);		
	#Get the GT value
	my $gt  = select_element_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_genotype_sample_table'},$cfg_hash->{'db_genotype_gt'}
					,\@fields,\@values);
	#Get the PL value
	my $pl  = select_element_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_genotype_sample_table'},$cfg_hash->{'db_genotype_pl'}
					,\@fields,\@values);
	#Get the varids
	#my $varids  = select_element_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
	#				$cfg_hash->{'db_pass'},$cfg_hash->{'db_genotype_sample_table'},$cfg_hash->{'db_var_ids'}
	#				,\@fields,\@values);
															
	#print_and_log("Varids: $varids sampid= $sample_id varpos = $var_pos. GT= $gt, PL= $pl..\n ",$log_file);			
	
	#THe new genotype
	my $new_gt = "";
			
	if ( (defined $pl) and ($pl ne '') and ($gt ne '') and ($pl ne '.')){
		#Get the separator of the two elements of the genotype
		my $sep = "";
		if ( $gt =~ /\|/){
			$sep = "|";
		}else{$sep= "/";}
		#Get the elements
		my @g_els = split($sep,$gt);
		#Get the phred scaled likelihoods
		my @pl_vals = split(",",$pl);
		
		#Get the location where to start to take the pl values
		#(subtract 1 because it starts from 1)
		$var_pos = ($var_pos-1) * 3;
		#Get the PL values associated with this var_pos
		my $pl0 = $pl_vals[$var_pos];
		my $pl1 = $pl_vals[$var_pos+1];
		my $pl2 = $pl_vals[$var_pos+2];
		
		#Get the minimum among the three that will specify the new genotype
		my $min = 100000000000;

		if ($pl0 < $min) {$min = $pl0;$new_gt = "0/0";}
		if ($pl1 < $min) {$min = $pl1;$new_gt = "0/1";}
		if ($pl2 < $min) {$min = $pl2;$new_gt = "1/1";}
		
		#Return a genotype only if the minimum reached is lower than a threshold
		if ($min > $pl_thr){
			$new_gt = "./.";
		}		
	}else{
			$new_gt = $gt;
	}								
	#print_and_log("Returning with new GT= $new_gt..\n ",$log_file);
	return $new_gt;
}



=head2 get_genotype_info_vcf

 Title   : get_genotype_info_vcf
 Usage   : get_genotype_info_vcf( -database => 'name of the database,
                               );

 Function: Gets information from the form_val related with the genotype for a given
						sample and a given variant
					
					If format_field = GT
						For multiallelic variant  (i.e. there are multiple variant 
						ids in varids) it calls the function get_gt_from_pl
						to infer the biallelic genotype
						
						Puts a tag to describe the zygosity
						If alleles are different the sample is heterozygous	
						If the alleles are the same the sample is homozygous
					If format_field = AD
						returns REF_reads, ALT reads, tot reads and percentage of ALT reads
					If format_field = GQ or DP
							returns the value of form_val
 Returns : a string with 
						if a value is missing will add $vargenius_empty_val

=cut
sub get_genotype_info_vcf {
	my $cfg_hash = shift;
	my $format_field = shift;
	my $form_val = shift;
	my $alt_ind = shift;
	my $mult_flag = shift;
	my $sep = shift;
	my $log_file = shift;
		
	
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	my $out_str = "";#Final output string
	
	$alt_ind = $alt_ind + 1;
	###################AD##########################
	#If the field is AD insert the first element, the element corrisponding
	#to the alternative allele and the last element
	# number of reference reads,number of alterate reads, total reads number, percentage alterate reads
	if ( $format_field eq 'AD'  ){	
		#print_and_log("AD is: $elem - varpos = $var_pos\t",$log_file);
		if (defined $form_val  and ($form_val ne $vargenius_empty_val) and ($form_val ne '.') and ($form_val ne '.,.')){	
			#print_and_log("Form_val is: $form_val\t",$log_file);
			my @all_dpts = split(",",$form_val);
			#Compute the percentage of alterate reads respect to the sum
			#of the number of altered reads and non altered
			my $perc_alterate = 0;
			my $readsum = $all_dpts[0] + $all_dpts[$alt_ind];
			if ( ($readsum ne 0) ){
				$perc_alterate = $all_dpts[$alt_ind]/$readsum;
			}else{
				$perc_alterate = $vargenius_empty_val;
			}
			$out_str .=  $all_dpts[0].$sep.$all_dpts[$alt_ind].$sep.$readsum.$sep.$perc_alterate.$sep;
		}else{
			$out_str .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
		}
	}
	
	###################GT##########################
	#If the field is GT than convert any number greater than zero to 1
	#because that number indicates the allele position but we have only one position
	#in the database
	if ( $format_field eq 'GT' ){
		#If this is a multiallelic variant then call the function get_gt_from_pl
		#to infer the biallelic genotype
		if( $mult_flag == 1 ){
				#print_and_log("Var:$varid Varids: ".$res[0]." applying get_gtfrompl..\n ",$log_file);
				my @gt_pl = split(/\Q$mult_ann_sep\E/,$form_val);
				$form_val = get_gt_from_pl_vcf($cfg_hash,$gt_pl[0],$gt_pl[1],$alt_ind,$log_file);
		}
		#print_and_log("GT is: $elem\t",$log_file);
		if (defined $form_val  and ($form_val ne $vargenius_empty_val) and ($form_val ne '.')){
			#Whatever is related with an additional alternative change to 1
			$form_val =~ s/[2-9]/1/g;
			$out_str .=  $form_val.$sep;
			#Now determine the Zygosity. Separate the GT using the slash or |
			my @gens = split("[\|,\/]",$form_val);
			#If the alleles are the same the sample is homozygous
			if ( $gens[0] eq $gens[1]){
				#0/0
				if ( $gens[0] eq 0){
					$out_str .=  "HOMREF".$sep;
				}#Situation with ./.
				elsif ( $gens[0] eq '.' ) {
					$out_str .=  $vargenius_empty_val.$sep;
				}#1/1
				else{
					$out_str .=  "HOM".$sep;
				}
			#If alleles are different the sample is heterozygous	
			}else{#0/1 or 1/0 or 0|1 or 1|0
				$out_str .=  "HET".$sep;
			}
		}
		#If no info is present, add empty values for GT and ZYG
		else{
			$out_str .=  $vargenius_empty_val.$sep.$vargenius_empty_val.$sep;
		}
	}

	###################DP and GQ##########################
	#If the field is DP or GQ just insert the value and the sep
	if ( $format_field eq 'GQ'  or  $format_field eq 'DP'  ){	
		#print_and_log("DP or GQ is: $elem\t",$log_file);
		if (defined $form_val  and ($form_val ne $vargenius_empty_val) and ($form_val ne '.')){	
			$out_str .=  $form_val.$sep;
		}else{
			$out_str .=  $vargenius_empty_val.$sep;
		}
	}
	chop($out_str);
	return $out_str;
}



=head2 get_gene_info

 Title   : get_gene_info
 Usage   : get_gene_info( -database => 'name of the database,
                               );

 Function: Gets information from the database related with a gene symbol
						HPO and RefSeq
 
 Returns : a string with HPO_id\tHPO_desc\tRefSeqID\t
						if a value is missing will add '-'

=cut
sub get_gene_info {
	my $cfg_hash = shift;
	my $dbh = shift;
	my $log_file = shift;
	my $genename = shift;
	
	#Get Columns number for fields in the genes table
	my $hpoids_ind = $cfg_hash->{'hpoids_ind'};
	my $refseq_ids_ind = $cfg_hash->{'refseq_ids_ind'};

	my $gdi_score_ind = $cfg_hash->{'gdi_score_ind'};
	my $rvis_score_ind = $cfg_hash->{'rvis_score_ind'};
	my $rvis_perc_ind = $cfg_hash->{'rvis_perc_ind'};
	my $entrez_id_ind = $cfg_hash->{'entrez_id_ind'};
	my $omim_id_ind = $cfg_hash->{'omim_id_ind'};
	
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	#Needed if the query fails or there is no value in input
	my $all_emtpy = "$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val".
	"\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t";#They are 8 '-' for 7 gene descriptions (1 additional is from refseq)
	
	#Get the fields to select from the configuration file
	#N.B. The header depends by the gene annotation fields where I do not include the HPO desc (because 
	#I fetch fields of the genes_table. If you change something here change also into get_all_genes_annotation()	
	my @sel_fields = split (",",$cfg_hash->{'gene_annotation_fields'});
	my $sel_f_str = "";
	foreach my $sel_fields (@sel_fields){
			$sel_f_str .= $cfg_hash->{$sel_fields}.",";
	}
	chop($sel_f_str);
	
	my $out_str = "";
	if ( $genename ne $vargenius_empty_val and $genename ne '.'){

		my $var_query = "SELECT $sel_f_str FROM genes WHERE genename='$genename';";	
		my @res = ();
		do_fetch_row_array($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$var_query,\@res);
		if (scalar (@res) > 0){
			my $fld_num = 1;
			foreach my $elem (@res){
				###################HPOS IDS##########################
				#If the field is hpoids
				if ( $fld_num == $hpoids_ind   ){
					if (defined $elem and ($elem ne $vargenius_empty_val)) {
						#Get the hpo ids singularly
						my @hpo_dbids = split(",",$elem);
						if ( scalar(@hpo_dbids) > 0 ){
							foreach my $hpo_dbid ( @hpo_dbids){
									#Obtain the hpo_id from the database given the hpo db id
									my $hpo_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_hpo_id'}
																	,$cfg_hash->{'db_phenotype_id'},$hpo_dbid);
									#The id is a string or -1 if not found
									if ( length $hpo_id > 2){
										$out_str .= $hpo_id.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
							#Put the HPO description
							foreach my $hpo_dbid ( @hpo_dbids){
									#Obtain the hpo_id from the database given the hpo db id
									my $hpo_desc = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_hpo_name'}
																	,$cfg_hash->{'db_phenotype_id'},$hpo_dbid);
									#Means the return is -1, element exists
									if ( length $hpo_desc > 2){
										$out_str .= $hpo_desc.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
						}else{
								$out_str .= "-\t-\t";
						}	
					}else{
							$out_str .= "-\t-\t";
					}

				}
				
				###################REFSEQ IDS##########################
				#If the field is refseq_ids_ind
				if ( $fld_num == $refseq_ids_ind ){
					if (defined $elem  and ($elem ne $vargenius_empty_val)){
						#Get the refseq ids singularly
						my @refseq_dbids = split(",",$elem);
						if ( scalar(@refseq_dbids) > 0 ){
							foreach my $refseq_dbid ( @refseq_dbids){
									#Obtain the refseq_id from the database given the refseq db id
									my $refseq_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_refseq_id'}
																	,$cfg_hash->{'db_transcript_id'},$refseq_dbid);
									#The id is a string or -1 if not found
									if ( length $refseq_id > 2){
										$out_str .= $refseq_id.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
						}else{
							$out_str .= "$vargenius_empty_val\t";
						}						
					}else{
						$out_str .= "$vargenius_empty_val\t";
					}
				}
				###################OTHER scores and IDENTIFIERS##########################
				#If the field is refseq_ids_ind
				if ( $fld_num == $gdi_score_ind or  $fld_num == $rvis_score_ind or
				    $fld_num == $rvis_perc_ind or  $fld_num == $entrez_id_ind or
				    $fld_num == $omim_id_ind){
					if (defined $elem  and ($elem ne $vargenius_empty_val)){
						$out_str .= "$elem\t";
					}else{
						$out_str .= "$vargenius_empty_val\t";
					}
				}
				$fld_num++;
			}
		}else{
				print_and_log( "WARNING: No result from query: $var_query..\n",$log_file);
				$out_str .= $all_emtpy;
		}	
	}else{$out_str .= $all_emtpy;}#They are 8 '-' for 7 gene descriptions (1 additional is from refseq)
	
	#remove the last tab
	chop($out_str);
	return $out_str;
}

=head2 get_allgene_info

 Title   : get_allgene_info
 Usage   : get_allgene_info( -database => 'name of the database,
                               );

 Function: Prints a table with information about all gene symbols
						HPO and RefSeq
 
 Returns : a string with HPO_id\tHPO_desc\tRefSeqID\t
						if a value is missing will add '-'

=cut
sub get_allgene_info {
	my $cfg_hash = shift;
	my $gen_ann_file = shift;
	my $log_file = shift;
	
	#Get Columns number for fields in the genes table
	my $gene_name_ind = 0;
	my $hpoids_ind = $cfg_hash->{'hpoids_ind'};
	my $refseq_ids_ind = $cfg_hash->{'refseq_ids_ind'};

	my $gdi_score_ind = $cfg_hash->{'gdi_score_ind'};
	my $rvis_score_ind = $cfg_hash->{'rvis_score_ind'};
	my $rvis_perc_ind = $cfg_hash->{'rvis_perc_ind'};
	my $entrez_id_ind = $cfg_hash->{'entrez_id_ind'};
	my $omim_id_ind = $cfg_hash->{'omim_id_ind'};
	
	my $tab = "\t";
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	#Needed if the query fails or there is no value in input
	my $all_emtpy = "$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val".
	"\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t$vargenius_empty_val\t";#They are 8 '-' for 7 gene descriptions (1 additional is from refseq)
	
	#Get the fields to select from the configuration file
	my @sel_fields = split (",",$cfg_hash->{'gene_annotation_fields'});
	my $sel_f_str = $cfg_hash->{'db_genes_name'}.",";
	foreach my $sel_fields (@sel_fields){
			$sel_f_str .= $cfg_hash->{$sel_fields}.",";
	}
	chop($sel_f_str);

	my $var_query = "SELECT $sel_f_str FROM genes;";
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
	$cfg_hash->{'db_pass'},$var_query);
	
	#Writes into a table
	open (GENE_ANN,">$gen_ann_file") or die "ERROR: Cannot open $gen_ann_file. The program will exit..\n";
	
	#To print the header I must add the hpo_desc that is not in @sel_fields
	#hence here I remove the hpo_id and write both
	print GENE_ANN $cfg_hash->{'db_genes_name'}.$tab.($cfg_hash->{shift(@sel_fields)}).$tab.$cfg_hash->{'db_hpo_name'};
	foreach my $field (@sel_fields){
			print GENE_ANN $tab.$cfg_hash->{$field};
	}
	print GENE_ANN "\n";
	
	#For each result print all the needed fields for a VCF
	foreach my $row (@$res) {
		
		my @arr = @$row;
		my $out_str = "";
		my $fld_num = 0;
		 #print join($tab, @$row), "\t.\n";
		foreach  my $elem (@arr) {
				###################HPOS IDS##########################
				#If the field is hpoids
				if ( $fld_num == $hpoids_ind   ){
					if (defined $elem and ($elem ne $vargenius_empty_val)) {
						#Get the hpo ids singularly
						my @hpo_dbids = split(",",$elem);
						if ( scalar(@hpo_dbids) > 0 ){
							foreach my $hpo_dbid ( @hpo_dbids){
									#Obtain the hpo_id from the database given the hpo db id
									my $hpo_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																	$cfg_hash->{'db_pass'},$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_hpo_id'}
																	,$cfg_hash->{'db_phenotype_id'},$hpo_dbid);
									#The id is a string or -1 if not found
									if ( length $hpo_id > 2){
										$out_str .= $hpo_id.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
							#Put the HPO description
							foreach my $hpo_dbid ( @hpo_dbids){
									#Obtain the hpo_id from the database given the hpo db id
									my $hpo_desc = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																	$cfg_hash->{'db_pass'},$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_hpo_name'}
																	,$cfg_hash->{'db_phenotype_id'},$hpo_dbid);
									#Means the return is -1, element exists
									if ( length $hpo_desc > 2){
										$out_str .= $hpo_desc.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
						}else{
								$out_str .= "-\t-\t";
						}	
					}else{
							$out_str .= "-\t-\t";
					}

				}
				
				###################REFSEQ IDS##########################
				#If the field is refseq_ids_ind
				if ( $fld_num == $refseq_ids_ind ){
					if (defined $elem  and ($elem ne $vargenius_empty_val)){
						#Get the refseq_id singularly
						my @refseq_dbids = split(",",$elem);
						if ( scalar(@refseq_dbids) > 0 ){
							foreach my $refseq_dbid ( @refseq_dbids){
									#Obtain the refseq_id from the database given the hpo db id
									my $refseq_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																	$cfg_hash->{'db_pass'},$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_refseq_id'}
																	,$cfg_hash->{'db_transcript_id'},$refseq_dbid);
									#The id is a string or -1 if not found
									if ( length $refseq_id > 2){
										$out_str .= $refseq_id.$cfg_hash->{'mult_ann_sep'};
									}
							}
							#Remove the last separator
							for (my $i = 0; $i < length($cfg_hash->{'mult_ann_sep'}); $i++ ){
								chop($out_str);
							}
							$out_str .= "\t";
						}else{
							$out_str .= "$vargenius_empty_val\t";
						}						
					}else{
						$out_str .= "$vargenius_empty_val\t";
					}
				}
				###################OTHER scores and IDENTIFIERS##########################
				#If the field is other
				if ( $fld_num == $gdi_score_ind or  $fld_num == $rvis_score_ind or
				    $fld_num == $rvis_perc_ind or  $fld_num == $entrez_id_ind or
				    $fld_num == $omim_id_ind or  $fld_num == $gene_name_ind ){
					if (defined $elem  and ($elem ne $vargenius_empty_val)){
						$out_str .= "$elem\t";
					}else{
						$out_str .= "$vargenius_empty_val\t";
					}
				}
				$fld_num++;
		} 
		#print GENE_ANN $arr[1].$tab.$arr[2].$tab.$arr[3].$tab.$arr[4].$tab.$arr[5].$tab.".".$tab.".".$tab.$arr[0]."\n";
		print GENE_ANN $out_str."\n";
	}
	close(GENE_ANN);	
}


		


=head2 import_annovarout_2_db_fast

 Title   : import_annovarout_2_db_fast
 Usage   : import_annovarout_2_db_fast( -database => 'name of the database,
                               );

 Function: Inserts in the database selected results from annotation with Annovar
			It first deletes all the information contained into the annotation table so that
			it does not do any check for items existence and simply inserts information.
			
					This function also adds the transcripts information by cutting the $tr_ann_field
					into different parts sepearated by comma. The first part will be separated in 5 fields
					and added to the db table.
					
 Returns : 

=cut
sub import_annovarout_2_db_fast {
	my $cfg_hash = shift;
	my $ann_out = shift;
	my $prog_used = shift;
	my $analysis_id = shift;
	my $log_file = shift;

	#Separators
	my $sep = "\t";
	my $db_sep = ",";
	my $tab = "\t";
	my $comma = ",";
	my $mult_ann_sep  = $cfg_hash->{'mult_ann_sep'};
	my $empty_val = $cfg_hash->{'vargenius_empty_val'};
	my @header_fields = ();
	
	#First thing that we do is to remove all elements from the annotation table
	delete_from_db ($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_var_annotations_table'});
	
	#Open output from the annotation
	open (ANN_OUT,"<".$ann_out) or die "ERROR: Cannot open $ann_out. The program will exit..\n";

	#Open a file that will be the table to import
	my $table_2_imp = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'db_var_annotations_table'};
	open (ANN_IN, ">".$table_2_imp)
		or die "ERROR: Cannot open $table_2_imp. The program will exit..\n";

	#Make a query to get all the variants compid so that you can match them from an hash 
	my $compids_hash;
	#to verify if it must be inserted
	my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_var_compid'}." FROM  ".$cfg_hash->{'db_variants_table'}.";";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'},$query);
	foreach my $row (@$res) {
		#Get the row values into an array
		my @arr = @$row;
		my $fld_num = 0;
		my $compid = "";
		my $varid = -1;
		#Each element of the array is now a field of the table
		foreach  my $elem (@arr) {
			if (defined $elem and $fld_num == 0) {
				$varid = $elem;
			}
			if (defined $elem and $fld_num == 1) {
				$compid = $elem;
			}	
			$fld_num++;		
		}
		#Put the elements into an hash
		$compids_hash->{$compid} = $varid;
		#print_and_log( "Putting  $varid into hash for: $compid\n",$log_file);#DEBUGCODE
	}
		
	
	#Numerates the line just to know what is the header of the file
	my $num_line = 0;
	my $tot_cols = 0;
	
	#Print the header of the table that will be imported
	my $fields_2_imp = $cfg_hash->{'db_var_id'}.$comma.$cfg_hash->{'fields_to_import'}.$comma.$cfg_hash->{'freq_fields_4_db'};
	print ANN_IN $fields_2_imp.",".$cfg_hash->{'pred_to_transform'}."\n";
	
	#The array of field names
	my @head_flds = ();
	#Needed generic fields to simply pick and add. So far I only added frequencies
	my @needed_flds = split(",",$cfg_hash->{'freq_fields_4_db'});
	my @needed_flds_inds = ();#indexes found 
	
	#Needed fields from gene annotation
	my @func_cols = split(",",$cfg_hash->{'func_cols'});#columns from the output table describing the variant function
	my @func_cols_inds = ();#indexes found for func columns
	my @exonic_cols = split(",",$cfg_hash->{'exonic_cols'});#descriptions for exonic variants into the func column

	#PREDICTIONS: this values need to be rearranged, hence need separate code
	my @pred_2_transf_cols = split(",",$cfg_hash->{'pred_to_transform'});
	my @pred_2_transf_cols_inds = ();
	
	
	#OtherInfo column used to define where the variants coordinates start
	my $var_info_start_col = get_col_index($ann_out,$cfg_hash->{'other_info_col'}) + $cfg_hash->{'start_info_ann_wo_db'};
	print_and_log("Parsing Annovar annotation into $ann_out and writing a table ($table_2_imp) to upload into database\n ",$log_file);
	#Go through each line of a TXT tab separated table
	while (my $line = <ANN_OUT>){
		#Remove last \n
		chop($line);
		#Split the line to get the field values
		my @fields = split($sep,$line);
			
		##########HEADER###########
		#If we are in the header
		if ($num_line == 0){

			$tot_cols = scalar(@fields);
			#Number of column will be inserted in the flds_touse
			my $num_col = 0;
			#Keep in count if a field has been already inserted
			my @insd_flds = ();
			#Get each field in the fields variable
			foreach my $head_fld (@fields){
				
				#print_and_log( "Headfield:$head_fld \n",$log_file);
				#Convert all - and . into _
				$head_fld =~ s/[\-\.\+]/\_/g;
				#Convert upper to lowercase
				$head_fld = lc( $head_fld);

				if ( grep {/\b$head_fld\b/} @needed_flds ) {
					#Add the index where to find the func cols
					push(@needed_flds_inds,$num_col);	
				}	
										
				if ( grep {/\b$head_fld\b/} @func_cols ) {
					#Add the index where to find the func cols
					push(@func_cols_inds,$num_col);	
				}											
				$num_col++;
			}
			################################
			#Get the columns names and select columns for predictions
			################################
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];	
				#Put all the field names in an array
				push(@header_fields,$field);
				#Store the indexes of the predictions to transform
				if ( grep {/\b$field\b/} @pred_2_transf_cols ) {
					#print "Putting the col $field ($col)in the array..\n";
					push(@pred_2_transf_cols_inds,$col);
				}
			}
		}##########HEADER###########		
		
	###########DATA############
		else{

			#Get the compid using the chr_pos_alt_ref
			my $compid = $fields[$var_info_start_col+$cfg_hash->{'vcf_chr_ind'}]."_".$fields[$var_info_start_col+$cfg_hash->{'vcf_pos_ind'}]."_".$fields[$var_info_start_col+$cfg_hash->{'vcf_ref_ind'}]."_".$fields[$var_info_start_col+$cfg_hash->{'vcf_alt_ind'}];
			my $varid = $compids_hash->{$compid};
			my $values = $varid.$comma;

			#print_and_log(" $varid ($compid:$num_line) - ",$log_file);#DEBUGCODE
			#print_and_log(" $varid ($num_line) - ",$log_file);#DEBUGCODE
			
			#Insert now the new field between the previous two
			for ( my $col=0; $col< scalar(@fields); $col++){

				my $field = $fields[$col];					
				#################################################
				######REARRANGE FUNCTION DETAILS COLUMNS#########
				#################################################
				if ( grep {/\b$col\b/} @func_cols_inds ) {
											

					my $func_column = $field;
					my $gene_column = $fields[$col+1];
					my $detail_column = $fields[$col+2];
					my $exonfunc_column = $fields[$col+3];
					my $aachange_column = $fields[$col+4];
					
					$values .= "$gene_column$comma";
					#If func is "exonic" add also the genedetail
					if ( (grep {/\b$func_column\b/} @exonic_cols))  {
						if ( $exonfunc_column ne $empty_val ){
							$values .= "$func_column;$exonfunc_column$comma";
						}	
					}else{
							$values .= "$func_column$comma";
					}
					my $aachange_first_field  = (split($comma,$aachange_column))[0];
					my @aachange_fields = split(":",$aachange_first_field);
					
					my $transcript = $empty_val;
					my $exon = $empty_val;
					my $nt_change = $empty_val;
					my $aa_change = $empty_val;
						
					if (defined $aachange_fields[1] ){
						$transcript = $aachange_fields[1];
					}
					if (defined $aachange_fields[2] ){
						$exon = $aachange_fields[2];
					}
					if (defined $aachange_fields[3] ){
						$nt_change = $aachange_fields[3];
					}
					if (defined $aachange_fields[4] ){
						$aa_change = $aachange_fields[4];
					}					
					$values .= "$transcript$comma$exon$comma$nt_change$comma$aa_change$comma";						

				}				
				#######################################
				######CHANGE PREDICTION SYMBOLS#########
				########################################
				elsif (grep {/\b$col\b/} @pred_2_transf_cols_inds) {
						#my $i = 0;
						#foreach my $h (@header_fields){
							#print "$i $h\n";
							#$i++;
						#}
						
						#print "\n\nselected:\n";
						#print join("\n", @pred_2_transf_cols_inds);
						if ( $field ne $empty_val){

							#print "Converting prediction at row $num_row col $col (".$header_fields[$col]."), field $field\n";#DEBUGCODE
							$values .= convert_prediction($cfg_hash,$header_fields[$col],$field).$comma;
						}else{
							$values .= $empty_val.$comma;
						}
				}
				
				#######################################
				######OTHER FIELDS NEEDED     #########
				#######################################
				elsif (grep {/\b$col\b/} @needed_flds_inds) {

						#print "\n\nselected:\n";
						#print join("\n", @needed_flds_inds);
						if ( $field ne $empty_val){

							#print "Adding field $field\n";#DEBUGCODE
							$values .= $field.$comma;
						}else{
							$values .= $empty_val.$comma;
						}
				}
								
			}

			#Check if the varid is numeric, the annotation could contain a row with
			#additional information
			if ( correct_type($varid,"positiveint") ){
				chop($values);
				print ANN_IN $values."\n";		
			}
			
		}
		$num_line++;	
	}	
	close(ANN_OUT);
	close(ANN_IN);
  
  #Finally load to db
  print_and_log(" Loading $table_2_imp to ".$cfg_hash->{'db_name'}.":".$cfg_hash->{'db_var_annotations_table'}."\n",$log_file);	
  load_csv_2_db($cfg_hash->{'db_name'},$cfg_hash->{'db_host'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
					$cfg_hash->{'db_var_annotations_table'},$table_2_imp,"YES","NO");	
}





=head2 get_non_covered_regions

 Title   : get_non_covered_regions
 Usage   : get_non_covered_regions(   );

 Function: 
#Using a specific dimension for a region, check
#any region to see if the sum of reads which coverage its bases
#is less of some threshold that defines if the region is covered or not

					
 Returns : nothing
=cut
sub get_non_cov_reg_from_perbase_coverage{
 my $cfg_hash = shift;
 my $filename = shift;
 my $region_len = shift;
 my $min_cov = shift;
 my $out_table = shift;
 my $log_file = shift;
 
 #Indexes
 my $chrom_col = $cfg_hash->{'cov_chrom_col'};
	my $start_col = $cfg_hash->{'cov_start_col'};
	my $end_col = $cfg_hash->{'cov_end_col'};
	my $code_col = $cfg_hash->{'cov_code_col'};
	my $base_num_col = $cfg_hash->{'cov_base_num_col'};
	my $read_l_col = $cfg_hash->{'cov_read_l_col'};
	
 #Get the lines of the file
	open (DATAFILE, $filename);
	my @lines = <DATAFILE>;
	close(DATAFILE);
	
	open (OUT, ">".$out_table);
	print OUT "CHROM\tSTART\tEND\tCOV\n";
	print_and_log( "Getting non covered regions...\n",$log_file);
	 print "(".scalar(@lines)." lines) lines count: ";
	 my $thresh = 100000;
		#Depending from the number of sequences the current state will be printed
		my $totlines = scalar(@lines);
		if ( $totlines < $thresh){
			$thresh = 100;
			if ( $totlines < $thresh){
				$thresh = 10;
				if ( $totlines < $thresh){
					$thresh = 1;
				}
			}
		}	
		
	#Loop through all the file with coverage information base per base	
	for (my $line=0; $line < scalar(@lines)-$region_len;$line++){

		#Prints the number of sequences evaluated
		if (($line >= 0) and ($line % $thresh == 0) ){
			print $line." - ";
		}
		
		#print ".";
		my $reads_sum = 0;
		#get the start and end of the region and check if it is the
		#same chromosome
		my $startBase = $line;
		my $endBase = $region_len+$line;
		my @fieldsSta = split ("\t",$lines[$startBase]);
		my @fieldsEnd = split ("\t",$lines[$endBase]);
		if ($fieldsSta[$chrom_col] eq $fieldsEnd[$chrom_col]){
			#Go line by line as this is the region base per base
			for (my $base = $line; $base <= $endBase; $base++){
				my @fields = split ("\t",$lines[$base]);
				#Sum the number of reads for this base
				chop($fields[$read_l_col]);
				$reads_sum = $reads_sum + $fields[$read_l_col];
			}
			#If the mean coverage is smaller than the minimum coverage, get the start
			#and end and the coverage
			my $mean_cov = $reads_sum/$region_len;
			if ( $mean_cov < $min_cov){
				#Separate the fields for the first base of the region and for the last
				my @fields1 = split("\t",$lines[$line]) ;
				my @fields2 = split("\t",$lines[$endBase]);
				#Print the start, end chromosome and coverage of this region of N nucleotides
				print OUT $fields1[$chrom_col]."\t".$fields1[$base_num_col]."\t".$fields2[$base_num_col]."\t".$mean_cov."\n";
			}			
		}
	}
	
	close(OUT);
}


#=head2 get_non_cov_reg_with_bedtools

 #Title   : get_non_cov_reg_with_bedtools
 #Usage   : get_non_cov_reg_with_bedtools(   );

 #Function: 
					#Gets non-covered regions using bedtools command from Aaron Quinlan
					#(https://twitter.com/aaronquinlan/status/421786507511205888)
					
					#Uses genomecov to get a BedGraph for the input bam
					#then awk to get only those lines where the coverage is less than the threshold
					#Then it merges the regions with bedtools merge
					#Finally intersects the regions with the target
					#For IntersectBed
					#Parameters used are:
						#-wb which allows to write the original entry in B for each overlap. 
							#Useful for knowing what A overlaps. Restricted by -f and -r.
						#-r is used to restrict the intersection of the target feature for at least the 95% of the
							#interval given
					
 #Returns : nothing
#=cut
#sub get_non_cov_reg_with_bedtoolsOLD{
	#my $cfg_hash = shift;
	#my $input_bam = shift;
	#my $noncovgenes_out = shift;
	#my $analysis_id = shift;
	#my $log_file = shift;
	
	##Variables for computing execution times
	#my $globalStart = time;
	#my $partTime  = time;
	#my $partDuration = undef;
	
	##Number of reads to be defined as LOW coverage
	#my $nc_threshold = $cfg_hash->{'noncov_threshold'};
	#my $nc_recip_overlap = $cfg_hash->{'noncov_recip_overlap'};
	#my $nc_cover_col = $cfg_hash->{'noncov_cover_col'};
	
	##UCSC genes bed file
	#my $ucsc_genes = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_genes_bed'};
	##Getting the target bed file using the folder for targets and the name contained in the database
	#my $target_bed = $cfg_hash->{'target_reg_f'}."/".
				#get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
				#$cfg_hash->{'db_analysis_id'},$analysis_id);
	#if (file_not_present($target_bed) ){ die "Cannot proceed Check: $target_bed.\n";}
	
	##my $out_reg = $out_path."_".extract_name($input_bam,1).".noncov$low_coverage";
		
	##Execute the command
	#print_and_log( "...getting bedgraph (GenomeCov)...",$log_file);#DEBUGCODE
	#my $bedgraph = extract_name($input_bam,'noext').".".$cfg_hash->{'bedgraph_ext'};
	#run_BEDTOOLS_genomeCoverageBed($cfg_hash,$input_bam,$bedgraph,$log_file);######
	##Execute AWK
	#print_and_log( "...filtering regions covered lt $nc_threshold (AWK)...",$log_file);#DEBUGCODE
	#my $bedgraph_nc = extract_name($input_bam,'noext')."_nc.".$cfg_hash->{'bedgraph_ext'};
	#my $command =	" awk "."'".'$4<'.$nc_threshold."' $bedgraph > $bedgraph_nc";
	##print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	#try_exec_command($command) or die "Unable to execute command: $command\n";#####
	##Execute the command
	##Merge intervals
	##Parameters for merge: take the fourth column (coverage) and do average 
	#print_and_log( "...merging adjacent regions (mergeBed)...",$log_file);#DEBUGCODE
	#my $merge_params = " -c $nc_cover_col -o min,max,mean";
	#my $merged_int = extract_name($bedgraph_nc,'noext')."_".$cfg_hash->{'bedgraph_ext'}."_merged.".$cfg_hash->{'bed_ext'};
	#run_BEDTOOLS_mergeBed($cfg_hash,$bedgraph_nc,$merged_int,$merge_params,$log_file);#######
	
	##Intersect the target with the non covered region and print both a and b entries 
	##to see how much is the overlap and also when overlap is missing	
	##Parameters used are:
	 ##  -a is the target file
	  ## -b is the merged bed file
		## -wao which allows to write the original entry in B and A . 
	#my $cov_out = extract_name($merged_int,'noext')."_"."noncov$nc_threshold.".$cfg_hash->{'bed_ext'};
	#my $inters_params = " -wao ";
	#print_and_log( "...intersecting with target (intersectBed)...",$log_file);#DEBUGCODE
	#run_BEDTOOLS_intersectBed($cfg_hash,$target_bed,$merged_int,$cov_out,$inters_params,$log_file);#######
		
	##Remove the column containing target region information
	#my @cols_to_remove = (3);
	##print_and_log( "\nRemoving column 3 (starting from 0) from $cov_out ..\n",$log_file);#DEBUGCODE
	#print_and_log( "...deleting target region info...",$log_file);#DEBUGCODE
	#delete_columns($cov_out,\@cols_to_remove);#######
	
	
	##Delete rows containing  -1\t-1. Rows of the target without intersection
	#print_and_log( "...removing rows without intersection...",$log_file);#DEBUGCODE
	#delete_rows_containing($cov_out,"-1\t-1",$cov_out.".temp");#####
	
	###Adding gene symbols
	##Intersect with ucsc genes regions to get the gene name	
	#my $nc_genes_out = extract_name($merged_int,'noext')."_genes.".$cfg_hash->{'bed_ext'};
	#print_and_log( "...intersection with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	#run_BEDTOOLS_intersectBed($cfg_hash,$cov_out.".temp",$ucsc_genes,$nc_genes_out,$inters_params,$log_file);
	
	#@cols_to_remove = (3,4,5);
	#print_and_log( "...deleting ucsc region start-end, sort and get unique lines...",$log_file);#DEBUGCODE
	#delete_columns($nc_genes_out,\@cols_to_remove);
	#$command =	" sort $nc_genes_out | uniq > $nc_genes_out.temp";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	#try_exec_command($command) or die "Unable to execute command: $command\n";	
		
	##print_and_log( "Moving $nc_genes_out.temp to $nc_genes_out..\n",$log_file);#DEBUGCODE	
	#move($nc_genes_out.".temp",$nc_genes_out);	
	#print_and_log( "...DONE!",$log_file);#DEBUGCODE
	#if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		#delete_file($bedgraph);
		#delete_file($merged_int);
		#delete_file($cov_out);
		#delete_file($cov_out.".temp");
		#delete_file($nc_genes_out.".temp");
	#}
#}


=head2 get_non_cov_reg_with_bedtools

 Title   : get_non_cov_reg_with_bedtools
 Usage   : get_non_cov_reg_with_bedtools(   );

 Function: 
					Gets non-covered regions using bedtools command from Aaron Quinlan
					(https://twitter.com/aaronquinlan/status/421786507511205888)
					
					Uses genomecov to get a BedGraph for the input bam
					then awk to get only those lines where the coverage is less than the threshold
					Then it merges the regions with bedtools merge
					Finally intersects the regions with the target
					For IntersectBed
					Parameters used are:
						-wb which allows to write the original entry in B for each overlap. 
							Useful for knowing what A overlaps. Restricted by -f and -r.
						-r is used to restrict the intersection of the target feature for at least the 95% of the
							interval given
					
 Returns : nothing
=cut
sub get_non_cov_reg_with_bedtools{
	my $cfg_hash = shift;
	my $input_bam = shift;
	my $noncovgenes_out = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Number of reads to be defined as LOW coverage
	my $nc_threshold = $cfg_hash->{'noncov_threshold'};
	my $nc_recip_overlap = $cfg_hash->{'noncov_recip_overlap'};
	my $nc_cover_col = $cfg_hash->{'noncov_cover_col'};
	
	#UCSC genes bed file
	my $ucsc_genes = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_genes_bed'};
	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = $cfg_hash->{'target_reg_f'}."/".
				get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
				$cfg_hash->{'db_analysis_id'},$analysis_id);
	if (file_not_present($target_bed) ){ die "Cannot proceed Check: $target_bed.\n";}
	
		
	#1.Generate the BEDGRAPH
	print_and_log( "...getting bedgraph (GenomeCov)...",$log_file);#DEBUGCODE
	my $bedgraph = extract_name($input_bam,'noext').".".$cfg_hash->{'bedgraph_ext'};
	run_BEDTOOLS_genomeCoverageBed($cfg_hash,$input_bam,$bedgraph,$log_file);######
	
	#2. Filter regions covered less than a threshold
	print_and_log( "...filtering regions covered lt $nc_threshold (AWK)...",$log_file);#DEBUGCODE
	my $bedgraph_nc = extract_name($input_bam,'noext')."_nc.".$cfg_hash->{'bedgraph_ext'};
	my $command =	" awk "."'".'$4<'.$nc_threshold."' $bedgraph > $bedgraph_nc";
	print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";#####

	#3. Merge intervals
	#Parameters for merge: take the fourth column (coverage) and do average 
	print_and_log( "...merging adjacent regions (mergeBed)...",$log_file);#DEBUGCODE
	my $merge_params = " -c $nc_cover_col -o min,max,mean";
	my $merged_int = extract_name($bedgraph_nc,'noext')."_".$cfg_hash->{'bedgraph_ext'}."_merged.".$cfg_hash->{'bed_ext'};
	run_BEDTOOLS_mergeBed($cfg_hash,$bedgraph_nc,$merged_int,$merge_params,$log_file);#######
	
	#4. Intersect the target with the non covered region and print both a and b entries 
	#to see how much is the overlap and also when overlap is missing	
	#Parameters used are:
	 #  -a is the target file
	  # -b is the merged bed file
		# -wao which allows to write the original entry in B and A . 
	my $cov_out = extract_name($merged_int,'noext')."_"."noncov$nc_threshold.".$cfg_hash->{'bed_ext'};
	my $inters_params = " -wao ";
	print_and_log( "...intersecting with target (intersectBed)...",$log_file);#DEBUGCODE
	run_BEDTOOLS_intersectBed($cfg_hash,$target_bed,$merged_int,$cov_out,$inters_params,$log_file);#######
		
	#5. Remove the column containing target region information
	my @cols_to_remove = (3);
	#print_and_log( "\nRemoving column 3 (starting from 0) from $cov_out ..\n",$log_file);#DEBUGCODE
	print_and_log( "...deleting target region info...",$log_file);#DEBUGCODE
	delete_columns($cov_out,\@cols_to_remove,$cov_out.".noTinfo");#######
	
	
	#6. Delete rows containing  -1\t-1. Rows of the target without intersection
	print_and_log( "...removing rows without intersection...",$log_file);#DEBUGCODE
	delete_rows_containing($cov_out.".noTinfo","-1\t-1",$cov_out.".noTinfo.no0rows");#####
	
	#7. Adding gene symbols
	#Intersect with ucsc genes regions to get the gene name	
	my $nc_genes_out = extract_name($cov_out,'noext')."_genes.".$cfg_hash->{'bed_ext'};
	print_and_log( "...intersection with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	run_BEDTOOLS_intersectBed($cfg_hash,$cov_out.".noTinfo.no0rows",$ucsc_genes,$nc_genes_out,$inters_params,$log_file);
	
	@cols_to_remove = (3,4,5);
	print_and_log( "...deleting ucsc region start-end, sort and get unique lines...",$log_file);#DEBUGCODE
	delete_columns($nc_genes_out,\@cols_to_remove,$nc_genes_out.".noUCSCinfo");
	$command =	" sort $nc_genes_out.noUCSCinfo | uniq > $nc_genes_out.noUCSCinfo.sort";
	print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";	
		
	#print_and_log( "Moving $nc_genes_out.temp to $nc_genes_out..\n",$log_file);#DEBUGCODE	
	move($nc_genes_out.".noUCSCinfo.sort",$nc_genes_out."noUCSCinfo.sort.bed");	
	print_and_log( "...DONE!",$log_file);#DEBUGCODE
	if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		delete_file($bedgraph);
		delete_file($merged_int);
		delete_file($cov_out);
		delete_file($cov_out.".temp");
		delete_file($nc_genes_out.".temp");
	}
}


=head2 get_non_cov_exons

 Title   : get_non_cov_exons
 Usage   : get_non_cov_exons(   );

 Function: 
					Gets non-covered exons using bedtools 
	
	After each interval in A, bedtools coverage will report:

	The number of features in B (reads) that overlapped (by at least one base pair) the A interval.
	The number of bases in A that had non-zero coverage from features in B.
	The length of the entry in A.
	The fraction of bases in A that had non-zero coverage from features in B.
					
 Returns : nothing
=cut
sub get_non_cov_exons{
	my $cfg_hash = shift;
	my $input_bam = shift;
	my $noncovexons_out = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = $cfg_hash->{'target_reg_f'}."/".
				get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
				$cfg_hash->{'db_analysis_id'},$analysis_id);
	if (file_not_present($target_bed) ){ die "Cannot proceed Check: $target_bed.\n";}

	#UCSC genes bed file
	my $ucsc_genes = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_genes_bed'};
			
	#Prepare Target Exons:
	my $ucsc_exons = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_exons_bed'};
	my $target_exons = extract_name($target_bed,"noext")."_".$cfg_hash->{'ucsc_exons_bed'};
	
	#if ( !(-e $target_exons) or (-z $target_exons) ){
		##1. Intersect the exons with the target so that we have the exons on-target
		##to see how much is the overlap and also when overlap is missing	
		##Parameters used are:
		 ##  -a the exons
		  ## -b the target file
			## -wa which allows to write the original entry in A . 
		#my $inters_params = " -wa ";
		#print_and_log( "...intersecting UCSC exons with target (intersectBed)...",$log_file);#DEBUGCODE
		#run_BEDTOOLS_intersectBed($cfg_hash,$ucsc_exons,$target_bed,$target_exons.".temp",$inters_params,$log_file);#######	
		
		##Get only the chrom and interval, sort and get uniq
		#print_and_log( "Get only the chrom and interval, sort and get unique lines...",$log_file);#DEBUGCODE
		#my $command =	"  cut -f1-3 $target_exons.temp | sort | uniq > $target_exons";
		#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
		#try_exec_command($command) or die "Unable to execute command: $command\n";			
	#}
	if ( -e $target_exons and !(-z $target_exons) ){
		#2: Get exon coverage
		my $target_exons_cov = extract_name($input_bam,"noext")."_".extract_name($target_exons,1).".".$cfg_hash->{'bed_ext'};
		#I give in input a genome file that has been obtained with 
		# samtools faidx ucsc.hg19.fa that produces ucsc.hg19.fa.fai
		#and then cut -f1,2  ucsc.hg19.fa.fai >  ucsc.hg19.genomefile
		#and with the -sorted parameter the computation will be very fast
		#
		#Target exons were previously sorted with
		#sortBed -i Agilent_SureSelect_CCP_v1_exons_UCSC.bed -faidx names.txt  > Agilent_SureSelect_CCP_v1_exons_UCSC.sort.bed
		#where names.txt is a file:
		#chrM
		#chr1
		#chr2
		#..
		#chrY
		my $covbedparams = " -g ".extract_name($cfg_hash->{'hum_ref'},"noext").".genomefile -sorted ";
		run_BEDTOOLS_coverageBed_gen($cfg_hash,$target_exons,$input_bam,$target_exons_cov,$covbedparams,$log_file);
		
		
		#3: Restrict to ipotetically low coverage exons only
		my $nc_threshold = 0.2;#percent of interval covered
		my $target_exons_nc = extract_name($target_exons_cov,"noext")."_nc.".$cfg_hash->{'bed_ext'};
		print_and_log( "Get only the chrom and interval, sort and get unique lines...",$log_file);#DEBUGCODE
		my $command =	" awk "."'".'$7<'.$nc_threshold."' $target_exons_cov > $target_exons_nc";
		print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
		try_exec_command($command) or die "Unable to execute command: $command\n";	

		#4: intersection with ucsc genes (intersectBed)
		#I use -wo because I want the dimension of the intersection and the Gene name.
		#But I don't want the regions that do not have intersection with Genes because I will not use them
		my $inters_params = " -wo ";
		print_and_log( "...intersecting UCSC exons with Genes...",$log_file);#DEBUGCODE
		run_BEDTOOLS_intersectBed($cfg_hash,$target_exons_nc,$ucsc_genes,$noncovexons_out,$inters_params,$log_file);#######	

		print_and_log( "...DONE!",$log_file);#DEBUGCODE
		if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			delete_file($target_exons_cov);
			delete_file($target_exons_nc);
		}
		
	}else{
		print_and_log( "ERROR: $target_exons file does not exist...",$log_file);#DEBUGCODE
	}

}
=head2 get_gene_4_interval

 Title   : get_gene_4_interval
 Usage   : get_gene_4_interval(   );

 Function: Obtains the GENE for each interval of the file given in input
					
 Returns : a BED file in which for each interval there is the GENE from UCSC
 
=cut
sub get_gene_4_interval {
	my $cfg_hash = shift;
	my $bed_input = shift;
	my $bed_out = shift;
	my $log_file = shift;
	
	my $inters_params = " -wao ";
	#UCSC genes bed file
	my $ucsc_genes = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_genes_bed'};
		
	#Intersect with ucsc genes regions to get the gene name	
	print_and_log( "Intersection  $bed_input with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	run_BEDTOOLS_intersectBed($cfg_hash,$bed_input,$ucsc_genes,$bed_out,$inters_params,$log_file);
	
	#delete useless fields
	my @cols_to_remove = (10,11,12);
	print_and_log( "...deleting ucsc region start-end, sort and get unique lines...",$log_file);#DEBUGCODE
	delete_columns($bed_out,\@cols_to_remove);
	my $command =	" sort $bed_out | uniq > $bed_out.temp";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#print_and_log( "Moving $nc_genes_out.temp to $nc_genes_out..\n",$log_file);#DEBUGCODE	
	move($bed_out.".temp",$bed_out) or die "Unable to move $bed_out.temp to $bed_out\n";	
	print_and_log( "...DONE!",$log_file);#DEBUGCODE
	if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		delete_file($bed_out.".temp");
	}	
		
}

=head2 get_all_samples_reads_number

 Title   : get_all_samples_reads_number
 Usage   : get_all_samples_reads_number(   );

 Function: Obtains the samples reads number (total, after duplicate removal and correctly aligned pairs against the reference
					from the stats folder (FlagStat Alignment Statistics BAM file aligned and sorted)
					and creates a plot for all samples in the same category (e.g. same project) identifiable with the same userid
					and the same target file
					
 Returns : a table with:
						SampleName TotalReads TotalRemoved ProperlyPaired
						
						moreover returns three variable: success, the total reads, the total removed and the total correctly paired for
						 the samples in this group
=cut
sub get_all_samples_reads_number {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $out_file = shift;
	my $steps_array = shift;
	my $log_file = shift;
	
	my $success = 0;

	#Get all the group ids from the same category of analyses
	my  @analysis_ids = ();
	#Obtain the userid from the database given the group id
	my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Obtain the target file name from the database given the group id
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
															
	#collects all the group ids related with this userid which could represent a project and that use the same target file
	my $query = "SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".
										$cfg_hash->{'db_analysis_userid'}." = $userid AND ".$cfg_hash->{'db_targetbed'}." = '$targetbed'".
										" AND ".$cfg_hash->{'db_analysis_infreq'}." = 1 ;";	
										
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $analysis_ids_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_analysis_id'});
	#Now put them into a normal array
	foreach my $all_analysis_id(keys %{$analysis_ids_fetch}){
		push(@analysis_ids,$all_analysis_id);
	}
	#Variable to store information for the subject group
	my $tot_reads_subject = "";
	my $prop_paired_subject = "";
	my $reads_removed_subject = "";
	
	####################Finished getting the analysis ids
	#Used to get the number of samples for which we have the kinship
	my $kinship_num = 0;
	my $all_samples_num = 0;
	if ( scalar(@analysis_ids) > 0 ){

		#Open the final table to write all the results
		open (OUT_FILE, ">$out_file") or print ("ERROR: Cannot open $out_file\n");
		#Open a file to store this information for the current analysis
		my $current_anal_stat_f = "$out_file.$analysis_id.totab";
		#open (SAMPLE_STATS, ">$current_anal_stat_f") or die ("Cannot open $current_anal_stat_f\n");
		#Print the header of this table
		print OUT_FILE "sample_name\tkinship\t".$cfg_hash->{'flagstat_totreads_field'}."\t".$cfg_hash->{'flagstat_readsremoved_field'}."\t".$cfg_hash->{'flagstat_propaired_field'}."\n";
		#print SAMPLE_STATS "sample_name\tkinship\t".$cfg_hash->{'flagstat_totreads_field'}."\t".$cfg_hash->{'flagstat_readsremoved_field'}."\t".$cfg_hash->{'flagstat_propaired_field'}."\n";
		
		#For each group id go into the stats folder, pick the sample_cumulative_coverage_proportions if exists
		#and get the columns given in @coverage_fields
		foreach my $all_analysis_id (@analysis_ids){
			#Obtain the group name from the database given the group id
			my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$all_analysis_id);	
																
			#Then for each sample of the group search the flagstat output
			#collects all the group ids related with this userid which could represent a project
			my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".
												$cfg_hash->{'db_analysis_id'}." = $all_analysis_id;";	
			#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
			my $sample_ids_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
			my @sample_ids = ();
			#Now put them into a normal array
			foreach my $sample_id (keys %{$sample_ids_fetch}){
				push(@sample_ids,$sample_id);
			}
			if ( scalar(@sample_ids) > 0 ){
				
				foreach my $sample_id (@sample_ids){
					my $at_least_one = 0;
					$all_samples_num++;
					#Obtain the group name from the database given the group id
					my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
																$cfg_hash->{'db_sample_id'},$sample_id);	

					#Verifying that there has been an execution of multiple samples
					my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    					
					#Here we get both the out of flagstat for the bam file after BWA and the same file
					#obtained after the duplicate removal
					
					#Aligned BAM
					my $flagstat_aln = "";

					#If there are more read_files per sample than search the merged output
					if ( $multiple == 1){					
						$flagstat_aln = $cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{$analysis_id.'_stats_out_f'}."/$sample_name\_",
												$cfg_hash->{'align_step'}.".".$cfg_hash->{'flagstat_ext'};
					}else{
						#The sample has a single read file. In this case we can access its name from the 
						#readf_table with the sample id that we are carrying 
						#print_and_log( "\n----There are no multiple samples..\n",$log_file);#DEBUGCODE
						my $params;
						my $main_name = $cfg_hash->{'db_readf_name'};
						my $step_needed = $cfg_hash->{'sort_idx_step'};
						getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},$sample_id);
						$flagstat_aln = $cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{$analysis_id.'_stats_out_f'}."/".build_input_name_from_executed($params,$step_needed,
											$params->{$main_name},$steps_array)."_".$step_needed.".".$cfg_hash->{'flagstat_ext'};		
					}
					#print_and_log("flagstat input: $flagstat_aln\n",$log_file);			
					
					my $flagstat_aln_sort_mrdup = $cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{$analysis_id.'_stats_out_f'}."/$sample_name\_"
												.$cfg_hash->{'mergebam_step'}."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'flagstat_ext'};
												
					#print_and_log( "Picking $flagstat_aln and $flagstat_aln_sort_mrdup ..\n",$log_file);#DEBUGCODE
				  
			
					##############################
					#Here I get a table with information from flagstat

					#Obtain the kinship from the database given the sample name
					my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
															$cfg_hash->{'db_analysis_id'}.','.$cfg_hash->{'db_sample_name'},$analysis_id.",'".$sample_name."'");							
					my $kinship_types = join("",split(",",$cfg_hash->{'kinship_types'}));
					my $new_row = "";
					if ( $kinship =~ /[$kinship_types]/){
						$kinship_num++;
						$new_row .= "$sample_name\t$kinship";
					}else{
						$new_row .= "$sample_name\t-";
					}					
					my $tot_reads_align = -1;
					my $tot_reads_mrdup = -1;
					my $reads_removed = -1;
					my $properly_paired = -1;
					
					#Get the number of total reads and properly paired from the file just after BWA
					if (-e $flagstat_aln and (!-z $flagstat_aln) ){
						#From the output of flagstat take the total read number and the properly paired
						#Get only the first column
						my $flagstat_firstcol = $flagstat_aln.".temp";
						extract_colnum_from_file_linux($flagstat_aln,1,$flagstat_firstcol,' ');
						#Put into an array and pick the columns related with the information wanted
						my @stats_list = list_to_array($flagstat_firstcol,'NO_NEW_LINE');
						$tot_reads_align = $stats_list[$cfg_hash->{'flagstat_totreads_row'}];
						$properly_paired = $stats_list[$cfg_hash->{'flagstat_propaired_row'}];
						 $at_least_one = 1;
					}else{$tot_reads_align = "-";}

					#Get the number of total reads and properly paired from the file after the duplicates removal
				  #The higher level of BAM used will be picked for the statistics
				  if (-e $flagstat_aln_sort_mrdup and (!-z $flagstat_aln_sort_mrdup)  ){
						
						#From the output of flagstat take the total read number and the properly paired
						#Get only the first column
						my $flagstat_firstcol = $flagstat_aln_sort_mrdup.".temp";
						extract_colnum_from_file_linux($flagstat_aln_sort_mrdup,1,$flagstat_firstcol,' ');
						#Put into an array and pick the columns related with the information wanted
						my @stats_list = list_to_array($flagstat_firstcol,'NO_NEW_LINE');
						$tot_reads_mrdup = $stats_list[$cfg_hash->{'flagstat_totreads_row'}];	
						$properly_paired = $stats_list[$cfg_hash->{'flagstat_propaired_row'}];
						 $at_least_one = 1;		
					}else{$tot_reads_mrdup = "-";}			
					
					#Compute the number of removed reads
					if ( $tot_reads_align ne '-' and $tot_reads_mrdup ne '-'){
						$reads_removed = $tot_reads_align - $tot_reads_mrdup;
					}
					
					#At least one, means that at least one among flagstat after BWA ore duplicate removal is found
					if ( $at_least_one > 0  ){		
						#print_and_log( "Opening : $flagstat_to_use\n",$log_file);#DEBUGCODE
						$success = 1;
											
						#Print into the output file
						print OUT_FILE "$new_row\t$tot_reads_align\t$reads_removed\t$properly_paired\n";
						##Store the information for the current analysis
						#if ( $all_analysis_id == $analysis_id){
							#print SAMPLE_STATS "$new_row\t$tot_reads_align\t$reads_removed\t$properly_paired\n";
						#}						
						#If this is the group id subject then stores the data in the given variables
						$tot_reads_subject .= "$tot_reads_align,";
						$reads_removed_subject .= "$reads_removed,";
						$prop_paired_subject .= "$properly_paired,";	
					}else{
						print_and_log( "WARNING : Cannot find any flagstat output to use for sample $sample_name\n",$log_file);
					}
				}
			}else{print_and_log( "WARNING : Cannot find any sample for group $all_analysis_id\n",$log_file);}

		}	
		close (OUT_FILE);
		
		#Sample Reads Alignment statistics
		#close (SAMPLE_STATS);
		#Obtain the analysis name of the current analysis from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);			
		#append_str_2_file_if_path_notexist($cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{'finalout_out_f'}."/".$cfg_hash->{'statstables_file'},"Alignment Statistics\t".$current_anal_stat_f);
	
	}
	#Remove last comma
	if ( $tot_reads_subject ne '' and $prop_paired_subject ne '' and $reads_removed_subject ne ''){
		chop($tot_reads_subject );
		chop($prop_paired_subject );
		chop($reads_removed_subject);
	}
	#If all samples have a kinship info return 2 instead of 1
	if ( $kinship_num == $all_samples_num){
		$success = 2;
	}	
	
	return $success,$tot_reads_subject,$reads_removed_subject,$prop_paired_subject;
}

#=head2 get_statistics_comparison

 #Title   : get_statistics_comparison
 #Usage   : get_statistics_comparison(   );

 #Function: Compares the statistics obtained for the given analysis
					#with the statistics present into the database:
					#- target coverage
					#- duplicate removed
					
 #Returns : a table where the comparison is printed
#=cut
#sub get_statistics_comparison {
	#my $cfg_hash = shift;
	#my $analysis_id = shift;
	#my $outFolder = shift;
	#my $log_file = shift;

	#my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						#$cfg_hash->{'db_analysis_id'},$analysis_id);
																	
	##Get the average number of duplicates
	#my $var_query = "SELECT AVG (".$cfg_hash->{'db_sampstats_duplicates'}.") FROM ".$cfg_hash->{'db_sampstats_duplicates'}.
							#" WHERE ".$cfg_hash->{'db_analyses_table'}." IN ".
							#" ( SELECT ".$cfg_hash->{'db_analyses_id'}." FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".
							#" ".$cfg_hash->{'db_targetbed'}." = $targetbed);";	
							
	#my @res = ();
	#do_fetch_row_array($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$var_query,\@res);
	#if (scalar (@res) > 0){
		#my $fld_num = 1;
		#foreach my $elem (@res){
			
		#}
	#}
	
#}

=head2 get_statistics_tables

 Title   : get_statistics_tables
 Usage   : get_statistics_tables(   );

 Function: Obtains tables into tab separated file using statistics obtained during
					the analyses
					 -  sample coverage table (from DepthOfCoverage .sample_summary file
					
 Returns : a file with a list of paths of tables
=cut
sub get_statistics_tables {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $outFolder = shift;
	my $log_file = shift;

	my $tables_list_file = $outFolder."/".$cfg_hash->{'statstables_file'};
	#Obtain the analysis name from the database given the group id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);		
	#Get the file with the summary of samples coverage
	my $sam_cov_summ_f = $cfg_hash->{$analysis_id.'_stats_out_f'}."/$analysis_name\_DOC.".$cfg_hash->{'DOC_sam_sum'};
	#Define the name for the table to be printed
	my $sam_cov_summ_totab = $sam_cov_summ_f.".totab";

	if ( -e $sam_cov_summ_f) {
		#print_and_log( "Extracting a table to be shown in HTML from $sam_cov_summ_f..\n",$log_file);
		#print_and_log( "Whose path will be printed in  $tables_list_file..\n",$log_file);#DEBUGCODE

		#From the table remove useless lines
		my $sam_cov_summ_f_temp = $sam_cov_summ_f.".temp";
		delete_rows_containing($sam_cov_summ_f,"^Total",$sam_cov_summ_f_temp);
		
		#ELIMINARE							
		###From the sample summary file remove the last row starting with "Total"
		##my $command = 'grep -v '."'".'^Total'."' ".$sam_cov_summ_f.' > '.$sam_cov_summ_f.".temp";
		###print_and_log( "Executing: $command..\n",$log_file);
		#try_exec_command($command) or die "Unable to execute command: $command\n";	
		
		#Now get only the 0,1,2,6,7,8,9,10 columns
		extract_columns_from_file($sam_cov_summ_f_temp,"0,1,2,6,7,8,9,10",$sam_cov_summ_totab);
		
		#Add the output to the file list
		append_str_2_file_if_path_notexist($tables_list_file,$cfg_hash->{'DOC_sam_sum'}."\t".$sam_cov_summ_totab);
	}else{
		print_and_log( "WARNING: Cannot find: $sam_cov_summ_f. get_statistics_tables will not be executed.\n",$log_file);
	}
	
	
	###############################################
	#Put coverage statistics into the database
	#############################
	#From the sample summary file remove the last row starting with "sample_id"
	my $command = 'grep -v '."'".'^sample_id'."' ".$sam_cov_summ_f.'.temp > '.$sam_cov_summ_f.".temp2";
	#print_and_log( "Executing: $command..\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
	my $cov_stat_f = $sam_cov_summ_f.".cov.temp";
	extract_columns_from_file($sam_cov_summ_f.".temp2","0,6",$cov_stat_f);
	my $cov_stat_hash;
	twocolfile_to_hash($cov_stat_f,\$cov_stat_hash);
	
	foreach my $samplename (keys %$cov_stat_hash){
		#get the sampleid
		my $sampleid =  get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_name'},"$analysis_id,'$samplename'");		
		#add element tp statistics table if the element do not exists, otherwise update it
		my $id = -1;
		my $fields =  $cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_sampstats_avgcov'}.",".$cfg_hash->{'db_analysis_id'};
		my $values = "$sampleid,".$cov_stat_hash->{$samplename}.",$analysis_id";
		my ($new_id,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_statistics'},$cfg_hash->{'db_sample_id'},
								$fields,$values,$cfg_hash->{'db_sample_id'},$sampleid);

		#If the sampleid exists insert that, otherwise the new one
		if ( $new_id >= 0){
			print_and_log( "Inserted the target coverage (".$cov_stat_hash->{$samplename}.") into ".
				$cfg_hash->{'db_sample_statistics'}." for sample $sampleid\n",$log_file);#DEBUGCODE
		}elsif ( $exist_id >= 0) {
			update_table_2($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
						$cfg_hash->{'db_sample_statistics'},$fields,$values,$cfg_hash->{'db_sample_id'},$sampleid);
		}else{
			log_and_exit("ERROR: during the import of the target coverage for $samplename ($sampleid : analysisid=$analysis_id). Cannot understand what is happen with the db...\n",$log_file);
		}				
	}

	if ( $cfg_hash->{'remove_temp'} eq 'YES' ){ 
		#Remove the temporary files
		delete_file($sam_cov_summ_f.".temp");
		delete_file($sam_cov_summ_f.".temp2");
		delete_file($sam_cov_summ_f.".cov.temp");
	}
}

=head2 get_all_samples_coverage_table

 Title   : get_all_samples_coverage_table
 Usage   : get_all_samples_coverage_table(   );

 Function: Obtains the samples coverage levels from the stats folder (sample_cumulative_coverage_proportions)
					and creates a plot for all samples in the same category (e.g. same project) identifiable with the same userid
					and the same target file
					
 Returns : a table with coverage levels given in the configuration file
=cut
sub get_all_samples_coverage_table {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $out_file = shift;
	my $input_suffix = shift;
	my $log_file = shift;
	
	my $success = 0;
	#The needed coverage level fields into the DepthOfCoverage result
	my @sam_cum_cov_prop_fields = split(",",$cfg_hash->{'DOC_sam_cum_cov_prop_fields'});
		
	#Get all the group ids from the same category of analyses
	my  @analysis_ids = ();
	#Obtain the userid from the database given the group id
	my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Obtain the target file name from the database given the group id
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
															
	#collects all the group ids related with this userid which could represent a project and that use the same target file
	my $query = "SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".
										$cfg_hash->{'db_analysis_userid'}." = $userid AND ".$cfg_hash->{'db_targetbed'}." = '$targetbed'  ".
										" AND ".$cfg_hash->{'db_analysis_infreq'}." = 1 ;";	
										
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $analysis_ids_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_analysis_id'});
	#Now put them into a normal array
	foreach my $analysis_id (keys %{$analysis_ids_fetch}){
		push(@analysis_ids,$analysis_id);
	}
	####################Finished getting the group ids
	#Used to get the number of samples for which we have the kinship
	my $kinship_num = 0;
	my $all_samples_num = 0;
	if ( scalar(@analysis_ids) > 0 ){
		
		#Open the final table to write all the results
		open (OUT_FILE, ">$out_file") or warn ("WARNING: Cannot open $out_file\n");
		
		#Print the header of this table
		print OUT_FILE "sample_name\tkinship\t".join("\t", @sam_cum_cov_prop_fields)."\n";
		
		#For each group id go into the stats folder, pick the sample_cumulative_coverage_proportions if exists
		#and get the columns given in @coverage_fields
		foreach my $analysis_id (@analysis_ids){
			#Obtain the group name from the database given the group id
			my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);	
			my $sam_cum_cov_prop_f = $cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{'stats_out_f'}."/$analysis_name\_DOC.".$input_suffix;
			
			
			if (-e $sam_cum_cov_prop_f){
				my @sam_cum_cov_prop_f_inds = ();
				#print_and_log( "Opening : $sam_cum_cov_prop_f\n",$log_file);#DEBUGCODE
				$success = 1;
				open (DATA, "<$sam_cum_cov_prop_f") or die ("Cannot open $sam_cum_cov_prop_f\n");
				my $num_row =1;
				#Write a new file with the given columns only
				while ( my $row = <DATA> ) {
					my $new_row = "";
					chop($row);
					my @fields = split /\t/, $row;
					#Go in the header and get the indexes of the columns needed
					if ( $num_row == 1){
						
						for ( my $col=0; $col< scalar(@fields); $col++){
							#Get the field name
							my $field = $fields[$col];
							#Put the index of the columns in an array
							#The first field in the header of the sam_cum_cov_prop_f is the sample name but is empty and leads
							#to error. So I avoid it here
							if ( $field ne ''  and grep {/\b$field\b/} @sam_cum_cov_prop_fields ) {
								#Add the index where to find the coverage levels needed
								push(@sam_cum_cov_prop_f_inds,$col);
							}
						}
						#print_and_log( "Num of columns to keep : ".scalar(@sam_cum_cov_prop_f_inds)." \n",$log_file);#DEBUGCODE
					}
					#We are now in the data
					else{
						#From the first field get the sample name and get the kinship
						my $sample_name = $fields[0];
						
						$all_samples_num++;
						#Obtain the kinship from the database given the sample name
						my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
									$cfg_hash->{'db_analysis_id'}.','.$cfg_hash->{'db_sample_name'},$analysis_id.",'".$sample_name."'");							
						my $kinship_types = join("",split(",",$cfg_hash->{'kinship_types'}));
						#print_and_log( "Sample name : $sample_name - Kinship: $kinship\n",$log_file);#DEBUGCODE
						if ( $kinship =~ /[$kinship_types]/){
							$kinship_num++;
							$new_row .= "$sample_name\t$kinship\t";
						}else{
							$new_row .= "$sample_name\t-\t";
						}
						
						#Insert the data given the needed indexes
						for ( my $col=0; $col< scalar(@fields); $col++){
							my $field = $fields[$col];
							if ( grep {/\b$col\b/} @sam_cum_cov_prop_f_inds ) {
								$new_row .= $field."\t";
							}
						}
						chop($new_row);#Remove last tab
						print OUT_FILE "$new_row\n";
					}
					$num_row++;
				}
				close(DATA);			
			}
			
		}	
		close (OUT_FILE);
	}
	#If all samples have a kinship info return 2 instead of 1
	if ( $kinship_num == $all_samples_num){
		$success = 2;
	}
	return $success;
}



=head2 get_gene_coverage_above

 Title   : get_gene_coverage_above
 Usage   : get_gene_coverage_above(   );

 Function: Given the analysisid and a file with a list of selected genes
					picks from the output of DepthOfCoverage the columns of the samples
					which indicate the coverage above $percentage_above and generates
					a table with all these genes and columns. If $sel_genes_f= ALL_GENES
					retains all the genes
				
 Returns : a file path
=cut
sub get_gene_coverage_above {
		my $cfg_hash = shift;
		my $analysis_id = shift;
		my $sel_genes_f = shift;
		my $percentage_above = shift;
		my $out_genecov_above = shift;
		my $statsFolder = shift;
		my $log_file = shift;
		
		#Check if selected genes are inputed
		if ( $sel_genes_f eq ''){
				$sel_genes_f  = "ALL_GENES";
		}
		
		#Get the gene summary file name
		#Obtain the analysis name from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
		my @doc_outs = split(",",$cfg_hash->{'DOC_outputs'});	
		my $gene_summary = $statsFolder."/".$analysis_name."_DOC.".$doc_outs[0];
		
		#If the gene summary exists..
		if ( -e $gene_summary  and !(-z $gene_summary)){
			#Getting the samples list associated to the groupid
			my $samples_h;#Samples hash to be reordered
			my @sort_samples = ();#The array with samples id reordered
			#Get the samples names involved for the group
			my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
			#print "Executing: $query\n";#DEBUGCODE
			my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
			#Print the header for information about the genotype for the sample
			my $sample_gen_fields1 = "";
			my $sample_gen_fields2 = "";
			#Get the kinship, to make the resorting
			foreach my $sample_name (keys %{$group_sam}){  #Print version

				#Obtain the kinship from the database given the sample id
				my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
																$cfg_hash->{'db_sample_name'},"'".$sample_name."'");	
				#Build an hash to reorder
				$samples_h->{$sample_name}->{'id'} = $sample_name;
				$samples_h->{$sample_name}->{'k'} = $kinship;
			}
			#Sort the samples as in samples_order parameter
			sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);	
			
			my $samples_l = join(",",@sort_samples);
			
			###############################
			# Rscripts
			###############################

			my $RLogPath = $statsFolder."/log/$analysis_name\_genecovabove$percentage_above".$cfg_hash->{'R_log_file'};
			my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_coverage_stats_script'};

			my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args TABLECOVABOVE $gene_summary $sel_genes_f $samples_l above_$percentage_above $out_genecov_above' $R_script $RLogPath";
		 
			print "The command is:\n $command\n";
			try_exec_command($command) or R_die($command,$R_script);			
		}else{
			print_and_log("WARNING: Cannot get gene coverage above $percentage_above. $gene_summary empty or non-existent... \n",$log_file);
		}

		
		#return $group_genecove_above;
}



=head2 add_genes_to_refseq_ids

 Title   : add_genes_to_refseq_ids
 Usage   : add_genes_to_refseq_ids(   );

 Function: Fetches from the database the genes associated to the refseq ids
					in the table. The refseq id is cut by using a split on the '_'
					and taking the first two elements
				
					Depending by the target file the information contained at the specified
					field (noncov_targ_name_col) will be different.
					Hence, I decided to get the gene only when the RefSeq transcript id
					is present. In some target the gene name is present as an information
					thus this step is not needed.
					
 Returns : nothing
=cut
sub add_genes_to_refseq_ids{
	my $cfg_hash = shift;
	my $noncov_tab = shift;
	my $out_table = shift;
	my $log_file = shift;	
	
	my $nc_targ_chr_col  = $cfg_hash->{'noncov_targ_chr_col'};
	my $nc_targ_name_col = $cfg_hash->{'noncov_targ_name_col'};
	my $nc_targ_start_col = $cfg_hash->{'noncov_targ_start_col'};
	my $nc_targ_end_col = $cfg_hash->{'noncov_targ_end_col'};
	my $nc_min_cov_col = $cfg_hash->{'noncov_min_cov_col'};
	my $nc_max_cov_col = $cfg_hash->{'noncov_max_cov_col'};
	my $nc_avg_cov_col = $cfg_hash->{'noncov_avg_cov_col'};
	my $nc_overlap_col = $cfg_hash->{'noncov_overlap_col'};
	
	
	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
	#The header
	my $header = "CHR\tSTART\tEND\tLocus\tMinCov\tMaxCov\tAverageCov\tGeneSymbol\n";
	###Inserting gene symbols using  a table with transcripts
	print "Reading from non covered regions file $noncov_tab ..\n";
	open(NONCOV,"<".$noncov_tab) or die "ERROR: Cannot open file $noncov_tab\n";
	open(OUT,">".$out_table) or die "ERROR: Cannot open file $out_table\n";
	print OUT $header;
	while  (my $row = <NONCOV>){
		chop($row);
		my @fields = split("\t",$row);
		my $genes = "-";
		if ( $fields[$nc_targ_name_col] ne '-'){
			#First separate using the comma. If there are more than two elements in the result do not
			#pick the gene name
			my @fields_in_field = split (",",$fields[$nc_targ_name_col]);
			if ( scalar (@fields_in_field) == 1){
				#Separate the field
				my @parts = split ("_",$fields[$nc_targ_name_col]);
				#If the field contains a RefSeq transcript id (e.g. NM_152486_cds_5_0_chr1_874420_f)
				#it has the NM_XXX value that we can find with the _ and the NM
				if (scalar (@parts) > 2){
					if ( $parts[0] eq 'NM'){
						#Get the refseq id
						my $rs_id = uc($parts[0]."_".$parts[1]);
												
						#Obtain the gene from the database given the refseq db id
						$genes = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_tr_genes_ids'},
															$cfg_hash->{'db_refseq_id'},"'".$rs_id."'");
																	
						if ( (length $genes) <= 2 ){
							$genes= "-";
						}	
					}
				}
				#if there is only one element in the field but it is not refseq. We just write that as gene information.
				else{
					$genes = $fields[$nc_targ_name_col];
				}
			}#Field contains more transcript ids or gene ids
			else{
					$genes = $fields[$nc_targ_name_col];
			}
		}#Field is empty
		else{
					$genes = "-";
		}
		#Print table info								
		print OUT $fields[$nc_targ_chr_col]."\t".$fields[$nc_targ_start_col]."\t".$fields[$nc_targ_end_col]."\t".$fields[$nc_targ_name_col]."\t".
							$fields[$nc_min_cov_col]."\t".$fields[$nc_max_cov_col]."\t".$fields[$nc_avg_cov_col]."\t$genes\n";	
	}
	
	#Disconnect db
  $dbh->disconnect(); 		
	
	close(NONCOV);
	close(OUT);
}

################################FUNCTIONS USED TO MODIFY THE FINAL OUTPUT

=head2 splicing_var_dist_2_final_out

 Title  : splicing_var_dist_2_final_out
 Usage  : splicing_var_dist_2_final_out( - filePath => 'the path to the file');
																			func_column => the column with the type of variant
																			detail_col => the column with the detail of the variant
 
 Function: starting from the  final output looks in func_column if the variant
					is splicing and gets the distance of the variant from the exon
  
 Returns : a new output file with one more column

=cut		
sub splicing_var_dist_2_final_out {
	my $cfg_hash = shift;
	my $file = shift;
	my $outFile = shift;
	
	my $func_column = $cfg_hash->{'spl_func_cols'};
	my $detail_col = $cfg_hash->{'spl_detail_col'};
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	open (OUT_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");

	
	my $func_column_ind = -1;
	my $detail_col_ind = -1;
	my $num_row = 1;
	my @splicing_cols = split($cfg_hash->{'parameters_sep'},$cfg_hash->{'splicing_cols'});
	my @func_cols = split($cfg_hash->{'parameters_sep'},$cfg_hash->{'spl_func_cols'});
	my @func_cols_inds = ();
	
	my $semicolon = ":";
	#Write a new file with the given columns only
	while ( my $row = <DATA> ) {
		my $new_row = "";
		chomp($row);
		my @fields = split /\t/, $row;
		#Go in the header and get the two indexes of the columns needed
		if ( $num_row == 1){
			#my $col_num = 0;
			#Insert now the new field between the previous two
			for ( my $col=0; $col< scalar(@fields); $col++){
				my $field = $fields[$col];
				$new_row .= $field."\t";
				#Put the index of the func col in an array
				if ( grep {/\b$field\b/} @func_cols ) {
					#Add the index where to find the func cols
					push(@func_cols_inds,$col);
					#Print both the header for the detail and for distance
					$col++;
					$new_row .= $fields[$col]."\t".$fields[$col].$cfg_hash->{'dist_from_exon_suff'}."\t";
				}
			}			
		}#Here we are in the data
		else{
			#my $col_num = 0;
			#Insert now the new field between the previous two
			for ( my $col=0; $col< scalar(@fields); $col++){
				my $field = $fields[$col];
				$new_row .= $field."\t";
				#If this index is among the indexes selected of func
				if ( grep {/\b$col\b/} @func_cols_inds ) {
					my $func_column = $field;
					my $detail_column = $fields[$col+1];
					$new_row .= $detail_column."\t";
					#Go in the next column and see if there is the splicing info
					if ( (grep {/\b$func_column\b/} @splicing_cols) and ( $detail_column ne $vargenius_empty_val) 
							and ( $detail_column ne $vargenius_empty_val) ){
						#Get the distance from the first result 
						#NM_001269053:exon7:c.505-27A>C]---[NM_138812:exon8:c.505-27A>C
						my @results = split(/\Q$mult_ann_sep\E/,$detail_column);
						#NM_001269053:exon7:c.505-27A>C
						my $first_res = shift @results;
						my @parts = split(/\Q$semicolon\E/,$first_res);
						#c.505-27A>C
						my $last = $parts[2];
						#Remove the c.505 => -27A>C
						$last =~  s/c\.\d+//g;
						#Get only the distance -27A>C => -27
						$last =~ /([+-]\d+)/;
						$new_row .= $1."\t";
					}else{
						$new_row .= $vargenius_empty_val."\t";
					}					
					#jump the detail column
					$col++;
				}
			}				
		}
		chop($new_row);
		print OUT_FILE "$new_row\n";
		$num_row++;
	}

	close(DATA);
	close(OUT_FILE);	
}



###########################################################################
###################REARRANGEMENT OF THE OUTPUT FILE#######################
##########################################################################

=head2 get_disease_gene_association

 Title  : get_disease_gene_association
 Usage  : get_disease_gene_association( 
 
 Function: 				
				
 Returns : returns an hash in which, for each gene there will be an association 
					to panels of diseases

=cut		
sub get_disease_gene_association {
	my $cfg_hash = shift;
	my ($gene_hash) = shift;
	my $log_file = shift;

	#Put into an hash of arrays the genes of each panel
	my @sel_panels_l = split(",",$cfg_hash->{'disease_genes_lists_f'});

	#The folder where are located all the gene panels
	my $gene_panels_fold = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'gene_panels_fold'};
	my %panels_genes;

	print_and_log("Getting disease_gene_association for panels in $gene_panels_fold... \n",$log_file);
	foreach my $gene_panel ( @sel_panels_l ){		
		my @panel_genes = list_to_array($gene_panels_fold."/".$gene_panel,'NO_NEW_LINE');
		$panels_genes{$gene_panel} = \@panel_genes;
		#print_and_log("panel $gene_panel contains: ".scalar(@panel_genes).". .. \n",$log_file);
	}

	#Get all the genes from the database
	my $genes;
	my $query = "";
	$query = "SELECT ".$cfg_hash->{'db_genes_name'}." FROM  ".$cfg_hash->{'db_genes_table'}.";";	
	print "Executing: $query\n";
	$genes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
						$query,$cfg_hash->{'db_genes_name'});			
	#For each gene
	foreach my $gene_name (keys %{$genes}){
		foreach my $gene_panel ( @sel_panels_l ){
			my $gp = $panels_genes{$gene_panel};
			if ( grep {/\b$gene_name\b/} @$gp ){
				
				if (! defined $$gene_hash->{$gene_name}) {
					$$gene_hash->{$gene_name} = $gene_panel;
				}else{
					$$gene_hash->{$gene_name} .= $cfg_hash->{'mult_ann_sep'}.$gene_panel;	
				}
				#print_and_log("Gene: $gene_name is in panel $gene_panel.. \n",$log_file);
			}
		}
	}
}



=head2 rearrange_rawtabular_out

 Title  : rearrange_rawtabular_out
 Usage  : rearrange_rawtabular_out( - filePath => 'the path to the file');
																			func_column => the column with the type of variant
																			detail_col => the column with the detail of the variant
 
 Function: 				
					- when it finds the gene symbol knows that after that there is the function details
						that will be rearranged by dividing the first result in different columns 
						In some cases there will be splicing variants for which the distance from exon will be 
						computed and added in the last column of this set.
					
					- when it finds HGMD values calls get_rearranged_hgmd_fields that remove some and merges other
					
					- when it finds EXAC columns gets the max among all and adds the new column
					
					- when it finds the columns with the predictions it changes the predictions values using
						the function convert_prediction
						
					- when it finds the 	gnomad colum
 Returns : a new output file with one more column

=cut		
sub rearrange_rawtabular_out {
	my $cfg_hash = shift;
	my $file = shift;
	my $outFile = shift;#output2rearrange
	#my $fb_vars = shift;
	my $log_file = shift;
	
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	#Open input and output
	open (OUT_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");

	##Genomic information inizialization
	my $func_column_ind = -1;
	my $detail_col_ind = -1;
	my $num_row = 1;
	my @exonic_cols = split(",",$cfg_hash->{'exonic_cols'});#descriptions for exonic variants into the func column
	my @splicing_cols = split(",",$cfg_hash->{'splicing_cols'});#descriptions for splicing variants into the func column
	my @func_cols = split(",",$cfg_hash->{'func_cols'});#columns from the output table describing the variant function
	my @func_cols_inds = ();#indexes found for func columns
	my @genename_cols_inds = ();#indexes found for genename columns
	my $num_refseq_fields = $cfg_hash->{'num_refseq_fields'};

	#Genomic annotation fields to use (RefSeq, Gencode, etc)
	my $genomic_ann_fields = $cfg_hash->{'genomic_ann_fields'};
	#transform the fields from annovar into all lowercase and without symbols
	$genomic_ann_fields =~ s/[\-\.\+]/\_/g;
	$genomic_ann_fields = lc($genomic_ann_fields);	
	my $func_col = $cfg_hash->{'func_col'};#describing the variant function column field
	my @genomic_ann_fields = split(",",$genomic_ann_fields);
	my @compl_gen_ann_fields = ();#Concatenate func_col with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_gen_ann_fields,$func_col."_".$genomic_ann_field);
	}
	my $gene_col = $cfg_hash->{'gene_col'};#describing the variant function column field
	my @compl_genename_fields = ();#Concatenate gene_col with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_genename_fields,$gene_col."_".$genomic_ann_field);
	}
		
	##FILTER column
	#my $filter_col = $cfg_hash->{'db_vcf_filter_str'};
	#my $filter_col_ind = -1;
	#PREDICTIONS
	my @pred_2_transf_cols = split(",",$cfg_hash->{'pred_to_transform'});
	my @pred_2_transf_cols_inds = ();
	#HGMD
	my $hgmd_start_col = $cfg_hash->{'hgmd_start_col'};
	my $hgmd_start_col_ind = -1;
	my @hgmd_cols_to_keep = split(",",$cfg_hash->{'hgmd_cols_to_keep'});
	my @hgmd_cols_to_group = split(",",$cfg_hash->{'hgmd_cols_to_group'});
	my $hgmd_group_col = $cfg_hash->{'hgmd_group_col'};
	#EXAC
	my $exac_all_col = $cfg_hash->{'exac_all_col'};
	my $exac_all_col_ind = -1;	
	my $exac_start_col = $cfg_hash->{'exac_start_col'};
	my $exac_start_col_ind = -1;	
	my $exac_nat_cols = $cfg_hash->{'exac_nat_cols'};	
	my $exac_empty_val = $cfg_hash->{'exac_empty_val'};

	#1000 genomes 
	my $_1000g_all_col = $cfg_hash->{'1000g_all_col'};
	my $_1000g_all_col_ind = -1;	
	my $_1000g_start_col = $cfg_hash->{'1000g_start_col'};
	my $_1000g_start_col_ind = -1;	
	my $_1000g_nat_cols = $cfg_hash->{'1000g_nat_cols'};		
	my $_1000g_empty_val = $cfg_hash->{'1000g_empty_val'};
	#gnomad exome 
	my $gnomad_ex_all_col = $cfg_hash->{'gnomad_ex_all_col'};
	my $gnomad_ex_all_col_ind = -1;	
	my $gnomad_ex_start_col = $cfg_hash->{'gnomad_ex_start_col'};
	my $gnomad_ex_start_col_ind = -1;	
	my $gnomad_ex_nat_cols = $cfg_hash->{'gnomad_ex_nat_cols'};		
	my $gnomad_ex_empty_val = $cfg_hash->{'gnomad_ex_empty_val'};
	#gnomad genome
	my $gnomad_gen_all_col = $cfg_hash->{'gnomad_gen_all_col'};
	my $gnomad_gen_all_col_ind = -1;	
	my $gnomad_gen_start_col = $cfg_hash->{'gnomad_gen_start_col'};
	my $gnomad_gen_start_col_ind = -1;	
	my $gnomad_gen_nat_cols = $cfg_hash->{'gnomad_gen_nat_cols'};	
	my $gnomad_gen_empty_val = $cfg_hash->{'gnomad_gen_empty_val'};

	#Get the disease gene association hash
	my $gene_hash;
	get_disease_gene_association($cfg_hash,\$gene_hash,$log_file);
	
	###########
	#Reads the first version of the VarGenius output line per line		
	my @header_fields = ();
	my $semicolon = ":";
	#Write a new file with the given columns only
	while ( my $row = <DATA> ) {
		my $new_row = "";
		chomp($row);
		my @fields = split /\t/, $row;
		
		#############################
		#HEADER: the first set of instruction si to go in the header and 
		#get the two indexes of the columns needed
		if ( $num_row == 1){
			#my $col_num = 0;
			#########################################
			#Print the new header jumping the func details columns
			#########################################
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];
				
				if ( grep {/\b$field\b/} @compl_genename_fields ){
					#Before the genename we put the gene panels where it is present
					#The field name is obtained by removing $gene_col from the column name
					my $panelcolname = $field;
					$panelcolname =~ s/$gene_col/panels/;
					$new_row .= "$panelcolname\t";
					push(@genename_cols_inds,$col);
				}
				#Put the index of the func col in an array and construct the new header:
				#gene_XX,func_XX,genedetail_XX,exonicfunc_XX,XX_symbol,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange
				#becomes:
				#gene_XX,func_XX;exonicfunc_XX,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange,genedetail_XX,splvar_dist_from_exon
				#Because exonic_func is absorbed into func and XX_symbol will be removed since there is already gene_XX
				if ( grep {/\b$field\b/} @compl_gen_ann_fields ) {
					$new_row .= "$field\t";
					#Add the index where to find the func cols
					push(@func_cols_inds,$col);					
					#Print both the header for the detail and for distance
					my $detail_column = $fields[$col+1];
					my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					
					#Get the genomic annotation name removing the aachange field
					my $gen_ann_name = $field;
					$gen_ann_name =~ s/^$func_col\_//;		
					$new_row .= $gen_ann_cols."\t".$detail_column."\t$gen_ann_name\_splvar_dist_from_exon\t";
					$col = $col + $num_refseq_fields;
				}		
				#Generate the new header for HGMD info and get the column index for later
				elsif ( $field eq $hgmd_start_col){
						#Get the index of the info col
						$hgmd_start_col_ind = $col;
						$new_row .= join("_hgmd\t",@hgmd_cols_to_keep)."_hgmd\t";
						$new_row .= $hgmd_group_col."\t";
						$col++;
				}
				#get the column index for EXAC_All
				elsif ( $field eq $exac_all_col){
						#Get the index of the info col
						$exac_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for EXAC info and get the column index for later
				elsif ( $field eq $exac_start_col){
						#Get the index of the info col
						$exac_start_col_ind = $col;
						$new_row .= "exac_max\t";
						$col = $col + $exac_nat_cols;
				}
				#get the column index for gnomad_exome_all
				elsif ( $field eq $gnomad_ex_all_col){
						#Get the index of the info col
						$gnomad_ex_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $gnomad_ex_start_col){
						#Get the index of the info col
						$gnomad_ex_start_col_ind = $col;
						$new_row .= "gnomad_ex_max\t";
						$col = $col + $gnomad_ex_nat_cols;
				}
				#get the column index for gnomad_genome_all
				elsif ( $field eq $gnomad_gen_all_col){
						#Get the index of the info col
						$gnomad_gen_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $gnomad_gen_start_col){
						#Get the index of the info col
						$gnomad_gen_start_col_ind = $col;
						$new_row .= "gnomad_gen_max\t";
						$col = $col + $gnomad_gen_nat_cols;
				}
				#get the column index for gnomad_genome_all
				elsif ( $field eq $_1000g_all_col){
						#Get the index of the info col
						$_1000g_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $_1000g_start_col){
						#Get the index of the info col
						$_1000g_start_col_ind = $col;
						$new_row .= "db1000g_max\t";
						$col = $col + $_1000g_nat_cols;
				}
				#elsif ( $field eq $filter_col){
						#$new_row .= $field."\tvc_prog\t";
						#$filter_col_ind = $col;
				#}
											
				#All other fields
				else{
					$new_row .= $field."\t";
				}
			}
			################################
			#Get the columns names and select columns for predictions
			################################
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];	
				#Put all the field names in an array
				push(@header_fields,$field);
				#Store the indexes of the predictions to transform
				if ( grep {/\b$field\b/} @pred_2_transf_cols ) {
					#print "Putting the col $field ($col)in the array..\n";
					push(@pred_2_transf_cols_inds,$col);
				}
			}
						
		}
		#########################################
		#DATA: The second set of instructions is for when we are in the data
		#Hence, we write the output file modifying the information located
		#in those columns we obtained in the HEADER part
		else{
			
			#Insert now the new field between the previous two
			for ( my $col=0; $col< scalar(@fields); $col++){
				my $field = $fields[$col];
				my $refseq_info = "";
				
				################################################
				##########GET THE GENE PANELS FOR THE GENE######
				################################################
				if ( grep {/\b$col\b/} @genename_cols_inds ) {
					my $genename = $field;
					if ( defined $gene_hash->{$genename} ){
						$new_row .= $gene_hash->{$genename}."\t$genename\t";
					}else{
						$new_row .= "$vargenius_empty_val\t$genename\t";
					}
				}
				#################################################
				######REARRANGE FUNCTION DETAILS COLUMNS#########
				#################################################
				elsif ( grep {/\b$col\b/} @func_cols_inds ) {
					my $func_column = $field;
					my $detail_column = $fields[$col+1];
					my $exonfunc_column = $fields[$col+2];
					my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					my $aachange_col = $fields[$col+8];
					#$new_row .= $detail_column."\t";
					######GET THE SPLICING VARIANT DISTANCE#########
					if ( (grep {/\b$func_column\b/} @splicing_cols) and ( $detail_column ne $vargenius_empty_val) ){
						#Get the distance from the first result 
						my $distance = "-";
						#NM_001269053:exon7:c.505-27A>C]---[NM_138812:exon8:c.505-27A>C
						my @results = split(/\Q$mult_ann_sep\E/,$detail_column);
						#NM_001269053:exon7:c.505-27A>C
						my $first_res = shift @results;
						my @parts = split(/\Q$semicolon\E/,$first_res);
						#c.505-27A>C
						my $last = $parts[2];
						#Remove the c.505 => -27A>C
						$last =~  s/c\.\d+//g;
						#Get only the distance -27A>C => -27 (if exists)
						if ( $last =~ /[+-]\d+/ ){
							$last =~ /([+-]\d+)/;				
							$distance = $1;			
						}
						$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$distance."\t";
					}#else{
						#$new_row .= $vargenius_empty_val."\t";
					#}		
					######GROUP INFORMATION FOR EXONIC FUNCTION#########
					elsif ( (grep {/\b$func_column\b/} @exonic_cols))  {
						if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$aachange_col."\t$vargenius_empty_val\t";
					}else{
						if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$vargenius_empty_val."\t";
					}				
					#jump all the refseq cols
					$col = $col + $num_refseq_fields;
					$new_row .= $refseq_info;
				}
				########################################
				######CHANGE PREDICTION SYMBOLS#########
				########################################
				elsif (grep {/\b$col\b/} @pred_2_transf_cols_inds) {
						#my $i = 0;
						#foreach my $h (@header_fields){
							#print "$i $h\n";
							#$i++;
						#}
						
						#print "\n\nselected:\n";
						#print join("\n", @pred_2_transf_cols_inds);
						if ( $field ne $vargenius_empty_val){

							#print "Converting prediction at row $num_row col $col (".$header_fields[$col]."), field $field\n";#DEBUGCODE
							$new_row .= convert_prediction($cfg_hash,$header_fields[$col],$field)."\t";
						}else{
							$new_row .= $vargenius_empty_val."\t";
						}
				}
				######################################## 
				######REARRANGE HGMD INFORMATION########
				########################################
				elsif ($col == $hgmd_start_col_ind) {
					$col++;
					my $info_field = $fields[$col];
					$new_row .= get_rearranged_hgmd_fields($info_field,\@hgmd_cols_to_keep,\@hgmd_cols_to_group)."\t";
				}
				######################################## 
				######REARRANGE EXAC NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $exac_all_col_ind) {
					
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$exac_empty_val\t";
					}
				}
				#######################################
				###EXAC: put only the maximum frequency among nationalities and ExAC_all
				########################################
				elsif ($col == $exac_start_col_ind) {

					#print "Rearranging from $col to ".($col+$exac_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @exac_fields = @fields[($col)..($col+$exac_nat_cols)];
					$new_row .= get_max_val(\@exac_fields , $vargenius_empty_val,$exac_empty_val)."\t";
					$col = $col + $exac_nat_cols;
				}
				######################################## 
				######REARRANGE gnomad exome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $gnomad_ex_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$gnomad_ex_empty_val\t";
					}
				}
				#######################################
				###gnomad exome: put only the maximum frequency among nationalities and gnomad_exome_all
				########################################
				elsif ($col == $gnomad_ex_start_col_ind) {
					#print "Rearranging from $col to ".($col+$gnomad_ex_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @gnomad_ex_fields = @fields[($col)..($col+$gnomad_ex_nat_cols)];
					$new_row .= get_max_val(\@gnomad_ex_fields,$vargenius_empty_val,$gnomad_ex_empty_val)."\t";
					$col = $col + $gnomad_ex_nat_cols;
				}
				######################################## 
				######REARRANGE gnomad genome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $gnomad_gen_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$gnomad_gen_empty_val\t";
					}
				}
				#######################################
				###gnomad genome: put only the maximum frequency among nationalities and gnomad_genome_all
				########################################
				elsif ($col == $gnomad_gen_start_col_ind) {
					#print "Rearranging from $col to ".($col+$gnomad_gen_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @gnomad_gen_fields = @fields[($col)..($col+$gnomad_gen_nat_cols)];
					$new_row .= get_max_val(\@gnomad_gen_fields,$vargenius_empty_val,$gnomad_gen_empty_val)."\t";
					$col = $col + $gnomad_gen_nat_cols;
				}
				######################################## 
				######REARRANGE 1000 genome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $_1000g_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$_1000g_empty_val\t";
					}
				}
				#######################################
				###1000 genome: put only the maximum frequency among nationalities and 1000 genome_all
				########################################
				elsif ($col == $_1000g_start_col_ind) {
					#print "Rearranging from $col to ".($col+$_1000g_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @_1000g_fields = @fields[($col)..($col+$_1000g_nat_cols)];
					$new_row .= get_max_val(\@_1000g_fields,$vargenius_empty_val,$_1000g_empty_val)."\t";
					$col = $col + $_1000g_nat_cols;
				}	
				##########################
				######QUESTA PARTE MI AGGIUNGEVA UN FLAG SE LA VARIANTE ERA DI FREEBAYES ORA L'
				##The variant caller
				#elsif ($col == $filter_col_ind) {
					##Get the compid
					#my $compid = $fields[0]."_".$fields[1]."_".$fields[3]."_".$fields[4];
					##print_and_log( "Checking if $compid is present into fb_vars (".$$fb_vars[0].")\n",$log_file);#DEBUGCODE
					#my @vars = @$fb_vars;
					#if ( grep {/\b$compid\b/} @vars ){ 
						#$new_row .= $field."\tFB\t";					
					#}else{
						#$new_row .= $field."\t$vargenius_empty_val\t";	
					#}
				#}						
				########################
				#Write all other needed columns as they are
				else{
					$new_row .= $field."\t";
				}	
			}				
		}
		chop($new_row);
		print OUT_FILE "$new_row\n";
		$num_row++;
	}

	close(DATA);
	close(OUT_FILE);	
}



=head2 rearrange_rawtabular_outOLD

 Title  : rearrange_rawtabular_out
 Usage  : rearrange_rawtabular_out( - filePath => 'the path to the file');
																			func_column => the column with the type of variant
																			detail_col => the column with the detail of the variant
 
 Function: 				
					- when it finds the gene symbol knows that after that there is the function details
						that will be rearranged by dividing the first result in different columns 
						In some cases there will be splicing variants for which the distance from exon will be 
						computed and added in the last column of this set.
					
					- when it finds HGMD values calls get_rearranged_hgmd_fields that remove some and merges other
					
					- when it finds EXAC columns gets the max among all and adds the new column
					
					- when it finds the columns with the predictions it changes the predictions values using
						the function convert_prediction
						
					- when it finds the 	gnomad colum
 Returns : a new output file with one more column

=cut		
sub rearrange_rawtabular_outOLD {
	my $cfg_hash = shift;
	my $file = shift;
	my $outFile = shift;#output2rearrange
	my $fb_vars = shift;
	my $log_file = shift;
	
	my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	#Open input and output
	open (OUT_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");

	##Genomic information inizialization
	my $func_column_ind = -1;
	my $detail_col_ind = -1;
	my $num_row = 1;
	my @exonic_cols = split(",",$cfg_hash->{'exonic_cols'});#descriptions for exonic variants into the func column
	my @splicing_cols = split(",",$cfg_hash->{'splicing_cols'});#descriptions for splicing variants into the func column
	my @func_cols = split(",",$cfg_hash->{'func_cols'});#columns from the output table describing the variant function
	my @func_cols_inds = ();#indexes found for func columns
	my @genename_cols_inds = ();#indexes found for genename columns
	my $num_refseq_fields = $cfg_hash->{'num_refseq_fields'};

	#Genomic annotation fields to use (RefSeq, Gencode, etc)
	my $genomic_ann_fields = $cfg_hash->{'genomic_ann_fields'};
	#transform the fields from annovar into all lowercase and without symbols
	$genomic_ann_fields =~ s/[\-\.\+]/\_/g;
	$genomic_ann_fields = lc($genomic_ann_fields);	
	my $func_col = $cfg_hash->{'func_col'};#describing the variant function column field
	my @genomic_ann_fields = split(",",$genomic_ann_fields);
	my @compl_gen_ann_fields = ();#Concatenate func_col with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_gen_ann_fields,$func_col."_".$genomic_ann_field);
	}
	my $gene_col = $cfg_hash->{'gene_col'};#describing the variant function column field
	my @compl_genename_fields = ();#Concatenate gene_col with RefSeq and wgEncodeGencode to get the exact fields name
	foreach my $genomic_ann_field (@genomic_ann_fields){
		push(@compl_genename_fields,$gene_col."_".$genomic_ann_field);
	}
		
	#FILTER column
	my $filter_col = $cfg_hash->{'db_vcf_filter_str'};
	my $filter_col_ind = -1;
	#PREDICTIONS
	my @pred_2_transf_cols = split(",",$cfg_hash->{'pred_to_transform'});
	my @pred_2_transf_cols_inds = ();
	#HGMD
	my $hgmd_start_col = $cfg_hash->{'hgmd_start_col'};
	my $hgmd_start_col_ind = -1;
	my @hgmd_cols_to_keep = split(",",$cfg_hash->{'hgmd_cols_to_keep'});
	my @hgmd_cols_to_group = split(",",$cfg_hash->{'hgmd_cols_to_group'});
	my $hgmd_group_col = $cfg_hash->{'hgmd_group_col'};
	#EXAC
	my $exac_all_col = $cfg_hash->{'exac_all_col'};
	my $exac_all_col_ind = -1;	
	my $exac_start_col = $cfg_hash->{'exac_start_col'};
	my $exac_start_col_ind = -1;	
	my $exac_nat_cols = $cfg_hash->{'exac_nat_cols'};	
	my $exac_empty_val = $cfg_hash->{'exac_empty_val'};

	#1000 genomes 
	my $_1000g_all_col = $cfg_hash->{'1000g_all_col'};
	my $_1000g_all_col_ind = -1;	
	my $_1000g_start_col = $cfg_hash->{'1000g_start_col'};
	my $_1000g_start_col_ind = -1;	
	my $_1000g_nat_cols = $cfg_hash->{'1000g_nat_cols'};		
	my $_1000g_empty_val = $cfg_hash->{'1000g_empty_val'};
	#gnomad exome 
	my $gnomad_ex_all_col = $cfg_hash->{'gnomad_ex_all_col'};
	my $gnomad_ex_all_col_ind = -1;	
	my $gnomad_ex_start_col = $cfg_hash->{'gnomad_ex_start_col'};
	my $gnomad_ex_start_col_ind = -1;	
	my $gnomad_ex_nat_cols = $cfg_hash->{'gnomad_ex_nat_cols'};		
	my $gnomad_ex_empty_val = $cfg_hash->{'gnomad_ex_empty_val'};
	#gnomad genome
	my $gnomad_gen_all_col = $cfg_hash->{'gnomad_gen_all_col'};
	my $gnomad_gen_all_col_ind = -1;	
	my $gnomad_gen_start_col = $cfg_hash->{'gnomad_gen_start_col'};
	my $gnomad_gen_start_col_ind = -1;	
	my $gnomad_gen_nat_cols = $cfg_hash->{'gnomad_gen_nat_cols'};	
	my $gnomad_gen_empty_val = $cfg_hash->{'gnomad_gen_empty_val'};

	#Get the disease gene association hash
	my $gene_hash;
	get_disease_gene_association($cfg_hash,\$gene_hash,$log_file);
	
	###########
	#Reads the first version of the VarGenius output line per line		
	my @header_fields = ();
	my $semicolon = ":";
	#Write a new file with the given columns only
	while ( my $row = <DATA> ) {
		my $new_row = "";
		chomp($row);
		my @fields = split /\t/, $row;
		
		#############################
		#HEADER: the first set of instruction si to go in the header and 
		#get the two indexes of the columns needed
		if ( $num_row == 1){
			#my $col_num = 0;
			#########################################
			#Print the new header jumping the func details columns
			#########################################
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];
				
				if ( grep {/\b$field\b/} @compl_genename_fields ){
					#Before the genename we put the gene panels where it is present
					$new_row .= "panels\t";
					push(@genename_cols_inds,$col);
				}
				#Put the index of the func col in an array and construct the new header:
				#gene_XX,func_XX,genedetail_XX,exonicfunc_XX,XX_symbol,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange
				#becomes:
				#gene_XX,func_XX;exonicfunc_XX,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange,genedetail_XX,splvar_dist_from_exon
				#Because exonic_func is absorbed into func and XX_symbol will be removed since there is already gene_XX
				if ( grep {/\b$field\b/} @compl_gen_ann_fields ) {
					$new_row .= "$field\t";
					#Add the index where to find the func cols
					push(@func_cols_inds,$col);					
					#Print both the header for the detail and for distance
					my $detail_column = $fields[$col+1];
					my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					
					#Get the genomic annotation name removing the aachange field
					my $gen_ann_name = $field;
					$gen_ann_name =~ s/^$func_col\_//;		
					$new_row .= $gen_ann_cols."\t".$detail_column."\t$gen_ann_name\_splvar_dist_from_exon\t";
					$col = $col + $num_refseq_fields;
				}		
				#Generate the new header for HGMD info and get the column index for later
				elsif ( $field eq $hgmd_start_col){
						#Get the index of the info col
						$hgmd_start_col_ind = $col;
						$new_row .= join("_hgmd\t",@hgmd_cols_to_keep)."_hgmd\t";
						$new_row .= $hgmd_group_col."\t";
						$col++;
				}
				#get the column index for EXAC_All
				elsif ( $field eq $exac_all_col){
						#Get the index of the info col
						$exac_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for EXAC info and get the column index for later
				elsif ( $field eq $exac_start_col){
						#Get the index of the info col
						$exac_start_col_ind = $col;
						$new_row .= "exac_max\t";
						$col = $col + $exac_nat_cols;
				}
				#get the column index for gnomad_exome_all
				elsif ( $field eq $gnomad_ex_all_col){
						#Get the index of the info col
						$gnomad_ex_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $gnomad_ex_start_col){
						#Get the index of the info col
						$gnomad_ex_start_col_ind = $col;
						$new_row .= "gnomad_ex_max\t";
						$col = $col + $gnomad_ex_nat_cols;
				}
				#get the column index for gnomad_genome_all
				elsif ( $field eq $gnomad_gen_all_col){
						#Get the index of the info col
						$gnomad_gen_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $gnomad_gen_start_col){
						#Get the index of the info col
						$gnomad_gen_start_col_ind = $col;
						$new_row .= "gnomad_gen_max\t";
						$col = $col + $gnomad_gen_nat_cols;
				}
				#get the column index for gnomad_genome_all
				elsif ( $field eq $_1000g_all_col){
						#Get the index of the info col
						$_1000g_all_col_ind = $col;
						$new_row .= $field."\t";
						#$col++;
				}
				#Generate the new header for gnomad info and get the column index for later
				elsif ( $field eq $_1000g_start_col){
						#Get the index of the info col
						$_1000g_start_col_ind = $col;
						$new_row .= "db1000g_max\t";
						$col = $col + $_1000g_nat_cols;
				}
				elsif ( $field eq $filter_col){
						$new_row .= $field."\tvc_prog\t";
						$filter_col_ind = $col;
				}
											
				#All other fields
				else{
					$new_row .= $field."\t";
				}
			}
			################################
			#Get the columns names and select columns for predictions
			################################
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];	
				#Put all the field names in an array
				push(@header_fields,$field);
				#Store the indexes of the predictions to transform
				if ( grep {/\b$field\b/} @pred_2_transf_cols ) {
					#print "Putting the col $field ($col)in the array..\n";
					push(@pred_2_transf_cols_inds,$col);
				}
			}
						
		}
		#########################################
		#DATA: The second set of instructions is for when we are in the data
		#Hence, we write the output file modifying the information located
		#in those columns we obtained in the HEADER part
		else{
			
			#Insert now the new field between the previous two
			for ( my $col=0; $col< scalar(@fields); $col++){
				my $field = $fields[$col];
				my $refseq_info = "";
				
				################################################
				##########GET THE GENE PANELS FOR THE GENE######
				################################################
				if ( grep {/\b$col\b/} @genename_cols_inds ) {
					my $genename = $field;
					if ( defined $gene_hash->{$genename} ){
						$new_row .= $gene_hash->{$genename}."\t$genename\t";
					}else{
						$new_row .= "$vargenius_empty_val\t$genename\t";
					}
				}
				#################################################
				######REARRANGE FUNCTION DETAILS COLUMNS#########
				#################################################
				elsif ( grep {/\b$col\b/} @func_cols_inds ) {
					my $func_column = $field;
					my $detail_column = $fields[$col+1];
					my $exonfunc_column = $fields[$col+2];
					my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					my $aachange_col = $fields[$col+8];
					#$new_row .= $detail_column."\t";
					######GET THE SPLICING VARIANT DISTANCE#########
					if ( (grep {/\b$func_column\b/} @splicing_cols) and ( $detail_column ne $vargenius_empty_val) ){
						#Get the distance from the first result 
						my $distance = "-";
						#NM_001269053:exon7:c.505-27A>C]---[NM_138812:exon8:c.505-27A>C
						my @results = split(/\Q$mult_ann_sep\E/,$detail_column);
						#NM_001269053:exon7:c.505-27A>C
						my $first_res = shift @results;
						my @parts = split(/\Q$semicolon\E/,$first_res);
						#c.505-27A>C
						my $last = $parts[2];
						#Remove the c.505 => -27A>C
						$last =~  s/c\.\d+//g;
						#Get only the distance -27A>C => -27 (if exists)
						if ( $last =~ /[+-]\d+/ ){
							$last =~ /([+-]\d+)/;				
							$distance = $1;			
						}
						$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$distance."\t";
					}#else{
						#$new_row .= $vargenius_empty_val."\t";
					#}		
					######GROUP INFORMATION FOR EXONIC FUNCTION#########
					elsif ( (grep {/\b$func_column\b/} @exonic_cols))  {
						if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$aachange_col."\t$vargenius_empty_val\t";
					}else{
						if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$vargenius_empty_val."\t";
					}				
					#jump all the refseq cols
					$col = $col + $num_refseq_fields;
					$new_row .= $refseq_info;
				}
				########################################
				######CHANGE PREDICTION SYMBOLS#########
				########################################
				elsif (grep {/\b$col\b/} @pred_2_transf_cols_inds) {
						#my $i = 0;
						#foreach my $h (@header_fields){
							#print "$i $h\n";
							#$i++;
						#}
						
						#print "\n\nselected:\n";
						#print join("\n", @pred_2_transf_cols_inds);
						if ( $field ne $vargenius_empty_val){

							#print "Converting prediction at row $num_row col $col (".$header_fields[$col]."), field $field\n";#DEBUGCODE
							$new_row .= convert_prediction($cfg_hash,$header_fields[$col],$field)."\t";
						}else{
							$new_row .= $vargenius_empty_val."\t";
						}
				}
				######################################## 
				######REARRANGE HGMD INFORMATION########
				########################################
				elsif ($col == $hgmd_start_col_ind) {
					$col++;
					my $info_field = $fields[$col];
					$new_row .= get_rearranged_hgmd_fields($info_field,\@hgmd_cols_to_keep,\@hgmd_cols_to_group)."\t";
				}
				######################################## 
				######REARRANGE EXAC NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $exac_all_col_ind) {
					
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$exac_empty_val\t";
					}
				}
				#######################################
				###EXAC: put only the maximum frequency among nationalities and ExAC_all
				########################################
				elsif ($col == $exac_start_col_ind) {

					#print "Rearranging from $col to ".($col+$exac_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @exac_fields = @fields[($col)..($col+$exac_nat_cols)];
					$new_row .= get_max_val(\@exac_fields , $vargenius_empty_val,$exac_empty_val)."\t";
					$col = $col + $exac_nat_cols;
				}
				######################################## 
				######REARRANGE gnomad exome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $gnomad_ex_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$gnomad_ex_empty_val\t";
					}
				}
				#######################################
				###gnomad exome: put only the maximum frequency among nationalities and gnomad_exome_all
				########################################
				elsif ($col == $gnomad_ex_start_col_ind) {
					#print "Rearranging from $col to ".($col+$gnomad_ex_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @gnomad_ex_fields = @fields[($col)..($col+$gnomad_ex_nat_cols)];
					$new_row .= get_max_val(\@gnomad_ex_fields,$vargenius_empty_val,$gnomad_ex_empty_val)."\t";
					$col = $col + $gnomad_ex_nat_cols;
				}
				######################################## 
				######REARRANGE gnomad genome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $gnomad_gen_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$gnomad_gen_empty_val\t";
					}
				}
				#######################################
				###gnomad genome: put only the maximum frequency among nationalities and gnomad_genome_all
				########################################
				elsif ($col == $gnomad_gen_start_col_ind) {
					#print "Rearranging from $col to ".($col+$gnomad_gen_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @gnomad_gen_fields = @fields[($col)..($col+$gnomad_gen_nat_cols)];
					$new_row .= get_max_val(\@gnomad_gen_fields,$vargenius_empty_val,$gnomad_gen_empty_val)."\t";
					$col = $col + $gnomad_gen_nat_cols;
				}
				######################################## 
				######REARRANGE 1000 genome NO INFORMATION########
				########################################
				#Put the new empty value where you find the vargenius empty val
				elsif ($col == $_1000g_all_col_ind) {
					if ( $field ne $vargenius_empty_val){
						$new_row .= "$field\t";
					}else{
						$new_row .= "$_1000g_empty_val\t";
					}
				}
				#######################################
				###1000 genome: put only the maximum frequency among nationalities and 1000 genome_all
				########################################
				elsif ($col == $_1000g_start_col_ind) {
					#print "Rearranging from $col to ".($col+$_1000g_nat_cols)." $field..\n";
					#Get the maximum from all the columns for the different nationalities
					my @_1000g_fields = @fields[($col)..($col+$_1000g_nat_cols)];
					$new_row .= get_max_val(\@_1000g_fields,$vargenius_empty_val,$_1000g_empty_val)."\t";
					$col = $col + $_1000g_nat_cols;
				}	
				##########################
				#The variant caller
				elsif ($col == $filter_col_ind) {
					#Get the compid
					my $compid = $fields[0]."_".$fields[1]."_".$fields[3]."_".$fields[4];
					#print_and_log( "Checking if $compid is present into fb_vars (".$$fb_vars[0].")\n",$log_file);#DEBUGCODE
					my @vars = @$fb_vars;
					if ( grep {/\b$compid\b/} @vars ){ 
						$new_row .= $field."\tFB\t";					
					}else{
						$new_row .= $field."\t$vargenius_empty_val\t";	
					}
				}						
				########################
				#Write all other needed columns as they are
				else{
					$new_row .= $field."\t";
				}	
			}				
		}
		chop($new_row);
		print OUT_FILE "$new_row\n";
		$num_row++;
	}

	close(DATA);
	close(OUT_FILE);	
}


#=head2 rearrange_rawtabular_out

 #Title  : rearrange_rawtabular_out
 #Usage  : rearrange_rawtabular_out( - filePath => 'the path to the file');
																			#func_column => the column with the type of variant
																			#detail_col => the column with the detail of the variant
 
 #Function: 				
					#- when it finds the gene symbol knows that after that there is the function details
						#that will be rearranged by dividing the first result in different columns 
						#In some cases there will be splicing variants for which the distance from exon will be 
						#computed and added in the last column of this set.
					
					#- when it finds HGMD values calls get_rearranged_hgmd_fields that remove some and merges other
					
					#- when it finds EXAC columns gets the max among all and adds the new column
					
					#- when it finds the columns with the predictions it changes the predictions values using
						#the function convert_prediction
						
					#- when it finds the 	gnomad colum
 #Returns : a new output file with one more column

#=cut		
#sub rearrange_rawtabular_out {
	#my $cfg_hash = shift;
	#my $file = shift;
	#my $outFile = shift;#output2rearrange

	
	#my $mult_ann_sep = $cfg_hash->{'mult_ann_sep'};
	#my $vargenius_empty_val = $cfg_hash->{'vargenius_empty_val'};
	
	##Open input and output
	#open (OUT_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	#open (DATA, "<$file") or die ("Cannot open $file\n");

	###Genomic information inizialization
	#my $func_column_ind = -1;
	#my $detail_col_ind = -1;
	#my $num_row = 1;
	#my @exonic_cols = split(",",$cfg_hash->{'exonic_cols'});#descriptions for exonic variants into the func column
	#my @splicing_cols = split(",",$cfg_hash->{'splicing_cols'});#descriptions for splicing variants into the func column
	#my @func_cols = split(",",$cfg_hash->{'func_cols'});#columns from the output table describing the variant function
	#my @func_cols_inds = ();#indexes found for func columns
	#my $num_refseq_fields = $cfg_hash->{'num_refseq_fields'};

	##Genomic annotation fields to use (RefSeq, Gencode, etc)
	#my $genomic_ann_fields = $cfg_hash->{'genomic_ann_fields'};
	##transform the fields from annovar into all lowercase and without symbols
	#$genomic_ann_fields =~ s/[\-\.\+]/\_/g;
	#$genomic_ann_fields = lc($genomic_ann_fields);	
	#my $func_col = $cfg_hash->{'func_col'};#describing the variant function column field
	#my @genomic_ann_fields = split(",",$genomic_ann_fields);
	#my @compl_gen_ann_fields = ();#Concatenate func_col with RefSeq and wgEncodeGencode to get the exact fields name
	#foreach my $genomic_ann_field (@genomic_ann_fields){
		#push(@compl_gen_ann_fields,$func_col."_".$genomic_ann_field);
	#}
	
	##PREDICTIONS
	#my @pred_2_transf_cols = split(",",$cfg_hash->{'pred_to_transform'});
	#my @pred_2_transf_cols_inds = ();
	##HGMD
	#my $hgmd_start_col = $cfg_hash->{'hgmd_start_col'};
	#my $hgmd_start_col_ind = -1;
	#my @hgmd_cols_to_keep = split(",",$cfg_hash->{'hgmd_cols_to_keep'});
	#my @hgmd_cols_to_group = split(",",$cfg_hash->{'hgmd_cols_to_group'});
	#my $hgmd_group_col = $cfg_hash->{'hgmd_group_col'};
	##EXAC
	#my $exac_all_col = $cfg_hash->{'exac_all_col'};
	#my $exac_all_col_ind = -1;	
	#my $exac_start_col = $cfg_hash->{'exac_start_col'};
	#my $exac_start_col_ind = -1;	
	#my $exac_nat_cols = $cfg_hash->{'exac_nat_cols'};	
	#my $exac_empty_val = $cfg_hash->{'exac_empty_val'};

	##1000 genomes 
	#my $_1000g_all_col = $cfg_hash->{'1000g_all_col'};
	#my $_1000g_all_col_ind = -1;	
	#my $_1000g_start_col = $cfg_hash->{'1000g_start_col'};
	#my $_1000g_start_col_ind = -1;	
	#my $_1000g_nat_cols = $cfg_hash->{'1000g_nat_cols'};		
	#my $_1000g_empty_val = $cfg_hash->{'1000g_empty_val'};
	##gnomad exome 
	#my $gnomad_ex_all_col = $cfg_hash->{'gnomad_ex_all_col'};
	#my $gnomad_ex_all_col_ind = -1;	
	#my $gnomad_ex_start_col = $cfg_hash->{'gnomad_ex_start_col'};
	#my $gnomad_ex_start_col_ind = -1;	
	#my $gnomad_ex_nat_cols = $cfg_hash->{'gnomad_ex_nat_cols'};		
	#my $gnomad_ex_empty_val = $cfg_hash->{'gnomad_ex_empty_val'};
	##gnomad genome
	#my $gnomad_gen_all_col = $cfg_hash->{'gnomad_gen_all_col'};
	#my $gnomad_gen_all_col_ind = -1;	
	#my $gnomad_gen_start_col = $cfg_hash->{'gnomad_gen_start_col'};
	#my $gnomad_gen_start_col_ind = -1;	
	#my $gnomad_gen_nat_cols = $cfg_hash->{'gnomad_gen_nat_cols'};	
	#my $gnomad_gen_empty_val = $cfg_hash->{'gnomad_gen_empty_val'};
	
	
	############
	##Reads the first version of the VarGenius output line per line		
	#my @header_fields = ();
	#my $semicolon = ":";
	##Write a new file with the given columns only
	#while ( my $row = <DATA> ) {
		#my $new_row = "";
		#chomp($row);
		#my @fields = split /\t/, $row;
		
		##############################
		##HEADER: the first set of instruction si to go in the header and 
		##get the two indexes of the columns needed
		#if ( $num_row == 1){
			##my $col_num = 0;
			##########################################
			##Print the new header jumping the func details columns
			##########################################
			#for ( my $col=0; $col< scalar(@fields); $col++){
				##Get the field name
				#my $field = $fields[$col];
				
				##Put the index of the func col in an array and construct the new header:
				##gene_XX,func_XX,genedetail_XX,exonicfunc_XX,XX_symbol,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange
				##becomes:
				##gene_XX,func_XX;exonicfunc_XX,XX_refgene,XX_exon_number,XX_nucl_change,XX_aachange,genedetail_XX,splvar_dist_from_exon
				##Because exonic_func is absorbed into func and XX_symbol will be removed since there is already gene_XX
				#if ( grep {/\b$field\b/} @compl_gen_ann_fields ) {
					#$new_row .= $field."\t";
					##Add the index where to find the func cols
					#push(@func_cols_inds,$col);					
					##Print both the header for the detail and for distance
					#my $detail_column = $fields[$col+1];
					#my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					
					##Get the genomic annotation name removing the aachange field
					#my $gen_ann_name = $field;
					#$gen_ann_name =~ s/^$func_col\_//;		
					#$new_row .= $gen_ann_cols."\t".$detail_column."\t$gen_ann_name\_splvar_dist_from_exon\t";
					#$col = $col + $num_refseq_fields;
				#}		
				##Generate the new header for HGMD info and get the column index for later
				#elsif ( $field eq $hgmd_start_col){
						##Get the index of the info col
						#$hgmd_start_col_ind = $col;
						#$new_row .= join("_hgmd\t",@hgmd_cols_to_keep)."_hgmd\t";
						#$new_row .= $hgmd_group_col."\t";
						#$col++;
				#}
				##get the column index for EXAC_All
				#elsif ( $field eq $exac_all_col){
						##Get the index of the info col
						#$exac_all_col_ind = $col;
						#$new_row .= $field."\t";
						##$col++;
				#}
				##Generate the new header for EXAC info and get the column index for later
				#elsif ( $field eq $exac_start_col){
						##Get the index of the info col
						#$exac_start_col_ind = $col;
						#$new_row .= "exac_max\t";
						#$col = $col + $exac_nat_cols;
				#}
				##get the column index for gnomad_exome_all
				#elsif ( $field eq $gnomad_ex_all_col){
						##Get the index of the info col
						#$gnomad_ex_all_col_ind = $col;
						#$new_row .= $field."\t";
						##$col++;
				#}
				##Generate the new header for gnomad info and get the column index for later
				#elsif ( $field eq $gnomad_ex_start_col){
						##Get the index of the info col
						#$gnomad_ex_start_col_ind = $col;
						#$new_row .= "gnomad_ex_max\t";
						#$col = $col + $gnomad_ex_nat_cols;
				#}
				##get the column index for gnomad_genome_all
				#elsif ( $field eq $gnomad_gen_all_col){
						##Get the index of the info col
						#$gnomad_gen_all_col_ind = $col;
						#$new_row .= $field."\t";
						##$col++;
				#}
				##Generate the new header for gnomad info and get the column index for later
				#elsif ( $field eq $gnomad_gen_start_col){
						##Get the index of the info col
						#$gnomad_gen_start_col_ind = $col;
						#$new_row .= "gnomad_gen_max\t";
						#$col = $col + $gnomad_gen_nat_cols;
				#}
				##get the column index for gnomad_genome_all
				#elsif ( $field eq $_1000g_all_col){
						##Get the index of the info col
						#$_1000g_all_col_ind = $col;
						#$new_row .= $field."\t";
						##$col++;
				#}
				##Generate the new header for gnomad info and get the column index for later
				#elsif ( $field eq $_1000g_start_col){
						##Get the index of the info col
						#$_1000g_start_col_ind = $col;
						#$new_row .= "db1000g_max\t";
						#$col = $col + $_1000g_nat_cols;
				#}							
				##All other fields
				#else{
					#$new_row .= $field."\t";
				#}
			#}
			#################################
			##Get the columns names and select columns for predictions
			#################################
			#for ( my $col=0; $col< scalar(@fields); $col++){
				##Get the field name
				#my $field = $fields[$col];	
				##Put all the field names in an array
				#push(@header_fields,$field);
				##Store the indexes of the predictions to transform
				#if ( grep {/\b$field\b/} @pred_2_transf_cols ) {
					##print "Putting the col $field ($col)in the array..\n";
					#push(@pred_2_transf_cols_inds,$col);
				#}
			#}
						
		#}
		##########################################
		##DATA: The second set of instructions is for when we are in the data
		##Hence, we write the output file modifying the information located
		##in those columns we obtained in the HEADER part
		#else{
		
			##Insert now the new field between the previous two
			#for ( my $col=0; $col< scalar(@fields); $col++){
				#my $field = $fields[$col];
				#my $refseq_info = "";
				
				##################################################
				#######REARRANGE FUNCTION DETAILS COLUMNS#########
				##################################################
				#if ( grep {/\b$col\b/} @func_cols_inds ) {
					#my $func_column = $field;
					#my $detail_column = $fields[$col+1];
					#my $exonfunc_column = $fields[$col+2];
					#my $gen_ann_cols = join("\t",@fields[($col+4)..($col+7)]);
					#my $aachange_col = $fields[$col+8];
					##$new_row .= $detail_column."\t";
					#######GET THE SPLICING VARIANT DISTANCE#########
					#if ( (grep {/\b$func_column\b/} @splicing_cols) and ( $detail_column ne $vargenius_empty_val) ){
						##Get the distance from the first result 
						#my $distance = "-";
						##NM_001269053:exon7:c.505-27A>C]---[NM_138812:exon8:c.505-27A>C
						#my @results = split(/\Q$mult_ann_sep\E/,$detail_column);
						##NM_001269053:exon7:c.505-27A>C
						#my $first_res = shift @results;
						#my @parts = split(/\Q$semicolon\E/,$first_res);
						##c.505-27A>C
						#my $last = $parts[2];
						##Remove the c.505 => -27A>C
						#$last =~  s/c\.\d+//g;
						##Get only the distance -27A>C => -27 (if exists)
						#if ( $last =~ /[+-]\d+/ ){
							#$last =~ /([+-]\d+)/;				
							#$distance = $1;			
						#}
						#$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$distance."\t";
					#}#else{
						##$new_row .= $vargenius_empty_val."\t";
					##}		
					#######GROUP INFORMATION FOR EXONIC FUNCTION#########
					#elsif ( (grep {/\b$func_column\b/} @exonic_cols))  {
						#if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							#$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$aachange_col."\t$vargenius_empty_val\t";
					#}else{
						#if ( $exonfunc_column ne $vargenius_empty_val){ $func_column .= ";".$exonfunc_column;}
							#$refseq_info .= $func_column."\t".$gen_ann_cols."\t".$detail_column."\t".$vargenius_empty_val."\t";
					#}				
					##jump all the refseq cols
					#$col = $col + $num_refseq_fields;
					#$new_row .= $refseq_info;
				#}
				#########################################
				#######CHANGE PREDICTION SYMBOLS#########
				#########################################
				#elsif (grep {/\b$col\b/} @pred_2_transf_cols_inds) {
						##my $i = 0;
						##foreach my $h (@header_fields){
							##print "$i $h\n";
							##$i++;
						##}
						
						##print "\n\nselected:\n";
						##print join("\n", @pred_2_transf_cols_inds);
						#if ( $field ne $vargenius_empty_val){

							##print "Converting prediction at row $num_row col $col (".$header_fields[$col]."), field $field\n";#DEBUGCODE
							#$new_row .= convert_prediction($cfg_hash,$header_fields[$col],$field)."\t";
						#}else{
							#$new_row .= $vargenius_empty_val."\t";
						#}
				#}
				######################################### 
				#######REARRANGE HGMD INFORMATION########
				#########################################
				#elsif ($col == $hgmd_start_col_ind) {
					#$col++;
					#my $info_field = $fields[$col];
					#$new_row .= get_rearranged_hgmd_fields($info_field,\@hgmd_cols_to_keep,\@hgmd_cols_to_group)."\t";
				#}
				######################################### 
				#######REARRANGE EXAC NO INFORMATION########
				#########################################
				##Put the new empty value where you find the vargenius empty val
				#elsif ($col == $exac_all_col_ind) {
					
					#if ( $field ne $vargenius_empty_val){
						#$new_row .= "$field\t";
					#}else{
						#$new_row .= "$exac_empty_val\t";
					#}
				#}
				########################################
				####EXAC: put only the maximum frequency among nationalities and ExAC_all
				#########################################
				#elsif ($col == $exac_start_col_ind) {

					##print "Rearranging from $col to ".($col+$exac_nat_cols)." $field..\n";
					##Get the maximum from all the columns for the different nationalities
					#my @exac_fields = @fields[($col)..($col+$exac_nat_cols)];
					#$new_row .= get_max_val(\@exac_fields , $vargenius_empty_val,$exac_empty_val)."\t";
					#$col = $col + $exac_nat_cols;
				#}
				######################################### 
				#######REARRANGE gnomad exome NO INFORMATION########
				#########################################
				##Put the new empty value where you find the vargenius empty val
				#elsif ($col == $gnomad_ex_all_col_ind) {
					#if ( $field ne $vargenius_empty_val){
						#$new_row .= "$field\t";
					#}else{
						#$new_row .= "$gnomad_ex_empty_val\t";
					#}
				#}
				########################################
				####gnomad exome: put only the maximum frequency among nationalities and gnomad_exome_all
				#########################################
				#elsif ($col == $gnomad_ex_start_col_ind) {
					##print "Rearranging from $col to ".($col+$gnomad_ex_nat_cols)." $field..\n";
					##Get the maximum from all the columns for the different nationalities
					#my @gnomad_ex_fields = @fields[($col)..($col+$gnomad_ex_nat_cols)];
					#$new_row .= get_max_val(\@gnomad_ex_fields,$vargenius_empty_val,$gnomad_ex_empty_val)."\t";
					#$col = $col + $gnomad_ex_nat_cols;
				#}
				######################################### 
				#######REARRANGE gnomad genome NO INFORMATION########
				#########################################
				##Put the new empty value where you find the vargenius empty val
				#elsif ($col == $gnomad_gen_all_col_ind) {
					#if ( $field ne $vargenius_empty_val){
						#$new_row .= "$field\t";
					#}else{
						#$new_row .= "$gnomad_gen_empty_val\t";
					#}
				#}
				########################################
				####gnomad genome: put only the maximum frequency among nationalities and gnomad_genome_all
				#########################################
				#elsif ($col == $gnomad_gen_start_col_ind) {
					##print "Rearranging from $col to ".($col+$gnomad_gen_nat_cols)." $field..\n";
					##Get the maximum from all the columns for the different nationalities
					#my @gnomad_gen_fields = @fields[($col)..($col+$gnomad_gen_nat_cols)];
					#$new_row .= get_max_val(\@gnomad_gen_fields,$vargenius_empty_val,$gnomad_gen_empty_val)."\t";
					#$col = $col + $gnomad_gen_nat_cols;
				#}
				######################################### 
				#######REARRANGE 1000 genome NO INFORMATION########
				#########################################
				##Put the new empty value where you find the vargenius empty val
				#elsif ($col == $_1000g_all_col_ind) {
					#if ( $field ne $vargenius_empty_val){
						#$new_row .= "$field\t";
					#}else{
						#$new_row .= "$_1000g_empty_val\t";
					#}
				#}
				########################################
				####1000 genome: put only the maximum frequency among nationalities and 1000 genome_all
				#########################################
				#elsif ($col == $_1000g_start_col_ind) {
					##print "Rearranging from $col to ".($col+$_1000g_nat_cols)." $field..\n";
					##Get the maximum from all the columns for the different nationalities
					#my @_1000g_fields = @fields[($col)..($col+$_1000g_nat_cols)];
					#$new_row .= get_max_val(\@_1000g_fields,$vargenius_empty_val,$_1000g_empty_val)."\t";
					#$col = $col + $_1000g_nat_cols;
				#}							
				#########################
				##Write all other needed columns as they are
				#else{
					#$new_row .= $field."\t";
				#}	
			#}				
		#}
		#chop($new_row);
		#print OUT_FILE "$new_row\n";
		#$num_row++;
	#}

	#close(DATA);
	#close(OUT_FILE);	
#}


=head2 print_variants_stats

 Title   : print_variants_stats
 Usage   : print_variants_stats(  config_file => the config hash
								);

 Function: Uses some functions to print stats in an HTML web page
 
 
 Returns :
 
=cut
sub print_variants_stats{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $vargenius_out = shift;
	my $task = shift;
	my $log_file = shift;
	my $stats_out = shift;
	
	my $R_stats_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_stats_script'};
	
	my $comma = ",";
	#Get samples names involved for this analysis
	#Getting the samples list associated to the groupid
  my $samples_h;#Samples hash to be reordered
	my @sort_samples = ();#The array with samples id reordered
	#Get the samples ids involved for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	#Print the header for information about the genotype for the sample
	my $sample_gen_fields1 = "";
	my $sample_gen_fields2 = "";
	#Get the kinship, to make the resorting
	foreach my $sample_id (keys %{$group_sam}){
		#Obtain the kinship from the database given the sample id
		my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
														$cfg_hash->{'db_sample_id'},$sample_id);	
		#Build an hash to reorder
		$samples_h->{$sample_id}->{'id'} = $sample_id;
		$samples_h->{$sample_id}->{'k'} = $kinship;
	}
	#Sort the samples as in samples_order parameter
	sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	my $sort_sam_names_str = "";
	
	my $sample_name ="";
	foreach my $sort_sample (@sort_samples){
			#Getting the sample name
			$sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sort_sample);	
			$sort_sam_names_str .= $sample_name.$comma;
	}
	#Remove last comma
	chop($sort_sam_names_str);
	
	#Variables needed into the R script	 
	#Name of the study
	my @sams = split("_", $sample_name);
	pop(@sams);
	my $analysis = join("_",@sams);
	#print_and_log(  "The study name is: $analysis (from $sample_name)\n",$log_file);
	#DbSNP version
	my $dbsnp_fld = $cfg_hash->{'dbSNPversion'};
	#CADD version
	my $cadd_phred = $cfg_hash->{'cadd_phred_version'}; #= "cadd13_phred"
	
	#Zygosity stats
	my $RLogPath = $cfg_hash->{$analysis_id.'_'.$task.'_log_fold'}."/".$cfg_hash->{'R_log_file'};
	 
  my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $vargenius_out $sort_sam_names_str $stats_out $analysis $dbsnp_fld $cadd_phred' $R_stats_script $RLogPath";
   
  print_and_log(  "The command is:\n $command\n",$log_file);
  try_exec_command($command) or R_die($command,$R_stats_script);
}

=head2 generate_igvsession_xml

 Title  : generate_igvsession_xml
 Usage  : generate_igvsession_xml( - filePath => 'the path to the file');
																			
 
 Function: Generates an XML file which contains the paths to VCF and the BAM files
					and that can be used in IGV
  
 Returns : a new output file with one more column

=cut		
sub generate_igvsession_xml {
		my $cfg_hash = shift;
		my $paths_hash = shift;
		my $xml_igv_session_out = shift;
		my $log_file = shift;
		
		#Needed variables
		my $covpage = $cfg_hash->{"outlist_cov_page"};
		my $outpage = $cfg_hash->{"outlist_output_page"};
		#Since January 2019 I am now using the recalibrated BAM output from BaseRecalibrator
		#before I used the aligned,sorted and eventually duplicate removed BAM (outlist_sortedbam_desc)
		my $bam_desc = $cfg_hash->{"outlist_recalibratedbam_desc"};
		my $vcf_desc = $cfg_hash->{"outlist_final_vcf"};
		my $build_ver = $cfg_hash->{"ann_build_ver"};
		my $host = $cfg_hash->{"html_host"};
		my $working_folder = $cfg_hash->{"work_fold"};
		
		print_and_log("Printing the IGV session file $xml_igv_session_out...\n",$log_file);
		#Open XML file		
		open (IGV,">$xml_igv_session_out") or die "ERROR[?]: Cannot open $xml_igv_session_out\n";
 		#Print the header of the XML giving indication of the genome
		print IGV '<Session genome="'.$build_ver.'" hasGeneTrack="true" hasSequenceTrack="true" version="8">'."\n";
		print IGV "\t".'<Resources>'."\n";
		
		#Print VCF path
		my $vcf_path = $paths_hash->{$outpage}->{'-'}->{'-'}->{$vcf_desc}->{'path'};
		$vcf_path =~ s/$working_folder//;
		my $vcf_link = $host.$vcf_path;
		print_and_log("Using vcf: $vcf_link ..\n",$log_file);
		print IGV "\t\t".'<Resource path="'.$vcf_link.'"/>'."\n";
		#Print BAM paths				
		while (my ($sample_name, $value) = each %{ $paths_hash->{$covpage} } ) {
			#print_and_log("Reading $covpage - $sample_name - $value - $bam_desc\n",$log_file);#DEBUGCODE
			my $bam_path = $paths_hash->{$covpage}->{$sample_name}->{'-'}->{$bam_desc}->{'path'};
			print_and_log("Bampath $bam_path\n",$log_file);#DEBUGCODE
			if ( defined $bam_path) {
				if ( -e $bam_path){
					$bam_path =~ s/$working_folder//;
					my $bam_link = $host.$bam_path;
					print_and_log("Using BAM: $bam_link ..\n",$log_file);
					print IGV "\t\t".'<Resource path="'.$bam_link.'"/>'."\n";				
				}else{
					print_and_log("ERROR during the generation of the IGV session file. File $bam_path does not exist..\n",$log_file);
				}				#code
			}else{
					print_and_log("ERROR Cannot find any BAM ($bam_desc) for $sample_name..\n",$log_file);
			}
			
		}
		
		print IGV "\t".'</Resources>'."\n";
		print IGV '</Session>'."\n";
		close(IGV);
		print_and_log("Finished printing  the IGV session file...\n",$log_file);
}

=head2 get_gene_annotation

 Title  : get_gene_annotation
 Usage  : get_gene_annotation( - filePath => 'the path to the file');
																			func_column => the column with the type of variant
																			detail_col => the column with the detail of the variant
 
 Function: Prints the gene annotation for the output of VarGenius given in input
  
 Returns : a new output file with one more column

=cut		
sub get_gene_annotation {
	my $cfg_hash = shift;
	my $file = shift;
	my $gene_annotation_out = shift;
	my $log_file = shift;
		
	my $gene_field = $cfg_hash->{"gene_field"};
	my $gene_field_ind = -1;
	#Naming the gene annotation file
	#my $gene_ann = $file.".".$cfg_hash->{"gene_annotation_ext"};
	open (DATA, "<$file") or die ("Cannot open $file\n");
	open (GENE_INFO, ">$gene_annotation_out") or die ("Cannot open $gene_annotation_out\n");
	
	my $num_row  = 1;
	my @genes_used = ();
	my $sep = ",";

	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
	#print_and_log( "Starting gene annotation using $file..\n ",$log_file);#DEBUGCODE

	#Print the header
	#Get the gene annotation fields to be used (e.g. db_hpo_ids,db_gdi_score,db_rvis_score,db_rvis_perc)
	#and print the header
	#N.B. The header depends by the gene annotation fields where I do not include the HPO desc (because 
	#I fetch fields of the genes_table. If you change something here change also into get_gene_info()
	my @gene_annotation_fields = split(",",$cfg_hash->{'gene_annotation_fields'});
	my $header = $cfg_hash->{'db_genes_name'}."\t";
	foreach my $gen_ann_field (@gene_annotation_fields){
		 $header .=  $cfg_hash->{$gen_ann_field}."\t";
		 #For HPO put also the desc
		 if ( $gen_ann_field eq 'db_hpo_ids'){
			$gen_ann_field =~ s/\_ids/\_name/; 
			$header .=  $cfg_hash->{$gen_ann_field}."\t"; 
		}
	}
	chop($header);#remove last tab
	print GENE_INFO "$header\n";
						
	#Write a new file with the given columns only
	while ( my $line = <DATA> ) {
		chomp($line);

		my @fields = split /\t/, $line;
		#Go in the header and get the two indexes of the columns needed
		if ( $num_row == 1){
			
			#print_and_log( "I'm in the header..\n ",$log_file);#DEBUGCODE
			#Search the gene name field number
			for ( my $col=0; $col< scalar(@fields); $col++){
				#Get the field name
				my $field = $fields[$col];
				
				#Save the index of the gene name
				if (  $field eq $gene_field ) {
					$gene_field_ind = $col;
					#print_and_log( "Found $gene_field at position $gene_field_ind..\n ",$log_file);#DEBUGCODE
				}
			}
		}else{
			#Now add the gene information in a different table. The gene symbol field may contain
			#more than one symbols (separated by semicolon) and we take the first
			my $gene_field_val = $fields[$gene_field_ind];
			#The gene field could have more gene names separated by ',' or ';'
			#Hence I change here all the ; in , because the ; are less in general
			$gene_field_val =~ s/;/,/g;
			my @gene_syms = split($sep,$gene_field_val);
			#Pick only the first gene name
			my $gene_sym = shift @gene_syms;
			
			#If this gene is not present in the table of GENE info
			if ( not (grep {/\b$gene_sym\b/} @genes_used) ) {
				#print_and_log( "$gene_sym has not been previuosly take, taking it..\n ",$log_file);#DEBUGCODE
				#Get the information for the Gene
				my $gene_info = get_gene_info($cfg_hash,$dbh,$log_file,$gene_sym);
				print GENE_INFO "$gene_sym\t$gene_info\n";
				push(@genes_used,$gene_sym);
			}#else{print_and_log( "The gene_sym  previuosly take $file..\n ",$log_file);#DEBUGCODE}
		}
		$num_row++;
	}
	#Disconnect db
  $dbh->disconnect(); 	
	
	close(DATA);
	close(GENE_INFO);
}


=head2 get_all_genes_annotation

 Title  : get_all_genes_annotation
 Usage  : get_all_genes_annotation( - filePath => 'the path to the file');
																			func_column => the column with the type of variant
																			detail_col => the column with the detail of the variant
 
 Function: Prints the gene annotation for all genes into the database
 
 Returns : a new output file with one more column

=cut		
sub get_all_genes_annotation {
	my $cfg_hash = shift;
	my $gene_annotation_out = shift;
	my $log_file = shift;
		

	#Naming the gene annotation file
	#my $gene_ann = $file.".".$cfg_hash->{"gene_annotation_ext"};
	open (GENE_INFO, ">$gene_annotation_out") or die ("Cannot open $gene_annotation_out\n");

	my $sep = ",";

	my $query = "SELECT ".$cfg_hash->{'db_genes_name'}." FROM  ".$cfg_hash->{'db_genes_table'}.";";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $genes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_genes_name'});
		
	print_and_log( "Starting all gene annotation using..\n ",$log_file);

	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
	#Print the header
	#Get the gene annotation fields to be used (e.g. db_hpo_ids,db_gdi_score,db_rvis_score,db_rvis_perc)
	#and print the header
	my @gene_annotation_fields = split(",",$cfg_hash->{'gene_annotation_fields'});
	my $header = $cfg_hash->{'db_genes_name'}."\t";
	foreach my $gen_ann_field (@gene_annotation_fields){
		 $header .=  $cfg_hash->{$gen_ann_field}."\t";
		 #For HPO put also the desc
		 if ( $gen_ann_field eq 'db_hpo_ids'){
			$gen_ann_field =~ s/\_ids/\_desc/; 
			$header .=  $cfg_hash->{$gen_ann_field}."\t"; 
		}
	}
	chop($header);#remove last tab
	print GENE_INFO "$header\n";		
	
	#Get all genes information 
	foreach my $gene (keys %{$genes}){
		#Get the information for the Gene
		my $gene_info = get_gene_info($cfg_hash,$dbh,$log_file,$gene);
		print GENE_INFO "$gene\t$gene_info\n";	
	}
	#Disconnect db
  $dbh->disconnect(); 	
		
	close(GENE_INFO);
}


=head2 get_max_val

 Title  : get_max_val
 Usage  : get_max_val( 
 
 Function: gets the max value from an array
					
  
 Returns : a number being the max

=cut		
sub get_max_val{
	my $array = shift;
	my $vargenius_empty_val = shift;
	my $exac_empty_val = shift;
	
	my $max = -100;
	my $empty_vals = 0;
	#print "Comparing ".scalar(@$array)." numbers\n";
	foreach my $el (@$array){
			if ( $el ne $vargenius_empty_val){
					if ( $el > $max){
						#print "$el > $max\n ";
						$max = $el;
					}
			}else{
					$empty_vals++;
			}
	}
	if ( $empty_vals == scalar(@$array)){
			$max =  $exac_empty_val;
	}
	return $max;
}

=head2 get_rearranged_hgmd_fields

 Title  : get_rearranged_hgmd_fields
 Usage  : get_rearranged_hgmd_fields( 
 
 Function: rearranges hgmd fields
					
  
 Returns : a string with rearranged fields

=cut		
sub get_rearranged_hgmd_fields{
	my $field = shift;
	my $hgmd_cols_to_keep = shift;
	my $hgmd_cols_to_group = shift;
	my $hgmd_info_sep = ";"; 
	my $hash_vals;
	my $vargenius_empty_val = "-";
	
	my $tab = "\t";
	my $out_field = "";
	if ( $field ne $vargenius_empty_val){
		#print "Rearranging $field\n";
		
		#Split the fields from HGMD info to get an hash with the values
		my @els = split($hgmd_info_sep,$field);
		foreach my $el (@els){
			my @pieces = split("=",$el);
			$hash_vals->{$pieces[0]} = $pieces[1]; 
		}
		#Two new variables to construct the output string
		my $diff_cols = "";
		my $grouped_cols = "";

		#Put in different columns the fields to keep in the order given
		foreach my $tokeep (@$hgmd_cols_to_keep){
			if (defined $hash_vals->{$tokeep}){
				$diff_cols .= $hash_vals->{$tokeep}.$tab;
			}else{
				$diff_cols .= $vargenius_empty_val.$tab;
			}
		}
		chop($diff_cols);
		#Put in the same column the fields to group in the order given
		if (scalar(@$hgmd_cols_to_group) > 0){
			foreach my $togroup (@$hgmd_cols_to_group){
				if (defined $hash_vals->{$togroup}){
					$grouped_cols .= $togroup."=".$hash_vals->{$togroup}.$hgmd_info_sep;
				}
			}
			chop($grouped_cols);
		}else{
			$grouped_cols .= $vargenius_empty_val;
		}
		
	$out_field = join($tab,$diff_cols,$grouped_cols);

	}else{
		#Put empty columns for each column to keep
		foreach my $tokeep (@$hgmd_cols_to_keep){
			$out_field .= $vargenius_empty_val.$tab;
		}
		#Insert one last column for the informative grouping field
		$out_field .= $vargenius_empty_val;
	}
		#print "Rearranged : $out_field\n";	
	return $out_field;
}

=head2 convert_prediction

 Title  : convert_prediction
 Usage  : convert_prediction( 
 
 Function: being one of the predictions changes the symbol according to my 
					
  
 Returns : a new symbol for the one given in input

=cut		
sub convert_prediction{
	my $cfg_hash = shift;
	my $score_name = shift;
	my $pred = shift;
	
	my $pred_hash;
	my $pred_t = $pred;
	
	#print "Getting hash of predictions from  ".$cfg_hash->{$score_name."_old"}." using ".$cfg_hash->{$score_name."_new"}."\n";
	if ( defined $cfg_hash->{$score_name."_old"} and defined $cfg_hash->{$score_name."_new"} ){
			get_preds_hash($cfg_hash->{$score_name."_old"},$cfg_hash->{$score_name."_new"},\$pred_hash);
	}
	
	
	my $chars = join "", keys %{$pred_hash};
	my %dr_pred_hash =  %$pred_hash;
	#print "Converting $pred_t using $chars\n";
	$pred_t  =~ s/([$chars])/$dr_pred_hash{$1}/g;
	
	#print "$pred becomes $pred_t\n";#DEBUGCODE
	return $pred_t;
}

=head2 get_preds_hash

 Title  : get_preds_hash
 Usage  : get_preds_hash( 
 
 Function: being one of the predictions changes the symbol according to my 
					
  
 Returns : a new symbol for the one given in input

=cut		
sub get_preds_hash{
	my $old = shift;
	my $new = shift;
	my ($preds) = shift;
	
	my @old_s = split(",",$old);
	my @new_s = split(",",$new);
	
	if (scalar(@old_s) == scalar(@new_s)){
		for (my $i = 0; $i< scalar(@old_s); $i++){
			$$preds->{$old_s[$i]} = $new_s[$i];
		}
	}
}

#####################SEND EMAIL

=head2 send_email

 Title   : send_email
 Usage   : send_email( -database => 'name of the database,
                               );

 Function: Sends an email to dests
 Returns : 

=cut			
sub send_email {
	my $from = shift;
	my $to = shift;
	my $subject = shift;
	my $message = shift;
	my $log_file= shift;
	
	#print_and_log( "Sending an email with sendmail...\n",$log_file);
	#Recipients must be separated by comma
	open(MAIL, "|/usr/sbin/sendmail -t");

	# Email Header
	print MAIL "To: $to\n";
	print MAIL "From: $from\n";
	print MAIL "Subject: $subject\n\n";
	# Email Body
	print MAIL $message;

	close(MAIL);
	print "Email Sent Successfully!\n";
	print_and_log( "Email Sent Successfully!\n",$log_file);
}



=head2 evaluate_quality_of_sequencing

 Title   : evaluate_quality_of_sequencing
 Usage   : evaluate_quality_of_sequencing(  config_file => the config hash
								);

 Function: Compares some statistics related with the samples and with the analysis
			with other sequencing
						
			Uses the userid and the targetbed file name to distinguish this analysis
			into a localized set.
					
						
						
 
 Returns :
 
=cut
sub evaluate_quality_of_sequencing{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $tables_list_file = shift;
	my $log_file = shift;

	#Obtain the userid from the database given the group id
	my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Obtain the target file name from the database given the group id
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Obtain the db_analysis_name from the database given the group id
	my $analysisname = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
														
	#############################
	#STATISTICS ON VARIANTS
	#############################								
	#Needed values
	my $minvars = 1000000000;
	my $maxvars = -1;
	my $meanvars = 0;
	
	my $mintrtv = 1000000000;
	my $maxtrtv = -1;
	my $meantrtv = 0;
	
	my $mindenovop = 1000000000;
	my $maxdenovop = -1;
	my $meandenovop = 0;
	
	my $tot_results = 0;
	
	#Get totvariants, trtvratio and percentage of denovo from similar analyses
	print_and_log( "Get totvariants, trtvratio and percentage of denovo from similar analyses\n",$log_file);
	my $query = "SELECT ".$cfg_hash->{'db_analstats_totvariants'}.",".$cfg_hash->{'db_analstats_tr_tvratio'}.
				",".$cfg_hash->{'db_analstats_denovoperc'}." FROM ".$cfg_hash->{'db_analysis_statistics'}."  WHERE ".
				$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}." FROM ".
					$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_userid'}." = $userid AND ".
					$cfg_hash->{'db_targetbed'}." = '$targetbed');";							
	print_and_log("Executing: $query\n",$log_file);
	#fetch_all_rows
	my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
								$cfg_hash->{'db_pass'},$query);
	
	print_and_log( "Getting average data for variants\n",$log_file);#DEBUGCODE
	foreach my $result (@$res) {
		my @arr = @$result;
		#Get the needed fields
		my $totvar = $arr[0];
		my $trtv = $arr[1];
		my $denovop = $arr[2];

		
		#Totvariants
		if ( $totvar ne $cfg_hash->{'vargenius_empty_val'} and $totvar > 0){
			if ( $totvar < $minvars) {$minvars = $totvar;}
			if ( $totvar > $maxvars) {$maxvars = $totvar;}
			$meanvars	= $meanvars + $totvar;			
		}
		
		#TrTV
		if ( $trtv ne $cfg_hash->{'vargenius_empty_val'} and $trtv > 0){		
			if ( $trtv < $mintrtv) {$mintrtv = $trtv;}
			if ( $trtv > $maxtrtv) {$maxtrtv = $trtv;}
			$meantrtv	= $meantrtv + $trtv;
		}
		
		#denovo
		if ( $denovop ne $cfg_hash->{'vargenius_empty_val'}  and $denovop > 0){	
			if ( $denovop < $mindenovop) {$mindenovop = $denovop;}
			if ( $denovop > $maxdenovop) {$maxdenovop = $denovop;}
			$meandenovop = $meandenovop + $denovop;
		}				
		$tot_results++;
		
	}
	
	#Obtain means
	$meanvars = $meanvars / $tot_results;
	$meantrtv = $meantrtv / $tot_results;
	$meandenovop = $meandenovop / $tot_results;
	
	my $anal_totvariants = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analysis_statistics'},$cfg_hash->{'db_analstats_totvariants'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $anal_tr_tvratio = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analysis_statistics'},$cfg_hash->{'db_analstats_tr_tvratio'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);	
	my $anal_denovoperc = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analysis_statistics'},$cfg_hash->{'db_analstats_denovoperc'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);	
	
	#Evaluate quality of the three parameters
	my $anal_totvariants_qual = "NA";
	if ( $anal_totvariants >= $minvars and $anal_totvariants <= $maxvars){$anal_totvariants_qual = "HIGH";}elsif ($anal_totvariants >= 0){$anal_totvariants = "LOW";}
	my $anal_tr_tvratio_qual = "NA";
	if ( $anal_tr_tvratio >= $mintrtv and $anal_tr_tvratio <= $maxtrtv){$anal_tr_tvratio_qual = "HIGH";}elsif ($anal_tr_tvratio >= 0){$anal_tr_tvratio = "LOW";}
	my $anal_denovoperc_qual = "NA";
	if ( $anal_denovoperc >= $mindenovop and $anal_denovoperc <= $maxdenovop){$anal_denovoperc_qual = "HIGH";}elsif ($anal_denovoperc >= 0 ){$anal_denovoperc = "LOW";}
			
	print_and_log( "Printing a table for variants\n",$log_file);#DEBUGCODE	
	#Finalization: print a table on a file that will be added to the list of tables to print in HTML
	my $varstats = extract_name($tables_list_file,"noext")."_VARSTATS.totab";
	open (VAR_ST,">$varstats") or die "ERROR: Cannot open $varstats . The program will exit..\n";
	print VAR_ST "-\tMIN\tMAX\tAVERAGE\t$analysisname\tQUALITY\n".
				$cfg_hash->{'db_analstats_totvariants'}."\t$minvars\t$maxvars\t$meanvars\t$anal_totvariants\t$anal_totvariants_qual\n".
				$cfg_hash->{'db_analstats_tr_tvratio'}."\t$mintrtv\t$maxtrtv\t$meantrtv\t$anal_tr_tvratio\t$anal_tr_tvratio_qual\n".
				$cfg_hash->{'db_analstats_denovoperc'}."\t$mindenovop\t$maxdenovop\t$meandenovop\t$anal_denovoperc\t$anal_denovoperc_qual\n";
	close(VAR_ST);

	#Add the output to the file list
	append_str_2_file_if_path_notexist($tables_list_file,$cfg_hash->{'VARSTATSCOMP_table'}."\t".$varstats);		
			
	#############################
	#STATISTICS ON ALIGNMENT
	#############################
	#Needed values
	my $minreads = 1000000000;
	my $maxreads = -1;
	my $meanreads = 0;
	
	my $minmapped = 1000000000;
	my $maxmapped = -1;
	my $meanmapped = 0;
	
	my $minduplicates = 1000000000;
	my $maxduplicates = -1;
	my $meanduplicates = 0;
	
	$tot_results = 0;
	
	#Get number of aligned and duplicate reads
	print_and_log( "Get number of aligned and duplicate reads\n",$log_file);#DEBUGCODE
	$query = "SELECT ".$cfg_hash->{'db_sampstats_totreads'}.",".$cfg_hash->{'db_sampstats_mappedreads'}.
				",".$cfg_hash->{'db_sampstats_duplicates'}." FROM ".$cfg_hash->{'db_sample_statistics'}."  WHERE ".
				$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}." FROM ".
					$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_userid'}." = $userid AND ".
					$cfg_hash->{'db_targetbed'}." = '$targetbed');";							
	print_and_log("Executing: $query\n",$log_file);
	#fetch_all_rows
	$res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
								$cfg_hash->{'db_pass'},$query);

	print_and_log( "Computing average data for samples \n",$log_file);#DEBUGCODE		
	foreach my $result (@$res) {
		my @arr = @$result;
		#Get the database id of the variant
		my $totreads = $arr[0];
		my $mapped = $arr[1];
		my $duplicates = $arr[2];
		
		#Totvariants
		if ( $totreads ne $cfg_hash->{'vargenius_empty_val'} and $totreads > 0){			
			if ( $totreads < $minreads) {$minreads = $totreads;}
			if ( $totreads > $maxreads) {$maxreads = $totreads;}
			$meanreads	= $meanreads + $totreads;
		}
		
		#TrTV
		if ( $mapped ne $cfg_hash->{'vargenius_empty_val'} and $mapped > 0){			
			if ( $mapped < $minmapped) {$minmapped = $mapped;}
			if ( $mapped > $maxmapped) {$maxmapped = $mapped;}
			$meanmapped	= $meanmapped + $mapped;
		}		
		
		#denovo
		if ( $duplicates ne $cfg_hash->{'vargenius_empty_val'} and $duplicates > 0){			
			if ( $duplicates < $minduplicates) {$minduplicates = $duplicates;}
			if ( $duplicates > $maxduplicates) {$maxduplicates = $duplicates;}
			$meanduplicates	= $meanduplicates + $duplicates;
		}						
		$tot_results++;
	}		

	#Obtain means
	if ( $meanreads > 0){
		$meanreads = $meanreads / $tot_results;
	}
	if ( $meanmapped > 0){
		$meanmapped = $meanmapped / $tot_results;
	}
	if ( $meanduplicates > 0){
		$meanduplicates = $meanduplicates / $tot_results;
	}		

	#Get the samples ids involved for the analysis
	$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $res_group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});


	#Finalization: print a table on a file that will be added to the list of tables to print in HTML
	my $alnstats = extract_name($tables_list_file,"noext")."_ALNSTATS.totab";
	open (ALN_ST,">$alnstats") or die "ERROR: Cannot open $alnstats. The program will exit..\n";
	
	my $head_str = "";
	my $totr_str = "";
	my $map_str = "";
	my $dup_str = "";
	
	foreach my $sample_id (keys %{$res_group_sam}){
		print_and_log( "Getiing data for sample $sample_id \n",$log_file);#DEBUGCODE		

		my $anal_totreads = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_statistics'},$cfg_hash->{'db_sampstats_totreads'},
														$cfg_hash->{'db_sample_id'},$sample_id);
		my $anal_mapped = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_statistics'},$cfg_hash->{'db_sampstats_mappedreads'},
															$cfg_hash->{'db_sample_id'},$sample_id);	
		my $anal_duplicates = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_statistics'},$cfg_hash->{'db_sampstats_duplicates'},
														$cfg_hash->{'db_sample_id'},$sample_id);			
	
		#Obtain the db_analysis_name from the database given the analysis id
		my $samplename = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
														$cfg_hash->{'db_sample_id'},$sample_id);
															
		print_and_log( "Compiling table for sample $sample_id \n",$log_file);#DEBUGCODE		

		#Evaluate quality of the values
		my $anal_totreads_qual = "NA";
		if ( $anal_totreads >= $minreads and $anal_totreads <= $maxreads){$anal_totreads_qual = "HIGH";}elsif ($anal_totreads >= 0 ){$anal_totreads_qual = "LOW";}
		my $anal_mapped_qual = "NA";
		if ( $anal_mapped >= $minmapped and $anal_mapped <= $maxmapped){$anal_mapped_qual = "HIGH";}elsif ($anal_mapped >= 0){$anal_mapped = "LOW";}
		my $anal_duplicates_qual = "NA";
		if ( $anal_duplicates >= $minduplicates and $anal_duplicates <= $maxduplicates){$anal_duplicates = "HIGH";}elsif ($anal_duplicates >= 0){$anal_totreads_qual = "LOW";}
	
		$head_str .= "$samplename\t$samplename\_QUAL\t";
		$totr_str .= "$anal_totreads\t$anal_totreads_qual\t";
	  $map_str .= "$anal_mapped\t$anal_mapped_qual\t";
	  $dup_str .= "$anal_duplicates\t$anal_duplicates_qual\t";
	  
	}
	chop($head_str);
	chop($totr_str);
	chop($map_str);
	chop($dup_str);
	
	print ALN_ST "-\tMIN\tMAX\tAVERAGE\t$head_str\n".
					$cfg_hash->{'db_sampstats_totreads'}."\t$minreads\t$maxreads\t$meanreads\t$totr_str\n".
					$cfg_hash->{'db_sampstats_mappedreads'}."\t$minmapped\t$maxmapped\t$meanmapped\t$map_str\n".
					$cfg_hash->{'db_sampstats_duplicates'}."\t$minduplicates\t$maxduplicates\t$meanduplicates\t$dup_str\n";	
	#Close file
	close(ALN_ST);	
	#Add the output to the file list
	append_str_2_file_if_path_notexist($tables_list_file,$cfg_hash->{'ALNSTATSCOMP_table'}."\t".$alnstats);	
			
}
	
=head2 print_alignment_statistics_table

 Title   : print_alignment_statistics_table
 Usage   : print_alignment_statistics_table(  config_file => the config hash
								);

 Function: The coverage analysis implies different steps. In particular this function
					is executed per-sample.
						
 
 Returns :
 
=cut
sub print_alignment_statistics_table{
	my $cfg_hash = shift;
	my $flagstat_aln = shift;
	my $flagstat_aln_sort_mrdup = shift;
	my $reads_stats_tab = shift;
	my $analysis_id = shift;
	my $sample_name = shift;
	my $log_file = shift;
	
	##############################
	#Here I get a table with information from flagstat

	my $sampleid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
											$cfg_hash->{'db_analysis_id'}.','.$cfg_hash->{'db_sample_name'},$analysis_id.",'".$sample_name."'");	
												
	#Obtain the kinship from the database given the sample name
	my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
											$cfg_hash->{'db_analysis_id'}.','.$cfg_hash->{'db_sample_name'},$analysis_id.",'".$sample_name."'");							
	my $kinship_types = join("",split(",",$cfg_hash->{'kinship_types'}));
	my $new_row = "";
	if ( $kinship =~ /[$kinship_types]/){
		$new_row .= "$sample_name\t$kinship";
	}else{
		$new_row .= "$sample_name\t-";
	}
	my $sample_present = 0;
	
	#Open the table file and concatenate if it is already there, otherwise create a new one
	if ( (-e $reads_stats_tab) and (!-z $reads_stats_tab) ){
		#Check if the sample is not present into the file
		extract_columns_from_file($reads_stats_tab,"0",$reads_stats_tab.".temp");
		my @written_s =  list_to_array($reads_stats_tab.".temp",'NO_NEW_LINE');
		if ( grep {/\b$sample_name\b/} @written_s){$sample_present = 1;}
 
		open (SAMPLE_STATS, ">>$reads_stats_tab") or die ("Cannot open $reads_stats_tab\n");
  	print_and_log( "Opened an old $reads_stats_tab...",$log_file);#DEBUGCODE

	}else{
		open (SAMPLE_STATS, ">$reads_stats_tab") or die ("Cannot open $reads_stats_tab\n");
  	print_and_log( "Opened a new $reads_stats_tab...",$log_file);#DEBUGCODE
		
		#Print the header of this table
		print SAMPLE_STATS "sample_name\tkinship\t".$cfg_hash->{'flagstat_totreads_field'}."\t".
			$cfg_hash->{'flagstat_readsmapped_field'}."\t".$cfg_hash->{'flagstat_readsremoved_field'}."\t".$cfg_hash->{'flagstat_propaired_field'}."\n";
	}
	#Needed variables
	
	my $tot_reads = -1;
	my $tot_mapped = -1;
	my $duplicate_reads = -1;
	my $properly_paired = -1;
	my $at_least_one = 0;
	my $success = 0;
	
	if ( $sample_present  == 0 ){
		
		#Get the number of total reads and properly paired from the file just after BWA
		if (-e $flagstat_aln and (!-z $flagstat_aln) ){
			#print_and_log( "Get the number of total reads and properly paired from the file $flagstat_aln just after BWA...",$log_file);#DEBUGCODE
			#From the output of flagstat take the total read number and the properly paired
			#Get only the first column
			my $flagstat_firstcol = $flagstat_aln.".temp";
			extract_colnum_from_file_linux($flagstat_aln,1,$flagstat_firstcol,' ');
			#Put into an array and pick the columns related with the information wanted
			my @stats_list = list_to_array($flagstat_firstcol,'NO_NEW_LINE');
			$tot_reads = $stats_list[$cfg_hash->{'flagstat_totreads_row'}];
			$tot_mapped = $stats_list[$cfg_hash->{'flagstat_readsmapped_row'}];
			$properly_paired = $stats_list[$cfg_hash->{'flagstat_propaired_row'}];
			$at_least_one++;
		}else{$tot_reads = "-";}

		#Get the number of total reads and properly paired from the file after the duplicates removal
	  #The higher level of BAM used will be picked for the statistics
	  if (-e $flagstat_aln_sort_mrdup and (!-z $flagstat_aln_sort_mrdup)  ){
			#print_and_log( "Get the number of total reads and properly paired from the file $flagstat_aln_sort_mrdup after the duplicates removal...",$log_file);#DEBUGCODE
			#From the output of flagstat take the total read number and the properly paired
			#Get only the first column
			my $flagstat_firstcol = $flagstat_aln_sort_mrdup.".temp";
			extract_colnum_from_file_linux($flagstat_aln_sort_mrdup,1,$flagstat_firstcol,' ');
			#Put into an array and pick the columns related with the information wanted
			my @stats_list = list_to_array($flagstat_firstcol,'NO_NEW_LINE');
			
			$tot_mapped = $stats_list[$cfg_hash->{'flagstat_readsmapped_row'}];
			$duplicate_reads = $stats_list[$cfg_hash->{'flagstat_readsremoved_row'}];	
			$properly_paired = $stats_list[$cfg_hash->{'flagstat_propaired_row'}];
			#print_and_log( "Get the removed duplicates from file $flagstat_firstcol row:".$cfg_hash->{'flagstat_readsremoved_row'}.":$duplicate_reads..",$log_file);#DEBUGCODE

			$at_least_one++;		
		}else{$duplicate_reads = "-";}			
		
		
		#At least one, means that at least one among flagstat after BWA ore duplicate removal is found
		if ( $at_least_one > 0 ){		
			
			$success = 1;
			#Store the information for the current analysis
			print_and_log( "Writing: $new_row\t$tot_reads\t$tot_mapped\t$duplicate_reads\t$properly_paired",$log_file);#DEBUGCODE
			print SAMPLE_STATS "$new_row\t$tot_reads\t$tot_mapped\t$duplicate_reads\t$properly_paired\n";
			
			#If the removed duplicates have been calculated put them into the db
			if ( $at_least_one > 1 ){	
				#Update statistics table if the element do not exists, otherwise update it
				my $id = -1;
				my $fields =  $cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_sampstats_totreads'}.",".$cfg_hash->{'db_sampstats_mappedreads'}
								.",".$cfg_hash->{'db_sampstats_duplicates'}.",".$cfg_hash->{'db_analysis_id'};
				my $values = "$sampleid,$tot_reads,$tot_mapped,$duplicate_reads,$analysis_id";
				my ($new_id,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_statistics'},$cfg_hash->{'db_sample_id'},
										$fields,$values,$cfg_hash->{'db_sample_id'},$sampleid);

				#If the sampleid exists insert that, otherwise the new one
				if ( $new_id >= 0){
					print_and_log( "Inserted the reads removed ($duplicate_reads) into ".$cfg_hash->{'db_sample_statistics'}." for sample $sampleid\n",$log_file);#DEBUGCODE
				}elsif ( $exist_id >= 0) {
					update_table_2($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
								$cfg_hash->{'db_sample_statistics'},$fields,$values,$cfg_hash->{'db_sample_id'},$sampleid);
				}else{
					log_and_exit("ERROR: during the import of the statistics status for $sample_name ($sampleid : analysisid=$analysis_id). Cannot understand what is happen with the db...\n",$log_file);
				}		
			}
		}else{
			print_and_log( "WARNING : Cannot find any flagstat output to use for sample $sample_name\n",$log_file);
		}
	}
	
	#Sample Reads Alignment statistics
	close (SAMPLE_STATS);				
}

1;
