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
    
    
package LIB::programs_management;

## programs_management.pm
#Permits the management of an input configuration file whose variables
#will be put inside an hash

#BEGIN { print (("  " x $main::x++) . "Beginning progmanagement compile\n") }
BEGIN
{
	require Exporter;
	use vars qw(@ISA @EXPORT);
	@ISA = qw(Exporter);
	@EXPORT_OK = qw( initialize_folders checkConfigVariables checkVariable
									 configFile2Hash check_input_fasta_file	correct_type 
									 get_chosen_session hash2ConfigFile
									 try_exec_command try_exec_job print_and_log log_and_exit
									 execute_threads build_input_name_from_executed separate_input_ids
									 kill_job remove_sample_sheet show_analyses_from_sample_sheet get_ped_file
									 sort_samples remove_analyses clean_folders get_sorted_sample_names
									 save_analyses R_die separate_bed_perchr show_analysesids_from_analysesnames
									 get_ped_file_for_enlarged_analysis show_samples_from_group
									 check_config_parameters check_genedb_links get_flowcell_and_lane
									 get_number_of_variants good_reads_quality execute_jobs
									 JOB_get_jobid_form_jobname alter_job JOB_get_info_from_jobname
									 store_and_compress_analyses outlist_to_hash_2
									 exec_store_and_compress_analyses parse_cand_genes_2_db
									 check_datasets_links load_table_into_db get_output_name_from_executed
									 fields_2_db run_FREEBAYES_vc get_samples_id_from_analysis
									 show_samplesids_from_samplesnames get_extended_regions
									 get_zygosity compress_data_to_store update_dependencies_hash
									 JOB_get_status_field JOB_is_active_queue JOB_get_dependencies);
}
             
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Cwd;#To change work directory
use Parallel::ForkManager;#For parallel executions
use File::Copy;#To manage files
use LWP::Simple;#For simple URL check (head)

#Using a library for database management
use LIB::db_management qw(get_id_if_exists_from_db db_print_selected_rows update_table 
													do_query_select_all delete_from_db get_id_if_exists_from_db_woconn
													insert_into_table_locked_woconn insert_only_if_not_exists_locked
													update_table update_table_2);

use LIB::files_management qw( extract_col_from_file check_presence check_URL my_extract_any_file
														extract_archive compress_folder delete_directory tar_folder_sys
														gzip_folder extract_name file_not_present
														list_to_array delete_file is_folder_empty
														create_folder save_hash load_hash
														copy_folder );
			
use LIB::html_management qw( html_page paths_hash_to_html);			
			
			
			

		
			
=head2 get_chosen_session

 Title   : get_chosen_session
 Usage   : get_chosen_session( -workingFolder => 'the folder where you search'
                               );

 Function: Opens the working directory and gives to the user the possibility to choose one

 Returns : nothing

=cut
sub get_chosen_session {
        my $workingFolder = shift;

        my $chosenFolder = "";

        if ( !(is_folder_empty($workingFolder)) ){
                my $workDir = getcwd;
                chdir $workingFolder;
                (system("ls -d */")) == 0 or die "Unable to list directory $workingFolder!\n";
                print "Choose your folder: ";
                $chosenFolder = <STDIN>;
                chomp $chosenFolder;
               	$chosenFolder = $1 if($chosenFolder=~/(.*)\/$/);#Cut last letter if it is /
                #print $myFolder;
                #my $validSession = 0;
                #if (valid_session($myFolder) == 1){
                                #$validSession = 1;
                #}
                while ( !(-d $chosenFolder)){
                        (system("ls -d */")) == 0 or die "Unable to list directory $chosenFolder!\n";
                       	print "Choose your folder: ";
                       	$chosenFolder = <STDIN>;
                        chomp $chosenFolder;
                        $chosenFolder = $1 if ($chosenFolder=~/(.*)\/$/);#Cut last letter if it is /
                        #if (valid_session($myFolder) == 1){
                                #$validSession = 1;
                        #}
                }
                $chosenFolder =  $workingFolder."/".$chosenFolder;
                chdir $workDir;
        }else{
                print "There are not existing sessions in $workingFolder. Please restart the program and create a new one!\n";
                exit 1;
        }

        return $chosenFolder;
}

=head2 load_table_into_db

 Title   : load_table_into_db
 Usage   : load_table_into_db( -database => 'name of the database,
                               );

 Function: Takes in input a table and for generates a table to be imported into the database
						getting a record for each column and for each row
					
						$id1 -> the  identifier of the columns to get the dbid 
						$id2 -> the identifier of the row used to get the dbid
						$id3 -> the identifier for the values that also is a table
						
 Returns : returns a table
=cut
sub load_table_into_db{
	my $cfg_hash = shift;
	my $input_tab = shift;
	my $analysisid = shift;
	my $id1 = shift;
	my $id2 = shift;
	my $id3 = shift;
	my $log_file = shift;
	#my $out_tab = shift;
	
	my @header = ();
	
	#If the input table exists
	if( -e $input_tab){
		#Connect to database
		my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
		# PERL DBI CONNECT
		my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
		#Open both input and output
		open(IN,"<".$input_tab) or die "ERROR: Cannot open file $input_tab\n";
		#open(OUT,">".$out_tab) or die "ERROR: Cannot open file $out_tab\n";
		
		my $num_row = 1;
		#For each line of the input table
		while  (my $row = <IN>){
			chop($row);
			my @fields = split("\t",$row);
			#From the header get the db identifiers for all columns but not the first
			#that is the element that is the rowname
			if ($num_row == 1){
				my $num_col = 0;
				foreach my $field (@fields){
					if( $num_col > 0 ){
						my $dbid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{"db_$id1\_table"},$cfg_hash->{"db_$id1\_id"},
										$cfg_hash->{"db_$id1\_name"}.",".$cfg_hash->{"db_analysis_id"},"'".$field."',".$analysisid);
						if ($dbid > 0){
							push(@header,$dbid);
					  }else{
								die "ERROR: no results from gettin id from the database.".
									"table: ".$cfg_hash->{"db_$id1\_table"}." fields: ".$cfg_hash->{"db_$id1\_name"}.",".$cfg_hash->{"db_analysis_id"}." values : $field, $analysisid";
						}
					}
					$num_col++;
				} 
			}else{
				my $num_col = 0;
				my $dbid = -1;
				foreach my $field (@fields){
					#Get the DB id from the first column
					if( $num_col == 0 ){
						$dbid = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{"db_$id2\_table"},$cfg_hash->{"db_$id2\_id"},
						$cfg_hash->{"db_$id2\_name"},"'".$field."'");
						#If the dbid is not present insert the value and get a new id
						if ($dbid < 0){
							#INSERT the value not present in table id2
							my $fields = $cfg_hash->{"db_$id2\_name"};
							my $values = "'".$field."'";											
							#Inserts the  id into the table only if that gene does not exists
							$dbid = insert_into_table_locked_woconn($dbh,$cfg_hash->{"db_$id2\_table"},$fields,$values);
						}									
					}
					#For each remaining column put into the table
					else{
						#print OUT $header[$num_col-1]."\t$dbid\t$field\n";
						my $fields = $cfg_hash->{"db_$id1\_id"}.",".$cfg_hash->{"db_$id2\_id"}.",".$cfg_hash->{"db_$id3"};
						my $values = $header[$num_col-1].",$dbid,$field";
						#print_and_log("Inserting $fields -> $values into ".$cfg_hash->{"db_$id3\_table"},$log_file);
						insert_into_table_locked_woconn($dbh,$cfg_hash->{"db_$id3\_table"},$fields,$values);
					}
					$num_col++;
				}				
			}

			$num_row++;
		}
		close(IN);
		#close(OUT);
		#Disconnect db
		$dbh->disconnect(); 
	}
	
}	

	
=head2 get_samples_id_from_analysis

 Title   : get_samples_id_from_analysis
 Usage   : get_samples_id_from_analysis( -database => 'name of the database,
                               );

 Function: Makes a query to the database to select all samples associated with an analysis id
						Needs in input the analysis id for which you want the samples 

 Returns : returns an array with sampleids

=cut
sub get_samples_id_from_analysis{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	
	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	
	#print_and_log( "$analysis_indication Executing: $query\n",$log_file);#DEBUGCODE
	my $samples_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	my @idsarray = ();
	foreach my $sampleid (keys %{$samples_ids}){
		push(@idsarray,$sampleid );
	}
}
	
=head2 get_sorted_sample_names

 Title   : get_sorted_sample_names
 Usage   : get_sorted_sample_names( -database => 'name of the database,
                               );

 Function: Sorts the samples belonging to the given group_id

 Returns : returns an array with sample_names sorted as for the kinship given
						in the configuration file. Uses the function sort_samples

=cut
sub get_sorted_sample_names{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	
	my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	print_and_log( "Executing: $query\n",$log_file);
	my $analysis_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});


	#Getting the samples list associated to the groupid
	my $samples_h;#Samples hash to be reordered
	my @sort_samples = ();#The array with samples id reordered
	
	#Get the kinship, to make the resorting
	foreach my $sample_name (keys %{$analysis_sam}){
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
		$samples_h->{$sample_name}->{'id'} = $sample_name;
		$samples_h->{$sample_name}->{'k'} = $kinship;
	}
	#Sort the samples as in samples_order parameter
	sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);
	
	return @sort_samples;
}
	
	
=head2 rank

 Title   : rank
 Usage   : rank(  );

 Function:  A function to make a comparison between elements using a
						order given in input
						
 Returns : returns the position of the element

=cut
sub rank {
	
	my $l = shift;
	my $s_order = shift;
	
	if ($l eq $s_order->[0]){return 0;}
	if ($l eq $s_order->[1]){return 1;}
	if ($l eq $s_order->[2]){return 2;}
	if ($l eq $s_order->[3]){return 3;}
	if ($l eq $s_order->[4]){return 4;}

}

=head2 sort_samples

 Title   : sort_samples
 Usage   : sort_samples(  );

 Function:  This subroutine takes in input an hash with sample information
						$s->{sample_id}->{name} and $s->{sample_id}->{kinship}
						then sorts the samples based on the kinship order in 
						program_config, parameter samples_order
					
						Exploits the function rank to sort based on the given order
						
 Returns : writes in an array given in input the samples ids sorted (or names if names are given)

=cut
sub sort_samples{
	my $samples = shift;
	my $order = shift;
	my $sort_samples = shift;
		
	my @s_order = split(",",$order);	

 #foreach my $key ( #
  #sort  { rank($samples->{$a}->{'k'},\@s_order) <=> rank( $samples->{$b}->{'k'},\@s_order) }
        #keys %{$samples})
	 #{	#print $samples->{$key}->{'n'}."\n";
			#push(@$sort_samples,$samples->{$key}->{'id'});
	 #}

 foreach my $key ( #
  sort  { rank($samples->{$a}->{'k'},\@s_order) <=> rank( $samples->{$b}->{'k'},\@s_order)  or 
					$samples->{$a}->{'id'} cmp $samples->{$b}->{'id'}
					}
        keys %{$samples})
	 {	#print $samples->{$key}->{'n'}."\n";
			push(@$sort_samples,$samples->{$key}->{'id'});
	 }
}

=head2 clean_folders

 Title   : clean_folders
 Usage   : clean_folders(  );

 Function: Clean the folders of the groups given in input
						- Remove all SAM files generated by BWA
						- Remove all bam files not sorted and not indexed
						- Remove all fastq useless files
						- Remove all the content of the refine, varcall and genotype folders
 Returns : 

=cut
sub clean_folders{
	my $cfg_hash = shift;
	my $group_ids = shift;
	my $work_folder = shift;
	my $log_file = shift;
	
	#my @group_ids = split(",",$groups_l);
	#foreach my $group_id (@$group_ids){
	
	for ( my $i = 0; $i < scalar(@$group_ids); $i++){
		#Obtain the group name from the database given the group id
		my $group_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$group_ids->[$i]);

		my $outFolder = $work_folder."/$group_name/".$cfg_hash->{'align_out_f'};					
		print_and_log("\n\n#########Cleaning folder $outFolder for group ".$group_ids->[$i].":  ..\n",$log_file);
							
		#if (!is_folder_empty($outFolder)){print "Removing from $outFolder..\n";}#DEBUGCODE
		
		#Remove all SAM files generated by BWA
		print_and_log("Removing all SAM files from $outFolder:".$outFolder."/*.".$cfg_hash->{'sam_ext'}."...\n",$log_file);
		delete_file($outFolder."/*.".$cfg_hash->{'sam_ext'});
		
		#Remove all bam files not sorted and not indexed
		print_and_log("Removing all bam files not sorted and not indexed from $outFolder:".
			      $outFolder."/*_".$cfg_hash->{'align_step'}."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'}.
			      "...\n",$log_file);
		delete_file($outFolder."/*_".$cfg_hash->{'align_step'}."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'});
		#Remove all bam files obtained on read files by lane
		print_and_log("Removing all bam files obtained from lanes from $outFolder:".
			      $outFolder."/*_L00*.".$cfg_hash->{'bam_ext'}."...\n",$log_file);
		delete_file($outFolder."/*_L00*.".$cfg_hash->{'bam_ext'});
		
		#Remove all fastq useless files
		print_and_log("Removing all fastq useless files from $outFolder:".
			      $outFolder."/*_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'}.
			      "...\n",$log_file);
		delete_file($outFolder."/*_".$cfg_hash->{'mergepe_step'}.".".$cfg_hash->{'fastq_ext'});
		
		#Refine folder
		$outFolder = $work_folder."/$group_name/".$cfg_hash->{'refine_out_f'};
		#Remove all the content of the refine folder
		print_and_log("Removing all the content of the folder: $outFolder\n",$log_file);
		delete_file($outFolder."/*");
		
		#varcall folder
		$outFolder = $work_folder."/$group_name/".$cfg_hash->{'varcall_out_f'};
		#Remove all the content of the varcall_out folder
		print_and_log("Removing all the content of the folder: $outFolder\n",$log_file);
		delete_file($outFolder."/*");

		#genotype folder
		$outFolder = $work_folder."/$group_name/".$cfg_hash->{'genot_out_f'};
		#Remove all the content of the genotype folder
		print_and_log("Removing all the content of the folder: $outFolder\n",$log_file);
		delete_file($outFolder."/*");
	
		print_and_log("#########$work_folder/$group_name/ cleaned!\n\n",$log_file);
	}
}

=head2 get_extended_regions

 Title : get_extended_regions
 Usage : get_extended_regions( );

 Function: given a bed file in input and the dimension. 

 Returns : writes a new bed file where each region is extended of the given
						number of nucleotides before and after

=cut 
# Given a bed file with the first three lines giving the chromosome, the start
# and the end. If there are two identical consecutive intervals keeps only the first one
sub get_extended_regions {
	my $bed = shift;
	my $outname = shift;
	my $dim = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname") or die "Cannot open $outname\n";
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines and printing extended bed into $outname\n";
	
	my $prevchr = "";
	my $prevend = -1;
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		if ( $lines[$i] =~ /^chr/){
			chomp($lines[$i]);
			chomp($lines[$i+1]);
			my @fields1 = split("\t",$lines[$i]);
			my @fields2 = split("\t",$lines[$i+1]);
			
			my $thischr = $fields1[0];
			my $thisstart = $fields1[1];
			my $thisend = $fields1[2];
			my $nextchr = $fields2[0];
			my $nextstart = $fields2[1];
			my $desc = "";
			if (defined $fields1[3]){
				$desc = $fields1[3];
			}
			
			#newstart is $dim nucleotides before and newend $dim after 
			my $newstart = $thisstart - $dim;
			my $newend = $thisend + $dim;
			
			#if now newstart is negative, set to 0
			if ( $newstart < 0 ){
					$newstart = 0;
			}
			#and if the newstart overlaps the previous end:
			if ( $newstart <= $prevend  and $thischr eq $prevchr){
					$newstart = $prevend + 1;
			}
			#or if the newend overlaps the next start:
			if ( $newend >= $nextstart  and $thischr eq $nextchr){
					$newend = $nextstart - 1;
			}
			#Set the new previous details
			$prevchr = $thischr;
			$prevend = $newend;
			print NEWFILE "$thischr\t$newstart\t$newend";
			if ( $desc ne '' ){
				print NEWFILE "\t$desc";
			}
			print NEWFILE "\n";	
		}
		
	}
}

=head2 get_number_of_variants

 Title   : get_number_of_variants
 Usage   : get_number_of_variants(  );

 Function: Given a VCF file simply executes grep -v '^#' file.vcf | wc -l
					to count the number of variants
						
 Returns : the number of variants

=cut
sub get_number_of_variants{
	my $vcf_file = shift;
	
	my $var_num = `grep -v '^#' $vcf_file | wc -l`;
	
	return $var_num;
}


=head2 get_ped_file_for_enlarged_analysis

 Title   : get_ped_file_for_enlarged_analysis
 Usage   : get_ped_file_for_enlarged_analysis(  );

 Function: Obtains a unique pedigree file for a large set of samples.
					Takes a sample and removes the last _X char to obtain the family name. 
					Then searches for the _P and _M (mother) and _F( father) and writes 
					the ped file
						
 Returns :

=cut
sub get_ped_file_for_enlarged_analysis{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $out_ped = shift;
	my $log_file = shift;
	
	my $sep = "\t";
	#Get the samples ids involved for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".
							$cfg_hash->{'db_analysis_id'}."= $analysis_id;";	
	#print_and_log( "Executing: $query\n",$log_file);
	my $analysis_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
	
	my $num_samples = scalar(keys %{$analysis_sam});

	#Prepare the analysisid as field of the query	
	my $fields = $cfg_hash->{'db_analysis_id'};
	my $values = $analysis_id;
	
	#Indicates if there are families or just have to write info for single samples
	my $nofamilies = 0;
	
	my @uniq_fam_names = ();
	my @single_samples_names = ();
	#Get a list of unique family names
	foreach my $sample_name (keys %{$analysis_sam}){
		my @sname_fields = split("_",$sample_name);
		pop(@sname_fields);
		#Get the family name removing the last _X
		my $family_name = join("_",@sname_fields);
		
		#Check that is actually a family by verifying that _P or _P1 are present in the db
		my $proband1 = $family_name."_P";
		my $proband2 = $family_name."_P1";		
		my $proband_id1 = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband1."'");
		my $proband_id2 = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband2."'");
		if ( $proband_id1 > 0 or  $proband_id2 > 0 ){																	
			if ( not grep {/\b$family_name\b/} @uniq_fam_names){
				push(@uniq_fam_names,$family_name);
			}
		}else{
			push(@single_samples_names,$sample_name);
		}
	}	
	#Open the ped file for writing
	open(PED_F,">".$out_ped) or die "Cannot open file $out_ped\n";
			

	#For each family add its ped file into the global one
	foreach my $family (@uniq_fam_names){
			
			#Get the mother
			my $mother_name = $family."_M";
			#Look if the sample exists in the database given the group id
			my $mother_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$mother_name."'");	
			#Get the father
			my $father_name = $family."_F";
			#Look if the sample exists in the database given the group id
			my $father_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$father_name."'");			
			#Print mother and father informaton
			print PED_F $family.$sep.$father_name.$sep."0".$sep."0".$sep."1".$sep."1\n";
			print PED_F $family.$sep.$mother_name.$sep."0".$sep."0".$sep."2".$sep."1\n";
			#Get the probands
			#If there is a single proband we have family_name_P
			my $proband_name = $family."_P";		
			#Look if the sample exists in the database given the group id
			my $proband_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");
			#But if we cannot find that name, we can go on to see if there are more probands and write all of them
			if ( $proband_id < 0){
				my $num_probands = 1;	
				$proband_name = $family."_P$num_probands";		
				#Look if the sample exists in the database given the group id
				my $proband_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");				

				#Look for other probands _P2,_P3 ...etc
				while ( $proband_id > 0){
					my $gender = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_gender'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");
					my $gender_id = -1;
					if (defined $gender){
						if ($gender eq 'M'){ $gender_id = 1;}else{$gender_id = 2;}
					}
					print PED_F $family.$sep.$proband_name.$sep.$father_name.$sep.$mother_name.$sep.$gender_id.$sep."2\n";
					$num_probands++;	
					$proband_name = $family."_P$num_probands";		
					#Look if the sample exists in the database given the group id
					$proband_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");								
				}				
			}else{
				my $gender = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_gender'},
																$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");
				my $gender_id = -1;
				if (defined $gender){
					if ($gender eq 'M'){ $gender_id = 1;}else{$gender_id = 2;}
				}
				print PED_F $family.$sep.$proband_name.$sep.$father_name.$sep.$mother_name.$sep.$gender_id.$sep."2\n";
			}					
	}

	#For each single samples  add it into the pedfile
	foreach my $proband_name (@single_samples_names){
		my $gender = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_gender'},
														$fields.",".$cfg_hash->{'db_sample_name'},$values.",'".$proband_name."'");
		my $gender_id = -1;
		if (defined $gender){
			if ($gender eq 'M'){ $gender_id = 1;}else{$gender_id = 2;}
		}
		print PED_F $proband_name.$sep.$proband_name.$sep."0".$sep."0".$sep.$gender_id.$sep."2\n";
	}
	close(PED_F);	
}


=head2 get_zygosity

 Title   : get_zygosity
 Usage   : get_zygosity( -database => 'name of the database,
                               );

 Function: determines from the VCF genotype if the variant is
			HOM, HOMREF, HET
 
 Returns : a string with a genotype (HOM, HOMREF, HET)

=cut
sub get_zygosity{
	my $gt = shift;

	my $out_str = "";
	#print "Getting zygosity for $gt\n";#DEBUGCODE
	if (defined $gt and $gt =~ /[\|,\/]/){
		#Now determine the Zygosity. Separate the GT using the slash or |
		my @gens = split("[\|,\/]",$gt);
		
		#If the alleles are the same the sample is homozygous
		if ( $gens[0] eq $gens[1]){
			#0/0
			if ( $gens[0] eq 0){
				$out_str =  "HOMREF";
			}#Situation with ./.
			elsif ( $gens[0] eq '.' ) {
				$out_str =  "UNDEF";
			}#1/1
			else{
				$out_str =  "HOM";
			}
		#If alleles are different the sample is heterozygous	
		}else{#0/1 or 1/0 or 0|1 or 1|0
			$out_str =  "HET";
		}		
	}

	return $out_str;
}


=head2 get_ped_file

 Title   : get_ped_file
 Usage   : get_ped_file(  );

 Function:  creates a PED file for any analysis.
						returns 2 if the analysis identifier given in input is a trio.
						To check this accesses the database and looks if there are 
						three samples and there are M,F and P (Mother, Father and Proband)
						
 Returns : 2 if is a trio, 1 if is not a trio, 0 if the pedfile cannot be created

=cut
sub get_ped_file{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $out_ped = shift;
	my $log_file = shift;
		
	my $retval = 0;
	my $sep = "\t";
	my $analysis_indication = "[an:$analysis_id ](".scalar(localtime)."): ";
	
	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".
							$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "Executing: $query\n",$log_file);
	my $analysis_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	my $num_samples = scalar(keys %{$analysis_sam});
	
	my $num_probs = 0;
	my $genders;
	
	#The pedfile is created if there is at least one sample
	if( $num_samples > 0){
		my $moth_name = "-";
		my $fath_name  ="-";
		my $prob_name  ="-";
		#Used to check if at least one kinship has been defined
		my $defined_kinships = 0;
		#Print the header for information about the genotype for the sample
		foreach my $sample_id (keys %{$analysis_sam}){
			#Obtain the kinship
			my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'}
														,$cfg_hash->{'db_sample_id'},$sample_id);
			
			if (defined $kinship){
				$defined_kinships++;
				#Found the father!
				if( $kinship eq 'F' ){
					#Getting the sample name
					$fath_name  = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
												$cfg_hash->{'db_sample_id'},$sample_id);		
				}
				#Found the mother!
				elsif( $kinship eq 'M' ){
					#Getting the sample name
					$moth_name  = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
												$cfg_hash->{'db_sample_id'},$sample_id);						
				}
				#Found a proband!
				elsif($kinship =~ 'P' ){
					$num_probs++;
					#Getting the sample name
					my $prob_name  = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
												$cfg_hash->{'db_sample_id'},$sample_id);					
					my $gender = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_gender'},
																$cfg_hash->{'db_sample_id'},$sample_id);
					#Use gender if present
					if (defined $gender){
						if ($gender eq 'M'){ $genders->{$prob_name} = 1;}else{$genders->{$prob_name} = 2;}
					}else{
						$genders->{$prob_name} = "-";
					}
				}
				##TEMPORANEO
				#If there is another family member, the analysis is not of a trio
				else{
						$retval = 1;
				}				
			}else{
				print_and_log("$analysis_indication Kinship is not defined for sample $sample_id..\n",$log_file);
			}
		}
		#Only if at least one kinship is defined, write a ped
		if ($defined_kinships > 0){
			#If all F,M and P have been found set return accordingly
			if( $fath_name ne '-' and $moth_name ne '-' and $num_probs == 1 and $retval == 0){
					$retval = 2;
			}else{
					$retval = 1;
			}
			
			#Obtain the analysis name from the database given the analysis id
			my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
															$cfg_hash->{'db_analysis_id'},$analysis_id);

			#Print the lines of the file
			open(PED_F,">".$out_ped) or die "Cannot open file $out_ped\n";
			#Print parents if exist
			if( $fath_name ne '-'){print PED_F $analysis_name.$sep.$fath_name.$sep."0".$sep."0".$sep."1".$sep."1\n"; }
			if( $moth_name ne '-'){print PED_F $analysis_name.$sep.$moth_name.$sep."0".$sep."0".$sep."2".$sep."1\n"; }
			#Print probands			
			foreach my $prob_name ( keys %$genders){
				print PED_F $analysis_name.$sep.$prob_name.$sep.$fath_name.$sep.$moth_name.$sep.$genders->{$prob_name}.$sep."2\n";
			}
			close(PED_F);
			
			print_and_log("$analysis_indication Pedfile written at: $out_ped\n",$log_file);
			#Add pedfile name in the database
			#UPDATE fields = values
			my $fields = $cfg_hash->{'db_pedfile'};
			my $values = "'".$analysis_name.".".$cfg_hash->{'ped_ext'}."'";
			#WHERE 
			my $fields2 = $cfg_hash->{'db_analysis_id'};
			my $values2 = $analysis_id;

			#Update
			update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);
			#print_and_log("$analysis_indication Updated table ".$cfg_hash->{'db_analyses_table'}." with the pedfile name\n",$log_file);#DEBUGCODE
		}			
	}

	return $retval;
}

=head2 remove_sample_sheet

 Title   : remove_sample_sheet
 Usage   : remove_sample_sheet(  );

 Function:  removes a sample sheet from the database given a list of analysis ids
						or given a list of sample sheet names
						
 Returns : nothing

=cut
sub remove_sample_sheet{
	my $cfg_hash = shift;
	my $log_file = shift;
	my $ss_to_rem = shift;
	
	print_and_log("Removing a sample sheet from the database..\n",$log_file);
	my $sep = $cfg_hash->{'parameters_sep'};


	#If a string has given in input..
	if ( !correct_type($ss_to_rem,"positiveint") ){
		#Obtain the group name from the database given the group id
		my $ss_to_del = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_id'}
										,$cfg_hash->{'db_sample_sheet_name'},"'".$ss_to_rem."'");
		if ( correct_type($ss_to_del,"positiveint") ){
			print_and_log("Removing sample sheet $ss_to_del from database:\n",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
			print_and_log($cfg_hash->{'db_readf_table'}."...",$log_file);					
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
			print_and_log($cfg_hash->{'db_sample_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
			print_and_log($cfg_hash->{'db_analyses_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
			print_and_log($cfg_hash->{'db_sample_sheets_table'}."...\n",$log_file);
		}else{print_and_log("Sample Sheet $ss_to_del is not present in db ".$cfg_hash->{'db_name'}."\n",$log_file);}
	}else{
		my @ssheets_to_del = split($sep,separate_input_ids($ss_to_rem,$sep));
		foreach my $ss_to_del (@ssheets_to_del){
			#If the parameter is a number then use db_sample_sheet_id, otherwise db_sample_sheet_name

			if ( correct_type($ss_to_del,"positiveint") ){					
				print_and_log("Removing Sample Sheet $ss_to_del from database:\n",$log_file);
				delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
				print_and_log($cfg_hash->{'db_readf_table'}."...",$log_file);					
				delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
				print_and_log($cfg_hash->{'db_sample_table'}."...",$log_file);
				delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
				print_and_log($cfg_hash->{'db_analyses_table'}."...",$log_file);
				delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_id'},$ss_to_del);
				print_and_log($cfg_hash->{'db_sample_sheets_table'}."...\n",$log_file);
			}else{print_and_log("Sample Sheet $ss_to_del is not present in db ".$cfg_hash->{'db_name'}."\n",$log_file);}
		}
	}
}

=head2 save_analyses

 Title   : save_analyses
 Usage   : save_analyses(  );

 Function:  Given a set of analysis ids moves all the output files needed into a folder named
			as the analysisname. 
			The files to move will be picked from the outfile_list file in the DATA folder
						
 Returns : nothing

=cut
sub save_analyses{
	my $cfg_hash = shift;
	my $working_folder = shift;
	my $ans_to_store_str = shift;
	my $force_run = shift;
	my $log_file = shift;	
	
	my $storage_dir = $cfg_hash->{'storage_f'};
	print_and_log("Storing files for the analysis into $storage_dir..\n",$log_file);
	#Get an array of group ids to store
	my $sep = $cfg_hash->{'parameters_sep'};
	my @ans_to_store = split($sep,separate_input_ids($ans_to_store_str,$sep));

	#Execute for each analysisid to store	 	
	foreach my $analysisid (@ans_to_store){
		
		
		#Obtain the stored variable
		my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_stored'},
											$cfg_hash->{'db_analysis_id'},$analysisid);	
		
		#If the analysis has not been already stored or the user wants to force the run
		if ( $stored == 0 or ($force_run) ){
			#If the parameter is a number then proceed
			if ( correct_type($analysisid,"positiveint") ){					
				
				#Obtain the analysis name from the database given the analysis id
				my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																				$cfg_hash->{'db_analysis_id'},$analysisid);											
				#Create a directory for the analysis into the storage dir
				my $analysis_fold_st = $storage_dir."/".$analysis_name;

				#Check if analysis_fold_st exists, otherwise it creates it
				unless(-d $analysis_fold_st){
					print_and_log("Creating folder $analysis_fold_st...\n",$log_file);
					mkdir $analysis_fold_st or die "ERROR: can't create folder $analysis_fold_st. Check permissions. \n";
				}		
											
				#Get the paths from the outlist file
				my $outlist_f = $working_folder."/".$analysis_name."/".$cfg_hash->{'data_fold'}."/".$analysis_name."_".$cfg_hash->{'outlist_file'};
				print_and_log("Getting the paths to move from $outlist_f..\n",$log_file);
				my $only_paths_f = $outlist_f.".onlypaths";
				extract_col_from_file($outlist_f,$cfg_hash->{'outlist_path_pos'},$only_paths_f);
				my @paths = list_to_array($only_paths_f,'NO_NEW_LINE');
				
				my $only_descs_f = $outlist_f.".onlydescs";
				extract_col_from_file($outlist_f,$cfg_hash->{'outlist_desc_pos'},$only_descs_f);
				my @descs = list_to_array($only_descs_f,'NO_NEW_LINE');	 

								
				#For each file path into the outlist file... 
				my $i = 0;
				foreach my $path ( @paths ){
					my $desc = $descs[$i];
					print_and_log("Checking if $desc matches one among IGV, SAVE, WEB.\n",$log_file); 
					#The description must have one among the following keywords
					if ( $desc =~ /WEB/ or $desc =~ /IGV/ or $desc =~ /SAVE/ ){
						print_and_log("Copying $path to  $analysis_fold_st..\n",$log_file);
						if ( -d $path ){
							copy_folder($path,$analysis_fold_st);			 
						}
						#Move the file into the storage directory
						copy($path,$analysis_fold_st) or print "ERROR: unable to copy $path in $analysis_fold_st\n";
					}else{
						print_and_log(" $path will not be copied to $analysis_fold_st..\n",$log_file);
					}
					$i++;
				}
				##################################################
				## UPDATE STORAGE STATUS
				##################################################
				
				#Update the analysis table with the fact that the analysis was stored (stored=1). 
				my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'},
									$cfg_hash->{'db_analysis_id'},$analysisid);

				my $fields = $cfg_hash->{'db_analysis_stored'};
				my $values = 1;
				my $fields2 = $cfg_hash->{'db_analysis_id'};
				my $values2 = $analysisid;
				update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
				##################################################		
						
			}else{print_and_log("Analsysis $analysisid is not a db id\n",$log_file);}			
		}else{print_and_log("The analysis $analysisid has been already stored. Use --force_run if you know what you are doing..\n",$log_file);}																	

	}
}	

#=head2 save_analyses

 #Title   : save_analyses
 #Usage   : save_analyses(  );

 #Function:  Given a set of analysis ids moves all the output files needed into a folder named
					#as the analysisname. 
					#The files to move will be picked from the outfile_list file in the DATA folder
					#The files moved will be substituted with symbolik links using the 'ln -s' command
						
 #Returns : nothing

#=cut
#sub save_analysesOLD{
	#my $cfg_hash = shift;
	#my $working_folder = shift;
	#my $ans_to_store_str = shift;
	#my $force_run = shift;
	#my $log_file = shift;	
	
	#my $storage_dir = $cfg_hash->{'storage_f'};
	#print_and_log("Storing files for the analysis into $storage_dir..\n",$log_file);
	##Get an array of group ids to store
	#my $sep = $cfg_hash->{'parameters_sep'};
	#my @ans_to_store = split($sep,separate_input_ids($ans_to_store_str,$sep));

	##Execute for each analysisid to store	 	
	#foreach my $analysisid (@ans_to_store){
		
		
		##Obtain the stored variable
		#my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_stored'},
											#$cfg_hash->{'db_analysis_id'},$analysisid);	
		
		##If the analysis has not been already stored or the user wants to force the run
		#if ( $stored == 0 or ($force_run) ){
			##If the parameter is a number then proceed
			#if ( correct_type($analysisid,"positiveint") ){					
				
				##Obtain the analysis name from the database given the group id
				#my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																				#$cfg_hash->{'db_analysis_id'},$analysisid);											
				##Create a directory for the analysis into the storage dir
				#my $analysis_fold_st = $storage_dir."/".$analysis_name;
				##Check if directory exists, otherwise creates it
				#if (! (-d $analysis_fold_st) ){

									
					##Get the paths from the outlist file
					#my $outlist_f = $working_folder."/".$analysis_name."/".$cfg_hash->{'data_fold'}."/".$analysis_name."_".$cfg_hash->{'outlist_file'};
					#print_and_log("Getting the paths to move from $outlist_f..\n",$log_file);
					#my $only_paths_f = $outlist_f.".onlypaths";
					#extract_col_from_file($outlist_f,$cfg_hash->{'outlist_path_pos'},$only_paths_f);
					#my @paths = list_to_array($only_paths_f,'NO_NEW_LINE');

					##Create a new folder
					#print_and_log("Creating folder $analysis_fold_st...\n",$log_file);
					#mkdir $analysis_fold_st or die "ERROR: can't create folder $analysis_fold_st. Check permissions. \n";
									
					##For each file path into the outlist file... 
					#foreach my $path ( @paths ){
						#print_and_log("Moving $path to  $analysis_fold_st..\n",$log_file);
						##Check if it is not a symbolink link
						#if ( ! (-l $path) ){
							##Move the file into the storage directory
							#move($path,$analysis_fold_st) or print "ERROR: unable to move $path in $analysis_fold_st\n";
						#}else{
						#print_and_log("ERROR: The file $path is already a symbolik link! Transfer was not performed...\n",$log_file);
						#}
						##delete_file($tempFile);
						##Create a symbolink link using the old name to the current path	
						#my $file_name = extract_name($path,0);
						#print_and_log("Creating symbolic link $analysis_fold_st/$file_name..\n",$log_file);
						#symlink $analysis_fold_st."/".$file_name, $path;
					#}
					###################################################
					### UPDATE STORAGE STATUS
					###################################################
					
					##Update the analysis table with the fact that the analysis was stored (stored=1). 
					#my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'},
										#$cfg_hash->{'db_analysis_id'},$analysisid);

					#my $fields = $cfg_hash->{'db_analysis_stored'};
					#my $values = 1;
					#my $fields2 = $cfg_hash->{'db_analysis_id'};
					#my $values2 = $analysisid;
					#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
					###################################################		
						
				#}else{print_and_log("WARNING: The folder $analysis_fold_st already exists. Check if the transfer has already done...\n",$log_file);}
			#}else{print_and_log("Analsysis $analysisid is not a db id\n",$log_file);}			
		#}else{print_and_log("The analysis $analysisid has been already stored. Use --force_run if you know what you are doing..\n",$log_file);}																	

	#}
#}	
	

=head2 store_and_compress_analyses

 Title   : store_and_compress_analyses
 Usage   : store_and_compress_analyses(  );

 Function:  Given a set of group ids moves all the output files needed into a folder named
					as the analysisname. 
					The files to move will be picked from the outfile_list file in the DATA folder
					There will be a single final folder that will be compressed and the HTML file 
					will now contain only the raw links
						
 Returns : nothing

=cut
sub exec_store_and_compress_analyses{
	my $cfg_hash = shift;
	my $ans_to_store_str = shift;
	my $force_run = shift;
	my $log_file = shift;
	
	my $workDir = getcwd;
	chdir($cfg_hash->{'work_fold'});
	my $cmd = $cfg_hash->{'program_folder'}."/".$cfg_hash->{'program_name'}." -c user_config.txt -vgstore $ans_to_store_str ";
	if ($force_run){$cmd .= " --force_run";}

	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};	
	my $task = "store";
	my $job_name = $task;
  print_and_log("Executing storage of analyses $ans_to_store_str \n",$log_file);# I'll write a log in $log.\n"; 	
	print_and_log("\nLaunching a job $job_name : ",$log_file);
	my $job_log = $cfg_hash->{'work_fold'}."/".$cfg_hash->{'logFolder'}."/$task.log";
	#Call the program using qsub and giving it parameters for the job and
	my $env_vars = "COMMAND='$cmd' ";
	my $job_id = try_exec_job( $cfg_hash,$env_vars,$task,$exec_job_script,$job_name,"no_depend",$job_log,$job_log,$log_file);	
	#return back
	chdir($workDir);
}


=head2 store_and_compress_analyses

 Title   : store_and_compress_analyses
 Usage   : store_and_compress_analyses(  );

 Function:  Given a set of group ids moves all the output files needed into a folder named
					as the analysisname. 
					The files to move will be picked from the outfile_list file in the DATA folder
					There will be a single final folder that will be compressed and the HTML file 
					will now contain only the raw links
						
 Returns : nothing

=cut
sub store_and_compress_analyses{
	my $cfg_hash = shift;
	my $working_folder = shift;
	my $ans_to_store_str = shift;
	my $force_run = shift;
	my $log_file = shift;	
	
	my $storage_dir = $cfg_hash->{'storage_f'};
	print_and_log("Storing files for the analysis into $storage_dir..\n",$log_file);
	#Get an array of analysis ids to store
	my $sep = $cfg_hash->{'parameters_sep'};
	my @ans_to_store = split($sep,separate_input_ids($ans_to_store_str,$sep));

	#Execute for each analysisid to store	 	
	foreach my $analysisid (@ans_to_store){
		
		
		#Obtain the group name from the database given the analysis id
		my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_stored'},
																			$cfg_hash->{'db_analysis_id'},$analysisid);	
		
		#If the analysis has not been already stored or the user wants to force the run
		if ( $stored == 0 or ($force_run) ){
			#If the parameter is a number then proceed
			if ( correct_type($analysisid,"positiveint") ){					
				
				#Obtain the analysis name from the database given the analysis id
				my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																				$cfg_hash->{'db_analysis_id'},$analysisid);	
			 print_and_log("\nStoring the analysis $analysis_name (analysis_id: $analysisid)..\n",$log_file);										
				#Create a directory for the analysis into the storage dir
				my $analysis_fold_st = $storage_dir."/".$analysis_name;
				#Check if directory exists and gives error if it does unless you do force the run
				if (! (-d $analysis_fold_st) or ($force_run)){

					my $paths_hash;				
					#Get the paths from the outlist file
					my $outlist_f = $working_folder."/".$analysis_name."/".$cfg_hash->{'data_fold'}."/".$analysis_name."_".$cfg_hash->{'outlist_file'};
					print_and_log("Getting the paths to move from $outlist_f..\n",$log_file);
					my $only_paths_f = $outlist_f.".onlypaths";
					extract_col_from_file($outlist_f,$cfg_hash->{'outlist_path_pos'},$only_paths_f);
					my @paths = list_to_array($only_paths_f,'NO_NEW_LINE');

					#Create a new folder
					print_and_log("Creating folder $analysis_fold_st...\n",$log_file);
					#mkdir $analysis_fold_st or die "ERROR: can't create folder $analysis_fold_st. Check permissions. \n";

					#Check if directory exists, otherwise it creates it
					unless(-d $analysis_fold_st){
						print_and_log("Creating folder $analysis_fold_st...\n",$log_file);
						mkdir $analysis_fold_st or die "ERROR: can't create folder $analysis_fold_st. Check permissions. \n";
					}	
														
					#For each file path into the outlist file... 
					foreach my $path ( @paths ){
						print_and_log("Moving $path to $analysis_fold_st..\n",$log_file);
						if ( -e $path){
							#Check if it is not a symbolink link and if extists
							if ( ! (-l $path) ){
								#Move the file into the storage directory
								move($path,$analysis_fold_st) or print "ERROR: unable to move $path in $analysis_fold_st\n";
							}else{
								print_and_log("ERROR: The file $path is already a symbolink link! Transfer was not performed...\n",$log_file);
							}
						}else{
								print_and_log("ERROR: The file $path does not exist..\n",$log_file);
							}
					}
					
					############################
					##HTML PAGE PRINT
					############################
					#Put it into an hash and change the folder where they are located
					outlist_to_hash_2($outlist_f,\$paths_hash,$analysis_fold_st);
			
					#Create a single html file containing all the links (OUT,QC, COVERAGE)
					print_and_log("Creating a new HTML page $analysis_fold_st/".$cfg_hash->{'html_index_page'}."..\n",$log_file);
					my $html_string = paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_output_page'},$paths_hash,$cfg_hash->{"work_fold"}); 
					$html_string .= paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_qc_page'},$paths_hash,$cfg_hash->{"work_fold"}); 
					$html_string .= paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_cov_page'},$paths_hash,$cfg_hash->{"work_fold"}); 
					#finally create an html file with this page
					#Arrays for info of pages
					my @descText = ($cfg_hash->{'html_index_desc'});
					my @links = ($cfg_hash->{'html_index_page'});
					html_page($cfg_hash,'none',\@descText,0,\@links,undef,$html_string,$analysis_fold_st);
								
					#Compress and remove
					print_and_log("Compressing the new folder $analysis_fold_st..\n",$log_file);
				#	compress_folder($storage_dir,$analysis_name,$analysis_name.".".$cfg_hash->{'targz_ext'});		
					#gzip_folder($storage_dir,$analysis_name,$analysis_name.".".$cfg_hash->{'zip_ext'});		
					print_and_log("Removing the folder $analysis_fold_st..\n",$log_file);	
					#delete_directory($analysis_fold_st);
					print_and_log("Removing the folder ".$working_folder."/".$analysis_name."..\n",$log_file);	
					delete_directory($working_folder."/".$analysis_name);
										
					##################################################
					## UPDATE STORAGE STATUS
					##################################################
					
					#Update the analysis table with the fact that the analysis was stored (stored=1). 
					my $stored = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'},
										$cfg_hash->{'db_analysis_id'},$analysisid);

					my $fields = $cfg_hash->{'db_analysis_stored'};
					my $values = 1;
					my $fields2 = $cfg_hash->{'db_analysis_id'};
					my $values2 = $analysisid;
					update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
					##################################################		
						
				}else{print_and_log("WARNING: The folder $analysis_fold_st already exists. Check if the transfer has been already done...\n",$log_file);}
			}else{print_and_log("Analsysis $analysisid is not a db id\n",$log_file);}			
		}else{print_and_log("The analysis $analysisid has been already stored. Use --force_run if you know what you are doing..\n",$log_file);}																	

	}
}	

=head2 outlist_to_hash

 Title   : outlist_to_hash
 Usage   : outlist_to_hash( config_file => the config hash );

 Function: Reads the outlist file and builds an hash for straightforward access
						The hash structure is: PAGE->SAMPLE->READFILE->DESCRIPTION->PATH
					
					So that, later, for each page, for each sample and for each readfile can be print a link with a description	
					
  
 Returns : the hash given in input is filled..
 
=cut
sub outlist_to_hash_2{
	my $outlist_f = shift;
	my ($paths_hash) = shift;
	my $new_folder = shift;
	
	
	my $samid_fld_ind = 0;
	my $rdf_fld_ind = 1;
	my $page_fld_ind = 2;
	my $desc_fld_ind = 3;
	my $path_fld_ind = 4;
	
	my $final_path = "";		
 
 #Createthe hash with the information contained
 open (FILE,"<$outlist_f") or die "ERROR [$!]: can't open file: $outlist_f. Check permissions.\n"; 
 while (my $line = <FILE>){
	chop($line);
	my @fields = split("\t",$line);

	if (defined $new_folder){
		$final_path =  $new_folder."/".extract_name($fields[$path_fld_ind],0);
	}else{
		$final_path = $fields[$path_fld_ind];
	}
	#The hash contains data separated by page
	$$paths_hash->{$fields[$page_fld_ind]}->{$fields[$samid_fld_ind]}->{$fields[$rdf_fld_ind]}->{$fields[$desc_fld_ind]}->{'path'} = $final_path;

 }
 close FILE;
}



=head2 compress_data_to_store

 Title   : compress_data_to_store
 Usage   : compress_data_to_store(  );

 Function:  compress_data_to_store
						
 Returns : nothing

=cut
sub compress_data_to_store{
	my $paths_hash = shift;
	my $infolder = shift;
	my $log_file = shift;
	
	my $storefolder = "/compressed";
	create_folder($infolder."/".$storefolder);
	
	while (my $page = each %{ $paths_hash } ) {
		while (my ($sample_name, $value) = each %{ $paths_hash->{$page} } ) {
			while (my ($readf_name, $value2) = each %{ $paths_hash->{$page}->{$sample_name}} ) {
				while (my ($desc, $value3) = each %{ $paths_hash->{$page}->{$sample_name}->{$readf_name}} ) {
					#the description indicates if the file must be printed or not [WEB]
					print_and_log( "Checking $desc if contains [STORE]\n ",$log_file);#DEBUGCODE
					if ( $desc =~ /STORE/ ){
						#Extract the path and put the file into the folder
						my $f_path = $paths_hash->{$page}->{$sample_name}->{$readf_name}->{$desc}->{'path'};
						if (defined $f_path) {
							#If it is a directory use copy_folder
							if ( -d $f_path ){
								print_and_log( "Copying folder $f_path to $infolder/$storefolder ",$log_file);#DEBUGCODE								
								copy_folder($f_path,$infolder."/".$storefolder);			 
							}else{
								#Move file to storefolder
								print_and_log( "Copying $f_path to $infolder/$storefolder ",$log_file);#DEBUGCODE
								copy ($f_path,$infolder."/".$storefolder) or print "ERROR: Unable to copy $f_path to $infolder/$storefolder ";
							}
						}	
					}

				}
			}
		}	
	
	}
	#Tar all the files into the folder using th esystem command
	tar_folder_sys($infolder,$storefolder,"$storefolder.tar");
	return "$infolder/$storefolder.tar"
}


=head2 get_flowcell_and_lane

 Title   : get_flowcell_and_lane
 Usage   : get_flowcell_and_lane(  );

 Function:  get_flowcell_and_lane
						
 Returns : nothing

=cut
sub get_flowcell_and_lane{
	my $cfg_hash = shift;
	my $fq = shift;
	my $log_file = shift;
	
	my $infos = "";
	 
	if ( $fq =~ '\.gz'){
		$infos = `zcat $fq |  head -n1 | sed 's/^@//' | cut -f1 -d' '`;
	}else{
		$infos = `head -n1 $fq | sed 's/^@//' | cut -f1 -d' '`;
	}
	my @fields = split (":",$infos);
	
	
	my $flowcell = $fields[$cfg_hash->{'fq_flowcell_pos'}];
	my $lane = $fields[$cfg_hash->{'fq_lane_pos'}];
	my $instrument = $fields[$cfg_hash->{'fq_instrument_pos'}];
	
	return $instrument,$flowcell,$lane;
}
	
=head2 remove_analyses

 Title   : remove_analyses
 Usage   : remove_analyses(  );

 Function:  removes an analysis from the database given a list of group ids
						
 Returns : nothing

=cut
sub remove_analyses{
	my $cfg_hash = shift;
	my $log_file = shift;
	my $an_to_rem = shift;
	
	print_and_log("Removing an analysis from the database..\n",$log_file);
	my $sep = $cfg_hash->{'parameters_sep'};
	my @ans_to_rem = split($sep,separate_input_ids($an_to_rem,$sep));

	foreach my $an_to_rem (@ans_to_rem){
		#If the parameter is a number then proceed
		if ( correct_type($an_to_rem,"positiveint") ){					
			#Obtain the group name from the database given the group id
			my $group_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'}
																			,$cfg_hash->{'db_analysis_id'},$an_to_rem);
			#Obtain the run id from the database given the group id
			my $analysis_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sample_sheet_id'}
																			,$cfg_hash->{'db_analysis_id'},$an_to_rem);
											
			print_and_log("Removing analysis $group_name from database (id: $an_to_rem):\n",$log_file);
			#delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				#					$cfg_hash->{'db_pass'},$cfg_hash->{'db_samples_hpo_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			#print_and_log($cfg_hash->{'db_samples_hpo_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_genotype_sample_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_genotype_sample_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_var_statistics_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_var_statistics_table'}."...\n",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_readf_table'}."...",$log_file);					
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_sample_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_analyses_table'}."...",$log_file);
			delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$an_to_rem);
			print_and_log($cfg_hash->{'db_analyses_table'}."...",$log_file);
			
			#If the run was only from this analysis remove also the sample_sheet name
			#Get the group ids involved for the analysis
			my $query = "SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".
									$cfg_hash->{'db_sample_sheet_id'}."=$analysis_id;";	
			print_and_log( "Executing: $query\n",$log_file);
			my $group_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_analysis_id'});
			
			my $num_groups = scalar(keys %{$group_ids});
			if ( $num_groups == 0){
				print_and_log("Last ".$cfg_hash->{'db_analysis_id'}." for ".$cfg_hash->{'db_sample_sheet_id'}.": $analysis_id. Removing from ".$cfg_hash->{'db_sample_sheets_table'}."...\n",$log_file);				
				delete_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_id'},$analysis_id);
			}
		}else{print_and_log($cfg_hash->{'db_analysis_id'}." $an_to_rem is not present in db ".$cfg_hash->{'db_name'}."\n",$log_file);}
	}
}

=head2 show_analyses_from_sample_sheet

 Title   : show_analyses_from_sample_sheet
 Usage   : show_analyses_from_sample_sheet(  );

 Function:  shows groups related to a specific analysis name
						
 Returns : nothing

=cut
sub show_analyses_from_sample_sheet{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	
	print_and_log("Showing groups from analysis $analysis_id..\n",$log_file);
	my @res = ();

	if ( !correct_type($analysis_id,"positiveint") ){
		#Obtain the group name from the database given the group id
		$analysis_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_sheets_table'},$cfg_hash->{'db_sample_sheet_id'}
								,$cfg_hash->{'db_sample_sheet_name'},"'".$analysis_id."'");
	}
	my $fields = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_analysis_name'};
	
	if ( correct_type($analysis_id,"positiveint") ){
		#Prints the selected fields of the groups table using the analysis id as a key value
		db_print_selected_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$fields,$cfg_hash->{'db_analyses_table'},
									$cfg_hash->{'db_sample_sheet_id'},$analysis_id);
	}else{
		log_and_exit("ERROR: Cannot find any analysis with your input ...\n",$log_file);	
	}
}

=head2 show_samples_from_group

 Title   : show_samples_from_group
 Usage   : show_samples_from_group(  );

 Function:  shows samples related to a specific group name
						
 Returns : nothing

=cut
sub show_samples_from_group{
	my $cfg_hash = shift;
	my $group_id = shift;
	my $log_file = shift;
		
	print_and_log("Showing samples from group $group_id..\n",$log_file);
	my @res = ();

	if ( !correct_type($group_id,"positiveint") ){
		#Obtain the group name from the database given the group id
		$group_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'}
								,$cfg_hash->{'db_analysis_name'},"'".$group_id."'");
	}
	my $fields = $cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_sample_name'};
	
	if ( correct_type($group_id,"positiveint") ){
		#Prints the selected fields of the groups table using the analysis id as a key value
		db_print_selected_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$fields,$cfg_hash->{'db_sample_table'},
									$cfg_hash->{'db_analysis_id'},$group_id);
	}else{
		log_and_exit("ERROR: Cannot find any analysis with your input ...\n",$log_file);	
	}
}

=head2 fields_2_db

 Title   : fields_2_db
 Usage   : fields_2_db(  );

 Function:  given a file with the first row containing the fields of
						the table and the second row contains the values puts this 
						values into the given table. Uses an id which can be an analysisid
						
 Returns : returns 1 else 0. -1 in case of error.

=cut
sub fields_2_db{
	my $cfg_hash = shift;
	my $dbidfield = shift;
	my $dbid = shift;
	my $file = shift;
	my $log_file = shift;
	
	#print_and_log("Checking $file.... \n",$log_file);#DEBUGCODE
	open (FILE,"<$file") or die "ERROR [$!]: can't open file: $file. Check permissions.\n"; 
	my $numRow = -1;
	 
	my @hFields = ();
	my @vFields = ();
	my $dbtable = "";
	while (my $line = <FILE>){
		chop($line);
		#A block of information starts with a set of # and a string
		if ($line =~ /#+DB/){
			$numRow = 0;
			my @fields = split(" ",$line);
			#Get the table involved
			$dbtable = $cfg_hash->{$fields[1]};	
			print_and_log("The table to use is $dbtable\n",$log_file);						
			
		}
		elsif ($line =~ /#+/) {
			#If a block is just finished, fill the database
			if ($numRow == 2){
				#Check the dimensions that must be the same
				if ( scalar(@hFields) ==  scalar(@vFields) ){
					#With a for loop put into a string all the fields and values from the file
					my $numField = 0;
					my $fields = "$dbidfield,";
					my $values = "$dbid,";
					foreach my $hField (@hFields){
						
						$fields .= "$hField,";
						$values .= $vFields[$numField].",";
						$numField++;
					}
					chop($fields);
					chop($values);
					
					print_and_log("Inserting into $dbtable. dbidfield: $dbidfield, fields: $fields, values:$values, dbidfield: $dbidfield, dbid: $dbid\n",$log_file);
					#insert into table if the element do not exists, otherwise update it			
					my ($new_id,$exist_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$dbtable,$dbidfield,$fields,$values,$dbidfield,$dbid);

					#If the id is new ok, otherwise update
					if ( $new_id >= 0){
						print_and_log( "Inserted $fields -> $values into $dbtable for $dbidfield = $dbid\n",$log_file);#DEBUGCODE
					}elsif ( $exist_id >= 0) {
						update_table_2($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
												$dbtable,$fields,$values,$dbidfield,$dbid);
					}else{
						log_and_exit("ERROR: during the import of information into $dbtable for $dbidfield : $dbid. Cannot understand what is happen with the db...\n",$log_file);
					}		

				}else{
					print_and_log("ERROR: rows in file $file differ. Different number of fields ".scalar(@hFields)." ".scalar(@vFields)."\n",$log_file);						
				}
			}
			$numRow = -1;
		}
		else{
			#If we are in this else means that we are in the block
			if ($numRow == 0 or $numRow == 1){
				#Get the header fields
				if ($numRow == 0){
					@hFields = split("\t",$line);		
				}
				#Get the values
				if ($numRow == 1){
					@vFields = split("\t",$line);		
				}
				$numRow++;				
			}
		}
	}
	close FILE; 
}

=head2 good_reads_quality

 Title   : good_reads_quality
 Usage   : good_reads_quality(  );

 Function:  extracts the report from the zipped file from FastQC and checks the content
						of the fastqc_data.txt file. 
						Then basically gets the average of the mean sequence qualities and if it is
						higher than a threshold, returns 1 else 0. -1 in case of error.
						
 Returns : returns 1 else 0. -1 in case of error.

=cut
sub good_reads_quality{
	my $cfg_hash = shift;
	my $qc_fq1_arch = shift;
	my $qc_fq2_arch = shift;
	my $out_folder = shift;
	my $log_file = shift;
	
	my $good = 0;
	
	#Extract the files
	extract_archive($qc_fq1_arch,$out_folder);
	
	my $arch_fold = extract_name($qc_fq1_arch,1);
	my $report_fq = $out_folder."/".$arch_fold."/".$cfg_hash->{'fastqc_report'};
	
 #Create
 #print_and_log("Checking $report_fq.... \n",$log_file);#DEBUGCODE
 open (FILE,"<$report_fq") or die "ERROR [$!]: can't open file: $report_fq. Check permissions.\n"; 
 while (my $line = <FILE>){
	chop($line);
	if ($line =~ />>Per base sequence quality/){
		my @fields = split("\t",$line);
		#print_and_log("quality check: ".$fields[1].".... \n",$log_file);#DEBUGCODE
		if ( $fields[1] eq 'pass'){
				$good = 0.5;
		}		
		last;
	}

 }
 close FILE; 
 
 #If fastq1 is good check fastq2
 if ($good > 0){
		#Extract the files
		extract_archive($qc_fq2_arch,$out_folder);
		$arch_fold = extract_name($qc_fq2_arch,1);
		$report_fq = $out_folder."/".$arch_fold."/".$cfg_hash->{'fastqc_report'};
		#print_and_log("Checking $report_fq.... \n",$log_file);#DEBUGCODE
	 #Create
	 open (FILE,"<$report_fq") or die "ERROR [$!]: can't open file: $report_fq. Check permissions.\n"; 
	 while (my $line = <FILE>){
		chop($line);
		if ($line =~ />>Per base sequence quality/){
			my @fields = split("\t",$line);
			#print_and_log("quality check: ".$fields[1].".... \n",$log_file);#DEBUGCODE
			if ( $fields[1] eq 'pass'){
					$good = 1;
			}				
		}
	}
	close FILE;
 }
 
 return $good;
}

=head2 show_analysesids_from_analysesnames

 Title   : show_analysesids_from_analysesnames
 Usage   : show_analysesids_from_analysesnames(  );

 Function:  shows groups ids using the groupnames
						
 Returns : nothing

=cut
sub show_analysesids_from_analysesnames{
	my $cfg_hash = shift;
	my $analysis_names = shift;
	my $log_file = shift;
		
	print "Showing analysis identifiers for: $analysis_names..\n";
	my @analyses_names = split(",",$analysis_names);

	foreach my $analysis_name (@analyses_names){
		#Obtain the group name from the database given the group id
		my $analysis_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
								$cfg_hash->{'db_analysis_name'},"'".$analysis_name."'");		
		
		if ( $analysis_id >= 0){
			print "$analysis_name : $analysis_id\n";
		}else{
			print "Mmmmh... could not find any analysis named $analysis_name...\n";
		}
	}
}

=head2 show_samplesids_from_samplesnames

 Title   : show_samplesids_from_samplesnames
 Usage   : show_samplesids_from_samplesnames(  );

 Function:  shows sample ids using the sample names
						
 Returns : nothing

=cut
sub show_samplesids_from_samplesnames{
	my $cfg_hash = shift;
	my $samples_names = shift;
	my $log_file = shift;
		
	print "Showing sample identifiers for: $samples_names..\n";
	my @samples_names = split(",",$samples_names);

	foreach my $sample_name (@samples_names){
		#Obtain the sample id from the database given the sample name id
		my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
								$cfg_hash->{'db_sample_name'},"'".$sample_name."'");		
		
		if ( $sample_id >= 0){
			print "$sample_name : $sample_id\n";
		}else{
			print "Mmmmh... could not find any sample named $sample_name...\n";
		}
	}
}


=head2 separate_bed_perchr

 Title   : separate_bed_perchr
 Usage   : separate_bed_perchr(  );

 Function:  this subroutine takes in input the target file given by the user and
						generates a file for each chromosome 
						
 Returns : nothing

=cut
sub separate_bed_perchr{
	my $target_bed = shift;
	my $dest_folder = shift;

	my $workDir = getcwd;
	chdir $dest_folder;
	
	my $command = "awk '".'{print $0 >> $1".bed"}'."' $target_bed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#Return to previous directory
	chdir $workDir;
}



    
=head2 initialize_folders

 Title   : initialize_folders
 Usage   : initialize_folders(  );

 Function:  this subroutine initializes the basic folders. If the folders.txt file is not present
            means that the program has not been installed in the folder used.
 Returns : nothing

=cut
sub initialize_folders{
  my $foldersFile = shift;

	#print "Opening $foldersFile..\n";
  #If the file with folders references is installed the program can start
  open (FOLD, $foldersFile) or die " ERROR: opening $foldersFile. Please install the program before.\n";
  my $line = <FOLD>;

	#print "Extracting folders path..\n";
  #Extracts folders paths
  my @folders = split(" ", $line);

  #The working folder is the first path
  my $workingFolder = $folders[0];

  #Cancel the final slash..if it is there it is removed
  if($workingFolder =~ /\/$/){
    chop($workingFolder);
  }

  my $programFolder = $folders[1];
  #Cancel the final slash..if it is there it is removed
  if($programFolder =~ /\/$/){
    chop($programFolder);
  }
  close(FOLD);

  return $workingFolder,$programFolder;
}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configFile -> file with the user configuration
                              - variablesFile -> the path to a file with all variables written
                              - lineToCheck -> line in the variables file to be used
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.

 Returns : nothing

=cut
sub checkConfigVariables {
  my $configFile = shift;
  my $variablesFile = shift;
  my $lineToCheck = shift;

  my $hashToCheck;

  if (! open(VARF,"<$variablesFile")){ die "ERROR: Failure opening '$variablesFile'. Your program version is corrupted - $!";}
  if (! open(CFILE,"<$configFile")){ die "ERROR: Cannot find '$configFile' - Your program version is corrupted - $!";}

  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CFILE>){
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashToCheck->{$1} = "OK";
    }
  }
  close(CFILE);

  my @confVars = ();

  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $lines=0;
  #For each line of the variables file
  while (my $line = <VARF>){
    $line =~ s/\n//;#Remove \n in the end of the line

    #get the variables in the line
    my @variables = split (/;/,$line);

    $lines++;
    #For each of the variables in the variables file
    foreach my $var (@variables){
      if ($lines == $lineToCheck){
        #print "program Variable: $var - - value: ".$hashToCheck->{$var}."\n";#DEBUGCODE
        #put the variable inside an array for program config
        push (@confVars, $var);
        if( !(defined($hashToCheck->{$var})) ){

          die "ERROR: in $configFile variable $var is missing. Please check the file. Closing...\n ";
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      }
    }
  }

  #print_array(\@allVars);
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashToCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @confVars )){
       die "ERROR: Variable $key is in the config file $configFile and not in $variablesFile file. This is completely wrong. Program will not work...\n ";
     }
  }
  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}



=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                             - configHash = the pointer to the hash to be filled
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
 Returns : nothing

=cut
sub configFile2Hash{
  my $configFilePath=shift;
  my ($configHash) = shift;

  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or die "ERROR: The file $configFilePath doesn't exists. The program will exit..\n";
  while (my $line = <configFile>){
		#This IF is useful if you want to put a description above in the text file. Delimit it with a set of hashtags
    if ($line =~ /#########/){
      $start = 1;
    }
   # if( ($line =~ /(\S+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /#/) ){
   if( ($line =~ /(\w+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /^#/) ){
      if ( $2 ne ''){
			  $$configHash->{$1} = $2;
				#annoPrint ("$1 = $2\n") ;#DEBUGCODE
			}
    #}elsif ( $line =~ /(\S+)\s*=/ ){
    }elsif ( $line =~ /(\w+)\s*=$/ ){
			delete($$configHash->{$1});
		}
  }
  close(configFile);
  #print Dumper\$configHash; #DEBUGCODE
}


=head2 hash2ConfigFile

 Title   : hash2ConfigFile
 Usage   : hash2ConfigFile( - configFilePath = path of the config file
                             - configHash = the pointer to the hash to be filled
                               );

 Function: puts variables from an hash inside a file to make a configuration file
 Returns : nothing
=cut
sub hash2ConfigFile{
  my $configFilePath=shift;
  my ($configHash) = shift;
  my $start_string = "################\n";
  my $start = 0;

  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,">$configFilePath") or die "ERROR: The file $configFilePath doesn't exists. The program will exit..\n";
  print configFile $start_string;
  foreach my $key ( keys %{$$configHash}){
    print configFile "$key = ".$$configHash->{$key}."\n";
  }
  close(configFile);
        #print Dumper\$configHash; #DEBUGCODE
}

=head2 checkVariable

 Title   : checkVariable
 Usage   : checkVariable(  - var -> value of the variable to check
                           - name -> name of the variable
                           - sentence -> something to write to the user in case of error)

 Function:  this subroutine checks if a variable has the YES or NO value.
						It takes in input a "yes and a "no sentence" which is used whether 
						the value corresponds

 Returns : nothing

=cut
sub checkVariable {
  my $var = shift;
  my $name = shift;
  my $yesSentence = shift;
  my $log_file = shift;

  if ( ($var ne 'YES') and ($var ne 'NO')){
    log_and_exit("ERROR: Check variable $name in the config file. $var is a wrong value!\n",$log_file);
  }elsif ($var eq 'YES'){
    print_and_log($yesSentence,$log_file);
  }
}


=head2 correct_type
Title  : correct_type
 Usage  : correct_type( -number => 'a number to check',
                      -typeWanted => 'the expected type');

 Function:       Check if the number is of the type given in input

  Returns 1 if true

=cut
sub correct_type {
  my $number = shift;
  my $typeWanted = shift;

  my $type = "";
  my $ret = 0;

  if ($typeWanted eq "positiveint"){
    #contains onli digits
    if ( $number =~ m/^\d+$/) {
      $ret=1;
    }
  }elsif ($typeWanted eq "real"){
    #Contains digits starting with (+-), separated by '.' and can be added the 'E-x'
    if ( $number =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ){
        $ret=1;
    }
  }
  return $ret;
}

=head2 check_input_fasta_file

 Title   : check_input_fasta_file
 Usage   : check_input_fasta_file(  - fastaSeqs = the name of the fasta
                                                                                                                                       $
                                                                                                                                       $
                                                                                                                                       $
                                                                                                                                       $
                               );

 Function:   check if the query FASTA file is well written and is in the working folder. Behaviour of the function is different
          if user is calling a saved session. The file is moved soon in the session folder.


 Returns : In the config hash will be returned:
                                               	- $$configHash->{$name.'Present'} => YES if the file is ok, NO otherwise
                                                - $$configHash->{$name.'Alphabet'} => alphabet type
                                                - $$configHash->{$name.'totSequences'} => total number of sequences in the file

=cut
sub check_input_fasta_file {

#CONTROLS IF THE QUERY EXISTS AND THE NAME IS .fa or .fasta
  my $fastaSeqs = shift;
  my ($configHash) = shift;#An hash with variables to be used in the main program
  my $newSession = shift;#Indicates if the session is loaded
  my $workingFolder = shift;#The folder where the user is working
  my $sessionFolder = shift;#The folder used for the output

  my $name = '';

  #The name has to be without dots inside. Only one dot can separate the extension from the name
  my $dots=0;
  my @line = split (//,$fastaSeqs);
  foreach my $char(@line){
      if ($char eq '.'){
          $dots++;
      }
  }
  if ($dots>1){
      die "Please change the name of your file removing all internal dots (not the one used for the extension). Exiting...\n";
  }
                                                                                                                                        
  #annoPrint ("Fasta name read: ".$fastaSeqs."\n";
  my @name = split(/\./,$fastaSeqs);
  #annoPrint ("Fasta name read: ".$name[0]."\n";
  #Here is created the name to use for all the files of output by extracting it from the fasta file name
  #It cannot be longer than 50 chars
  if (@name > 1){
    #$$configHash->{'name'} = extract_name($$configHash->{'sessionFolder'},0);
    $name = $name[0];
    #Checks for permitted characters
    if ( (length($name) > 50) # or !($$configHash->{'fastaSeqs'} =~ /(\.fasta?)|(\.fa?)/)
            or ($fastaSeqs !~ /^[A-za-z0-9\_\-]+\.(fa|fasta)$/i)){
            die "$fastaSeqs is not a correct FASTA file name. Allowed characters [A-Z,a-z,0-9,_,-]. Allowed extensions [.fasta,.fa]. Max length: 50 chars\n";
    }else{
      #If session is a new one
      #Fasta file must stay only in the working folder and it will be shifted. If it is not there,then only db creation can be done or $
      if ( $newSession == 1 ){
       my $fastaPath = $workingFolder.'/'.$fastaSeqs;
       print "Path to control for fasta: $fastaPath\n";
        #The Sequences file can also stay already in the session folder created by the user in a previous execution
        if (-e $fastaPath) {
          print "Sequences found in your working folder. Checking...\n";
          #After the name all the fasta is checked to see if it respects the standards
                                        my $warnings = checkFastaFormat($fastaPath,$$configHash->{'maxFastaSeqLen'},$$configHash->{'nuclIUPAC'},$$configHash->{'protIUPAC'});
                                        if($warnings > 0){
                                               	annoDie("Fasta check failed! The transcriptome you are using have the listed problems ".
                                                " which can cause stochastic error in BLAST programs execution. ".
                                                "Please try to correct the problems. Exiting...\n");
                                       }
          print "Copying the fasta in session folder...\n";
          #move($fastaPath,$sessionFolder);
          copy($fastaPath,$sessionFolder) or annoDie("Cannot copy $fastaPath in $sessionFolder");
          $$configHash->{$name.'Present'} = 'YES';
        }#else{
          ##This else happens when sequences file is not there. Then, if the user chose to execute analysis or to print output
          ## the program has to die.
          #if ( ($$configHash->{'doExecutePrograms'} eq 'YES') or ($$configHash->{'doBuildOutput'} eq 'YES')
            #or ($$configHash->{'extractStatistics'} eq 'YES')){
              #die "No sequences file found in $fastaPath. You'll cannot execute analysis...\n";
            #}
          #print "No sequences file found in  $fastaPath. But you are only creating a DB...\n";
          #$$configHash->{$name.'Present'} = 'NO';
        #}
      }else{#If session is loaded then the fasta must stay only in the session folder. No other fasta can be shifted there. Die otherwi$
        if (-e $sessionFolder.'/'.$fastaSeqs) {
          print "Sequences found in your session folder. Checking it...\n";
          $$configHash->{$name.'Present'} = 'YES';
                                #After the name all the fasta is checked to see if it respects the standards
					my $warnings = checkFastaFormat($sessionFolder.'/'.$$configHash->{'fastaSeqs'},$$configHash->{'maxFastaSeqLen'},$$configHash->{'nuclIUPAC'},$$configHash->{'protIUPAC'});
                                        if($warnings > 0){
                                                die "Fasta check failed! The transcriptome you are using have the listed problems ".
                                               	" which can cause stochastic error in BLAST programs execution. ".
                                               	"Please try to correct the problems. Exiting...\n";
                                        }
        }else{
          die "Unable to find ".$$configHash->{'fastaSeqs'}." in $sessionFolder. The program will exit and you can figure out why.".
                    " \nThere are 2 possible reasons: \n\t 1. you written a wrong fasta name in the configuration file.".
                    "\n\t 2. the name in the config file is right but not the folder you chose\n Bye!\n ";
        }
      }

      if ($$configHash->{$name.'Present'} eq 'YES'){

        #Here we control what kind of alphabet has the sequence and store the value in a variable
        $$configHash->{$name.'Alphabet'} = detect_fasta_type($sessionFolder.'/'.$fastaSeqs);
        print "The program detected that your sequences are ";
        if ( ($$configHash->{$name.'Alphabet'} eq 'dna') or ($$configHash->{$name.'Alphabet'} eq 'rna') ){
            print "transcripts\n";
        }else{
            print $$configHash->{$name.'Alphabet'}."\n";
         }
        #Here we count and store the number of sequences in the fasta file
        my $seqsPath = $sessionFolder.'/'.$fastaSeqs;
        my $totSequences = count_sequences($seqsPath); #Extract number of sequences from the query file
        print "Number of sequences: $totSequences";
        $$configHash->{$name.'totSequences'} = $totSequences;
      }
    }
  }else {
    die "The name of file with transcripts needs to have a .fa or .fasta extension. Exiting...\n";
  }
}

=head2 try_exec_command

 Title   : try_exec_command
 Usage   : try_exec_command( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.

 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_command{
    my $command = shift;

    my $maxTimes = 5;
    my $success = -1;
    my $timesCount = 0;

    while ($success == -1 and $timesCount < $maxTimes){
        if ( (system $command) == 0) {
          $success = 1;
        }
        else{
         if ($? == -1) {
              print "failed to execute: $!\n";
          }
          elsif ($? & 127) {
              printf "child died with signal %d, %s coredump\n",
                  ($? & 127),  ($? & 128) ? 'with' : 'without';
          }
          elsif ($? == 35 ) {
              printf "child exited normally\n",
              $success = 1;   
          }          
          else {
              printf "child exited with value %d\n", $? >> 8;
          }
         $timesCount++;
        }
    }
    return $success;
}


=head2 build_input_name_from_executed

 Title   : build_input_name_from_executed
 Usage   : build_input_name_from_executed(   );

 Function: Takes in input the hash with parameters from the sample table
					and the step in which we are and uses the array steps_array
					to construct the name depending by the mode.
					NB: The last step will not be included in the name
					
 Returns : a string with the name of the file to be used as input in some program
=cut
sub build_input_name_from_executed{
	my $params = shift;
	my $step = shift;
	my $in = shift;
	my $steps_array = shift;
	
	my $out = $in;
	my $stop = 0;
	for ( my $i = 0; $i < scalar(@$steps_array) and $stop == 0; $i++){
		#If from the database we see that the specific step has been run add its suffix
		if ( $steps_array->[$i] eq $step ) {
			$stop = 1;
		}elsif ( $params->{$steps_array->[$i]} == 1){
			$out .= '_'.$steps_array->[$i];
		}
	}
	return $out;
}

=head2 get_output_name_from_executed

 Title   : get_output_name_from_executed
 Usage   : get_output_name_from_executed(   );

 Function: Takes in input the hash with parameters from the sample table
					and the step in which we are and uses the array steps_array
					to obtain the name of an output file checking which steps has been executed
					NB: The last step will be included in the name
					
					Returns : a string with the name of the file to be used as input in some program
=cut
sub get_output_name_from_executed{
	my $params = shift;
	my $last_step = shift;
	my $in = shift;
	my $steps_array = shift;
	
	my $out = $in;
	my $stop = 0;
	for ( my $i = 0; $i < scalar(@$steps_array) and $stop == 0; $i++){
		#If from the database we see that the specific last_step has been run add its suffix
		if ( $steps_array->[$i] eq $last_step ) {
			$stop = 1;
		}
		if ( $params->{$steps_array->[$i]} == 1){
			$out .= '_'.$steps_array->[$i];
		}
	}
	return $out;
}

=head2 alter_job

 Title   : alter_job
 Usage   : alter_job(   );

 Function:  uses qalter to modify a job given in input
 
 Returns : nothing

=cut
sub alter_job{
	my $job_id = shift;
	my $string = shift;
	my $log_file = shift;
	
	my $command = "qalter $string $job_id";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}

=head2 JOB_is_active_queue

 Title   : JOB_is_active_queue
 Usage   : JOB_is_active_queue(   );

 Function: goes iteratively through a job queue and checks if the 
			queue is completely active, hence there are no jobs waiting
			for exited jobs (zombie status)
 
 Returns : 1 if the queue of waiters is active
 
=cut
sub JOB_is_active_queue{
	my $jobid = shift;
	
	my $retval = 1;
	#Get the job dependencies. If it has no jobs to wait return back
	my $jobdeps = JOB_get_dependencies($jobid);
	if ( $jobdeps !~ /finished/ ){
		#Separate the type of waiting
		my @waits =  split(",",$jobdeps);
		foreach my $wait ( @waits){
			#Separate the elements using :
			my @elems = split(":",$jobdeps);
			my $wtype = shift(@elems);
			#If there is afterok means there is a list to check
			if ( $wtype eq 'afterok'){
				#until there is a job to wait in the tree, go ahed and check that job
				foreach my $job2wait (@elems){
					if ( JOB_is_active_queue($job2wait) == 0){
						$retval = 0;
					}
				} 
			}
		}
	}else{$retval = 0;}
	
	return $retval;
}

=head2 JOB_get_dependencies

 Title   : JOB_get_dependencies
 Usage   : JOB_get_dependencies(   );

 Function:  uses qstat -f to get dependencies of a job given in input
 
 Returns : a string with the list of dependencies
 
=cut
sub JOB_get_dependencies{
	my $job_id = shift;
	my $log_file = shift;
	
	#Here we take the current time 
	my $time = scalar(localtime);
	$time =~ s/ /_/g;
	$time =~ s/:/-/g;

	#Printing qstat result in a file
	my $qstat_out = "qstat_out_$time.txt";
	#my $deps = `qstat -f $job_id | grep 'depend =' | cut -f2 -d'=' | sed 's/ //' `;
	
	my $command = "qstat -f $job_id > $qstat_out";
	try_exec_command($command) or die "Unable to execute command: $command\n";
		
	open (QSTAT_OUT,"<$qstat_out") or die "ERROR: Cannot open $qstat_out. The program will exit..\n";
	my $start = 0;
	my $deps = "";
	while (my $line = <QSTAT_OUT>){
		chop($line);
		#If the depend was found..
		if ($start == 1){
			if ( defined $log_file){
				print_and_log( "reading nex line of depends $line...\n",$log_file);
			}
			#print $line."\n";
			#Be careful dependences are not finished with a new element
			if( $line != /=/ ){
				
				$line =~ s/\s//g;
				$deps .= $line;
				if ( defined $log_file){
					print_and_log( "adding:'line'. Obtaning:'$deps'\n",$log_file);
				}
			}else{
				$start = 0;
			}
		}
		#Find the line with depend = 
		if( $line =~ /depend =/ ){
			if ( defined $log_file){
				print_and_log( "depends are found...\n",$log_file);
			}	
			#print $line."\n";
			my @fields = split("=",$line);
			my $stdeps = $fields[1];
			$stdeps =~ s/ //g;
			$deps .= $stdeps;
			$start = 1;
		}

	}
	close(QSTAT_OUT);
		
	return $deps;
}

=head2 JOB_get_status_field

 Title   : JOB_get_status_field
 Usage   : JOB_get_status_field(   );

 Function:  uses qstat to get status of a job given in input
 
 Returns : a string with status which can be
		
			qstat: 120294.node001 Job has finished, use -x or -H to obtain historical job information
			if the job has finished
=cut
sub JOB_get_status_field{
	my $job_id = shift;
	my $status_field = shift;
	my $log_file = shift;
	

	my $status = "";
	#Here we take the current time 
	my $time = scalar(localtime);
	$time =~ s/ /_/g;
	$time =~ s/:/-/g;

	#Printing qstat result in a file
	my $qstat_out = "qstat_out_$time.txt";
	my $command = "qstat -n1 -t -w $job_id  > $qstat_out";
	try_exec_command($command) or die "Unable to execute command: $command\n";

	my $start =0 ;
	#print "Looping...\n";
	open (QSTAT_OUT,"<$qstat_out") or die "ERROR: Cannot open $qstat_out. The program will exit..\n";
	while (my $line = <QSTAT_OUT>){
			if( $start == 1 ){
				#print $line."\n";
				my @fields = split(/\s+/,$line);
				$status = $fields[$status_field];
			}
			if ($line =~ /---------------/){
				$start = 1;
			}
			#Job could be finished
			if ($line =~ /finished/){
				$status = "FINISHED";	
			}
	}
	close(QSTAT_OUT);	
	delete_file($qstat_out);
	return $status;
}

=head2 run_FREEBAYES_vc

 Title   : run_FREEBAYES_vc
 Usage   : run_FREEBAYES_vc(   );

 Function: Takes in input a bam file and an output VCF and executes freebayes
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_FREEBAYES_vc{
	my $cfg_hash = shift;
	my $bam_inputs = shift;
	my $outFile = shift;
	my $analysis_id = shift;	
	my $log_file = shift;
	
	
	#Prepare parameters
	my $param_str = " ";	
	
	if ( defined $cfg_hash->{'FB_min_mapping_quality'} ){
					$param_str .= " --min-mapping-quality ".$cfg_hash->{'FB_min_mapping_quality'};
	}
	
	if ( defined $cfg_hash->{'FB_min_base_quality'} ){
					$param_str .= " --min-base-quality ".$cfg_hash->{'FB_min_base_quality'};
	}
	
	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
					$cfg_hash->{'db_analysis_id'},$analysis_id);
	if (file_not_present($target_bed) > 0 ){ print_and_log("WARNING no target bed file was usesìd using in freebayes!\n",$log_file);}
			
	if ( $target_bed ne '-'){ 
		#Get the target without the header
		my $command = " grep '^chr' $target_bed > $target_bed.temp";
		print_and_log("The command is: $command\n",$log_file);
		try_exec_command($command) or die "ERROR: cannot execute $command\n";	
					
		$param_str .= " --targets $target_bed.temp ";	
	}
	
	my @bam_files = split(",",$bam_inputs);
	$bam_inputs = "";
	foreach my $bam_file ( @bam_files){
			$bam_inputs .= " --bam $bam_file "; 
	}
		
	#Execute the command
	my $command = $cfg_hash->{'freebayes_path'}." -f ".$cfg_hash->{'hum_ref'}." $param_str $bam_inputs  > $outFile";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}


=head2 update_dependency_hash

 Title   : update_dependency_hash
 Usage   : update_dependency_hash(   );

 Function: updates the dependencies hash that is routinely loaded by the 
			job cleaner system. When the alter command is called this 
			subroutine should always be called because dependencies are changed
 
 Returns : nothing

=cut
sub update_dependencies_hash{
	my $deps_f = shift;
	my $job_name = shift;
	my $job_id = shift;
	my $new_depend = shift;
	my $log_file = shift;

	#Load the hash of dependencies
	my $deps = load_hash($deps_f);
	
	#updates the hash with the new dependencies
	$deps->{$job_name}->{'id'} = $job_id;
	$deps->{$job_name}->{'deps'} = $new_depend;
	
	#Saving the hash with the ne dependencies
	save_hash(\$deps,$deps_f);
}


=head2 separate_input_ids

 Title   : separate_input_ids
 Usage   : separate_input_ids(   );

 Function:  given ids in the form 1,2,3 or 1-4,5,6 will separate this ids
						and will return a list separated by comma
 
 Returns : nothing

=cut
sub separate_input_ids{
	my $ids = shift;
	my $sep = shift;
	
	my $final = '';
	
	my @ids = split($sep,$ids);
	#print "I am separating: $ids\n";#DEBUG_CODE
	foreach my $id (@ids){
		if ($id =~ /\-/){
				#print "There is sep with -\n";#DEBUG_CODE
				my @nums = split("-",$id);
				#print "From: ".$nums[0]." to ".$nums[1]."\n";#DEBUG_CODE
				for (my $i = $nums[0]; $i <= $nums[1]; $i++){
						$final .= "$i,";
				}
		}else{
				$final .= "$id,";
		}
	}
	chop($final);
	#print "Returning $final\n";#DEBUG_CODE
	return $final;
}


=head2 try_exec_job

 Title   : try_exec_job
 Usage   : try_exec_job(    $qsub_account => account to use for execution;
																$qsub_queue => queue to use;
																$env_vars =s> variables to pass to subject program
																$program_to_run => program to be run;
																$lnodes => number of nodes needed;
																$ppn => number of cpus needed;
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.
						
						If the variable qsub_param_for_all=YES then it takes the parameters for the resources
						request from the program_config. This means that all the tasks will ask for the same
						number of resources. Otherwise the suffix of the variable will change with the task name 
						(qc,align, refine, etc..) and the resources will be ask depending by the task 
						from the user_config.
						The following are the parameters for qsub that can be changed
							qsub_mem = 120GB
							qsub_ncpus = 12
							qsub_nodes =
							qsub_select = 1
							
							
 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_job{
	my $cfg_hash = shift;
	my $env_vars =shift;
	my $task = shift;
	my $program_path = shift;
	my $job_name = shift;
	my $dependencies = shift;
	my $std_out = shift;
	my $std_err = shift;
	my $log_file = shift;
	
	my $qsub_command;
	#Account to be used to launch job on the cluster
	if (defined $cfg_hash->{'qsub_cmd'} ){
		$qsub_command = $cfg_hash->{'qsub_cmd'};
	}	
	#A string with parameters for the job taken from the program config file
	my $qsub_params = "";
	#Account to be used to launch job on the cluster
	if (defined $cfg_hash->{'qsub_account'} ){
		$qsub_params .= " -A ".$cfg_hash->{'qsub_account'};
	}


	#Active the restartable job
	if ( $cfg_hash->{'qsub_restartable'} eq 'YES'){
		$qsub_params .= " -r y";
	}
	######Given in input
	#Define the jobname
	if ($job_name ne 'null' ){
		$qsub_params .= " -N ".$job_name;
	}

	if (defined $cfg_hash->{'qsub_account'} ){
		$qsub_params .= " -W group_list=".$cfg_hash->{'qsub_account'}." ";
	}
	#Define the standard output
	if ($std_out ne 'null' ){
		$qsub_params .= " -o ".$std_out;
	}
	#Define the standard error
	if ($std_err ne 'null' ){
		$qsub_params .= " -e ".$std_err;
	}	
	#From the user configuration file (dependent by program)
	my $qsub_resources = "";
	#Define the resources needed
	my $suff ="";
	if ( $cfg_hash->{'qsub_param_for_all'} eq 'YES' ){
		$suff = $qsub_command;
	}else{
		$suff = $task;
	}
	
	#The queue to use
	if (defined $cfg_hash->{$suff.'_queue'}){
		$qsub_params .= " -q ".$cfg_hash->{$suff.'_queue'};
	}else{	
		if (defined $cfg_hash->{'qsub_queue'} ){
			$qsub_params .= " -q ".$cfg_hash->{'qsub_queue'};
		}
	}
		
	#QSUB RESOURCES
	if (defined $cfg_hash->{$suff.'_select'} or defined $cfg_hash->{$suff.'_cpus'} 
				or defined $cfg_hash->{$suff.'_nodes'} or defined $cfg_hash->{$suff.'_job_memory'}
				or defined $cfg_hash->{$suff.'_ncpus'}){
		$qsub_resources .= " -l ";
	}
	if (defined $cfg_hash->{$suff.'_select'}){
		$qsub_resources .= "select=".$cfg_hash->{$suff.'_select'}.":";
	}
	if (defined $cfg_hash->{$suff.'_ncpus'}){
		$qsub_resources .= "ncpus=".$cfg_hash->{$suff.'_ncpus'}.":";
	}
	if (defined $cfg_hash->{$suff.'_nodes'}){
		$qsub_resources .= "nodes=".$cfg_hash->{$suff.'_nodes'}.":";
	}
	if (defined $cfg_hash->{$suff.'_mem'}){
		$qsub_resources .= "mem=".$cfg_hash->{$suff.'_mem'}.":";
	}

	chop($qsub_resources);	
			
	#The walltime to use
	if (defined $cfg_hash->{$suff.'_walltime'} ){
		$qsub_params .= " -l walltime=".$cfg_hash->{$suff.'_walltime'};
	}
	#Use dependencies
	if ($dependencies ne 'no_depend'){
		$qsub_params .= " -W depend=".$cfg_hash->{$suff.'_depend'}.":$dependencies ";
	}		
	
	#Add environmental variables
	$qsub_params .= " $qsub_resources ";
	if ( $env_vars ne 'NONE'){
		 $qsub_params .= " -v $env_vars ";
	}
	my $command = $qsub_command.$qsub_params.$program_path;
	print_and_log( "Executing command: $command\n",$log_file);
	my $qsub_id = `$command`;
	chomp($qsub_id);
	if ($qsub_id ne ''){
		#print "Output from qsub: $qsub_id\n";;#DEBUGCODE
		#The qsub id should come out with number.host. Hence I remove the host
		$qsub_id =~ s/\..*//;
		#print "qsub_id: $qsub_id\n";#DEBUGCODE
	}else{
		print "ERROR: Unable to start $program_path. Please check error from qsub..\n ";
	}
	return $qsub_id;
	#die "Unable to execute: $command\n";
}



=head2 execute_threads

 Title  : execute_threads
 Usage  : execute_threads( - commands => 'the commands to be executed',
                              );

 Function: for each command it will launch a thread performing that operation
			
 Returns : nothing

=cut	 
sub execute_threads{
  my $commands = shift;#reference to array of commands
  
  my $processes = scalar(@$commands);
  print "Executing $processes commands \n";# I'll write a log in $log.\n"; 
  #Creates an object to fork processes in parallel
  my $pm = new Parallel::ForkManager($processes);

  #You can define a subroutine which is called when a child is terminated. It is called in the parent process.
  #  - pid of the process, which is terminated
  # - exit code of the program
  $pm->run_on_finish(
    sub {
      my($pid,$exit_code) = @_;
      #print "** Just got out of the pool with PID $pid and exit code: $exit_code\n";#DEBUCODE
    }
  );

	foreach my $cmd (@$commands){
    # Forks and returns the pid for the child:
    my $pid = $pm->start and next;
   
    # Here is the parallelized block
    # -----------
    #print "$pid ---> running \n";
    print  "Launching $cmd\n";
    try_exec_command($cmd);
    
    # Terminates the child process
    $pm->finish;
    print "Thread $pid finished. \n";
	}

  $pm->wait_all_children; 
}

=head2 execute_jobs

 Title  : execute_jobs
 Usage  : execute_jobs( - commands => 'the commands to be executed executing a job',
                              );

 Function: for each command it will launch a thread performing that operation
			If the command contains ' char (like R commands) it will not be executed 
			by qsub since the ENV contains the same char.
			Hence write the command into a file		
			
 Returns : nothing

=cut	 
sub execute_jobs{
  my $cfg_hash = shift;
  my $analysis_id = shift;
  my $commands = shift;#reference to array of commands
  my $program_call = shift;
  my $dependency = shift;
  my $fold_suff = shift;
  my $name = shift;
  my $log_file = shift;
  
  my $jobs_num = scalar(@$commands);
  print_and_log("Executing $jobs_num jobs \n",$log_file);# I'll write a log in $log.\n"; 

	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	my $job_ids = "";
	
	my $num_cmd = 1;					    
	foreach my $cmd (@$commands){
		#print_and_log("$analysis_indication This job will wait for these jobs: $dependency\n",$log_file);#DEBUGCODE
		my $job_name = "$name\_a$analysis_id\_$num_cmd";
		print_and_log("\nLaunching a job $job_name : ",$log_file);
		my $job_log = $cfg_hash->{$analysis_id.'_'.$fold_suff.'_log_fold'}."/$analysis_name\_$name\_$num_cmd.log";

		#If the command contains ' char it will not be executed by qsub since the ENV contains the same char.
		#Hence write the command into a file
		if ( $cmd =~ /\'/){
			open(CMD_F,">$job_log.tmpcmd");
			print CMD_F $cmd;
			close(CMD_F);
			
			#The command then is the path to the file
			$cmd = "$job_log.tmpcmd"; 
		}
		#Call the program using qsub and giving it parameters for the job and
		my $env_vars = "COMMAND='$cmd' ";
		my $job_id = try_exec_job( $cfg_hash,$env_vars,$fold_suff,$program_call,$job_name,$dependency,$job_log,$job_log,$log_file);	
		$num_cmd++;
		$job_ids .= $job_id.":";
		
		#Update the hash with dependencies
		my $deps_f = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'jc_dep_hash_f'};
		update_dependencies_hash($deps_f,$job_name,$job_id,$dependency,$log_file);
	}
	chop($job_ids);
	return $job_ids;
}


=head2 kill_job

 Title   : kill_job
 Usage   : kill_job( 
							- job_id -> the job id to be removed
							);
							
 Function: will kill bill
 
 Returns : nothing

=cut
sub kill_job{
	my $bill = shift;
	
	my $command = "qdel $bill";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	print "Job $bill killed!\n";#DEBUGCODE
}

###############################LOG FUNCTIONS#########################################
=head2 print_and_log

 Title   : print_and_log
 Usage   : print_and_log( - string -> the sentence that have to be print_and_log 
											- onlyLog -> a number);

 Function: will print (the string in input always in the log file and on the STDOUT
					if onlyLog is used then the print will be only in the log file
 Returns : nothing

=cut
sub print_and_log{
  my $string = shift;    
  my $logFile = shift; 
  my $onlyLog = shift;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  if ( defined $onlyLog){
    print LOG $string;
  }else{
    my $STDOUT = *STDOUT;
    my $LOG = *LOG;
    #Prints on both OUT
    for ($LOG, $STDOUT) { print $_ $string; }
  }
  #Close the log
	close(LOG)
}

=head2 log_and_exit

 Title   : log_and_exit
 Usage   : log_and_exit( - string -> the sentence that have to be print in the log

 Function: will print the string and die
 
 Returns : nothing

=cut
sub log_and_exit{
  my $string = shift; 
  my $logFile = shift;     
  
  my $STDOUT = *STDOUT;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  my $LOG = *LOG;
  #Print on both the log and the STDOUT the error
  for ($LOG, $STDOUT) { print $_ $string; }
  #if ( -d $logFolder){
    ##Moving the log file to the log folder before to exit the program
    #move($workingFolder.'/'.$logFile, $logFolder."/".$logFile) or die "Cannot move $workingFolder/$logFile to $logFolder/$logFile\n";
  #}
  #Close the log
	close(LOG);
  exit;
}


=head2 R_die

 Title : R_die
 Usage : R_die( );

 Function: die with R script and says to see the R log file. It moves the R log file in the log folder.

 Returns : nothing.. it dies..

=cut 
sub R_die{
  my $command = shift;
  my $R_script = shift;
  my $logFolder = shift;
  
  print "The program will die for problems with R. Command: $command\n";
 
  if ( -d $logFolder){
   # print "Trying to move ".$R_plot_script."out in ".$logFolder."/".$R_plot_script."out\n";
    my $RlogName = extract_name($R_script,0);
   #Moving the log file to the log folder before to exit 
    move($RlogName."out", $logFolder."/".$RlogName."out") 
      or die "Cannot move ".$R_script."out in ".$logFolder."/".$R_script."out";
    print "Please check the R log file at $logFolder/".$R_script."out\n";
  }else{
    print "The program tried without succeed to move the log file ".$R_script."out: $logFolder does not exists.\n";
  }
  die;
}  
####################################################################################################

=head2 JOB_get_info_from_jobname

 Title   : JOB_get_info_from_jobname
 Usage   : JOB_get_info_from_jobname( -database => 'name of the database,
                               );

 Function: Given the job name, this subroutine resolves if it is
					related with a sample, readfile or analysis and extracts
					its name.

 Returns : returns three variables: the task, the type of element, the name

=cut
sub JOB_get_info_from_jobname{
	my $cfg_hash = shift;
	my $jobname = shift;
	
	
	#Get the information about the job
	my @fields = split("_",$jobname);
	my $task = $fields[0];
	#Extract the type of element and its id
	$fields[1] =~ /([a-z])(\d+)/;
	
	#Get the db fields to search depending by the info
	my $table_f ="";
	my $name_f = "";
	my $id_f = "";
	my $type = "";
	if( $1 eq 'a'){
		$table_f = $cfg_hash->{'db_analyses_table'};
		$name_f = $cfg_hash->{'db_analysis_name'};
		$id_f = $cfg_hash->{'db_analysis_id'};
		$type = "analysis";
	}
	if ( $1 eq 's'){
		$table_f = $cfg_hash->{'db_sample_table'};
		$name_f = $cfg_hash->{'db_sample_name'};
		$id_f = $cfg_hash->{'db_sample_id'};
		$type = "sample";
	}
	if ( $1 eq 'r'){
		$table_f = $cfg_hash->{'db_readf_table'};
		$name_f = $cfg_hash->{'db_readf_name'};
		$id_f = $cfg_hash->{'db_readf_id'};
		$type = "readfile";
	}	
	my $name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																				$cfg_hash->{'db_pass'},$table_f,$name_f,$id_f,$2);
	
	return $task,$type,$name;
}

=head2 JOB_get_jobid_form_jobname

 Title : JOB_get_jobid_form_jobname
 Usage : JOB_get_jobid_form_jobname( );

 Function:  gets a job id using the job name 
 
 Returns : the job id

=cut 
sub JOB_get_jobid_form_jobname{
	my $jobname = shift;
	my $qsub_username  =shift;
	
	my $jobid = -1;
	my $id_field = 0;
	my $name_field = 3;
	
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
					
					if ( $fields[$name_field] eq $jobname){
						$jobid = $name_parts[0];
						print "jobid of  ".$name_parts[0].": $jobid \n" ;#DEBUGCODE
					}

					
				}
				if ($line =~ /---------------/){
					$start = 1;
				}
		}
		close(QSTAT_OUT);
	return 	$jobid;
}



=head2 parse_cand_genes_2_db

 Title   : parse_cand_genes_2_db
 Usage   : parse_cand_genes_2_db( -database => 'name of the database,
                               );

 Function: Parses a table like the following
					analysis_name	model	vmodel	gene	key	fmax
					UD_NA009_TRIO	new_ar	HTm,known,rare	RARS2	chr6_88299677_T_C	0.005


 Returns : nothing

=cut 
sub parse_cand_genes_2_db {
	my $cfg_hash = shift;
	my $file = shift;

	#Used to remember analyses already inserted in db
	my $analyses_inserted;
	
	#Temporary hash for the header
	my $head_hash;
	my $sep = "\t";
	print "Opening info file $file\n";
	open (FILE, "<$file") or die "Cannot open $file\n";
		
		#Go line by line and get the field values for each sample
		my $row_num = 0;
		while (my $line = <FILE>){
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
					my @fields = split($sep,$line);
					my $n_field = 0;
					my $ss_hash;	
					#Use an hash (ss_hash) to keep all the values for this iteration
					#put the value where the field-index says at the location of the sample id given
					foreach my $field (@fields){
							#print $head_hash->{$n_field}.":$field ";
							$ss_hash->{$head_hash->{$n_field}} = $field;#if sample id from db is used
							$n_field++;
					}

					###############################
					#INSERT ANALYSIS INFO					#
					###############################
					#
					my $analysisname = $ss_hash->{$cfg_hash->{'db_analysis_name'}};
					print "Considering $analysisname  and variant ".$ss_hash->{$cfg_hash->{'db_var_compid'}}."\n";
						#Get the analysisid
							my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},
																					$cfg_hash->{'db_analysis_name'},"'".$analysisname."'");

						#Get the analysisid
							my $varid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																					$cfg_hash->{'db_pass'},$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'},
																					$cfg_hash->{'db_var_compid'},"'".$ss_hash->{$cfg_hash->{'db_var_compid'}}."'");

							#Now try to insert the value
							if ( $analysisid > 0 and $varid > 0){
								my $fields = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'};
								my $values = "$analysisid,$varid";		
								#Model
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_model'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_model'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_model'}}."'";
								}	
								#Vmodel
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_vmodel'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_vmodel'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_vmodel'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_vmodel'}}."'";
								}									
																
								#causative
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_causative'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_causative'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_causative'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_causative'}}."'";
								}									
								#scriptver
								if ( $ss_hash->{$cfg_hash->{'db_progver'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_progver'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_progver'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_progver'}}."'";
								}	
								#validation
								if ( $ss_hash->{$cfg_hash->{'db_genes_name'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_genes_name'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_genes_name'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_genes_name'}}."'";
								}	
								#validation
								if ( $ss_hash->{$cfg_hash->{'db_freqmax'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_freqmax'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_freqmax'};
										$values .= ",".$ss_hash->{$cfg_hash->{'db_freqmax'}};
								}								
								#notes
								if ( $ss_hash->{$cfg_hash->{'db_cand_vars_notes'}} ne '-' and  $ss_hash->{$cfg_hash->{'db_cand_vars_notes'}} ne '' ){
										$fields .= ",".$cfg_hash->{'db_cand_vars_notes'};
										$values .= ",'".$ss_hash->{$cfg_hash->{'db_cand_vars_notes'}}."'";
								}																	

								#The second values are used to evaluate if the db_analysis_id is present
								my $fields2 = $cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_var_id'};
								my $values2 = "$analysisid,$varid";								
								#print "Inserting $fields = $values in table ".$cfg_hash->{'db_cand_vars_table'}."\n";#DEBUGCODE
								my ($new_group_id,$exist_group_id) = insert_only_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_cand_vars_table'},$cfg_hash->{'db_analysis_id'},
														$fields,$values,$fields2,$values2);														
							}else{
								print "Cannot find analysisid: $analysisname or varid: ".$ss_hash->{$cfg_hash->{'db_var_compid'}}.". Line will not be inserted\n";#DEBUGCODE
							}

				}#ELSE
		$row_num++;
	}
	close(FILE);

}
#######CHECK PARAMETERS OF THE CONFIGURATION FILE ################################

=head2 check_genedb_links

 Title : check_genedb_links
 Usage : check_genedb_links( );

 Function: this function checks the links to the files needed to build
 tables for gene annotation

 Returns : nothing

=cut 
sub check_genedb_links {
	my $cfg_hash = shift;
	my $log_file = shift;
	
	##PROGRAM CONFIG #THIS CHECK SHOULD BE PERFORMED ONLY THE FIRST TIME THAT THE PROGRAM IS RAN
	unless (head($cfg_hash->{'hpo_table_link'}) ) {
		unless (head($cfg_hash->{'hpo_table_link'}) ) {
				#rint_and_log ("Please check if [".$cfg_hash->{'hpo_table_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}
	unless (head($cfg_hash->{'gdi_table_link'}) ){
		unless (head($cfg_hash->{'gdi_table_link'}) ) {
				print_and_log ("Please check if [".$cfg_hash->{'gdi_table_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}
	unless (head($cfg_hash->{'rvis_table_link'}) ){
		unless (head($cfg_hash->{'rvis_table_link'}) ){
				print_and_log ("Please check if [".$cfg_hash->{'rvis_table_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}		
	unless (head($cfg_hash->{'omim_table_link'}) ){
		unless (head($cfg_hash->{'omim_table_link'}) ){
				print_and_log ("Please check if [".$cfg_hash->{'omim_table_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}
	unless (head($cfg_hash->{'refseq2genes_link'}) ){
		unless (head($cfg_hash->{'refseq2genes_link'}) ){
				print_and_log ("Please check if [".$cfg_hash->{'refseq2genes_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}		
}


=head2 check_datasets_links

 Title : check_datasets_links
 Usage : check_datasets_links( );

 Function: this function checks the links to the files needed to run the analyses
 
 Returns : nothing

=cut 
sub check_datasets_links {
	my $cfg_hash = shift;
	my $log_file = shift;
	
	##PROGRAM CONFIG #THIS CHECK SHOULD BE PERFORMED ONLY THE FIRST TIME THAT THE PROGRAM IS RAN
	unless (head($cfg_hash->{'genome_2bit_link'}) ) {
		unless (head($cfg_hash->{'genome_2bit_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'genome_2bit_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}
	unless (head($cfg_hash->{'bitToFasta_link'}) ){
		unless (head($cfg_hash->{'bitToFasta_link'}) ) {
				print_and_log ("Please check if [".$cfg_hash->{'bitToFasta_link'}."] is a correct URL.  will continue...\n",$log_file); 
		}
	}
	
	
	#unless (head($cfg_hash->{'known_sites_mills_link'}) ) {
		#unless (head($cfg_hash->{'known_sites_mills_link'}) ) {
				##print_and_log ("Please check if [".$cfg_hash->{'known_sites_mills_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}
	#unless (head($cfg_hash->{'known_db_snp_link'}) ){
		#unless (head($cfg_hash->{'known_db_snp_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_db_snp_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}
	#unless (head($cfg_hash->{'known_hapmap_link'}) ) {
		#unless (head($cfg_hash->{'known_hapmap_link'}) ) {
				##print_and_log ("Please check if [".$cfg_hash->{'known_hapmap_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}
	#unless (head($cfg_hash->{'known_1000g_link'}) ){
		#unless (head($cfg_hash->{'known_1000g_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_1000g_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}
	#unless (head($cfg_hash->{'known_omni_link'}) ) {
		#unless (head($cfg_hash->{'known_omni_link'}) ) {
				##print_and_log ("Please check if [".$cfg_hash->{'known_omni_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}
	#unless (head($cfg_hash->{'known_sites_mills_link'}) ){
		#unless (head($cfg_hash->{'known_sites_mills_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_sites_mills_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}			
	#unless (head($cfg_hash->{'known_sites_db_snp_link'}) ){
		#unless (head($cfg_hash->{'known_sites_db_snp_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_sites_db_snp_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}	
	#unless (head($cfg_hash->{'known_sites_hapmap_link'}) ){
		#unless (head($cfg_hash->{'known_sites_hapmap_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_sites_hapmap_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}			
	#unless (head($cfg_hash->{'known_sites_1000g_link'}) ){
		#unless (head($cfg_hash->{'known_sites_1000g_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_sites_1000g_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}		
	#unless (head($cfg_hash->{'known_sites_omni_link'}) ){
		#unless (head($cfg_hash->{'known_sites_omni_link'}) ) {
				#print_and_log ("Please check if [".$cfg_hash->{'known_sites_omni_link'}."] is a correct URL.  will continue...\n",$log_file); 
		#}
	#}				
	#if ( defined $cfg_hash->{'known_sites_1000gP3_link'}){
		#unless (head($cfg_hash->{'known_sites_1000gP3_link'}) ){
			#unless (head($cfg_hash->{'known_sites_1000gP3_link'}) ) {
					#print_and_log ("Please check if [".$cfg_hash->{'known_sites_1000gP3_link'}."] is a correct URL.  will continue...\n",$log_file); 
			#}
		#}		
	#}
	
}


=head2 check_config_parameters

 Title : check_config_parameters
 Usage : check_config_parameters( );

 Function: this function checks the parametersof the configuration file.
					It is dependent by the program we are running.

 Returns : nothing

=cut 
sub check_config_parameters {
	my $cfg_hash = shift;
	my $log_file = shift;
	
	#All the parameters from the configuration file are checked
	
	###PIPELINE STEPS EXECUTION WITH PRINTS USE FOR DEBUGCODE
	##QUALITY CHECK AND TRIMMING
	#checkVariable($cfg_hash->{'qc_exec'},'qc_exec',$cfg_hash->{'qc_prog'}." of raw reads will be performed.\n",$log_file);
	#checkVariable($cfg_hash->{'trimming'},'trimming',$cfg_hash->{'trim_prog'}." of raw reads will be performed.\n",$log_file);
	#checkVariable($cfg_hash->{'qc_after_trim'},'qc_after_trim',$cfg_hash->{'qc_prog'}." of trimmed reads will be performed.\n",$log_file);
	#
	##ALIGNMENT
	#checkVariable($cfg_hash->{'convert_scores'},'convert_scores',"Base scores will be converted in Phred using:".$cfg_hash->{'seqtk_seq_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'merge_pairs'},'merge_pairs',"Read pairs will be merged in a single fastq:".$cfg_hash->{'seqtk_mergepe_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'align_prog_exec'},'align_prog_exec',"Read will be aligned against ref genome with:".$cfg_hash->{'alignment_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'mark_rem_dup_exec'},'mark_rem_dup_exec',"Duplicate reads will be rearranged with:".$cfg_hash->{'mark_rem_dup_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'map_alt_exec'},'map_alt_exec',"ALT sequences will be found with:".$cfg_hash->{'map_alt_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'sam_sort_idx'},'sam_sort_idx',"The output from the alignment will be sorted and indexed with:".
	#								"samtools ".$cfg_hash->{'samview_prog'}.",".$cfg_hash->{'sort_prog'}." and ".$cfg_hash->{'index_prog'}.".\n",$log_file);
	#
	##REFINEMENT (BestPractices GATK)
	#checkVariable($cfg_hash->{'realign_target_creator'},'realign_target_creator',"GATK: Target regions for indel realignment will be found with:".$cfg_hash->{'realign_target_creator_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'indel_realigner'},'indel_realigner',"GATK: INDELs will be realigned with:".$cfg_hash->{'indel_realigner_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'base_recalibrator'},'base_recalibrator',"GATK: A file to use for base recalibration is created with:".$cfg_hash->{'base_recalibrator_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'print_reads'},'print_reads',"GATK: recalibration of bases will be performed with:".$cfg_hash->{'print_reads_prog'}.".\n",$log_file);
	#
	##MERGE READFILES
	#checkVariable($cfg_hash->{'mergesam'},'mergesam',"VarGenius will verify if is necessary to run:".$cfg_hash->{'mergebed_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'mrdup_groups'},'mrdup_groups',"If read files are merged , duplicates will be removed with:".$cfg_hash->{'mrdupbam_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'index_after_merge'},'index_after_merge',"After merge the BAM file is re-indexed with: samtools".$cfg_hash->{'index_prog'}.".\n",$log_file);
	#
	##VARIANT CALLING
	#checkVariable($cfg_hash->{'varcall'},'varcall',"GATK: Variants will be called with:".$cfg_hash->{'varcall_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'catvar'},'catvar',"GATK: If each sample has different VCFs. They will be merged with:".$cfg_hash->{'catvar_prog'}.".\n",$log_file);
	#
	##GENOTYPING
	#checkVariable($cfg_hash->{'genot'},'genot',"GATK: Genotypes will be inferred with:".$cfg_hash->{'genot_prog'}.".\n",$log_file);
	#
	##VARIANT FILTERING
	#checkVariable($cfg_hash->{'varfilt'},'varfilt',"GATK: Variants will be filtered with:".$cfg_hash->{'varfilt_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'varrecal'},'varrecal',"GATK: Variants will be recalibrated with:".$cfg_hash->{'varrecal_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'apprecal'},'apprecal',"GATK: Recalibration is applied with:".$cfg_hash->{'apprecal_prog'}.".\n",$log_file);
	#
	##PHASING
	#checkVariable($cfg_hash->{'phasing'},'phasing',"GATK: Haplotype will be inferred with:".$cfg_hash->{'phasing_prog'}.".\n",$log_file);
	#
	##STATISTICS
	#checkVariable($cfg_hash->{'mergebam'},'mergebam',"Statistics: BAM for reads files will be merged with".$cfg_hash->{'mergebed_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'coverbed'},'coverbed',"Statistics: Coverbed is deprecated!!".$cfg_hash->{'mergebed_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'noncovered_regions'},'noncovered_regions',"Statistics: Non covered regions will be found\n",$log_file);
	#checkVariable($cfg_hash->{'flagstat_sorted'},'flagstat_sorted',"Statistics: alignment statistics for sorted BAM with ".$cfg_hash->{'flagstat_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'flagstat_rmdup'},'flagstat_sorted',"Statistics: alignment statistics for BAM without duplicates with ".$cfg_hash->{'flagstat_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'flagstat_inters_targ'},'flagstat_inters_targ',"Statistics: alignment statistics for BAM covering the target with ".$cfg_hash->{'flagstat_prog'}.".\n",$log_file);
	#checkVariable($cfg_hash->{'coverage_analysis'},'coverage_analysis',"Statistics: Coverage of gene, intervals and samples will be performed with ".$cfg_hash->{'depth_coverage_prog'}.".\n",$log_file);
	#
	##OUTPUT MANAGEMENT
	#checkVariable($cfg_hash->{'create_vcf_to_ann'},'create_vcf_to_ann',"Output: VCF file will be prepared for annotation\n",$log_file);
	#checkVariable($cfg_hash->{'annov_db_dl'},'annov_db_dl',"Output: Databases for Annovar will be downloaded\n",$log_file);SB_ID_A551_TRY.hg19.final_out_resgr14.txt
	#checkVariable($cfg_hash->{'annot_all'},'annot_all',"Output: Variants will be annotated with Annovar\n",$log_file);
	#checkVariable($cfg_hash->{'write_output'},'write_output',"Output: The output file will be created using the VCF and the annotation\n",$log_file);
	#checkVariable($cfg_hash->{'rearrange_output'},'rearrange_output',"Output: The columns of the annotation will be rearranged\n",$log_file);
	#checkVariable($cfg_hash->{'get_coverage_plots'},'get_coverage_plots',"Output: Coverage plots will be generated\n",$log_file);
	#checkVariable($cfg_hash->{'update_overall_samples_plots'},'update_overall_samples_plots',"Output: Global coverage plots will be generated\n",$log_file);
	#checkVariable($cfg_hash->{'generate_html'},'generate_html',"Output: An HTML website will be created to show results\n",$log_file);
	#checkVariable($cfg_hash->{'gene_annotation'},'gene_annotation',"Output: A table for gene annotation will be created\n",$log_file);
	#checkVariable($cfg_hash->{'sep_joint_analyses'},'sep_joint_analyses',"Output: The output from a joint analysis will be separated per-sample\n",$log_file);
	#checkVariable($cfg_hash->{'supplementary_output'},'supplementary_output',"Output: A supplementary output will be created\n",$log_file);
	#checkVariable($cfg_hash->{'send_email'},'send_email',"Output: An email will be sent when the process finishes\n",$log_file);
	#
	##IMPORT VARIANTS AND ANNOTATIONS
	#checkVariable($cfg_hash->{'vcfimport'},'vcfimport',"Output: Variants will be imported into the database\n",$log_file);
	#
	##Remove all the temporary files
	#checkVariable($cfg_hash->{'remove_temp'},'remove_temp',"Temporary files will be removed\n",$log_file);
	#
	##Block database for exclusive use
	#checkVariable($cfg_hash->{'block_db'},'block_db',"Each analysis will use exclusively the database\n",$log_file);
	#
	#
	print_and_log("Checking configuration variables... ",$log_file);#DEBUGCODE
	##PIPELINE STEPS EXECUTION WITHOUT PRINT
	#QUALITY CHECK AND TRIMMING
	checkVariable($cfg_hash->{'qc_exec'},'qc_exec',"",$log_file);
	checkVariable($cfg_hash->{'trimming'},'trimming',"",$log_file);
	checkVariable($cfg_hash->{'qc_after_trim'},'qc_after_trim',"",$log_file);
	
	#ALIGNMENT
	checkVariable($cfg_hash->{'convert_scores'},'convert_scores',"",$log_file);
	checkVariable($cfg_hash->{'merge_pairs'},'merge_pairs',"",$log_file);
	checkVariable($cfg_hash->{'align_prog_exec'},'align_prog_exec',"",$log_file);
	checkVariable($cfg_hash->{'mark_rem_dup_exec'},'mark_rem_dup_exec',"",$log_file);
	checkVariable($cfg_hash->{'map_alt_exec'},'map_alt_exec',"",$log_file);
	checkVariable($cfg_hash->{'sam_sort_idx'},'sam_sort_idx',"",$log_file);
	
	#print_and_log("Checking REFINEMENT variables..\n",$log_file);#DEBUGCODE
	#REFINEMENT (BestPractices GATK)
	checkVariable($cfg_hash->{'realign_target_creator'},'realign_target_creator',"",$log_file);
	checkVariable($cfg_hash->{'indel_realigner'},'indel_realigner',"",$log_file);
	checkVariable($cfg_hash->{'base_recalibrator'},'base_recalibrator',"",$log_file);
	checkVariable($cfg_hash->{'print_reads'},'print_reads',"",$log_file);

	#MERGE READFILES
	checkVariable($cfg_hash->{'mergesam'},'mergesam',"",$log_file);
	checkVariable($cfg_hash->{'mrdup_groups'},'mrdup_groups',"",$log_file);
	checkVariable($cfg_hash->{'index_after_merge'},'index_after_merge',"",$log_file);

	#VARIANT CALLING
	checkVariable($cfg_hash->{'varcall'},'varcall',"",$log_file);
	checkVariable($cfg_hash->{'catvar'},'catvar',"",$log_file);

	#GENOTYPING
	checkVariable($cfg_hash->{'genot'},'genot',"",$log_file);
	
	#print_and_log("Checking VARIANT FILTERING variables..\n",$log_file);#DEBUGCODE
	#VARIANT FILTERING
	checkVariable($cfg_hash->{'varfilt'},'varfilt',"",$log_file);
	checkVariable($cfg_hash->{'varrecal'},'varrecal',"",$log_file);
	checkVariable($cfg_hash->{'apprecal'},'apprecal',"",$log_file);

	#PHASING
	checkVariable($cfg_hash->{'phasing'},'phasing',"",$log_file);

	#STATISTICS
	checkVariable($cfg_hash->{'mergebam'},'mergebam',"",$log_file);
	checkVariable($cfg_hash->{'coverbed'},'coverbed',"",$log_file);
	checkVariable($cfg_hash->{'noncovered_regions'},'noncovered_regions',"",$log_file);
	checkVariable($cfg_hash->{'flagstat_sorted'},'flagstat_sorted',"",$log_file);
	checkVariable($cfg_hash->{'flagstat_baserecal'},'flagstat_baserecal',"",$log_file);
	checkVariable($cfg_hash->{'flagstat_inters_targ'},'flagstat_inters_targ',"",$log_file);
	checkVariable($cfg_hash->{'coverage_analysis'},'coverage_analysis',"",$log_file);

	#OUTPUT MANAGEMENT
	checkVariable($cfg_hash->{'create_vcf_to_ann'},'create_vcf_to_ann',"",$log_file);
	checkVariable($cfg_hash->{'annov_db_dl'},'annov_db_dl',"",$log_file);
	checkVariable($cfg_hash->{'annot_all'},'annot_all',"",$log_file);
	checkVariable($cfg_hash->{'write_output'},'write_output',"",$log_file);
	checkVariable($cfg_hash->{'rearrange_output'},'rearrange_output',"",$log_file);
	checkVariable($cfg_hash->{'get_coverage_plots'},'get_coverage_plots',"",$log_file);
	checkVariable($cfg_hash->{'update_overall_samples_plots'},'update_overall_samples_plots',"",$log_file);
	checkVariable($cfg_hash->{'generate_html'},'generate_html',"",$log_file);
	checkVariable($cfg_hash->{'gene_annotation'},'gene_annotation',"",$log_file);
	checkVariable($cfg_hash->{'sep_joint_analyses'},'sep_joint_analyses',"",$log_file);
	checkVariable($cfg_hash->{'supplementary_output'},'supplementary_output',"",$log_file);
	checkVariable($cfg_hash->{'send_email'},'send_email',"",$log_file);

	#IMPORT VARIANTS AND ANNOTATIONS
	checkVariable($cfg_hash->{'vcfimport'},'vcfimport',"",$log_file);

	#Remove all the temporary files
	checkVariable($cfg_hash->{'remove_temp'},'remove_temp',"",$log_file);

	#Block database for exclusive use
	checkVariable($cfg_hash->{'block_db'},'block_db',"",$log_file);
	
	##*************Group 2*******************
	##PROGRAM PATHS
	
	#print_and_log("CheckingPATHS TO PROGRAMS ..\n",$log_file);#DEBUGCODE
  #CHECK THE PATHS TO PROGRAMS
  my @prog_paths = ($cfg_hash->{'fastqc_path'},$cfg_hash->{'java_path'},$cfg_hash->{'trimmomatic_path'},
										$cfg_hash->{'picard_path'},$cfg_hash->{'samtools_path'},$cfg_hash->{'gatk_path'},
										$cfg_hash->{'bwa_path'},$cfg_hash->{'bedtools_path'},
										$cfg_hash->{'annovar_path'},$cfg_hash->{'R_path'});
	print_and_log("Checking program paths... ",$log_file);#DEBUGCODE
	#	print "Checking fastqc_path: ".$cfg_hash->{'fastqc_path'}." java_path: ".$cfg_hash->{'java_path'}."..";
	#						print " trimmomatic_path: ".$cfg_hash->{'trimmomatic_path'}." picard_path: ".$cfg_hash->{'picard_path'};
		#					print " samtools_path: ".$cfg_hash->{'samtools_path'}." gatk_path: ".$cfg_hash->{'gatk_path'};
		#					print " bwa_path: ".$cfg_hash->{'bwa_path'}." bedtools_path: ".$cfg_hash->{'bedtools_path'};
		#					print " annovar_path: ".$cfg_hash->{'annovar_path'}." R_path: ".$cfg_hash->{'R_path'};
  foreach my $prog_path ( @prog_paths){

		if (!(check_presence($prog_path))){
			log_and_exit("$prog_path is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
		} 
	}
	
	#REFERENCE FILES
	#print_and_log("Checking REFERENCE FILES ..\n",$log_file);#DEBUGCODE
	#Reference Genome
	if (!(check_presence($cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'}))){
			log_and_exit($cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'}." is missing. Please check the links in the configuration file. Exiting...\n",$log_file);
	}

	if (!(-d $cfg_hash->{'target_reg_f'} )){
			log_and_exit($cfg_hash->{'target_reg_f'}." is not a correct folder. Please check into the configuration file. Exiting...\n",$log_file);
	}
  ##CHECK GATK REFERENCES
  ##print_and_log("Checking GATK REFERENCE FILES ..\n",$log_file);#DEBUGCODE
	#if (!(-d $cfg_hash->{'gatk_ref_f'} )){
			#log_and_exit($cfg_hash->{'gatk_ref_f'}." is not a correct folder. Please check into the configuration file. Exiting...\n",$log_file);
	#}
  #my @gatk_ref_paths = ($cfg_hash->{'known_mills'},$cfg_hash->{'known_db_snp'},$cfg_hash->{'known_hapmap'},
										#$cfg_hash->{'known_1000g'},$cfg_hash->{'known_omni'},$cfg_hash->{'known_sites_mills'},
										#$cfg_hash->{'known_sites_db_snp'},$cfg_hash->{'known_sites_hapmap'},$cfg_hash->{'known_sites_1000g'},
										#$cfg_hash->{'known_sites_omni'});	
  #foreach my $gatk_ref_path ( @gatk_ref_paths){
		#if ( defined $gatk_ref_path){
			#if (!(check_presence($cfg_hash->{'gatk_ref_f'}."/".$gatk_ref_path))){
				#log_and_exit("$gatk_ref_path is missing. Please check the into the configuration file or downlad it from GATK bundle. Exiting...\n",$log_file);
			#} 			
		#}

	#}
	
	#Gene panels files
	if ( defined $cfg_hash->{'disease_genes_lists_f'} ){
		my @disease_genes_lists_files = split(",",$cfg_hash->{'disease_genes_lists_f'});
		foreach my $disease_genes_lists_f ( @disease_genes_lists_files){
			if (!(-e $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'gene_panels_fold'}."/".$disease_genes_lists_f )){
					print_and_log("WARNING: cannot find $disease_genes_lists_f. This means you chose a disease list gene file that is not".
					" in ".$cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'gene_panels_fold'}.". Please check 'disease_genes_lists_f' into the configuration file and add any ".
					" new gene list to that folder. Exiting...\n",$log_file);
			}	
		}		
	}
	
		
	#print_and_log("Checking Folders ..\n",$log_file);#DEBUGCODE 
	
	if (!(-d $cfg_hash->{'scratch_f'} )){
			log_and_exit($cfg_hash->{'scratch_f'}." is not a correct folder. Please check into the configuration file. Exiting...\n",$log_file);
	}		
	if (!(-d $cfg_hash->{'storage_f'} )){
			log_and_exit($cfg_hash->{'storage_f'}." is not a correct folder. Please check into the configuration file. Exiting...\n",$log_file);
	}	
	if (!(-d $cfg_hash->{'annovar_db_f'} )){
			log_and_exit($cfg_hash->{'annovar_db_f'}." is not a correct folder. Please check into the configuration file. Exiting...\n",$log_file);
	}

	
	#If the html_host is not a directory, check the URL with head
	if (!-d ($cfg_hash->{'html_host'})){
					print_and_log ("WARNING: [".$cfg_hash->{'html_host'}."] is not a directory. Please be sure that is a correct URL. VarGenius will not check this...\n",$log_file); 
	}

		
	#Data for email sending
  if ( !($cfg_hash->{'email_author'} =~ /@/) or !($cfg_hash->{'email_author'} =~ /\.[a-zA-Z]+?/) ){
        log_and_exit("ERROR : ".$cfg_hash->{'email_author'}." is not a correct email address. Please check it in config file \n",$log_file);
  }
	my @recipients = split(",",$cfg_hash->{'email_recipients'});
  foreach my $recipient ( @recipients){
		if ( !($recipient =~ /@/) or !($recipient =~ /\.[a-zA-Z]+?/) ){
			log_and_exit("$recipient is not a correct email address. Exiting...\n",$log_file);
		} 
	}
	#*************Group 3*******************
	###JOBS RESOURCES REQUESTS	
	#print_and_log("Checking JOBS RESOURCES REQUEST  ..\n",$log_file);#DEBUGCODE 
	checkVariable($cfg_hash->{'qsub_param_for_all'},'qsub_param_for_all',"Resource parameters for QSUB will be picked from the qsub_ parameters for all programs\n",$log_file);
  if ( !($cfg_hash->{'qsub_cmd'} eq 'qsub') ){
        log_and_exit("ERROR : qsub_cmd must be 'qsub'. Please check it in config file \n");
  }
  
  if ( !($cfg_hash->{'qsub_walltime'} =~ /^\d+\:\d+:\d+$/) ){
        log_and_exit("ERROR : ".$cfg_hash->{'qsub_walltime'}." qsub_walltime must be set in this form hours:minutes:seconds. Please check it in config file \n",$log_file);
  }  
	checkVariable($cfg_hash->{'qsub_restartable'},'qsub_restartable',"",$log_file);
  
	my @tasks = ('qsub','qc','align','refine','MS','varcall','genot','varfilt','phasing','stats',
			'anstats','finalout','fastfout','jobclean');
	my @qsub_params = ("nodes","select","ncpus","threads");
  foreach my $task ( @tasks){
		if (defined $cfg_hash->{$task.'_ncpus'}){
			if ( !($cfg_hash->{$task.'_mem'} =~ /^\d+GB$/) ){
						log_and_exit("ERROR : ".$cfg_hash->{'qsub_mem'}." qsub_mem must be set in this form NumberGB (e.g. 60GB). Please check it in config file \n",$log_file);
			}  
		}
		#qqsub parameters that are positive integers
		foreach my $qsub_param (@qsub_params){
			if (defined $cfg_hash->{$task.'_'.$qsub_param }){
				if ( $cfg_hash->{$task.'_'.$qsub_param } ne ' ') {
					
					if ( !correct_type($cfg_hash->{$task.'_'.$qsub_param },"positiveint") ){
						log_and_exit("ERROR ".$cfg_hash->{$task.'_'.$qsub_param }.": $task\_$qsub_param must be a positive number. Please check it in config file \n",$log_file);	
					}	
				}
			}		
		}
	}
	
	#Other qsub parameters
	if ( defined $cfg_hash->{'align_threads_mem'} ){
		if ( !($cfg_hash->{'align_threads_mem'} =~ /^\d+GB$/) and !($cfg_hash->{'align_threads_mem'} =~ /^\d+G$/)  ){
					log_and_exit("ERROR : ".$cfg_hash->{'align_threads_mem'}." align_threads_mem must be set in this form NumberGB (e.g. 60GB). Please check it in config file \n",$log_file);
		}	
	}
	
	#Java maximum memory for GATK Programs
	my @java_progs = ("RTC","IR","BR","PR","MD","MSF","HAPC","CV","GGVCF","VR","VF","PBT","DOC","SV","CV","VA","CGP");
  foreach my $java_prog ( @java_progs){
		#Other qsub parameters
		
		if ( !($cfg_hash->{'javamem_'.$java_prog} =~ /^Xmx\d+[k,m,g]$/) ){
					log_and_exit("ERROR : ".$cfg_hash->{'javamem_'.$java_prog}." javamem_$java_prog must be set in this form Xmx[number][k,m,g]. Please check it in config file \n",$log_file);
		}			
	}	
	
	#*************Group 4*******************
	###SOFTWARE PARAMETERS
	#print_and_log("Checking SOFTWARE PARAMETERS  ..\n",$log_file);#DEBUGCODE 
	if (defined $cfg_hash->{'trimg_min_qual'}){
		if ( !correct_type($cfg_hash->{'trimg_min_qual'},"positiveint") ){
		 log_and_exit("ERROR : ".$cfg_hash->{'trimg_min_qual'}." must be a positive number. Please check it in config file \n",$log_file);	
		}
	}		
	if (defined $cfg_hash->{'clip_R1'}){
		if ( !correct_type($cfg_hash->{'clip_R1'},"positiveint") ){
		 log_and_exit("ERROR : clip_R1 must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	if (defined $cfg_hash->{'clip_R2'}){
		if ( !correct_type($cfg_hash->{'clip_R2'},"positiveint") ){
		 log_and_exit("ERROR : clip_R2 must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	#Samtools parameters
	if (defined $cfg_hash->{'samsort_threads'}){
		if ( !correct_type($cfg_hash->{'samsort_threads'},"positiveint") ){
		 log_and_exit("ERROR : samsort_threads must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	if ( defined $cfg_hash->{'samsort_mem'}){
		if ( !($cfg_hash->{'samsort_mem'} =~ /^\d+GB$/) and !($cfg_hash->{'samsort_mem'} =~ /^\d+G$/) ){
					log_and_exit("ERROR : samsort_mem must be set in this form NumberGB (e.g. 60GB). Please check it in config file \n",$log_file);
		}			
	}


	my @gatk_known = ("known_mills","known_db_snp","known_hapmap","known_1000g","known_omni","known_sites_mills",
								"known_sites_db_snp","known_sites_hapmap","known_sites_1000g","known_sites_omni");
	#RealignerTargetCreator
	my $var = $cfg_hash->{'known_RTC'}; 
  if (! grep {/\b$var\b/} @gatk_known ){
		 log_and_exit("ERROR : known_RTC must be one among [...]. Please check it in config file \n",$log_file);
	}
	#Checks if the value of the variable is appropriate
	checkVariable($cfg_hash->{'filter_mismatching_base_and_quals_RTC'},'filter_mismatching_base_and_quals_RTC',"",$log_file);
	checkVariable($cfg_hash->{'filter_mismatching_base_and_quals_BR'},'filter_mismatching_base_and_quals_BR',"",$log_file);
	checkVariable($cfg_hash->{'fix_misencoded_quality_scores'},'fix_misencoded_quality_scores',"",$log_file);
	checkVariable($cfg_hash->{'allow_potentially_misencoded_quality_scores'},'allow_potentially_misencoded_quality_scores',"",$log_file);
	if (defined $cfg_hash->{'nt_RTC'}){
		if ( !correct_type($cfg_hash->{'nt_RTC'},"positiveint") ){
		 log_and_exit("ERROR : nt_RTC must be a positive number. Please check it in config file \n",$log_file);	
		}
	}		

	#IndelRealigner parameters
	$var = $cfg_hash->{'known_IR'}; 
  if (! grep {/\b$var\b/} @gatk_known ){
		 log_and_exit("ERROR : known_IR must be one among [...]. Please check it in config file \n",$log_file);
	}
	checkVariable($cfg_hash->{'filter_bases_not_stored'},'filter_bases_not_stored',"",$log_file);


#BaseRecalibrator parameters
	foreach my $var_br (split(",",$cfg_hash->{'known_BR'})){
		if (! grep {/\b$var_br\b/} @gatk_known ){
			 log_and_exit("ERROR : known_BR must be one among [...]. Please check it in config file \n",$log_file);
		}		
	} 
	if ( defined $cfg_hash->{'covariates_BR'} ){
		if ( $cfg_hash->{'covariates_BR'}  !~ /^[A-za-z0-9\_\-\,]*$/i ){
			log_and_exit("ERROR : covariates_BR must be formed by only chars,numbers and '-', '_'. Please check it in config file \n",$log_file);
		}
	}
	if (defined $cfg_hash->{'nct_BR'}){
		if ( !correct_type($cfg_hash->{'nct_BR'},"positiveint") ){
		 log_and_exit("ERROR : nct_BR must be a positive number. Please check it in config file \n",$log_file);	
		}
	}
	checkVariable($cfg_hash->{'BQSR_uniq'},'BQSR_uniq',"",$log_file);

	#PrintReads parameters
	if (defined $cfg_hash->{'baq_PR'}){
		if ( $cfg_hash->{'baq_PR'}  !~ /^[A-za-z0-9\_\-]*$/i ){
			log_and_exit("ERROR : baq_PR must be formed by only chars,numbers and '-', '_'. Please check it in config file \n",$log_file);
		}
	}
	if (defined $cfg_hash->{'nct_PR'}){
		if ( !correct_type($cfg_hash->{'nct_PR'},"positiveint") ){
		 log_and_exit("ERROR : nct_PR must be a positive number. Please check it in config file \n",$log_file);	
		}
	}

	if (defined $cfg_hash->{'nct_HAPC'}){
		if ( !correct_type($cfg_hash->{'nct_HAPC'},"positiveint") ){
		 log_and_exit("ERROR : nct_HAPC must be a positive number. Please check it in config file \n",$log_file);	
		}
	}
	checkVariable($cfg_hash->{'do_genotypeGVCF_4_single_sample'},'do_genotypeGVCF_4_single_sample',"",$log_file);
	
	#Annovar annov_splicing_threshold
	if (defined $cfg_hash->{'nct_HAPC'}){
		if ( !correct_type($cfg_hash->{'annov_splicing_threshold'},"positiveint") ){
		 log_and_exit("ERROR : annov_splicing_threshold must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	

	#Non-covered regions parameters
	if (defined $cfg_hash->{'noncov_threshold'}){
		if ( !correct_type($cfg_hash->{'noncov_threshold'},"positiveint") ){
		 log_and_exit("ERROR : noncov_threshold must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	if (defined $cfg_hash->{'noncov_recip_overlap'}){
		if ( !correct_type($cfg_hash->{'noncov_recip_overlap'},"real") ){
		 log_and_exit("ERROR : noncov_recip_overlap must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	
	#Plots
	if (defined $cfg_hash->{'min_gene_coverage'}){
		if ( !correct_type($cfg_hash->{'min_gene_coverage'},"positiveint") ){
		 log_and_exit("ERROR : min_gene_coverage must be a positive number. Please check it in config file \n",$log_file);	
		}
	}	
	print_and_log("...DONE!\n",$log_file);#DEBUGCODE 
}

###############SUB ROUTINES FOR THE CHECK OF THE VARIABLES



#BEGIN { print (("  " x $main::x++) . "finish progmanagement compile\n") }

1;

