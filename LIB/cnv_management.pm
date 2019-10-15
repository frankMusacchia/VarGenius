#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2019>  <Francesco Musacchia>

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
    
    
package LIB::cnv_management;

## cnv_management.pm
#Permits the management of programs and pipelines for the detection and annotation of 
#CNV Copy-Number Variations

#BEGIN { print (("  " x $main::x++) . "Beginning cnvmanagement compile\n") }
BEGIN
{
	require Exporter;
	use vars qw(@ISA @EXPORT);
	@ISA = qw(Exporter);
	@EXPORT_OK = qw( run_exome_depth calculate_gc_content run_xhmm_pipeline
						ED_reallocate_output ED_output_preprocessing
						XHMM_output_preprocessing rearrange_cnv_output);
}
             
use strict;
use warnings;

#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked get_id_if_exists_from_db 
						do_query_select_all update_table get_id_if_exists_from_db_woconn );

#Using a library to manage files
use LIB::files_management qw( append_str_2_file extract_name delete_file append_hash_to_file
								insert_col_in_file_table delete_rows_containing append_file_2_file
								get_rows_containing count_lines_file insert_col_in_file
								list_to_array);
				
			
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(try_exec_command R_die print_and_log build_input_name_from_executed
								execute_jobs JOB_get_jobid_form_jobname alter_job try_exec_job);

#Using a library to manage gatk libraries
use LIB::gatk_management qw(run_GATK_GCContentByInterval run_GATK_DepthOfCoverage_gen);
				
										
############################################
##EXOME DEPTH
############################################
=head2 run_exome_depth

 Title   : run_exome_depth
 Usage   : run_exome_depth(  config_file => the config hash
								);

 Function: CNV pipeline
 
			This subroutine obtains BAM files from samples sequenced with the same 
			target kit and launches ExomeDepth
 
 Returns :
 
=cut
sub run_exome_depth {
	my $cfg_hash = shift;
	my $outFolder = shift;
	my $analysis_name = shift;
	my $analysis_id = shift;
	my $steps_array = shift;
	my $task = shift;	
	my $log_file = shift;
	
	#Use a maximum of 100 bam files randomly chosen
	my $max_bams = $cfg_hash->{'cnv_max_bam_to_use'};
	#The step at which the BAM should be collected: the one before the variant calling hence after Base Recalibration	
	my $step = $cfg_hash->{'varcall_step'};
	my $in_fold_suff = $cfg_hash->{'refine_step'};

	#Obtain the userid from the database given the group id
	my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
															
	#Get the target name
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    
	#Get the samples ids involved which are those obtained with the same kit
	#As from instructions of ExomeDepth the samples should be unrelated. Hence I am getting only probands
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." IN ".
	" (SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_targetbed'}." = '$targetbed' AND ".
	" ".$cfg_hash->{'db_analysis_infreq'}." = 1 AND ".$cfg_hash->{'db_analysis_userid'}." = $userid ".
	"AND ".$cfg_hash->{'db_sample_kinship'}." = 'P' ) ORDER BY RANDOM() LIMIT $max_bams;";
		
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $similar_anal_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#Add the current samples to the set of samples to consider
	print_and_log("Adding samples from the subject analysis..\n",$log_file);
	$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." = $analysis_id";
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $anal_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	foreach my $sample_id (keys %{$anal_samples}){
		#The content is useles while the index only will be used
		$similar_anal_sam->{$sample_id} = $sample_id;
	}
	my @curr_samples_edout = ();
	
	print_and_log("Getting BAM files for selected samples..\n",$log_file);
	#Create a file to write all the ExomeDepht that will be produced for which analysis
	#they are so that later they can be reallocated
	my $ed_out_paths_f = $outFolder."/". $cfg_hash->{'ED_cnvout_f'};
	#Get the black list of analysis names to be avoided
	#April 2019: I added this black list since ExomeDepth ends with weird errors while using some bam files
	#If the file does not exist the array will be empty
	my @ed_bam_blacklist_sex = list_to_array($cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ED_bamblacklist_f'}."_sex.txt","NO_NEW_LINE");
	my @ed_bam_blacklist_auto = list_to_array($cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ED_bamblacklist_f'}."_auto.txt","NO_NEW_LINE");
	#Just open the file to overwrite it
	open(EOP,">$ed_out_paths_f") or die "ERROR: Could not open $ed_out_paths_f\n";
	
					
	#The following array is needed to save all the bam file names
	#to give in input to the script
	my $bam_paths_f = $outFolder."/". $cfg_hash->{'ED_bamlist_f'};
	#Just open the file to overwrite it
	open(BP,">$bam_paths_f") or die "ERROR: Could not open $bam_paths_f\n";
	close(BP);
	foreach my $sample_id (keys %{$similar_anal_sam}){
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
		#When doing trials in trial folder use the commented variables
		#my $temp = "/pico/work/TELET_TIGEM/vargenius_analyses/";
		#my $inFolder = $temp."/".$curr_anal_name."/".$cfg_hash->{$in_fold_suff.'_out_f'};
		#my $cnvoutFolder = $temp."/".$analysis_name."/".$cfg_hash->{'fastfout_out_f'}."/".$cfg_hash->{'cnv_out_f'};		

		my $inFolder = $cfg_hash->{'work_fold'}."/".$curr_anal_name."/".$cfg_hash->{$in_fold_suff.'_out_f'};
		#Write the name of the ExomeDepth output here for further check (this is where ExomeDepth saves the outputs)
		my $cnvoutFolder = $cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{'fastfout_out_f'}."/".$cfg_hash->{'cnv_out_f'};		

		my $paramsIn;
		getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
		my $bam_path = $inFolder."/".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array).".".$cfg_hash->{'bam_ext'};
		
		#The sample name must not be in the blacklist
		if ( (not grep {/\b$curr_sample_name\b/} @ed_bam_blacklist_sex) and (not grep {/\b$curr_sample_name\b/} @ed_bam_blacklist_auto) ){
			#If exists and is non-zero
			if (-e $bam_path and !(-z $bam_path) ){
				print_and_log("Appending $bam_path..\n",$log_file);

				#Add the path to the list of BAM files
				append_str_2_file($bam_paths_f,$bam_path);	
				
				my $auto_o = $cnvoutFolder."/auto/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array)
				.".CNVs.tsv"; 
				my $sex_o = $cnvoutFolder."/sex/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array)
				.".CNVs.tsv";
				#If the sampleid is among those of the current analysis, put in the array	
				if ( grep {/\b$sample_id\b/} keys %{$anal_samples} ){
					push(@curr_samples_edout,$auto_o);
					push(@curr_samples_edout,$sex_o);	
				}
				
				#Record the output names for this file. The output later will be saved in each analysis folder
				#$similar_anal_sam->{$sample_id}	= "$curr_anal_name".$cfg_hash->{'mult_ann_sep'}."$auto_o,$sex_o";#REMOVE
				print EOP "$curr_sample_name,$curr_anal_name,$auto_o,$sex_o\n";
			}else{
				print_and_log("$bam_path cannot be used: non-existent or empty..\n",$log_file);			
			}			
		}else{
			print_and_log("$bam_path cannot be used since $curr_sample_name has been blacklisted..\n",$log_file);			

		}


	}
	close(EOP);
	
	#Check if the output from ExomeDepth is already present do not start it
	my $num_present = 0;
	foreach my $curr_sample_edout (@curr_samples_edout){
		if (-e $curr_sample_edout and ! (-z $curr_sample_edout)){
				$num_present = $num_present +1;
		}
	}
	my $exomedepthjobids = -1;
	#If at least one is missing for all samples and sex/auto, start ExomeDepth
	if ( $num_present != scalar(keys %{$anal_samples})){
		#ExomeDepth will be executed independently using Autosomes and Sex chromosomes
		my @commands_to_exec = ();
		
		#AUTO
		my $RLogPath = $outFolder."/log/$analysis_name\_exomedepth_auto".$cfg_hash->{'R_log_file'};
		my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_exome_depth_script'};
									
		my $args =  " $analysis_name $bam_paths_f  ".$cfg_hash->{'hum_ref'}." ".$cfg_hash->{'target_bed_auto'}." ".$cfg_hash->{'cnv_auto_f'}." $outFolder/".$cfg_hash->{'cnv_auto_f'}."/ "; #Analysis name.		
		#Execute ExomeDepth.R
		my $command = $cfg_hash->{'R_path2'}." CMD BATCH --no-save --no-restore '--args $args'"." $R_script $RLogPath";
		print_and_log("The command is: $command\n",$log_file);
		#Add the command to the list of commands
		push(@commands_to_exec,$command);				
		
		##SEX
		$RLogPath = $outFolder."/log/$analysis_name\_exomedepth_sex".$cfg_hash->{'R_log_file'};
		$args =  " $analysis_name $bam_paths_f  ".$cfg_hash->{'hum_ref'}." ".$cfg_hash->{'target_bed_sex'}." ".$cfg_hash->{'cnv_sex_f'}." $outFolder/".$cfg_hash->{'cnv_sex_f'}."/ "; #Analysis name.		
		#Execute ExomeDepth.R
		$command = $cfg_hash->{'R_path2'}." CMD BATCH --no-save --no-restore '--args $args'"." $R_script $RLogPath";
		print_and_log("The command is: $command\n",$log_file);
		#Add the command to the list of commands
		push(@commands_to_exec,$command);				
		
		#Execute jobs for DoC for Autosomes and Sex chromosomes
		#Path to the script to execute multiple jobs using an array of commands
		my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
		$exomedepthjobids = execute_jobs($cfg_hash,$analysis_id,\@commands_to_exec,$exec_job_script,"no_depend",$task,$task,$log_file);		
	}
	

	#I return the job ids so that the final output will wait them
	return $exomedepthjobids;
}



=head2 ED_output_preprocessing

 Title   : ED_output_preprocessing
 Usage   : ED_output_preprocessing(  config_file => the config hash
								);

 Function: Pre-process the output from ExomeDepht prior to merge with 
			other CNV output
			1. The output from ExomeDepth contains CNVs found in the single samples  
			hence for ExomeDepth we take the ouput for all samples for the analysis and concatenate them together.
			to use this file later we need to add a column specifing which sample is it.
 Returns :
 
=cut
sub ED_output_preprocessing {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $cnvoutFolder = shift;
	my $steps_array = shift;
	my $ed_out = shift;
	my $log_file = shift;
	
	#Operation 1:
	#In Exome Depth output, for each sample of the analysis find the output, add the column with the sample name and concatenate all the files
	#But if the user wants to have the complete output all the files will be used
	my $ed_out_temp = $ed_out.".temp";	
	#Clean ed_out_temp
	open(EOT,">$ed_out_temp");
	close(EOT);
	
	#Get sample ids for the analysis
	my $query = "";
	#If only the current analyis shold be used...
	if ( $cfg_hash->{'cnv_use_subject_only'} eq 'YES') {
		$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." = $analysis_id";
	}
	else{
		#Get the samples names reading all the available output of ExomeDepht from
		my $ed_out_paths_f = $cnvoutFolder."/". $cfg_hash->{'ED_cnvout_f'};
		my $sample_names = "";
		open (ED_OUTS,"<$ed_out_paths_f") or die "ERROR: Cannot open $ed_out_paths_f\n";
		while (my $row = <ED_OUTS>){
			chop($row);
			my @fields = split(",",$row);
			$sample_names .= "'".$fields[0]."',";

		}
		close(ED_OUTS);
		chop ($sample_names);
		$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_sample_name'}." IN ($sample_names);";
	}
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $anal_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	
	print_and_log( "Generating the complete ExomeDepth output in $ed_out using $ed_out_temp\n",$log_file);#DEBUGCODE

	my $auto_o = "";
	my $sex_o = "";
	my $headerok = 0;
	foreach my $sample_id (keys %{$anal_samples}){
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
						    
		#Get the ExomeDepth output
		#The step at which the BAM should be collected: the one before the variant calling hence after Base Recalibration	
		my $step = $cfg_hash->{'varcall_step'};
		my $paramsIn;
		getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		
		$auto_o = $cnvoutFolder."/auto/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array).".".$cfg_hash->{'bam_ext'}.".CNVs.tsv"; 
		$sex_o = $cnvoutFolder."/sex/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array).".".$cfg_hash->{'bam_ext'}.".CNVs.tsv";

		print_and_log( "Using $auto_o and $sex_o ..\n",$log_file);#DEBUGCODE
		
		#For both auto and sex do the following:
		#	
		#1. Add the column for the sample name (create an array to insert)
		#2. Remove the header 
		#3. Concatenate in a single ExomeDepth output		
		#Auto
		if ( -e $auto_o and !(-z $auto_o)){
			my @sample_name_col_auto = ($cfg_hash->{'db_sample_name'});
			my $num_lines = count_lines_file($auto_o);
			for (my $i = 1; $i < $num_lines; $i++){
				push(@sample_name_col_auto,$curr_sample_name);
			}			
			my $outfile_auto =  $auto_o.".temp1";
			insert_col_in_file($auto_o,\@sample_name_col_auto,0,"b",$outfile_auto);
			delete_rows_containing($outfile_auto,"start.p",$auto_o.".temp2");
			print_and_log( "appending $auto_o.temp2 to $ed_out_temp ..\n",$log_file);#DEBUGCODE			
			append_file_2_file($auto_o.".temp2",$ed_out_temp);			
			#At the first good output, put the header (I have to use .temp1 since it has the new header)
			if ( $headerok == 0){
				my $command = " head -n1 $auto_o.temp1 > $ed_out "; 
				(system ($command)) == 0 or die "ERROR: Cannot execute $command ";
				$headerok++;
			}
		}else{
			print_and_log( "WARNING: Could not find $auto_o ..\n",$log_file);#DEBUGCODE
		}
		#Sex
		if ( -e $sex_o and !(-z $sex_o)){
			my @sample_name_col_sex = ($cfg_hash->{'db_sample_name'});
			my $num_lines = count_lines_file($sex_o);
			for (my $i = 1; $i < $num_lines; $i++){
				push(@sample_name_col_sex,$curr_sample_name);
			}
			my $outfile_sex =  $sex_o.".temp1";
			insert_col_in_file($sex_o,\@sample_name_col_sex,0,"b",$outfile_sex);
			delete_rows_containing($outfile_sex,"start.p",$sex_o.".temp2");
			print_and_log( "appending $sex_o.temp2 to $ed_out_temp ..\n",$log_file);#DEBUGCODE
			append_file_2_file("$sex_o.temp2",$ed_out_temp);
		}else{
			print_and_log( "WARNING: Could not find $sex_o ..\n",$log_file);#DEBUGCODE
		}
	}

	#Concatenate the two files into a single ExomeDepth output
	append_file_2_file($ed_out_temp,$ed_out);
}	

=head2 rearrange_cnv_output

 Title   : rearrange_cnv_output
 Usage   : rearrange_cnv_output( -database => 'name of the database,
                               );

 Function:  
						
							
 Returns : 

=cut
sub rearrange_cnv_output {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $cnv_out = shift;
	my $cnv_rearranged_out = shift;
	my $cnv_geneann = shift;
	my $log_file = shift;

	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
	my $analysis_name = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
													$cfg_hash->{'db_analysis_id'},$analysis_id);
														
	#First add the omim annotation to the cnv_gene annotation file
	my $omim_info_hash;
	my $cnv_gene_filt_ann = "$cnv_geneann.ann";
	open (CNV_GA,"<$cnv_geneann") or die "ERROR cannot open $cnv_geneann";
	open (CNV_GAOUT,">$cnv_gene_filt_ann") or die "ERROR cannot open $cnv_gene_filt_ann";
	print CNV_GAOUT "compid\tgene\t".$cfg_hash->{'db_omim_ids'}."\t".$cfg_hash->{'db_omim_genedesc'}."\t".$cfg_hash->{'db_omim_phenotype'}."\n";
	while (my $row = <CNV_GA>){
		next if ($row =~ 'compid');
		chop($row);
		my @fields = split("\t",$row);
		my $genes = $fields[1];
		
		my @genes = split(",",$genes);
		my $db_omim_ids_str = "";
		my $db_omim_genedesc_str = "";
		my $db_omim_phenotype_str = "";
		foreach my $gene (@genes){
			my $db_omim_ids = "";
			my $db_omim_genedesc = "";
			my $db_omim_phenotype = "";
			if (! defined $omim_info_hash->{$gene}){
				$db_omim_ids = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_omim_ids'},
								$cfg_hash->{'db_genes_name'},"'".$gene."'");
				if ( length($db_omim_ids) <= 2 or $db_omim_ids eq ''){$db_omim_ids ="-";}
				$omim_info_hash->{$gene}->{"db_omim_ids"} = $db_omim_ids;
				$db_omim_genedesc = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_omim_genedesc'},
								$cfg_hash->{'db_genes_name'},"'".$gene."'");
				if ( length ($db_omim_genedesc) <= 2  or $db_omim_genedesc eq ''){$db_omim_genedesc = "-";}			
				$omim_info_hash->{$gene}->{"db_omim_genedesc"} = $db_omim_genedesc;
				$db_omim_phenotype = get_id_if_exists_from_db_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_omim_phenotype'},
								$cfg_hash->{'db_genes_name'},"'".$gene."'");
				if ( length($db_omim_phenotype) <= 2  or $db_omim_phenotype eq ''){$db_omim_phenotype = "-";}		
				$omim_info_hash->{$gene}->{"db_omim_phenotype"} = $db_omim_phenotype;						    				
			}
			$db_omim_ids_str .= $omim_info_hash->{$gene}->{"db_omim_ids"}.";";
			$db_omim_genedesc_str .= $omim_info_hash->{$gene}->{"db_omim_genedesc"}.";";
			$db_omim_phenotype_str .= $omim_info_hash->{$gene}->{"db_omim_phenotype"}.";";							    			    
		}
		chop($db_omim_ids_str);
		chop($db_omim_genedesc_str);
		chop($db_omim_phenotype_str);
		print CNV_GAOUT "$row\t$db_omim_ids_str\t$db_omim_genedesc_str\t$db_omim_phenotype_str\n";
	} 
	close (CNV_GA);
	close (CNV_GAOUT);
	
	#Disconnect db
	$dbh->disconnect();
  
	#Now use R merge to merge the CNV output and the OMIM annotation using the compid
	print_and_log(" Merging CNV output with OMIM annotation for $analysis_name..\n",$log_file);		
	my $cnvoutFolder  = $cfg_hash->{$analysis_id.'_cnv_out_f'};
	
	#compid is used as field to use for the merge. x: means that only the first dataset must include all rows
	my $RLogPath = $cnvoutFolder."/".$cfg_hash->{'log_fold'}."/$analysis_name\_merge_omimann".$cfg_hash->{'R_log_file'};
	my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_utils'};	
	my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args MERGE2 $cnv_out $cnv_gene_filt_ann compid x $cnv_rearranged_out' $R_script $RLogPath";
	
	print_and_log("The command is: $command\n",$log_file);
	try_exec_command($command) or R_die($command,$R_script);	
	
}
	
=head2 ED_reallocate_output

 Title   : ED_reallocate_output
 Usage   : ED_reallocate_output(  config_file => the config hash
								);

 Function: Reallocate output from ExomeDepth printed in a single folder into 
			each analaysis folder.
			Opens the file at ED_cnvout_f (program_config.txt) and reads
			the analysis name  and puts it into finalout_out/cnv folders
			
 Returns :
 
=cut
sub ED_reallocate_output {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	
	
	my $outFolder = $cfg_hash->{$analysis_id.'_cnv_out_f'};
	my $work_fold = $cfg_hash->{'work_fold'};#"/pico/work/TELET_TIGEM/vargenius_analyses/";
	
	#############
	##ExomeDepth output management
	##1. Put all the output from ExomeDepth for the different samples into the indicated analysis folder
	my $ed_out_paths_f = $outFolder."/". $cfg_hash->{'ED_cnvout_f'};
	#The file format is the following:
	#sampleid = analysisname]--[exomedepth.path.auto,exomedepth.path.sex
	open (ED_OUTS,"<$ed_out_paths_f") or die "ERROR: Cannot open $ed_out_paths_f\n";
	while (my $row = <ED_OUTS>){
		chop($row);
		my @fields = split(",",$row);
		my $curr_sample_name = $fields[0];
		my $curr_anal_name = $fields[1];
		my $auto_o = $fields[2];
		my $sex_o = $fields[3];
		
		#If exists, put the given ExomeDepht output into the finalout_out/cnv/sex|auto folder
		my $cnvoutFolder = $work_fold."/".$curr_anal_name."/".$cfg_hash->{'fastfout_out_f'}."/".$cfg_hash->{'cnv_out_f'};
		if ( -e $auto_o and ! (-z $auto_o)){
			move($auto_o,$cnvoutFolder."/".$cfg_hash->{'cnv_auto_f'}."/") or die "ERROR: Cannot move $auto_o to $cnvoutFolder/".$cfg_hash->{'cnv_auto_f'}."/";
		}
		if ( -e $sex_o and ! (-z $sex_o)){
			move($sex_o,$cnvoutFolder."/".$cfg_hash->{'cnv_sex_f'}."/") or die "ERROR: Cannot move $sex_o to $cnvoutFolder/".$cfg_hash->{'cnv_auto_f'}."/";
		}		
	}
	close (ED_OUTS);
}

############################################
##EXOME DEPTH
############################################







###########################################
#XHMM
############################################

=head2 XHMM_output_preprocessing

 Title   : XHMM_output_preprocessing
 Usage   : XHMM_output_preprocessing(  config_file => the config hash
								);

 Function: Pre-process the output from XHMM prior to merge with 
			other CNV output
		2. XHMM's contains CNVs for all samples analyzed, hence we need to grep results only for the samples of the subject analysis

 Returns :
 
=cut
sub XHMM_output_preprocessing {
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $cnvoutFolder = shift;
	my $xhmm_out = shift; 
	my $log_file = shift;
	#Get from the XHMM output only those lines concerning this analysis, hence starting with one of the subject samples names
	#

	#Get sample ids for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." = $analysis_id";
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $anal_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	
	my $auto_o = $cnvoutFolder."/auto/".$cfg_hash->{'xhmm_xcnv_out'};
	my $sex_o = $cnvoutFolder."/sex/".$cfg_hash->{'xhmm_xcnv_out'};
	#Put first the header
	(system ( " head -n1 $auto_o > $xhmm_out " )) == 0 or die "ERROR: Cannot execute head -n1 $auto_o > $xhmm_out ";	
	foreach my $sample_id (keys %{$anal_samples}){
		my $curr_sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    		
		get_rows_containing($auto_o,$curr_sample_name,$auto_o.".temp");
		get_rows_containing($sex_o,$curr_sample_name,$sex_o.".temp");
		print_and_log( "Appending $auto_o.temp and $sex_o.temp to $xhmm_out\n",$log_file);#DEBUGCODE		
		append_file_2_file($auto_o.".temp",$xhmm_out);
		append_file_2_file($sex_o.".temp",$xhmm_out);
	}
}

=head2 calculate_gc_content

 Title   : calculate_gc_content
 Usage   : calculate_gc_content(  config_file => the config hash
								);

 Function: CNV pipeline
 
			Using GATK returns a file containing regions with high GC content
			It is useful to perform upfront filtering of exome targets expected to pose more challenges
			to detecting CNV due to their more extreme properties (e.g., high and low GC content
			targets tend to be captured less well on exome hybridization arrays).
 
 Returns :
 
=cut
sub calculate_gc_content {
	my $cfg_hash = shift;
	my $targetbed = shift;
	my $log_file = shift;
	
	my $gc_content_locus = $cfg_hash->{'target_reg_f'}."/".$cfg_hash->{'xhmm_gc_content_locus_f'};
	my $target_bed_path = $cfg_hash->{'target_reg_f'}."/".$targetbed;
	my $extreme_gc_f = $cfg_hash->{'target_reg_f'}."/".extract_name($targetbed,"noext")."_".$cfg_hash->{'xhmm_extreme_gc_f'};

	print_and_log("Checking if $extreme_gc_f exists and is non-empty\n",$log_file);
	if ( !(-e $extreme_gc_f) or -z $extreme_gc_f ){
		run_GATK_GCContentByInterval($cfg_hash,$target_bed_path,$gc_content_locus,$log_file);

		my $command = "cat  $gc_content_locus | awk '".'{if ($2 < 0.1 || $2 > 0.9) print $1}'."' > $extreme_gc_f ";
		print_and_log("Executing $command\n",$log_file);
		try_exec_command($command) > 0 
		 or die "ERROR [$?]: Cannot generate extreme gc content file $extreme_gc_f (Command: $command)\n";		
	}else{
		print_and_log("$extreme_gc_f exists and will not be created...\n",$log_file);
	}

	#Update table
	my $fields = $cfg_hash->{'db_target_gccontentf'};
	my $values = "1";
	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
	$cfg_hash->{'db_pass'},$cfg_hash->{'db_target_files_table'},$fields,$values,$cfg_hash->{'db_targetname'},"'".$targetbed."'");
    
}
	
	

=head2 run_xhmm_pipeline

 Title   : run_xhmm_pipeline
 Usage   : run_xhmm_pipeline(  config_file => the config hash
								);

 Function: XHMM CNV pipeline
 	
 
 Returns :
 
=cut
sub run_xhmm_pipeline {
	
	#Launch a job for each BAM file executing DepthOfCoverage
	my $cfg_hash = shift;
	my $outFolder = shift;
	my $analysis_name = shift;
	my $analysis_id = shift;
	my $steps_array = shift;
	my $task = shift;
	my $log_file = shift;
	
	my $max_bams = $cfg_hash->{'cnv_max_bam_to_use'};
	#The step at which the BAM should be collected: the one before the variant calling hence after Base Recalibration	
	my $step = $cfg_hash->{'varcall_step'};
	my $in_fold_suff = $cfg_hash->{'refine_step'};
	
	my @chr_types = ("auto","sex");
	
	#Obtain the userid from the database given the analysis id
	my $userid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
														
	#Get the target name
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    
	#Get the samples ids involved which are those obtained with the same kit
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." IN ".
	" (SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_targetbed'}." = '$targetbed' AND ".
	" ".$cfg_hash->{'db_analysis_infreq'}." = 1 AND ".$cfg_hash->{'db_sample_kinship'}." = 'P' AND ".
	$cfg_hash->{'db_analysis_userid'}." = $userid) ORDER BY RANDOM() LIMIT $max_bams;";	
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $similar_anal_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#Add the subject samples to the set of samples to consider
	$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." = $analysis_id";
	print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	my $anal_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	foreach my $sample_id (keys %{$anal_samples}){
		#The content is useles while the index only will be used
		$similar_anal_sam->{$sample_id} = $sample_id;
	}
	
	
	print_and_log("Getting BAM files for selected samples..\n",$log_file);
	
	#Parameters for DepthOfCoverage
	my $param_str = $cfg_hash->{'xhmm_DOC_params'};
		
	#The following array is needed to save all the bam file names
	#to give in input to the script
	my @commands_to_exec = ();
	foreach my $sample_id (keys %{$similar_anal_sam}){
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
		#When doing trials in trial folder use the commented variables
		#my $temp = "/pico/work/TELET_TIGEM/vargenius_analyses/";
		#my $inFolder = $temp."/".$curr_anal_name."/".$cfg_hash->{$in_fold_suff.'_out_f'};
		#my $statsFolder = $temp."/".$curr_anal_name."/".$cfg_hash->{'stats_out_f'};
		#
		my $inFolder = $cfg_hash->{'work_fold'}."/".$curr_anal_name."/".$cfg_hash->{$in_fold_suff.'_out_f'};
		my $statsFolder = $cfg_hash->{'work_fold'}."/".$curr_anal_name."/".$cfg_hash->{'stats_out_f'};
		
		
		my $paramsIn;
		getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
		my $bam_path = $inFolder."/".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,$steps_array).".".$cfg_hash->{'bam_ext'};
		
		my $get_doc_command = "";
		#If exists and is non-zero
		if ( -e $bam_path and !(-z $bam_path) ){
			print_and_log("Appending $bam_path..\n",$log_file);
			#AUTO
			my $out_suffix = "$statsFolder/$curr_sample_name.XHMM.AUTO.DOC";
			#Check first if it exists already
			if ( !(-e $out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'}) or -z $out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'} ){
				$get_doc_command = run_GATK_DepthOfCoverage_gen($cfg_hash,$bam_path,$param_str,$cfg_hash->{'target_bed_auto'},$out_suffix,$curr_anal_id,$log_file,"GET");
				#print_and_log("Putting $get_doc_command into the list..\n",$log_file);
				#Add the command to the list of commands
				push(@commands_to_exec,$get_doc_command);
			}else{
				print_and_log($out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'}." already present!..\n",$log_file);			
			}
						
			#SEX
			$out_suffix = "$statsFolder/$curr_sample_name.XHMM.SEX.DOC";
			#Check first if it exists already
			if ( !(-e $out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'}) or -z $out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'} ){
				$get_doc_command = run_GATK_DepthOfCoverage_gen($cfg_hash,$bam_path,$param_str,$cfg_hash->{'target_bed_sex'},$out_suffix,$curr_anal_id,$log_file,"GET");
				print_and_log("Putting $get_doc_command into the list..\n",$log_file);#DEBUGCODE		
				#Add the command to the list of commands
				push(@commands_to_exec,$get_doc_command);
			}else{
				print_and_log($out_suffix.".".$cfg_hash->{'DOC_sam_interval_sum'}." already present!..\n",$log_file);			
			}
		}else{
			#print_and_log("$bam_path cannot be used: non-existent or empty..\n",$log_file);			
		}
	}
	
	print_and_log("Execute jobs for DoC for Autosomes and Sex chromosomes..\n",$log_file);
	#Execute jobs for DoC for Autosomes and Sex chromosomes
	#Path to the script to execute multiple jobs using an array of commands
	my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};
	my $doc_jobs = execute_jobs($cfg_hash,$analysis_id,\@commands_to_exec,$exec_job_script,"no_depend",$task,"XDOC",$log_file);
	
	
	###	
	##Now launch a job that waits for doc_jobs and executes the steps of the pipeline for XHMM
	#each command of XHMM is obtained and printed into an .sh file that is finally launched
	#as a single job
	
	
	#Get a string with al the samples required:
	my $samples_str = "";
	foreach my $sample_id (keys %{$similar_anal_sam}){
		$samples_str .= "$sample_id-";
	}	
	chop($samples_str);
	
	my $out_f_with_paths = "$outFolder/$analysis_name\_XHMM_BAMS_2_merge.txt";
	my $params = "'samples=$samples_str;out_f=$out_f_with_paths'";
	my $func = "CHECK_CNV_INPUT";
	my $env_vars = "CONFIG_FILE=".$cfg_hash->{'cfg_hash_f'}.",FUNC=$func,PARAMS=$params,LOG_FILE=$log_file ";
	$exec_job_script = $cfg_hash->{'lib_folder'}."/exec_job_PERL_func.pl";
	my $job_name = "CNV_CHECK";
	my $job_log = "$outFolder/log/$analysis_name\_XHMM_$func.log";
	my $job_err = "$outFolder/log/$analysis_name\_XHMM_$func.err";
	my $check_bam_job_id = try_exec_job( $cfg_hash,$env_vars,"cnvo",$exec_job_script,$job_name,$doc_jobs,$job_log,$job_err,$log_file);
		
	my $xhmm_job_2_wait = "";
	#For each chromosome type (autosomes and sex)	
	foreach my $chr_type (@chr_types){
		print_and_log("Getting the commands for the pipeline for $chr_type chromosome..\n",$log_file);

		my @commands_4_job = ();
		#Merge
		my $parameters = "";
		#Get output from DoC for all samples
		#if ($chr_type eq 'AUTO'){
			#$parameters .= " $samp_int_sums_paths_auto ";
		#}else{
			#$parameters .= " $samp_int_sums_paths_sex ";
		#}
		#
		my $var_command = "";
		if ($chr_type eq 'auto'){
			$var_command = 'samp_int_sums_paths_auto=$(sed '."'"."1q;d"."'"." $out_f_with_paths)";
			$parameters .= ' $samp_int_sums_paths_auto ';
		}else{
			$var_command = 'samp_int_sums_paths_sex=$(sed '."'"."2q;d"."'"." $out_f_with_paths)";
			$parameters .= ' $samp_int_sums_paths_sex ';
		}
		push(@commands_4_job,$var_command."\n");
		
		#MERGE COMMAND
		my $merge_out = $outFolder."/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_merge_out"}.".txt";
		my $xhmm_command = run_xhmm($cfg_hash,"NONE","--mergeGATKdepths",$parameters,$merge_out,"GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");


		#STEP1: Filter samples and targets and prepare for normalization
		my $RD_filtered_targets = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_RD_filtered_targets"}.".txt";
		my $RD_filtered_samples = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_RD_filtered_samples"}.".txt";	
		
		$parameters = "--outputExcludedTargets $RD_filtered_targets ".
						"--outputExcludedSamples $RD_filtered_samples ";
						
		#Check whether exclude targets file is present, add it
		my $target_bed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
		my $exclude_gc_f = $cfg_hash->{'target_reg_f'}."/".extract_name($targetbed,"noext")."_".$cfg_hash->{'xhmm_extreme_gc_f'};
		if ( -e $exclude_gc_f ){
			$parameters .= "--excludeTargets $exclude_gc_f ";
		}
		#Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step1'}){
			$parameters .= $cfg_hash->{'xhmm_params_step1'};
		}							
		my $filtered_centered_RD = $outFolder."/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_filtered_centered_RD"}.".txt";
		$xhmm_command = run_xhmm($cfg_hash,$merge_out,"--matrix",$parameters,$filtered_centered_RD,"GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");
		
		
		#STEP 2 Run principal component analysis (PCA) on mean-centered data
		my $pca_out = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_pca_out"};
		$parameters = " --PCAfiles $pca_out ";
		#Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step2'}){
			$parameters .= $cfg_hash->{'xhmm_params_step2'};
		}
		$xhmm_command = run_xhmm($cfg_hash,$filtered_centered_RD,"--PCA",$parameters,"NONE","GET",$log_file);	
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");	
		
		#STEP 3 Normalize mean-centered data using PCA information
		my $PCA_normalized = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_PCA_normalized"}.".txt";
		$parameters = " --PCAfiles $pca_out ".
		" --normalizeOutput $PCA_normalized ";
		#Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step3'}){
			$parameters .= $cfg_hash->{'xhmm_params_step3'};
		}	
		$xhmm_command = run_xhmm($cfg_hash,$filtered_centered_RD,"--normalize",$parameters,"NONE","GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");
		
		#STEP 4 Filter and calculate z-scores for the data
		my $filtered_sample_zscores_RD = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_filtered_sample_zscores_RD"}.".txt ";
		my $zscores_RD_filtered_targets = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_zscores_RD_filtered_targets"}.".txt";
		my $zscores_RD_filtered_samples = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_zscores_RD_filtered_samples"}.".txt";
		$parameters = " --outputExcludedTargets $zscores_RD_filtered_targets ".
		" --outputExcludedSamples $zscores_RD_filtered_samples ";
		#Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step4'}){
			$parameters .= $cfg_hash->{'xhmm_params_step4'};
		}
		$xhmm_command = run_xhmm($cfg_hash,$PCA_normalized,"--matrix",$parameters,$filtered_sample_zscores_RD,"GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");	
		
		
		#STEP 5 Filter original read-depth data to restrict to same samples and targets as filtered, normalized data
		my $same_filtered_RD = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_same_filtered_RD"}.".txt ";;
		$parameters = " --excludeTargets $RD_filtered_targets ".
		" --excludeTargets $zscores_RD_filtered_targets ".
		" --excludeSamples $RD_filtered_samples ".
		" --excludeSamples $zscores_RD_filtered_samples ";
		#Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step5'}){
			$parameters .= $cfg_hash->{'xhmm_params_step5'};
		}	
		$xhmm_command = run_xhmm($cfg_hash,$merge_out,"--matrix",$parameters,$same_filtered_RD,"GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");	
		

		#STEP 6 Call CNVs in normalized data 
		$parameters = " -p ".$cfg_hash->{'main_data_folder'}."/".$cfg_hash->{"xhmm_params_f"}.
		 " -R $same_filtered_RD ".
		 " -c $outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_xcnv_out"}.
		 " -a $outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_xcnvaux_out"}.
		 " -s $outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_out"};
		 #Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step6'}){
			$parameters .= $cfg_hash->{'xhmm_params_step6'};
		}
		$xhmm_command = run_xhmm($cfg_hash,$filtered_sample_zscores_RD,"--discover",$parameters,"NONE","GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");	
		
		
		#STEP 7 Genotype
		$parameters = " -p ".$cfg_hash->{'main_data_folder'}."/".$cfg_hash->{"xhmm_params_f"}.
		 " -R $same_filtered_RD ".
		 " -g $outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_xcnv_out"}.
		 " -F ".$cfg_hash->{'hum_ref'}." ".#Fasta reference
		 " -v $outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_vcf_out"};
		 #Add parameters from config file (if exist)
		if ( defined $cfg_hash->{'xhmm_params_step7'}){
			$parameters .= $cfg_hash->{'xhmm_params_step7'};
		}
		$xhmm_command = run_xhmm($cfg_hash,$filtered_sample_zscores_RD,"--genotype",$parameters,"NONE","GET",$log_file);
		#print_and_log("Adding: $xhmm_command..\n",$log_file);#DEBUGCODE
		push(@commands_4_job,$xhmm_command."\n");
		
		
		##RUN XHMM JOB
		#Now that all the commands have been registered as strings, write a BASH script to execute as a job
		my $bash_script = "$outFolder/".$cfg_hash->{"cnv_$chr_type\_f"}."/".$cfg_hash->{"xhmm_script"};
		print_and_log("Saving commands for XHMM pipeline into $bash_script..\n",$log_file);

		open(XHMM_S,">$bash_script") or die "ERROR: Cannot open $bash_script. Check permissions...\n";
		print XHMM_S @commands_4_job;
		close(XHMM_S);
		#Execute the pipeline XHMM into a job
		my $job_name = "XHMM_a$analysis_id\_$chr_type";
		my $job_log = "$outFolder/log/$analysis_name\_XHMM".$cfg_hash->{"cnv_$chr_type\_f"}.".log";
		my $job_err = "$outFolder/log/$analysis_name\_XHMM".$cfg_hash->{"cnv_$chr_type\_f"}.".err";
		print_and_log("Executing job $job_name for XHMM pipeline ..\n",$log_file);
		$xhmm_job_2_wait .= try_exec_job( $cfg_hash,"NONE","XDOC",$bash_script,$job_name,$doc_jobs.":$check_bam_job_id",$job_log,$job_err,$log_file).":";						
		
	}
	#remove last :
	chop($xhmm_job_2_wait);
	return $xhmm_job_2_wait;

}	


=head2 run_xhmm

 Title   : run_xhmm
 Usage   : run_xhmm(  config_file => the config hash
								);

 Function: Executes XHMM. The user can choose 
			function, input and output
 	
 
 Returns :
 
=cut
sub run_xhmm {
	my $cfg_hash = shift;
	my $input = shift;
	my $function = shift;
	my $parameters = shift;
	my $output = shift;
	my $get_string = shift;
	my $log_file = shift;
	
	#Set the XHMM call
	my $xhmm_call = $cfg_hash->{'xhmm_path'};
	
	if ( defined $function ) {
		$xhmm_call .= " $function ";
	}else{
		die "ERROR: cannot execute XHMM without the function\n";
	}
	
	my $params = "";

	#All the times that the -r parameter is used
	if ( $input ne 'NONE'){
		$params .= " -r $input ";
	}
		
	if ( $output ne 'NONE'){
		$params .= " -o $output ";
	}
	
	if ( $parameters ne 'NONE'){
		$params .= " $parameters ";	
	}else{
		die "ERROR: No params in input. All functions of XHMM need parameters. There's something strange here...\n";
	}
	
	my $command = "$xhmm_call $params ";
	#If you ask to get the command for further execution it will not launch here the command			 
	if ($get_string ne 'GET' ){
		print_and_log("Executing command: $command\n",$log_file);
		try_exec_command($command) or die "Unable to execute command: $command\n";
	}
	
	return $command;
}
###########################################
#XHMM
############################################
1;
