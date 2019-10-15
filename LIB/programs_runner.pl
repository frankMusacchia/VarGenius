#!/usr/bin/perl

#vargenius.pl path
use lib '/pico/work/TELET_TIGEM/VarGeniusBeta/';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

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
    


#To use it you need:
#Libraries
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use IO::Handle;#To immediately print with autoflush 
use File::Copy;#To manage files
 



#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash correct_type try_exec_command
				print_and_log log_and_exit  execute_threads build_input_name_from_executed
				sort_samples get_flowcell_and_lane 
				get_ped_file_for_enlarged_analysis get_number_of_variants separate_input_ids
				good_reads_quality  get_jobid_form_jobname alter_job execute_jobs
				parse_cand_genes_2_db load_table_into_db get_output_name_from_executed
				fields_2_db try_exec_job run_FREEBAYES_vc get_samples_id_from_analysis
				compress_data_to_store update_dependencies_hash);


#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present
				download_file dl_and_extract delete_file extract_columns_from_file
				append_str_2_file check_presence file_num_rows 
				invert_cols_position get_col_index insert_col_in_file_table delete_columns
				extract_col_from_file file_name_wostrange_chars append_hash_to_file
				delete_rows_containing file_list_to_array extract_colnum_from_file_linux
				shift_column compress_folder extract_name append_str_2_file_if_path_notexist
				join_files_with_cat separate_elements_on_col join_vcf_files_with_cat
				vcf_insert_string_into_field vcf_filter_identical_calls
				overwrite_str_in_file_if_exists tsv_2_xls list_to_array);
				
#Using a library for database management
use LIB::db_management qw(getSampleConfiguration_locked update_analysis_status_locked
			select_distinct_samples get_id_if_exists_from_db insert_into_table_locked
			parse_HPO_genes_to_phenotypes parse_GDI_scores parse_RVIS_scores parse_RefSeq_2_Genes
			parse_genes2RefSeq parse_OMIM_2_Genes do_query_select_all update_table get_count_from_selected
			fetch_all_rows update_table_woconn createHPOTable);
			
#Using a library for gatk functions management
use LIB::gatk_management qw( run_GATK_RealignerTargetCreator run_GATK_IndelRealigner 
            run_GATK_BaseRecalibrator_perlane run_GATK_HaplotypeCaller
            run_GATK_CatVariants run_GATK_GenotypeGVCFs run_GATK_VariantFiltration
            run_GATK_VariantRecalibration run_GATK_ApplyRecalibration run_GATK_PhaseByTransmission
            run_GATK_DepthOfCoverage run_GATK_BaseRecalibrator_persample run_GATK_PrintReads_persample
            run_GATK_SelectVariants run_GATK_CombineVariants run_GATK_VariantAnnotator 
            run_GATK_CalculateGenotypePosteriors run_GATK_SplitNCigarReads
            run_GATK_GenotypeGVCFs_parallel);

#Using a library for gatk4 functions management
use LIB::gatk4_management qw(run_GATK4_RealignerTargetCreator run_GATK4_IndelRealigner 
            run_GATK4_BaseRecalibrator_perlane run_GATK4_HaplotypeCaller
            run_GATK4_CatVariants run_GATK4_GenotypeGVCFs run_GATK4_VariantFiltration
            run_GATK4_VariantRecalibration run_GATK4_ApplyRecalibration 
            run_GATK4_PhaseByTransmission run_GATK4_DepthOfCoverage
            run_GATK4_BaseRecalibrator_persample run_GATK4_PrintReads_persample
            run_GATK4_SelectVariants run_GATK4_CombineVariants
            run_GATK4_VariantAnnotator run_GATK4_CalculateGenotypePosteriors 
            run_GATK4_SplitNCigarReads run_GATK4_GenotypeGVCFs_parallel
            run_GATK4_ApplyBQSR run_GATK_MergeVcfs);
            
#Using a library for picard functions management
use LIB::picard_management qw( run_PICARD_AddOrReplaceReadGroups run_PICARD_MarkDuplicates_gen
				run_PICARD_MarkDuplicates_chrom
				run_PICARD_MergeSamFiles_gen run_PICARD_MergeSamFiles_chrom
				run_PICARD_BuildBamIndex run_PICARD_SortSam run_PICARD_AddOrReplaceReadGroups
				run_PICARD_CreateSequenceDictionary run_PICARD_ReorderSam
				run_PICARD_AddOrReplaceReadGroups_gen run_PICARD_SortVcf) ;
            
#Using a library for samtools functions management
use LIB::samtools_management qw(  run_SAMTOOLS_Index run_SAMTOOLS_Sort run_SAMTOOLS_View
									run_SAMTOOLS_View_gen run_SAMTOOLS_Sort_gen run_SAMTOOLS_Index_gen
									run_SAMTOOLS_flagstat run_SAMTOOLS_mpileup);

#Using a library for gatk functions management
use LIB::bwakit_management qw( run_BWAKIT_altmapping run_BWAKIT_samblaster 
													run_BWAKIT_bwa run_BWAKIT_seqtk_seq run_BWAKIT_seqtk_mergepe
													run_STAR);

#Using a library for bedtools functions management
use LIB::bedtools_management qw( run_BEDTOOLS_coverageBed run_BEDTOOLS_intersectBed
																		run_BEDTOOLS_subtractBed);

#Using a library for annotation functions management
use LIB::annotation_management qw( run_ANNOVAR_download_db run_ANNOVAR_convert2annovar
																run_ANNOVAR_annotate_variation run_ANNOVAR_table_annovar
																run_ANNOVAR_download_WG_fasta run_ANNOVAR_retrieveseq
																all_variants_2_VCF checkIncludedAnnovarDBs);
													
#Using a library for managing the output (VCF,..)
use LIB::output_management qw( import_vcf_2_db import_ann_out_2_db get_non_covered_regions
																export_results get_non_cov_reg_with_bedtools
																splicing_var_dist_2_final_out rearrange_rawtabular_out
																build_output_from_annotation 
																generate_vcf_to_annotate output_from_annotation_without_db
																get_gene_annotation print_variants_stats generate_igvsession_xml
																get_allgene_info separate_joint_outputs send_email
																get_gene_coverage_above get_all_samples_coverage_table
																get_all_samples_reads_number get_statistics_tables
																get_all_genes_annotation separate_analysis_from_joint
																tag_panel_genes import_annovarout_2_db import_vgout_2_db
																import_vcf_2_db_fast import_annovarout_2_db_fast
																print_alignment_statistics_table evaluate_quality_of_sequencing);
																
#Using a library for bedtools functions management
use LIB::html_management qw( file_2_html html_page paths_hash_to_html var_stats_2_html
															table_2_html table_list_2_html);

#Using a library for standard utilities								
use LIB::std_lib qw(print_array);

#Using a library for fastq quality management								
use LIB::fastq_quality_management qw(run_trimmomatic  run_fastqc run_trimgalore);

#Using a library for bcftools				
use LIB::bcftools_management qw(run_BCFTOOLS_call run_BCFTOOLS_filter);

#Using a library for vcftools				
use LIB::vcftools_management qw(run_VCFTOOLS_sort);

#Using a library for managing CNVs
use LIB::cnv_management qw( run_exome_depth run_xhmm_pipeline ED_reallocate_output 
							ED_output_preprocessing XHMM_output_preprocessing);

								
#Get environmental variables
#Get the file with the parameters for the module
my $config_file = $ENV{CONFIG_FILE};
#Get the task name to run
my $task = $ENV{TASK};
#Get the program to run
#my $prog_used = $ENV{PROG_USED};
#Gets the job id to be run if it is defined
my $sample_id = $ENV{SAMPLE_ID};
#Gets here the group id
my $analysis_id = $ENV{ANALYSIS_ID};
#Gets the log file
my $log_file = $ENV{LOG_FILE};
#Gets the type of pipeline that is going to be executed
my $pipeline = $ENV{PIPE};

STDOUT->autoflush(1);#This makes STDOUT "hot" in the sense that everything will be print immediately
#print_and_log ("At this time the job id is: $sample_id\n",$log_file);#DEBUGCODE

my $program_name = "programs_runner.pl";

#The steps array is used both to construct the input and output file names
#and to have a reference for the type of pipeline to use

my @steps_array;
my @steps_array2;
my @steps_array3;
my @stats_steps_array;




#Fill the steps arrays using the $pipeline mode given in input
fill_steps_arrays();

if ($task eq 'qc'){
	#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
	run_quality_check_steps($config_file);
}
if ($task eq 'trim'){
	#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
	run_trimming_steps($config_file);
}
#Selection of program to execute
if ($task eq 'align'){
	my $cfg_hash = load_hash($config_file);

	# get the sequencing type
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);		

	#For RNA we have a different pipeline for the alignment task
	if ( $sequencingtype eq 'rna' ){
		#print_and_log( "Running the $task for the RNA sample $sample_id...\n",$log_file);#DEBUGCODE
		run_rna_alignment_steps($config_file);
	}else{
		#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
		run_alignment_steps($config_file);		
	}
}
if ($task eq 'IR'){
					#print_and_log("Steps to run in array1: ",$log_file);#DEBUGCODE
					#foreach my $a (@steps_array){
						#print_and_log("$a - ",$log_file);#DEBUGCODE
						#}
					#print_and_log("\n",$log_file);#DEBUGCODE	
	#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
	#Changing the task name to use the same folder later
	$task = 'refine';
	run_indel_realignment_steps($config_file);
}
if ($task eq 'BR'){
	#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
	#Changing the task name to use the same folder later
	$task = 'refine';
	run_base_recalibration_steps($config_file);
}
if ($task eq 'MS'){
	my $cfg_hash = load_hash($config_file);

	# the sequencing performed must be also  the targeted
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);		
	if ( $pipeline eq 'pipeline'){
		
		#For Amplicons panels the base recalibration has not been ran, so we need to merge and remove duplicates
		if ( $sequencingtype eq 'targeted' ){
			#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
			run_merge_samples_steps($config_file);	
		}
	}
}
if ($task eq 'varcall'){
	#print_and_log( "Running the $task for the sample $sample_id...\n",$log_file);#DEBUGCODE
	run_variant_calling_steps($config_file);
}
if ($task eq 'genot'){
	#print_and_log( "Running the $task for the group $sample_id...\n",$log_file);#DEBUGCODE
	run_genotyping_steps($config_file);
}
#The variant filtering step can be executed in two different ways
# 1. VQSR (Variant Quality Score Recalibration) if:
			# - the sequencing is WES
      #	- there are more than 30 exomes 
      # - the GenotypeGVCFs returned more than 100K variants
# 2. Hard Filtering with VariantFiltration in all other cases 
if ($task eq 'varfilt'){
	my $cfg_hash = load_hash($config_file);
	#Get the type of sequencing for this analysis
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Get the count of distinct samples where the genotype has been predicted (must be in use for frequency calculation and for this type of sequencing)
	my $query_num_samples =  "SELECT ".$cfg_hash->{'db_sample_id'}." FROM ".$cfg_hash->{'db_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id";
	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
										 	 
	my $num_samples = get_count_from_selected($dbh,$query_num_samples," ",$log_file);	
	#Disconnect db
  $dbh->disconnect();
	#################
	  	
	if ( ($sequencingtype eq 'exome') and ($num_samples >= $cfg_hash->{'min_samples_4_joint'}) ){
		print_and_log("Analysis $analysis_id contains $num_samples exomes. Trying VQSR...\n",$log_file);				
		#print_and_log("Running the $task (VARIANT QUALITY SCORE RECALIBRATION) for the group $sample_id...\n",$log_file);#DEBUGCODE
		run_variant_quality_score_recalibration_steps($config_file);	
		run_genotype_refinement_steps($config_file);				
	}else{
		#print_and_log( "Running the $task (HARD FILTERING) for the group $sample_id...\n",$log_file);#DEBUGCODE
		run_var_filt_INDEL_SNP_steps($config_file);		
	}	
}

if ($task eq 'phasing'){
	#print_and_log( "Running the $task for the group $sample_id...\n",$log_file);#DEBUGCODE
	run_phasing_steps($config_file);
}
if ($task eq 'stats'){
	#print_and_log( "Running the $task for the group $sample_id...\n",$log_file);#DEBUGCODE
	#run_tigem_coverage_steps($config_file);
	#run_tigem_flagstats_steps($config_file);
	run_coverage_steps($config_file);
}
if ($task eq 'anstats'){
	#print_and_log( "Running the $task for the group $sample_id...\n",$log_file);#DEBUGCODE
	run_coverage_per_analysis_steps($config_file);
}
if ($task eq 'annotate'){
	#print_and_log( "Running the $task ...\n",$log_file);#DEBUGCODE
	run_annotation_steps($config_file);
}
if ($task eq 'finalout'){
	#Changing the task name to use the same folder later
	#$task = 'finalout';	
	print_and_log( "Running the $task ...\n",$log_file);#DEBUGCODE
	run_import_2_db_steps($config_file);
}

if ($task eq 'fastfout'){
	#Changing the task name to use the same folder later
	$task = 'finalout';
	#print_and_log( "Running the $task ...\n",$log_file);#DEBUGCODE
	run_fast_final_out_steps($config_file);
}

#######CNV Calling
if ($task eq 'cnv'){
	print_and_log( "Running the $task ...\n",$log_file);#DEBUGCODE
	#Changing the task name to use the same folder later
	#$task = 'finalout';
	run_cnv_calling_pipeline($config_file);
}
if ($task eq 'cnvo'){
	print_and_log( "Running the $task ...\n",$log_file);#DEBUGCODE
	#Changing the task name to use the same folder later
	$task = 'cnv';
	run_cnv_finalout_steps($config_file);
}

#Consensus approach tasks
if ($task eq 'VCFB'){
	#Changing the task name to use the same folder later
	$task = 'varfilt';	
	print_and_log( "Running the variant calling with Freebayes ...\n",$log_file);
	run_freebayes_pipeline($config_file);	
}

if ($task eq 'VCSAM'){
	#Changing the task name to use the same folder later
	$task = 'varfilt';	
	print_and_log( "Running the variant calling with Samtools ...\n",$log_file);
	run_samtools_vc_pipeline($config_file);
}


#Independent tasks

if ($task eq 'genedb'){
	print_and_log( "Running the $task ...\n",$log_file);
	run_genedb_creation_steps($config_file);
	#This task uses the main config file. At the end remove it
	delete_file($config_file);		
}
if ($task eq 'geneann'){
	print_and_log( "Running the $task ...\n",$log_file);
	run_allgenes_annotation_steps($config_file);
	#This task uses the main config file. At the end remove it
	delete_file($config_file);		
}
if ($task eq 'reannote'){
	print_and_log( "Running the $task ...\n",$log_file);
	run_reannotate_variants($config_file);
	#This task uses the main config file. At the end remove it
	delete_file($config_file);	
}

if ($task eq 'update_freqs'){
	print_and_log( "Running the $task ...\n",$log_file);
	run_recompute_variants_frequencies($config_file);
	#This task uses the main config file. At the end remove it
	delete_file($config_file);	
}

if ($task eq 'dldatasets'){
	
	print_and_log( "Running the $task ...\n",$log_file);
	run_download_datasets($config_file);
	#This task uses the main config file. At the end remove it
	delete_file($config_file);		
}


=head2 run_freebayes_pipeline

 Title   : run_freebayes_pipeline
 Usage   : run_freebayes_pipeline(  config_file => the config hash
								);

 Function: Executes the variant calling based on the usage of FREEBAYES
 
 Returns :
 
=cut
sub run_freebayes_pipeline{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);

	my $analysis_indication = "\n[an:$analysis_id ](".scalar(localtime)."): ";
	
	#Obtain the analysis name of the current analysis from the database given the group id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);	

									
	my $in_fold_suff = $cfg_hash->{'align_task'};
	my $out_fold_suff = $cfg_hash->{$task.'_task'};
			
	#InFolder is the alignment one at this stage
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	
	#The VCF output from Freebayes
	my $vcfFile	= "$outFolder/$analysis_name\_freebayes.vcf";


	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "$analysis_indication Executing: $query\n",$log_file);#DEBUGCODE
	my $samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#The following array is needed to save all the bam file names
	my @paths_to_write = ();

	
	print_and_log("$analysis_indication [FREEBAYES-pipeline] Getting the BAM files to be used for the calling\n",$log_file);
	foreach my $sample_id (keys %{$samples}){
		#Verifying that there has been an execution of multiple samples
		my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
					    $cfg_hash->{'db_sample_id'},$sample_id);		
		
		my $steps_arr1_in_use;
		my $main_name = "";
		my $params;
		
		if ( $multiple == 1){
			$steps_arr1_in_use = \@stats_steps_array;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);	
			$main_name = $cfg_hash->{'db_sample_name'};		
			#print_and_log("$analysis_indication There are multiple readfiles: $main_name - ".@$steps_arr1_in_use[0]."\n",$log_file);#DEBUGCODE
		}else{
			$steps_arr1_in_use = \@steps_array;	
			my $readfid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
					    $cfg_hash->{'db_sample_id'},$sample_id);				
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readfid);	
			$main_name = $cfg_hash->{'db_readf_name'};		
			#print_and_log("$analysis_indication [FREEBAYES-pipeline] There are not multiple readfiles: $main_name - ".@$steps_arr1_in_use[0]."\n",$log_file);#DEBUGCODE
		}
				
		my $bam_input = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'sort_idx_step'},
					$params->{$main_name},$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'};		
							
		#print_and_log("$analysis_indication [FREEBAYES-pipeline] Using input file: $bam_input\n",$log_file);#DEBUGCODE
		if ( -e $bam_input) {
				
			my $bam_temp = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'sort_idx_step'},
					$params->{$main_name},$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'}.".rg";					
			#Each BAM file needs that the Read Group is set with precision to be used in freebayes
			my $prog_used = $cfg_hash->{'picard_arg_prog'}; 	
		
			print_and_log("$analysis_indication [FREEBAYES-pipeline] Adding Read Group to $bam_input into $bam_temp\n",$log_file);#DEBUGCODE
			run_PICARD_AddOrReplaceReadGroups_gen($cfg_hash,$prog_used,$sample_id,$analysis_id,$bam_input,$bam_temp,$log_file);	
			move($bam_temp, $bam_input) or die "$analysis_indication ERROR: Cannot move $bam_temp in $bam_input\n";
			#Create index for bam file
			print_and_log("$analysis_indication [FREEBAYES-pipeline] Indexing $bam_input into $bam_temp\n",$log_file);#DEBUGCODE	
			run_SAMTOOLS_Index_gen($cfg_hash,$bam_input,$log_file);
			
			push(@paths_to_write,$bam_input);
		}else{
			die "$analysis_indication [FREEBAYES-pipeline] ERROR: Cannot find file $bam_input that is needed for $task. Please run alignment step..\n";
		}
	}		
			
	my $bam_inputs = join(",",@paths_to_write);
	
	###################VARIANT CALLING
	print_and_log("$analysis_indication [FREEBAYES-pipeline] Running FreeBayes using these BAM files: $bam_inputs\n",$log_file);
	run_FREEBAYES_vc($cfg_hash,$bam_inputs,$vcfFile,$analysis_id,$log_file);	

	
	###################VARIANT FILTERING
	#Before to do variant filtering I re-sort the VCF file since with the output from freebayes I always get the following error:
	#https://software.broadinstitute.org/gatk/documentation/article.php?id=1328
	my $prog_used = $cfg_hash->{'picard_sortvcf_prog'};
	print_and_log("$analysis_indication [FREEBAYES-pipeline] Running $prog_used using $vcfFile\n",$log_file);
	my $sortvcf = extract_name($vcfFile,"noext").".sort.".$cfg_hash->{'vcf_ext'};
	run_PICARD_SortVcf($cfg_hash,$prog_used,$vcfFile,$sortvcf,$log_file);
	
	#Finally I want to filter the VCF with the same filters used with GATK 
	#I have to use the -U parameter because of the same problem of before. 
	my $filters = "7";
	my $param_str = "  -U ALLOW_SEQ_DICT_INCOMPATIBILITY ";
	$prog_used = $cfg_hash->{'varfilt_prog'};
	print_and_log("$analysis_indication [FREEBAYES-pipeline] Executing  $prog_used for filtering VCF file $sortvcf\n",$log_file);

	#The output name is again only with freebayes key
	run_GATK_VariantFiltration($cfg_hash,$prog_used,$sortvcf,$vcfFile,$param_str,$filters,$log_file);
	
	#At the end of each variant calling pipeline add the resulting VCF into the configHash
	if ( defined $cfg_hash->{'vc_output_paths'}){
		$cfg_hash->{'vc_output_paths'} .= ",".$vcfFile;
		$cfg_hash->{'vc_algorithms'} .= ",freebayes";
	}else{
		$cfg_hash->{'vc_output_paths'} = $vcfFile;
		$cfg_hash->{'vc_algorithms'} = "freebayes";
	}
	
	if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		print_and_log("$analysis_indication [FREEBAYES-pipeline]  Removing the following temporary files: $sortvcf..\n",$log_file);	
		delete_file($sortvcf);	
	}
}


=head2 run_samtools_vc_pipeline

 Title   : run_samtools_vc_pipeline
 Usage   : run_samtools_vc_pipeline(  
							config_file => the config hash
								);

 Function: Executes the variant calling based on the usage of SAMTOOLS
 
 Returns :
 
=cut
sub run_samtools_vc_pipeline{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);

	#Obtain the analysis name of the current analysis from the database given the group id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);
																
	my $analysis_indication = "\n[an:$analysis_id ](".scalar(localtime)."): ";

	my $in_fold_suff = $cfg_hash->{'align_task'};
	my $out_fold_suff = $cfg_hash->{$task.'_task'};
			
	#InFolder is the alignment one at this stage
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	
	print_and_log("$analysis_indication [SAMTOOLS-pipeline] Getting the BAM files to be used for the calling...\n",$log_file);
		
	#Get the samples ids involved for the analysis
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "$analysis_indication Executing: $query\n",$log_file);#DEBUGCODE
	my $samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#The following array is needed to save all the bam file names
	my @paths_to_write = ();
	
	foreach my $sample_id (keys %{$samples}){
		#Verifying that there has been an execution of multiple samples
		my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
					    $cfg_hash->{'db_sample_id'},$sample_id);		
		
		my $steps_arr1_in_use;
		my $main_name = "";
		my $params;
		
		if ( $multiple == 1){
			$steps_arr1_in_use = \@stats_steps_array;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);	
			$main_name = $cfg_hash->{'db_sample_name'};		
			#print_and_log("$analysis_indication There are multiple readfiles: $main_name - ".@$steps_arr1_in_use[0]."\n",$log_file);#DEBUGCODE
		}else{
			$steps_arr1_in_use = \@steps_array;	
			my $readfid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
					    $cfg_hash->{'db_sample_id'},$sample_id);				
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$readfid);	
			$main_name = $cfg_hash->{'db_readf_name'};		
			#print_and_log("$analysis_indication There are not multiple readfiles: $main_name - ".@$steps_arr1_in_use[0]."\n",$log_file);#DEBUGCODE
		}
				
		my $bam_input = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'sort_idx_step'},
					$params->{$main_name},$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'};		
							
		#print_and_log("$analysis_indication Using input file: $bam_input\n",$log_file);#DEBUGCODE
		if ( -e $bam_input) {
			push(@paths_to_write,$bam_input);
		}else{
			die "$analysis_indication ERROR: Cannot find file $bam_input that is needed for $task. Please run alignment step..\n";
		}
	}		
			
	my $bam_inputs = join(" ",@paths_to_write);
	#print_and_log("$analysis_indication [SAMTOOLS-pipeline] ...DONE!: $bam_inputs\n",$log_file);#DEBUGCODE

	#RUN Samtools MPileup
	print_and_log("$analysis_indication [SAMTOOLS-pipeline] Executing Samtools variant calling using BAM files: $bam_inputs \n",$log_file);
	my $bcf_input	= "$outFolder/$analysis_name\_mpileup.bcf";
	run_SAMTOOLS_mpileup($cfg_hash,$bam_inputs,$bcf_input,$log_file);	

	
	#Run filtering with bcftools
	print_and_log("$analysis_indication [SAMTOOLS-pipeline] Executing bcftools call\n",$log_file);
	my $parameters = " -vmO z ";
	my $bcf_output = "$outFolder/$analysis_name\_bcfcall.bcf";
	run_BCFTOOLS_call($cfg_hash,$bcf_input,$bcf_output,$parameters,$log_file);
	
	print_and_log("$analysis_indication [SAMTOOLS-pipeline] Executing bcftools filter\n",$log_file);
	$parameters = "  ";
	my $filters = "8";
	#Getting the target bed file using the folder for targets and the name contained in the database
	my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
					$cfg_hash->{'db_analysis_id'},$analysis_id);
						
	my $vcf_output = "$outFolder/$analysis_name\_mpileup.vcf.gz";
	run_BCFTOOLS_filter($cfg_hash,$bcf_input,$vcf_output,$parameters,$filters,$target_bed,$log_file);
}


=head2 run_download_datasets

 Title   : run_download_datasets
 Usage   : run_download_datasets(   );

 Function: Fills the steps array using the type of pipeline that must be 
					executed
					
 Returns : 
=cut
sub run_download_datasets{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
		
	################################OLD GENOME DOWNLOAD
	##Download the 2bit format of the genome
	#my $genome_2bit_link = $cfg_hash->{'genome_2bit_link'};
	#print_and_log("Downloading the 2bit format of the genome from $genome_2bit_link into ".$cfg_hash->{'ref_genome_f'}." \n",$log_file);
	##download_file($genome_2bit_link,$cfg_hash->{'ref_genome_f'});
	#my $genome_2bit = extract_name($genome_2bit_link,0);
	
	##Download the 2bittofasta program for the conversion to fasta
	#print_and_log("Downloading 2bittofasta program for the conversion to fasta.\n",$log_file);
	#my $twoBitToFasta_link = $cfg_hash->{'bitToFasta_link'};
	##download_file($twoBitToFasta_link,$cfg_hash->{'main_data_folder'});
	#my $twoBitToFasta = extract_name($twoBitToFasta_link,0);
	
	##Executing the conversion from  2bit format of the Genome to Fasta format
	#chmod 0755, $cfg_hash->{'main_data_folder'}."/".$twoBitToFasta ;
	#my $command = $cfg_hash->{'main_data_folder'}."/".$twoBitToFasta." ".$cfg_hash->{'ref_genome_f'}."/".$genome_2bit.
							#" ".$cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'};
	#print_and_log("Executing $command \n",$log_file);
	##try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;
	##delete_file($cfg_hash->{'main_data_folder'}."/".$twoBitToFasta);
	##delete_file($cfg_hash->{'ref_genome_f'}."/".$genome_2bit);
	
	#print_and_log("Indexing the genome for GATK.\n",$log_file);
	#$command = $cfg_hash->{'samtools_path'}." faidx ".$cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'};
	#try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;	

	#print_and_log("Generate the sequence dictionary for GATK.\n",$log_file);
	#run_PICARD_CreateSequenceDictionary($cfg_hash,$cfg_hash->{'picard_seqdict_prog'},
	#$cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'},
	#$cfg_hash->{'ref_genome_f'}."/".extract_name($cfg_hash->{'hum_gen_align'},'noext').".dict",$log_file);
	
	############################### GENOME DOWNLOAD FROM GATK BUNDLE
	#I download the reference genome and the dictionary from GATK bundle because I was having an error 
	#while running GATK :(Errors about contigs in BAM or VCF files not being properly ordered or sorted)
	#https://software.broadinstitute.org/gatk/documentation/article.php?id=1328. This error disappears
	#when using the genome from the bundle with its  dictionary and when using the PICARD ReorderSam 
	#on the BAM files
	
	my $ref_genome_f = $cfg_hash->{'ref_genome_f'};
	my $gatk_ref_f = $cfg_hash->{'gatk_ref_f'};

	#Download the genome from GATK bundle
	
	my $gatk_genome_link = $cfg_hash->{'gatk_genome_link'};
	my $gatk_genome_com = $ref_genome_f."/".extract_name($gatk_genome_link,'0');
	my $gatk_genome_unc = $ref_genome_f."/".extract_name($gatk_genome_link,'gz');
	print_and_log("Downloading the genome from $gatk_genome_link into ".$cfg_hash->{'ref_genome_f'}." \n",$log_file);
	dl_and_extract($gatk_genome_link,$gatk_genome_com,$gatk_genome_unc,$cfg_hash->{'ref_genome_f'},$log_file);
	
	##Download the index file for the genome
	my $gatk_genome_ind_link = extract_name($gatk_genome_link,'noext').".".$cfg_hash->{'fai_ext'}.".".$cfg_hash->{'gz_ext'};
	my $gatk_genome_ind_com = $ref_genome_f."/".extract_name($gatk_genome_ind_link,'0');
	my $gatk_genome_ind_unc = $ref_genome_f."/".extract_name($gatk_genome_ind_link,'gz');	
	print_and_log("Downloading the genome index from $gatk_genome_ind_link into $ref_genome_f \n",$log_file);
	print_and_log("gatk_genome_ind_link $gatk_genome_ind_link - gatk_genome_ind_com $gatk_genome_ind_com, gatk_genome_ind_unc $gatk_genome_ind_unc, ref_genome_f $ref_genome_f \n",$log_file);
	dl_and_extract($gatk_genome_ind_link,$gatk_genome_ind_com,$gatk_genome_ind_unc,$ref_genome_f,$log_file);	
	
	#Download the dictionary
	my $gatk_genome_dict_link = extract_name($gatk_genome_link,'no2ext').".".$cfg_hash->{'dict_ext'}.".".$cfg_hash->{'gz_ext'};
	my $gatk_genome_dict_com = $ref_genome_f."/".extract_name($gatk_genome_dict_link,'0');
	my $gatk_genome_dict_unc = $ref_genome_f."/".extract_name($gatk_genome_dict_link,'gz');	
	print_and_log("Downloading the genome dictionary from $gatk_genome_dict_link into ".$cfg_hash->{'ref_genome_f'}." \n",$log_file);
	dl_and_extract($gatk_genome_dict_link,$gatk_genome_dict_com,$gatk_genome_dict_unc,$cfg_hash->{'ref_genome_f'},$log_file);		
	
	print_and_log("Indexing the genome for BWA.\n",$log_file);
	my $command = $cfg_hash->{'bwa_path'}." index ".$cfg_hash->{'ref_genome_f'}."/".$cfg_hash->{'hum_gen_align'};
	try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;
	
	
	
	
	#Downloading the datasets to be used in GATK
	print_and_log("Downloading the datasets to be used in GATK...\n",$log_file);
	print_and_log($cfg_hash->{'known_sites_mills_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_mills_link'}){
		my $link = $cfg_hash->{'known_sites_mills_link'};

		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');	
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to  $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($cfg_hash->{'known_sites_mills_link'},'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);		
																		
		#download_file($cfg_hash->{'known_sites_mills_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_mills_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
	}
	print_and_log($cfg_hash->{'known_db_snp_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_db_snp_link'}){
		my $link = $cfg_hash->{'known_db_snp_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($cfg_hash->{'known_db_snp_link'},'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);		
		
				
		#download_file($cfg_hash->{'known_db_snp_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_db_snp_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
	}
	print_and_log($cfg_hash->{'known_hapmap_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_hapmap_link'}){
		my $link = $cfg_hash->{'known_hapmap_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);	
		
				
		#download_file($cfg_hash->{'known_hapmap_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_hapmap_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_1000g_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_1000g_link'}){
		my $link = $cfg_hash->{'known_1000g_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);	
				
		#download_file($cfg_hash->{'known_1000g_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_1000g_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_omni_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_omni_link'}){
		my $link = $cfg_hash->{'known_omni_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);	
		
				
		#download_file($cfg_hash->{'known_omni_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_omni_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_mills_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_mills_link'}){
		my $link = $cfg_hash->{'known_sites_mills_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);
				
		#download_file($cfg_hash->{'known_sites_mills_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_mills_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_db_snp_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_db_snp_link'}){
		my $link = $cfg_hash->{'known_sites_db_snp_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);
						
		#download_file($cfg_hash->{'known_sites_db_snp_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_db_snp_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_hapmap_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_hapmap_link'}){
		my $link = $cfg_hash->{'known_sites_hapmap_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);
				
		#download_file($cfg_hash->{'known_sites_hapmap_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_hapmap_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_1000g_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_1000g_link'}){
		my $link = $cfg_hash->{'known_sites_1000g_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);		
		
		#download_file($cfg_hash->{'known_sites_1000g_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_1000g_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_omni_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_omni_link'}){
		my $link = $cfg_hash->{'known_sites_omni_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);	
				
		#download_file($cfg_hash->{'known_sites_omni_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_omni_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	print_and_log($cfg_hash->{'known_sites_1000gP3_link'}."...",$log_file);
	if (defined $cfg_hash->{'known_sites_1000gP3_link'}){

		my $link = $cfg_hash->{'known_sites_1000gP3_link'};
		my $compr = $gatk_ref_f."/".extract_name($link,'0');
		my $uncompr = $gatk_ref_f."/".extract_name($link,'gz');
 
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);

		$link = extract_name($link,'noext').".".$cfg_hash->{'idxgz_ext'};
		$compr = $gatk_ref_f."/".extract_name($link,'0');
		$uncompr = $gatk_ref_f."/".extract_name($link,'gz');
		#print_and_log("Downloading $compr from $link and extracting into $uncompr to $gatk_ref_f. \n",$log_file);		
		dl_and_extract($link, $compr, $uncompr,$gatk_ref_f,$log_file);	
				
		#download_file($cfg_hash->{'known_sites_1000gP3_link'},$cfg_hash->{'gatk_ref_f'});
		#download_file(extract_name($cfg_hash->{'known_sites_1000gP3_link'},'noext').".".$cfg_hash->{'idxgz_ext'},$cfg_hash->{'gatk_ref_f'});
		}
	
	print_and_log("Datasets have been downladed and genome has been indexed! This step is complete!\n",$log_file);
}


=head2 fill_steps_arrays

 Title   : fill_steps_arrays
 Usage   : fill_steps_arrays(   );

 Function: Fills the steps array using the type of pipeline that must be 
					executed
					
 Returns : 
=cut
sub fill_steps_arrays{
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	my $comma = ",";
	my $pipe_ind = 1;

	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
										
	#Change the index of pipeline
	if ($sequencingtype eq 'exome' or $sequencingtype eq 'targeted'){
		$pipe_ind = 1;
	}
	if ($sequencingtype eq 'rna' ) {
		$pipe_ind = 2;
	}
	#Split the variables from the configuration file and fill the arrays
	@steps_array = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st1'});
	@steps_array2 = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st2'});
	@steps_array3 = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_st3'});
	@stats_steps_array = split($comma,$cfg_hash->{'pipe'.$pipe_ind.'_stats'});
}




##################################################


#################################################

=head2 run_quality_check_steps

 Title   : run_quality_check_steps
 Usage   : run_quality_check_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the quality check
 
 Returns : the path of the last output file
 
=cut
sub run_quality_check_steps{
	my $config_file = shift;
	
	 #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	#READFID is the sampleid
	my $readf_id = $sample_id;


	#Get also the sample_id
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'}
						    ,$cfg_hash->{'db_readf_id'},$readf_id);

	#Get also the sample_name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'}
						    ,$cfg_hash->{'db_sample_id'},$sample_id);
	
	my $analysis_indication = "\n[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
						    						    
	#Get infos for subject read file
	my $params;
	my $readfname = $cfg_hash->{'db_readf_name'};
	#print_and_log("$analysis_indication Getting information for read file id $readf_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
	
	##Here I get information about the flowcell for each fastq
	#my ($instrument,$flowcell, $lane) = get_flowcell_and_lane($cfg_hash,$params->{'fqdir'}."/".$params->{'fq1'},$log_file);
	##..and put them into the db
	##print_and_log("$analysis_indication Updating ".$cfg_hash->{'db_readf_table'}." putting flowcell: $flowcell and lane $lane for $readf_id ..\n",$log_file);#DEBUGCODE
	#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							#$cfg_hash->{'db_readf_instrument'},"'".$instrument."'",$cfg_hash->{'db_readf_id'},$readf_id);
	#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							#$cfg_hash->{'db_readf_flowcell'},"'".$flowcell."'",$cfg_hash->{'db_readf_id'},$readf_id);
	#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							#$cfg_hash->{'db_readf_lane'},"'".$lane."'",$cfg_hash->{'db_readf_id'},$readf_id);
	
	
		
	########################################
	##########FASTQC   						##########
	########################################		
	#Run quality check with one of the available software
	if ($cfg_hash->{'qc_exec'} eq 'YES' ){
		my $prog_used = $cfg_hash->{'qc_prog'};
		my $step = $cfg_hash->{'qc_step'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
		if ($prog_used eq 'fastqc'){
			my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
			my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
			run_fastqc($config_file,$task,$analysis_id,$sample_id,$prog_used,$fq1,$fq2,$log_file);
		}
		#print_and_log("$analysis_indication Updating analysis status for read file $readf_id in ".$cfg_hash-> {'db_name'}." $step=1 ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
				$cfg_hash-> {'db_readf_id'},$readf_id,$step);
		$partDuration = time - $partTime;
		print_and_log("$analysis_indication Finished\t$prog_used\tread file\t$readf_id\t$partDuration\tseconds\n",$log_file);
		$partTime = time;		
	
		##########FILLING OUTLIST FILE##########
		##########RAW QC files#####
	
		#build the names of the html output from FastQC
		my $qc_fq1 = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name($params->{'fq1'},1)."_".$cfg_hash->{'fastqc_html_ext'};
		my $qc_fq2 = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name($params->{'fq2'},1)."_".$cfg_hash->{'fastqc_html_ext'};

		#Check if both the files for read1 and 2 are present
		if ( check_presence($qc_fq1) and check_presence($qc_fq2) ){
			print_and_log( "$analysis_indication Appending $qc_fq2 and $qc_fq2 paths to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);	
			my $line = "$sample_name\t".$params->{$readfname}."_1\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qc_desc'}."\t".$qc_fq1;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
			$line = "$sample_name\t".$params->{$readfname}."_2\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qc_desc'}."\t".$qc_fq2;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
		}else{
			print_and_log( "$analysis_indication WARNING: $qc_fq2 or $qc_fq2 are not present. Please check what is happen. \n",$log_file);	
		}		
		
		##########If the quality is good then the trimming can be avoided.
		#Here I call a function that, given the zip file containing the fastqc report, returns a flag saying if
		#the trimming is suggested or not
		#build the names of the html output from FastQC
		my $qc_fq1_arch = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name($params->{'fq1'},1)."_".$cfg_hash->{'fastqc_zip_ext'};
		my $qc_fq2_arch = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name($params->{'fq2'},1)."_".$cfg_hash->{'fastqc_zip_ext'};	
		print_and_log( "$analysis_indication Using $qc_fq1_arch and $qc_fq2_arch to check if trimming is needed from the quality of reads.... \n",$log_file);
		if ( good_reads_quality($cfg_hash,$qc_fq1_arch,$qc_fq2_arch,$cfg_hash->{$analysis_id.'_qc_out_f'},$log_file) == 1 ){
			update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							$cfg_hash->{'db_readf_goodquality'},"1",$cfg_hash->{'db_readf_id'},$readf_id);			
			print_and_log("$analysis_indication trimming is not needed for this read file.... \n",$log_file);
		}else{
			print_and_log("$analysis_indication trimming is needed for this read file.... \n",$log_file);
		}
			
	}else{print_and_log("$analysis_indication No quality check today..\n",$log_file);}

	

}

=head2 run_trimming_steps

 Title   : run_trimming_steps
 Usage   : run_trimming_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the trimming
 
 Returns : the path of the last output file
 
=cut
sub run_trimming_steps{
	my $config_file = shift;
	
	 #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	#READFID is the sampleid
	my $readf_id = $sample_id;

		
	#Get also the sample_id
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_readf_id'},$readf_id);

	#Get also the sample_name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
	
	my $analysis_indication = "\n[an:$analysis_id sam: $sample_id rf: $readf_id](".scalar(localtime)."): ";
						    						    
	#Get infos for subject read file
	my $params;
	my $readfname = $cfg_hash->{'db_readf_name'};
	#print_and_log("$analysis_indication Getting information for read file id $readf_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);

	#Check the quality if is good to avoid the trimming 
	my $goodquality = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_goodquality'},
						$cfg_hash->{'db_readf_id'},$readf_id);
	if ( $goodquality == 1){
		$cfg_hash->{'trimming'} = "NO";
	}
########################################
	##########TRIMMING						##########
	########################################	
	#Run trimming with one of the available software
	if ($cfg_hash->{'trimming'} eq 'YES' ){
		my $prog_used = $cfg_hash->{'trim_prog'};
		my $step = $cfg_hash->{'trim_step'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);

		###########################
		##In case of Trimgalore
		if ($prog_used eq 'trimgalore'){
			run_trimgalore($config_file,$task,$analysis_id,$readf_id,$prog_used,$log_file);
			##########FILLING OUTLIST FILE##########
			########## Trimmed QC files#####
			#Build the names of the quality check of the two trimmed files
			my $qc_trimmed1 = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name( $params->{'fq1'},'fqgz');
			$qc_trimmed1 .= $cfg_hash->{'trim_out_suff'}."1.".$cfg_hash->{'fastq_ext'}."_".$cfg_hash->{'fastqc_html_ext'};
			my $qc_trimmed2 = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".extract_name( $params->{'fq2'},'fqgz');
			$qc_trimmed2 .= $cfg_hash->{'trim_out_suff'}."2.".$cfg_hash->{'fastq_ext'}."_".$cfg_hash->{'fastqc_html_ext'};

			#Check if both the files for read1 and 2 are present
			if ( check_presence($qc_trimmed1) and check_presence($qc_trimmed2) ){
				print_and_log( "$analysis_indication Appending $qc_trimmed1 and $qc_trimmed2 paths to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
				my $line = "$sample_name\t".$params->{$readfname}."_1\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qctrim_desc'}." \t".$qc_trimmed1;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
				$line = "$sample_name\t".$params->{$readfname}."_2\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qctrim_desc'}."\t".$qc_trimmed2;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
			}else{
				print_and_log( "$analysis_indication WARNING: $qc_trimmed1 or $qc_trimmed2 are not present. Please check what is happen. The next analysis will not be run\n",$log_file);	
			}
		}
		###########################
		##In case of Trimmomatic
		elsif ($prog_used eq 'trimmomatic'){
			#Sample sheet variables (These I will take from the DB)
			my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
			my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
			
			#Output folder
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			#The ouput names that is giving VarGenius	
			my $trimmed1 = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
			my $trimmed2 = $outFolder."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2_".$cfg_hash->{'trim_step'}.$cfg_hash->{'trim_out_ext'};
				
			run_trimmomatic($config_file,$task,$sample_id,$prog_used,$fq1,$fq2,$trimmed1,$trimmed2,$log_file);
			#Run FastQC on the trimmed sequences
			$prog_used = $cfg_hash->{'qc_prog'}; 
			run_fastqc($config_file,$task,$analysis_id,$sample_id,$prog_used,$trimmed1,$trimmed2,$log_file);
			
			#Check if both the files for read1 and 2 are present
			if ( check_presence($trimmed1) and check_presence($trimmed2) ){
				print_and_log( "$analysis_indication Appending $trimmed1 and $trimmed2 paths to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
				#Here sample_id is readf_id
				my $line = "$sample_name\t".$params->{$readfname}."_1\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qctrim_desc'}."\t".$trimmed1;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
				$line = "$sample_name\t".$params->{$readfname}."_2\t".$cfg_hash->{'outlist_qc_page'}."\t".$cfg_hash->{'outlist_qctrim_desc'}." \t".$trimmed2;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
			}else{
				print_and_log( "$analysis_indication WARNING: $trimmed1 or $trimmed2 are not present. Please check what is happen. The next analysis will not be run\n",$log_file);	
			}
		}
		else{
				log_and_exit("$analysis_indication Trimming program is not correct in user configuration file trim_prog = $prog_used. Pleas
						use one of those available!\n",$log_file);
		}
		
		#print_and_log("$analysis_indication Updating analysis status for read file $readf_id  in ".$cfg_hash-> {'db_name'}." $step=1 ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
					$cfg_hash-> {'db_readf_id'},$readf_id,$step);
		$partDuration = time - $partTime;
		print_and_log("$analysis_indication Finished\t$prog_used\tread file\t$readf_id\t$partDuration\tseconds\n",$log_file);
		$partTime = time;
	}
	#If the user does not want to execute the trimming and it has not been already executed
	#VarGenius will copy the fastq files into the qc_out directory
	else{
		print_and_log("$analysis_indication No trimming today.. \n",$log_file);
		my $trim_step = $cfg_hash->{'trim_step'};
		if ($params->{$trim_step} == 0){
			my $fq1 = $params->{'fqdir'}."/".$params->{'fq1'};
			my $fq2 = $params->{'fqdir'}."/".$params->{'fq2'};
			my $fq1_new = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".$params->{$cfg_hash->{'db_readf_name'}}."_R1".$cfg_hash->{'trim_out_ext'};
			my $fq2_new = $cfg_hash->{$analysis_id.'_qc_out_f'}."/".$params->{$cfg_hash->{'db_readf_name'}}."_R2".$cfg_hash->{'trim_out_ext'};
			print_and_log("$analysis_indication Copying $fq1 to $fq1_new.. \n",$log_file);
			copy($fq1,$fq1_new) or die "Cannot copy $fq1 in $fq1_new";
			print_and_log("$analysis_indication Copying $fq2 to $fq2_new.. \n",$log_file);
			copy($fq2,$fq2_new) or die "Cannot copy $fq2 in $fq2_new";
		}
	}

	##################################################
	## UPDATE STATUS
	##################################################
	
	#Update the group table with the status. 
	#At least one of the samples has finished the refinement step. Some other
	#sample could be in a previous task hence here I check what is the value of status
	my $status = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	#the status should be lower that the one we are supposed to put
	if ( $status <	$cfg_hash->{'analysis_qc'}){
		my $fields = $cfg_hash->{'db_analysis_status'};
		my $values = $cfg_hash->{'analysis_alnref'};
		my $fields2 = $cfg_hash->{'db_analysis_id'};
		my $values2 = $analysis_id;
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
	}
	##################################################	
}


=head2 run_rna_alignment_steps

 Title   : run_rna_alignment_steps
 Usage   : run_rna_alignment_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the alignment
 
 Returns : the path of the last output file
 
=cut
sub run_rna_alignment_steps{
	my $config_file = shift;
	
	 #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	#Depending by the steps used this will be the last file path
	my $last_output = "";
	

	#readfid is the sampleid
	my $readf_id = $sample_id;
	
	#Get also the sample_id
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'}
						    ,$cfg_hash->{'db_readf_id'},$readf_id);
						    	
	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id rf:$readf_id](".scalar(localtime)."): ";
	
	my $dataFolder = $cfg_hash->{$analysis_id."_f"}."/".$cfg_hash->{'data_fold'};
	
	#There are two types of alignment BWA and BWA-KIT.
	#For the first case is used BWA,SAMTOOLS and PiCARD MarkDuplicates
	#You can define wha to use in align_prog
	
	#print_and_log("Running alignment steps.\n",$log_file);
	

	#################################################
	#						GATK BEST PRACTICES
	#Alignment strictly following the GATK BestPractices
	#
	#################################################
	
	#elsif ( $cfg_hash->{'align_mode'} eq 'gatk_BP'){
	if ( $cfg_hash->{'align_mode'} eq 'gatk_BP'){
	
		#Launching STAR
		if ($cfg_hash->{'align_prog_exec'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'alignment_prog'};
			my $step = $cfg_hash->{'align_step'};
			my $runMode = "NONE";
			print_and_log("\n$analysis_indication RUNNING $prog_used - First Alignment: ",$log_file);
			#First step of star: first alignment
			my $param_in = "";
			#Set as GenomeDir a directory into the main data folder
			$param_in .= " --genomeDir ".$cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'STAR_genome_f'};
			my $out_suffix = run_STAR($cfg_hash,$prog_used,$readf_id,$analysis_id,$log_file,$task,\@steps_array,$step,$param_in,$runMode);

			
			#Second step of stare: create an indexing using the splice junction information
			print_and_log("\n$analysis_indication RUNNING $prog_used - Indexing based on splice junction information: ",$log_file);

			$runMode = "genomeGenerate";
			#Put the new indexing into the analysis data folder
			#Check if directory exists, otherwise it creates it
			$param_in = $dataFolder."/".$cfg_hash->{'STAR_genome_f'};
			unless(-d $param_in){
				print_and_log("Creating folder $param_in...\n",$log_file);
				mkdir $param_in or die "ERROR: can't create folder $param_in. Check permissions. \n";
			}	
			$param_in = " --genomeDir $param_in ";
			$param_in .= " --sjdbOverhang ".$cfg_hash->{'STAR_sjdbOverhang'}." --sjdbFileChrStartEnd $out_suffix".$cfg_hash->{'STAR_sjout'}." ";
			run_STAR($cfg_hash,$prog_used,$readf_id,$analysis_id,$log_file,$task,\@steps_array,$step,$param_in,$runMode);						
			
						
			#Second step: redo alignment using the new indexing
			print_and_log("\n$analysis_indication RUNNING $prog_used - Alignment based on the new indexing: ",$log_file);
			$runMode = "NONE";
			$param_in = " --genomeDir ".$dataFolder."/".$cfg_hash->{'STAR_genome_f'};
			run_STAR($cfg_hash,$prog_used,$readf_id,$analysis_id,$log_file,$task,\@steps_array,$step,$param_in,$runMode);	
			
			
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							$cfg_hash->{'db_readf_id'},$readf_id,$step);			
			$param_in = "";
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		}
		
		##Launching PICARD  SortSam
		if ($cfg_hash->{'sam_sort_idx'} eq 'YES' ){
				
			#PICARD Sort execution
			print_and_log("$analysis_indication Starting PICARD to Add Read Group, sort and do the bam file \n",$log_file);
			my $prog_used = $cfg_hash->{'picard_arg_prog'};
			my $step = $cfg_hash->{'sort_idx_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
						$cfg_hash-> {'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $input_file = "";#Input file
			
			#print_and_log("Building name using $step..\n",$log_file);#DEBUGCODE		 						
			$input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},\@steps_array).".".$cfg_hash->{'sam_ext'};
	
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);	#DEBUGCODE
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$step.".".$cfg_hash->{'bam_ext'};
			#print_and_log("Output file will be $outFile..\n",$log_file);	#DEBUGCODE

			#print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);			
			#run_PICARD_AddOrReplaceReadGroups($cfg_hash,$prog_used,$readf_id,$analysis_id,$input_file,$outFile,$log_file);
			
			#update status
			#print_and_log("$analysis_indication Updating analysis status for readf_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		}		
		
		#Launching PICARD MarkDuplicates
		if ($cfg_hash->{'mark_rem_dup_exec'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'mark_rem_dup_prog'};
			my $step = $cfg_hash->{'mark_rem_dup_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
						$cfg_hash-> {'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $input_file = "";#Input file
			#print_and_log("Building name using ".$cfg_hash->{'mark_rem_dup_step'}."..\n",$log_file);	#DEBUGCODE		 						

			$input_file = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
							$params->{$main_name},\@steps_array).".".$cfg_hash->{'bam_ext'};
	
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);#DEBUGCODE	
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication  Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'};
			#print_and_log("Output file will be $outFile..\n",$log_file);	#DEBUGCODE		 						
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);		
		  #Execute the generic function to run Picard MarkDuplicates
			run_PICARD_MarkDuplicates_gen($cfg_hash,$prog_used,$input_file,$outFile,$log_file);
						
			#print_and_log("$analysis_indication Updating analysis status for readf_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		
		}
		
		#Launching PICARD BAM indexing
		if ($cfg_hash->{'split_ncigar_exec'} eq 'YES' ){
			
			#PICARD Sort execution
			my $prog_used = $cfg_hash->{'gatk_splitncig_prog'};
			#print_and_log("Starting PICARD $prog_used\n",$log_file);#DEBUGCODE
			my $step = $cfg_hash->{'split_n_cigar_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
						$cfg_hash-> {'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			
			
			#Prepare input file name
			my $input_file = "";#Input file
			#print_and_log("Building name using $step..\n",$log_file);#DEBUGCODE			 						
			$input_file = $outFolder."/".build_input_name_from_executed($params,$step,$params->{$main_name},\@steps_array).".".$cfg_hash->{'bam_ext'};;
			##For the specific case of the indexing we do not have a step more so we 
			##then add into the name the mark_rem_dup_step if it has been executed				
			#if ($cfg_hash->{'mark_rem_dup_exec'} eq 'YES' ){
				#$input_file .= "_".$cfg_hash->{'mark_rem_dup_step'};
			#}
			#Complete the name construction with the extension
			#$input_file .= ".".$cfg_hash->{'bam_ext'};

			#Set output name
			my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
				$params->{$main_name},\@steps_array)."_".$step.".".$cfg_hash->{'bam_ext'};
																	
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);#DEBUGCODE	
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);			
			run_GATK_SplitNCigarReads($cfg_hash,$prog_used,$input_file,$out_bam,$log_file);
			#print_and_log("$analysis_indication Updating analysis status for readf_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id,$step);			
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	

		}
				
	}
	
	save_hash(\$cfg_hash,$config_file);
	return $last_output;
}

=head2 run_alignment_steps

 Title   : run_alignment_steps
 Usage   : run_alignment_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the alignment
 
 Returns : the path of the last output file
 
=cut
sub run_alignment_steps{
	my $config_file = shift;
	
	 #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	#Depending by the steps used this will be the last file path
	my $last_output = "";
	

	#readfid is the sampleid
	my $readf_id = $sample_id;
	
	#Get also the sample_id
	my $sample_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'}
						    ,$cfg_hash->{'db_readf_id'},$readf_id);
						    	
	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id rf:$readf_id](".scalar(localtime)."): ";
	
	#There are two types of alignment BWA and BWA-KIT.
	#For the first case is used BWA,SAMTOOLS and PiCARD MarkDuplicates
	#You can define wha to use in align_prog
	
	#print_and_log("Running alignment steps.\n",$log_file);
	

	#################################################
	#						GATK BEST PRACTICES
	#Alignment strictly following the GATK BestPractices
	#
	#################################################
	
	#elsif ( $cfg_hash->{'align_mode'} eq 'gatk_BP'){
	if ( $cfg_hash->{'align_mode'} eq 'gatk_BP'){
	
		#Launching BWA
		if ($cfg_hash->{'align_prog_exec'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'alignment_prog'};
			my $step = $cfg_hash->{'align_step'};
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);			
			run_BWAKIT_bwa($cfg_hash,$prog_used,$readf_id,$analysis_id,$log_file,$task,\@steps_array,$step);
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
							$cfg_hash->{'db_readf_id'},$readf_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		}
		##Launching PICARD  SortSam
		if ($cfg_hash->{'sam_sort_idx'} eq 'YES' ){
				
			#PICARD Sort execution
			print_and_log("$analysis_indication Starting PICARD Sorting\n",$log_file);
			my $prog_used = $cfg_hash->{'picard_sort_prog'};
			my $step = $cfg_hash->{'sort_idx_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
						$cfg_hash-> {'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $input_file = "";#Input file
			
			#print_and_log("Building name using $step..\n",$log_file);#DEBUGCODE		 						
			$input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},\@steps_array).".".$cfg_hash->{'sam_ext'};
	
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);	#DEBUGCODE
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$step.".".$cfg_hash->{'bam_ext'};
			#print_and_log("Output file will be $outFile..\n",$log_file);	#DEBUGCODE

			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);			
			run_PICARD_SortSam($cfg_hash,$prog_used,$input_file,$outFile,$log_file);
			#update status
			#print_and_log("$analysis_indication Updating analysis status for readf_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		}		
		
		#Launching PICARD MarkDuplicates
		if ($cfg_hash->{'mark_rem_dup_exec'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'mark_rem_dup_prog'};
			my $step = $cfg_hash->{'mark_rem_dup_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},
						$cfg_hash-> {'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			my $input_file = "";#Input file
			#print_and_log("Building name using ".$cfg_hash->{'mark_rem_dup_step'}."..\n",$log_file);	#DEBUGCODE		 						

			$input_file = $outFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
							$params->{$main_name},\@steps_array).".".$cfg_hash->{'bam_ext'};
	
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);#DEBUGCODE	
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication  Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			#Prepare output name
			my $outFile = $outFolder."/".extract_name($input_file,1)."_".$cfg_hash->{'mark_rem_dup_step'}.".".$cfg_hash->{'bam_ext'};
			#print_and_log("Output file will be $outFile..\n",$log_file);	#DEBUGCODE		 						
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);		
		  #Execute the generic function to run Picard MarkDuplicates
			run_PICARD_MarkDuplicates_gen($cfg_hash,$prog_used,$input_file,$outFile,$log_file);
						
			#print_and_log("$analysis_indication Updating analysis status for readf_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		
		}
		
		#Launching PICARD BAM indexing
		if ($cfg_hash->{'sam_sort_idx'} eq 'YES' ){

			#Which program will be used..
			my $prog_used = "";
			#depends on if the reorder is needed or not
			if ( $cfg_hash->{'reorder_sam'} eq 'YES'){
				$prog_used = $cfg_hash->{'picard_reorder_prog'};
			}else{
				$prog_used = $cfg_hash->{'picard_index_prog'};
			}	
			#print_and_log("Starting PICARD $prog_used\n",$log_file);#DEBUGCODE

			my $step = $cfg_hash->{'mark_rem_dup_step'};
			
			#Set the input file name
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$readf_id);
			my $main_name = $cfg_hash->{'db_readf_name'};			
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
			
			
			#Prepare input file name
			my $input_file = "";#Input file
			#print_and_log("Building name using $step..\n",$log_file);#DEBUGCODE			 						
			$input_file = $outFolder."/".build_input_name_from_executed($params,$step,
							$params->{$main_name},\@steps_array);
			#For the specific case of the indexing we do not have a step more so we 
			#then add into the name the mark_rem_dup_step if it has been executed				
			if ($cfg_hash->{'mark_rem_dup_exec'} eq 'YES' ){
				$input_file .= "_".$cfg_hash->{'mark_rem_dup_step'};
			}
			#Complete the name contstruction with the extension
			$input_file .= ".".$cfg_hash->{'bam_ext'};
			
			#print_and_log("Starting $step with $prog_used with $input_file picked from $outFolder..\n",$log_file);#DEBUGCODE	
			if (file_not_present($input_file) > 0 ){ die "$analysis_indication Cannot proceed with $prog_used! Check: $input_file.\n";}
			
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);			
			

			#In some cases I had errors with GATK tools because the BAM file was not indexed in the same way of the 
			#reference FASTA file. Hence I got the suggestion in the forum to use PICARD ReorderSam to reorder the BAM
			#This completely solved and I decided to use it as default behaviour. But still in some cases ReorderSam
			#ends up with errors. Thus I introduced (January 2019) the possibility to choose if to run it or 
			#simply generate the index file with BuildBamIndex
			if ( $cfg_hash->{'reorder_sam'} eq 'YES'){
				#Finally reorder the BAM file to have the same order of the reference file
				#ReorderSam also makes the index
				my $out_file = $input_file.".temp";
				run_PICARD_ReorderSam($cfg_hash,$prog_used,$input_file,$cfg_hash->{'hum_ref'},$out_file,$log_file);
				move($out_file, $input_file) or die "$analysis_indication ERROR: Cannot move $out_file in $input_file\n";
				move($out_file.".bai", $input_file.".bai") or die "$analysis_indication ERROR: Cannot move $out_file in $input_file\n";				
			}else{
				#PICARD BuildBamIndex execution				
				#Just generates the bai index file
				run_PICARD_BuildBamIndex($cfg_hash,$prog_used,$input_file,"NONE",$log_file);				
			}
						
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
			

		}
				
	}
	
	save_hash(\$cfg_hash,$config_file);
	return $last_output;
}


=head2 run_indel_realignment_steps

 Title   : run_indel_realignment_steps
 Usage   : run_indel_realignment_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the indel realignment
 
 Returns : the path of the last output file
 
=cut
sub run_indel_realignment_steps{
	my $config_file = shift;
	my $in_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id](".scalar(localtime)."): ";

	#Check that is not being using the GATK>3  version
	if ( $cfg_hash->{'gatk_ver'} >= 4){
			print_and_log("ERROR: This is weird. You cannot do IndelRealignment with GATK version". $cfg_hash->{'gatk_ver'}." ...\n",$log_file); 
  }	
  
  
	#Creating targets for the realignment
	
	if ($cfg_hash->{'realign_target_creator'} eq 'YES' ){
		my $step = $cfg_hash->{'real_targ_step'};
		my $prog_used = $cfg_hash->{'realign_target_creator_prog'};	
		print_and_log("$analysis_indication  RUNNING $prog_used : ",$log_file);
		run_GATK_RealignerTargetCreator($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,\@steps_array,$step);
		#Save status in the db by searching for the specific sample name
		print_and_log("Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
					$cfg_hash->{'db_readf_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}
	#Creating targets for the realignment
	if ($cfg_hash->{'indel_realigner'} eq 'YES' ){
		my $step = $cfg_hash->{'indel_real_step'};
		my $prog_used = $cfg_hash->{'indel_realigner_prog'};
		print_and_log("$analysis_indication  RUNNING $prog_used : ",$log_file);
		run_GATK_IndelRealigner($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,\@steps_array);
		#Save status in the db by searching for the specific sample name
		print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
					$cfg_hash->{'db_readf_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}
}



=head2 run_base_recalibration_steps

 Title   : run_base_recalibration_steps
 Usage   : run_base_recalibration_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the refinement
					I divided the indel realignment and base recalibration steps in 
					september 2016 because with OPBG we decided to make the base recalibration
					per sample and not anymore per-lane. 
					Hence the refinement must be divided in two different steps
					the indel realingment works per-lane and the base recalibration
					waits indel relaingment for each lane then starts per-sample
					
 Returns : the path of the last output file
 
=cut
sub run_base_recalibration_steps{
	my $config_file = shift;
	my $in_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);

	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id] (".scalar(localtime)."): ";
			
	my $last_step = "";
	my $last_task = "";
	#Check if the refinement has not been executed. In that case, take the output from the alignment
	#Getting the indel_real_step status
	my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'},
																			$cfg_hash->{'db_sample_id'},$sample_id);	

	#Get also the sample_name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						      
  #If indel realignment steps was not been executed previously and not now. Start from the aligned BAM. 
  #Furthermore, without the indel realignment there will be problems in executing the refinement step
  #per-chromosome, hence update the perchrom mode to not execute it
	if ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO' and $indel_real_step==0){
		$last_step	= $cfg_hash->{'sort_idx_step'};
		$last_task	= $cfg_hash->{'align_task'};
		#Set perchrom=0
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},"0",$cfg_hash->{'db_analysis_id'},$analysis_id);
	}#Otherwise take the input from indel Realigner
	else{
		$last_step	= $cfg_hash->{'indel_real_step'};
		$last_task	= $cfg_hash->{'refine_task'};			
	}
	
	#CBase recalibration step
	if ($cfg_hash->{'base_recalibrator'} eq 'YES' ){
		#my $step_needed  = $cfg_hash->{'indel_real_step'};
		#We check that the indel_real_step was performed for all the reads
		#file and then we update the status of indel_real_step=1 in the 'samples' table
		update_sample_status($cfg_hash,$sample_id,$last_step);
		my $step = $cfg_hash->{'base_recal_step'};
		my $prog_used = $cfg_hash->{'base_recalibrator_prog'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);

		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_BaseRecalibrator_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,$last_step,\@steps_array,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
		  run_GATK4_BaseRecalibrator_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,$last_step,\@steps_array,\@steps_array2);
	  }
		

		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}		
	
	#Print reads using the recalibrated bases
	if ($cfg_hash->{'print_reads'} eq 'YES' ){
		#my $step_needed  = $cfg_hash->{'base_recal_step'};
		#We check that the indel_real_step was performed for all the reads
		#file and then we update the status of indel_real_step=1 in the 'samples' table
		#update_sample_status($cfg_hash,$sample_id,$step_needed);
		my $step = $cfg_hash->{'print_reads_step'};
		my $prog_used = $cfg_hash->{'print_reads_prog'};
		print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
		
		
		####################N.B.
		###################HERE I AM WRITING THE CODE TO ADD THE BAM OUTPUT FILE TO THE OUTLIST FILE
		###################BUT I DID NOT CONSIDER THE POSSIBILITY THAT THE PROCESS IS RUN PER-CHROMOSOME!!
		#Out folder
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
		my $params;
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
	
		#Prepare output name for files for data recalibration
		my $out_bam = $outFolder."/".build_input_name_from_executed($params,$step,
						$params->{$cfg_hash->{'db_sample_name'}},\@steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};
		
							
		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_PrintReads_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$out_bam,$log_file,$last_task,$task,\@steps_array,\@steps_array2);
		  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
			run_GATK4_ApplyBQSR($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,\@steps_array,\@steps_array2);
		  }

		########################################
		##########FILLING OUTLIST FILE##########
		########################################		
		#Putting the path to the BAM RECALIRATED into the outlist file
		if ( check_presence($out_bam) ){
			print_and_log( "$analysis_indication Appending $out_bam path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
			my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_recalibratedbam_desc'}."\t".$out_bam;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
		}else{
			print_and_log( "$analysis_indication WARNING: $out_bam is not present. Please check what is happen. \n",$log_file);	
		}
			  		
		#Save status in the for all read files for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for all read files in sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		set_readf_status($cfg_hash,$sample_id,$step);
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log("$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}	

	##################################################
	## UPDATE STATUS
	##################################################
	
	#Update the group table with the status. 
	#At least one of the samples has finished the refinement step. Some other
	#sample could be in a previous task hence here I check what is the value of status
	my $status = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	#the status should be lower that the one we are supposed to put
	if ( $status <	$cfg_hash->{'analysis_alnref'}){
		my $fields = $cfg_hash->{'db_analysis_status'};
		my $values = $cfg_hash->{'analysis_alnref'};
		my $fields2 = $cfg_hash->{'db_analysis_id'};
		my $values2 = $analysis_id;
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
	}
	##################################################	

}


=head2 run_base_recalibration_steps

 Title   : run_base_recalibration_steps
 Usage   : run_base_recalibration_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the refinement
					I divided the indel realignment and base recalibration steps in 
					september 2016 because with OPBG we decided to make the base recalibration
					per sample and not anymore per-lane. 
					Hence the refinement must be divided in two different steps
					the indel realingment works per-lane and the base recalibration
					waits indel relaingment for each lane then starts per-sample
					
 Returns : the path of the last output file
 
=cut
sub run_base_recalibration_stepsOLD{
	my $config_file = shift;
	my $in_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);

	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id] (".scalar(localtime)."): ";
			
	my $last_step = "";
	my $last_task = "";
	#Check if the refinement has not been executed. In that case, take the output from the alignment
	#Getting the indel_real_step status
	my $indel_real_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'indel_real_step'},
																			$cfg_hash->{'db_sample_id'},$sample_id);	
  #If indel realignment steps was not been executed previously and not now. Start from the aligned BAM. 
  #Furthermore, without the indel realignment there will be problems in executing the refinement step
  #per-chromosome, hence update the perchrom mode to not execute it
	if ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO' and $indel_real_step==0){
				$last_step	= $cfg_hash->{'sort_idx_step'};
				$last_task	= $cfg_hash->{'align_task'};
				#Set perchrom=0
				update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},"0",$cfg_hash->{'db_analysis_id'},$analysis_id);
	}#Otherwise take the input from indel Realigner
	else{
				$last_step	= $cfg_hash->{'indel_real_step'};
				$last_task	= $cfg_hash->{'refine_task'};			
	}
	
	#CBase recalibration step
	if ($cfg_hash->{'base_recalibrator'} eq 'YES' ){
		#my $step_needed  = $cfg_hash->{'indel_real_step'};
		#We check that the indel_real_step was performed for all the reads
		#file and then we update the status of indel_real_step=1 in the 'samples' table
		update_sample_status($cfg_hash,$sample_id,$last_step);
		my $step = $cfg_hash->{'base_recal_step'};
		my $prog_used = $cfg_hash->{'base_recalibrator_prog'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);

		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_BaseRecalibrator_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,$last_step,\@steps_array,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
		  run_GATK4_BaseRecalibrator_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,$last_step,\@steps_array,\@steps_array2);
	  }
		

		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}		
	
	#Print reads using the recalibrated bases
	if ($cfg_hash->{'print_reads'} eq 'YES' ){
		#my $step_needed  = $cfg_hash->{'base_recal_step'};
		#We check that the indel_real_step was performed for all the reads
		#file and then we update the status of indel_real_step=1 in the 'samples' table
		#update_sample_status($cfg_hash,$sample_id,$step_needed);
		my $step = $cfg_hash->{'print_reads_step'};
		my $prog_used = $cfg_hash->{'print_reads_prog'};
		print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
		
		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_PrintReads_persample($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,\@steps_array,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
			run_GATK4_ApplyBQSR($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$last_task,$task,\@steps_array,\@steps_array2);
	  }
	  		
		#Save status in the for all read files for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for all read files in sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		set_readf_status($cfg_hash,$sample_id,$step);
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step);
		$partDuration = time - $partTime;
		print_and_log("$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;	
	}	

	##################################################
	## UPDATE STATUS
	##################################################
	
	#Update the group table with the status. 
	#At least one of the samples has finished the refinement step. Some other
	#sample could be in a previous task hence here I check what is the value of status
	my $status = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_status'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	#the status should be lower that the one we are supposed to put
	if ( $status <	$cfg_hash->{'analysis_alnref'}){
		my $fields = $cfg_hash->{'db_analysis_status'};
		my $values = $cfg_hash->{'analysis_alnref'};
		my $fields2 = $cfg_hash->{'db_analysis_id'};
		my $values2 = $analysis_id;
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);			
	}
	##################################################	

}


=head2 set_readf_status

 Title   : set_readf_status
 Usage   : set_readf_status(  config_file => the config hash
								);

 Function: Given a sample id updates the given field of all the read
 files associated to that sample
 
 Returns :
 
=cut
sub set_readf_status{
	my $cfg_hash = shift;
	my $sample_id = shift;
	my $step_needed = shift;

	#Here we get the read file ids associated with that sample (string)
	my $distinct_readf = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $params;
	#If there are more files associated to this sample
	if ( (scalar keys %{$distinct_readf}) > 0){
		#$step_needed  = $cfg_hash->{'mrdup_groups_step'};
		#Check from the DB if the step has been performed for all the samples
		foreach my $dist_readf ( keys %{$distinct_readf}){
			#Save status in the db by searching for the specific sample name
			print_and_log("Updating status for readfid $dist_readf in ".$cfg_hash-> {'db_name'}." ($step_needed=1) ...\n",$log_file);
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
						$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},
						$cfg_hash->{'db_readf_id'},$dist_readf,$step_needed);
		}

	}
}

=head2 update_sample_status

 Title   : update_sample_status
 Usage   : update_sample_status(  config_file => the config hash
								);

 Function: Checks that the given step has finished for all the files and 
					updates the samples table with the final step
 
 Returns : 
 
=cut
sub update_sample_status{
	my $cfg_hash = shift;
	my $sample_id = shift;
	my $step_needed = shift;

	my $analysis_indication = "\n[an:$analysis_id ] (".scalar(localtime)."): [sam: $sample_id]";
	
	#Here we get the read file ids associated with that sample (string)
	my $distinct_readfs = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
	my $main_name = $cfg_hash->{'db_readf_name'};
	my $params;
	#If there are more files associated to this sample
	if ( (scalar keys %{$distinct_readfs}) > 0){
		#$step_needed  = $cfg_hash->{'mrdup_groups_step'};
		#Check from the DB if the step has been performed for all the samples
		foreach my $distinct_readf ( keys %{$distinct_readfs}){
			#print_and_log("Checking status for readfile $dist_sample...\n",$log_file);#DEBUGCODE
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$distinct_readf);
			log_and_exit( "$analysis_indication ERROR: The step $step_needed is needed for readfile ".
					$params->{$main_name}." and has not been completed (".$cfg_hash->{'db_readf_id'}.
					" = $distinct_readf)\n", $log_file) unless  ($params->{$step_needed} == 1) ;
		}
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step_needed=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step_needed);
	}
}

=head2 update_group_status

 Title   : update_group_status
 Usage   : update_group_status(  config_file => the config hash
								);

 Function: Checks that the given step has finished for all the samples and 
					updates the groups table with the final step
 
 Returns : 
 
=cut
sub update_group_status{
	my $cfg_hash = shift;
	my $analysis_id = shift;
	my $step_needed = shift;
	
	my $main_name = $cfg_hash->{'db_sample_name'};
	#Here we get the read file ids associated with that sample (string)
	my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $params;
	#If there are more files associated to this sample
	#if ( (scalar keys %{$distinct_samples}) > 1){
	if ( (scalar keys %{$distinct_samples}) > 0){
		#$step_needed  = $cfg_hash->{'mrdup_groups_step'};
		#Check from the DB if the MergeGroups step has been performed for all the samples
		foreach my $dist_sample ( keys %{$distinct_samples}){
			
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$dist_sample);
			log_and_exit( "ERROR: The step $step_needed is needed for sample ".
					$params->{$main_name}." and has not been completed (".$params->{$cfg_hash->{'db_sample_id'}}.
					" = $dist_sample)", $log_file) unless  ($params->{$step_needed} == 1) ;
		}
		#Save status in the db by searching for the specific group
		print_and_log("Updating analysis status for group $analysis_id in ".$cfg_hash-> {'db_name'}." ($step_needed=1) ...\n",$log_file);
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step_needed);
	}
}

=head2 run_mark_duplicates_steps

 Title   : run_mark_duplicates_steps
 Usage   : run_mark_duplicates_steps(  config_file => the config hash
								);

 Function: THIS PROCEDURE BELONGS TO PIPELINE 2
						Where when you have the same sample divided in different lanes
						then to merge the lanes together we use the PrintReads step in GATK
						and here we have only to mark and remove duplicates.
						
						It also checks if the files from the refinement step were subdivided
						per-chromosome and if this happens the multiple files will be
						merged per each chromosome sub-division.
						
						Finally the bam files are indexed with samtools index
						
						NEEDS print_reads_step TO BE EXECUTED
						
 Returns : the path of the last output file
 
=cut
sub run_mark_duplicates_steps{
	my $config_file = shift;

	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'}
						    ,$cfg_hash->{'db_sample_id'},$sample_id);
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);
	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'}
						    ,$cfg_hash->{'db_sample_id'},$sample_id);		
	
	my $last_refine_step = "";
	#If the base recalibration has not been performed during the refine step
	#then the last step will be the Indel Realignment that provides BAM files as well
	# the sequencing performed must be also  the targeted
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'}
									,$cfg_hash->{'db_analysis_id'},$analysis_id);		
	if ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
			$last_refine_step	= $cfg_hash->{'indel_real_step'};
	}else{
			$last_refine_step	= $cfg_hash->{'print_reads_step'};
	}
		
	
	if ( $multiple == 1){
	#Here we need the print reads step has been done. 
		my $step_needed = $last_refine_step;
		#We check that the print reads step was performed for all the reads
		#file and then we update the status of precalread=1 in the 'samples' table
		#update_sample_status($cfg_hash,$sample_id,$step_needed);
		my $step  = "";
		my $in_fold_suff = 'refine';#The folder that will be used to take the input
		my $out_fold_suff = 'refine';#The folder that will be used for the output
			

		#Applying MarkDuplicates
		if ($cfg_hash->{'mrdup_groups'} eq 'YES' ){	
			my $prog_used = $cfg_hash->{'mrdup_groups_prog'};
			$step = $cfg_hash->{'mrdup_groups_step'};
			my $previous_task = 'refine';#Used to know the folder where to take the input bam files	
			my $prev_step = $last_refine_step;

			print_and_log("[an:$analysis_id] (".scalar(localtime)."): [sam: $sample_id] RUNNING $prog_used..\n",$log_file);
			print_and_log("Starting $prog_used\n",$log_file);
			
			#Here I give in input the stepsArray2 for samples for both the input and output
			run_PICARD_MarkDuplicates_chrom($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,$step,$prev_step,\@steps_array2,\@steps_array2,$perchrom);
			#Save status in the db by searching for the specific sample name
			print_and_log("Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
			$partDuration = time - $partTime;
			print_and_log( "Finished\t$prog_used\tsample\t$sample_id\t$partDuration\tseconds\n",$log_file);
			$partTime = time;
		}
		
		#INDEXING of the BAM file
		if ($cfg_hash->{'index_after_merge'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'sort_prog'};
			print_and_log("[an:$analysis_id] (".scalar(localtime)."): [sam: $sample_id] RUNNING $prog_used..\n",$log_file);
			print_and_log("Starting $prog_used\n",$log_file);	
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
							
			#inFolder for the bam files
			my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
			#If it was performed per chromosome
			my $numExec = 1;
			if ( $perchrom == 1 ){
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
				#for (my $exec = 1; $exec <= $numExec; $exec++){
					#Get a single bam file from separated
					my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
										$sample_name,\@steps_array2)."_".$step."_".$exec.".".$cfg_hash->{'bam_ext'};
					print "Starting SAMTOOLS index for $input_file\n";		
					#Checks if the SAM output exists and is non-zero
					if ( file_not_present($input_file) > 0){ log_and_exit( "ERROR: Cannot find $input_file. Please check merge_groups task\n", $log_file)}
					run_SAMTOOLS_Index_gen($cfg_hash,$input_file,$log_file);
				}
			}#Otherwise, if it has not been performed per-chromosome
			else{
				#Get a single bam file from separated
				my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
									$sample_name,\@steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};
										
				#Checks if the SAM output exists and is non-zero
				if ( file_not_present($input_file)  > 0){ log_and_exit( "ERROR: Cannot find $input_file. Please check merge_groups task\n", $log_file)}
				run_SAMTOOLS_Index_gen($cfg_hash,$input_file,$log_file);
			}
		}
	}else{ print_and_log( "Step ".$cfg_hash->{'mrdup_groups'}." for sample $sample_id is not needed\n",$log_file);}					    
	
}


=head2 run_merge_samples_steps

 Title   : run_merge_samples_steps
 Usage   : run_merge_samples_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the merge
						This step is called when one sample file has been divided
						in more files. For example when you have the same sample with
						different lanes. 
						With this program all the read files coming from the same sample 
						will be merged in a single one.
						
						It also checks if the files from the refinement step were subdivided
						per-chromosome and if this happens the multiple files will be
						merged for each chromosome sub-division.
						
						Both input and output folder will be refined
												
						Needs print_reads_step to be ran
						
 Returns : the path of the last output file
 
=cut
sub run_merge_samples_steps{
	my $config_file = shift;

	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id] (".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	#Getting the sample name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);		
	
	my $last_step = "";
	my $last_task = "";
	#If the base recalibration has not been performed during the refine step
	#then the last step will be the Indel Realignment that provides BAM files as well
	# the sequencing performed must be also  the targeted
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'}
									,$cfg_hash->{'db_analysis_id'},$analysis_id);		
	
	#The refinement has not been executed
	if (( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' )){
				$last_step	= $cfg_hash->{'sort_idx_step'};
				$last_task	= $cfg_hash->{'align_task'};
	}
	#Only Indel Realignment
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
			$last_step	= $cfg_hash->{'indel_real_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#Only Base Recalibrator 
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#All refinement steps are executed
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}else{
		log_and_exit("$analysis_indication ERROR: There is something wrong with your setting. This refinement is not allowed..\n",$log_file);
	}
		
	
	if ( $multiple == 1){
	#Here we need the print reads step has been done. 
		my $step_needed = $last_step;
		#We check that the print reads step was performed for all the reads
		#file and then we update the status of precalread=1 in the 'samples' table
		update_sample_status($cfg_hash,$sample_id,$step_needed);
		my $step  = "";
		my $in_fold_suff = $last_task;#The folder that will be used to take the input
		my $out_fold_suff = $last_task;#The folder that will be used for the output
		#print_and_log("$analysis_indication multiple: $multiple and mergesam=".$cfg_hash->{'mergesam'} eq 'YES'."\n",$log_file);#DEBUGCODE
		#Creating targets for the realignment
		if ($cfg_hash->{'mergesam'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'mergesam_prog'};
			$step = $cfg_hash->{'mergesam_step'};
			my $prev_step = $last_step;
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			run_PICARD_MergeSamFiles_chrom($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,
							$step,$prev_step,\@steps_array,\@steps_array2,$perchrom);
			#Save status in the db by searching for the specific sample name
			print_and_log("$analysis_indication Updating sample status in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;
		}
		##Creating targets for the realignment
		#if ($cfg_hash->{'mrdup_groups'} eq 'YES' ){	
			#my $prog_used = $cfg_hash->{'mrdup_groups_prog'};
			#$step = $cfg_hash->{'mrdup_groups_step'};
			#my $previous_task = $last_task;#Used to know the folder where to take the input bam files	
			#my $prev_step = $cfg_hash->{'mergesam_step'};

			#print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			#run_PICARD_MarkDuplicates_chrom($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,$step,$prev_step,\@steps_array2,\@steps_array2,$perchrom);
			##Save status in the db by searching for the specific sample name
			##print_and_log("Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			#update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					#$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
					#$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
			#$partDuration = time - $partTime;
			#print_and_log( "Finished $prog_used $partDuration seconds\n",$log_file);
			#$partTime = time;
		#}
		#INDEXING of the BAM file
		if ($cfg_hash->{'index_after_merge'} eq 'YES' ){
			my $prog_used = $cfg_hash->{'sort_prog'};
			print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
							
			#inFolder for the bam files
			my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
			#If it was performed per chromosome
			my $numExec = 1;
			if ( $perchrom == 1 ){
				##Using also Chromosomex X and Y that Romans did not use!!!
				$numExec = $cfg_hash->{'chromosomes_num'};
				my @chroms = 1..$numExec;
				push(@chroms,split(",",$cfg_hash->{'sex_chromosomes'}));
							
				#An array will contain all the commands to be executed
				my @commands = ();
				
				#Executing one GATK command for each target chromosome using bed files
				foreach my $exec (@chroms){
					#Get a single bam file from separated
					my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
										$sample_name,\@steps_array2)."_".$step."_".$exec.".".$cfg_hash->{'bam_ext'};
					#print "Starting SAMTOOLS index for $input_file\n";		
					#Checks if the SAM output exists and is non-zero
					if ( file_not_present($input_file) > 0){ log_and_exit( "$analysis_indication ERROR: Cannot find $input_file. Please check merge_groups task\n", $log_file)}
					run_SAMTOOLS_Index_gen($cfg_hash,$input_file,$log_file);
				}
			}#Otherwise, if it has not been performed per-chromosome
			else{
				#Get a single bam file from separated
				my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,
									$sample_name,\@steps_array2)."_".$step.".".$cfg_hash->{'bam_ext'};
										
				#Checks if the SAM output exists and is non-zero
				if ( file_not_present($input_file)  > 0){ log_and_exit( "$analysis_indication ERROR: Cannot find $input_file. Please check merge_groups task\n", $log_file)}
				run_SAMTOOLS_Index_gen($cfg_hash,$input_file,$log_file);
			}
		}
	}else{ print_and_log( "$analysis_indication Step mergesam for sample $sample_id is not needed\n",$log_file);}					    
	
}

	
	
=head2 run_variant_calling_steps

 Title   : run_variant_calling_steps
 Usage   : run_variant_calling_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the variant calling
 
 Returns : the path of the last output file
 
=cut
sub run_variant_calling_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed = "";
	#my $previous_task = "refine";


	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id](".scalar(localtime)."): ";

	#Obtain the analysis name of the current analysis from the database given the group id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);	
																		
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Starting variant calling task \n",$log_file);
	#print_and_log("####################################################...\n",$log_file);
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
						    $cfg_hash->{'db_sample_id'},$sample_id);

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	
	
	my $last_step = "";
	my $last_task = "";
		

	#Get sequencing type
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
				$cfg_hash->{'db_analysis_id'},$analysis_id);		

	#The refinement has not been executed
	#the sequencing performed must be also  the targeted to go on
	if ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
				$last_step	= $cfg_hash->{'sort_idx_step'};
				$last_task	= $cfg_hash->{'align_task'};
	}
	#Only Indel Realignment
	#If the base recalibration has not been performed during the refine step
	#then the last step will be the Indel Realignment that provides BAM files as well
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
			$last_step	= $cfg_hash->{'indel_real_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#Only Base Recalibrator 
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#All refinement steps are executed
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}else{
		log_and_exit("$analysis_indication ERROR: There is something wrong with your setting. This refinement is not allowed:\n".
			     "realign_target_creator: ".$cfg_hash->{'realign_target_creator'}." ".
			     "indel_realigner: ".$cfg_hash->{'indel_realigner'}." ".
			     "base_recalibrator: ".$cfg_hash->{'base_recalibrator'}." ".
			     "print_reads: ".$cfg_hash->{'print_reads'}.
			"seqtype: $sequencingtype..\n",$log_file);
	}
	#print_and_log("$analysis_indication As from the pipeline you want to execute: last step: $last_step - last_task: $last_task - seqtype: $sequencingtype\n",$log_file);#DEBUGCODE
	
	#If there are multiple read files we need a step of merge					    						    
	if ( $multiple == 1){
		$step_needed  = $cfg_hash->{'mrdup_groups_step'};
		print_and_log("$analysis_indication Multiple samples are present. Taking the input from $last_task folder...\n",$log_file);
	}else {
		$step_needed  = $last_step;
		update_sample_status($cfg_hash,$sample_id,$step_needed);
		print_and_log("$analysis_indication Multiple samples are not present. Taking the input from $last_task folder...\n",$log_file);
	}

	#Getting the flag indicating if the joint genotype has to be run
	my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					       $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $out_ext;
	#If it has been run, take the output form the genotype, else from the varcalling
	if ($dogenot){
		$out_ext = $cfg_hash->{'gvcf_ext'};
		print_and_log("$analysis_indication The genotype step has to be run. The output extension will be: $out_ext...\n",$log_file);
	}else{
		$out_ext = $cfg_hash->{'vcf_ext'};
		print_and_log("$analysis_indication The genotype will not be run. The output extension will be: $out_ext...\n",$log_file);
	}
				
	#Detecting the variants with HaplotypeCaller
	if ($cfg_hash->{'varcall'} eq 'YES' ){

		
		my $prog_used = $cfg_hash->{'varcall_prog'};	
		my $step = $cfg_hash->{'varcall_step'};
		my $in_fold_suff = $last_task;
		my $out_fold_suff = $task;
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);					

		#Construct input and output file names     
		my $steps_arr1_in_use;
		my $paramsIn;
		my $paramsOut;

		#Case IR only or no refinement step: there is a single file and only IndelRealignment or no refinement step was executed: 
		#construct the name using readfile level fields
		if ( $step_needed eq $cfg_hash->{'indel_real_step'} or $step_needed eq $cfg_hash->{'sort_idx_step'} ){
			$steps_arr1_in_use = \@steps_array;
			my $distinct_samples = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
									$cfg_hash->{'db_analysis_id'},$analysis_id,$cfg_hash->{'db_sample_id'},$sample_id);
			foreach my $dist_sample ( keys %{$distinct_samples}){
				getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},$dist_sample);
			}
		}	
		#CASE: IR-BR or IR-BR-M or IR-M: there are multiple read files or the last step is BR
		#construct the name using sample level fields
		elsif ($step_needed eq $cfg_hash->{'mrdup_groups_step'} or $step_needed eq $cfg_hash->{'print_reads_step'}){
			$steps_arr1_in_use = \@steps_array2;
			getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);
		}
		
		#Get the sample parameters from the the db for the output file
		getSampleConfiguration_locked(\$paramsOut,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

																		
		#InFolder is the alignment one at this stage
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
		
		#Getting the sample name
		my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
									$cfg_hash->{'db_sample_id'},$sample_id);		
																			
		#Set the input and output file names
		my $input_file = $inFolder."/".build_input_name_from_executed($paramsIn,$step,$sample_name,$steps_arr1_in_use).".".$cfg_hash->{'bam_ext'};
		my $out_file = $outFolder."/".build_input_name_from_executed($paramsOut,$step,$sample_name,\@steps_array2)."_".$step.".".$out_ext;
					
		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
				run_GATK_HaplotypeCaller($cfg_hash,$prog_used,$sample_id,$analysis_id,$input_file,$out_file,$inFolder,$outFolder,$log_file);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				run_GATK4_HaplotypeCaller($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,$step,$step_needed,\@steps_array,\@steps_array2);
	  }
	  		
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
				$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;
	}
	
	#Merging gvcf files if they have been separated per-chromosome
	if ( ($cfg_hash->{'catvar'} eq 'YES') and ( $perchrom == 1) ){	
		my $prog_used = $cfg_hash->{'catvar_prog'};
		my $step = $cfg_hash->{'catvar_step'};
		print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
		
		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_CatVariants($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
			run_GATK4_CatVariants($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,\@steps_array2);
	  }		
		
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
				$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;
	}

	##At this point VarGenius may remove all the files from different read files if samples are separated
	#Haplotype caller will pick the BAM files from the alignment folder only if the refinement steps have not been 
	#run. Anyway there is no difference in removing them now or before...
	if ( $multiple == 1){
		my $step_needed = $cfg_hash->{'sort_idx_step'};
		my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};		
		
		#Here get the samples names to be merged
		my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$step_needed,$inFolder,\@steps_array,$cfg_hash->{'bam_ext'});			
		#HERE I CAN REMOVE THE SAM FILES BECAUSE THEY HAVE BEEN ALREADY USED
		#and we are sure that the bam_input is present!
		#After generating the merged bam remove the sam file 
		if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			#Check if this is the last task using the BAM files 
			my $bamfilesused = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
			if ( $bamfilesused ==  ($cfg_hash->{'tasks_using_bam'} - 1) )	{
				foreach my $bam_path (@paths){
					print_and_log("$analysis_indication Removing the alignment $bam_path file..\n",$log_file);	
					delete_file($bam_path);				
					delete_file($bam_path.".".$cfg_hash->{'bai_ext'});				
				}
				$bamfilesused++;
				update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
								$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},$bamfilesused,$cfg_hash->{'db_sample_id'},$sample_id);								
			}
		}
		
		
	}
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Variant calling task completed.\n",$log_file);
	#print_and_log("####################################################...\n",$log_file);	
}


=head2 run_variant_calling_steps

 Title   : run_variant_calling_steps
 Usage   : run_variant_calling_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the variant calling
 
 Returns : the path of the last output file
 
=cut
sub run_variant_calling_stepsOLD{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed = "";
	#my $previous_task = "refine";


	my $analysis_indication = "\n[an:$analysis_id sam:$sample_id](".scalar(localtime)."): ";
		
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Starting variant calling task \n",$log_file);
	#print_and_log("####################################################...\n",$log_file);
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'}
						    ,$cfg_hash->{'db_sample_id'},$sample_id);

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'}
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);
	
	
	my $last_step = "";
	my $last_task = "";
		

	#Get sequencing type
	my $sequencingtype = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
				$cfg_hash->{'db_analysis_id'},$analysis_id);		

	#The refinement has not been executed
	#the sequencing performed must be also  the targeted to go on
	if ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
				$last_step	= $cfg_hash->{'sort_idx_step'};
				$last_task	= $cfg_hash->{'align_task'};
	}
	#Only Indel Realignment
	#If the base recalibration has not been performed during the refine step
	#then the last step will be the Indel Realignment that provides BAM files as well
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'NO' and $cfg_hash->{'print_reads'} eq 'NO'
			and $sequencingtype eq 'targeted' ){
			$last_step	= $cfg_hash->{'indel_real_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#Only Base Recalibrator 
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'NO' and $cfg_hash->{'indel_realigner'} eq 'NO'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}
	#All refinement steps are executed
	elsif ( $cfg_hash->{'realign_target_creator'} eq 'YES' and $cfg_hash->{'indel_realigner'} eq 'YES'
			and $cfg_hash->{'base_recalibrator'} eq 'YES' and $cfg_hash->{'print_reads'} eq 'YES'){
			$last_step	= $cfg_hash->{'print_reads_step'};
			$last_task	= $cfg_hash->{'refine_task'};
	}else{
		log_and_exit("$analysis_indication ERROR: There is something wrong with your setting. This refinement is not allowed:\n".
			     "realign_target_creator: ".$cfg_hash->{'realign_target_creator'}." ".
			     "indel_realigner: ".$cfg_hash->{'indel_realigner'}." ".
			     "base_recalibrator: ".$cfg_hash->{'base_recalibrator'}." ".
			     "print_reads: ".$cfg_hash->{'print_reads'}.
			"seqtype: $sequencingtype..\n",$log_file);
	}
	#print_and_log("$analysis_indication As from the pipeline you want to execute: last step: $last_step - last_task: $last_task - seqtype: $sequencingtype\n",$log_file);#DEBUGCODE
	
	#If there are multiple read files we need a step of merge					    						    
	if ( $multiple == 1){
		$step_needed  = $cfg_hash->{'mrdup_groups_step'};
		#This code was used when I still used to make the merge of samples with samtools merge.
		#But now this is performed by printRecalibratedReads...
		#if ( $pipeline eq 'pipeline') {
			##Verifying that there has been an execution of merge samples
			#my $merge_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mergesam_step'}
							    #,$cfg_hash->{'db_sample_id'},$sample_id);		
			#if ( $merge_executed == 0 ){
				#log_and_exit("ERROR: You need to execute $step_needed before the variant calling..\n",$log_file);
			#}
		#}
		print_and_log("$analysis_indication Multiple samples are present. Taking the input from $last_task folder...\n",$log_file);
	}else {
		$step_needed  = $last_step;
		update_sample_status($cfg_hash,$sample_id,$step_needed);
		print_and_log("$analysis_indication Multiple samples are not present. Taking the input from $last_task folder...\n",$log_file);
	}

	#Getting the flag indicating if the joint genotype has to be run
	my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					       $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
							$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $out_ext;
	#If it has been run, take the output form the genotype, else from the varcalling
	if ($dogenot){
		$out_ext = $cfg_hash->{'gvcf_ext'};
		print_and_log("$analysis_indication The genotype step has to be run. The output extension will be: $out_ext...\n",$log_file);
	}else{
		$out_ext = $cfg_hash->{'vcf_ext'};
		print_and_log("$analysis_indication The genotype will not be run. The output extension will be: $out_ext...\n",$log_file);
	}
				
	#Detecting the variants with HaplotypeCaller
	if ($cfg_hash->{'varcall'} eq 'YES' ){
		my $prog_used = $cfg_hash->{'varcall_prog'};	
		my $step = $cfg_hash->{'varcall_step'};
		my $in_fold_suff = $last_task;
		my $out_fold_suff = $task;
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);					

		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
				run_GATK_HaplotypeCaller($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,$step,$step_needed,\@steps_array,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				run_GATK4_HaplotypeCaller($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$in_fold_suff,$out_fold_suff,$step,$step_needed,\@steps_array,\@steps_array2);
	  }
	  		
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
				$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;
	}
	
	#Merging gvcf files if they have been separated per-chromosome
	if ( ($cfg_hash->{'catvar'} eq 'YES') and ( $perchrom == 1) ){	
		my $prog_used = $cfg_hash->{'catvar_prog'};
		my $step = $cfg_hash->{'catvar_step'};
		print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
		
		#Select GATK depending by the version
		if ( $cfg_hash->{'gatk_ver'} < 4){
			run_GATK_CatVariants($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,\@steps_array2);
	  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
			run_GATK4_CatVariants($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,\@steps_array2);
	  }		
		
		#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
				$cfg_hash->{'db_sample_id'},$sample_id,$step,$cfg_hash->{'db_analysis_id'},$analysis_id);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;
	}

	##At this point VarGenius may remove all the files from different read files if samples are separated
	#Haplotype caller will pick the BAM files from the alignment folder only if the refinement steps have not been 
	#run. Anyway there is no difference in removing them now or before...
	if ( $multiple == 1){
		my $step_needed = $cfg_hash->{'sort_idx_step'};
		my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};		
		
		#Here get the samples names to be merged
		my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$step_needed,$inFolder,\@steps_array,$cfg_hash->{'bam_ext'});			
		#HERE I CAN REMOVE THE SAM FILES BECAUSE THEY HAVE BEEN ALREADY USED
		#and we are sure that the bam_input is present!
		#After generating the merged bam remove the sam file 
		if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			#Check if this is the last task using the BAM files 
			my $bamfilesused = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
			if ( $bamfilesused ==  ($cfg_hash->{'tasks_using_bam'} - 1) )	{
				foreach my $bam_path (@paths){
					print_and_log("$analysis_indication Removing the alignment $bam_path file..\n",$log_file);	
					delete_file($bam_path);				
					delete_file($bam_path.".".$cfg_hash->{'bai_ext'});				
				}
				$bamfilesused++;
				update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
								$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},$bamfilesused,$cfg_hash->{'db_sample_id'},$sample_id);								
			}
		}
		
		
	}
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Variant calling task completed.\n",$log_file);
	#print_and_log("####################################################...\n",$log_file);	
}

=head2 run_genotyping_steps

 Title   : run_genotyping_steps
 Usage   : run_genotyping_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the variant calling
 
 Returns : the path of the last output file
 
=cut
sub run_genotyping_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed = "";

	#Obtain the analysis name of the current analysis from the database given the group id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);	
																	
	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
	
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Starting genotyping task...\n",$log_file);
	#print_and_log("####################################################...\n",$log_file);
	

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    
	#If there are more files associated to this sample
	if ( $perchrom == 1 ){
		$step_needed = $cfg_hash->{'catvar_step'};
	}else{
		$step_needed = $cfg_hash->{'varcall_step'};
	}
	
	#Checks that step_needed is executed for all the samples and updates that step
	#in the group table
	update_group_status($cfg_hash,$analysis_id,$step_needed);

		
	#Creating targets for the realignment
	if ($cfg_hash->{'genot'} eq 'YES' ){
		my $prog_used = $cfg_hash->{'genot_prog'};	
		my $step = $cfg_hash->{'genot_step'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);

		#Get the count of distinct samples where the genotype has been predicted (must be in use for frequency calculation and for this type of sequencing)
		my $query_num_samples =  "SELECT ".$cfg_hash->{'db_sample_id'}." FROM ".$cfg_hash->{'db_sample_table'}.
												 " WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id";
		#################DB
		#Connect to database
		my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
		# PERL DBI CONNECT
		my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
												 
		my $num_samples = get_count_from_selected($dbh,$query_num_samples," ",$log_file);	
		#Disconnect db
		$dbh->disconnect();
		#################
			
		#If the number of samples is huge
		if ( $num_samples >= $cfg_hash->{'min_samples_4_joint'} ){		
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			
			my ($jobs_to_wait,$vcfs,$final_vcf);
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				($jobs_to_wait,$vcfs,$final_vcf) = run_GATK_GenotypeGVCFs_parallel($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,$step_needed,\@steps_array2,\@steps_array3);
		  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				($jobs_to_wait,$vcfs,$final_vcf) = run_GATK4_GenotypeGVCFs_parallel($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,$step_needed,\@steps_array2,\@steps_array3);
		  }
		  
			$prog_used = $cfg_hash->{'combvar_prog'};
			
			#Launch a job to combine the VCF obtained but get the command before
			print_and_log("$analysis_indication Launch a job to combine the VCF obtained:\n ",$log_file);
			my $combine_cmd = "";
			
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
			  $combine_cmd = run_GATK_CombineVariants($cfg_hash,$prog_used,$vcfs,$final_vcf,$log_file,"YES");
		  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
			  $combine_cmd = run_GATK4_CombineVariants($cfg_hash,$prog_used,$vcfs,$final_vcf,$log_file,"YES");
		  }
			
			#Run a job for CombineVariants
			my @commands = ($combine_cmd);
			#Path to the script to execute multiple jobs using an array of commands
			my $exec_job_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'exec_job_script'};

			my $combvar_job = execute_jobs($cfg_hash,$analysis_id,\@commands,$exec_job_script,$jobs_to_wait,$task,$task,$log_file);

			#Now alter the dependencies of the variant filtering job: it has to wait the combine variant job
			my $new_depend = "-W depend=afterok:$combvar_job ";
			my $job_name = "varfilt_a".$analysis_id;
			my $varfilt_jobid = get_jobid_form_jobname($job_name,$cfg_hash->{'qsub_username'});
			print_and_log("$analysis_indication Altering dependencies of "."varfilt_a".$analysis_id." ($varfilt_jobid):\n ",$log_file);
			alter_job($varfilt_jobid,$new_depend,$log_file);

			#Update the hash with dependencies
			my $deps_f = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'jc_dep_hash_f'};
			update_dependencies_hash($deps_f,$job_name,$varfilt_jobid,$new_depend,$log_file);
		}else{
			if ( $cfg_hash->{'gatk_ver'} < 4){
				#Here I give two steps array corresponding to input and output locations				
				run_GATK_GenotypeGVCFs($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,$step_needed,\@steps_array2,\@steps_array3);
		  }elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				run_GATK4_GenotypeGVCFs($cfg_hash,$prog_used,$sample_id,$analysis_id,$log_file,$task,$step,$step_needed,\@steps_array2,\@steps_array3);
		  }
		}
		
			#Save status in the db by searching for the specific sample name
		#print_and_log("$analysis_indication Updating analysis status for group $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
				$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
		$partTime = time;
	}

	##################################################
	## UPDATE STATUS
	##################################################
	
	#Update the group table with the status
	#The end of the genotype step means that the variant calling is finished
	#for the complete analysis
	my $fields = $cfg_hash->{'db_analysis_status'};
	my $values = $cfg_hash->{'analysis_varcall'};
	my $fields2 = $cfg_hash->{'db_analysis_id'};
	my $values2 = $analysis_id;
	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);
	##################################################
	
	
	#print_and_log("\n####################################################...\n",$log_file);
	print_and_log("$analysis_indication Genotyping task completed.\n",$log_file);
	#print_and_log("####################################################...\n",$log_file);	
}



=head2 run_var_filt_INDEL_SNP_steps

 Title   : run_var_filt_INDEL_SNP_steps
 Usage   : run_var_filt_INDEL_SNP_steps(  config_file => the config hash
								);

 Function: Executes step by step the programs needed for the variant filtering
					but then it selects the INDELs and SNPs put them in two different
					files and filters them with a different strong filter
					This because the HARD-FILTERING for SNPs could be too high for
					INDELs.
 
 Returns : the path of the last output file
 
=cut
sub run_var_filt_INDEL_SNP_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	#(here sample_id is the group name, since we bring this value from vargenius.pl)
	my $analysis_name = $sample_id;

	my $analysis_indication = "\n[an:$analysis_id ($analysis_name)](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed;
	my $extension;
	my @steps_arr = ();
	my $infold_step;
	#Parameters needed to construct the name
	my $table_4_name;
	my $id_4_name;
	my $sample_or_gr_id;
	
	#Getting the flag indicating if the joint genotype has to be run
	my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    

	#If it has been run, take the output form the genotype, else from the varcalling
	if ($dogenot){
		$step_needed = $cfg_hash->{'genot_step'}; 
		$extension = $cfg_hash->{'gvcf_ext'};
		@steps_arr = @steps_array3;
		$infold_step = $step_needed;
		#Set variable to construct the name
		$table_4_name = $cfg_hash->{'db_analyses_table'};
		$id_4_name = $cfg_hash->{'db_analysis_id'};
		$sample_or_gr_id = $analysis_id;
	}#If joint genotype has not been executed means that there is a single sample
	else{
		#Get the sample id from the group id
		my $sample_id_for_group = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    		
		#If there are more files associated to this sample
		if ( $perchrom == 1 ){
			#Checks that step_needed is executed for all the samples and updates that step
			#in the group table
			$step_needed  = $cfg_hash->{'catvar_step'};
			update_group_status($cfg_hash,$analysis_id,$cfg_hash->{'varcall_step'});
			
		}else{
			$step_needed  = $cfg_hash->{'varcall_step'};
		}
		$extension = $cfg_hash->{'vcf_ext'};
		@steps_arr = @steps_array2;
		$infold_step = $cfg_hash->{'varcall_step'};
		#Set variable to construct the name
		$table_4_name = $cfg_hash->{'db_sample_table'};
		$id_4_name = $cfg_hash->{'db_sample_id'};
		$sample_or_gr_id = $sample_id_for_group;
		
		#Checks that step_needed is executed for all the samples and updates that step
		#in the group table
		update_group_status($cfg_hash,$analysis_id,$step_needed);
	}


		
	#Verifying that step has been executed for this group
	my $db_step_needed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$step_needed,
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
		    	
	#Executing variantFiltration
	if ($cfg_hash->{'varfilt'} eq 'YES' ){
		if ( $db_step_needed == 1){
			my $prog_used = $cfg_hash->{'varfilt_prog'};	
			my $step = $cfg_hash->{'varfilt_step'};
			
			#input and output folder setting
			my $inFolder = $cfg_hash->{$analysis_id.'_'.$infold_step.'_out_f'};
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$step.'_out_f'};
			
			#Get infos for subject sample
			my $params_in;

			#print_and_log("$analysis_indication Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#
			getSampleConfiguration_locked(\$params_in,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
								$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$table_4_name,$id_4_name,$sample_or_gr_id);
			#Setting the input file name (could be gvcf or vcf)
			my $input_file = $inFolder."/".build_input_name_from_executed($params_in,$step,
										$analysis_name,\@steps_arr).".$extension ";
			
			#Set the unique output name (always VCF)
			my $params_out;
			#print_and_log("$analysis_indication Getting information for $analysis_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
			getSampleConfiguration_locked(\$params_out,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
								$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
			my $outFile = $outFolder."/".build_input_name_from_executed($params_out,$step,
							$analysis_name,\@steps_array3)."_".$step.".".$cfg_hash->{'vcf_ext'};
			
			#print_and_log( "$analysis_indication Configuration finished input: $input_file and output: $outFile \n",$log_file);#DEBUGCODE
			
			#Before to filter extract a file for INDELs and one for SNPs
			#At each step check if the input from the previous step is there and is not zero
			$prog_used = $cfg_hash->{'selvar_prog'};
			my $vartype = "INDEL";
			#Output indels called file
			my $indelsFile = $outFolder."/".build_input_name_from_executed($params_out,$step,
							$analysis_name,\@steps_array3)."_".$step."_$vartype.".$cfg_hash->{'vcf_ext'};
			print_and_log("$analysis_indication  RUNNING $prog_used for $vartype : ",$log_file);
			
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				#check input
				my $selvar_param_str = " -selectType $vartype ";	
				run_GATK_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$input_file,$indelsFile,$log_file);
			}elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				#check input
				my $selvar_param_str = " --select-type-to-include $vartype ";	
				run_GATK4_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$input_file,$indelsFile,$log_file);
			}	

			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used [$vartype] $partDuration seconds\n",$log_file);
			$partTime = time;
						
			#Now filter INDELs
			$prog_used = $cfg_hash->{'varfilt_prog'};
			my $filtIndOut = $outFolder."/".build_input_name_from_executed($params_out,$step,
							$analysis_name,\@steps_array3)."_".$step."_$vartype\_filt.".$cfg_hash->{'vcf_ext'};
			#check input
			print_and_log("$analysis_indication RUNNING $prog_used for $vartype : ",$log_file);
	
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				#Set the filters by searching in the user file
				my $param_str = " ";	
				#Cluster Window size
				if (defined $cfg_hash->{'VF_clusterWindowSize'}){$param_str .= " --clusterWindowSize ".$cfg_hash->{'VF_clusterWindowSize'}." ";}
				if (defined $cfg_hash->{'VF_clusterSize'}){$param_str .= " --clusterSize ".$cfg_hash->{'VF_clusterSize'}." ";}			
				run_GATK_VariantFiltration($cfg_hash,$prog_used,$indelsFile,$filtIndOut,$param_str,$cfg_hash->{'indel_filters'},$log_file);
			}elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				#Set the filters by searching in the user file
				my $param_str = " ";	
				#Cluster Window size
				if (defined $cfg_hash->{'VF_clusterWindowSize'}){$param_str .= " --cluster-window-size ".$cfg_hash->{'VF_clusterWindowSize'}." ";	}
				if (defined $cfg_hash->{'VF_clusterSize'}){	$param_str .= " --cluster-size ".$cfg_hash->{'VF_clusterSize'}." ";}	
				run_GATK4_VariantFiltration($cfg_hash,$prog_used,$sample_id,$analysis_id,$indelsFile,$filtIndOut,$param_str,$cfg_hash->{'indel_filters'},$log_file);
			}	
			
			
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used [$vartype] $partDuration seconds\n",$log_file);
			$partTime = time;
							
			#Select SNPs
			$vartype = "SNP";
			$prog_used = $cfg_hash->{'selvar_prog'};
			my $snpsFile = $outFolder."/".build_input_name_from_executed($params_out,$step,
							$analysis_name,\@steps_array3)."_".$step."_$vartype.".$cfg_hash->{'vcf_ext'};
			print_and_log("$analysis_indication RUNNING $prog_used for $vartype : ",$log_file);
			
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				#check input
				my $selvar_param_str = " -selectType $vartype ";
				run_GATK_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$input_file,$snpsFile,$log_file);			
			}elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				#check input
				my $selvar_param_str = " --select-type-to-include $vartype ";
				run_GATK4_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$input_file,$snpsFile,$log_file);
			}	
			
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used [$vartype] $partDuration seconds\n",$log_file);
			$partTime = time;
						
			#Now filter SNPs											
			$prog_used = $cfg_hash->{'varfilt_prog'};
			my $filtSnpOut = $outFolder."/".build_input_name_from_executed($params_out,$step,
							$analysis_name,\@steps_array3)."_".$step."_$vartype\_filt.".$cfg_hash->{'vcf_ext'};
			print_and_log("$analysis_indication RUNNING $prog_used for $vartype : ",$log_file);
			
			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				#Set the filters by searching in the user file
				my $param_str = " ";	
				#Cluster Window size
				if (defined $cfg_hash->{'VF_clusterWindowSize'}){$param_str .= " --clusterWindowSize ".$cfg_hash->{'VF_clusterWindowSize'}." ";}
				if (defined $cfg_hash->{'VF_clusterSize'}){$param_str .= " --clusterSize ".$cfg_hash->{'VF_clusterSize'}." ";}							
				run_GATK_VariantFiltration($cfg_hash,$prog_used,$snpsFile,$filtSnpOut,$param_str,$cfg_hash->{'snp_filters'},$log_file);	
			}elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				#Set the filters by searching in the user file
				my $param_str = " ";	
				#Cluster Window size
				if (defined $cfg_hash->{'VF_clusterWindowSize'}){$param_str .= " --cluster-window-size ".$cfg_hash->{'VF_clusterWindowSize'}." ";}
				if (defined $cfg_hash->{'VF_clusterSize'}){$param_str .= " --cluster-size ".$cfg_hash->{'VF_clusterSize'}." ";}							
				run_GATK4_VariantFiltration($cfg_hash,$prog_used,$sample_id,$analysis_id,$snpsFile,$filtSnpOut,$param_str,$cfg_hash->{'snp_filters'},$log_file);
			}	
			
					
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used [$vartype] $partDuration seconds\n",$log_file);
			$partTime = time;		
					
			#Combine VCFs into the same file
			$prog_used = $cfg_hash->{'combvar_prog'};
			#Check filtered SNP and INDELs outputs
			if (file_not_present($filtSnpOut) ){ log_and_exit( "Cannot proceed with $prog_used [$vartype]! Check: $filtSnpOut.\n",$log_file);}
			if (file_not_present($filtIndOut) ){ log_and_exit( "Cannot proceed with $prog_used [$vartype]! Check: $filtIndOut.\n",$log_file);}

			#Select GATK depending by the version
			if ( $cfg_hash->{'gatk_ver'} < 4){
				run_GATK_CombineVariants($cfg_hash,$prog_used,$filtSnpOut.",".$filtIndOut,$outFile,$log_file,"NO");
			}elsif ( $cfg_hash->{'gatk_ver'} >= 4){ 
				run_GATK_MergeVcfs($cfg_hash,$prog_used,$filtSnpOut.",".$filtIndOut,$outFile,$log_file);
			}
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;
						
			#Check that the output has been created
			if (file_not_present($outFile)  > 0){ log_and_exit( "ERROR: Something did not work with $prog_used! Check: $outFile.\n",$log_file);}
			
			#Save status in the db by searching for the specific sample name
			#print_and_log("$analysis_indication Updating analysis status for analysis $analysis_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step);

			
			#Delete all temporary files
			if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
				print_and_log("$analysis_indication Removing the following temporary files: $snpsFile\n$indelsFile\n$filtSnpOut\n$filtIndOut..\n",$log_file);	
				delete_file($snpsFile);
				delete_file($indelsFile);
				delete_file($filtSnpOut);
				delete_file($filtIndOut);
			}
		}else{
				print_and_log( "$analysis_indication WARNING: $step_needed has not been run. Execute this before varfilt step.. \n",$log_file);
		}
	}
	
}



=head2 run_variant_quality_score_recalibration_steps

 Title   : run_variant_quality_score_recalibration_steps
 Usage   : run_variant_quality_score_recalibration_steps(  config_file => the config hash
								);

 Function: Executes the variant_quality_score_recalibration steps given in
						GATK best practices
 
 Returns : nothing
 
=cut
sub run_variant_quality_score_recalibration_steps{

	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#(here sample_id is the group name, since we bring this value from vargenius.pl)
	my $group_name = $sample_id;

	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed;
	my $extension;
	my @steps_arr = ();
	#The last task is needed to identify where to pick the input file/s
	my $last_task;
	#Parameters needed to construct the name
	my $table_4_name;
	my $id_4_name;
	my $sample_or_gr_id;
	
	#Getting the flag indicating if the joint genotype has to be run
	my $dogenot = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_dogenot'},
						$cfg_hash->{'db_analysis_id'},$analysis_id);

	#Verifying that there has been an execution per chromosome
	my $perchrom = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_perchrom'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    

	#If it has been run, take the output form the genotype, otherwise from the varcalling
	if ($dogenot){
		$step_needed = $cfg_hash->{'genot_step'}; 
		$extension = $cfg_hash->{'gvcf_ext'};
		@steps_arr = @steps_array3;
		$last_task = $step_needed;
		#Set variable to construct the name
		$table_4_name = $cfg_hash->{'db_analyses_table'};
		$id_4_name = $cfg_hash->{'db_analysis_id'};
		$sample_or_gr_id = $analysis_id;
	}#If joint genotype has not been executed means that there is a single sample
	else{
		#Get the sample id from the group id
		my $sample_id_for_group = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
						    		
		#If there are more files associated to this sample catVariants is needed. We cannot proceed 
		#with multiple samples at this step
		if ( $perchrom == 1 ){
			#Checks that step_needed is executed for all the samples and updates that step
			#in the group table
			$step_needed  = $cfg_hash->{'catvar_step'};
			update_group_status($cfg_hash,$analysis_id,$cfg_hash->{'varcall_step'});
			
		}else{
			$step_needed  = $cfg_hash->{'varcall_step'};
		}
		$extension = $cfg_hash->{'vcf_ext'};
		@steps_arr = @steps_array2;
		$last_task = $cfg_hash->{'varcall_step'};
		#Set variable to construct the name
		$table_4_name = $cfg_hash->{'db_sample_table'};
		$id_4_name = $cfg_hash->{'db_sample_id'};
		$sample_or_gr_id = $sample_id_for_group;
		
		#Checks that step_needed is executed for all the samples and updates that step
		#in the group table
		#print_and_log("$analysis_indication Updating analysis status for analysis $analysis_id in ".$cfg_hash-> {'db_name'}." ($step_needed=1) ...\n",$log_file);	#DEBUGCODE
		update_group_status($cfg_hash,$analysis_id,$step_needed);
	}
	
	
	#Verifying that step has been executed for this group
	my $db_step_needed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$step_needed
						    ,$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $step = $cfg_hash->{'varrecal_step'};
		
	#input and output folder setting
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$last_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	
	#Get infos for subject sample
	my $params_in;

	#print_and_log("Getting information for $id_4_name: $sample_or_gr_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params_in,$cfg_hash->{'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$table_4_name,$id_4_name,$sample_or_gr_id);
	#Setting the input file name (could be gvcf or vcf)
	my $input_file = $inFolder."/".build_input_name_from_executed($params_in,$step,$group_name,\@steps_arr).".$extension ";
			
	#Set the unique output name (always VCF)
	my $params_out;
	#print_and_log("$analysis_indication Getting information for ".$cfg_hash->{'db_analysis_id'}.":$analysis_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params_out,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);

	#Before doing the analysis alert the user if there are less than X variants
	my $var_num = get_number_of_variants($input_file);
	if ( $var_num < 100000){
		print_and_log("$analysis_indication WARNING: VCF file $input_file has $var_num variants. (Consider that VQSR could stop with less than 100K)\n",$log_file);
	}
	
	#VARIANT RECALIBRATION
	#Call Variant Recalibration for both indels and snps
	if ( ($cfg_hash->{'varrecal'} eq 'YES')  ){	
		my $prog_used = $cfg_hash->{'varrecal_prog'};
		my $step = $cfg_hash->{'varrecal_step'};	

		#SNP Output files
		my $mode = "SNP";
		my $recal_file_snp = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.recal";
		my $tranches_file_snp = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.tranches";
		my $rscript_file_snp = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.r";				
		#Setting the first output file name (will be a vcf with only SNPs recalibrated and raw INDELs)
		my $recal_snps_raw_indels = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_$step\_$mode.".$cfg_hash->{'vcf_ext'};
		my $recal_snps_and_indels = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_$step.".$cfg_hash->{'vcf_ext'};
		#print_and_log( "Configuration finished input: $input_file , recal_file: $recal_file_snp,tranches_file= $tranches_file_snp, rscript_file= $rscript_file_snp\n",$log_file);#DEBUGCODE

		print_and_log( "$analysis_indication Doing VQSR for $mode...\n",$log_file);#DEBUGCODE
	
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
		run_GATK_VariantRecalibration($cfg_hash,$prog_used,$sample_id,$analysis_id,$mode,$input_file,$recal_file_snp,$tranches_file_snp,$rscript_file_snp,$log_file);

		$prog_used = $cfg_hash->{'apprecal_prog'};	
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
		run_GATK_ApplyRecalibration($cfg_hash,$prog_used,$sample_id,$analysis_id,$mode,$input_file,$recal_file_snp,$tranches_file_snp,$rscript_file_snp,$recal_snps_raw_indels,$log_file);


		#INDEL Execution (uses as input the first output obtained with SNPs)
		$mode = "INDEL";
		my $recal_file_indel = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.recal";
		my $tranches_file_indel = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.tranches";
		my $rscript_file_indel = $outFolder."/".build_input_name_from_executed($params_out,$step,$group_name,\@steps_array3)."_".$step."_$mode.r";				
		#print_and_log( "Configuration finished input: $input_file , recal_file: $recal_file_indel,tranches_file= $tranches_file_indel, rscript_file= $rscript_file_indel\n",$log_file);#DEBUGCODE
		
		print_and_log( "Doing now VQSR for $mode...\n",$log_file);#DEBUGCODE
		$prog_used = $cfg_hash->{'varrecal_prog'};		
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
		run_GATK_VariantRecalibration($cfg_hash,$prog_used,$sample_id,$analysis_id,$mode,$recal_snps_raw_indels,$recal_file_indel,$tranches_file_indel,$rscript_file_indel,$log_file);
		
		$prog_used = $cfg_hash->{'apprecal_prog'};
		print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
		run_GATK_ApplyRecalibration($cfg_hash,$prog_used,$sample_id,$analysis_id,$mode,$recal_snps_raw_indels,$recal_file_indel,$tranches_file_indel,$rscript_file_indel,$recal_snps_and_indels,$log_file);
		
		
		#Save status in the db 
		#print_and_log("$analysis_indication Updating analysis status for group $analysis_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
		update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
				$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
		$partDuration = time - $partTime;
		print_and_log( "Finished\t$prog_used\tsample\t$sample_id\t$partDuration\tseconds\n",$log_file);
		$partTime = time;
	}
	
}
	
	
=head2 run_genotype_refinement_steps

 Title   : run_genotype_refinement_steps
 Usage   : run_genotype_refinement_steps(  config_file => the config hash
								);

 Function: Executes the genotype refinement workflow as in GATK best practices.
					This works only when doing Joint Genotyping and VQSR hence just for more than
					30 exome samples.
 
 Returns : nothing
 
=cut
sub run_genotype_refinement_steps{
		my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed = $cfg_hash->{'varfilt_step'} ;
	my $step = $cfg_hash->{'genotref_step'};
	
	my $group_name = $sample_id;

	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Genotype refinement workflow execution
	if ($cfg_hash->{'genotref'} eq 'YES' ){
		
			#input folder setting
			my $inFolder = $cfg_hash->{$analysis_id.'_'.$step_needed.'_out_f'};
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$step_needed.'_out_f'};
			print_and_log("$analysis_indication STARTING GENOTYPE REFINEMENT... \n",$log_file);
			##Getting the pedigree file
			
			my $all_ped_file = $cfg_hash->{$analysis_id."_data_fold"}."/$group_name.".$cfg_hash->{'ped_ext'};	
			#print_and_log("$analysis_indication Obtaining a complete PED file : $all_ped_file \n",$log_file);
			#get_ped_file_for_enlarged_analysis($cfg_hash,$analysis_id,$all_ped_file,$log_file);
			
			my $params;
			#print_and_log("Getting information for ".$cfg_hash->{'db_analysis_id'}.":$analysis_id from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE#
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
						$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
			
			#Get the input file obtained after the variant quality score recalibration
			my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,$group_name,\@steps_array3).".".$cfg_hash->{'vcf_ext'};
			#Set the unique output name
			my $out_vcf = $outFolder."/".build_input_name_from_executed($params,$step,$group_name,\@steps_array3)."_$step.".$cfg_hash->{'vcf_ext'};
			
			my $prog_used = $cfg_hash->{'calc_gen_post_prog'};
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			my $vcf_posteriors = 	$input_file."_pos".".".$cfg_hash->{'vcf_ext'};
			my $param_str = " --supporting ".$cfg_hash->{'gatk_ref_f'}."/".$cfg_hash->{$cfg_hash->{'GRW_CGP_supporting'}}." ";		
			run_GATK_CalculateGenotypePosteriors($cfg_hash,$prog_used,$sample_id,$analysis_id,$input_file,$all_ped_file,$vcf_posteriors,$param_str,$log_file);
		
			$prog_used = $cfg_hash->{'varfilt_prog'};	
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			my $vcf_gen_filtered = 	$input_file."_gFiltered.".$cfg_hash->{'vcf_ext'};	
			#Set the filters by searching in the user file
			$param_str = " ";			
			run_GATK_VariantFiltration($cfg_hash,$prog_used,$vcf_posteriors,$vcf_gen_filtered,$param_str,$cfg_hash->{'GRW_genotype_filters'},$log_file);
		
			$prog_used = $cfg_hash->{'varann_prog'};	
			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			my $vcf_gfilt_denovo = 	$input_file."_gFiltered_denovo.".$cfg_hash->{'vcf_ext'};	
			run_GATK_VariantAnnotator($cfg_hash,$prog_used,$sample_id,$analysis_id,$input_file,$out_vcf,$all_ped_file,$log_file);
								
			#Save status in the db by searching for the specific sample name
			#print_and_log("Updating analysis status for group $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;				
	}	
}


=head2 run_phasing_steps

 Title   : run_phasing_steps
 Usage   : run_phasing_steps(  config_file => the config hash
								);

 Function: Executes the phasing based on the usage of the pedigree file
 
 Returns :
 
=cut
sub run_phasing_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	my $step_needed = $cfg_hash->{'varfilt_step'} ;
	my $step = $task;
	
		
	#(here sample_id is the group name, since we bring this value from vargenius.pl)
	my $analysis_name = $sample_id;

	my $analysis_indication = "\n[an:$analysis_id ($analysis_name)](".scalar(localtime)."): ";
	
	#Verifying that there has been an execution of multiple samples
	my $filter_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$step_needed,
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	
	#Phase by transmission execution
	if ($cfg_hash->{'phasing'} eq 'YES' ){
		#checking that filtering was performed
		if ( $filter_step == 1 ){
			my $prog_used = $cfg_hash->{'phasing_prog'};	
			my $step = $cfg_hash->{'phasing_step'};
			my $params;

			#input folder setting
			my $inFolder = $cfg_hash->{$analysis_id.'_'.$step_needed.'_out_f'};
			my $outFolder = $cfg_hash->{$analysis_id.'_'.$step_needed.'_out_f'};

			#Getting the pedigree file
			my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
										$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
										$cfg_hash->{'db_analysis_id'},$analysis_id);	
			if ( $ped_file eq 'none' ){ 
				log_and_exit("$analysis_indication WARNING: the analysis $analysis_id does not have a pedigree file. $prog_used will not be run\n",$log_file);
			}
			#The pedigree file is in the data folder of the group folder
			my $ped_f_path = $cfg_hash->{$analysis_id.'_data_fold'}."/$ped_file";
			
			#Get infos for subject sample
			#print_and_log("$analysis_indication Getting information for $group_name from database ".$cfg_hash-> {'db_name'}."...\n",$log_file);#DEBUGCODE
			getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
								$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);

			my $input_file = $inFolder."/".build_input_name_from_executed($params,$step,$analysis_name,\@steps_array3).".".$cfg_hash->{'vcf_ext'};
			#Set the unique output name
			my $out_vcf = $outFolder."/".build_input_name_from_executed($params,$step,$analysis_name,\@steps_array3)."_$step.".$cfg_hash->{'vcf_ext'};
			

			print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);						
			run_GATK_PhaseByTransmission($cfg_hash,$prog_used,$analysis_name,$analysis_id,$outFolder,$input_file,$out_vcf,$ped_f_path,$log_file);

			#Use the mendelian violations output to get the number of MV and give it in the html output
			my $mendel_viols_out = extract_name($out_vcf,"noext").".".$cfg_hash->{'mendel_viols_ext'};
			#Save the path for the future :-) (I will print the number of mendelian violations into the output
			$cfg_hash->{'mendel_viols_out'} = $mendel_viols_out;
			if ( -e $mendel_viols_out){
				########################################
				##########FILLING OUTLIST FILE##########
				########################################	

				#Put the path to the output in the outlist file
				if ( check_presence($mendel_viols_out) ){
					print_and_log( "$analysis_indication Appending $mendel_viols_out path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
					my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_mendviol'}."\t".$mendel_viols_out;
					append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
				}else{
					print_and_log( "$analysis_indication WARNING: $mendel_viols_out is not present. Please check what is happen. \n",$log_file);	
				}		
				########################################
				##########FILLING OUTSTATS FILE##########
				########################################					
						
				#MENDELIAN VIOLATIONS Add number of mendelian violations into the output
				#counting the number of lines of the output from PhaseByTransmission and removing lines with at least one undefined genotype (./.)
				if ( !-z $mendel_viols_out){
					my $command = 'grep -v '."'".'\./\.'."' ".$mendel_viols_out.' | wc -l';
					print_and_log( "$analysis_indication Executing: $command..\n",$log_file);
					my $num_mend_viol =  `$command`;
					if (defined $num_mend_viol){
						chomp($num_mend_viol);
						my $line = "Mendelian Violations (lines with at least one undefined genotype (./.) are removed) : $num_mend_viol\t[MendelianViolations]";
						print_and_log( "$analysis_indication Appending number of mendelian violation to ".$cfg_hash->{$analysis_id."_outstats_file"}." file... \n",$log_file);
						overwrite_str_in_file_if_exists($cfg_hash->{$analysis_id."_outstats_file"},$line);
						#Update statistics table
						update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
													$cfg_hash->{'db_analysis_statistics'},$cfg_hash->{'db_analstats_mendelviol'},$num_mend_viol,$cfg_hash->{'db_analysis_id'},$analysis_id);
					}else{
						print_and_log( "ERROR: with command: $command... \n",$log_file);
					}
				}
			}
			
			
			#Save status in the db by searching for the specific sample name
			#print_and_log("$analysis_indication Updating analysis status for analysis $analysis_name in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
			$partTime = time;		
		
		
		}else{
			print_and_log( "$analysis_indication WARNING: Cannot run ".$cfg_hash->{'phasing_step'}." step because $step_needed has not been run...\n",$log_file);
		}
		
	}
}



=head2 run_reannotate_variants

 Title   : run_reannotate_variants
 Usage   : run_reannotate_variants(  config_file => the config hash
								);

 Function: This module is the reannotation module. It fetches all the variants from the database
					reannotates them and imports the new annotation into the database.
 
 Returns : 
 
=cut
sub run_reannotate_variants{
	my $config_file = shift;

	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
		
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	my $prog_used = $cfg_hash->{'table_annov_prog'};
	
	#Get a VCF with all the variants of the database
	my $vcf_to_ann = $cfg_hash->{'main_data_folder'}."/all_variants.vcf";
	print_and_log("Generating a VCF file from the database with all variants: $vcf_to_ann ...\n",$log_file);
	####all_variants_2_VCF($cfg_hash,$vcf_to_ann ,$log_file);
	
	my $in_annov = $cfg_hash->{'main_data_folder'}."/all_variants.annov";
	#Converts the VCF to a format readable by Annovar
	print_and_log("Convert VCF format to ANNOVAR format in $in_annov ...\n",$log_file);									
	####run_ANNOVAR_convert2annovar($cfg_hash,$prog_used ,$log_file,$vcf_to_ann,$in_annov);

	#Annotate missing variants
	my $out_file_suff = $cfg_hash->{'main_data_folder'}."/all_variants";	
	#Annotate all with all the databases and settings from the user configuration
	print_and_log("Annotating $in_annov...\n",$log_file);
	
	#Reset the configuration variables according with the annotation
	$cfg_hash->{'annov_g_dbs'} = $cfg_hash->{'annov_g_2db'};
	$cfg_hash->{'annov_f_dbs'} = $cfg_hash->{'annov_f_2db'};
	$cfg_hash->{'annov_r_dbs'} = $cfg_hash->{'annov_r_2db'};
			
	####run_ANNOVAR_table_annovar($cfg_hash,$prog_used,$log_file,$in_annov,$out_file_suff);
	print_and_log("Completed the annotation of $in_annov. Output file is: $out_file_suff...\n",$log_file);	
	
	#Imports into the db
	my $ann_out =  $cfg_hash->{'main_data_folder'}."/all_variants.".$cfg_hash->{'ann_build_ver'}."_".$cfg_hash->{'annovar_out_file'};
	if ( -e $ann_out) {
		my $prog_used = "VarGenius";	

		print_and_log("[GLOBAL] (".scalar(localtime)."): RUNNING $prog_used variants reannotation..\n",$log_file);	

		import_annovarout_2_db_fast($cfg_hash,$ann_out,$prog_used,$analysis_id,$log_file);
		
		#Save status in the db by searching for the specific sample name
		#print_and_log("Updating analysis status for group $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
		#update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
				#$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
				#$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
		#$partDuration = time - $partTime;
		#print_and_log( "Finished\t$prog_used\tsample\t$sample_id\t$partDuration\tseconds\n",$log_file);
		#$partTime = time;	
	}
	
	#Here I print a log with a specific format for GLOBAL tasks
	$partDuration = time - $partTime;
	print_and_log( "Finished [GLOBAL]\t$prog_used\t$partDuration\tseconds\n",$log_file);
	$partTime = time;
}
	



=head2 run_recompute_variants_frequencies

 Title   : run_recompute_variants_frequencies
 Usage   : run_recompute_variants_frequencies(  config_file => the config hash
								);

 Function: This module permits to recompute all the frequencies of the variants 
						into the database.
						
						To compute the allelic frequency of a given variant we first compute its:
							- number of times the variant is found heterozygous (V_het)
							- number of times the variant is found homozygous (V_hom)
							- number of samples for which for the variant has been computed a genotype (Tot_gen)
							
						Then the formula to obain the allelic frequency of the variant V is
						
						V_all_f = V_het + (V_hom * 2) / (Tot_gen * 2)
						
						Where we multiply V_hom by two because the homozygous variant has been found in both alleles
						and Tot_gen because the samples have two alleles. We use only those variant for which GATK
						is able to get a genotype.
						
 
 Returns : updates the allelefreq field of the variants table
 
=cut
sub run_recompute_variants_frequencies{
	my $config_file = shift;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
		
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	my $step = $cfg_hash->{'update_freqs_step'};

	#Get all sequencing type and for each type do the update
	print_and_log( "Get all sequencing types\n",$log_file);
	my $query = "SELECT DISTINCT ".$cfg_hash->{'db_sequencingtype'}." FROM  ".$cfg_hash->{'db_analyses_table'}.";";	
	print_and_log( "Executing: $query\n",$log_file);
	my $seqtypes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sequencingtype'});
				
	#For each sequencing type execute an update process of the frequencies
	if (scalar(keys %{$seqtypes}) > 1){
		foreach my $seqtype (keys %{$seqtypes}){
			print_and_log("Starting the update process of variants frequencies for $seqtype\n",$log_file);
			
			#################DB
			#Connect to database
			my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
			# PERL DBI CONNECT
			my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
			
				
			#Get all the variant ids
			print_and_log( "Get all the variant ids\n",$log_file);
			my $query = "SELECT ".$cfg_hash->{'db_var_id'}." FROM ".$cfg_hash->{'db_variants_table'}.";";							
			print_and_log("Executing: $query\n",$log_file);
			#fetch_all_rows
			my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
										$cfg_hash->{'db_pass'},$query);
	

			#Table of all variants in heterozygosity (must be in use for frequency calculation and for this type  of sequencing)
			my $het_vars_query = "SELECT ".$cfg_hash->{'db_var_ids'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype')".
											 " AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values('1|0'),('1/0'),('0|1'),('0/1'))";

			#Table of all variants in homozygosity	(must be in use for frequency calculation and for this type  of sequencing)							 
			my $hom_vars_query = "SELECT ".$cfg_hash->{'db_var_ids'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype')".
											 " AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values ('1|1'),('1/1'))";
																					
			#Get the count of distinct samples where the genotype has been predicted (must be in use for frequency calculation and for this type of sequencing)
			my $distinct_gen_samples_infreq =  "SELECT DISTINCT ".$cfg_hash->{'db_sample_id'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}.
											 " WHERE ".$cfg_hash->{'db_analysis_id'}." IN ( SELECT ".$cfg_hash->{'db_analysis_id'}.
											 " FROM ".$cfg_hash->{'db_analyses_table'}." WHERE ".$cfg_hash->{'db_analysis_infreq'}."=1 AND ".$cfg_hash->{'db_sequencingtype'}."='$seqtype')".
											 " AND  ".$cfg_hash->{'db_genotype_gt'}." LIKE ANY ( values ('1|1'),('1/1'),('0|0'),('0/0'),('1|0'),('1/0'),('0|1'),('0/1') )";	 
			my $tot_infreq_gen_sam = get_count_from_selected($dbh,$distinct_gen_samples_infreq," ",$log_file);
			
			#If distinct genotyped samples are present...
			if ( defined $tot_infreq_gen_sam){
				#Get two lists corresponding to the variants in heterozygosity and those in homozygosity
				my $het_vars_lists = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
											$cfg_hash->{'db_pass'},$het_vars_query);
				my $hom_vars_lists = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, 
											$cfg_hash->{'db_pass'},$hom_vars_query);
				
				#The variants will result as lists separated by comma. I want here a full list with all the variants
				#in heterozygosity and homozygsity so that later I use the grep on it
				my $het_vars_hash;
				foreach my $het_vars_list (@$het_vars_lists){
					my @fields = @$het_vars_list;
					foreach (split(",",$fields[0])){
						$het_vars_hash->{$_}++;
					}
					#print_and_log("Splitting and saving ".$fields[0]."\n",$log_file);#DEBUGCODE
				}
				my $hom_vars_hash;
				foreach my $hom_vars_list (@$hom_vars_lists){
					my @fields = @$hom_vars_list;
					foreach (split(",",$fields[0])){
						$hom_vars_hash->{$_}++;
					}
				}
				print_and_log("Obtained ".scalar(keys(%$het_vars_hash))." het_vars and ".scalar(keys(%$hom_vars_hash))." hom vars\n",$log_file);	
				
				#Print some statistics									
				print_and_log("het_vars_query: $het_vars_query\n",$log_file); 
				print_and_log("hom_vars_query: $hom_vars_query\n",$log_file); 
				print_and_log("distinct_gen_samples_infreq: $distinct_gen_samples_infreq\n",$log_file); 
				print_and_log("Total Genotyped $seqtype samples: $tot_infreq_gen_sam\n",$log_file);
				print_and_log("Out of ".scalar(@$res)." $seqtype variants, frequencies updated:\n",$log_file);
				
				#Using some variables to print at any 10000 variants the status of the update
				my $countThresh = 10000;									
				my $var_count = 0;
				
				#For each variant get the count of occurrences from all the samples in the database
				#avoiding those whose field infreq=0. And update the variant frequency table
				foreach my $var_field (@$res) {
					$var_count++;
					#Prints the number of variants evaluated
					if (($var_count >= 0) and ($var_count % $countThresh == 0) ){
						print_and_log(" $var_count - ",$log_file);
					}
					my @arr = @$var_field;
					#Get the database id of the variant
					my $varid = $arr[0];
					
					#print_and_log("\nComputing the frequency for $varid\n",$log_file);
					
					#Compute the number of times the variants is found hetorozygous
					my $het_var_count = "";
					if (defined $het_vars_hash->{$varid}){
						$het_var_count = $het_vars_hash->{$varid};
					}else{$het_var_count = 0;}
					#Compute the number of times the variants is found homozygous
					my $hom_var_count = "";
					if (defined  $hom_vars_hash->{$varid}){
						$hom_var_count = $hom_vars_hash->{$varid};
					}else{$hom_var_count = 0;}

																				
					#Divide heterozygosity and homozygosity by the total number of samples
					#Check if the denominator is greater than zero
					my $sc_var_freq = -1;
					
					if ( $tot_infreq_gen_sam > 0 ){
						$sc_var_freq = ($het_var_count + ($hom_var_count*2))/($tot_infreq_gen_sam*2);
					}
					##Update table
					#Truncate too much long values coming out from the division				
					my $var_freq = 	sprintf("%.10f",$sc_var_freq);			
					#Build the freq factors field that contains all the factors used in the division
					my $freq_factors = "$het_var_count".$cfg_hash->{'mult_ann_sep'}."$hom_var_count".$cfg_hash->{'mult_ann_sep'}."$tot_infreq_gen_sam";
					#Build fields and values to put into the db
					my $fields = $cfg_hash->{"db_$seqtype\_allele_freq"}.";".$cfg_hash->{"db_$seqtype\_freq_factors"};
					my $values = $var_freq.";'$freq_factors'";
					#Query to update
					update_table_woconn($dbh,$cfg_hash->{'db_variants_table'},$fields,$values,$cfg_hash->{'db_var_id'},$varid);
					
					#print_and_log("Frequency of $varid: ($het_var_count + ($hom_var_count*2))/($tot_infreq_gen_sam*2)\n",$log_file);#DEBUGCODE
				}
				#Disconnect db
				$dbh->disconnect(); 		
				#Update the variants frequency update history table with the date
				my $fields = $cfg_hash->{'db_sequencingtype'}.",".$cfg_hash->{'db_var_freq_update'};
				my ($d,$m,$y) = (localtime)[3,4,5];
				my $year = $y+1900; # In PERL year
				$m = $m + 1; #Month in perl starts from 0 (0:january,11:december)
				my $values = "'$seqtype','$year\_$m\_$d'";
				insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_var_freq_update_hist_table'},$fields,$values);
									
				$partDuration = time - $partTime;
				print_and_log("\nFinished [GLOBAL]\trecomputing variants frequencies\t$seqtype\t$partDuration\tseconds\n",$log_file);
				$partTime = time;					
			}
		}
	}	
}


=head2 run_allgenes_annotation_steps

 Title   : run_allgenes_annotation_steps
 Usage   : run_allgenes_annotation_steps(  config_file => the config hash
								);

 Function: This module launches a subroutine that generates a table with
					the annotation of all genes
 
 Returns :
 
=cut
sub run_allgenes_annotation_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	my $prog_used = "VarGenius";	
	print_and_log("[an:$analysis_id] (".scalar(localtime)."):  RUNNING $prog_used all_genes_annotation..\n",$log_file);	
	print_and_log( "Executing the  annotation of all genes...\n",$log_file);
  my $ann_file = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'db_name'}."_all_genes.". $cfg_hash->{'gene_annotation_ext'};		
	get_all_genes_annotation($cfg_hash,$ann_file,$log_file);	

	#Here I print a log with a specific format for GLOBAL tasks
	$partDuration = time - $partTime;
	print_and_log( "Finished [GLOBAL]\t$prog_used\t$partDuration\tseconds\n",$log_file);
	$partTime = time;	
}



=head2 run_genedb_creation_steps

 Title   : run_genedb_creation_steps
 Usage   : run_genedb_creation_steps(  config_file => the config hash
								);

 Function: This module launches the creation of parts of the database needed for further
					analyses
 
 Returns :
 
=cut
sub run_genedb_creation_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	my $gene_db_fold = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'gene_db_fold'};

	#An hash that will contain the db_id of the HPO terms
	my $hpoTableHash;
	#HPO
	my $create_hpo_table = 0;
	
	print_and_log( "\n--Filling the tables for genes and phenotypes: \n",$log_file);
	
	if ( $cfg_hash->{'genes2hpo_parsing'} eq 'YES'){
		$create_hpo_table = 1;
	}
	
	if ( $cfg_hash->{'hpo_parsing'} eq 'YES'){
		$create_hpo_table = 1;
		my $prog_used = 'hpo_parsing';
		my $hpo_table_link = $cfg_hash->{'hpo_table_link'};
		my $hpo_table =  extract_name($hpo_table_link,0);#Get the name with extension
		print_and_log( "\n----Starting hpo_parsing ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$hpo_table) ){
			print_and_log(  "\n Downloading File: ".$hpo_table_link."\n",$log_file);
			download_file($hpo_table_link,$gene_db_fold);
		}else{
			print_and_log(  "File $hpo_table_link already downloaded...\n",$log_file);
		}
		#Fill table with HPO phenotypes and its associated Genes
		createHPOTable($cfg_hash,$gene_db_fold."/".$hpo_table,\$hpoTableHash,$gene_db_fold,$create_hpo_table,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}
		
	
	
	#HPO
	if ( $cfg_hash->{'genes2hpo_parsing'} eq 'YES'){
		my $prog_used = 'genes2hpo_parsing';
		my $genes_2_hpo_table_link = $cfg_hash->{'genes_2_hpo_table_link'};
		my $genes_2_hpo_table =  extract_name($genes_2_hpo_table_link,0);#Get the name with extension
		print_and_log( "\n----Starting genes to phenotypes parsing ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$genes_2_hpo_table) ){
			print_and_log(  "\n Downloading File: ".$genes_2_hpo_table_link."\n",$log_file);
			download_file($genes_2_hpo_table_link,$gene_db_fold);
		}else{
			print_and_log( "File $genes_2_hpo_table_link already downloaded...\n",$log_file);
		}
		#Fill table with HPO phenotypes and its associated Genes
		parse_HPO_genes_to_phenotypes($cfg_hash,\$hpoTableHash,$gene_db_fold."/".$genes_2_hpo_table,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}

	
	#GDI
	if ( $cfg_hash->{'gdi_parsing'} eq 'YES'){
		my $prog_used = 'gdi_parsing';
		my $gdi_table_link = $cfg_hash->{'gdi_table_link'};
		my $gdi_table = extract_name($gdi_table_link,0);#Get the name with extension
		print_and_log( "\n----Starting $prog_used ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$gdi_table) ){
			print_and_log( "\n Downloading File: $gdi_table_link \n",$log_file);
			download_file($gdi_table_link,$gene_db_fold);
		}else{
			print_and_log( "File $gdi_table_link already downloaded...\n",$log_file);
		}
		#Fill table with HPO phenotypes and its associated Genes
		parse_GDI_scores($cfg_hash,$gene_db_fold."/".$gdi_table,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}

	#RVIS
	if ( $cfg_hash->{'rvis_parsing'} eq 'YES'){
		my $prog_used = 'rvis_parsing';
		my $rvis_table_link = $cfg_hash->{'rvis_table_link'};
		my $rvis_table = extract_name($rvis_table_link,0);#Get the name with extension
		print_and_log( "\n----Starting $prog_used ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$rvis_table) ){
			#print "\n Downloading File: ".$rvis_table_link."\n";
			download_file($rvis_table_link,$gene_db_fold);
		}else{
			print_and_log( "File $rvis_table_link already downloaded...\n",$log_file);
		}
		#Fill table with HPO phenotypes and its associated Genes
		parse_RVIS_scores($cfg_hash,$gene_db_fold."/".$rvis_table,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}

	#RefSeq 2 genes parsing
	if ( $cfg_hash->{'refseq2genes_parsing'} eq 'YES'){
		my $prog_used = 'refseq2genes_parsing';
		my $refseq2genes_table_link = $cfg_hash->{'refseq2genes_link'};
		my $refseq2genes_compres = $gene_db_fold."/".extract_name($refseq2genes_table_link,0);#Get the name with extension
		my $refseq2genes_uncompres = $gene_db_fold."/".extract_name($refseq2genes_table_link,1);#Get the name without
		print_and_log( "\n----Starting $prog_used ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$refseq2genes_uncompres) ){
			#print "\n Downloading File: ".$rvis_table_link."\n";
			dl_and_extract($refseq2genes_table_link,$refseq2genes_compres,$refseq2genes_uncompres,$gene_db_fold,$log_file);
		}else{
			print_and_log(  "File $refseq2genes_compres already downloaded...\n",$log_file);
		}
		#Map gene symbols to RefSeq ids
		##Get RefSeq ids and Entrez id associated to gene symbols
		parse_genes2RefSeq($cfg_hash,$refseq2genes_uncompres,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}
	
	
	#OMIM
	if ( $cfg_hash->{'omim_parsing'} eq 'YES'){
		my $prog_used = 'omim_parsing';
		my $omim_table_link = $cfg_hash->{'omim_table_link'};
		my $omim_table = extract_name($omim_table_link,0);#Get the name with extension
		print_and_log( "\n----Starting omim_parsing ..\n",$log_file);
		#Download file if it doesn't exists
		if (not (-e $gene_db_fold."/".$omim_table) ){
			#print "\n Downloading File: ".$omim_table_link."\n";
			download_file($omim_table_link,$gene_db_fold);
		}else{
			print_and_log(  "File $omim_table_link already downloaded...\n",$log_file);
		}
		#Fill table with OMIM ids associated to Genes
		parse_OMIM_2_Genes($cfg_hash,$gene_db_fold."/".$omim_table,$log_file);
		$partDuration = time - $partTime;
		print_and_log( "Finished[GLOBAL]\t$prog_used in\t".$partDuration." seconds\n",$log_file);
		$partTime = time;
	}
	
	
		
}


=head2 get_readf_paths_using_step

 Title   : get_readf_paths_using_step
 Usage   : get_readf_paths_using_step(  config_file => the config hash
								);

 Function: Get the path of the files using the step
					Takes in input the config hash, the sample id to which the read files
					are associated, the step with which it has to construct the file name
					the steps array to use and the folder to put in front of the 
					file name.
					
 Returns : an array with paths
 
=cut
sub get_readf_paths_using_step{
	my $cfg_hash = shift;
	my $sample_id = shift;
	my $log_file = shift;
	my $step = shift;
	my $inFolder = shift;
	my $steps_array = shift;
	my $extension = shift;
	
	my $main_name = $cfg_hash->{'db_readf_name'};
	#Here we get the read file ids associated with this sample. They will be used for the input file names
	my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
					$cfg_hash->{'db_sample_id'},$sample_id);
							
	#Here get the samples names to be merged
	my @paths  = ();
	#I defined params here because I need it in the construction of output file name
	#there we need just one among all
	my $params;
	foreach my $readf_id (
		sort keys %{$distinct_readf_ids}
		)
	{
		#print_and_log("Sample id is: $readf_id..\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
							$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
		 #Build the input file name using as stop point the base recalibration step
		 my $path = $inFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array)."_".$step.".".$extension;
		push(@paths,$path);
	}
	return @paths;
}

=head2 get_readf_paths_and_names_using_step

 Title   : get_readf_paths_and_names_using_step
 Usage   : get_readf_paths_and_names_using_step(  config_file => the config hash
								);

 Function: Get the path of the files and their readfile names using the step
					Takes in input the config hash, the sample id to which the read files
					are associated, the step with which it has to construct the file name
					the steps array to use and the folder to put in front of the 
					file name.
					
 Returns : fills the two arrays given in input paths and names
 
=cut
sub get_readf_paths_and_names_using_step{
	my $cfg_hash = shift;
	my $sample_id = shift;
	my $log_file = shift;
	my $step = shift;
	my $inFolder = shift;
	my $steps_array = shift;
	my $extension = shift;
	my $paths  = shift;
	my $names  = shift;
	
	my $main_name = $cfg_hash->{'db_readf_name'};
	#Here we get the read file ids associated with this sample. They will be used for the input file names
	my $distinct_readf_ids = select_distinct_samples($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_readf_id'},
					$cfg_hash->{'db_sample_id'},$sample_id);
							

	#I defined params here because I need it in the construction of output file name
	#there we need just one among all
	my $params;
	foreach my $readf_id (
		sort keys %{$distinct_readf_ids}
		)
	{
		print_and_log("Sample id is: $readf_id..\n",$log_file);#DEBUGCODE
		getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},$cfg_hash-> {'db_user'},
							$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_readf_table'},$cfg_hash-> {'db_readf_id'},$readf_id);
		 #Build the input file name using as stop point the base recalibration step
		 my $path = $inFolder."/".build_input_name_from_executed($params,$step,
						$params->{$main_name},$steps_array)."_".$step.".".$extension;
		push(@$paths,$path);
		push(@$paths,$params->{$main_name});
	}
}




=head2 run_coverage_stepsOLD

 Title   : run_coverage_steps
 Usage   : run_coverage_steps(  config_file => the config hash
								);

 Function: The coverage analysis implies different steps. In particular this function
					is executed per-sample.
					
					Two main files are used to get the statistics:	
						- the BAM file just after the BWA alignment
						- the BAM file where the duplicates have been removed.
					If the DuplicateRemoval was not executed (e.g. targeted sequencing), then
					the second BAM will not be present.
					
					STEP 1: Construct the name of the last BAM (bam_input)
						IF sample is divided into multiple read files	:					
							In case there are multiple read files for the sample, these must be merged into a single
							and the duplicate removal must be performed again.
							The input BAM files per-lane can be removed after this step
							(if the variant calling is already finished)
							Before to remove I check if this task is the last one that is using the lane files
							(The variant calling step, for example, uses them.)
						IF the sample bam file is unique the name is reconstucted using the database
					
					STEP2: get non covered regions.
						Using $bam_input a procedure to get the non-covered regions is ran
						
					
					STEP3: Flagstat Statistics
						Aligned:
						Flagstat is executed first with the bam file just after the alignment. 
						In this case we do not have a bam file but the SAM file after the alignment
						So we merge these files into a single BAM
					
						With duplicate Removed:
						Flagstat is also executed on the file with duplicate removed that is the
						$bed_merged_rmdup previously constructed using $cfg_hash->{'mark_rem_dup_step'}
						parameter. If it exists, of course.
						
						Intersected with the target file:
						Last alignment statistics are performed on the bam file intersected with the target
						$bam_input is used for the intersection.			
 
 Returns :
 
=cut
sub run_coverage_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	
	my $analysis_indication = "\n[an:$analysis_id sam: $sample_id] (".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("\n$analysis_indication Starting alignment statistics with Flagstat and non-covered region analysis [per-sample]..\n",$log_file);
	
	
	my $step_needed = $cfg_hash->{'sort_idx_step'};
	#Checks if sort and indexing has been performed for each read file
	update_sample_status($cfg_hash,$sample_id,$cfg_hash->{'sort_idx_step'});
	
	my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
	my $out_fold_suff = $cfg_hash->{'stats_step'};#The folder that will be used for the output is the one for statistics

	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	
	#Get also the sample_name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	#Check if markDuplicate after the alignment was executed
	#Update the sample status
	update_sample_status($cfg_hash,$sample_id,$cfg_hash->{'mark_rem_dup_step'});
	my $mark_rem_dup_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mark_rem_dup_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
	#Checks if sort and indexing has been performed for each read file
	my $mergebam_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mergebam_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);			
	my $mergealn_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mergealn_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    						    
	###GET THE NAMES OF THE FILES
	my $bam_input = "";
	my $bam_aligned = "";

	#Out folder
	my $refinetask = $cfg_hash->{'refine_step'};
	my $refineOutFolder = $cfg_hash->{$analysis_id.'_'.$refinetask.'_out_f'};
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash-> {'db_name'},$cfg_hash-> {'db_dsn'},
					$cfg_hash-> {'db_user'},$cfg_hash-> {'db_pass'},$cfg_hash-> {'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

	#get bam recalibrated
	my $bam_recalibrated = $refineOutFolder."/".get_output_name_from_executed($params,$cfg_hash->{'print_reads_step'},
					$params->{$cfg_hash->{'db_sample_name'}},\@steps_array2).".".$cfg_hash->{'bam_ext'};
					
					 
	#Names of the flagstat outputs
	my $flagstat_aln = "";
	my $flagstat_recalibrated = "";
	
	#########
	####STEP 1: check if there are multiple read files. All the steps will use a unique bam file for the sample.		    
	#########
	
	#If there are more read_files per sample than the merge is needed merge them						        						    
	if ( $multiple == 1){
		print_and_log( "\n$analysis_indication Sample $sample_name is sub-divided in multiple read files. A merge is needed..\n",$log_file);
		
		my $params;
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		#The bam file merged with duplicates removed and sorted
		$bam_input = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'mergebam_step'}."_".$cfg_hash->{'sort_idx_step'}.".".$cfg_hash->{'bam_ext'};

		#The merged file from only aligned
		$bam_aligned = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'bam_ext'};
								
								
		print_and_log("\n$analysis_indication Input BAM (eventually duplicates removed): $bam_input - BAM aligned only: $bam_aligned..\n",$log_file);#DEBUGCODE
		
		#If readfiles are multiple and one among the flagstat and the non-covered region step
		#must be executed, they must be merged
		if ($cfg_hash->{'flagstat_mrdup'} eq 'YES' 
					or $cfg_hash->{'noncovered_regions'} eq 'YES' 
						or $cfg_hash->{'consensus_vc_approach'} ne 'NO' ){
			#If the merged BAM file has not been already producedor size is zero			
			 if ( !-e $bam_input or -z  $bam_input){
					if ( $cfg_hash->{'mergebam'} eq 'YES' ){
						#Here get the samples to be merged
						my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$cfg_hash->{'mark_rem_dup_step'},$inFolder,\@steps_array,$cfg_hash->{'bam_ext'});		

						print_and_log( "\n$analysis_indication Starting MergeBam to merge multiple samples..\n",$log_file);

						my $samples_to_merge  = " I=".join(" I=",@paths)." ";
						#Create the name for the bam file of the merge from the lane files
						my $bam_input_merged = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'mergebam_step'}.".".$cfg_hash->{'bam_ext'};						
						my $prog_used = $cfg_hash->{'mergebam_prog'};
						print_and_log( "RUNNING $prog_used : ",$log_file);#DEBUGCODE
						run_PICARD_MergeSamFiles_gen($cfg_hash,$prog_used,$samples_to_merge,$bam_input_merged,$log_file);
						
						my $bam_input_merge_sort = $inFolder."/".extract_name($bam_input_merged,1)."_".$cfg_hash->{'sort_idx_step'}.".".$cfg_hash->{'bam_ext'};
						print_and_log("$analysis_indication RUNNING Sort : ",$log_file);
						run_SAMTOOLS_Sort_gen($cfg_hash,$bam_input_merged,$bam_input_merge_sort,$log_file);
						print_and_log("$analysis_indication RUNNING Index : ",$log_file);
						run_SAMTOOLS_Index_gen($cfg_hash,$bam_input_merge_sort,$log_file);						

						#Save status in the db by searching for the specific sample name
						#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
						update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
							$cfg_hash->{'db_sample_id'},$sample_id,$cfg_hash->{'mergebam_step'},$cfg_hash->{'db_analysis_id'},$analysis_id);
						
						#Remove the mergebam file before the sortidx
						if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
							delete_file($bam_input_merged);	
						}
						
						##################################	
						#READ FILES REMOVAL
						#################################
						#Before to remove I check if this task is the last one that is using the lane files
						#The variant calling step, for example, uses them. 
						if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
							print_and_log( "\n$analysis_indication Removing read files\n",$log_file);
							#Check if this is the last task using the BAM files 
							my $bamfilesused = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},
												$cfg_hash->{'db_sample_id'},$sample_id);	
							if ( $bamfilesused == ($cfg_hash->{'tasks_using_bam'} - 1) )	{
								foreach my $bam_path (@paths){
									print_and_log("$analysis_indication Removing the alignment $bam_path file..\n",$log_file);	
									delete_file($bam_path);		
									delete_file($bam_path.".".$cfg_hash->{'bai_ext'});		
								}
								$bamfilesused++;
								update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
												$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},$bamfilesused,$cfg_hash->{'db_sample_id'},$sample_id);								
							}else{
								print_and_log( "\n$analysis_indication Read files will not be removed. Still ".($cfg_hash->{'tasks_using_bam'} - 1) - $bamfilesused." steps have to use them..\n",$log_file);
							}
						}							
						#Print time needed
						$partDuration = time - $partTime;
						print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
						$partTime = time;					
					}
			 }		
		}
		
		
		#Now check if flagstat for aligned read file must be executed and in case do the merge
		if ($cfg_hash->{'flagstat_sorted'} eq 'YES' and $cfg_hash->{'mergebam'} eq 'YES'){
			 #If the merged aligned BAM file has not been already produced				
			 if ( !-e $bam_aligned ){
				print_and_log("\n$analysis_indication Merging all the SAM files in a unique BAM..\n",$log_file);
				my $main_name = $cfg_hash->{'db_readf_name'};
				
				#Here get the read files names to be merged using the sample name and the step needed (align)
				my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$cfg_hash->{'align_step'},$inFolder,\@steps_array,$cfg_hash->{'sam_ext'});
				
				my $samples_to_merge  = " I=".join(" I=",@paths)." ";
				print_and_log("$analysis_indication Samples: $samples_to_merge..\n",$log_file);
				my $prog_used = $cfg_hash->{'mergesam_prog'};	
				print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
				#Merge all the bam files from the alignment in a single bam 
				run_PICARD_MergeSamFiles_gen($cfg_hash,$prog_used,$samples_to_merge,$bam_aligned,$log_file);####COMMENTTOREMOVE###

				update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
							$cfg_hash->{'db_sample_id'},$sample_id,$cfg_hash->{'mergealn_step'},$cfg_hash->{'db_analysis_id'},$analysis_id);
																
				#Remove the SAM files of the aligned and sorted 
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					foreach my $path (@paths){
						print_and_log("$analysis_indication Removing the alignment $path file..\n",$log_file);	
						delete_file($path);						
					}
				}
			}	
		}	
	}else{
	
		#The sample has a single read file. In this case we can access its name from the 
		#readf_table with the sample id that we are carrying 
		print_and_log( "\n$analysis_indication There are no multiple samples..\n",$log_file);
		my $params;
		my $main_name = $cfg_hash->{'db_readf_name'};
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		#If there are no multiple read files, the BAM will not be merged and take the name that they already have
		$bam_input = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@steps_array).".".$cfg_hash->{'bam_ext'};
		
		$bam_aligned = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'sort_idx_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@steps_array).".".$cfg_hash->{'bam_ext'};
		print_and_log("$analysis_indication bam input: $bam_input\n",$log_file);

		#If we have to execute flagstat on the bam aligned and sorted we need to 
		# get the bam aligned. If only the SAM file is there then convert with view
		if (! (-e $bam_aligned) or -z  $bam_aligned ){
			#Construct the same name with SAM extension
			my $sam_aligned = extract_name($bam_aligned,"noext").".".$cfg_hash->{'sam_ext'};
			if (-e $sam_aligned){
				print_and_log("$analysis_indication RUNNING  View : ",$log_file);
				#If only the SAM file exists then convert in BAM, sort and index it
				my $bam_out = $inFolder."/".extract_name($sam_aligned,1).".".$cfg_hash->{'bam_ext'};
				run_SAMTOOLS_View_gen($cfg_hash,$sam_aligned,$bam_out,$log_file);####COMMENTTOREMOVE###
				print_and_log("$analysis_indication RUNNING  Sort : ",$log_file);
				run_SAMTOOLS_Sort_gen($cfg_hash,$bam_out,$bam_aligned,$log_file);####COMMENTTOREMOVE###
				print_and_log("$analysis_indication RUNNING  Index : ",$log_file);
				run_SAMTOOLS_Index_gen($cfg_hash,$bam_aligned,$log_file);####COMMENTTOREMOVE###
				#Remove the BAM files not sorted
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					delete_file($bam_out);						
				}						
			}else{
					log_and_exit("$analysis_indication File $sam_aligned does not exists. Please run the alignment...\n",$log_file);
			}
		}		
	}
	#Finally get the name of the flagstat output
	$flagstat_recalibrated = $outFolder."/".extract_name($bam_recalibrated,1).".".$cfg_hash->{'flagstat_ext'};
	$flagstat_aln = $outFolder."/".extract_name($bam_aligned,1).".".$cfg_hash->{'flagstat_ext'};
	
	
	
	#########	
	#STEP2  #######GET NON-COVERED REGIONS SECTION
	#########
	#BEDTOOLS_coverageBed execution
	my $task_name = 'noncovered_regions';
	my $noncovgenes_out = $outFolder."/".extract_name($bam_input,1).".".$cfg_hash->{'noncov_ext'}.$cfg_hash->{'noncov_threshold'};
	if ($cfg_hash->{'noncovered_regions'} eq 'YES' ){
			#METHOD WITH GENOMECOV. IT IS EXECUTED AUTOMATICALLY
			if ( -e  $bam_recalibrated or -z  $bam_recalibrated){
				print_and_log( "\n$analysis_indication Starting the step to get $task_name : ",$log_file);
				#Get non covered regions
				get_non_cov_reg_with_bedtools($cfg_hash,$bam_input,$noncovgenes_out,$analysis_id,$log_file);
				#print_and_log( "\n$analysis_indication The output file is: $noncovgenes_out..\n",$log_file);#DEBUGCODE
				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished get $task_name $partDuration seconds\n",$log_file);
				$partTime = time;		
			}else{print_and_log( "$analysis_indication ERROR: Cannot get $task_name for sample $sample_id. Input file $bam_input does not exists.\n",$log_file);}
	}

	########################################
	##########FILLING OUTLIST FILE##########
	########################################		
	#Putting the path to the BAM sorted and eventually merged for multiple samples
	#if ( check_presence($bam_input) ){
	#	print_and_log( "$analysis_indication Appending $bam_input path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
	#	my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_sortedbam_desc'}."\t".$bam_input;
	#	append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	#}else{
	#	print_and_log( "$analysis_indication WARNING: $bam_input is not present. Please check what is happen. \n",$log_file);	
	#}		
	
	#Put the path to the non covered regions file in the outlist file
	if ( check_presence($noncovgenes_out) ){
		print_and_log( "$analysis_indication Appending $noncovgenes_out path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_final_noncov'}."\t".$noncovgenes_out;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	}else{
		print_and_log( "$analysis_indication WARNING: $noncovgenes_out is not present. Please check what is happen. \n",$log_file);	
	}		



	#########	
	#STEP3#######FLAGSTATS SECTION
	##################

	#############################
	#SAMTOOLS_flagstat execution using the only aligned file (without remove dup)
	#############################
	if ($cfg_hash->{'flagstat_sorted'} eq 'YES' ){

		#Now, using the just constructed aln_input, run flagstat
		my $prog_used = $cfg_hash->{'flagstat_prog'};	
		print_and_log("$analysis_indication RUNNING $prog_used : \n",$log_file);		
		$flagstat_aln = $outFolder."/".extract_name($bam_aligned,1).".".$cfg_hash->{'flagstat_ext'};	
		run_SAMTOOLS_flagstat($cfg_hash,$bam_aligned,$flagstat_aln,$log_file);####COMMENTTOREMOVE###

		########################################
		##########FILLING OUTLIST FILE##########
		########################################	
		#Put the path to the flagstat output in the outlist file
		if ( check_presence($flagstat_aln) ){
			print_and_log( "$analysis_indication Appending $flagstat_aln path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
			my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_FS_alignsort_desc'}."\t".$flagstat_aln;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
		}else{
			print_and_log( "$analysis_indication WARNING: $flagstat_aln is not present. Please check what is happen. \n",$log_file);	
		}		

		#Remove the SAM files of the aligned and sorted used for the statistics
		if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			print_and_log("$analysis_indication Removing the alignment $bam_aligned file..\n",$log_file);	
			delete_file($bam_aligned);						
		}	
		
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used (flagstat_sorted) $partDuration seconds\n",$log_file);
		$partTime = time;		
	}#END FLAGSTAT ON ALIGNED (NO RMDUP) BAM FILE
	
	
	
	#############################
	#SAMTOOLS_flagstat execution using the alignment and removed duplicates (if present)
	############################
	if ($cfg_hash->{'flagstat_baserecal'} eq 'YES' ){
			#Then the BAM file should be there
			if (-e $bam_recalibrated){ 
				#Flagstat
				$flagstat_recalibrated = $outFolder."/".extract_name($bam_input,1).".".$cfg_hash->{'flagstat_ext'};	
				my $prog_used = $cfg_hash->{'flagstat_prog'};	
				print_and_log("$analysis_indication SAMTOOLS_flagstat execution using to get statistics for the file with removed duplicates..\n",$log_file);	
				run_SAMTOOLS_flagstat($cfg_hash,$bam_input,$flagstat_recalibrated,$log_file);####COMMENTTOREMOVE###
					
				#Put the path to the flagstat output in the outlist file
				if ( check_presence($flagstat_recalibrated) ){
					print_and_log("$analysis_indication Appending $flagstat_recalibrated path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
					my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_FS_remdup_desc'}."\t".$flagstat_recalibrated;
					append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
				}else{
					print_and_log("$analysis_indication WARNING: $flagstat_recalibrated is not present. Please check what is happen. \n",$log_file);	
				}		
					
				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished\t$prog_used (flagstat_rmdup)\tsample\t$sample_id\t$partDuration\tseconds\n",$log_file);
				$partTime = time;			
			}else{
				print_and_log( "$analysis_indication WARNING: Cannot start flagstat_rmdup: $bam_input is not present. Please check what is happen. \n",$log_file);
			}
	}
	
	##############################
	# Make a table with alignment statistics using the results from flagstat
	# and add the path to statstables_file
	##############################
	
	my $plot_name = $cfg_hash->{'ALIGNMENTSTATS_name'};
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	my $reads_stats_tab = $outFolder."/$analysis_name\_$plot_name.totab";	
	#If at least one of the two has been ran
	if ( ! ($flagstat_aln eq ''  and $flagstat_recalibrated eq '') ){
		print_and_log( "$analysis_indication Writing a table $reads_stats_tab with statistics from flagstat $flagstat_aln and $flagstat_recalibrated ... \n",$log_file);	
		#Here I want to fill a table using the data extracted from the flagstat output
		print_alignment_statistics_table($cfg_hash,$flagstat_aln,$flagstat_recalibrated,$reads_stats_tab,$analysis_id,$sample_name,$log_file );	
		#Obtain the analysis name of the current analysis from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);			
		#append the path to statstables_file
		append_str_2_file_if_path_notexist($cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{'finalout_out_f'}."/".$cfg_hash->{'statstables_file'},"Alignment Statistics\t".$reads_stats_tab);

	}else{
		print_and_log( "$analysis_indication WARNING: Both $flagstat_aln and $flagstat_recalibrated are empty. Please check what is happen. \n",$log_file);	
	}
}

=head2 run_coverage_stepsOLD

 Title   : run_coverage_steps
 Usage   : run_coverage_steps(  config_file => the config hash
								);

 Function: The coverage analysis implies different steps. In particular this function
					is executed per-sample.
					
					Two main files are used to get the statistics:	
						- the BAM file just after the BWA alignment
						- the BAM file where the duplicates have been removed.
					If the DuplicateRemoval was not executed (e.g. targeted sequencing), then
					the second BAM will not be present.
					
					STEP 1: Construct the name of the last BAM (bam_input)
						IF sample is divided into multiple read files	:					
							In case there are multiple read files for the sample, these must be merged into a single
							and the duplicate removal must be performed again.
							The input BAM files per-lane can be removed after this step
							(if the variant calling is already finished)
							Before to remove I check if this task is the last one that is using the lane files
							(The variant calling step, for example, uses them.)
						IF the sample bam file is unique the name is reconstucted using the database
					
					STEP2: get non covered regions.
						Using $bam_input a procedure to get the non-covered regions is ran
						
					
					STEP3: Flagstat Statistics
						Aligned:
						Flagstat is executed first with the bam file just after the alignment. 
						In this case we do not have a bam file but the SAM file after the alignment
						So we merge these files into a single BAM
					
						With duplicate Removed:
						Flagstat is also executed on the file with duplicate removed that is the
						$bed_merged_rmdup previously constructed using $cfg_hash->{'mark_rem_dup_step'}
						parameter. If it exists, of course.
						
						Intersected with the target file:
						Last alignment statistics are performed on the bam file intersected with the target
						$bam_input is used for the intersection.			
 
 Returns :
 
=cut
sub run_coverage_stepsOLD{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	
	my $analysis_indication = "\n[an:$analysis_id sam: $sample_id] (".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("\n$analysis_indication Starting alignment statistics with Flagstat and non-covered region analysis [per-sample]..\n",$log_file);
	
	
	my $step_needed = $cfg_hash->{'sort_idx_step'};
	#Checks if sort and indexing has been performed for each read file
	update_sample_status($cfg_hash,$sample_id,$cfg_hash->{'sort_idx_step'});
	
	my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
	my $out_fold_suff = $cfg_hash->{'stats_step'};#The folder that will be used for the output is the one for statistics

	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	
	#Get also the sample_name
	my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	#Verifying that there has been an execution of multiple samples
	my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    
	#Check if markDuplicate after the alignment was executed
	#Update the sample status
	update_sample_status($cfg_hash,$sample_id,$cfg_hash->{'mark_rem_dup_step'});
	my $mark_rem_dup_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mark_rem_dup_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
	#Checks if sort and indexing has been performed for each read file
	my $mergebam_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mergebam_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);			
	my $mergealn_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mergealn_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);
						    						    
	###GET THE NAMES OF THE FILES
	my $bam_input = "";
	my $bam_aligned = "";
 
	#Names of the flagstat outputs
	my $flagstat_aln = "";
	my $flagstat_aln_sort_mrdup = "";
	my $flagstat_aln_mrdup_tinters = "";
	
	#########
	####STEP 1: check if there are multiple read files. All the steps will use a unique bam file for the sample.		    
	#########
	
	#If there are more read_files per sample than the merge is needed merge them						        						    
	if ( $multiple == 1){
		print_and_log( "\n$analysis_indication Sample $sample_name is sub-divided in multiple read files. A merge is needed..\n",$log_file);
		
		my $params;
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		#The bam file merged with duplicates removed and sorted
		$bam_input = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'mergebam_step'}."_".$cfg_hash->{'sort_idx_step'}.".".$cfg_hash->{'bam_ext'};

		#The merged file from only aligned
		$bam_aligned = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'align_step'}.".".$cfg_hash->{'bam_ext'};
								
								
		print_and_log("\n$analysis_indication Input BAM: $bam_input - BAM aligned only: $bam_aligned..\n",$log_file);#DEBUGCODE
		#If readfiles are multiple and one among the flagstat and the non-covered region step
		#must be executed, they must be merged
		if ($cfg_hash->{'flagstat_rmdup'} eq 'YES' 
					or $cfg_hash->{'noncovered_regions'} eq 'YES' 
						or $cfg_hash->{'flagstat_inters_targ'} eq 'YES' ){
			#If the merged BAM file has not been already produced				
			 if ( !-e $bam_input ){
					if ( $cfg_hash->{'mergebam'} eq 'YES' ){
						#Here get the samples names to be merged
						my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$cfg_hash->{'mark_rem_dup_step'},$inFolder,\@steps_array,$cfg_hash->{'bam_ext'});		

						print_and_log( "\n$analysis_indication Starting MergeBam to merge multiple samples..\n",$log_file);

						my $samples_to_merge  = " I=".join(" I=",@paths)." ";
						#Create the name for the bam file of the merge from the lane files
						my $bam_input_merged = $inFolder."/".build_input_name_from_executed($params,$cfg_hash->{'mergebam_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@stats_steps_array)."_".$cfg_hash->{'mergebam_step'}.".".$cfg_hash->{'bam_ext'};						
						my $prog_used = $cfg_hash->{'mergebam_prog'};
						print_and_log( "RUNNING $prog_used : ",$log_file);#DEBUGCODE
						run_PICARD_MergeSamFiles_gen($cfg_hash,$prog_used,$samples_to_merge,$bam_input_merged,$log_file);
						
						my $bam_input_merge_sort = $inFolder."/".extract_name($bam_input_merged,1)."_".$cfg_hash->{'sort_idx_step'}.".".$cfg_hash->{'bam_ext'};
						print_and_log("$analysis_indication RUNNING Sort : ",$log_file);
						run_SAMTOOLS_Sort_gen($cfg_hash,$bam_input_merged,$bam_input_merge_sort,$log_file);
						print_and_log("$analysis_indication RUNNING Index : ",$log_file);
						run_SAMTOOLS_Index_gen($cfg_hash,$bam_input_merge_sort,$log_file);						

						#Save status in the db by searching for the specific sample name
						#print_and_log("$analysis_indication Updating analysis status for sample $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
						update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
							$cfg_hash->{'db_sample_id'},$sample_id,$cfg_hash->{'mergebam_step'},$cfg_hash->{'db_analysis_id'},$analysis_id);
						
						#Remove the mergebam file before the sortidx
						if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
							delete_file($bam_input_merged);	
						}
						
						##################################	
						#READ FILES REMOVAL
						#################################
						#Before to remove I check if this task is the last one that is using the lane files
						#The variant calling step, for example, uses them. 
						if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
							print_and_log( "\n$analysis_indication Removing read files\n",$log_file);
							#Check if this is the last task using the BAM files 
							my $bamfilesused = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
												$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},
												$cfg_hash->{'db_sample_id'},$sample_id);	
							if ( $bamfilesused == ($cfg_hash->{'tasks_using_bam'} - 1) )	{
								foreach my $bam_path (@paths){
									print_and_log("$analysis_indication Removing the alignment $bam_path file..\n",$log_file);	
									delete_file($bam_path);		
									delete_file($bam_path.".".$cfg_hash->{'bai_ext'});		
								}
								$bamfilesused++;
								update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},
												$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_bamfilesused'},$bamfilesused,$cfg_hash->{'db_sample_id'},$sample_id);								
							}else{
								print_and_log( "\n$analysis_indication Read files will not be removed. Still ".($cfg_hash->{'tasks_using_bam'} - 1) - $bamfilesused." steps have to use them..\n",$log_file);
							}
						}							
						#Print time needed
						$partDuration = time - $partTime;
						print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
						$partTime = time;					
					}
			 }		
		}
		
		
		#Now check if flagstat for aligned read file must be executed and in case do the merge
		if ($cfg_hash->{'flagstat_sorted'} eq 'YES' and $cfg_hash->{'mergebam'} eq 'YES'){
			 #If the merged aligned BAM file has not been already produced				
			 if ( !-e $bam_aligned ){
				print_and_log("\n$analysis_indication Merging all the SAM files in a unique BAM..\n",$log_file);
				my $main_name = $cfg_hash->{'db_readf_name'};
				
				#Here get the read files names to be merged using the sample name and the step needed (align)
				my @paths = get_readf_paths_using_step($cfg_hash,$sample_id,$log_file,$cfg_hash->{'align_step'},$inFolder,\@steps_array,$cfg_hash->{'sam_ext'});
				
				my $samples_to_merge  = " I=".join(" I=",@paths)." ";
				print_and_log("$analysis_indication Samples: $samples_to_merge..\n",$log_file);
				my $prog_used = $cfg_hash->{'mergesam_prog'};	
				print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
				#Merge all the bam files from the alignment in a single bam 
				run_PICARD_MergeSamFiles_gen($cfg_hash,$prog_used,$samples_to_merge,$bam_aligned,$log_file);####COMMENTTOREMOVE###

				update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},
							$cfg_hash->{'db_sample_id'},$sample_id,$cfg_hash->{'mergealn_step'},$cfg_hash->{'db_analysis_id'},$analysis_id);
																
				#Remove the SAM files of the aligned and sorted 
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					foreach my $path (@paths){
						print_and_log("$analysis_indication Removing the alignment $path file..\n",$log_file);	
						delete_file($path);						
					}
				}
			}	
		}	
	}else{
	
		#The sample has a single read file. In this case we can access its name from the 
		#readf_table with the sample id that we are carrying 
		print_and_log( "\n$analysis_indication There are no multiple samples..\n",$log_file);
		my $params;
		my $main_name = $cfg_hash->{'db_readf_name'};
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		#If there are no multiple read files, the BAM will not be merged and take the name that they already have
		$bam_input = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@steps_array).".".$cfg_hash->{'bam_ext'};
		
		$bam_aligned = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'sort_idx_step'},
								$params->{$cfg_hash->{'db_sample_name'}},\@steps_array).".".$cfg_hash->{'bam_ext'};
		print_and_log("$analysis_indication bam input: $bam_input\n",$log_file);

		#If we have to execute flagstat on the bam aligned and sorted we need to 
		# get the bam aligned. If only the SAM file is there then convert with view
		if (! (-e $bam_aligned) ){
			#Construct the same name with SAM extension
			my $sam_aligned = extract_name($bam_aligned,"noext").".".$cfg_hash->{'sam_ext'};
			if (-e $sam_aligned){
				print_and_log("$analysis_indication RUNNING  View : ",$log_file);
				#If only the SAM file exists then convert in BAM, sort and index it
				my $bam_out = $inFolder."/".extract_name($sam_aligned,1).".".$cfg_hash->{'bam_ext'};
				run_SAMTOOLS_View_gen($cfg_hash,$sam_aligned,$bam_out,$log_file);####COMMENTTOREMOVE###
				print_and_log("$analysis_indication RUNNING  Sort : ",$log_file);
				run_SAMTOOLS_Sort_gen($cfg_hash,$bam_out,$bam_aligned,$log_file);####COMMENTTOREMOVE###
				print_and_log("$analysis_indication RUNNING  Index : ",$log_file);
				run_SAMTOOLS_Index_gen($cfg_hash,$bam_aligned,$log_file);####COMMENTTOREMOVE###
				#Remove the BAM files not sorted
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					delete_file($bam_out);						
				}						
			}else{
					log_and_exit("$analysis_indication File $sam_aligned does not exists. Please run the alignment...\n",$log_file);
			}
		}		
	}
	#Finally get the name of the flagstat output
	$flagstat_aln_sort_mrdup = $outFolder."/".extract_name($bam_input,1).".".$cfg_hash->{'flagstat_ext'};
	$flagstat_aln = $outFolder."/".extract_name($bam_aligned,1).".".$cfg_hash->{'flagstat_ext'};
	
	
	
	#########	
	#STEP2  #######GET NON-COVERED REGIONS SECTION
	#########
	#BEDTOOLS_coverageBed execution
	my $task_name = 'noncovered_regions';
	my $noncovgenes_out = $outFolder."/".extract_name($bam_input,1).".".$cfg_hash->{'noncov_ext'}.$cfg_hash->{'noncov_threshold'};
	if ($cfg_hash->{'noncovered_regions'} eq 'YES' ){
			#METHOD WITH GENOMECOV. IT IS EXECUTED AUTOMATICALLY
			if ( -e  $bam_input){
				print_and_log( "\n$analysis_indication Starting the step to get $task_name : ",$log_file);
				#Get non covered regions
				get_non_cov_reg_with_bedtools($cfg_hash,$bam_input,$noncovgenes_out,$analysis_id,$log_file);
				#print_and_log( "\n$analysis_indication The output file is: $noncovgenes_out..\n",$log_file);#DEBUGCODE
				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished get $task_name $partDuration seconds\n",$log_file);
				$partTime = time;		
			}else{print_and_log( "$analysis_indication ERROR: Cannot get $task_name for sample $sample_id. Input file $bam_input does not exists.\n",$log_file);}
	}

	########################################
	##########FILLING OUTLIST FILE##########
	########################################		
	#Putting the path to the BAM sorted and eventually merged for multiple samples
	#if ( check_presence($bam_input) ){
	#	print_and_log( "$analysis_indication Appending $bam_input path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
	#	my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_sortedbam_desc'}."\t".$bam_input;
	#	append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	#}else{
	#	print_and_log( "$analysis_indication WARNING: $bam_input is not present. Please check what is happen. \n",$log_file);	
	#}		
	
	#Put the path to the non covered regions file in the outlist file
	if ( check_presence($noncovgenes_out) ){
		print_and_log( "$analysis_indication Appending $noncovgenes_out path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_final_noncov'}."\t".$noncovgenes_out;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	}else{
		print_and_log( "$analysis_indication WARNING: $noncovgenes_out is not present. Please check what is happen. \n",$log_file);	
	}		



	#########	
	#STEP3#######FLAGSTATS SECTION
	##################

	#############################
	#SAMTOOLS_flagstat execution using the only aligned file (without remove dup)
	#############################
	if ($cfg_hash->{'flagstat_sorted'} eq 'YES' ){

		#Now, using the just constructed aln_input, run flagstat
		my $prog_used = $cfg_hash->{'flagstat_prog'};	
		print_and_log("$analysis_indication RUNNING $prog_used : \n",$log_file);		
		$flagstat_aln = $outFolder."/".extract_name($bam_aligned,1).".".$cfg_hash->{'flagstat_ext'};	
		run_SAMTOOLS_flagstat($cfg_hash,$bam_aligned,$flagstat_aln,$log_file);####COMMENTTOREMOVE###

		########################################
		##########FILLING OUTLIST FILE##########
		########################################	
		#Put the path to the flagstat output in the outlist file
		if ( check_presence($flagstat_aln) ){
			print_and_log( "$analysis_indication Appending $flagstat_aln path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
			my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_FS_alignsort_desc'}."\t".$flagstat_aln;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
		}else{
			print_and_log( "$analysis_indication WARNING: $flagstat_aln is not present. Please check what is happen. \n",$log_file);	
		}		

		#Remove the SAM files of the aligned and sorted used for the statistics
		if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			print_and_log("$analysis_indication Removing the alignment $bam_aligned file..\n",$log_file);	
			delete_file($bam_aligned);						
		}	
		
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished $prog_used (flagstat_sorted) $partDuration seconds\n",$log_file);
		$partTime = time;		
	}#END FLAGSTAT ON ALIGNED (NO RMDUP) BAM FILE
	
	
	
	#############################
	#SAMTOOLS_flagstat execution using the alignment and removed duplicates (if present)
	############################
	if ($cfg_hash->{'flagstat_rmdup'} eq 'YES' ){
		#Check if the MarkDuplicates step after the alignment was performed and if the bam file exists
		if ( $mark_rem_dup_executed  ){
			#Then the BAM file should be there
			if (-e $bam_input){ 
				#Flagstat
				$flagstat_aln_sort_mrdup = $outFolder."/".extract_name($bam_input,1).".".$cfg_hash->{'flagstat_ext'};	
				my $prog_used = $cfg_hash->{'flagstat_prog'};	
				print_and_log("$analysis_indication SAMTOOLS_flagstat execution using to get statistics for the file with removed duplicates..\n",$log_file);	
				run_SAMTOOLS_flagstat($cfg_hash,$bam_input,$flagstat_aln_sort_mrdup,$log_file);####COMMENTTOREMOVE###
					
				#Put the path to the flagstat output in the outlist file
				if ( check_presence($flagstat_aln_sort_mrdup) ){
					print_and_log("$analysis_indication Appending $flagstat_aln_sort_mrdup path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
					my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_FS_remdup_desc'}."\t".$flagstat_aln_sort_mrdup;
					append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
				}else{
					print_and_log("$analysis_indication WARNING: $flagstat_aln_sort_mrdup is not present. Please check what is happen. \n",$log_file);	
				}		
					
				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished\t$prog_used (flagstat_rmdup)\tsample\t$sample_id\t$partDuration\tseconds\n",$log_file);
				$partTime = time;			
			}else{
				print_and_log( "$analysis_indication WARNING: Cannot start flagstat_rmdup: $bam_input is not present. Please check what is happen. \n",$log_file);
			}
		}else{
			print_and_log( "$analysis_indication I will not start flagstat_rmdup because duplicates are not removed...\n",$log_file);
		}
	}
	
	
	##############################
	##SAMTOOLS_flagstat execution using the aligned against the reference and intersected with the target
	##This step is executed on the last file obtained (bam_input)
	##############################
	#if ($cfg_hash->{'flagstat_inters_targ'} eq 'YES' ){
		#if ( -e $bam_input){
			#print_and_log("$analysis_indication SAMTOOLS_flagstat execution using the aligned and intersected with the target..\n",$log_file);	
			#my $main_name = $cfg_hash->{'db_sample_name'};
			#my $step_needed = $cfg_hash->{'sort_idx_step'};
			#my $prog_used;
			#my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
			#my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};

			##Intersect with the target file
			##Getting the target bed file using the folder for targets and the name contained in the database
			#my $target_bed = $cfg_hash->{'target_reg_f'}."/".get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},	$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'},
											#$cfg_hash->{'db_analysis_id'},$analysis_id);
			##Check if it exists								
			#if (file_not_present($target_bed)  > 0){ die "Cannot proceed with $prog_used! Check: $target_bed.\n";}
			
			#$prog_used = $cfg_hash->{'intersectbed_prog'};	
			#print_and_log("$analysis_indication RUNNING $prog_used : ",$log_file);
			#my $interTargCodBed = $inFolder."/".extract_name($bam_input,1)."_inters_".extract_name($target_bed,1).".".$cfg_hash->{'bed_ext'};	
			##Set parameters for intersect bed. giving None so that reports the shared interval between the two overlapping features
			#my $inters_params = " ";
			#run_BEDTOOLS_intersectBed($cfg_hash,$bam_input,$target_bed,$interTargCodBed,$inters_params,$log_file);####COMMENTTOREMOVE###
			#my $flagstat_aln_mrdup_tinters = $outFolder."/".extract_name($interTargCodBed,1).".".$cfg_hash->{'flagstat_ext'};	
			
			##Flagstat
			#$prog_used = $cfg_hash->{'flagstat_prog'};	
			#print_and_log("$analysis_indication  RUNNING $prog_used..\n",$log_file);	
			#run_SAMTOOLS_flagstat($cfg_hash,$interTargCodBed,$flagstat_aln_mrdup_tinters,$log_file);####COMMENTTOREMOVE###

			#########################################
			###########FILLING OUTLIST FILE##########
			#########################################	
			##Put the path to the flagstat output in the outlist file

			##Delete the BAM file intersected with the target
			#if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
				#print_and_log("$analysis_indication Removing  $interTargCodBed..\n",$log_file);	
				#delete_file($interTargCodBed);
			#}
						
			##Put the path to the flagstat output in the outlist file
			#if ( check_presence($flagstat_aln_mrdup_tinters) ){
				#print_and_log( "$analysis_indication Appending $flagstat_aln_mrdup_tinters path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
				#my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_FS_intertarg_desc'}."\t".$flagstat_aln_mrdup_tinters;
				#append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
			#}else{
				#print_and_log( "$analysis_indication WARNING: $flagstat_aln_mrdup_tinters is not present. Please check what is happen. \n",$log_file);	
			#}		
			
			#$partDuration = time - $partTime;
			#print_and_log( "$analysis_indication Finished  $prog_used (flagstat_inters_targ) $partDuration seconds\n",$log_file);
			#$partTime = time;		
		#}	
	#}	
	
	##############################
	# Make a table with alignment statistics using the results from flagstat
	# and add the path to statstables_file
	##############################
	
	my $plot_name = $cfg_hash->{'ALIGNMENTSTATS_name'};
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);
	my $reads_stats_tab = $outFolder."/$analysis_name\_$plot_name.totab";	
	#If at least one of the two has been ran
	if ( ! ($flagstat_aln eq ''  and $flagstat_aln_sort_mrdup eq '') ){
		print_and_log( "$analysis_indication Writing a table $reads_stats_tab with statistics from flagstat $flagstat_aln and $flagstat_aln_sort_mrdup ... \n",$log_file);	
		#Here I want to fill a table using the data extracted from the flagstat output
		print_alignment_statistics_table($cfg_hash,$flagstat_aln,$flagstat_aln_sort_mrdup,$reads_stats_tab,$analysis_id,$sample_name,$log_file );	
		#Obtain the analysis name of the current analysis from the database given the group id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
																$cfg_hash->{'db_analysis_id'},$analysis_id);			
		#append the path to statstables_file
		append_str_2_file_if_path_notexist($cfg_hash->{'work_fold'}."/".$analysis_name."/".$cfg_hash->{'finalout_out_f'}."/".$cfg_hash->{'statstables_file'},"Alignment Statistics\t".$reads_stats_tab);

	}else{
		print_and_log( "$analysis_indication WARNING: Both $flagstat_aln and $flagstat_aln_sort_mrdup are empty. Please check what is happen. \n",$log_file);	
	}
}

=head2 run_cnv_finalout_steps

 Title   : run_cnv_finalout_steps
 Usage   : run_cnv_finalout_steps(  config_file => the config hash
								);

 Function: CNV pipeline output
 
			Starting from the output of the CNV detection programs,
			parses and obtains a tabular output containing the final CNVs
			
			ExomeDepht: ExomeDepth prints an output for each of the BAM used for 
			the analysis. The path to the BAMs are contained into a file called 
			(ED_bamlist_f from config file). 
			This file can be used also to build the ExomeDepth output names.
			 1- each output from ExomeDepht must go into the own analysis folder
			 2. only the output files for this anaysis remain in the folder
 
 Returns :
 
 
=cut
sub run_cnv_finalout_steps{
	
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the analysis name from the initial sample_id
	my $analysis_name = $sample_id;
	
	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("$analysis_indication Starting CNV pipeline..\n",$log_file);

	#The input folder is given by the working folder  (vargenius_analyses) and the analysis name
	my $temp = "/pico/work/TELET_TIGEM/vargenius_analyses/";
			
	my $outFolder = $cfg_hash->{$analysis_id.'_cnv_out_f'};
	
	
	#############
	##ExomeDepth output management
	##1. Put all the output from ExomeDepth for the different samples into the indicated analysis folder
	ED_reallocate_output($cfg_hash,$analysis_id);
	
	#my $ed_out_paths_f = $outFolder."/". $cfg_hash->{'ED_cnvout_f'};
	##The file format is the following:
	##sampleid = analysisname]--[exomedepth.path.auto,exomedepth.path.sex
	#open (ED_OUTS,"<$ed_out_paths_f") or die "ERROR: Cannot open $ed_out_paths_f\n";
	#while (my $row = <ED_OUTS>){
		#chop($row);
		#my @fields = split(" = ",$row);
		#my $sampleid = $fields[0];
		
		#my @fields2 = split($cfg_hash->{'mult_ann_sep'},$fields[1]);
		#my $curr_anal_name = $fields2[0];
		
		#my @fields3 = split(",",$fields2[1]);
		#my $auto_o = $fields3[0];
		#my $sex_o = $fields3[1];
		
		##If exists, put the given ExomeDepht output into the finalout_out/cnv/sex|auto folder
		#my $cnvoutFolder = $temp."/".$curr_anal_name."/".$cfg_hash->{'fastfout_out_f'}."/".$cfg_hash->{'cnv_out_f'};
		#if ( -e $auto_o and ! (-z $auto_o)){
			#move($auto_o,$cnvoutFolder."/".$cfg_hash->{'cnv_auto_f'}."/") or die "ERROR: Cannot move $auto_o to $cnvoutFolder/".$cfg_hash->{'cnv_auto_f'}."/";
		#}
		#if ( -e $sex_o and ! (-z $sex_o)){
			#move($sex_o,$cnvoutFolder."/".$cfg_hash->{'cnv_sex_f'}."/") or die "ERROR: Cannot move $sex_o to $cnvoutFolder/".$cfg_hash->{'cnv_auto_f'}."/";
		#}		
	#}
	#close (ED_OUTS);
	
	######
	#MERGE OUTPUT FROM EXOMEDEPTH AND XHMM
	###
	#To merge the results we need a preprocessing:
	#	1. The output from ExomeDepth contains CNVs found in the single samples  hence for ExomeDepth we take the ouput for all samples for the analysis and concatenate them together.
	#		to use this file later we need to add a column specifing which sample is it.
	#	2. XHMM's contains CNVs for all samples analyzed, hence we need to grep results only for the samples of the subject analysis
	#	3. We will use R merge to merge the ouputs and obtain a starting table. But we need a common column to use. The "id" column in exome depth
	#		and INTERVAL in XHMM are exactly what we need (chromosome and coordinates chr1:69092-70008) but we need to add if it is a DEL or a DUP.
	 
	my $cnvoutFolder = $cfg_hash->{$analysis_id.'_cnv_out_f'};
	 
	#Operation 1:
	#In Exome Depth output, for each sample of the analysis find the output, add the column with the sample name and concatenate all the files
	###
	my $ed_out = $cnvoutFolder."exomedepth.txt";
	ED_output_preprocessing	($cfg_hash,$analysis_id,$cnvoutFolder,\@steps_array2,$ed_out,$log_file);
	
	#my $ed_out_temp = $cnvoutFolder."exomedepth.txt.temp";	
	
	##Get sample ids for the analysis
	#my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}." = $analysis_id";
	#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	#my $anal_samples = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
	#my $auto_o = "";
	#my $sex_o = "";
	#foreach my $sample_id (keys %{$anal_samples}){
		##Get the analysis id and name for the given sample
		#my $curr_anal_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},
						    #$cfg_hash->{'db_sample_id'},$sample_id);
		#my $curr_anal_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    #$cfg_hash->{'db_analysis_id'},$curr_anal_id);	
		#my $curr_sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    #$cfg_hash->{'db_sample_id'},$sample_id);
						    
		##Get the ExomeDepth output
		##The step at which the BAM should be collected: the one before the variant calling hence after Base Recalibration	
		#my $step = $cfg_hash->{'varcall_step'};
		#my $paramsIn;
		#getSampleConfiguration_locked(\$paramsIn,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						#$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

		
		#$auto_o = $cnvoutFolder."/auto/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,\@steps_array2)
			#.".CNVs.tsv"; 
		#$sex_o = $cnvoutFolder."/sex/exomedepth.".build_input_name_from_executed($paramsIn,$step,$curr_sample_name,\@steps_array2)
			#.".CNVs.tsv";
			
		##Add the column for the sample name
		##First create an array to add as the first column
		#my @sample_name_col = ($cfg_hash->{'db_sample_name'});
		#my $num_lines = `wc -l $auto_o`;
		#for (my $i = 1; $i < $num_lines; $i++){
			#push(@sample_name_col,$curr_sample_name);
		#}
		##Add to both sex and auto outputs
		#insert_col_in_file_table($auto_o,\@sample_name_col,"start.p","i");
		#insert_col_in_file_table($sex_o,\@sample_name_col,"start.p","i");
		##Remove the headers
		#delete_rows_containing($auto_o,"start.p",$auto_o.".temp");
		#delete_rows_containing($sex_o,"start.p",$sex_o.".temp");
		
		##Concatenate the two files into a single ExomeDepth output		
		#append_file_2_file($auto_o.".temp",$ed_out_temp);
		#append_file_2_file($sex_o.".temp",$ed_out_temp);
	#}

	##Put the header
	#system ( " head -n1 $auto_o > $ed_out " ) or die "ERROR: Cannot execute head -n1";
	##Concatenate the two files into a single ExomeDepth output		
	#append_file_2_file($ed_out_temp,$ed_out);
	
	#Operation 2:
	#Get from the XHMM output only those lines concerning this analysis, hence starting with one of the subject samples names
	#
	my $xhmm_out = $cnvoutFolder."/xhmm.txt";
	XHMM_output_preprocessing($cfg_hash,$analysis_id,$cnvoutFolder,$xhmm_out,$log_file);


	#$auto_o = $cnvoutFolder."/auto/".$cfg_hash->{'xhmm_xcnv_out'};
	#$sex_o = $cnvoutFolder."/sex/".$cfg_hash->{'xhmm_xcnv_out'};
	##Put first the header
	#system ( " head -n1 $auto_o > $xhmm_out " ) or die "ERROR: Cannot execute head -n1";	
	#foreach my $sample_id (keys %{$anal_samples}){
		#my $curr_sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    #$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    #$cfg_hash->{'db_sample_id'},$sample_id);
						    		
		#get_rows_containing($auto_o,$curr_sample_name,$auto_o.".temp");
		#get_rows_containing($sex_o,$curr_sample_name,$sex_o.".temp");
		
		#append_file_2_file($auto_o.".temp",$xhmm_out);
		#append_file_2_file($sex_o.".temp",$xhmm_out);
	#}
	
	#Operation 3. Merge the output from XHMM and ExomeDepth
	my $cnv_out = $cnvoutFolder."/vargenius_cnv.txt";
	my $RLogPath = $cnvoutFolder."/".$cfg_hash->{'log_fold'}."/$analysis_name\_MERGE_CNV_OUTS".$cfg_hash->{'R_log_file'};
	my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_utils'};	
	my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args MERGE_CNV_OUTS $xhmm_out $ed_out $cnv_out' $R_script $RLogPath";
	
	print_and_log("The command is: $command\n",$log_file);
	try_exec_command($command) or R_die($command,$R_script);	
	
	#ANNOTATION
	#Finally annotate using OMIM		
			
}

=head2 run_cnv_calling_pipeline

 Title   : run_cnv_calling_pipeline
 Usage   : run_cnv_calling_pipeline(  config_file => the config hash
								);

 Function: CNV pipeline
 
			This subroutine obtains BAM files from samples sequenced with the same 
			target kit and launches ExomeDepth and XHMM.
 
 Returns :
 
=cut
sub run_cnv_calling_pipeline{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the analysis name from the initial sample_id
	my $analysis_name = $sample_id;
	
	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("$analysis_indication Starting CNV pipeline..\n",$log_file);

	#my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
	#my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};

	my $out_fold_suff = $cfg_hash->{'finalout_step'};#The folder that will be used for the output is CNV into final output folder
	my $outFolder = $cfg_hash->{$analysis_id.'_cnv_out_f'};
	
	#Run ExomeDepth
	print_and_log("$analysis_indication Running ExomeDepth..\n",$log_file);
	my $edjobids = "";
	$edjobids = run_exome_depth($cfg_hash,$outFolder,$analysis_name,$analysis_id,\@steps_array2,$task,$log_file);
	
	#Run XHMM
	print_and_log("$analysis_indication Running XHMM..\n",$log_file);
	my $xhmmjobids = run_xhmm_pipeline($cfg_hash,$outFolder,$analysis_name,$analysis_id,\@steps_array2,$task,$log_file);
	
	#The final output is parsed into another job for which dependencies are altered.
	#Now alter the dependencies of the variant filtering job: it has to wait the cnv detection job
	my $new_depend = "-W depend=afterok:";
	if ( $edjobids ne ''){ $new_depend .= "$edjobids";}
	if ( $xhmmjobids ne ''){ $new_depend .= ":$xhmmjobids";}
	my $job_name = $cfg_hash->{'cnv_finalout_step'}."_a".$analysis_id;
	my $cnvfout_jobid = get_jobid_form_jobname($job_name,$cfg_hash->{'qsub_username'});
	if ( $cnvfout_jobid > 0 ){ 
		print_and_log("Altering dependencies of ".$cfg_hash->{'cnv_finalout_step'}."_a".$analysis_id." ($cnvfout_jobid):\n ",$log_file);
		alter_job($cnvfout_jobid,$new_depend,$log_file);
		
		#Update the hash with dependencies
		my $deps_f = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'jc_dep_hash_f'};
		update_dependencies_hash($deps_f,$job_name,$cnvfout_jobid,$new_depend,$log_file);
	}
}

=head2 run_coverage_per_analysis_steps

 Title   : run_coverage_per_analysis_steps
 Usage   : run_coverage_per_analysis_steps(  config_file => the config hash
								);

 Function: Coverage per analysis using DepthOfCoverage
						 I see if there are multiple samples 
						#because in that case we need that the mergesam step has been executed 
						#and it is executed during the stats step. Hence it will block you and 
						#ask to run that before
						
						This subroutine launches DepthOfCoverage
 
 Returns :
 
=cut
sub run_coverage_per_analysis_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the analysis name from the initial sample_id
	my $analysis_name = $sample_id;
	
	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("$analysis_indication Starting coverage study using GATK-DepthOfCoverage [per-analysis]..\n",$log_file);

	my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
	my $out_fold_suff = $cfg_hash->{'stats_step'};#The folder that will be used for the output is the one for statistics

	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	
	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "$analysis_indication Executing: $query\n",$log_file);#DEBUGCODE
	my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#The following array is needed to save all the bam file names
	#to use in DepthOfCoverage
	my @paths_to_write = ();
	foreach my $sample_id (keys %{$group_sam}){
		#print_and_log("Verifying sample $sample_id: ",$log_file);#DEBUGCODE
			
			#I use getSampleConfiguration_locked just because I need the sample name
			#Out folder
			my $refinetask = "refine";
			my $refineOutFolder = $cfg_hash->{$analysis_id.'_'.$refinetask.'_out_f'};
			my $params;
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
							$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_id'},$sample_id);

			#get bam recalibrated
			my $bam_input = $refineOutFolder."/".get_output_name_from_executed($params,$cfg_hash->{'print_reads_step'},
					$params->{$cfg_hash->{'db_sample_name'}},\@steps_array2).".".$cfg_hash->{'bam_ext'};	
								
			print_and_log("no! -> samplename: ".$params->{$cfg_hash->{'db_sample_name'}}.". Using input file: $bam_input\n",$log_file);
			if ( -e $bam_input) {
				push(@paths_to_write,$bam_input);
			}else{
				die "$analysis_indication ERROR: Cannot find file $bam_input that is needed for $task. Please run alignment step..\n";
			}
	}
	
	#For DepthOfCoverage we can obtain a file with all the coverages per group
	#Now print a file with bam files paths from the different samples
	my $bam_paths_f = $outFolder."/$analysis_name\_bam_files_paths.".$cfg_hash->{'bamlist_ext'}; 
	#Open output from the annotation
	open (PATHS,">".$bam_paths_f) or die "$analysis_indication ERROR: Cannot open $bam_paths_f. The program will exit..\n";
	foreach my $path (@paths_to_write){
		print PATHS $path."\n";
	}
	close(PATHS);

	#########GET COVERAGE USING GATK DEPTHOFCOVERAGE
	my $task_name = 'coverage_analysis';
	my $covgenes_suff = $outFolder."/$analysis_name\_DOC";
	if ($cfg_hash->{$task_name} eq 'YES' ){
		if ( -e  $bam_paths_f){
			print_and_log( "$analysis_indication Starting the step to get $task_name : ",$log_file);
			run_GATK_DepthOfCoverage($cfg_hash,$bam_paths_f,$covgenes_suff,$analysis_id,$log_file);
			print_and_log( "$analysis_indication The output suffix is: $covgenes_suff..\n",$log_file);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $task_name $partDuration seconds\n",$log_file);
			$partTime = time;		
		}else{print_and_log( "$analysis_indication ERROR: Cannot get $task_name for analysis $analysis_name. File of paths $bam_paths_f does not exists.\n",$log_file);}
	}


	########################################
	##########FILLING OUTLIST FILE##########
	##########DepthOfCoverage files#####
	########################################	
	#Put the path to the coverage analysis files in the outlist file
	my @doc_outs = split(",",$cfg_hash->{'DOC_outputs'});
	foreach my $doc_out (@doc_outs){
		my $out_f = $covgenes_suff.".".$doc_out;
		if ( -e $out_f){
			#$doc_out =~ s/\_/ /g;
			
			print_and_log( "$analysis_indication Appending $out_f path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);	
			my $line = "-\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_cov_'.$doc_out}."\t".$out_f;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
		}else{
			print_and_log( "$analysis_indication WARNING: $out_f is not present. Please check what is happen. \n",$log_file);	
		}	
	}

}

=head2 run_coverage_per_analysis_steps

 Title   : run_coverage_per_analysis_steps
 Usage   : run_coverage_per_analysis_steps(  config_file => the config hash
								);

 Function: Coverage per analysis using DepthOfCoverage
						 I see if there are multiple samples 
						#because in that case we need that the mergesam step has been executed 
						#and it is executed during the stats step. Hence it will block you and 
						#ask to run that before
						
						This subroutine launches DepthOfCoverage
 
 Returns :
 
=cut
sub run_coverage_per_analysis_steps_ALN{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Get the analysis name from the initial sample_id
	my $analysis_name = $sample_id;
	
	my $analysis_indication = "\n[an:$analysis_id](".scalar(localtime)."): ";
		
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	print_and_log("$analysis_indication Starting coverage study using GATK-DepthOfCoverage [per-analysis]..\n",$log_file);

	my $in_fold_suff = $cfg_hash->{'align_step'};#The folder that will be used for the input is the one for alignments
	my $out_fold_suff = $cfg_hash->{'stats_step'};#The folder that will be used for the output is the one for statistics
	#my $statsFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$out_fold_suff.'_out_f'};
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$in_fold_suff.'_out_f'};
	
	#Get the samples ids involved for the group
	my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
	#print_and_log( "$analysis_indication Executing: $query\n",$log_file);#DEBUGCODE
	my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

	#With this bunch of code I see if there are multiple samples 
	#because in that case we need that the mergesam step has been executed 
	#and it is executed during the stats step. Hence it will block you and 
	#ask to run that before
	print_and_log("$analysis_indication Checking if there are multiple read file for each sample. Because in that case".
								"a step of merge should have been ran...\n",$log_file);
	#The following array is needed to save all the bam file names
	#to use in DepthOfCoverage
	my @paths_to_write = ();
	foreach my $sample_id (keys %{$group_sam}){
		#print_and_log("Verifying sample $sample_id: ",$log_file);#DEBUGCODE
		#Check if the markDuplicate after the readfiles alignment was executed
		my $mark_rem_dup_step = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mark_rem_dup_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
						    		
		#Check if markDuplicate after the readfiles merge was executed (This is optionally made after Base Recalibration)
		my $mrdup_groups_step_executed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'mrdup_groups_step'},
						    $cfg_hash->{'db_sample_id'},$sample_id);	
						    		
		#Verifying that there has been an execution of multiple samples
		my $multiple = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_multiple'},
					    $cfg_hash->{'db_sample_id'},$sample_id);		
		#If there are multiple read files than it should have been run the step for merge
		if ( $multiple == 1) {
			my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
						    $cfg_hash->{'db_sample_id'},$sample_id);			
			my $bed_merged = $inFolder."/".$sample_name."_".$cfg_hash->{'mergebam_step'};#NEWCOVERAGE
			
			#In the case we do not use the mark and remove duplicates
			if ( $mrdup_groups_step_executed == 1 ){
					$bed_merged .= "_".$cfg_hash->{'mark_rem_dup_step'};
			}
			#Add finally sortidx tag
			$bed_merged .= "\_".$cfg_hash->{'sort_idx_step'}.".".$cfg_hash->{'bam_ext'};
			
			#print_and_log("$analysis_indication yes! -> samplename: $sample_name. Using input file: $bed_merged\n",$log_file);#DEBUGCODE
			if ( -e $bed_merged) {
				push(@paths_to_write,$bed_merged);
			}else{
				die "$analysis_indication ERROR: Cannot find file $bed_merged that is needed for $task. Please run stats step..\n";
			}
		}#Otherwise go on with the result from the alignment
			#The sample has a single read file. In this case we can access its name from the 
			#readf_table with the sample id that we are carrying
		else{ 
			
			
			#I use getSampleConfiguration_locked just because I need the sample name
			my $params;
			my $main_name = $cfg_hash->{'db_readf_name'};
			getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_readf_table'},$cfg_hash->{'db_sample_id'},$sample_id);
					
			my $bam_input = $inFolder."/".get_output_name_from_executed($params,$cfg_hash->{'mark_rem_dup_step'},
						$params->{$main_name},\@steps_array).".".$cfg_hash->{'bam_ext'};		
								
			print_and_log("no! -> samplename: ".$params->{$main_name}.". Using input file: $bam_input\n",$log_file);
			if ( -e $bam_input) {
				push(@paths_to_write,$bam_input);
			}else{
				die "$analysis_indication ERROR: Cannot find file $bam_input that is needed for $task. Please run alignment step..\n";
			}
		}
	}
	#For DepthOfCoverage we can obtain a file with all the coverages per group
	#Now print a file with bam files paths from the different samples
	my $bam_paths_f = $outFolder."/$analysis_name\_bam_files_paths.".$cfg_hash->{'bamlist_ext'}; 
	#Open output from the annotation
	open (PATHS,">".$bam_paths_f) or die "$analysis_indication ERROR: Cannot open $bam_paths_f. The program will exit..\n";
	foreach my $path (@paths_to_write){
		print PATHS $path."\n";
	}
	close(PATHS);

	#########GET COVERAGE USING GATK DEPTHOFCOVERAGE
	my $task_name = 'coverage_analysis';
	my $covgenes_suff = $outFolder."/$analysis_name\_DOC";
	if ($cfg_hash->{$task_name} eq 'YES' ){
		if ( -e  $bam_paths_f){
			print_and_log( "$analysis_indication Starting the step to get $task_name : ",$log_file);
			run_GATK_DepthOfCoverage($cfg_hash,$bam_paths_f,$covgenes_suff,$analysis_id,$log_file);
			print_and_log( "$analysis_indication The output suffix is: $covgenes_suff..\n",$log_file);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $task_name $partDuration seconds\n",$log_file);
			$partTime = time;		
		}else{print_and_log( "$analysis_indication ERROR: Cannot get $task_name for analysis $analysis_name. File of paths $bam_paths_f does not exists.\n",$log_file);}
	}


	###Generate a table with gene coverage
	
	
	########################################
	##########FILLING OUTLIST FILE##########
	##########DepthOfCoverage files#####
	########################################	
	#Put the path to the coverage analysis files in the outlist file
	my @doc_outs = split(",",$cfg_hash->{'DOC_outputs'});
	foreach my $doc_out (@doc_outs){
		my $out_f = $covgenes_suff.".".$doc_out;
		if ( -e $out_f){
			#$doc_out =~ s/\_/ /g;
			
			print_and_log( "$analysis_indication Appending $out_f path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);	
			my $line = "-\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_cov_'.$doc_out}."\t".$out_f;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line);
		}else{
			print_and_log( "$analysis_indication WARNING: $out_f is not present. Please check what is happen. \n",$log_file);	
		}	
	}

}



=head2 run_fast_final_out_steps

 Title   : run_fast_final_out_steps
 Usage   : run_fast_final_out_steps(  config_file => the config hash
								);

 Function: Uses the final VCF:
						- extracts an input for Annovar from the VCF
						- installs Annovar databases if needed
						- runs Annovar to get variants annotated
						- writes a final output merging the information of the genotype from the VCF and the annotation from
							the output of Annovar
						- rearranges the final output
						- executes two R scripts to prioritize variants
						- prints coverage plots 
						- prints plots of statistics related with all same samples from the same target
						- gets the gene annotation
						- prints an HTML Web page to show results
						- sends an email to claim results are completed
 
					JOINT ANALYSIS
					  In case of a Joint analysis of more than 30 samples, VQSR will be executed	
						then we select from the VCF only those samples that we are analyzing in this run. 
						I check the genotype refinement step into the db
						
						Hence, the user may want to execute the pipeline until the VCF filtering and later
						execute the final output step for each of the families
 Returns :
 
=cut
sub run_fast_final_out_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	
	my $analysis_indication = "\n[an:$analysis_id ] (".scalar(localtime)."):";
	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	
	##################################################
	######Before to start take posses of the database db_busy=1
	print_and_log("$analysis_indication Setting analysis_id $analysis_id to be the last waiter of the db for the final output ...\n",$log_file);

	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_busy'},"1",$cfg_hash->{'db_info_row'},"1");	
	#################################################
	


	#Set this variable to get the last output
	my $previous_task = $cfg_hash->{'varfilt_step'};	
	my $step = $cfg_hash->{'finalout_step'};
	
	
	#Input and output folders
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	
	#Getting the sample name
	my $db_group_name = $cfg_hash->{'db_analysis_name'};
	

	#Getting the analysis name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
											$cfg_hash->{'db_analysis_id'},$analysis_id);
												
	my $params;
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $vcf_file = $inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$db_group_name},\@steps_array3).".".$cfg_hash->{'vcf_ext'};
	
	

	######CONSENSUS VARIANT CALLING
	###NB. CURRENTLY ONLY FREEBAYES IS LAUNCHED SINCE SAMTOOLS TAKES TOO MUCH TIME FOR THE EXECUTION
	#In this block of instructions we check if the user wants to execute the consensus approach
	#if YES VarGenius runs two jobs for Freebayes and Samtools-mpileup pipelines
	#and changes the dependencies of the fout task that will wait for these two jobs
	#Getting the consensus flag
	my $consensus = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_consensus'},
											$cfg_hash->{'db_analysis_id'},$analysis_id);
		
	#If this hash variable is set means that another calling program has been executed.
	if ( $consensus) {
		
		#Search if there are VCF files out of GATK
		my @vc_algorithms = split(",",$cfg_hash->{'vc_algorithms'});
		
		if (! defined $cfg_hash->{'vc_output_paths'}){
			$cfg_hash->{'vc_output_paths'} = "";
			foreach my $vc_algorithm (@vc_algorithms){
				my $vcfSupplOut	= "$inFolder/$analysis_name\_$vc_algorithm.".$cfg_hash->{'vcf_ext'};
				#If exists and size is non zero
				if (-e $vcfSupplOut ){
					if ( -s $vcfSupplOut){
						$cfg_hash->{'vc_output_paths'} .= $vcfSupplOut;
					}
				}
			}			
		}


		
		my $vc_merge_performed = 0;
		#Check the variant calling approaches executed which VCF file have returned
		if ( $cfg_hash->{'vc_output_paths'} ne ''){
			
			my @vc_output_paths = split(",",$cfg_hash->{'vc_output_paths'});
			
			#All the VCF files must have a field in the INFO indicating the algorithm used.
			#Thus I execute here a script that adds this field for each calling program
			my $final_vcf = "";
			my $vc_alg_num = 0;
			foreach my $vc_algorithm (@vc_algorithms){
				vcf_insert_string_into_field($vc_output_paths[$vc_alg_num],"VC_ALGORITHM=".$vc_algorithms[$vc_alg_num],uc($cfg_hash->{'db_vcf_info_str'}),$vc_output_paths[$vc_alg_num].".temp");
				$vc_alg_num++;
			}
			
			#Depending by user choice "consensus_vc_approach" here the VCF files will be
			#unified (UNION), intersected (INTERSECT) or unique variants in all approaches will be added to GATK result
			
			##UNION
			if ($cfg_hash->{'consensus_vc_approach'} eq 'UNION'){
				
				#The final VCF brings the name of all var call approaches
				$final_vcf = "$inFolder/$analysis_name"."_".$cfg_hash->{'consensus_vc_approach'}."_".join("_",@vc_algorithms).".".$cfg_hash->{'vcf_ext'};
				print_and_log("$analysis_indication The ".$cfg_hash->{'consensus_vc_approach'}." of VCF outputs will be performed and $final_vcf generated..\n",$log_file);	
				
				#Join the VCF files generated from all the variant calling software
				my $vc_alg_num = 0;
				my $first_vcf = $vcf_file;
				foreach my $vc_algorithm (@vc_algorithms){
					print_and_log("$analysis_indication Joining $first_vcf and ".$vc_output_paths[$vc_alg_num].".temp into $final_vcf.temp..\n",$log_file);			
					join_vcf_files_with_cat($first_vcf,$vc_output_paths[$vc_alg_num].".temp",$final_vcf.".temp");
					$first_vcf = $final_vcf.".temp";
					$vc_alg_num++;
				}
				$vc_merge_performed = 1;
				
				#After the join, we need to sort the VCF file and remove all the duplicate calls (where all samples have the same
				#genotype)
				print_and_log("$analysis_indication Sorting $final_vcf.temp into $final_vcf.sort.temp..\n",$log_file);			
				run_VCFTOOLS_sort($cfg_hash,$final_vcf.".temp",$final_vcf.".sort.temp","",$log_file);
				print_and_log("$analysis_indication Removing identical calls from into $final_vcf.sort.temp and saving into $final_vcf..\n",$log_file);
				vcf_filter_identical_calls($cfg_hash,$final_vcf.".sort.temp",$final_vcf,$log_file);
				
				$vcf_file = $final_vcf;

				#Delete all temporary files
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					print_and_log("$analysis_indication Removing the following temporary files: $final_vcf.temp..\n",$log_file);	
					delete_file($final_vcf.".temp");
					delete_file($final_vcf.".sort.temp");
				}		
				
			}
			
			##Addition of Unique variants to GATK result
			elsif ($cfg_hash->{'consensus_vc_approach'} eq 'ADDUNIQ'){
				
				#The final VCF has the name of all approaches
				$final_vcf = "$inFolder/$analysis_name"."_".$cfg_hash->{'consensus_vc_approach'}."_".join("_",@vc_algorithms);
				
				#Join the two VCF files
				my $vc_alg_num = 0;
				my $first_vcf = $vcf_file;
				foreach my $vc_algorithm (@vc_algorithms){
					
					my $params  ="";
		
					#I have to get the variants from the common that are not in the results from  GATK
					my $current_vcf = $vc_output_paths[$vc_alg_num].".temp";
					print_and_log("$analysis_indication Get the variants from $current_vcf that are not in GATK $vcf_file..\n",$log_file);	
					my $vcf_uniq = "$inFolder/$analysis_name\_$vc_algorithm\_uniq.vcf";
					run_BEDTOOLS_subtractBed($cfg_hash,$current_vcf,$vcf_file,$vcf_uniq,$params,$log_file);
					
					#The output of subtractBed is without the VCF header. Hence I just concatenate the results and then order
					#Finally I combine the results of GATK with those that are new to GATK just extracted	
					#Join the two VCF files
					print_and_log("$analysis_indication Joining $first_vcf and $vcf_uniq into $final_vcf.temp..\n",$log_file);			
				
					join_vcf_files_with_cat($first_vcf,$vcf_uniq,$final_vcf.".temp");
					$first_vcf = $final_vcf.".temp";
					$vc_alg_num++;
					
					#Delete all temporary files
					if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
						print_and_log("$analysis_indication Removing the following temporary files: $vcf_uniq..\n",$log_file);	
						delete_file($vcf_uniq);
					}	
				}
				$vc_merge_performed = 1;
				
				#Sort the obtained VCF and it will be the VCF to annotate
				print_and_log("$analysis_indication Sorting $final_vcf.temp into $final_vcf..\n",$log_file);			
				run_VCFTOOLS_sort($cfg_hash,$final_vcf.".temp",$final_vcf,"",$log_file);			
				$vcf_file = $final_vcf;

				#Delete all temporary files
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					print_and_log("$analysis_indication Removing the following temporary files: $final_vcf.temp..\n",$log_file);	
					delete_file($final_vcf.".temp");
				}						
				
			}
			#This is the real consensus approach because choses all the calls only that are called from all the callers
			elsif ($cfg_hash->{'consensus_vc_approach'} eq 'INTERSECT'){
					######WRITE CODE
			}else{
				print_and_log("$analysis_indication No consensus approach has been used (".$cfg_hash->{'consensus_vc_approach'}."). The result will be only from GATK.\n",$log_file);
			}
						
		}
	}
	###########################

	#######CONSENSUS VARIANT CALLING
	####NB. CURRENTLY ONLY FREEBAYES IS LAUNCHED SINCE SAMTOOLS TAKES TOO MUCH TIME FOR THE EXECUTION
	##In this block of instructions we check if the user wants to execute the consensus approach
	##if YES VarGenius runs two jobs for Freebayes and Samtools-mpileup pipelines
	##and changes the dependencies of the fout task that will wait for these two jobs
	##Getting the consensus flag
	#my $consensus = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											#$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_consensus'},
											#$cfg_hash->{'db_analysis_id'},$analysis_id);
											
	##When using more variant calling methods I use an hash that for each compid (variant) 
	##records which program obtained it. The GATK unique result is not recorded
	#my @fb_vars = ();
	#if ( $consensus ){
		
		##In this snip of code I am adding the variants extracted with Freebayes and Samtools to 
		##those extracted by GATK. In particular:
		#my $fb_vcf_output = "$inFolder/$analysis_name\_freebayes.vcf";
		##my $sam_vcf_output = "$inFolder/$analysis_name\_mpileup.vcf.gz";
		#my $params  ="";
		
		###First the variants from Freebayes and mpileup are intersected to pick the common
		##my $vcf_intersect = $fb_vcf_output;#"$inFolder/$analysis_name\_fb-sam.vcf.gz";
		##run_BEDTOOLS_intersectBed($cfg_hash,$fb_vcf_output,$sam_vcf_output,$vcf_intersect,$params,$log_file);#CON SAMTOOLS
		
		##Then I have to get the variants from the common that are not in the results from  GATK
		#print_and_log("$analysis_indication Get the variants from $fb_vcf_output that are not in GATK $vcf_file..\n",$log_file);	
		#my $vcf_fbuniq = "$inFolder/$analysis_name\_fbuniq.vcf";
		#run_BEDTOOLS_subtractBed($cfg_hash,$fb_vcf_output,$vcf_file,$vcf_fbuniq,$params,$log_file);
		
		##The output of subtractBed is without the VCF header. Hence I just concatenate the results and then order
		##Finally I combine the results of GATK with those that are new to GATK just extracted
		#my $vcfs = "$vcf_file,$vcf_fbuniq";
		#my $final_vcf = "$inFolder/$analysis_name\_gatk_fb.vcf";		
		##Join the two VCF files
		#print_and_log("$analysis_indication Joining $vcf_file and $vcf_fbuniq into $final_vcf.temp..\n",$log_file);			
		#join_vcf_files_with_cat($vcf_file,$fb_vcf_output,$final_vcf.".temp");
		##join_files_with_cat($vcf_file,$vcf_fbuniq,$final_vcf.".temp");
		##Sort the obtained VCF
		#print_and_log("$analysis_indication Sorting $final_vcf.temp into $final_vcf..\n",$log_file);			
		#run_VCFTOOLS_sort($cfg_hash,$final_vcf.".temp",$final_vcf,"",$log_file);
		##...and the final VCF name changes
		#$vcf_file = $final_vcf;
		
		###Furthermore I want to take the variants here that are uniq to GATK, just to know which one has been added
		##my $vcf_gatk_uniq = "$inFolder/$analysis_name\_gatk_uniq.vcf.gz";
		##print_and_log("$analysis_indication Getting uniq to GATK into $vcf_gatk_uniq..\n",$log_file);	
		##run_BEDTOOLS_subtractBed($cfg_hash,$vcf_file,$vcf_fbuniq,$vcf_gatk_uniq,$params,$log_file);

		##Register the variants uniq to freebayes into an array
		#print_and_log("$analysis_indication Recording the variants uniq to freebayes into an array..\n",$log_file);			
	  #extract_columns_from_file($vcf_fbuniq,"0,1,3,4",$vcf_fbuniq.".temp");
	  #separate_elements_on_col($vcf_fbuniq.".temp",$vcf_fbuniq.".sepvar.temp",4);
	  #open(FB_V ,"<$vcf_fbuniq.sepvar.temp") or die "ERROR: cannot open $vcf_fbuniq.sepvar.temp\n";
	  #while ( my $row = <FB_V> ) {
			#chop($row);
			#my @fields = split("\t",$row); 
			#my $compid = join("_",@fields);
			#push(@fb_vars,$compid);
		#}
		#close(FB_V);
		#print_and_log("$analysis_indication Elements of the array:\n",$log_file);	
		#print_array(\@fb_vars);
		
		##Delete all temporary files
		#if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
			#print_and_log("$analysis_indication Removing the following temporary files: $vcf_fbuniq\n$final_vcf.temp..\n",$log_file);	
			#delete_file($vcf_fbuniq);
			#delete_file($final_vcf.".temp");
			#delete_file($vcf_fbuniq.".temp");
			#delete_file($vcf_fbuniq.".sepvar.temp");
		#}				
	#}
	#######CONSENSUS VARIANT CALLING
	
	
	#File name for the VCF with variants that must be annotated (missing)
	my $vcf_to_ann = $inFolder."/".extract_name($vcf_file,1).".".$cfg_hash->{'miss_ann_ext'};
	#File VCF with all variants
	my $vcf_all_var = $inFolder."/".extract_name($vcf_file,1).".".$cfg_hash->{'allvar_ext'};
			
																								
	#VCF import. Reads the  variants from the VCF file and generates a VCF with all the missing	variants and another with all the variants
	if ($cfg_hash->{'create_vcf_to_ann'} eq 'YES'){
		
	#################VQSR Considerations
  #If VQSR was executed	then we select from the VCF only those samples that we are
  #analyzing in this run. I check the genotype refinement step into the db
 	my $genotref = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_genotref'},
									$cfg_hash->{'db_analysis_id'},$analysis_id); 
	#print_and_log( "genotref = $genotref and runsamples: ".$cfg_hash->{'run_samples'}." \n",$log_file);									
  if ($genotref == 1){
		if ($cfg_hash->{'run_samples'} ne 'ALL'){

			my $sep = ",";
			my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},$sep);
			#print "Obtained: $samples_ids\n ";#DEBUGCODE
			my @samples_ids = split($sep,$samples_ids);
			my $selvar_param_str = "";
			foreach my $sel_sample_id (@samples_ids){
					my $sel_sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
									$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'},"$analysis_id,$sel_sample_id"); 
				$selvar_param_str .= " --sample_name $sel_sample_name ";		
			}
			my $temp_out = $vcf_file.".temp";
			my $prog_used = $cfg_hash->{'selvar_prog'};
			print_and_log( "$analysis_indication Since the VQSR has been executed we use a VCF with only selected samples:\n ".
			"Executing $prog_used with $selvar_param_str\n",$log_file);
			#If the file was created before means that the vcf_file has no more all the samples
			if ( -e $vcf_file.".all"){
				run_GATK_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$vcf_file.".all",$temp_out,$log_file);				
			}else{
				run_GATK_SelectVariants($cfg_hash,$prog_used,$analysis_id,$selvar_param_str,$vcf_file,$temp_out,$log_file);				
				#Change the name of the complete VCF file adding .all and the .temp is the new VCF file
				print_and_log( "$analysis_indication Moving $vcf_file to $vcf_file.all \n",$log_file);			
				move($vcf_file, $vcf_file.".all") or die "$analysis_indication ERROR: Cannot move $vcf_file in $vcf_file.all\n";
			}
			print_and_log( "$analysis_indication Moving $temp_out in $vcf_file \n",$log_file);	
			move($temp_out, $vcf_file) or die "$analysis_indication ERROR: Cannot move $temp_out in $vcf_file\n";      
		}
		#Otherwise, if the samples must all be used, then just check if the .all file is present.
		#In that case remove the .all before to proceed
		else{
			if ( -e $vcf_file.".all"){
				move($vcf_file.".all",$vcf_file) or die "ERROR: Cannot move $vcf_file.all in $vcf_file\n";	
			}
		}
	}
	#################END VQSR Considerations	
	
	
		my $prog_used = "VarGenius_create_vcf_to_ann";
		if ( -e $vcf_file){
			
			#Import if the same analysis has not been inserted yet
			my $step = $cfg_hash->{'vcfimport_step'};
			print_and_log( "$analysis_indication Executing create_vcf_to_ann. Two tables will be generated from $vcf_file".
			" $vcf_all_var including variants info and $vcf_all_var.info for genotype information\n",$log_file);
			generate_vcf_to_annotate($cfg_hash,$prog_used,$analysis_id,$log_file,$vcf_file,$vcf_all_var);
			#Save status in the db by searching for the specific sample name
			#print_and_log("$analysis_indication Updating analysis status for group $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used $partDuration seconds\n",$log_file);
			$partTime = time;	
		}else{
				log_and_exit("$analysis_indication ERROR: Cannot start $prog_used for analysis $sample_id. Missing the VCF output: $vcf_file\n",$log_file);
		}
	}
	

	#Check existence of Annovar databases	
	if ($cfg_hash->{'annov_db_dl'} eq 'YES'){
		my $prog_used = $cfg_hash->{'annot_var_prog'};
		print_and_log("$analysis_indication Checking the existence of Annovar database...\n",$log_file);
		#Get the annotation types wanted
		my @ann_types = split($cfg_hash->{'word_sep'},$cfg_hash->{'annov_ann_types'});
		#For each annotation type get its databases
		foreach my $ann_type (@ann_types){
			if ( defined $cfg_hash->{'annov_'.$ann_type.'_dbs'}) {
				my @databases = split($cfg_hash->{'word_sep'},$cfg_hash->{'annov_'.$ann_type.'_dbs'});
				print_and_log("$analysis_indication Installing ".scalar(@databases)." databases for Annovar ($ann_type) using $prog_used \n",$log_file);
				print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);	
				#Download each database if it does not exists yet
				foreach my $db_2_dl (@databases){
						print_and_log( "$analysis_indication Checking database: $db_2_dl\n",$log_file);#DEBUGCODE
						#Obtain the database id if the database is already installed
						#picks from the annovar_path the last part of the path (annovar db version) and uses it for the query
						my $annovar_ver = extract_name($cfg_hash->{'annovar_db_f'},0);
						my $db_present = "";
						#print_and_log( "$analysis_indication Searching $db_2_dl and $annovar_ver into database: $db_2_dl\n",$log_file);#DEBUGCODE
						$db_present = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
											$cfg_hash->{'db_pass'},$cfg_hash->{'db_annovar_dbs_table'},$cfg_hash->{'db_annovar_dbs_name'},
											$cfg_hash->{'db_annovar_dbs_name'}.",".$cfg_hash->{'db_annovar_dbs_ver'},"'".$db_2_dl."','".$annovar_ver."'");
						#Means the return is -1, db does not exist
						#Otherwise the return is the name of the database
						if ( length $db_present <= 2 ){
							#The suffix _EXT next to the database names, indicates that
							#the db type is generic. The suffix _BED means that the db is a bed file.
							#In those cases the user needs to download
							#it manually and put in path written at annovar_db_f (program_config)
							if( $db_2_dl !~ /\_EXT$/  and $db_2_dl !~ /\_BED$/){							
								#For GenCode database we have a different procedure
								if( $db_2_dl =~ /Gencode/ ){
									my $fasta_seqs = "seq";
									run_ANNOVAR_download_WG_fasta($cfg_hash,$prog_used,$log_file,$fasta_seqs);					
									print_and_log( "$analysis_indication Installing database $db_2_dl for Annovar : ",$log_file);
									run_ANNOVAR_download_db($cfg_hash,$prog_used,$log_file,$db_2_dl,"r",'null');
									run_ANNOVAR_retrieveseq($cfg_hash,$prog_used,$log_file,$db_2_dl);
									
								}else{
									my $source;
									#Save the db2dl as it is because if it is ucsc I will drop the 'ucsc'
									my $temp = $db_2_dl;
									#If it is not EXT and not Gencode, check if user wants to download from UCSC or from Annovar
									if( $db_2_dl =~ /\_UCSC$/ ){
										$source = "UCSC";
										$db_2_dl =~ s/\_UCSC//;
									}else{
										$source = "ANN";
									}
									print_and_log( "$analysis_indication Installing database $db_2_dl for Annovar : ",$log_file);
									run_ANNOVAR_download_db($cfg_hash,$prog_used,$log_file,$db_2_dl,$ann_type,$source);
									#restore db2dl
									$db_2_dl = $temp;
								}		
							}else{
								print_and_log( "$analysis_indication Database $db_2_dl is EXT or BED \n",$log_file);#DEBUGCODE
								my $path = "";
								my $db_name = $db_2_dl;
								if ($db_name =~ /\_EXT$/){
									checkIncludedAnnovarDBs($cfg_hash,$db_name,$cfg_hash->{'txt_ext'},$log_file);
									#Check if the EXT database exists in the folder. The name must be HGVERSION_DBNAME.txt
									$db_name =~ s/\_EXT//;
									$path = $cfg_hash->{'annovar_db_f'}."/".$cfg_hash->{'ann_build_ver'}."_".$db_name.".".$cfg_hash->{'txt_ext'};
								}
								if ($db_name =~ /\_BED$/){
									checkIncludedAnnovarDBs($cfg_hash,$db_name,$cfg_hash->{'bed_ext'},$log_file);
									#Check if the BED database exists in the folder. The name must be HGVERSION_DBNAME.bed
									$db_name =~ s/\_BED//;
									$path = $cfg_hash->{'annovar_db_f'}."/".$cfg_hash->{'ann_build_ver'}."_".$db_name.".".$cfg_hash->{'bed_ext'};
								}
								#my $path = $cfg_hash->{'annovar_db_f'}."/".$cfg_hash->{$db_2_dl.'_f'};
								
								if ( ! -e  $path ){
									log_and_exit( "$analysis_indication WARNING: Database $path for Annovar does not exist in ".$cfg_hash->{'annovar_db_f'}.
									" To use it you must download from Annovar website and insert there. Now it won't be used. \n",$log_file);
								}
							}
							
							#Save in the db the specific database name
							my $fields = " ".$cfg_hash->{'db_annovar_dbs_name'}.",".$cfg_hash->{'db_annovar_dbs_type'}.",".$cfg_hash->{'db_annovar_dbs_ver'};
							my $values = " '".$db_2_dl."','".$ann_type."','".$annovar_ver."'";
							print_and_log("$analysis_indication Updating table with database $db_2_dl in ".$cfg_hash-> {'db_name'}." ...\n",$log_file);
							insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
									$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_annovar_dbs_table'},
									$fields,$values);
						}else{
							print_and_log( "$analysis_indication Database $db_2_dl for Annovar already exists... \n",$log_file);
						}	
				}#code
			}	
		}

		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
		$partTime = time;	
	}
	


	###################################################
	############### COMPLETE ANNOTATION 		###########
	###################################################
	
	#Since the import-export operation using the database is very computationally intensive
	#I decided (Sept 2016) to implement the annotation of the complete VCF
	my $allvar_out_suff = $outFolder."/".$params->{$db_group_name}."_".$cfg_hash->{'allvar_ext'};
	if ($cfg_hash->{'annot_all'} eq 'YES'){

		my $prog_used = $cfg_hash->{'table_annov_prog'};
		print_and_log("$analysis_indication Annotate all variants in VCF file $vcf_all_var ...\n",$log_file);
		print_and_log("$analysis_indication  RUNNING $prog_used : ",$log_file);	
		#inFolder for the VCF files
		my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
		my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
		
		my $params;
		getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
		#Getting the sample name
		my $db_group_name = $cfg_hash->{'db_analysis_name'}; 
			
		my $in_annov = $outFolder."/".extract_name($vcf_file,1)."_".$cfg_hash->{'allvar_ext'}.".".$cfg_hash->{'annov_in_ext'};
		
		if (file_num_rows($vcf_all_var) > 1) {
			#Convert VCF format to ANNOVAR format					
			print_and_log("$analysis_indication Convert VCF format file $vcf_all_var to ANNOVAR format in $in_annov : ",$log_file);									
			run_ANNOVAR_convert2annovar($cfg_hash,$prog_used ,$log_file,$vcf_all_var,$in_annov);
			#Annotate all with all the databases and settings from the user configuration
			print_and_log("$analysis_indication Annotating $in_annov : ",$log_file);			
			run_ANNOVAR_table_annovar($cfg_hash,$prog_used,$log_file,$in_annov,$allvar_out_suff);
			print_and_log("$analysis_indication Completed the annotation of $in_annov. Output file is: $allvar_out_suff...\n",$log_file);
			
			#Delete the input for Annovar
			if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
				print_and_log("$analysis_indication Removing  $in_annov..\n",$log_file);	
				delete_file($in_annov);
				#Delete the invalid_input files
				delete_file("*.invalid_input");
			}
			
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
			$partTime = time;	
		}else{
			log_and_exit("$analysis_indication WARNING: Cannot annotate  $in_annov. File empty ...\n",$log_file);
		}	
		
	}

	#Export annotation from the database using fields from the
	#configuration file
	my $compl_annovar_out =  $allvar_out_suff.".".$cfg_hash->{'ann_build_ver'}."_".$cfg_hash->{'annovar_out_file'};#Annovar Output	
	my $vargenius_out = $allvar_out_suff.".".$cfg_hash->{'ann_build_ver'}.".final_out.txt";#First VarGenius Output
	#The file with INFO and FORMAT fields from the VCF. The format information is concatenated for each sample
	my $other_info_f = $vcf_all_var.".info";
	if ($cfg_hash->{'write_output'} eq 'YES'){
		#Both the VCF and the annotated vcf must exist to run
		if ( -e $compl_annovar_out and !(-z $compl_annovar_out) ){
			if ( -e $vcf_file  and   !(-z $vcf_file) ){
				my $prog_used = "VarGenius-write_output";	
				print_and_log( "$analysis_indication Extracting the final output using the annotation ($compl_annovar_out), $vcf_all_var and $vcf_all_var.info into ($vargenius_out) ..\n",$log_file);
				print_and_log("$analysis_indication  RUNNING $prog_used : ",$log_file);	
				output_from_annotation_without_db($cfg_hash,$prog_used,$analysis_id,$log_file,$vargenius_out,$compl_annovar_out,$other_info_f,$vcf_file);
				
				###############################################
				#To allow the use of the rearrange_output function. Here I consider that the 
				#func column is before the gene name column usually in the annotation.
				#while the function considers this two columns inverted. Hence I invert them here
				my @func_cols = split(",",$cfg_hash->{'func_cols'});
				foreach my $func_col (@func_cols){
					my $func_col_ind = get_col_index($vargenius_out,$func_col);
					#print_and_log("$analysis_indication Inverting columns $func_col_ind and ".($func_col_ind+1)." to shift $func_col ",$log_file);#DEBUGCODE
					#print_and_log("(This operation is required for rearrange_output)\n",$log_file);#DEBUGCODE
					invert_cols_position($vargenius_out,$func_col_ind+1,$func_col_ind+2);					
				}

				###############################################
				
				################################################
				#Shift the interval frequency fields after the 1000genomes frequency
				#Get the position of the ivf field to remove it later
				my $ivf_position = get_col_index($vargenius_out,$cfg_hash->{'intern_var_freq_fld'});
				
				#print_and_log( "$analysis_indication Shifting columns ".$cfg_hash->{'intern_var_freq_fld'}." and ".$cfg_hash->{'var_frequency_factors_fld'}." after ".$cfg_hash->{'last_1000g_field'}."\n",$log_file);#DEBUGCODE
				#Set the file names for the two columns which will be extracted from the file
				my $ivf_col_f = $vargenius_out.".temp_ivf";
				my $ffact_col_f = $vargenius_out.".temp_ffact";
				#Extract the columns into the two files
				#print_and_log( "Extracting column ".$cfg_hash->{'intern_var_freq_fld'}." from $vargenius_out to $ivf_col_f\n",$log_file);#DEBUGCODE
				#print_and_log( "Extracting column ".$cfg_hash->{'var_frequency_factors_fld'}." from $vargenius_out to $ffact_col_f\n",$log_file);#DEBUGCODE
				extract_col_from_file($vargenius_out,$cfg_hash->{'intern_var_freq_fld'},$ivf_col_f);
				extract_col_from_file($vargenius_out,$cfg_hash->{'var_frequency_factors_fld'},$ffact_col_f);
				#Set the two arrays that will be the columns to insert and insert them
				my @ivf_list = list_to_array($ivf_col_f,'NO_NEW_LINE');
				my @ffact_list = list_to_array($ffact_col_f,'NO_NEW_LINE');
				#	push the header into the lists
				unshift @ivf_list, $cfg_hash->{'intern_var_freq_fld'};
				unshift @ffact_list, $cfg_hash->{'var_frequency_factors_fld'};
				
				#Choose the point where to insert and insert it
				#print_and_log( "Inserting column ".$cfg_hash->{'var_frequency_factors_fld'}." from file $ivf_col_f into $vargenius_out \n",$log_file);
				my $last1000g_info_ind = get_col_index($vargenius_out,$cfg_hash->{'last_1000g_field'});
				#print_and_log( "ivf_list array has ".scalar(@ivf_list)." elements and the first is ".$ivf_list[0]."\n",$log_file);#DEBUGCODE
				#print_and_log( "ffact_list array has ".scalar(@ffact_list)." elements and the first is ".$ffact_list[0]."\n",$log_file);#DEBUGCODE
				insert_col_in_file_table($vargenius_out,\@ivf_list,$cfg_hash->{'last_1000g_field'},'i');
				insert_col_in_file_table($vargenius_out,\@ffact_list,$cfg_hash->{'intern_var_freq_fld'},'i');
				#Delete the previous columns using the previously obtained index $ivf_position
				#If you put them before the previous position, now that position needs you sum 2
				my @cols_to_remove_later = ();
				if ( $last1000g_info_ind < $ivf_position){
					@cols_to_remove_later = (($ivf_position+2),($ivf_position+3));
				}else{
					@cols_to_remove_later = (($ivf_position),($ivf_position+1));
				}
				#print_and_log( "Deleting ".scalar(@cols_to_remove_later)." [".$cols_to_remove_later[0]." and ".$cols_to_remove_later[1]."] columns from $vargenius_out\n",$log_file);#DEBUGCODE
				delete_columns($vargenius_out,\@cols_to_remove_later);
	
				#Delete the two column-files
				if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
					print_and_log("$analysis_indication Removing $ivf_col_f and $ffact_col_f..\n",$log_file);	
					delete_file($ivf_col_f);
					delete_file($ffact_col_f);
				}			
				#############################################################
				

				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
				$partTime = time;		
			}else{
				log_and_exit("$analysis_indication WARNING: Cannot write the output $vcf_file is missing or empty ...\n",$log_file);
			}	
		}else{
			log_and_exit("$analysis_indication WARNING: Cannot write the output $compl_annovar_out is missing or empty..\n",$log_file);
		}	

	}	
	
	
	#Here I change the columns disposition. This change is based on the use of the
	#profiles and is specifically performed when a research group asks to do it
	#hence we need to get the research group before so that the output is
	#associated with their identifier
	
	#The quick output writing (write_output) has priority
	my $output_2_rearrange = $vargenius_out;
	
	#Research group specific requests
	my $resgroupid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_resgroupid'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
	
	my $rearr_ann_out =  $outFolder."/".$params->{$db_group_name}.".".$cfg_hash->{'ann_build_ver'}.".final_out_resgr$resgroupid.txt";
	my $rearr_ann_out_xls =  $outFolder."/".$params->{$db_group_name}.".".$cfg_hash->{'ann_build_ver'}.".final_out_resgr$resgroupid.xlsx";
	if ($cfg_hash->{'rearrange_output'} eq 'YES'){
		my $prog_used = "VarGenius-rearrange_output";	
		if ( -e $output_2_rearrange and !(-z $output_2_rearrange)) {
			if ( correct_type($resgroupid,"positiveint")){							
				
				print_and_log("$analysis_indication  RUNNING $prog_used rearrange_output..\n",$log_file);	
				#print_and_log( "Request from research group: $resgroupid\n",$log_file);
				#configuration file#I seguenti tre commenti possono essere elimnati
				#my $temp_file = $outFolder."/".$params->{$db_group_name}.".".$cfg_hash->{'ann_build_ver'}.".final_out_resgr$resgroupid.temp";
				#print_and_log( "Repositioning columns as from output_new_cols variable...\n",$log_file);
				#extract_columns_from_file($output_2_rearrange,$cfg_hash->{'output_new_cols'},$temp_file);
				
				print_and_log( "$analysis_indication Rearranging the output in $output_2_rearrange making $rearr_ann_out...\n",$log_file);
				rearrange_rawtabular_out($cfg_hash,$output_2_rearrange,$rearr_ann_out,$log_file);
				
				#rearrange_rawtabular_out($cfg_hash,$output_2_rearrange,$rearr_ann_out,\@fb_vars,$log_file);
				
				#Convert the tsv file to xls in the same location
				print_and_log( "$analysis_indication Converting the TSV file ($rearr_ann_out) to XLS format ($rearr_ann_out_xls)...\n",$log_file);
				tsv_2_xls($rearr_ann_out,$rearr_ann_out_xls);
				
			
				$partDuration = time - $partTime;
				print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
				$partTime = time;		
			}else{
				log_and_exit( "$analysis_indication ERROR: Cannot execute rearrange_output for research group: $resgroupid\n",$log_file);
			}
		}else{
			log_and_exit("$analysis_indication WARNING: Cannot find $output_2_rearrange. $prog_used will not be executed...\n",$log_file);
		}
	}

				
	#If Joint Analysis is used then create different outputs from the final
	if ($cfg_hash->{'sep_joint_analyses'} eq 'YES'){
		print_and_log("$analysis_indication   RUNNING sep_joint_analyses..\n",$log_file);	
		#Separate the joint analyses in the unique file
		my @out_paths = ();
		separate_joint_outputs($cfg_hash,$rearr_ann_out,\@out_paths,$log_file);
		
		#Now for each path convert in XLS
		foreach my $out_path (@out_paths){
			#Convert the tsv file to xls in the same location
			my $out_path_xls = extract_name($out_path,'noext').".xlsx";
			print_and_log( "$analysis_indication Converting the TSV file ($out_path) to XLS format ($out_path_xls)...\n",$log_file);
			tsv_2_xls($out_path,$out_path_xls);				

			#Put the path to the Final tabular output in the outlist file
			if ( check_presence($out_path) ){
				print_and_log( "$analysis_indication Appending $out_path path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
				my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_out'}." ".extract_name($out_path,0)."\t".$out_path;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
			}else{
				print_and_log( "$analysis_indication WARNING: $rearr_ann_out is not present. Please check what is happen. \n",$log_file);	
			}		
			
			#Put the path to the XLS version of the Final tabular output in the outlist file
			if ( check_presence($out_path_xls) ){
				print_and_log( "$analysis_indication Appending $out_path_xls path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
				my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_xls'}." ".extract_name($out_path_xls,0)."\t".$out_path_xls;
				append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
			}else{
				print_and_log( "$analysis_indication WARNING: $rearr_ann_out_xls is not present. Please check what is happen. \n",$log_file);	
			}	
		}
	}	


	#If Joint Analysis is used and samples are selected, then fetch the needed samples from the output
	#and generate a new final output that will be used through the rest of the code
	if ($cfg_hash->{'separate_analysis_from_joint'} eq 'YES'){
		
		print_and_log("$analysis_indication  RUNNING separate_analysis_from_joint..\n",$log_file);	
		
		#Separate the joint analyses in the unique file
		my $sep = ",";
		my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},$sep);
		my @samples_ids = split($sep,$samples_ids);
		my $sample_names = "";
		#print_and_log( "$analysis_indication Picking informations only for $sample_names..\n",$log_file);#DEBUGCODE
		foreach my $sel_sample_id (@samples_ids){
			#print_and_log( "Picking informations only for $sample_names [$sel_sample_id]..\n",$log_file);#DEBUGCODE
			my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
							$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_id'},"$analysis_id,$sel_sample_id"); 
			$sample_names .= $sample_name.",";		
		}
		chop($sample_names);
		
		#Define output names
		my $out_path = extract_name($rearr_ann_out,'noext')."_".$sample_names.".".$cfg_hash->{'txt_ext'};
		#substitute commas with underscores
		$out_path =~ s/,/_/g;
		my $out_path_xls = extract_name($out_path,'noext').".".$cfg_hash->{'xlsx_ext'};
		
		print_and_log( "$analysis_indication Writing information only for $sample_names into file $out_path...\n",$log_file);
		separate_analysis_from_joint($cfg_hash,$cfg_hash->{'run_samples'},$analysis_id,$rearr_ann_out,$out_path,$log_file);
		
		if ( check_presence($out_path) ){
			#The output is now separated for only the samples in run_samples but still it contains all the variants
			#coming out from the joint analysis. We need to filter it removing those which are not variants at all for the 
			#family. Hence we use th function filter_no_variant which removes all those cases
			#where must happens that all the samples of the family have GT=0/0 or 0|0 (R script)
			print_and_log( "$analysis_indication Filtering variants for $sample_names (remove GT=0/0 or 0|0 for all samples)...\n",$log_file);
			my $RLogPath = $outFolder."/log/$analysis_name\_filter_no_variants".$cfg_hash->{'R_log_file'};
			my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_utils'};
					
			my $command = $cfg_hash->{'R2_path'}." CMD BATCH --no-save --no-restore '--args FILTER_NO_VARIANTS $out_path $out_path.temp'  $R_script $RLogPath";
			print_and_log("$analysis_indication The command is: $command\n",$log_file);
			try_exec_command($command) or R_die($command,$R_script);			
			
			#Change the temp filename and remove it
			move("$out_path.temp",$out_path);
			if ( $cfg_hash->{'remove_temp'} eq 'YES'){
					delete_file("$out_path.temp");
			}	
			
			#Now convert in XLSX
			#Convert the tsv file to xls in the same location
			print_and_log( "$analysis_indication Converting the TSV file ($out_path) to XLS format ($out_path_xls)...\n",$log_file);
			tsv_2_xls($out_path,$out_path_xls);				
		}else{
			print_and_log( "$analysis_indication WARNING: $out_path is not present. Please check what is happen. \n",$log_file);	
		}
		
		#Put the path to the Final tabular output in the outlist file
		if ( check_presence($out_path) ){
			print_and_log( "$analysis_indication Appending $out_path path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
			my $line = "$sample_names\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_out'}." ".extract_name($out_path,0)."\t".$out_path;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
		}else{
			print_and_log( "$analysis_indication WARNING: $out_path is not present. Please check what is happen. \n",$log_file);	
		}		
		
		#Put the path to the XLS version of the Final tabular output in the outlist file
		if ( check_presence($out_path_xls) ){
			print_and_log( "$analysis_indication Appending $out_path_xls path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
			my $line = "$sample_names\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_xls'}." ".extract_name($out_path_xls,0)."\t".$out_path_xls;
			append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
		}else{
			print_and_log( "$analysis_indication WARNING: $out_path_xls is not present. Please check what is happen. \n",$log_file);	
		}	
		
		#Now change the name of the rearranged output to this that has just been changed
		$rearr_ann_out = $out_path;
		
	}	
	
	############################################
	################PRIORITIZATION
	#Supplementary annotation wrote by Margherita and Micheles Prioritization
	#############################################
	
	my $out2nd_lev = extract_name($rearr_ann_out,1)."_2ndlev";
	if ($cfg_hash->{'supplementary_output'} eq 'YES'){
		my $prog_used = "VarGenius-supplementary_output";	
		if ( -e $rearr_ann_out and !(-z $rearr_ann_out) ){
			#Working directory
			my $priorwdir = "$outFolder/prioritization";
			
			#The script needs in input the list of samples in a correct order (probands,mother, father)
			#so we need here to now if the analysis is a joint one or it is a simple single, Trio or Quartet analysis
			#we can use the parameter $cfg_hash->{'run_samples'}
			###############Variables for SORT SAMPLES
			my @samples_names = ();#The samples names

			if ($cfg_hash->{'run_samples'} ne 'ALL') {
				
				my $sep = ",";
				my $samples_ids = separate_input_ids($cfg_hash->{'run_samples'},$sep);
				#print "Obtained: $samples_ids\n ";#DEBUGCODE		
				my @samples_ids = split($sep,$samples_ids);
				
				#Now get the names
				foreach my $sample_id (@samples_ids){
					my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
													$cfg_hash->{'db_sample_id'}.",".$cfg_hash->{'db_analysis_id'},$sample_id.",".$analysis_id);		
					push(@samples_names,$sample_name);		
				}
			}else{
				#Get the samples names involved for the group
				my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
				#print "Executing: $query\n";#DEBUGCODE
				my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
				foreach my $sample_name (keys %{$group_sam}){ 
					push(@samples_names,$sample_name);
				}	
			}
			###############SORT SAMPLES
			my $samples_h;#Samples hash to be reordered
			my @sort_samples = ();#The array with samples id reordered
			my $kinship_ok = 1;			
			#Get the kinship, to make the resorting
			foreach my $sample_name (@samples_names){
				
				
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
				#Get the order from the db
				@sort_samples = @samples_names;
			}
			
			#print_and_log("$analysis_indication Samples header: " ,$log_file);#DEBUGCODE
			#foreach my $sample (@sort_samples){	
				#print_and_log(" $sample " ,$log_file);#DEBUGCODE
			#}print_and_log("\n" ,$log_file);#DEBUGCODE
			############END SORT
			
			#Put samples sorted into a string separated by space
			my $sort_samples = join(" ",@sort_samples);
			#Change output names if samples to analyze are restricted
			if ($cfg_hash->{'run_samples'} ne 'ALL') {
				#Change also supplementary ouytput name and..
				$out2nd_lev .= "_".join("_",@sort_samples);
				#Change working dir for prioritization name
				$priorwdir .= "_".join("_",@sort_samples);
			}		

			#Getting the pedigree file 
			my $ped_file = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_pedfile'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);		
			#Getting the pedigree file 
			my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						    $cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
						    $cfg_hash->{'db_analysis_id'},$analysis_id);								    
			

			##################
			#Supplementary Second level output (Margherita Mutarelli)
			#################
			print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
			my $RLogPath = $outFolder."/log/$analysis_name\_supplementary_output".$cfg_hash->{'R_log_file'};
			my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_suppl_out_script'};
			
			my $rprep_file = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/".$cfg_hash->{'prepare_rdata'};
			
			my $command = $cfg_hash->{'R2_path'}." CMD BATCH --no-save --no-restore '--args $analysis_name ".
					$cfg_hash->{'db_host'}." ".$cfg_hash->{'db_name'}." ".$cfg_hash->{'db_user'}." ".$cfg_hash->{'db_pass'}.
					" $outFolder ".extract_name($rearr_ann_out,0)." ".$outFolder."/".$out2nd_lev." $rprep_file'  $R_script $RLogPath";
		 
			print_and_log("$analysis_indication The command is: $command\n",$log_file);
			try_exec_command($command) or R_die($command,$R_script);		

			my $out2nd_lev_fin = $outFolder."/".$out2nd_lev.".".$cfg_hash->{'txt_ext'}; 
			my $out2nd_lev_fin_xls = $outFolder."/".$out2nd_lev.".".$cfg_hash->{'xlsx_ext'}; 
			
			#Convert the tsv file to xls in the same location
			print_and_log( "$analysis_indication Converting the TSV file ($out2nd_lev_fin) to XLS format ($out2nd_lev_fin_xls)...\n",$log_file);
			tsv_2_xls($out2nd_lev_fin,$out2nd_lev_fin_xls);	
			
						
			#Without pedfile the prioritization will not be run!
			if ($ped_file ne 'none'){				
				#################
				##Prioritization (Michele Pinelli)
				#################
				my $prog_used = "VarGenius-prioritization";	

				print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);
				#Parameters for script
				$RLogPath = $outFolder."/log/$analysis_name\_".$cfg_hash->{'prioritization_fold'}.$cfg_hash->{'R_log_file'};
				#Here we use a script that starts another one. The first one contains a configuration
				$R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_start_prioritization_script'};
				my $gene_table = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/".$cfg_hash->{'gene_full_table_rdata'};
				
				my $dataFolder = $cfg_hash->{$analysis_id."_f"}."/".$cfg_hash->{'data_fold'};
				my $candidate_panels = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/gene_annotations";
				my $annotation_panels = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/gene_annotations";
				my $candidate_variants = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/variant_annotations";
				my $annotation_variants = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'prioritization_fold'}."/variant_annotations";

				
				my $out2nd_lev_rdata = $outFolder."/".$out2nd_lev.".".$cfg_hash->{'rdata_ext'};
				
				#Check if directory exists, otherwise it is created
				unless(-d $priorwdir){
					print_and_log("$analysis_indication Creating folder $priorwdir...\n",$log_file);
					mkdir $priorwdir or die "$analysis_indication ERROR: can't create folder $priorwdir. Check permissions. \n";
				}
				`chmod 777 $priorwdir`;
				copy($out2nd_lev_rdata, $priorwdir) or die "Cannot copy $out2nd_lev_rdata in $priorwdir\n";
				
				my $date = `date +"%d-%m-%Y_%H-%M-%S"`;
				copy( $R_script, $priorwdir."/script_$date") or die "Cannot copy $R_script in $priorwdir/script_$date";
				
				my $args =  $analysis_name." ". #Analysis name.
									$priorwdir." ".# working dir, where input and output files are placed
								 $out2nd_lev.".".$cfg_hash->{'rdata_ext'}." ".#Output file in Rdata format
								 $dataFolder."/".$ped_file." ". #Pedfile full path 
								 $gene_table." ". #Gene annotation Rdata file
								 $candidate_panels." ".$annotation_panels." ".$candidate_variants." ".$annotation_variants." ".#annotation and panels  absolute paths
								 $cfg_hash->{'db_vcf_chr_str'}." ".$cfg_hash->{'db_vcf_pos_str'}." ". #Variants VCF columns
								 $cfg_hash->{'db_vcf_ref_str'}." ".$cfg_hash->{'db_vcf_alt_str'}." ".
								 $cfg_hash->{'dbSNPversion'}." ". #Dbsnp field
								 $cfg_hash->{'gene_field'}." ".$cfg_hash->{'nucl_change'}." ".$cfg_hash->{'aa_change'}." ". #gene annotation
								$cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_prioritization_script'}." "; #script path
				
				$command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $args'"." $R_script $RLogPath";
				print_and_log("$analysis_indication The command is: $command\n",$log_file);
				try_exec_command($command) or R_die($command,$R_script);					
			
				#Update the candidate variants table
				#my $cand_vars_f = $priorwdir."/".$group_name."_".$cfg_hash->{'prior_cand_var_f'};
				#parse_cand_genes_2_db($cfg_hash,$cand_vars_f);
			}		
		}else{
				print_and_log( "$analysis_indication ERROR: cannot start $prog_used. File $rearr_ann_out does not exists... \n",$log_file);
		}		
	}
	$out2nd_lev = "$outFolder/$out2nd_lev.".$cfg_hash->{'xlsx_ext'};	
	
	#######################################
	#
	#Print a set of statistics and plots
	if ($cfg_hash->{'get_coverage_plots'} eq 'YES'){	
		my $prog_used = "VarGenius-get_coverage_plots";
		print_and_log("$analysis_indication RUNNING  $prog_used..\n",$log_file);	
		#########################
		#Set and get needed data
		#########################
		#Some statistics will be take from the stats folder
		my $statsFolder = $cfg_hash->{$analysis_id.'_'.$cfg_hash->{'stats_step'}.'_out_f'};	
		#All the output plots and other data generated here will go in the finalout/img folder
		#print_and_log( "Setting the img folder using $analysis_id, $task, $analysis_id\_$task\_img_fold... \n",$log_file);
		my $imgFolder = $cfg_hash->{$analysis_id.'_'.$task.'_img_fold'};
		#print_and_log( "img folder :$imgFolder... \n",$log_file);
		
		
		#SAMPLES LIST: Getting the samples list associated with the analysisid
		print_and_log("$analysis_indication Getting the samples list associated with the analysisid... \n",$log_file);
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
		
		#DISEASE LIST FILES: Get the list of diseases files
		my @sel_genes_f_l = split(",",$cfg_hash->{'disease_genes_lists_f'});
		my $disease_l = join(",",@sel_genes_f_l);
		
		#Coverage-above file (All genes, all samples
		#Pick the extensions of the output from DepthOfCoverage
		my $gene_summary_ext = (split(",",$cfg_hash->{'DOC_outputs'}))[0];	
		#Get a table with gene coverage above $percentage_above for each sample
		my $sel_genes_f = "ALL_GENES";
		my $percentage_above = $cfg_hash->{'min_gene_coverage'};
		
		#Obtain a file containing for all the samples of the analysis and for all the genes the gene coverage above $percentage_above
		print_and_log("$analysis_indication Getting the gene coverage above $percentage_above... \n",$log_file);
		my $out_genecov_above = $statsFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'}."_".$sel_genes_f."_".$gene_summary_ext."_above$percentage_above.".$cfg_hash->{'txt_ext'};
		get_gene_coverage_above($cfg_hash,$analysis_id,$sel_genes_f,$percentage_above,$out_genecov_above,$statsFolder,$log_file);
		
		#The script which contains all the needed functionalities for the coverage plots
		my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_coverage_stats_script'};
		
		#The folder where are located all the gene panels
		my $gene_panels_fold = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'gene_panels_fold'};
		
		
		####################################
		#Get plots of coverage per-disease
		####################################
		my $plot_name = $cfg_hash->{'DISEASECOVBOXPLOT_name'};
		if ( -e $out_genecov_above){
			print_and_log("$analysis_indication Getting plots of coverage per-disease... \n",$log_file);
			my $RLogPath = $statsFolder."/log/$analysis_name\_$plot_name$percentage_above".$cfg_hash->{'R_log_file'};
			
			my $outsuffix = $imgFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'}."_".$sel_genes_f."_above$percentage_above";
			
			my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $out_genecov_above $gene_panels_fold $disease_l $outsuffix' $R_script $RLogPath";

			print "The command is:\n $command\n";
			try_exec_command($command) or R_die($command,$R_script);
		}else{
			print_and_log( "$analysis_indication WARNING: cannot start $plot_name. File $out_genecov_above does not exists... \n",$log_file);
		}	
				
			
			
		##############################		
		#Get tables for all genes and all disease with coverage less than 99%
		##############################
		$plot_name = $cfg_hash->{'COVLT99_name'};
		if ( -e $out_genecov_above){
			
			print_and_log("$analysis_indication Getting tables for all genes and all disease with coverage less than 99%... \n",$log_file);
							
			my $RLogPath = $statsFolder."/log/$analysis_name\_$plot_name".$cfg_hash->{'R_log_file'};
			my $outsuffix = $imgFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_tab_suff'}.$sel_genes_f."_$plot_name";
			
			my $command = $cfg_hash->{'R_path'}."  CMD BATCH --no-save --no-restore '--args $plot_name $out_genecov_above $samples_l $gene_panels_fold $disease_l $outsuffix' $R_script $RLogPath";
		 
			print_and_log("$analysis_indication The command is:\n $command\n",$log_file);
			try_exec_command($command) or R_die($command,$R_script);
		
			#These tables should be printed as links in the HTML page. Hence I add here this links to the outlist file
			#There is a table for each sample containing the data of all the diseases
			foreach my $sample_name (keys %{$group_sam}){  
				#Put the path to the VCF file in the outlist file
				if ( check_presence($outsuffix."_".$sample_name) ){
					print_and_log( "$analysis_indication Appending $outsuffix\_$sample_name path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
					my $line = "$sample_name\t-\t".$cfg_hash->{'outlist_cov_page'}."\t".$cfg_hash->{'outlist_cov_'.$plot_name}."\t".$outsuffix."_".$sample_name;
					append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
				}else{
					print_and_log( "$analysis_indication WARNING: $outsuffix\_$sample_name is not present. Please check what is happen. \n",$log_file);	
				}	
			}
		}else{
				print_and_log( "$analysis_indication WARNING: cannot start $plot_name. File $out_genecov_above does not exists... \n",$log_file);
		}	
		
		
		
		################################
		# High Quality Variants Plots (Dispersion)
		################################
		#Input: uses the rearranged output from VarGenius
		$plot_name = $cfg_hash->{'HQVARPLOT_name'};
		if ( -e $rearr_ann_out and !(-z $rearr_ann_out)){
			
			print_and_log("$analysis_indication Getting High Quality Variants Plots (Dispersion) : ",$log_file);	
			my $RLogPath = $statsFolder."/log/$analysis_name\_$plot_name".$cfg_hash->{'R_log_file'};
			my $outsuffix = $imgFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'};
			
			my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $rearr_ann_out $samples_l $outsuffix' $R_script $RLogPath";
		 
			print_and_log("$analysis_indication The command is: $command\n",$log_file);
			try_exec_command($command) or R_die($command,$R_script);
		}else{
				print_and_log( "$analysis_indication ERROR: cannot start $plot_name. File $rearr_ann_out does not exists... \n",$log_file);
		}	
	}
		
	#Gene annotation. A new file with .ann extension that for each gene into the output table 
	#has the gene annotation associated
	my $gene_annotation_out = $rearr_ann_out.".".$cfg_hash->{"gene_annotation_ext"};
	if ($cfg_hash->{'gene_annotation'} eq 'YES'){			
		if ( -e $rearr_ann_out and !(-z $rearr_ann_out)){
			my $prog_used = "VarGenius";	
			print_and_log("$analysis_indication RUNNING $prog_used gene_annotation..\n",$log_file);	

			print_and_log( "$analysis_indication Executing gene annotation of $rearr_ann_out...\n",$log_file);
			get_gene_annotation($cfg_hash,$rearr_ann_out,$gene_annotation_out,$log_file);
			
			#get_all_genes_annotation($cfg_hash,$gene_annotation_out.".ALLGENES",$log_file);#CAN BE REMOVED
			$partDuration = time - $partTime;
			print_and_log( "$analysis_indication Finished $prog_used-gene_annotation $partDuration seconds\n",$log_file);
			$partTime = time;					
		}else{
				print_and_log( "$analysis_indication WARNING: cannot start gene_annotation. File $rearr_ann_out does not exists... \n",$log_file);
		}	

	}		
		
	################################		
	#Ad HOC modifies for given users
	if ( $cfg_hash->{'do_users_modifies'} eq 'YES' ){
			my $user_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_userid'},
															$cfg_hash->{'db_analysis_id'},$analysis_id);
			
			#TAG Panel Genes: asked by Sandro Banfi, this snip of code uses data from a main data folder ad hoc for the user
			#This folder is DATA/userid_data and must contain a text file with the name of the panel (e.g. retinopathy) and 
			#a text file with the name of the candidates panel (retinopathy_cand)
			#The variable panel_to_tag is only declared in the user_config but defined into the userid_config.txt.
			#Any user can define a different panel with the same rules as just declared.
			if ( (defined $cfg_hash->{'panel_to_tag'}) and length($cfg_hash->{'panel_to_tag'}) > 2 ){
					my $panel = $cfg_hash->{'panel_to_tag'};
					
					print_and_log( "$analysis_indication Generating an additional output (XLSX) adding a flag for each $panel gene...\n",$log_file);
					my $out_path = tag_panel_genes($cfg_hash,$analysis_id,$rearr_ann_out,$panel,$user_id,$outFolder,$log_file);
					my $out_path_xls = extract_name($out_path,"noext").".".$cfg_hash->{'xlsx_ext'};
					
					#Shifting the column after the gene name. The column has the same name of the panel
					shift_column($out_path,$panel,$cfg_hash->{'gene_field'},$log_file);
					#Convert to xlsx
					print_and_log( "$analysis_indication Converting the TSV file ($out_path) to XLS format ($out_path_xls)...\n",$log_file);
					tsv_2_xls($out_path,$out_path_xls);								
					#Remove the txt
					if ( $cfg_hash->{'remove_temp'} eq 'YES'){
						delete_file($out_path);
					}
					#Add to outlist file
					if ( check_presence($out_path_xls) ){
						print_and_log( "$analysis_indication Appending $out_path_xls path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
						my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t$panel filter: ".extract_name($out_path_xls,1)."\t".$out_path_xls;
						append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
					}else{
						print_and_log( "$analysis_indication WARNING: $out_path_xls is not present. Please check what is happen. \n",$log_file);	
					}			
			}			
	}
	
	########################################
	##########FILLING OUTLIST FILE##########
	########################################
	
	#Put the path to the VCF file in the outlist file
	if ( check_presence($vcf_file) ){
		print_and_log( "$analysis_indication Appending $vcf_file path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_vcf'}."\t".$vcf_file;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	}else{
		print_and_log( "$analysis_indication WARNING: $vcf_file is not present. Please check what is happen. \n",$log_file);	
	}		

	#Put the path to the Final tabular output in the outlist file
	if ( check_presence($rearr_ann_out) ){
		print_and_log( "$analysis_indication Appending $rearr_ann_out path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_out'}."\t".$rearr_ann_out;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	}else{
		print_and_log( "$analysis_indication WARNING: $rearr_ann_out is not present. Please check what is happen. \n",$log_file);	
	}		

	#Put the path to the XLS version of the Final tabular output in the outlist file
	if ( check_presence($rearr_ann_out_xls) ){
		print_and_log( "$analysis_indication Appending $rearr_ann_out_xls path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_xls'}."\t".$rearr_ann_out_xls;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );
	}else{
		print_and_log( "$analysis_indication WARNING: $rearr_ann_out_xls is not present. Please check what is happen. \n",$log_file);	
	}		

	#Put the supplementary  output into the outlist file if present
	if ( check_presence($out2nd_lev) ){
		print_and_log( "$analysis_indication Appending $out2nd_lev path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_2ndlev'}."\t".$out2nd_lev;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );		
	}else{
		print_and_log( "$analysis_indication WARNING: $out2nd_lev is not present. Please check what is happen. \n",$log_file);	
	}

	#Put the supplementary  output into the outlist file if present
	if ( check_presence($gene_annotation_out) ){
		print_and_log( "$analysis_indication Appending $gene_annotation_out path to ".$cfg_hash->{$analysis_id."_outlist_file"}." file... \n",$log_file);
		my $line = "-\t-\t".$cfg_hash->{'outlist_output_page'}."\t".$cfg_hash->{'outlist_final_geneann'}."\t".$gene_annotation_out;
		append_str_2_file_if_path_notexist($cfg_hash->{$analysis_id."_outlist_file"},$line );		
	}else{
		print_and_log( "$analysis_indication WARNING: $gene_annotation_out is not present. Please check what is happen. \n",$log_file);	
	}

	#Plots describing statistics about the overall samples into a given project or set of analyses
	if ($cfg_hash->{'update_overall_samples_plots'} eq 'YES'){
		print_and_log("$analysis_indication RUNNING  update_overall_samples_plots..\n",$log_file);	
		#########################
		#Set and get needed data
		#########################
		#Some statistics will be take from the stats folder
		my $statsFolder = $cfg_hash->{$analysis_id.'_'.$cfg_hash->{'stats_step'}.'_out_f'};	
		#All the output plots and other data generated here will go in the finalout/img folder
		#print_and_log( "Setting the img folder using $analysis_id, $task, $analysis_id\_$task\_img_fold... \n",$log_file);
		my $imgFolder = $cfg_hash->{$analysis_id.'_'.$task.'_img_fold'};
			
		#The script which contains all the needed functionalities for the coverage plots
		my $R_script = $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_coverage_stats_script'};
			
		##################################
		#All samples coverage plot proportion
		#################################
		#Uses the sample cumulative coverage proportion output form DepthOfCoverage from any file belonging to the 
		#same user id and same target file and makes a plot describing what is the sample coverage at different levels for each sample
		#.sample_cumulative_coverage_proportions header is the following:
		#			SampleName gte_0 gte_1 .. gte_500
		#This plot contains information for all the samples
		print_and_log("$analysis_indication Getting all samples coverage plot (percentage of locus)... \n",$log_file);	
		my $plot_name = $cfg_hash->{'ALLSAMPLESCOVPLOT_name'};
		my $all_sam_cov_tab = $statsFolder."/$analysis_name\_$plot_name\_PROPS.txt";
		my $DOC_sam_cov_type = $cfg_hash->{'DOC_sam_cum_prop'};
		
		my $success1 = get_all_samples_coverage_table($cfg_hash,$analysis_id,$all_sam_cov_tab,$DOC_sam_cov_type,$log_file);
		#my $success = 2;#DEBUGCODE
		if ( -e $all_sam_cov_tab and !(-z $all_sam_cov_tab)){
			my $RLogPath = $statsFolder."/log/$analysis_name\_$plot_name".$cfg_hash->{'R_log_file'};

			my $outsuffix = $imgFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'};
			my $command = "";
			my $plot_title = $cfg_hash->{'DOC_sam_cum_prop_plot_title'};
			#Success = 2 means the kinship can be used
			if ( $success1 > 0){
				if ($success1 == 1){
					$command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $all_sam_cov_tab 0 ".
											$cfg_hash->{'DOC_sam_cum_cov_prop_fields'}." $plot_title $outsuffix' $R_script $RLogPath";				
				}
				if ($success1 == 2){
					$command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $all_sam_cov_tab 1 ".
											$cfg_hash->{'DOC_sam_cum_cov_prop_fields'}." $plot_title $outsuffix' $R_script $RLogPath";								
				}
				print_and_log( "$analysis_indication The command is: $command\n",$log_file);
				try_exec_command($command) or R_die($command,$R_script);						
			}else {
				print_and_log( "$analysis_indication WARNING: cannot start $plot_name. Bad result from get_all_samples_coverage_table... \n",$log_file);
			}

		}else {
				print_and_log( "$analysis_indication WARNING: cannot start $plot_name. Empty result from get_all_samples_coverage_table... \n",$log_file);
		}		

		##################################
		#All samples coverage plot count of locus
		#################################
		#Uses the sample cumulative coverage count output form DepthOfCoverage from any file belonging to the 
		#same user id and makes a  plot describing what is the number of locus of the sample coverad at different levels for each sample
		print_and_log("$analysis_indication Getting all samples coverage plot (number of locus)... \n",$log_file);	
		$plot_name = $cfg_hash->{'ALLSAMPLESCOVPLOT_name'};
		my $all_sam_cov_counts_tab = $statsFolder."/$analysis_name\_$plot_name\_COUNTS.txt";
		$DOC_sam_cov_type = $cfg_hash->{'DOC_sam_cum_counts'};
		my $success3 = get_all_samples_coverage_table($cfg_hash,$analysis_id,$all_sam_cov_counts_tab,$DOC_sam_cov_type,$log_file);

		
		##################################
		#All samples reads plot
		#################################
		#Uses the output from FLAGSTAT of all the samples to print a plot with points for both
		#the total number of reads and the properly paired		
		
		$plot_name = $cfg_hash->{'ALLSAMPLESREADSPLOT_name'};
		print_and_log("$analysis_indication Getting all samples reads number plot ($plot_name) ... \n",$log_file);
		my $all_sam_reads_tab = $statsFolder."/$analysis_name\_$plot_name.txt";
		#Gets a table for all samples from same user id and target giving: SampleName TotalReads TotalRemoved ProperlyPaired
		#In addition we get the values for the current analysis
		my ($success2,$tot_reads_subject,$reads_removed_subject,$prop_paired_subject) = get_all_samples_reads_number($cfg_hash,$analysis_id,$all_sam_reads_tab,\@steps_array,$log_file);
		
		#If the table was succesfully created
		if ( -e $all_sam_reads_tab and !(-z $all_sam_reads_tab)){
			#From the table remove lines containing absent values
			my $all_sam_reads_tab_temp = $all_sam_reads_tab.".temp";
			delete_rows_containing($all_sam_reads_tab,"-",$all_sam_reads_tab_temp);
		
			my $RLogPath = $statsFolder."/log/$analysis_name\_$plot_name".$cfg_hash->{'R_log_file'};

			my $outsuffix = $imgFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'};
			my $command = "";
		
			#Success = 2 means the kinship can be used
			if ( $success2 > 0){
				if ($success2 == 1){
					$command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $all_sam_reads_tab_temp 0 ".
											$cfg_hash->{'flagstat_totreads_field'}.",".$cfg_hash->{'flagstat_propaired_field'}." $outsuffix' $R_script $RLogPath";				
				}
				if ($success2 == 2){
					$command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $plot_name $all_sam_reads_tab_temp 1 ".
											$cfg_hash->{'flagstat_totreads_field'}.",".$cfg_hash->{'flagstat_propaired_field'}." $outsuffix' $R_script $RLogPath";								
				}
				print_and_log( "$analysis_indication The command is: $command\n",$log_file);
				try_exec_command($command) or R_die($command,$R_script);						
			}else {
				print_and_log( "$analysis_indication WARNING: cannot start $plot_name. Bad result from get_all_samples_reads_number... \n",$log_file);
			}
			#delete_file($all_sam_reads_tab_temp);
		}else {
			print_and_log( "$analysis_indication WARNING: cannot start $plot_name. Empty result from get_all_samples_reads_number... \n",$log_file);
		}				
	}		
		
	
	#Print an HTML informative website about the results obtained
	if ($cfg_hash->{'generate_html'} eq 'YES'){			
		my $imgFolder = $cfg_hash->{$analysis_id.'_'.$task.'_img_fold'};
		
		my $prog_used = "VarGenius-generate_html";	
		print_and_log("$analysis_indication RUNNING $prog_used..\n",$log_file);	

		print_and_log( "$analysis_indication Executing HTML productions of $rearr_ann_out...\n",$log_file);
		
		##Write the informations from the file into an hash so that we can extract
		#the info in the correct order: QC, COVERAGE, OUTPUT and write them in HTML
		my $paths_hash;
		my $outlist_f = $cfg_hash->{$analysis_id."_outlist_file"};
		
		if ( !-e $outlist_f ) {log_and_exit("$analysis_indication ERROR: Cannot proceed with generate_html. $outlist_f missing.\n",$log_file);}
		if ( !-e $rearr_ann_out) {log_and_exit("$analysis_indication ERROR: Cannot proceed with generate_html. $rearr_ann_out missing.\n",$log_file);}

		print_and_log("$analysis_indication Reading $outlist_f to generate an hash...\n",$log_file);
		outlist_to_hash($outlist_f,\$paths_hash);
		
		#################################
		#Print R Statistics
		#################################
		#Print a text file with variants stats
		my $stats_out = $rearr_ann_out.".".$cfg_hash->{'Rstats_ext'};
		print_and_log( "$analysis_indication Printing statistics about the results obtained in $stats_out...\n",$log_file);
		print_variants_stats($cfg_hash,$analysis_id,$rearr_ann_out,$task,$log_file,$stats_out);
		#Update table with statistics from this file 
		print_and_log( "$analysis_indication Loading statistics into the database...\n",$log_file);
		fields_2_db($cfg_hash,$cfg_hash->{'db_analysis_id'},$analysis_id,$stats_out,$log_file);
		
	
		#################################
		###IGV SESSION
		#################################
		#Print the IGV session file using the paths_hash
		my $xml_igv_session_out = $outFolder."/".$analysis_name."_".$cfg_hash->{'igv_sess_ext'};
		print_and_log( "$analysis_indication Generating the IGV session file in $xml_igv_session_out...\n",$log_file);
		generate_igvsession_xml($paths_hash,$cfg_hash,$xml_igv_session_out,$log_file);
		#Add the IGV session file to the paths hash
		$paths_hash->{ $cfg_hash->{"outlist_cov_page"}}->{'-'}->{'-'}->{$cfg_hash->{'outlist_igvsession_xml'}}->{'path'} = $xml_igv_session_out;

		
		#################################
		##HTML WEB PAGES
		#################################
		print_and_log( "$analysis_indication Generating the HTML website in $outFolder...\n",$log_file);

		################################
		#DATA COMPRESSION
		#Compress data and add path to the paths hash (first page)
		print_and_log( "$analysis_indication Compressing data to store into $outFolder...\n",$log_file);
		my $compressed_data_path = compress_data_to_store($paths_hash,$outFolder,$log_file);
		$paths_hash->{ $cfg_hash->{"outlist_output_page"}}->{'-'}->{'-'}->{$cfg_hash->{'outlist_storage_compressed'}}->{'path'} = $compressed_data_path;
		
		#Arrays for info of pages
		my @descText = ($cfg_hash->{'html_index_desc'},$cfg_hash->{'html_qc_desc'},$cfg_hash->{'html_cov_desc'});
		my @links = ($cfg_hash->{'html_index_page'},$cfg_hash->{'html_qc_page'},$cfg_hash->{'html_cov_page'});				
		
		#INDEX PAGE: REPORT AND DATA
		#Reorder the hash and convert in a HTML table OUT
		my $html_string = paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_output_page'},$paths_hash,$cfg_hash->{"html_host"},$cfg_hash->{"work_fold"}); 
		
		#########STATISTICS THROUGH THE ANALYSIS
		#Add statistics obtained during the analysis using the outstats file
		print_and_log("$analysis_indication STATISTICS DURING THE ANALYSIS:\n",$log_file);
		if ( -e $cfg_hash->{$analysis_id."_outstats_file"}){
			print_and_log("$analysis_indication Printing ".$cfg_hash->{$analysis_id."_outstats_file"}." to HTML...\n",$log_file);
			$html_string .= file_2_html($cfg_hash->{$analysis_id."_outstats_file"});
		}else{print_and_log("$analysis_indication ".$cfg_hash->{$analysis_id."_outstats_file"}." file does not exist...\n",$log_file);}


		my $tables_list_file = $outFolder."/".$cfg_hash->{'statstables_file'};
		################################
		#QUALITY STATISTICS COMPARISON
		#Get the quality statistics compared with the database average to print in the web page		
		#Prints two tables (statistics for variants and for the alignment) and adds the paths to the statstables_file
		evaluate_quality_of_sequencing($cfg_hash,$analysis_id,$tables_list_file,$log_file);
		
		########STATISTICS FROM TABLES
		#Currently I'm extracting the following two tables
		# 1. sample_summary table (Statistics of coverage for the samples (into get_statistics_tables))
		# 2. alignment statistics using flagstat output (into the function run_coverage_steps)
		print_and_log("$analysis_indication STATISTICS FROM TABLES:\n",$log_file);

		print_and_log( "$analysis_indication Extract and print statistics from tables (from $tables_list_file)...\n",$log_file);
		#print_and_log( "$analysis_indication The files will be printed into $outFolder...\n",$log_file);#DEBUGCODE
		#Extract additional statistics from data obtained during the analysis and list them into $tables_list_file file	
		get_statistics_tables($cfg_hash,$analysis_id,$outFolder,$log_file);
		#Add to the HTML page the tables obtained from stats_table
		if (-e $tables_list_file ) {
			print_and_log( "$analysis_indication $tables_list_file exists. Using table2html...\n",$log_file);#DEBUGCODE
			$html_string .= table_list_2_html($tables_list_file,$cfg_hash->{'remove_temp'});
		}else{print_and_log("$analysis_indication WARNING: $tables_list_file has not been found..\n",$log_file);}
		
				
		#########VARIANTS STATISTICS		
		#Add to the HTML page statistics obtained with a R script to the html string
		$html_string .= var_stats_2_html($stats_out,$cfg_hash->{'mult_ann_sep'});


		##########MAKE WEBSITE
		#Print the HTML page using the string created
		html_page($cfg_hash,'none',\@descText,0,\@links,undef,$html_string,$outFolder);

		#Reorder the hash and convert in a HTML table QC
		$html_string = paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_qc_page'},$paths_hash,$cfg_hash->{"html_host"},$cfg_hash->{"work_fold"}); 
		#Print the HTML page using the string created
		html_page($cfg_hash,'none',\@descText,1,\@links,undef,$html_string,$outFolder);
		
		#Reorder the hash and convert in a HTML table COVERAGE
		$html_string = paths_hash_to_html($cfg_hash,$cfg_hash->{'outlist_cov_page'},$paths_hash,$cfg_hash->{"html_host"},$cfg_hash->{"work_fold"}); 
		#The Coverage page will contain also some images which are plots
		my $keyword = $cfg_hash->{'html_cov_img_suff'};
		#Print the HTML page using the string created
		html_page($cfg_hash,$imgFolder,\@descText,2,\@links,$keyword,$html_string,$outFolder); 								 
		
		#print_and_log( "The string is: $string\n",$log_file);
		
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
		$partTime = time;		
	}
	

	
	
	#Send email when the analysis is completed
	if ($cfg_hash->{'send_email'} eq 'YES'){			
		my $prog_used = "VarGenius";	
		print_and_log("$analysis_indication   RUNNING $prog_used send_email..\n",$log_file);	
		
		
		print_and_log( "$analysis_indication Sending an email to ".$cfg_hash->{'email_recipients'}."...\n",$log_file);
		my $subject = $cfg_hash->{'email_subject'}.": $analysis_name";
		my $message = $analysis_name." ".$cfg_hash->{'email_message'};
		send_email($cfg_hash->{'email_author'},$cfg_hash->{'email_recipients'},$subject,$message,$log_file);
				
		$partDuration = time - $partTime;
		print_and_log( "$analysis_indication Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
		$partTime = time;		
	}

	##################################################
	##INSERT INFORMATIONS INTO THE DB FOR THE WEBSITE
	##################################################
	
	#Update the analyses table with the date and VarGenius version
	my $fields = $cfg_hash->{'db_analysis_date'}.";".$cfg_hash->{'db_progver'};
	my ($d,$m,$y) = (localtime)[3,4,5];
	my $year = $y+1900; # In PERL year
	$m = $m + 1; #Month in perl starts from 0 (0:january,11:december)
	my $values = "'$year\_$m\_$d'".";".$cfg_hash->{'VarGeniusVer'};
	my $fields2 = $cfg_hash->{'db_analysis_id'};
	my $values2 = $analysis_id;
	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);
					
	#Insert some information related with the analysis into the INFO field 	for the group
	$fields = $cfg_hash->{'db_analysis_info'};
	my $targetbed = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_targetbed'}
						,$cfg_hash->{'db_analysis_id'},$analysis_id);
	$values = "'".$targetbed."'";
	$fields2 = $cfg_hash->{'db_analysis_id'};
	$values2 = $analysis_id;
	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);	
	
	##################################################

	##################################################
	## UPDATE STATUS
	##################################################
	
	#Update the group table with the final status.
	#At this point all the analyses per sample are finished and
	#all the analysis is finished
	$fields = $cfg_hash->{'db_analysis_status'};
	$values = $cfg_hash->{'analysis_fout'};
	$fields2 = $cfg_hash->{'db_analysis_id'};
	$values2 = $analysis_id;
	
	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$fields,$values,$fields2,$values2);
	##################################################
	
	
	##################################################
	## FREE THE DATABASE (IF USED THE EXCLUSIVE USAGE)
	##################################################	
	
	#Before	to leave, check if this is the last analysis using the db
	#and if yes leave the database: db_busy=0
	#Get the id of the job that is using the db
	my $analysisin = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
		$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_analysis_id'},$cfg_hash->{'db_info_row'},"1");
	if ( $analysisin eq $analysis_id){
		print_and_log("$analysis_indication Freeing the db (db_busy=0) since $analysis_id is the last using it ...\n",$log_file);
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_busy'},"0",$cfg_hash->{'db_info_row'},"1");	
	}else{
		print_and_log("$analysis_indication Analysis $analysisin will use the database. This analysis ($analysis_id) does not have to free it up with db_busy=0...\n",$log_file);
	}	


	##################################################
	## 	Save the config user into the log file
	##################################################
	#print_and_log("\n\nHere it is the user configuration file you used: $config_file...\n",$log_file);
	#append_hash_to_file($log_file,$cfg_hash);
	
}

=head2 paths_hash_to_text

 Title   : paths_hash_to_text
 Usage   : paths_hash_to_text( config_file => the config hash );

 Function: Reads the paths hash content  and writes a file ready for the
					HTML conversion
  
 Returns : the hash given in input is filled..
 
=cut
sub paths_hash_to_text{
	my $cfg_hash = shift;
	my $paths_hash = shift;
	#my $pages = shift;
	my $out  = shift;
	
	open (FILE,">$out") or die "ERROR [$!]: can't open file: $out. Check permissions.\n";
	print FILE "Quality Check and Trimming\n";
	
	my @pages = ($cfg_hash->{'outlist_qc_page'},$cfg_hash->{'outlist_cov_page'},$cfg_hash->{'outlist_output_page'});
	
	
	foreach my $page (@pages){
		print FILE $page."\n";
	#Quality Check Page
		while (my ($sample_name, $value) = each %{ $paths_hash->{$page} } ) {
			if ($sample_name ne '-'){print FILE $sample_name."\n"}
			while (my ($readf_name, $value2) = each %{ $paths_hash->{$page}->{$sample_name}} ) {
				if ($readf_name ne '-'){print FILE $readf_name."\n"}
				while (my ($desc, $value3) = each %{ $paths_hash->{$page}->{$sample_name}->{$readf_name}} ) {
					#if ($desc ne '-'){print FILE $desc."\n"}
					print FILE $desc."\t".$paths_hash->{$page}->{$sample_name}->{$readf_name}->{$desc}->{'path'}."\n";
				}
			}
		}
	}
	
	close(FILE);

}
	
=head2 outlist_to_hash

 Title   : outlist_to_hash
 Usage   : outlist_to_hash( config_file => the config hash );

 Function: Reads the outlist file and builds an hash for straightforward access
						The hash structure is: PAGE->SAMPLE->READFILE->DESCRIPTION->PATH
					
					So that, later, for each page, for each sample and for each readfile can be print a link with a description	
					
  
 Returns : the hash given in input is filled..
 
=cut
sub outlist_to_hash{
	my $outlist_f = shift;
	my ($paths_hash) = shift;

	my $samid_fld_ind = 0;
	my $rdf_fld_ind = 1;
	my $page_fld_ind = 2;
	my $desc_fld_ind = 3;
	my $path_fld_ind = 4;
			
 #Create
 open (FILE,"<$outlist_f") or die "ERROR [$!]: can't open file: $outlist_f. Check permissions.\n"; 
 while (my $line = <FILE>){
	chop($line);
	my @fields = split("\t",$line);
	
	#The hash contains data separated by page
	$$paths_hash->{$fields[$page_fld_ind]}->{$fields[$samid_fld_ind]}->{$fields[$rdf_fld_ind]}->{$fields[$desc_fld_ind]}->{'path'} = $fields[$path_fld_ind];
 }
 close FILE;
}





	
=head2 run_import_2_db_steps

 Title   : run_import_2_db_steps
 Usage   : run_import_2_db_steps(  config_file => the config hash
								);

 Function: Executes the importation of data  into the database 
						- Variants from the VCF 
						- annotations from the Annovar output
 
 Returns :
 
=cut
sub run_import_2_db_steps{
	my $config_file = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

	print_and_log("IMPORT VCF function\n",$log_file);#DEBUGCODE

	#Get the parameters from input
	my $cfg_hash = load_hash($config_file);
	

	##################################################
	######Before	to start take posses of the database db_busy=1
	print_and_log("Setting analysis $analysis_id to be the last waiter of the db for the final output ...\n",$log_file);

	update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
			$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_busy'},"1",$cfg_hash->{'db_info_row'},"1");	
	#################################################
	
	
	#Used to take the last input
	my $previous_task = $cfg_hash->{'varfilt_step'};	
	my $step = $cfg_hash->{'vcfimport_step'};
	
	#Input and output folders
	my $inFolder = $cfg_hash->{$analysis_id.'_'.$previous_task.'_out_f'};
	my $outFolder = $cfg_hash->{$analysis_id.'_'.$task.'_out_f'};
	
	#Getting the sample name
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $params;
	print_and_log("Getting sample configuration for analysisid $analysis_id .\n",$log_file);#DEBUGCODE
	getSampleConfiguration_locked(\$params,$cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $vcf_file = $inFolder."/".build_input_name_from_executed($params,$step,
								$params->{$cfg_hash->{'db_analysis_name'}},\@steps_array3).".".$cfg_hash->{'vcf_ext'};
								
	#File name for the VCF with variants that must be annotated (missing)
	print_and_log("Getting filenames and flags\n",$log_file);#DEBUGCODE
	my $vcf_to_ann = $inFolder."/".extract_name($vcf_file,1).".".$cfg_hash->{'miss_ann_ext'};
	#File VCF with all variants
	my $vcf_all_var = $inFolder."/".extract_name($vcf_file,1).".".$cfg_hash->{'allvar_ext'};
		
	my $vcf_import = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'vcfimport_step'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);

	my $import_missing = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'import_missing_step'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);


		print_and_log("Checking if vcfimport is active ".$cfg_hash->{'vcfimport'}."\n",$log_file);#DEBUGCODE
																							
	#VCF import. Reads the  variants from the VCF file filling the statistics and genotype_sample table
	#and generates a VCF with all the missing	variants and another with all the variants
	if ($cfg_hash->{'vcfimport'} eq 'YES'){

		#Check if variants for this analysis were already imported
		my $query = "SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_genotype_sample_table'}." WHERE ".
														$cfg_hash->{'db_analysis_id'}." = $analysis_id;";	
		#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
		my $anal_ids_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_analysis_id'});
		#Now check if they are present
		if ( scalar(keys %{$anal_ids_fetch}) > 0 ){
			log_and_exit("ERROR: Cannot import variants from analysis  $analysis_id. I found ".scalar(keys %{$anal_ids_fetch}).
			" genotypes into the ".$cfg_hash->{'db_genotype_sample_table'}." table. Please use vcfimport=NO if you don't want to insert variants. Exiting...\n",$log_file);
		}				
		
		my $prog_used = "VarGenius";
		#Import if the same analysis has not been inserted yet
		#if ( $vcf_import == 0){
			my $step = $cfg_hash->{'vcfimport_step'};
			print_and_log( "Executing $step. The VCF file $vcf_file will be imported in the db".
			" and another VCF ($vcf_to_ann), including all the non-existent variants, will be created for the".
			" annotation\n",$log_file);
			#import_vcf_2_db($cfg_hash,$prog_used,$analysis_id,$log_file,$vcf_file,$vcf_to_ann,$vcf_all_var);
			import_vcf_2_db_fast($cfg_hash,$prog_used,$analysis_id,$log_file,$vcf_file,$vcf_to_ann,$vcf_all_var);
			#Save status in the db by searching for the specific sample name
			#print_and_log("Updating analysis status for group $sample_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);#DEBUGCODE
			update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
			$partDuration = time - $partTime;
			print_and_log( "Finished\t$prog_used\tanalysis\t$analysis_name\t$partDuration\tseconds\n",$log_file);
			$partTime = time;	
		#}else{
		#	print_and_log( "This analysis results ($analysis_id) have been already imported in the database... \n",$log_file);
		#}
	}


	#After importing variants, we recompute the frequencies
	if ($cfg_hash->{'update_frequencies'} eq 'YES'){

		my $step = "update_frequencies";
		print_and_log( "Executing $step. The variants frequencies for samples from all sequencing types will be update.\n",$log_file);
		run_recompute_variants_frequencies($config_file);
	}
	

	#Research group specific requests
	my $resgroupid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_resgroupid'},
									$cfg_hash->{'db_analysis_id'},$analysis_id);
	
	#Get the name of the re-arranged output
	my $rearr_ann_out =  $outFolder."/".$params->{$cfg_hash->{'db_analysis_name'}}.".".$cfg_hash->{'ann_build_ver'}.".final_out_resgr$resgroupid.txt";
	#QUESTO CODICE PROBABILMENTE NON E' PIU' UTILE SE USO SEMPRE IL REANNOTE GLOBALE
	#if ($cfg_hash->{'import_missing'} eq 'YES'){
		#if ( -e $rearr_ann_out) {
		#my $prog_used = "VarGenius-import-Annotation";	
		#my $step = $cfg_hash->{'import_missing_step'};
		#print_and_log("[an:$analysis_id] (".scalar(localtime)."): RUNNING $prog_used import missing..\n",$log_file);	


			#import_vgout_2_db($cfg_hash,$rearr_ann_out,$prog_used,$analysis_id,$log_file);
			
			##Save status in the db by searching for the specific sample name
			#print_and_log("Updating analysis status for analysis $analysis_id in ".$cfg_hash-> {'db_name'}." ($step=1) ...\n",$log_file);
			#update_analysis_status_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},
					#$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},
					#$cfg_hash->{'db_analysis_id'},$analysis_id,$step);
			#$partDuration = time - $partTime;
			#print_and_log( "Finished\t$prog_used for analysis $analysis_id in\t".$partDuration." seconds\n",$log_file);
			#$partTime = time;	

		#}else{
			#print_and_log("WARNING: Cannot find $rearr_ann_out. Annotation import will not be performed...\n",$log_file);
		#}
	#}
	
	##Import data about the coverage
	if ($cfg_hash->{'import_coverage_data'} eq 'YES'){
		my $prog_used = "VarGenius-import_coverage_data";
		print_and_log("[an:$analysis_id] (".scalar(localtime)."): RUNNING $prog_used ..\n",$log_file);	
	
		#Check the coverage of any of the samples in this analysis has already been inserted
		#Get samples ids
		my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
		print_and_log( "Executing: $query\n",$log_file);
		my $sample_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

		#For each sample check its presence
		foreach my $sample_id (keys %{$sample_ids}){
			#Check if coverage data for this analysis were already imported
			$query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_genecoverage_table'}." WHERE ".$cfg_hash->{'db_sample_id'}." = $analysis_id;";	
			#print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
			my $sample_ids_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});
			#Now check if they are present
			if ( scalar(keys %{$sample_ids_fetch}) > 0 ){
				log_and_exit("ERROR: Cannot import variants from analysis  $analysis_id. I found ".scalar(keys %{$sample_ids_fetch}).
				" genes coverage values for sampleid $sample_id into the ".$cfg_hash->{'db_genecoverage_table'}." table. Please use import_coverage_data=NO if you don't want to insert coverage data. Exiting...\n",$log_file);
			}		
			
		}
			
		my $statsFolder = $cfg_hash->{$analysis_id.'_stats_out_f'};
		my $gene_summary_ext = (split(",",$cfg_hash->{'DOC_outputs'}))[0];
		my $percentage_above = $cfg_hash->{'min_gene_coverage'};
		
		my $out_genecov_above = $statsFolder."/".$analysis_name."_".$cfg_hash->{'html_cov_img_suff'}."_ALL_GENES_".$gene_summary_ext."_above$percentage_above.".$cfg_hash->{'txt_ext'};
	
		if ( -e $out_genecov_above){
			load_table_into_db($cfg_hash ,$out_genecov_above,$analysis_id,"sample","genes","genecoverage",$log_file);
		}else{
			print_and_log("WARNING: Cannot find $out_genecov_above. $prog_used  will not be performed...\n",$log_file);
		}
	}
	
		
	#Before	to leave, check if this is the last group using the db
	#and if yes leave the database: db_busy=0
	#Get the id of the job that is using the db
	my $analysisin = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
		$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_analysis_id'},$cfg_hash->{'db_info_row'},"1");
	if ( $analysisin eq $analysis_id){
		print_and_log("Freeing the db (db_busy=0) since $analysis_id is the last using it ...\n",$log_file);
		update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
				$cfg_hash->{'db_pass'},$cfg_hash->{'db_info_table'},$cfg_hash->{'db_busy'},"0",$cfg_hash->{'db_info_row'},"1");	
	}else{
		print_and_log("Analysis $analysisin will use the database. This analysis ($analysis_id) do not have to free it up with db_busy=0...\n",$log_file);
	}
}

