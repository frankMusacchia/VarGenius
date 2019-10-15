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
    
package LIB::samtools_management;
## finalout_management.pm
#Author: Francesco Musacchia  (2019)
#Permits the management of the subroutines used to manage the output of VarGenius
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(  );
}

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
				compress_data_to_store);


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


								


=head2 run_fast_final_out_steps

 Title   : run_fast_final_out_steps
 Usage   : run_fast_final_out_steps(  config_file => the config hash
								);

 Function: Uses the final VCF:
						- extracts an input for Annovar from the VCF

 Returns :
 
=cut
sub fout_create_vcf_to_ann{

			
}
