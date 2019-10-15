#VarGenius - History

#Spring 2019

On spring 2019 we have completed a series of important updates:

1. Since the approval of the GDPR (Global Data Proctection Regulation) we decided that the results of the analysis must be stored into an area of the cluster
	protected from external access. Hence, we also reduced the output present into the HTML file.
	To achieve this objective we inserted an additional keyword to the outlist file containing the paths to the output files.
	The keyword will be useful to understand if the file must be:
	[STORE]: stored into a TAR file
	[WEB]: used into the HTML web page
	[IGV]: used as an output for IGV visualization
	The files with the [STORE] keyword can be downloaded as a TAR file from the HTML web page. The [WEB] files are uncompressed.
	The BAM files, with the [IGV] keyword cannot be downloaded but only viewed in IGV.
	
	1a. At the end of the fast_final_output_steps function now we have the call to save_analysis that copies all files indicated in the outlist file with the keywords [IGV] or [WEB].
	1b. the HTML file contains now only data that can be shown for privacy. While other supplementary information are stored into a tar file.

2. The consensus variant calling developed during the summer 2018 has been updated. Now different calling software can be used.
	If a job for a new calling software is used this strip of code must be added:
	
	#At the end of each variant calling pipeline add the resulting VCF into the configHash
	if ( defined $cfg_hash->{'vc_output_paths'}){
		$cfg_hash->{'vc_output_paths'} .= ",".$vcfFile;
		$cfg_hash->{'vc_algorithms'} .= ",freebayes";
	}else{
		$cfg_hash->{'vc_output_paths'} = $vcfFile;
		$cfg_hash->{'vc_algorithms'} = "freebayes";
	}

	Hence config variables are used to instruct the fast_final_output_steps to use more VCF files. 
	A new field will be added to the INFO field of each VCF file found but not for GATK. GATK calling will be always used as a primary tool.

	Currently two approaches are implemented to merge the results from different callers:
		UNION: In this case, after merge, the complete VCF is sorted and duplicate calls (where all samples have the same genotype) are removed.
			The GATK variant calling is always kept!
		ADDUNIQ: The variants that are present in any of the callers and are not in the result of GATK are only added to the result of GATK.

3. The union of variants from different callers increased the times for the computation of the final ouput of *VarGenius*. We reduced the time of the production of the output using the Annovar and VCF output.

4. The BAM files used in IGV are now those after the GATK Base Recalibrator step. This is the last file before the calling and is the one that must be used for the first check in IGV and the one that will be used for coverage statistics and for the CNV calling

5. VarGenius can now call CNVs using two different tools: ExomeDepth and XHMM. A new task is added to run this pipeline that can be automatically executed after the refinement task. 
	- ExomeDepth requires BAM files from as much as possible samples and will be run independently for sex and autosomes using two different jobs.
	- XHMM requires GATK DepthOfCoverage output that unfortunately will need different parameters from what we used for the coverage analysis, hence, for each BAM file it must be executed using a different job since it is computationally expensive. XHMM will be launched for sex and autosomes separately, too.

	The output from both CNV callers is preprocessed and merged. 

	a. Since the callers are used separately for sex chrosomes and autosomes new targets file must be produced. Hence, if the input parameter **--cnv_detection** is used, VarGenius will immediately produce two separated target files that will be named **TargetName_auto.bed** and **TargetName_sex.bed** and saved in the **target_reg_f** folder.

6. Given the panels of genes into the DATA/gene_panels/ folders the VarGenius output will indicate if the variated gene belongs to one of more panels
	The new function **get_disease_gene_association** returns an hash in which, for each gene there will be an association to panels of diseases.
	This function is called only once in **rearrange_rawtabular_out** (where the output is rearranged for simplicity of use) and the hash used through the annotation process.
	The panels to use must be written into **user_config.txt** (**disease_genes_lists_f** parameter).

7. Target Extension: The target file related with the enrichment kit may be enlarged. The user can choose how many nucleotides from the 3' an 5' using the parameter:
	*target_extension = YES* to say that the target must be enlarged 
	*target_extens_dim = 200* to say how many nucleotides

8. Automatic quality check. VarGenius compares now the coverage, the number of variants with the one obtained in experiment with the same kit. If these values are far from the mean the user is alerted in the result webpage.

9. Modified the generation of the outlist file. Now most of the output files to save are identified and written in the outlist.txt file even if the program is not executed. This permits to generate on demand the outlist file with all the requested output to save and store.

N.  Minor errors have been corrected: 
	- corrected a bug for which the IGV file was print without the BAM files paths.
	- paths_hash_to_html does not take in input the HOST name now. It is taken from the configuration hash.
	- corrected a bug for which the vcfimport was not started automatically
	- corrected a bug for which the non-covered region were not being created
	- corrected all the log file names for the jobs for all tasks that now are divided by log and error as requested for qsub

#April 2017

GATK updates to version 3.7
- variant_index_parameter and variant_index_type is no longer necessary when producing GVCF files
https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_engine_CommandLineGATK.php#--	variant_index_parameter

- Indel Realigner is no more useful
