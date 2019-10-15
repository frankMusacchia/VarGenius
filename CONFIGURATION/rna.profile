#Parameters for RNA sequencing Variant Call
#These are based on the Best Practices for RNA Seq from GATK
#https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

######################
alignment_prog = STAR

#TrimGalore
clip_R1 =
clip_R2 =

#ALIGNMENT
#Merges read1 and 2 in a single file
merge_pairs = NO
split_ncigar_exec = YES

#REFINEMENT (BestPractices GATK)
realign_target_creator = YES
indel_realigner = YES
base_recalibrator = YES
print_reads = YES

#VARIANT CALLING STEPS
#http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
#VQSR won't work well on gene panels because there's just not enough sites to build a stable model. 
varrecal = NO
apprecal = NO


#PICARD TOOLS PARAMETERS
#AddOrReplaceReadGroups parameters
picard_RGID = YES
picard_RGLB = YES
picard_RGPL = illumina
picard_RGPU = YES
picard_RGSM = YES
picard_RGCN = TIGEM
picard_SORT_ORDER = coordinate
picard_USE_THREADING = true

#MarkDuplicates
#MarkDuplicates Parameters
#http://gatkforums.broadinstitute.org/gatk/discussion/2799
#markduplicate is sufficient; HC will automatically apply filters to ignore dupes 
picard_CREATE_INDEX = true
picard_VALIDATION_STRINGENCY = SILENT
picard_REMOVE_DUPLICATES =
picard_METRICS_FILE = duplicate.metrics
picard_ASSUME_SORTED =
picard_MAX_FILE_HANDLES_FOR_READ_ENDS_MAP =

#SplitNCigarReads parameters
SNCR_read_filter = ReassignOneMappingQuality
SNCR_reassign_mapping_quality_from = 255
SNCR_reassign_mapping_quality_to = 60
SNCR_unsafe = ALLOW_N_CIGAR_READS

#GATK HaplotypeCaller (Old pipeline call  1000.0 and emit 200.0
HC_dontUseSoftClippedBases = YES
HC_stand_call_conf = 20.0
hapc_annotation =
#GenotypeGVCFs parameters
ggvcf_single_sam_stand_call_conf = 30

#VariantFiltration parameters
VF_clusterWindowSize = 35
VF_clusterSize = 3

#JOBS parameters
#Alignment (this task is ran per read file)
#mode: you can choose 'bwa' or 'bwa-kit'
align_mem = 80GB
align_walltime = 10:00:00

#Refinement (this task is ran per read file)
refine_ncpus = 8
refine_mem = 80GB

