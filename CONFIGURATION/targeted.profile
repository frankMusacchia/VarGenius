#Parameters for panels of amplicons

######################

#Steps that must be avoided
mark_rem_dup_exec = NO
#Refinement
#When using Haloplex system I had problems
#with BaseRecalibrator, hence I am jumping this step
#http://gatkforums.broadinstitute.org/gatk/discussion/4272/targeted-sequencing-appropriate-to-use-baserecalibrator-bqsr-on-150m-bases-over-small-intervals
realign_target_creator = NO
indel_realigner = NO
base_recalibrator = NO
print_reads = NO
#MERGE READFILES 
mergesam = YES
mrdup_groups = NO
#VARIANT CALLING STEPS
#http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
#VQSR won't work well on gene panels because there's just not enough sites to build a stable model. 
varrecal = NO
apprecal = NO


#GATK HaplotypeCaller (Old pipeline call  1000.0 and emit 200.0
HC_stand_call_conf = 30
hapc_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,HaplotypeScore,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest,SpanningDeletions

#GenotypeGVCFs parameters
ggvcf_single_sam_stand_call_conf = 30

#GATK Variant Recalibrator
vr_max_gaussians = 8

#Refinement (this task is ran per read file)
refine_ncpus = 8
refine_mem = 80GB
