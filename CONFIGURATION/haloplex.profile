#Parameters for Haloplex

######################



qsub_walltime = 10:00:00
#Steps that must be avoided
mark_rem_dup_exec = NO
map_alt_exec = NO
sam_sort_idx = YES
#Refinement
#When using Haloplex system I had problems
#with BaseRecalibrator, hence I am jumping this step
#http://gatkforums.broadinstitute.org/gatk/discussion/4272/targeted-sequencing-appropriate-to-use-baserecalibrator-bqsr-on-150m-bases-over-small-intervals
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
HC_stand_emit_conf = 10
hapc_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,HaplotypeScore,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest,SpanningDeletions

#GenotypeGVCFs parameters
ggvcf_single_sam_stand_call_conf = 30
ggvcf_single_sam_stand_emit_conf = 10

#GATK Variant Recalibrator
vr_max_gaussians = 8

#This is a link to the working directory
html_host = http://UDProject:udpNGSd4t4@130.186.13.104/tigem

#EMAIL TO SEND WHEN ANALYSIS IS COMPLETED
email_author = VarGenius@tigem.it
email_recipients = francescomusacchia@gmail.com,annalaura.torella@gmail.com,msavarese@tigem.it
email_subject = VarGenius analysis is completed

#Refinement (this task is ran per read file)
refine_ncpus = 8
refine_mem = 80GB
