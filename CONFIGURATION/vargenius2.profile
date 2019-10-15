###################
#This is the configuration we used for VarGenius2.0
#If you want to use freebayes, samtools and mix results
#Indicate the type of merge: UNION, INTERSECT, ADDUNIQ (this last means GATK+Unique in all other approches)
#or NO if you do not want to use multiple variant calling programs
consensus_vc_approach = UNION
gatk_ver = 4
target_extension = YES
target_extens_dim = 100

VCF_filter_vcf = NO

gatk_path = /pico/work/TELET_TIGEM/bin/gatk-4.1.3.0/gatk
#GATK special parameters
ggvcf_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,InbreedingCoeff,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest

hapc_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,InbreedingCoeff,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest
#When using a VCF from old version of GATK in GATK4, it should allow old MQ scores (DEPRECATED)
ggvcf_allow_old_rms = YES

#Genotyping (this task is ran per-group)
genot_mem = 80GB
genot_walltime = 100:00:00
qsub_walltime_max = 100:00:00

#CNV discovery
cnv_queue =

#GATK4 programs new names
apprecal_prog = ApplyVQSR
print_reads_prog = ApplyBQSR

#Needed sites for Variant Quality Score Recalibration
known_sites_mills = Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
known_sites_hapmap = hapmap_3.3.hg19.sites.vcf
known_sites_1000g = 1000G_phase1.snps.high_confidence.hg19.sites.vcf
known_sites_omni = 1000G_phase1.indels.hg19.sites.vcf
known_sites_1000gP3 = 1000G_phase3_v4_20130502.VG.sites.vcf.gz
#DBSNP The original file contains chromosomes number without "chr" string. I added it here
known_sites_db_snp = dbsnp_150_last.hg19.vcf.gz
