#Parameters to change to have the most recently update
#GATK Best Practices Pipeline in GATK4
#GATK uses Picard to align, sort and mark duplicates after the alignment
###################

#In GATK4 annotations for HaplotypeCaller and GenotypeGVCFs change. 
#HaplotypeScore and SpanningDeletion are now not applicable
hapc_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,InbreedingCoeff,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest
ggvcf_annotation = Coverage,FisherStrand,BaseQualityRankSumTest,InbreedingCoeff,MappingQualityRankSumTest,MappingQualityZero,QualByDepth,RMSMappingQuality,ReadPosRankSumTest

#BaseRecalibrator parameters
known_BR = known_db_snp

known_db_snp = dbsnp_138.hg19.vcf
known_sites_db_snp = dbsnp_138.hg19.vcf

print_reads_prog = ApplyBQSR
combvar_prog = MergeVcfs
baq_PR =
