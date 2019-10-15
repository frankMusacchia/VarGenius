#Parameters for Whole Genome Analysis

######################

#Steps of VarGenius2
#Freebayes is still not ready
consensus_vc_approach = NO

#Quality Check (this task is ran per read file)
qc_mem = 30GB
qc_walltime = 10:00:00
#Alignment (this task is ran per read file)
align_mem = 30GB
align_walltime = 10:00:00
#Refinement (this task is ran per read file)
refine_mem = 70GB
refine_walltime = 24:00:00

#Java Memory is the standard used for an whole exome
#Java maximum memory for GATK Programs
javamem_SNCR = Xmx25g
javamem_RTC = Xmx4g
javamem_IR = Xmx4g
javamem_BR = Xmx8g
javamem_PR = Xmx4g
javamem_MD = Xmx4g
javamem_MSF = Xmx4g
javamem_HAPC = Xmx64g
javamem_CV = Xmx100g
javamem_GGVCF = Xmx25g
javamem_VR = Xmx64g
javamem_VF = Xmx64g
javamem_PBT = Xmx8g
javamem_DOC = Xmx64g
javamem_SV = Xmx8g
javamem_CV = Xmx8g
javamem_VA = Xmx8g
javamem_CGP = Xmx8g
javamem_SS = Xmx4g
javamem_BBI = Xmx4g
