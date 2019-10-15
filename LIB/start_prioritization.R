#Configuration script for the prioritization


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

#Analysis name
analysis_name <- args[1] #"UD_NA001" # 
analysis_name
ver = "v1"

# working dir, where input and output files are placed
wdir = args[2] #Dir where the games are played #'/pico/work/TELET_UDP/PROVE/prior_michele/' # 
wdir
# it may be Tab-separated (.tsv), or R-data (.Rdata )
datafile = args[3] #Output file in Rdata format #"UD_NA001_TRIO_2ndlev.RData" 
datafile
#Optional. comment if you want to avoid use
pedfile = args[4]#Pedfile path  # # 
pedfile
gene_table_file =  args[5] #Gene annotation Rdata file  # "/pico/work/TELET_UDP/PROVE/prior_michele/UDP_gene_full_table.ver02.RData" #
gene_table_file
#annotation and panels  absolute paths
candidate_gene_dir = args[6] #paste(startdir,"gene_annotations",sep="/") #
candidate_gene_dir
annotation_gene_dir = args[7] #paste(startdir,"gene_annotations",sep="/") #
annotation_gene_dir
candidate_variants_dir = args[8] # paste(startdir,"variant_annotations",sep="/") #
candidate_variants_dir
annotation_variants_dir = args[9] # paste(startdir,"variant_annotations",sep="/") #
annotation_variants_dir

#columns from the output files
chromosome_cname = args[10] # 'chrom' #
chromosome_cname
position_cname = args[11] # 'pos' #
position_cname
ref_cname = args[12] #'ref' #
ref_cname
alt_cname = args[13] #'alt' #
alt_cname

rs_cname = args[14] #'avsnp147' #
rs_cname
gene_cname = args[15] #'gene_refgene' #
gene_cname
cds = args[16] # "refgene_nucl_change" #
cds
aa = args[17] # "refgene_aa_change" #
aa
# qual_cname = args[19] #'QUAL'
# depth_cname = args[20] #'DEPTH'

#source script 
prioritization_script = args[18] #  "/pico/work/TELET_UDP/PROVE/prior_michele/prioritization.1.12.r"
prioritization_script

################################
#THESE VARIABLES MUST BE MODIFIED INTO THIS CONFIGURATION FILE
id_cname='key'
maxfrequency_cname='Freq_Max'
varianteffect_cname='GeneFunctionClass'

#R objects
effects<-c('CodingChange', 'LossOfFunction')   # ==
cdiseasevar=c('ClinVar_CLASS','HGMD_CLASS')  # grep
diseasevar_grep=list()
diseasevar_grep[['ClinVar_CLASS']]=c('*patho*')
diseasevar_grep[['HGMD_CLASS']]=c('*DM*')
phenotype=c("hgmd_phenotype")
pubmed_terms=list(c("intellectual disability"))

#FRANK: All the following could be inserted into the VarGenius Config
GTsuff='_GT'
# filtering options
Fmax=0.001
#Transmission_model=c("Dominant","Recessive","X-linked")[1]
#Penetrance=1
penetranceTollerance_dominantControls=0
penetranceTollerance_xlinkedControls=0
penetranceTollerance_recessiveControls=0

# Summary
summarycolumns<-c("affected","fathers","mothers","controls","HGVS",aa,maxfrequency_cname,rs_cname,"freq_factors","vmodel")

# whether consider the failed genotype of controls as homozygous reference
if_lax_segregation=FALSE  

#To be executed:
if_bestHit_gene_report=TRUE
if_gene_report = TRUE
if_write_vartable = TRUE
if_pubmed_besthits = TRUE
if_pubmed_others = TRUE
if_search_panel = TRUE
if_write_00 = TRUE
if_write_R = TRUE
#Leave false the following
if_panel_app=FALSE

######################
#Start prioritization
source(prioritization_script)



