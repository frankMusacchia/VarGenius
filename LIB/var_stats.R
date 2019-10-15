
#R CMD BATCH --no-save --no-restore '--args <output_from_vargenius> <samples_separate_by_comma>  <out_file_name> ' var_stats.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)


#Output from VarGenius
var_in <- args[1] #"NA01/finalout_out/NA01.hg19.final_out_resgr14.txt"#
var_in 
#List of samples separated by comma
samples_str <- args[2] #c("UD15","UD16","UD17") #
samples_str
#Output file name
out_name <- args[3] #"chr.buttami"#
out_name
#Analysis
analysis <- args[4]
analysis
#DBSNP field in the output
dbsnp_fld <- args[5]#'avsnp150'
dbsnp_fld
#CADD version
cadd_phred <- args[6] #= "cadd13_phred"
#
hgmdgene_h <- "GENE_hgmd"
exac_max_h <- "exac_max"
exac_all_h <- "exac_all"

#CADD threshold
	#Usually, a scaled CADD score of 20 means that a variant is amongst the top 1% of deleterious variants 
	#in the human genome. A scaled CADD score of 30 means that the variant is in the top 0.1% and so forth.
cadd_thr = 20

#Fields to access
tot_reads_fld = "tot_reads"
zyg_fld = "ZYG"
vg_separator = "_"
perc_alt_reads_fld = "PercALT_reads"



#Get the output
v<-read.delim(var_in,sep="\t",head=T,stringsAsFactors=F,check.names=FALSE)
 
#Split string with sample names 
samples_l <- unlist(strsplit(samples_str,",",fixed=T))
samples_l[[1]] 
#Filter the X chromosome variants
v_X<-v[v$chrom == 'chrX',]
dim(v_X)
#Filter the Y chromosome variants
v_Y<-v[v$chrom == 'chrY',]
dim(v_Y)
#Number of probands is the number of samples minus the parents
prob_num <- length(samples_l)-2 


#####################################################
#################NEEDED FUNCTIONS####################
#####################################################

#This says if the nucleotide change is a transition or a trasversion
snv_type = function (nucl_1,nucl_2){
	nucl_1 = toupper(nucl_1)
	nucl_2 = toupper(nucl_2)
	
	type = "-"
	
	if ( nucl_1 == nucl_2){
			type = "not_a_SNV"
	}
	#Transversion: A,G <-> T,C 
	if (  ( (nucl_1 == 'A' | nucl_1 == 'G') & (nucl_2 == 'C' | nucl_2 == 'T') ) |
				( (nucl_1 == 'T' | nucl_1 == 'C') & (nucl_2 == 'A' | nucl_2 == 'G') ) ){
		type = "transversion"		
	}
	
	#Transversion: A <-> G or T <-> C
	if (  ( (nucl_1 == 'A' & nucl_2 == 'G') | (nucl_1 == 'G' & nucl_2 == 'A') ) |
				( (nucl_1 == 'T' & nucl_2 == 'C') | (nucl_1 == 'C' & nucl_2 == 'T') ) ){
		type = "transition"		
	}
	type
}


analysis_stats <- function(vgout) {
	n_res <- nrow(vgout)#num resuts
	#Set all non present values to NA
	vgout$qual[vgout$qual=='-']<- NA
	mean_q <- mean(na.omit(as.numeric(vgout$qual))) #mean quality
		 
	#synonymous and nonsynonymous SNV and INDELs
	snp <- vgout[nchar(vgout$ref) == 1 & nchar(vgout$alt) == 1,]
	snp_n <- nrow(snp)
	snp_perc <- (snp_n/n_res)*100
  
	indel <- vgout[nchar(vgout$ref) > 1 | nchar(vgout$alt) > 1,]	
	indel_n <- nrow(indel)
	indel_perc <- (indel_n/n_res)*100
	
	syn <- grepl('synonymous SNV',vgout$func_refgene)
	syn_n <- length(syn[syn==TRUE])
	syn_perc <- (syn_n/n_res)*100
	
	nonsyn <- grepl('nonsynonymous SNV',vgout$func_refgene)
	nonsyn_n <- length(nonsyn[nonsyn==TRUE])
	nonsyn_perc <- (nonsyn_n/n_res)*100
	ns_syn_ratio <- -1
	if (nonsyn_n > 0 && syn_n > 0){
		ns_syn_ratio <- nonsyn_n/syn_n
	}

	#Check for denovo variants if proband parents are both present
	p_zyg = paste(analysis,"P","ZYG",sep="_")
	m_zyg = paste(analysis,"M","ZYG",sep="_")
	f_zyg = paste(analysis,"F","ZYG",sep="_")	
	if ( m_zyg %in% colnames(vgout) &  f_zyg %in% colnames(vgout)) {
		#Get the DENOVO
		denovo <- nrow(vgout[vgout[[p_zyg]] !='HOMREF' & vgout[[p_zyg]] != '-'  & vgout[[m_zyg]]  =='HOMREF' & vgout[[f_zyg]]  =='HOMREF',])
		denovoperc <- (denovo/n_res)*100
	}else{
		denovo <- '-1'
		denovoperc <- '-1' 
	}	
		
	#Transitions and transversions among SNVs
	snp$snptype = mapply(snv_type,snp$ref,snp$alt)
	transit <- grepl('transition',snp$snptype)
	transit_n <- length(transit[transit==TRUE])
	transit_perc <- (transit_n/n_res)*100
	transv <- grepl('transversion',snp$snptype)
	transv_n <- length(transv[transv==TRUE])
	transv_perc <- (transv_n/n_res)*100
	tr_tv_ratio <- -1
	if (transit_n > 0 & transv_n > 0){
		tr_tv_ratio <- transit_n/transv_n
	}
		
	#Frequencies and dbsnp
	#Here I remove both those <0 because VarGenius uses -1 when data is missing and also I remove the NA
	vgout$exac_max[vgout$exac_max=='-']<- NA
	exac_max_mean <- mean(as.numeric(vgout$exac_max[as.numeric(vgout$exac_max)>0]),na.rm=TRUE)
	if ( is.na(exac_max_mean) ){
		exac_max_mean <- -1
	}
	vgout$exac_all[vgout$exac_all=='-']<- NA
	exac_all_mean <- mean(as.numeric(vgout$exac_all[as.numeric(vgout$exac_all)>0]),na.rm=TRUE)
	if ( is.na(exac_all_mean) ){
		exac_all_mean <- -1
	}	
	n_dbsnp <- nrow(vgout[vgout[[dbsnp_fld]] != '-',])
	n_dbsnp_perc <- (n_dbsnp/n_res)*100
	n_notdbsnp <- nrow(vgout[vgout[[dbsnp_fld]] == '-',])
	n_notdbsnp_perc <- (n_notdbsnp/n_res)*100
	#HGMD
	n_hgmd <- nrow(vgout[vgout[[hgmdgene_h]] != '-',])
	n_hgmd_perc <- (n_hgmd/n_res)*100
	#CADD PHRED SCORE > Threshold 
	#Usually, a scaled CADD score of 20 means that a variant is amongst the top 1% of deleterious variants 
	#in the human genome. A scaled CADD score of 30 means that the variant is in the top 0.1% and so forth.
	cadd_del_n <- nrow(na.omit(vgout[as.numeric(vgout[[cadd_phred]])>cadd_thr,]))
	cadd_del_n_perc <-(cadd_del_n/n_res)*100
	#Splicing
	spl_var <-	sum(!is.na(as.numeric(vgout$refgene_splvar_dist_from_exon))) #splicing variants
	spl_var_perc <- (spl_var/n_res)*100			
	
  output <- list(n_res,snp_n,snp_perc,indel_n,indel_perc,mean_q,denovo,denovoperc,syn_n,syn_perc,nonsyn_n,nonsyn_perc,ns_syn_ratio,transit_n,transit_perc,
  transv_n,transv_perc,tr_tv_ratio,exac_max_mean,exac_all_mean,n_dbsnp,n_dbsnp_perc,n_notdbsnp,n_notdbsnp_perc,n_hgmd,n_hgmd_perc,cadd_del_n,
  cadd_del_n_perc,spl_var,spl_var_perc)
  
  return(output)
}

##################################################
## BLOCK1: Reads on X and Y chromosomes
##################################################

#Delimitation of the region with X Y reads count
delim_1 = "##############Total Reads on X and Y chromosomes\n"
cat(delim_1,file=out_name)

header = c("SampleName","tot_reads_X","tot_reads_Y","perc_reads_X","perc_reads_Y")

tab = "\t"
cat(header,file=out_name,sep=tab,append = TRUE)

#Analyze and print statistics about chromosomes
for ( sample in samples_l ){
	#sam <- samples_l[1]
	#sam
	cat("\n",file=out_name,append = TRUE)
	#paste(sam,tot_reads_fld,sep=vg_separator)
	#Total reads on X chromosome 
	v_X_f <- as.numeric(v_X[[paste(sample,tot_reads_fld,sep=vg_separator)]])
	v_X_f[is.na(v_X_f)] <- 0
	X_sum <- sum( v_X_f )
	
	#Total reads on Y chromosome 
	v_Y_f <- as.numeric(v_Y[[paste(sample,tot_reads_fld,sep=vg_separator)]])
	v_Y_f[is.na(v_Y_f)] <- 0
	Y_sum <- sum( v_Y_f )		
	
	sum_xy = X_sum + Y_sum
	
	cat(sample,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)
	cat(X_sum,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)
	cat(Y_sum,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)
	cat(paste(round(X_sum/sum_xy, 4)*100, "%", sep=""),file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)
	cat(paste(round(Y_sum/sum_xy, 4)*100, "%", sep=""),file=out_name,append = TRUE)	
	cat("\n",file=out_name,append = TRUE)
}
end_del = "##############\n"
cat(end_del,file=out_name,append = TRUE)


##################################################
## BLOCK2: Zygosity Frequencies with X chromosome
##################################################

#Delimitation of the region with X Y reads count
delim_2 = "##############Zygosity Frequencies with X chromosome\n"
cat(delim_2,file=out_name,append = TRUE)
#Analyze and print statistics about chromosomes
for ( sample in samples_l ){
#	cat(paste("Frequencies for:",sample,sep=" "),file=out_name,append = TRUE)
#	cat("\n",file=out_name,append = TRUE)
	
	#Frequency of each zygosity in chromosome X
	cat(paste("Chromosome X:",sample,sep=" "),file=out_name,append = TRUE)
	cat("\n",file=out_name,append = TRUE)
	X_freq <- ftable(v_X[[paste(sample,zyg_fld,sep=vg_separator)]])
	write.table(as.data.frame(X_freq),file=out_name,append = TRUE,quote=F,row.names=F,col.names=F)
	cat("\n",file=out_name,append = TRUE)
	
#	#Frequency of each zygosity in chromosome Y
#	cat(paste("Chromosome Y:",sample,sep=" "),file=out_name,append = TRUE)
#	cat("\n",file=out_name,append = TRUE)
#	Y_freq <- ftable(v_Y[[paste(sample,zyg_fld,sep=vg_separator)]])
#	write.table(as.data.frame(Y_freq),file=out_name,append = TRUE,quote=F,row.names=F,col.names=F)
#	cat("\n",file=out_name,append = TRUE)
#	cat("\n",file=out_name,append = TRUE)
}
cat(end_del,file=out_name,append = TRUE)
min_gq = 50
altp_thr1 = 0.25
altp_thr2 = 0.8
tab = "\t"


##################################################
##BLOCK3:  SEGREGATION
##################################################
#Prints the table of segregation	
print_seg_table <- function (mat,dimx,dimy,zygosities,comp1,comp2){
	
	#Print the header
	cat(paste(comp1,comp2,sep="/"),file=out_name,append = TRUE)
	cat("\n",file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)		
	for (z in zygosities){
		cat(z,file=out_name,append = TRUE)
		cat(tab,file=out_name,append = TRUE)		
	}
	cat("\n",file=out_name,append = TRUE)	
	
	#Print the matrix
	for ( x in 1:dimx){
		#Print Row names
		cat(zygosities[x],file=out_name,append = TRUE)
		cat(tab,file=out_name,append = TRUE)
		for (y in 1:dimy){
			cat(mat[x,y],file=out_name,append = TRUE)
			cat(tab,file=out_name,append = TRUE)			
		}
		cat("\n",file=out_name,append = TRUE)		
	}
}	


	#Delimitation of the region with Segregation
	delim = "##############Segregation\n"
	cat(delim,file=out_name,append = TRUE)

	moth_in = prob_num + 1#The mother is the first one after the probands
	fath_in = prob_num + 2#The father is the second one after the probands

	#Delimitation of the region with Segregation
	delim = "Segregation test tables\n"
	cat(delim,file=out_name,append = TRUE)

	zygosities = c("HOMREF","HET","HOM","-")
	#Define a matrix to fill
	mat <- matrix(NA, nrow = length(zygosities), ncol = length(zygosities))

#Do the segregation only if there is a trio or quartet
if (length(samples_l)==3 | length(samples_l)==4){		
	#for each proband print the segregation table with mother and father	
	for (prob_in in 1:prob_num){

		prob_zyg = paste(samples_l[prob_in],"ZYG",sep="_")
		moth_zyg = paste(samples_l[moth_in],"ZYG",sep="_")
		fath_zyg = paste(samples_l[fath_in],"ZYG",sep="_")
		
	#	##########################
	#	#Mother and Proband Table
	#	##########################
		rowin <- 1
		colin <- 1
		for (z1 in zygosities){
			for (z2 in zygosities){
				mat[rowin,colin] <- nrow(v[v[[prob_zyg]] == z1 & v[[moth_zyg]] == z2,])
				colin = colin +1
			}
			rowin = rowin +1
			colin = 1
		}
		print_seg_table(mat,length(zygosities),length(zygosities),zygosities,samples_l[prob_in],samples_l[moth_in])
		
	#	##########################
	#	#Father and Proband Table
	#	##########################	
		rowin <- 1
		colin <- 1
		for (z1 in zygosities){
			for (z2 in zygosities){
				mat[rowin,colin] <- nrow(v[v[[prob_zyg]] == z1 & v[[fath_zyg]] == z2,])
				colin = colin +1
			}
			rowin = rowin +1
			colin = 1
		}
		print_seg_table(mat,length(zygosities),length(zygosities),zygosities,samples_l[prob_in],samples_l[fath_in])

	}

	#	##########################
	#	#Mother and Father Table
	#	##########################	
		rowin <- 1
		colin <- 1
		for (z1 in zygosities){
			for (z2 in zygosities){
				mat[rowin,colin] <- nrow(v[v[[moth_zyg]] == z1 & v[[fath_zyg]] == z2,])
				colin = colin +1
			}
			rowin = rowin +1
			colin = 1
		}
		print_seg_table(mat,length(zygosities),length(zygosities),zygosities,samples_l[moth_in],samples_l[fath_in])
}else{
	#Delimitation of the region with Segregation
	delim = "Segregation test has not been performed. VarGenius computes it only for TRIOs or QUARTETs\n"
	cat(delim,file=out_name,append = TRUE)
}

end_del = "##############\n"
cat(end_del,file=out_name,append = TRUE)




############################################
# BLOCK4: PRINT NUMBER OF EXISTING VARIANTS IN DBSNP#
############################################
#Delimitation of the region with dbSNP existing variants number
delim = "##############dbSNP Existent variants\n"
cat(delim,file=out_name,append = TRUE)


dbsnp_exist <- nrow(v[v[[dbsnp_fld]] != '-',])
dbsnp_nonexist <- nrow(v[v[[dbsnp_fld]] == '-',])

cat("dbSNP existing variants: ",file=out_name,append = TRUE)
cat(dbsnp_exist,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)

cat("dbSNP non existing variants: ",file=out_name,append = TRUE)
cat(dbsnp_nonexist,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)

#Optional to be commented..
sample = samples_l[[1]]
palt_fld = paste(sample,"PercALT_reads",sep="_")
refr_fld = paste(sample,"REF_reads",sep="_")
zyg_fld = paste(sample,"ZYG",sep="_")
gq_fld = paste(sample,"GQ",sep="_")

dbsnp_nonexist <- v[v[[dbsnp_fld]] == '-',]	
dsnp_ne_gq <- dbsnp_nonexist[dbsnp_nonexist[[gq_fld]] != '-' ,]
a <- mean(as.numeric(dsnp_ne_gq[[gq_fld]]))
dsnp_ne_rr <- dbsnp_nonexist[dbsnp_nonexist[[refr_fld]] != '-' ,]
b <- mean(as.numeric(dsnp_ne_rr[[refr_fld]]))
dsnp_ne_ar <- dbsnp_nonexist[dbsnp_nonexist[[palt_fld]] != '-' ,]
c <- mean(as.numeric(dsnp_ne_ar[[palt_fld]]))
dsnp_ne_ex <- dbsnp_nonexist[dbsnp_nonexist[['exac_all']] != '-1' & dbsnp_nonexist[['exac_all']] != '-' ,]
d <- mean(as.numeric(dsnp_ne_ex[['exac_all']]))

ne_dbSNP_stats_str = paste("Among the non-existing in dbSNP:",paste("Mean GQ:",a,"Mean PercREFReads:",b,"Mean PercALTReads:",c,"Mean ExAC:",d,sep=" "),sep="\n")
cat(ne_dbSNP_stats_str,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)

############################################
# BLOCK5: ZYGOSITY CONCORDANCE
############################################

#Delimitation of the region with Zygosity concordance
delim_3 = "##############Zygosity concordance\n"
cat(delim_3,file=out_name,append = TRUE)

header2 = paste("Sample","THE","SBHE","ONHE","THOR","SBHOR","ONHOR","THO","SBHO","ONHO","FHOHE","OFHOHE",sep=tab)
cat(header2,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)

#Print the statistics for each sample
for ( sample in samples_l ){
	cat(sample,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)
	palt_fld = paste(sample,"PercALT_reads",sep="_")
	palt_fld
	zyg_fld = paste(sample,"ZYG",sep="_")
	zyg_fld
	gq_fld = paste(sample,"GQ",sep="_")
	gq_fld
	
	#####################HET
	#TRUE HET:
	true_het <- nrow(v[v[[palt_fld]]>altp_thr1 & v[[palt_fld]]<altp_thr2 & v[[zyg_fld]]=='HET',c(palt_fld,zyg_fld)])
	cat(true_het,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#SHOULD BE HET
	shbe_het <- nrow(v[v[[palt_fld]]>altp_thr1 & v[[palt_fld]]<altp_thr2 & v[[zyg_fld]]!='HET',c(palt_fld,zyg_fld)])
	cat(shbe_het,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#OK_NOT_HET
	oknot_het <- nrow(v[v[[palt_fld]]>altp_thr1 & v[[palt_fld]]<altp_thr2 & v[[zyg_fld]]!='HET' & v[[gq_fld]]<min_gq,c(palt_fld,zyg_fld)])
	cat(oknot_het,file=out_name,append = TRUE)
	cat(" ",file=out_name,append = TRUE)
	cat(paste(round(oknot_het/shbe_het, 2)*100, "%", sep=""),file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#####################HOMREF
	#TRUE HOMREF:
	true_homref <- nrow(v[v[[palt_fld]]<altp_thr1 & v[[zyg_fld]]=='HOMREF',c(palt_fld,zyg_fld)])
	cat(true_homref,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#SHOULD BE HOMREF:
	shbe_homref <- nrow(v[v[[palt_fld]]<altp_thr1 & v[[zyg_fld]]!='HOMREF',c(palt_fld,zyg_fld)])
	cat(shbe_homref,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#OK NOT HOMREF:
	oknot_homref <- nrow(v[v[[palt_fld]]<altp_thr1 & v[[zyg_fld]]!='HOMREF' & v[[gq_fld]]<min_gq,c(palt_fld,zyg_fld)])
	cat(oknot_homref,file=out_name,append = TRUE)
	cat(" ",file=out_name,append = TRUE)
	cat(paste(round(oknot_homref/shbe_homref, 2)*100, "%", sep=""),file=out_name,append = TRUE)	
	cat(tab,file=out_name,append = TRUE)


	######################HOM
	#TRUE HOM:
	true_hom <- nrow(v[v[[palt_fld]]>altp_thr2 & v[[zyg_fld]]=='HOM',c(palt_fld,zyg_fld)])
	cat(true_hom,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#SHOULD BE HOM:
	shbe_hom <- nrow(v[v[[palt_fld]]>altp_thr2 & v[[zyg_fld]]!='HOM',c(palt_fld,zyg_fld)])
	cat(shbe_hom,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#OK NOT HOM:
	oknot_hom <- nrow(v[v[[palt_fld]]>altp_thr2 & v[[zyg_fld]]!='HOM' & v[[gq_fld]]<min_gq,c(palt_fld,zyg_fld)])
	cat(oknot_hom,file=out_name,append = TRUE)
	cat(" ",file=out_name,append = TRUE)
	cat(paste(round(oknot_hom/shbe_hom, 2)*100, "%", sep=""),file=out_name,append = TRUE)	
	cat(tab,file=out_name,append = TRUE)


	##########FALSE HOM AND HET
	#FALSE HOM or HET:
	false_homhet <- nrow(v[v[[palt_fld]]==0 & (v[[zyg_fld]]=='HOM' | v[[zyg_fld]]=='HET') ,c(palt_fld,zyg_fld)])
	cat(false_homhet,file=out_name,append = TRUE)
	cat(tab,file=out_name,append = TRUE)

	#OK FALSE HOM or HET:
	ok_false_homhet <- nrow(v[v[[palt_fld]]==0 & (v[[zyg_fld]]=='HOM' | v[[zyg_fld]]=='HET') & v[[gq_fld]]<min_gq ,c(palt_fld,zyg_fld)])
	cat(ok_false_homhet,file=out_name,append = TRUE)
	cat(" ",file=out_name,append = TRUE)
	cat(paste(round(ok_false_homhet/false_homhet, 2)*100, "%", sep=""),file=out_name,append = TRUE)	
	cat("\n",file=out_name,append = TRUE)
}
cat(end_del,file=out_name,append = TRUE)


############################################
# BLOCK6: ZYGOSITIT CONCORDANCE legendS
############################################


#CONCORDANCE BETWEEN READS NUMBER AND GENOTYPE CALLED BY GATK
#Delimitation of the region with X Y reads count
delim = "##############Zygosity concordance (legenda)\n"
cat(delim,file=out_name,append = TRUE)

####################LEGEND
#####################HET
true_het_desc = paste("(THE) TRUE HET:Variants with perc of alt reads between",altp_thr1,"and", altp_thr2,"and are claimed as HET",sep=" ")
cat(true_het_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
shbe_het_desc = paste("(SBHE) SHOULD BE HET: Variants with perc of alt reads between",altp_thr1,"and", altp_thr2,"and are NOT claimed as HET",sep=" ")
cat(shbe_het_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
oknot_het_desc = paste("(ONHE) OK_NOT_HET: Variants with perc of alt reads between",altp_thr1,", ", altp_thr2,"NOT claimed as HET but with GQ<",min_gq,sep=" ")
cat(oknot_het_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
#####################HOMREF
true_homref_desc = paste("THOR) TRUE HOMREF:Variants with perc of alt reads lower than",altp_thr1,"claimed as HOMREF",sep=" ")
cat(true_homref_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
shbe_homref_desc = paste("(SBHOR) SHOULD BE HOMREF:Variants with perc of alt reads lower than",altp_thr1,"NOT claimed as HOMREF",sep=" ")
cat(shbe_homref_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
oknot_homref_desc = paste("(ONHOR) OK NOT HOMREF:Variants with perc of alt reads lower than",altp_thr1,"NOT HOMREF but with GQ<",min_gq,sep=" ")
cat(oknot_homref_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
######################HOM
true_hom_desc = paste(" (THO) TRUE HOM:Variants with perc of alt reads greater than",altp_thr2,"claimed as HOM",sep=" ")
cat(true_hom_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
shbe_hom_desc = paste("(SBHO) SHOULD BE HOM:Variants with perc of alt reads greater than",altp_thr2,"NOT claimed as HOM",sep=" ")
cat(shbe_hom_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
oknot_hom_desc = paste("(ONHO) OK NOT HOM:Variants with perc of alt reads greater than",altp_thr2,"NOT HOM but with GQ<",min_gq,sep=" ")
cat(oknot_hom_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
##########FALSE HOM AND HET
false_homhet_desc = paste("(FHOHE) FALSE HOM or HET:Variants with perc of alt = 0 claimed as HOM or HET",sep=" ")
cat(false_homhet_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
ok_homhet_desc = paste("(OFHOHE) OK FALSE HOM or HET:Variants with perc of alt = 0 claimed as HOM or HET but with GQ<",min_gq,sep=" ")
cat(ok_homhet_desc,file=out_name,append = TRUE)
cat("\n",file=out_name,append = TRUE)
cat(end_del,file=out_name,append = TRUE)


############################################
# BLOCK7: Variants statistics
# This block of information is highlighted with ##DB since it 
# can be loaded into the database into the table with the name following ##DB 
# hence analysis_statistics.
# The fields name are the same in the table of db_schema.sql Hence to modify here
# means also to modify there.
############################################
delim = "##############DB db_analysis_statistics\n"
cat(delim,file=out_name,append = TRUE)

#STATISTICS PER-ANALYSIS
#N.B. THIS FIELDS MUST BE IDENTICAL TO FIELDS OF THE analysis_statistics TABLE 
#INTO THE SCRIPT db_schema.sql
header_a = c("totvariants","totsnvs","snvsperc","totindels","indelsperc","meanqual","denovo","denovoperc","synsnv","synsnvperc","nonsynsnv","nonsynsnvperc","ns_sratio","snvtr","snvtrperc","snvtv","snvtvperc","tr_tvratio","meanexacmax","meanexacall","dbsnpe","dbsnpeperc","dbsnpne","dbsnpneperc","inhgmd","inhgmdperc","caddgt20","caddgt20perc","splvar","splvarperc")
cat(header_a,file=out_name,sep=tab,append = TRUE)
cat("\n",file=out_name,sep=tab,append = TRUE)
					
a_stats <- analysis_stats(v)
cat(unlist(a_stats),file=out_name,append = TRUE,sep=tab)

cat("\n",file=out_name,append = TRUE)
cat(end_del,file=out_name,append = TRUE)
	

