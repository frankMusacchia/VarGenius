#R CMD BATCH --no-save --no-restore '--args <vgout1_f> <vgout2_f> <analysis>' r_utils.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

length(args)
if (length(args) < 5){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <table1> <table2> <selcol> <outname>' r_utils.R")
}

###################COMPARISON OF TWO OUTPUT FROM VARGENIUS
#The samples names must be the analysis name concatenating P, M and F
#This part after the header so that we have a unique table!!!
#SAMPLE INFO TO MODIFY


anal_type = "TRIO"
analysis = "UD_MO017_TRIO"
analysis2 = "UD_WGS_MO017_TRIO"
vgout1_f = paste("/pico/work/TELET_TIGEM/vargenius_analyses",analysis,"finalout_out",paste(analysis,"hg19.final_out_resgr14.txt",sep="."),sep="/");
vgout1_f
vgout2_f = paste("/pico/work/TELET_TIGEM/vargenius_analyses",analysis2,"finalout_out",paste(analysis2,"hg19.final_out_resgr14.txt",sep="."),sep="/");
vgout2_f


vgout1_f = args[1]
vgout1_f
vgout2_f = args[2] 
vgout2_f
analysis = args[3]
analysis2 = args[4]
anal_type = args[5]
selected_genes_f = args[6]

proband = "P"

	
dbsnp_h <- "avsnp150"
hgmdgene_h <- "GENE_hgmd"
exac_max_h <- "exac_max"
exac_all_h <- "exac_all"
cadd_phred = "cadd13_phred"
cadd_thr = 20

#HEADER
out_name_a = paste(analysis,"comp.txt",sep="_")
out_name_s = paste(analysis,"samples_comp.txt",sep="_")

  
header_a_OLD = c("Analysis","TotVar","TotSNVs","SNVsPerc","TotINDELs","INDELsPerc","Mean_qual","SynSNV","SynSNVPerc","NonSynSNV","NonSynSNVPerc",
					"NS/S","SNVTr","SNVTrPerc","SNVTv","SNVTvPerc","Tr/Tv","MeanExACMax","MeanExACAll","DBSnp_E","DBSnp_EPerc","DBSNP_NE","DBSNP_NEPerc","in_hgmd","in_hgmdPerc",paste("cadd",cadd_thr,sep=">"),paste("cadd",cadd_thr,"Perc",sep="_"),"SplVar","SplVarPerc")

header_a = c("Analysis","TotVar","TotSNVs","SNVsPerc","TotINDELs","INDELsPerc","Mean_qual","SynSNV","SynSNVPerc","NonSynSNV","NonSynSNVPerc",
					"NS/S","SNVTr","SNVTrPerc","SNVTv","SNVTvPerc","Tr/Tv","MeanExACMax","MeanExACAll","DBSnp_E","DBSnp_EPerc","DBSNP_NE","DBSNP_NEPerc","in_hgmd","in_hgmdPerc",paste("cadd",cadd_thr,sep=">"),paste("cadd",cadd_thr,"Perc",sep="_"),"SplVar","SplVarPerc","stopgain","stopgain_perc","stoploss","stoploss_perc","frameshift","frameshift_perc")

header_s = c("Sample","Mean_GQ","TotReads")

tab = "\t"
cat(header_a,file=out_name_a,sep=tab,append = FALSE)
cat("\n",file=out_name_a,sep=tab,append = TRUE)
cat(header_s,file=out_name_s,sep=tab,append = FALSE)
cat("\n",file=out_name_s,sep=tab,append = TRUE)


#Open first result and create the compid column
vgout1 <- read.delim(vgout1_f,stringsAsFactors=F,sep="\t",head=T)	
vgout1$compid = paste(vgout1$chrom,vgout1$pos,vgout1$ref,vgout1$alt,sep="_")

#Open the second result
vgout2 <- read.delim(vgout2_f,stringsAsFactors=F,sep="\t",head=T)
#Open Tigem result and create the compid column
vgout2$compid = paste(vgout2$chrom,vgout2$pos,vgout2$ref,vgout2$alt,sep="_")


#If the analysis is a trio, then the samples are three, otherwise uses a single sample
if (anal_type == "TRIO"){
	#Get the samples name concatenating P, M and F for the first out
	samples_str1 = paste(paste(analysis,proband,sep="_") ,paste(analysis,"M",sep="_"),paste(analysis,"F",sep="_"),sep=",")
	#Get the samples name concatenating P, M and F for the second out
	samples_str2 = paste(paste(analysis2,proband,sep="_") ,paste(analysis2,"M",sep="_"),paste(analysis2,"F",sep="_"),sep=",")
}else{
	samples_str1 = analysis
	samples_str2 = analysis2
}


#Split string with sample names 
samples_l1 <- unlist(strsplit(samples_str1,",",fixed=T))
samples_l2 <- unlist(strsplit(samples_str2,",",fixed=T))

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


sample_stats <- function(vgout,sample_name) {
	refreads_f = paste(sample_name,"REF_reads",sep="_")
	altreads_f = paste(sample_name,"ALT_reads",sep="_")
	totreads_f = paste(sample_name,"tot_reads",sep="_")
	paltreads_f = paste(sample_name,"PercALT_reads",sep="_") 
	gt_f = paste(sample_name,"GT",sep="_")
	zyg_f = paste(sample_name,"ZYG",sep="_")	
	gq_f = paste(sample_name,"GQ",sep="_")			 
	
	mean_gq <- mean(na.omit(as.numeric(vgout[[gq_f]])))
	totreads <- mean(na.omit(as.numeric(vgout[[totreads_f]])))
				
	output <- list(mean_gq,totreads)
  
	return(output)
}




#In this analysis statistics I added also the stopgain, stoploss and frameshift number
#Nonsynonymous=missens
#Stop= nonsense	
analysis_stats <- function(vgout) {
	
	n_res <- nrow(vgout)#num results
	mean_q <- mean(na.omit(vgout$qual)) #mean quality
		 
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
	ns_syn_ratio <- nonsyn_n/syn_n
	
	stopgain <- grepl('exonic;stopgain',vgout$func_refgene)
	stopgain_n <- length(stopgain[stopgain==TRUE])
	stopgain_perc <- (stopgain_n/n_res)*100

	stoploss <- grepl('exonic;stoploss',vgout$func_refgene)
	stoploss_n <- length(stoploss[stoploss==TRUE])
	stoploss_perc <- (stoploss_n/n_res)*100

	frameshift <- grepl('exonic;frameshift',vgout$func_refgene)
	frameshift_n <- length(frameshift[frameshift==TRUE])
	frameshift_perc <- (frameshift_n/n_res)*100
	
	#Transitions and transversions among SNVs
	snp$snptype = mapply(snv_type,snp$ref,snp$alt)
	transit <- grepl('transition',snp$snptype)
	transit_n <- length(transit[transit==TRUE])
	transit_perc <- (transit_n/n_res)*100
	transv <- grepl('transversion',snp$snptype)
	transv_n <- length(transv[transv==TRUE])
	transv_perc <- (transv_n/n_res)*100
	tr_tv_ratio <- transit_n/transv_n
	
	#Frequencies and dbsnp
	#Here I remove both those <0 because VarGenius uses -1 when data is missing and also I remove the NA
	exac_max_mean <- mean(na.omit(vgout$exac_max[vgout$exac_max>0]))
	exac_all_mean <- mean(na.omit(vgout$exac_all[vgout$exac_all>0]))
	n_dbsnp <- nrow(vgout[vgout[[dbsnp_h]] != '-',])
	n_dbsnp_perc <- (n_dbsnp/n_res)*100
	n_notdbsnp <- nrow(vgout[vgout[[dbsnp_h]] == '-',])
	n_notdbsnp_perc <- (n_notdbsnp/n_res)*100
	#HGMD
	n_hgmd <- nrow(vgout[vgout[[hgmdgene_h]] != '-',])
	n_hgmd_perc <- (n_hgmd/n_res)*100
	#CADD PHRED SCORE > Threshold 
	cadd_del_n <- nrow(na.omit(vgout[as.numeric(vgout[[cadd_phred]])>cadd_thr,]))
	cadd_del_n_perc <-(cadd_del_n/n_res)*100
	#Splicing
	spl_var <-	sum(!is.na(as.numeric(vgout$refgene_splvar_dist_from_exon))) #splicing variants
	spl_var_perc <- (spl_var/n_res)*100			
	
  output <- list(n_res,snp_n,snp_perc,indel_n,indel_perc,mean_q,syn_n,syn_perc,nonsyn_n,nonsyn_perc,ns_syn_ratio,transit_n,transit_perc,
  transv_n,transv_perc,tr_tv_ratio,exac_max_mean,exac_all_mean,n_dbsnp,n_dbsnp_perc,n_notdbsnp,n_notdbsnp_perc,n_hgmd,n_hgmd_perc,cadd_del_n,
  cadd_del_n_perc,spl_var,spl_var_perc,stopgain_n,stopgain_perc,stoploss_n,stoploss_perc,frameshift_n,frameshift_perc)
  
  return(output)
}
########################################PRINT STATISTICS
print_statistics <- function(vgout1,vgout2,out_name_a,out_name_s,suffix) {

	#STATISTICS PER-ANALYSIS
	a_stats_1 <- analysis_stats(vgout1)

	cat(paste(analysis,suffix,sep="_"),file=out_name_a,append = TRUE)#Sample name
	for (i in a_stats_1) {
		cat(tab,file=out_name_a,append = TRUE)
		cat(i,file=out_name_a,append = TRUE)
	}
	cat("\n",file=out_name_a,append = TRUE)


	a_stats_2 <- analysis_stats(vgout2)
	
	cat(paste(analysis2,suffix,sep="_"),file=out_name_a,append = TRUE)#Sample name
	for (i in a_stats_2) {
		cat(tab,file=out_name_a,append = TRUE)
		cat(i,file=out_name_a,append = TRUE)
	}
	cat("\n",file=out_name_a,append = TRUE)

		
	#STATISTICS PER-SAMPLE	
	#For each sample..
	sample_num = 1
	for ( sample_name in samples_l1 ){
		#Get samples statistics for analysis1
		s_stats_1 <- sample_stats(vgout1,sample_name)
		
		cat(paste(sample_name,"1",suffix,sep="_"),file=out_name_s,append = TRUE)#Sample name
		for (i in s_stats_1) {
			cat(tab,file=out_name_s,append = TRUE)
			cat(i,file=out_name_s,append = TRUE)
	
		}
		cat("\n",file=out_name_s,append = TRUE)
		
		#For the second analysis use the sample name picked from the second list
		s_stats_2 <- sample_stats(vgout2,samples_l2[sample_num])
		cat(paste(samples_l2[sample_num],"2",suffix,sep="_"),file=out_name_s,append = TRUE)#Sample name
		for (i in s_stats_2) {
			cat(tab,file=out_name_s,append = TRUE)
			cat(i,file=out_name_s,append = TRUE)
		}
		cat("\n",file=out_name_s,append = TRUE)	
		
		sample_num = sample_num + 1
	}
}

#######################MAIN

###
#Here I filter those variant which have all HOMREF along all the samples
a1_zyg_p <- paste(analysis,"P","ZYG",sep="_")
a1_zyg_m <- paste(analysis,"M","ZYG",sep="_")
a1_zyg_f <- paste(analysis,"F","ZYG",sep="_")
vgout1_filt <- vgout1[ !(vgout1[[a1_zyg_p]] == 'HOMREF' & vgout1[[a1_zyg_m]] == 'HOMREF' & vgout1[[a1_zyg_f]] == 'HOMREF'  ),]
vgout1 <- vgout1_filt

a2_zyg_p <- paste(analysis2,"P","ZYG",sep="_")
a2_zyg_m <- paste(analysis2,"M","ZYG",sep="_")
a2_zyg_f <- paste(analysis2,"F","ZYG",sep="_")
vgout2_filt <- vgout2[ !(vgout2[[a2_zyg_p]] == 'HOMREF' & vgout2[[a2_zyg_m]] == 'HOMREF' & vgout2[[a2_zyg_f]] == 'HOMREF'  ),]
vgout2 <- vgout2_filt

	
	
#Analysis on the complete set of variants
print_statistics(vgout1,vgout2,out_name_a,out_name_s,"ALL")

#Analysis on unique variants
#Unique to vgout1
vgout1_uniq <- setdiff(vgout1$compid,vgout2$compid)
vgout1_uniq_d = merge (vgout1,as.data.frame(vgout1_uniq),by.x="compid",by.y="vgout1_uniq",all.x=F,all.y=T)
vgout1_uniq_p <- na.omit(vgout1_uniq_d)
#Unique to vgout2
vgout2_uniq <- setdiff(vgout2$compid,vgout1$compid)
vgout2_uniq_d = merge (vgout2,as.data.frame(vgout2_uniq),by.x="compid",by.y="vgout2_uniq",all.x=F,all.y=T)
vgout2_uniq_p <- na.omit(vgout2_uniq_d)
print_statistics(vgout1_uniq_p,vgout2_uniq_p,out_name_a,out_name_s,"UNIQ")

#print tables for uniq variants if wanted
write.table(vgout1_uniq_p,file=paste(analysis,"uniq",sep="."),sep="\t",row.names=F,quote=F )
write.table(vgout2_uniq_d,file=paste(analysis2,"uniq",sep="."),sep="\t",row.names=F,quote=F )

if (length(args) > 5){
	#Analysis on selected genes (matched with RefSeq)
	sel_g <- read.table(selected_genes_f,stringsAsFactors=F,sep="\t",head=T)
	vgout1_selg <- merge (vgout1,sel_g,by.x="gene_refgene",by.y="Gene",all.x=F,all.y=T)
	vgout1_selg_p <- na.omit(vgout1_selg)
	vgout2_selg <- merge (vgout2,sel_g,by.x="gene_refgene",by.y="Gene",all.x=F,all.y=T)
	vgout2_selg_p <- na.omit(vgout2_selg)
	print_statistics(vgout1_selg_p,vgout2_selg_p,out_name_a,out_name_s,"SELGMOTOR")

	#Analysis on unique variants with selected genes
	#Unique to vgout1
	vgout1_selguniq <- setdiff(vgout1_selg_p$compid,vgout2_selg_p$compid)
	vgout1_selguniq_d = merge (vgout1_selg_p,as.data.frame(vgout1_selguniq),by.x="compid",by.y="vgout1_selguniq",all.x=F,all.y=T)
	vgout1_selguniq_p <- na.omit(vgout1_selguniq_d)
	#Unique to vgout2
	vgout2_selguniq <- setdiff(vgout2_selg_p$compid,vgout1_selg_p$compid)
	vgout2_selguniq_d = merge (vgout2_selg_p,as.data.frame(vgout2_selguniq),by.x="compid",by.y="vgout2_selguniq",all.x=F,all.y=T)
	vgout2_selguniq_p <- na.omit(vgout2_selguniq_d)
	print_statistics(vgout1_selguniq_p,vgout2_selguniq_p,out_name_a,out_name_s,"SELGUNIQ")

}
	


#####################################LEGENDA

#TotVar: totale varianti
#TotSNVs: totale SNV
#SNVsPerc: percentuale di SN
#TotINDELs: totale INDEL
#INDELsPerc: percentuale di indel
#Mean_qual: Qualità media della chiamata
#SynSNV: SNV trovate sinonime
#SynSNVPerc: percentuale di SNV trovate sinonime
#NonSynSNV: SNV trovate non sinonime
#NonSynSNVPerc: percentuale di SNV trovate non sinonime
#NS/S: Non sinonime/sinonime Ratio
#SNVTr: Transitions
#SNVTrPerc: Transitions Percentage
#SNVTv: Transversions
#SNVTvPerc: Transversions PErcentage
#Tr/Tv: Transitions/Transversions Ratio
#MeanExACMax: Media tra tutti gli ExAC massimi
#MeanExACAll: Media tra tutti gli ExAC_all
#DBSnp_E: Presenti in DBSNP
#DBSnp_EPerc: Percentuale di presenti in DBSNP
#DBSNP_NE: Non presenti in DBSNP
#DBSNP_NEPerc: Percentuale di non presenti in DBSNP
#in_hgmd: Presenti in HGMD
#in_hgmdPerc: Percentuale di non presenti in HGMD
#cadd>20: Numero di varianti con CADD Phred Score > 20 (Probabile deleterious)
#cadd_20_Perc: Percentuale Numero di varianti con CADD Phred Score > 20
#SplVar: Putative Splicing
#SplVarPerc: Percentuale di Putative Splicing
