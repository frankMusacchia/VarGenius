##########Script to find common variants from multiple VarGenius output and filter them
#User Guide: 
#1. When you are searching for rare variants in TRIOS you can use functions to check for AR, denovo,
#HETCOMP and Xlinked
#2. When you are searching for "modifiers", hence variations that are common but apparently could not have
#a visible disruptive function you can use the functions to get common LoF or missense variants 
#3. If you want to just get specific type of variants in a subset of genes in a set of TRIO analyses, use this option	

#R CMD BATCH --no-save --no-restore '--args <output_from_vargenius> <samples_separate_by_comma>  <out_file_name> ' var_stats.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

length(args)
if (length(args) < 1){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <function>' multiple_samples_prioritization.R")
}

func = args[1]


out_folder = "." #Per lavorare in locale
suffix = "hg19.final_out_resgr14.txt"

compid_all <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
min_vc_qual = 50 #Minimum variant calling quality 
gene_selection = 0 
rare_variants = 0
#Thresholds for common variants
min_maf = 0.1
max_maf = 0.15
#Threshold for rare variants 
maf=0.02
min_var_in_same_gene = 2
gene_name_f = "gene_wgencodegencodebasicv19"
gene_func_f  = "func_wgencodegencodebasicv19"


##LOWE DISEASE May 2019 for De Matteis         
#out_folder = "/pico/work/TELET_TIGEM/FRANCESCO/PROVE/DEMATTEIS" 

#selgenes = "" #"genes_contr.txt" #"selected_genes_LOWE.txt"  #"trafficking.genes" #"selected_genes_DENT.txt"#
#analyses <- c('ADM_10','ADM_127','ADM_130','ADM_167','ADM_174','ADM_197','ADM_20','ADM_230','ADM_236','ADM_24','ADM_269','ADM_49','ADM_51','ADM_76','ADM_86','ADM_8','ADM_LW15','ADM_LW7','ADM_DN2','ADM_DN3','ADM_242')
#cases <- c('ADM_10','ADM_127','ADM_130','ADM_167','ADM_174','ADM_197','ADM_20','ADM_230','ADM_236','ADM_24','ADM_269','ADM_49','ADM_51','ADM_76','ADM_86','ADM_8','ADM_LW15','ADM_LW7')
#controls <- c('ADM_DN2','ADM_DN3','ADM_242')

#trio=0
#analysis = "ADM_DN2"



#Neural Tube defects for Marcello from Gaslini
#out_folder = "/pico/work/TELET_TIGEM/FRANCESCO/PROVE/GASLINI"
#analyses <- c('UD_GE001','UD_GE005', 'UD_GE025', 'UD_GE032','UD_GE046')
#trio=1
#analysis = 'UD_GE001'
#trio=0


####UDP Project - Groups
#analyses <- c('UD_BES017','UD_BES023','UD_GE004','UD_GE012','UD_GE015','UD_GE019','UD_GE023','UD_GE027','UD_GE039','UD_GE042','UD_MO009','UD_MO016','UD_MO019','UD_MO021','UD_MO030','UD_MO032','UD_MO033','UD_MO040','UD_MO049','UD_MO055','UD_MO076','UD_MO090','UD_ZOL001')

#Gruppo per NA105
analyses <- c('UD_GE059','UD_MO073','UD_MO092','UD_MO093','UD_MO119','UD_NA026','UD_NA068','UD_NA079','UD_NA087','UD_NA088','UD_NA105','UD_RM004','UD_RM1005','UD_RM1010','UD_RM1017','UD_RM1018','UD_RM1021','UD_RM1024')

#Given a singleton analysis
extract_var_counts <- function (analyses,sel_genes_f,trio,min_vc_qual,min_maf,max_maf,maf,out_name,vars_type,rare_variants,filter){
		
	#	ar_num <- 0
	#	denovo_num <- 0
	#	hetcomp_num <- 0
	#	xl_num <- 0
	#vars_type = "missense" #"lof" #
	savedcolnames = ""
	vgout_f_num <-0
	vgout_sel_num <- 0
	result_present <- 0
	
	#Change the output name accordingly with what you want to do
	if (rare_variants == 1){
		out_name = paste(out_name,"rare_vars",filter,sep="_")
	}else{
		out_name = paste(out_name,"common",vars_type,sep="_")
	}
	
	for (analysis in analyses){
		if (trio == 1){
			analysis_name = paste(analysis,"TRIO",sep="_")
			prob = paste(analysis,"P",sep="_")
			mother = paste(analysis,"M",sep="_")
			father = paste(analysis,"F",sep="_")
			gq = paste(prob,"GQ",sep="_")
			dp = paste(prob,"DP",sep="_")
			zyg = paste(prob,"ZYG",sep="_")
			f_zyg = paste(father,"ZYG",sep="_")
			m_zyg = paste(mother,"ZYG",sep="_")
		}else{
			analysis_name = analysis
			prob = analysis
			zyg = paste(prob,"ZYG",sep="_")
		}

		vgout_name <- paste(analysis_name,"finalout_out",paste(analysis_name,suffix,sep="."),sep="/")
		#vgout_name <- paste(analysis_name,suffix,sep=".") #IN LOCALE
		vgout_name
		vgout <- read.delim(vgout_name,stringsAsFactors=F,head=T,quote="",sep="\t")
		vgout$compid <- paste(vgout$chrom,vgout$pos,vgout$ref,vgout$alt,sep="_")
		#Filter relevant variants
		a<-vgout
		#There is an error in some output of VarGenius that will not permit to consider a lot of rows
		a[a$IVF=='-',]<--1


		 #If we have the trio, the variant must be present in the proband
		if (trio == 1){

		}else{ #Single sample


	  #Filter with a given list of genes
	  if ( sel_genes_f != '' ){
		  
		#File with genes to filter 
		genes<-read.delim(sel_genes_f,stringsAsFactors=F,head=T,quote="",sep="\t")
		  
		if (rare_variants == 1){
			
			 #Search RARE VARIANTS 			
			if (filter == 'strict'){
				#Strict Filter
				vgout_f <-a[ (a[[gene_name_f]] %in% genes$gene_name ) & 
				(a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
				grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
			}else{
				vgout_f <-a[ (a[[gene_name_f]] %in% genes$gene_name ) & 
				(a[[zyg]] == 'HOM' | a[[zyg]] == 'HET')  &
				!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
				a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
			}			 			
		}else{
			#Search COMMON VARIANTS WHICH COULD BE MODIFIERS
			if (vars_type == "lof"){
				#Lof in HET  EXCLUDED##| grepl('splicing',a[[gene_func_f]])
				vgout_f <-a[ a$qual > min_vc_qual & a$VC_ALGORITHM=='-' & (a[[zyg]] == 'HET')  & (a[[gene_name_f]] %in% genes$gene_name | a$gene_wgencodegencodebasicv19 %in% genes$gene_name)  &  (grepl('stop',a[[gene_func_f]]) | grepl(';frameshift',a[[gene_func_f]]) ) & as.numeric(a$gnomad_ex_max) < max_maf & as.numeric(a$db1000g_max) < max_maf & as.numeric(a$IVF) < max_maf, ]
			}else{
				#Missense
				vgout_f <-a[ a$qual > min_vc_qual & a$VC_ALGORITHM=='-' & (a[[zyg]] == 'HET')  & (a[[gene_name_f]] %in% genes$gene_name | a$gene_wgencodegencodebasicv19 %in% genes$gene_name)  &  grepl('exonic;nonsyn',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < max_maf & as.numeric(a$db1000g_max) < max_maf & as.numeric(a$IVF) < max_maf, ]
			}
		}
	  }
	  #NO gene list. Use all genes
	  else{
		 if (rare_variants == 1){
			 #Search RARE VARIANTS 			
			if (filter == 'strict'){
				#Strict Filter
				vgout_f <-a[ (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a$filter == 'PASS' & 
				!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
				grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
			}else{
				vgout_f <-a[ (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET')  &
				!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
				a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
			}
			#Search RARE VARIANTS 
			#Less stringent
			vgout_f <-a[ a$qual > min_vc_qual & a$VC_ALGORITHM=='-' & (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & !(grepl('ncRNA',a[[gene_func_f]])) & !(grepl('synonymous SNV',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]])) & as.numeric(a$exac_max)<maf & as.numeric(a$gnomad_ex_max)<maf & as.numeric(a$db1000g_max)<maf & as.numeric(a$IVF)<maf, ]
		}
		else{	
			#Search COMMON VARIANTS WHICH COULD BE MODIFIERS
			#Lof in HET EXCLUDED:| grepl('splicing',a[[gene_func_f]])
			if (vars_type == "lof"){
				vgout_f <-a[ a$qual > min_vc_qual & a$VC_ALGORITHM=='-' & (a[[zyg]] == 'HET') & (grepl('stop',a[[gene_func_f]]) | grepl(';frameshift',a[[gene_func_f]])  ) &  as.numeric(a$gnomad_ex_max) < max_maf & as.numeric(a$db1000g_max) < max_maf & as.numeric(a$IVF) < max_maf, ]
			}else{
				#Missense
				vgout_f <-a[ a$qual > min_vc_qual & a$VC_ALGORITHM=='-' & (a[[zyg]] == 'HET')  &  grepl('exonic;nonsyn',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < max_maf & as.numeric(a$db1000g_max) < max_maf & as.numeric(a$IVF) < max_maf, ]
			}
		}
	  }

	  #Generate a table with all the genes in the rows and the samples in columns
	  #containing the number of rare variants found in each sample for each gene
	  #The table is printed out of the for loop
	  #vgout_f_all <- data.frame(stringsAsFactors=FALSE)
	  if (nrow(vgout_f) > 0){
		  if (vgout_f_num == 0){
		   vgout_f_all<-as.data.frame(ftable(vgout_f$gene_refgene))
		   colnames(vgout_f_all) = c("gene",analysis_name)
		   vgout_f_num = vgout_f_num + 1
		   
		   if ( sel_genes_f != '' ){
			vgout_vars_all<-as.data.frame(vgout_f)
			vgout_vars_all$analysisname <- analysis
			savedcolnames <- colnames(vgout_vars_all)	
			colnames(vgout_vars_all)<-c(1:ncol(vgout_vars_all))	
		   }
		  }else{
		   curr<-as.data.frame(ftable(vgout_f$gene_refgene))
		   colnames(curr) = c("gene",analysis_name)
		   vgout_f_all<-merge(vgout_f_all,curr,by="gene",all=T)
		   vgout_f_num = vgout_f_num + 1
		  
		   if ( sel_genes_f != '' ){
			   curr<-as.data.frame(vgout_f)
			   curr$analysisname <- analysis
			   colnames(curr)<-c(1:ncol(curr))
			   vgout_vars_all<-rbind(vgout_vars_all,curr)			   
		   }
		   		   
		  }
	 }
	
	###END ELSE
 }
	}

	#If the  gene file is given, print the VarGenius output for all the variants in the selected genes
	if ( sel_genes_f != '' & length(savedcolnames) > 1 ){	
		colnames(vgout_vars_all) <- savedcolnames
		write.table(file=paste(out_folder,paste(out_name,"variants_in_genes.txt",sep="_"),sep="/"),vgout_vars_all,col.names=T,sep="\t",quote=F,row.names=F)
	}
	
	#Return the table with all the counts
	if (result_present == 1){
		vgout_f_all
	}
}

#Selects variated genes based on a ratio of presence of variants in cases and controls
get_top_variated_genes <- function (freq_table_all,analyses,outpath){
	num_analyses = ncol(freq_table_all)-1 #do not count the genename as an analysis
	#Add columns to count the NA
	freq_table_all$na_case <- apply(freq_table_all, 1, function(x) sum(is.na(x)))
	#Get ratio
	freq_table_all$ratio_case <- (num_analyses-freq_table_all$na_case)/length(analyses)
	#Mutated count
	freq_table_all$mutated <- num_analyses-freq_table_all$na_case
	#Get the genes where:
	#the ratio of mutated cases is gt the threshold
	#genes_case <- as.data.frame(case[case$ratio_case >= top_var_genes_thr ,c("gene")])
	genes_case <- as.data.frame(freq_table_all[freq_table_all$mutated >= min_var_in_same_gene ,c("gene")])
	
	results = nrow(genes_case )
	if ( results > 0){
		colnames(genes_case)=c("gene_name")
		write.table(file=outpath,genes_case,col.names=T,sep="\t",quote=F,row.names=F)
	}
	
	results	
}


#Selects variated genes based on a ratio of presence of variants in cases and controls
get_variants_on_sel_genes <- function (vgout_f_all){

	#Analysis of the table with variants count
	#Add columns to count the NA
	case <- vgout_f_all[,c("gene",cases)]
	contr <- vgout_f_all[,c("gene",controls)]
	case$na_case <- apply(case, 1, function(x) sum(is.na(x)))
	contr$na_contr <- apply(contr, 1, function(x) sum(is.na(x)))
	case$ratio_case <- (length(cases)-case$na_case)/length(cases)
	contr$ratio_contr <- (length(controls)-contr$na_contr)/length(controls)
	cas_contr <- cbind(case,contr)
	cas_contr$diff <- case$ratio_case - contr$ratio_contr
	cas_contr_s <- cas_contr[order(cas_contr$diff),]
	#Get the genes where:
	#the difference in ratio is gt half 0.5 
	# or in the first set maximum 1/3 are variated and only one is variated in the second set
	genes_case <- as.data.frame(cas_contr_s[cas_contr_s$diff >= 0.5 | (cas_contr_s$na_case < length(cases)/3  & cas_contr_s$na_contr == 1),c("gene")])
	genes_contr <- as.data.frame(cas_contr_s[cas_contr_s$diff <= -0.5 | (cas_contr_s$na_case < length(controls)/3  & cas_contr_s$na_contr == 1),c("gene")])
	colnames(genes_case)=c("gene_name")
	colnames(genes_contr)=c("gene_name")
	#Print tables for the list of genes relevant for cases and controls
	write.table(file=paste(out_folder,"genes_case.txt",sep="/"),genes_case,col.names=T,sep="\t",quote=F,row.names=F)
	write.table(file=paste(out_folder,"genes_contr.txt",sep="/"),genes_contr,col.names=T,sep="\t",quote=F,row.names=F)
	

}


select_suspect_genes <- function (analyses,count_table){
	numanalyses <- ncol(count_table) - 1 #Avoid the Gene column
	case <- count_table
	case$na_case <- apply(count_table, 1, function(x) sum(is.na(x)))
	case$ratio_case<-(numanalyses-case$na_case)/numanalyses
	#Get the genes where:
	#the ratio is at least > 0.1. Means that at least 10% of the samples have a variant in that gene
	genes_case <- as.data.frame(case[case$ratio_case >= 0.1 ,c("gene")])
	colnames(genes_case)=c("gene_name")

	write.table(file=paste(out_folder,"genes_case.txt",sep="/"),genes_case,col.names=T,sep="\t",quote=F,row.names=F)
	get_variants_types_count_on_multiple_trios(analyses,"genes_case.txt","denovo",out_folder)
	get_variants_types_count_on_multiple_trios(analyses,"genes_case.txt","ar",out_folder)
	get_variants_types_count_on_multiple_trios(analyses,"genes_case.txt","hetcomp",out_folder)
	get_variants_types_count_on_multiple_trios(analyses,"genes_case.txt","xl",out_folder)
}







#Given a list of analyses (TRIOS) searches for common AR, denovo, HETCOMP and Xlinked variants in genes.
#The protocol is the following:
#For each type of variation (denovo, hetcomp, ar, xl)
#1. Run it for all the genes (selgenes="") and a
#2. run get_top_variated_genes to obtain a selection of genes
#3. run again get_variants_types_count_on_multiple_trios with a list of selgenes for each type of variant

#STRICT FILTERING
#a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
#grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf

#OPEN FILTER
#!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
#a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
get_variants_types_count_on_multiple_trios <- function (analyses,sel_genes_f,vartype,filter,out_folder){

	#Depending by the filter set the MAF
	if (filter == 'strict'){
		#Strict Filter
		maf = 0.01
	}else{
		maf = 0.1	
	}
	
	ar_num <- 0
	denovo_num <- 0
	hetcomp_num <- 0
	xl_num <- 0
	vgout_f_num <-0
	vgout_sel_num <- 0
	
	#Initialize as empty data frame the complete vgout for the filtering
	vgout_all <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
 
 #analysis = 'UD_BES017','UD_BES023','UD_GE004','UD_GE012','UD_GE015'
                
	for (analysis in analyses){
		analysis_name = paste(analysis,"TRIO",sep="_")
		print (analysis)
		prob = paste(analysis,"P",sep="_")
		mother = paste(analysis,"M",sep="_")
		father = paste(analysis,"F",sep="_")
		gq = paste(prob,"GQ",sep="_")
		dp = paste(prob,"DP",sep="_")
		zyg = paste(prob,"ZYG",sep="_")
		f_zyg = paste(father,"ZYG",sep="_")
		m_zyg = paste(mother,"ZYG",sep="_")	  
		
		vgout_name <- paste(analysis_name,"finalout_out",paste(analysis_name,suffix,sep="."),sep="/")
		#vgout_name <- paste(analysis_name,suffix,sep=".") #IN LOCALE
		vgout_name
		if (file.info(vgout_name)$size > 0 & file.exists(vgout_name)){
			vgout <- read.delim(vgout_name,stringsAsFactors=F,head=T,quote="",sep="\t")
			vgout$compid <- paste(vgout$chrom,vgout$pos,vgout$ref,vgout$alt,sep="_")
			#Filter relevant variants
			a<-vgout
			#There is an error in some output of VarGenius that will not permit to consider a lot of rows
			a[is.na(a$IVF),] <- -1
			a[a$IVF=='-',] <- -1
			

			#Filter with a given list of genes
			if ( sel_genes_f != '' ){
			  
				#File with genes to filter 
				genes<-read.delim(sel_genes_f,stringsAsFactors=F,head=T,quote="",sep="\t")	
			
				if (vartype == 'denovo'){
					if (filter == 'strict'){
						#Strict Filter
						denovos <-a[ (a[[gene_name_f]] %in% genes$gene_name ) & 
						(a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HOMREF' &
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						denovos <-a[ (a[[gene_name_f]] %in% genes$gene_name ) & 
						(a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HOMREF' &
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
					vgout_f <- denovos[,c(1:vg_cols_num)]
				}
				
				#Find Autosomal Recessive
				if (vartype == 'ar'){
					if (filter == 'strict'){				
						#Strict Filter
						ar <-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a[[zyg]] == 'HOM' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HET' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{	
						#Permissive filter			
						ar <-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a[[zyg]] == 'HOM' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HET' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
					vgout_f <- ar[,c(1:vg_cols_num)]
				}
				
				#Find HETCOMP
				if (vartype == 'hetcomp'){
					if (filter == 'strict'){	
						#Strict Filter
						aHetP<-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a[[zyg]] == 'HET' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HOMREF' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
						
						aHetM<-a[(a[[gene_name_f]] %in% genes$gene_name ) &
						 a[[zyg]] == 'HET' &  a[[m_zyg]] == 'HET' & a[[f_zyg]] == 'HOMREF' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						#Permissive filter
						aHetP<-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a[[zyg]] == 'HET' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HOMREF' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
						
						aHetM<-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a[[zyg]] == 'HET' &  a[[m_zyg]] == 'HET' & a[[f_zyg]] == 'HOMREF' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
					
					genes_het <- intersect(aHetP$gene_refgene,aHetM$gene_refgene)
					aHetP_f<-aHetP[aHetP$gene_refgene %in% genes_het,]
					aHetM_f<-aHetM[aHetM$gene_refgene %in% genes_het,]
					hetcomp <- rbind(aHetP_f,aHetM_f)	
					vgout_f <- hetcomp[,c(1:vg_cols_num)]
				}

				#Find X-Linked Recessive - Mother is carrier and proband inherits. If the proband is male will manifest symptoms. If female will be carrier as well
				if (vartype == 'xl'){
					if (filter == 'strict'){				
						#Strict Filter
						xl <-a[ (a[[gene_name_f]] %in% genes$gene_name ) & 
						a$chrom == 'chrX' & (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HET' &
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						#Permissive filter
						xl <-a[ (a[[gene_name_f]] %in% genes$gene_name ) &
						a$chrom == 'chrX' & (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HET' &
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}		
					vgout_f <- xl[,c(1:vg_cols_num)]
				}
				
				#Append all results to a uniq table
				if (nrow(vgout_f) > 0){
					if (vgout_f_num == 0){
						vgout_all<-as.data.frame(vgout_f)
						vgout_all$analysisname <- analysis
						savedcolnames <- colnames(vgout_all)	
						vgout_all<-vgout_all[,c("analysisname", setdiff(names(vgout_all), "analysisname"))]
						colnames(vgout_all)<-c(1:ncol(vgout_all))
						
					}else{	
						curr<-as.data.frame(vgout_f)
					   curr$analysisname <- analysis
					   curr<-curr[,c("analysisname", setdiff(names(curr), "analysisname"))]
					   colnames(curr)<-c(1:ncol(curr))
					   vgout_all<-rbind(vgout_all,curr)
					}
					vgout_f_num = vgout_f_num + 1
				}
				
			}
			#######
			#No gene list used. The operation is performed for all genes
			else{
				#Find denovos
				if (vartype == 'denovo'){
					if (filter == 'strict'){
						#Strict Filter
						vgout_f_filt <-a[ (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HOMREF' &
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						vgout_f_filt <-a[ (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HOMREF' &
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
				}
				
				#Find Autosomal Recessive
				if (vartype == 'ar'){
					if (filter == 'strict'){				
						#Strict Filter
						vgout_f_filt <-a[ a[[zyg]] == 'HOM' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HET' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{	
						#Permissive filter			
						vgout_f_filt <-a[ a[[zyg]] == 'HOM' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HET' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
				}
				
				#Find HETCOMP
				if (vartype == 'hetcomp'){
					if (filter == 'strict'){	
						#Strict Filter
						aHetP<-a[ a[[zyg]] == 'HET' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HOMREF' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
						
						aHetM<-a[ a[[zyg]] == 'HET' &  a[[m_zyg]] == 'HET' & a[[f_zyg]] == 'HOMREF' & 
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						#Permissive filter
						aHetP<-a[ a[[zyg]] == 'HET' &  a[[f_zyg]] == 'HET' & a[[m_zyg]] == 'HOMREF' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
						
						aHetM<-a[ a[[zyg]] == 'HET' &  a[[m_zyg]] == 'HET' & a[[f_zyg]] == 'HOMREF' & 
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}
					
					genes_het <- intersect(aHetP$gene_refgene,aHetM$gene_refgene)
					aHetP_f<-aHetP[aHetP$gene_refgene %in% genes_het,]
					aHetM_f<-aHetM[aHetM$gene_refgene %in% genes_het,]
					vgout_f_filt <- rbind(aHetP_f,aHetM_f)
				}

				#Find X-Linked Recessive - Mother is carrier and proband inherits. If the proband is male will manifest symptoms. If female will be carrier as well
				if (vartype == 'xl'){
					if (filter == 'strict'){					
						#Strict Filter
						maf = 0.01				
						vgout_f_filt <-a[ a$chrom == 'chrX' & (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HET' &	
						a$filter == 'PASS' & !(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) &
						grepl('exonic',a[[gene_func_f]]) & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf,]
					}else{
						#Permissive filter
						vgout_f_filt <-a[ a$chrom == 'chrX' & (a[[zyg]] == 'HOM' | a[[zyg]] == 'HET') & a[[f_zyg]] == 'HOMREF' & a[[m_zyg]] == 'HET' &	
						!(grepl('synonymous SNV',a[[gene_func_f]])) & !(grepl('ncRNA',a[[gene_func_f]])) & (grepl('exonic',a[[gene_func_f]]) |
						a[[gene_func_f]] == 'splicing') & as.numeric(a$gnomad_ex_max) < maf & as.numeric(a$db1000g_max) < maf & as.numeric(a$IVF)< maf ,]
					}				
				}
				
				#Define the complete file if it is empty
				if (nrow(vgout_f_filt) > 0){
					
					
						#Initialize as empty data frame the complete vgout for the filtering
		a <- data.frame(Date=as.Date(character()),
					 File=character(), 
					 User=character(), 
					 stringsAsFactors=FALSE) 
			  colnames(a) = c(analysis_name,"gene_wgencodegencodebasicv19")       
					 
					if (vgout_f_num == 0){
						freq_table_all <-as.data.frame(ftable(vgout_f_filt$gene_wgencodegencodebasicv19))
						colnames(freq_table_all) = c("gene",analysis_name)
						vgout_f_num = vgout_f_num+1					
					}else{
						#Merge with the all file
						curr<-as.data.frame(ftable(vgout_f_filt$gene_wgencodegencodebasicv19))
						colnames(curr) = c("gene",analysis_name)
						freq_table_all<-merge(freq_table_all,curr,by="gene",all=T)
						vgout_f_num = vgout_f_num+1				
					}
				}	
			}			
		}
	}
	
	#If the  gene file is given, print the VarGenius output for all the variants in the selected genes
	if ( sel_genes_f != ''  ){
		if ( nrow(vgout_all) > 0){
			out_name = analyses[1]
			colnames(vgout_all) <- savedcolnames
			write.table(file=paste(out_folder,paste("group",out_name,vartype,filter,"variants_in_genes.txt",sep="_"),sep="/"),vgout_all,col.names=T,sep="\t",quote=F,row.names=F)
			#Return the filtered output
			vgout_all
		}
	}else{
		freq_table_all
	}
}


###############################################MAIN


if (func == 1){
	####################################

	#Given a set of Trios extract first a table with the number of variants of each "vartype" (denovo,ar,etc..)
	#for all genes
	#Then the table of occurrencies is filtered to obtain a list of most variated genes
	#Finally the most variated genes list is used to obtain the VarGenius annotation (columns from 1: vg_cols_num)
	#
	filepath = args[2]
	filepath
	f = file(filepath, "r")
	while ( TRUE ) {
		line = readLines(f, n = 1)
		if ( length(line) == 0 ) {	
			break
		}
		analyses = unlist(strsplit(line, ","))
		
		vg_cols_num <- 70 #Columns of VarGenius output to take into the result
		filters <- c("strict","none")
		vartypes <- c("denovo","ar","hetcomp","xl")
		
#		filter="strict"
#		vartype="denovo"
		for (filter in filters) {
			for (vartype in vartypes) {
				#Get the table with count of variations per-gene for all genes for a type of variants
				sel_genes_f = ""
				freq_table_all <- get_variants_types_count_on_multiple_trios(analyses,sel_genes_f,vartype,filter,out_folder)
				
				#Get the TOP-Variated gene list
				sel_genes_f = paste(out_folder,paste(vartype,"genes_case.txt",sep="_"),sep="/")
				num_results <- get_top_variated_genes(freq_table_all,analyses,sel_genes_f)
				if (num_results > 0 ){
					#Get the VarGenius ouput for this type of variants and for only the selected genes
					get_variants_types_count_on_multiple_trios(analyses,sel_genes_f,vartype,filter,out_folder)
				}
			}
		}
	}

  close(f)
	########################################
}

#With this method I search within a set of cases and a set of controls the genes which could have
#rare variants. The variants are compared and genes are selected.
#After gene selection the variants are again fetched and the complete VarGenius output is produced for all samples
#and all suspect variants
if (func == 2){
	rare_variants = 1
	vartype = "" #"lof" #
	trio=0
	filter ="strict"
	###################
	########Anaysis of single samples with cases and controls
	#First extracts the table with variants count for each gene
	vgout_f_all <- extract_var_counts(analyses,"",trio,min_vc_qual,min_maf,max_maf,maf,out_name,vartype,rare_variants)

	#Compares cases and controls and prints the VarGenius output for a set of relevant genes in both cases
	get_variants_on_sel_genes(vgout_f_all)
		
	#Use the gene sets to get the variants associated
	#extract_var_counts(cases,trio,min_vc_qual,min_maf,max_maf,maf,"cases")
	#extract_var_counts(controls,trio,min_vc_qual,min_maf,max_maf,maf,"controls")
	
	#Now that you have relevant genes, reexecute the extraction of variants to get only variants on the selected genes
	#and print the VarGenius output
	genes_f = paste(out_folder,"genes_case.txt",sep="/")
	extract_var_counts(cases,genes_f,trio,min_vc_qual,min_maf,max_maf,maf,"cases",vartype,rare_variants,filter)
	genes_f = paste(out_folder,"genes_contr.txt",sep="/")
	extract_var_counts(controls,"genes_contr.txt",trio,min_vc_qual,min_maf,max_maf,maf,"controls",vartype,rare_variants,filter)
	

	if (trio == 1){

	}else{
		#Print the table with the occurrences of variants per gene
		if (nrow(vgout_f_all) > 0){
			
			#Print it
			write.table(file=paste(out_folder,paste("variants_in_all_genes.txt",sep="_"),sep="/"),vgout_f_all,col.names=T,sep=" ",quote=F,row.names=F)
		}
	}
	####################################
}

### If you want to just get specific type of variants in a subset of genes in a set of TRIO analyses, use this option	
# Input: - a file with the list of analyses to use. They must be all TRIOs and all the analysis names must be XXX_TRIO, YYY_TRIO, etc.
#				while the list must contain only XXX,YYY
# 		 - a file with a list of genes. Column name must be "gene_name"	 
if (func == 3){
	
	filepath = args[2]
	print(filepath)
	sel_genes_f = args[3]
	print(sel_genes_f)
	out_folder = args[4]
	
	filters <- c("strict","none")
	vartypes <- c("denovo","ar","hetcomp","xl")#all
	#vartypes <- c("ar","hetcomp")#only biallelic
	
	analyses <- read.delim(filepath,stringsAsFactors=F,head=F,quote="",sep="\t")
	vg_cols_num <- 70 #Columns of VarGenius output to take into the result

	for (filter in filters) {
		for (vartype in vartypes) {
			#Get the VarGenius ouput for this type of variants and for only the selected genes
			get_variants_types_count_on_multiple_trios(analyses$V1,sel_genes_f,vartype,filter,out_folder)
		}
	}
}

### SEARCHING FOR MODIFIERS IN SUSPECT GENES in single samples

if (func == 4){
	#Thresholds for common variants
	min_maf = 0.1
	max_maf = 0.16
	trio=0
	rare_variants = 0
	filter=""
#	cases_f="/pico/work/TELET_TIGEM/prove_vargenius/ADN_DENT_samples.txt"
#	controls_f="/pico/work/TELET_TIGEM/prove_vargenius/ADN_Lowe_samples.txt"
#	sel_genes_f="/pico/work/TELET_TIGEM/vargenius_analyses/DENT-LOWE_genes.txt"
	cases_f = args[2]#The names of samples
	print(cases_f)
	controls_f = args[3]#The names of samples
	print(controls_f)
	sel_genes_f = args[4]
	print(sel_genes_f)
	
	cases <- read.delim(cases_f,stringsAsFactors=F,head=F,quote="",sep="\t")
	controls <- read.delim(controls_f,stringsAsFactors=F,head=F,quote="",sep="\t")
				
	#vg_cols_num <- 70 #Columns of VarGenius output to take into the result
	vartypes <- c("lof","missense")#only biallelic

	for (vartype in vartypes) {
		extract_var_counts(cases$V1,sel_genes_f,trio,min_vc_qual,min_maf,max_maf,maf,"cases",vartype,rare_variants,filter)
		extract_var_counts(controls$V1,sel_genes_f,trio,min_vc_qual,min_maf,max_maf,maf,"controls",vartype,rare_variants,filter)
	}

}


if (func == 5){
	#Thresholds for common variants
	maf=0.02
	trio=0
	rare_variants = 1
	filter="strict"
#	cases_f="/pico/work/TELET_TIGEM/prove_vargenius/ADN_DENT_samples.txt"
#	controls_f="/pico/work/TELET_TIGEM/prove_vargenius/ADN_Lowe_samples.txt"
#	sel_genes_f="/pico/work/TELET_TIGEM/vargenius_analyses/DENT-LOWE_genes.txt"
	cases_f = args[2]#The names of samples
	print(cases_f)
	controls_f = args[3]#The names of samples
	print(controls_f)
	sel_genes_f = args[4]
	print(sel_genes_f)
	
	cases <- read.delim(cases_f,stringsAsFactors=F,head=F,quote="",sep="\t")
	controls <- read.delim(controls_f,stringsAsFactors=F,head=F,quote="",sep="\t")
				
	#vg_cols_num <- 70 #Columns of VarGenius output to take into the result
	#vartypes <- c("lof","missense")#only biallelic

	#for (vartype in vartypes) {
		extract_var_counts(cases$V1,sel_genes_f,trio,min_vc_qual,min_maf,max_maf,maf,"cases",vartype,rare_variants,filter)
		extract_var_counts(controls$V1,sel_genes_f,trio,min_vc_qual,min_maf,max_maf,maf,"controls",vartype,rare_variants,filter)
	#}

}
