#R CMD BATCH --no-save --no-restore '--args <table1> <table2> <outname>' r_utils.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

job_type = args[1] 

job_type

out_format = "png"

#Dimensions and quality of images
pdf_width = 20
pdf_height = 20
pdf_unit = "cm"
image_width = 10
image_height = 10
image_dpi = 150 #I used 75 because of problems with R


#JOB SELECTION
if (job_type == "HQVARPLOT"){
	library(ggplot2)
	hq_gq = 50 #High quality GQ value
	hq_varfilter = "PASS" #High quality variants filter
	
	#This function takes in input the output table from VarGenius and builds a table with high quality variants
	#for each sample of the analysis. Then, it filters from that 6 tables and for each of that 
	#generates a plot
	#  - table with only HET SNPs
	#  - table with only HOM SNPs
	#  - table with only HOMREF SNPs
	#  - table with only HET INDELs
	#  - table with only HOM INDELs
	#  - table with only HOMREF INDELs
	
	length(args)
	if (length(args) < 4){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <VarGenius output> <samples_list> <out_prefix>' r_utils.R")
	}
	vg_out_f = args[2] #VarGenius output
	vg_out_f
	samples_str = args[3] #Samples list
	samples_str
	out_prefix = args[4] #Output prefix for the images
	out_prefix
	
	#Open the table
	vg_out <- read.table(vg_out_f,stringsAsFactors=F,head=T,sep="\t")
	
	#Split string with sample names 
	samples_l <- unlist(strsplit(samples_str,",",fixed=T))

	#For each sample..
	for ( sample in samples_l ){
		#Filter the variants given in the table using GQ and FILTER values
		hq_vars <- vg_out[as.numeric(vg_out[[paste(sample,"GQ",sep="_")]]) > hq_gq & vg_out$filter == hq_varfilter,]

		#Put a new column which gives the proportion of ALT reads respect with the total reads number
		hq_vars$ALTReadsProp <- hq_vars[[paste(sample,"ALT_reads",sep="_")]] / (hq_vars[[paste(sample,"ALT_reads",sep="_")]] + hq_vars[[paste(sample,"REF_reads",sep="_")]])
		
	  #Filter to select only the needed columns
	  refreads_f = paste(sample,"REF_reads",sep="_")
	  altreads_f = paste(sample,"ALT_reads",sep="_")
	  totreads_f = paste(sample,"tot_reads",sep="_")
	  paltreads_f = paste(sample,"PercALT_reads",sep="_") 
	  gt_f = paste(sample,"GT",sep="_")
	  zyg_f = paste(sample,"ZYG",sep="_")
	  gq_f = paste(sample,"GQ",sep="_")
	  altreadsp_f = paste(sample,"GQ",sep="_") 
	  hq_vars_filt <- hq_vars[,c('chrom','pos','ref','alt','qual','filter',refreads_f,altreads_f,totreads_f,paltreads_f,gt_f,zyg_f,gq_f,'ALTReadsProp')]

		
		#Now separate the table in two tables containing one the SNPs and another the INDELs
		hq_snps <- hq_vars_filt[nchar(hq_vars_filt$ref) == 1 & nchar(hq_vars_filt$alt) == 1,]
		hq_indels <- hq_vars_filt[nchar(hq_vars_filt$ref) > 1 | nchar(hq_vars_filt$alt) > 1,]
		
		#######################
		#Get plots for SNPs
		######################
		
		#Heterozygous SNPs	
		snps_het <- hq_snps[hq_snps[[zyg_f]] == 'HET',]
		snps_het[[paltreads_f]] <- as.numeric(snps_het[[paltreads_f]])
		snps_het[[totreads_f]] <- as.numeric(snps_het[[totreads_f]])
		#Omit NA rows
		snps_het_nona <- na.omit(snps_het)
		p_het <- ggplot(snps_het_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="blue")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HET SNPs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HET_SNPs_dispersion",sep="_"),out_format,sep="."), width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HET_SNPs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}
		
		#Homozygous SNPs	
		snps_hom <- hq_snps[hq_snps[[zyg_f]] == 'HOM',]
		snps_hom[[paltreads_f]] <- as.numeric(snps_hom[[paltreads_f]])
		snps_hom[[totreads_f]] <- as.numeric(snps_hom[[totreads_f]])
		#Omit NA rows
		snps_hom_nona <- na.omit(snps_hom)
		p_hom <- ggplot(snps_hom_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="green")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HOM SNPs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#ggsave(paste(paste(out_prefix,sample,"HOM_SNPs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HOM_SNPs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HOM_SNPs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}		
				
		#Homozygous Reference SNPs	
		snps_homref <- hq_snps[hq_snps[[zyg_f]] == 'HOMREF',]
		snps_homref[[paltreads_f]] <- as.numeric(snps_homref[[paltreads_f]])
		snps_homref[[totreads_f]] <- as.numeric(snps_homref[[totreads_f]])
		#Omit NA rows
		snps_homref_nona <- na.omit(snps_homref)
		p_homref <- ggplot(snps_homref_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="red")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HOMREF SNPs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#ggsave(paste(paste(out_prefix,sample,"HOMREF_SNPs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HOMREF_SNPs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HOMREF_SNPs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}	
				
		######################
		#Get plots for INDELs
		######################
		
		#Heterozygous INDELs	
		indels_het <- hq_indels[hq_indels[[zyg_f]] == 'HET',]
		indels_het[[paltreads_f]] <- as.numeric(indels_het[[paltreads_f]])
		indels_het[[totreads_f]] <- as.numeric(indels_het[[totreads_f]])
		#Omit NA rows
		indels_het_nona <- na.omit(indels_het)
		p_het <- ggplot(indels_het_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="blue")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HET INDELs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#ggsave(paste(paste(out_prefix,sample,"HET_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HET_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HET_INDELs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}	

		#Homozygous indels	
		indels_hom <- hq_indels[hq_indels[[zyg_f]] == 'HOM',]
		indels_hom[[paltreads_f]] <- as.numeric(indels_hom[[paltreads_f]])
		indels_hom[[totreads_f]] <- as.numeric(indels_hom[[totreads_f]])
		#Omit NA rows
		indels_hom_nona <- na.omit(indels_hom)
		p_hom <- ggplot(indels_hom_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="green")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HOM INDELs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#ggsave(paste(paste(out_prefix,sample,"HOM_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HOM_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HOM_INDELs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}
						
		#Homozygous Reference indels	
		indels_homref <- hq_indels[hq_indels[[zyg_f]] == 'HOMREF',]
		indels_homref[[paltreads_f]] <- as.numeric(indels_homref[[paltreads_f]])
		indels_homref[[totreads_f]] <- as.numeric(indels_homref[[totreads_f]])
		#Omit NA rows
		indels_homref_nona <- na.omit(indels_homref)
		p_homref <- ggplot(indels_homref_nona,aes_string(y=totreads_f,x=paltreads_f)) + geom_point(size=0.1, colour="red")+
		ggtitle("TotReads compared with \nPercentage of ALT alleles in HOMREF INDELs") + theme(axis.text.x=element_text(angle=60, hjust=1))
		#ggsave(paste(paste(out_prefix,sample,"HOMREF_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)		
		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(out_prefix,sample,"HOMREF_INDELs_dispersion",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(out_prefix,sample,"HOMREF_INDELs_dispersion",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}
			
	}
}

##################################################DISEASECOVBOXPLOT

#Input: needs in input
#				- the >10X coverage data for all genes for a single analysis
#				- a set of files with gene lists for each disease wanted. The file name will be used in the plot
if (identical(job_type,"DISEASECOVBOXPLOT")){
	length(args)
	if (length(args) < 5){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <VarGenius output> <samples_list> <out_prefix>' r_utils.R")
	}
	cov_data_f = args[2] 
	cov_data_f
	disease_folder = args[3] #Diseases files
	disease_folder
	disease_f_str = args[4] #Diseases files
	disease_f_str
	outsuffix = args[5]
	outsuffix
	
	library(ggplot2)
	library(reshape2)#Needed for 'melt'

	#Split string with disease file names
	disease_l <- unlist(strsplit(disease_f_str,",",fixed=T))

	cov_t <- read.table(cov_data_f,stringsAsFactors=F,head=T,sep="\t")
	
	#Rearrange the table putting all the values in a single column
	cov_t_rearr <- melt(cov_t)
	
	#
	complete_table <- data.frame(Gene=character(),Sample=character(), Coverage=double(),Panel =character(), stringsAsFactors=FALSE) 
	#For each disease
	for ( disease in disease_l ){
		#Select those genes belonging to the given disease
		sel_genes <- read.table(paste(disease_folder,disease,sep="/"),stringsAsFactors=F,sep="\t",quote="")
		#Filter the coverage table picking only the selected genes
		cov_t_sel_genes <- merge(cov_t_rearr,sel_genes,by.x="Gene",by.y="V1",all.y=T)	
		#Add a column with the disease description
		cov_t_sel_genes$Panel <- disease

		colnames(cov_t_sel_genes) = c("Gene","Sample","Coverage","Panel")
		#Convert values to log10
		cov_t_sel_genes$Coverage <- log10(as.numeric(cov_t_sel_genes$Coverage))
		#Bind this table to rows of the previous
		complete_table <- rbind(complete_table,cov_t_sel_genes)

	}
		#Do a boxlot
		ggplot(data = complete_table, aes(x=Panel, y=Coverage)) + geom_boxplot(aes(fill=Sample))+
			facet_wrap( ~ Panel, scales="free") + ggtitle("Coverage per panel")

		#Different mode of save depending by the format
		if ( out_format == 'pdf'){ 
			ggsave(paste(paste(outsuffix,"coverage_per_panel",sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit)
		}else{
			ggsave(paste(paste(outsuffix,"coverage_per_panel",sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi)
		}
}

########################################COVLT99
#Obtains a tables for all genes and all disease with coverage less than 99%
# e.g.: Gene Coverage DiseasePanel
if (job_type == "COVLT99"){
	length(args)
	if (length(args) < 6){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <VarGenius output> <samples_list> <out_prefix>' r_utils.R")
	}
	cov_data_f = args[2] 
	cov_data_f
	samples_str = args[3] #Samples list
	samples_str
	disease_folder = args[4] #Diseases files
	disease_folder
	disease_f_str = args[5] #disease list
	disease_f_str
	outsuffix = args[6]
	outsuffix
		
	cov_p_field = "above_10"

	cov_data <- read.table(cov_data_f,stringsAsFactors=F,head=T,sep="\t",quote="")
	
	#Split string with sample names 
	samples_l <- unlist(strsplit(samples_str,",",fixed=T))

	#Split string with disease file names
	disease_l <- unlist(strsplit(disease_f_str,",",fixed=T))
			
	#For each sample..
	for ( sample in samples_l ){
		#sample<-samples_l[1]#
		#sample
		field = paste(sample,cov_p_field,sep="_._")#"above_10"
		field
		cov_data[field]
		lt99 <- cov_data[cov_data[[field]] < 99,c("Gene",field)]

		complete_table <- data.frame(Gene=character(),Coverage=double(),Panel =character(), stringsAsFactors=FALSE) 

		#For each disease
		for ( disease in disease_l ){
			#disease = disease_l[1]#DEBUGCODE
			dis_genes <- read.table(paste(disease_folder,disease,sep="/"),stringsAsFactors=F,sep="\t",quote="")
			lt99_dis <- merge(lt99,dis_genes,by.x="Gene",by.y="V1",all.y=T)
			names <- unlist(strsplit(disease,"/"))
			lt99_dis$Panel <- names[length(names)]
			lt99_dis_nona <- na.omit(lt99_dis)
			complete_table <- rbind(complete_table,lt99_dis_nona)
		}
		write.table(complete_table,file=paste(outsuffix,sample,sep="_"),quote=F,row.names=F,sep="\t")
	}
	
}

########################################TABLECOVABOVE
##INPUT: 
			# - The table from DepthOfCoverage "sample_gene_summary" giving the coverage at different level for all samples in the analysis
			# - list of selected genes to pick for the table (use ALL_GENES to pick all genes)
			# - list of samples of the analysis
			# - the coverage percentage field to use (above_10, above_100, etc)
##OUTPUT
			# Returns a table 
if (job_type == "TABLECOVABOVE"){

	length(args)
	if (length(args) < 6){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <in_table> <sel_genes> <samples_str> <cov_p_field> <out_name>' r_utils.R")
	}
	#Output file name
	in_table <- args[2] #"chr.buttami"#
	in_table
	#Output from VarGenius
	sel_genes <- args[3] 
	sel_genes 
	#List of samples separated by comma
	samples_str <- args[4] #c("UD15","UD16","UD17") #
	samples_str

	#The coverage percentage field to use
	cov_p_field <- args[5]
	cov_p_field

	#Output name
	out_name <- args[6]
	out_name

	#open the sample_gene_summary table from DOC
	gc<-read.table(in_table,stringsAsFactors=F,head=T,sep="\t",quote="")
	
	#Filter the table with the list of genes, if needed
	if (sel_genes != "ALL_GENES"){
		sel <- read.table(sel_genes,stringsAsFactors=F,sep="\t",quote="")
		gc_sel <- merge(gc,sel,by.x="Gene",by.y="V1",all.y=T)
	}else{
		gc_sel <- gc
	}

	#Split string with sample names 
	samples_l <- unlist(strsplit(samples_str,",",fixed=T))

	#Get the list of genes into the result
	numcol = 1
	df <- gc_sel[1]
	
	#bind the column with the selected coverage associated to the any sample into the final table
	for ( sample in samples_l ){
		field = paste(sample,cov_p_field,sep="_._")#"above_10"
		df <- cbind(df,gc_sel[field])
		numcol = numcol + 1
	}
  colnames(df) = append(c("Gene"),samples_l )
	#Write output
	write.table(df,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)
}


#################################################JOB ALLSAMPLESCOVPLOT
#INPUT:
# - a table containing for all samples considered selected columns from sample_cumulative_coverage_proportions
		#.sample_cumulative_coverage_proportions header is the following:
		#			SampleName gte_0 gte_1 .. gte_500
# - a flag indicating if the information about the kinship is present (used to give colors)
# - the list of selected coverage levels (e.g. 1X, 10x, 20x -> gte_1,gte_10,gte_20
# - the title of the plot
# - the name of the output file
##OUTPUT:
# Prints using ggplot as many plots as the selected coverage levels. Hence we have global plot sdescribing
# what is for each sample the coverage level
if (job_type == "ALLSAMPLESCOVPLOT"){

	library(ggplot2)
	length(args)
	if (length(args) < 5){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <task> <in_table> <kinship> <covlev_str> <out_name>' r_utils.R")
	}
	#Output file name
	in_table <- args[2] #"chr.buttami"#
	in_table
	#Flag to indicate presence of kinship information
	kinship <- args[3] 
	kinship 
	#List of coverage levels separated by comma
	covlev_str <- args[4] 
	covlev_str

	#Plot title
	title <-  args[5] 
	title
	
	#Output name
	out_name <- args[6]
	out_name

	out_format = 'png'
	#Split string with sample names 
	covlev_l <- unlist(strsplit(covlev_str,",",fixed=T))

	cov_t <- read.table(in_table,stringsAsFactors=F,head=T,sep="\t",quote="")
	#Print a plot for each coverage level
	for ( covlev in covlev_l ){
	  cov_t$logcl <- log10(cov_t[[covlev]])
		if ( kinship == 1){
			#Plot the coverage levels for all the families using the color for the kinship
			ggplot(cov_t,aes_string(x="sample_name",y="logcl"))+
			geom_point(aes(color=kinship))+theme(axis.text.x=element_text(angle=60, hjust=1))	+
			scale_y_continuous(labels=function(x){format(round(10^x, 2), nsmall = 2)}) +#Changing the labels restoring the 0-1 value
			labs(x="Samples",y="Percentage of locus") + 
			ggtitle(paste(title,covlev,"X",sep=" "))
			
		}else{
			#Plot the coverage levels for all the families
			ggplot(cov_t,aes_string(x="sample_name",y="logcl"))+geom_point()+theme(axis.text.x=element_text(angle=60, hjust=1))+
			scale_y_continuous(labels=function(x){format(round(10^x, 2), nsmall = 2)}) + #Changing the labels restoring the 0-1 value
			labs(x="Samples",y="Percentage of locus") + 
			ggtitle(paste(title,covlev,"X",sep=" "))		
		}

		#Print output
		#The plot width will be extended for each additional sample considered. Will add 0.1 for PNG format and 0.2 for PDF
		pdf_width = nrow(cov_t) * 0.5
		image_width = nrow(cov_t) * 0.2
		if ( out_format == 'pdf'){
			ggsave(paste(paste(out_name,"ALLSAMPLESCOVPLOT",covlev,sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit, limitsize = FALSE)
		}else{
			ggsave(paste(paste(out_name,"ALLSAMPLESCOVPLOT",covlev,sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi, limitsize = FALSE)
		}		
	}

}



############################################JOB ALLSAMPLESREADSPLOT
##INPUT:
# 	- a table containing for all samples considered the following values: SampleName TotalReads TotalRemoved ProperlyPaired
# 	-	a flag saying if the kinship values is present
# 	- list of columns names to be plotted into different plots
##OUTPUT:
# Prints using ggplot as many dot plots as the columns names. Hence, we can have a global plot describing
# the ovarall amount of reads, the overall removed reads and the properly paired
if (job_type == "ALLSAMPLESREADSPLOT"){

	library(ggplot2)
	length(args)
	if (length(args) < 5){
		stop("USAGE: R CMD BATCH --no-save --no-restore '--args <task> <in_table> <kinship> <colnames_str> <out_name>' r_utils.R")
	}
	#Output file name
	in_table_f <- args[2] #"chr.buttami"#
	in_table_f
	#Flag to indicate presence of kinship information
	kinship <- args[3] 
	kinship 
	#List of columns names to be plotted
	colnames_str <- args[4] 
	colnames_str

	#Output name
	out_name <- args[5]
	out_name

	out_format = 'png'
	#Split string with sample names 
	colnames_l <- unlist(strsplit(colnames_str,",",fixed=T))

	in_table <- read.table(in_table_f,stringsAsFactors=F,head=T,sep="\t",quote="")
	#Print a plot for each column name
	for ( colname in colnames_l ){
	  #Get the logarithm here
	  in_table$logcl <- log10(in_table[[colname]])
		if ( kinship == 1){
			#Plot the data for all the rows
			#ggplot(in_table,aes_string(x="sample_name",y=colname))+
			ggplot(in_table,aes_string(x="sample_name",y="logcl"))+
			geom_point(aes(color=kinship))+theme(axis.text.x=element_text(angle=60, hjust=1))	+
			ggtitle(paste("Global ",colname,sep=" "))
		}else{
			#Plot the coverage levels for all the families
			#ggplot(in_table,aes_string(x="sample_name",y=colname))+geom_point()+theme(axis.text.x=element_text(angle=60, hjust=1))+
			ggplot(in_table,aes_string(x="sample_name",y="logcl"))+geom_point()+theme(axis.text.x=element_text(angle=60, hjust=1))+
			ggtitle(paste("Global ",colname,sep=" "))		
		}

		#Print output
		#The plot width will be extended for each additional sample considered. Will add 0.2 for PNG format and 0.5 for PDF
		pdf_width = nrow(in_table) * 0.5
		image_width = nrow(in_table) * 0.2
		if ( out_format == 'pdf'){
			ggsave(paste(paste(out_name,"ALLSAMPLESREADSPLOT",colname,sep="_"),out_format,sep="."),  width = pdf_width, height = pdf_height, units = pdf_unit, limitsize = FALSE)
		}else{
			ggsave(paste(paste(out_name,"ALLSAMPLESREADSPLOT",colname,sep="_"),out_format,sep="."), width = image_width, height = image_height, dpi = image_dpi, limitsize = FALSE)
		}		
	}

}
