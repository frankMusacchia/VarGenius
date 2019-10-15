#This function takes in input:
	#the gene summary table, the selected genes, the samples names, the output names
# and builds a table with only the column of gene coverage above 10 for each sample
# and for the selected genes

#R CMD BATCH --no-save --no-restore '--args <in_table> <sel_genes> <samples_str> <out_name>' r_utils.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

#Output file name
in_table <- args[1] #"chr.buttami"#
in_table
#Output from VarGenius
sel_genes <- args[2] 
sel_genes 
#List of samples separated by comma
samples_str <- args[3] #c("UD15","UD16","UD17") #
samples_str

#The coverage percentage field to use
cov_p_field <- args[4]
cov_p_field

#Output name
out_name <- args[5]
out_name

length(args)
if (length(args) < 5){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <in_table> <sel_genes> <samples_str> <cov_p_field> <out_name>' r_utils.R")
}



gc<-read.table(in_table,stringsAsFactors=F,head=T,sep="\t",quote="", check.names=FALSE)

if (sel_genes != "ALL_GENES"){
	sel <- read.table(sel_genes,stringsAsFactors=F,sep="\t",quote="")
	gc_sel <- merge(gc,sel,by.x="Gene",by.y="V1",all.y=T)
}else{
	gc_sel <- gc
}

#Split string with sample names 
samples_l <- unlist(strsplit(samples_str,",",fixed=T))


numcol = 1
#Put the Column with GeneName
df <- gc_sel[1]

for ( sample in samples_l ){
	field = paste(sample,cov_p_field,sep="_%_")#"above_10"
	field
	df <- cbind(df,gc_sel[field])
	numcol = numcol + 1
}

#Write output
write.table(df,file=out_name,sep="\t",row.names=F,col.names=T,quote=F)



