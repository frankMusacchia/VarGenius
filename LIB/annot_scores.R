#VarGenius: 
# Author: Francesco Musacchia (2016)
#
#INPUT: the sample_file with the grouping of replicates
#       the filtered table with transcripts counts 
#				(We need filtered file of counts since we make the analysis
#				 only on the significant transcripts)
#	
#R CMD BATCH --no-save --no-restore '--args <input_folder>  <counts_filtered> <organism_name> <target_file> <minFDR> <out_suffix> <analysis_type> <annocript_out>' DE_analysis_generic.R

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples considered are 4 - please comment whenever you will not use all
length(args)
if (length(args) < 2){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <input_folder> <counts_filtered> <organism_name> <target_file> <minFDR> <out_suffix>  <analysis_type> <annocript_out> ' DE_analysis_generic.R")
}
###############################################METTERE IL FOLD CHANGE SOGLIA COME PARAMETRO IN INPUT

#Folder where all is taken and output goes
ann_f = args[1]
ann_f
col_mapper_ann = args[2]#Raw counts filtered before (not CPM)!!


mapping_file = args[3]
mapping_file
col_mapper_map = args[4]
col_mapper_map

ann_f = "Trio_VN_CGH06/annotate_out/Trio_VN_CGH06.hg19_multianno.csv"
col_mapper_ann = "Gene.refGene"
mapping_file = ""


ann <- read.delim(file=ann_f, head=T, comment.char='', sep=',')



