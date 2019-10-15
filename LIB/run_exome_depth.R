library(ExomeDepth)
library(GenomicRanges)
data(exons.hg19)
data(exons.hg19.X)
data(Conrad.hg19)

#R CMD BATCH --no-save --no-restore '--args <bam_paths> <fasta> <bedfile> ' var_stats.R


#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

analysis <- args[1] #analysis <- "SB_ID_A551_TRY"
analysis 
#File with list of BAM files to use
bam_paths <- args[2] #bam_paths<-"/pico/work/TELET_TIGEM/prove_vargenius/SB_ID_A551_TRY/finalout_out/cnv/exome_depth_bam_list.txt"
bam_paths 
#Fasta file for the Reference Genome
fasta <- args[3] #fasta<-"/pico/work/TELET_TIGEM/ngsworkspace/references/genomes//ucsc.hg19.fa"
fasta
#target BED file
bedfile <- args[4] #bedfile<-"/pico/work/TELET_TIGEM/ngsworkspace/references/targetbed//Agilent_ClearSeq_inherited_diseases_sex.bed"
bedfile
#Chr_type
chr_type <- args[5] #chr_type<-"sex"
chr_type
#Output folder
outfold <- args[6] #outfold<-"/pico/work/TELET_TIGEM/prove_vargenius/SB_ID_A551_TRY/finalout_out/cnv/sex/"
outfold

#bam_paths = "/pico/work/TELET_UDP/PROVE/exomedepth/probs_clinex1_bam_paths.txt"
#bam_paths = "/pico/work/TELET_UDP/PROVE/exomedepth/mothers_clinex1_bam_paths.txt"
#bam_paths = "/pico/work/TELET_UDP/PROVE/exomedepth/fathers_clinex1_bam_paths.txt"

#Read the bam paths list
dat <- read.table(file = bam_paths, header = FALSE)
my.bam <- as.vector(dat$V1)

#fasta = "/pico/work/TELET_UDP/home/shared/references/genome/ucsc.hg19.fa"
#fasta
#bedfile = "/pico/work/TELET_UDP/PROVE/exomedepth/clinical_exome_clean.bed"
#bedfile

#Open the BED file
bedfile.text <- read.table(file = bedfile, header = FALSE, fill = TRUE)

#Generate an object with ranges
exons.hg19.GRanges <- GRanges(seqnames = bedfile.text[,1],IRanges(start=bedfile.text[,2],end=bedfile.text[,3]),names = bedfile.text[,4])
colnames(bedfile.text)<-c("chromosome","start","end","names")
bedfile.text$chromosome <- as.character(bedfile.text$chromosome) 

#Exome depth command
if (chr_type == 'auto' ){
	my.counts <- getBamCounts(bed.frame = exons.hg19,bam.files = my.bam,include.chr = TRUE,referenceFasta = fasta)
}
if (chr_type == 'sex' ){
	my.counts <- getBamCounts(bed.frame = exons.hg19.X,bam.files = my.bam,include.chr = TRUE,referenceFasta = fasta)
}
#my.counts <- getBamCounts(bed.frame = bedfile.text,bam.files = my.bam,include.chr = FALSE,referenceFasta = fasta)


Complete.df <- as(my.counts[, colnames(my.counts)],'data.frame')
ExomeCount.mat <- as.matrix(Complete.df[, grep(names(Complete.df),pattern = '*.bam')])
nsamples <- ncol(ExomeCount.mat)

message('Now Looping through samples')

for (i in 1:nsamples) {
  #Set the name of the output file
  output.file <- paste0(outfold,'exomedepth.',colnames(ExomeCount.mat)[i],'.CNVs.tsv')
  
  #### Create the aggregate reference set for this sample (DEFAULT)
  my.choice <- select.reference.set (
  test.counts = ExomeCount.mat[,i],
  reference.counts = ExomeCount.mat[,-i],
  bin.length = (Complete.df$end - Complete.df$start)/1000,
    n.bins.reduced = 10000)
    
    
  col <- my.choice$reference.choice
  my.reference.selected <- apply(X = ExomeCount.mat[, col, drop = FALSE], MAR = 1, FUN = sum)
  
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth', test = ExomeCount.mat[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  #### Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = Complete.df$space,
                        start = Complete.df$start,
                        end = Complete.df$end,
                        name = Complete.df$names)
                        
                        
  #### Get rid of pesky "chr" that prevents subsequent intersect
  id <- gsub("chr","",all.exons@CNV.calls$id)
  all.exons@CNV.calls$id<-id
  chromosome <- gsub("chr","",all.exons@CNV.calls$chromosome)
  all.exons@CNV.calls$chromosome<- chromosome

  # To prevent stopping the whole Loop if only one bam does not produce any CNVs
  if(length(all.exons@CNV.calls)==2) {
    cat ("No CNVs detected at ",output.file)
    next
    }

  #### Annotate the ExomeDepth object
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.5,
                             column.name = 'Conrad.hg19.0.5')
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.0001,
                             column.name = 'Conrad.hg19.0.0001')
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.hg19.GRanges,
min.overlap = 0.0001,
                             column.name = 'exons.hg19')
  #### Write to TSV file
  cat ("Writing ",output.file)
  write.table(file = output.file, x = all.exons@CNV.calls,
              row.names = FALSE, quote = FALSE, sep = "\t")
}

