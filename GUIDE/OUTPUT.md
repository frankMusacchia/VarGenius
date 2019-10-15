# VarGenius - Output

## The sample sheet

The sample sheet can be created manually or with the PERL script get\_sample\_sheet contained in the *VarGenius* folder that basically uses the sample folder to fill the sample sheet with all the needed information. It has many option that you can read typing 

```
perl  VarGenius/get_sample_sheet.pl --help 
```

Here we see an example of command that you can use to obtain your first sample sheet. 
Go into the working directory: cd vargenius_analyses and test the command.

```
perl /home/users/frank/VarGenius/get_sample_sheet.pl 
-h /home/users/frank/VarGenius/CONFIGURATION/ss_head_file.txt 
-dir /home/users/frank/samples/ 
-m exome 
-t target.bed 
-uid 1 
-rgid 1
-cgn MY_ANALYSIS
-o /home/users/frank/vargenius_analyses/ss_first_exome_run.tsv 
```

With this command we take the sample sheet header file from -h, the samples folders from -dir, we write that all the samples are from exome sequencing, the target file is clinical_exome.bed, the user and research group identifiers are both 1, the output sample sheet path name is -o, the analysis name is indicated with -cgn and is MY_ANALYSIS.
Please use only these parameters for just one analysis for your first run of *VarGenius*. 
We reccomend that the sample name is identical to the folder name containing the fastq files. Since we are using Illumina technology *VarGenius* expects that fastq files contain the string L001,L002, etc.. to indicate the lanes.
While using this script you will understand how to build sample_sheets more complex to run hundreds of analysis with a single command.

After the sample sheet is created you need to read it using a Linux editor (like nano or VI) and verify that the data is correct.

The sample sheet has the following mandatory columns:
 - heritability: if the disease is heritable or not (0 or 1);
 - gender: gender of the sample (M,F);
 - diseaseid: identifier of the disease (not available in this version, use '-');
 - kinship: kinship of the sample (P: proband, M: mother, F: father, A: affected non proband (brother, syster or cousin of proband), R: any relative inserted for any reason);
 - userid: identifier of the user that is running the software (a number chosen by you);
 - resgroupid: identifier of the research group that is running vargenius (a number chosen by you);
 - fqdir: full path of the folder where the fastq is located;
 - fq1: fastq file name of read 1;
 - fq2: fastq file name of read 2;
 - readfname: a read file name to use into the DB;
 - samplename: a sample name to use into the DB ;
 - groupname: a name for the analysis to use into the DB (groupnames must be uniq!);
 - reference: this software can be run only with "hg19";
 - targetbed: name of the BED file with the target regions (this BED files must be all located in the same folder SEE BELOW);
 - seqtype: sequencing type ('targeted' or 'exome')
 - perchrom: if the analysis must be run per-chromosome or not (0 or 1, default 0);
 - convscores: if the base scores must be converted. Implies the use of **seqtk seq** (0 or 1). Use only if you have to convert scores to Phred.
 - infreq: if this analysis should be kept in count for the variant frequencies calculation.



## How *VarGenius* manages folders and files

Once you have installed *VarGenius* it will save all the analyses output into your preferred working directory (vargenius_analyses). Each analysis will have a folder with the same name of the analysis. The analysis folder contains the output from each of the tasks that *VarGenius* executes:
 - qc_out
 - alignment_out
 - refine_out
 - varcall_out
 - genotype_out
 - varfilt_out
 - phasing_out
 - finalout_out

N.B. We call **tasks** the main tasks of the pipeline and **steps** the steps contained in each task. We use the word **job** to indicate the jobs that are executed from the cluster and that in *VarGenius* execute a single task.

To the aim of having a software  which can execute different pipelines, we designed an algorithm which builds the output names using the steps names included in the tasks. Steps names are human readable so that is very simple to understand which step was ran to obtain a certain output. Moreover file names change depending by the fact that the task is executed on the *read files*, on *samples* or on the complete *group* of samples of the analysis.

Suppose you have 
For example the lane 1 raw reads fastq file of the sample **SampleA** has been trimmed, paired ends merged and finally aligned with BWA-mem. We obtain the following name:

vargenius_analyses/SampleA/alignment_out/SampleA_L001_trim_mergepe_align.bam

In a successive task the reads coming from different lanes are merged into a single file and GATK HaplotypeCaller is used for variant calling:

vargenius_analyses/SampleA/varcall_out/SampleA_mergegrp_varcall.g.vcf

Then, the sample VCF is genotyped with GenotypeGVCF and filtered with GATK VariantFiltering

vargenius_analyses/SampleA/varfilt_out/SampleA_genot_varfilt.vcf

Finally the output is produced and saved into:

vargenius_analyses/SampleA/finalout_out/SampleA.hg19.final_out_resgr14.txt


The **finalout_out** folder contains all the files that are used for the visualization. In particular we have
 -  the files used for the website:
   - index.html, coverage.html, qc.html

 - SampleA_allvar.hg19.final_out.txt: the raw output. It will be rearranged and saved into the next file
 - SampleA.hg19.final_out_resgr14.txt: the final tabular output from VarGenius. The resgr14 suffix is useful to indicate the research group that make the analysis.
 - SampleA.hg19.final_out_resgr14.xlsx: the XLS version of the final tabular output
 - SampleA_allvar.hg19_multianno.txt: the output of the annotation as it comes from Annovar
 - SampleA_igv_session.xml: IGV session file

The other useful outputs are
 - FastQC output files contained in **qc_out** folder
  - fastqc out of reads files finishing with fastqc.html
  - fastqc out of trimmed files finishing with val_1.fq_fastqc.html
 - The BAM files of the alignment of the reads against the genome (useful to check alignments in IGV). They are in the **alignment_out** folder.
  - For example a file finishing in align_mrdup_sortidx.bam represents the bam aligned with duplicates removed and reads sorted and indexed
 - The VCF files contained in the **varfilt_out** folder. They are the following:
   - SampleA_genot_varfilt_INDEL_filt.vcf: the filtered VCF for INDELs
   - SampleA_genot_varfilt_SNP_filt.vcf: the filtered VCF for SNPs
   - SampleA_genot_varfilt.vcf: the filtered VCF obtained merging filtered INDELs and SNPs.

VCF and BAM files paths are written into the IGV session file that can be open automatically in IGV (read IGV user guide).


The **stats_out** folder contains the output of statistics that have been computed on the BAM files. Specifically we use:
- **samtools flagstat** to compute the statistics about the alignment for:
 - the BAM file of the basic alignment against the genome
 - the BAM file of the reads aligned, duplicates removed and  sorted
 - the BAM file of the reads aligned, duplicates removed and sorted and intersected against the target file
- The pipeline designed by Aaron Quinlan [here](https://twitter.com/aaronquinlan/status/421786507511205888) to get the non-covered regions into the file: SampleA_XXX.noncov20 (where 20 means that the regions are covered with less than 20 reads and the parameter can be changed in the configuration file with *noncov_threshold*)
- **GATK DepthOfCoverage** to get the coverage per-gene (.sample_gene_summary) and per-interval (.UD\_MO021\_TRIO\_DOC.sample\_interval\_summary). Go the the [DepthOfCoverage page](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php) to obtain a description of the command and to read about its useful output.


## Log files

Log files are of two types. One is general and contains a description of the overall execution. All the jobs running for a given analysis write into this file.
The second type of log file is the one associated to the task (job log). That file contains the logs from the programs executed at each step of a given task from the job.

The general log is contained into the folder vargenius\_analyses/Sample1/data/. It has the name of the analysis id from the database and the date of execution (e.g. 616\_Thu\_2016\_XX\_XX.log ).
The job logs are contained into the task folders in a specific log folder which contains a log file for each of the jobs executed per-readfile, per-sample or per-analysis. Job logs do not contain the date into the name so that it is not puzzling to understand what happens after several runs, while the general log contains the data so that one is able to keep in count different logs from different runs of the pipeline.
Including the data into the log file is pertinent with the possibility to execute different tasks for the same analysis in different times.




## What does *VarGenius* gives in the tabular output


*VarGenius* generates an output in text format and one in XLS format. Both of them contain tab separated values. The final output that we describe here is the we cited as the *rearranged* output (SampleA.hg19.final_out_resgr14.txt). It contains three sections of fields/results. The first section are the basic information about the variant (chrom, pos,ref, alt etc..), the second is the genotype information and the third the annotation. Th first and second are both extracted from the VCF file. The thirs is obtained from the output of Annovar.

**Variant information**

Chrom: chromosome where the variant has been discovered
Pos: position in the chromosome
Ref: reference nucleotide or set of reference nucleotides
Alt: alternative nucleotide or set (if GATK found more alternative alleles this have been separated)
Filter:  PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. (SEE BELOW)
Qual: phred-scaled quality score for the assertion made in ALT. High QUAL scores indicate high confidence calls.
MQ (RMS Mapping Quality): this annotation provides an estimation of the overall mapping quality of reads supporting a variant call. 

**Samples genotype information**


Sample_XX_GT: genotype
Sample_XX_ZYG: zygosity
Sample_XX_REF_reads: allele depth for reference allele. The number refers to the total number of reads with each allele (including reads not used for calling due to low quality or whatever)
Sample_XX_ALT_reads: allele depth for alternative allele. The number refers to the total number of reads with each allele (including reads not used for calling due to low quality or whatever)
Sample_XX_DP: read depth at this position for this sample. I.e. the number of reads used to call the variant.

N.B. *Multiallelic variants*
For those variants which result to be multi allelic we used a different approach to show the genotype. Actually the genotype given in the FORMAT field of the VCF is not informative and must be inferred from the PL field. Even though the PL field is not inserted in the final output we have stored it into the database. 
It is possible to define a new GT field using the three values from PL. Those values score three possible genotypes: 0/0, 0/1 and 1/1. Once that the genotype is inferred, the phase will be lost, if present.

**Annotation**

The annotation is separated in: gene based, filter based, region based annotation as in Annovar.

*Gene Based annotation (http://annovar.openbioinformatics.org/en/latest/user-guide/gene/)*
Gene_refgene: RefSeq gene symbol

Func_refgene:  type of mutation (possible values are: nonsynonymous SNV, synonymous SNV, frameshift insertion, frameshift deletion, nonframeshift insertion, nonframeshift deletion, frameshift block substitution, nonframshift block substitution)

In this point Annovar gives a detail of the transcripts where the variant was found (NM_152486:exon8:c.707-25A>G). This information contains which transcripts includes the variant and which exon. Most probable result is the first that we decompose in multiple fields separating using the colon to the aim of the downstream filtering and analysis. Learn more [here](http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-2-refseq-gene-annotation). The following are the field that we obtain:
 - refgene\_refseq: RefSeq identifier of the transcript involved;
 - refgene\_exon_number: Exon number
 - refgene\_nucl\_change: nucleotidic change;
 - refgene\_aachange:the aminoacidic change in the gene and in the transcript only for exonic variants;
 - genedetail\_refgene: contains all the remaining predicted transcripts information as they are given from Annovar;
 - splvar\_dist\_from_exon: in those cases where the variant external or close to the exon it could be a splicing variant and this field indicates the distance from the exon (negative if it is before, positive if after)


*Filter Based annotation (http://annovar.openbioinformatics.org/en/latest/user-guide/filter/)*
This annotation permits to compare the variants in our experiment with those in other public (or not) databases.

exac\_XXX: latest Exome Aggregation Consortium dataset with allele frequencies for different populations (ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian))

esp6500siv2\_all:  ESP is a NHLBI funded exome sequencing project aiming to identify genetic variants in exonic regions from over 6000 individuals, including healthy ones as well as subjects with different diseases. 

1000g2015aug\_all: latest 1000 Genomes Project dataset with allele frequencies in six populations including ALL, AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian)

IVF: This is the internal variants frequency computed by *VarGenius* keeping in count separately the frequency of variants for whole exome and targeted sequencing. This frequency is computed for all the analysis for which the **infreq** parameters was set to 1 into the sample sheet and the variants have been imported into the database.

freq_factors: are the factors used for the IVF frequency calculation separated with an unique separator: number of variants found heterozygous, number found homozygous and total number of variants for which a genotype has been calculated. See below (Variants allelic frequency).

*Clinvar Annotation*: NCBI has three relatively new online resources for information about genetic tests, genetic conditions, and genetic variations:
 - MedGen – a medical genetics portal that focuses on information about medical conditions with a genetic component
 - ClinVar – an archival database that contains reported assertions about the relationship between genetic variations and phenotypes
 - The Genetic Testing Registry, or GTR – a registry of genetic tests for heritable and somatic changes in humans

clinsig: refers to Variant Clinical Significance, including unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response,histocompatibility, other. 
clndbn: refers to Variant disease name
clnacc: refers to Variant Accession and Versions
Clndsdb: database included in ClinVar 
Clndsdbid: database id

*Region Based annotation (http://annovar.openbioinformatics.org/en/latest/user-guide/region/)*

The classification of the variant based on the score is present only in tha final output of *VarGenius*. Since these scores are often different, *VarGenius* makes a "translation" of the classification symbol from the one present to one among those in the [ACMG-AMP classification guideline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908185/): P: pathogenic, LP:likely pathogenic, U: unknown, LB: likely benign, B:benign.

Following are all the included fields of annotation from Annovar output. 

Functional prediction of variants in whole-exome data with dbnsfp which contains data from: SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores, but ONLY on coding variants.

**N.B. The scores have been adapted to have a uniq format: P:pathogenic; LP: likely pathogenic; B: Benign; LB: Likely benign; U: unknown**


*sift\_pred*: SIFT predicts whether an amino acid substitution affects protein function. SIFT prediction is based on the degree of conservation of amino acid residues in sequence alignments derived from closely related sequences, collected through PSI-BLAST. SIFT can be applied to naturally occurring nonsynonymous polymorphisms or laboratory-induced missense mutations. (http://sift.jcvi.org/)
sift_pred

polyphen2\_hdiv\_pred (Polymorphism Phenotyping v2) is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. (http://genetics.bwh.harvard.edu/pph2/)
There are two databases for PolyPhen2: HVAR and HDIV. They are explained below (http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#summary-of-databases):
ljb2\_pp2hvar should be used for diagnostics of Mendelian diseases, which requires distinguishing mutations with drastic effects from all the remaining human variation, including abundant mildly deleterious alleles.The authors recommend calling "probably damaging" if the score is between 0.909 and 1, and "possibly damaging" if the score is between 0.447 and 0.908, and "benign" is the score is between 0 and 0.446.
ljb2\_pp2hdiv should be used when evaluating rare alleles at loci potentially involved in complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural selection from sequence data. The authors recommend calling "probably damaging" if the score is between 0.957 and 1, and "possibly damaging" if the score is between 0.453 and 0.956, and "benign" is the score is between 0 and 0.452.

lrt_score (log likelihood ratio) Using a comparative genomics data set of 32 vertebrate species we show that a likelihood ratio test (LRT) can accurately identify a subset of deleterious mutations that disrupt highly conserved amino acids within protein-coding sequences, which are likely to be unconditionally deleterious. The LRT is also able to identify known human disease alleles and performs as well as two commonly used heuristic methods, SIFT and PolyPhen. Details [here](http://www.ncbi.nlm.nih.gov/pubmed/19602639)

mutation\_assessor: Mutation assessor uses a multiple sequence alignment (MSA), partitioned to reflect functional specificity, and generates conservation scores for each column to represent the functional impact of a missense variant. 

fathmm_pred: A high-throughput web-server capable of predicting the functional consequences of both coding variants, i.e. non-synonymous single nucleotide variants (nsSNVs), and non-coding variants (http://fathmm.biocompute.org.uk/).
If a score is smaller than -1.5 the corresponding NS is predicted as "D(AMAGING)"; otherwise it is predicted as "T(OLERATED)". All the four missense variants were predicted as tolerated by FATHMM. The method is less well known compared to SIFT and PolyPhen, but in our experience it works really well and better than SIFT/PolyPhen (annovar documentation).

provean\_pred: PROVEAN (Protein Variation Effect Analyzer) is a software tool which predicts whether an amino acid substitution or indel has an impact on the biological function of a protein.
PROVEAN is useful for filtering sequence variants to identify nonsynonymous or indel variants that are predicted to be functionally important. Details [here](http://provean.jcvi.org/index.php).

vest3_score: VEST (Variant Effect Scoring Tool): is a machine learning method that predicts the functional significance of missense mutations observed through genome sequencing, allowing mutations to be prioritized in subsequent functional studies, based on the probability that they impair protein activity. Details [here](http://karchinlab.org/apps/appVest.html).

cadd13\_rawscore and cadd13\_phred and CADD is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome Details [here](http://cadd.gs.washington.edu/).

dann: score same as CADD but uses Neural Networks.

metasvm\_pred: METASVM: Accurate deleteriousness prediction of nonsynonymous variants for distinguishing pathogenic mutations from background polymorphisms in whole exome sequencing (WES) studies.

gerp\_elem: GERP: identifies constrained elements in multiple alignments by quantifying substitution deficits. These deficits represent substitutions that would have occurred if the element were neutral DNA, but did not occur because the element has been under functional constraint. We refer to these deficits as "Rejected Substitutions". Rejected substitutions are a natural measure of constraint that reflects the strength of past purifying selection on the element.

Avsnp150: this are reformatted dbsnp rs identifiers.

Cytoband: Giemsa-stained chromosomes bands. column in the output file below represent cytogenetfreq_factorsic bands. For example: 1p36.33.

Genomicsuperdups: describe variants in segmental duplications. The "Name" field in output represents the other "matching" segment in genome (which is located in the same chromosome at chr1:13142561). The "Score" field is the sequence identity with indels between two genomic segments (fracMatchIndel in the UCSC database table).(Example. Score=0.99448;Name=chr1:323738-chr5:180750353)

cpgislandext
rmsk
simplerepeat

Dgvmerged:  annotate deletions and duplications and compare them to previously published variants in Database of Genomic Variants (DGV) (Example: Name=dgv2n100-esv2758911-nsv945706-nsv469759)


Spidex
Spidex is a tool that first examines RefSeq transcripts to locate exons on the reference genome (hg19) whose inclusion levels may be affected. For each exon, the tool extracts 1393 features from proximal DNA sequence and uses a computational model to predict the percent of transcripts with the exon spliced in (PSI) for each of 16 human tissues, using both the wildtype (reference genome) and mutated sequences. 
The tool reports:
Dpsi_max_tissue: the maximum mutation-induced change in PSI across 16 tissues
Dpsi_zscore: a regulatory score obtained by aggregating the change in PSI across the tissues


## Gene Annotation

*VarGenius* associates gene information by using the gene symbol in refSeq. More than one gene symbol may be present separated by semicolon. Only for the first will be associated an information. To get this table the parameter gene_annotation must be set to YES.

The following are the annotations of the genes:

**GDI**: [REF](http://lab.rockefeller.edu/casanova/GDI)
The gene damage index (GDI) is the accumulated mutational damage of each human gene in healthy human population, based on the 1000 Genomes Project database (Phase 3) gene variations of healthy individuals and of the CADD score for calculating impact. They have shown that highly damaged human genes are unlikely to be disease-causing and yet they generate a big proportion of false postive variants harbored in such genes. GDI is very effective to filter out variants harbored in highly damaged (high GDI) genes that are unlikely to be disease-causing.

**RVIS**: [REF](http://genic-intolerance.org/)
The Residual Variation Intolerance Score measures genetic intolerance of genes to functional mutations. This scoring system predicts the expected amount of common functional variation based on the total amount of variation in each gene. The intolerance score itself is a measure of the deviation from this prediction.


**OMIM**: [REF](http://www.omim.org/)
The mendelian Inheritance in Man (OMIM) is a compendium of human genes and phenotypes that is freely available. OMIM contanis information about all known mendelian disorders an over 15000 gnes and focuses on the relatinship between phenotype and genotype.


# VariantFiltration parameters

The parameters used in VariantFiltration are the same suggested in the GATK best practices [here](http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)

We set 4 filterExpressions and then associated the expression to the SNP or INDELS using the two parameters snp\_filters and indel\_filters.
As we can see only the HARD\_TO\_VALIDATE filter changes between the two variants types.

snp\_filters = 1,2,3
indel\_filters = 4,2,3
clusterWindowSize = 10
filterName1 = HARD\_TO\_VALIDATE
filterExpression1 = **QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0**
filterName2 = VeryLowQual
filterExpression2 = **QUAL < 30.0**
filterName3 = LowQual
filterExpression3 = **QUAL > 30.0 && QUAL < 100.0**
filterName4 = HARD\_TO\_VALIDATE
filterExpression4 = **QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0**

Where:
QD: Variant Confidence/Quality by Depth
FS: Phred-scaled p-value using Fisher's exact test to detect strand bias
MQ: RMS Mapping Quality
MQRankSum: Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
ReadPosRankSum: Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias


### HTML output with statistics

*VarGenius* gives in output a little website composed by few pages showing statistics. There are three pages:
- **index.html** contains the links to the final output and the VCF file. Moreover it shows statistics coming from the analysis:
 - Total Reads on X and Y chromosomes
 - Zygosity Frequencies with X chromosome
 - A table comparing the number of homozygous and heterozygous variants in the different samples when the joint analysis is used. This information is useful when lookin for the segregation in family analysis.
 - number of existing variants in dbSNP. This value is useful to understand if there are errors in the call.
 - Zygosity Concordance: used to understand if the zygosity call is optimally correlated with the amount of ALTernative reads
- **qc.html** contains the links to all the HTML files output  of FastQC. They are used to check if the sequencing was correct before the analysis.
- **coverage.html** Contains the links to files contained in *stats_out* folder. They are the output of programs to compute the coverage (see below) and the output of the samtools flagstats program that computes statistics about the alignment.

## PLOTS

Into the **coverage.html** page *VarGenius* prints many plot showing statistics about the coverage. Some plots are pertinent with the current analysis while others come from the comparison with all the analyses performed with the same kit by the same profile.
The configuration file permits to select one or both of them using the parameters:

get_coverage_plots and update_overall_samples_plots


# Global coverage plots

The "global" plots show statistics about the coverage for all the samples sequenced with the same kit and the same profile.
The HTML page shows the following plots for global statistics:
 - Global samples coverage: at different levels of coverage shows the coveage level of each sample. It is useful to verify the coverage of different sequencing runs performed with the same kit.
 - Global samples ProperlyPaired: it is the plot showing the number of reads properly paired and aligned for each sample. This statistic is obtained using samtools flagstat.
 - Global samples TotReads: it is the plot showing the total number of reads for each sampe. This statistic is obtained using samtools flagstat.

# Plots per-analysis

Plots per analysis show statistics picking the percentage of locus, for each gene, which is covered more than a given threshold given in the configuration file (**min_gene_coverage** parameter, default: 10X). 

- Coverage per disease gene panel: this plot shows a boxplot for each sample into the analysis. Each boxplot has a value for each gene of a specific diseases panel and shows the percentage of gene covered for more than the **min_gene_coverage**. It is useful to check if there are genes related with a specific pathology with low coverage. This plot is correlated with the table present into the **Coverage** HTML page named **Disease Genes with Coverage<99**. This table has three columns: gene name, panel and the percentage of the gene which is covered less than 99%. You can use this table to check if there are specific genes with a particularly low coverage.

- Plots of TotReads VS  Percentage of Alternative Alleles: this plot is shown for homozygous, homozygous reference and heterozygous SNPs (and INDELs). They correlate the number of alternative alleles found with the genotype identified in GATK. Hence, they are useful to check the correctness of the genotype inferred.
For homozygous non-reference genotypes the cloud of point should be flatten towards the right of the plot, while for reference ones the cloud should be flatten towards the left.
For heterozygous genotypes the cloud of points has a gaussian shape in the middle of the plot.
Generally, INDELs are less than the SNPs, hence comparing the clouds they have same shapes but different density.



## Variants allelic frequency

*VarGenius* computes the allelic frequency of all variants present into the database, generating an Internal Variant Frequency (IVF). The command performing the calculation (recompute\_freqs) must be manually run every time a set of samples have been analyzed. For the calculation is used a query to get the number of times that the variant is found heterozygous and homozygous.
Then the formula to obain the allelic frequency of the variant V is
	
```
V\_all\_f = V\_het + (V\_hom * 2) / (Tot\_gen * 2)
```
						
Where we multiply V\_hom by two because the homozygous variant was found in both alleles and Tot\_gen because the samples have two alleles. We use only those variant for which GATK is able to compute the genotype.

## Coverage

*VarGenius* gives many informations about the coverage.

# Sample, Intervals and Genes coverage with DepthOfCoverage

DepthOfCoverage is used to analyze the coverage of the samples, the intervals of the target file and of the genes covered by the target. The output files from DepthOfCoverage are written into the **stats** folder of your analysis and start with *AnalysisName_DOC_XXX*.

Inside the **stats** folder you will find:
  - *Analysis_DOC_.sample_XXX_summary*: for each region of the target and for each sample of the analysis shows the percentage of reads above a certain percentage of coverage. This statistic is present for **sample**, **interval** and **gene** coverage.
  - *Analysis_DOC.sample_cumulative_coverage_proportions*: for each sample shows the percentage of the regions of the samples covered with more than Y reads (gte_Y).

Please refer to the [GATK DepthOfCoverage](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php) page for supplementary information.

Coverage thresholds to be used in GATK DepthOfCoverage can be changed with the parameter: **DOC_coverage_thresholds** 
using comma separated values.

# Non-covered regions 

Another useful information provided by *VarGenius* is a table showing the regions of the target with poor coverage.
We define as poor covered regions those which are covered by less than 20 reads. This threshold can be changed into the configuration file.

The table is written into the **stats** folder for each sample and named as *SampleName_XXX.noncov20*. 
It contains the chromosome, start, end transcript (or gene name).
For each line we have a region that is composed of all the adjacent regions (without gaps).
Additional information in the table is:
 - MinCov: is the minimum coverage among all the adjacent merged regions of the interval;
 - MaxCov: is the maximum coverage among all the adjacent merged regions of the interval;
 - AverageCov: is the average coverage among all the adjacent merged regions of the interval;
 - GeneSymbol: is the gene containing the transcript name shown in the Locus (where present).

```
CHR	START	END	Locus	MinCov	MaxCov	AverageCov	GeneSymbol
chr1	1535211	1535212	NM_001242659_cds_0_0_chr1_1534715_r	19	19	19	C1orf233
```
The method to extract the non covered regions is the one described by Aaron Quinlan [HERE](https://twitter.com/aaronquinlan/status/421786507511205888).
Pseudo code:
 - Use genomecov to get a BedGraph for the input bam
 - Use AWK to get only those lines where the coverage is less than the threshold
 - Merges the regions with bedtools merge
 - Finally intersects the regions with the target

You can change the coverage level threshold using the parameter: **noncov_threshold**.

###  Database Creation

*VarGenius* creates a PosgreSQL database the first time you run it. Three sets of tables are present: 1. those regarding the management of the analyses; 2. the tables with variants information and 3. the tables with the gene's information. 
While the first two are filled during your analyses, the latter is constructed the first time you run *VarGenius* and it allows to get a link between the variant detected and the genes. Genes information are downloaded from multiple sources (RefSeq, HPO, GDI, RVIS) and stored into tables using the command -genedb. You can keep this tables up to date by running again the database generation task.


---------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/vargenius
