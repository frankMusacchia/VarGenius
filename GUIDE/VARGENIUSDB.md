# VarGenius - VarGenius Database


## Introduction
*VarGenius* creates a database at its first run. It allows to keep memory of information about the samples, the analyses execution, variants and their genotypes in each analysis and contains tables for the annotation of the genes.

Here is a brief explanation of each table that is created into the database:

  - **readfiles**: is the table containing information about the single unit of fastq files. If there are N fastq (for example from different lanes) than there will be N readfiles. If the samples is composed by a single fastq there is only 1 readfile. Each readfile has an identifier (readfid). This table contains the exact location of the fastq file in the system (fqdir) and the file names (fq1 and fq2). Moreover, contains flags indicating if a particular step has been run using this readfile.
  - **samples**: is the table containing information for each sample. Contains the information about the samples (gender, kinship, affected or not). Contains also the fields to indicate the steps executed.
  - **analyses**: contains information about all the analyses. It contains also the name of the pedfile, the target BED. Contains fields to indicate the steps executed and if the analysis is included into the database for the frequency calculation.
  - **sample_sheets**: contains the names of the sample sheet used for a specific run of the pipeline.
While the following are the tables used to store variants information. 
  - **target_files**: information about the target files used for the analyses.
  - **variants**: contains chromosome, position, reference  and alternative nucleotide
  - **var_freq_update_hist**: contains history of any update of variants frequency computation performed.
  - **statistics**: this table contains the scores obtained from GATK which define the quality of the variant. In particular contains the QUAL, INFO and FILTER fields from the VCF file
  - **genotype_sample**: is a table that for each sample from the analysis contains the genotype obtained for a specific variant. Hence, the format field is included (GT, AD, DP, GQ, PL, TP, PGT, PID).
  - **annotations**: is the table of the annotations of variants obtained with Annovar (almost 100 fields). The use of this table is OPTIONAL. It could be useful whenever someone is interested in doing queries per population. A specific step has been included in the pipeline to annotate all available variants and store the result into this table.
Finally we have the gene information:
  - **genes**: is the table containing information about the genes (gene name, hpo identifiers associated, GDI score, RVIS score, RefSeq, OMIM and Entrez identifiers of the gene.
  - **phenotypes**: the *genes* table contains the Human Phenotype Ontology (HPO) identifiers field with the identifiers of the this table. Each record contains the RefSeq identifier of the HPO phenotype and its description.
  - **transcripts**:  Each record contains the RefSeq identifier of the transcript and its description. The *genes* table contains the refseqids field with the identifiers of the this table.
  - **info**: this table keeps few information about the database. It is used in *VarGenius* for exclusive access of the database.

The genes information are downloaded from known sources:

The table with correspondences of the HPO phenotypes with genes_http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt

Accumulated mutational damage of genes is given by the Gene Damage Index (GDI) score. The table with values is at:

http://lab.rockefeller.edu/casanova/assets/file/GDI_full_10282015.txt

The variation intolerance score of genes is given by the Residual Variation Intolerance Score (RVIS). The table is downloaded from:

http://genic-intolerance.org/data/GenicIntolerance_v3_12Mar16.txt

Phenotypes associated to genes are taken by the Online Mendelian Inheritance in Man (OMIM) database. We used the one from:

http://omim.org/static/omim/data/mim2gene.txt

Transcripts associated to genes are obtained from RefSeq:

ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz 


#query_2_db.pl

The script query_2_db.pl has been created to allow users to execute automated queries to the db. While the VarGenius database can always be queried using SQL language, we have developed PERL functions into the module db_management.pm allowing to fetch data programmatically using PERL.

The script query_2_db can be used giving the argument -f and executing the following functions:

	"VARIANTS (given a list of variants returns the samples in which the variant is present with their genotype"
	" To use this function you need to give in input:"
	"-i: a file with a list of compid"
          					
	"GENE_PHEN (Given a list of sample names, returns, for each gene, the number of phenotypes which are involved with the gene."
	" To use this function you need to give in input:"
	" -i: a file with a list of samples for which you want to investigate the phenotypes"
								
	"FREQUENCIES (permits to obtain the frequencies of the variants into the database adjusting for relatedness or by subgroup of samples use also -i."
	" To use this function you need to give in input:\n"
	" -o: the output file"
	" -log: a log file"
        " -parameters: the filter that you want to use"
          
        "SAMPLES_GENE_COVERAGE: Given a list of sample identifiers and a gene name, returns for each sample the coverage of that gene"
        " To use this function you need to give in input:".
	"-i: a list of comma separated sample identifiers (you can use vargenius.pl -c XXX -sh\_s ANALYSIS_ID\n"
        " -p: the name of the gene for which you want to get the coverage for the samples in -i"
          
        "VARS_ON_GENE: Given a list of gene names, extracts a file with all the variants present on that gene"
        " To use this function you need to give in input:"
        " -p: the name of the gene for which you want to get the variants\n"


## Useful SQL queries

This queries can be executed from the Linux terminal with the psql command:

**Get all the analysis identifiers of the analyses for which the variants have not been inserted:**

```
psql -h postgresql.pico.cineca.it -U ttudp001 telethon -c "\COPY ( select analysisid,analysisname from analyses where analysisid in (select distinct analysisid from analyses where infreq=1 and analysisid not in (select distinct analysisid from genotype\_sample)) ) TO 'all\_anal\_to\_import.csv' with csv header"
```

**Get all phenotypes for all samples:**

```
psql -h postgresql.pico.cineca.it -U ttudp001 telethon -c "\COPY ( select hpoid,hponame,hpodef from phenotypes where phenid in (select distinct phenid from samples_hpo)) TO 'all\_phenotypes.txt' with csv header"
```

**Get the parents of a list of probands sequenced with a specific kit and with a specific owner:**

```
psql -h postgresql.pico.cineca.it -U ttudp001 telethon -c "\COPY ( select sampleid,samplename from samples where analysisid in (select analysisid from samples where sampleid in (4258,4259,4260,4262,4267)) and (kinship = 'M' or kinship = 'F') and analysisid in (select analysisid from analyses where userid=1 and infreq=1 and targetbed like 'XXX') order by samplename) TO '/../parents_of_prob_ids.txt' with csv header"
```



---------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/vargenius
