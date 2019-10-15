# VarGenius - Advanced Usage

## VarGenius pipelineS

In *VarGenius* we implemented the possibility to avoid given steps for different reasons. For example of targeted sequencing removing duplicates is not suggested or the VQSR steps of GATK are to be executed only for joint analysis.

Use the configuration file to set the step you want to execute
Mandatory steps are the following: 

In QUALITY CHECK AND TRIMMING
qc_exec
In ALIGNMENT
align_prog_exec and sam_sort_idx
In REFINEMENT
none is mandatory
In MERGE READFILES
steps are automatically executed if there are multiple read files. You don't need to take care of them
In VARIANT CALLING
varcall
In GENOTYPING
none is mandatory
In VARIANT FILTERING
varfilt
In PHASING
phasing is automatically executed if there is a proband and the two parents
In STATISTICS, OUTPUT MANAGEMENT and IMPORT VARIANTS AND ANNOTATIONS all the steps are optional

You can optionally run: 
In QUALITY CHECK AND TRIMMING
trimming and qc\_after_trim 
In ALIGNMENT
convert_scores, merge_pairs and mark_rem_dup_exec
In REFINEMENT
all steps are optional
In VARIANT CALLING
catvar is automatically executed and varrecal, apprecal are optional for the GATK VQSR procedure


N.B. *VarGenius* will not work if you ask only for non-consecutive tasks (e.g. "-qc -vc"  or "-aln -vf" ).  
## DEFAULT BEHAVIOURS
 1. If there are more than one analysis, VarGenius executes the update of frequencies only after the last analysis. This is hard-coded (vargenius.pl->pipeline_execution().
 2. If there are more than user_config.txt->min_samples_4_joint samples, the analysis id considered a joint analysis and a single PED file is created for all the samples. 
 3. If the analysis is of a TRIO (mother, father and a single proband) the phasing will be performed automatically. 


## HOW TO CREATE THE SAMPLE SHEET

The sample sheet can be created manually or with the PERL script get\_sample\_sheet contained in the *VarGenius* folder that basically uses the sample folder to fill the sample sheet with all the needed information. 

Here we see an example of command that you can use to obtain your first sample sheet. 
Go into the work directory (vargenius_analyses) and test the command:

```
perl /home/users/frank/VarGenius/get_sample_sheet.pl 
-dir /home/users/frank/samples/ 
-m targeted 
-t target.bed 
-uid 1 
-rgid 1
-o /home/users/frank/vargenius_analyses/sample\_sheet\_test.tsv 
```

After the sample sheet is created you need to read it using a Linux editor (like nano or VI) and verify that the data is correct.
For more details about the sample sheet read the [OUTPUT](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/OUTPUT.md) page.



## HOW TO SET Group3 OPTIONAL parameters

Optional parameters into the configuration file can be changed according with your system and substantially there are two groups of variables that you may want to modify.

The first group is **JOB RESOURCES REQUESTS** that are the resources requested by each job. Depending by your cluster you may want to modify RAM memory and maximum time (wall_time)
requested by QSUB to accelerate the jobs scheduling for each task.
To modify these parameters please consider the syntax of the configuration file:
 - put one space after '=' symbol
 - the ram request is always in the form: numberGB (e.g. 120GB)
 - the walltime request is always in the form: HOURS:MINUTES:SECONDS (e.g. 10:30:00)
 
 
The second group of variables indicate the folders to use for reference databases. You could have your own databases in specific folders and you can use them here.
If these variables are set, the default databases will not be used.

You can modify: 
 - **scratch_f**: The scratch area is a folder were temporary files are written and deleted at the end of the process (GATK uses it).
 - **storage_f**: The storage area is a location where you can store the old results of *VarGenius* in a folder different from the one that is used to make the analysis.  You must have privileges to read and write there.
 - **target_reg_f**: Together with the enrichment kit for the creation of the library you should find a BED file containing the target regions with which the reads will be aligned. All the target files that you use for your analyses must be put inside the same folder. Then tell to *VarGenius* what is the folder by setting:
 - **gatk_ref_f**: GATK uses known indels and snps to work with well assessed variants. VarGenius will search VCF files inside this folder

Moreover, here we have:
 - **html_host** if you want to use a web server to provide results to users, you can set this variable to the URL where the results are located (vargenius_analysis folder)
 - **disease_genes_lists_f**: these are predefined text files with lists of genes. each file is related with a specific disease or group of diseases. You may want to add your own panel for the downstream statistics.
 - **hum_gen_align**: the fasta file to use as a reference genome.



## INSTRUCTIONS TO GET THE FASTA REFERENCE HG19

You can download the Hg19 fasta reference at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
then you can follow the instructions to convert the 2bit format to fasta here (https://genome.ucsc.edu/goldenpath/help/twoBit.html). Downlad 2bitToFasta program here (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa)
```

Put the reference hg19 fasta into a newly created subfolder /home/reference/genome/
Then set the folder and the file with the two parameters:

```
ref_genome_f = /home/reference/genome/
hum_gen_align = ucsc.hg19.fa
```

```
## INSTRUCTIONS TO GET THE INDEX FILE FOR BWA AND GATK

To use the reference genome in BWA and GATK you need to go into the folder where the Hg19 fasta is located and do the following:
 1. Generate the BWA index: **bwa index reference.fa**
 2. Generate the fasta file index by running the following SAMtools command: **samtools faidx reference.fa**
 3. Generate the sequence dictionary by running the following Picard command: 



```	
java -jar picard CreateSequenceDictionary.jar \\
  REFERENCE=reference.fa \
  OUTPUT=reference.dict
```



## Indel realignment and base recalibrator
The following files are needed to run refinement steps in GATK.
Please download the following files (your preferred versions) from the **bundle** (ftp://ftp.broadinstitute.org/bundle/hg19/) into the just created folder gatk\_ref\_f and set the variable names correspondingly. 
You can use either .vcf or .vcf.gz format.

```
known_mills = Mills\_and\_1000G\_gold\_standard.indels.hg19.vcf.gz
known_db_snp = dbsnp\_147\_hg19.vcf.gz
known_hapmap = hapmap\_3.3.hg19.vcf
known_1000g = 1000G\_phase1.indels.hg19.vcf.gz
known_omni = 1000G\_omni2.5.hg19.vcf.gz
```

## Variant Quality Recalibration Steps*
The following files are needed to run the variant quality recalibration steps in GATK.
Please download the following files from the **bundle** into the just created folder gatk\_ref\_f and set the variable names correspondingly. 
N.B. From the bundle you will also need the corresponding *.idx* index files.

```
known_sites_mills = Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
known_sites_db_snp = dbsnp_147_hg19.vcf.gz
known_sites_hapmap = hapmap_3.3.hg19.sites.vcf.gz
known_sites_1000g = 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
known_sites_omni = 1000G_phase1.indels.hg19.sites.vcf.gz
```


## HOW TO RUN VARGENIUS USING A FILTERED VCF

Let's suppose you have obtained a VCF file (AnalysisX.vcf) for an analysis using a different program and now you want to  annotate it
with *VarGenius* :
- Create a sample sheet for this new analysis (AnalysisX). To do this, you have to create a fake sample folder and fake fastq files
 - mkdir /path/to/dir/AnalysisX
 - touch /path/to/dir/AnalysisX/AnalysisX_R1.fastq.gz
 - touch /path/to/dir/AnalysisX/AnalysisX_R2.fastq.gz

N.B. If you have more than a sample (for example a VCF for a trio analysis) you will need to create fake folders and fastq files for each sample and

N.B.2 Get the folder name from the sample name found inside the VCF
 
Now get the sample sheet file
 perl VarGenius/get_sample_sheet.pl -dir /path/to/dir/AnalysisX  -m exome -o ss_AnalysisX.tsv

Load the sample sheet into the database to create the analysis and its folder:

perl VarGenius/vargenius.pl -c user_config.txt -ss ss_AnalysisX.tsv

Put the VCF in the varfilt\_out folder of the analysis folder that has just been created. 

```
mv AnalysisX.vcf /path/to/dir/AnalysisX/varfilt_out/
```

Get the analysis id to use for the run

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt -aid AnalysisX
AnalysisX : 2251
```

Finally run the final output task and wait for the results:

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt -ra 2251 -fout
```

For the BAM files you should use those aligned, sorted and with removed duplicates:
- insert both BAM and BAI files in alignment_out folder. 
- Change the name accordingly with VarGenius pre-defined names (usually sample\_name\_sortidx\_mergebam\_rmdup.bam).
- run only gr_st because other statistics need you run the alignment steps.

## HOW TO RUN VARGENIUS IN JOINT GENOTYPE MODE

The Joint Genotype mode is always used in VarGenius for a group. In fact, if you have a family it will be genotyped together using GenotypeGVCF. Nevertheless the method can be used for hundreds of samples together as stated in this article (http://gatkforums.broadinstitute.org/gatk/discussion/3686/why-do-joint-calling-rather-than-single-sample-calling).

Joint genotype for multiple samples may be used to increase the probability to call a rare variant in a sample by using data from multiple.
To use this method for many different samples you need to produce a sample sheet with the samples to process all belonging to the same groupname which could be “Joint_analysis”.
In this case a new folder will be created to contain all the samples to be analyzed. Once that the Joint genotype has been called and the output created, VarGenius has compiled a table containing the data about the variants and the genotype information for each sample.

To separate the Genotype information  per-sample we can use the option in the user_config.txt file 

sep_\joint\_analyses = YES

which means that, starting from the "common" output a new output will be created for each sample.



## RUNNING MULTIPLE ANALYSES TOGETHER

*VarGenius* allows to start multiple analyses with a single command using one single sample sheet. The sample sheet may be obtained with the provided script "get_sample_sheet.pl" but many times it needs to be modified manually. For example if some samples are to be analyzed together they must have the same groupname. The groupname is the way *VarGenius* identifies a single analysis.

Please consider that often PostgreSQL databases have a limit of accesses and you must think that the final task of output production involves thousands of db accesses. Hence if you want to run hundreds of analyses together you may want to use the parameter **block_db = YES** so that each final task will access in an exclusive way.




## HOW TO ADD A SAMPLE TO THE JOINT ANALYSIS

To add an analysis to the joint sample safely and without re-running all the variant call for all the contained samples you must do two different steps.

You generate a sample sheet with the samples to insert belonging to the joining group name. Then you execute *VarGenius* without execution parameters but just the sample sheet. You must add the --force_run command, so that *VarGenius* won't block you for using an existing analysis name.

```
perl VarGenius -ss sample_sheet.tsv -c user_config.txt --force_run
```

By doing that, you obtain that a new sample has been added to the group and finally you can run vargenius by executing the analyses only for the just added samples. The analysis will be executed per-sample until the VariantCalling step. Then, starting from the genotyping the analysis will use all the output already present from the variant calling. This will overwrite the current GenotypeGVCFs output with an updated one.

To get what are those samples you can execute:

```
perl vargenius.pl -c user_config.txt -sh_s groupname 
```

Finally, start *VarGenius* using only the samples you just inserted. 

```
perl VarGenius -ra <analysis_id> --run_analyses <samples_ids> -c user_config.txt -qc -aln -ref -vc -vf -fout
```

## Variants Frequencies computation

*VarGenius* stores information about any variant detected into your analyses. Go to the [VARGENIUSDB](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/VARGENIUSDB.md) page for a survey of *VarGenius* database. Go to the [OUTPUT](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/OUTPUT.md) page to understand how allelelic frequency is computed.

The genotype_sample table is the table used to calculate the frequency of the variant because it contains the genotype information for each variant.

Variants are not automatically imported using the **--start** command because we observed that when running many analyses together the continuous parallel access of many jobs to the database may cause problems such as crossing the limit of Postgres accesses (that in many cases is under the administrator control).
The **--vcfimport** command must be used instead.

To import the variants you can decide to just import them or import and update the frequencies. To modify this behaviour change to YES or NOT the following parameters
into the configuration file:

```
# IMPORT VARIANTS AND ANNOTATIONS
vcfimport = YES
update_frequencies = YES
```

To only update frequencies (both exome and targeted) you can use also a single command:

```
perl VarGenius -c user_config.txt  --update_freqs
```

Our suggestion is to use the --vcfimport after the complete analysis (--start) if you are executing a single analsysis whilst you should use the following parameters
if you are executing many analyses togheter:

```
# IMPORT VARIANTS AND ANNOTATIONS
vcfimport = YES
update_frequencies = NO

# Block database for exclusive use
block_db = YES
```

The **block_db** parameter gives to the vcfimport process an exclusive access to the database. Other processes doing the same will wait it to complete.
When you see that all import jobs are completed you can finally run the **--update_freqs** command.


## MYSQL USEFUL QUERIES

# Get the analyses ids for which the frequencies have not been inserted into the database

The following SQL query permits in the PostgreSQL database to obtain the identifiers and names of the analyses that have not been inserted into the database. You need to specify the userid that in *VarGenius* corresponds with a set of similar analyses from the same person and with the same target file.

```
select groupid,groupname from groups where groupid in (select groupid from groups where userid=145 and infreq=1 EXCEPT select distinct groupid from genotype_sample) order by groupid;
```


## INTERSECTING TARGET FILE WITH CODING REGIONS

We have wrote a script to intersect the target bed file with coding regions downloaded from Genome Browser. 
In some cases for WES we decided to use only the coding regions to avoid puzzling variants.

This script is included into the *VarGenius* directory 'utilities' and its name is **intersect_target_with_cod_regions.pl**

You must use it from the vargenius working folder with:

```
perl [VarGeniusPath]/utilities/intersect_target_with_cod_regions.pl [target_file] [coding_regions.bed] [output_file]
```

Remember that the output folder MUST be the one you set into the configuration file with **target_reg_f**.

## FIX AN ANALYSIS FOR WHICH ONLY SOME SAMPLES ARE FINISHED WITH AN ERROR

If you executed a joint analysis and something went wrong with a task only for a specific sample X, you may want to run that specific step only for the sample X, right? Without wasting your time and cpus for the same process which were good before, right?

You can do this with *VarGenius* by using the --analyze_samples parameter. For example if for the sample X of the group G the refine step had some error you can just run

```
perl VarGenius -c user_config.txt  -ra G --analyze_sample X --refine
```

and it's all done!


## STATISTICS USING GENE PANELS

*VarGenius* may produce some plots and statistics related with specific panels of genes. Specifically we have a plot describing what is the coverage (at least 20X) of the genes of a specific panel (e.g. related with a specific disease).

The package auto-installs some gene panels into the vargenius_analyses/data/gene_panels/ folder and uses them if specified in the user configuration file with the variable disease_genes_lists_f.

If you want to add your own gene panel for the automatic construction of a new plot, you need to write a pure text file with a list of gene names separated by newline and give it the name of the disease (or whatever). The name will be used into the plot as it is hence we suggest to give it a name without extension (e.g. Epilessia).
Then, put it into the vargenius_analyses/data/gene_panels/ and add it to the list in disease_genes_lists_f (comma separated values).


## GATK EXECUTION PER-CHROMOSOME

GATK can be executed per-chromosome. This means that some of the operations performed by the kit can be executed using a single chromosome during the run. The result will be identical while the overall process can be decomposed in multiple processes which exploit the multiple processors of an HPC machine.

*VarGenius* will automatically separate your target file per-chromosome if you set perchrom=1 into the sample_sheet that you create for your analysis.


## STORE ANALYSIS INTO A STORAGE FOLDER

If your working area has limited space you may want to store the useful output files generated into a specific folder. 
First set the folder full path to use with the **storage_f** parameter into the configuration file. 
Second use the *VarGenius* **store_analyses** command to store a given set of analyses using their database identifiers.
You can separate intervals of identifiers with '-' and atomic with comma. (e.g --store_analyses 1-10,14,18,23-34)

All the files linked into the HTML web site will be stored into the given folder and the corresponding original files will be transformed in symbolink links to those.

You may want to check if an analysis is already stored into the database using the **stored** field of the **groups** table.

## PROFILES


To avoid to change continuosly a single configuration file and to avoid mistakes we created the **profiles**. Profiles are text files that basically re-set some of the existing configuration parameters. They allow to change any parameter at once.
They must have a specific name and must be located into the folder CONFIGURATION of *VarGenius*. 
The name of this files must be *user_XX.profile* where *user* means that a specific user may want different types of analyses. The *XX* number is the user id is associated with each sample and analysis and must be set into the sample sheet.

If you do not create it, do not care, nothing will happen and default parameters will be used.


## HOW TO USE THE LOG FILES TO FIND THE ERRORS DURING THE EXECUTION

*VarGenius* gives the possibility to monitor what happens step by step to your analysis by writing two types of log files.

We write one main log file into the ANALYSIS_NAME/data folder whose name depends by the database analysis identifier (numeric) and the current date and time.
In addition, all the tasks folder contain a LOG folder which contains the log file for the job executed for the given step.



## HOW TO REDO AN ANALYSIS OR DOING IT MULTIPLE TIMES WITH DIFFERENT PIPELINES

You can do an analysis more than once (e.g. different parameters or pipelines, single analysis or family) by just giving it a different analysis name into the sample sheet. 

N.B. You must be aware that  if you run two analyses of the same sample(s) the variants imported into the database will be almost identical and will confuse the variants allelic frequency. To this aim we use the parameter **infreq** into the sample sheet to inform *VarGenius* if the variants of a given analysis must be kept in count for the computation of the frequency. Use infreq=0 if you do not want to consider the current analysis.

## HOW TO EXECUTE AN ANALYSIS WHERE ONE OR MORE SAMPLES ARE DIFFERENT SEQUENCING TYPES (E.G. SampleA is exome, SampleB is amplicons panel)

Whether you have an analysis with two or more samples and they are made with different library construction (eg. exome and amplicons) then VarGenius cannot execute the same analysis with two pipelines. 
But there is an hack you can do!

You first set the analysis being an "exome" and run VarGenius tasks: -quality_check, -trimming, -alignment -refinement. In this way for all samples will be executed the normal pipeline for exomes. Once completed, you must use SQL commands and do the following for only the samples which are amplicons based:

 - table **readfiles** set mrdup=0: *update readfiles set mrdup=0 where sampleid=XX;*
 - table **samples** set precalread=0 and baserecal=0: *update samples set precalread=0,baserecal=0 where sampleid=XX;*
 - table **analyses** set seqtype=targeted :  **update  analyses set seqtype ='targeted' where analysisid=YY;*
 - execute variant calling step in VarGenius for only that samples: *perl /pico/work/TELET_UDP/VarGeniusBeta/vargenius.pl -c user_config.txt -ra YY -rs XX -vc*

The variant calling will be executed only for those samples and the genotype step will merge all variant calling files. Once the variant calling finishes:
 - set again the analysis to be for exome:**update  analyses set seqtype ='exome' where analysisid=YY;*
 - execute VarGenius final steps for the complete analysis:  *perl /pico/work/TELET_UDP/VarGeniusBeta/vargenius.pl -c user_config.txt -ra YY -vf -fout*

## RAM Memory requests

This part of the configuration file contains all such variables used to ask resources for each task.
Without going crazy with these parameters, just set the RAM memory and number of CPUs requests according to RAM and CPUs limits of the cluster nodes.
For example if the nodes of your cluster have a RAM memory of 80GB then you set all the XXX\_mem variables greater than 80GB to that limit. And if the number of CPUs in the nodes is 8, then you set all XXX\_ncpus variables greater than 8 to that limit.

Therefore, *VarGenius* can be programmed in such a way it can use HPC cluster resources in an efficient way.

```
qc_nodes =
qc_select = 1
qc_ncpus = 5
qc_threads = 5
qc_mem = 10GB
```

Go through these parameters for all the tasks and change the memory and cpu accordingly.
The same considerations are to be done for the javamem\_YYY variables. Please change the values greater than your limits accordingly.

```
javamem_RTC = Xmx4g
javamem_IR = Xmx4g
javamem_BR = Xmx8g
javamem_PR = Xmx4g
```

---------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/vargenius
