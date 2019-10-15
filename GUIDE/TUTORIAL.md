# VarGenius - Configure and run!

## Before to start

In this tutorial we will call *task* all the tasks that are executed by *VarGenius* (quality check, alignment, refinement, etc) and *steps* all the specific steps contained into the tasks (quality check: quality check, trimming, quality check after trimming; etc).

To use *VarGenius* you should have:
 - the fastq files of your sequenced samples;
 - the target .BED file included into the enrichment kit;
 - a configuration file;
 - a sample sheet with all information about the samples;

Yet, we provide sample data for this tutorial. Please follow carefully this page of tutorial. It looks long but after doing this steps the first time you will always use almost the same configuration and thank God that *VarGenius* exists!

## Configuration File

Before you run *VarGenius* you must set some parameters of the user_config.txt file related with your system and configuration. Without setting all of them *VarGenius* will not work.

The user_config.txt file has been copied during the installation in your working folder. Please open it with a text editor (nano, vi, etc..) and change the following parameters accordingly with your installation and configuration.

		---------------- TERMINAL ------------------ 
		frank@compaq2:~$ cd vargenius_analyses 
		frank@compaq2:~$ nano user_config.txt
		---------------- TERMINAL ------------------ 

### Group 1: Pipeline Steps Execution

The first set of parameters concerns the pipeline steps to execute. Default steps to use are present. Do not modify them now.
Once you will be expert with *VarGenius* you could use these parameters to include/exclude steps for the analyses.


### Group 2: NEEDED TO RUN VARGENIUS

**Programs Paths**

Full paths of the installed programs: folders or executable paths. Please set all of them accordingly whit what the installer gave you before or with your manual installation. 
N.B. We reccomend the manual download of GATK3.8!!

```
Needs path to the executable:
R_path = /.../bin/R
java_path = /.../java
fastqc_path = /.../fastqc
trimmomatic_path = /.../trimmomatic/0.33/binary/trimmomatic-0.33.jar
samtools_path = /.../samtools
bwa_path = /.../bwa_folder/bin/bwa
picard_path = /.../picard.jar
gatk_path =  /.../GenomeAnalysisTK.jar

Needs path to the folder
bedtools_path = /.../bedtools/bin/
annovar_path = /.../annovar_072017/
```

... if still your tick head is telling you to use GATK4 please modify the last parameter of group2 config variables to
gatk_ver = 4

*Hint: you may want to use the linux **which** command to get the full paths of the programs.* 

### Database Configuration
**Instruction on how to create an user account with PostgreSQL**

You need system administrator or sudo permissions to create a database and a PostgreSQL account. 

First you need to create the database: 

              ---------------- TERMINAL ------------------ 
              > CREATE DATABASE vargenius_db;
              ---------------- TERMINAL ------------------

Then, you add a new user:

              ---------------- TERMINAL ------------------ 
              > CREATE USER vargenius
              ---------------- TERMINAL ------------------ 
              
Then you have to grant the user the privileges to access and modify the database:

              ---------------- TERMINAL ------------------ 
              > GRANT ALL PRIVILEGES ON DATABASE vargenius_db to vargenius;
              ---------------- TERMINAL ------------------ 


If you cannot either install a PostgreSQL database or create an user with all privileges, please ask to your sys admin. Once that the sysadmin created it and gave you an user account, set the following variables. Change just the host and the port of db_dsn with yours.

```
DATABASE configuration
db_name = vargeniusdb
db_dsn = dbi:Pg:host=domain.yournetwork.com;port=5432;sslmode=allow;
db_host = domain.yournetwork.com
db_user = vargenius_user
db_pass = password
```
Notice that if you have root privileges you will be able to create the database automatically later with this settings. Otherwise, you will just create the tables for the database.
 
### QSUB Jobs Resources request

To execute automatically jobs with the cluster you need to set the appropriate parameters for qsub. Please ask your system administrator what are the parameters that you can use for account, queue, walltime, your username and maximum memory and cpu:

```
qsub_account = your_account_name
qsub_queue = queue
qsub_walltime = 5:00:00
qsub_username = frank
qsub_walltime_max = 24:00:00
qsub_mem = 120GB
qsub_ncpus = 12
```


### Data to send email

At the end of the analysis an email will be sent using a linux command. Set these parameters as preferred.
The mandatory one is the email_recipients which contains comma separated email addresses of users interested in the analysis.

```
email_author = VarGenius@tigem.it
email_recipients = youraddress@fakedomain.com
email_subject = VarGenius analysis is completed
email_message =  analysis is complete! Please check results in WebServer Path. This message is automatically generated. Please do not use the reply button. For any question send an email to youraddress@fakedomain.com.
```

## Group3: OPTIONAL

These parameters can be changed optionally. 
Please leave them as they are and read the [ADVANCED_USAGE](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/ADVANCED_USAGE.md) page if you want to modify them later.


## Group4: SOFTWARE PARAMETERS


There are many parameters of the software used in *VarGenius* that you may want to change in this section. Their names are the same as in the software you are using (BWA, GATK, etc.) but for this tutorial you do not change any of them. Just use default parameters. 

### RAM Memory requests

There is a part of the configuration file containing variables used to ask resources for each task. For this tutorial you do not change any of them. Just use default parameters.

If you have any problem with memory usage please read the [ADVANCED_USAGE](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/ADVANCED_USAGE.md) page.

And that's all for the configuration file! Let's fill the database!


# FILL THE DATABASE


## Create the database

The database must be created or filled the first time that you run *VarGenius* otherwise it will not work at all.

If an administrator made an user account for you and gave full privileges you can use the following command and automatically create the database and tables:


              ---------------- TERMINAL ------------------ 
		perl ../VarGenius/vargenius.pl -c user_config.txt --create_database
              ---------------- TERMINAL ------------------ 


If you created manually the database **vargenius_db**, now just run:

              ---------------- TERMINAL ------------------ 
		perl ../VarGenius/vargenius.pl -c user_config.txt --create_tables
              ---------------- TERMINAL ------------------ 


## Download references

Now you need to download the reference genome and index it for BWA and GATK. And download the databases to use in GATK.

Please run the command:

              ---------------- TERMINAL ------------------ 
		perl ../VarGenius/vargenius.pl -c user_config.txt --dldatasets
              ---------------- TERMINAL ------------------ 


A job will be launched to perform this process that will need some time. You can follow the job looking at the log file created: vargenius_analyses/LOG/dldatasets_DAY_MONTH_TIME_YEAR.log

In parallel, you can launch the job to fill the database with genes information.

## Genes information

*VarGenius* database has tables containing information about genes. Specifically: gene name, HPO identifiers associated, GDI score, RVIS score, RefSeq, OMIM and Entrez identifiers of the gene.
These information is automatically downloaded by *VarGenius* but you need to start the process to create it before you use this software.

Just run:

              ---------------- TERMINAL ------------------ 
		perl ../VarGenius/vargenius.pl -c user_config.txt --genedb
              ---------------- TERMINAL ------------------ 

This process will need some time. 

In the meanwhile you can go on and modify the example sample sheet for the first run...



# RUN VARGENIUS


### Spidex database

Before to run the first analysis please move the Spidex database that you downloaded previously into the folder:

```
/path/to/vargenius_analyses/DATA/references/annovar\_humandb\_2018/
```


## Sample sheet creation

The sample sheet collects all the information about the localization of the fastq files, samples gender and kinship info, sequencing type (WES or targeted) and other parameters. For this tutorial you can use the example files from the VarGenius/EXAMPLE_DATA folder:
 - target.bed
 - SampleA_R1_001.fastq.gz
 - SampleA_R2_001.fastq.gz
 - sample\_sheet\_test.tsv

We suggest to create a folder for your samples and samples sheets (outside of vargenius_analyses folder). Let's say: /home/users/frank/samples/ and /home/users/frank/sample\_sheets/
Then, copy there the EXAMPLE_DATA/SampleA

```
cp -r VarGenius/EXAMPLE_DATA/SampleA/ samples/
```

Put the target.bed file into your vargenius_analyses/DATA/references/targets folder.

```
cp -r VarGenius/EXAMPLE_DATA/target.bed vargenius_analyses/DATA/references/targets/
```

Copy the file   VarGenius/EXAMPLE_DATA/sample\_sheet\_test.tsv into the sample\_sheets folder. Then open it with an editor and modify the paths of the fastq files (fqdir column) to the one you have created above: /home/users/musacchia/samples/SampleA/

```
cp -r VarGenius/EXAMPLE_DATA/sample_sheet_test.tsv /home/users/frank/sample\_sheets/
nano /home/users/frank/sample\_sheets/sample_sheet_test.tsv
```


If you want to try our in-house script for the creation of the sample sheet please read the [ADVANCED_USAGE](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/ADVANCED_USAGE.md) section.



## Finally ...Run an analysis

Please wait that the **genedb** and **dldatasets** processes finish. You can look at the current job run status using the qstat command.

After you finished this configuration, at each analysis you will just need to launch two commands to run *VarGenius*: the sample sheet creation and the execution.

The perl script contains an helpful command to look specifically at all the functionalities. Please read it later because you will find that there are a set of useful operation that will make your life simpler.

To launch your first analysis you go in the vargenius\_analyses folder and type the command:


```
perl /home/frank/VarGenius/vargenius.pl -c user\_config.txt -ss /home/users/frank/sample\_sheets/sample\_sheet\_test.tsv --start\_all
```

This simple command will execute all the tasks we implemented in *VarGenius* with the default steps depending by the sequencing type (exome or targeted) that you selected into the sample sheet. Different analyses can be run using the same sample sheet.

There is not much more to say after this to run this program. If all stuffs work, enjoy *VarGenius*!!


### (HINTS) IF YOU HAVE PROBLEMS WITH THE FIRST RUN:
Once that the sample_sheet is loaded into the db, if you want to rerun the analysis, you must first get the db id of the analysis:

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt -aid SampleA 
10
```

and then use that id to run the analysis again:

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt -ra 10 --start_all
```

if you want to remove a sample sheet and its analyses please use 

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt -remss sample\_sheet\_test.tsv
```

if you want to remove an analysis use:

```
perl /pico/work/TELET_UDP/VarGenius/vargenius.pl -c user_config.txt  -reman 10
```

If you want to execute only some step from the pipeline please look at the *VarGenius pipelineS* section into the [ADVANCED_USAGE](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/ADVANCED_USAGE.md) page 
If you are sure that a task has completed all its steps for all the samples you can successfully run all the remaining tasks of the cascade. 



Now to understand what is the output and how to interpret it go to the [OUTPUT](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/OUTPUT.md) page.


---------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/vargenius

