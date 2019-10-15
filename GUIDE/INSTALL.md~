# VarGenius - Installation


## Install Dependencies  

Currently *VarGenius* can be used only within a HPC cluster infrastructure managed with PBS (Portable Batch System) and needs some software to be installed. 
We wrote a PERL script for automated installation of the tools and databases needed (GATK, BWA, and so on..) but not Annovar that you will need to download manually since it needs a registration. The script works only for Linux and needs that you have apt-get installed.

**If you do not have SUDO privileges** :please go into the file VarGenius/LIB/install_all.sh and comment the section where sudo privileges are needed and indicated with : "##COMMENT THIS SECTION IF YOU DON'T HAVE SUDO PRIVILEGES"

To install *VarGenius* you have to create a folder to save all the files and folders produced during the analyses and a folder for the dependencies. Then, you run the install.pl script from the *VarGenius* folder.

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ mkdir vargenius_analyses 
              frank@compaq2:~$ mkdir bin 
              frank@compaq2:~$ cd VarGenius 
              frank@compaq2:~$ perl install.pl 
              ---------------- TERMINAL ------------------
              


This will ask to you the full path to the working folder you just created and to a bin folder for the programs. Then all the needed tools will be installed. Please wait for the installer to finish.

              ---------------- TERMINAL ------------------ 
		This script will prepare vargenius.pl to work in your directory.
		Are you sure you want to enjoy this software?(y or n) y
		Write an existing complete path where you want to play with vargenius.pl (/home/username/analyses): /home/francesco/vargenius_analyses/
		Write an existing complete path where you want to install the programs (/home/username/bin/): /home/francesco/bin/
		Installing dependencies with  /home/francesco/VarGenius/LIB/install_all.sh

		[.....DEPENDENCIES ARE BEING INSTALLED...]

		##############################PROGRAM LINKS #################################
		Please use the following links for the programs into the configuration file:
		R_path = /home/francesco/R
		java_path = /home/francesco/java
		fastqc_path = /home/francesco/bin/FastQC/fastqc
		trimmomatic_path = /home/francesco/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
		samtools_path = /home/francesco/bin/samtools-1.8/samtools
		bwapath = /home/francesco/bin/bwa-0.7.17/bwa
		bedtools_path = /home/francesco/bin/bedtools2/bin/bedtools
		picard_path = /home/francesco/bin/picard.jar
		gatk_path = /home/francesco/bin/
		You will need to manually download Annovar and Spidex because they need a registration 
		GATK4 was automatically download but we still reccomend to manually download GATK3.8
		##############################PROGRAM LINKS #################################

		[...]

		Now you can work with VarGenius from this folder!
		Only the first time that you run VarGenius you must: 
		1. Create the tables of the database 
		2. Download the genes information 
		3. Download the genome fasta and the datasets for GATK

              ---------------- TERMINAL ------------------ 

As you can see *VarGenius* installer gives you the paths where the programs are installed. Please keep this paths at hand to insert them later during the Tutorial.

**Download GATK Annovar and Spidex**

*VarGenius* can run either GATK4 and GATK3. GATK4 is now an open source software and can be automatically downloaded by VarGenius.
Yet, we suggest to users to download the previous version (GATK3.8) that must be manually downloaded because previously GATK was not an open source.
Version 3.8 is still preferable for *VarGenius* users because there are still few suggestions on the Broad website on how to use GATK4 even though it is no more in beta version. But you can decide to use GATK4 at your own risk.

Annovar and Spidex need user registration hence they cannot be automatically downloaded by *VarGenius*. Yet, the Annovar databases will be automatically downloaded and installed by *VarGenius*!

Download Annovar from the website: http://annovar.openbioinformatics.org/en/latest/user-guide/download/
and Spidex [HERE](http://www.openbioinformatics.org/annovar/spidex_download_form.php). You need to fill the form and wait for the response from the administrator. 

Once you finished with Annovar, see how to use the first time *VarGenius*  in the [TUTORIAL](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/TUTORIAL.md)

If you could not install automatically the dependencies the following are some instructions to do it manually.







## Manual installation of dependencies

If you are not confident with Parallel Architecture management or you do not have sudo privileges you should ask her/him for the installation of the following software and get their complete paths. Otherwise please refer to this manual installation guide.

 
  - PostgreSQL database (https://www.postgresql.org/)
  - Perl: (tested with ver5.10) (http://www.perl.org/get.html) ;
  - BioPerl: (tested with ver1.6) (http://www.bioperl.org/wiki/Getting_BioPerl); 
  - R: (tested with ver3.2.3) (http://www.r-project.org/). 
  - Java: (tested with ver.1.8.0) 

The following are the supplementary Perl and R modules. If they are not included in the default installation, install them manually.

  - Supplementary Perl MODULES:
    - Parallel::ForkManager.pm (http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm)
    - Excel::Writer::XLSX (http://search.cpan.org/~jmcnamara/Excel-Writer-XLSX/lib/Excel/Writer/XLSX.pm)
    - Vcf.pm (http://search.cpan.org/~ajpage/Bio-Pipeline-Comparison-1.123050/lib/Vcf.pm) 

   You may want to use the following Linux commands to check if you already have these perl modules:

	- perl -MParallel::ForkManager -e 'print "Module installed.\n";'
	- perl -MExcel::Writer::XLSX -e 'print "Module installed.\n";'
	- perl -MVcf -e 'print "Module installed.\n";'

    To install them on a Debian derivate linux distribution use from a terminal: 
    - sudo apt-get install libparallel-forkmanager-perl 
    - sudo apt-get install libexcel-writer-xlsx-perl
    - sudo apt-get install libVcf
  *HINT: if you do not have system administrator permissions to install into the cluster you can download manually these packages. Put the needed .pm files into a directory that you choose (e.g. /home/perl_lib/) and add this path into the perl_lib_paths.txt file located into the VarGenius/CONFIGURATION folder. Later the install script will add the needed libraries into VarGenius.*

  - R needed libraries:
      - RPostgreSQL
      - ggplot2
      - reshape2
  Start an R session and use library("package") to check if the libraries are installed.

Then you have to download the programs that *VarGenius* uses: BWA, FastQC, Trimmomatic, TrimGalore, Picard, Samtools, GATK, bedtools, Annovar. The tested versions are in bracket. Different versions may not work with this *VarGenius* version.
Actually, *VarGenius* uses the full path to each of the software. Hence, it does not require installation of them. As described into the TUTORIAL page, you have to write their full paths into the configuration file.

Download them at: 

- BWA: http://bio-bwa.sourceforge.net/ (ver. 0.7.17)
- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/download.html (ver. 0.11.7)
- Trimmomatic: ttp://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip (ver. 0.36)
- Picard tools: https://broadinstitute.github.io/picard/ (ver. 2.18.2)
- Samtools: http://www.htslib.org/  (ver. 1.8)
- GATK: https://software.broadinstitute.org/gatk/download/ (gatk. 3.8)
- bedtools: http://bedtools.readthedocs.io/en/latest/content/installation.html (ver. 2.25)
- Annovar: http://annovar.openbioinformatics.org/en/latest/user-guide/download/ ver. 2017Jul)

*HINT: you may want to download the folders and put them all into a **bin** directory to avoid to have a puzzling situation with paths.*





------------------------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/VarGenius
