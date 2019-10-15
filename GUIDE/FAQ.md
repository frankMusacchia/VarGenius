# VarGenius - Frequently Asked Questions


Here you can find issues already solved at any step with *VarGenius*. Please search here before to post any question on the forum.
   

*Question*: *VarGenius* gave the following error:
```
Exception in thread "main" java.lang.UnsupportedClassVersionError: org/broadinstitute/gatk/engine/CommandLineGATK : Unsupported major.minor version 52.0
```

*Answer*: This is an error that happens with GATK when using an uncompatible java version. Please search the GATK website the best java version to use for GATK.


-------------------------

*Question*: *VarGenius* gave the following error:
```
Error: the required database file hg19_cadd13.txt does not exist
```

*Answer*: This is an error occurring when Annovar cannot download a database. You can try to re-run *VarGenius* using **annov_db_dl = YES** or try to install manually the database using **annotate_variation.pl -downdb  -webfrom annovar -buildver hg19 cadd13  /..pathto../annovar/humandb**

-------------------------


------------------------------------------------

*Question*: *VarGenius* gave the following error durin ExomeDepth execution:
	
```
It looks like the test samples has only 0 bins with more than 5 reads. The coverage is too small to perform any meaningful inference so no likelihood will be computed.
```

*Answer*: This is an error occurring when ExomeDepth could not analyze your BAM file. We could not find a solution to this problem hence we managed to write a bit of code to exclude samples giving errors from ExomeDepth computation. You can create the files exome_depth_bam_blacklist_sex.txt and exome_depth_bam_blacklist_auto.txt depending by which of the two runs gave you the error and simply add the analysis name into that file.

The exome_depth_bam_blacklist must be located into vargenius_analyses/DATA folder
-------------------------------------------------


If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/VarGenius
