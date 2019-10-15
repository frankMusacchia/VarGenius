# *VarGenius*
--------------------------------

## License
VarGenius - Variant Discovery and Annotation Tool
Copyright (C) <2018>  <Francesco Musacchia>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Mission
*VarGenius* is a pipeline for the detection and annotation of genetic variants. It can be currently used with WES and targeted sequencing.

## Installation and tutorial

To start using *VarGenius* you must follow the instructions into the **GUIDE** folder. There you will find additional information on the content of the User Guide.

- The **install and go** section is here: [INSTALL](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/INSTALL.md). There you will also learn how to install *VarGenius*;
- After you install you can try to run the pipeline for the first time. Few settings are needed before to run. Please follow instructions from the [TUTORIAL](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/TUTORIAL.md).


Please enjoy this product and send us questions and suggestions!

## System Requirements

To work, *VarGenius* needs:
- PostgreSQL database
- Perl
- BioPerl
- R 
- samtools
- picard
- GATK
- BWA
- TrimGalore
- FastQC
- bedtools
- Annovar

Check for the tested versions in the [INSTALL](https://github.com/frankMusacchia/VarGenius/blob/master/GUIDE/INSTALL.md) page.


## *VarGenius* times and space

*VarGenius* can be used only with a PBS cluster. It has been desigend to create jobs with the **qsub** command. The execution time depends strongly from the speed and memory of your machine. 

The pipeline needs around 5 hours to run the complete sequence of tasks. Depending by the queue and resources used by other users of the cluster, jobs will wait less or more time.


#### Space

The space needed to build the PosgreSQL database depends from the number of variants contained.  100 samples (WES and targeted) consisting of 300 thousands variants use around 4GB of disk space. 

The space needed for all the default databases of Annovar takes 1.2 Terabyte in our system

------------------------------------------------

If you get some error during the installation or the running of *VarGenius* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/VarGenius
