#!/bin/sh

#Installation script for sotware needed in VarGenius
##Go to the bin directory
bin_dir=$1
cd $bin_dir

##COMMENT THIS SECTION IF YOU DON'T HAVE SUDO PRIVILEGES
###########################COMMENT FROM HERE
# update & upgrade #
sudo apt-get update
sudo apt-get upgrade

echo "INSTALLING VARGENIUS PRE-REQUISITES"
#MySQL
echo "INSTALLING PostgreSQL"
sudo apt-get install postgresql
#PERL
echo "INSTALLING PERL"
sudo apt-get install -y perl
#BioPERL
echo "INSTALLING BIOPERL"
sudo apt-get install -y bioperl
#PERL Libraries
echo "INSTALLING Additional PERL libraries"
sudo apt-get install -y libexcel-writer-xlsx-perl
sudo apt-get install -y libparallel-forkmanager-perl
sudo apt-get install libVcf
sudo apt-get install libdbd-pg-perl
#R
echo "INSTALLING R"
sudo apt-get install -y r-base
#JAVA
echo "INSTALLING JAVA"
sudo apt-get install openjdk-8-jdk
###########################COMMENT TO HERE
#COMMENT THIS SECTION IF YOU DON'T HAVE SUDO PRIVILEGES


#TAR: use -xfj for tar.bz2 and -zxvf for tar.gz


bwa_link="https://github.com/lh3/bwa/archive/v0.7.17.tar.gz"
bwa_gz="v0.7.17.tar.gz"
bwa_f="bwa-0.7.17"

samtools_link="https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2"
samtools_gz="samtools-1.8.tar.bz2"
samtools_f="samtools-1.8"

fastqc_link="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip"
fastqc_zip="fastqc_v0.11.7.zip"
fastqc_f="FastQC"

trimgalore_link="https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip"
trimgalore_zip="0.4.5.zip"
trimgalore_f="TrimGalore-0.4.5"

trimmomatic_link="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip"
trimmomatic_zip="Trimmomatic-0.36.zip"
trimmomatice_f="Trimmomatic-0.36"
trimmmomatic_exec="trimmomatic-0.36.jar"

cutadapat_link="https://pypi.python.org/packages/68/73/2ae48245bbf6d84a24bdf29540ed01669f09c5d21c26258f2ce07e13c767/cutadapt-1.16.tar.gz#md5=222d062207e2a3a2a0760caa83bf14ff"
cutadapat_gz="cutadapt-1.16.tar.gz"
cutadapat_f="cutadapt-1.16"

bedtools_link="https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz"
bedtools_gz="bedtools-2.25.0.tar.gz"
bedtools_f="bedtools2"


picard_link="https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar"

gatk_ver="4.1.4.0"
gatk_link="https://github.com/broadinstitute/gatk/releases/download/${gatk_ver}/gatk-${gatk_ver}.zip"
gatk_zip="gatk-${gatk_ver}.zip"
gatk_f="gatk-${gatk_ver}/gatk"

vcfpm_link="https://github.com/vcftools/vcftools/archive/master.zip"
vcfpm_zip="master.zip"

freebayes_ver="v1.3.1"
freebayes_link="https://github.com/ekg/freebayes/releases/download/${freebayes_ver}/freebayes-${freebayes_ver}"
freebayes_f="freebayes-${freebayes_ver}"

xhmm_ver=""
xhmm_link="https://bitbucket.org/statgen/xhmm/get/cc14e528d909.zip"
xhmm_zip="statgen-xhmm-cc14e528d909.zip"
xhmm_f="statgen-xhmm-cc14e528d909"


#R packages
Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("ggplot2")'
Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("reshape2")'
Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("RPostgreSQL")'
Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("ExomeDepth")'

##Dowload Vcf library for PERL
echo "DOWNLOADING VCF LIBRARIES FOR PERL"
wget $vcfpm_link 
unzip $vcfpm_zip
rm $vcfpm_zip

#FASTQC
echo "DOWNLOADING FASTQC"
wget $fastqc_link
unzip $fastqc_zip
rm $fastqc_zip
chmod -R 755 $bin_dir/$fastqc_f

#TRIMMOMATIC
echo "DOWNLOADING TRIMMOMATIC"
wget $trimmomatic_link
unzip $trimmomatic_zip
rm $trimmomatic_zip

#BWA
echo "DOWNLOADING BWA"
wget $bwa_link 
tar -zxvf $bwa_gz
rm $bwa_gz
cd $bwa_f
./configure --prefix=$bin_dir
make
cd $bin_dir

#SAMTOOLS
echo "DOWNLOADING SAMTOOLS"
wget $samtools_link
tar -jxvf $samtools_gz
rm $samtools_gz
cd $samtools_f
./configure --prefix=$bin_dir
make
cd $bin_dir

#PICARD
echo "DOWNLOADING PICARD"
wget $picard_link

#GATK (has been copied into the bin folder)
echo "DOWNLOADING GATK"
wget $gatk_link
unzip $gatk_zip
rm $gatk_zip

#BEDTOOLS
echo "DOWNLOADING BEDTOOLS"
wget $bedtools_link 
tar -zxvf $bedtools_gz
cd $bedtools_f
./configure --prefix=$bin_dir
make
cd $bin_dir

#FREEBAYES
echo "DOWNLOADING FREEBAYES"
wget $freebayes_link

#XHMM
echo "DOWNLOADING XHMM"
wget $xhmm_link
unzip $xhmm_zip
cd $xhmm_f
make 
rm $xhmm_zip
cd $bin_dir

echo "##############################PROGRAM LINKS #################################"
echo "Please use the following links for the programs into the configuration file:"
eval "echo -n 'R_path = ';  which R"
eval "echo -n 'java_path = ';  which java"
echo fastqc_path = $bin_dir/$fastqc_f/fastqc
echo trimmomatic_path = $bin_dir/$trimmomatice_f/$trimmmomatic_exec
echo samtools_path = $bin_dir/$samtools_f/samtools
#which samtools
echo bwa_path = $bin_dir/$bwa_f/bwa
echo bedtools_path = $bin_dir/$bedtools_f/bin/bedtools
#which bedtools
echo picard_path = $bin_dir/picard.jar
echo gatk_path = $bin_dir/$gatk_f
echo freebayes_path = $bin_dir/$freebayes_f
echo "You will need to manually download Annovar and Spidex because they need a registration"
echo "GATK4 was automatically download but we still reccomend to manually download GATK3.8"
echo "##############################PROGRAM LINKS #################################"

exit 1
