
/*Creating a type enumerative for the type of analysis
P_  proband
A_ affected non proband (brother, syster or cousin of proband)
M_ 	mother
F_  father
R_  any relative inserted for any reason
*/
CREATE TYPE kinship_type AS ENUM ('P','F','M','A','R','-');

/*Creating a type enumerative for the type of kinship
The first level of the family has a single letter
P_  proband
B_  brother (B1,B2...)
Z_ 	sister (Z1,Z2,...)
F_  father
M_  mother
S_ son of the proband
D_ daughter of the proband
W_ wife

From the second level an additional letter indicates different relationships
always related with the proband
MZ_ sister of the mother
MB_ brother of the mother
MBS_ son of the brother of the mother (the male cousin)
FM_ mother of the father
MF_ father of the mother
FFB_ brother of the father of the father (old uncle, brother of the grandfather)
FSD_ daughter of a sister of the father (female cousin)
FSD2_ second daughter of a sister of the father
*/

/*Creating a type enumerative for the gender*/
CREATE TYPE gender_type AS ENUM ('M','F');

/*Creating a type enumerative for the inheritance models*/
CREATE TYPE inheritance_type AS ENUM ('DN','AD','HC','HO','XL','Other');

/*Creating a type enumerative for the type of validation with Sanger
	0: not considered for the validation
	1: validated and resulted positive
	2: validated and resulted negative
*/
CREATE TYPE validation_type AS ENUM ('0','1','2');

/*A table for information about the program and the execution*/
CREATE TABLE info (
	infoid serial PRIMARY KEY,
  inforow int DEFAULT 1,
  progver VARCHAR(50),
	dbversion VARCHAR(50),
	dbbusy int DEFAULT 0,
  jobin VARCHAR(100),
  analysisid int not null,
  ssid int not null
);

/*A table for target files information*/
CREATE TABLE target_files (
	targetid serial PRIMARY KEY,
	targetname  VARCHAR(100) not null,
	perchrdiv int DEFAULT 0,
	perchrdivdir VARCHAR(1000) DEFAULT 'none',
	xhmmGCcontentFile int DEFAULT 0,
	bedtointerval  int DEFAULT 0
);

/*A table for profiles information (users)*/
CREATE TABLE profiles (
	profid int not null,
	proffile  VARCHAR(1000),
	profdesc  VARCHAR(1000),
	proffolder  VARCHAR(1000)
);

/*A table for Annovar database installed*/
CREATE TABLE annovar_dbs (
	db_id serial PRIMARY KEY,
	ann_db_name VARCHAR(1000),
	ann_db_ver VARCHAR(100),
	ann_db_type VARCHAR(10)
);

/*A table for sample sheets paths*/
CREATE TABLE sample_sheets (
	ssid serial PRIMARY KEY,
	ss_name VARCHAR(1000),
	unique (ss_name)
);

/*Table with analyses information*/
CREATE TABLE analyses (
	analysisid serial PRIMARY KEY,
	userid int DEFAULT -1,
	resgroupid int DEFAULT -1,
  analysisname VARCHAR(50),
  ssid int not null,
  toremove  int DEFAULT 0,/*If the analysis must be removed*/
  stored  int DEFAULT 0,/*If the analysis has been stored into the storage area*/
  pedfile VARCHAR(150) DEFAULT 'none',
  targetbed VARCHAR(150) DEFAULT '-',/*Target file name*/
  perchrom int DEFAULT 0,/*is 1 if analysis has been ran per chromosome*/
  convscores int DEFAULT 0,/*is 1 if the scores must be converted*/
  dogenot int DEFAULT 0,/*is 1 if the analysis is composed of more samples to run the joint genotype*/
  seqtype VARCHAR(150) not null,
  infreq int DEFAULT 0,/*If the analysis must be kept in account for the frequency calculation*/
  consensus int DEFAULT 0,/*If the analysis was run with the consensus approach*/
  targetextended int DEFAULT 0,/*If the analysis was run extending the targetbed file*/
  cnvanalysis int DEFAULT 0,/*If the analysis of CNV was performed*/
  varcall int DEFAULT 0,
  catvar int DEFAULT 0,
  genot int DEFAULT 0,
	varfilt  int DEFAULT 0,
	varrecal int DEFAULT 0,
	apprecal int DEFAULT 0,	
	genotref int DEFAULT 0,
 	phasing int DEFAULT 0,
 	finalout int DEFAULT 0,
 	annotate int DEFAULT 0,
 	vcfimport int DEFAULT 0,
 	missimport int DEFAULT 0,
 	status int DEFAULT 0, /*Five different status are possible*/
 	analdate VARCHAR(150) DEFAULT 'none',
 	analinfo VARCHAR(500) DEFAULT 'none',
 	scriptver VARCHAR(1000) DEFAULT '-',/*Which version of VarGenius was used*/
  FOREIGN KEY(ssid) REFERENCES sample_sheets(ssid)
);

/*Table with samples information*/
CREATE TABLE samples (
	sampleid serial PRIMARY KEY,
	analysisid int not null,
	personid int DEFAULT -1,
	samplename VARCHAR(50),
	gender gender_type,
	kinship kinship_type,
	target_cov float DEFAULT -1,	
	ssid int not null,
	affected int DEFAULT 0,
	multiple int DEFAULT 0,
	sortidx int DEFAULT 0,
	indreal int DEFAULT 0,
	baserecal int DEFAULT 0, 
	precalread int DEFAULT 0,
	mergegrp int DEFAULT 0,
	mrdupgrp int DEFAULT 0,
	varcall int DEFAULT 0,
	catvar int DEFAULT 0,
	mergebam int DEFAULT 0,
	mergealn int DEFAULT 0,
	coverbed int DEFAULT 0,
	readfbamused int DEFAULT 0,/*This is a flag that can be incremented each time a process finishes to use the bam files of the reads*/
  FOREIGN KEY(ssid) REFERENCES sample_sheets(ssid),
  FOREIGN KEY(analysisid) REFERENCES analyses(analysisid)  
);

/*Table with files information*/
CREATE TABLE readfiles (
	readfid serial PRIMARY KEY,
	readfname VARCHAR(100), 
	sampleid int not null,
	analysisid int not null,
	ssid int not null,
	fqdir VARCHAR(1000), 
	fq1 VARCHAR (200), 
	fq2 VARCHAR(200),
	instrument VARCHAR(200),
	flowcell VARCHAR(200),
	lane VARCHAR(200),
	goodquality int DEFAULT 0,
	qc int DEFAULT 0,
	trim int DEFAULT 0,
	mergepe int DEFAULT 0,
	align int DEFAULT 0,
	mrdup int DEFAULT 0,
	splnc int DEFAULT 0,
	alt int DEFAULT 0,
	samview int DEFAULT 0,
	samsort int DEFAULT 0,
	sortidx int DEFAULT 0,
	realtar int DEFAULT 0,
	indreal int DEFAULT 0,
	baserecal int DEFAULT 0,
	precalread int DEFAULT 0,
	FOREIGN KEY(ssid) REFERENCES sample_sheets(ssid),
	FOREIGN KEY(analysisid) REFERENCES analyses(analysisid),
	FOREIGN KEY(sampleid) REFERENCES samples(sampleid)
);

/*Tables for variants*/

/*A table for Variants*/
CREATE TABLE variants (
	varid serial PRIMARY KEY,
	chrom VARCHAR(50),
	pos VARCHAR(50),
	id VARCHAR(10) DEFAULT '-',
	ref VARCHAR(1000),
	alt VARCHAR(1000),
	compid VARCHAR(1000),
	targeted_allele_freq decimal,
	targeted_freq_factors VARCHAR(1000),
	exome_allele_freq decimal,
	exome_freq_factors VARCHAR(1000),
	wgs_allele_freq decimal,
	wgs_freq_factors VARCHAR(1000),	
	rna_allele_freq decimal,
	rna_freq_factors VARCHAR(1000)
);

/*A table for Variant frequencies update history*/
CREATE TABLE var_freq_update_hist (
	insertid serial PRIMARY KEY,
	seqtype VARCHAR(150) not null,
	varfrequpdate VARCHAR(150) DEFAULT 'none',
	insertinfo VARCHAR(500) DEFAULT 'none'
);

/*A table for Statistics about Variants*/
CREATE TABLE statistics (
	varstid serial PRIMARY KEY,
	varid int not null,
	analysisid int not null,
	qual float  DEFAULT -1,
	filter VARCHAR(10000),
	DP float  DEFAULT -1,
	DB boolean,
	AC float  DEFAULT -1,
	AF float  DEFAULT -1, 
	AN float  DEFAULT -1,
	ExcessHet float  DEFAULT -1,
	ClippingRankSum float  DEFAULT -1,
	SOR float  DEFAULT -1,
	FS float  DEFAULT -1, 
	MLEAC float DEFAULT -1,
	MLEAF float DEFAULT -1,
	MQ float DEFAULT -1,
	MQ0 float DEFAULT -1,
	QD float DEFAULT -1,
	BaseQRankSum float DEFAULT -1,
	MQRankSum float DEFAULT -1,
	ReadPosRankSum float DEFAULT -1,
	set VARCHAR(1000),
	FOREIGN KEY(analysisid) REFERENCES analyses(analysisid),
	FOREIGN KEY(varid) REFERENCES variants(varid)
);

/*A table for Genotypes associated to a given sample*/
CREATE TABLE genotype_sample (
	genid serial PRIMARY KEY,
	analysisid int not null,
	sampleid int not null,
	varids VARCHAR(1000),
  GT VARCHAR(10) DEFAULT '-',
  AD VARCHAR(1000) DEFAULT '-',
  DP VARCHAR(1000) DEFAULT '-',
  GQ VARCHAR(1000) DEFAULT '-',
  PL VARCHAR(1000) DEFAULT '-',
  TP VARCHAR(1000) DEFAULT '-',
  PGT VARCHAR(1000) DEFAULT '-',
  PID VARCHAR(1000) DEFAULT '-',
	FOREIGN KEY(analysisid) REFERENCES analyses(analysisid),
	FOREIGN KEY(sampleid) REFERENCES samples(sampleid)
);

/*A table for Annotations associated to a given variant*/
CREATE TABLE annotations (
	/*annid serial PRIMARY KEY,*/
	varid int not null,
	gene_refgene VARCHAR(10000) DEFAULT '-',
	func_refgene VARCHAR(10000) DEFAULT '-',
	refgene_transcript VARCHAR(10000) DEFAULT '-',/*Transcript from RefSeq*/
	refgene_exon VARCHAR(10000) DEFAULT '-',/*Exon*/
	refgene_nt_change  VARCHAR(10000) DEFAULT '-',/*nucleotidic change*/
	refgene_aa_change VARCHAR(10000) DEFAULT '-',/*Brings all the aminoacids changes*/
	gene_gencode VARCHAR(10000) DEFAULT '-',	
	func_gencode VARCHAR(10000) DEFAULT '-',
	gencode_transcript VARCHAR(10000) DEFAULT '-',/*Transcript from RefSeq*/
	gencode_exon VARCHAR(10000) DEFAULT '-',/*Exon*/
	gencode_nt_change  VARCHAR(10000) DEFAULT '-',/*nucleotidic change*/
	gencode_aa_change VARCHAR(10000) DEFAULT '-',/*Brings all the aminoacids changes*/
	exac_all VARCHAR(100) DEFAULT '-',
	gnomad_exome_all VARCHAR(100) DEFAULT '-',
	polyphen2_hdiv_pred VARCHAR(100) DEFAULT '-',
	polyphen2_hvar_pred VARCHAR(100) DEFAULT '-',
	mutationtaster_pred VARCHAR(100) DEFAULT '-',
	mutationassessor_pred VARCHAR(100) DEFAULT '-',
	lrt_pred VARCHAR(100) DEFAULT '-',
	provean_pred VARCHAR(100) DEFAULT '-',
	fathmm_pred VARCHAR(100) DEFAULT '-',
	sift_pred VARCHAR(100) DEFAULT '-',
	fathmm_mkl_coding_pred VARCHAR(100) DEFAULT '-',
	metasvm_pred VARCHAR(100) DEFAULT '-',
	metalr_pred VARCHAR(100) DEFAULT '-',
	FOREIGN KEY(varid) REFERENCES variants(varid)
);

/*A table for genes information (scores,names etc)*/
CREATE TABLE genes (
	geneid serial PRIMARY KEY,
	genename  VARCHAR(1000) ,
	hpoids  VARCHAR(100000),
	gdi float,
	rvisscore float,
	rvisperc float,
	refseqids  VARCHAR(100000),
	entrezid VARCHAR(100),
	omimids  VARCHAR(100000),
	omimgenedesc VARCHAR(100000),
	omimphenotype VARCHAR(100000)
);

/*A table for phenotypes associated to a given gene*/
CREATE TABLE phenotypes (
	phenid serial PRIMARY KEY,
	hpoid  VARCHAR(100) not null,
	hponame VARCHAR(10000) not null,
	hpodef VARCHAR(100000)
);

/*A table for transcripts associated to a given gene*/
CREATE TABLE transcripts (
	trid serial PRIMARY KEY,
	refseqid  VARCHAR(100) not null,
	refseqdesc VARCHAR(10000),
	geneids  VARCHAR(100000)
);


/*A table for information associated to a given sample*/
CREATE TABLE samples_info (
	personid serial PRIMARY KEY,
	samplename VARCHAR(50),
	familyid VARCHAR(1000) DEFAULT '-',
	dob VARCHAR(1000) DEFAULT '-',
	pob VARCHAR(10000) DEFAULT '-',
	gender VARCHAR(10) DEFAULT '-',
	kinship VARCHAR(10) DEFAULT '-',
	affected  int DEFAULT 0,
	diseaseid int DEFAULT 0,
	notes VARCHAR(100000) DEFAULT '-'
);


/*Table with readfiles statistics*/
CREATE TABLE readfile_statistics (
	readfid int not null,
	analysisid int not null,
	avgreadsqual float DEFAULT -1
);

/*Table with sample statistics*/
CREATE TABLE sample_statistics (
	sampleid int not null,
	analysisid int not null,
	avgreadsqual float DEFAULT -1,
	totreads float DEFAULT -1,
	mappedreads float DEFAULT -1,
	avgcov float DEFAULT -1,
	duplicates float DEFAULT -1
);

/*Table with analysis statistics
	N.B. This field names MUST be identical to the fields
	inserted in the header into script var_stats.R
	section: #STATISTICS PER-ANALYSIS
*/
CREATE TABLE analysis_statistics (
	analysisid int not null,
	totvariants int DEFAULT 0,
	totsnvs int DEFAULT 0,
	snvsperc float DEFAULT -1,
	totindels int DEFAULT 0,
	indelsperc float DEFAULT -1,
	meanqual float DEFAULT -1,
	denovo int DEFAULT 0,
	denovoperc float DEFAULT -1,	
	synsnv int DEFAULT 0,
	synsnvperc float DEFAULT -1,	
	nonsynsnv int DEFAULT 0,
	nonsynsnvperc float DEFAULT -1,
	ns_sratio float DEFAULT -1,
	snvtr int DEFAULT 0,
	snvtrperc float DEFAULT -1,
	snvtv int DEFAULT 0,
	snvtvperc float DEFAULT -1,		
	tr_tvratio float DEFAULT -1,
	meanexacmax float DEFAULT -1,		
	meanexacall float DEFAULT -1,			
	dbsnpe int DEFAULT 0,
	dbsnpeperc float DEFAULT -1,
	dbsnpne int DEFAULT 0,
	dbsnpneperc float DEFAULT -1,
	inhgmd int DEFAULT 0,
	inhgmdperc float DEFAULT -1,
	caddgt20 float DEFAULT -1,
	caddgt20perc float DEFAULT -1,
	splvar int DEFAULT 0,
	splvarperc float DEFAULT -1,
	mendelviol int DEFAULT 0
);


/*A table for genecoverage for each gene and for each sampleid*/
CREATE TABLE genecoverage (
	sampleid int not null,
	geneid int not null,
	genecoverage float DEFAULT -1
);

/*A table for HPO associated to a given sample*/
CREATE TABLE samples_hpo (
	personid int not null,
	phenid int not null,
	observed boolean,
	FOREIGN KEY(personid) REFERENCES samples_info(personid),
	FOREIGN KEY(phenid) REFERENCES phenotypes(phenid)
);

/*A table for candidate variants and the transmission*/
/*Here analysisid and varid are correspondingly correlated with the analysis
and the variant but the data will not be removed if an analysis is simply repeated. 
In that case the user should manually update the analysisid with the newer*/
CREATE TABLE candidate_variants_2 (
	candvarid serial PRIMARY KEY,
	varid int not null,
	compid VARCHAR(1000),	
	analysisid int not null,
	analysisname VARCHAR(50),
	avsnp VARCHAR(1000) DEFAULT '-',
	qual float  DEFAULT -1,
	uncertain_samples int,
	freq_max float DEFAULT -1,
	ivf float DEFAULT -1,
	freq_factors VARCHAR(1000) DEFAULT '-',
	exac_ac_het int,
	exac_ac_hom int,
	exac_ac_hemi int,
	genename VARCHAR(1000) DEFAULT '-',
	refgene_aa_change VARCHAR(10000) DEFAULT '-',
	known_variant VARCHAR(1000) DEFAULT '-',
	genefunctionclass VARCHAR(100) DEFAULT '-',
	selection_criteria VARCHAR(10000) DEFAULT '-',
	model VARCHAR(1000) DEFAULT '-',
	interpretation int DEFAULT -1,
	exac_pli float DEFAULT -1,
	exac_prec float DEFAULT -1,
	inpanel VARCHAR(1000) DEFAULT '-',
	causative boolean DEFAULT false,
	automatic_selection boolean DEFAULT false,
	manual_selection boolean DEFAULT false,
	ai_selection boolean DEFAULT false,
	manual_revision boolean DEFAULT false,
	gene_matching boolean DEFAULT false,
	validated int DEFAULT -1,
	notes VARCHAR(100000) DEFAULT '-'
);

/*A table for the status of any analysis that will be introduced into the db*/
CREATE TABLE gene_matcher (
	familyid VARCHAR(1000) DEFAULT '-',
	genename VARCHAR(1000) DEFAULT '-'
);

--/*A table for the status of any analysis that will be introduced into the db*/
--CREATE TABLE research_centers (
	--centerid serial PRIMARY KEY,
	--centername VARCHAR(10000) DEFAULT '-',
	--centerlocation VARCHAR(10000) DEFAULT '-'
--);

--/*A table for the status of any analysis that will be introduced into the db*/
--CREATE TABLE analysis_status (
	--caseid int not null, 
	--analysisid int not null, 
	--receipt date,	
	--sequencing date,	
	--analysis date,	
	--interpretation date,		
	--report date,		
	--notes VARCHAR(100000) DEFAULT '-'
	--/*FOREIGN KEY(analysisid) REFERENCES analyses(analysisid)*/
--);
