#!/usr/bin/perl
#To work at CINECA
##!/opt/software/perl/bin/perl
#To work at TIGEM

#Works at CINECA
use lib '/pico/work/TELET_UDP/VarGeniusBeta'; #PICO CINECA SHARED
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5/x86_64-linux-thread-multi';
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5';
use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

# Author: Francesco Musacchia (2018)
#
#This software is a pipeline to discover SNPs and Variants through clinical
#sequencing data and to update our internal databases


#Libraries
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs
use File::Copy;#To manage files
use Getopt::Long;
use File::Copy;#To manage files
use Cwd;

#Using a library to manage files
use LIB::bedtools_management qw( run_BEDTOOLS_mergeBed run_BEDTOOLS_intersectBed);

#Using a library to manage files
use LIB::files_management qw(  delete_file save_hash load_hash delete_columns
														insert_line_into_file separate_elements_on_col);

#Files manipulation
use LIB::files_manipulation qw( extract_name );

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(hash2ConfigFile configFile2Hash try_exec_command 
			try_exec_job initialize_folders print_and_log log_and_exit
			separate_input_ids R_die sort_samples correct_type);

#DB management
use LIB::db_management qw(get_id_if_exists_from_db do_query_select_all);

use LIB::output_management qw( tag_panel_genes);

my $program_name  = "coverage_per_gene.pl";			
my $program_version = "0.1";
my $program_config = "program_config.txt";#The path to the configuration file for the program
my $config_folder = "CONFIGURATION";
my $dataFolder = "DATA";

my $working_folder = '';#Folder where to work
my $session_folder = '';#Folder where all the analyses will be written
my $program_folder = '';#Folder of the program


#Configuration hash
my $cfg_hash;

my $log_file = "";
my $foldersFile = "folders.txt";


#Selected genes
#my $sel_genes_f = $ARGV[0];
#my $given_samples = $ARGV[1];
#my $user_config = $ARGV[0];

#Input variables
my $user_config = "";   #It's mandatory
my $sel_genes_f = "";
my $run_analyses = "";
my $percentage_above = "";
my $user_id = "";
my $target = "";
my $input_f = "";
my $panel = "";
#Operations
my $gene_coverage = "";
my $cnv_analysis = "";

my $sep = ",";

parse_command_line_args();
if ( -e $user_config and $user_config ne ''){
	print "Configuring..\n";
	configureScript();
	#Set the main data folder
	 $cfg_hash->{'main_data_folder'} = $working_folder."/$dataFolder/";
	 $cfg_hash->{'lib_folder'} = $program_folder."/".$cfg_hash->{'lib_folder'};
}else{
	die  "Cannot open configuration file. Please give a correct parameter for --user_config! Exiting...\n";
}

#Get the names of the files output from DepthOfCoverage
my @doc_outs = split(",",$cfg_hash->{'DOC_outputs'});		

#Gete coverage per gene
if ( $gene_coverage ne '' ){
	get_gene_coverage();
}

#Get coverage per gene
if ( $cnv_analysis ne '' ){
	my @analyses_needed = ();
	#If used get the analyses needed	
	if ( $run_analyses eq "" ){
		print "You did not use run_analyses. \n";#,$log_file);
		#Get from the input
		#@analyses_needed = @analysis_ids;
	}else{
		@analyses_needed = split($sep,separate_input_ids($run_analyses,$sep));
	}
	
	foreach my $analysis_id (@analyses_needed){
		#Obtain the analysis_name from the database given the analyisis id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);

		my $proband_name = 	get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
														$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_kinship'},$analysis_id.",'P'");
		
		my $ed_out =  "$input_f/$proband_name\_ExomeDepth.csv";
		print "Analysis $analysis_name using $ed_out ...\n";#DEBUGCODE
		if ( -e $ed_out){
			my $log_file =  $input_f.".log";
			if ($panel ne ''){
				coverage_analysis_from_exomedepth($cfg_hash,$ed_out,$analysis_id,$log_file,$panel);	
			}else{
				coverage_analysis_from_exomedepth($cfg_hash,$ed_out,$analysis_id,$log_file);
			}
						
		}else{
			print "WARINING: Cannot find $ed_out ...\n";#DEBUGCODE
			}
	
	}	

	
}	



			
			
=head2 get_non_cov_reg_with_bedtools

 Title   : coverage_analysis_from_exomedepth
 Usage   : coverage_analysis_from_exomedepth(   );

 Function: 
					Uses the output from ExomeDepth to identify if there are relevant 
					CNVs
					Assumes that ExomeDepth is run for each analysis and .csv files are
					updated
					
 Returns : nothing
=cut
sub coverage_analysis_from_exomedepth{
	my $cfg_hash = shift;
	my $ed_out = shift;
	my $analysis_id = shift;
	my $log_file = shift;
	my $panel = shift;
	
	#Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Number of reads to be defined as LOW coverage
	my $deletion_thr = 0.1;#$cfg_hash->{'noncov_threshold'};

	#Check that also the user id is specified and correct
	if ( !( -e $ed_out ) ){
		die "File $ed_out does not exists! Exiting...\n";
	}

	#Obtain the analysis_name from the database given the analyisis id
	my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
	my $proband_name = 	get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
														$cfg_hash->{'db_analysis_id'}.",".$cfg_hash->{'db_sample_kinship'},$analysis_id.",'P'");
															
	#UCSC genes bed file
	my $ucsc_genes =  $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'ucsc_genes_bed'};
	my $inters_params = " -wao ";

	###################ANALYSIS OF HOMOZYGOUS DELETION
	#Execute AWK
	#print_and_log( "...filtering regions with reads ratio lt $deletion_thr (AWK)...",$log_file);#DEBUGCODE
	print "...filtering regions with reads ratio lt $deletion_thr (AWK)...";#DEBUGCODE
	my $ed_out_del = extract_name($ed_out,'noext')."_homdel.".$cfg_hash->{'bed_ext'};

	my $adapt_com = '{ print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7 }';
	my $command =	"cut -f5-7,9-12  $ed_out | awk "."'".'$7<'.$deletion_thr." $adapt_com' > $ed_out_del";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	#print "\nExecuting command: $command\n";#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";#####

	##Adding gene symbols
	#Intersect with ucsc genes regions to get the gene name	
	my $ed_out_del_genes = extract_name($ed_out_del,'noext')."_genes.".$cfg_hash->{'bed_ext'};
	#print_and_log( "...intersection with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	print "...intersection with ucsc genes $ucsc_genes (intersectBed)...";#DEBUGCODE
	run_BEDTOOLS_intersectBed($cfg_hash,$ed_out_del,$ucsc_genes,$ed_out_del_genes,$inters_params);
	
	my @cols_to_remove = (7,8,9);
	#print_and_log( "...deleting ucsc region start-end, sort and get unique lines...",$log_file);#DEBUGCODE
	print "...deleting ucsc region start-end, sort and get unique lines...";#DEBUGCODE
	delete_columns($ed_out_del_genes,\@cols_to_remove);
	$command =	" sort $ed_out_del_genes | uniq > $ed_out_del_genes.temp";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	#print "\nExecuting command: $command\n";#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";	

	my $new_header = "chr\tstart\tend\tPB\texpected\tobserved\tratio\t".$cfg_hash->{'gene_field'}."\tgenecov\n";
	insert_line_into_file($ed_out_del_genes.".temp",$new_header,0);
	
	#Separate multiple genes indication
	print "\nSeparating $ed_out_del_genes.temp into $ed_out_del_genes.temp2\n";
	separate_elements_on_col($ed_out_del_genes.".temp",$ed_out_del_genes.".temp2",8);
			
	#print_and_log( "Moving $ed_out_del_genes.temp2 to $ed_out_del_genes..\n",$log_file);#DEBUGCODE	
	move($ed_out_del_genes.".temp2",$ed_out_del_genes);	
	##print_and_log( "...DONE!",$log_file);#DEBUGCODE
	#print "...DONE!\n";#DEBUGCODE

	my $outFolder = "";
	my $out_path = "";
	if (defined $panel){
		#Now tag the panel genes and get a new output 
		$user_id = 147;
		$outFolder = getcwd;
		print "Now tagging panel genes into $ed_out_del_genes!\n";#DEBUGCODE
		$out_path = tag_panel_genes($cfg_hash,$analysis_id,$ed_out_del_genes,$panel,$user_id,$outFolder,$log_file);		
	}

	##Delete all temp files
	if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		delete_file($ed_out_del);
		delete_file($ed_out_del_genes.".temp");
		delete_file($ed_out_del_genes);		
	}


	###################ANALYSIS OF HETEROZYGOUS DELETIONS
	print "\n\nNow analyzying the heterozygous deletions...";#DEBUGCODE						
	my $cons_het_snps = $working_folder."/$analysis_name/finalout_out/$analysis_name\_het_cons.bed";
	my $rearr_ann_out =  $working_folder."/$analysis_name/finalout_out/$analysis_name.".$cfg_hash->{'ann_build_ver'}.".final_out_resgr14.txt";

	#Getting the regions with consecutive snps
	my $min_het_snps = 10;
	get_cons_het_snps_2_bed($rearr_ann_out,$cons_het_snps,$proband_name);
	my $cons_het_snps_hs = extract_name($cons_het_snps,'noext')."_hs.".$cfg_hash->{'bed_ext'};
	print "...filtering consecutive snps taking only those with at least $min_het_snps snps (AWK)...";#DEBUGCODE
	$command =	"awk "."'".'$4>'.$min_het_snps."'"." $cons_het_snps > $cons_het_snps_hs";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	print "\nExecuting command: $command\n";#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";#####
	
		
	#print_and_log( "...filtering regions with reads ratio lt $deletion_thr (AWK)...",$log_file);#DEBUGCODE
	my $het_del_min = 0.1;
	my $het_del_max = 1;
	print "...filtering het deletions with reads ratio lt $het_del_max and gt $het_del_min (AWK)...";#DEBUGCODE
	$ed_out_del = extract_name($ed_out,'noext')."_hetdel.".$cfg_hash->{'bed_ext'};
	$adapt_com = '{ print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7 }';
	$command =	"cut -f5-7,9-12  $ed_out | awk "."'".'$7<'.$het_del_max." && ".'$7>'.$het_del_min." $adapt_com' > $ed_out_del";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	print "\nExecuting command: $command\n";#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";#####

  #Intersect with regions with consecutive het snps
	my $ed_out_del_conssnp = extract_name($ed_out_del,'noext')."_conssnp.".$cfg_hash->{'bed_ext'};
	#print_and_log( "...intersection with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	print "...intersection with consecutive het snps $cons_het_snps_hs (intersectBed)...";#DEBUGCODE
	$inters_params = " -f 0.50 ";#overlap dimension
	run_BEDTOOLS_intersectBed($cfg_hash,$ed_out_del,$cons_het_snps_hs,$ed_out_del_conssnp,$inters_params);
	
	##Adding gene symbols
	#Intersect with ucsc genes regions to get the gene name	
	$ed_out_del_genes = extract_name($ed_out_del,'noext')."_genes.".$cfg_hash->{'bed_ext'};
	#print_and_log( "...intersection with ucsc genes (intersectBed)...",$log_file);#DEBUGCODE
	print "...intersection with ucsc genes $ucsc_genes (intersectBed)...";#DEBUGCODE
	$inters_params = " -wao ";
	run_BEDTOOLS_intersectBed($cfg_hash,$ed_out_del_conssnp,$ucsc_genes,$ed_out_del_genes,$inters_params);
	
	@cols_to_remove = (7,8,9);
	#print_and_log( "...deleting ucsc region start-end, sort and get unique lines...",$log_file);#DEBUGCODE
	print "...deleting ucsc region start-end, sort and get unique lines...";#DEBUGCODE
	delete_columns($ed_out_del_genes,\@cols_to_remove);
	$command =	" sort $ed_out_del_genes | uniq > $ed_out_del_genes.temp";
	#print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	print "\nExecuting command: $command\n";#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";	

	$new_header = "chr\tstart\tend\tPB\texpected\tobserved\tratio\t".$cfg_hash->{'gene_field'}."\tgenecov\n";
	insert_line_into_file($ed_out_del_genes.".temp",$new_header,0);
	
	#Separate multiple genes indication
	print "\nSeparating $ed_out_del_genes.temp into $ed_out_del_genes.temp2\n";
	separate_elements_on_col($ed_out_del_genes.".temp",$ed_out_del_genes.".temp2",8);
			
	#print_and_log( "Moving $ed_out_del_genes.temp2 to $ed_out_del_genes..\n",$log_file);#DEBUGCODE	
	move($ed_out_del_genes.".temp2",$ed_out_del_genes);	
	##print_and_log( "...DONE!",$log_file);#DEBUGCODE
	#print "...DONE!\n";#DEBUGCODE

	if (defined $panel){
		#Now tag the panel genes and get a new output 
		
		$user_id = 147;
		$outFolder = getcwd;
		print "Now tagging panel genes into $ed_out_del_genes!\n";#DEBUGCODE
		$out_path = tag_panel_genes($cfg_hash,$analysis_id,$ed_out_del_genes,$panel,$user_id,$outFolder,$log_file);		
	}


	###Delete all temp files
	#if ( $cfg_hash->{'remove_temp'} eq 'YES' ){
		#delete_file($ed_out_del);
		#delete_file($ed_out_del_genes.".temp");
		#delete_file($ed_out_del_genes);		
	#}
	
	
}

###Given the output from VarGenius searches
#the intervals where 2 or more het snps are consecutives
#takes the start and the end of each of this intervals
sub get_cons_het_snps_2_bed{
	my $vgout = shift;
	my $outname = shift;
	my $sample_name = shift;
	
	print "Writing the bed for consecutive het snps into $outname\n";
	
	my $chr_ind = 0;
	my $pos_ind = 1;
	my $gt_field = $sample_name."_GT";
	#print "GT field $gt_field\n";
	my $gt_field_ind = -1;
	####################
	open (FILE,"<$vgout") or die "Cannot open $vgout\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	my $gt_count = 0;
	my $last_het_pos = -1;
	my $chr2keep = "";
	my $lineclosed = 0;
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		if ($numline == 0){
			my $numcol = 0;
			foreach my $field (@fields){
					if ( $field eq $gt_field){
						$gt_field_ind = $numcol;
					}
					$numcol++;
			}
		}else{

			my $gt_val = $fields[$gt_field_ind];
			#print "Checking field $gt_field_ind :$gt_val\n";
			#In case of HET
			if ( $gt_val eq '0/1' or  $gt_val eq '1/0' or $gt_val eq '0|1' or  $gt_val eq '1|0'){
				#print "Het SNP at line $numline:$gt_val\n";
				$gt_count++;
				#Initialize the chromosome when is the first HET
				if ($gt_count == 1){
					$lineclosed = 0;
					$chr2keep = $fields[$chr_ind];
				}
				#If this HET is the same chromosome of the one before
				if ( $fields[$chr_ind] eq $chr2keep){
					#If it is the first HET, start to write position
					if ($gt_count == 1){
						$chr2keep = $fields[$chr_ind];
						print NEWFILE $fields[$chr_ind]."\t".$fields[$pos_ind]."\t";
					}					
				}
				#If this HET is on a new chromosome..
				else{
					print NEWFILE $last_het_pos."\t".$gt_count."\n".$fields[$chr_ind]."\t".$fields[$pos_ind]."\t";
					$chr2keep = $fields[$chr_ind];
					$gt_count = 1;
					$lineclosed = 1;
				}
				#Mark the last HET position
				$last_het_pos = $fields[$pos_ind];
			}
			#In case of HOM
			else{
				#If there is at least one het snp, print the end
				if ($gt_count >= 1){
						print NEWFILE $last_het_pos."\t".$gt_count."\n";
						$gt_count = 0;
						$lineclosed = 1;
				}				
			}
		}
		$numline++;
	}
	#If the line wasn't closed, do it
	if ( $lineclosed == 0 ){ print NEWFILE $last_het_pos."\t".$gt_count."\n";}
	close(FILE);
	close(NEWFILE);	
}

		
sub get_gene_coverage {
	
	print "Parameters needed:\n";
	print "percentage_above=$percentage_above\n";
	print "user_id=$user_id\n";
	print "sel_genes_f=$sel_genes_f\n";
	print "target=$target\n";
		
	#Check that also the percentage above is specified and correct
	if ( !( correct_type($percentage_above,"positiveint") and $percentage_above ne '') ){
		die "Percentage above for the coverage must be given as a positive integer with --percentage_above! Exiting...\n";
	}
	#Check that also the user id is specified and correct
	if ( !( correct_type($user_id,"positiveint") and $user_id ne '') ){
		die "User id must be given as a positive integer with --user_id! Exiting...\n";
	}
	#Check that also the user id is specified and correct
	if ( $target eq '' ){
		die "Target file name must be given with --target! Exiting...\n";
	}
	#Get the names of the files output from DepthOfCoverage
	my @doc_outs = split(",",$cfg_hash->{'DOC_outputs'});	
	#Check if selected genes are inputed
	if ( $sel_genes_f eq ''){
			$sel_genes_f  = "ALL_GENES";
	}
	#Set the name of the complete table for all samples from the database
	my $genecove_above_all = $working_folder."/DATA/$user_id\_$target\_".extract_name($sel_genes_f,1)."_".$doc_outs[0]."_above$percentage_above.".$cfg_hash->{'txt_ext'};

	#Groups for which we want to do the job
	my @analysis_ids = ();
	my @analyses_needed = ();															
	#collects all the group ids related with this userid which could represent a project and that use the same target file
	my $query = "SELECT ".$cfg_hash->{'db_analysis_id'}." FROM  ".$cfg_hash->{'db_analyses_table'}." WHERE ".
										$cfg_hash->{'db_analysis_userid'}." = $user_id AND ".$cfg_hash->{'db_targetbed'}." = '".$target."'".
										" AND ".$cfg_hash->{'db_analysis_infreq'}." = 1;";										
	print "Executing: $query\n";#,$log_file);
	my $analysis_id_fetch = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_analysis_id'});

	#Now put into a normal array
	foreach my $analysis_id (keys %{$analysis_id_fetch}){
		push(@analysis_ids,$analysis_id);
	}
		
	#If used get the analyses needed	
	if ( $run_analyses eq "" ){
		print "You did not use run_analyses. All the analyses for the given user id and target will be used\n";#,$log_file);
		#Get from the input
		@analyses_needed = @analysis_ids;
	}else{
		@analyses_needed = split($sep,separate_input_ids($run_analyses,$sep));
	}

	my $groups_count = 1;

	#if ( scalar(keys %{$analysis_ids}) > 0){
	if ( scalar(@analysis_ids) > 0){

		#Open and close just to get an empty file
		open (GC_ALL,">$genecove_above_all ");close(GC_ALL);
		#foreach my $group_id (keys %{$group_ids}){
		foreach my $analysis_id (@analysis_ids){
			if ( grep {/\b$analysis_id\b/} @analyses_needed) {
				print "Getting gene coverage for analysis $analysis_id\n";
				my $analysis_genecove_above = get_gene_coverage_above( $cfg_hash,$analysis_id,$sel_genes_f,$percentage_above);
				
				#If this is the first file then copy to the all file
				if ( $groups_count == 1 ){
					print "Copying $analysis_genecove_above to $genecove_above_all\n";
					copy($analysis_genecove_above,$genecove_above_all) or die "Cannot copy $analysis_genecove_above in $genecove_above_all";		
				}else{
					my $R_utils =  $cfg_hash->{'lib_folder'}."/".$cfg_hash->{'R_utils'};
					
					my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args MERGE $genecove_above_all $analysis_genecove_above Gene $genecove_above_all' $R_utils";
				 
					print "The command is:\n $command\n";
					try_exec_command($command) or R_die($command,$R_utils);			
				}
				#Increase groups counter
				$groups_count++;			
			}
		}
		print "Finished writing the table: $genecove_above_all\n";
	}
}			


sub get_gene_coverage_above {
		my $cfg_hash = shift;
		my $analysis_id = shift;
		my $sel_genes_f = shift;
		my $percentage_above = shift;
		
		
		#Get the gene summary file name
		#Obtain the analysis_name from the database given the analyisis id
		my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
														$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
														$cfg_hash->{'db_analysis_id'},$analysis_id);
													
		my $gene_summary = $working_folder."/".$analysis_name."/".$cfg_hash->{'stats_out_f'}."/".$analysis_name."_DOC.".$doc_outs[0];
		
		#Getting the samples list associated to the analysisid
		my $samples_h;#Samples hash to be reordered
		my @sort_samples = ();#The array with samples id reordered
		#Get the samples names involved for the analysis
		my $query = "SELECT ".$cfg_hash->{'db_sample_name'}." FROM  ".$cfg_hash->{'db_sample_table'}." WHERE ".$cfg_hash->{'db_analysis_id'}."=$analysis_id;";	
		print "Executing: $query\n";#DEBUGCODE
		my $group_sam = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_name'});
		#Print the header for information about the genotype for the sample
		my $sample_gen_fields1 = "";
		my $sample_gen_fields2 = "";
		#Get the kinship, to make the resorting
		foreach my $sample_name (keys %{$group_sam}){  #Print version

			#Obtain the kinship from the database given the sample id
			my $kinship = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_kinship'},
															$cfg_hash->{'db_sample_name'},"'".$sample_name."'");	
			#Build an hash to reorder
			$samples_h->{$sample_name}->{'id'} = $sample_name;
			$samples_h->{$sample_name}->{'k'} = $kinship;
		}
		#Sort the samples as in samples_order parameter
		sort_samples($samples_h,$cfg_hash->{'samples_order'},\@sort_samples);	
		
		my $samples_l = join(",",@sort_samples);
		
		###############################
		# Rscripts
		###############################
		
		#Now call the R script that will write a new table
		my $analysis_genecove_above = $working_folder."/".$analysis_name."/".$cfg_hash->{'stats_out_f'}."/".$analysis_name."_DOC_".extract_name($sel_genes_f,1)."_".$doc_outs[0]."_above$percentage_above.".$cfg_hash->{'txt_ext'};
		
		#
		my $RLogPath = $working_folder."/".$analysis_name."/".$cfg_hash->{'stats_out_f'}."/log/$analysis_name\_genecovabove$percentage_above".$cfg_hash->{'R_log_file'};
		my $R_script = $cfg_hash->{'lib_folder'}."/coverage_per_gene.R";
		
		my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $gene_summary $sel_genes_f $samples_l above_$percentage_above $analysis_genecove_above' $R_script $RLogPath";
   
		print "The command is:\n $command\n";
		try_exec_command($command) or R_die($command,$R_script);
		
		return $analysis_genecove_above;
}


sub configureScript {
	#print "Inizializing folders with $foldersFile...\n";
	($working_folder,$program_folder) = initialize_folders($foldersFile);
	$config_folder = $program_folder."/".$config_folder;
	$program_config = $config_folder."/".$program_config;


	#loads the parameters in the configuration hash as from the configuration files
	print "Loading the parameters for the run as from $program_config...\n";	
	configFile2Hash($program_config,\$cfg_hash);
	
	print "Loading the parameters for the run as from $user_config...\n";	
	configFile2Hash($user_config,\$cfg_hash);
	
	print "Configuration completed...\n";	
}	


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
 sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

	my $howToUse = "Use with: \nperl $program_name \n\n".
	"\t-c|--user_config: The configuration file from the user.\n".
	"\t-s|--sel_genes_f: A text file with a list of needed genes to filter.\n".
	"\t-p|--percentage_above: The percentage to use to filter the gene coverage percentage table.\n".
	"\t-ra|--run_analyses: Analysis ids to consider for the table\n".
	"\t-uid|--user_id: user id to consider for the table\n".
	"\t-t|--target: target file to consider for the table\n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "version" => \$VERSION,
           "c|user_config=s" => \$user_config,   #It's mandatory
           "gene_cov|gene_coverage" => \$gene_coverage,
           "cnv|cnv_analysis" => \$cnv_analysis,
           "s|sel_genes_f=s" => \$sel_genes_f,
           "i|input_f=s" => \$input_f,
           "perc|percentage_above=s" => \$percentage_above,
           "pan|panel=s" => \$panel,
           "ra|run_analyses=s" => \$run_analyses,
           "uid|user_id=s" => \$user_id,
           "t|target=s" => \$target
            );

  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "version: $program_version \n";
    exit;
  }

}
