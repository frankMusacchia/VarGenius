#!/usr/bin/perl
#To work at CINECA
##!/opt/software/perl/bin/perl
#To work at TIGEM


#GLI HASHTAG #TOCHECK E .. SONO DA CONTROLLARE PER CODICE DA AGGIUNGERE

#Works at CINECA
#use lib '/pico/home/userexternal/fmusacch/VarGenius/'; #PICO CINECA HOME
use lib '/pico/work/TELET_UDP/VarGenius'; #PICO CINECA SHARED
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5/x86_64-linux-thread-multi';
use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5';
#Parallel::ForkManager
#use lib qw(/pico/home/userexternal/fmusacch/bin/perl_lib/Parallel-ForkManager-1.17/lib/);
use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

#TIGEM
use lib '/home/users/musacchia/VarGenius';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED


# Package VarGenius - Variant Discovery and Annotation Tool
# Author: Francesco Musacchia (2017)
#
#This software is a pipeline to discover SNPs and Variants through clinical
#sequencing data and to update our internal databases
#This tool is needed to intersect the target regions chosen for a particular
#Whole Exome analysis with the coding regions to obtain a target only for
#coding regions


#Libraries
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Time::HiRes qw( time ); #To compute the running time of jobs

#Using a library to manage files
use LIB::bedtools_management qw( run_BEDTOOLS_mergeBed run_BEDTOOLS_intersectBed);

#Using a library to manage files
use LIB::files_management qw(  delete_file save_hash load_hash);

#Files manipulation
use LIB::files_manipulation qw( extract_name);

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(hash2ConfigFile configFile2Hash try_exec_command 
			try_exec_job initialize_folders print_and_log log_and_exit
			separate_input_ids);
			
my $program_config = "program_config.txt";#The path to the configuration file for the program
my $config_folder = "CONFIGURATION";

my $working_folder = '';#Folder where to work
my $session_folder = '';#Folder where all the analyses will be written
my $program_folder = '';#Folder of the program

my $task = $ARGV[0];
#Configuration hash
my $cfg_hash;

my $log_file = "";
my $foldersFile = "folders.txt";



sub configureScript {
	#print "Inizializing folders with $foldersFile...\n";
	($working_folder,$program_folder) = initialize_folders($foldersFile);
	$config_folder = $program_folder."/".$config_folder;
	$program_config = $config_folder."/".$program_config;


	#loads the parameters in the configuration hash as from the configuration files
	print "Loading the parameters for the run as from $program_config...\n";	
	configFile2Hash($program_config,\$cfg_hash);
}	

=head2 add_col_with_chrom_code

 Title   : add_col_with_chrom_code
 Usage   : add_col_with_chrom_code(   );

 Function: Adds a column to a bed file where there is a chromosome code
				formed by the chromosome number-start-end-strand
				
 Returns : nothing
=cut
sub add_col_with_chrom_code{
	
	my $infile = shift;
	my $outfile = shift;
	
	my $chr_num_ind = 0;
	my $start_ind = 1;
	my $end_ind = 2;
	my $strand_ind = 5;
	
	open (IN,"<$infile") or die "Cannot open $infile\n";
	open (OUT,">$outfile") or die "Cannot open $outfile\n";
	while (my $row = <IN>){
			chomp($row);
			my @fields = split ("\t",$row);
			my $code = $fields[$chr_num_ind]."-".$fields[$start_ind]."-".$fields[$end_ind]."-".$fields[$strand_ind];
			 
			print OUT "$row\t$code\n";
			#print "$num\t$row";
	}
	close(IN);
	close(OUT);
	
}


=head2 sort_by_chrom_and_startpos

 Title   : sort_by_chrom_and_startpos
 Usage   : sort_by_chrom_and_startpos(   );

 Function: Runs get_file_with_chrom_num
 Returns : nothing
=cut
sub sort_by_chrom_and_startpos{
	my $infile = shift;
	my $outfile = shift;
	
	
	my $command = "sort -k1,1 -k2,2n $infile > $outfile";
	
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	
}

=head2 get_unique_lines

 Title   : get_unique_lines
 Usage   : get_unique_lines(   );

 Function: Returns a file with unique rows where the identfiers to be unique
					are in the ind_col column
 Returns : nothing
=cut
sub get_unique_lines{
	my $infile = shift;
	my $ind_col = shift;
	my $outfile = shift;
	
	my $command = "uniq -f $ind_col $infile > $outfile";
	
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 add_chromnum_to_bed

 Title   : add_chromnum_to_bed
 Usage   : add_chromnum_to_bed(   );

 Function: Puts the chromosome number as first element of the bed file
 
 Returns : nothing
=cut
sub add_chromnum_to_bed{
	
	my $infile = shift;
	my $outfile = shift;
	
	open (IN,"<$infile") or die "Cannot open $infile\n";
	open (OUT,">$outfile") or die "Cannot open $outfile\n";
	while (my $row = <IN>){
			my @fields = split ("\t",$row);
			my $num = $fields[0];
			$num =~ s/chr//;
			 
			print OUT "$num\t$row";
			#print "$num\t$row";
	}
	close(IN);
	close(OUT);
	
}

=head2 intersect_target_with_cod_regions

 Title   : intersect_target_with_cod_regions
 Usage   : intersect_target_with_cod_regions(   );

 Function: Intersects a target file in input with given coding regions
						Both files must be BED format.
						
						
 Returns : nothing
=cut
sub intersect_target_with_cod_regions{
	my $target_bed = shift;
	my $cod_regions = shift;
	my $out_target = shift;

	print "Starting with target $target_bed and coding reg: $cod_regions and log:$log_file...\n";
	
	#Put the chromosome id in front of the target file
	my $sortedTargetBed = "sortedTarget.bed";
	sort_by_chrom_and_startpos($target_bed,$sortedTargetBed);

	#combines overlapping features in an interval file in a single feature
	my $mergedSortTargBed = "mergedSortTarg.bed";
	run_BEDTOOLS_mergeBed($cfg_hash,$sortedTargetBed,$mergedSortTargBed,"",$log_file);
	
	#intersect the target with the coding regions given in input
	my $interTargCodBed = "interTargCod.bed";
	run_BEDTOOLS_intersectBed($cfg_hash,$cod_regions,$mergedSortTargBed,$interTargCodBed,"",$log_file);
	
	#sorts the output file per chromosome id
	my $interTargCodBedSort = "interTargCodSort.bed";
	sort_by_chrom_and_startpos($interTargCodBed,$interTargCodBedSort);

	#Adds a column with a chromosome id		
	my $interTargCodBedChromid = "interTargCodChromid.bed";
	add_col_with_chrom_code($interTargCodBedSort,$interTargCodBedChromid);
	
	#Gets unique lines using the chromosome id (column 6)
	#my $interTargCodBedChromidUniq = "interTargCodChromidUniq.bed";
	get_unique_lines($interTargCodBedChromid,6,$out_target);
	
	#Remove all the temporary files
	print "Removing temporary files...\n";
	delete_file($sortedTargetBed);
	delete_file($mergedSortTargBed);
	delete_file($interTargCodBed);
	delete_file($interTargCodBedSort);
	delete_file($interTargCodBedChromid);
}

print "Configuring..\n";
configureScript();
my @tasks = ('TARG_COD_INTERS');
my $tasks_st = 'TARG_COD_INTERS';
if ( !(grep {$task} @tasks) ){ die "Unknown task used: $task. Please select one among $tasks_st\n";}

if ($task eq 'TARG_COD_INTERS'){
	die "Needed arguments: target_bed cod_regions out_target" unless scalar(@ARGV == 4);
	print "Running $task\n";
	my $target_bed =  $ARGV[1];
	my $cod_regions = $ARGV[2];
	my $out_target = $ARGV[3];
	
	$log_file = "TARG_COD_INTERS.log";
	intersect_target_with_cod_regions($target_bed,$cod_regions,$out_target);
}
