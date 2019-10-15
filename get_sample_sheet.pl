#!/opt/software/perl/bin/perl
#To work at TIGEM
#!/usr/bin/perl
#To work at CINECA

use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED


#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2017>  <Francesco Musacchia>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    
# Package VarGenius - Variant Discovery and Annotation Tool
# Author: Francesco Musacchia (2017)
#
#This software is a pipeline to discover SNPs and Variants through clinical
#sequencing data and to update our internal databases


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;#Used to write an usage

use FindBin;#To search libraries
use lib $FindBin::Bin;#To search libraries

#Using a library to manage files 
use LIB::files_management qw( join_files_with_zcat );

use LIB::programs_management qw( initialize_folders );

my $program_name = "get_sample_sheet.pl";
my $program_version = "1.0";

#Input variables
my $target_bed = "";
my $intern_fold = "";
my $perchrom = "";
my $reference = "";
my $convscores = "";
my $user_id = "";
my $research_group_id = "";
my $add_to_analysisname = "";
my $change_analysisname = "";
my $infreq = "";
my $set_kinship = "";
my $output = "";
my $head_file = "";
my $dir = "";
my $mode = "";
my $mult_lanes = "";
my $fileList = "";
my $r1_suff = "";
my $r2_suff = "";
my $lanes_list = "";


my $foldersFile = "folders.txt";
my $ss_head_file = "CONFIGURATION/ss_head_file.txt ";

           
#Standard values for parameters
my $motor_target = "motor7_targeted.bed";
my $exome_target = "clinical_exome_cod.bed";
my $rna_target = "rna.bed";
my $base_perchrom = 0;
my $base_reference = "hg19";
my $base_convscores = 0;#Convert scores
#Users and research group ids
my $base_user_id = 142;#Motor 142 UDP 133
my $base_research_group_id = 14;

my $base_add_to_analysisname = "";#"";#
my $base_change_analysisname = "";#If this is not '-' the analysis name will be identical for all
#The analysis must be kept in count for the frequencies computation
my $base_infreq = 0;
#Folder that is inside each sample folder before to arrive to fastq files
my $base_intern_fold = "";#"/Files/";

#Reads 1 and 2 suffix
my $base_r1_suff = '_R1_';
my $base_r2_suff = '_R2_';

#Indicator of lanes
my @lanes = ("_L001_","_L002_","_L003_","_L004_","_L005_","_L006_","_L007_","_L008_","_L009_");
####################################MAIN

#Start getting the input parameters
parse_command_line_args();


#Get the working folder and vargenius	
my ($working_folder,$program_folder) = initialize_folders($foldersFile);
my $start = 1;	


#Check mandatory values
#if ( $head_file eq ''){
		#die  "Header file is needed. Please use --head_file\n";
#}
$head_file = $program_folder."/".$ss_head_file;
if ( $dir eq ''){
		die  "Directory where to search samples is needed. Please use -d-ir\n";
}
if ( $mode eq ''){
		die  "mode is needed [exome|targeted|rna]. Please use --mode\n";
}



#Set variables whether they are inputed or not
unless ( $intern_fold ne "" ) {$intern_fold = $base_intern_fold; }
unless ( $perchrom ne "" ) {$perchrom = $base_perchrom;}
unless ( $reference ne "" ) {$reference = $base_reference;}
unless ( $convscores ne "" ) {$convscores = $base_convscores; }
unless ( $user_id ne "" ) {$user_id = $base_user_id; }
unless ( $research_group_id ne "" ) {$research_group_id = $base_research_group_id; }
unless ( $r1_suff ne "" ) {$r1_suff = $base_r1_suff; }
unless ( $r2_suff ne "" ) {$r2_suff = $base_r2_suff; }
unless ( $lanes_list eq "" ) {@lanes = split(",",$lanes_list); }


my $ext_to_search = "gz";
my $exts_to_search = "fq.gz|fastq.gz|fastq|fq;";#('.fq.gz','.fastq.gz','.fastq','fq');

#Groupname	
unless ( $add_to_analysisname ne "" ) {$add_to_analysisname = $base_add_to_analysisname; }
unless ( $change_analysisname ne "" ) {$change_analysisname = $base_change_analysisname; }
#The analysis must be kept in count for the frequencies computation
unless ( $infreq ne "" ) {$infreq = $base_infreq;}


#Print what will be used (if the output goes somewhere else...)
if ($output ne ''){
	
if ( $output eq ''){
		print  "You did not use the file path.. ( --output)\n";
}
if ( $mult_lanes eq ''){
		print  "You did not use the  mult_lanes parameters...\n";
}
if ( $fileList eq ''){
		print  "You did not use a file with the list of sample names [-fileList]...\n";
}

	
	print "head_file: $head_file\n";
	print "dir: $dir\n";
	print "mode: $mode\n";
	print "mult_lanes: $mult_lanes\n";
	print "fileList: $fileList\n";
	print "Internal folder: $intern_fold\n";
	print "perchrom : $perchrom\n"; 
	print "reference: $reference\n";	
	print "convscores: $convscores\n";
	print "user_id: $user_id\n";
	print "research_group_id: $research_group_id\n";
	print "add_to_analysisname: $add_to_analysisname\n";
	print "change_analysisname: $change_analysisname\n";
	print "infreq: $infreq\n";
	print "output: $output\n";
	print "Read1 suffix: $r1_suff\n";
	print "Read2 suffix: $r2_suff\n";
	print "lanes list: $lanes_list\n";
}
	
if ( $start){
	#write_sample_sheet_both_selected();
	write_sample_sheet();
}


=head2 get_gender_hash

 Title  : get_gender_hash
 Usage  : get_gender_hash( - dirname => 'the folder to check',
                      );

 Function:                       
         If the user gives in input a file with: sample_name\tgender
         the gender is picked from the file and inserted into an hash.
         Then this function fills the list of sample names
         
 
 Returns : imputed hash and array are filled

=cut
sub get_gender_hash {
	my ($gender_hash) = shift;
	my ($list) = shift;
	my $fileList = shift;
	
	my $retval = 0;
	open(LIST,"<$fileList") or die "Cannot open $fileList..\n";

	while ( my $row = <LIST>){
		chomp($row);
		my @fields = split("\t",$row);
		if (defined $fields[1]){
			$retval = 1;
			$$gender_hash->{$fields[0]} = $fields[1];					
		}
		push(@$list,$fields[0]);
	}		
	close(LIST);
	
	return $retval;
}

=head2 get_sample_info_hash

 Title  : get_sample_info_hash
 Usage  : get_sample_info_hash( - dirname => 'the folder to check',
                      );

 Function:                       
         This files contains fields with samples info: 
         sample_name\tgender\tdob\tpob\tkinship\taffected\n
         
         infos are picked from the file and inserted into an hash.
         Then this function fills the list of sample names
         
 
 Returns : imputed hash and array are filled

=cut
sub get_sample_info_hash {
	my ($sample_info_hash) = shift;
	my ($list) = shift;
	my $fileList = shift;
	
	my $retval = 0;
	open(LIST,"<$fileList") or die "Cannot open $fileList..\n";

	while ( my $row = <LIST>){
		chomp($row);
		my @fields = split("\t",$row);
		if (defined $fields[1]){
			$retval = 1;
			#Gender
			$$sample_info_hash->{$fields[0]}->{'gender'} = $fields[1];
			$$sample_info_hash->{$fields[0]}->{'dob'} = $fields[2];
			$$sample_info_hash->{$fields[0]}->{'pob'} = $fields[3];
			$$sample_info_hash->{$fields[0]}->{'affected'} = $fields[4];					
		}
		push(@$list,$fields[0]);
	}		
	close(LIST);
	
	return $retval;
}
=head2 write_sample_sheet_both_selected

 Title  : write_sample_sheet_both_selected
 Usage  : write_sample_sheet_both_selected( - dirname => 'the folder to check',
                      );
                      
              dir: source directory
              mode: motor or exome or exome_gr
							mult_lanes: can be 'lanes' or something other
							
							Use exome_gr if you know there are multiple samples inside a folder
							
							set the intern_fold variable if needed
							
 Function: Writes a sample sheet by using the folder of the fastq files.
					MotorV7: motor7_haloplex.bed
					MotorV8: MotorPlex_v8_TuttiGeni_HPnormale_Covered.bed
					
 Returns : nothing

=cut
sub write_sample_sheet_both_selected{

	my $sep = "\t";
	my $no_val = "-";
	

	#my $intern_fold = "";#"/Files/";

	#Put the names from the file in an array
	my @list = ();
	my @samples_used = ();
	my $gender_hash;
	my $gender_taken = 0;
	
	if (defined $fileList and $fileList ne 'none' and $fileList ne ''){
		if ( get_gender_hash(\$gender_hash,\@list,$fileList) == 1){
			$gender_taken = 1;
		}
		else{
			@list = list_to_array($fileList,'NO_NEW_LINE');
		}
	}
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs= readdir DIR;#All the files are in this list

	
	my $final_info = "";
	if ($mode eq "targeted"){
		$perchrom = 0;
		$final_info = "targeted".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		#$target_bed = "MotorPlex_v8_TuttiGeni_HPnormale_Covered.bed";
		unless ( $target_bed ne "" ) {$target_bed = $motor_target };
	}elsif ($mode eq "exome"){#Exome single sample
		#$perchrom = 1;
		$final_info = "exome".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $exome_target;}
	}elsif ($mode eq "exome_gr"){#Exomes
		$final_info = "exome".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $exome_target;}
	}elsif ($mode eq "rna"){#RNA
		$final_info = "rna".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $rna_target;}
	}else{
			die "Need a mode to insert the parameters! [exome, motor, exome_gr, rna]\n";
	}

	
	#Start to write the file 
	if ($output ne ''){
		open (OUT,">$output") or die "ERROR: Cannot open $output. The program will exit..\n";
	}
	#print "There are ".scalar(@samples_dirs)." elements in $dir\n";#DEBUGCODE
	my $couple;
	#Get and print header
	open (HEAD,"<$head_file") or die "ERROR: Cannot open $head_file. The program will exit..\n";
	my $header = <HEAD>;
	close(HEAD);
	#Substitute commas in the header with the separator used here
	$header =~ s/,/$sep/g;
	if ($output ne ''){
		print OUT $header;
	}else{
		print $header;
	}
	#A control starts on each sample among all the samples        
	my $numFile= 1;
	my $fastq1 = "";
	my $fastq2 = "";
	
	#Flag used to check if the sample is to take or not
	my $get_sel_only = 0;
	my $analyze_sample = 1;

	
	#Go through the different samples directories
	foreach my $sample_dir (@samples_dirs){
			#print "Considering sample $sample_dir \n";#DEBUGCODE
			
			my $analysis_name = $sample_dir;			
			##########CONSIDER SAMPLELIST
			#If the file list is defined set the analyze_sample flag to NO
			if (defined $fileList and $fileList ne 'none'  and $fileList ne ''){
						$get_sel_only = 1;
						if ( not grep {/\b$sample_dir\b/} @list ){
							$analyze_sample = 0;
							#print "The sample $sample_dir will not be analyzed..\n";#DEBUGCODE
						}else{
							#print "The sample $sample_dir will be analyzed..\n";#DEBUGCODE
						}
			}
			
			#Get all the samples unless the user wants to get only those selected and this sample is ok
			unless ( $get_sel_only == 1 and $analyze_sample == 0){
			
				#If the folder name does not start with a dot and is a  directory
				if ( ($sample_dir !~ /^\..*/) and (-d $dir.$sample_dir)){

					#print "$sample_dir will be analyzed.\n";#DEBUGCODE
					my @lanes_used = ();
					my $kinship = "-";
					my $gender = "-";

					############SET KINSIP
					#Kinship is associated when the last letter is given
					#More than one proband are specified with the number after the _P
					#The groupname is picked from the proband sample name
					if ( $set_kinship){
						if ($sample_dir =~ /_P\d*$/ ) {
							$kinship = "P";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_P\d*$//;
						}elsif ($sample_dir =~ /_M$/) {
							$kinship = "M";
							$gender  ="F";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_M$//;
						}elsif ($sample_dir =~ /_F$/) {
							$kinship = "F";
							$gender  ="M";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_F$//;
						}else{
							$kinship = "P";
							$gender  = "-";
							#pick the analysis name from the sample name 
							$analysis_name = $sample_dir;
						}
					}
					
					##############SET GENDER
					#If the gender has been put into 
					if ( defined $gender_hash->{$sample_dir}){
						$gender  = $gender_hash->{$sample_dir};
					}
					
					
					############ADD SUFFIX TO GROUPNAME
					#If the group name is not specified than it is the same as the sample dir
					if ( ($mode eq 'targeted') or ($mode eq 'exome') or ($mode eq 'rna')){
						#If change_analysisname is set, then it is used. Otherwise add_to is verified 
						if ($change_analysisname ne "" ){
							$analysis_name = $change_analysisname;
						}else{
							#If I already changed the groupname in set_kinship
							if ($set_kinship){
								$analysis_name = $analysis_name.$add_to_analysisname;
							}else{
								$analysis_name = $sample_dir.$add_to_analysisname;
							}
						}
					}	

					#print "Entering $sample_dir...\n";
					my $dir1 = $dir."/".$sample_dir;
					#print $dir1."\t";
					$dir1 .= "/$intern_fold/";
					opendir DIR1, $dir1 or die "ERROR: cannot open dir $dir1";
					#my @files= readdir DIR1;#All the files are in this list
					my @files = grep { /\.($exts_to_search)$/ } readdir DIR1;
					die "ERROR: Cannot find files with extensions $exts_to_search in $dir1\n".
							"Maybe there is an internal folder? (intern_fold: $intern_fold)" unless @files > 0;
				
					#A control starts on each file of the folder        
							
					############MULTIPLE LANES
					############					
					if ( $mult_lanes ){
						my $ext_found = 0;
						foreach my $file (@files){
							#print "Using $file..\n";#DEBUGCODE
								$ext_found = 1;
								foreach my $lane (@lanes){
									#print "Checking if $lane is contained...\n";
									if ( $file =~ /$lane/ ){
										#print "$file contains $lane...\n";
										#Check if is the read 1 or two
										if ( $file =~ /$r1_suff/ ){
											$couple->{$lane}->{1} = $file;
											push(@lanes_used,$lane);
										}else{
											$couple->{$lane}->{2} = $file;
										}
									}
								}								
						}
						
						#Error not found fq
						die "ERROR: Cannot find fastq files ($ext_to_search) in $dir1. Exiting\n"
							unless ($ext_found);
													
						#Print both the files for each lane used
						#print "for folder $sample_dir used ".scalar(@lanes_used)."\n";
						foreach my $lane (@lanes_used){
							#heritability,gender,diseaseID,kinship are NOVAL
							if ($output ne ''){
								print OUT $no_val.$sep.$gender.$sep.$no_val.$sep.$kinship.$sep.$user_id.$sep.$research_group_id.$sep;
								print OUT $dir1.$sep;
							}else{
								print $no_val.$sep.$gender.$sep.$no_val.$sep.$kinship.$sep.$user_id.$sep.$research_group_id.$sep;
								print $dir1.$sep;								
							}
							#Add read 1 file name
							if (defined $couple->{$lane}->{1} ) {
								if ($output ne ''){							
									print  OUT $couple->{$lane}->{1}.$sep;
								}else{
									print  $couple->{$lane}->{1}.$sep;
								}
							}
							#Add read 2 file name
							if (defined $couple->{$lane}->{2} ) {
								if ($output ne ''){		
								print OUT $couple->{$lane}->{2}.$sep;
								}else{
								print $couple->{$lane}->{2}.$sep;	
								}
							}
							#The name of the read file must contain the lane number and 
							#must be preceeded by an underscore an should contain nothing after
							#the number
							my $lane_num = $lane;
							#if present, remove the last '_'
							if ($lane_num =~ /\_$/){chop($lane_num);}
							#if not present, add an '_' at the start
							if ($lane_num !~ /^_/){$lane_num = "_".$lane_num;}
														
							#If theere is a single read file just use the sample name for the readfname
							if ( scalar(@lanes_used) == 1){
									$lane_num = "";
							}
							#ReadFName is obtained by adding the lane number to the sample name
							#If the user specified an output file
							if ($output ne ''){	
								print OUT $sample_dir.$lane_num.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";
							}#else she/he didn't specify the out file
							else{
								print $sample_dir.$lane_num.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";									
							}
							$numFile++;
							
						}
					}
					############NO MULTIPLE LANES
					############
					else{
						my $fq1_found = 0;
						my $fq2_found = 0;
						
						#Check if there are more than two files in the folder. Could mean the user used lanes
						if ( scalar(@files) > 2 ){
							print "WARNING: there are multiple files in $dir1. If multiple lanes fastq are present use -ml\n"
						}
						#Expects only two files inside the folder with _R1_ and _R2_
						foreach my $file (@files){

							#print "Using $file..\n";#DEBUGCODE

								#Search the files with _R1_ and R2
								if ( $file =~ /$r1_suff/ ){
									$fastq1 = $file;
									$fq1_found = 1;
								}elsif ( $file =~ /$r2_suff/ ){
									$fastq2 = $file;
									$fq2_found = 1;
								}						
						}
							#Error not found fq1
							die "ERROR: Cannot find fastq with extension $r1_suff. Please use --r1_suff to change how to search.\n"
								unless ($fq1_found);
							#Error not found fq1
							die "ERROR: Cannot find fastq with extension $r2_suff. Please use --r2_suff to change how to search.\n"
								unless ($fq2_found);

							
							#print "for folder $sample_dir used ".scalar(@lanes_used)."\n";
							if ($output ne ''){	
								print OUT $no_val.$sep.$gender.$sep.$no_val.$sep.$kinship.$sep.$user_id.$sep.$research_group_id.$sep;
								print OUT $dir1.$sep;
								print OUT $fastq1.$sep.$fastq2.$sep;
								print OUT $sample_dir.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";		
							}else{
								print $no_val.$sep.$gender.$sep.$no_val.$sep.$kinship.$sep.$user_id.$sep.$research_group_id.$sep;
								print $dir1.$sep;
								print $fastq1.$sep.$fastq2.$sep;
								print $sample_dir.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";		
							}							
						$numFile++;						
					}	
					push(@samples_used,$sample_dir);
				}
			}
			#Reset the flag to 1 so that the next sample is to analyze
			#unless it is set again to 0
			$analyze_sample = 1			
	}
	if (scalar(@samples_used) > 0){
		print "Used ".scalar(@samples_used)." samples \n";
	}else{
		print "WARNING: no sample has been found in $dir \n";
	}
	close(OUT);
}


=head2 write_sample_sheet_both_selected

 Title  : write_sample_sheet_both_selected
 Usage  : write_sample_sheet_both_selected( - dirname => 'the folder to check',
                      );
                      
		Using a file in input with information about the samples 
		generates a text file with the data that will be imported into the db 
              
              dir: source directory
              mode: motor or exome or exome_gr
							mult_lanes: can be 'lanes' or something other
							
							Use exome_gr if you know there are multiple samples inside a folder
							
							set the intern_fold variable if needed
										
 Function: 
					
 Returns : nothing

=cut
sub write_sample_sheet{

	my $sep = "\t";
	my $no_val = "-";
	

	#my $intern_fold = "";#"/Files/";

	#Put the names from the file in an array
	my @list = ();
	my @samples_used = ();
	my $sample_info_hash;
	my $gender_taken = 0;
	
	if (defined $fileList and $fileList ne 'none' and $fileList ne ''){
		if ( get_sample_info_hash(\$sample_info_hash,\@list,$fileList) == 1){
			$gender_taken = 1;
		}
		else{
			@list = list_to_array($fileList,'NO_NEW_LINE');
		}
	}
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs= readdir DIR;#All the files are in this list

	
	my $final_info = "";
	if ($mode eq "targeted"){
		$perchrom = 0;
		$final_info = "targeted".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		#$target_bed = "MotorPlex_v8_TuttiGeni_HPnormale_Covered.bed";
		unless ( $target_bed ne "" ) {$target_bed = $motor_target };
	}elsif ($mode eq "exome"){#Exome single sample
		#$perchrom = 1;
		$final_info = "exome".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $exome_target;}
	}elsif ($mode eq "exome_gr"){#Exomes
		$final_info = "exome".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $exome_target;}
	}elsif ($mode eq "rna"){#RNA
		$final_info = "rna".$sep.$perchrom.$sep.$convscores.$sep.$infreq;
		unless ( $target_bed ne "" ) {$target_bed = $rna_target;}
	}else{
			die "Need a mode to insert the parameters! [exome, motor, exome_gr, rna]\n";
	}

	
	#Start to write the file 
	if ($output ne ''){
		open (OUT,">$output") or die "ERROR: Cannot open $output. The program will exit..\n";
	}
	#print "There are ".scalar(@samples_dirs)." elements in $dir\n";#DEBUGCODE
	my $couple;
	#Get and print header
	open (HEAD,"<$head_file") or die "ERROR: Cannot open $head_file. The program will exit..\n";
	my $header = <HEAD>;
	close(HEAD);
	#Substitute commas in the header with the separator used here
	$header =~ s/,/$sep/g;
	if ($output ne ''){
		print OUT $header;
	}else{
		print $header;
	}
	#A control starts on each sample among all the samples        
	my $numFile= 1;
	my $fastq1 = "";
	my $fastq2 = "";
	
	#Flag used to check if the sample is to take or not
	my $get_sel_only = 0;
	my $analyze_sample = 1;

	
	#Go through the different samples directories
	foreach my $sample_dir (@samples_dirs){
			#print "Considering sample $sample_dir \n";#DEBUGCODE
			
			my $analysis_name = $sample_dir;			
			##########CONSIDER SAMPLELIST
			#If the file list is defined set the analyze_sample flag to NO
			if (defined $fileList and $fileList ne 'none'  and $fileList ne ''){
						$get_sel_only = 1;
						if ( not grep {/\b$sample_dir\b/} @list ){
							$analyze_sample = 0;
							#print "The sample $sample_dir will not be analyzed..\n";#DEBUGCODE
						}else{
							#print "The sample $sample_dir will be analyzed..\n";#DEBUGCODE
						}
			}
			
			#Get all the samples unless the user wants to get only those selected and this sample is ok
			unless ( $get_sel_only == 1 and $analyze_sample == 0){
			
				#If the folder name does not start with a dot and is a  directory
				if ( ($sample_dir !~ /^\..*/) and (-d $dir.$sample_dir)){

					#print "$sample_dir will be analyzed.\n";#DEBUGCODE
					my @lanes_used = ();
					my $kinship = "-";
					my $gender = "-";
					my $dob = "-";
					my $pob = "-";
					my $affected = '-';
					
					############SET KINSIP
					#Kinship is associated when the last letter is given
					#More than one proband are specified with the number after the _P
					#The groupname is picked from the proband sample name
					if ( $set_kinship){
						if ($sample_dir =~ /_P\d*$/ ) {
							$kinship = "P";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_P\d*$//;
						}elsif ($sample_dir =~ /_M$/) {
							$kinship = "M";
							$gender  ="F";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_M$//;
						}elsif ($sample_dir =~ /_F$/) {
							$kinship = "F";
							$gender  ="M";
							#pick the analysis name from the sample name removing the kinship
							$analysis_name = $sample_dir;
							$analysis_name =~ s/_F$//;
						}else{
							$kinship = "P";
							$gender  = "-";
							#pick the analysis name from the sample name 
							$analysis_name = $sample_dir;
						}
					}
					
					##############SET GENDER
					#If the gender has been put into 
					if ( defined $sample_info_hash->{$sample_dir}->{'gender'}){
						$gender  = $sample_info_hash->{$sample_dir}->{'gender'};
					}
					#If the dob has been put into 
					if ( defined $sample_info_hash->{$sample_dir}->{'dob'}){
						$dob  = $sample_info_hash->{$sample_dir}->{'dob'};
					}
					#If the pob has been put into 
					if ( defined $sample_info_hash->{$sample_dir}->{'pob'}){
						$pob  = $sample_info_hash->{$sample_dir}->{'pob'};
					}
					#If the affected has been put into 
					if ( defined $sample_info_hash->{$sample_dir}->{'affected'}){
						$affected  = $sample_info_hash->{$sample_dir}->{'affected'};
					}															
					
					############ADD SUFFIX TO GROUPNAME
					#If the group name is not specified than it is the same as the sample dir
					if ( ($mode eq 'targeted') or ($mode eq 'exome') or ($mode eq 'rna')){
						#If change_analysisname is set, then it is used. Otherwise add_to is verified 
						if ($change_analysisname ne "" ){
							$analysis_name = $change_analysisname;
						}else{
							#If I already changed the groupname in set_kinship
							if ($set_kinship){
								$analysis_name = $analysis_name.$add_to_analysisname;
							}else{
								$analysis_name = $sample_dir.$add_to_analysisname;
							}
						}
					}	

					#print "Entering $sample_dir...\n";
					my $dir1 = $dir."/".$sample_dir;
					#print $dir1."\t";
					$dir1 .= "/$intern_fold/";
					opendir DIR1, $dir1 or die "ERROR: cannot open dir $dir1";
					#my @files= readdir DIR1;#All the files are in this list
					my @files = grep { /\.($exts_to_search)$/ } readdir DIR1;
					die "ERROR: Cannot find files with extensions $exts_to_search in $dir1\n".
							"Maybe there is an internal folder? (intern_fold: $intern_fold)" unless @files > 0;
				
					#A control starts on each file of the folder        
							
					############MULTIPLE LANES
					############					
					if ( $mult_lanes ){
						my $ext_found = 0;
						foreach my $file (@files){
							#print "Using $file..\n";#DEBUGCODE
								$ext_found = 1;
								foreach my $lane (@lanes){
									#print "Checking if $lane is contained...\n";
									if ( $file =~ /$lane/ ){
										#print "$file contains $lane...\n";
										#Check if is the read 1 or two
										if ( $file =~ /$r1_suff/ ){
											$couple->{$lane}->{1} = $file;
											push(@lanes_used,$lane);
										}else{
											$couple->{$lane}->{2} = $file;
										}
									}
								}								
						}
						
						#Error not found fq
						die "ERROR: Cannot find fastq files ($ext_to_search) in $dir1. Exiting\n"
							unless ($ext_found);
													
						#Print both the files for each lane used
						#print "for folder $sample_dir used ".scalar(@lanes_used)."\n";
						foreach my $lane (@lanes_used){
							#gender,kinship
							if ($output ne ''){			 
								print OUT $gender.$sep.$kinship.$sep.$dob.$sep.$pob.$sep.$affected.$sep.$user_id.$sep.$research_group_id.$sep;
								print OUT $dir1.$sep;
							}else{
								print $gender.$sep.$kinship.$sep.$dob.$sep.$pob.$sep.$affected.$sep.$user_id.$sep.$research_group_id.$sep;
								print $dir1.$sep;								
							}
							#Add read 1 file name
							if (defined $couple->{$lane}->{1} ) {
								if ($output ne ''){							
									print  OUT $couple->{$lane}->{1}.$sep;
								}else{
									print  $couple->{$lane}->{1}.$sep;
								}
							}
							#Add read 2 file name
							if (defined $couple->{$lane}->{2} ) {
								if ($output ne ''){		
								print OUT $couple->{$lane}->{2}.$sep;
								}else{
								print $couple->{$lane}->{2}.$sep;	
								}
							}
							#The name of the read file must contain the lane number and 
							#must be preceeded by an underscore an should contain nothing after
							#the number
							my $lane_num = $lane;
							#if present, remove the last '_'
							if ($lane_num =~ /\_$/){chop($lane_num);}
							#if not present, add an '_' at the start
							if ($lane_num !~ /^_/){$lane_num = "_".$lane_num;}
														
							#If theere is a single read file just use the sample name for the readfname
							if ( scalar(@lanes_used) == 1){
									$lane_num = "";
							}
							#ReadFName is obtained by adding the lane number to the sample name
							#If the user specified an output file
							if ($output ne ''){	
								print OUT $sample_dir.$lane_num.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";
							}#else she/he didn't specify the out file
							else{
								print $sample_dir.$lane_num.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";									
							}
							$numFile++;
							
						}
					}
					############NO MULTIPLE LANES
					############
					else{
						my $fq1_found = 0;
						my $fq2_found = 0;
						
						#Check if there are more than two files in the folder. Could mean the user used lanes
						if ( scalar(@files) > 2 ){
							print "WARNING: there are multiple files in $dir1. If multiple lanes fastq are present use -ml\n"
						}
						#Expects only two files inside the folder with _R1_ and _R2_
						foreach my $file (@files){

							#print "Using $file..\n";#DEBUGCODE

								#Search the files with _R1_ and R2
								if ( $file =~ /$r1_suff/ ){
									$fastq1 = $file;
									$fq1_found = 1;
								}elsif ( $file =~ /$r2_suff/ ){
									$fastq2 = $file;
									$fq2_found = 1;
								}						
						}
							#Error not found fq1
							die "ERROR: Cannot find fastq with extension $r1_suff. Please use --r1_suff to change how to search.\n"
								unless ($fq1_found);
							#Error not found fq1
							die "ERROR: Cannot find fastq with extension $r2_suff. Please use --r2_suff to change how to search.\n"
								unless ($fq2_found);

							
							#print "for folder $sample_dir used ".scalar(@lanes_used)."\n";
							if ($output ne ''){	
								print OUT $gender.$sep.$kinship.$sep.$dob.$sep.$pob.$sep.$affected.$sep.$user_id.$sep.$research_group_id.$sep;
								print OUT $dir1.$sep;
								print OUT $fastq1.$sep.$fastq2.$sep;
								print OUT $sample_dir.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";		
							}else{
								print $gender.$sep.$kinship.$sep.$dob.$sep.$pob.$sep.$affected.$sep.$user_id.$sep.$research_group_id.$sep;
								print $dir1.$sep;
								print $fastq1.$sep.$fastq2.$sep;
								print $sample_dir.$sep.$sample_dir.$sep.$analysis_name.$sep.$reference.$sep.$target_bed.$sep.$final_info."\n";		
							}							
						$numFile++;						
					}	
					push(@samples_used,$sample_dir);
				}
			}
			#Reset the flag to 1 so that the next sample is to analyze
			#unless it is set again to 0
			$analyze_sample = 1			
	}
	if (scalar(@samples_used) > 0){
		print "Used ".scalar(@samples_used)." samples \n";
	}else{
		print "WARNING: no sample has been found in $dir \n";
	}
	close(OUT);
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
	"\t-h|--head_file: The header file is the file with the header of the sample sheet. It is name ss_head_file.txt and ".
									"located into the VarGenius/CONFIGURATION folder. Please use its full path.\n".
	"\t-d|--dir: Full path of the directory where the samples folders are located.\n".
	"\t-m|--mode: [exome|targeted|rna]\n".
	"\t-ml|--mult_lanes: use this parameter if into the folders you find different files for different samples. ".
										" The file names must contain _L001_,_L002_,_L003_,_L004_ and so on.\n".
  "\t-r1|--r1_suff: set this parameter to set how to find the fastq read 1. Default value is _R1_\n".
  "\t-r2|--r2_suff: set this parameter to set how to find the fastq read 2. Default value is _R2_\n".	
	"\t-fl|--fileList: if you want to use selected samples in the folder for this sample sheet, set this parameter ".
									"giving in input a text file with a list of sample names separated by newline.\n".
									"This file could also be used to give the gender of selected samples. Please separate gender and sample name with a TAB.
									(samplename	F). Allowed genders: F and M".
	"\t-lanes|--lanes_list: set this parameter to set how to find the lane names. Give it a list of lanes names ".
	" (e.g. lane1,lane2,lane3,lane4). Default value is L001,L002,L003,L004,L005,L007,L008\n".
  "\t-t|--target_bed: the target BED file name (just the name, not the path)\n".
  "\t-if|--intern_fold: if the fastq are contained into a subdirectory of the sample folder specify the subfolder name here.\n".
  "\t-pchr|--perchrom: just 0 or 1 if you want to do the analysis per-chromosome or not. Default: 1\n".
  "-r|--reference: reference genome identifier[hg19]\n".
  "\t-conv|--convscores: if scores must be converted to Phred [0 or 1]\n".
  "\t-uid|--user_id: identifier of the user that is running vargenius (a number that is choose by you);\n".
  "\t-rgid|--research_group_id: identifier of the research group that is running vargenius (a number that is choose by you);\n".
  "\t-agn|--add_to_analysisname: string to add to the analysis name. Default: both the sample name and analysis name are the sample folder name. \n".
  "\t-cgn|--change_analysisname: set and assign same analysis name to each file of the sample sheet.\n".
  "\t-infr|--infreq: if this analysis should be kept in count for the variant frequencies calculation [0|1]\n".
  "\t-o|--output: output sample sheet path and name\n".
  "\t-k|--set_kinship: This is a flag just to search for kinship at the end of the name. The name must end with one among: ".
										"P:proband,M:mother,F:father,A:affected non proband (brother, syster or cousin of proband),".
										"R:any relative inserted for any reason\n\n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "version" => \$VERSION,
           "h|head_file=s" => \$head_file,   #It's mandatory
           "d|dir=s" => \$dir,   #It's mandator
           "m|mode=s" => \$mode,   #It's mandatory
           "ml|mult_lanes" => \$mult_lanes,   
					 "r1|r1_suff=s" => \$r1_suff,
					 "r2|r2_suff=s" => \$r2_suff,  
					 "lanes|lanes_list=s" => \$lanes_list,
           "fl|fileList=s" => \$fileList,   
           "t|target_bed=s" => \$target_bed,   #It's mandatory
           "o|output=s" => \$output,
           "if|intern_fold=s" => \$intern_fold,
           "pchr|perchrom=s" => \$perchrom, 
           "r|reference=s" => \$reference,
           "conv|convscores=s" => \$convscores,
           "uid|user_id=s" => \$user_id,
           "rgid|research_group_id=s" => \$research_group_id,
           "aan|add_to_analysisname=s" => \$add_to_analysisname,
           "can|change_analysisname=s" => \$change_analysisname,
           "infr|infreq=s" => \$infreq,
           "k|set_kinship" => \$set_kinship,#Flag just to search for kinship at the end of the name
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




#############################UTILITIES

=head2 list_to_array

 Title   : list_to_array
 Usage   : list_to_array( listFile = file path
													new_line = take or not the new line
										);

 Function:  puts lines of a file inside an array and returns the array.
					You can use it taking or not the new line character '\n' at the end of eac line
					by using the parameter new_line (NO_NEW_LINE, to take or nothing or NO to
					not take.
					DEFAULT: will take
 Returns : an array with the lines

=cut
sub list_to_array{ 
	my $listFile = shift;
	my $new_line = shift;
	
	my @list = ();
	
	#you can do it using the newline or not
	if ($new_line eq 'NO_NEW_LINE'){
		open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
		while ( my $row = <LIST>){
				chomp($row);
				push(@list,$row);
		}
		close(LIST);
	}else{
		#Reading all in one with <FILE> command
		open (FILE,"<$listFile") or die "Cannot open $listFile\n";
		@list = <FILE>;
		close(FILE);
	}
	return @list;
}
