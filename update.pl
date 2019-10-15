#!/usr/bin/perl

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
    

use strict;
use warnings;
use Cwd;
use File::Copy;

my $program_name = "vargenius.pl";
my $configUserFile = "user_config.txt";
#Files included in VarGenius to be inserted during the installation
my $gene_tracks = "geneTrack1xy_sort.txt";
my $gene_panels_fold = "gene_panels";
my $canonSplSites = "hg19_canonSplSites.bed";

my $data_folder = "DATA";
my $repository_folder = "REPOSITORY";
my $config_folder = "CONFIGURATION";
my $log_folder = "LOG";
my $perl_libs_paths_f = "perl_libs_paths.txt";#Into config folder

print "This script will prepare $program_name to work in your directory.\n";


my $sure= getInput("Are you sure you want to enjoy this software?(y or n)",'^[YynN]$');

my @main_pls = ("vargenius.pl","programs_runner.pl");

if ( ($sure eq "y") or ($sure eq "Y")){ 
  
  #Asks to the user what is the directory where he wants to work
  my $workDir = '';    
  while ( !(-e $workDir) ){
    print "\nWrite an existing complete path where you want to play with $program_name (/home/username/analyses): ";
    $workDir= <STDIN>;
    chomp $workDir; #to cancel the return 
  }
  
  #print "working dir: $workDir\n";
  #I take the folder path of the program
  my $progDir = getcwd;
  chmod 0755,$program_name;
  
  #Writing a file with the paths to the current folder and the one
  open (FOLD, ">$workDir/folders.txt") or die "Cannot create file folders.txt. Exiting.. \n";
  print FOLD "$workDir $progDir";
  close(FOLD);

  #Moving the User configuration file in the working directory
  chmod 0755,$progDir."/$configUserFile";
  copy($progDir."/$config_folder/$configUserFile",$workDir."/$configUserFile")
		or die "Unable to copy ".$progDir."/$config_folder/$configUserFile in $workDir/$configUserFile\n";

  
  #Creating the DATA folder into the working directory
  my $main_data_folder = $workDir."/$data_folder";
	#Check if directory exists, otherwise it creates it
	unless(-d $main_data_folder){
		print "Creating folder $main_data_folder...\n";
		mkdir $main_data_folder or die "ERROR: can't create folder $main_data_folder. Check permissions. \n";
	}

	#Copying the sorted GENE TRACKS (downloaded from GenomeBrowser and rearranged) into the DATA folder
  chmod 0755,$progDir."/$configUserFile";
  print "Copying gene tracks (downloaded from GenomeBrowser and rearranged) for DepthOfCoverage into $workDir/$data_folder...\n";
  print "(See section 'Coverage' in the OUTPUT page of the user guide )\n";
  copy($progDir."/$repository_folder/$gene_tracks",$workDir."/$data_folder/")
		or die "Unable to copy ".$progDir."/$repository_folder/$gene_tracks in $workDir/$data_folder/\n";	
  copy($progDir."/$repository_folder/$gene_tracks.idx",$workDir."/$data_folder/")
		or die "Unable to copy ".$progDir."/$repository_folder/$gene_tracks in $workDir/$data_folder/\n";	
			
	#Copying the gene_panels folder into the DATA folder
	print "Copying gene panels into  $workDir/$data_folder...\n";
  (system "cp -r $progDir/$repository_folder/$gene_panels_fold $workDir/$data_folder" ) == 0 
		or die "Unable to copy ".$repository_folder."/$gene_panels_fold in $workDir/$data_folder\n";	

	#Copying CanonSplSites into the DATA folder
	print "Copying Canonic Splicing Sites for Annovar into $workDir/$data_folder...\n";
  (system "cp -r $progDir/$repository_folder/$canonSplSites $workDir/$data_folder" ) == 0 
		or die "Unable to copy $repository_folder/$canonSplSites in $workDir/$data_folder\n";
					
  #Creating the LOG folder into the working directory
  my $main_log_folder = $workDir."/$log_folder";
	#Check if directory exists, otherwise it creates it
	unless(-d $main_log_folder){
		print "Creating folder $main_log_folder...\n";
		mkdir $main_log_folder or die "ERROR: can't create folder $main_log_folder. Check permissions. \n";
	}
	
	#Now change the starting part of the .pl scripts in VarGenius accordingly to
	#system settings
	print "Updating the perl scripts with specific platform settings...\n";
	my $start_string = "####PLATFORM_SPECIFIC_SETTINGS_TERMINATED";
	my $perl_location = `which perl`;
	my $vargenius_lib = "use lib '".$progDir."/';";
	chop($perl_location);
	my $new_start_string = "#!".$perl_location."\n\n#$program_name path\n".$vargenius_lib."\n";
	
	#Open the file with libraries if it exists
	my $perl_libs_paths_f_p = "$progDir/$config_folder/$perl_libs_paths_f";
	my @perl_libs_paths = ();
	if ( -e $perl_libs_paths_f_p and ! (-z $perl_libs_paths_f_p) ){
		@perl_libs_paths = file_list_to_array($perl_libs_paths_f_p);
		foreach  my $perl_libs_path (@perl_libs_paths){ 
			$new_start_string .= "use lib '".$perl_libs_path."';\n";
		}
	}
	
	foreach my $main_pl (@main_pls){
		print "$main_pl...";
		substitute_start_of_file($progDir."/".$main_pl,$start_string,$new_start_string);			
	}

	print "Done! Now you can start from this folder!\n";	 
}else {die "Change directory and start again the install script.\n";}
  





=head2 getInput

 Title   : getInput
 Usage   : getInput( - sentence: a sentence that will printed in input to ask something to the user;
					 - regex: what the answer of the user have to respect
                               );

 Function: Takes in input a sentence and a regex. Asks to the user the sentence and controls its input with regex
 
 Returns : input given by the user

=cut
sub getInput{
  my $sentence = shift;
  my $regex = shift;
  my $input='';
		
  while (!($input =~ /$regex/) or ($input eq "")){  
    print $sentence." ";
    $input = <STDIN>;
    chomp $input;
    print "\n";
  }	
  return $input;
}


=head2 substitute_start_of_file

 Title   : substitute_start_of_file
 Usage   : substitute_start_of_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: substitutes the start of the file given in input.
					It searches for the string_to_find, and until that string
					writes only the string_to_ins in the new file. Once that the string
					is found it copies all the rest of the input file into the output
					
 Returns : nothing

=cut
sub substitute_start_of_file{
  my $file = shift;
  my $string_to_find = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	
	#Print the input string in the new file
	print TEMP $string_to_ins."\n";
	
	# print the lines before the change
	while( my $line = <FILE> ) {
    chop($line);
    last if $line eq $string_to_find; # chosen string to find found
   }

	print TEMP $string_to_find."\n";
	# print the rest of the lines
	while( my $line = <FILE> )  {
     print TEMP $line;
   }
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	unlink($temp_f);
}

=head2 file_list_to_array

 Title   : file_list_to_array
 Usage   : file_list_to_array(    );

 Function:  puts lines of a file inside an array and returns the array
 Returns : nothing

=cut
sub file_list_to_array{ 
	my $listFile = shift;
	
	my @list = ();
	
	open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
	while ( my $row = <LIST>){
			chomp($row);
			if (-e $row ){
				push(@list,$row);
			}else{
				die "ERROR in $listFile: $row does not exist.\n ";
			}
	}
	close(LIST);
	return @list;
}

