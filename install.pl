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
use Getopt::Long;


my $update = "";
my $program_config = "CONFIGURATION/program_config.txt";
my $variables_file = "CONFIGURATION/variables.txt";
				
my $program_name = "vargenius.pl";
my $configUserFile = "user_config.txt";

print "This script will prepare $program_name to work in your directory.\n";

#Get the confighash and pre-defined variables
my $cfg_hash;
checkConfigVariables($program_config,$variables_file,1);
configFile2Hash($program_config,\$cfg_hash);
my $data_folder = $cfg_hash->{'data_folder'};
my $lib_folder = $cfg_hash->{'lib_folder'};#"LIB";
my $config_folder = $cfg_hash->{'config_folder'};#"CONFIGURATION";
my $log_folder = $cfg_hash->{'log_folder'};#"LOG";
my $repository_folder = $cfg_hash->{'repository_folder'};#"REPOSITORY";
#Files included in VarGenius to be inserted during the installation
my $gene_tracks = $cfg_hash->{'DOC_genelist'};#"geneTrack1xy_sort.txt";
my $gene_panels_fold = $cfg_hash->{'gene_panels_fold'};#"gene_panels";
my $canonSplSites = $cfg_hash->{'annov_ext_canonsplsites'};#"hg19_canonSplSites.bed";
my $lastGATK3vers = "GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2";
my $perl_libs_paths_f = $cfg_hash->{'perl_libs_paths_f'};# "perl_libs_paths.txt";#Into config folder
my $install_dep_script =  $cfg_hash->{'install_dep_script'};# "install_all.sh";
my $folders_file = $cfg_hash->{'folders_file'};

my $sure= getInput("Are you sure you want to enjoy this software?(y or n)",'^[YynN]$');

my @main_pls = ("vargenius.pl","programs_runner.pl");

parse_command_line_args();

if ( ($sure eq "y") or ($sure eq "Y")){ 
  
  #Asks to the user what is the directory where he wants to work
  my $workDir = '';    
  while ( !(-e $workDir) ){
    print "\nWrite an existing complete path where you want to play with $program_name (/home/username/analyses): ";
    $workDir= <STDIN>;
    chomp $workDir; #to cancel the return 
  }

	#Asks to the user what is the directory where he wants to install the programs
	my $binDir = '';    
	while ( !(-e $binDir) ){
		print "\nWrite an existing complete path where you want to install the programs (/home/username/bin/): ";
		$binDir= <STDIN>;
		chomp $binDir; #to cancel the return 
	}
		  
  #print "working dir: $workDir\n";
  #I take the folder path of the program
  my $progDir = getcwd;
  chmod 0755,$program_name;
  
  #Writing a file with the paths to the current folder and the one
  open (FOLD, ">$workDir/$folders_file") or die "Cannot create file $folders_file. Exiting.. \n";
  print FOLD "$workDir $progDir";
  close(FOLD);
	
 #print "UPDATE: $update\n";
 
  #Set the file with libraries if it exists
	my $perl_libs_paths_f_p = "$progDir/$config_folder/$perl_libs_paths_f";
	
	#################INSTALL DEPENDENCIES
	##This part of code allows the automated installation of all dependencies
	if ( $update > 0 ){ 
		print "The program will only be updated. Dependencies will not be installed.\n";
	}else{
		##Copying the last version of GATK3 This is illegal!
		#print "Copying the last version of GATK3  into $binDir...\n";
		#(system "cp -r $progDir/$repository_folder/$lastGATK3vers $binDir" ) == 0 
			#or die "Unable to copy $progDir/$repository_folder/$lastGATK3vers into $binDir\n";
			
					
		#Now go to the $lib_folder folder
		chdir $binDir;
		#Install dependencies
		my $install_script = "$progDir/$lib_folder/$install_dep_script";
		print "\nInstalling dependencies with  $install_script\n";
		chmod 0755,"$install_script";
		system("$install_script $binDir/" ) == 1 or print "DONE..\n";
		#Go back into the program folder
		chdir $progDir;  
	
		#Put the Vcf.pm library into the path of external libraries
		append_str_2_file_if_path_notexist($perl_libs_paths_f_p,$binDir."/vcftools-master/src/perl/");	
		
	}
	#################INSTALL DEPENDENCIES
	
	
		
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

	#Create the references folder as a subfolder of the main data folder
	my $references_fold = $main_data_folder."/".$cfg_hash->{'references_fold'};
	unless(-d $references_fold){
		print "Creating folder $references_fold...\n";
		mkdir $references_fold or die "ERROR: can't create folder $references_fold. Check permissions. \n";
	 }
	$cfg_hash->{'annovar_db_f'} = $references_fold."/".$cfg_hash->{'annovar_db_f_def'};
	#Create the folders for annovar databases
	my $annovar_db_f = $cfg_hash->{'annovar_db_f'};
	unless(-d $annovar_db_f){
		print "Creating folder $annovar_db_f...\n";
		mkdir $annovar_db_f or die "ERROR: can't create folder $annovar_db_f. Check permissions. \n";
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

	##Copying CanonSplSites into the DATA folder
	#print "Copying Canonic Splicing Sites for Annovar into $workDir/$data_folder...\n";
  #(system "cp -r $progDir/$repository_folder/$canonSplSites $workDir/$data_folder" ) == 0 
		#or die "Unable to copy $repository_folder/$canonSplSites in $workDir/$data_folder\n";

	##Copying CanonSplSites into the DATA folder
	print "Copying Canonic Splicing Sites for Annovar into ".$cfg_hash->{'annovar_db_f'}."...\n";
  (system "cp -r $progDir/$repository_folder/$canonSplSites ".$cfg_hash->{'annovar_db_f'} ) == 0 
		or die "Unable to copy $repository_folder/$canonSplSites in ".$cfg_hash->{'annovar_db_f'}."\n";
		

					
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

	
	print "Done! Now you can work with VarGenius from this folder!\n";
	print "Only the first time that you run VarGenius you must: ";
	print "\n\t1. Create the tables of the database";
	print "\n\t2. Download the genes information ";
	print "\n\t3. Download the genome fasta and the datasets for GATK";	 
	print "\nPlease read carefully the user guide (TUTORIAL page) to do this\n";	
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


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.

	my $howToUse = "Use with: \nperl complete_install.pl \n\n".
	"\t-u|--update: Use this command if you just want to update VarGenius with the same dependencies you installed before\n";
 
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "u|update" => \$update );
	#print "UPDATE: $update\n";
	if ($update){$update=1;}else{$update=-1;}
	#print "UPDATE: $update\n";
  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

}



=head2 append_str_2_file_if_path_notexist

 Title   : append_str_2_file_if_path_notexist
 Usage   : append_str_2_file_if_path_notexist( - path -> the file where to write,
                                - string => the string to be written)

 Function: 
						will append the string at the file path and put a newline
						Both string and file list must be this format: 
							AnyInformationAlsoMultipleColumns\tpathToFile
						but only if the file does not exist or there is not the title
						
						NOTICE: The path MUST be always the last element in the row!

 Returns : nothing

=cut
sub append_str_2_file_if_path_notexist{
  my $file = shift;
  my $string = shift;
  
  #print "String: $string will be appended in $file..\n";
  #Check if the file exists
  if (-e $file){
		#Get all the paths of the files in an array
		open(DATA,"<$file") or die "Could not open file $file";
		my @pres = ();
		while (my $row = <DATA>){
			chomp($row);
			my @fields = split("\t",$row);
			#The path must be always the last element in the row!
			my $last_el_ind = scalar(@fields)-1;
			push(@pres,$fields[$last_el_ind]);
		}
		close(DATA);
		
		#print_array(\@pres);#DEBUGCODE
		#Insert the element if it does not exists
		#Get the file name from the string in input
		my @in_str_flds = split("\t",$string);
		my $last_el_ind = scalar(@in_str_flds)-1;
		my $filename = $in_str_flds[$last_el_ind];
		
		#Append to file only if not present
		if ( !(grep {/\Q$filename\E/} @pres) ){
			open(DATA,">>$file") or die "Couldn't open file $file";
			#print DATA "index: $last_el_ind $filename \n";#DEBUGCODE
			print DATA $string."\n";
			close(DATA);			
		}
	}else{
		open(DATA,">$file") or die "Couldn't open file $file";
		print DATA $string."\n";
		close(DATA);
	}
}



=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                             - configHash = the pointer to the hash to be filled
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
 Returns : nothing

=cut
sub configFile2Hash{
  my $configFilePath=shift;
  my ($configHash) = shift;

  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or die "ERROR: The file $configFilePath doesn't exists. The program will exit..\n";
  while (my $line = <configFile>){
		#This IF is useful if you want to put a description above in the text file. Delimit it with a set of hashtags
    if ($line =~ /#########/){
      $start = 1;
    }
   # if( ($line =~ /(\S+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /#/) ){
   if( ($line =~ /(\w+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /^#/) ){
      if ( $2 ne ''){
			  $$configHash->{$1} = $2;
				#annoPrint ("$1 = $2\n") ;#DEBUGCODE
			}
    #}elsif ( $line =~ /(\S+)\s*=/ ){
    }elsif ( $line =~ /(\w+)\s*=$/ ){
			delete($$configHash->{$1});
		}
  }
  close(configFile);
  #print Dumper\$configHash; #DEBUGCODE
}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configFile -> file with the user configuration
                              - variablesFile -> the path to a file with all variables written
                              - lineToCheck -> line in the variables file to be used
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.

 Returns : nothing

=cut
sub checkConfigVariables {
  my $configFile = shift;
  my $variablesFile = shift;
  my $lineToCheck = shift;

  my $hashToCheck;

  if (! open(VARF,"<$variablesFile")){ die "ERROR: Failure opening '$variablesFile'. Your program version is corrupted - $!";}
  if (! open(CFILE,"<$configFile")){ die "ERROR: Cannot find '$configFile' - Your program version is corrupted - $!";}

  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CFILE>){
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashToCheck->{$1} = "OK";
    }
  }
  close(CFILE);

  my @confVars = ();

  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $lines=0;
  #For each line of the variables file
  while (my $line = <VARF>){
    $line =~ s/\n//;#Remove \n in the end of the line

    #get the variables in the line
    my @variables = split (/;/,$line);

    $lines++;
    #For each of the variables in the variables file
    foreach my $var (@variables){
      if ($lines == $lineToCheck){
        #print "program Variable: $var - - value: ".$hashToCheck->{$var}."\n";#DEBUGCODE
        #put the variable inside an array for program config
        push (@confVars, $var);
        if( !(defined($hashToCheck->{$var})) ){

          die "ERROR: in $configFile variable $var is missing. Please check the file. Closing...\n ";
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      }
    }
  }

  #print_array(\@allVars);
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashToCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @confVars )){
       die "ERROR: Variable $key is in the config file $configFile and not in $variablesFile file. This is completely wrong. Program will not work...\n ";
     }
  }
  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}
