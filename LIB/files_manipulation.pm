
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
    

package LIB::files_manipulation;

## files_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of an files with many manipulation functions
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(  delete_file is_folder_empty tsv_2_xls copy_2_location 
									append_file_2_file append_str_2_file_if_path_notexist
									list_to_array extract_name );
}
									
use strict;
use warnings;
#FILE MANAGEMENT MODULES
use File::Path qw(make_path remove_tree ) ;
use File::Fetch; #TO use FTP download functions
use File::Copy; #Used to move files
use List::Util qw(first);#to get the index of a string in an array

use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;
use Archive::Extract;
use IO::Compress::Gzip qw(gzip $GzipError);
use Data::Dumper;#To print the hashes
use LIB::utilities qw(  extract_fasta_from_fasta );
use Storable;#To store data structures
use Excel::Writer::XLSX;#For conversion of tab sep in xls format

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(  try_exec_command  );
use LIB::std_lib qw( print_array correct_type);


#MOVE, DELETE

=head2 delete_file

 Title   : delete_file
 Usage   : delete_file( -filePath => 'a path to a file to delete'
			       );

 Function:  Delete a file in a given location. It can erase also more 
						than one file starting with the same name. This
            can be done by using the usual '*'. When it is found the 
            'glob' keyword will be used.
 Returns : Error code: 1 -> correctly deleted; -1 -> error with perl function; -2 file does not exist

=cut
sub delete_file{
  my $filePath = shift;
 
  my $retVal = 0;
  
	#With this first IF we can delete a set of file with the same name
	if ($filePath =~ /\*/){
		$retVal = 1;
		unlink glob $filePath or $retVal = -1;
	 }elsif ( -e $filePath ){
			if  ( unlink($filePath) == 1) {
				$retVal = 1;
				#deleted successfully
			}else{
			 #not deleted for some problems with unlink subroutine or the file does not exist anymore
				$retVal = -1;
			}
	}else{
	 #does not exist
		$retVal = -2;
	}
    
  return $retVal;
}

########################FILE READING
#=head2 list_to_array

 #Title   : list_to_array
 #Usage   : list_to_array( listFile = file path
													#new_line = take or not the new line
										#);

 #Function:  puts lines of a file inside an array and returns the array.
					#You can use it taking or not the new line character '\n' at the end of eac line
					#by using the parameter new_line (NO_NEW_LINE, to take or nothing or NO to
					#not take.
					#DEFAULT: will take
 #Returns : an array with the lines

#=cut
#sub list_to_array{ 
	#my $listFile = shift;
	#my $new_line = shift;
	
	#my @list = ();
	
	##you can do it using the newline or not
	#if ($new_line eq 'NO_NEW_LINE'){
		#open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
		#while ( my $row = <LIST>){
				#chomp($row);
				#push(@list,$row);
		#}
		#close(LIST);
	#}else{
		##Reading all in one with <FILE> command
		#open (FILE,"<$listFile") or die "Cannot open $listFile\n";
		#@list = <FILE>;
		#close(FILE);
	#}
	#return @list;
#}


########################FILE WRITING







  
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
  my $path = shift;
  my $string = shift;
  
  #print "String: $string will be appended in $path..\n";
  #Check if the field exists
  if (-e $path){
		#Get all the paths of the files in an array
		open(DATA,"<$path") or die "Could not open file $path";
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
			open(DATA,">>$path") or die "Couldn't open file $path";
			#print DATA "index: $last_el_ind $filename \n";#DEBUGCODE
			print DATA $string."\n";
			close(DATA);			
		}
	}else{
		open(DATA,">>$path") or die "Couldn't open file $path";
		print DATA $string."\n";
		close(DATA);
	}

}


############################CONVERSION OF FORMAT



=head2 copy_2_location
 Title   : copy_2_location
 Usage   : copy_2_location(  - hash_file -> the complete path of the file where the hash is saved
                      );

 Function:  Copies a file to some location given the ip
 
 Returns : the inputed output file is written

=cut
sub copy_2_location{
	my $infile = shift;
	my $dest_path = shift;
	
	my $command = "scp $infile $dest_path";
	try_exec_command($command) or die "Unable to execute command: $command\n";		
	
}
############################CHECK FUNCTIONS

=head2 is_folder_empty

 Title  : is_folder_empty
 Usage  : is_folder_empty( - dirname => 'the folder to check',
                      );

 Function: Checks if a folder is empty

 Returns : nothing

=cut
sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

######################OTHER
=head2 extract_name

 Title   : extract_name
 Usage   : extract_name( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: extract the name from a complete path of a file. Even the file name only with the extension
              0: the complete name with the extension
              1: the name only
              2: the first two words joined by a dot
              noext: just remove the last extension from the path
              no2ext: remove the last 2 extensions from the path
              gz: the name from the .gz
              targz: the name from the .tar.gz
              zip: the name from the .zip
              tar: the name from the .tar
              fqgz: the name before two extensions

 Returns : the name only

=cut
sub extract_name {
  my $filePath = shift;#Path to the file
  my $type = shift;#The type of file
  
  #Separate the path in pieces using slashes
  my @list = split("/",$filePath);
  my $complName = pop(@list);
  
  
  my @nameElements = split (/\./, $complName);
  my $name;
  if ($type eq "0"){ $name = $complName;}
  elsif ($type eq "1"){ $name = $nameElements[0];}
  elsif ($type eq "2"){ $name = $nameElements[0].'.'.$nameElements[1];}
  elsif ($type eq 'noext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq 'no2ext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq "gz"){ $complName =~ /(\S+).gz/;
                 $name= $1;##Name for uncompressed file
                 }
  elsif ($type eq "targz"){$complName =~ /(\S+).tar.gz/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "zip"){$complName =~ /(\S+).zip/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "tar"){$complName =~ /(\S+).tar/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "fqgz"){
		$complName =~ /(\S+)\.\S+\.\S+/;
		
                 $name= $1;##Name for uncompressed file
                 }
   else { die	"ERROR [$?]: $type is not a valid input extracting a name: ?\n";}
  return $name;
  
}

1;
