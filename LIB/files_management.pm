
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
    
package LIB::files_management;
## files_management.pm
#Permits the management of an files with many manipulation functions
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( read_and_print check_presence 
								compare_dimensions delete_directory 
								dl_and_extract download_file extract_archive
								extract_special_chars join_files_with_cat match_array_with_string
								my_extract_any_file my_extract_file_sys print_file check_presence
								extract_file_folder merge_columns extract_col_from_file 
								create_folder is_folder_empty file_list_to_array
								search_file_with_pattern file_not_present
								join_files_with_zcat save_hash load_hash extract_columns_from_file
								append_str_2_file file_num_rows invert_cols_position get_col_index
								delete_file insert_col_in_file_table delete_columns
								extract_colnum_from_file_linux file_name_wostrange_chars
								append_hash_to_file append_string_at_line_to_file 
								substitute_start_of_file check_URL delete_rows_containing
								shift_column compress_folder list_to_array insert_line_into_file
								separate_elements_on_col tar_folder_sys gzip_folder extract_name
								append_str_2_file_if_path_notexist list_to_hash
								twocolfile_to_hash join_vcf_files_with_cat vcf_insert_string_into_field
								vcf_filter_identical_calls overwrite_str_in_file_if_exists
								append_file_2_file tsv_2_xls separate_bed_auto_and_sex_chr
								get_only_first_on_col_from_BED get_rows_containing
								count_lines_file insert_col_in_file compid_intervals_2_bed
								copy_folder merge_rows_with_same_col split_bedfile);
}
use strict;
use warnings;
#FILE MANAGEMENT MODULES
use File::Path qw(make_path remove_tree ) ;
use File::Fetch; #TO use FTP download functions
use File::Copy; #Used to move files
use List::Util qw(first);#to get the index of a string in an array
use Vcf;#Use vcf_tools library

use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;
use Archive::Extract;
use Archive::Tar;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use IO::Compress::Gzip qw(gzip $GzipError);
use Data::Dumper;#To print the hashes
#use LIB::utilities qw(  extract_fasta_from_fasta );
use Storable;#To store data structures

#HTTP FTP CHECKS
use LWP::UserAgent;
use LWP::Simple;
use HTTP::Request;
use HTTP::Response;
use HTTP::Date;
use Net::FTP;

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(try_exec_command correct_type print_and_log get_zygosity);


#################GENERIC OPERATIONS ON FILES

=head2 overwrite_str_in_file_if_exists

 Title   : overwrite_str_in_file_if_exists
 Usage   : overwrite_str_in_file_if_exists( - path -> the file where to write,
                                - string => the string to be written)

 Function: 
						will append the string at the file path and put a newline
						Both string and file list must be this format: 
							AnyInformationAlsoMultipleColumns\tdescriptor
						If the descriptor exists the existing string will be overwrite
						
						NB: The descriptor MUST be always the last element in the row!

 Returns : nothing

=cut
sub overwrite_str_in_file_if_exists{
  my $file = shift;
  my $string = shift;
  
  #print "String: $string will be appended in $file..\n";
  #Check if the file exists
  if (-e $file){
	my $lines;
		#Get all the paths of the files in an array
		open(DATA,"<$file") or die "Could not open file $file";
		my @pres = ();
		while (my $row = <DATA>){
			chomp($row);
			my @fields = split("\t",$row);
			#The path must be always the last element in the row!
			my $last_el_ind = scalar(@fields)-1;
			$lines->{$fields[$last_el_ind]} = $row;
		}
		close(DATA);
		
		#print_array(\@pres);#DEBUGCODE
		#Insert the element if it does not exists
		#Get the file name from the string in input
		my @in_str_flds = split("\t",$string);
		my $last_el_ind = scalar(@in_str_flds)-1;
		my $filename = $in_str_flds[$last_el_ind];
		
		#Append to file only if not present
		if ( ! defined $lines->{$in_str_flds[$last_el_ind]} ){
			open(DATA,">>$file") or die "Couldn't open file $file";
			#print DATA "index: $last_el_ind $filename \n";#DEBUGCODE
			print DATA $string."\n";
			close(DATA);			
		}else{
			open(DATA,">$file") or die "Couldn't open file $file";
			#Go through the hash and print all the lines but not the existing one..
			foreach my $key (keys %{$lines}){
					if ( $key ne $in_str_flds[$last_el_ind]){
						print DATA $lines->{$key}."\t$key\n";
					}
			}
			#..update that
			print DATA $string."\n";
			close(DATA);
		}
	}else{
		open(DATA,">$file") or die "Couldn't open file $file";
		print DATA $string."\n";
		close(DATA);
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

=head2 copy_folder

 Title   : copy_folder
 Usage   : copy_folder( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function:Copies the content of a folder using the copy command on all the content of the folder
			and creating the folder before
 Returns : 

=cut
sub copy_folder{
	my $input_f = shift;
	my $destination =  shift;
	
	
	#Check if directory exists, otherwise it creates it
	my $folder = $destination."/".extract_name($input_f,0); 
	unless(-d $folder){
		print "Creating folder $folder...\n";
		mkdir $folder or die "ERROR: can't create folder $folder. Check permissions. \n";
	}
	print "Searching into $input_f/*...\n";
	my @files = glob "$input_f/*";
	foreach my $file (@files){
		print "Coying  $file into $folder..\n";
		copy($file,$folder) or print "ERROR: unable to copy $input_f/$file in $folder\n";	
	}
}


=head2 get_only_first_on_col_from_BED

 Title   : get_only_first_on_col_from_BED
 Usage   : get_only_first_on_col_from_BED( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: on a given column separates the fields and remains only the first element
		e.g. you give in input
		chr1 12 14 ens|ENST00000483767,mRNA|JX093079
		with: col2sep_ind = 4
		 it returns
		chr1 12 14 ens|ENST00000483767,
			NB: Column index given in input must start from 1!

 Returns : 

=cut
sub get_only_first_on_col_from_BED{
	
	my $bed = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	my $retval = 1;
	my $sep = ",";
	####################SEPARATE the rows based on one column with multiple elemnts
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	foreach my $line (<FILE>){
		#Jump initial description
		if ($line !~ /^chr/){
			next;
		}
		chop($line);
		my @fields = split("\t",$line);
		
		#Check if the number of columns is greater than the column to sep
		if ( $col2Sep_ind > scalar(@fields) ){$retval = -1;last;}
		
		#Get the column to sep
		my $col2Sep = $fields[$col2Sep_ind - 1];
		my @elems = split(/$sep/,$col2Sep);	
		my $newline = "";

		for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
			$newline .= $fields[$col]."\t";
		}
		$newline .= $elems[0]."\t";
		for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
			$newline .= $fields[$col]."\t";
		}			
		chop($newline);
		$newline .= "\n";
		
		#print $newline;
		print NEWFILE $newline;
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
	return $retval;
}


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
		
		
=head2 shift_column

 Title  : shift_column
 Usage  : shift_column ();
 
 Function: Shifts the column "col2shift" after the column "afterthis" into the file "intable".
					col2shift and afterthis must be the strings 
  
 Returns : the final output

=cut		
sub shift_column {
	my $intable = shift;
	my $col2shift = shift;
	my $afterthis = shift;
	my $log_file = shift;
	
	#Get the position of the ivf field to remove it later
	my $col2shift_ind = get_col_index($intable,$col2shift);	
	#Set the file names for the column which will be extracted from the file
	my $col_f = $intable.".temp_col2shift";
	
	#Extract the columns into the file
	print_and_log( "Extracting column $col2shift from $intable to $col_f\n",$log_file);#DEBUGCODE
	extract_col_from_file($intable,$col2shift,$col_f);
	
	#Set the array that will be the column to insert and insert it
	my @list = list_to_array($col_f,'NO_NEW_LINE');			
	#	push the header into the list
	unshift @list, $col2shift;
	
	#Choose the point where to insert and insert it

	my $afterthis_ind = get_col_index($intable,$afterthis);
	print_and_log( "Inserting column $col2shift from file $col_f into $intable after $afterthis [$afterthis_ind] \n",$log_file);
	print_and_log( "list array has ".scalar(@list)." elements and the first is ".$list[0]."\n",$log_file);#DEBUGCODE
	insert_col_in_file_table($intable,\@list,$afterthis,'i');
	#Delete the previous columns using the previously obtained index $col2shift_ind
	#If you put them before the previous position, now that position needs you sum 1
	my @cols_to_remove = ();
	if ( $afterthis_ind < $col2shift_ind){
		@cols_to_remove = (($col2shift_ind + 1));
	}else{
		@cols_to_remove = (($col2shift_ind));
	}
	print_and_log( "Deleting ".scalar(@cols_to_remove)." [".$cols_to_remove[0]."] columns from $intable\n",$log_file);#DEBUGCODE
	delete_columns($intable,\@cols_to_remove);				
		
				#################################################
				##Shift the interval frequency fields after the 1000genomes frequency
				##Get the position of the ivf field to remove it later
				#my $ivf_position = get_col_index($vargenius_out,$cfg_hash->{'intern_var_freq_fld'});
				
				##print_and_log( "$analysis_indication Shifting columns ".$cfg_hash->{'intern_var_freq_fld'}." and ".$cfg_hash->{'var_frequency_factors_fld'}." after ".$cfg_hash->{'last_1000g_field'}."\n",$log_file);#DEBUGCODE
				##Set the file names for the two columns which will be extracted from the file
				#my $ivf_col_f = $vargenius_out.".temp_ivf";
				#my $ffact_col_f = $vargenius_out.".temp_ffact";
				##Extract the columns into the two files
				##print_and_log( "Extracting column ".$cfg_hash->{'intern_var_freq_fld'}." from $vargenius_out to $ivf_col_f\n",$log_file);#DEBUGCODE
				##print_and_log( "Extracting column ".$cfg_hash->{'var_frequency_factors_fld'}." from $vargenius_out to $ffact_col_f\n",$log_file);#DEBUGCODE
				#extract_col_from_file($vargenius_out,$cfg_hash->{'intern_var_freq_fld'},$ivf_col_f);
				#extract_col_from_file($vargenius_out,$cfg_hash->{'var_frequency_factors_fld'},$ffact_col_f);
				##Set the two arrays that will be the columns to insert and insert them
				#my @ivf_list = list_to_array($ivf_col_f,'NO_NEW_LINE');
				#my @ffact_list = list_to_array($ffact_col_f,'NO_NEW_LINE');
				##	push the header into the lists
				#unshift @ivf_list, $cfg_hash->{'intern_var_freq_fld'};
				#unshift @ffact_list, $cfg_hash->{'var_frequency_factors_fld'};
				
				##Choose the point where to insert and insert it
				##print_and_log( "Inserting column ".$cfg_hash->{'var_frequency_factors_fld'}." from file $ivf_col_f into $vargenius_out \n",$log_file);
				#my $last1000g_info_ind = get_col_index($vargenius_out,$cfg_hash->{'last_1000g_field'});
				##print_and_log( "ivf_list array has ".scalar(@ivf_list)." elements and the first is ".$ivf_list[0]."\n",$log_file);#DEBUGCODE
				##print_and_log( "ffact_list array has ".scalar(@ffact_list)." elements and the first is ".$ffact_list[0]."\n",$log_file);#DEBUGCODE
				#insert_col_in_file_table($vargenius_out,\@ivf_list,$cfg_hash->{'last_1000g_field'},'i');
				#insert_col_in_file_table($vargenius_out,\@ffact_list,$cfg_hash->{'intern_var_freq_fld'},'i');
				##Delete the previous columns using the previously obtained index $ivf_position
				##If you put them before the previous position, now that position needs you sum 2
				#my @cols_to_remove_later = ();
				#if ( $last1000g_info_ind < $ivf_position){
					#@cols_to_remove_later = (($ivf_position+2),($ivf_position+3));
				#}else{
					#@cols_to_remove_later = (($ivf_position),($ivf_position+1));
				#}
				##print_and_log( "Deleting ".scalar(@cols_to_remove_later)." [".$cols_to_remove_later[0]." and ".$cols_to_remove_later[1]."] columns from $vargenius_out\n",$log_file);#DEBUGCODE
				#delete_columns($vargenius_out,\@cols_to_remove_later);
				##############################################################
}	



=head2 file_name_wostrange_chars

 Title   : file_name_wostrange_chars
 Usage   : file_name_wostrange_chars(  );

 Function:  this subroutine check if a filename with a specific extension
					does not have dots inside and if contains only chars, numbers or _ and -
            Needs that the file in input has an extension!
            
 Returns : 1 if check is positive, 0 otherwise

=cut
sub file_name_wostrange_chars{
	my $file_name = shift;
	my $extension = shift;
	
	my $retVal = 1;
	
 #The name has to be without dots inside. Only one dot can separate the extension from the name
  my $dots=0;
  my @line = split (//,$file_name);
  foreach my $char(@line){
      if ($char eq '.'){
          $dots++;
      }
  }
  if ($dots>1){
    die "Please change the file name: $file_name by removing all internal dots (not the one used for the extension)...\n";
		$retVal = 0;
  }
  
 
 #print "File name: ".$file_name."\n";
  my @name_parts = split(/\./,$file_name);

  #It cannot be longer than 50 chars
  if (@name_parts > 1){
		my $name = $name_parts[0];
    #Checks for permitted characters
    if ( length($name) > 50 or ($file_name !~ /^[A-za-z0-9\_\-]+\.($extension)$/i) ){
       die "$name is not a correct file name. Allowed characters [A-Z,a-z,0-9,_,-]. Allowed extensions [.$extension]. Max length: 50 chars\n";
			$retVal = 0;
    }
	}
	return $retVal;
}
 
 
 
=head2 extract_colnum_from_file_linux

 Title  : extract_colnum_from_file_linux
 Usage  : extract_colnum_from_file_linux( - filePath => 'the path to the file');

 Function: given a file and a name of a column (or the index), fetches all the values
					under the selected column
  
  Returns : a file with only one column

=cut		
sub extract_colnum_from_file_linux {
	my $file = shift;
	my $col_num = shift;
	my $retFile = shift;
	my $sep = shift;
	
	my $command = "cut -f".$col_num." -d '".$sep."' $file > $retFile";
	system($command) == 0 or die "Unable to execute command: $command\n";
}


=head2 read_and_print

 Title   : read_and_print
 Usage   : read_and_print( -file => 'a path to a file to print'
			       );

 Function: Reads a file using the more command

 Returns : nothing

=cut
sub read_and_print {
		my $file = shift;
		my $command = "more $file";
		try_exec_command($command) or die "Unable to execute command: $command\n";
}


=head2 tsv_2_xls
 Title   : tsv_2_xls
 Usage   : tsv_2_xls(  - hash_file -> the complete path of the file where the hash is saved
                      );

 Function:  Converts a file in TSV format to XLS using the Spreadsheet::WriteExcel
 PERL library
 
 Returns : the inputed output file is written

=cut
sub tsv_2_xls{
	my $in_tsv = shift;
	my $out_xls = shift;
	
	# Open the tab delimited file
	open (TABFILE, $in_tsv) or die "$in_tsv: $!";

	# Create a new Excel workbook
	my $workbook  = Excel::Writer::XLSX->new($out_xls);
	my $worksheet = $workbook->add_worksheet();

	# Row and column are zero indexed
	my $row = 0;

	while (<TABFILE>) {
			chomp;
			# Split on single tab
			my @Fld = split('\t', $_);

			my $col = 0;
			foreach my $token (@Fld) {
					$worksheet->write($row, $col, $token);
					#$worksheet->write_string($row, $col, $token);
					$col++;
			}
			$row++;
	}
	close(TABFILE);
}

=head2 file_list_to_array

 Title   : file_list_to_array
 Usage   : file_list_to_array(    );

 Function:  A list of file is written into the input file.
			All the paths are inserted into an array if the file exists
			
			
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


=head2 list_to_hash

 Title   : list_to_hash
 Usage   : list_to_hash(    );

 Function:  puts lines of a file inside an hash that is automatically modified
 Returns : nothing

=cut
sub list_to_hash{ 
	my $listFile = shift;
	my ($hash) = shift;
	
	
	open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
	while ( my $row = <LIST>){
		
			chop($row);
			#print "Inserting $row...\n";
			$$hash->{$row} = $row;
	}
	close(LIST);
}

=head2 twocolfile_to_hash

 Title   : twocolfile_to_hash
 Usage   : twocolfile_to_hash(    );

 Function:  given a file with two columns separated by tab 
					put the second value into an hash 
						where the key is the first value
 Returns : nothing

=cut
sub twocolfile_to_hash{ 
	my $listFile = shift;
	my ($hash) = shift;
	
	
	open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
	while ( my $row = <LIST>){
		
			chop($row);
			my @pieces = split("\t",$row);
			#print "Inserting $row...\n";
			$$hash->{$pieces[0]} = $pieces[1];
	}
	close(LIST);
}

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

=head2 search_file_with_pattern

 Title   : search_file_with_pattern
 Usage   : search_file_with_pattern( -pattern => 'a pattern to search'
													-folder => the folder where the file with pattern is to search
			       );

 Function:  Opens a directory and searches a filename containing the pattern
				given in input. Gives an error if more than one file is found. In that
				case the first found is given
 Returns : the file name identified or an empty string

=cut
sub search_file_with_pattern {
 my $pattern = shift;
 my $folder = shift;
 
 opendir DIR, $folder or die "ERROR [$!]: cannot open dir $folder\n";
 my @files= readdir DIR;
 closedir DIR;
 
 #Search in folder
 my $ret_file = "";
 my @results = grep {/$pattern/} @files;
  if ( scalar(@results) == 1 ){
		$ret_file = $results[0];
	}elsif ( scalar(@results) > 1 ){
			print "WARNING: too many files matching $pattern in $folder\n";
		}else{
			print "Can't find a file matching $pattern in $folder..\n";
			}
			
 return $ret_file;			
}


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
			# not deleted for some problems with unlink subroutine or the file does not exist anymore
				$retVal = -1;
			}
	}else{
	# does not exist
		$retVal = -2;
	}
    
  return $retVal;
}

=head2 delete_directory

 Title   : delete_directory
 Usage   : delete_directory( -filePath => 'a path to a directory to delete'
			       );perl  VarGenius/get_sample_sheet.pl --help

 Function:  Delete a directory with all its content in a given location
 Returns : nothing

=cut
sub delete_directory{
  my $dirPath = shift;
  
  remove_tree($dirPath);
}  


=head2 dl_and_extract

 Title   : dl_and_extract
 Usage   : dl_and_extract( -link = the link to the file on the web
                            - compressed = the file name compressed
                             - uncompressed = the file name uncompressed
                              - user = the username to access some ftp (optional)
                               - pass = the password to access some ftp (optional)
                               );

 Function:  takes in input a link to a compressed file and downloads it if the uncompressed o the compressed
            are not in the folder. The user and password informations need only in some cases. So you can
            give to the function these two last parameters too.

 Returns : nothing


=cut
sub dl_and_extract{
  my $link = shift;
  my $compressed = shift;
  my $uncompressed = shift;
  my $out_folder = shift;
  my $log_file = shift;
  
  my $user = shift;
  my $pass = shift;
  my $originalFileDim = shift;
  
 # print_and_log( "Checking File: ".$link."\n",$log_file);#DEBUGCODE
  #Download the file 
 unless (-e $uncompressed){#Doesn't exists the  file, proceed
		unless (-e $compressed) {#Doesn't exists the compressed, download it
			 #Download
			 #	print_and_log( "Downloading File: ".$link."\n",$log_file);#DEBUGCODE

       print "\n Downloading File: ".$link."\n";
        if (defined ($user) and defined ($pass)){
           download_file($link,$out_folder,$user,$pass);
        }else{
           download_file($link,$out_folder);
        }
        #Check if the dimensions are  equal to those taken during the link checking
        #If original is 0 it means that the check has been incorrect
        if ( defined $originalFileDim ){
					if ( $originalFileDim > 0){
						if (compare_dimensions($originalFileDim,$compressed) <= 0){
								die "Error: ".$compressed." has not been well downloaded. Delete it before to start again Annocript.\n";
								#delete_file($dbDataFolder."/".$compressed);
						}
					}
				}
		 }else{ print_and_log( "File ".$compressed." already downloaded...\n",$log_file);#DEBUGCODE 
			 print "File ".$compressed." already downloaded...\n"; }
		#Uncompression
			print_and_log( "\n Uncompress File: $compressed --> $uncompressed\n",$log_file);#DEBUGCODE
    #print "\n Uncompress File: ".$out_folder."/".$compressed."-->".$out_folder."/".$uncompressed."\n";
    my_extract_any_file($compressed,$out_folder);
 }else{
	 print_and_log( "File ".$uncompressed." already present...\n",$log_file);#DEBUGCODE
	 print "File ".$uncompressed." already present...\n";}#Exists the file 
}



=head2 download_file

 Title   : download_file
 Usage   : download_file( -fileAdd => is the path of the file
                          -folderName => is the folder where we have to download
                          -user => username froruniprot
                            -pass => password for uniprot
                               );

 Function:  Given in input a path to a file, it download it in the data folder. It assigns as name of the file the last
#element of the array composed by the split of the string path.
The function tries more than one download types. The first one is with the GET of LWP. If it does not work, Annocript tries again with LWP Get
and if this again doesn't then it uses the Fetch::File class that calls various methods.
If this kind of download still does not work and we are trying to download Uniprot databases, then other domains are used because it can happen
that one of them is not working. The access is tried on Uniprot.org, EBI.org, Expasy.org

 Returns : it gives 1 if download succed, 0 if it is not, 2 if succeed with a change of link. In this case the new_link is valued

=cut
sub download_file{
  my $fileAdd = shift;
  #print "file name: ".$fileAdd."\n";
  my $folderName = shift;
  my $user = shift;
  my $pass =shift;
  my $originalFileDim = shift;
  
  my $new_link = '';
  my $retVal = 0;

  print "Downloading: $fileAdd...";
  
  my $ua = LWP::UserAgent->new;
  my $req = HTTP::Request->new (GET => $fileAdd);
  
  if (defined($user) and defined ($pass)){
    $req->authorization_basic($user, $pass);
  }
  my $fileToDl = extract_name($fileAdd,"0");
  my $res = $ua->request($req, $folderName."/".$fileToDl);
	 print "checking result!\n";#DEBUGCODE
  if ($res->is_success) {
     print "..completed!\n";
     $retVal = 1; 
  }
  else {
	  #Here we add a specific control and re-download because sometimes for problems of networks
	  #the file is not downloadable but it becomes after some seconds.
     print "Some problem occurred while downloading $fileAdd. Trying again...\n";
	   my $res2 = $ua->request($req, $folderName."/".$fileToDl);
	   if ($res2->is_success) {
			print "..completed!\n";
      $retVal = 1;
      }
      #Commented because I do not need to download at this time. Uncomment if needed
      else{
     
        #This command will force the download with different methods
        if (ftp_fetch_file($fileAdd,$folderName,$user,$pass) == 1){
            print "..completed!\n";
            $retVal = 1;
        #If it does not work and we are searching Uniprot, then access different servers
        }#elsif (match_array_with_string($fileAdd, \@uni_domains) == 1){
          #$new_link  = download_from_diff_sources($fileAdd,$folderName,$user,$pass);
          #$retVal = 2;
        #}
      }
  }
  if ($retVal == 0){
    die "Unable to download file $fileAdd. Restart again later. Program will close...\n";
  }
  return $retVal, $new_link;
}


=head2 match_array_with_string 

 Title   : match_array_with_string
 Usage   : match_array_with_string( -string => the string to search in
                                    - @array => the array with elements to search
                               );

 Function: Searches if at least one element of an array is present in a string
 
 Returns : 1 if succeed, 0 otherwise

=cut
sub match_array_with_string {
  my $string = shift;
  my $array = shift;

  my $found = 0;
  
  foreach my $el (@$array){
    if ( $string =~ /\Q$el\E/){
      $found = 1;
    }
  }
  return $found;
}


=head2 compare_dimensions 

 Title   : compare_dimensions
 Usage   : compare_dimensions( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.
 
 Returns : 1 if succeed, -1 otherwise

=cut
sub compare_dimensions{
  my $originalDim = shift;
  my $filePath = shift;
  
  my $retVal = -1;
  
  if ( (-s $filePath) == ($originalDim/1024) ) {
      $retVal = 1;
  }
  
  return $retVal;
}

=head2 compress_folder
 Title   : compress_folder
 Usage   : compress_folder( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of a file and that of the final compressed
            and compresses it as a tar

 Returns : nothing
=cut
sub compress_folder{
	my $dir = shift;
	my $infile = shift;
	my $outfile = shift;
	
	# Create a new tar object:
	my $tar = Archive::Tar->new();
	
	#Go to the input dir path
	chdir $dir;
	
	#Compress
	$tar->create_archive( $outfile,COMPRESS_GZIP,glob $infile."/*");
}

=head2 gzip_folder
 Title   : gzip_folder
 Usage   : gzip_folder( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of a file and that of the final compressed
            and compresses it as a tar

 Returns : nothing
=cut
sub gzip_folder{
	my $dir = shift;
	my $infile = shift;
	my $outfile = shift;
	
	my $in = $dir."/".$infile;
	my $out = $dir."/".$outfile;
	
	my $zip = Archive::Zip->new();

  # Add a directory
  my $dir_member = $zip->addDirectory( $in);
 
  # Save the Zip file
  unless ( $zip->writeToFileNamed($out) == AZ_OK ) {
    die 'write error';
  } 
}

=head2 tar_folder_sys
 Title   : tar_folder_sys
 Usage   : tar_folder_sys( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of a file and that of the final compressed
            and compresses it as a tar with the linux command

 Returns : nothing
=cut
sub tar_folder_sys{
	my $dir = shift;
	my $infold = shift;
	my $outfold = shift;
	
	
	#Use the system command TAR to tar the archive
	my $command = "tar -czvf $dir/$outfold $dir/$infold";
	#print_and_log( "Executing: $command..\n",$log_file);
	print "Executing: $command..\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 extract_archive
 Title   : extract_archive
 Usage   : extract_archive( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen

 Returns : nothing
=cut
sub extract_archive{
  my $input = shift;
  my $outDir = shift;
	
  my $path = '';
  my $ae = Archive::Extract->new( archive => $input);
  $Archive::Extract::PREFER_BIN = 1;
  
  #Before to extract save a copy of the ZIP file
  copy($input,$input.".temp") or die "Cannot copy $input to $input.temp";

  my $resp = $ae->extract( to => $outDir );
  $path = File::Spec->catfile( $outDir, $ae->files->[0] );
  
  #If there are some problems delete the file and die
  if ( !$resp ){
     delete_file($path);
     die "ERROR [$?]: ".$ae->error;
  } 
  
  #Move the copy to the original name
  move($input.".temp", $input) or die "ERROR: Cannot move $input.temp to $input\n";

  
  return $path;
}

=head2 my_extract_any_file
 Title   : my_extract_any_file
 Usage   : my_extract_any_file( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen.
            It can uncompress .tar.gz calling the extract_archive subroutint
            .gz, .bz2 and .zip calling anyuncompress subroutine.
            Thus, it uncompress always using PERL modules.

 Returns : nothing
=cut
sub my_extract_any_file{
  my $input = shift;
  my $outDir = shift;
  
  print "\n Uncompressing $input...\n";
  my $command = '';
  my $outName = '';
  if ($input =~ /\.tar.gz$/){
    extract_archive($input,$outDir);
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        #$command = "gunzip -c $input > $outName";
        anyuncompress $input => $outName
          or die " gz uncompress failed: $AnyUncompressError\n";
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        anyuncompress $input => $outName
          or die "zip uncompress failed: $AnyUncompressError\n";
      }elsif ($input =~ /\.bz2$/){
        $outName = extract_name($input,"bz2");
        $outName = $outDir."/".$outName;
        anyuncompress $input => $outName
          or die "bz2 uncompress failed: $AnyUncompressError\n";
        }  
  #print "Uncompression of ".$input." finished\n";
  print "...completed!\n";
}

=head2 my_extract_file_sys
 Title   : my_extract_file_sys
 Usage   : my_extract_file_sys( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen
                It extracts by using system calls.

 Returns : nothing
=cut
sub my_extract_file_sys{
  my $input = shift;
  my $outDir = shift;
  
  my $command = '';
  my $outName = '';
  
  my $removeText = 'Please go in the folder and remove all the files coming out from this erroneous execution and try again.'.
  'If you get again an error please consider to unzip manually the files and leave them inside the folder. Then restart Annocript.';
  if ($input =~ /\.tar.gz$/){
      $command = "tar -zxf $input -C $outDir";
      try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText \n" ;
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }else{
        extract_archive($input,$outDir);
      } 
      
  printf "Decompression of ".$input." finished\n";
}


=head2 ftp_fetch_file

 Title   : ftp_fetch_file
 Usage   : ftp_fetch_file( -fileAdd => is the path of the file
                          -folderName => is the folder where we have to download
                          -user => username froruniprot
                            -pass => password for uniprot
                               );

 Function:  Given in input a path to a file, it download it in a specific folder. 

 Returns : nothing

=cut
sub ftp_fetch_file{
  my $fileAdd = shift;
  my $folderName = shift;
  my $user = shift;
  my $pass =shift;

  my $result = 1;
  #print "Trying download with ftp_fetch_file: ".$fileAdd."\n";
  
  my $maxTimes = 5;#MAx number of times to execute the command
  my $success = -1;#If it is success finish trying
  my $timesCount = 0;#Counter of faults
   
  #attempts maxTimes to access and download the file  
  while ( ($success == -1) and ($timesCount < $maxTimes) ){
   ## build a File::Fetch object ###
   my $ff = File::Fetch->new(uri => $fileAdd);
   $File::Fetch::WARN = 0;
      # fetch the uri to folderName ###
      if ( $ff->fetch(to =>$folderName) ){
        $success = 1;
      }
      else{
        print "ERROR [$?]: an error occurred while fetching file $fileAdd. Trying again...\n"; 
        $timesCount++;
      }
  }
  #print $ff->error();
  return $success;
}

=head2 extract_file_folder

 Title   : extract_file_folder
 Usage   : extract_file_folder( -filePath => 'complete path of the file',
			       );

 Function: extracts the folder of a file from its complete path
              

 Returns : the name only

=cut
sub extract_file_folder {
  my $filePath = shift;#Path to the file
  
  #Separate the path in pieces using slashes
  my @list = split("/",$filePath);
  pop(@list);
  my $folder = join("/",@list);
  return $folder;  
}



=head2 file_num_rows

 Title   : file_num_rows
 Usage   : file_num_rows( - path1 -> the file to be used)

 Function: returns number of rows of a file
 
 Returns : nothing

=cut
sub file_num_rows{
  my $path1 = shift;

 my $row_num = 0;
  open(DATA1, "<$path1") or die "Couldn't open file $path1";
    
  while(<DATA1>){
   $row_num++;
  }
  close(DATA1);
return $row_num;
}


=head2 invert_cols_position

 Title   : invert_cols_position
 Usage   : invert_cols_position( - file -> the file to write,
                                - colA => the column index A
                                - colB => the column index B 
                                )

 Function: Inverts the position of the column A and column B.
					If colA=10 and colB=11, The column 11 will be the column 10 and
					viceversa.
						 
 Returns : the same file is modified

=cut
sub invert_cols_position{
  my $file = shift;
  my $colA = shift;
  my $colB = shift;
  
  my $tempFile = $file."_TEMP";
  my $command = 'awk -F"\t" '."'".'BEGIN{OFS=FS;}{t=$i;$i=$j;$j=t;}1'."'"." i=$colA j=$colB $file > $tempFile";
	print "Executing command: $command\n";
	system ($command) == 0 or die "Unable to execute command: $command\n";
					
	move($tempFile,$file) or print "ERROR: unable to move $tempFile in $file\n";
	delete_file($tempFile);
}

=head2 append_file_2_file

 Title   : append_file_2_file
 Usage   : append_file_2_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append the file on the path1 to the one on the path2
 
 Returns : nothing

=cut
sub append_file_2_file{
  my $path1 = shift;
  my $path2 = shift;
  my $noheader = shift;
  
  open(DATA1, "<$path1") or die "Couldn't open file $path1";
  open(DATA2,">>$path2") or die "Couldn't open file $path2";
  
  if (! (defined $noheader) ){$noheader = "YES";}
  my $numline = 1;
  while(<DATA1>){
	if ( $noheader eq "NOHEAD"){
		if ($numline > 1){
			print DATA2 "$_";
		}
	}else{
		print DATA2 "$_";
	}
	$numline++;
  }
  close(DATA1);
  close(DATA2);
}


=head2 compid_intervals_2_bed

 Title   : compid_intervals_2_bed
 Usage   : compid_intervals_2_bed( - path -> the file where to write,
                                - string => the string to be written)

 Function: # Given the output from ExomeDepth and XHMM 
#It extracts a BED structure
# From  chr5:175391955-175477735_DUP
#chr5    175391955   175477735 chr5:175391955-175477735_DUP
# ....
 
 Returns : nothing

=cut
sub compid_intervals_2_bed {
	my $input_f = shift;
	my $out_bed = shift;
		
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	#Generate a new bed file
	open(NEWFILE,">$out_bed");
				
	#Read all the file
	foreach my $line (<FILE>){
		chop($line);
		#First separate chr
		my @fields1 = split(":",$line);
		my $chr = $fields1[0];
		my @fields2 = split("-",$fields1[1]);
		my $start = $fields2[0];
		my @fields3 = split("_",$fields2[1]);
		my $end = $fields3[0];	
		
		print NEWFILE "$chr\t$start\t$end\t$line\n";
	}
	close(NEWFILE);	
	close(FILE);
}
  
=head2 append_str_2_file

 Title   : append_str_2_file
 Usage   : append_str_2_file( - path -> the file where to write,
                                - string => the string to be written)

 Function: 
						will append the string at the file path and put a newline.
					
 
 Returns : nothing

=cut
sub append_str_2_file{
  my $path = shift;
  my $string = shift;
  
  print "String: $string will be appended in $path..\n";
  
  #Check if the field exists
	open(DATA,">>$path") or die "Couldn't open file $path";
	print DATA $string."\n";
	close(DATA);
}


=head2 separate_bed_auto_and_sex_chr

 Title   : separate_bed_auto_and_sex_chr
 Usage   : separate_bed_auto_and_sex_chr( - path -> the file where to write,
                                - string => the string to be written)

 Function:  separates a BED with intervals in two files: autosomes and sex chromosomes
			#I am using this function for use in CNV calling with ExomeDepth and XHMM
			#If the description of the BED is in the header ExomeDepth will fail
 
 Returns : nothing

=cut
sub separate_bed_auto_and_sex_chr{
	
	my $bed = shift;
	my $sexbed = shift;
	my $autobed = shift;
	
	#SEX Chromosomes
	#Put chrX
	my $command = "grep '^chrX' $bed > $sexbed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	#append chrY
	$command = "grep '^chrY' $bed >> $sexbed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#Autosomes
	#Remove chrX
	$command = "grep -v '^chrX' $bed > $autobed.temp";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";

	#Remove chrY
	$command = "grep -v '^chrY' $autobed.temp > $autobed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#Remove descriptions
	$command = "grep  '^chr' $autobed > $autobed.temp";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
		
	#Move the copy to the original name
	move($autobed.".temp", $autobed) or die "ERROR: Cannot move $autobed.temp to $autobed\n";
	 
}



=head2 count_lines_file

 Title   : count_lines_file
 Usage   : count_lines_file( - path1 -> the file to be copied,
                             - path2 => the file to be written)

 Function: counts lines in a file
 
 Returns : nothing

=cut
sub count_lines_file{
	my $file = shift;
	
	open(FILE,"<$file") or die "ERROR: Cannot open $file\n";
	my @lines = <FILE>;
	close(FILE);

	return scalar( @lines );
}

=split_bedfile

Title   : split bedfile
 Usage   : split bed file(  - bed_file -> the entire file for chr1
                      );
Function:  this subroutine splits a BED file in two bed files
 
 Returns : nothing

=cut
sub split_bedfile {
	my $filein = shift;
	my $output_1 = shift;
	my $output_2 = shift;

	my @array_chr1 = () ;
	open(FILEIN, "<", $filein);
	while(<FILEIN>) {
		push @array_chr1, $_;
		chomp;
	}
	close(FILEIN);

	open (OUTPUT1,">$output_1") or die "ERROR[?]: Cannot open $output_1\n";
	open (OUTPUT2,">$output_2") or die "ERROR[?]: Cannot open $output_2\n";
	#splits array in half and print in two bed files
	for (my $i = 0; $i < (scalar(@array_chr1)/2); $i++) {
		print OUTPUT1 $array_chr1[$i];
	}
	for (my $z = (scalar(@array_chr1)/2)+1; $z < (scalar(@array_chr1)-1); $z++){
		print OUTPUT2 $array_chr1[$z];
	}
	close(OUTPUT1);
	close(OUTPUT2);
}

=head2 append_hash_to_file

 Title   : append_hash_to_file
 Usage   : append_hash_to_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append the hash in input on the file
 
 Returns : nothing

=cut
sub append_hash_to_file{
  my $file = shift;
  my $hash = shift;
  
  open(FILE, ">>$file") or die "Couldn't open file $file";
  
  foreach my $key (keys %$hash){
		print FILE $key." = ".$hash->{$key}."\n";
	}
  close(FILE);
}

=head2 append_string_at_line_to_file

 Title   : append_string_at_line_to_file
 Usage   : append_string_at_line_to_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append a string after a specific line (detected with a string)
					into a file. 
 
 Returns : nothing

=cut
sub append_string_at_line_to_file{
  my $file = shift;
  my $string_to_find = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	
	# print the lines before the change
	while( my $line = <FILE> ) {
    chop($line);
    print TEMP $line."\n";
    last if $line eq $string_to_find; # chosen string to find found
   }
   #Print the input line
	print TEMP $string_to_ins."\n";
	# print the rest of the lines
	while( my $line = <FILE> )  {
     print TEMP $line."\n";
   }
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
}

=head2 insert_line_into_file

 Title   : insert_line_into_file
 Usage   : insert_line_into_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: inserts an input line into a file given an index.
						The line is inserted at the index
					
 Returns : nothing

=cut
sub insert_line_into_file{
  my $file = shift;
  my $string_to_ins = shift;
  my $index = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	

	my $num_row = 0;
	# print the lines before the change
	while( my $line = <FILE> ) {
		if ( $num_row == $index){
			#Print the input string in the new file
			print TEMP $string_to_ins;		
		}
		#The current line is print after the new one
    print TEMP $line;
    $num_row++;
   }
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
	close(TEMP);
	close(FILE);
}



###SEPARATE the rows based on one column with multiple elemnts
##Example:
## chr1 9292 299923 GeneA,GeneB,GeneC
#becomes:
## chr1 9292 299923 GeneA
## chr1 9292 299923 GeneB
## chr1 9292 299923 GeneC
sub separate_elements_on_col{
	
	my $bed = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	####################SEPARATE the rows based on one column with multiple elemnts
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	#print "Executing gene separation on file $bed to fie $outname..\n";#DEBUGCODE
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		#print $line;#DEBUGCODE
		#Get the column to sep
		my $col2Sep = $fields[$col2Sep_ind - 1];
		my @elems = split(",",$col2Sep);	
		my $newline = "";
		foreach my $elem (@elems){
			#if (scalar(@elems) > 1){print "$elem \t";}
			for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
				$newline .= $fields[$col]."\t";
			}
			$newline .= $elem."\t";
			for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
				$newline .= $fields[$col]."\t";
			}			
			chop($newline);
			$newline .= "\n";
		}
		print NEWFILE $newline;
		#if (scalar(@elems) > 1){print "\n";}
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
}


=head2 merge_rows_with_same_col

 Title   : merge_rows_with_same_col
 Usage   : merge_rows_with_same_col( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: given a file and an index of a column, merges the rows where the first 
			column is the same and separates with a comma the "index" columns
			The 	field2merge is always added with comma at the end of line	
 Returns : nothing

=cut
sub merge_rows_with_same_col{
	my $input = shift;
	my $outname = shift;
	my $fields = shift;
	
	my @fields = split(",",$fields);
	my $field4comp = $fields[0];
	my $field2merge = $fields[1];
	
	open (FILE,"<$input") or die "Cannot open $input\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	#print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $compid1 = $fields1[$field4comp];
		my $compid2 = $fields2[$field4comp];

		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $lines[$i];
		while( $compid1 eq $compid2 and ($i < ($totlines-2)) ){
			
			#Add the column at the end of lne
			$row_2_write .= ",".$fields2[$field2merge];
			
			#Go ahed with next line
			$i++;				
			#print "Line $i: $compid1 and $compid2 is same information.\n";
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$compid2 = $fields2[$field4comp];
			

		}
		$row_2_write .= "\n";
		print NEWFILE $row_2_write
		
	}
	close(FILE);
	close(NEWFILE);
}

#=head2 separate_elements_on_col

 #Title   : separate_elements_on_col
 #Usage   : separate_elements_on_col( - path1 -> the file to be copied,
                                #- path2 => the file to be written)

 #Function: given a file and an index of a column, separates the elements
						#separated with comma and writes a line for each element
					
 #Returns : nothing

#=cut
#sub separate_elements_on_col{
	#my $file = shift;
	#my $outname = shift;
	#my $col2Sep_ind = shift;
	
	#print  "Separating $file into $outname\n";
	#####################SEPARATE elements separated with comma
	#open (FILE,"<".$file) or die "Cannot open $file\n";
	#open(NEWFILE,">".$outname) or die "Cannot open $outname\n";
	#print  "Column $col2Sep_ind\n";
	
	#while (my $line = <FILE>){
		#print  "Column 2 $col2Sep_ind: $line\n";
		#chop($line);
		#my @fields = split("\t",$line);
		
		##If the field contains elements separated by comma
		#if ($fields[$col2Sep_ind] =~ /,/){
			
			##Get index of the column wanted
			#my @elements = split(",",$fields[$col2Sep_ind]);	
							
			#foreach my $element  (@elements){
				#my $newline = "";
				#for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
					#$newline .= $fields[$col]."\t";
				#}
				#$newline .= $element."\t";
				#for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
					#$newline .= $fields[$col]."\t";
				#}			
				#chop($newline);
				#$newline .= "\n";
				#print NEWFILE $newline;
			#}			
		#}else{
				#print NEWFILE $line."\n";
		#}
	#}
	#print "Exiting\n";
	#close(FILE);
	#close(NEWFILE);	
#}


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
	delete_file($temp_f);
}
    
    

=head2 print_file

 Title   : print_file
 Usage   : print_file( - path1 -> the file to be print
                        )

 Function: will print the file on the screen
 
 Returns : nothing

=cut
sub print_file{
  my $path1 = shift;
  
  open(DATA1, "<$path1") or die "Couldn't open file $path1";
    
  while(<DATA1>){
   print "$_";
  }
  close(DATA1);
}


=head2 check_presence
 Title  : check_presence
 Usage  : check_presence(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : 1 if the file is ok, otherwise 0

=cut
sub check_presence{
 my $fileToCheck = shift;
 
 my $retVal = 1;
 #print $fileToCheck." \n";
 if ( $fileToCheck ne ''){
   if(-z "$fileToCheck"){
    $retVal = 0;
    print "$fileToCheck is empty...!\n";
   }
   if(!-e "$fileToCheck"){
    $retVal = 0;
    print " $fileToCheck does not exists..\n";
   }
 }else{
	 print ("The file to be checked is an empty string. Please check the code...\n");
	 $retVal = 0;
	}
 
 return $retVal;
}

=head2 file_not_present
 Title  : file_not_present
 Usage  : file_not_present(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : 1 if the filename is an empty string
					 2 if the file is not present
					 3 if the file is present but the size is 0

=cut
sub file_not_present{
 my $fileToCheck = shift;
 
 my $retVal = 0;
 #print $fileToCheck." \n";
 if ( $fileToCheck ne ''){
   if(-z $fileToCheck){
    $retVal = 1;
    print "$fileToCheck is present but empty...!\n";
   }
   #-e could return false if the directory is not readable by the program giving a misleading behaviour
	
   if(! (-e $fileToCheck) ){
    $retVal = 2;
		#print "$fileToCheck does not exists..\n";
   }
 }else{
	 print "The file to be checked is an empty string. Please check the code...\n";
	 $retVal = 3;
	}
 
 return $retVal;
}

=head2 vcf_filter_identical_calls

 Title   : vcf_filter_identical_calls
 Usage   : vcf_filter_identical_calls( -database => 'name of the database,
                               );

 Function:
		#takes in input a VCF that is supposed to be created
		#joining different programs calls hence there could be
		#many identical calls. This subroutine removes all the 
		#calls from the 2nd..ith program where all samples have the same
		#genotype of the first program called.
		#The 2nd...ith programs are recognizable from the VC_ALGORITHM 
		#field of the INFO colum
					
		#NB The VCF file must be sorted!	
 Returns : 

=cut			
sub vcf_filter_identical_calls {
	my $cfg_hash = shift;
	my $vcf_file = shift;
	my $out_vcf = shift;
	my $log_file = shift;

	my $gt_str = $cfg_hash->{'gt_str'};

	#Field names and indexes
	my $chr_str = $cfg_hash->{'db_vcf_chr_str'};
	my $chr_ind = $cfg_hash->{'vcf_chr_ind'};
	my $pos_str = $cfg_hash->{'db_vcf_pos_str'};
	my $pos_ind = $cfg_hash->{'vcf_pos_ind'};
	my $id_str = $cfg_hash->{'db_vcf_id_str'};
	my $id_ind = $cfg_hash->{'vcf_id_ind'};
	my $ref_str = $cfg_hash->{'db_vcf_ref_str'};
	my $ref_ind = $cfg_hash->{'vcf_ref_ind'};
	my $alt_str = $cfg_hash->{'db_vcf_alt_str'};
	my $alt_ind = $cfg_hash->{'vcf_alt_ind'};
	my $qual_str = $cfg_hash->{'db_vcf_qual_str'};
	my $qual_ind = $cfg_hash->{'vcf_qual_ind'};
	my $filter_str = $cfg_hash->{'db_vcf_filter_str'};
	my $filter_ind = $cfg_hash->{'vcf_filter_ind'};
	my $info_str = $cfg_hash->{'db_vcf_info_str'};
	my $info_ind = $cfg_hash->{'vcf_info_ind'};
	my $genot_ind = $cfg_hash->{'vcf_genot_ind'};
	my $prob_gen_ind = $cfg_hash->{'vcf_prob_gen_ind'};
	my $comp_id_str = $cfg_hash->{'db_vcf_comp_id_str'};
	#Separators
	my $sep = ",";
	my $tab = "\t";
	my $vcf_info_sep = $cfg_hash->{'vcf_info_sep'};
	my $vcf_gen_sep = $cfg_hash->{'vcf_gen_sep'};
	my $semicolon = ";";
	
	#A little hash to keep the number of saved calls from each algorithm
	my $saved_calls;
	my $identical_calls = 0;
	
	  #Open the ouput file that will be write
   open (OUT,">$out_vcf") or die "ERROR: Cannot open $out_vcf\n";
   
  #First print the header of the VCF 
  open (VCF, "<$vcf_file") or die "ERROR: Cannot open $vcf_file\n";
  while (my $row = <VCF>){
		if ($row =~ /^#/){
			print OUT $row;
		}else{
			last;
		}
	}
  close(VCF);

  #print_and_log( "Opening VCF file: $vcf_file\n",$log_file);#DEBUGCODE
  #Loads the VCF
	my $vcf = Vcf->new(file=>$vcf_file);
	#Reads the header
	$vcf->parse_header();
	
	#Get the list of samples for the VCF
	my (@samples) = $vcf->get_samples();
	
	#print "There are ".scalar(@samples)." samples in this VCF: ".$samples[0]."\n";#DEBUGCODE
	#Get all VCF lines into an array
	my @lines = ();
	while (my $row = $vcf->next_data_array()){
		push(@lines,$row);
	}
	my $totlines = scalar(@lines);

	print_and_log("Analyzing all calls in $vcf_file... ",$log_file);
	
	#Start the process going through all the lines of the input VCF file 
	#This loop gets couples of consecutive lines and matches the compids (Chr_pos_ref_alt)
	#for equal calls from different algorithms and if they do not call the same genotype the line is kept,
	#otherwise removed
	for (my $line = 0; $line < ($totlines-1); $line++){
		my $row1 = $lines[$line];
		my $row2 = $lines[$line+1];
		#Get compids
		my $compid1 = $$row1[$chr_ind]."_".$$row1[$pos_ind]."_".$$row1[$ref_ind]."_".$$row1[$alt_ind];
		my $compid2 = $$row2[$chr_ind]."_".$$row2[$pos_ind]."_".$$row2[$ref_ind]."_".$$row2[$alt_ind];
		
		#####Take all the consecutive equal calls#####
		#Save into an hash (row_hash) all the rows containing equal calls
		my $eq_rows = 0;
		my $row_hash;
		$row_hash->{$eq_rows} = $row1;
		while( $compid1 eq $compid2 and $eq_rows < ($totlines-1) ){
			#print "On line $line : $compid1 eq $compid2\n";#DEBUGCODE
			$eq_rows++;
			$row_hash->{$eq_rows} = $row2;
			$line++;#Go ahed with line
			$row2 = $lines[$line+1];
			$compid2 = $$row2[$chr_ind]."_".$$row2[$pos_ind]."_".$$row2[$ref_ind]."_".$$row2[$alt_ind];
		}
		##########
		
		#Work on the hash containing all equal rows
		if ( $eq_rows > 0 ){
			
			$identical_calls = $identical_calls + $eq_rows;
			#print "There are $eq_rows equal variants to the primary..\n";#DEBUGCODE			
			
			######Set the primary call
			#Get the VCF call from the primary algorithm (without the VC_ALGORITHM info field)
			my $primary_vc = -1;
			for (my $i = 0; $i < $eq_rows; $i++){
				my $subjrow = $row_hash->{$i};
				if ( ! defined ($vcf->get_info_field($$subjrow[$info_ind],"VC_ALGORITHM")) ){
					$primary_vc = $i;
				}
			}
			#If for this variant there is no value for the main caller, the primary is the first one
			if ( $primary_vc < 0 ){
				$primary_vc = 1;
			}
			my $primary_row = $row_hash->{$primary_vc};
			
			#print "The VCF has ".scalar(@$primary_row)." columns\n";#DEBUGCODE
			
			# print the primary row
			print  OUT $$primary_row[$chr_ind]."\t".$$primary_row[$pos_ind]."\t".$$primary_row[$id_ind]."\t".$$primary_row[$ref_ind]."\t".$$primary_row[$alt_ind].
			"\t".$$primary_row[$qual_ind]."\t".$$primary_row[$filter_ind]."\t".$$primary_row[$info_ind]."\t".$$primary_row[$genot_ind];
			for (my $i=1; $i <= scalar(@samples); $i++){
					print OUT "\t".$$primary_row[$genot_ind + $i];
			}
			print OUT "\n";	
			############
			
			
			#Get the items for each sample for the primary row
			#print "Get the items for each sample for the primary row: ";#DEBUGCODE			
			my @prim_items = ();
			foreach my $sample (@samples){
				my $sample_col = $vcf->get_column($primary_row, $sample);
				#Get the index of the specific FORMAT information
				my $form_idx = $vcf->get_tag_index($$primary_row[$genot_ind],$gt_str,':');
				#Get the value for the specified sample using the index
				my $val = $vcf->get_field($sample_col,$form_idx);
				push(@prim_items,$val);
				#print "$val ";#DEBUGCODE
			}
			
			#Go through the row_hash and print all the lines having at least one genotype information 
			#for one of the samples different from the primary algorithm call
			for (my $i = 0; $i <= $eq_rows; $i++){
				#Just avoid the primary call
				if ( $i ne $primary_vc ){
					#print "Checking $i th equal row\n";#DEBUGCODE
					my $subjrow = $row_hash->{$i};

					#Count how many values are equal
					my $eq_vals = 0;
					my $samplenum = 0;
					foreach my $sample (@samples){
						#print "Checking for $sample\n";#DEBUGCODE

						my $sample_col = $vcf->get_column($subjrow, $sample);
						#Get the index of the specific FORMAT information
						my $form_idx = $vcf->get_tag_index($$subjrow[$genot_ind],$gt_str,':');
						#Get the value for the specified sample using the index
						my $formval = $vcf->get_field($sample_col,$form_idx);
						if (defined $formval){
							#Get the zigosities
							my $val = get_zygosity($formval);
							my $prim_val = get_zygosity($prim_items[$samplenum]);
							#print "Checking if $val = $prim_val\n";#DEBUGCODE
							#Match zygosities
							if ($val eq $prim_val ){
								$eq_vals++;
							}								
						}
						$samplenum++;
					}
					
					#print "There are $eq_vals identical values for the ".scalar(@samples)." samples..\n";#DEBUGCODE			
					my $alg = "GATK";
					if (defined  $vcf->get_info_field($$subjrow[$info_ind],"VC_ALGORITHM")){
							$alg = $vcf->get_info_field($$subjrow[$info_ind],"VC_ALGORITHM");
					}
						
					#If the values are not all equal along the samples, write the line
					if ( $eq_vals < scalar(@samples)){
						#print "Keeping genotype from $alg !\n";#DEBUGCODE
						print OUT $$subjrow[$chr_ind]."\t".$$subjrow[$pos_ind]."\t".$$subjrow[$id_ind]."\t".$$subjrow[$ref_ind]."\t".$$subjrow[$alt_ind].
						"\t".$$subjrow[$qual_ind]."\t".$$subjrow[$filter_ind]."\t".$$subjrow[$info_ind]."\t".$$subjrow[$genot_ind];
						for (my $i=1; $i <= scalar(@samples); $i++){
								print OUT "\t".$$subjrow[$genot_ind + $i];
						}
						print OUT "\n";
						
						#Update counts
						if (not defined $saved_calls->{$alg} ){
							$saved_calls->{$alg} = 1; 
						}else{
							$saved_calls->{$alg} = $saved_calls->{$alg} + 1;
						}
					}
					#else{
						#print "I won't keep genotype from $alg since are equal!\n";#DEBUGCODE
					#}
				}
			}			
		}else{
			#Just print the row
			print OUT $$row1[$chr_ind]."\t".$$row1[$pos_ind]."\t".$$row1[$id_ind]."\t".$$row1[$ref_ind]."\t".$$row1[$alt_ind].
			"\t".$$row1[$qual_ind]."\t".$$row1[$filter_ind]."\t".$$row1[$info_ind]."\t".$$row1[$genot_ind];
			for (my $i=1; $i <= scalar(@samples); $i++){
				print OUT "\t".$$row1[$genot_ind + $i];
			}
			print  OUT "\n";	
		}
	}	
	close(OUT);
	
	print_and_log("..Process finished! Identical calls: $identical_calls. Calls saved from different programs:",$log_file);
	foreach my $a ( keys %{$saved_calls}){
		print_and_log(" $a (".$saved_calls->{$a}.") ",$log_file);
	}

}


=head2 join_vcf_files_with_cat

 Title   : join_vcf_files_with_cat
 Usage   : join_vcf_files_with_cat( -file1 = first file
                              - $file2 = second file
                              - outFile = otuput file
                               );

 Function:  Join two files in input and write them on a third file using the system call CAT

 Returns : nothing

=cut
sub join_vcf_files_with_cat{
 my $file1 = shift;
 my $file2 = shift;
 my $outFile = shift;
  
 #If the file is already present and its dimensions equals the sum of the two file it will not create again  
 unless ( (-e $outFile) and ( (-s $outFile) == (-s $file1) + (-s $file2) ) ){
    print "Concatenating VCFs ".$file1." and ".$file2." in ".$outFile."...";
       
   my $command="cat $file1 > $outFile";
   
   
   try_exec_command($command) > 0 
     or die "ERROR [$?]: The cat of $file1 failed: ?\n";
     
     $command="grep -v '^#' $file2 >> $outFile"; 
     
     try_exec_command($command) > 0 
     or die "ERROR [$?]: The grep and append of $file1 and $file2 failed: ?\n";
     
   print " ..completed. \n";
  }else{print "$outFile already present. I will not create..\n";} 

}

=head2 join_files_with_cat

 Title   : join_files_with_cat
 Usage   : join_files_with_cat( -file1 = first file
                              - $file2 = second file
                              - outFile = otuput file
                               );

 Function:  Join two files in input and write them on a third file using the system call CAT

 Returns : nothing

=cut
sub join_files_with_cat{
 my $file1 = shift;
 my $file2 = shift;
 my $outFile = shift;
  
 #If the file is already present and its dimensions equals the sum of the two file it will not create again  
 unless ( (-e $outFile) and ( (-s $outFile) == (-s $file1) + (-s $file2) ) ){
    print "Concatenating ".$file1." and ".$file2." in ".$outFile."...";
   my $command="cat ".$file1." ".$file2." > ".$outFile;
    
   try_exec_command($command) > 0 
     or die "ERROR [$?]: The concatenation of $file1 and $file2 failed: ?\n";
     
   print " ..completed. \n";
  }else{print "$outFile already present. I will not create..\n";} 

 #close($outFile);
}



=head2 join_files_with_zcat

 Title   : join_files_with_zcat
 Usage   : join_files_with_zcat( -folder = folder with the files
                              - $out_zip = otuput file
                              - compr_f_list = list of compressed files
                               );

 Function:  Joins a list of compressed files in input in a single compressed
			file

 Returns : nothing

=cut
sub join_files_with_zcat{
	my $folder = shift;
	my $out_zip = shift;
	my $compr_f_list = shift;
	
	#$out_zip = $folder."/".$out_zip;
	if ( ! -e $out_zip ){
	my $temp_folder = extract_file_folder($out_zip);
		my $temp_out = $temp_folder."/temp.out";
		foreach my $comp_file (@$compr_f_list){
			#chomp($data_file);
			
			my $command =  "zcat ".$folder."/".$comp_file." >> $temp_out";
			print "Executing: $command\n";
			if (try_exec_command($command) < 1){
				die "Cannot execute $command\n";
			}
		}
		print "Now compressing with gzip to $out_zip \n";
		#Compress with gzip
		gzip $temp_out => $out_zip;
		#Remove the temp file
		delete_file($temp_out);
		print $temp_out." removed..\n";
	}else{
			die "ERROR: file $out_zip already exists. Please delete that before to start...\n";
	}
}

=head2 check_joined_with_zcat

 Title   : check_joined_with_zcat
 Usage   : check_joined_with_zcat( -folder = folder with the files
                              - $out_zip = otuput file
                              - compr_f_list = list of compressed files
                              - fusion_file = zipped file with the merge
                               );

 Function:  Joins a list of compressed files in input in a single compressed
			file

 Returns : nothing

=cut
sub check_joined_with_zcat{
	my $folder = shift;
	my $out_zip = shift;
	my $compr_f_list = shift;
	my $fusion_file = shift;
	
	
	$out_zip = $folder."/".$out_zip;
	if ( ! -e $out_zip ){
	
		my $temp_out = $folder."/temp.out";
		foreach my $comp_file (@$compr_f_list){
			#chomp($data_file);
			my $command =  "zcat ".$folder."/".$comp_file." >> $temp_out";
			print "Executing: $command\n";
			if (try_exec_command($command) < 1){
				die "Cannot execute $command\n";
			}
		}
		print "Now compressing with gzip to $out_zip \n";
		#Compress with gzip
		gzip $temp_out => $out_zip;
		#Remove the temp file
		delete_file($temp_out);
		print $temp_out." removed..\n";
	}else{
			die "ERROR: file $out_zip already exists. Please delete that before to start...\n";
	}
}






=head2 merge_columns

 Title  : merge_columns
 Usage  : merge_columns( - filePath => 'the path to the file',
												- outFile => 'output file',
												- separator => 'to indicate which char separates columns'
												);
 Function: all the columns are separated and merged row by row
  
  Returns : a string

=cut		
sub merge_columns {
	my $file = shift;
	my $outFile = shift;
	my $in_sep = shift;
	
	my $separator;
	
	print "Merging columns from $file...\n";
	if (!(defined $in_sep)){
		$separator = "\t";
	}else{
		$separator = $in_sep;
	}
	open (FILE,"<$file") or die "Cannot open $file\n";
	open (OUT,">$outFile") or die "Cannot open $outFile\n";
	
	while (my $row = <FILE>){
			my @pieces = split($separator,$row);

			print OUT join("\n",@pieces);
	}
	
	close(OUT);
	close(FILE);
}

=head2 extract_col_from_file

 Title  : extract_col_from_file
 Usage  : extract_col_from_file( - filePath => 'the path to the file');

 Function: given a file and a name of a column (or the index), fetches all the values
					under the selected column
  
  Returns : a file with only one column

=cut		
sub extract_col_from_file {
	my $file = shift;
	my $colName = shift;
	my $outFile = shift;
	my $sep_c = shift;
	
	open (COL_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");
	
	my $header = <DATA>;
	
	my $sep = "\t"; 
	if (defined $sep_c){
		$sep = $sep_c;
	}	
	my @header_fields = split (/$sep/, $header);

	
	my $num_col = 0;
	
	#If the colName is the index of the column takes that
	if ( $colName =~ m/^\d+$/ ){
		$num_col = $colName;
	}
	#Otherwise it is the column name
	else{
		for my $field (@header_fields) {
			last if $field =~ m/$colName/;
			$num_col++;
		}		
	}
	
	while ( my $row = <DATA> ) {
		last unless $row =~ /\S/;
		chomp($row);
		my @fields = split /$sep/, $row;
		print COL_FILE $fields[$num_col]."\n";
	}

	close(DATA);
	close(COL_FILE);	
}




=head2 get_col_index

 Title  : get_col_index
 Usage  : get_col_index( - filePath => 'the path to the file');

 Function: given a file and a name of a column, searches the index
						of the column in input
  
  Returns : the index of the column as integer (starts from zero)
  
=cut		
sub get_col_index {
	my $file = shift;
	my $colName = shift;
	
	open (FILE, "<$file") or die ("Cannot open $file\n");
	
	my $header = <FILE>;
	
	my @header_fields = split (/\t/, $header);
	my $num_col = -1;
	
	for my $field (@header_fields) {
		$num_col++;
	  last if $field =~ m/$colName/;
	}
	
	close(FILE);	
	return $num_col;
}


=head2 extract_columns_from_file

 Title  : extract_columns_from_file
 Usage  : extract_columns_from_file( - filePath => 'the path to the file');
																			columns => a string with numbers separated by comma
 Function: given a file and a list of a columns indexes, returns only the 
					selected column. This script allows also to resort a table with the
					given order
					N.B. Starts from 0 as they are indexes of a PERL array
					
  Returns : a file with only the needed columns in the order as given in input

=cut		
sub extract_columns_from_file {
	my $file = shift;
	my $columns = shift;#a string with numbers separated by comma
	my $outFile = shift;
	
	open (COL_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");
	#Get the columns in an array
	my @cols = split(",",$columns );
	
	my $sep = "\t";
	#print "Getting columns $columns from $file..\n";
	#Write a new file with the given columns only
	while ( my $row = <DATA> ) {
		#print "$row";
		my $newline = "";
		#last unless $row =~ /\S/;
		chop($row);
		my @fields = split ("\t", $row);
		if (scalar(@cols) > scalar(@fields) ){ die "ERROR in extract_columns_from_file: More columns ids (".scalar(@cols).") than columns (".scalar(@fields).")\n";}
		foreach my $col (@cols){
			$newline .=  $fields[$col]."\t";
		}
		#Remove the last \t
		chop($newline);
		#print newline
		print COL_FILE $newline."\n";
	}

	close(DATA);
	close(COL_FILE);	
}




=head2 extract_special_chars

 Title  : extract_special_chars
 Usage  : extract_special_chars( - filePath => 'the path to the file');

 Function: this function searches special chars inside a file and returns a string with them
  
  Returns : a string

=cut		
sub extract_special_chars {
  my $filePath = shift;
  my $goodLetters = shift;
      
  my @collected = ();  
  open(FILE,"<$filePath");
        
  while(<FILE>){
    my $line = $_;
    $line =~ s/[$goodLetters]//;
    $line =~ s/\s//; 
  } 
  close(FILE);
} 

=head2 create_folder

 Title  : create_folder
 Usage  : create_folder( - folderPath => 'the path to the folder');

 Function: creates the folder if it does not exist
  
  Returns : nothing
  
=cut		
sub create_folder {
	my $folderPath = shift;
	
	#Creation
	unless(-d $folderPath){
		print "$folderPath doesn't exists. Creating folder...\n";
		mkdir $folderPath or die "ERROR: can't create folder $folderPath\n ";
	}
}

=head2 vcf_insert_string_into_field
 Title   : vcf_insert_string_into_field
 Usage   : vcf_insert_string_into_field(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - an array of values in input
						# - a column name after which you want to add the array as the new column
						#You can insert a string at the end of the any of the fields of a VCF file.	
						#Once that it finds the field then for each line adds at the end of the field
						#the separator and the strin given
 Returns : nothing

=cut
sub vcf_insert_string_into_field {
	my $vcffile = shift;
	my $string = shift;#the string to insert
	my $colname = shift;#string representing the column
	my $outvcf = shift;#output name
	my $sep = shift;#separator to use, default: \t
	
	#If a separator is defined, use that, else it is ;
	unless (defined $sep) {$sep = ";";}
	
	my $tab = "\t";
	#Open input and output file
	open (IN,"<$vcffile") or die "Cannot open $vcffile\n";
	open (OUT,">$outvcf") or die "Cannot open $outvcf\n";
	
	my $row = "";
	my $totCols = 0;
	my $index = -1;
	while ( $row  = <IN> ) {
			
		#Jump (but print) all lines contaning comments
		if ($row =~ /^##/){
			print OUT $row;
			next;
		}
		
		#Once arrived at the header,get it
		if ( $row =~ /^#/){
			print OUT $row;
			chomp($row);
			my @header = split($tab, $row);

			
			#This function will take the index of the first element in the header
			#which matches the string col
			$index = first {$header[$_] eq $colname} 0..$#header;
			print "$colname is at position $index\n";

			$totCols = scalar(@header);
		}else{
			chomp($row);
			my @fields = split($tab,$row );
			#Take all the columns until the selected index
			my @toPrint = @fields[0 .. ($index-1)];
			#Print them separated with tab
			print OUT join($tab,@toPrint);
			#Then appends the string at the column wanted
			print OUT $tab.$fields[$index].$sep.$string;#
			
			#Then finish to print the rest of the line
			if ($index + 1 <= $totCols) {
				print OUT $tab;
				@toPrint = @fields[$index+1 .. $#fields];
				print OUT join($tab,@toPrint);
			}
			print OUT "\n";
		}
	}
	close(IN);
	close(OUT);
	
}

=head2 insert_col_in_file
 Title   : insert_col_in_file
 Usage   : insert_col_in_file(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
			# - a file made as a table
			# - an array of values in input
			# - a column name after which you want to add the array as the new column
				if it is a number i it will find the ith column
			# - job: (i: insertion; s: substitution) what to do, insert or substitute?
			# - position: b: before; a: after	
			# opens the file, reads the header and finds the index in which should be 
			# inserted the column. Then it goes through all the rest of the file copying
			# the content in a new file and at that index it always insert the corresponding
			# value from the array set.
						
 Returns : nothing

=cut
sub insert_col_in_file {
	my $fileTable = shift;
	my ($newCol) = shift;#the vector to insert
	my $index = shift;#number of the column
	my $pos = shift;#b: before ; a: after
	my $outfile = shift;
	my $sep = shift;#separator to use, default: \t

	#If a separator is defined, use that, else it is tab
	unless (defined $sep) {$sep = "\t";}

	#Open input and output file
	open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	open (OUT,">$outfile") or die "Cannot open $outfile\n";

	my $row_ind = 0;
	my $isfirst = 0;

	#Look at the position (before or after)
	if ( $pos eq 'b'){
		$index--;
		if ( $index < 0 ){
			$isfirst = 1;
			$index = 1;
		}
	}else{
		$index++;	
	}

	#For each line insert the given column
	while ( my $row = <IN> ){
		chomp($row);
		my @fields = split($sep,$row );
		
		if ( $isfirst ) {
			#Then prints immediately the column wanted
			print OUT $$newCol[$row_ind].$sep;#$newCol->[$row_ind]	
		}
		
		#Take all the columns until the selected index
		my @toPrint = @fields[0 .. ($index-1)];
		#Print them separated with tab
		print OUT join($sep,@toPrint);
		if ( $isfirst == 0){
			#Then prints the column wanted
			print OUT $sep.$$newCol[$row_ind];#$newCol->[$row_ind]
		}

		#Then finish to print the rest of the line
		if ($index + 1 <= scalar (@fields)) {
			print OUT $sep;
			@toPrint = @fields[$index .. $#fields];
			print OUT join($sep,@toPrint);
		}
		print OUT "\n";
		
		$row_ind++;
	}

	close (IN);
	close (OUT);
 
}

=head2 insert_col_in_file_table
 Title   : insert_col_in_file_table
 Usage   : insert_col_in_file_table(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
			# - a file made as a table
			# - an array of values in input
			# - a column name after which you want to add the array as the new column
			# opens the file, reads the header and finds the index in which should be 
			# inserted the column. Then it goes through all the rest of the file copying
			# the content in a new file and at that index it always insert the corresponding
			# value from the array set.
			
			#To put as first colum use $colNum = -1
			
 Returns : nothing

=cut
sub insert_col_in_file_table {
	my $fileTable = shift;
	my ($newCol) = shift;#the vector to insert
	my $colNum = shift;#string representing the column
	my $job = shift;#Kind of job to exectue: substitution or insertion
	my $sep = shift;#separator to use, default: \t

	#If a separator is defined, use that, else it is tab
	unless (defined $sep) {$sep = "\t";}
 my $newFile = "newfile.temp";
 
 #Open input and output file
 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
 open (OUT,">$newFile") or die "Cannot open $newFile\n";
 
 my $row = <IN>;
 chomp($row);
 my @header = split($sep,$row );
 
 
 #This function will take the index of the first element in the header
 #which matches the string col
 my $index;
 if ( correct_type($colNum,"positiveint") ){
	$index = $colNum;
 }else{
  $index = first {$header[$_] eq $colNum} 0..$#header;
  print "$colNum is at position $index\n";
 }
 
 my $totCols = scalar(@header);	
 my $row_ind = 0;
 
 #Here print the header
 if ($job eq 's'){
	print OUT join($sep,@header);
	print OUT "\n";
	
 }elsif ($job eq 'i'){
	 close (IN);
	 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	 $index++;
	}
 
 #For each line insert the given column
 while ( $row = <IN>){
		chomp($row);
		my @fields = split($sep,$row );
		#Take all the columns until the selected index
		my @toPrint = @fields[0 .. ($index-1)];
		#Print them separated with tab
		print OUT join($sep,@toPrint);
		#Then prints the column wanted
		print OUT $sep.$$newCol[$row_ind];#$newCol->[$row_ind]
		
		my $next_ind;
		if ($job eq 's'){
			$next_ind = $index+1;
		}elsif ($job eq 'i'){
			$next_ind = $index;
		}
		
		#Then finish to print the rest of the line
		if ($next_ind+1 <= $totCols) {
			print OUT $sep;
			@toPrint = @fields[$next_ind .. $#fields];
			print OUT join($sep,@toPrint);
		}
		print OUT "\n";
		
		$row_ind++;
	}
	
 close (IN);
 close (OUT);
 
	move($newFile,$fileTable) or print "ERROR: unable to move $newFile in $fileTable\n";
	delete_file($newFile);
}


=head2 delete_columns
 Title   : delete_columns
 Usage   : delete_columns(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - a set of column index that you want to remove (starting from 0)
						# opens the file, goes through all the file copying
						# the content in a new file and it will jump those indexes
 
 Returns : nothing

=cut
sub delete_columns {
	my $fileTable = shift;
	my $cols_to_rem = shift;#the vector to insert
	my $newFile = shift;
	
	my $sep = "\t";
	
	#If an output is defined, use that, else it is tab
	unless (defined $newFile) {$newFile = "newfile.temp";}
 
	 #Open input and output file
	 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	 open (OUT,">$newFile") or die "Cannot open $newFile\n";
 

	#For each line insert the given column
	while ( my $row = <IN>){
		chomp($row);
		my @fields = split($sep,$row );
		my $col_num = 0;
		my $new_row = "";
		foreach my $field (@fields){
			if ( ! grep {/\b$col_num\b/} @$cols_to_rem){
				#Then prints the column wanted
				$new_row = $new_row.$field.$sep;
			}
			$col_num++;
		}
		chop($new_row);
		print OUT $new_row."\n";
	}

	close (IN);
	close (OUT);
	
	unless (defined $newFile) {
		move($newFile,$fileTable) or print "ERROR: unable to move $newFile in $fileTable\n";
		delete_file($newFile);
	}
	
}


=head2 delete_rows_containing
 Title   : delete_rows_containing
 Usage   : delete_rows_containing(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - a symbol
						# opens the file, goes through all the file removing
						#all lines containing a given symbol in input
 
 Returns : nothing

=cut
sub delete_rows_containing {
	my $fileTable = shift;
	my $symbol = shift;#something to grep
	my $newFile = shift;
		
	#Open input and output file
	open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	open (OUT,">$newFile") or die "Cannot open $newFile\n";

	#For each line insert the given column
	while ( my $row = <IN>){
	 if ( $row !~ /$symbol/){
		print OUT $row;		 
	 }
	}
	close (IN);
	close (OUT);
}


=head2 get_rows_containing
 Title   : get_rows_containing
 Usage   : get_rows_containing(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - a symbol
						# opens the file, goes through all the file getting
						#all lines containing a given symbol in input
 
 Returns : nothing

=cut
sub get_rows_containing {
	my $fileTable = shift;
	my $symbol = shift;#something to grep
	my $newFile = shift;
		
	#Open input and output file
	open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	open (OUT,">$newFile") or die "Cannot open $newFile\n";

	#For each line insert the given column
	while ( my $row = <IN>){
	 if ( $row =~ /$symbol/){
		print OUT $row."\n";		 
	 }
	}
	close (IN);
	close (OUT);
}


=head2 save_hash
 Title   : save_hash
 Usage   : save_hash(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  this subroutine saves an hash given in input into a file
					using the package Storable
 
 Returns : nothing

=cut
sub save_hash {
	my ($hash) = shift;
	my $hash_file = shift;
	#print "Storing hash..";
	#print $$hash->{'B'}."\n";#DEBUGCODE
	use Storable;
	store $$hash, $hash_file;
}

=head2 load_hash
 Title   : load_hash
 Usage   : load_hash(  - hash_file -> the complete path of the file where the hash is saved
                      );

 Function:  this subroutine loads an hash from a file
 
 Returns : the hash loaded

=cut
sub load_hash{
	my $hash_file = shift;
	return retrieve($hash_file);
}

=head2 check_URL

 Title   : check_URL
 Usage   : check_URL( - link => the link to check
                    );

 Function: It checks if an URL to external resource is working. If it is, it returns 1. If it is not then it will check
          what kind of error is coming out and alert the user with a message.

 Returns : 1 if link is good, -1 otherwise

    
=cut
sub check_URL{
    my $link = shift;
    
    my $agent = LWP::UserAgent->new();
    my $res = $agent->get($link);

    my $retVal = -1;
    
    if($res->is_success) {
        #print $res->as_string();
        $retVal = 1;
    } else {
        #if ( (split(/ /,$res->status_line))[0] eq '404'){
          #print "ERROR: ".$res->status_line.".\n";
        #}
        print "ERROR: ".$res->status_line."\n";
        $retVal = -1;
    }
    return $retVal;
}

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
	if ( -e $listFile ){
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
	}

	return @list;
}



1;
