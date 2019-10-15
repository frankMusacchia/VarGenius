
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
    
package LIB::html_management;
## html_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of HTML files for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( file_2_html html_page paths_hash_to_html var_stats_2_html
								table_2_html table_list_2_html);
}
use strict;
use warnings;
use CGI; #load cgi module to produce HTML output

#Using a library to manage files
use LIB::files_management qw(delete_file list_to_array extract_name);

use LIB::files_manipulation qw( );

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(print_and_log log_and_exit );

###########################################

=head2 paths_hash_to_html

 Title   : paths_hash_to_html
 Usage   : paths_hash_to_html( config_file => the config hash );

 Function: Reads the paths hash obtained from the outlist file and writes a file in HTML
					with for each PAGE and for each SAMPLE and for each READFILE the description linking to 
					the given path into the webserver.
					The variable $host is used to change the host name 
  
 Returns : a string that is actually an HTML page
 
=cut
sub paths_hash_to_html{
	my $cfg_hash = shift;
	my $page = shift;
	my $paths_hash = shift;
	my $workingFolder = shift;
	
	#The place where the analysis folders are located
	my $host = $cfg_hash->{"html_host"}; 
	my $cgi = CGI->new; #Instantiate a CGI class
	my $string = "<center>";

	#
	while (my ($sample_name, $value) = each %{ $paths_hash->{$page} } ) {
		#Print the sample name if exists
		if ($sample_name ne '-'){
			$string = $string."Sample: $sample_name<br>";
		}
		while (my ($readf_name, $value2) = each %{ $paths_hash->{$page}->{$sample_name}} ) {
			while (my ($desc, $value3) = each %{ $paths_hash->{$page}->{$sample_name}->{$readf_name}} ) {
				#the description indicates if the file must be printed or not [WEB]
				if ( $desc =~ /\[WEB\]/ ){
					#Extract only the file name to make the relative link
					my $fname = extract_name($paths_hash->{$page}->{$sample_name}->{$readf_name}->{$desc}->{'path'},0);
					
					if (defined $fname) {
						#$from_wf =~ s/$workingFolder//;
						my $link = $fname;#$host.$from_wf;
						#Concatenate description with sample and read file name
						if ($sample_name ne '-'){$desc = "$desc $sample_name ";}
						if ($readf_name ne '-'){$desc = "$desc $readf_name ";}
						#Concatenate the link
						$string = $string.$cgi->a({href=>$link},$desc)."<br>";						
						#code
					}	
				}

			}
		}
	}

 $string .= "</center>";
 #print_and_log( "String is (inside) $string ..\n",$log_file);
 return $string;
}


=head2 paths_hash_to_htmlOLD

 Title   : paths_hash_to_html
 Usage   : paths_hash_to_html( config_file => the config hash );

 Function: Reads the paths hash obtained from the outlist file and writes a file in HTML
					with for each PAGE and for each SAMPLE and for each READFILE the description linking to 
					the given path into the webserver.
					The variable $host is used to change the host name 
  
 Returns : a string that is actually an HTML page
 
=cut
sub paths_hash_to_htmlOLD{
	my $cfg_hash = shift;
	my $page = shift;
	my $paths_hash = shift;
	my $host = shift;
	my $workingFolder = shift;
		
	my $cgi = CGI->new; #Instantiate a CGI class
	
	my $string = "<center>";

		#
		while (my ($sample_name, $value) = each %{ $paths_hash->{$page} } ) {
			#Print the sample name if exists
			if ($sample_name ne '-'){
				$string = $string."Sample: $sample_name<br>";
			}
			while (my ($readf_name, $value2) = each %{ $paths_hash->{$page}->{$sample_name}} ) {
				while (my ($desc, $value3) = each %{ $paths_hash->{$page}->{$sample_name}->{$readf_name}} ) {
					#if ($desc ne '-'){print FILE $desc."\n"}
					#Extract the path and make the link
					my $from_wf = $paths_hash->{$page}->{$sample_name}->{$readf_name}->{$desc}->{'path'};
					if (defined $from_wf) {
						$from_wf =~ s/$workingFolder//;
						my $link = $host.$from_wf;
						#Concatenate description with sample and read file name
						if ($sample_name ne '-'){$desc = "$desc $sample_name ";}
						if ($readf_name ne '-'){$desc = "$desc $readf_name ";}
						#Concatenate the link
						$string = $string.$cgi->a({href=>$link},$desc)."<br>";						#code
					}
					

				}
			}
		}

 $string .= "</center>";
 #print_and_log( "String is (inside) $string ..\n",$log_file);
 return $string;
}

=head2 file_2_html

 Title   : file_2_html
 Usage   : file_2_html(   );

 Function: Reads a file line by line and transforms it in a TAGged HTML text
					It is supposed that each line is a distinct paragraph
 Returns : a string with the HTML code
 
=cut
sub file_2_html{
	my $inFile = shift;

 #Create
 #print_and_log( "Opening $inFile ..\n",$log_file);
 open(INFILE, "<$inFile") or warn "WARNING [$!]: can't open file: $inFile. Check permissions.\n"; 
 my $string = "<center>";
 while (my $line = <INFILE>){
	chop($line);
	#Use <p> to separate paragraphs
	$string = $string."<p>".$line."</p>\n";
	#print_and_log( "Line read $line ..\n",$log_file);
 }
 close(INFILE);
 $string .= "</center>";
 #print_and_log( "String is (inside) $string ..\n",$log_file);
 return $string;
}


=head2 table_list_2_html

 Title   : table_list_2_html
 Usage   : table_list_2_html(   );

 Function: Transforms a list of tables in HTML. The input is a file with a list of
						TableName TablePath
						And the tables are text file tab separated. The function table_2_html is
						iteratively called.
					
 Returns : a string with the HTML code
 
=cut
sub table_list_2_html{
	my $tables_list_file = shift;
	my $remove_temp = shift;

	my $html_string = "";

	
	if (-e $tables_list_file) {
		#print "Opening $tables_list_file ..\n";
		open(INFILE, "<$tables_list_file") or print "WARNING [$!]: can't open file: $tables_list_file. Check permissions.\n"; 
		
		#Transform all the tables in html
		my @tables_paths = list_to_array($tables_list_file,'NO_NEW_LINE');
		foreach my $table_path (@tables_paths){
			my @fields = split("\t",$table_path);
				if ( -e $fields[1]){
				$html_string .= "<p>".$fields[0];
				$html_string .= table_2_html($fields[1]);
				$html_string .= "</p>\n";
				if ( $remove_temp eq 'YES' ){
					delete_file($fields[1]);
				}
		 }
		}
		close(INFILE);		
	}#else{print "WARNING [$!]: file $tables_list_file does not exist.\n";}

	return $html_string;
}


=head2 table_2_html

 Title   : table_2_html
 Usage   : table_2_html(   );

 Function: Transforms a file with a table (columns tab separated)
					in a TAGged HTML text
					N.B. First line contains the header of the table! columns are TAB separated!!
					
 Returns : a string with the HTML code
 
=cut
sub table_2_html{
	my $inFile = shift;

	#print_and_log( "Opening $inFile ..\n",$log_file);
	open(INFILE, "<$inFile") or die "ERROR [$!]: can't open file: $inFile. Check permissions.\n"; 
	my $string = '<center><table border="1" style="width:50%">';

	my $num_row = 0;
	while (my $line = <INFILE>){
		chop($line);
		my @fields = split("\t",$line);
		$string .= "\t<tr>\n";
		#Header
		if ($num_row == 0){
			foreach my $field (@fields){
					$string .= "\t\t<th>$field</th>\n";
			}
		}else{
			foreach my $field (@fields){
					$string .= "\t\t<td>$field</td>\n";
			}			
		}
		$string .= "\t</tr>\n";
		$num_row++;
		#print_and_log( "Line read $line ..\n",$log_file);
	}
	close(INFILE);
	$string .= "</table></center>\n";
	#print_and_log( "String is (inside) $string ..\n",$log_file);
	return $string;
}


=head2 var_stats_2_html

 Title   : var_stats_2_html
 Usage   : var_stats_2_html(  config_file => the config hash
								);

 Function: Prints an HTML page using the statistics generated with R
					for the result of VarGenius.
 
 Returns :
 
=cut
sub var_stats_2_html{
	my $var_stats_f = shift;
	
	my $datasep = shift;
	
	my $totreads_block_ind = 1;
	my $freqs_block_ind = 2;
	my $segregation_block_ind = 3;
	my $dbsnpexist_block_ind = 4;
	my $zygosity_block_ind = 5;
	my $zyglegenda_block_ind = 6;
	my $db_analysis_statistics_block_ind = 7;
	
	#print_and_log( "Opening $inFile ..\n",$log_file);
	 open(INFILE, "<$var_stats_f") or die "ERROR [$!]: can't open file: $var_stats_f. Check permissions.\n"; 
	 my $string = '<center>';
	 my $start = 0;
	 my $tab_line = 0;
	 while (my $line = <INFILE>){
			chop($line);
			#A block of information starts with a set of # and a string
			if ($line =~ /#+\w+/){
				$start++;
				$line =~ s/#//g;
				$string .= "<br>".$line."<br>";
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind
					or $start == $db_analysis_statistics_block_ind ) {
					$tab_line = 0;
					$string .= '<table border="1" style="width:50%">';
				}
				#$line = <INFILE>; 
			}#and ends with only a set of #
			elsif ($line =~ /#+/) {
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind
				or $start == $db_analysis_statistics_block_ind ){
					$string .= '</table>';
				}
				$string .= "<br>";
			}else{
				#This instructions are to read and build tables
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind
				or $start == $db_analysis_statistics_block_ind ){
					$string .= "<tr>";
					my $tag_s;
					my $tag_e;
					my @parts = split("\t",$line);		
					if ($tab_line == 0){
						$tag_s = "<th>";
						$tag_e = 	"</th>";
						$tab_line++;
					}else{
						$tag_s = "<td>";
						$tag_e = 	"</td>";					
					}
					for (my $i = 0; $i <scalar(@parts) ; $i++){
							$string .= $tag_s.$parts[$i].$tag_e;
					}
					$string .= "</tr>";	
				}

				#Second table
				if ( $start == $freqs_block_ind){
					#if ($line =~ /^Freq/){
						#$string .= "<br><br>".$line."<br>";
						##$line = <INFILE>; 
					#}els
					if ($line =~ /^Chrom/){
						$string .= "<br>".$line." : ";
						#$line = <INFILE>; 
					}elsif ($line ne ''){
						$string .= $line.$datasep;
					}
				}				
	
				#Legenda
				if ( $start == $zyglegenda_block_ind or $start == $dbsnpexist_block_ind){
						$string = $string."<h6>".$line."</h6>";
				}					
			}
			#print_and_log( "Line read $line ..\n",$log_file);
	 }
	 close(INFILE);
	 $string .= "</center>";
	 return $string;
}


=head2 html_page

 Title : html_page
 Usage : html_page(  imgDir -> #Directory where the images are stored	
					 htmlFile -> Path to the html file to create
					 ($descrText) -> array pointer descriptive texts of the links to the other pages
					  $descrTextIndex -> Index of the descriptive text of this page
					 ($links) -> pointer to array with the links at the other pages
					 ($imgsKeywords)-> pointer to array of keywords to use to search the images inside the folder
		);

 Function: creates an html page using a directory where are searched the images. Prints the header using $descrTextIndex
			from ($descrText). Puts the logo from the path. A loop searches for all the images in $imgDir
			with title comprising a word inside ($imgsKeywords). Print a series of links to the other pages using a loop
			searching in ($links) whose names are in descrText. Finally prints information at the bottom of the page.

 Returns : nothing

=cut 
sub html_page{
 my $cfg_hash = shift;
 my $imgPath = shift;#Directory where the images are stored	
 my ($descrText) = shift;#descriptive texts of the links to the other pages
 my $descrTextIndex = shift;#Index of the descriptive text of this page
 my ($links) = shift;#array with the links at the other pages
 my $imgsKeywords = shift;#keywords to use to search the images inside the folder
 my $infoText = shift;#A string to print as informations in the BODY of the page
 my $html_folder = shift;
 
 my $cgi = CGI->new; #Instantiate a CGI class
 
 #Open a new html file
 open (OUT,">$html_folder/".$$links[$descrTextIndex]) || die "ERROR [$!]: can't open file ".$$links[$descrTextIndex]."\n";
  
 print OUT $cgi->start_html("VarGenius ".$cfg_hash->{'VarGeniusVer'}."- ".$$descrText[$descrTextIndex]), # Begin HTML page
       #$cgi->header, # First we create a header
       $cgi->center; 
 if ( $imgPath ne 'none'){
	print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$cfg_hash->{'VarGeniusLogo'}."' ALT='image' ALIGN=TOP> <br>";
       #$cgi->h1('<i>Annocript '.$configHash->{'AnnocriptVer'}.'</i>'), # create the Tag <h1>TEXT</h1>
 }
 print OUT $cgi->h3($$descrText[$descrTextIndex]."\n");
 #BODY
 print OUT $cgi->p($infoText);
 
 #Open images folder if it exists
 if ( $imgPath ne 'none'){
	 opendir IMG, $imgPath or die "ERROR [$!]: cannot open dir $imgPath.\n";
	 my @images= readdir IMG;
	 my @imgsSorted = sort {$a cmp $b} @images;#order alfabetically the images names
	closedir IMG;
	 #Search and write HTML code to show images if it exitsts
	 #print OUT "</center>";
	 if(defined($imgsKeywords)){
		 foreach my $img (@imgsSorted){
			#if (scalar(grep {/\b$img\b/} @$imgsKeywords) == 1){
			if ($img =~ /$imgsKeywords/){
				print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$img."' ALT='image' ALIGN=TOP> <br>";
			} 
		 }
	 }
 }
 

	#FOOTER
 #print the website map
 for ( my $link=0; $link<scalar(@$links);$link++ ){
	print OUT $cgi->a({href=>$$links[$link]},$$descrText[$link]);
	print OUT " | ";
 }

 #Print information
 print OUT $cgi->h5('VarGenius '.$cfg_hash->{'VarGeniusVer'}.' - Copyright of Tigem'),
      $cgi->p(scalar localtime); # create <p>Text</p>
 print OUT "</center>";
 print OUT $cgi->end_html; # end of HTML
 close(OUT);
}


1;
