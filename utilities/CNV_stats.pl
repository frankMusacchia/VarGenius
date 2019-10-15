#!/usr/bin/perl

#Script to intersect CNVs found with various methods and print the length 
#of the intersection. I used it on 201908 to verify the different CNVs found with 
#different methods

#Author: Francesco Musacchia
#HOW TO USE:
#perl CNV_stats.pl cnvout_ED,cnvout_XHMM,cnvout_VG ED,XHMM,VG
#The table must contain a column with the sample and a column with the compid (eg. chrX:31500000-37600000_DUP)

use strict;
use warnings;
use File::Copy;#To manage files


####
#Change only bedtools_path and comment uncomment paths to programs output
#Comment the common CNVs path if you dont want to use them!
#CINECA
my $bedtools_path = "/pico/work/TELET_TIGEM/bin/bedtools2/bin/bedtools";


my $intersectbed_prog = "intersect";
my $subtractbed_prog  = "subtract";
my $exomeDepth_out_f = "";
my $xhmm_out = "";
my $conifer_out = "";
my $commonbed = "";

#An hash to store families information
my $family_h;

my $remove_common_cnvs = 0;

my @chr_types = ("auto","sex");
my $symb_gt = ">";
my $symb_eq = "==";
my $outfolder="";
my @overlaps = ("0.0001","0.30","0.50","0.90");
my $cnvouts = $ARGV[0];
my $algos = $ARGV[1];
my @cnvouts = split(",",$cnvouts);
my @algos = split(",",$algos);
my $log_file = "CVN_comparisons.log";

my $inters_hash; 
my $all_cnv_hash;

#The table must contain a column with the sample and a column with the compid (eg. chrX:31500000-37600000_DUP)
for (my $i = 0; $i < scalar(@algos)-1; $i++){
	for (my $j = $i+1; $j < scalar(@algos); $j++){
		compare_CNV_results($algos[$i],$algos[$j],$cnvouts[$i],$cnvouts[$j],\$inters_hash,\$all_cnv_hash);		
	}	
}
#compare_CNV_results("ED","ACGH",$cnvouts[0],$cnvouts[0],\$inters_hash,\$all_cnv_hash,17,1);
#compare_CNV_results("ED","VG2",$cnvouts[0],$cnvouts[0],\$inters_hash,\$all_cnv_hash,17,1);
#compare_CNV_results("XHMM","ACGH",$cnvouts[0],$cnvouts[0],\$inters_hash,\$all_cnv_hash,15,1);
#compare_CNV_results("XHMM","VG2",$cnvouts[0],$cnvouts[0],\$inters_hash,\$all_cnv_hash,15,1);
#compare_CNV_results("ACGH","VG2",$cnvouts[0],$cnvouts[0],\$inters_hash,\$all_cnv_hash,1,1);

#my @algos = ("ED","XHMM","ACGH","VG2");
print "sample_cnv\tLength\t";
foreach my $algo (@algos){
	print "\t".$algo;
}
for (my $i = 0; $i < scalar(@algos); $i++){
	for (my $j = $i+1; $j < scalar(@algos); $j++){
		print "\t".$algos[$i]."^".$algos[$j];
	}
}
print "\n";

#Go thorgh the all_cnv_hash containing all the compid encountered during the process
#writes 1 for the algorithm if it was found there.
#Later it adds the intersection length if there was an intersection between the two algorithms
foreach my $sample_cnv (keys %{$all_cnv_hash}){
	my $newline = $sample_cnv."\t";
	#Get and print CNV lenght
	my @parts = split(":",$sample_cnv);
	my @ints = split("-",$parts[1]);
	my $len = $ints[1]-$ints[0];
	$newline .= $len."\t";
	#For each algorithm print 1 if it found the cnv
	for (my $i = 0; $i < scalar(@algos); $i++){
		if ( defined $all_cnv_hash->{$sample_cnv}->{$algos[$i]}){
			$newline .= $all_cnv_hash->{$sample_cnv}->{$algos[$i]}."\t";
		}else{
			$newline .= "0\t";
		}
	}
	#For each intersection
	for (my $i = 0; $i < scalar(@algos); $i++){
		for (my $j = $i+1; $j < scalar(@algos); $j++){
			if ( defined $inters_hash->{$sample_cnv}->{$algos[$i]."-".$algos[$j]}){
				$newline .= $inters_hash->{$sample_cnv}->{$algos[$i]."-".$algos[$j]}."\t";
			}else{
				$newline .= "0\t";
			}
		}
	}	
	
	chop($newline);
	print $newline."\n";
}



#This function opens the two files with CNV results from ED and XHMM and 
#builds an hash containing for each sample name a list of CNVs.
#
#Whether both XHMM and ED have a sample in the hash, the corresponding CNV
#list will be transformed in two BED files which can be intersected with bedtools
sub compare_CNV_results {
	my $algo1 = shift;
	my $algo2 = shift;
	
	my $cnvout1 = shift;
	my $cnvout2 = shift;

	#Hash with intersected CNV
	my ($inters_hash) = shift; 
	#Hash with uniq CNVs
	my ($all_cnv_hash) = shift;

	#XHMM
	my $compid_field1 = 1; 
	my $compid_field2 = 1;
			
	my $cnv_h1;
	my $cnv_h2;
	#ED
	my $sn_field1 = 0;
	my $sn_field2 = 0;

	
	#print "Starting intersectbed using $cnvout1 and $cnvout2. \n";#DEBUGCODE
	
	#ED
	open (FILE,"<$cnvout1") or die "ERROR: Cannot open $cnvout2\n ";
	my $numline = 0;
	while (my $row = <FILE>){
		
		if( $numline > 0 ) {
			chop($row);
			my @fields = split("\t",$row);
			#print $row."-".$fields[$compid_field1]." fieldnum: $compid_field1\n";
			my $compid = $fields[$compid_field1];
			if ( $compid ne 'NA'){
				#Get CNV type
				my $cnv_type = "";
				if ( $compid =~ /DEL/){
					$cnv_type = "DEL";
					$compid =~ s/\_DEL//;
				}else{
					$cnv_type = "DUP";
					$compid =~ s/\_DUP//;
				}
				my $sample_cnv = $fields[$sn_field1]."_".$cnv_type;
				if (! defined $cnv_h1->{$sample_cnv}){
					$cnv_h1->{$sample_cnv} = $compid;
				}else{
					$cnv_h1->{$sample_cnv} .= ",".$compid;
				}
				
				#Fill the hash with all CNVs
				my $compl_id = "$sample_cnv\_$compid";
				#print "Inserting $compl_id for $algo1\n";
				if (! defined ($$all_cnv_hash->{$compl_id}->{$algo1}) ){
					$$all_cnv_hash->{$compl_id}->{$algo1} = 1;
				}			
			}

		}
		$numline++;
	}
	close(FILE);

	#XHMM
	open (FILE,"<$cnvout2") or die "ERROR: Cannot open $cnvout2\n ";
	$numline = 0;
	while (my $row = <FILE>){
		chop($row);
		if( $numline > 0 ){
			my @fields = split("\t",$row);
			my $compid = $fields[$compid_field2];
			if ( $compid ne 'NA'){
				#Get CNV type
				my $cnv_type = "";
				if ( $compid =~ /DEL/){
					$cnv_type = "DEL";
					$compid =~ s/\_DEL//;
				}else{
					$cnv_type = "DUP";
					$compid =~ s/\_DUP//;
				}
				my $sample_cnv = $fields[$sn_field2]."_".$cnv_type;				
				if ( ! defined $cnv_h2->{$sample_cnv}){
					$cnv_h2->{$sample_cnv} = $compid;
				}else{
					$cnv_h2->{$sample_cnv} .= ",".$compid;
				}
				#Fill the hash with all CNVs
				my $compl_id = "$sample_cnv\_$compid";
				#print "Inserting $compl_id for $algo2\n";
				if (! defined $$all_cnv_hash->{$compl_id}->{$algo2}){
					$$all_cnv_hash->{$compl_id}->{$algo2} = 1;
				}		
			}

		}
		$numline++;
	}
	close(FILE);
	 

	#Set parameters for bedtools
	my $params = " -wo ";#-f $overlap ";
	#For each sample for which CNV have been found, if it has been found also
	#in the second case, create  BED files and 
	foreach my $sample_cnv (keys %{$cnv_h1}){
		if (defined $cnv_h2->{$sample_cnv}){
			my $bed1 = 	intervals_list_compid_2_bed($sample_cnv."_1",$cnv_h1->{$sample_cnv});
			my $bed2 = 	intervals_list_compid_2_bed($sample_cnv."_2",$cnv_h2->{$sample_cnv});
			
			my $outbed = "$sample_cnv\_intersect.bed";
			#print "Starting intersectbed using $bed1 and $bed2. Outbed in: \n";#DEBUGCODE	
			run_BEDTOOLS_intersectBed($bedtools_path,$bed1,$bed2,$outbed,$params,$log_file);
			#Now count lines of intersection 	
			#my $wc = (split(" ",`wc -l $outbed`))[0];	
			
			#Insert all the intersected region into an hash including the length of the intersection
			#for the two algorithms used
			open (FILE,"<$outbed") or die "ERROR: Cannot open $outbed\n";
			my $numline = 0;
			while (my $row = <FILE>){
				chop($row);
				my @pieces = split("\t",$row);
				my $compid1 = $sample_cnv."_".$pieces[0].":".$pieces[1]."-".$pieces[2];#."_".$pieces[3];
				my $compid2 = $sample_cnv."_".$pieces[3].":".$pieces[4]."-".$pieces[5];#."_".$pieces[7];
				$$inters_hash->{$compid1}->{$algo1."-".$algo2} = $pieces[6];
				$$inters_hash->{$compid2}->{$algo1."-".$algo2} = $pieces[6];
			}
			close(FILE);						
		}
	}
}


# Given a table where into a specific field there is an interval with a chromosome 
#It extracts from the id field a BED file
# From  1:65510-65625_DEL to:
#chr1    65510   65625
# ....
sub intervals_list_compid_2_bed {
	my $sample = shift;
	my $compid_list = shift;

	my $outbed = "$sample.bed";
	
	#Generate a new bed file
	open(NEWFILE,">$outbed");
	#if ($compid_list =~ /,/){chop($compid_list);}
	my @compids = split(",",$compid_list);	
	
	foreach my $compid (@compids){
		my @pieces = split(":",$compid);
		my $chr = $pieces[0];
		my @ints = split("-",$pieces[1]);
		#my @cnvt = split("_",$ints[1]);
		print NEWFILE $chr."\t".$ints[0]."\t".$ints[1]."\n";
	}
	close(NEWFILE);	
	
	return $outbed;
}


###############################################

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
		print OUT $row."\n";		 
	 }
	}
 close (IN);
 close (OUT);

}

=head2 extract_rows_from_table_where_col_contains_AWK

 Title   : extract_rows_from_table_where_col_contains_AWK
 Usage   : extract_rows_from_table_where_col_contains_AWK(   );

 Function: One liner using AWK filtering a table using a specific value that must be present
			on a determined column
			 
 Returns : nothing
=cut
sub extract_rows_from_table_where_col_contains_AWK {
	my $in_file = shift;
	my $col_num = shift;
	my $value = shift;
	my $out_file = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = 'awk -F"\t" '."'".'NR==1 || $'.$col_num.' == "'.$value.'" { print $0 }'."'"." $in_file > $out_file";
	#if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	#else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";	
	
	
}

=head2 extract_rows_from_table_where_AWK

 Title   : extract_rows_from_table_where_AWK
 Usage   : extract_rows_from_table_where_AWK(   );

 Function: One liner using AWK filtering a table using a specific value and a symbol
			for comparison 	on a determined column
 
 Returns : nothing
=cut
sub extract_rows_from_table_where_AWK {
	my $in_file = shift;
	my $col_num = shift;
	my $value = shift;
	my $symbol = shift;
	my $out_file = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = 'awk -F"\t" '."'".'NR==1 || $'.$col_num.' '.$symbol.' "'.$value.'" { print $0 }'."'"." $in_file > $out_file";
	#if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	#else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";	
	
	
}

=head2 run_BEDTOOLS_intersectBed

 Title   : run_BEDTOOLS_intersectBed
 Usage   : run_BEDTOOLS_intersectBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_intersectBed{
	my $bedtools_path = shift;
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = "";
	if ( $bed1 =~ /\.bam$/ and $bed2 =~ /\.bed$/ ){
		#$params .= " -bed ";
		$files_str = " -abam $bed1 -b $bed2 ";
	}elsif ($bed2 =~ /\.bam$/ and $bed1 =~ /\.bed$/){
		$files_str = " -abam $bed2 -b $bed1 ";
	}else{
		$files_str = " -a $bed1 -b $bed2 ";
	}
		
	#Execute the command
	my $command = "$bedtools_path $intersectbed_prog $files_str $params > $outFile";
	#print "Executing command: $command ";
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

=head2 run_BEDTOOLS_subtractBed

 Title   : run_BEDTOOLS_subtractBed
 Usage   : run_BEDTOOLS_subtractBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools subtract
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_subtractBed{
	my $bedtools_path = shift;
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = " -a $bed1 -b $bed2 ";
		
	#Execute the command
	my $command = " $bedtools_path $subtractbed_prog $files_str $params > $outFile";
	#if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	#else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}

# Given a table where into a specific field there is an interval with a chromosome 
#It extracts from the id field a BED file
# From  1:65510-65625 to:
#chr1    65510   65625
# ....
sub intervals_2_bed {
	my $input_f = shift;
	my $column = shift;
	my $outname = shift;
	
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	
	my $numline = 0;
	my $idfield = -1;
	
	#Generate a new bed file
	open(NEWFILE,">$outname");
				
	#Read all the zscores file
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		#Get the intervals from the header into the array
		if ($numline == 0){
			for (my $i = 1; $i < scalar(@fields) ; $i++){
				if ($fields[$i] eq $column){
						$idfield = $i;
				}
			}
		}
		else{
			my $interval = $fields[$idfield];
			#Print the header
			my @pieces = split(":",$interval);
			my $chr = $pieces[0];
			my @ints = split("-",$pieces[1]);
			print NEWFILE "chr".$chr."\t".$ints[0]."\t".$ints[1]."\n";
		}
		$numline++;
	}
	close(NEWFILE);	
	close(FILE);
}

# Given the output from ExomeDepth 
#It extracts from the id field a BED structure
# From  1:65510-65625 to:
#chr1    65510   65625
# ....
sub ExomeDepth_intervals_2_bed {
	my $input_f = shift;
	
	my $outname = shift;
	
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	
	my $numline = 0;
	my $idfield = -1;
	
	#Generate a new bed file
	open(NEWFILE,">$outname");
				
	#Read all the zscores file
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		#Get the intervals from the header into the array
		if ($numline == 0){
			for (my $i = 1; $i < scalar(@fields) ; $i++){
				if ($fields[$i] eq 'id'){
						$idfield = $i;
				}
			}
		}
		else{
			my $interval = $fields[$idfield];
			#Print the header
			my @pieces = split(":",$interval);
			my $chr = $pieces[0];
			my @ints = split("-",$pieces[1]);
			print NEWFILE "chr".$chr."\t".$ints[0]."\t".$ints[1]."\n";
		}
		$numline++;
	}
	close(NEWFILE);	
	close(FILE);
}

=head2 print_and_log

 Title   : print_and_log
 Usage   : print_and_log( - string -> the sentence that have to be print_and_log 
											- onlyLog -> a number);

 Function: will print (the string in input always in the log file and on the STDOUT
					if onlyLog is used then the print will be only in the log file
 Returns : nothing

=cut
sub print_and_log{
  my $string = shift;    
  my $logFile = shift; 
  my $onlyLog = shift;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  if ( defined $onlyLog){
    print LOG $string;
  }else{
    my $STDOUT = *STDOUT;
    my $LOG = *LOG;
    #Prints on both OUT
    for ($LOG, $STDOUT) { print $_ $string; }
  }
  #Close the log
	close(LOG)
}



=head2 try_exec_command

 Title   : try_exec_command
 Usage   : try_exec_command( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.

 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_command{
    my $command = shift;

    my $maxTimes = 5;
    my $success = -1;
    my $timesCount = 0;

    while ($success == -1 and $timesCount < $maxTimes){
        if ( (system $command) == 0) {
          $success = 1;
        }
        else{
         if ($? == -1) {
              print "failed to execute: $!\n";
          }
          elsif ($? & 127) {
              printf "child died with signal %d, %s coredump\n",
                  ($? & 127),  ($? & 128) ? 'with' : 'without';
          }
          else {
              printf "child exited with value %d\n", $? >> 8;
              
          }
         $timesCount++;
        }
    }
    return $success;
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
