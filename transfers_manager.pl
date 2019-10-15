#!/opt/software/perl/bin/perl
#To work at TIGEM
#!/usr/bin/perl
#To work at CINECA

 
# Package VarGenius - Variant Discovery and Annotation Tool
#
# Author: Francesco Musacchia (2016)
#
#This software is a pipeline to discover SNPs and Variants through clinical
#sequencing data and to update our internal databases

#GLI HASHTAG #TOCHECK E .. SONO DA CONTROLLARE PER CODICE DA AGGIUNGERE

#Works at CINECA
#use lib '/pico/home/userexternal/fmusacch/VarGenius/'; #PICO CINECA
#use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5/x86_64-linux-thread-multi';
#use lib '/cineca/prod/libraries/bioperl/1.6.923/gnu--4.8.3/lib/perl5';
#Parallel::ForkManager
#use lib qw(/pico/home/userexternal/fmusacch/bin/perl_lib/Parallel-ForkManager-1.17/lib/);
use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';

#TIGEM
#use lib '/home/users/musacchia/VarGenius';
###########################################################

use strict;
use warnings;

use FindBin;#To search libraries
use lib $FindBin::Bin;#To search libraries

download_selected($ARGV[0],$ARGV[1],$ARGV[2]);
#check_md5sum_selected($ARGV[0],$ARGV[1],$ARGV[2]);

#put_in_ASPERA($ARGV[0],$ARGV[1]);

=head2 copy_selected

 Title  : copy_selected
 Usage  : copy_selected( 
									- dirname => 'the folder where to copy',
                  - fileList => list of samples to download
									- destDir => destination folder
                      );

 Function: Copies samples folders from a given folder and changes names 
					accordingly to those given into a file

 Returns :

=cut
sub copy_selected {
	my $dir = shift;
	my $destDir = shift;
	my $fileList = shift;

	my $intern_fold = "";#"";
		
	my $correspondences;
	#The file list has two columns, one with the original name and another with the new name
	#Reading all in one with <FILE> command
	open (FILE,"<$fileList") or die "Cannot open $fileList\n";
	my @lines = <FILE>;
	close(FILE);
	
	#Fill an hash with the correpsondences
	foreach my $line (@lines){
		chop($line);
			my @fields = split("\t",$line);
			$correspondences->{"Sample_".$fields[0]} = $fields[1];
			print "The correspondence of Sample_".$fields[0]." will be ".$fields[1]."\n";
	}
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs = readdir DIR;#All the files are in this list
	print "Copying from: $dir (containing ".scalar(@samples_dirs)." directories)\n";	
	
	foreach my $sample ( @samples_dirs ){
		print "Checking $sample \n";
		if ( defined $correspondences->{$sample} ){
			my $outFolder = $destDir."/".$correspondences->{$sample};
			print "Copying $sample to..$outFolder\n";
				if (! (-d $outFolder) ){
					#Check if directory exists, otherwise it creates it
					unless(-d $outFolder){
						print "Creating folder $outFolder...\n";
						mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
					 }
					my $command = "cp -Lr $dir/$sample".$intern_fold."/* $outFolder";
					print "Executing: $command\n";
					try_exec_command($command);
				}else{
					print "Unable to copy $outFolder Folder exists...\n";
				}
		}
	}
			
}
	

=head2 move_selected

 Title  : move_selected
 Usage  : move_selected( 
									- dirname => 'the folder where to copy',
                  - fileList => list of samples to download
									- destDir => destination folder
                      );

 Function: Moves folders in $dir to $destDir changing the name if needed
					new names must be in fileList in format:
					 oldname\tnewname

 Returns :

=cut
sub move_selected {
	my $dir = shift;
	my $destDir = shift;
	my $fileList = shift;
	my $suffix = shift;
	
	my $intern_fold = "";#"";
	#Here if you put IF undef(x) it does not work
	if (! defined $suffix){
			$suffix = "";
	}	
	#print "suffix:$suffix-\n";
	my $correspondences;
	#The file list has two columns, one with the original name and another with the new name
	#Reading all in one with <FILE> command
	open (FILE,"<$fileList") or die "Cannot open $fileList\n";
	my @lines = <FILE>;
	close(FILE);
	
	#Fill an hash with the correpsondences
	foreach my $line (@lines){
		chop($line);
			my @fields = split("\t",$line);
			$correspondences->{"$suffix".$fields[0]} = $fields[1];
			#print "The correspondence of $suffix".$fields[0]." will be ".$fields[1]."\n";
	}
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs = readdir DIR;#All the files are in this list
	#print "Copying from: $dir (containing ".scalar(@samples_dirs)." directories)\n";	
	
	foreach my $sample ( @samples_dirs ){
		#print "Checking $sample \n";
		if ( defined $correspondences->{$sample} ){
			my $outFolder = $destDir."/".$correspondences->{$sample};
			print "Moving $sample to $outFolder\n";
				if (! (-d $outFolder) ){
					my $command = "mv $dir/$sample $outFolder";
					print "Executing: $command\n";
					try_exec_command($command);
				}else{
					print "Unable to copy $outFolder Folder exists...\n";
				}
		}
	}
}	
	

####################put_in_ASPERA
#Performs the command to put files using Aspera. 
#(Used for uploads in Rd-Connect (SolveRD)
#################
sub put_in_ASPERA {
	my $infolder = shift;
	my $samples_f = shift;	
	
	open (FILE,"<$samples_f") or die "Cannot open $samples_f\n";
	my @lines = <FILE>;
	close(FILE);
	
	
	foreach my $line ( @lines ){
		chop($line);
		
		my @fields = split("\t",$line);
		my @files = split(",",$fields[1]);
		
		my $insert = 1;
		my $paths_str = "";
		my $folder = $infolder."/".$fields[0]."/";
		foreach my $file (@files) {
			if ( -z $folder."/".$file  or  !(-e $folder."/".$file )){
				print "WARNING: $folder/$file  is not available\n";	
				$insert = 0;
			}else {
				$paths_str .= "$folder/$file ";
			}	
		}
	
		if ( $insert == 1 ){
			#Execute the command for the pipeline using qsub
			my $command = "mlia shares repository upload --url=https://fileshare.rediris.es/ ".
			'--username=f.musacchia@tigem.it --password="`hYB/VJ]mX9M3B!=" '.
			'--to-folder=/CNAG_CRG_PLATFORM/PLATAFORMA_CNAG_CRG/RDCONNECT_U04 --sources=@args '.$paths_str;
			
			print "\nExecuting: $command\n";
			 (system $command)== 0	or die "ERROR [$?]: an error occurred while executing $command. \n"; 
		}else{
				print 
		}
				 
	}
}	

=head2 check_md5sum_selected

 Title  : check_md5sum_selected
 Usage  : check_md5sum_selected( 
									- dirname => 'the folder where to copy',
                  - fileList => list of samples to download
									- destDir => destination folder
                      );

 Function: Given an origin folder, a file with sample names and a destination folder
					checks whether the md5sums are equal and reports errors
					

 Returns : just copies to the destination folder

=cut
sub check_md5sum_selected {
	my $dir = shift;
	my $destDir = shift;
	my $fileList = shift;
	
	#Put the names from the file in an array
	my @list = ();
	my $intern_fold = "/Files/";#"";
	
	if (defined $fileList){
		@list = list_to_array($fileList,'NO_NEW_LINE');
	}
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs = readdir DIR;#All the files are in this list
	print "Checking from: $dir\n";
	
	foreach my $sample ( @samples_dirs ){
			if ( $sample !~ /^\..*/ ){
				#print "Checking $sample..\n";
				if (defined $fileList){
					if ( grep {/$sample/} @list){
						print "analyzing $sample..\n";
						my $outFolder = $destDir."/".$sample;

						#Get Md5sum and substitute the two separator spaces with a single tab
						print "Executing md5sum for origin folder...\n";
						my $md5sum_origin = `md5sum $dir/$sample$intern_fold/* | sed 's/ \+/\t/g'`;
						#Get Md5sum and substitute the two separator spaces with a single tab
						print "Executing md5sum for destination folder...\n";
						my $md5sum_dest = `md5sum $outFolder/* | sed 's/ \+/\t/g'`;
						#Check all the fastq transferred
						my $res = compare_md5sums($md5sum_origin,$md5sum_dest);
						if ( length($res) > 0){
								print "ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n";
						}else{
								print "...DONE! All is ok!";
						}
					}					
				}else{

					print "You did not give a list. All samples will be checked..\n";
					my $outFolder = $destDir."/".$sample;

					print "Executing md5sum for origin folder...\n";
					my $md5sum_origin = `md5sum $dir/$sample$intern_fold/* | sed 's/ \+/\t/g'`;
					#Get Md5sum and substitute the two separator spaces with a single tab
					print "Executing md5sum for destination folder...\n";
					my $md5sum_dest = `md5sum $outFolder/* | sed 's/ \+/\t/g'`;
					#Check all the fastq transferred
					my $res= compare_md5sums($md5sum_origin,$md5sum_dest);
					if ( length($res) > 0){
							print "ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n";
					}else{
								print "...DONE! All is ok!";
						}
				}
			}
	}
}



=head2 download_selected

 Title  : download_selected
 Usage  : download_selected( 
									- dirname => 'the folder where to copy',
                  - fileList => list of samples to download
									- destDir => destination folder
                      );

 Function: This subroutine I used to copy files from basespace following
					the link and using a list of MotorPlex identifiers since in Basespace
					I found the overall set of samples generated since the times..
					I consider only the number to be the sample but the folder contains
					Motor_NUMBER
					

 Returns : just copies to the destination folder

=cut
sub download_selected {
	my $dir = shift;
	my $destDir = shift;
	my $fileList = shift;
	
	#Put the names from the file in an array
	my @list = ();
	my $intern_fold = "/Files/";#"";
	
	if (defined $fileList){
		@list = list_to_array($fileList,'NO_NEW_LINE');
	}
	
	#Create also a log file
	my $log_file = "$fileList.log";
	open(FILE,"<$log_file") or die "Cannot open $log_file..\n";

	
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs = readdir DIR;#All the files are in this list
	 
	print_and_log ("Downloading from: $dir\n",$log_file);
	foreach my $sample ( @samples_dirs ){
			if ( $sample !~ /^\..*/ ){
				#$sample =~ s/Motor\_//;
				#print "Found $sample..\n";
				if (defined $fileList){
					if ( grep {/$sample/} @list){
						print_and_log ("Using $sample that is in the list..\n",$log_file);
						#print "Using $sample that is in the list..\n";
						my $outFolder = $destDir."/".$sample;
							if (! (-d $outFolder) ){
								
								#Check if directory exists, otherwise it creates it
								unless(-d $outFolder){
									#print "Creating folder $outFolder...\n";
									print_and_log ("Creating folder $outFolder...\n",$log_file);
									mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
								 }
								my $command = "cp -Lr $dir/$sample".$intern_fold."/* $outFolder";
								#print "Executing: $command\n";
								print_and_log ("Executing: $command\n",$log_file);
								try_exec_command($command);
								
								#Get Md5sum and substitute the two separator spaces with a single tab
								my $md5sum_origin = `md5sum $dir/$sample$intern_fold/* | sed 's/ \+/\t/g'`;
								#Get Md5sum and substitute the two separator spaces with a single tab
								my $md5sum_dest = `md5sum $outFolder/* | sed 's/ \+/\t/g'`;
								#Check all the fastq transferred

								my $res = compare_md5sums($md5sum_origin,$md5sum_dest);
								if ( length($res) > 0){
									print_and_log ("ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n",$log_file);
									#print "ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n";
								}else{
									print_and_log ("All samples correctly transferred into $outFolder!\n",$log_file);
									#print "All samples correctly transferred into $outFolder!\n";
								}
								
							}else{
								print_and_log ( "Unable to copy! $outFolder folder exists...\n",$log_file);
								#print "Unable to copy! $outFolder folder exists...\n";
								}
					}					
				}else{
					print_and_log ( "You did not give a list. All samples will be downloaded..\n",$log_file);
					#print "You did not give a list. All samples will be downloaded..\n";
					my $outFolder = $destDir."/".$sample;
						if ( ! (-d $outFolder) ){

						#Check if directory exists, otherwise it creates it
						unless(-d $outFolder){
							print_and_log ("Creating folder $outFolder...\n",$log_file);
							#print "Creating folder $outFolder...\n";
							mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
						 }
							my $command = "cp -Lr $dir/$sample".$intern_fold." $outFolder";
							#print "Executing2: $command\n";
							#print_and_log ("Executing2: $command\n",$log_file);
							try_exec_command($command);

							#Get Md5sum and substitute the two separator spaces with a single tab
							my $md5sum_origin = `md5sum $dir/$sample.$intern_fold/* | sed 's/ \+/\t/g'`;
							#Get Md5sum and substitute the two separator spaces with a single tab
							my $md5sum_dest = `md5sum $outFolder/* | sed 's/ \+/\t/g'`;
							#Check all the fastq transferred
							my $res= compare_md5sums($md5sum_origin,$md5sum_dest);
							if ( length($res) > 0){
								print_and_log ("ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n",$log_file);
									#print "ERROR: while copying $sample to $outFolder. The following md5sums do not match:\n$res\n";
							}else{
								print_and_log ("All samples correctly transferred into $outFolder!\n",$log_file);
								#	print "All samples correctly transferred into $outFolder!\n";
								}
						}else{
							print_and_log ( "Unable to copy! $outFolder folder exists...\n",$log_file);
							#print "Unable to copy $outFolder Folder exists...\n";
						}	
				}
			}
	}
	close(FILE);
}

=head2 compare_md5sums

 Title  : compare_md5sums
 Usage  : compare_md5sums( - dirname => 'the folder to check',
                      );

 Function: compares md5sum into two variable containing a list of 
						md5sum\tfilepath\nmd5sum2\tfilepath2\n
 
 
 Returns : a string of paths of files with non matching md5sums

=cut
sub compare_md5sums{
	my $md5sum_origin = shift;
	my $md5sum_dest = shift;
	
	my $res = "";
	#Two hashes to make the match faster
	my $md5_h_o;
	my $md5_h_d;
	#print "ORIGIN: $md5sum_origin \n DEST: $md5sum_dest\n";
	#Extract lines
	my @md5lines_o = split("\n",$md5sum_origin);
	my @md5lines_d = split("\n",$md5sum_dest);
	
	#Put the mds5sums in hashes indexed with paths
	foreach my $md5line_o (@md5lines_o){
		#$md5line_o =~ s/ \+/\t/g;
		my @md5line_o_fields = split("  ",$md5line_o);
		#print $md5line_o_fields[0] ." e ".$md5line_o_fields[1]."  \n";
		$md5_h_o->{extract_name($md5line_o_fields[1],0)} = $md5line_o_fields[0];
	}

	foreach my $md5line_d (@md5lines_d){
		#$md5line_d =~ s/ \+/\t/g;
		my @md5line_d_fields = split("  ",$md5line_d);
		#print $md5line_d_fields[0] ." e ".$md5line_d_fields[1]."  \n";
		$md5_h_d->{extract_name($md5line_d_fields[1],0)} = $md5line_d_fields[0];
	}	
	#Compares all the md5sums using an hash and saves the filename if
	#the dimensions are not equal
	foreach my $key ( keys %$md5_h_o){
		print "Comparing MD5: ".$md5_h_o->{$key}." and ".$md5_h_d->{$key}." ($key)\n";
		if ($md5_h_o->{$key} ne $md5_h_d->{$key} ){
				$res .= " ".$key;
		}
	}
	
	return $res;
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


=head2 write_old_sample_sheet

 Title  : write_old_sample_sheet
 Usage  : write_old_sample_sheet( - dirname => 'the folder to check',
                      );

 Function: Writes a sample sheet by using the folder of the fastq files.

 Returns : nothing

=cut
sub write_old_sample_sheet{
	my $head_file = shift;
	my $dir = shift;
	my $mode = shift;
	
	my $sep = ",";
	my $no_val = "";
	my $target_bed = "";
	my $r1 = '_R1_';#'pe_1';
	my $intern_fold = "";#"/Files/";
	my $group_name = "";
	
	#Opens the dir with all the samples
	opendir DIR, $dir or die "ERROR : cannot open dir $dir";
	my @samples_dirs= readdir DIR;#All the files are in this list
	my $final_info = "";
	if ($mode eq "motor"){
		$final_info = "4".$sep."2".$sep."5000".$sep."4".$sep."0.12".$sep."30".$sep."0".$sep."generic".$sep."patient";
		$target_bed = "motor7_haloplex.bed";
	}elsif ($mode eq "exome"){#Exomes
		$final_info = "4".$sep."2".$sep."250".$sep."6".$sep."0.03".$sep."30".$sep."1".$sep."generic".$sep."patient";
		$target_bed = "clinical_exome.bed";
	}else{
			die "Need a mode to insert the parameters! [exome, motor]\n";
	}
	
	print "There are ".scalar(@samples_dirs)." in $dir..\n";
	my $couple;
	#Get and print header
	open (HEAD,"<$head_file") or die "ERROR: Cannot open $head_file. The program will exit..\n";
	my $header = <HEAD>;
	close(HEAD);
	#print "S.No.,heritability,gender,diseaseID,userID,diseaseVal,freqCalc,sampleID,fqDir,fq1,fq2,sampleName,groupName,groupViewName,sampleRef,target_bed,GroupID,ploidy,dcov,maxg,percBed,minQS,rmdup,exometype,analysistype\n";
	print $header;
	#A control starts on each sample among all the samples        
	my $numFile= 1;
	 my $fastq1 = "";
	 my $fastq2 = "";
	foreach my $sample_dir (@samples_dirs){
		if ($mode eq 'motor'){
			$group_name = $sample_dir;
		}
		my @lanes_used = ();
		#If the folder name does not start with a dot and is a  directory
		if ( ($sample_dir !~ /^\..*/) and (-d $dir.$sample_dir)){
			
			#print "Entering $sample_dir...\n";
			#my $dir1 = "/home/ngsworkspace/runs/ MotorPlex_v7_18122015/Samples/".$file
			my $dir1 = $dir.$sample_dir;
			#print $dir1."\t";
			$dir1 .= "$intern_fold";
			opendir DIR1, $dir1 or die "ERROR : cannot open dir $dir1";
			my @files= readdir DIR1;#All the files are in this list
			#A control starts on each file of the folder        

			foreach my $file (@files){
				#print "Using $file..\n";
				
				if ( $file =~ /$r1/ ){
					$fastq1 = $file;
				}else{
					$fastq2 = $file;
				}
			}
			#print "for folder $sample_dir used ".scalar(@lanes_used)."\n";
				print $numFile.$sep.$no_val.$sep.$no_val.$sep.$no_val.$sep."14".$sep."FALSE".$sep."FALSE".$sep.$sample_dir.$sep;
				print $dir1.$sep;
				print $fastq1.$sep.$fastq2.$sep;
				print $sample_dir.$sep.$group_name.$sep.$group_name.$sep."hg19".$sep.$target_bed.$sep.$final_info."\n";
				$numFile++;
		}
	}
}



=head2 find_and_merge_from_basespace_dl

 Title  : find_and_merge_from_basespace_dl
 Usage  : find_and_merge_from_basespace_dl( - dirname => 'the folder to check',
                      );

 Function: Given a folder searches files for Read 1 (R1) and read 2 (R2).
			Then, specifically this script works for those cases when
			there are more fastq for the same sample coming from different lanes.
			Hence, more files with the same name but distinguished by _L00X_
			could be present.
			This script first get all the files for R1 and R2,
			sorts them and finally calls a function to merge this files using zcat. 
			If the fileList is given with the names of the motor, then only those are
			considered for the merge

			Use rem_suffix to say if there is or not the suffix (e.g. Motor_) in the fileList
			
 Returns : nothing

=cut
sub find_and_merge_from_basespace_dl {
	my $dir = shift;
	my $out_dir = shift;
	my $fileList = shift;	
	
	my $suffix = "Motor";
	my $rem_suffix = 'NO';
	
	print "Starting find and merge procedure!\n";
	#Put the names from the file in an array
	my @list = ();
	if (defined $fileList){
		print "I will use only those in $fileList..\n";
	  @list = list_to_array($fileList,'NO_NEW_LINE');
	}
	
	my $bsUselessFold = "Files";
	opendir DIR, $dir or die "ERROR : cannot open $dir \n";
	my @folders= readdir DIR;#All the files are in this list
	
	#For each folder in the directory makes a new identical one in 
	#the output directory and calls the merge subroutine
	print "Found ".scalar(@folders)." at : $dir\n";
	foreach my $fold (@folders){
		print "Folder: $fold\n";
		if ($fold !~ /^\..*/) {
			#print "Using folder $fold\n";
			#If the list of folders is given than check that it exists
			if ( defined $fileList){
					print "Searching and merging only those in $fileList\n";
					my $outFolder;
					if ( $rem_suffix eq 'YES'){
						$fold =~ s/Motor\_//;#Remove the word Motor_
						$outFolder = $out_dir."/Motor\_$fold";
					}else{
						$outFolder = $out_dir."/$fold";
					}
					if ( grep {/$fold/} @list){
						print "\n\n---------------------------------------\n";
						print "Input: $fold - Output: $outFolder\n";
						#Check if directory exists, otherwise it creates it
						unless(-d $outFolder){
							print "Creating folder $outFolder...\n";
							mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
						 }
						if ( $rem_suffix eq 'YES'){
							merge_fastq($dir."/Motor\_".$fold."/$bsUselessFold",$outFolder);	
						}else{
							merge_fastq($dir."/".$fold."/$bsUselessFold",$outFolder);	
						}
						
					}
			}else{
					print "\n\n---------------------------------------\n";
					print "Searching and merging from folder: $fold...\n";
					my $outFolder = $out_dir."/$fold";
					#Check if directory exists, otherwise it creates it
					unless(-d $outFolder){
						print "Creating folder $outFolder...\n";
						mkdir $outFolder or die "ERROR: can't create folder $outFolder. Check permissions. \n";
					 }
					merge_fastq($dir."/".$fold."/$bsUselessFold",$outFolder);	
			}
		}
	}
}

=head2 merge_fastq

 Title  : merge_fastq
 Usage  : merge_fastq( - dirname => 'the folder to check',
                      );

 Function: Given a folder searches files for Read 1 (R1) and read 2 (R2).
			Then, specifically this script works for those cases when
			there are more fastq for the same sample coming from different lanes.
			Hence, more files with the same name but distinguished by _L00X_
			could be present.
			This script first get all the files for R1 and R2,
			sorts them and finally calls a function to merge this files using zcat. 
			

 Returns : nothing

=cut
sub merge_fastq {
	my $dir = shift;
	my $out_dir = shift;
	
	opendir DIR, $dir or die "ERROR : cannot open $dir \n";
	my @files= readdir DIR;#All the files are in this list
	
	#Gets the files for both Reads1 and reads2 usin grep
	my $r1 = "_R1_";
	my $r2 = "_R2_";
	my @r1_files = grep {/$r1/} @files;#Get all the R1 files
	my @r2_files = grep {/$r2/} @files;#Get all the R2 files


	#Sort the filenames in a new array
	#print "Results: ".scalar(@r1_files)."\n";
	my @r1_sorted = sort  @r1_files;
	my @r2_sorted = sort  @r2_files;
	
	#Do it for read1
	if ( scalar(@r1_files) ){
		#Get the name for the output by removing the lane number from the first file
		my $out_file_r1 = $r1_sorted[0];
		$out_file_r1 =~ s/\_L001//;
		
		print "Putting: \n";
		foreach my $comp (@r1_sorted){
				print $comp."\n";
		}
		print "inside one single file: $out_file_r1\n";	
		
		#Merge files together
		if ( ! -e $out_dir."/".$out_file_r1){
			print "Merging the files together\n";
			join_files_with_zcat($dir,$out_dir."/".$out_file_r1,\@r1_sorted);
			#Check if the number of reads in the merged file is equal to the sum of sorted reads
			check_merged_from_basespace_dl($dir,$out_dir."/".$out_file_r1,\@r1_sorted);
		}else{
			print "ERROR: cannot merge the files. ".$out_dir."/".$out_file_r1." exists..\n";
		}
	}else{ print "No file with $r1 present..\n";}

	#Do it for read2
	if ( scalar(@r2_files) ){
		my $out_file_r2 = $r2_sorted[0];
		$out_file_r2 =~ s/\_L001//;

		print "Putting: \n";
		foreach my $comp (@r2_sorted){
				print $comp."\n";
		}
		print "inside one single file: $out_file_r2\n";	
		
		if ( ! -e $out_dir."/".$out_file_r2){	
			join_files_with_zcat($dir,$out_dir."/".$out_file_r2,\@r2_sorted);
			#Check if the number of reads in the merged file is equal to the sum of sorted reads
			check_merged_from_basespace_dl($dir,$out_dir."/".$out_file_r2,\@r2_sorted);
		}else{
			print "ERROR: cannot merge the files. ".$out_dir."/".$out_file_r2." exists..\n";
		}
	}else{ print "No file with $r2 present..\n";}

}



=head2 check_merged_from_basespace_dl

 Title  : check_merged_from_basespace_dl
 Usage  : check_merged_from_basespace_dl( - in_dir => 'the folder to check',
																		-file_merge => the merged file,
																		- sorted => the array of file names
                      );

 Function: Given a folder and an array of fastq files sorted this subroutine
					counts with grep -c '^@' the number of reads present summing for each file
					and then will compared with the merged file in input $file_merge
 
 Returns : prints an error if there is a mismatch

=cut
sub check_merged_from_basespace_dl {
	my $in_dir = shift;
	my $file_merge = shift;
	my $file_names = shift;
	
	my $total_reads = 0;
		print "Checking $file_merge..\n";
	foreach my $file (@$file_names){
		
			my $command = 'zcat '.$in_dir.'/'.$file.' | grep -c "^@" ';
			#print "Executing: $command..\n";
			my $read_num =  `$command`;
			chomp($read_num);
			#print $read_num;
			$total_reads = $total_reads + $read_num;
			#print "Reads in $file: $read_num\n";
	}
	print "Total reads in $in_dir: $total_reads\n";
	
	my $command = 'zcat '.$file_merge.' | grep -c "^@" ';
	#print "Executing: $command..\n";
	my $read_merge =  `$command`;
	chomp($read_merge);
	print "Merged reads in $file_merge: $read_merge\n";
	
	if( $total_reads != $read_merge ){
			print "ERROR: number of reads in the merged file is $read_merge while total are  $total_reads\n";
	}else{
			print "$file_merge is ok!\n";
	}
	
}



######################################################################
####################UTILITIES#########################################

###############################################################




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
