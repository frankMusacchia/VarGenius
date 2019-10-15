#!/usr/bin/perl

# INPUT:
# A file with a list of:
# chr	pos
# The number of bases before and after

# OUTPUT: 
# Returns the output of GATK HaplotypeCaller -bamout
# and gives the bam realigned by this tool

# Author: Francesco Musacchia (2017)

# NEEDED LIBRARIES
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;#Used to write an usage

my $bam_file = "";
my $positions_list = "";
my $interval_size = "";
my $vcf_file = "";
my $out_file = "";
my $forceActive = "";
my $program_name = "get_bamout.pl";
my $gatk_path = "/cineca/prod/applications/gatk/3.8/jre--1.8.0_73/GenomeAnalysisTK.jar";
my $java_path = "/cineca/prod/compilers/jre/1.8.0_73/none/bin/java";
my $reference = "/pico/work/TELET_TIGEM/ngsworkspace/references/genomes/ucsc.hg19.fa";

parse_command_line_args();

#Open file and save list of lines
open (FILE,"<$positions_list") or die "Cannot open $positions_list\n";
my @lines = <FILE>;
close(FILE);
	
#For each line if the position and chromosome are correct
foreach my $line ( @lines ){
	chop($line);
	my @fields = split("\t",$line);

	
	if ( !correct_type($fields[1],"positiveint") ){
		die "ERROR : ".$fields[1]." must be a positive number. \n";	
	}

	if ($fields[0] !~ /chr/ ){
		die "ERROR : ".$fields[0]." must be a correct chromosome indication (chr1,chr2,...) \n";	
	}
	
	#Create the interval string
	if ( $interval_size eq ''){
		$interval_size = 1000;
	}
	my $interval;
	$interval = $fields[0].":".($fields[1]-$interval_size)."-".($fields[1]+$interval_size);
	
	#Execute the command for the pipeline using qsub
	my $bamout = $out_file.".".$interval.".bam";
  $bamout =~ s/:/_/;
	
	my $params = "";
	if ($forceActive ){
		$params = " -forceActive -disableOptimizations ";
	}
  
	my $command = "$java_path -Xmx4g -jar $gatk_path -T HaplotypeCaller -R $reference -I $bam_file -o $vcf_file -L $interval $params -bamout $bamout";
	
	print "\nExecuting: $command\n";
	 (system $command)== 0	or die "ERROR [$?]: an error occurred while executing $command. \n";	
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

	
	my $howToUse = "\nUse with: perl $program_name \n".
  "-p|--positions_list: A file with a list of positions for which you want to get the bamout from GATK (chromosome\tposition).\n".
  "-s|--interval_size: number of bases to use for left and right interval size\n".
  "-b|--bam_file: path to the bam file\n".
  "-v|--vcf: OUTPUT VCF file\n".
  "-force|--forceActive: Force GATK to print also unactive regions\n".
  "-o|--out: output suffix (use a path and the name of the analysis only, position will be added)";
  
  #  Parse options
  GetOptions(
           "help" =>        \$HELP,
           "version" =>      \$VERSION,
           "p|positions_list=s" => \$positions_list,   #It's mandatory
           "s|interval_size=s" => \$interval_size,
           "b|bam_file=s" => \$bam_file,
           "v|vcf=s" => \$vcf_file,
           "force|forceActive" => \$forceActive,
           "o|out=s" => \$out_file
            );

  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "Version: $program_name \n";
    exit;
  }
 
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


=head2 correct_type
Title  : correct_type
 Usage  : correct_type( -number => 'a number to check',
                      -typeWanted => 'the expected type');

 Function: 	 Check if the number is of the type given in input
  
  Returns 1 if true

=cut
sub correct_type {
  my $number = shift;
	my $typeWanted = shift;
	
  my $type = "";
	my $ret = 0;

  if ($typeWanted eq "positiveint"){
    #contains onli digits
    if ( $number=~ m/^\d+$/) {
      $ret=1;
    }
  }elsif ($typeWanted eq "real"){
    #Contains digits starting with (+-), separated by '.' and can be added the 'E-x'
    if ( $number =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ){
        $ret=1;
    }
  }
  return $ret;
}
