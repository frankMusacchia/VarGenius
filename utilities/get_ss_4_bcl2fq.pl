#!/usr/bin/perl

#VarGenius utilities
#This script generates the samples sheet  for bcl2fastq
#taking in input the samples and their index codes


use strict;
use warnings;


my $samples_tab = $ARGV[0];					
my $investigator = $ARGV[1];
my $expname = $ARGV[2];
my $date = $ARGV[3];
my $rc_active = $ARGV[4];


my $ss_out = "ss_$samples_tab.csv";
my $index_map_f = "index_map_f.txt";

#Generate an hash with all the indexes associtated with their codes
my $index_map;
open (MAPF,"<$index_map_f") or die "Cannot open $index_map_f\n"; 
while (my $row = <MAPF>){
	chop($row);
	my @fields = split("\t",$row);
	$index_map->{$fields[0]} = $fields[1];
}
close (MAPF);

################HEader section
my $header = "[Header],,,,\n".
"IEMFileVersion,4,,,\n".
"Investigator Name,$investigator,,,\n".
"Experiment Name,$expname,,,\n".
"Date,$date,,,\n".
"Workflow,GenerateFASTQ,,,\n".
"Application,FASTQ Only,,,\n".
",,,,\n".
",,,,\n".
"[Settings],,,,\n".
",,,,\n".
",,,,\n".
"[Data],,,,\n";


################Samples and indexes
my @used_ind_couples = ();
open (IN,"<$samples_tab") or die "Cannot open $samples_tab\n"; 
open (OUT,">$ss_out") or die "Cannot open $ss_out\n"; 
#Print the header
print OUT $header;
my $num_line = 0;
my $lane_present = 0;
while (my $row = <IN>){
	chop($row);
	my @fields = split("\t",$row);
	
	if ($num_line == 0){
		if ( $row =~ /Lane|lane/){
			$lane_present = 1;
			print OUT "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project\n";
		}else{
			print OUT "Sample_ID,Sample_Name,index,index2,Sample_Project\n"
		}
	}else{
		my $lane = -1;
		if ($lane_present == 1){
			$lane = $fields[$lane_present];
		}
		my $index1 = $fields[$lane_present + 1];
		my $index2 = $fields[$lane_present + 2];

		
		my $comb = "$lane\_$index1\_$index2";
		#print "ind1 : $index1, ind2 : $index2 , comb : $comb\n";
		if ( not grep {/\b$comb\b/} @used_ind_couples){
			if ($lane_present == 1){
				print OUT "$lane,";
			}
			my $dna_index2 = $index_map->{$index2};
			#If the Rev Compl flag is active, do it
			if ( $rc_active eq "RC"){
				$dna_index2 = rev_compl($dna_index2);
			}			
			#Sample_ID,Sample_Name,index,index2,Sample_Project
			#print $fields[0].",".$fields[0].",".$index_map->{$index1}.",".$index_map->{$index2}.",".$fields[0]."\n";#DEBUGCODE
			print OUT $fields[0].",".$fields[0].",".$index_map->{$index1}.",".$dna_index2.",".$fields[0]."\n";
			push (@used_ind_couples,$comb);
		}else{
			print "WARNING $comb already used!\n";
		}		
	}

	$num_line++;
}
close (IN);
close(OUT);


sub rev_compl {
	my $dna = shift;
	
	my $rc = reverse $dna;
	$rc =~ tr/ACGTacgt/TGCAtgca/;
	return $rc;
}
