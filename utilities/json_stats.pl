#!/usr/bin/perl
#To work at CINECA

use lib '/cineca/prod/compilers/perl/5.x/none/lib/perl5';
use Data::Dumper;#To print the hashes

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2017>  <Francesco Musacchia>

use strict;
use warnings;
use Getopt::Long;

binmode STDOUT, ":utf8";
use utf8;
 
use JSON;

my $file = $ARGV[0];

my $json;
{
  local $/; #Enable 'slurp' mode
  open my $fh, "<", $file;
  $json = <$fh>;
  close $fh;
}

my $data = decode_json($json);
#print Dumper $data;

for my $report ( @{$data} ) {
   my $samplename = $report->{external_id};
   my $sex =  $report->{'sex'}."\t";
   my $feats = $report->{'features'};
   $samplename =~ s/-/\_/;
    $samplename =~ s/\_P\d+//;  
	
	if ( $samplename !~ /00000/ and $samplename !~ /NAN0/ ){
   for my $feat (@$feats) {
	  print uc($samplename)."\t".$feat->{'id'}."\t".$feat->{'observed'}."\n"; 
	}
}
}



#my $studies_hash;
##############################CODE TO USE FOR NIG JSON OUTPUT
#foreach my $SampleId (@$DemuxResults){
	#print "\n\n$SampleId";
	##print "\n\n".$arr_el->{'attributes'}->{'accession'}." ".$arr_el->{'attributes'}->{'name'}."\n"; 
	##$studies_hash->{$arr_el->{'attributes'}->{'name'}} = $arr_el->{'attributes'}->{'accession'};
	##foreach my $key (keys %$arr_el){
	 ##print $key."\n"; 
	##}
#}
