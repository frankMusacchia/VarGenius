#!/usr/bin/perl

#Given a list of compid returns the information about what are the samples continaing it and the sequencing type

#vargenius.pl path
use lib '/pico/work/TELET_UDP/VarGenius/';

#Using a library for database management
use LIB::db_management qw(get_id_if_exists_from_db do_query_select_all );
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(configFile2Hash);
					
my $program_config = $ARGV[0];
my $user_config = $ARGV[1];
my $input_f = $ARGV[2];

my $cfg_hash;

configFile2Hash($program_config,\$cfg_hash);

configFile2Hash($user_config,\$cfg_hash);

print "Getting info for all variants of $input_f\n";
open (FILE, "<$input_f") or die "Cannot open $input_F\n";
print "Variant\tAnalysis\tSample\tGT\tSequencingType\n";
while (my $row = <FILE>){
		chop($row);
		#print "Analysing $row...\n";
		#Get the variant id
		my $varid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_variants_table'},$cfg_hash->{'db_var_id'},
															$cfg_hash->{'db_var_compid'},"'".$row."'");	

		#Get all samples and analyses with the given variant id
		my $query = "SELECT ".$cfg_hash->{'db_sample_id'}." FROM  ".$cfg_hash->{'db_genotype_sample_table'}.
		" WHERE ".$cfg_hash->{'db_var_ids'}." LIKE ANY ( values('$varid'),('$varid,%'),('%,$varid'),('%,$varid,%') ) ;";	
		#print "Executing: $query\n";#DEBUGCODE
		my $sample_ids = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$query,$cfg_hash->{'db_sample_id'});

		foreach my $sample_id (keys %{$sample_ids}){
				my $analysisid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_analysis_id'},
															$cfg_hash->{'db_sample_id'},$sample_id);	
				my $sample_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_sample_table'},$cfg_hash->{'db_sample_name'},
															$cfg_hash->{'db_sample_id'},$sample_id);
				my $seq_type = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_sequencingtype'},
															$cfg_hash->{'db_analysis_id'},$analysisid);	
				my $analysis_name = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
															$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_analysis_name'},
															$cfg_hash->{'db_analysis_id'},$analysisid);				
				#Get the genotype info
				my $var_query = "SELECT ".$cfg_hash->{'db_genotype_gt'}." FROM ".$cfg_hash->{'db_genotype_sample_table'}." WHERE ".$cfg_hash->{'db_var_ids'}." LIKE ".
			"ANY ( values('$varid'),('$varid,%'),('%,$varid'),('%,$varid,%') ) AND sampleid=$sample_id ;";	

				my $genotypes = do_query_select_all($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},$cfg_hash->{'db_pass'},$var_query,$cfg_hash->{'db_genotype_gt'});
				my $genotype = "";	
				
				foreach my $gen (keys %{$genotypes}){
						$genotype = $gen;
						#print $gen."\n";#DEBUGCODE
				}
				print "$row\t$analysis_name\t$sample_name\t".$genotype."\t$seq_type\n"									
		}	
}
close(FILE);
