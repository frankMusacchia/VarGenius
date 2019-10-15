#VarGenius - Variant Discovery and Annotation Tool
#Copyright (C) <2018>  <Francesco Musacchia>

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
    
    
package LIB::db_management;
## db_management.pm
#Permits the management of the database for the pipeline

#BEGIN { print (("  " x $main::x++) . "Beginning dbmanagement compile\n") }#DEBUGCODE

BEGIN
{
require Exporter;
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
@EXPORT_OK = qw( createDBAndTables getSampleConfiguration copyTable2DB getSampleConfiguration_locked
										update_analysis_status_locked exists_ss_locked insert_into_table_locked
										get_id_if_exists_from_db select_all_from_db db_present
										select_distinct_samples createTables insert_if_not_exists_locked
										parse_HPO_genes_to_phenotypes parse_GDI_scores parse_RVIS_scores
										do_query_select_all download_RefSeq_2_Genes do_query
										parse_RefSeq_2_Genes parse_genes2RefSeq parse_OMIM_2_Genes
										do_fetch_row_array db_print_selected_rows update_table
										select_element_from_db insert_only_if_not_exists_locked fetch_all_rows
										delete_file get_count_from_selected get_dbelement_ifexists_woconn
										update_table_woconn insert_only_if_not_exists_locked_woconn
										insert_into_table_locked_woconn delete_from_db do_fetch_row_array_woconn
										get_id_if_exists_from_db_like createHPOTable insert_into_table_woconn
										update_table_2 get_id_if_exists_from_db_woconn load_csv_2_db
										do_query_select_all_woconn fetch_all_rows_woconn);
}
use strict;
use warnings;

use DBI;#DB library


#Using a library to manage files
use LIB::files_management qw( download_file extract_columns_from_file
																delete_file extract_name);

#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(print_and_log log_and_exit try_exec_command);

#BEGIN { print (("  " x $main::x) . "Just used programs_management in dbmanagement compile\n") }#DEBUGCODE

#Standard useful functions for use everyday
use LIB::std_lib qw(correct_type );

#Pre-defined tables to be set here
my $sample_conf_table = "sample_conf";
my $sample_id_field = "sampleid";


=head2 delete_from_db

 Title   : delete_from_db
 Usage   : delete_from_db( -database => 'name of the database,
                               );

 Function: deletes elements from a table into the database.
			Needs database, dsn, user, pass, and table name
			if fields and values are qualified a WHERE clause is built
			while if they are not provided the items into the table
			are all deleted.

 Returns : returns nothing

=cut
sub delete_from_db{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;

	
	#Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
		my $whe_clause = "";
		if (defined $fields and defined $values){
			#Prepare the WHERE clause
			my @fields_arr = split(",",$fields);
			my @values_arr = split(",",$values);
			
			if ( scalar(@fields_arr) == scalar(@values_arr) ) {
				for ( my $i = 0; $i < scalar(@fields_arr); $i++){
					$whe_clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
				}
				#Remove the last AND
				chop($whe_clause);
				chop($whe_clause);
				chop($whe_clause);
			}				
		}
		
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "DELETE FROM $table ";
			if ( $whe_clause ne '' ){
					$query .= "WHERE $whe_clause";
			}
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}

=head2 getSampleConfiguration

 Title   : getSampleConfiguration
 Usage   : getSampleConfiguration( -database => 'name of the database,
                               );

 Function: gets in input an hash for the sample and returns all the
						variables given in the table

 Returns : returns the inputed hash

=cut
sub getSampleConfiguration{
	my ($in_hash) = shift;
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $sample_id = shift;
	
  #Opening database if needed 
	my $dbh;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT AND CREATE TABLE
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open an access to the database to collect informations
		my $select_handle;
  
		#Extracts information from database defined.
    #fetch informations from the database
    my $query = "SELECT * FROM $sample_conf_table WHERE $sample_id_field='".$sample_id."'";
    $select_handle = $dbh->prepare($query);
    die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
    $select_handle->execute() or die $select_handle->errstr;
    my $res = $select_handle->fetchrow_hashref;
    if (defined $res){
      #print "There is no result for query: $query\n";
       #Copy the output hash to the inputed one
			foreach my $elem (keys %$res){
					$$in_hash->{$elem} = $res->{$elem};
			}
    }
   
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}	



	
=head2 update_analysis_status_locked

 Title   : update_analysis_status_locked
 Usage   : update_analysis_status_locked( -database => 'name of the database,
                               );

 Function: gets in input 

 Returns : returns the inputed hash

=cut
sub update_analysis_status_locked{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $field1 = shift;
	my $value1 = shift;
	my $analysis = shift;
	#OPTIONAL
	my $field2 = shift;
	my $value2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");#LOCK
			#SET client_min_messages=WARNING;
			#Extracts information from database defined.
			#fetch informations from the database
			$query = "UPDATE $table SET $analysis=1 WHERE $field1=$value1";
			#Add second conditon if exists
			if (defined $field2 and defined $value2){
				$query .= " AND $field2='".$value2."'";
			}
			print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			print "Executing query: $query\n";
			$select_handle->execute() or die $select_handle->errstr;
		 
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}



=head2 insert_into_table_locked

 Title   : insert_into_table_locked
 Usage   : insert_into_table_locked( -database => 'name of the database,
                               );

 Function: makes a simple insert into the database using the LOCK 
					values may be multiple by separating them with comma
					This function does not activate the database connection hence it needs the 
					database reference is given in input
					
 Returns : returns the id of the sample inserted

=cut
sub insert_into_table_locked_woconn{
	my $dbh = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;

	
	#Opening database if needed 
	my $query;
	my $last_id = -1;
	if (defined $dbh){

		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) VALUES ($values)";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Undefined database connection!\n";}
	
	return $last_id;
}


=head2 insert_into_table_locked

 Title   : insert_into_table_locked
 Usage   : insert_into_table_locked( -database => 'name of the database,
                               );

 Function:  makes a simple insert into the database using the LOCK 
					values may be multiple by separating them with comma


 Returns : returns the id of the sample inserted

=cut
sub insert_into_table_locked{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;

	
	#Opening database if needed 
	my $dbh;
	my $query;
	my $last_id = -1;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) VALUES ($values)";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $last_id;
}



=head2 load_csv_2_db

 Title   : load_csv_2_db
 Usage   : load_csv_2_db( -database => 'name of the database,
                               );

 Function: uses the \COPY psql to load the content of a file
			into the database
			
			any user can run it
 Returns : returns nothing

=cut
sub load_csv_2_db{
	my $database = shift;
	my $host = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $file = shift;
	my $header = shift;
	my $delimiter = shift;
	
	#Opening database if needed 
	my $dbh;
	my $command = "PGPASSWORD='$pass' psql -h $host -U $user $database -c ";
	
	$command .= '"\COPY '."$table FROM '$file'";
	if ($delimiter ne 'NO'){
		$command .= " DELIMITER $delimiter ";	
	}
		$command .= " CSV ";
	if ($header eq 'YES'){
		$command .= " HEADER ";	
	}		
	$command .= '"';
	
	print "Executing command: $command\n";
	try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;
}


=head2 load_csv_2_db_sudo

 Title   : load_csv_2_db_sudo
 Usage   : load_csv_2_db_sudo( -database => 'name of the database,
                               );

 Function: uses the COPY function to load the content of a file
			into the database
			needs you have the privileges 
			
			NEVER TESTED!
			
 Returns : returns nothing

=cut
sub load_csv_2_db_sudo{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $file = shift;
	my $header = shift;
	my $delimiter = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#construct the query with COPY and DELIMITER and HEADER if needed
			$query = '\copy'." $table FROM '$file' ";
			if ($delimiter ne 'NO'){
				$query .= " DELIMITER $delimiter ";	
			}
			$query .= " CSV ";
			if ($header eq 'YES'){
				$query .= " HEADER ";	
			}			
			
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}	

}

=head2 insert_into_table_woconn

 Title   : insert_into_table_woconn
 Usage   : insert_into_table_woconn( -database => 'name of the database,
                               );

 Function: makes a simple insert into the database using the LOCK 
					values may be multiple by separating them with comma
					This function does not activate the database connection hence it needs the 
					database reference is given in input
					
 Returns : returns the id of the sample inserted

=cut
sub insert_into_table_woconn{
	my $dbh = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;

	
	#Opening database if needed 
	my $query;
	my $last_id = -1;
	if (defined $dbh){

		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) VALUES ($values)";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Undefined database connection!\n";}
	
	return $last_id;
}


=head2 update_table_2

 Title   : update_table_2
 Usage   : update_table_2( -database => 'name of the database,
                               );

 Function: Update the table $table the fields $fields and values $values
					 whether the $fields and $values do not exist.
					 In this last case the record is updated
					
						Uses experimentally the comma as a separator
 Returns : 

=cut
sub update_table_2{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;
	#OPTIONAL
	my $fields2 = shift;
	my $values2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		#Prepare the UPDATE clause
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		my $upd_clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$upd_clause .= " ".$fields_arr[$i]."=".$values_arr[$i].",";
			}
			#Remove the last comma
			chop($upd_clause);
		}else{die "ERROR [$?]: Values in fields and values differ in number: \nfields:$fields\nvalues:$values\n" ;}
		#Prepare an additional 
		my $whe_clause = "";
		if (defined $fields2 and defined 	$values2){
			#Prepare the WHERE clause
			@fields_arr = split(",",$fields2);
			@values_arr = split(",",$values2);
			
			if ( scalar(@fields_arr) == scalar(@values_arr) ) {
				for ( my $i = 0; $i < scalar(@fields_arr); $i++){
					$whe_clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
				}
				#Remove the last AND
				chop($whe_clause);
				chop($whe_clause);
				chop($whe_clause);
			}		
		}
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");#LOCK

			#Updates the database
			$query = "UPDATE $table SET $upd_clause ";
			#Add the where clause if needed
			if ( $whe_clause ne '' ){
				$query .= " WHERE $whe_clause";
			}
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
		 
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
}



=head2 update_table

 Title   : update_table
 Usage   : update_table( -database => 'name of the database,
                               );

 Function: Update the table $table the fields $fields and values $values
					 whether the $fields and $values do not exist.
					 In this last case the record is updated

 Returns : 

=cut
sub update_table{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;
	#OPTIONAL
	my $fields2 = shift;
	my $values2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		#Prepare the UPDATE clause
		my @fields_arr = split(";",$fields);
		my @values_arr = split(";",$values);
		my $upd_clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$upd_clause .= " ".$fields_arr[$i]."=".$values_arr[$i].",";
			}
			#Remove the last comma
			chop($upd_clause);
		}else{die "ERROR [$?]: Values in fields and values differ in number: \nfields:$fields\nvalues:$values\n" ;}
		#Prepare an additional 
		my $whe_clause = "";
		if (defined $fields2 and defined 	$values2){
			#Prepare the WHERE clause
			@fields_arr = split(";",$fields2);
			@values_arr = split(";",$values2);
			
			if ( scalar(@fields_arr) == scalar(@values_arr) ) {
				for ( my $i = 0; $i < scalar(@fields_arr); $i++){
					$whe_clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
				}
				#Remove the last AND
				chop($whe_clause);
				chop($whe_clause);
				chop($whe_clause);
			}		
		}
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");#LOCK

			#Updates the database
			$query = "UPDATE $table SET $upd_clause ";
			#Add the where clause if needed
			if ( $whe_clause ne '' ){
				$query .= " WHERE $whe_clause";
			}
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
		 
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
}

=head2 update_table_woconn

 Title   : update_table_woconn
 Usage   : update_table_woconn( -database => 'name of the database,
                               );

 Function: Update the table $table the fields $fields and values $values
					 whether the $fields and $values do not exist.
					 In this last case the record is updated

 Returns : 

=cut
sub update_table_woconn{
	my $dbh = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;
	#OPTIONAL
	my $fields2 = shift;
	my $values2 = shift;
	
	#Opening database if needed 
	my $query;
	if (defined $dbh){
		#Prepare the UPDATE clause
		my @fields_arr = split(";",$fields);
		my @values_arr = split(";",$values);
		my $upd_clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$upd_clause .= " ".$fields_arr[$i]."=".$values_arr[$i].",";
			}
			#Remove the last comma
			chop($upd_clause);
		}else{
			print "ERROR (update_table_woconn): fields and values are of different cardinalities\n";
		}		
				
		#Prepare the WHERE clause
		@fields_arr = split(";",$fields2);
		@values_arr = split(";",$values2);
		my $whe_clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$whe_clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
			}
			#Remove the last AND
			chop($whe_clause);
			chop($whe_clause);
			chop($whe_clause);
		}else{
			print "ERROR (update_table_woconn): fields2 and values2 are of different cardinalities\n";
		}		
			
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");
			#SET client_min_messages=WARNING;
			#Extracts information from database defined.
			#fetch informations from the database
			$query = "UPDATE $table SET $upd_clause WHERE $whe_clause";
			print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			print "Executing query: $query\n";
			$select_handle->execute() or die $select_handle->errstr;
		 
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Undefined database connection!\n";}
	
}

=head2 insert_on_duplicate

 Title   : insert_on_duplicate
 Usage   : insert_on_duplicate( -database => 'name of the database,
                               );

 Function: Inserts into table $table the fields $fields and values $values
					 whether the $fields and $values do not exist.
					 In this last case the record is updated

 Returns : 

=cut
sub insert_on_duplicate{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $fields = shift;
	my $values = shift;
	my $fields2 = shift;
	my $values2 = shift;
		
	#Opening database if needed 
	my $dbh;
	my $query;
	#my $last_id = -1;
	if (defined $database){
				$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 0, AutoCommit => 1 });
		
		#Prepare the WHERE clause
		my @fields_arr = split(",",$fields2);
		my @values_arr = split(",",$values2);
		my $update_str = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$update_str .= " ".$fields_arr[$i]."=VALUES(".$values_arr[$i]."),";
			}
			#Remove the last comma
			chop($update_str);
		}
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) VALUES ($values) ".
					"ON DUPLICATE KEY UPDATE $update_str";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			print "Executing query: $query\n";#DEBUGCODE
			#$select_handle->execute() or die $select_handle->errstr;
			###$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
	$dbh->disconnect();
	#if ( ! defined $last_id ) {
		#$last_id = get_id_if_exists_from_db($database,$dsn,$user,$pass,$table,$id_field,$fields2,$values2);
	#}
	
	}else {print "ERROR: Undefined database $database!\n";}
	
	#return $last_id;
	
}

=head2 insert_if_not_exists_locked

 Title   : insert_if_not_exists_locked
 Usage   : insert_if_not_exists_locked( -database => 'name of the database,
                               );

 Function: Inserts into table $table the fields $fields and values $values
					 only if $fields2 and $values2 do not exist. Then returns the last id.
					 If they exist instead, the $id_field is extracted from the database calling the subroutine
					 get_dbelement_ifexists.

 Returns : returns the id of the sample inserted othervise returns the id of
					the element already present. -1 on error

=cut
sub insert_only_if_not_exists_locked{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $id_field = shift;
	my $fields = shift;
	my $values = shift;
	my $fields2 = shift;
	my $values2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query = "";
	my $last_id = -1;
	my $existing_id = -1;
	if (defined $database){
		#Check fields and values equality
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		if ( scalar(@fields_arr) != scalar(@values_arr) ) {
			print "ERROR in insert_if_not_exists_locked: You are using ".scalar(@fields_arr)." fields and ".scalar(@values_arr)." values ?\n" ;
		}
		#Prepare the WHERE clause
		my @fields_arr2 = split(",",$fields2);
		my @values_arr2 = split(",",$values2);
		my $clause = "";
		if ( scalar(@fields_arr2) == scalar(@values_arr2) ) {
			for ( my $i = 0; $i < scalar(@fields_arr2); $i++){
				$clause .= " ".$fields_arr2[$i]."=".$values_arr2[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 0, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) SELECT * FROM ( SELECT $values) AS tmp ".
					"WHERE NOT EXISTS (SELECT * FROM $table WHERE $clause) RETURNING $id_field";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			#$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$last_id = $res->{$id_field};
					#print "Last inserted id is $last_id\n";
			}
			
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	$dbh->disconnect();
	#Get the ID inserted or if it was already there the clause wil be the same
	if ( $last_id == -1) {
		$existing_id = get_id_if_exists_from_db($database,$dsn,$user,$pass,$table,$id_field,$fields2,$values2);
	}
	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $last_id,$existing_id;
}


=head2 insert_only_if_not_exists_locked_woconn

 Title   : insert_only_if_not_exists_locked_woconn
 Usage   : insert_only_if_not_exists_locked_woconn( -database => 'name of the database,
                               );

 Function: Inserts into table $table the fields $fields and values $values
					 only if $fields2 and $values2 do not exist. Then returns the last id.
					 If they exist instead, the $id_field is extracted from the database calling the subroutine
					 get_dbelement_ifexists_woconn.

					This function does not activate the database connection hence it needs the 
					database reference is given in input
					
 Returns : returns the id of the sample inserted othervise returns the id of
					the element already present. -1 on error

=cut
sub insert_only_if_not_exists_locked_woconn{
	my $dbh = shift;
	my $table = shift;
	my $id_field = shift;
	my $fields = shift;
	my $values = shift;
	my $fields2 = shift;
	my $values2 = shift;
	
	my $query = "";
	my $last_id = -1;
	my $existing_id = -1;
	
	if (defined $dbh){
		#Check fields and values equality
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		if ( scalar(@fields_arr) != scalar(@values_arr) ) {
			print "ERROR in insert_if_not_exists_locked: You are using ".scalar(@fields_arr)." fields and ".scalar(@values_arr)." values ?\n" ;
		}
		#Prepare the WHERE clause
		my @fields_arr2 = split(",",$fields2);
		my @values_arr2 = split(",",$values2);
		my $clause = "";
		if ( scalar(@fields_arr2) == scalar(@values_arr2) ) {
			for ( my $i = 0; $i < scalar(@fields_arr2); $i++){
				$clause .= " ".$fields_arr2[$i]."=".$values_arr2[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}

		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) SELECT * FROM ( SELECT $values) AS tmp ".
					"WHERE NOT EXISTS (SELECT * FROM $table WHERE $clause) RETURNING $id_field";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			#$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$last_id = $res->{$id_field};
					#print "Last inserted id is $last_id\n";
			}
			
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;

	#Get the ID inserted or if it was already there the clause will be the same
	if ( $last_id == -1) {
		$existing_id = get_dbelement_ifexists_woconn($dbh,$table,$id_field,$fields2,$values2);
	}
	
	}else {print "ERROR: Undefined database connection!\n";}
	
	return $last_id,$existing_id;
}

=head2 insert_if_not_exists_locked

 Title   : insert_if_not_exists_locked
 Usage   : insert_if_not_exists_locked( -database => 'name of the database,
                               );

 Function: Inserts into table $table the fields $fields and values $values
					 whether the $fields2 and $values2 do not exist.
					 In this last case the $id_field is extracted from the database

 Returns : returns the id of the sample inserted othervise returns the id of
					the element already present. -1 on error

=cut
sub insert_if_not_exists_locked{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $id_field = shift;
	my $fields = shift;
	my $values = shift;
	my $fields2 = shift;
	my $values2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query = "";
	my $last_id = -1;
	if (defined $database){
		#Check fields and values equality
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		if ( scalar(@fields_arr) != scalar(@values_arr) ) {
			print "ERROR in insert_if_not_exists_locked: You are using ".scalar(@fields_arr)." fields and ".scalar(@values_arr)." values ?\n" ;
		}
		#Prepare the WHERE clause
		my @fields_arr2 = split(",",$fields2);
		my @values_arr2 = split(",",$values2);
		my $clause = "";
		if ( scalar(@fields_arr2) == scalar(@values_arr2) ) {
			for ( my $i = 0; $i < scalar(@fields_arr2); $i++){
				$clause .= " ".$fields_arr2[$i]."=".$values_arr2[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 0, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			$dbh->do("LOCK TABLE $table  IN SHARE MODE");

			#Extracts information from database defined.
			#fetch informations from the database
			$query = "INSERT INTO $table ($fields) SELECT * FROM ( SELECT $values) AS tmp ".
					"WHERE NOT EXISTS (SELECT * FROM $table WHERE $clause)";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			print "Executing query: $query\n";#DEBUGCODE
			$select_handle->execute() or die $select_handle->errstr;
			#$last_id = $dbh->last_insert_id( undef, undef, $table, undef );
			#print "Last inserted id is $last_id\n";
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	$dbh->disconnect();
	#Get the ID inserted or if it was already there the clause wil be the same
	#if ( ! defined $last_id ) {
	$last_id = get_id_if_exists_from_db($database,$dsn,$user,$pass,$table,$id_field,$fields2,$values2);
	#}
	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $last_id;
}

=head2 get_id_if_exists_from_db_like

 Title   : get_id_if_exists_from_db_like
 Usage   : get_id_if_exists_from_db_like( -database => 'name of the database,
                               );

 Function: selects a single field given a condition in a given table
						Uses LIKE instead that =
 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub get_id_if_exists_from_db_like{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $sel_field = shift;
	my $fields = shift;
	my $values = shift;
		
	#Opening database if needed 
	my $dbh;
	my $query;
	my $exists = -1;
	if (defined $database){
		#Prepare the WHERE clause
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		my $clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$clause .= " ".$fields_arr[$i]." LIKE ".$values_arr[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}else{die "ERROR: fields number differs from values number ?\n"; }
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $clause";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$exists = $res->{$sel_field};
			}
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $exists;
}

=head2 get_id_if_exists_from_db

 Title   : get_id_if_exists_from_db
 Usage   : get_id_if_exists_from_db( -database => 'name of the database,
                               );

 Function: selects a single field given a condition in a given table

 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub get_id_if_exists_from_db_woconn{
	my $dbh = shift;
	my $table = shift;
	my $sel_field = shift;
	my $fields = shift;
	my $values = shift;
		
	#Opening database if needed 
	my $query;
	my $exists = -1;
	if (defined $dbh){
		#Prepare the WHERE clause
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		my $clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}else{die "ERROR: fields number differs from values number ?\n"; }

		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $clause";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res->{$sel_field}){
					$exists = $res->{$sel_field};
			}
		 $dbh->do("COMMIT WORK");  
		};
		warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Undefined database !\n";}
	
	return $exists;
}


#=head2 get_hash_from_table

 #Title   : get_hash_from_table
 #Usage   : get_hash_from_table( -database => 'name of the database,
                               #);

 #Function: obtains an hash from a table using a select where the first
		#field is the value and the second field the key

 #Returns : returns the filled hash

#=cut
#sub get_hash_from_table{
	#my 

	##to verify if it must be inserted
	#my $query = "SELECT ".$cfg_hash->{'db_var_id'}.",".$cfg_hash->{'db_var_compid'}." FROM  ".$cfg_hash->{'db_variants_table'}.";";	
	##print_and_log( "Executing: $query\n",$log_file);#DEBUGCODE
	#my $res = fetch_all_rows($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'},$query);
	#foreach my $row (@$res) {
		##Get the row values into an array
		#my @arr = @$row;
		#my $fld_num = 0;
		#my $compid = "";
		#my $varid = -1;
		##Each element of the array is now a field of the table
		#foreach  my $elem (@arr) {
			#if (defined $elem and $fld_num == 0) {
				#$varid = $elem;
			#}
			#if (defined $elem and $fld_num == 1) {
				#$compid = $elem;
			#}	
			#$fld_num++;		
		#}
		##Put the elements into an hash
		#$compids_hash->{$compid} = $varid;
		##print_and_log( "Putting  $varid into hash for: $compid\n",$log_file);#DEBUGCODE
	#}
	
#}
	

=head2 get_id_if_exists_from_db

 Title   : get_id_if_exists_from_db
 Usage   : get_id_if_exists_from_db( -database => 'name of the database,
                               );

 Function: selects a single field given a condition in a given table

 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub get_id_if_exists_from_db{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $sel_field = shift;
	my $fields = shift;
	my $values = shift;
		
	#Opening database if needed 
	my $dbh;
	my $query;
	my $exists = -1;
	if (defined $database){
		#Prepare the WHERE clause
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		my $clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}else{die "ERROR: fields number differs from values number ?\n"; }
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $clause";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$exists = $res->{$sel_field};
			}
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $exists;
}



=head2 get_dbelement_ifexists_woconn

 Title   : get_dbelement_ifexists_woconn
 Usage   : get_dbelement_ifexists_woconn( 
 -database => 'name of the database,
                               );

 Function: selects a single field from the database given a condition in a given table.
					This function does not activate the database connection hence it needs the 
					database reference is given in input
					

 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub get_dbelement_ifexists_woconn{
	my $dbh = shift;
	my $table = shift;
	my $sel_field = shift;
	my $fields = shift;
	my $values = shift;

	my $query = "";
	my $exists = -1;
	if (defined $dbh){
		#Prepare the WHERE clause
		my @fields_arr = split(",",$fields);
		my @values_arr = split(",",$values);
		my $clause = "";
		if ( scalar(@fields_arr) == scalar(@values_arr) ) {
			for ( my $i = 0; $i < scalar(@fields_arr); $i++){
				$clause .= " ".$fields_arr[$i]."=".$values_arr[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}else{die "ERROR: fields number differs from values number ?\n"; }

		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE ACCESS SHARE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $clause";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$exists = $res->{$sel_field};
			}
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;

	}else {print "ERROR: Connection to database looks undefined!\n";}
	
	return $exists;
}


=head2 select_element_from_db

 Title   : select_element_from_db
 Usage   : select_element_from_db( -database => 'name of the database,
                               );

 Function: equal to get_id_if_exists_from_db, selects a single field given 
					fields and values given in a reference to array

 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub select_element_from_db{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $sel_field = shift;
	my $fields = shift;
	my $values = shift;
		
	#Opening database if needed 
	my $dbh;
	my $query;
	my $exists = -1;
	if (defined $database){
		#Prepare the WHERE clause
		#my @fields_arr = split(",",$fields);
		#my @values_arr = split(",",$values);
		my $clause = "";
		if ( scalar(@$fields) == scalar(@$values) ) {
			for ( my $i = 0; $i < scalar(@$fields); $i++){
				$clause .= " ".$fields->[$i]."=".$values->[$i]." AND";
			}
			#Remove the last AND
			chop($clause);
			chop($clause);
			chop($clause);
		}else{die "ERROR: fields number differs from values number ?\n"; }
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $clause";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$exists = $res->{$sel_field};
			}
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $exists;
}


=head2 get_count_from_selected

 Title   : get_count_from_selected
 Usage   : get_count_from_selected( -database => 'name of the database,
                               );

 Function: this function uses the sql function count(*) to obtain a count
						of the results obtained from a SELECT 
					
 Returns : returns the id of the sample fetched, otherwise -1

=cut
sub get_count_from_selected{
	my $dbh = shift;
	my $sel_table = shift;#
	my $where_clause = shift;
	my $log_file = shift;
	
	my $query;
	my $count = 0;
	if (defined $dbh){
		#Here we open a locked access to the database 
		#my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");
			
			#fetch informations from the database
			$query = "SELECT count(*) FROM ( $sel_table	) AS res $where_clause";
			#print "Executing: $query\n";
			#print_and_log("Executing: $query\n",$log_file);
			$count = $dbh->selectrow_array($query, undef);

		  $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	# $select_handle->finish;
	}else {print "ERROR: Undefined database connection!\n";}
	
	return $count;
}

=head2 select_all_from_db

 Title   : select_all_from_db
 Usage   : select_all_from_db( -database => 'name of the database,
                               );

 Function: selects all the sel_field which are output from the query 
					SELECT $sel_field FROM $table WHERE $field='".$value."'"
					Uses fetchall_hashref
					
					Useful when you want to get all the samples belonging to a given group
					Use the reference hash in output with a FOR construct to iter them
					
					//foreach my $group_id ( keys %{$groups_associated}){
					
 Returns : returns an hash reference to the results

=cut
sub select_all_from_db{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $sel_field = shift;
	my $field = shift;
	my $value = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;

	my $res; #Response to use as output
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT $sel_field FROM $table WHERE $field='".$value."'";
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			$res = $select_handle->fetchall_hashref($sel_field);
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $res;
}


=head2 do_query_select_all_woconn

 Title   : do_query_select_all_woconn
 Usage   : do_query_select_all_woconn( -database => 'name of the database,
                               );

 Function: executes a query given in input

 Returns : returns an hash reference to the results

=cut
sub do_query_select_all_woconn{
	my $dbh = shift;
	my $query = shift;
	my $field = shift;


	my $res; #Response to use as output
	if (defined $dbh){
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			$res = $select_handle->fetchall_hashref($field);
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Undefined database!\n";}
	
	return $res;
}


=head2 do_query_select_all

 Title   : do_query_select_all
 Usage   : do_query_select_all( -database => 'name of the database,
                               );

 Function: fetches an entire row from a table and returns an hash containing the requested fields

 Returns : returns an hash reference to the results

=cut
sub do_query_select_all{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $query = shift;
	my $field = shift;
	
	#Opening database if needed 
	my $dbh;

	my $res; #Response to use as output
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			$res = $select_handle->fetchall_hashref($field);
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $res;
}


=head2 db_print_selected_rows

 Title   : db_print_selected_rows
 Usage   : db_print_selected_rows( -database => 'name of the database,
                               );

 Function: Uses selectall_hashref to print a table of results using a primary key id
	that is the first element of $fields. $fields must be separated with comma.
	$fields2 and Values2 are used for the WHERE clause
	
	The function goes through the hash reference of results obtained (rows of the table)
	and prints all the columns contained in $fields

 Returns : returns an array reference to the results

=cut
sub db_print_selected_rows{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $fields = shift;
	my $table = shift;
	my $fields2 = shift;
	my $values2  =shift;
	
	my $results;
	
	
	#Opening database if needed 
	my $dbh;

	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	my $query = "SELECT $fields  FROM $table WHERE  $fields2 = $values2";
		my $select_handle;
		eval{
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			#Separating the fields so that we can take the first to be the primary key for
			#select_all_hashref. The remaining elements will be used to print results
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			my @arr_flds = split(",",$fields);
			#Shift is used to give the first element as the primary key and leaves only the
			#the remaining fields in the array
			$results = $dbh->selectall_hashref($query,shift @arr_flds);
			if ( scalar(keys %{$results}) > 0  ) {
				foreach my $id (keys %$results) {
					#If more fields are used to print it will print all
					if (scalar (@arr_flds) > 0){
						foreach my $fld (@arr_flds){
							print "$id\t".$results->{$id}->{$fld}."\n";	
						}	
					}else{#otherwise only the first will be print
						print "$id\n";
					}
				}
			}else{
					print "WARNING: no result for query $query\n";
			}
			
		};
		warn "Query $query aborted with error: $@\n" if $@;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}

=head2 do_fetch_row_array

 Title   : do_fetch_row_array
 Usage   : do_fetch_row_array( -database => 'name of the database,
                               );

 Function: Executes fetch_row_array 
					It is useful when you want to get two or more fields of a row in a table
					
					Needs you give in input a reference to an array. Set it outside as
					my @res = ();
					Then call the function giving n input the reference \@res
					Finally you will get the fields in the order you gave in the query
					that you must write before to call the function.

 Returns : returns an array reference to the results

=cut
sub do_fetch_row_array{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $query = shift;
	my $ref = shift;
	
	
	#Opening database if needed 
	my $dbh;

	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			@$ref = $select_handle->fetchrow_array;
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}

=head2 do_fetch_row_array_woconn

 Title   : do_fetch_row_array_woconn
 Usage   : do_fetch_row_array_woconn( -database => 'name of the database,
                               );

 Function: Executes fetch_row_array 
					It is useful when you want to get two or more fields of a row in a table
					
					Needs you give in input a reference to an array. Set it outside as
					my @res = ();
					Then call the function giving n input the reference \@res
					Finally you will get the fields in the order you gave in the query
					that you must write before to call the function.

 Returns : returns an array reference to the results

=cut
sub do_fetch_row_array_woconn{
	my $dbh = shift;
	my $query = shift;
	my $ref = shift;
	
	if (defined $dbh){
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			@$ref = $select_handle->fetchrow_array;
			$dbh->do("COMMIT WORK");  
		};
		warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
	}else {print "ERROR: Connection to database looks undefined!\n";}
}

=head2 do_query

 Title   : do_query
 Usage   : do_query( -database => 'name of the database,
                               );

 Function: executes a query given in input

 Returns : returns an hash reference to the results

=cut
sub do_query{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $query = shift;
	
	#Opening database if needed 
	my $dbh;

	my $res; #Response to use as output
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $res;
}

=head2 select_distinct_samples

 Title   : select_distinct_samples
 Usage   : select_distinct_samples( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: selects all the distinct samplesid for a given group

 Returns : returns an hash reference to the results

=cut
sub select_distinct_samples{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $sel_field = shift;
	my $field = shift;
	my $value = shift;
	#OPTIONAL
	my $field2 = shift;
	my $value2 = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;

	my $res; #Response to use as output
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT DISTINCT $sel_field FROM $table WHERE $field=$value";
			#Add second conditon if exists
			if (defined $field2 and defined $value2){
					$query .= " AND $field2=$value2";
			}
			#print "Executing query: $query\n";#DEBUGCODE
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			$res = $select_handle->fetchall_hashref($sel_field);
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted with error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $res;
}

=head2 exists_ss_locked

 Title   : exists_ss_locked
 Usage   : exists_ss_locked( -database => 'name of the database,
                               );

 Function: select

 Returns : returns the id of the sample inserted

=cut
sub exists_ss_locked{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $field = shift;
	my $value = shift;
	
	#Opening database if needed 
	my $dbh;
	my $query;
	my $exists = -1;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");

			#fetch informations from the database
			$query = "SELECT * FROM $table WHERE $field='".$value."'";
			#print "Executing query: $query\n";
			$select_handle = $dbh->prepare($query);
			die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
			$select_handle->execute() or die $select_handle->errstr;
			my $res = $select_handle->fetchrow_hashref;
			if (defined $res){
					$exists = 1;
			}
		 $dbh->do("COMMIT WORK");  
		};
	warn "Query $query aborted wth error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $exists;
}


=head2 getSampleConfiguration_locked

 Title   : getSampleConfiguration_locked
 Usage   : getSampleConfiguration_locked( -database => 'name of the database,
                               );

 Function: gets in input an hash for the sample and returns all the
						variables given in the table

 Returns : returns the inputed hash

=cut
sub getSampleConfiguration_locked{
	my ($in_hash) = shift;
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $table = shift;
	my $field = shift;
	my $sample_id = shift;
	
  #Opening database if needed 
	my $dbh;
	my $query;
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open an access to the database to collect informations
		my $select_handle;
		eval{
		$dbh->do("BEGIN WORK");
		#$dbh->do("LOCK TABLE $table");
		
		#Extracts information from database defined.
    #fetch informations from the database
    $query = "SELECT * FROM $table WHERE $field=$sample_id";
    $select_handle = $dbh->prepare($query);
    die "ERROR [$?]: cannot prepare queries; aborting: ?\n" unless defined $select_handle;
    #print "Executing query: $query\n";
    $select_handle->execute() or die $select_handle->errstr;
    my $res = $select_handle->fetchrow_hashref;
    if (defined $res){
       #Copy the output hash to the inputed one
			foreach my $elem (keys %$res){
					$$in_hash->{$elem} = $res->{$elem};
					#print "$elem: ".$$in_hash->{$elem}."\n";
			}
    }else{ print "There is no result for query: $query\n";
     }
	 
	 $dbh->do("COMMIT WORK");  
	};
	warn "Query $query aborted by error: $@\n" if $@;
	
	 $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
}	


=head2 copyTable2DB

 Title   : copyTable2DB
 Usage   : copyTable2DB( -database => 'name of the database,
                               );

 Function: copies a file table inside a destination table of the database

 Returns : nothing

=cut
sub copyTable2DB{
	  # CONFIG VARIABLE
  my $file_table = shift;
  my $database = shift;
  my $table_name = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $host = shift;
	
	$dsn .= "dbname=$database;";
	 my $dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1 } )
          or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
          
	#Then it will be created
	my $query = '\COPY'." $table_name FROM '$file_table'";
	print("Trying to execute: $query\n");
	#my $rc = $dbh->do($query);
	
	#Executing SQL code from the .sql file
	my $command = qq(PGPASSWORD=$pass psql --username $user --dbname $database --host $host -c '$query');
	print "Executing: $command\n";
	 (system  $command)== 0	or die "ERROR [$?]: an error occurred while executing $query. \n"; 
	if ($dbh->err()) { die "$DBI::errstr\n"; }
	$dbh->disconnect(); 
}



=head2 createDBAndTables

 Title   : createDBAndTables
 Usage   : createDBAndTables( -database => 'name of the database,
                               );

 Function:  creates the database with all the tables

 Returns : nothing

=cut
sub createDBAndTables{
  # CONFIG VARIABLE
  my $database = shift;
	my $sqlFile= shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $host= shift;
	
	my $basic_db = 'postgres';
	
  if ( -e $sqlFile){
    #First the database with that name will be erased
    destroy_db($database,$dsn, $user, $pass);

    print "Creating db ".$database." with user:".$user." and pwd: ".$pass."\n";
    # DATA SOURCE NAME
    #my $dsn = "dbi:$platform:host=$host;port=$port;sslmode=allow;dbname=$dbname;";
    $dsn .= "dbname=$basic_db;";
    my $dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1 } )
              or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
        
    #Then it will be created
    my $rc = $dbh->do("create database $database");
     
    if ($dbh->err()) { die "$DBI::errstr\n"; }
    $dbh->disconnect();
    print "Uploading ".$sqlFile." into ".$database."\n"; 
    
    #Executing SQL code from the .sql file
    (system qq(PGPASSWORD=$pass psql --username $user --dbname $database --host $host --file $sqlFile) )== 0
			or die "ERROR [$?]: an error occurred while importing $sqlFile in $database. \n";
     
    print "DB and tables created!\n";
    }else {die "ERROR [$?]: the file $sqlFile with the SQL commands to create tables is not present.Exiting...\n";}

}

=head2 createTables

 Title   : createTables
 Usage   : createTables( -database => 'name of the database,
                               );

 Function:  creates the tables from an SQL file

 Returns : nothing

=cut
sub createTables{
  # CONFIG VARIABLE
  my $database = shift;
	my $sqlFile= shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $host= shift;

  if ( -e $sqlFile){
     print "Uploading ".$sqlFile." into ".$database."\n"; 
    $dsn .= "dbname=$database;";
    #Executing SQL code from the .sql file
    (system qq(PGPASSWORD=$pass psql --username $user --dbname $database --host $host --file $sqlFile) )== 0
			or die "ERROR [$?]: an error occurred while importing $sqlFile in $database. \n";
     
    print "DB and tables created!\n";
    }else {die "ERROR [$?]: the file $sqlFile with the SQL commands to create tables is not present.Exiting...\n";}
 
  
}

=head2 destroy_db

 Title   : destroy_db
 Usage   : destroy_db( -database => 'name of the database,
                               );

 Function:  destroy a database with all the tables. Here I always use the connection to another database to drop that. 

 Returns : nothing

=cut
sub destroy_db{
    my $database = shift;
    my $dsn = shift;
    my $user = shift;
    my $pass = shift;
    
    my $basic_db = 'postgres';
    
    # DATA SOURCE NAME
    #my $dsn = "dbi:$platform:information_schema:$host:$port";
    $dsn .= "dbname=$basic_db;";
    my $dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1 } )
              or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
    
    #First the database with that name will be erased
    $dbh->do("DROP DATABASE IF EXISTS $database") or die "ERROR: unable to cancel database $database! ?"; 
}


=head2 db_present

 Title   : db_present
 Usage   : db_present( -dbToSearch => 'name of the database);

 Function:  checks if a database exists.

 Returns : 1 if the database exists, 0 otherwise

=cut
sub db_present {
  my $dbToSearch = shift;
  my $dsn = shift;
  my $user = shift;
  my $pass = shift;
  
  #my $basic_db = 'postgres';
  
  $dsn .= "dbname=$dbToSearch;";
  #my $dsn = "dbi:$platform:information_schema:$host:$port";
  print "Searching database $dbToSearch\n";#DEBUGCODE          
  my $dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 0} );
  #my $databases = $dbh->selectcol_arrayref("show databases like '".$dbToSearch."'");#MySQL
  #my $numKeys = @$databases; #MYSQL
  #return $numKeys;
  my $retVal = 0;
	
	if (defined $dbh){
		$retVal =1 ;
		#print "Database $dbToSearch is present\n";#DEBUGCODE
	}
  return $retVal;
}


#####################################HPO

=head2 parse_RVIS_scores

 Title   : parse_RVIS_scores
 Usage   : parse_RVIS_scores( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: parses the given RVIS file and gives this scores to genes in the database.
					Some values may have the 'NA' values. In that case the gene will not be taken.

 Returns : 

=cut
sub parse_RVIS_scores{
	my $cfg_hash = shift;
	my $rvis_file = shift;
	my $log_file = shift;
	
	my $db_genes_table = $cfg_hash->{'db_genes_table'};
	
	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $gene_ind = 0;
	my $rvis_score_ind = 3;
	my $rvis_perc_ind = 4;

	#Check if the file with the mapping was downloaded
  if (not (-e $rvis_file) ){
    die "ERROR: $rvis_file not present...\n";
  }
 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });


	open(GENES,"<".$rvis_file) or die "Cannot open file $rvis_file\n";
	my $not_first_row = 0;
	while  (my $row = <GENES>){
		chomp($row);
		#Header line
		if ($not_first_row == 0){
			my @fld_names = split($sep,$row);
			$not_first_row++;
		}else{
			my @fields = split($sep,$row);
			my $gene = $fields[$gene_ind];
			#Remove all chars which are not:
			$gene =~ s/[^A-Za-z0-9\-\_]//g;
			my $rvis_score = $fields[$rvis_score_ind];
			my $rvis_perc = $fields[$rvis_perc_ind];
			#They should be two real number. No "NA" and no other stupid characters
			#Genes having NA will not be considered
			if (correct_type($rvis_score,"real") and  correct_type($rvis_perc,"real") ) {
							#Get the gene id if exists
			#print_and_log( "Searching for gene:$gene..\n ",$log_file);
			my $gene_id =	get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																				$cfg_hash->{'db_genes_name'},"'".$gene."'");
			#If the gene is already present use an update																	
			if ($gene_id >= 0){
					#UPDATE fields = values
					my $fields = $cfg_hash->{'db_rvis_score'}.";".$cfg_hash->{'db_rvis_perc'};
					my $values = $rvis_score.";".$rvis_perc;
					#WHERE genename= gene 
					my $fields2 = $cfg_hash->{'db_genes_name'};
					my $values2 = "'".$gene."'";
					
					#print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
					#Update
					update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);		

				}
			#If the gene is not present insert it..	
				else{
					#print_and_log("$gene is not present. Inserting...\n",$log_file);
					#INSERT genename,gdiscore VALUES gene,gdiscore
					my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_gdi_score'}.",".$cfg_hash->{'db_rvis_perc'};
					my $values = "'".$gene."',$rvis_score,$rvis_perc";
																	
					#Inserts the gene id into the table only if that gene does not exists
					insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);
				}
			}
		}
	}	
	close(GENES);
	#Disconnect db
  $dbh->disconnect(); 		
}


=head2 download_RefSeq_2_Genes

 Title   : download_RefSeq_2_Genes
 Usage   : download_RefSeq_2_Genes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: parses the given GDI file and gives this scores to genes in the database

 Returns : 

=cut
sub download_RefSeq_2_Genes{
	my $cfg_hash = shift;
	my $log_file = shift;
	
	#The file will be downloaded in the main DATA folder
	my $ref2gene_f = $cfg_hash->{'main_data_folder'}."/".$cfg_hash->{'refseq2genes_f'};
	my $ucsc_user = $cfg_hash->{'ucsc_user'};
	my $ucsc_host = $cfg_hash->{'ucsc_host'};
	my $ann_build_ver = $cfg_hash->{'ann_build_ver'};
	my $ucsc_refseq2genes_query =  $cfg_hash->{'ucsc_refseq2genes_query'};
	
	#Executing a query to UCSC genome browser
	my $command = "mysql --user=$ucsc_user -N --host=$ucsc_host -A -D $ann_build_ver -e '$ucsc_refseq2genes_query' > $ref2gene_f";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
}



=head2 parse_GDI_scores

 Title   : parse_GDI_scores
 Usage   : parse_GDI_scores( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: parses the given GDI file and gives this scores to genes in the database

 Returns : 

=cut
sub parse_GDI_scores{
	my $cfg_hash = shift;
	my $gdi_file = shift;
	my $log_file = shift;
	
	my $db_genes_table = $cfg_hash->{'db_genes_table'};
	#field into the database for the GDI score
	my $gdi_field = $cfg_hash->{'db_gdi_score'};
	
	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $gene_ind = 0;
	my $gdi_score_ind = 1;
	my $gdi_phred_ind = 2;
	
	#Check if the file with the mapping was downloaded
  if (not (-e $gdi_file) ){
    die "ERROR: $gdi_file not present...\n";
  }

 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });


	open(GENES,"<".$gdi_file) or die "Cannot open file $gdi_file\n";
	while  (my $row = <GENES>){
		chomp($row);
		#Header line
		if ($row =~ /^Gene/){
			my @fld_names = split($sep,$row);
		}else{
			my @fields = split($sep,$row);
			my $gene = $fields[$gene_ind];
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
			my $gdi_score = $fields[$gdi_score_ind];
			
			#Get the gene id if exists
			#print_and_log( "Searching for gene:$gene..\n ",$log_file);

			my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
										$cfg_hash->{'db_genes_name'},"'".$gene."'");
			#If the gene is already present use an update																	
			if ($gene_id >= 0){
					#UPDATE fields = values
					my $fields = $cfg_hash->{'db_gdi_score'};
					my $values = $gdi_score;
					#WHERE genename= gene 
					my $fields2 = $cfg_hash->{'db_genes_name'};
					my $values2 = "'".$gene."'";
					
					#print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
					#Update
					#update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							#$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
					update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
				}
			#If the gene is not present insert it..	
				else{
					#print_and_log("$gene is not presente. Inserting...\n",$log_file);#DEBUGCODE
					#INSERT genename,gdiscore VALUES gene,gdiscore
					my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_gdi_score'};
					my $values = "'".$gene."',$gdi_score";
																	
					#Inserts the gene id into the table only if that gene does not exists
					#insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
							#$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$fields,$values);
					insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);			
				}
		}
	}	
	close(GENES);
	
	#Disconnect db
  $dbh->disconnect(); 
}



=head2 parse_genes2RefSeq

 Title   : parse_genes2RefSeq
 Usage   : parse_genes2RefSeq( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function:  grep '^9606' | cut -f1,2,4,16
					grep takes the human species
					cut takes the taxonomy id, entrez id, refseq transcript, gene symbol
					
					constructs both a table with trnascripts associated to genes and
					also associates the entrez id to the genes.
					The gene hash will be filled with the concatenation of transcripts
					and with the entrez id (if it has not been inserted yet)

 Returns : returns an hash reference to the results

=cut
sub parse_genes2RefSeq{
	my $cfg_hash = shift;
	my $refseq2genes_f = shift;
	my $log_file = shift;

	#Separator to use for the table
	my $sep  = "\t";

	my $tax_col = $cfg_hash->{'g2ref_tax_col'};
	my $entrez_col = $cfg_hash->{'g2ref_entrez_col'};
	my $refseq_col = $cfg_hash->{'g2ref_refseq_col'};
	my $gensymb_col = $cfg_hash->{'g2ref_gensymb_col'};

	#Check if the file with the mapping was downloaded
  if (not (-e $refseq2genes_f) ){
    die "ERROR: $refseq2genes_f not present...\n";
  }

	my $human_only_f = $refseq2genes_f.".human";
	#########For debug comment the following lines and..
	#Apply a grep to select only those lines related with human
	my $command = 'grep '."'^".$cfg_hash->{'human_tax_id'}.'\s'."' $refseq2genes_f > $human_only_f";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
				
	#Now filter the file to get one only with the needed columns
	#Entrez id\tTranscript description\tGene\n
	my $filtered_f = $human_only_f.".filt";
	my $needed_cols = "$entrez_col,$refseq_col,$gensymb_col";
	extract_columns_from_file($human_only_f,$needed_cols,$filtered_f);
	
	#...uncomment these for debug of less...
	###############################################
	#my $filtered_f = $human_only_f.".filt_RED";#DEBUGCODE
	###############################################
		
	#Change now the indexes because I just filtered the file 
	$entrez_col = 0;
	$refseq_col = 1;
	$gensymb_col = 2;
		
	#An hash for the genes
	my $genes_hash;
	my $tr_hash;
	
 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });

	
	###Inserting refseq ids in a table for transcripts
	print_and_log("Reading from filtered mapping file $filtered_f ..\n",$log_file);
	open(REF2GEN,"<".$filtered_f) or die "ERROR: Cannot open file $filtered_f\n";
	my $fields1 = $cfg_hash->{'db_refseq_id'};
	my $fields2 = $cfg_hash->{'db_refseq_id'};
	my $line_num=1;
	while  (my $row = <REF2GEN>){
		chop($row);
		#print_and_log("Row $row\n",$log_file);
		my @fields = split($sep,$row);
		#If refseqid exists (empty value is '-') ..
		if ( length($fields[$refseq_col]) > 1 ){
			#RefSeq ids could contain a '.1': NM_006983.1
			my $refseq_id = (split(/\./,$fields[$refseq_col]))[0];
			my $values1 = "'".$refseq_id."'";
			my $values2 = "'".$refseq_id."'";
			#Put the refseq information in the transcripts table if the id is not present
			#If it is already there, it takes its identifier
			my ($exist_tr_id, $new_tr_id) = insert_only_if_not_exists_locked_woconn($dbh,$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_transcript_id'},
								$fields1,$values1,$fields2,$values2); 	
			#I will insert the RefSeq id only if it is not present in the DB
			#print_and_log("Row $line_num: Transcript ".$fields[$refseq_col]." with id  $transcript_id...\n",$log_file);#DEBUGCODE
			
			#Get the transcript id selecting from the one that has been filled from the insert function
			my $transcript_id = -1;
			if ($exist_tr_id >=0 ){
				$transcript_id = $exist_tr_id;
			}elsif ( $new_tr_id >= 0){
				$transcript_id = $new_tr_id;
			}else{
				print_and_log("WARNING: cannot find $refseq_id in ".$cfg_hash->{'db_transcripts_table'}."..\n",$log_file);
			}
				
			#Put the RefSeq symbols
			if ( !(defined($genes_hash->{$fields[$gensymb_col]}->{'rs'})) ){
				$genes_hash->{$fields[$gensymb_col]}->{'rs'} = $transcript_id;
				$genes_hash->{$fields[$gensymb_col]}->{'rs_str'} = $fields[$refseq_col];
				#print_and_log("New value for ".$fields[$gensymb_col]." \n",$log_file);#DEBUGCODE
			}else{
				my $rs_temp = $fields[$refseq_col];
				if (scalar( grep {/\b$rs_temp\b/} split(",",$genes_hash->{$fields[$gensymb_col]}->{'rs_str'}) ) == 0){
					#print_and_log("I can add $transcript_id to ".$fields[$gensymb_col]." . $rs_temp was not added..\n",$log_file);#DEBUGCODE
					$genes_hash->{$fields[$gensymb_col]}->{'rs'} .= ",".$transcript_id;
					$genes_hash->{$fields[$gensymb_col]}->{'rs_str'} .= ",".$fields[$refseq_col];
				}
			}					
			#Put the genes associated to the transcript
			#GENETICS: Each transcripts is present only in one gene
			if ( $transcript_id >= 0 ){
				my $gene_temp = $fields[$gensymb_col];
				
				if ( !(defined($tr_hash->{$transcript_id})) ){
					$tr_hash->{$transcript_id} = $gene_temp;
				}else{
					#If the gene has not been concatenated yet..
					if (scalar( grep {/\b$gene_temp\b/} split(",",$tr_hash->{$transcript_id}) ) == 0){
						#print_and_log("I can add gene to $transcript_id . $gene_temp was not added..\n",$log_file);#DEBUGCODE
						$tr_hash->{$transcript_id} .= ",".$gene_temp;
					}#else{print_and_log("$gene_temp already present among ".$tr_hash->{$transcript_id}."..\n",$log_file);}#DEBUGCODE
				}				
			}
		}
		#Put the Entrez id, if present
		if ( (! defined $genes_hash->{$fields[$gensymb_col]}->{'ent'}) and ($fields[$entrez_col] ne '-') ){
				#print_and_log("Adding also entrez".$fields[$entrez_col]."..\n",$log_file);#DEBUGCODE
				$genes_hash->{$fields[$gensymb_col]}->{'ent'} = $fields[$entrez_col];
		}
		$line_num++;
	}
	close(REF2GEN);
	#Remove the temporary filtered files
	#delete_file($filtered_f);
	#delete_file($human_only_f);
	
	#Update transcripts table with Genes associated
	print_and_log("Finished copying. Now updating table ".$cfg_hash->{'db_transcripts_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	foreach my $tr_id (keys %{$tr_hash} ){
		if ( $tr_id >= 0 ){
			my $fields = $cfg_hash->{'db_tr_genes_ids'}; 
			my $values = "'".$tr_hash->{$tr_id}."'";
			#WHERE genename = gene
			my $fields2 = $cfg_hash->{'db_transcript_id'};
			my $values2 = $tr_id;
			#Update
			update_table_woconn($dbh,$cfg_hash->{'db_transcripts_table'},$fields,$values,$fields2,$values2);		
		}
	}
	
	print_and_log("Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id, the refseq ids and the entrez id in the genes table
	my $genes_added = 0;
	my $genes_updated = 0;
	
	foreach my $gene (keys %{$genes_hash} ){
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");	
		#If there is at least an information for the gene
		if ( (defined ($genes_hash->{$gene}->{'rs'})) or (defined ($genes_hash->{$gene}->{'ent'})) ){
			#If the gene is already present use an update																	
			if ($gene_id >= 0){
					#UPDATE fields = values
					my $fields = ""; 
					my $values = "";
					if (defined ($genes_hash->{$gene}->{'rs'}) ){
						$fields .= $cfg_hash->{'db_refseq_ids'}.";";
						$values .= "'".$genes_hash->{$gene}->{'rs'}."';";
					}
					if (defined ($genes_hash->{$gene}->{'ent'}) ){
						$fields .= $cfg_hash->{'db_entrez_id'}.";";
						$values .= "'".$genes_hash->{$gene}->{'ent'}."';";
					}
					#Remove last comma
					chop($fields);
					chop($values);
					#WHERE genename = gene
					my $fields2 = $cfg_hash->{'db_genes_id'};
					my $values2 = $gene_id;
					#print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
					#Update
					update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
					$genes_updated++;
			}
			#If the gene is not present insert it..	
			else{
					#INSERT genename,refseqids,entrezid VALUES gene,ids
					my $fields = $cfg_hash->{'db_genes_name'}.",";
					my $values = "'".$gene."',";

					if (defined $genes_hash->{$gene}->{'rs'}){
						$fields .= $cfg_hash->{'db_refseq_ids'}.",";
						$values .= "'".$genes_hash->{$gene}->{'rs'}."',";
					}
					if (defined $genes_hash->{$gene}->{'ent'}){
						$fields .= $cfg_hash->{'db_entrez_id'}.",";
						$values .= "'".$genes_hash->{$gene}->{'ent'}."',";
					}
					#Remove last comma
					chop($fields);
					chop($values);																				
					#Inserts the gene id into the table only if that gene does not exists
					insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);
					$genes_added++;
				}
		}
	}
	print_and_log("Finished updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	print_and_log("Added  $genes_added genes and updated $genes_updated...\n",$log_file);
}


=head2 parse_genes2RefSeq

 Title   : parse_genes2RefSeq
 Usage   : parse_genes2RefSeq( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function:  grep '^9606' | cut -f1,2,4,16
					grep takes the human species
					cut takes the taxonomy id, entrez id, refseq transcript, gene symbol
					
					constructs both a table with trnascripts associated to genes and
					also associates the entrez id to the genes.
					The gene hash will be filled with the concatenation of transcripts
					and with the entrez id (if it has not been inserted yet)

 Returns : returns an hash reference to the results

=cut
sub parse_genes2RefSeq_wo_tr{
	my $cfg_hash = shift;
	my $refseq2genes_f = shift;
	my $log_file = shift;

	#Separator to use for the table
	my $sep  = "\t";

	my $tax_col = $cfg_hash->{'g2ref_tax_col'};
	my $entrez_col = $cfg_hash->{'g2ref_entrez_col'};
	my $refseq_col = $cfg_hash->{'g2ref_refseq_col'};
	my $gensymb_col = $cfg_hash->{'g2ref_gensymb_col'};

	#Check if the file with the mapping was downloaded
  if (not (-e $refseq2genes_f) ){
    die "ERROR: $refseq2genes_f not present...\n";
  }

	my $human_only_f = $refseq2genes_f.".human";
	#Apply a grep to select only those lines related with human
	my $command = 'grep '."'^".$cfg_hash->{'human_tax_id'}.'\s'."' $refseq2genes_f > $human_only_f";
	print_and_log( "Executing command: $command\n",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";
				
	#Now filter the file to get one only with the needed columns
	my $filtered_f = $human_only_f.".filt";
	my $needed_cols = "$entrez_col,$refseq_col,$gensymb_col";
	extract_columns_from_file($human_only_f,$needed_cols,$filtered_f);
	
	#New indexes
	$entrez_col = 0;
	$refseq_col = 1;
	$gensymb_col = 2;
		
	#An hash for the genes
	my $genes_hash;
	
	###Inserting refseq ids in a table for transcripts
	print_and_log("Reading from filtered mapping file $filtered_f ..\n",$log_file);
	open(REF2GEN,"<".$filtered_f) or die "ERROR: Cannot open file $filtered_f\n";
	my $fields1 = $cfg_hash->{'db_refseq_id'};
	my $fields2 = $cfg_hash->{'db_refseq_id'};
	my $line_num=1;
	while  (my $row = <REF2GEN>){
		chop($row);
		my @fields = split($sep,$row);
		#If refseqid exists..
		if ( length($fields[$refseq_col]) > 1 ){
			my $refseq_id = $fields[$refseq_col];
			my $values1 = "'".$refseq_id."'";
			my $values2 = "'".$refseq_id."'";
			#Put the refseq information in the transcripts table if the id is not present
			#If it is already there, it takes its identifier
			my $transcript_id = insert_if_not_exists_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
								$cfg_hash->{'db_pass'},$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_transcript_id'}
								,$fields1,$values1,$fields2,$values2); 	
			#I will insert the RefSeq id only if it is not present in the DB
			#print_and_log("Row $line_num: Transcript ".$fields[$refseq_col]." with id  $transcript_id...\n",$log_file);
				
				#Put the RefSeq symbols
				if ( !(defined($genes_hash->{$fields[$gensymb_col]}->{'rs'})) ){
					$genes_hash->{$fields[$gensymb_col]}->{'rs'} = $transcript_id;
					$genes_hash->{$fields[$gensymb_col]}->{'rs_str'} = $fields[$refseq_col];
				}else{
					my $rs_temp = $fields[$refseq_col];
					if (scalar( grep {/\b$rs_temp\b/} split(",",$genes_hash->{$fields[$gensymb_col]}->{'rs_str'}) ) == 0){
						$genes_hash->{$fields[$gensymb_col]}->{'rs'} .= ",".$transcript_id;
						$genes_hash->{$fields[$gensymb_col]}->{'rs_str'} .= ",".$fields[$refseq_col];
					}
				}								
		}
		#Put the Entrez id, if present
		if ( (! defined $genes_hash->{$fields[$gensymb_col]}->{'ent'}) and ($fields[$entrez_col] ne '-') ){
			#	print_and_log("Adding also entrez".$fields[$entrez_col]."..\n",$log_file);
				$genes_hash->{$fields[$gensymb_col]}->{'ent'} = $fields[$entrez_col];
		}
		$line_num++;
	}
	close(REF2GEN);
	#Remove the temporary filtered files
	#delete_file($filtered_f);
	#delete_file($human_only_f);
	
	print_and_log("Finished copying. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id, the refseq ids and the entrez id in the genes table
	my $genes_added = 0;
	my $genes_updated = 0;
	
	foreach my $gene (keys %{$genes_hash} ){
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		#Get the gene id if exists
		my $gene_id = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
																			$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");	
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = ""; 
				my $values = "";
				if (defined $genes_hash->{$gene}->{'rs'}){
					$fields .= $cfg_hash->{'db_refseq_ids'}.";";
					$values .= "'".$genes_hash->{$gene}->{'rs'}."';";
				}
				if (defined $genes_hash->{$gene}->{'ent'}){
					$fields .= $cfg_hash->{'db_entrez_id'}.";";
					$values .= "'".$genes_hash->{$gene}->{'ent'}."';";
				}
				#Remove last comma
				chop($fields);
				chop($values);
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_id'};
				my $values2 = $gene_id;
				#print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
				#Update
				update_table($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
				$genes_updated++;
		}
		#If the gene is not present insert it..	
			else{
				#INSERT genename,refseqids,entrezid VALUES gene,ids
				my $fields = $cfg_hash->{'db_genes_name'}.",";
				my $values = "'".$gene."',";

				if (defined $genes_hash->{$gene}->{'rs'}){
					$fields .= $cfg_hash->{'db_refseq_ids'}.",";
					$values .= "'".$genes_hash->{$gene}->{'rs'}."',";
				}
				if (defined $genes_hash->{$gene}->{'ent'}){
					$fields .= $cfg_hash->{'db_entrez_id'}.",";
					$values .= "'".$genes_hash->{$gene}->{'ent'}."',";
				}
				#Remove last comma
				chop($fields);
				chop($values);																				
				#Inserts the gene id into the table only if that gene does not exists
				insert_into_table_locked($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
						$cfg_hash->{'db_pass'},$cfg_hash->{'db_genes_table'},$fields,$values);
				$genes_added++;
			}
	}
	print_and_log("Finished updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	print_and_log("Added  $genes_added genes and updated $genes_updated...\n",$log_file);
}


=head2 parse_RefSeq_2_Genes

 Title   : parse_RefSeq_2_Genes
 Usage   : parse_RefSeq_2_Genes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: downloads the subjext file using the link given and parses it to 
					create the genes and phenotypes table

 Returns : returns an hash reference to the results

=cut
sub parse_RefSeq_2_Genes{
	my $cfg_hash = shift;
	my $refseq2genes_f = shift;
	my $log_file = shift;

	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $refseq_ind = 0;
	my $gene_ind = 1;
	
	
	#Check if the file with the mapping was downloaded
  if (not (-e $refseq2genes_f) ){
    die "ERROR: $refseq2genes_f not present...\n";
  }

 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });

	#An hash for the genes
	my $genes_hash;

	#An array to store the inserted phenotypes
	print_and_log("Reading from $refseq2genes_f ..\n",$log_file);

	open(REF2GEN,"<".$refseq2genes_f) or die "Cannot open file $refseq2genes_f\n";
	my $fields1 = $cfg_hash->{'db_refseq_id'};
	my $fields2 = $cfg_hash->{'db_refseq_id'};
	while  (my $row = <REF2GEN>){
	 chomp($row);
			my @fields = split($sep,$row);
			my $refseq_id = $fields[$refseq_ind];
			my $values1 = "'".$refseq_id."'";
			my $values2 = "'".$refseq_id."'";
			#Put the refseq information in the transcripts table if the id is not present
			my $refseq_db_id = insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_transcripts_table'},$cfg_hash->{'db_refseq_id'}
								,$fields1,$values1,$fields2,$values2); 										
			if (! defined $genes_hash->{$fields[$gene_ind]}){
				$genes_hash->{$fields[$gene_ind]} = $refseq_db_id;
			}else{
				$genes_hash->{$fields[$gene_ind]} .= ",".$refseq_db_id;
			}
	}
	close(REF2GEN);

	print_and_log("Finished copying. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id and the corresponding set of HPO ids in the
	#genes table

	foreach my $gene (keys %{$genes_hash} ){
			#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																		$cfg_hash->{'db_genes_name'},"'".$gene."'");	
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = $cfg_hash->{'db_refseq_ids'};
				my $values = "'".$genes_hash->{$gene}."'";
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_name'};
				my $values2 = "'".$gene."'";
				#print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
				#Update
				update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);

			}
		#If the gene is not present insert it..	
			else{
				#INSERT genename,hpoids VALUES gene,ids
				my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_refseq_ids'};
				my $values = "'".$gene."','".$genes_hash->{$gene}."'";
																
				#Inserts the gene id into the table only if that gene does not exists
				insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);
			}
	}
	#Disconnect db
  $dbh->disconnect(); 			
}


=head2 parse_OMIM_2_Genes

 Title   : parse_OMIM_2_Genes
 Usage   : parse_OMIM_2_Genes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: downloads the subjext file using the link given and parses it to 
					create the genes and phenotypes table

 Returns : returns an hash reference to the results

=cut
sub parse_OMIM_2_Genes{
	my $cfg_hash = shift;
	my $omim2genes_f = shift;
	my $log_file = shift;

	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $omim_ind = 5;
	my $mim_entry_ind = 1;
	my $gene_ind = 6;
	my $omim_gene_desc_ind = 7;
	my $omim_phenotype_ind = 12;
	
	#Check if the file with the mapping was downloaded
  if (not (-e $omim2genes_f) ){
	die "ERROR: $omim2genes_f not present...\n";
  }

	#An hash for the genes
	my $genes_hash;

	#Get and save the associated ids
	print_and_log("Reading from $omim2genes_f ..\n",$log_file);

	open(OMIM2GEN,"<".$omim2genes_f) or die "Cannot open file $omim2genes_f\n";
	while  (my $row = <OMIM2GEN>){
		#Jump all lines staring with #
		next unless $row !~ /^#/;
		chop($row); 
		my @fields = split($sep,$row);
		#Read only the rows with the 'gene' indication as MIM entry type
			if ( defined $fields[$gene_ind] ){
				if (  $fields[$gene_ind] ne '') {
					my $gene_field = $fields[$gene_ind];
					my $omim_id = "-";
					my $omim_gene_desc = "-";
					my $omim_phenotype = "-"; 
					if (defined $fields[$omim_ind] and $fields[$omim_ind] ne '') {$omim_id = $fields[$omim_ind];}
					if (defined $fields[$omim_gene_desc_ind] and $fields[$omim_gene_desc_ind] ne '') {$omim_gene_desc = $fields[$omim_gene_desc_ind];}
					if (defined $fields[$omim_phenotype_ind] and $fields[$omim_phenotype_ind] ne '') {$omim_phenotype = $fields[$omim_phenotype_ind];}
					#Remove all semicolon from description which can corrupt the query
					$omim_phenotype =~ s/;//g;
					$omim_gene_desc =~ s/;//g;
					
					my @genes = split(",",$gene_field);
					foreach my $gene (@genes){
						if (! defined $genes_hash->{$gene}){
							$genes_hash->{$gene}->{'omim_id'} = $omim_id;
							$genes_hash->{$gene}->{'omim_gene_desc'} = $omim_gene_desc;
							$genes_hash->{$gene}->{'omim_phenotype'} = $omim_phenotype;
						}else{
							$genes_hash->{$gene}->{'omim_id'} .= ",".$omim_id;
							$genes_hash->{$gene}->{'omim_gene_desc'} .= ",".$omim_gene_desc;
							$genes_hash->{$gene}->{'omim_phenotype'} .= ",".$omim_phenotype;
						}						
					}

				}				
			}
	}
	close(OMIM2GEN);

	print_and_log("Finished parsing and storing information in memory. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	
	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
	#Now insert the gene id and the corresponding set of ids in the
	#genes table
	my $genes_added = 0;
	my $genes_updated = 0;
	foreach my $gene (keys %{$genes_hash} ){
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");																				
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = $cfg_hash->{'db_omim_ids'}.";".$cfg_hash->{'db_omim_genedesc'}.";".$cfg_hash->{'db_omim_phenotype'};
				my $values = "'".$genes_hash->{$gene}->{'omim_id'}."';'".$genes_hash->{$gene}->{'omim_gene_desc'}."';'".$genes_hash->{$gene}->{'omim_phenotype'}."'";
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_name'};
				my $values2 = "'".$gene."'";
				print_and_log("$gene already present. Updating with fields = $fields and values= $values...\n",$log_file);#DEBUGCODE
				#Update
				update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);		
				$genes_updated++;
		}
		#If the gene is not present insert it..	
		else{
			#INSERT genename,hpoids VALUES gene,ids
			my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_omim_ids'}.",".$cfg_hash->{'db_omim_genedesc'}.",".$cfg_hash->{'db_omim_phenotype'};
			my $values = "'".$gene."','".$genes_hash->{$gene}->{'omim_id'}."','".$genes_hash->{$gene}->{'omim_gene_desc'}."','".$genes_hash->{$gene}->{'omim_phenotype'}."'";
															
			#Inserts the gene id into the table only if that gene does not exists
			insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);		
			$genes_added++;
		}
	}
	print_and_log("Finished updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	print_and_log("Added  $genes_added genes and updated $genes_updated...\n",$log_file);

	#Disconnect db
  $dbh->disconnect(); 	
}

=head2 parse_OMIM_2_Genes

 Title   : parse_OMIM_2_Genes
 Usage   : parse_OMIM_2_Genes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: downloads the subjext file using the link given and parses it to 
					create the genes and phenotypes table

 Returns : returns an hash reference to the results

=cut
sub parse_OMIM_2_GenesOLD{
	my $cfg_hash = shift;
	my $omim2genes_f = shift;
	my $log_file = shift;

	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $omim_ind = 0;
	my $mim_entry_ind = 1;
	my $gene_ind = 3;
	
	
	#Check if the file with the mapping was downloaded
  if (not (-e $omim2genes_f) ){
	die "ERROR: $omim2genes_f not present...\n";
  }

	#An hash for the genes
	my $genes_hash;

 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });


	#Get and save the associated ids
	print_and_log("Reading from $omim2genes_f ..\n",$log_file);

	open(OMIM2GEN,"<".$omim2genes_f) or die "Cannot open file $omim2genes_f\n";
	while  (my $row = <OMIM2GEN>){
		#Jump all lines staring with #
		next unless $row !~ /^#/;
		chop($row); 
		my @fields = split($sep,$row);
		#Read only the rows with the 'gene' indication as MIM entry type
		if ( $fields[$mim_entry_ind] =~ /gene/) {
			if ( defined $fields[$gene_ind] ){
				if (  $fields[$gene_ind] ne '') {
					my $omim_id = $fields[$omim_ind];
					if (! defined $genes_hash->{$fields[$gene_ind]}){
						$genes_hash->{$fields[$gene_ind]} = $omim_id;
					}else{
						$genes_hash->{$fields[$gene_ind]} .= ",".$omim_id;
					}
				}				
			}

		}
	}
	close(OMIM2GEN);

	print_and_log("Finished copying. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id and the corresponding set of ids in the
	#genes table
	my $genes_added = 0;
	my $genes_updated = 0;
	foreach my $gene (keys %{$genes_hash} ){
			#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");																				
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = $cfg_hash->{'db_omim_ids'};
				my $values = "'".$genes_hash->{$gene}."'";
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_name'};
				my $values2 = "'".$gene."'";
				print_and_log("$gene already present. Updating...\n",$log_file);#DEBUGCODE
				#Update
				update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);		
				$genes_updated++;
			}
		#If the gene is not present insert it..	
			else{
				#INSERT genename,hpoids VALUES gene,ids
				my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_omim_ids'};
				my $values = "'".$gene."','".$genes_hash->{$gene}."'";
																
				#Inserts the gene id into the table only if that gene does not exists
				insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);		
				$genes_added++;
			}
	}
	print_and_log("Finished updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	print_and_log("Added  $genes_added genes and updated $genes_updated...\n",$log_file);

	#Disconnect db
  $dbh->disconnect(); 	
}


=head2 createGOTable

 Title   : createGOTable
 Usage   : createGOTable( -GOTermsLink => link to a file with GO terms and descriptions,
                          -uniprotWebUser, uniprotWebPass  =>  user-id and password to give to Uniprot to downloadd,
                          -tableHash => reference to an hash to fill with GO ids from the database
                          - go_create => use this variable to ask to the function if the table should be created or not 
			       );

 Function: creates the GOTable with each of the GO terms, definitions and divisions. It  parses a file with lines
          like this: GO:0000050\tGO:0006594\sGO:0006871\t urea cycle\tP\n


 Returns : nothing

=cut
sub createHPOTable{
  #Input Variables
  my $cfg_hash = shift;
  my $termsFileLink  = shift;
  my ($tableHash) =shift;
	my $dbDataFolder = shift;
  my $create = shift;
  my $log_file = shift;
  
  my $termsFile = extract_name($termsFileLink,"0");
  my $termstable = $cfg_hash->{'db_phenotypes_table'}; #Name of the goTable
  my $termstablePath =  $dbDataFolder."/".$termstable;
    
  #Download file if it doesn't exists
  if (not (-e $dbDataFolder."/".$termsFile) ){
    #print "\n Downloading File: ".$termsFileLink."\n";
    download_file($termsFileLink,$dbDataFolder);
  }else{
    print_and_log( "File ".$dbDataFolder."/".$termsFile." already downloaded...\n",$log_file);
  }
    
  open(FILE,$dbDataFolder."/".$termsFile);
  if ($create == 1){
    open(TABLE,">$termstablePath");
  }
  	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });
	
	   
  my $sep = ": ";
  #print "...creating a file and an hash with GO information...";
  my $term_id = "";
  my $termname = "";
  my $def = "";
  my @alt_ids = ();

  my $getflag = 0;
  print_and_log( "Starting reading ".$dbDataFolder."/".$termsFile." ...\n",$log_file);
  #read the   file
  while(my $row = <FILE> ) {
	 chop($row);
    #The first lines are descriptions and they start with '?'
    if ($row =~ /^\[Term/){
      #If the flag was active before, then print info
      if ( $getflag == 1 ){

				#One term is always there so... take it!
        if ($create == 1){
					#Remove any ' or " from the name
					$termname =~ s/[\'\"\,]//g;
					$def =~ s/[\'\"\,]//g;
					my $fields = $cfg_hash->{'db_hpo_id'}.",".$cfg_hash->{'db_hpo_name'}.",".$cfg_hash->{'db_hpo_def'};
					my $values = "'$term_id','$termname','$def'";												
					#print_and_log( "Inserting $fields = $values\n",$log_file);#DEBUGCODE
					##Inserts the gene id into the table only if that gene does not exists
					my $phenid = insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$fields,$values);
					if ($phenid >= 0 ){
							#Put the db id into the hash
							$$tableHash->{$term_id} = $phenid;
					}	
        }else{
					# get the db id of the term type
					my $phenid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
									$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_phenotypes_table'},
									$cfg_hash->{'db_hpo_id'},$term_id);			
					#Put the db id into the hash
					$$tableHash->{$term_id} = $phenid;						
				} 
								
				if ( scalar(@alt_ids) > 0){
						 foreach my $alt_id (@alt_ids){
							 #print " - ".$key." - ";#DEBUGCODE
							 if ($create == 1){
								#Remove any ' or " from the name
								$termname =~ s/[\'\"\,]//g;
								$def =~ s/[\'\"\,]//g;
								my $fields = $cfg_hash->{'db_hpo_id'}.",".$cfg_hash->{'db_hpo_name'}.",".$cfg_hash->{'db_hpo_def'};
								my $values = "'$term_id','$termname','$def'";		
								#print_and_log( "Inserting $fields = $values\n",$log_file);#DEBUGCODE										
								##Inserts the gene id into the table only if that gene does not exists
								my $phenid = insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$fields,$values);
								if ($phenid >=0 ){
										#Put the db id into the hash
										$$tableHash->{$term_id} = $phenid;
								}	
							 }else{
									# get the db id of the term type
									my $phenid = get_id_if_exists_from_db($cfg_hash->{'db_name'},$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
													$cfg_hash->{'db_pass'},$cfg_hash->{'db_analyses_table'},$cfg_hash->{'db_phenotypes_table'},
													$cfg_hash->{'db_hpo_id'},$term_id);			
									#Put the db id into the hash
									$$tableHash->{$term_id} = $phenid;								 
							 }
							 #print "\n";#DEBUGCODE
							}			
				}        
				#print "\n";#DEBUGCODE  
				#Re-initialize all.. you never know...
				$term_id = "";
				$termname = "";
				$def = "";
				@alt_ids = ();    

				$getflag = 0;
			}
			$getflag = 1;
    }
    #When the getflag is 1, means the row [Term] has been met
    if ($getflag == 1) {
			if ($row =~ /^id:/){
				my @strings = split ($sep, $row);
				$term_id = $strings[1];
			}

			if ($row =~ /^name:/){
				my @strings = split ($sep, $row);
				$termname = $strings[1];
			}	

			if ($row =~ /^def:/){
				my @strings = split ($sep, $row);
				$def = $strings[1];		
			}
			
			if ($row =~ /^alt_id:/){
				my @strings = split ($sep, $row);
				push(@alt_ids,$strings[1]);
			}
		}
	
  }
    
  close(FILE);
  if ($create == 1){
    close(TABLE);
  }
  
 	$dbh->disconnect();  

  if ($cfg_hash->{'remove_temp'} eq "YES"){
    delete_file($dbDataFolder."/".$termsFile);
  }
  #print "....DONE!(TABLE Upload completed)\n ";
}

=head2 parse_HPO_genes_to_phenotypes

 Title   : parse_HPO_genes_to_phenotypes
 Usage   : parse_HPO_genes_to_phenotypes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: downloads the subjext file using the link given and parses it to 
					create the genes and phenotypes table

 Returns : returns an hash reference to the results

=cut
sub parse_HPO_genes_to_phenotypes{
	my $cfg_hash = shift;
	my ($hpoTableHash) = shift;
	my $gen2phen_file = shift;
	my $log_file = shift;
	
	#A table and a file for phenotypes
	my $ph_table = $cfg_hash->{'db_phenotypes_table'};#"phenotypes";	
	my $db_genes_table = $cfg_hash->{'db_genes_table'};
	#field into the database for the HPO terms
	my $hpos_field = $cfg_hash->{'db_hpo_ids'};
	

	#Separator in the header
	my $h_sep = "<tab>";
	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $gene_ind = $cfg_hash->{'hpo2gen_gene_col'};
	my $hpo_name_ind = $cfg_hash->{'hpo2gen_hponame_col'};
	my $hpo_id_ind = $cfg_hash->{'hpo2gen_hpoid_col'};
	
	#Check if the file with the mapping was downloaded
  if (not (-e $gen2phen_file) ){
    die "ERROR: $gen2phen_file not present...\n";
  }
	#An hash for the HPO ids
	my $hpo_hash;
	#An hash for the genes
	my $genes_hash;

 # print join "\n", keys %INC;#DEBUGCODE
 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });

	#An array to store the inserted phenotypes
	print_and_log("Reading from $gen2phen_file and getting phenotypes associated with genes..\n",$log_file);
	
	my @ph_inserted = ();
	open(GEN2PHEN,"<".$gen2phen_file) or die "Cannot open file $gen2phen_file\n";
	open(PHEN,">".$ph_table) or die "Cannot open file $ph_table\n";
	my $fields1 = $cfg_hash->{'db_hpo_id'}.",".$cfg_hash->{'db_hpo_name'};
	my $fields2 = $cfg_hash->{'db_hpo_id'};
	while  (my $row = <GEN2PHEN>){
	 chomp($row);
	 if ($row =~ /^#Format/){
			$row =~ s/^#Format//;
			my @fld_names = split($h_sep,$row);
	 }else{
			my @fields = split($sep,$row);
			#Get the hpo id
			my $hpo_id = $fields[$hpo_id_ind];
			
			#Get the phenid from the hash
			my $phen_id = -1;
			if ( defined $$hpoTableHash->{$hpo_id} ){
				$phen_id = $$hpoTableHash->{$hpo_id};
			}else{
				print_and_log("WARNING: cannot find $hpo_id in ".$cfg_hash->{'db_phenotypes_table'}."..\n",$log_file);
			}
			if( $phen_id >= 0){
				if (! defined $genes_hash->{$fields[$gene_ind]}){
					#CLPP->HP:0007360
					$genes_hash->{$fields[$gene_ind]} = $phen_id;
				}else{
					#CLPP->HP:0007360,HP:0000130,HP:0001284
					$genes_hash->{$fields[$gene_ind]} .= ",".$phen_id;
				}
			}
	 }
	}
	close(GEN2PHEN);
	close(PHEN);
	
	print_and_log("Finished parsing. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id and the corresponding set of HPO ids in the
	#genes table

	foreach my $gene (keys %{$genes_hash} ){
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");	
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = $cfg_hash->{'db_hpo_ids'};
				my $values = "'".$genes_hash->{$gene}."'";
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_name'};
				my $values2 = "'".$gene."'";
				print_and_log("$gene already present. Updating...\n",$log_file);
				#Update
				update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
			}
		#If the gene is not present insert it..	
			else{
				#INSERT genename,hpoids VALUES gene,ids
				my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_hpo_ids'};
				my $values = "'".$gene."','".$genes_hash->{$gene}."'";
																
				##Inserts the gene id into the table only if that gene does not exists
				insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);		
			}
	}
	#Disconnect db
  $dbh->disconnect(); 
	
}


=head2 parse_HPO_genes_to_phenotypes

 Title   : parse_HPO_genes_to_phenotypes
 Usage   : parse_HPO_genes_to_phenotypes( -database => 'name of the database,
																		-group => 'group wanted'
                               );

 Function: downloads the subjext file using the link given and parses it to 
					create the genes and phenotypes table

 Returns : returns an hash reference to the results

=cut
sub parse_HPO_genes_to_phenotypesOLD{
	my $cfg_hash = shift;
	my $gen2phen_file = shift;
	my $log_file = shift;
	
	#A table and a file for phenotypes
	my $ph_table = $cfg_hash->{'db_phenotypes_table'};#"phenotypes";	
	my $db_genes_table = $cfg_hash->{'db_genes_table'};
	#field into the database for the HPO terms
	my $hpos_field = $cfg_hash->{'db_hpo_ids'};
	

	#Separator in the header
	my $h_sep = "<tab>";
	#Separator to use for the table
	my $sep  = "\t";
	#Indexes of the fields
	my $gene_ind = $cfg_hash->{'hpo2gen_gene_col'};
	my $hpo_name_ind = $cfg_hash->{'hpo2gen_hponame_col'};
	my $hpo_id_ind = $cfg_hash->{'hpo2gen_hpoid_col'};
	
	#Check if the file with the mapping was downloaded
  if (not (-e $gen2phen_file) ){
    die "ERROR: $gen2phen_file not present...\n";
  }
	#An hash for the HPO ids
	my $hpo_hash;
	#An hash for the genes
	my $genes_hash;

 # print join "\n", keys %INC;#DEBUGCODE
 	#################DB
	#Connect to database
	my $dsn = $cfg_hash->{'db_dsn'}."dbname=".$cfg_hash->{'db_name'}.";";
	# PERL DBI CONNECT
	my $dbh = DBI->connect($dsn, $cfg_hash->{'db_user'}, $cfg_hash->{'db_pass'}, { RaiseError => 1, AutoCommit => 1 });

	#An array to store the inserted phenotypes
	print_and_log("Reading from $gen2phen_file and filling table $ph_table..\n",$log_file);
	
	my @ph_inserted = ();
	open(GEN2PHEN,"<".$gen2phen_file) or die "Cannot open file $gen2phen_file\n";
	open(PHEN,">".$ph_table) or die "Cannot open file $ph_table\n";
	my $fields1 = $cfg_hash->{'db_hpo_id'}.",".$cfg_hash->{'db_hpo_name'};
	my $fields2 = $cfg_hash->{'db_hpo_id'};
	while  (my $row = <GEN2PHEN>){
	 chomp($row);
	 if ($row =~ /^#Format/){
			$row =~ s/^#Format//;
			my @fld_names = split($h_sep,$row);
	 }else{
			my @fields = split($sep,$row);
			#Get the hpo id
			my $hpo_id = $fields[$hpo_id_ind];
			#Remove any ' or " from the name
			$fields[$hpo_name_ind] =~ s/[\'\"\,]//g;
			my $values1 = "'".$hpo_id."','".$fields[$hpo_name_ind]."'";
			my $values2 = "'".$hpo_id."'";
			
			#Take the db phenotype id both if exists or if the element will be new
			my ($exist_phen_id, $new_phen_id) = insert_only_if_not_exists_locked_woconn($dbh,$cfg_hash->{'db_phenotypes_table'},$cfg_hash->{'db_phenotype_id'},
								$fields1,$values1,$fields2,$values2); 
			#HP:0007360->Aplasia/Hypoplasia of the cerebellum
			#$hpo_hash->{$fields[$hpo_id_ind]} = $fields[$hpo_name_ind];
			my $phen_id = -1;
			if ($exist_phen_id >=0 ){
				$phen_id = $exist_phen_id;
			}elsif ( $new_phen_id >= 0){
				$phen_id = $new_phen_id;
			}else{
				print_and_log("WARNING: cannot find $hpo_id in ".$cfg_hash->{'db_phenotypes_table'}."..\n",$log_file);
			}
			if( $phen_id >= 0){
				if (! defined $genes_hash->{$fields[$gene_ind]}){
					#CLPP->HP:0007360
					$genes_hash->{$fields[$gene_ind]} = $phen_id;
				}else{
					#CLPP->HP:0007360,HP:0000130,HP:0001284
					$genes_hash->{$fields[$gene_ind]} .= ",".$phen_id;
				}
			}
	 }
	}
	close(GEN2PHEN);
	close(PHEN);
	
	#print_and_log("Finished parsing. Now copying table $ph_table_f in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Load the table PHEN into the database
	#copyTable2DB($ph_table_f,$cfg_hash->{'db_name'},$ph_table,$cfg_hash->{'db_dsn'},$cfg_hash->{'db_user'},
					#$cfg_hash->{'db_pass'},$cfg_hash->{'host'});

	print_and_log("Finished copying. Now updating table ".$cfg_hash->{'db_genes_table'}." in db ".$cfg_hash->{'db_name'}."...\n",$log_file);
	#Now insert the gene id and the corresponding set of HPO ids in the
	#genes table

	foreach my $gene (keys %{$genes_hash} ){
		#Remove all chars which are not:
		$gene =~ s/[^A-Za-z0-9\-\_]//g;
		#Get the gene id if exists
		my $gene_id = get_dbelement_ifexists_woconn($dbh,$cfg_hash->{'db_genes_table'},$cfg_hash->{'db_genes_id'},
																			$cfg_hash->{'db_genes_name'},"'".$gene."'");	
		#If the gene is already present use an update																	
		if ($gene_id >= 0){
				#UPDATE fields = values
				my $fields = $cfg_hash->{'db_hpo_ids'};
				my $values = "'".$genes_hash->{$gene}."'";
				#WHERE genename = gene
				my $fields2 = $cfg_hash->{'db_genes_name'};
				my $values2 = "'".$gene."'";
				print_and_log("$gene already present. Updating...\n",$log_file);
				#Update
				update_table_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values,$fields2,$values2);
			}
		#If the gene is not present insert it..	
			else{
				#INSERT genename,hpoids VALUES gene,ids
				my $fields = $cfg_hash->{'db_genes_name'}.",".$cfg_hash->{'db_hpo_ids'};
				my $values = "'".$gene."','".$genes_hash->{$gene}."'";
																
				##Inserts the gene id into the table only if that gene does not exists
				insert_into_table_locked_woconn($dbh,$cfg_hash->{'db_genes_table'},$fields,$values);		
			}
	}
	#Disconnect db
  $dbh->disconnect(); 
	
}


=head2 fetch_all_rows_woconn

 Title   : fetch_all_rows_woconn
 Usage   : fetch_all_rows_woconn( -database => 'name of the database,
                               );

 Function: gets all rows from a table. Useful when you want to get all rows
					output of a given select.

 Returns : returns an hash reference to the results

=cut
sub fetch_all_rows_woconn{
	my $dbh = shift;
	my $query = shift;
	
	my $res; #Response to use as output
	if (defined $dbh){
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");
			#print  "Executing: $query\n";
			$res = $dbh->selectall_arrayref($query);
			$dbh->do("COMMIT WORK");  
		};
		warn "Query $query aborted with error: $@\n" if $@;
	}else {print "ERROR: Undefined database !\n";}
	
	return $res;
}

=head2 fetch_all_rows

 Title   : fetch_all_rows
 Usage   : fetch_all_rows( -database => 'name of the database,
                               );

 Function: gets all rows from a table. Useful when you want to get all rows
					output of a given select.

 Returns : returns an hash reference to the results

=cut
sub fetch_all_rows{
	my $database = shift;
	my $dsn = shift;
	my $user = shift;
	my $pass = shift;
	my $query = shift;
	
	#Opening database if needed 
	my $dbh;

	my $res; #Response to use as output
	if (defined $database){
		$dsn .= "dbname=$database;";
		# PERL DBI CONNECT
		$dbh = DBI->connect($dsn, $user, $pass, { RaiseError => 1, AutoCommit => 1 });
	
		#Here we open a locked access to the database 
		my $select_handle;
		eval{
			$dbh->do("BEGIN WORK");
			#$dbh->do("LOCK TABLE $table");
			#print  "Executing: $query\n";
			$res = $dbh->selectall_arrayref($query);
			$dbh->do("COMMIT WORK");  
		};
		warn "Query $query aborted with error: $@\n" if $@;
	
	# $select_handle->finish;
   $dbh->disconnect(); 	
	}else {print "ERROR: Undefined database $database!\n";}
	
	return $res;
}
#BEGIN { print (("  " x $main::x++) . "finish dbmanagement compile\n") }#DEBUGCODE

1;
