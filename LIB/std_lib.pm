
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
    
    
#######################utilities.pm###########################################################
# std_lib.pm - A module that contains a series of utils subroutine                              #
################################################################################################

package LIB::std_lib;

BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( print_array print_hash correct_type);
}

#GENERAL UTILITY MODULES
use strict;
use warnings;



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

=head2 print_array

 Title   : print_array
 Usage   : print_array(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints simply the values inside the array in input. Give it as \@array. USes the join function.
 Returns : nothing

=cut
sub print_array{
  my $array = shift;

  print join ' ', @$array, "\n";
  
 }

=head2 print_hash

 Title   : print_hash
 Usage   : print_hash(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints the first two levels of an hash in input. 
 Returns : nothing

=cut
sub print_hashOLD{
  my $hash = shift;
  
	foreach my $key ( keys %{$hash}){
		#print $key.": ".$hash->{$key}."\n";			
		while (my ($intkey, $value) = each %{ $hash->{$key} } ) {
			print "$key->$intkey = ".$hash->{$key}->{$intkey}." \n";
		}
	}
}


=head2 print_hash

 Title   : print_hash
 Usage   : print_hash(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints the first three level of an hash in input. 
 Returns : nothing

=cut
sub print_hash{
  my $hash = shift;
  
	foreach my $key ( keys %{$hash}){
		#print $key.": ".$hash->{$key}."\n";			
		while (my ($intkey, $value) = each %{ $hash->{$key} } ) {
			if (ref($value) eq 'HASH'){
					while (my ($intkey2, $value2) = each %{ $hash->{$key}->{$intkey}} ) {
						print "$key->$intkey-> $intkey2 = ".$hash->{$key}->{$intkey}->{$intkey2}." \n";
					}
			}else{
				print "$key->$intkey = ".$hash->{$key}->{$intkey}." \n";
			}
		}
	}
}

=head2 rev_compl

 Title   : rev_compl
 Usage   : rev_compl(  - dna = the dna sequence to reverse complement
                               );

 Function:  makes reverse complement
 Returns : returns the reverse complement of the dna seq given in input

=cut
sub rev_compl {
	my $dna = shift;
	
	my $rc = reverse $dna;
	$rc =~ tr/ACGTacgt/TGCAtgca/;
	return $rc;
}
1;
