#!/usr/bin/perl

#Libraries
use strict;
use warnings;
####################
#
#This script is used to launch a command. The command can be given in input
#as a variable, or in input could have been given a text file containing the command
#For example for R commands I prefer to write the R CMD BATCH into a file because
# the ' char corrupts the passage of parameters
# 
#N.B. The file should contain a uniq command on a single line!

#Get the file with the parameters for the module
my $command = $ENV{COMMAND};

#Check if the command is a path to a file
my $file = "";
if ( -e $command ){
	open(FILE,"<$command");
	$command = <FILE>;
	close(FILE);
}
#Execute the command
unless ( (system $command) == 0) {
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
}



