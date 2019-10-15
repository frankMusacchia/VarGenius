#!/usr/bin/perl
#To work at CINECA

my $job1 = $ARGV[0];
my $jobn = $ARGV[1];

for (my $i = $job1; $i <= $jobn; $i++){
	print "Deleting job $i..\n";	
	system ("qdel $i");
}
