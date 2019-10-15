#!/opt/software/perl/bin/perl
#To work at TIGEM


use strict;
use warnings;

my $usage = "\n\tUSAGE: perl $0 [function_to_Execute] [bed] [parameter]\n\t ".
						"\n\t1. get start and end	\n".
						"\n\t2. separate genes\n".
						"\n\t3. get start and end exons\n".
						"\n\t4. remove duplicates based on coverage\n".
						"\n\t5. remove duplicates based on intervals\n".
						"\n\t6. separate the rows based on one column with multiple elements\n".
						"\n\t7. do a CUT on a given column. (Change the separator in the code. Def: ,)\n".
						"\n\t8. get consecutive heterozygous snps using VarGenius output\n".
						"\n\t9. remove duplicates rows when the field you give in input is identical\n".
						"\n\t10. remove variants closer than the threshold you give in input\n".
						"\n\t11. put on same line the output of VARIANTSONGENE if is the same variant and same analysis\n".
						"\n\t12. separate_ALT_variants\n".
						"\n\t13. Get extended regions (give: bedfile dimension) \n".
						"\n\t14. Generate BEDGRAPH files from the output of XHMM (zscores) \n".
						"\n\t15. Generate BED files from the output of ExomeDepth\n".
						"\n\t16. Separate Autosomes and Sex Chromosomes\n ".
						"\n\t17. Separate fields of a file containing a single column chr19:54801928-54803753_DUP into a BED\n ".
						"\n\t18. Merge duplicate rows based on X,Y where X is the field to compare and Y the field to add on end of line\n ";
						
						
my $function = $ARGV[0];#function to execute						
my $bed = $ARGV[1];#BED file is the first argument
my $param1 = $ARGV[2];#A parameter give in input

die $usage."\ntoo less parameters\n" unless scalar(@ARGV) >= 2;#If arguments are less than 2, die
die $usage."\nbed $bed does not exist\n" unless (-e $bed);#If bed file does not exist
die $usage."\nnot valid option\n" unless ($function < 19 );#Check the option

#Build a name with 'slice' in the end
my $outname = "$bed";
$outname =~ s/\..+$//;#erase the .bed extension from the name


if ($function eq '1'){
	$outname .= '_st-end.bed';#
	get_start_end($bed,$outname,$param1);
}
if ($function eq '2'){
	$outname .= '_sepgenes.bed';#add _slice
	separate_genes($bed,$outname);
}

if ($function eq '3'){
	$outname .= '_st-end.bed';#
	get_start_end_exons($bed,$outname);
}


if ( $function eq '4' ){
	$outname .= '_uniq.bed';#
	remove_duplicates_based_on_coverage($bed,$outname,$param1);
}

if ( $function eq '5' ){
	$outname .= '_uniq.bed';#
	remove_duplicates_intervals($bed,$outname);
}

if ( $function eq '6' ){
	$outname .= '_seprows.bed';#
	separate_rows($bed,$outname,$param1);
}

if ( $function eq '7' ){
	$outname .= '_cleancol.bed';#
	operate_on_col($bed,$outname,$param1);
}

if ( $function eq '8' ){
	$outname .= '_het_cons.bed';#
	get_cons_het_snps_2_bed($bed,$outname,$param1);
}

if ( $function eq '9' ){
	$outname .= '_uniq.bed';#
	remove_duplicates_rows($bed,$outname,$param1);
}

if ( $function eq '10' ){
	$outname .= '_dist.bed';#
	remove_close_variants($bed,$outname,$param1);
}

if ( $function eq '11' ){
	$outname .= '_2analysis.txt';#
	who_has_variant_row_2_analysis($bed,$outname);
}


if ( $function eq '12' ){
	$outname .= '_sepvar.vcf';#
	separate_ALT_variants($bed,$outname);
}

if ( $function eq '13' ){
	$outname .= '_extended.bed';#
	get_extended_regions($bed,$outname,$param1);
}
if ( $function eq '14' ){
	XHMM_zscores_2_bedgraph($bed);
}
if ( $function eq '15' ){
	ExomeDepth_intervals_2_bed($bed);
}
if ( $function eq '16' ){
	separate_bed_auto_and_sex_chr($bed,$outname);
}
if ( $function eq '17' ){
	compid_intervals_2_bed($bed,$outname);
}
if ( $function eq '18' ){
	$outname .= '_merged.txt';
	merge_rows_with_same_col($bed,$outname,$param1);
}
if ( $function eq '19' ){
	cut_SAM_sequences_using_CIGAR($bed,$param1);
}

# For each interval into the bed file prints a new file with two intervals 
#  1. the one with the start of the interval (plus some nucleotides if added)
#  2. the one with the end of the interval (plus some nucleotides if added)

sub get_start_end {

	my $bed = shift;
	my $outname = shift;
	my $param1 = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		print NEWFILE $fields[0]."\t".($fields[1]-$param1)."\t".($fields[1]+$param1)."\t".$fields[3]."_st\n";
		print NEWFILE $fields[0]."\t".($fields[2]-$param1)."\t".($fields[2]+$param1)."\t".$fields[3]."_end\n";
	}

	close(FILE);
	close(NEWFILE);
	
}


# For each interval into the bed file prints a new file with two intervals 
# the one with only the start of the regions and the other with the end of the region
# 
sub get_start_end_exons {

	my $bed = shift;
	my $outname = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		print NEWFILE $fields[0]."\t".($fields[1]+1)."\t".($fields[1]+1)."\t".$fields[3]."_st\n";
		print NEWFILE $fields[0]."\t".($fields[2])."\t".($fields[2])."\t".$fields[3]."_end\n";
	}

	close(FILE);
	close(NEWFILE);
	
}


###SEPARATE the rows based on one column with multiple elemnts
##Example:
## chr1 9292 299923 GeneA,GeneB,GeneC
#becomes:
## chr1 9292 299923 GeneA
## chr1 9292 299923 GeneB
## chr1 9292 299923 GeneC
sub separate_rows{
	
	my $bed = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	####################SEPARATE the rows based on one column with multiple elemnts
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		#Get the column to sep
		my $col2Sep = $fields[$col2Sep_ind - 1];
		my @elems = split(",",$col2Sep);	
		my $newline = "";
		foreach my $elem (@elems){
			#if (scalar(@elems) > 1){print "$elem \t";}
			for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
				$newline .= $fields[$col]."\t";
			}
			$newline .= $elem."\t";
			for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
				$newline .= $fields[$col]."\t";
			}			
			chop($newline);
			$newline .= "\n";
		}
		print NEWFILE $newline;
		#if (scalar(@elems) > 1){print "\n";}
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
}


###Given the output from VarGenius searches
#the intervals where 2 or more het snps are consecutives
#takes the start and the end of each of this intervals
sub get_cons_het_snps_2_bed{
	my $vgout = shift;
	my $outname = shift;
	my $sample_name = shift;
	
	print "Writing the bed for consecutive het snps into $outname\n";
	
	my $chr_ind = 0;
	my $pos_ind = 1;
	my $gt_field = $sample_name."_GT";
	#print "GT field $gt_field\n";
	my $gt_field_ind = -1;
	####################
	open (FILE,"<$vgout") or die "Cannot open $vgout\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	my $gt_count = 0;
	my $last_het_pos = -1;
	my $chr2keep = "";
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		if ($numline == 0){
			my $numcol = 0;
			foreach my $field (@fields){
					if ( $field eq $gt_field){
						$gt_field_ind = $numcol;
					}
					$numcol++;
			}
		}else{

			my $gt_val = $fields[$gt_field_ind];
			#print "Checking field $gt_field_ind :$gt_val\n";
			#In case of HET
			if ( $gt_val eq '0/1' or  $gt_val eq '1/0' or $gt_val eq '0|1' or  $gt_val eq '1|0'){
				#print "Het SNP at line $numline:$gt_val\n";
				$gt_count++;
				#Initialize the chromosome when is the first HET
				if ($gt_count == 1){
					$chr2keep = $fields[$chr_ind];
				}
				#If this HET is the same chromosome of the one before
				if ( $fields[$chr_ind] eq $chr2keep){
					#If it is the first HET, start to write position
					if ($gt_count == 1){
						$chr2keep = $fields[$chr_ind];
						print NEWFILE $fields[$chr_ind]."\t".$fields[$pos_ind]."\t";
					}					
				}
				#If this HET is on a new chromosome..
				else{
					print NEWFILE $last_het_pos."\t".$gt_count."\n".$fields[$chr_ind]."\t".$fields[$pos_ind]."\t";
					$chr2keep = $fields[$chr_ind];
					$gt_count = 1;
				}
				#Mark the last HET position
				$last_het_pos = $fields[$pos_ind];
			}
			#In case of HOM
			else{
				#If there is at least one het snp, print the end
				if ($gt_count >= 1){
						print NEWFILE $last_het_pos."\t".$gt_count."\n";
						$gt_count = 0;
				}				
			}
		}
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
}


#Like a CUT on a defined column. The separator used for the cut is $sep
sub operate_on_col{
	
	my $bed = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	my $sep = ",";
	####################SEPARATE the rows based on one column with multiple elemnts
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		#Get the column to sep
		my $col2Sep = $fields[$col2Sep_ind - 1];
		my @elems = split(/$sep/,$col2Sep);	
		my $newline = "";

		for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
			$newline .= $fields[$col]."\t";
		}
		$newline .= $elems[0]."\t";
		for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
			$newline .= $fields[$col]."\t";
		}			
		chop($newline);
		$newline .= "\n";
		
		#print $newline;
		print NEWFILE $newline;
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
}

=head2 merge_rows_with_same_col

 Title   : merge_rows_with_same_col
 Usage   : merge_rows_with_same_col( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: given a file and an index of a column, merges the rows where the first 
			column is the same and separates with a comma the "index" columns
			The 	field2merge is always added with comma at the end of line	
 Returns : nothing

=cut
sub merge_rows_with_same_col{
	my $input = shift;
	my $outname = shift;
	my $fields = shift;
	
	my @fields = split(",",$fields);
	my $field4comp = $fields[0];
	my $field2merge = $fields[1];
	
	open (FILE,"<$input") or die "Cannot open $input\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $compid1 = $fields1[$field4comp];
		my $compid2 = $fields2[$field4comp];

		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $lines[$i];
		while( $compid1 eq $compid2 and ($i < ($totlines-2)) ){
			
			
			#Go ahed with next line
			$i++;				
			print "Line $i: $compid1 and $compid2 is same information.\n";
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$compid2 = $fields2[$field4comp];
			
			#Add the column at the end of lne
			$row_2_write .= ",".$fields2[$field2merge];
		}
		$row_2_write .= "\n";
		print NEWFILE $row_2_write
		
	}
	close(FILE);
	close(NEWFILE);
}

=head2 separate_elements_on_col

 Title   : separate_elements_on_col
 Usage   : separate_elements_on_col( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: given a file and an index of a column, separates the elements
						separated with comma and writes a line for each element
					
 Returns : nothing

=cut
sub separate_elements_on_col{
	my $file = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	print  "Separating $file into $outname\n";
	####################SEPARATE elements separated with comma
	open (FILE,"<".$file) or die "Cannot open $file\n";
	open(NEWFILE,">".$outname) or die "Cannot open $outname\n";
	print  "Column $col2Sep_ind\n";
	
	while (my $line = <FILE>){
		print  "Column 2 $col2Sep_ind: $line\n";
		chop($line);
		my @fields = split("\t",$line);
		
		#If the field contains elements separated by comma
		if ($fields[$col2Sep_ind] =~ /,/){
			
			#Get index of the column wanted
			my @elements = split(",",$fields[$col2Sep_ind]);	
							
			foreach my $element  (@elements){
				my $newline = "";
				for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
					$newline .= $fields[$col]."\t";
				}
				$newline .= $element."\t";
				for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
					$newline .= $fields[$col]."\t";
				}			
				chop($newline);
				$newline .= "\n";
				print NEWFILE $newline;
			}			
		}else{
				print NEWFILE $line."\n";
		}
	}
	print "Exiting\n";
	close(FILE);
	close(NEWFILE);	
}


sub separate_bed_auto_and_sex_chr{
	
	my $bed = shift;
	my $outname = shift;
	
	my $sexbed = $outname.".sex.bed";
	my $autobed = $outname.".auto.bed";
	
	#SEX Chromosomes
	#Put chrX
	my $command = "grep '^chrX' $bed > $sexbed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	#append chrY
	$command = "grep '^chrY' $bed >> $sexbed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
	
	#Autosomes
	#Remove chrX
	$command = "grep -v '^chrX' $bed > $autobed.temp";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";

	#Remove chrY
	$command = "grep -v '^chrY' $autobed.temp > $autobed";			 
	print "Executing command: $command\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";
		
}

#When in a BED file we have the interval and gene names falling in the same interval 
#separated by comma. this function can generate a new bed file with replicated
#intervals and single genes
sub separate_genes{
	
	my $bed = shift;
	my $outname = shift;
	my $gene_col = shift;
	####################SEPARATE GENES
	#The last column must be the one of the genes
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		#Get index of the last column, the one with genes indication
		my $lastcol = pop(@fields);
		my @genes = split(",",$lastcol);	
		my $newline = join("\t",@fields);
		
		foreach my $gene  (@genes){
			if (length($gene) > 1){
				print NEWFILE $newline."\t".$gene."\n";
		  }
		}

	}
	close(FILE);
	close(NEWFILE);	
}

sub separate_ALT_variants{
	
	my $input = shift;
	my $outname = shift;
	
	####################SEPARATE VARIANTS
	#The last column must be the one of the genes
	open (FILE,"<$input") or die "Cannot open $input\n";
	open(NEWFILE,">$outname");
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		
		#Get index of the last column, the one with genes indication
		my $lastcol = pop(@fields);
		my @alt_vars = split(",",$lastcol);	
		my $newline = join("\t",@fields);
		
		foreach my $alt_var  (@alt_vars){
			if (length($alt_var) > 0){
				print NEWFILE $newline."\t".$alt_var."\n";
		  }
		}

	}
	close(FILE);
	close(NEWFILE);	
}


# Given a bed file with the first three lines giving the chromosome, the start
# and the end. If there are identical consecutive intervals keeps only the one
# with a given field greater than the other
sub remove_duplicates_based_on_coverage {

	my $bed = shift;
	my $outname = shift;
	my $field4comp = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $compid1 = $fields1[0]."_".$fields1[1]."_".$fields1[2];
		my $compid2 = $fields2[0]."_".$fields2[1]."_".$fields2[2];

		my @exon_pieces1 = split("_",$fields1[3]);
		my $transcript1 = $exon_pieces1[0]."_".$exon_pieces1[1];
		my $exon1 = $exon_pieces1[2]."_".$exon_pieces1[3]."_".$exon_pieces1[4];
		my $stend1 = pop(@exon_pieces1);
		my $highercov = $fields1[$field4comp];
		my $gene1 = $fields1[4];
						
		my @exon_pieces2 = split("_",$fields2[3]);
		my $transcript2 = $exon_pieces2[0]."_".$exon_pieces2[1];
		my $exon2 = $exon_pieces2[2]."_".$exon_pieces2[3]."_".$exon_pieces2[4];
		my $stend2 = pop(@exon_pieces2);
		my $cov2 = $fields2[$field4comp];			
		my $gene2 = $fields2[4];
		my $coverage2 = $fields2[$field4comp];
				
		my $tableid1 = $compid1."_".$transcript1."_".$exon1;
		my $tableid2 = $compid2."_".$transcript2."_".$exon2;
		
		###
		my $row_2_write = $compid1."_$transcript1\_$exon1\_$stend1\_$gene1\t$highercov\n";
		while( $tableid1 eq $tableid2 ){
			print "Line 1 and Line2 with same information. Comparing :".$fields1[$field4comp]." and ".$fields2[$field4comp]."\n";
			#..then match the $field4comp and if the first is greater..
			if ( $highercov < $cov2 ){
				my $gene = $fields2[4];
				my $coverage = $fields2[$field4comp];
				$row_2_write = $compid2."_$transcript2\_$exon2\_$stend2\_$gene2\t$coverage2\n";
				$tableid1 = $tableid2;
				$highercov =  $cov2;
			}
			
			#Go ahed with next line
			$i++;				
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			
			$compid2 = $fields2[0]."_".$fields2[1]."_".$fields2[2];

			@exon_pieces2 = split("_",$fields2[3]);
			$transcript2 = $exon_pieces2[0]."_".$exon_pieces2[1];
			$exon2 = $exon_pieces2[2]."_".$exon_pieces2[3]."_".$exon_pieces2[4];
			$stend2 = pop(@exon_pieces2);
			$cov2 = $fields2[$field4comp];		
			
			$tableid2 = $compid2."_".$transcript2."_".$exon2;			
		}
		print NEWFILE $row_2_write
		
		##If chr start end are equal in this and the next line...
		#if ( $tableid1 eq $tableid2 ){
			
			#print "Line 1 and Line2 with same information. Comparing :".$fields1[$field4comp]." and ".$fields2[$field4comp]."\n";
			##..then match the $field4comp and if the first is greater..
			#if ($fields1[$field4comp] > $fields2[$field4comp]){
				#my $gene = $fields1[4];
				#my $coverage = $fields1[$field4comp];
				#print NEWFILE $compid1."_$transcript1\_$exon1\_$stend1\_$gene\t$coverage\n";
				##otherwise print the second line
			#}else{
				#my $gene = $fields2[4];
				#my $coverage = $fields2[$field4comp];
				#print NEWFILE $compid2."_$transcript2\_$exon2\_$stend2\_$gene\t$coverage\n";
			#}
			##Go ahed with next line
			#$i++;
		#}else{
				#my $gene = $fields1[4];
				#my $coverage = $fields1[$field4comp];
				#print NEWFILE $compid1."_$transcript1\_$exon1\_$stend1\_$gene\t$coverage\n";
		#}
	}
	close(FILE);
	close(NEWFILE);
}

# Given the output of query_2_db-VARIANTS. It puts on the same line
#the samples if the analysis name and the variant are the same
sub who_has_variant_row_2_analysis {

	my $file = shift;
	my $outname = shift;
	
	open (FILE,"<$file") or die "Cannot open $file\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $var1 = $fields1[0];
		my $anal1 = $fields1[1];
		my $var2 = $fields2[0];
		my $anal2 = $fields2[1];
		
		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $var1."\t$anal1\t".$fields1[2]."\t".$fields1[3]."\t";
		while( $var1 eq $var2  and $anal1 eq $anal2){
			#print " $var1 - $anal1 and $var2 - $anal2  is same analysis.\n";
			$row_2_write .= $fields2[2]."\t".$fields2[3]."\t";
			#Go ahed with next line
			$i++;				
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$var2 = $fields2[0];
			$anal2 = $fields2[1];			
			
		}
		chop($row_2_write);
		print NEWFILE $row_2_write."\n";
		
	}
	close(FILE);
	close(NEWFILE);
}



# Given a file with chromosome and position of a variant
# If there are two close variants depending on maxdist keeps only the first one
sub remove_close_variants {

	my $file = shift;
	my $outname = shift;
	my $maxdist = shift;
	
	open (FILE,"<$file") or die "Cannot open $file\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $chr1 = $fields1[0];
		my $pos1 = $fields1[1];
		my $chr2 = $fields2[0];
		my $pos2 = $fields2[1];
		
		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $lines[$i]."\n";
		while( $chr1 eq $chr2  and $pos2-$pos1 <=  $maxdist){
			#print "$compid1 and $compid2 is same information.\n";
			
			#Go ahed with next line
			$i++;				
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$chr2 = $fields2[0];
			$pos2 = $fields2[1];			
			
		}
		print NEWFILE $row_2_write
		
	}
	close(FILE);
	close(NEWFILE);
}

# Given the output from XHMM in zscores for N samples 
#It just extract each row to generate a bedgraph with this structure:

#track type=bedGraph name=sample-name description=center_label visibility=full graphType=points color=200,100,0 altColor=0,100,200
#chr1    65510   65625   0.19308252
#chr1    65832   65973   -1.35452024
# ....
sub XHMM_zscores_2_bedgraph {
	my $input_f = shift;
	
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	
	my $numline = 0;
	#Set an array with all the intervals found on first line
	my @intervals = ();
	#Read all the zscores file
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		#Get the intervals from the header into the array
		if ($numline == 0){
			for (my $i = 1; $i < scalar(@fields) ; $i++){
					push(@intervals,$fields[$i]);
			}
		}
		else{
			my $samplename = $fields[0];
			#Generate a new bedgraph file
			open(NEWFILE,">$samplename.bedgraph");
			print "Printing file $samplename.bedgraph\n";
			#Print the header
			print NEWFILE "track type=bedGraph name=$samplename description=center_label visibility=full graphType=points color=200,100,0 altColor=0,100,200\n";
			if ( scalar(@intervals) != scalar(@fields)-1){die "ERROR: number of intervals differs from number of scores at line $numline.\n";}
			my $int_num = 1;
			#For each interval read in the header print the information separated by tab
			foreach my $interval (@intervals){
				my @pieces = split(":",$interval);
				my $chr = $pieces[0];
				my @ints = split("-",$pieces[1]);
				print NEWFILE $chr."\t".$ints[0]."\t".$ints[1]."\t".$fields[$int_num]."\n";
				$int_num++;
			}
			close(NEWFILE);
		}
		$numline++;
	}
		
	close(FILE);
}


# Given the output from ExomeDepth and XHMM 
#It extracts a BED structure
# From  chr5:175391955-175477735_DUP
#chr5    175391955   175477735 chr5:175391955-175477735_DUP
# ....
sub compid_intervals_2_bed {
	my $input_f = shift;
	
	my $outname = extract_name($input_f,"noext").".bed";
	
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	#Generate a new bed file
	open(NEWFILE,">$outname");
				
	#Read all the file
	foreach my $line (<FILE>){
		chop($line);
		#First separate chr
		my @fields1 = split(":",$line);
		my $chr = $fields1[0];
		my @fields2 = split("-",$fields1[1]);
		my $start = $fields2[0];
		my @fields3 = split("_",$fields2[1]);
		my $end = $fields3[0];	
		
		print NEWFILE "$chr\t$start\t$end\t$line\n";
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
	
	my $outname = extract_name($input_f,"noext").".bed";
	
	open (FILE,"<$input_f") or die "Cannot open $input_f\n";
	
	
	my $numline = 0;
	my $idfield = -1;
	
	#Generate a new bed file
	open(NEWFILE,">$outname");
				
	#Read all the file
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

# This script is useful when we have a BAM file containing sequences which contain a subject sequence and 
#for the rest they align on the genome. The CIGAR string helps in finding which sequences match only partially.
#You should give in input a BAM file on which you already obtained sequences with a CIGAR containing S (that indicates the
#lenght of the query that do not map).
#Using the CIGAR string, this script will cut the sequence and generate a fasta file with all suspect sequences
#that can be used for downstream analyses
sub cut_SAM_sequences_using_CIGAR {
	my $sam = shift;
	my $fasta = shift;
		
	#Minuimun seq length to print in fasta
	my $min_seq_len = 5;	
	open (FASTA,">$fasta") or die "Cannot open $fasta\n";
	open (SAM,"<$bam") or die "Cannot open $bam\n";
	my $num_seq_l = 0;
	my $num_seq_r = 0;
	my $tot_seq = 0;
	foreach my $line (<SAM>){
			chop($line);
			print "Line $tot_seq: $line\n";
			my @fields = split("\t",$line);
			my $rname = $fields[0];
			my $cigar = $fields[5];
			my $seq = $fields[9];

			my $left_l = 0;
			my $right_l = 0;
			my $centr_l = 0;
			if ($cigar =~ /(\d+)S(\d+)M(\d+)S/){
					$left_l = $1;
					$centr_l = $2;
					$right_l = $3;
			}elsif ( $cigar =~ /(\d+)M(\d+)S/ ){
					$right_l = $2;
					}elsif ($cigar =~ /(\d+)S(\d+)M/){
							$left_l = $1;
							}else{
									print "ERROR: Cannot match any combination: SMS,SM,MS. Verify your SAM file..\n";
									}
			if ( $centr_l > 1){$suff="S";};#If there is a center and a right and left, this match is doubtful
			if ( $left_l > $min_seq_len){print FASTA ">seq_$tot_seq\_$left_l\_L-$rname\n".substr($seq,0,$left_l)."\n";}
			if ( $right_l > $min_seq_len){print FASTA ">seq_$tot_seq\_$right_l\_R-$rname\n".substr($seq,-($right_l))."\n";}

			$tot_seq++;
	}
	close(SAM);
	close(FASTA);
}



# Given a bed file with the first three lines giving the chromosome, the start
# and the end. If there are two identical consecutive intervals keeps only the first one
sub get_extended_regions {
	my $targetbed = shift;
	my $outbed = shift;
	my $dim = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines and printing extended bed into $outbed\n";
	
	my $prevchr = "";
	my $prevend = -1;
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		if ( $lines[$i] =~ /^chr/){
			chomp($lines[$i]);
			chomp($lines[$i+1]);
			my @fields1 = split("\t",$lines[$i]);
			my @fields2 = split("\t",$lines[$i+1]);
			
			my $thischr = $fields1[0];
			my $thisstart = $fields1[1];
			my $thisend = $fields1[2];
			my $nextchr = $fields2[0];
			my $nextstart = $fields2[1];
			
			
			#newstart is $dim nucleotides before and newend $dim after 
			my $newstart = $thisstart - $dim;
			my $newend = $thisend + $dim;
			
			#but if... now newstart is negative:
			if ( $newstart < 0 ){
					$newstart = 0;
			}
			#and if the newstart overlaps the previous end:
			if ( $newstart <= $prevend  and $thischr eq $prevchr){
					$newstart = $prevend + 1;
			}
			#or if the newend overlaps the next start:
			if ( $newend >= $nextstart  and $thischr eq $nextchr){
					$newend = $nextstart - 1;
			}
			#Set the new previous details
			$prevchr = $thischr;
			$prevend = $newend;
			print NEWFILE "$thischr\t$newstart\t$newend\n";	
		}
		
	}
}

# Given a bed file with the first three lines giving the chromosome, the start
# and the end. If there are two identical consecutive intervals keeps only the first one
sub remove_duplicates_intervals {

	my $bed = shift;
	my $outname = shift;
	my $field4comp = shift;
	
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $compid1 = $fields1[0]."_".$fields1[1]."_".$fields1[2];
		my $compid2 = $fields2[0]."_".$fields2[1]."_".$fields2[2];

		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $lines[$i]."\n";
		while( $compid1 eq $compid2 ){
			#print "$compid1 and $compid2 is same information.\n";
			
			#Go ahed with next line
			$i++;				
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$compid2 = $fields2[0]."_".$fields2[1]."_".$fields2[2];
		}
		print NEWFILE $row_2_write
		
	}
	close(FILE);
	close(NEWFILE);
}

# This function simply removes duplicate rows if the element
# indicated by $field4comp is identical. Keeps only the first that is met
sub remove_duplicates_rows {

	my $input = shift;
	my $outname = shift;
	my $field4comp = shift;
	
	open (FILE,"<$input") or die "Cannot open $input\n";
	open(NEWFILE,">$outname");
	#Read all file into an array
	my @lines = <FILE>;
	my $totlines = scalar(@lines);
	print "Analyzing $totlines lines\n";
	#Header
	#print NEWFILE "compid\tcov\n";
	
	#Go through all the array
	for (my $i=0; $i < ($totlines-1); $i++){
		#print "Line $i: ".$lines[$i]."\n";
		#print "Line ".($i+1).": ".$lines[$i+1]."\n";
		
		chop($lines[$i]);
		chop($lines[$i+1]);
		my @fields1 = split("\t",$lines[$i]);
		my @fields2 = split("\t",$lines[$i+1]);
		
		my $compid1 = $fields1[$field4comp];
		my $compid2 = $fields2[$field4comp];

		#While the current and next line have the same compid, go on with the index of row
		my $row_2_write = $lines[$i]."\n";
		while( $compid1 eq $compid2 ){
			#print "$compid1 and $compid2 is same information.\n";
			
			#Go ahed with next line
			$i++;				
			chop($lines[$i+1]);
			@fields2 = split("\t",$lines[$i+1]);
			$compid2 = $fields2[$field4comp];
		}
		print NEWFILE $row_2_write
		
	}
	close(FILE);
	close(NEWFILE);
}


##########utilities


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
