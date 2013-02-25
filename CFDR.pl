#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::Tools::GFF;

# r.farrer09@imperial.ac.uk

### Opening commands 

my $usage = "
Program:  Comparison of False Discovery Rates (CFDR)
Version:  0.11
Contact:  Rhys Farrer <r.farrer09\@imperial.ac.uk>

Usage:    $0 <commands>

Commands: -i\tFull list of introduced mutations
          -m\tSNP-calls in Variant Call Format (VCF)
          -p\tPileup of alignment
          
Optional: -c\tGFF/GTF (Takes longer)
          -f\tFeature in the GFF/GTF to CFDR (CDS, exon, mRNA etc.)
          -o\tOutput folder location. default: folder of found mutations
          -s\tSAM file of alignment (needed if -f and for correct true negatives)
";
our($opt_i, $opt_m, $opt_s, $opt_p, $opt_c, $opt_f, $opt_o);
getopt('imspcfo');
die $usage unless ($opt_i && $opt_m && $opt_p);  
die "Cannot open $opt_i : $!" unless (-e $opt_i);
die "Cannot open $opt_m : $!" unless (-e $opt_m);
die "Cannot open $opt_p : $!" unless (-e $opt_p);
if($opt_c) { die "Cannot open $opt_c : $!" unless (-e $opt_c); }
die "Specify both -c and -f or neither\n" if ((($opt_c) && (!defined $opt_f)) || (($opt_f) && (!defined $opt_c)));
if(defined $opt_o) { 
	if($opt_o !~ m/\/$/) { $opt_o .= '/'; }
	die "Output folder does not exist: $opt_o\n" unless (-d $opt_o); 
}
my ($TP_SNP, $FP_SNP, $FP_DEL)  = (0, 0, 0);
my ($TP_HET, $FP_HET) = (0,0);
my ($FN_no_cov, $FN_ref_base, $total_mutations_found) = (0, 0, 0);
my (%feature_positions, %contig_names);
my ($reference_seq, $number_of_introduced_mutations, $type_of_introduced_mutation) = make_ref($opt_i);

my $TN2_value = 10000000;
if($opt_c) {
	
	### Find the length of the reference sequence, and the names of valid supercontigs for feature file

	warn "Reading headers from $opt_s ...\n";

	open IN2, "<$opt_s" or die "Cannot open $opt_s: $!\n";
	READHEADERS: while (my $line = <IN2>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($header, $chr, $length) = @bits;
		if($header eq '@SQ') {
			$chr =~ s/^SN://;
			$length =~ s/^LN://;
			$contig_names{$chr} = $length;
			$TN2_value += $length;
		}
		last READHEADERS if($header !~ m/^@/); 
	}
	close IN2;
	
	### Go through GFF if specified to check for features
	
	warn "Saving locations of $opt_f in $opt_c...\n";
	
	### Check if it is a gff version 2
	
	my $feature_count = 0;
	my $gff_version = 2;
	$gff_version = &check_version_of_gff($opt_c, $gff_version, $opt_f);
	
	### Check if it is a gff version 3
	if($gff_version eq 3) {
		$gff_version = &check_version_of_gff($opt_c, $gff_version, $opt_f);
		if($gff_version eq 4) { die "No feature name for the first $opt_f in $opt_c\n"; }
	}
	
	my $gffio = Bio::Tools::GFF->new('-file' => "<$opt_c", -gff_version => $gff_version);
	while(my $feat = $gffio->next_feature()) {
		my $str = $feat->gff_string;
		my @bits = split /\t/, $str;
		my ($contig, $source, $feature, $from, $to, $score, $strand, $frame, $name) = @bits;
		
		if($opt_f eq $feature) {
			# Check the feature are in range of the contig and have a corresponding value in the genome
			die "$contig in $opt_c not found in $opt_s" if (!defined $contig_names{$contig});
			die "Feature $feature $contig $name out of range (<1) in genome" if(($from < 1) || ($to < 1));
			die "Feature $feature $contig $name out of range in genome" if(($from > $contig_names{$contig}) || ($to > $contig_names{$contig}));
			
			$feature_count++;
			my $temp_feature_positions = ($from . '-' . $to . " ");
			$feature_positions{$contig}{$name} .= $temp_feature_positions;
		}
        }
        $gffio->close();
	warn "Have found $feature_count $opt_f 's...\n";
}
%contig_names = ();

### Go through the file of polymorphic sites from the SNP-caller and find TP/FP rate

my (@FP_lines, @FN_lines, @FN_lines_no_cov);
my $location_of_genotype_in_format = 0;

open IN1, "<$opt_m";
warn "Cataloging SNPs in from $opt_m...\n";
FOUNDPOLY: while (my $line = <IN1>) {
   	chomp $line;
   	
   	### Read details from VCF
   	
   	next FOUNDPOLY if ($line =~ m/^\#/);
   	my @bits = split /\t/, $line;
   	next FOUNDPOLY if (@bits < 9);
   	die "$0 Currently not made to check more than a VCF for single isolate\n" if (@bits > 10);
   	my ($supercontig, $position, $id, $reference_VCF_format, $consensus_VCF_format, $cons_qual, $filter, $info, $format, $sample_info) = @bits;
	my @sample_info_parts = split /:/, $sample_info;
	my $alleles = $sample_info_parts[$location_of_genotype_in_format];
	
	# Check alleles are different
	if($alleles =~ m/\//) {
		my @allele_parts = split /\//, $alleles;
		if($allele_parts[0] eq $allele_parts[1]) { $alleles = $allele_parts[0]; }
	}

	### if the feature file is being used, check if position is over the wanted feature
	
	my $over_feature = 0;
	if($opt_c) {
		foreach my $feature_names(sort keys %{$feature_positions{$supercontig}}) {
			my $feature_positions = $feature_positions{$supercontig}{$feature_names};
			my @seperate_positions = split /\s/, $feature_positions;
			foreach(@seperate_positions) {
				my @terminals = split /-/, $_;
				if(($position >= $terminals[0]) && ($position <= $terminals[1])) { $over_feature = 1; }
			}
		}
		next FOUNDPOLY if ($over_feature eq 0);  
	}
	
	$$reference_seq{'seen'}{$supercontig}{$position} = 1;
	
	### If found mutation is over a saved position
	
	if(exists $$reference_seq{$supercontig}{$position}) {
		
		### Get saved positions
		
		my $corresponding_ref = $$reference_seq{$supercontig}{$position}{'ref_base'};
		my $corresponding_con = $$reference_seq{$supercontig}{$position}{'cons_base'};
		my $corresponding_con2;
		if(defined $$reference_seq{$supercontig}{$position}{'cons_base2'}) { $corresponding_con2 = $$reference_seq{$supercontig}{$position}{'cons_base2'}; }
		
		# Deletion
		if($info =~ m/INDEL/) { $FP_DEL++; }
		
		# Heterozygous position (bi-allelic) 
		elsif($alleles =~ m/(\d)\/(\d)/) { 
			
			if($type_of_introduced_mutation eq 'SNP') { $FP_HET++; }
			elsif($type_of_introduced_mutation eq 'HET') { 
				
				# Hets. neither reference
				if(defined $corresponding_con2) {
					
					if($consensus_VCF_format =~ m/,/) {
						my @consensi = split /,/, $consensus_VCF_format;
						
						# 1-way or the other way round
						if(($corresponding_con eq $consensi[0]) && ($corresponding_con2 eq $consensi[1])) { $TP_HET++; }
						elsif(($corresponding_con eq $consensi[1]) && ($corresponding_con2 eq $consensi[0])) { $TP_HET++; }
						else { $FP_HET++; }
					}
					
					# If it doesn;t have both, consider it a FP
					else { $FP_HET++; }
					
				} 
				# Het with reference
				else {
					
					if(($corresponding_ref eq $reference_VCF_format) && ($corresponding_con eq $consensus_VCF_format)) { $TP_HET++; }
					else {
						$FP_HET++;
						#warn "should be 0/1 with ref: $corresponding_ref $corresponding_con : $line\n";
					}
				}
			}
			
			#elsif($type_of_introduced_mutation eq 'HETTRIP') { $TP_HET++; }
		}
		# Homozygous base - Not Reference base (FP)
		elsif(($consensus_VCF_format ne $corresponding_ref) && ($alleles !~ m/\//)) { $FP_SNP++; push @FP_lines, $line; }
		# Homozygous SNP
		elsif(($consensus_VCF_format eq $corresponding_ref) && ($reference_VCF_format eq $corresponding_con)) { 
			if($type_of_introduced_mutation eq 'SNP')     { $TP_SNP++; }
			elsif($type_of_introduced_mutation eq 'HET') { $FP_SNP++; }
		}
		# Unrecognised
		else { warn "Unrecognised: $line found at correct location: ref base: $reference_VCF_format con: $consensus_VCF_format corr ref: $corresponding_ref corr con $corresponding_con\n"; }
	}
	
	### Not a saved position
	
	else {
		# Deletion
		if($info =~ m/INDEL/) { $FP_DEL++; }
		# Heterozygous base (FP)
		elsif($alleles =~ m/\//) { $FP_HET++; }
		# Homozygous base - Not Reference base (FP)
		elsif($alleles !~ m/\//) { $FP_SNP++; push @FP_lines, $line; }
		# Unrecognised
		else { warn "Unrecognised: $line ref base: $reference_VCF_format con: $consensus_VCF_format\n"; }
	}
}
close IN1;

### True Negatives

warn "Calculating true negatives (specified) from $opt_p ...\n";
my $genome_length_covered = `wc -l $opt_p`;
$genome_length_covered =~ s/^\s+//;
my @first_num = split /\s/, $genome_length_covered;
$first_num[0] =~ s/\s//g;
$genome_length_covered = $first_num[0];

my $TN_value = ($genome_length_covered - ($TP_SNP + $FP_SNP + $FP_HET + $FP_DEL));
$TN2_value -= $genome_length_covered;

### False negatives (covered)

open IN2, "<$opt_p";
warn "Identifying false negatives from $opt_p...\n";
my $num_positions_covered = 0;
PILEUP: while (my $line = <IN2>) {
   	chomp $line;
   	my @bits = split /\t/, $line;
   	next PILEUP if (@bits != 6);
   	
	### If it is a saved position of a mutation (introduced)
	if(exists $$reference_seq{$bits[0]}{$bits[1]}) { 
		$num_positions_covered++; 
		
		# But not a position found
		if(!exists $$reference_seq{'seen'}{$bits[0]}{$bits[1]}) { 
			$FN_ref_base++; push @FN_lines, $line; 
			$$reference_seq{'seen'}{$bits[0]}{$bits[1]} = 0;
		}
	}
}
close IN2;
$FN_no_cov = ($number_of_introduced_mutations - $num_positions_covered);

### False negatives (no coverage)

FNNOCOV: foreach my $contig(sort keys %$reference_seq) {
	next FNNOCOV if ($contig =~ m/^seen$/);
	foreach my $positions(sort {$a<=>$b} keys %{$$reference_seq{$contig}}) {
		if(!exists $$reference_seq{'seen'}{$contig}{$positions}) {
			push @FN_lines_no_cov, "$contig\t$positions";
		}
	}
}

### Calculate results

my ($FP_percent, $FN_percent, $TN2_percent, $FN2_percent, $TN_percent) = (0, 0, 0, 0, 0);
my $TP_percent = (($TP_SNP / $number_of_introduced_mutations) * 100);
$TP_percent = sprintf("%.2f", $TP_percent);
if(($TP_SNP ne 0) && ($FP_SNP ne 0)) { $FP_percent = (($FP_SNP / ($TP_SNP + $FP_SNP)) * 100); }
$FP_percent = sprintf("%.2f", $FP_percent);
if($TP_SNP ne 0) { $TN_percent = (($TN_value / ($genome_length_covered + $TP_SNP)) * 100); }
$TN_percent = sprintf("%.2f", $TN_percent);
if(($FN_ref_base ne 0) && ($FN_no_cov ne 0)) { 
	$FN_percent = ((($FN_ref_base + $FN_no_cov) / $number_of_introduced_mutations) * 100);
}
$FN_percent = sprintf("%.2f", $FN_percent);
if(defined $genome_length_covered) {
	$TN2_percent = (($TN2_value / $genome_length_covered) * 100);
}
$TN2_percent = sprintf("%.2f", $TN2_percent);
$FN2_percent = (($FN_no_cov / $number_of_introduced_mutations) * 100);
$FN2_percent = sprintf("%.2f", $FN2_percent);

### Add spaces to the beggining of percents and dashes for printing

$TP_percent  = (&add_chars((6 - length($TP_percent)), ' ')  . $TP_percent);
$FP_percent  = (&add_chars((6 - length($FP_percent)), ' ')  . $FP_percent);
$TN_percent  = (&add_chars((6 - length($TN_percent)), ' ')  . $TN_percent);
$FN_percent  = (&add_chars((6 - length($FN_percent)), ' ')  . $FN_percent);
$TN2_percent = (&add_chars((6 - length($TN2_percent)), ' ') . $TN2_percent);
$FN2_percent = (&add_chars((6 - length($FN2_percent)), ' ') . $FN2_percent);
my $TP_spaces = '';
if($TP_SNP > 10) { $TP_spaces = &add_chars((length($TP_SNP) - 1), ' '); }
my $short_line = '---------------';
my $line;
if($TP_SNP ne 0) { $line = &add_chars((length($TP_SNP) - 1), '-'); } 
$line .= '-----------------------------------------------------------';
my ($add_to_FP_line, $add_to_TP_line) = ('','');
if(length($TP_SNP) > length($FP_SNP)) { $add_to_FP_line = &add_chars((length($TP_SNP) - length($FP_SNP)), ' '); }
if(length($FP_SNP) > length($TP_SNP)) { $add_to_TP_line = &add_chars((length($FP_SNP) - length($TP_SNP)), ' '); }
my $add_to_TN_line = &add_chars((10 - length($TN_value)), ' ');
my $add_to_FN_line = &add_chars((10 - length($FN_ref_base)), ' ');
my $add_to_TN2_line = &add_chars((12 -  length($TN2_value)), ' ');
my $add_to_FN2_line = &add_chars((12 -  length($FN_no_cov)), ' ');

### Print results to file

my $OUTFIL = $opt_m;
if(defined $opt_o) {
	my @location_parts = split /\//, $OUTFIL;
	my $OUTFIL_name = $location_parts[(scalar(@location_parts) - 1)];
	$OUTFIL = ($opt_o . $OUTFIL_name);
}
my $outfile = $OUTFIL . '-CFDR-Summary';
if($opt_f) { $outfile .= ('-over_' . $opt_f); }
$outfile .= '.tab';
open OUT, ">$outfile";

### Print results (for homozygous)

if($type_of_introduced_mutation eq 'SNP') { 

	# outfile
	
	print OUT "True Positive Homozygous SNPs identified:\t$TP_SNP\t($TP_percent\%)\n";
	print OUT "False Positive Homozygous SNPs identified:\t$FP_SNP\t($FP_percent\%)\n";
	print OUT "True Negative Homozygous reference:\t$TN_value\t($TN_percent\%)\n";
	print OUT "False Negative Homozygous reference:\t$FN_ref_base\t($FN_percent\%)\n";
	print OUT "True Negative Homozygous (no coverage):\t$TN2_value\t($TN2_percent\%)\n";
	print OUT "False Negative Homozygous (no coverage):\t$FN_no_cov\t($FN2_percent\%)\n";
	print OUT "Heterozygous SNPs:\t$FP_HET\n";
	print OUT "Homo/hetero Deletions:\t$FP_DEL\n$short_line\n";
	
	print OUT "False Positive Homozygous SNPs\t($FP_SNP)\n$short_line\n";
	foreach (@FP_lines) { print OUT ($_ . "\n"); }
	print OUT "\n$short_line\nFalse Negatives reference ($FN_ref_base)\n$short_line\n";
	foreach (@FN_lines) { print OUT ($_ . "\n"); }
	print OUT "\n$short_line\nFalse Negatives no coverage ($FN_no_cov)\n$short_line\n";
	foreach (@FN_lines_no_cov) { print OUT ($_ . "\n"); }

	# terminal

	warn "\n\t\t\t\t\t\tHomozygous\n";
warn "         ------" . $line . "
         | Positive " . $TP_spaces . $add_to_TP_line . "   | Negative (reference) | Negative (no coverage)   |\n" .
$short_line . $line . "
| True   | " . $add_to_TP_line . "$TP_SNP ($TP_percent\%) | " . $add_to_TN_line . 
	"$TN_value ($TN_percent\%)   | " . $add_to_TN2_line . "$TN2_value ($TN2_percent\%) |\n" .
$short_line . $line . "
| False  | " . $add_to_FP_line . "$FP_SNP ($FP_percent\%) | " . $add_to_FN_line . 
	"$FN_ref_base ($FN_percent\%)   | " . $add_to_FN2_line . "$FN_no_cov ($FN2_percent\%) |\n" .
$short_line . $line . "\n\t\t\t\t\t\tOther:\n\n";

	warn "Heterozygous SNPs:     $FP_HET\n";
	warn "Homo/hetero Deletions: $FP_DEL\n\n";

}

### Print results (for heterozygous)

if($type_of_introduced_mutation eq 'HET') { 

	$TP_percent = (($TP_HET / $number_of_introduced_mutations) * 100);
	$TP_percent = sprintf("%.2f", $TP_percent);
	if(($TP_HET ne 0) && ($FP_HET ne 0)) { $FP_percent = (($FP_HET / ($TP_HET + $FP_HET)) * 100); }
	$FP_percent = sprintf("%.2f", $FP_percent);
	
	# outfile
	
	print OUT "True Positive Heterozygous SNPs identified:\t$TP_HET\t($TP_percent\%)\n";
	print OUT "False Positive Heterozygous SNPs identified:\t$FP_HET\t($FP_percent\%)\n";	
	
	# terminal
	
	warn "True Positive Heterozygous SNPs identified:\t$TP_HET\t($TP_percent\%)\n";
	warn "False Positive Heterozygous SNPs identified:\t$FP_HET\t($FP_percent\%)\n";	
}
close OUT;

sub make_ref {
	my $filename = shift;
	my (%found);
	my $count = 0;
	my $type;
	open IN, "<$filename";
	warn "Cataloging SNPs in from $filename...\n";
	DETAILS: while(my $line = <IN>) {
		chomp $line;
		if($line =~ m/^Contig\tPosition/) {
			if($line =~ m/SNP\t/) { $type = 'SNP'; }
			elsif($line =~ m/HET\t/) { $type = 'HET'; }
			elsif($line =~ m/HETTRIP\t/) { $type = 'HETTRIP'; }
			next DETAILS;
		}
		my @bits = split /\t/, $line;
		# The inserted mutations have their array position so need to +1
		my ($contig, $position, $nt, $introduced_mutation) = (@bits);
		die "Check if this is the details file: $filename\n" if ((!defined $contig) || (!defined $position) || (!defined $nt) || (!defined $introduced_mutation));
		die "Check if this is the details file (for mutation type): $filename\n" if (!defined $type);
		
		### If HET or HETTRIP are over a 'homologous sequence', 
		### should have been aligned to non-modified reference.
		if(($type eq 'HET') || ($type eq 'HETTRIP')) {
			if($contig =~ m/\_homologous/) {
				my @contig_parts = split /\_homologous/, $contig;
				$contig = $contig_parts[0];
				#$found_twice
			}
		}
		
		# Genotype 1/2 for HET
		if((defined $found{$contig}{($position + 1)}{'ref_base'}) || (defined $found{$contig}{($position + 1)}{'cons_base'})) {
			my $prev_found_ref_base = $found{$contig}{($position + 1)}{'ref_base'};
			#die "Ref base of homologous chromosomes don't match. Re-run analysis or report error on $contig, $position + 1 as $nt not equal $prev_found_ref_base\n" if($nt ne $found{$contig}{($position + 1)}{'ref_base'});
			$found{$contig}{($position + 1)}{'cons_base2'} = $introduced_mutation;
		}
		else {
			$found{$contig}{($position + 1)}{'ref_base'} = $nt;
			$found{$contig}{($position + 1)}{'cons_base'} = $introduced_mutation;
		}
		$count++;
	}
	close IN;
	die "Unrecognised details file. Currently needs to have either introduced SNP or HET" if(!defined $type);
	return (\%found, $count, $type);
}

sub add_chars {
	my ($num, $type) = @_;	
	my $string_spaces;
	for(my $i=0; $i<$num; $i++) { $string_spaces .= $type; }
	return $string_spaces;
}

sub check_version_of_gff {
	my ($gff_file, $gff_version, $feature_looking_for) = @_;
	my $gffio = Bio::Tools::GFF->new('-file' => "<$gff_file", -gff_version => $gff_version);
	TESTGFF: while(my $feat = $gffio->next_feature()) {
		my $str = $feat->gff_string;
		my @bits = split /\t/, $str;
		my ($feature, $name) = ($bits[2], $bits[8]);
		
		if($feature eq $feature_looking_for) {
			if(!defined $name) { $gff_version++; }
			last TESTGFF;
		}
	}
        $gffio->close();
        return ($gff_version);
}
