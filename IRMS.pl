#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use Bio::Tools::GFF;

# r.farrer09@imperial.ac.uk

### Opening commands 

my $usage = "
Program:  IRMS (Introduce Random Mutations into Sequence)
Version:  0.11
Contact:  Rhys Farrer <r.farrer09\@imperial.ac.uk>

Usage:    Introduce-Random-Mutations-into-Sequence.pl <commands>

Commands: -g\tReference sequence in FASTA format
          -n\tNumber of mutations to introduce
          
Optional: -t\tType of mutation to introduce (SNP/DEL/INS/HET/CNV/HETRIP) note. INS/DEL not supported for -c/-f (default SNP). CNV needs -c and -f
          -c\tGFF/GTF (Takes longer)
          -f\tFeature in the GFF/GTF to mutate (CDS, exon, mRNA etc.)
";

our($opt_g, $opt_c, $opt_n, $opt_t, $opt_f);
getopt('gcntf');
if(!defined $opt_t) { $opt_t = 'SNP'; }
if(!defined $opt_f) { $opt_f = 'Genome'; }
die $usage unless ($opt_g && $opt_n);     
die "Specify both -c and -f or neither -c specified ($opt_c)\n" if ((defined $opt_c) && (!defined $opt_f));
die "Specify both -c and -f or neither -f specified ($opt_f)\n" if (($opt_f ne 'Genome') && (!defined $opt_c));
die "Specify both -c and -f if -t is CNV\n" if (($opt_t eq 'CNV') && ((!defined $opt_c) && ($opt_f eq 'Genome')));
die "-t has to be SNP, DEL, INS or CNV\n" unless (($opt_t eq 'SNP') || ($opt_t eq 'DEL') || ($opt_t eq 'INS') || ($opt_t eq 'CNV') || ($opt_t eq 'HET') || ($opt_t eq 'HETTRIP'));
die "Cannot open $opt_g\n" unless (-e $opt_g);
die "Number of mutations to introduce must be greater than 1\n" if ($opt_n < 1);
die "Currently, -t has to be 'SNP' or 'HET' if specifying -c/-f\n" if ((($opt_t eq 'DEL') || ($opt_t eq 'INS')) && ($opt_c));
my (%feature_positions, %rand_numbers, %genome_fna, %genome_size, %feature_strand, %position_names_to_gene_names, %feature_seqs);
my ($range, $number_of_Ns_over_rand_numbers) = (0, 0);

### Specify output files

my $genome_with_identifies = ($opt_g . "-$opt_n-" . $opt_t . "-in-" . $opt_f);
my $mutated_genome_details = ($genome_with_identifies . ".details");
my $mutated_genome = ($genome_with_identifies . ".fasta");

### Save genome/sequence for mutating (make homologous version if heterozygous)

warn "Saving $opt_g...\n";
my $inseq = Bio::SeqIO->new('-file' => "<$opt_g",'-format' => 'fasta');
while (my $seq_obj = $inseq->next_seq ) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	$genome_fna{$id} = $seq;
	$genome_size{$id} = length($seq);
	if(!defined $opt_c) { $range += length($seq); }
	
	### For hets
	
	if (($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) {
		$genome_fna{($id . '_homologous')} = $seq;
		$genome_size{($id . '_homologous')} = length($seq);
		if(!defined $opt_c) { $range += length($seq); }
	}
	
	if ($opt_t eq 'HETTRIP') {
		$genome_fna{($id . '_homologous2')} = $seq;
		$genome_size{($id . '_homologous2')} = length($seq);
		if(!defined $opt_c) { $range += length($seq); }
	}
}

### Go through GFF if specified to check for features. 
### Ignore features that overlap (such as splice variants) 

if($opt_c) {
	warn "Checking $opt_c for $opt_f overlaps...\n";
	
	### Check if it is a gff version 2
	my $gff_version = 2;
	$gff_version = &check_version_of_gff($opt_c, $gff_version, $opt_f);
	
	### Check if it is a gff version 3
	if($gff_version eq 3) {
		$gff_version = &check_version_of_gff($opt_c, $gff_version, $opt_f);
		if($gff_version eq 4) { die "No feature name for the first $opt_f in $opt_c\n"; }
	}
	
	my $gffio = Bio::Tools::GFF->new('-file' => "<$opt_c", -gff_version => $gff_version);
	my @temporary_positions_for_checking_overlap;
	while(my $feat = $gffio->next_feature()) {
		my $str = $feat->gff_string;
		my @bits = split /\t/, $str;
		my ($contig, $source, $feature, $from, $to, $score, $strand, $frame, $name) = @bits;
		
		if($opt_f eq $feature) {
			
			### Check the feature are in range of the contig and have a corresponding value in the genome
			
			die "$contig in $opt_c not found in $opt_g" if (!defined $genome_fna{$contig});
			die "Feature $feature $contig $name out of range in genome" if(($from < 1) || ($to < 1) || ($from > $genome_size{$contig}) || ($to > $genome_size{$contig}));
			die "No name for feature $feature on $contig $from - $to\n" if(!defined $name);
			
			### Miss out any overlapping exons/splice variants
			
			my $ignore = 0;
			OVERLAPS: foreach my $positions(@temporary_positions_for_checking_overlap) {
				$positions =~ s/\s//;
				my @Start_Stop = split /-/, $positions;
				
				# Check if new feature overlaps the start (and all) 
				if(($from < $Start_Stop[0]) && ($to > $Start_Stop[0])) { $ignore = 1; last OVERLAPS; }
				# Check if new feature overlaps over the end (or middle)
				if(($from > $Start_Stop[0]) && ($from < $Start_Stop[1])) { $ignore = 1; last OVERLAPS; }
			}
			if($ignore eq 0) {
				if($opt_t eq 'CNV') { 
				
					if(!defined $feature_seqs{$contig}{$name}) {
						$range++; 
					}
					
					# Pull out feature seq if modifying CNV
					# +1 to turn array position to genome position
					my $feature_sequence = substr $genome_fna{$contig}, ($from - 1), ($to - $from + 1);
					#warn "my features sequence added to $name - $from - $to = $feature_sequence\n";
					$feature_seqs{$contig}{$name} .= $feature_sequence;
				}
				else { 
					$range += ($to - $from); 
					if ($opt_t eq 'HET') { $range += ($to - $from); }
					if ($opt_t eq 'HETTRIP') { $range += (($to - $from) * 2); }
				}
				
				#warn "$range from adding $contig $name\n";
				my $temp_feature_positions = ($from . '-' . $to . " ");
				push (@temporary_positions_for_checking_overlap, $temp_feature_positions);
				$feature_positions{$contig}{$name} .= $temp_feature_positions;
				
				if (($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) { 
					$feature_positions{($contig . '_homologous')}{$name} .= $temp_feature_positions;
				}
				if ($opt_t eq 'HETTRIP') { 
					$feature_positions{($contig . '_homologous2')}{$name} .= $temp_feature_positions;
				}
			}
		}
        }
        $gffio->close();
	
}

### Check range is bigger than number of mutations

if($opt_t eq 'CNV') { 
	die "Too many mutations required ($opt_n) for number of features ($opt_f) found for mutating ($range)\n" if($range < $opt_n);
} else {
	die "Too many mutations required ($opt_n) for amount of sequence (nt) found for mutating ($range)\n" if($range < $opt_n);
}

### Random positions generator

warn "Generating $opt_n / $range $opt_t...\n";
while(scalar(keys %rand_numbers) < $opt_n) {
	my $random_number = int(rand($range)) + 1; 
	$rand_numbers{$random_number} = 1;
}

### Modify anywhere in the genome

open OUT1, ">$mutated_genome_details";
print OUT1 "Contig\tPosition\tnt\tintroduced $opt_t\tFeature\n";
if(!defined $opt_c) {
	warn "Introducing mutations into $opt_g\n";
	my $seq_length_gone_through = 0;
	foreach my $id(sort keys %genome_fna) {
		my $seq = $genome_fna{$id};
		my @nts = split //, $seq;
		
		foreach my $numbers(sort {$a<=>$b} keys %rand_numbers) {
			my $random_number_in_id = ($numbers - $seq_length_gone_through);
			if(($random_number_in_id <= length($seq)) && ($random_number_in_id > 1)) {
				
				my $nt = $nts[$random_number_in_id];
				if($nt !~ m/[ACTGactg]/) { 
					
					$number_of_Ns_over_rand_numbers++; 
					my ($new_random_nucleotide, $random_number2);
					while(!defined $new_random_nucleotide) {
						$random_number2 = int(rand(scalar(@nts))) + 1; 
						if(!exists $rand_numbers{($seq_length_gone_through + $random_number2)}) {
							my $test_nt = $nts[$random_number2];
							if($test_nt =~ m/[ACTGactg]/) { $new_random_nucleotide = $test_nt; }
						}
					}
					$random_number_in_id = ($random_number2);
					$nt = $new_random_nucleotide;
				}
				my ($newbase);	
				if($opt_t eq 'DEL') { $newbase = ''; }
				if(($opt_t eq 'SNP') || ($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) { $newbase = randomnucleotide($nt); }
				if($opt_t eq 'INS') { $newbase = ($nt . randomnucleotide($nt)); }
				print OUT1 "$id\t$random_number_in_id\t$nt\t$newbase\t$opt_g\n";
				$nts[$random_number_in_id] = $newbase;
			}
		}
		$seq_length_gone_through += length($seq);
		$genome_fna{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t + $number_of_Ns_over_rand_numbers that were picked over Ns in first call\n";
}

### Modify specific places in the genome

if($opt_c) {
	warn "Introducing mutations into the $opt_f of $opt_g\n";
	my $random_numbers_gone_through = 0;
	# Go through each of the contigs
	foreach my $id(sort keys %genome_fna) {
		my $seq = $genome_fna{$id};
		my @nts = split //, $seq;
		# Go through each feature that did not have any overlaps with other features
		foreach my $feature_names(sort keys %{$feature_positions{$id}}) {
			
			if(($opt_t eq 'SNP') || ($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) { 
				my $feature_positions = $feature_positions{$id}{$feature_names};
				my @seperate_positions = split /\s/, $feature_positions;
				foreach(@seperate_positions) {
					my @terminals = split /-/, $_;
					for(my $i=$terminals[0]; $i<$terminals[1]; $i++) {
						$random_numbers_gone_through++;
						if(exists $rand_numbers{$random_numbers_gone_through}) {
							my $nt = $nts[$i];
							if($nt eq 'N') { $number_of_Ns_over_rand_numbers++; }
							else {
								my ($newbase);	
								$newbase = randomnucleotide($nt);
								print OUT1 "$id\t$i\t$nt\t$newbase\t$feature_names\n";
								$nts[$i] = $newbase;
							}
						}
					}
				}
			}
			if($opt_t eq 'CNV') { 
				$random_numbers_gone_through++;
				if(exists $rand_numbers{$random_numbers_gone_through}) {
					my @feature_nts = split //, $feature_seqs{$id}{$feature_names};
					
					#warn "my feature nts = @feature_nts \n";
					
					my $pre_CNV_contig_length = (scalar(@nts) + 1);
					push @nts, @feature_nts;
					my $post_CNV_contig_length = scalar(@nts);
					print OUT1 "$id\t$pre_CNV_contig_length-$post_CNV_contig_length\tCNV+\tCNV+\t$feature_names\n";
					
					#warn "pre = $pre_CNV_contig_length and post = $post_CNV_contig_length\n";
					
				}
			}
		}
		$genome_fna{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t ($number_of_Ns_over_rand_numbers Ns found but not changed)\n";
	warn "$random_numbers_gone_through Random numbers gone through in feature file\n";
}
close OUT1;

### Print modified genome

warn "Printing modified $opt_g...\n";
open OUT2, ">$mutated_genome";
foreach my $contig(sort keys %genome_fna) { 
	print OUT2 ">$contig\n"; 
	my $seq = $genome_fna{$contig};
	my @seq_bits = split //, $seq;
	my ($count, $contig_count) = (0, 0);
	foreach(@seq_bits) {
		$count++; $contig_count++;
		print OUT2 "$_";
		if(($count eq 60) && (defined $seq_bits[$contig_count + 1])) { print OUT2 "\n"; $count=0 }
	}
	print OUT2 "\n";
	$contig_count = 0;
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


sub randomnucleotide { 
	my $consbase = shift;
	my @nucs;
	if($consbase eq 'A') { @nucs = ('C','G','T'); }
	if($consbase eq 'C') { @nucs = ('A','G','T'); }
	if($consbase eq 'T') { @nucs = ('A','C','G'); }
	if($consbase eq 'G') { @nucs = ('A','C','T'); }
	return $nucs[rand @nucs];
}

