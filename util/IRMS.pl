#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_GFF;
use read_FASTA;
use File::Basename;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Program:  IRMS (Introduce Random Mutations into Sequence)
Usage: perl $0 -s <Reference FASTA> -n <Number of mutations to introduce>\n
Optional: -t\tType of mutation to introduce (SNP, DEL, INS, HET, CNV, HETRIP) [SNP]
          -g\tGFF/GTF [none]
          -f\tFeature to mutate (CDS, exon, mRNA, genome) [genome]
          -p\tSeperator in description for gene names (\" ; etc) [;]
          -c\tColumn that says which parent/gene it belongs to the description [0]
          -r\tRemove additional comments in column [Parent=]\n
Output:   -o\tOutput details [<opt_s>-n-<opt_n>-t-<opt_t>-f-<opt_f>.details]
          -a\tOutput FASTA [<opt_s>-n-<opt_n>-t-<opt_t>-f-<opt_f>.fasta]\n
Note: INS/DEL not supported for -c/-f (default SNP). CNV needs -c and -f\n";
our($opt_a, $opt_c, $opt_f, $opt_g, $opt_n, $opt_o, $opt_p, $opt_r, $opt_s, $opt_t);
getopt('acfgnoprst');
die $usage unless ($opt_s && $opt_n);     
die "Cannot open $opt_s\n" unless (-e $opt_s);
die "Number of mutations to introduce must be greater than 1\n" if ($opt_n < 1);
if(!defined $opt_t) { $opt_t = 'SNP'; }
if(!defined $opt_g) { $opt_g = 'none'; }
if(!defined $opt_f) { $opt_f = 'genome'; }
if(!defined $opt_p) { $opt_p = ';'; }
if(!defined $opt_c) { $opt_c = 0; }
if(!defined $opt_r) { $opt_r = 'Parent='; }
if(!defined $opt_o) { $opt_o = "$opt_s-n-$opt_n-t-$opt_t-f-$opt_f.details"; } 
if(!defined $opt_a) { $opt_a = "$opt_s-n-$opt_n-t-$opt_t-f-$opt_f.fasta"; } 
die "-t has to be SNP, DEL, INS or CNV\n" unless ($opt_t =~ m/SNP|DEL|INS|CNV|HET|HETTRIP/);
die "Specify both -c and -f or neither -c specified ($opt_c)\n" if ((defined $opt_c) && (!defined $opt_f));
die "Specify both -c and -f or neither -f specified ($opt_f)\n" if (($opt_f ne 'genome') && (!defined $opt_c));
die "Specify both -c and -f if -t is CNV\n" if (($opt_t eq 'CNV') && ((!defined $opt_c) && ($opt_f eq 'Genome')));
die "Currently, -t has to be 'SNP' or 'HET' if specifying -c/-f\n" if ((($opt_t eq 'DEL') || ($opt_t eq 'INS')) && ($opt_c));

# Specify output files
open my $ofh1, '>', $opt_o or die "Cannot open $opt_o : $!\n";
open my $ofh2, '>', $opt_a or die "Cannot open $opt_a : $!\n";

# Save genome/sequence to in silico mutate
my $fasta = fastafile::fasta_to_struct($opt_s);

# for heterozygous positions, need to add homologous genome
if($opt_t eq 'HET') { $fasta = fastafile::add_homologous_genome($fasta, '_homologous'); }
if($opt_t eq 'HETRIP') { $fasta = fastafile::add_homologous_genome($fasta, '_homologous2'); }

# range for mutating
my $genome_length = $$fasta{'total_length'};
my $range = $genome_length;

# Go through GFF (if specified) to save feature sequences and length
my $gff_struct;
if($opt_f ne 'genome') { 
	$gff_struct = gfffile::gff_to_struct($opt_g, $opt_f, $opt_p, $opt_c, $opt_r);
	($gff_struct, $range) = &extend_features_and_count_length($gff_struct, $opt_t);
}

# Check range is bigger than number of mutations
die "Error: Too many mutations required for range ($opt_n / $range) in feature $opt_f\n" if($range < $opt_n);

# Generate random numbers
my $rand_numbers = &generate_random_numbers($opt_n, $range);

# Modify anywhere in the genome
my $number_of_Ns_over_rand_numbers = 0;
print $ofh1 "Contig\tPosition\tnt\tintroduced $opt_t\tFeature\n";
if($opt_f eq 'genome') {
	warn "Introducing mutations into $opt_g\n";
	my $seq_length_gone_through = 0;
	foreach my $id(sort keys %{$$fasta{'seq'}}) {
		my $seq = $$fasta{'seq'}{$id};
		my @nts = split //, $seq;

		foreach my $numbers(sort {$a<=>$b} keys %{$rand_numbers}) {
			my $random_number_in_id = ($numbers - $seq_length_gone_through);

			next if($random_number_in_id < length($seq));
				
			my $nt = uc $nts[$random_number_in_id];

			# choose a new random number if previous one falls over an 'N' or similarly unsuitable base.
			if($nt !~ m/[ACTG]/) { 
				$number_of_Ns_over_rand_numbers++; 
				($nt, $random_number_in_id) = &find_new_random_number($rand_numbers, $seq_length_gone_through, \@nts);
			}

			my ($newbase);	
			if($opt_t eq 'DEL') { $newbase = ''; }
			if(($opt_t eq 'SNP') || ($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) { $newbase = randomnucleotide($nt); }
			if($opt_t eq 'INS') { $newbase = ($nt . randomnucleotide($nt)); }
			print $ofh1 "$id\t$random_number_in_id\t$nt\t$newbase\t$opt_g\n";
			$nts[$random_number_in_id] = $newbase;
		}
		$seq_length_gone_through += length($seq);
		$$fasta{'seq'}{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t + $number_of_Ns_over_rand_numbers that were picked over Ns in first call\n";
}

# Modify specific places in the genome
if(defined $opt_g) {
	warn "Introducing mutations into the $opt_f of $opt_g\n";
	my $random_numbers_gone_through = 0;

	# if it's CNV, i need the sequence of the cds
	my $cds_sequences;
	if($opt_t eq 'CNV') { $cds_sequences = &save_sequences_from_cds($fasta, $gff_struct); }

	# Go through each of the contigs
	foreach my $id(sort keys %{$$fasta{'seq'}}) {
		my $seq = $$fasta{'seq'}{$id};
		my @nts = split //, $seq;

		# Go through each feature that did not have any overlaps with other features
		foreach my $gene(sort keys %{$$gff_struct{'contig_to_gene'}{$id}}) {

			# SNPs, Hets and Het-triploids
			if($opt_t =~ m/SNP|HET/) { 
				
				my $cds = $$gff_struct{'gene_to_cds'}{$gene};
				my @exons = split /\s/, $cds;
				EXON: foreach my $exon(@exons) {
					my @terminals = split /-/, $exon;
					POSITION: for(my $i=$terminals[0]; $i<$terminals[1]; $i++) {
						$random_numbers_gone_through++;

						next POSITION if(!exists $$rand_numbers{$random_numbers_gone_through});

						my $nt = $nts[$i];

						# tally if random number falls over an 'N'
						if($nt eq 'N') { 
							$number_of_Ns_over_rand_numbers++; 
							next POSITION;
						}	

						my $newbase = &randomnucleotide($nt);
						print $ofh1 "$id\t$i\t$nt\t$newbase\t$gene\n";
						$nts[$i] = $newbase;
					}
				}
			}

			# CNV
			elsif($opt_t eq 'CNV') { 
				my $gene_seq = $$cds_sequences{$gene};
				my @feature_nts = split //, $gene_seq;

				# put the new gene right on the end of the genome!
				my $pre_CNV_contig_length = (scalar(@nts) + 1);
				push @nts, @feature_nts;
				my $post_CNV_contig_length = scalar(@nts);
				print $ofh1 "$id\t$pre_CNV_contig_length-$post_CNV_contig_length\tCNV+\tCNV+\t$gene\n";
			}

			else { die "Error: no know opt_t : $opt_t"; }
		}
		$$fasta{'seq'}{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t ($number_of_Ns_over_rand_numbers Ns found but not changed)\n";
	warn "Random numbers gone through in feature file: $random_numbers_gone_through\n";
}
close $ofh1;

# Print modified genome
fastafile::fasta_struct_print_simple_outfile($fasta);

sub randomnucleotide { 
	my $consbase = shift;
	my @nucs;
	if($consbase eq 'A') { @nucs = ('C','G','T'); }
	if($consbase eq 'C') { @nucs = ('A','G','T'); }
	if($consbase eq 'T') { @nucs = ('A','C','G'); }
	if($consbase eq 'G') { @nucs = ('A','C','T'); }
	return $nucs[rand @nucs];
}

sub generate_random_numbers {
	my ($number_wanted, $range) = @_;
	my %rand_numbers;
	warn "generate_random_numbers: $opt_n / $range\n";
	NUMBERS: while(scalar(keys %rand_numbers) < $number_wanted) {
		my $random_number = int(rand($range)) + 1;
		next NUMBERS if($random_number <= 1);
		$rand_numbers{$random_number} = 1;
	}
	return \%rand_numbers;
}

sub find_new_random_number {
	my ($rand_numbers, $seq_length_gone_through, $nts) = @_;
	my ($new_random_nucleotide, $random_number2);
	while(!defined $new_random_nucleotide) {
		$random_number2 = int(rand(scalar(@{$nts}))) + 1; 
		if(!exists $$rand_numbers{($seq_length_gone_through + $random_number2)}) {
			my $test_nt = $$nts[$random_number2];
			if($test_nt =~ m/[ACTGactg]/) { $new_random_nucleotide = $test_nt; }
		}
	}
	return ($new_random_nucleotide, $random_number2);
}

sub extend_features_and_count_length {
	my ($gff_struct, $mutation_type) = @_;
	my $feature_range = 0;

	warn "extend_features_and_count_length:";
	foreach my $gene(keys %{$$gff_struct{'gene_to_cds'}}) {
		my $contig = $$gff_struct{'gene_to_contig'}{$gene};
		my $cds = $$gff_struct{'gene_to_cds'}{$gene};
		my $strand = $$gff_struct{'gene_to_strand'}{$gene};
		my $cds_length = gfffile::CDS_to_CDS_length($cds);
		$feature_range += $cds_length;

		# extend feature range
		if($mutation_type eq 'HET') { 
			my $hom_gene_id = ($gene . '_homologous');
			$feature_range += $cds_length;
			$$gff_struct{'gene_to_cds'}{$hom_gene_id} = $cds;
			$$gff_struct{'gene_to_contig'}{$hom_gene_id} = $contig;
			$$gff_struct{'gene_to_strand'}{$hom_gene_id} = $strand;
			$$gff_struct{'contig_to_gene'}{$contig}{$hom_gene_id} = 1;
		}
		if($mutation_type eq 'HETRIP') {
			my $hom_gene_id = ($gene . '_homologous2');
			$feature_range += $cds_length;
			$$gff_struct{'gene_to_cds'}{$hom_gene_id} = $cds;
			$$gff_struct{'gene_to_contig'}{$hom_gene_id} = $contig;
			$$gff_struct{'gene_to_strand'}{$hom_gene_id} = $strand;
			$$gff_struct{'contig_to_gene'}{$contig}{$hom_gene_id} = 1;
		}
	}
	return ($gff_struct, $feature_range);
}

sub save_sequences_from_cds {
	my ($fasta, $gff) = @_;
	my %cds_sequences;

	foreach my $gene(sort keys %{$$gff_struct{'gene_to_cds'}}) {
		my $cds_coords = $$gff_struct{'gene_to_cds'}{$gene};
		my $contig = $$gff_struct{'gene_to_contig'}{$gene};
		my $strand = $$gff_struct{'gene_to_strand'}{$gene};
		my $contig_seq = $$fasta{'seq'}{$contig};

		$cds_sequences{$gene} = gfffile::extract_gene_from_coords_new_negative($contig_seq, $cds_coords, $strand); 
	}
	return \%cds_sequences;
}