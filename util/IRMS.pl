#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_GFF;
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

# Save genome/sequence for mutating (make homologous version if heterozygous)
warn "Saving $opt_s...\n";
my ($genome_fna, $genome_size, $genome_length) = &save_sequence_and_lengths_from_FASTA($opt_s, $opt_c, $opt_t);

# Go through GFF (if specified) to save feature sequneces and length
my  ($CDS_info, $strand_info, $feature_positions);
my $feature_range = 0;
if($opt_f ne 'genome') { 
	($CDS_info, $strand_info) = gfffile::gff_to_contig_parent_to_cds_hash($opt_g, $opt_f, $opt_p, $opt_c, $opt_r); 
	$feature_positions = gfffile::gff_to_contig_parent_to_start_stop_string($opt_g, $opt_f, $opt_p, $opt_c, $opt_r);
	foreach my $contig(keys %{$CDS_info}) { 
		foreach my $gene(keys %{$$CDS_info{$contig}}) {
			my $seq = $$CDS_info{$contig}{$gene};

			# Extend feature range
			$feature_range += length($seq);
			if($opt_t =~ m/CNV|HET/) { $feature_range += length($seq); }
			if($opt_t =~ m/HETTRIP/) { $feature_range += length($seq); }

			# Save extra feature positions
			if($opt_t =~ m/HET|HETTRIP/) { $$feature_positions{($contig . '_homologous')}{$gene} .= $$feature_positions{$contig}{$gene}; }
			if($opt_t eq 'HETTRIP') { $$feature_positions{($contig . '_homologous2')}{$gene} .= $$feature_positions{$contig}{$gene}; }
		}
	}
}

# Which range to use
my $range = $genome_length;
if($feature_range ne 0) { $range = $feature_range; }

# Check range is bigger than number of mutations
if($opt_t eq 'CNV') { die "Too many mutations required ($opt_n) for number of features ($opt_f) found for mutating ($feature_range)\n" if($feature_range < $opt_n); } 
else { die "Too many mutations required ($opt_n) for amount of sequence (nt) found for mutating ($genome_length)\n" if($genome_length < $opt_n); }

# Random positions generator
my %rand_numbers;
warn "Generating $opt_n / $range $opt_t...\n";
while(scalar(keys %rand_numbers) < $opt_n) {
	my $random_number = int(rand($range)) + 1; 
	$rand_numbers{$random_number} = 1;
}

# Modify anywhere in the genome
my $number_of_Ns_over_rand_numbers = 0;
print $ofh1 "Contig\tPosition\tnt\tintroduced $opt_t\tFeature\n";
if($opt_f eq 'genome') {
	warn "Introducing mutations into $opt_g\n";
	my $seq_length_gone_through = 0;
	foreach my $id(sort keys %{$genome_fna}) {
		my $seq = $$genome_fna{$id};
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
				print $ofh1 "$id\t$random_number_in_id\t$nt\t$newbase\t$opt_g\n";
				$nts[$random_number_in_id] = $newbase;
			}
		}
		$seq_length_gone_through += length($seq);
		$$genome_fna{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t + $number_of_Ns_over_rand_numbers that were picked over Ns in first call\n";
}

# Modify specific places in the genome
if($opt_c) {
	warn "Introducing mutations into the $opt_f of $opt_g\n";
	my $random_numbers_gone_through = 0;
	# Go through each of the contigs
	foreach my $id(sort keys %{$genome_fna}) {
		my $seq = $$genome_fna{$id};
		my @nts = split //, $seq;
		# Go through each feature that did not have any overlaps with other features
		foreach my $feature_names(sort keys %{$$feature_positions{$id}}) {
			
			if(($opt_t eq 'SNP') || ($opt_t eq 'HET') || ($opt_t eq 'HETTRIP')) { 
				my $feature_positions = $$feature_positions{$id}{$feature_names};
				my @seperate_positions = split /\s/, $feature_positions;
				foreach(@seperate_positions) {
					my @terminals = split /-/, $_;
					for(my $i=$terminals[0]; $i<$terminals[1]; $i++) {
						$random_numbers_gone_through++;
						if(exists $rand_numbers{$random_numbers_gone_through}) {
							my $nt = $nts[$i];
							if($nt eq 'N') { $number_of_Ns_over_rand_numbers++; }
							else {	
								my $newbase = &randomnucleotide($nt);
								print $ofh1 "$id\t$i\t$nt\t$newbase\t$feature_names\n";
								$nts[$i] = $newbase;
							}
						}
					}
				}
			}
			if($opt_t eq 'CNV') { 
				$random_numbers_gone_through++;
				if(exists $rand_numbers{$random_numbers_gone_through}) {
					my @feature_nts = split //, $$CDS_info{$id}{$feature_names};
					
					#warn "my feature nts = @feature_nts \n";
					
					my $pre_CNV_contig_length = (scalar(@nts) + 1);
					push @nts, @feature_nts;
					my $post_CNV_contig_length = scalar(@nts);
					print $ofh1 "$id\t$pre_CNV_contig_length-$post_CNV_contig_length\tCNV+\tCNV+\t$feature_names\n";
					
					#warn "pre = $pre_CNV_contig_length and post = $post_CNV_contig_length\n";	
				}
			}
		}
		$$genome_fna{$id} = join('',@nts);
	}
	my $number_of_successful_positions = ($opt_n - $number_of_Ns_over_rand_numbers);
	warn "Introduced $number_of_successful_positions $opt_t ($number_of_Ns_over_rand_numbers Ns found but not changed)\n";
	warn "$random_numbers_gone_through Random numbers gone through in feature file\n";
}
close $ofh1;

# Print modified genome
warn "Printing modified $opt_s...\n";
foreach my $contig(sort keys %{$genome_fna}) { 
	print $ofh2 ">$contig\n"; 
	my $seq = $$genome_fna{$contig};
	my @seq_bits = split //, $seq;
	my ($count, $contig_count) = (0, 0);
	foreach(@seq_bits) {
		$count++; $contig_count++;
		print $ofh2 "$_";
		if(($count eq 60) && (defined $seq_bits[$contig_count + 1])) { print $ofh2 "\n"; $count=0 }
	}
	print $ofh2 "\n";
	$contig_count = 0;
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

sub save_sequence_and_lengths_from_FASTA {
	my ($fasta, $gff, $mutation_type) = @_;

	# Save FASTA
	my (%genome_fna, %genome_size);
	my $total_length = 0;
	my $inseq = Bio::SeqIO->new('-file' => "<$fasta",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) {
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		$genome_fna{$id} = $seq;
		$genome_size{$id} = length($seq);
		$total_length += length($seq);

		# For hets
		if($mutation_type =~ m/HET|HETTRIP/) {
			$genome_fna{($id . '_homologous')} = $seq;
			$genome_size{($id . '_homologous')} = length($seq);
			$total_length += length($seq);
		}
		if($mutation_type eq 'HETTRIP') { 
			$genome_fna{($id . '_homologous2')} = $seq;
			$genome_size{($id . '_homologous2')} = length($seq);
			$total_length += length($seq);
		}
	}
	return (\%genome_fna, \%genome_size, $total_length);
}
