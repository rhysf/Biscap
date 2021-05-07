#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_Pileup;
use read_FASTA;
use read_VCF_lines;
use call_variants_using_binomial_probabilities;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Program:  Binomial SNP Caller from Pileup (BiSCaP)
Version: 0.12
Usage: perl $0 -p <mpileup> -r <reference.fasta>\n
Optional: -m\tMinimum read depth to be considered a mutation [4]
          -e\tProbability of error (e.g. 0.1 or 0.01) [0.1]
          -q\tRead Quality minimum cut-off (e.g. 10, 20...) [0]
          -n\tSample name for VCF [WGS]\n";
our($opt_p, $opt_r, $opt_m, $opt_e, $opt_q, $opt_n);
getopt('prmelsaqi');
die $usage unless ($opt_p && $opt_r);   
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_e) { $opt_e = 0.1; }
if(!defined $opt_q) { $opt_q = 0; }
if(!defined $opt_n) { $opt_n = 'WGS'; }

# Save Binomial probabilities and check all inputs are valid
my $Cumulative_Binoms = "$Bin/util/Cumulative_Binomial_Probabilities_for_$opt_e.tab";
foreach($opt_p, $opt_r, $Cumulative_Binoms) { die "Cannot open $_: $!\n" unless (-e $_); }
die "-m -e -q need to be numerical\n" if ((($opt_m) !~ /\d+/) || (($opt_e) !~ /\d+/) || (($opt_q) !~ /\d+/)); 
die "-e need to be numerical between 0 and 1\n" if (($opt_e > 1) || ($opt_e < 0));
my $settings = ('m-' . $opt_m . '-e-' . $opt_e . '-q-' . $opt_q);
&warn_settings($settings, $opt_p, $opt_r);
my ($probability, $max_depth) = binomial_probabilities::save_binomial_distributions($Cumulative_Binoms);

# Prepare output files
my $outfile1 = ($opt_p . "-Biscap-$settings.vcf");
my $outfile2 = ($opt_p . "-Biscap-$settings-summary.txt");
open my $ofh, '>', $outfile1 or die "Cannot open $outfile1 : $!\n";
open my $ofh2, '>', $outfile2 or die "Cannot open $outfile2 : $!\n";

# Save the length of reference for VCF header and print VCF header
my ($reference_length) = fastafile::fasta_id_to_seq_length_hash($opt_r);
my $print_header = &print_VCF_header($0, $opt_r, $reference_length, $settings, $opt_n);

# Read from Pileup and determine the state at each base
my ($tally_HOMSNP, $tally_HOMINS, $tally_HOMDEL, $tally_HETSNP, $tally_HETINS, $tally_HETDEL) = (0, 0, 0, 0, 0, 0);
my ($tally_OUTDIS, $tally_HOMAGR, $tally_MINDEP, $tally_UNCHAR, $removed_bases) = (0, 0, 0, 0, 0);
warn "Reading from file: $opt_p...\n";
open my $fh, '<', $opt_p or die "Can't open $opt_p : $!\n";
my @mapping_qualities;
PILEUP: while(my $line=<$fh>) {
	chomp $line;
	next PILEUP if($line =~ m/^\n$"/);

	my @bits = split /\t/, $line;
	my ($contig, $pos, $ref_base, $depth, $aligned, $read_qual) = @bits;
	die "Pileup line not separated by tab or missing column: $line\n" if (@bits != 6);
	$ref_base = uc $ref_base;

	my ($saved_depth, $saved_ref) = ($depth, $ref_base);

	# Save and remove the mapping qualities
	my ($aligned2, $end_of_reads, $map_qualities) = pileupfile::save_and_remove_mapping_quality($aligned);
	push (@mapping_qualities, @{$map_qualities});

	# Tally minimum depth
	if($depth < $opt_m) { $tally_MINDEP++; }

	if(($ref_base ne 'N') && ($depth >= $opt_m)) {

		# Count bases 
		my $bases = pileupfile::count_bases($ref_base, $aligned2);

		# Remove N's
		if(defined $$bases{'N'}) {
			$depth -= $$bases{'N'};
			delete $$bases{'N'};
		}
		my $variants = &find_most_common_bases($bases);

		# Outside of calculated distributions, use up to max depth as representative 
		if($depth > $max_depth) {
			$variants = &scale_counts_to_max_depth($variants, $depth, $max_depth);
			$depth = $max_depth;
		}

		# Compare probabilities
		my ($consensus, $genotype, $info) = binomial_probabilities::compare_probability($variants, $depth, $opt_e, $opt_m, $ref_base, $probability); 
		warn "(line, consensus and genotype): $line = $consensus and $genotype based on $variants, $depth, $opt_e, $opt_m, $ref_base, $probability\n";
		#warn Dumper($variants);

		# Calculate mean mapping quality, base quality and max mapping quality 
		my ($mean_mapping_quality, $mean_base_quality) = (0,0);
		if($opt_q) { ($mean_mapping_quality, $mean_base_quality) = &calculate_mean_mapping_qual_and_base_qual(\@mapping_qualities, $read_qual); }

		### Record the polymorphism
		# Reference base
		if($genotype eq 0) { 
			$tally_HOMAGR++; 
			next PILEUP; 
		}

		# base-change
		my $VCF_line = "$contig\t$pos\t.\t$ref_base\t$consensus\t$mean_base_quality\t";
		$VCF_line .= ".\t$info\tGT:DP:MMQ\t$genotype:$saved_depth:$mean_mapping_quality";
		print $ofh "$VCF_line\n";

		# Tally polymorphisms
		my $tally_print = &tally_and_print_variants($ref_base, $info, $consensus, $genotype, $line);
	}

	# Remove end of read mapping scores
	foreach(sort { $b <=> $a } @{$end_of_reads}) { 
		if(defined $mapping_qualities[$_]) {
			splice (@mapping_qualities, $_, 1); 
		}
	}
}
close $fh;

# Summary of SNP calling
my $warn_summary = &warn_summary($tally_HOMSNP, $tally_HOMINS, $tally_HOMDEL, 
	$tally_HOMAGR, $tally_HETSNP, $tally_HETINS, $tally_HETDEL, $tally_UNCHAR, 
	$tally_OUTDIS, $tally_MINDEP, $max_depth, $removed_bases, $opt_q, $opt_m);
print $ofh2 "$warn_summary\n";

# Close output files
close $ofh; 
close $ofh2; 

sub scale_counts_to_max_depth {
	my ($variants, $depth, $max_depth) = @_;
	my %new_variants;
	foreach my $order(sort keys %{$variants}) {
		foreach my $base(sort keys %{$$variants{$order}}) {
			my $tally = $$variants{$order}{$base};
			my $scale = sprintf "%.0f", (($tally / $depth) * $max_depth);
			$new_variants{$order}{$base} = $scale;
		}
	}
	return (\%new_variants);
}

sub find_most_common_bases {
	my ($bases) = $_[0];
	my %variants;

	# Save the order of the most common bases
	my $i=1;
	foreach my $base (sort { $$bases{$b} <=> $$bases{$a} } keys %{$bases}) {	
		$variants{$i}{$base} = $$bases{$base};
		warn "most common base: $i -> $base = $$bases{$base}\n";
		$i++;
	}
	return (\%variants);
}

sub calculate_mean_mapping_qual_and_base_qual {
	my ($mapping_qualities, $read_qual) = @_;
	my $max_mapping_quality = (sort { $b <=> $a } @{$mapping_qualities})[0];
	my ($sum_base_qual, $sum_map_qual) = (0, 0);
	foreach(@{$mapping_qualities}) { $sum_map_qual += $_; }
	
	# This might have to be split into an array first
	foreach(@{$read_qual}) { $sum_base_qual += $_; }

	my $mean_mapping_quality = int($sum_map_qual / scalar(@mapping_qualities));
	my $mean_base_quality = int($sum_base_qual / scalar(@{$read_qual}));
	return ($mean_mapping_quality, $mean_base_quality)
}

sub MAX {
	my (@values) = @_;
	my ($position_of_number, $big_number) = (0, 0);
	my @positions_of_numbers;
	for(my $i=0; $i<(scalar(@values) - 1); $i++) {
		if(($values[$i] > $values[$i+1]) && ($values[$i] > $big_number)) {  $position_of_number = $i; $big_number = $values[$i]; }
		elsif(($values[$i+1] > $values[$i]) && ($values[$i+1] > $big_number)) { $position_of_number = ($i+1); $big_number = $values[($i+1)]; }
	}
	return $position_of_number;
}

sub print_VCF_header {
	my ($source, $reference, $FASTA_length_ref, $settings, $sample_name) = @_;
	my @reference_name_parts = split /\//, $reference;
	my $reference_file_name = $reference_name_parts[(scalar(@reference_name_parts) - 1)];
	my @source_name_parts = split /\//, $source;
	my $source_file_name = $source_name_parts[(scalar(@source_name_parts) - 1)];
	print $ofh ('##fileformat=VCFv4.1' . "\n");
	print $ofh ('##source=' . "$source_file_name\n");
	print $ofh ('##' . "$settings\n");
	print $ofh ('##reference=' . "$reference_file_name\n");
	foreach my $contig_name(sort keys(%{$FASTA_length_ref})) {
		my $length = $$FASTA_length_ref{$contig_name};
		print $ofh ('##contig=<ID=' . $contig_name . ',length=' . $length . '>' . "\n");
	}
	print $ofh ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' . "\n");
	print $ofh ('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' . "\n");
	print $ofh ('##FORMAT=<ID=MMQ,Number=1,Type=Integer,Description="Mean mapping quality score">' . "\n");
	print $ofh ('##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">' . "\n");
	print $ofh ('#' . "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n");	
	return 1;
}

sub warn_settings {
	my ($param_full, $pileup_name, $ref) = @_;
	my @par = split /-/, $param_full;
	my ($min_depth, $error_prob, $qual) = ($par[1], $par[3], $par[5]);
	my $settings = "Settings:\n";
	$settings .= "Pileup:          $pileup_name\n";
	$settings .= "Reference FASTA: $ref\n";
	$settings .= "Min Read Depth:  $min_depth\n";
	$settings .= "Predicted Error: $error_prob\n";
	$settings .= "Quality Filter:  $qual\n";
	warn "$settings\n";
	return $settings;	
}

sub warn_summary { 
	my ($HOMSNP, $HOMINS, $HOMDEL, $HOMAGR, $HETSNP, $HETINS, $HETDEL, $UNCHAR, $OUTDIS, $MINDEP, $MAXDEP, $REMBAS, $optq, $optm) = @_;
	my $summary = "Finished\n\n-----------------------------------------\n";
	$summary .= "Homozygous SNPs identified: $HOMSNP\n";
	$summary .= "Homozygous INSs identified: $HOMINS\n";
	$summary .= "Homozygous DELs identified: $HOMDEL\n";
	$summary .= "Homozygous agree with ref : $HOMAGR\n";
	$summary .= "Heterozygous SNPs identified: $HETSNP\n"; 
	$summary .= "Heterozygous INSs identified: $HETINS\n"; 
	$summary .= "Heterozygous DELs identified: $HETDEL\n"; 
	$summary .= "-----------------------------------------\n";
	$summary .= "Unable to resolve (number of bases): $UNCHAR\n";
	
	if($tally_OUTDIS >= 1) { 
		$summary .= "Bases unresolved due to being greater than lookup table ($MAXDEP): $OUTDIS (increase by using GBiD.pl)\n"; 
	}
	$summary .= "\nIgnored bases lower than Quality $optq: $REMBAS\n";
	$summary .= "Ignored bases lower than min depth $optm: $MINDEP\n";
	warn "$summary\n";
	return 1;
}

sub tally_and_print_variants {
	my ($ref_base, $info, $consensus, $genotype, $line) = @_;
	# Insertions
	if(($info =~ m/INDEL/) && ($consensus =~ m/[ACTG][ACTG]/)) {
		if($genotype =~ m/\//) {
			$tally_HETINS++; 
		} else {
			$tally_HOMINS++; 
		}
	}
	# Deletions
	if(($info =~ m/INDEL/) && ($ref_base =~ m/[ACTG][ACTG]/)) {
		if($genotype =~ m/\//) {
			$tally_HETDEL++; 
		} else {
			$tally_HOMDEL++; 
		}
	}
	# SNP 
	if(($info !~ m/INDEL/) && ($consensus !~ m/[ACTG][ACTG]/) && ($ref_base !~ m/[ACTG][ACTG]/)) {
		if($genotype =~ m/\//) {
			$tally_HETSNP++; 
		} else {
			$tally_HOMSNP++; 
		}
	}
	return 1;
}
