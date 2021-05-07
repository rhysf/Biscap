#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::Tools::GFF;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF_lines;
use read_GFF;
use read_SAM;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Program:  Comparison of False Discovery Rates (CFDR)
Usage: perl $0 -i <Full list of introduced mutations (.details)> -m <VCF>\n";
our($opt_i, $opt_m);
getopt('im');
die $usage unless ($opt_i && $opt_m);
foreach($opt_i, $opt_m) { die "Cannot open $_: $!\n" unless(-e $_); }

# Output file
my $outfile = ($opt_m . '-CFDR-Summary.tab');

# save IRMS.details
my ($IRMS_hash, $number_of_introduced_mutations, $type_of_introduced_mutation) = &IRMS_details_to_contig_pos_type_to_base_hash($opt_i);

# Counts
my %accuracy_metric_counts;
my @FP_lines;

# Go through VCF and calculate TP/FP rate
my $last_contig;
open my $fh, '<', $opt_m or die "Cannot open $opt_m : $!\n";
warn "Cataloging variants in $opt_m...\n";
VCF: while (my $line = <$fh>) {
   	chomp $line;
	my ($VCF_line) = vcflines::read_VCF_lines($line);
	next VCF if($$VCF_line{'next'} eq 1);
	my $contig = $$VCF_line{'supercontig'};
	my $pos = $$VCF_line{'position'};
	my $ref = $$VCF_line{'reference_VCF_format'};
	my $cons =  $$VCF_line{'consensus_VCF_format'};
	my $base_type = $$VCF_line{'base_type0'};
	my $base1 = $$VCF_line{'0base1'};
	my $base2 = $$VCF_line{'0base2'};

	# Progress
	if((!defined $last_contig) || ($last_contig ne $contig)) {
		$last_contig = $contig;
		warn "\tprocessing $contig...\n";
	} 

	# Is variant over an introduced position
	my ($corresponding_ref, $corresponding_con, $corresponding_con2);
	if(defined $$IRMS_hash{$contig}{$pos}{'ref_base'}) { $corresponding_ref = $$IRMS_hash{$contig}{$pos}{'ref_base'}; }
	if(defined $$IRMS_hash{$contig}{$pos}{'cons_base'}) { $corresponding_con = $$IRMS_hash{$contig}{$pos}{'cons_base'}; }
	if(defined $$IRMS_hash{$contig}{$pos}{'cons_base2'}) { $corresponding_con2 = $$IRMS_hash{$contig}{$pos}{'cons_base2'}; }

	# Determine if true positive or false positive and base type

	# Ref True negative (TN)
	if(($base_type eq 'reference') && (!defined $corresponding_ref)) { $accuracy_metric_counts{'TN'}{$base_type}++; }

	# Variant FP
	elsif(($base_type ne 'reference') && (!defined $corresponding_ref)) {
		$accuracy_metric_counts{'FP'}{$base_type}++;
		push @FP_lines, $line;
	}

	# Variant FP (wrong base)
	elsif((($type_of_introduced_mutation =~ m/SNP/) && ($base_type ne 'snp')) ||
	     (($type_of_introduced_mutation =~ m/HET/) && ($base_type ne 'heterozygous')) ||
             (($type_of_introduced_mutation =~ m/INS/) && ($base_type ne 'deletion')) ||
	     (($type_of_introduced_mutation =~ m/DEL/) && ($base_type ne 'insertion'))) {
		$accuracy_metric_counts{'FP'}{$base_type}++;
		push @FP_lines, $line;
	}

	# SNP TP or FP
	elsif(($type_of_introduced_mutation =~ m/^SNP$/i) && ($base_type eq 'snp')) {
		if(($cons eq $corresponding_ref) && ($ref eq $corresponding_con)) { $accuracy_metric_counts{'TP'}{$base_type}++; }
		else { 
			$accuracy_metric_counts{'FP'}{$base_type}++;
			push @FP_lines, $line;
		}
	}

	# Heterozygous TP or FP
	elsif(($type_of_introduced_mutation =~ m/^HET$/) && ($base_type eq 'heterozygous')) {
		if($corresponding_con2) {
			if(($corresponding_con eq $base1) && ($corresponding_con2 eq $base2)) { $accuracy_metric_counts{'TP'}{$base_type}++; }
			elsif(($corresponding_con eq $base2) && ($corresponding_con2 eq $base1)) { $accuracy_metric_counts{'TP'}{$base_type}++; }
			else { 
				$accuracy_metric_counts{'FP'}{$base_type}++;
				#warn "FP het at position with introduced het: $line\n"; 
			}
		} else {
			if(!defined $corresponding_ref) {
				die "corresponding ref not found for $line\n";
			}

			if(($corresponding_ref eq $base1) && ($corresponding_con eq $base2)) { $accuracy_metric_counts{'TP'}{$base_type}++; }
			elsif(($corresponding_ref eq $base2) && ($corresponding_con eq $base1)) { $accuracy_metric_counts{'TP'}{$base_type}++; }
			else { 
				$accuracy_metric_counts{'FP'}{$base_type}++;
				#warn "FP het at position with introduced het: $line\n"; 
			}
		}
	}
	# Other
	elsif($type_of_introduced_mutation =~ m/^INS$/i) { 
		die "Error: Haven't coded for INS yet, but appears at a correct position: $line\n";
		#$TP_INS++; 
	}
	elsif($type_of_introduced_mutation =~ m/^DEL$/i) { 
		die "Error: Haven't coded for DEL yet, but appears at a correct position: $line\n";
		#$TP_INS++; 
	}
	elsif($type_of_introduced_mutation =~ m/^HETTRIP$/i) {  
		die "Error: Haven't coded for HETTRIP yet, but appears at a correct position: $line\n";
	}
	# Unrecognised
	else { die "Error: Unrecognised mutation type from details ($opt_i): $type_of_introduced_mutation\n"; }
}
close $fh;

# Calculate results
my $all_variant_count = 0;
foreach my $accuracy_type(keys %accuracy_metric_counts) {
	foreach my $base_type(keys %{$accuracy_metric_counts{$accuracy_type}}) {
		$all_variant_count += $accuracy_metric_counts{$accuracy_type}{$base_type};
	}
}

# Print results
open my $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
print $ofh "Type\tBase_Type\tTally\tPercent\n";
warn "Type\tBase_Type\tTally\tPercent\n";
foreach my $accuracy_type(keys %accuracy_metric_counts) {
	foreach my $base_type(keys %{$accuracy_metric_counts{$accuracy_type}}) {
		my $tally = $accuracy_metric_counts{$accuracy_type}{$base_type};
		my $percent;
		if($accuracy_type eq 'TP') { $percent = sprintf("%.2f", (($tally / $number_of_introduced_mutations) * 100)); } 
		else { $percent = 'N/A'; }
		print $ofh "$accuracy_type\t$base_type\t$tally\t$percent\n";
		warn "$accuracy_type\t$base_type\t$tally\t$percent\n";
	}
}
print $ofh "\nFalse Positive $type_of_introduced_mutation\n\n";
foreach my $line(@FP_lines) { print $ofh "$line\n"; }

sub IRMS_details_to_contig_pos_type_to_base_hash {
	my $filename = $_[0];
	my %found;
	my $count = 0;

	# Save mutation type from header
	my $type = `head -1 $filename | awk -F"\t| " '{print \$5}'`;
	if(! $type) { die "Cannot find the introduced muation type from file provided: $filename\n"; }

	# Introduced changes
	warn "Cataloging SNPs in from $filename...\n";
	open my $fh, '<', $filename or die "Cannot open $filename : $!\n";
	DETAILS: while(my $line = <$fh>) {
		chomp $line;
		next DETAILS if ($line =~ m/^Contig\tPosition/);
		my @bits = split /\t/, $line;
		die ("Incorrect number of columns found on line $line in $filename: " . scalar(@bits) . "\n") if (scalar(@bits) ne 5);
		my ($contig, $position, $nt, $introduced_mutation, $feature) = (@bits);
		
		# Heterozygous use an additional non-modified chromosome for simulating reads
		$contig =~ s/\_homologous//g;
		
		# All inserted mutations have their array position so need to +1
		$position++;
		# Insertions should be detected as a deletion so need a -1
		# check against gatk too!
		#if($type =~ m/^INS$/i) { $position -= 1; }

		# Genotype 1/2 for HET
		if(defined $found{$contig}{($position)}{'cons_base'}) {
			$found{$contig}{($position)}{'cons_base2'} = $introduced_mutation;
		}
		else {
			$found{$contig}{($position)}{'ref_base'} = $nt;
			$found{$contig}{($position)}{'cons_base'} = $introduced_mutation;
		}
		$count++;
	}
	close $fh;
	return (\%found, $count, $type);
}
