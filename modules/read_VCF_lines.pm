package vcflines;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub read_VCF_lines {
	my $VCF_line = $_[0];
	my @bits = split /\t/, $VCF_line;
	my %VCF_info;

	# Save header info including isolate names
	if($VCF_line =~ m/^\#/) { 
		my $VCF_struct = &VCF_header_to_struct($VCF_line, \%VCF_info); 
		return $VCF_struct; 
	} else {
		$VCF_info{'next'}=0;
		$VCF_info{'header'}='N';
	}

	# Initial quality check
	if(@bits < 9) {
		warn "$0: Bad VCF with < 9 columns: $VCF_line\n";
		$VCF_info{'next'}=1; 
		return \%VCF_info;
	}

	# Multi VCF
	if(@bits > 10) {
		for(my $i=9; $i < scalar(@bits); $i++) {
			my $sample_info_id = ('sample_info' . ($i - 9));
			#warn "$sample_info_id = $i = $bits[$i]\n";
			$VCF_info{$sample_info_id} = $bits[$i];
		}
	}

	# Parts continued
	$VCF_info{'supercontig'}          = $bits[0];
	$VCF_info{'position'}             = $bits[1];
	$VCF_info{'id'}                   = $bits[2];
	$VCF_info{'reference_VCF_format'} = $bits[3];
	$VCF_info{'consensus_VCF_format'} = $bits[4];
	$VCF_info{'cons_qual'}            = $bits[5];
	$VCF_info{'filter'}               = $bits[6];
	$VCF_info{'info'}                 = $bits[7];
	$VCF_info{'format'}               = $bits[8];

	# Split format parts, and save available names 
	my @format_parts = split /:/, $VCF_info{'format'};
	my %format_part_ids_available;
	foreach(@format_parts) { $format_part_ids_available{$_} = 1; }
	$VCF_info{'number_of_samples'} = (scalar(@bits) - 9);

	# Sample_info deals with Multi VCF as well
	for(my $i=9; $i < scalar(@bits); $i++) {

		# Save sample_info[isolate number]
		my $isolate_number = ($i - 9);
		my $sample_info_id = ('sample_info' . $isolate_number);
		$VCF_info{$sample_info_id} = $bits[$i];
		my @sample_info_parts = split /:/, $VCF_info{$sample_info_id};

		# Check format parts matches sample info parts. Error if less format parts than sample info parts
		die "$VCF_line: scalar(@sample_info_parts) < scalar(@format_parts) (sample info < format parts). Have not coded for this eventuality!\n" if(scalar(@format_parts) < scalar(@sample_info_parts));
		# Subset sample_info_parts by format_parts
		if(scalar(@format_parts) > scalar(@sample_info_parts)) {
			my @reduced_sample_info_parts;
			for(my $i=0; $i<scalar(@format_parts); $i++) {
				push @reduced_sample_info_parts, $sample_info_parts[$i];
			}
			@sample_info_parts = @reduced_sample_info_parts;
		}
		die "Format and Sample Info do not match: $VCF_line\n" if(scalar(@format_parts) ne scalar(@sample_info_parts));

		# Save genotype, depth etc.
		my ($GT_id, $DP_id, $base_type_ID, $amb_char_ID, $pid) = ("GT$isolate_number", "DP$isolate_number", "base_type$isolate_number", "amb_char$isolate_number", "PID$isolate_number");
		for(my $f=0; $f<scalar(@format_parts); $f++) {
			my $format_part = $format_parts[$f];
			my $sample_part = $sample_info_parts[$f];
			my $format_part_id = ($format_part . $isolate_number);
			die "$format_part_id overwrites known format id on $VCF_line\n" if(defined $format_part_ids_available{$format_part_id});
			$VCF_info{$format_part_id} = $sample_part;
		}
		die "Unable to find genotype ($GT_id): $VCF_line\n" if(!defined $VCF_info{$GT_id});
		if(!defined $VCF_info{$DP_id}) { $VCF_info{$DP_id} = '?'; }

		# Check alleles are different
		if($VCF_info{$GT_id} =~ m/\/|\|/) {
			my @allele_parts = split /\/|\|/, $VCF_info{$GT_id};
			die "Have not coded for multiple alleles for $VCF_line in check alleles are different\n" if(scalar(@allele_parts) ne 2);
			if($allele_parts[0] eq $allele_parts[1]) { 
				$VCF_info{$GT_id} = $allele_parts[0]; 
			}
		}

		# Determine base (base2 only if diploid)
		my ($base1, $base2, $base_type) = &VCF_struct_determine_bases_and_base_type(\%VCF_info, $GT_id);
		$VCF_info{$base_type_ID}= $base_type;
		$VCF_info{($isolate_number . 'base1')} = $base1;
		$VCF_info{($isolate_number . 'base2')} = $base2;
		
		# Ambiguity character
		if($VCF_info{$base_type_ID} eq 'heterozygous') { $VCF_info{$amb_char_ID} = &get_ambiguity_char($base1, $base2); }

		# GATK phasing
		if(defined $VCF_info{$pid}) {
			$VCF_info{'phased'} = 1;
			$VCF_info{'phase_group'} = $VCF_info{$pid}; # Not necessary. It's already saved. but for now, keep it so it works with my phasing
			#die "new code. I've found PID: $VCF_info{$pid} . All good?\n";
		}
	}	

	# Phased (by my scripts)
	if($VCF_info{'info'} =~ m/PHASE\=(.+)/) {
		$VCF_info{'phased'} = 1;
		my $phase_name = $1;
		$phase_name =~ s/\;//g;
		$VCF_info{'phase_group'} = $phase_name;
		#$VCF_info{'base_type'}='phased'; 
	}

	# Return
	return \%VCF_info;
}

sub get_ambiguity_char {
	my ($base1, $base2) = @_;
	my $ambiguity_char;
				
	# K
	if(($base1 eq 'T') && ($base2 eq 'G')) { $ambiguity_char = 'K'; }
	elsif(($base1 eq 'G') && ($base2 eq 'T')) { $ambiguity_char = 'K'; }
							
	# M
	elsif(($base1 eq 'A') && ($base2 eq 'C')) { $ambiguity_char = 'M'; }
	elsif(($base1 eq 'C') && ($base2 eq 'A')) { $ambiguity_char = 'M'; }
										
	# R
	elsif(($base1 eq 'A') && ($base2 eq 'G')) { $ambiguity_char = 'R'; }
	elsif(($base1 eq 'G') && ($base2 eq 'A')) { $ambiguity_char = 'R'; }
													
	# Y
	elsif(($base1 eq 'T') && ($base2 eq 'C')) { $ambiguity_char = 'Y'; }
	elsif(($base1 eq 'C') && ($base2 eq 'T')) { $ambiguity_char = 'Y'; }
																
	# S
	elsif(($base1 eq 'G') && ($base2 eq 'C')) { $ambiguity_char = 'S'; }
	elsif(($base1 eq 'C') && ($base2 eq 'G')) { $ambiguity_char = 'S'; }
																			
	# W
	elsif(($base1 eq 'A') && ($base2 eq 'T')) { $ambiguity_char = 'W'; }
	elsif(($base1 eq 'T') && ($base2 eq 'A')) { $ambiguity_char = 'W'; }
	return $ambiguity_char;
}

sub resolve_ambiguity_char {
	my $ambiguity_char = $_[0];
	my ($base1, $base2);
	if($ambiguity_char eq 'K') { $base1='T'; $base2='G'; }
	if($ambiguity_char eq 'M') { $base1='A'; $base2='C'; }
	if($ambiguity_char eq 'R') { $base1='A'; $base2='G'; }
	if($ambiguity_char eq 'Y') { $base1='T'; $base2='C'; }
	if($ambiguity_char eq 'S') { $base1='G'; $base2='C'; }
	if($ambiguity_char eq 'W') { $base1='A'; $base2='T'; }
	return ($base1, $base2);
}

sub transistion_or_transversion {
	my ($base1, $base2) = @_;
	my $type;

	# Transition 
	if(($base1 =~ m/G|g/) && ($base2 =~ m/A|a/)) { $type = 'Ts: G:C<A:T'; }
	if(($base1 =~ m/C|c/) && ($base2 =~ m/T|t/)) { $type = 'Ts: G:C<A:T'; }
	if(($base1 =~ m/A|a/) && ($base2 =~ m/G|g/)) { $type = 'Ts: A:T<G:C'; }
	if(($base1 =~ m/T|t/) && ($base2 =~ m/C|c/)) { $type = 'Ts: A:T<G:C'; }

	# Transversion
	if(($base1 =~ m/T|t/) && ($base2 =~ m/G|g/)) { $type = 'Tv: T:A<G:C'; }
	if(($base1 =~ m/A|a/) && ($base2 =~ m/C|c/)) { $type = 'Tv: T:A<G:C'; }
	if(($base1 =~ m/C|c/) && ($base2 =~ m/A|a/)) { $type = 'Tv: C:G<A:T'; }
	if(($base1 =~ m/G|g/) && ($base2 =~ m/T|t/)) { $type = 'Tv: C:G<A:T'; }

	if(($base1 =~ m/T|t/) && ($base2 =~ m/A|a/)) { $type = 'Tv: T:A<A:T'; }
	if(($base1 =~ m/A|a/) && ($base2 =~ m/T|t/)) { $type = 'Tv: T:A<A:T'; }
	if(($base1 =~ m/C|c/) && ($base2 =~ m/G|g/)) { $type = 'Tv: C:G<G:C'; }
	if(($base1 =~ m/G|g/) && ($base2 =~ m/C|c/)) { $type = 'Tv: C:G<G:C'; }

	if(!defined $type) { die "$base1 -> $base2 not defined as either Ts or Tv in sub routine\n"; }
	return $type;
}

# Replace with this!
sub get_phase_metrics_new {
	my ($phase_info, $current_phase_group, $current_position, $current_supercontig) = @_;

	# Entend the previous phase group
	if(($$phase_info{'previous_group'}) && ($$phase_info{'previous_group'} eq $current_phase_group)) {
		$$phase_info{'previous_last_position'} = $current_position;

		# Remove save
		if(defined $$phase_info{'save'}) { delete $$phase_info{'save'}; }
	}

	# New phase group
	else {
		# Save
		if(defined $$phase_info{'previous_group'}) {
			$$phase_info{'save'}{'group'} = $$phase_info{'previous_group'};
			$$phase_info{'save'}{'supercontig'} = $$phase_info{'previous_supercontig'};
			$$phase_info{'save'}{'first_position'} = $$phase_info{'previous_first_position'};
			$$phase_info{'save'}{'last_position'} = $$phase_info{'previous_last_position'};
		}

		# New
		$$phase_info{'previous_group'} = $current_phase_group;
		$$phase_info{'previous_supercontig'} = $current_supercontig;
		$$phase_info{'previous_first_position'} = $current_position;
		$$phase_info{'previous_last_position'} = $current_position;
	}
	return $phase_info;
}

sub get_phase_metrics {
	my ($current_phase_name, $position, $prev_phase_group, $first_pos, $last_pos) = @_;
	my ($length_of_last_phase_group) = (0,0);

	# Extend previous phase group
	if($current_phase_name eq $prev_phase_group) {
		#warn "phase group $current_phase_name extends to $position \n";
		$last_pos = $position;
	}

	# New phase group / Initialise. Save previous phase length
	else {
		#warn "new Group: $last_pos, $first_pos, $length_of_last_phase_group\n";
		$length_of_last_phase_group = ($last_pos - $first_pos);
		
		# change phase group to new one
		$prev_phase_group = $current_phase_name;
		$first_pos = $position;
		$last_pos = $position;
	}
	return ($prev_phase_group, $first_pos, $last_pos, $length_of_last_phase_group);
}

sub summarise_phases {
	my $phase_sizes = $_[0];
	my ($min, $max, $mean, $total, $num) = (0,0,0);
	foreach my $length(sort { $a <=> $b } keys %{$phase_sizes}) {
		my $tally = $$phase_sizes{$length};
		if(!defined $min) { $min = $length; }
		if($length > $max) { $max = $length; }
		$num += $tally;
		$total += ($length * $tally);
		#warn "in summarise_phases with $length and $tally\n";
	}
	$mean = ($total / $num);
	return ($min, $max, $mean, $total, $num);
}

sub print_VCF_header {
	my ($line, $printing_option) = @_; 
	print "$line\n" if($printing_option eq 'vcf');
	print "$line\n" if($printing_option eq 'header');
	return 1;
}

sub print_VCF_lines {
	my ($line, $printing_option, $feature, $base_type, $outfile) = @_;
	if(($feature eq 'all') || ($feature eq $base_type)) { 
		print "$line\n" if($printing_option eq 'vcf');
		print "$line\n" if($printing_option eq 'no_header');
		if($printing_option eq 'split_by_contig') { 
			open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n";
			print OUT "$line\n";
			close OUT;
		}
	}
	return 1;
}

# New modules for looping and saving info
sub parse_VCF {
	my ($vcf_file, $printing_option, $feature, $summarise_per_contig) = @_;
	my (%counts, %types_found);
	# Printing options (vcf, split_by_contig, no_header, header, opt_f, none) [none]\n
	
	# local (temp) phase info
	my ($previous_phase_group, $first_phase_group_position, $last_phase_group_position, $length_of_last_phase_group) = ('NA',0,0,0);

	# Go through the VCF 
	open IN1, "<$vcf_file" or die "Cannot open $vcf_file\n";
	warn "Reading $vcf_file...\n";
	VCF1: while (my $line = <IN1>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);
		
		# Print 
		# Header
		if($$VCF_line{'header'} eq 'Y') {
			&print_VCF_header($line, $printing_option);
			next VCF1;
		}
		die "update parse_VCF" if($$VCF_line{'number_of_samples'} > 1);

		# Print non-header VCF lines
		last VCF1 if(($printing_option eq 'header') && ($$VCF_line{'header'} eq 'N'));
		&print_VCF_lines($line, $printing_option, $feature, $$VCF_line{'base_type0'}, "$vcf_file-$$VCF_line{'supercontig'}");

		# Tally variants
		my $hash_contig_name = "all";
		if($summarise_per_contig eq 'y') { $hash_contig_name = $$VCF_line{'supercontig'}; }
		$counts{$vcf_file}{$hash_contig_name}{$$VCF_line{'base_type0'}}++; 
		$types_found{'Variant_types'}{$$VCF_line{'base_type0'}} = 1;

		# Special types
		if(($$VCF_line{'base_type0'} eq 'heterozygous') && ($$VCF_line{'info'} =~ m/TRIPROB|TRIALLELIC/)) { $counts{$vcf_file}{$hash_contig_name}{$$VCF_line{'info'}}++; } 
		if(($$VCF_line{'base_type0'} eq 'snp') && ($$VCF_line{'info'} eq 'phased')) { $counts{$vcf_file}{$hash_contig_name}{'Phased_SNP'}++; }
		if(($$VCF_line{'base_type0'} eq 'heterozygous') && ($$VCF_line{'info'} eq 'phased')) { $counts{$vcf_file}{$hash_contig_name}{'Phased_HET'}++; }
		if(defined $$VCF_line{'amb_char0'}) { $types_found{'Type_of_Het'}{$$VCF_line{'amb_char0'}}{$vcf_file}{$hash_contig_name}++; }
		next VCF1 if($$VCF_line{'base_type0'} eq 'ambigious');
		next VCF1 if($$VCF_line{'consensus_VCF_format'} eq '.');

		# Phase metrics
		if($$VCF_line{'phased'}) {
			($previous_phase_group, $first_phase_group_position, $last_phase_group_position, $length_of_last_phase_group) = &get_phase_metrics($$VCF_line{'phase_group'}, $$VCF_line{'position'}, $previous_phase_group, $first_phase_group_position, $last_phase_group_position);
			#warn "Phased line = $previous_phase_group, $first_phase_group_position, $last_phase_group_position, $length_of_last_phase_group\n";
			if($length_of_last_phase_group ne 0) { $types_found{'Phase_group_lengths'}{$vcf_file}{$hash_contig_name}{$length_of_last_phase_group}++; }
			$length_of_last_phase_group = 0;
		}

		# Ts/Tv
		if(($$VCF_line{'GT0'} eq 1) && ((length($$VCF_line{'reference_VCF_format'}) eq 1) && (length($$VCF_line{'consensus_VCF_format'}) eq 1))) {
			my $type = &transistion_or_transversion($$VCF_line{'reference_VCF_format'}, $$VCF_line{'consensus_VCF_format'});
			$types_found{'TsTv'}{$type}{$vcf_file}{$hash_contig_name}++;
		}
	}
	close IN1;
	return (\%counts, \%types_found);
}

sub save_VCF_to_hash {
	my ($file, $settings) = @_;
	warn "Saving VCF file: $file\n";
	my ($VCF_header);
	my %VCF;
	open IN1, "<$file" or die "Cannot open VCF: $file: $!";
	VCF1: while (my $line = <IN1>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Save the header of VCF. Ignore ambigious sites
		if($$VCF_line{'next'} eq 1) {
   			if ($line =~ m/^\#/) { $VCF_header .= "$line\n"; }
			next VCF1;
		}
		die "update save_VCF_to_hash" if($$VCF_line{'number_of_samples'} > 1);

		my $newline = $line;
		if($settings eq "for_phasing") {
   			# Save VCF (removing info column that is used for phasing), and only include alleles and depth for format
			$newline = "$$VCF_line{'supercontig'}\t$$VCF_line{'position'}\t$$VCF_line{'id'}\t$$VCF_line{'reference_VCF_format'}\t";
			$newline .= "$$VCF_line{'consensus_VCF_format'}\t$$VCF_line{'cons_qual'}\t$$VCF_line{'filter'}\t\tGT:DP\t$$VCF_line{'GT0'}:$$VCF_line{'DP0'}";
		} 
		elsif($settings eq "phased") {
			next VCF1 if(!defined $$VCF_line{'phased'});
		}
		elsif($settings eq "ignore_reference") {
			next VCF1 if($$VCF_line{'base_type0'} eq 'reference');
		}
		else { }
	   	$VCF{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}} = $newline;
	}
	close IN1;
	
	if($settings eq "for_phasing") {
		# Make new Phase block info
		$VCF_header .= "##INFO=<ID=PHASE,Number=1,Type=Flag,Description=\"Indicates a section that is phased\">\n";
		$VCF_header .= "##INFO=<ID=PHASEINFO,Number=1,Type=Flag,Description=\"Indicates why a section is not phased\">\n";
	}
	return (\%VCF, $VCF_header);
}

sub save_consensus_hash_from_vcf {
	my $file = $_[0];
	warn "Reading $file...\n";
	my %polymorphisms;
	open IN1, "<$file" or die "Cannot open $file\n";
	VCF1: while (my $line = <IN1>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);
		next VCF1 if($$VCF_line{'next'} eq 1);
		die "update save_consensus_hash_from_vcf\n" if($$VCF_line{'number_of_samples'} > 1);

		# Ignore indels
		next VCF1 if($$VCF_line{'base_type0'} =~ m/insertion|deletion/);

		# Hom and Het SNPs
		my $polymorphism;
		if($$VCF_line{'base_type0'} eq 'snp') { $polymorphism = $$VCF_line{'consensus_VCF_format'}; } 
		if(defined $$VCF_line{'amb_char0'}) { $polymorphism = $$VCF_line{'amb_char'}; }

		# Save
		$polymorphisms{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}} = $polymorphism;
	}
	close IN1;
	return \%polymorphisms;
}

############### Local subroutines
sub VCF_header_to_struct {
        my ($VCF_line, $VCF_struct) = @_;
        my @bits = split /\t/, $VCF_line;
        $$VCF_struct{'next'}=1; 
        $$VCF_struct{'header'}='Y';
        if($VCF_line =~ m/^\#CHROM\tPOS\tID\tREF/) {
                for(my $i=9; $i < scalar(@bits); $i++) {
                        $$VCF_struct{'isolate_names'}{($i - 9)} = $bits[$i];
                }
        }
        return $VCF_struct;
}

sub VCF_struct_determine_bases_and_base_type {
        my($VCF_struct, $GT_id) = @_;
        my ($base1, $base2, $base_type);
        $base2 = 'None';

        # ambiguous
        if(($$VCF_struct{'reference_VCF_format'} eq 'N') || ($$VCF_struct{'consensus_VCF_format'} eq 'N') || ($$VCF_struct{$GT_id} eq '.')) { 
                $base1 = 'N';
                $base_type = 'ambiguous';
                return ($base1, $base2, $base_type);
        }

        # Homozygous ref-calls
        if($$VCF_struct{$GT_id} eq 0) { 
                $base1 = $$VCF_struct{'reference_VCF_format'};
                $base_type = 'reference';
                return ($base1, $base2, $base_type);
        }

        # QC that GT matches a base for homozygous variants
        my @bases = split /,/, $$VCF_struct{'consensus_VCF_format'};
        if(($$VCF_struct{$GT_id} !~ m/(\d)([\/\|])(\d)/) && (!defined $bases[($$VCF_struct{$GT_id} - 1)])) {
                warn "Nothing found for this VCF entry:\n";
                warn Dumper($VCF_struct);
                $base1 = 'N';
                $base_type = 'ambiguous';
                return ($base1, $base2, $base_type);
        }

        # Not heterozygous
        if(($$VCF_struct{$GT_id} ne 0) && ($$VCF_struct{$GT_id} !~ m/(\d)([\/\|])(\d)/)) {
                my $consensus = $bases[($$VCF_struct{$GT_id} - 1)]; # won't be defined for heterozygous positions
        
                # Homozygous SNP
                if(length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) { 

                        # A SNP
                        if((length($$VCF_struct{'reference_VCF_format'}) eq 1) && (length($consensus) eq 1)) { 
                                $base1 = $consensus;
                                $base_type = 'snp';
                                return ($base1, $base2, $base_type);
                        }
        
                        # SNP(s) disguised as an indel
                        if((length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) && ($consensus !~ m/\./)) {
                                my @bases_reference = split //, $$VCF_struct{'reference_VCF_format'};
                                my @bases_consensus = split //, $consensus;
                                my $snp_count = 0;
                                my ($ref_base_saved, $cons_base_saved);
                                for(my $i=0; $i<scalar(@bases_reference); $i++) {
                                        my $ref_base = $bases_reference[$i];
                                        my $cons_base = $bases_consensus[$i];
                                        if($ref_base ne $cons_base) { 
                                                $ref_base_saved = $ref_base;
                                                $cons_base_saved = $cons_base;
                                                $snp_count++; 
                                        }
                                }
                                if($snp_count eq 0) {
                                        $base1 = $$VCF_struct{'reference_VCF_format'};;
                                        $base_type = 'reference';
                                        return ($base1, $base2, $base_type);
                                }
                                elsif($snp_count eq 1) {
                                        $base1 = $ref_base_saved;
                                        $base2 = $cons_base_saved;
                                        $base_type = 'snp';
                                        return ($base1, $base2, $base_type);
                                }
                                if($snp_count > 1) {
                                        $base1 = $consensus;
                                        $base_type = ('snp_multi' . $snp_count);
                                        return ($base1, $base2, $base_type);
                                }
                        }
        
                        # Ambiguous
                        warn "Nothing found for this apparant homozygous snp:\n";
                        warn Dumper($VCF_struct);
                        $base1 = 'N';
                        $base_type = 'ambiguous';
                        return ($base1, $base2, $base_type);
                }

                # Homozygous indel
                if(length($$VCF_struct{'reference_VCF_format'}) ne length($consensus)) {

                        # Deletion (maybe with snps in there too!)
                        if(length($$VCF_struct{'reference_VCF_format'}) > length($consensus)) { 
                                $base1 = $consensus;
                                $base_type = 'deletion';
                                return ($base1, $base2, $base_type);
                        }
                        if((length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) && ($consensus =~ m/^\./)) { 
                                $base1 = $consensus;
                                $base_type = 'deletion';
                                return ($base1, $base2, $base_type);
                        }       
        
                        # Insertion (maybe with snps in there too!)
                        if(length($$VCF_struct{'reference_VCF_format'}) < length($consensus)) { 
                                $base1 = $consensus;
                                $base_type = 'insertion';
                                return ($base1, $base2, $base_type);
                        }
                        
                        # Ambiguous
                        warn "Nothing found for this apparent homozygous indel:\n";
                        warn Dumper($VCF_struct);
                        $base1 = 'N';
                        $base_type = 'ambiguous';
                        return ($base1, $base2, $base_type);
                }
        }

        # Bi-allelic heterozygous positions & indels
        if($$VCF_struct{$GT_id} =~ m/(\d)([\/\|])(\d)/) { 
                $base_type = 'heterozygous';
                my @bases_het;
                if($$VCF_struct{'consensus_VCF_format'} =~ m/\,/) {
                        @bases_het = split /,/, $$VCF_struct{'consensus_VCF_format'};
                        foreach(@bases) {
                                if(length($_) > length($$VCF_struct{'reference_VCF_format'})) { $base_type = 'het_insertion'; }
                                if(length($_) < length($$VCF_struct{'reference_VCF_format'})) { $base_type = 'het_deletion'; }
                        }
                } else { 
                        push @bases_het, $$VCF_struct{'reference_VCF_format'};
                        push @bases_het, $$VCF_struct{'consensus_VCF_format'};
                        if(length($bases_het[1]) > length($bases_het[0])) { $base_type = 'het_insertion'; }
                        if(length($bases_het[1]) < length($bases_het[0])) { $base_type = 'het_deletion'; }
                }
                $base1 = $bases_het[0];
                $base2 = $bases_het[1];
                return ($base1, $base2, $base_type);
        }
        #return ($base1, $base2, $base_type);
}

1;
