package pileupfile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use ascii;
my %ascii = ascii::ascii();

### rfarrer@broadinstitute.org

# Find the ^ and $ in pileup aligned column
sub save_and_remove_mapping_quality {
	my ($aligned) = @_;
	my (@mapping_qualities, @positions_of_ends); 
	
	# Start of read segments are given along with Mapping Quality scores 
	# as ASCII values -33 if -s option is not applied.
	while($aligned =~ m/\^./ig) { 
		my @quality_parts = split //, $&;
		my $mapping_quality = ($ascii{$quality_parts[1]} -33);
		push (@mapping_qualities, $mapping_quality);
	}
	$aligned =~ s/\^.//g;

	# End of reads
	my @aligned_bits = split //, $aligned;
	my $position_in_line = 0;
	for(my $i=0; $i<scalar(@aligned_bits); $i++) {
		if($aligned_bits[$i] eq '$') {
			$position_in_line--;
			push (@positions_of_ends, $position_in_line); 
		}
		$position_in_line++;
	}
	$aligned =~ s/\$//ig;

	# Return the remaining bases, and details of the start- and ends- of reads
	return ($aligned, \@positions_of_ends, \@mapping_qualities);
}

sub filter_aligned_on_read_quality {
	my ($aligned, $qual, $min_qual) = @_;
	my ($new_aligned, $new_qual);
	my @qual_parts = split //, $qual;
	my @read_qual_parts;

	my ($offset, $removed_bases) = (0,0);
	for(my $i=0; $i<scalar(@qual_parts); $i++) {
		my $quality_value = ($ascii{$qual_parts[$i]});
		if($quality_value ne 0) {
			push (@read_qual_parts, $quality_value);
		}
		if($quality_value < $min_qual) { 
			# Remove quality score
			$new_qual = substr $qual, 0, ($i - $offset);
			$new_qual .= substr $qual, (($i+1)  - $offset), length($qual);
			$qual = $new_qual;
			
			# Remove corresponding read
			$new_aligned = substr $aligned, 0, ($i - $offset);
			$new_aligned .= substr $aligned, (($i+1)  - $offset), length($aligned);
			$aligned = $new_aligned;
			
			$offset++;
		}
	}
	return ($aligned, $removed_bases, \@read_qual_parts);
}

sub count_bases {
	my ($ref_base, $aligned) = @_;
	my (%bases);

	# Count insertions (deletions are for subsequent positions)
	while($aligned =~ m/([\+|\-])(\d+)([ACGTNacgtn]+)/ig) {
	#while($aligned =~ m/([\.\,ACTGNactgn])([\+|\-])(\d+)([ACGTNacgtn]+)/ig) {
		my $length_of_indel = $2;
		my $indel = substr $3, 0, $2;
		$indel =~ tr/a-z/A-Z/;
		warn "indel = $indel\n";
		#if($1 =~ m/\-/) { $bases{('-' . $indel)}++; }
		if($1 !~ m/\-/) { $bases{$indel}++; }
	}

	# Remove insertions and deletions (count '*' when they appear) to stop it being confused with mismatches
	while($aligned =~ m/[\-\+](\d+)([ACGTNacgtn]+)/ig) {
		my $length_of_insertion = $1;
		my $newaligned = substr $aligned, 0, length($`);
		my $offset = (length($`) + length($length_of_insertion) + $length_of_insertion + 1);
		$newaligned .= substr $aligned, $offset, (length($aligned) - $offset);
		$aligned = $newaligned;
	}

	# Start of read segments are given along with Quality scores as ASCII values -33 if -s option is not applied.
	while($aligned =~ m/\^./ig) { 
		my @quality_parts = split //, $&;
		my $quality = ($ascii{$quality_parts[1]} -33);
		$aligned =~ s/\^.//g;
	}
	my $end_of_read_segment = ($aligned =~ tr/$//);

	# Number of matches and mismatches
	$bases{$ref_base} = ($aligned =~ tr/\.|\,//);
	$bases{'A'} += ($aligned =~ tr/A|a//);
	$bases{'C'} += ($aligned =~ tr/C|c//);
	$bases{'T'} += ($aligned =~ tr/T|t//);
	$bases{'G'} += ($aligned =~ tr/G|g//);
	$bases{'N'} += ($aligned =~ tr/N|n//);
	#$bases{'-' . $ref_base} += ($aligned =~ tr/\*//);
	#foreach my $insertion(keys %insertions) { $bases{$insertion} = $insertions{$insertion}; }
	foreach my $base(keys %bases) {
		delete $bases{$base} if($bases{$base} eq 0);
	}
	return \%bases;
}

sub find_base {
	my ($ref_base, $depth, $aligned, $mindepth, $percent_consensus, $min_het_percent_consensus, $max_het_percent_consensus) = @_;
	my ($nucleotide, $sample_info);
	
	# Depths required for variant-calling. Return if too low
	if($depth < $mindepth) { return ($ref_base, 'NA'); }
	my $expected_depth_for_snp = ($depth * $percent_consensus);
	my $expected_min_depth_for_het = ($depth * $min_het_percent_consensus);
	my $expected_max_depth_for_het = ($depth * $max_het_percent_consensus);
	
	# Count bases and determine consensus based on expected depths
	my $base_counts = &count_bases($ref_base, $aligned);
	foreach my $base(keys %{$base_counts}) {
		my $tally = $$base_counts{$base};
		# Homozygous SNP or indel
		if($tally >= $expected_depth_for_snp) { $nucleotide = $base; }
		# Heterozygous position or indel
		elsif(($tally >= $expected_min_depth_for_het) && ($tally <= $expected_max_depth_for_het)) { $nucleotide .= ($base . ","); }
	}
	$nucleotide =~ s/,$//;

	# Identify genotype
	if(!defined $nucleotide) { 
		$nucleotide = '?';
		$sample_info = 'NA'; 
	}
	elsif($nucleotide eq $ref_base) { $sample_info = ('0:' . $depth); }
	elsif($nucleotide !~ m/,/) { $sample_info = ('1:' . $depth); }
	else {
		my @nt_parts = split /,/, $nucleotide;

		# if neither ref base: 1/2
		if(($nt_parts[0] ne $ref_base) && ($nt_parts[0] ne $ref_base)) { $sample_info = ('1/2:' . $depth); }
		# if one ref base
		elsif($nt_parts[0] eq $ref_base) {
			$nucleotide = $nt_parts[1];
			$sample_info = ('0/1:' . $depth); 
		}
		elsif($nt_parts[1] eq $ref_base) {
			$nucleotide = $nt_parts[0];
			$sample_info = ('0/1:' . $depth); 
		}
	}
	return ($nucleotide, $sample_info);
}

sub split_pileup_into_contigs {
	my ($pileup, $directory) = @_;

	# make directory
	my $mkdir = "mkdir $directory";
	system($mkdir);

	# Save the pileup into smaller parts
	my $lastcontig;
	warn "Writing data from $pileup into smaller parts...\n";
	open IN, "<$pileup" or die "Can't open $pileup : $!";
	while(my $line=<IN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($contig, $pos, $ref_base, $depth, $aligned, $qual) = (@bits);
		if(!defined $lastcontig) { 
			$lastcontig = $contig;
			my $outfile = ($directory . $contig);
			open OUTPUT, ">$outfile";
		}
		elsif($lastcontig ne $contig) {
			close OUTPUT;
			$lastcontig = $contig;
			my $outfile = ($directory . $contig);
			open OUTPUT, ">$outfile";
		}
		print OUTPUT "$line\n";
	}
	return 1;
}

sub remove_split_pileup {
	my $dir = $_[0];
	my $cmd = "rm -r $dir";
	system($cmd);
	return 1;
}

sub save_pileup_to_hash {
	my $pileup = $_[0];
	my (%pileupdata);
	if(! -e $pileup) {
		warn "Cannot open $pileup - either no coverage or the coords have IDs not matching the pileup's\n"; 
		return \%pileupdata;
	}
	
	open IN, "<$pileup";
	warn "Saving data from $pileup...\n";
	while(my $line=<IN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($contig, $pos, $ref_base, $depth, $aligned, $qual) = (@bits);
		$pileupdata{$contig}{$pos}{$ref_base}{$depth}{$aligned} = $qual;
	}
	close IN;
	return \%pileupdata;
}

sub query_pileup_from_nucleotide {
	my ($pileup_hash) = @_;
	my ($ref_base, $depth_of_coverage, $aligned_bases, $qual);	
	foreach my $ref(keys %{$pileup_hash}) { 
		$ref_base = $ref;
		foreach my $depth(keys %{$$pileup_hash{$ref}}) { 
			$depth_of_coverage = $depth;
			foreach my $aligned(keys %{$$pileup_hash{$ref}{$depth}}) {
				$aligned_bases = $aligned;
				$qual = $$pileup_hash{$ref}{$depth}{$aligned};
			}
		}
	}
	return ($ref_base, $depth_of_coverage, $aligned_bases, $qual);
}

1;
