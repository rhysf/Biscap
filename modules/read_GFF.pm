package gfffile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;

### rfarrer@broadinstitute.org

sub read_GFF_lines {
	my ($GFF_line, $feature_wanted, $desc_seperator, $desc_column, $desc_replace) = @_;
	my %GFF_info;
	$GFF_info{'next'} = 0;
	if($GFF_line =~ /^#/) {
		$GFF_info{'next'} = 1; 
		$GFF_info{'header'}='Y';
		if($GFF_line =~ m/^\#\#FASTA/) { $GFF_info{'fasta'} = 1; }
		return \%GFF_info;
	}
	if($GFF_line =~ /^$/) { 
		$GFF_info{'next'} = 1; 
		return \%GFF_info;
	}
		
	# initial quality check
	my @bits = split /\t/, $GFF_line;
	if(@bits < 8) {
		warn "Bad GFF with < 8 columns: $GFF_line\n";
		$GFF_info{'next'}=1; 
		return \%GFF_info;
	}

	# Parts continued
	$GFF_info{'contig'}  = $bits[0];
	$GFF_info{'source'}  = $bits[1];
	$GFF_info{'feature'} = $bits[2];
	$GFF_info{'start'}   = $bits[3];
	$GFF_info{'stop'}    = $bits[4];
	$GFF_info{'score'}   = $bits[5];
	$GFF_info{'strand'}  = $bits[6];
	$GFF_info{'frame'}   = $bits[7];
	$GFF_info{'description'} = $bits[8];
	$GFF_info{'CDS'} = "$GFF_info{'start'}-$GFF_info{'stop'}";
	
	# Find feature and parent name
	if(($GFF_info{'feature'} ne $feature_wanted) && ($feature_wanted ne 'all')) { $GFF_info{'feature_wanted'} = 'N'; }
	else {
		$GFF_info{'feature_wanted'} = 'Y';
		my @description_parts = split /$desc_seperator/, $GFF_info{'description'};
		$GFF_info{'feature_parent'} = $description_parts[$desc_column];
		die "read_GFF_lines: Feature parent could not be defined from part $desc_column of $GFF_info{'description'} split by $desc_seperator\n" if(!defined $GFF_info{'feature_parent'});
		$GFF_info{'feature_parent'} =~ s/$desc_replace//g;
	}

	return \%GFF_info;
}

sub gff_to_struct {
	my ($input, $feature_wanted, $desc_seperator, $desc_column, $desc_replace) = @_;

	my %gff_struct;

	# save GFF to memory
	my (%hash_info, %hash_strand);
	warn "gff_to_struct: saving $feature_wanted from $input (split col $desc_column by $desc_seperator and remove $desc_replace)...\n";
	open my $fh, '<', $input or die "Cannot open $input: $!\n";
	GFF: while(my $line = <$fh>) {
		chomp $line;
		my $GFF_info = &read_GFF_lines($line, $feature_wanted, $desc_seperator, $desc_column, $desc_replace);
		next GFF if($$GFF_info{'next'} eq 1);
		next GFF if($$GFF_info{'feature_wanted'} eq 'N');
		if(defined $$GFF_info{'fasta'}) {
			warn "FASTA found. Ending subroutine.\n";
			last;
		}
		my ($contig, $feature_parent, $cds, $strand) = ($$GFF_info{'contig'}, $$GFF_info{'feature_parent'}, $$GFF_info{'CDS'}, $$GFF_info{'strand'});

		# save contig -> gene = 1
		$gff_struct{'contig_to_gene'}{$contig}{$feature_parent} = 1;

		# save gene = contig
		if(defined $gff_struct{'gene_to_contig'}{$feature_parent}) {
			if($gff_struct{'gene_to_contig'}{$feature_parent} ne $contig) { die "gff_to_struct: ERROR: $feature_parent found on both $gff_struct{'gene_to_contig'}{$feature_parent} and $contig\n"; }
		}
		$gff_struct{'gene_to_contig'}{$feature_parent} = $contig;

		# save gene = strand
		if(defined $gff_struct{'gene_to_strand'}{$feature_parent}) {
			if($gff_struct{'gene_to_strand'}{$feature_parent} ne $strand) { 
				warn "gff_to_struct: WARNING: $feature_parent found on both $gff_struct{'gene_to_strand'}{$feature_parent} and $strand strands\n"; 
			}
		}
		$gff_struct{'gene_to_strand'}{$feature_parent} = $strand;

		# Save exons (avoid leaving a space at the start or end)
		if(!defined $gff_struct{'gene_to_cds'}{$feature_parent}) { $gff_struct{'gene_to_cds'}{$feature_parent} = $cds; }
		else { $gff_struct{'gene_to_cds'}{$feature_parent} .= " $cds"; }
	}
	close $fh;

	# check for and remove duplicates
	# and save range 
	warn "gff_to_struct: sort exons, and check and remove any duplicate exons...\n";
	my $count_changed = 0;
	my $count_not_changed = 0;
	foreach my $gene(sort keys(%{$gff_struct{'gene_to_cds'}})) {
		my $cds = $gff_struct{'gene_to_cds'}{$gene};

		# find and replace
		my $new_cds = &sort_and_check_exons_for_duplicates($cds);
		if($cds ne $new_cds) { 
			$count_changed++; 

			# replace
			$gff_struct{'gene_to_cds'}{$gene} = $new_cds;
		}
		else { $count_not_changed++; }

		# cds range
		my @cds_parts = split / /, $new_cds;
		$gff_struct{'gene_to_cds_range'}{$gene}{'start'} = $cds_parts[0];
		$gff_struct{'gene_to_cds_range'}{$gene}{'stop'} = $cds_parts[-1];
		$gff_struct{'gene_to_cds_range'}{$gene}{'range'} = ($cds_parts[0] . '-' . $cds_parts[-1]);
	}
	my $total_seen = ($count_changed + $count_not_changed);

	# warning
	if($count_changed > 0) { warn "gff_to_struct: $count_changed / $total_seen changed (order or duplicate exons)\n"; }
	else { warn "gff_to_struct: $total_seen entires saved\n"; }

	return (\%gff_struct);
}

sub CDS_to_CDS_length {
	my $cds = $_[0];
	my $length = 0;
	my @cds_parts = split / /, $cds;
	foreach my $exon(@cds_parts) {
		my @start_stop = split /-/, $exon;
		$length += ($start_stop[1] - $start_stop[0]);
	}
	return $length;
}

sub extract_gene_from_coords_new_negative {
	my ($seq, $exon_string, $strand) = @_;
	my ($extracted_gene);

	my @exons = split /\s/, $exon_string;
	if($strand eq '-') {
		# Alphanumeric for [\d+-\d+\s]+
		#my @exons_negative = sort { $b <=> $a || $b cmp $a } @exons;
		my @exons_negative = sort { $b cmp $a } @exons;
		@exons = @exons_negative;
	}
	foreach my $exon(@exons) {
		my @terminals = split /-/, $exon;
		for(my $i=0; $i< scalar(@terminals); $i+=2) { 
			# +1 to turn array position to genome position
			my $exon_sequence = substr $seq, ($terminals[$i] - 1), ($terminals[($i + 1)] - $terminals[$i] + 1);
			if($strand eq '-') { $exon_sequence = fastafile::reverse_compliment($exon_sequence); }
			$extracted_gene .= $exon_sequence;
		}
	}

	# New. Uppercase
	$extracted_gene = uc $extracted_gene;

	# Reverse compliment for negative strand
	#if($strand eq '-') { $extracted_gene = fastafile::reverse_compliment($extracted_gene); }
	return $extracted_gene;
}

# local sub routines

sub sort_and_check_exons_for_duplicates {
	my $exons = $_[0];

	# Check none of the exons are completely repeated
	my @unsorted_exons = split /\s/, $exons;
	my (%seen);
	my @unique = grep { ! $seen{ $_ }++ } @unsorted_exons;

	# Numerically sort
	my @terminals;
	foreach my $CDS(@unique) {
		my @temp_terminals = split /-|\s/, $CDS;
		push (@terminals, @temp_terminals);
	}
	my @sorted_terminals = sort { $a <=> $b } @terminals;

	# Join
	my $newline = '';	
	for(my $i=0; $i< scalar(@terminals); $i+=2) {
		$newline .= ($terminals[$i] . '-' . $terminals[($i + 1)] . ' '); 
	}
	$newline =~ s/\s$//;

	# Return
	return $newline;
}

1;