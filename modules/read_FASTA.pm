package fastafile;
use strict;
use Bio::SeqIO;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use color_space;

### rfarrer@broadinstitute.org

sub fasta_id_to_seq_hash {
	my $input = $_[0];
	my (%sequences, %descriptions);
	my @order;
	warn "fasta_id_to_seq_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		$sequences{$id}=$seq;
		$descriptions{$id}=$desc;
		push @order, $id;
	}
	return (\%sequences, \%descriptions, \@order);
}

sub fasta_id_to_seq_hash_1_entry {
	my ($input, $want) = @_;
	my (%sequences);
	warn "fasta_id_to_seq_hash_1_entry: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	FASTA: while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		next FASTA if ($want ne $id);
		$sequences{$id}=$seq;
	}
	return (\%sequences);
}

sub fasta_id_to_order_array {
	my $input = $_[0];
	my @order;
	warn "fasta_id_to_order_array: saving order from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		push @order, $id;
	}
	return (\@order);
}

sub fasta_id_to_seq_length_hash {
	my $input = $_[0];
	my (%lengths);
	warn "fasta_id_to_seq_length_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		$lengths{$id} = $length;
	}
	return (\%lengths);
}

sub fasta_id_to_seq_length_array {
	my $input = $_[0];
	my (@lengths);
	warn "fasta_id_to_seq_length_array: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		push @lengths, length($seq);	
	}
	return (\@lengths);
}

sub fasta_to_total_seq_length {
	my $input = $_[0];
	my ($lengths) = 0;
	warn "fasta_to_total_seq_length: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
    		my $seq = $seq_obj->seq;
		$lengths += length($seq);
	}
	return ($lengths);
}

sub fasta_hash_truncate {
	my ($fasta_hash, $truncate_from_5prime, $truncate_from_3prime, $minimum_length) = @_;
	my ($trimmed_count, $skipped_count) = (0, 0);
	my %sequences;
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		my $length = length($seq);

		# truncate 5' and 3'
		my $truncated_length = ($length - $truncate_from_3prime);
		my $truncated_seq = substr $seq, $truncate_from_5prime, $truncated_length;
		if($truncated_length < $minimum_length) {
			$skipped_count++;
			next FASTA;
		}
		$sequences{$id}=$truncated_seq;
		$trimmed_count += (length($seq) - length($truncated_seq));
	}
	warn "Trimmed $trimmed_count bases.\n";
	warn "Skipped $skipped_count entries for < $minimum_length\n";
	return (\%sequences);
}

sub fasta_hash_append_id {
	my ($fasta_hash, $append_id) = @_;
	warn "fasta_hash_append_id: $append_id\n";
	my %sequences;
	foreach my $id(keys %{$fasta_hash}) {
		my $newid = ($id . $append_id);
		$sequences{$newid} = $$fasta_hash{$id};
	}
	return (\%sequences);
}

sub fasta_order_append_id {
	my ($fasta_order, $append_id) = @_;
	warn "fasta_order_append_id: $append_id\n";
	my @order;
	foreach my $id(@{$fasta_order}) {
		$id .= $append_id;
		push @order, $id;
	}
	return (\@order);
}

sub fasta_hash_replace_id_with_desc {
	my ($fasta_hash, $descriptions, $replace_with_desc_part_number) = @_;
	warn "fasta_hash_replace_id_with_desc: $replace_with_desc_part_number\n";
	my %sequences;
	foreach my $id(keys %{$fasta_hash}) {
		my @desc_parts = split /\s/, $$descriptions{$id};
		my $new_id = $desc_parts[$replace_with_desc_part_number];
		$sequences{$new_id}=$$fasta_hash{$id};
	}
	return (\%sequences);
}

sub fasta_order_replace_id_with_desc {
	my ($fasta_order, $descriptions, $replace_with_desc_part_number) = @_;
	warn "fasta_order_replace_id_with_desc: $replace_with_desc_part_number\n";
	my @order;
	foreach my $id(@{$fasta_order}) {
		my @desc_parts = split /\s/, $$descriptions{$id};
		my $new_id = $desc_parts[$replace_with_desc_part_number];
		push @order, $new_id;
	}
	return (\@order);
}

sub fasta_order_remove_id {
	my ($fasta_order, $remove_id) = @_;
	warn "fasta_order_remove_id: $remove_id\n";
	my @order;
	foreach my $id(@{$fasta_order}) {
		if($id ne $remove_id) { push @order, $id; }
	}
	return (\@order);
}


sub fasta_order_remove_words_from_id {
	my ($fasta_order, $remove_words) = @_;
	warn "fasta_order_remove_words_from_id: removing $remove_words\n";
	my @words = split /,/, $remove_words;
	my @order;
	foreach my $id(@{$fasta_order}) {
		foreach(@words) {
			$id =~ s/$_//g;
		}
		push @order, $id;
	}
	return (\@order);
}

sub fasta_hash_remove_words_from_id {
	my ($fasta_hash, $remove_words) = @_;
	warn "fasta_hash_remove_words_from_id: removing $remove_words\n";
	my @words = split /,/, $remove_words;
	my %sequences;
	foreach my $id(keys %{$fasta_hash}) {
		my $old_id = $id;
		foreach(@words) {
			$id =~ s/$_//g;
		}
		$sequences{$id} = $$fasta_hash{$old_id};
	}
	return (\%sequences);
}

sub fasta_hash_reverse_compliment {
	my ($fasta_hash, $ids, $reverse) = @_;
	my %sequences = %{$$fasta_hash};
	if($reverse eq 'y') {
		foreach my $id(keys %{$fasta_hash}) {
			my @reverse_ids = split /\,/, $ids;
			foreach(@reverse_ids) {
				if($_ eq $id) {
					 $sequences{$id} = &reverse_compliment($$fasta_hash{$id});
				}
			}
		}	
	}
	return (\%sequences);
}

# Order by None/ignore (n), Biggest to smallest (b), IDs from opt_i (i), File from opt_k (k)
sub fasta_hash_order {
	my ($fasta_hash, $order_by, $ids, $additional_file) = @_;
	my @order;

	# Should not be in here (waste of time) but just in case!
	if($order_by eq 'n') {
		foreach my $id(keys %{$fasta_hash}) {
			push @order, $id;
		}
	}

	# Reorder according to length
	elsif($order_by eq 'b') { 
		my %lengths;
		# Save lengths
		foreach my $id(keys %{$fasta_hash}) {
			my $length = length($$fasta_hash{$id});
			LENGTHS: while(defined $lengths{$length}) { 
				#die "Two entries are same length and asked to order by lengths in $0\n"; 
				$length++;	
			}
			$lengths{$length} = $id;
		}
		# Reorder
		foreach my $length(sort {$b<=>$a} keys %lengths) {
			my $id = $lengths{$length};
			push @order, $id;
		}
	}

	# Reorder according to ids specificed at run time (opt_i)
	elsif($order_by eq 'i') {
		@order = split /,/, $ids;
	}

	# Reorder according to file of IDs
	elsif($order_by eq 'k') {
		open IN, "<$additional_file" or die "Cannot open $additional_file: $!\n";
		while(my $line=<IN>) {
			chomp $line;
			push @order, $line;
		}
		close IN;
	}
	else { die "Reorder according to what in $0?: $order_by\n"; }
	return \@order;
}

sub fasta_hash_print_simple {
	my $fasta_hash = $_[0];
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print ">$id\n$seq\n";
	}
	return;
}

sub fasta_hash_print_simple_outfile {
	my ($fasta_hash, $outfile) = @_;
	open OUT, ">$outfile" or die "Cannot open $outfile: $!\n";
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print OUT ">$id\n$seq\n";
	}
	close OUT;
	return;
}


sub fasta_hash_print {
	my ($fasta_hash, $original_FASTA_file, $output_format, $data_type, $order_array, $split_into_files) = @_;
	my ($sequence_length, $sequence_count, $split_count, $split_name) = (0, 0, 0, 0);
	my $output_name;

	warn "fasta_hash_print output=$output_format\n";
	FASTA: foreach my $id(@{$order_array}) {
		#next FASTA if(!defined $$fasta_hash{$id});
		die "$id is not found in fasta_hash\n" if(!defined $$fasta_hash{$id});
		my $seq = $$fasta_hash{$id};
		$sequence_length = length($seq);
		$sequence_count++;

		# Color-space FASTA
		if($output_format eq 'color') { $seq = color_space::base_to_color($seq); }

		# FASTA (base or color) output
		# Print with sequence lines up to 60 characters long
		if(($output_format eq 'fasta') || ($output_format eq 'color')) {
			$seq =~ s/(\S{60})/$1\n/g;
			my $entry = ">$id\n$seq\n";	

			# Print to standard output
			if($split_into_files eq 'n') { print "$entry"; }

			# Print to separate files
			else {
				if($split_count eq 0) {
					if($split_into_files eq 1) { 
						$output_name = ($original_FASTA_file . '-' . $id . '.fasta');
					}
					else {
					#if($split_into_files > 1) {
						$output_name = ($original_FASTA_file . '-split-into-' . $split_into_files . '-entries-per-file-pt-' . $split_name . '.fasta');
					}
					open OUTPUT, ">$output_name" or die "Cannot open $output_name: $!\n";		
				}
				print OUTPUT $entry;
				$split_count++;
				if($split_count eq $split_into_files) {
					$split_count = 0;
					close OUTPUT;
					$split_name++;
				}
			}
		}

		# Threaded-block aligner
		if($output_format eq 'tba') {
			my $new_id= ">$data_type:$id:1:+:$sequence_length";
			print "$new_id\n$seq\n";
		}
		# FASTQ
		if($output_format eq 'fastq') {
			my $quality;
			for(my $i=0; $i<$sequence_length; $i++) { $quality .= "I"; }
			print '@' . "$id\n$seq\n+\n$quality\n";
		}
	}

	# NEXUS output
	if($output_format eq 'nexus') { &print_nexus($fasta_hash, $sequence_count, $sequence_length, $data_type); }
	return 1;
}

sub fasta_summary {
	my ($input, $print_summary_of_full_file, $print_summary_for_each_entry) = @_;
	my %information;
	my @fasta_lengths;
	
	# Print headers for summary per entry
	if($print_summary_for_each_entry ne 'n') {
		# Long
		if($print_summary_for_each_entry eq 'l') { warn "Identity\tDescription\tSequence_length\tA\tC\tT\tG\tN\tGC\tGCpercent\n"; }
		# Short
		elsif($print_summary_for_each_entry eq 's') { warn "Identity\tSequence_length\n"; }
		# Unknown
		else { die "fasta_summary in read_FASTA.pm: $print_summary_for_each_entry is not l (long) or s (short)\n"; }
	}

	#warn "fasta_summary: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		if($desc eq '') { $desc = 'NA'; }
		my $length = length($seq);
		
		# Base counts per entry (not returned in hash)
		my %base_counts;
		$base_counts{'A'} = ($seq =~ tr/A|a//);
		$base_counts{'C'} = ($seq =~ tr/C|c//);
		$base_counts{'T'} = ($seq =~ tr/T|t//);
		$base_counts{'G'} = ($seq =~ tr/G|g//);
		$base_counts{'N'} = ($seq =~ tr/N|n//);
		$base_counts{'GC'} = ($seq =~ tr/G|g|C|c//);

		# Summary (returned in hash)
		$information{'description'}{$id}=$desc;
		$information{'full_length'} += $length;
		$information{'lengths'}{$id}=$length;
		$information{'n_of_seq'}++;
		$information{'full_length_minus_N'} += ($information{'lengths'}{$id} - $base_counts{'N'});
		foreach my $type(keys %base_counts) { $information{'base_counts'}{$type} += $base_counts{$type}; }
		push @fasta_lengths, $length;

		# Print summary per entry
		if($print_summary_for_each_entry eq 'l') { 
			my $GC_percent = sprintf("%.3f", (($base_counts{'GC'} / $length)*100));
			print "$id\t$desc\t$length\t$base_counts{'A'}\t$base_counts{'C'}\t$base_counts{'T'}\t$base_counts{'G'}\t$base_counts{'N'}\t$base_counts{'GC'}\t$GC_percent\n";
		}
		if($print_summary_for_each_entry eq 's') { print "$id\t$length\n"; }
	}
	my $n50_n90 = &calculate_n50_and_n90(\@fasta_lengths);
	%information = (%information, %{$n50_n90});

	# Print summary of whole file
	if($print_summary_of_full_file ne 'n') {
		warn "\nSummary of $input:\n";
		warn "Number of sequences counted:\t$information{'n_of_seq'}\n";
		warn "Length of sequences total:\t$information{'full_length'}\n";
		warn "Length of sequences minus Ns:\t$information{'full_length_minus_N'}\n";
		foreach my $type(keys %{$information{'base_counts'}}) {
			warn "Number of $type:\t$information{'base_counts'}{$type}\n";
		}
		warn "NMAX:\t$information{'NMAX'}\n";
		warn "N50:\t$information{'N50'}\n";
		if(defined $information{'N90'}) { warn "N90:\t$information{'N90'}\n"; }
	}
	return (\%information);
}

sub calculate_n50_and_n90 {
	my $lengths = $_[0];
	my %information;
	my @fasta_lengths=sort{$b<=>$a} @{$lengths}; 
	my $seq_length = 0;
	foreach(@fasta_lengths) {
		$seq_length += $_;
	}
	my ($count,$half)=(0,0);
	COUNT: for(my $i=0; $i<@fasta_lengths; $i++) {
		# Top entry
		if($count eq 0) { 
			$information{'NMAX'} = $fasta_lengths[$i]; 
		}
		$count += $fasta_lengths[$i];
		# Midway
		if(($count >= ($seq_length / 2)) && ($half == 0)) {
			$information{'N50'} = $fasta_lengths[$i];
			$half=$fasta_lengths[$i];
		}
		# 90% through
		elsif($count >= ($seq_length * 0.9)) {
			$information{'N90'} = $fasta_lengths[$i];
			last COUNT;
		}
	}
	return (\%information);
}

sub reverse_compliment { 
	my $sequence = $_[0];
	$sequence = uc ($sequence);

	# Standard
	$sequence =~ tr/[A,C,T,G]/[T,G,A,C]/;

	# IUPAC Ambiguity Codes
	$sequence =~ tr/[Y,R,K,M]/[R,Y,M,K]/;

	my @characters = split(//, $sequence);
	my $reverse_compliment = '';
	for(my $i=scalar(@characters); $i > 0; $i--) { $reverse_compliment .=  $characters[($i - 1)]; }
	return $reverse_compliment;	
}

sub print_nexus {
	my ($fasta_hash, $sequence_count, $sequence_length, $data_type) = @_;
	print "#NEXUS\n";
	print "begin data;\n";
	print "dimensions ntax=$sequence_count nchar=$sequence_length;\n";
	print "format datatype=$data_type interleave=no gap=-;\n";
	print "matrix\n";
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$id =~ s/-/_/g;
		print  "$id\n$seq\n";
	}
	print ";\n";
	print "end;\n";
}

1;
