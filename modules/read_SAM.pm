package samlines;
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

sub read_SAM_lines {
	my $SAM_line = $_[0];
	my %SAM_info;
	$SAM_info{'next'}=0;

	# Save headers
	if($SAM_line =~ m/^\@/) { 
		$SAM_info{'next'}=1; 
		$SAM_info{'header'}='Y';

		# Contig length
		if($SAM_line =~ m/^\@SQ/) {
			my @bits = split /\t/, $SAM_line;
			my ($header, $contig, $length) = @bits;
			$contig =~ s/^SN://;
			$length =~ s/^LN://;
			$SAM_info{'contigs'}{$contig} = $length;
		}
		return \%SAM_info; 
	}
	$SAM_info{'header'}='N';

	# Parts
	my @bits = split /\t/, $SAM_line;

	# initial quality check
	if(@bits < 11) {
		warn "$0: Bad SAM with < 11 columns: $SAM_line\n";
		$SAM_info{'next'}=1;
		return \%SAM_info;
	}

	# Parts continued
	$SAM_info{'Qname'}               = $bits[0];
	$SAM_info{'bitwiseFLAG'}         = $bits[1];
	$SAM_info{'Ref_seq_name'}        = $bits[2];
	$SAM_info{'leftmost_position'}   = $bits[3];
	$SAM_info{'mapping_quality'}     = $bits[4];
	$SAM_info{'CIGAR'}               = $bits[5];
	$SAM_info{'reference_name_mate'} = $bits[6];
	$SAM_info{'position_mate'}       = $bits[7];
	$SAM_info{'TLEN'}                = $bits[8];
	$SAM_info{'SEQ'}                 = $bits[9];
	$SAM_info{'QUAL'}                = $bits[10];
	$SAM_info{'Length'}              = length($SAM_info{'SEQ'});
	$SAM_info{'rightmost_position'}  = ($SAM_info{'leftmost_position'} + $SAM_info{'Length'});

	# Optional fields TAG:TYPE:VALUE (not interested in all of them!)
	# Types = A(character), B(generalarray), f(real number), H(hexadecimal array), i(integer), or Z(string)
	for(my $i=11; $i<scalar(@bits); $i++) {
		my @parts = split /:/, $bits[$i];
		my ($tag, $type, $value) = @parts;

		# SAM generic
		if(($tag eq 'AS') && ($type eq 'i')) { $SAM_info{'alignment_score'} = $value; }
		if(($tag eq 'MD') && ($type eq 'Z')) { $SAM_info{'mismatching_positions'} = $value; }

		# BWA specific fields (http://bio-bwa.sourceforge.net/bwa.shtml)
		if(($tag eq 'X0') && ($type eq 'i')) { $SAM_info{'number_of_best_hits'} = $value; }
		if(($tag eq 'X1') && ($type eq 'i')) { $SAM_info{'number_of_suboptimal_hits'} = $value; }
		
	}

	# Unaligned / how many times mapped?
	if(($SAM_info{'Ref_seq_name'} eq '*') || ($SAM_info{'bitwiseFLAG'} eq 4)) { $SAM_info{'aligned'} = 'N'; }
	else { $SAM_info{'aligned'} = 'Y'; }

	# Uniquely mapped? (either 256 or 0)
	if(($SAM_info{'bitwiseFLAG'} & 256) ne 0) { $SAM_info{'secondary'} = 'Y'; } 
	else { $SAM_info{'secondary'} = 'N'; }

	# Return
	return \%SAM_info;
}

sub SAM_to_contig_lengths {
	my $sam_file = $_[0];
	my %contigs;
	my $sum_length;
	warn "SAM_to_contig_lengths: $sam_file...\n";
	open my $fh, '<', $sam_file or die "Cannot open $sam_file: $!\n";
	READHEADERS: while (my $line = <$fh>) {
		chomp $line;
		#my $SAM_line = &read_SAM_lines($line);
		
		my @bits = split /\t/, $line;
		my ($header, $chr, $length) = @bits;
		if($header eq '@SQ') {
			$chr =~ s/^SN://;
			$length =~ s/^LN://;
			$contigs{$chr} = $length;
			$sum_length += $length;
		}
		last READHEADERS if($header !~ m/^@/); 
	}
	close $fh;
	return (\%contigs, $sum_length);
}

1;
