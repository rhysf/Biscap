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
use lib "$Bin";

### rfarrer@broadinstitute.org

sub init_fasta_struct {
	my $fasta_struct = $_[0];
	die "init_fasta_struct: fasta_struct format not recognised\n" if(!defined $$fasta_struct{'filename'});
	warn "init_fasta_struct...\n";
	my %new_struct;	
	$new_struct{'filename'} = $$fasta_struct{'filename'};
	$new_struct{'total_length'} = $$fasta_struct{'total_length'};
	return \%new_struct
}

sub fasta_to_struct {
	my $input = $_[0];
	my %struct;
	$struct{'filename'} = $input;
	warn "fasta_to_struct: saving from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		my $length = length($seq);
		
		# Save
		$struct{'seq'}{$id} = $seq;
		$struct{'desc'}{$id} = $desc;
		$struct{'seq_length'}{$id} = $length;
		$struct{'total_length'} += $length;
		push @{$struct{'order'}}, $id;
	}
	return \%struct;
}

sub add_homologous_genome {
	my ($struct, $append_word) = @_;
	warn "add_homologous_genome: $append_word";

	FASTA: foreach my $id(@{$$struct{'order'}}) {
		my $new_id = ($id . $append_word);
		$$struct{'seq'}{$new_id} = $$struct{'seq'}{$id};
		$$struct{'desc'}{$new_id} = $$struct{'desc'}{$id};
		$$struct{'seq_length'}{$new_id} = $$struct{'seq_length'}{$id};
		$$struct{'total_length'} += $$struct{'seq_length'}{$id};
		push @{$$struct{'order'}}, $new_id;
	}
	return $struct;
}

sub fasta_struct_print_simple_outfile {
	my ($fasta_struct, $outfile) = @_;
	warn "fasta_struct_print_simple_outfile...\n";
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!\n";
	FASTA: foreach my $id(@{$$fasta_struct{'order'}}) {
		my $seq = $$fasta_struct{'seq'}{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print $ofh ">$id\n$seq\n";
	}
	close $ofh;
	return;
}

1;