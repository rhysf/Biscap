package color_space;
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

sub color_to_base{
	my $seq = shift;
	
	my %colourspace;
	$colourspace{A0} = "A";
	$colourspace{C0} = "C";
	$colourspace{G0} = "G";
	$colourspace{T0} = "T";
	$colourspace{A1} = "C";
	$colourspace{C1} = "A";
	$colourspace{G2} = "A";
	$colourspace{A2} = "G";
	$colourspace{A3} = "T";
	$colourspace{T3} = "A";
	$colourspace{C2} = "T";
	$colourspace{T2} = "C";
	$colourspace{C3} = "G";
	$colourspace{G3} = "C";
	$colourspace{G1} = "T";
	$colourspace{T1} = "G";	
	
	my @letters = split //, $seq;
	my $first_base = $letters[0];
	
	for(my $i = 1; $i < (scalar(@letters) - 1); $i++) {
		my $colour = $letters[$i];
		my $encoding = $first_base.$colour;
		$first_base = $colourspace{ $encoding };
		$letters[$i] = $first_base;    
	}
	my $base_space_seq = join('',@letters);
	return $base_space_seq;
}

sub base_to_color{
	my $seq = shift;
	
	my %colourspace;
	$colourspace{A}{A} = 0;
	$colourspace{C}{C} = 0;
	$colourspace{T}{T} = 0;
	$colourspace{G}{G} = 0;
	$colourspace{A}{C} = 1;
	$colourspace{C}{A} = 1;
	$colourspace{G}{T} = 1;
	$colourspace{T}{G} = 1;
	$colourspace{A}{G} = 2;
	$colourspace{G}{A} = 2;
	$colourspace{T}{C} = 2;
	$colourspace{C}{T} = 2;	
	$colourspace{A}{T} = 3;
	$colourspace{T}{A} = 3;
	$colourspace{G}{C} = 3;
	$colourspace{C}{G} = 3;	
	
	my @letters = split //, $seq;
	my $color_space_seq = $letters[0];
	
	for(my $i = 1; $i < scalar(@letters); $i++) {
		my $first_base = $letters[($i - 1)];
		my $second_base = $letters[$i];
		my $colspace = $colourspace{$first_base}{$second_base};
		$color_space_seq .= $colspace;    
	}
	return $color_space_seq;
}

1;
