package binomial_probabilities;
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

sub save_binomial_distributions {
	my $file = $_[0];
	my %probability;
	my $max_Depth = 0;

	warn "Reading from file: $file\n";
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	PROB: while(my $line=<$fh>) {
		chomp $line;
		next PROB if($line =~ m/^Depth/);
		my @line_parts = split /\t/, $line;
		my ($depth, $agree, $lower_tail_ProbHom, $lower_or_upper_tail_ProbHet, $lower_or_upper_tail_ProbHet2, $lower_or_upper_tail_ProbHet3) = @line_parts;
		$probability{'hom'}{$depth}{$agree} = $lower_tail_ProbHom;
		$probability{'het1'}{$depth}{$agree} = $lower_or_upper_tail_ProbHet;
		$probability{'het2'}{$depth}{$agree} = $lower_or_upper_tail_ProbHet2;
		$probability{'het3'}{$depth}{$agree} = $lower_or_upper_tail_ProbHet3;
		if($depth > $max_Depth) { $max_Depth = $depth; }
	}
	close $fh;
	warn "Max depth to be considered is $max_Depth. Run GBiD.pl for greater depth. Using reads up to this depth to determine genotype\n";
	return (\%probability, $max_Depth);
}

sub compare_probability {
	my ($variant_counts, $depth, $error_prob, $min_depth_pass, $refbase, $probability) = @_;
	my ($consensus, $genotype, $info) = ('',0,'.');	

	# Pass min depth
	if($min_depth_pass eq 'N') { return ($consensus, $genotype, $info); }

	# Pull info from hash
	my ($type1, $type2, $type3) = ('', '', '');
	my ($num1, $num2, $num3) = (0, 0, 0);
	foreach my $prevelance(keys %{$variant_counts}) {
		foreach my $type(keys %{$$variant_counts{$prevelance}}) {
			my $count = $$variant_counts{$prevelance}{$type};
			if($prevelance eq 1) {
				$type1 = $type;
				$num1 = $count;
			}
			if($prevelance eq 2) {
				$type2 = $type;
				$num2 = $count;
			}
			if($prevelance eq 3) {
				$type3 = $type;
				$num3 = $count;
			}
		}
	}

	### Get probabilities 
	# hom
	my $bin_hom_agree = $$probability{'hom'}{$depth}{$num1};
	# dip. het
	my $bin_het_agree = $$probability{'het1'}{$depth}{$num1};
	my $bin_het_agree2 = $$probability{'het1'}{$depth}{$num2};
	# trip. het (33:33:33)
	my $bin_het_agree3 = $$probability{'het2'}{$depth}{$num1};
	my $bin_het_agree4 = $$probability{'het2'}{$depth}{$num2};
	my $bin_het_agree5 = $$probability{'het2'}{$depth}{$num3};
	# trip. het (66:33)
	my $bin_het_agree6 = $$probability{'het3'}{$depth}{$num1};
	my $result = &MAX($bin_hom_agree, $bin_het_agree, $bin_het_agree2, $bin_het_agree3, $bin_het_agree4, $bin_het_agree5, $bin_het_agree6);

	### Compare probabilities
	# Homozygous variant
	if($result eq 0) { 
		if($bin_hom_agree >= $error_prob) {
			#warn "homozygous variant or refbase\n";
			$consensus = ($type1);
			if(($type1 !~ m/-|\+/) && ($type1 eq $refbase)) { $genotype = '0'; }
			else { $genotype = '1'; }
		}
	}

	# Heterozygous variant (biallelic)
	elsif($result eq 1) { 
		if($bin_het_agree >= $error_prob) {
			
			# biallelic with other allele specifying no variant
			if($bin_het_agree2 < $error_prob) { 
				if(($type1 !~ m/-|\+/) && ($type1 eq $refbase)) { $genotype = '0'; }
				else {
					$consensus = ($type1);
					$genotype = '0/1'; 
				}
			}
			# biallelic with other allele specifying variant
			elsif($bin_het_agree2 >= $error_prob) { 
				if(($type1 !~ m/-|\+/) && ($type1 eq $refbase)) { 
					$consensus = ($type2);
					$genotype = '0/1'; 
				}
				elsif(($type1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type1);
					$genotype = '0/1'; 
				}
				else {
					$consensus = ($type1 . ',' . $type2);
					$genotype = '1/2'; 
				}
			}
		}
	}

	# Heterozygous variant (triallelic) (33:33:33)
	elsif($result eq 3) { 
		if($bin_het_agree3 >= $error_prob) {
			
			# 2nd and 3rd alleles not specified
			if($bin_het_agree4 < $error_prob) {
				$consensus = ($type1);
				$info = 'TRIPROB:0/1/?;';
				$genotype = '0/1'; 
			}
			# 3rd allele specifies no variant
			elsif(($bin_het_agree4 >= $error_prob) && ($bin_het_agree5 < $error_prob)) {
				$consensus = ($type1 . ',' . $type2);
				$info = 'TRIALLELIC:0/1/2;';
				$genotype = '0/1/2'; 
			}
			# 3rd allele specifies variant
			elsif(($bin_het_agree4 >= $error_prob) && ($bin_het_agree5 >= $error_prob)) {
				if(($type1 !~ m/-|\+/) && ($type1 eq $refbase)) {
					$consensus = ($type2 . ',' . $type3);
					$info = 'TRIALLELIC:0/1/2;';
					$genotype = '0/1/2';
				}
				elsif(($type1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type1 . ',' . $type3);
					$info = 'TRIALLELIC:0/1/2;';
					$genotype = '0/1/2';
				}
				else {
					$consensus = ($type1 . ',' . $type2 . ',' . $type3);
					$info = 'TRIALLELIC:1/2/3;';
					$genotype = '1/2/3'; 
				}
			}
		}
	}

	# Heterozygous variant (biallelic) (66:33)
	elsif($result eq 6) { 
		if($bin_het_agree6 >= $error_prob) {
			# 2nd allele specifies no variant
			if($bin_het_agree4 < $error_prob) {
				$consensus = ($type1);
				$info = 'TRIPROB:0/1/1;';
				$genotype = '0/1'; 	
			}
			# 2nd allele specifies variant
			elsif($bin_het_agree4 >= $error_prob) {
				if(($type1 !~ m/-|\+/) && ($type1 eq $refbase)) {
					$consensus = ($type2);
					$info = 'TRIPROB:0/0/1;';
					$genotype = '0/1';
				}
				elsif(($type1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type1);
					$info = 'TRIPROB:0/1/1;';
					$genotype = '0/1';
				}
				else {
					if($num2 != 0) {
						$consensus = ($type1 . ',' . $type2);
						$info = 'TRIPROB:1/1/2;';
						$genotype = '1/2';
					} else {
						$consensus = ($type1);
						$info = 'TRIPROB:1/1/2;';
						$genotype = '0/1';
					}
				}
			}
		}
	}

	if((defined $genotype) && ($type1 =~ m/-|\+/)) { 
		if($info eq '.') { $info = 'INDEL'; }
		else { $info .= 'INDEL'; }
	}

	if((!defined $consensus) || ($consensus eq '')) { $consensus = $refbase; }
	return ($consensus, $genotype, $info);
}

sub MAX {
	my (@values) = @_;
	my ($position_of_number, $big_number) = (0, 0);
	for(my $i=0; $i<(scalar(@values) - 1); $i++) {
		if(($values[$i] > $values[$i+1]) && ($values[$i] > $big_number)) {  $position_of_number = $i; $big_number = $values[$i]; }
		elsif(($values[$i+1] > $values[$i]) && ($values[$i+1] > $big_number)) { $position_of_number = ($i+1); $big_number = $values[($i+1)]; }
	}
	return $position_of_number;
}

1;
