#!/usr/bin/perl -w
use strict;
use lib "/home/unix/rfarrer/perl5/lib/perl5/";
use Statistics::R;
use Getopt::Std;
use FindBin qw($Bin);
#use lib "$Bin/../probabilities/";

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Program: Generate binomial distributions in R for SNP calling (GBiD)
Usage: perl $0 -x <Maximum depth to be calculated> -p <Probability of error e.g. 0.1 or 0.01>\n";
our($opt_x, $opt_p);
getopt('xp');
die $usage unless ($opt_x && $opt_p);    
die "-x and -p need to be numerical\n" if ((($opt_x) !~ /\d+/) || (($opt_p) !~ /\d+|\d/));
die "-x needs to be greater than 1\n" if ($opt_x < 1);
die "-p needs to be between 0 and 1\n" if (($opt_p > 1) || ($opt_p < 0));
my $prob_homozygous = (1 - $opt_p);
my $prob_heterozygous = (($prob_homozygous / 2) + ($opt_p * 0.5));
my $prob_heterozygous2 = (($prob_homozygous / 3) + ($opt_p * 0.25));
my $prob_heterozygous3 = ((($prob_homozygous / 3) * 2) + ($opt_p * 0.5));
warn "P(homozygous) = $prob_homozygous\n";
warn "P(heterozygous diploid) = $prob_heterozygous\n";
warn "P(heterozygous triploid) = $prob_heterozygous2 or $prob_heterozygous3\n";

# Input/Output
my $binomial_file = "$Bin/../probabilities/Binomial_Probabilities_for_$opt_p.tab";
my $cumulative_binomial_file = "$Bin/../Cumulative_Binomial_Probabilities_for_$opt_p.tab";

# 1) If binomial file already exists, save probabilities
my ($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3);
my $found_max_depth = 0;
if(-e $binomial_file) { ($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3, $found_max_depth) = &save_probabilities($binomial_file); }

# Print headers
else {
	open my $ofh1, '>', $binomial_file or die "Cannot open $binomial_file : $!\n";
	print $ofh1 "Depth\tAgree\tProb.Hom.\tProb.Het.Dip\tProb.Het.Trip1\tProb.Het.Trip2\n";
	close $ofh1;
}

# Calculate and print new binomial probabilities
if($opt_x <= $found_max_depth) { warn "Already calculated these probabilities\n"; }
else {
	# Calculate
	($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3) = &calculate_probabilities(($found_max_depth + 1), $opt_x, $prob_homozygous, $prob_heterozygous, $prob_heterozygous2, $prob_heterozygous3);

	# Print
	&print_binomial_probabilities($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3, $binomial_file, ($found_max_depth + 1), $opt_x);
}

# 2) If cumulative binomial file already exists, save max depth
my $found_max_depth2 = 0;
if (-e $cumulative_binomial_file) { $found_max_depth2 = &save_cumulative_probabilities($cumulative_binomial_file); }

# Print headers
else {
	open my $ofh2, '>', $cumulative_binomial_file or die "Cannot open $cumulative_binomial_file :$!\n";
	print $ofh2 "Depth\tAgree\tLower tail ProbHom\tLower or upper tail ProbHet\tLower or upper tail ProbHet2\tLower or upper tail ProbHet3\n";
	close $ofh2;
}

# Print new cumulative probabilities
if($opt_x <= $found_max_depth2) { die "Already calculated these probabilities. Done\n"; }
else {
	warn "Print cumulative probabilities between $found_max_depth2 and $opt_x\n";
	open my $ofh2, '>>', $cumulative_binomial_file or die "Cannot open $cumulative_binomial_file :$!\n";
	for(my $agree=0; $agree<=$opt_x; $agree++) {
		for(my $depth=$found_max_depth2; $depth<=$opt_x; $depth++) {
			if($agree > $depth) { $depth=$agree; }

			my $lower_hom_prob_hash1 = &tally_tail_probabilities($depth,$agree,$hom_prob_hash, 'l');
			my $lower_het_prob_hash1 = &tally_tail_probabilities($depth,$agree,$het_prob_hash1, 'l');
			my $upper_het_prob_hash1 = &tally_tail_probabilities($depth,$agree,$het_prob_hash1, 'u');
			my $lower_het_prob_hash2 = &tally_tail_probabilities($depth,$agree,$het_prob_hash2, 'l');
			my $upper_het_prob_hash2 = &tally_tail_probabilities($depth,$agree,$het_prob_hash2, 'u');
			my $lower_het_prob_hash3 = &tally_tail_probabilities($depth,$agree,$het_prob_hash3, 'l');
			my $upper_het_prob_hash3 = &tally_tail_probabilities($depth,$agree,$het_prob_hash3, 'u');

			## Expand
			my $lower_hom_prob_hash1e = &expand($lower_hom_prob_hash1);
			my $lower_het_prob_hash1e = &expand($lower_het_prob_hash1);
			my $upper_het_prob_hash1e = &expand($upper_het_prob_hash1);
			my $lower_het_prob_hash2e = &expand($lower_het_prob_hash2);
			my $upper_het_prob_hash2e = &expand($upper_het_prob_hash2);
			my $lower_het_prob_hash3e = &expand($lower_het_prob_hash3);
			my $upper_het_prob_hash3e = &expand($upper_het_prob_hash3);

			## Get the smallest probability for lower or upper tail
			my $prob_of_het1 = $upper_het_prob_hash1e; 
			my $prob_of_het2 = $upper_het_prob_hash2e; 
			my $prob_of_het3 = $upper_het_prob_hash3e; 
			if($lower_het_prob_hash1e < $upper_het_prob_hash1e) { $prob_of_het1 = $lower_het_prob_hash1e; }
			if($lower_het_prob_hash2e < $upper_het_prob_hash2e) { $prob_of_het2 = $lower_het_prob_hash2e; }
			if($lower_het_prob_hash3e < $upper_het_prob_hash3e) { $prob_of_het3 = $lower_het_prob_hash3e; }

			# Round to 10dp otherwise look-up tables get very large!
			my $rounded_homozygous = sprintf("%.10f", $lower_hom_prob_hash1e);
			my $rounded_heterozygous1 = sprintf("%.10f", $prob_of_het1);
			my $rounded_heterozygous2 = sprintf("%.10f", $prob_of_het2);
			my $rounded_heterozygous3 = sprintf("%.10f", $prob_of_het3);
			if($rounded_homozygous >= 1) { $rounded_homozygous = 1; }
			if($rounded_homozygous <= 0) { $rounded_homozygous = 0; }
			if($rounded_heterozygous1 >= 1) { $rounded_heterozygous1 = 1; }
			if($rounded_heterozygous1 <= 0) { $rounded_heterozygous1 = 0; }
			if($rounded_heterozygous2 >= 1) { $rounded_heterozygous2 = 1; }
			if($rounded_heterozygous2 <= 0) { $rounded_heterozygous2 = 0; }
			if($rounded_heterozygous3 >= 1) { $rounded_heterozygous3 = 1; }
			if($rounded_heterozygous3 <= 0) { $rounded_heterozygous3 = 0; }
		
			print $ofh2 "$depth\t$agree\t$rounded_homozygous\t$rounded_heterozygous1\t$rounded_heterozygous2\t$rounded_heterozygous3\n";
		}
	}
	close $ofh2;
}

sub tally_tail_probabilities {
	my ($depth, $agree, $homs_or_hets, $tail) = @_;
	
	my $lower_or_upper_tail = 0;
	if($tail eq 'l') {
		for(my $i=$agree; $i>=0; $i--) { 
			$lower_or_upper_tail += $$homs_or_hets{$depth}{$i};
		}
	}
	elsif($tail eq 'u') {
		for(my $i=$agree; $i<=$depth; $i++) { 
			$lower_or_upper_tail += $$homs_or_hets{$depth}{$i};
		}		
	}
	return $lower_or_upper_tail;
}


sub print_binomial_probabilities {
	my ($hom_prob_hash, $het_prob_hash, $het_prob_hash2, $het_prob_hash3, $OUTFILE, $min, $max) = @_;
	open OUT1, ">>$OUTFILE";
	for(my $agree=0; $agree<=$max; $agree++) {
		for(my $depth=$min; $depth<=$max; $depth++) {
			if($agree > $depth) { $depth=$agree; }
			my $prob_hom = $$hom_prob_hash{$depth}{$agree};
			my $prob_het = $$het_prob_hash{$depth}{$agree};
			my $prob_het2 = $$het_prob_hash2{$depth}{$agree};
			my $prob_het3 = $$het_prob_hash3{$depth}{$agree};
			
			if($prob_hom >= 1) { $prob_hom = 1; }
			if($prob_hom <= 0) { $prob_hom = 0; }
			if($prob_het >= 1) { $prob_het = 1; }
			if($prob_het <= 0) { $prob_het = 0; }
			if($prob_het2 >= 1) { $prob_het2 = 1; }
			if($prob_het2 <= 0) { $prob_het2 = 0; }
			if($prob_het3 >= 1) { $prob_het3 = 1; }
			if($prob_het3 <= 0) { $prob_het3 = 0; }

			print OUT1 "$depth\t$agree\t$prob_hom\t$prob_het\t$prob_het2\t$prob_het3\n";
		}
	}
	close OUT1;
	return 1;	
}

### Calculate the density binomial probability for chosen depth and probability of successes
sub calculate_probabilities {
	my ($min, $max, $prob_of_homozygous, $prob_of_heterozygous, $prob_of_heterozygous2, $prob_of_heterozygous3) = @_;
	my (%hom_prob_hash, %het_prob_hash, %het_prob_hash2, %het_prob_hash3);
	my $R = Statistics::R->new();
	$R->startR;
	warn "Generating binomial distributions for homozygous agree/SNP ($prob_of_homozygous), heterozygous diploid ($prob_of_heterozygous) and heterozygous triploid ($prob_of_heterozygous2 and $prob_of_heterozygous3)...\n";
	for(my $agree=0; $agree<=$max; $agree++) {
		for(my $depth=$min; $depth<=$max; $depth++) {
			if($agree > $depth) { $depth=$agree; }	
		
			# Homozygous agree/SNP
			$R->send(qq`z = dbinom($agree, $depth, $prob_of_homozygous) \n print(z)`);
			my $return1 = $R->read;
			my $homozygous = &expand($return1);		
		
			# Heterozygous SNP Diploid
			$R->send(qq`y = dbinom($agree, $depth, $prob_of_heterozygous) \n print(y)`);
			my $return2 = $R->read;
			my $heterozygous = &expand($return2);
			
			# Heterozygous SNP Triploid1
			$R->send(qq`a = dbinom($agree, $depth, $prob_of_heterozygous2) \n print(a)`);
			my $return3 = $R->read;
			my $heterozygous2 = &expand($return3);
			
			# Heterozygous SNP Triploid2
			$R->send(qq`b = dbinom($agree, $depth, $prob_of_heterozygous3) \n print(b)`);
			my $return4 = $R->read;
			my $heterozygous3 = &expand($return4);
		
			# Save values to calculate cumulative probabilities later
			# and round to 10dp otherwise look-up tables get very large!
			my $rounded_homozygous = sprintf("%.10f", $homozygous);
			my $rounded_heterozygous = sprintf("%.10f", $heterozygous);
			my $rounded_heterozygous2 = sprintf("%.10f", $heterozygous2);
			my $rounded_heterozygous3 = sprintf("%.10f", $heterozygous3);
			
			$hom_prob_hash{$depth}{$agree} = $rounded_homozygous;
			$het_prob_hash{$depth}{$agree} = $rounded_heterozygous;
			$het_prob_hash2{$depth}{$agree} = $rounded_heterozygous2;
			$het_prob_hash3{$depth}{$agree} = $rounded_heterozygous3;
		}
		warn "$agree / $max\n";
	}
	$R->stopR();	
	return(\%hom_prob_hash, \%het_prob_hash, \%het_prob_hash2, \%het_prob_hash3);
}

sub expand {
        my $n = shift;
        $n =~ s/\[1\]|\s//g;
        return $n unless $n =~ /^(.*)e([-+]?)(.*)$/;
        my ($num, $sign, $exp) = ($1, $2, $3);
        my $sig = $sign eq '-' ? "." . ($exp - 1 + length $num) : '';
        return sprintf "%${sig}f", $n;
}

sub save_probabilities {
	my $file = $_[0];
	my ($found_max_depth) = 0;
	my ($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3);
	warn "save_probabilities: $file\n";
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @line_parts = split /\t/, $line;
		my ($depth, $agree, $ProbHom, $ProbHet1, $ProbHet2, $ProbHet3) = @line_parts;

		next if($line =~ m/^Depth/);
		$$hom_prob_hash{$depth}{$agree} = $ProbHom;
		$$het_prob_hash1{$depth}{$agree} = $ProbHet1;
		$$het_prob_hash2{$depth}{$agree} = $ProbHet2;
		$$het_prob_hash3{$depth}{$agree} = $ProbHet3;
		if($depth > $found_max_depth) { $found_max_depth = $depth; }
	}
	close $fh;
	warn "save_probabilities: max depth = $found_max_depth\n";
	return($hom_prob_hash, $het_prob_hash1, $het_prob_hash2, $het_prob_hash3, $found_max_depth);
}

sub save_cumulative_probabilities {
	my $file = $_[0];
	my ($found_max_depth) = 0;
	warn "save_cumulative_probabilities: $file\n";
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @line_parts = split /\t/, $line;
		my ($depth, $agree, $lower_tail_ProbHom, $lower_or_upper_tail_ProbHet, $lower_or_upper_tail_ProbHet2, $lower_or_upper_tail_ProbHet3) = @line_parts;

		next if($line =~ m/^Depth/);
		if($depth > $found_max_depth) { $found_max_depth = $depth; }
	}
	close $fh;
	warn "save_cumulative_probabilities: max depth = $found_max_depth\n";
	return ($found_max_depth);
}
