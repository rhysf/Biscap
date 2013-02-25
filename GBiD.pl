#!/usr/bin/perl -w
use strict;
use Statistics::R;
use Getopt::Std;
use Cwd;
my $dir = getcwd;

# r.farrer09@imperial.ac.uk

### Opening commands 

my $usage = "
Program:  Generate binomial distributions in R for SNP calling (GBiD)
Version:  0.11
Contact:  Rhys Farrer <r.farrer09\@imperial.ac.uk>

Usage:    $0 <commands>

Commands: -x\tMaximum depth to be calculated
          -p\tProbability of error e.g. 0.1 or 0.01
";

our($opt_x, $opt_p);
getopt('xp');
die $usage unless ($opt_x && $opt_p);    
my $opt_m = 1;
die "-x and -p need to be numerical\n" if ((($opt_x) !~ /\d+/) || (($opt_p) !~ /\d+|\d/));
die "-x needs to be greater than 1\n" if ($opt_x < 1);
die "-p needs to be between 0 and 1\n" if (($opt_p > 1) || ($opt_p < 0));
my $prob_homozygous = (1 - $opt_p);
my $prob_heterozygous = (($prob_homozygous / 2) + ($opt_p * 0.5));
my $prob_heterozygous2 = (($prob_homozygous / 3) + ($opt_p * 0.25));
my $prob_heterozygous3 = ((($prob_homozygous / 3) * 2) + ($opt_p * 0.5));

### Check if the depths have already been calculated. 
### Save what has been calculated, and calculate anything that hasn't 

my ($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3);
my $Binom_OUT = ($dir . "/Binomial_Probabilities_for_$opt_p.tab");
if (-e $Binom_OUT) {
	# Read what has been calculated to find what else is needed and for cumulative probabilities
	my ($found_max_depth) = 0;
	warn "Reading from file: $Binom_OUT\n";
	open IN1, "<$Binom_OUT";
	while(my $line=<IN1>) {
		chomp $line;
		my @line_parts = split /\t/, $line;
		my ($depth, $agree, $ProbHom, $ProbHet1, $ProbHet2, $ProbHet3) = @line_parts;
		
		next if($line =~ m/^Depth/);
		$$homozygous_probability{$depth}{$agree} = $ProbHom;
		$$heterozygous_probability1{$depth}{$agree} = $ProbHet1;
		$$heterozygous_probability2{$depth}{$agree} = $ProbHet2;
		$$heterozygous_probability3{$depth}{$agree} = $ProbHet3;
		if($depth > $found_max_depth) { $found_max_depth = $depth; }
	}
	close IN1;
	warn "Max depth found in file is $found_max_depth\n";
	if($opt_x <= $found_max_depth) { warn "Already calculated these probabilities\n"; }
	elsif($opt_m eq ($found_max_depth + 1)) {
		warn "Calculating new values...\n";
		($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3) = &calculate_probabilities($opt_m, $opt_x, $prob_homozygous, $prob_heterozygous, $prob_heterozygous2, $prob_heterozygous3);
		&print_binomial_probabilities($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3, $Binom_OUT, $opt_m, $opt_x, 'F');
	}
	elsif(($opt_m > ($found_max_depth + 1)) || ($opt_m < ($found_max_depth + 1))) {
		warn "Calculating new values...\n";
		($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3) = &calculate_probabilities(($found_max_depth + 1), $opt_x, $prob_homozygous, $prob_heterozygous, $prob_heterozygous2, $prob_heterozygous3);
		&print_binomial_probabilities($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3, $Binom_OUT, ($found_max_depth + 1), $opt_x, 'F');
	}
	else { warn "May only have headers. Delete file and try and again!\n"; }	
}
else {
	($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3) = &calculate_probabilities(1, $opt_x, $prob_homozygous, $prob_heterozygous, $prob_heterozygous2, $prob_heterozygous3);
	&print_binomial_probabilities($homozygous_probability, $heterozygous_probability1, $heterozygous_probability2, $heterozygous_probability3, $Binom_OUT, 1, $opt_x, 'T');
}

### Next, calculate the cumulative probability of <= agreeing for homozygous and heterozygous 
### and >= for heterozygous and print them to a new file

my $find_min_cumulative_prob = (0);
my $Cumulative_Binom_OUT = ($dir . "/Cumulative_Binomial_Probabilities_for_$opt_p.tab");
if (-e $Cumulative_Binom_OUT) {
	my ($found_max_depth2) = 0;
	warn "Reading from file: $Cumulative_Binom_OUT\n";
	open IN2, "<$Cumulative_Binom_OUT";
	while(my $line=<IN2>) {
		chomp $line;
		my @line_parts = split /\t/, $line;
		my ($depth, $agree, $lower_tail_ProbHom, $lower_or_upper_tail_ProbHet, $lower_or_upper_tail_ProbHet2, $lower_or_upper_tail_ProbHet3) = @line_parts;
		
		next if($line =~ m/^Depth/);
		#warn "In the cumulative file, I've found $depth\t$agree\t$lower_tail_ProbHom\t$lower_or_upper_tail_ProbHet\n";
		if($depth > $found_max_depth2) { $found_max_depth2 = $depth; }
	}
	close IN2;
	warn "Max depth found in cumulative file is $found_max_depth2\n";
	if($opt_x <= $found_max_depth2) { die "Already calculated these probabilities. Done\n"; }
	elsif($opt_m eq ($found_max_depth2 + 1)) {
		warn "Calculated up to this value. Now Calculating new values\n";
		$find_min_cumulative_prob = $opt_m;
	}
	elsif(($opt_m > ($found_max_depth2 + 1)) || ($opt_m < ($found_max_depth2 + 1))) {
		warn ("Have not calculated between " . ($found_max_depth2 + 1) . " and $opt_x. Finding these first.\n");	
		$find_min_cumulative_prob = ($found_max_depth2 + 1);
	}
	else { warn "May only have headers. Delete file and try and again!\n"; }
}
else {
	if($opt_m ne 1) { warn ("Have not calculated between 1 and $opt_m. Finding these first.\n"); }
	$find_min_cumulative_prob = 1;
	warn "Making new file: $Cumulative_Binom_OUT\n";
	open OUT2, ">>$Cumulative_Binom_OUT";
	print OUT2 "Depth\tAgree\tLower tail ProbHom\tLower or upper tail ProbHet\tLower or upper tail ProbHet2\tLower or upper tail ProbHet3\n";
	close OUT2;
}

### Print cumulative probabilities

warn "Calculate cumulative probabilities between $find_min_cumulative_prob and $opt_x\n";
open OUT2, ">>$Cumulative_Binom_OUT";
for(my $agree=0; $agree<=$opt_x; $agree++) {
	for(my $depth=$find_min_cumulative_prob; $depth<=$opt_x; $depth++) {
		if($agree > $depth) { $depth=$agree; }
		
		my $lower_hom_prob1 = &tally_tail_probabilities($depth,$agree,$homozygous_probability, 'l');
		my $lower_het_prob1 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability1, 'l');
		my $upper_het_prob1 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability1, 'u');
		my $lower_het_prob2 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability2, 'l');
		my $upper_het_prob2 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability2, 'u');
		my $lower_het_prob3 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability3, 'l');
		my $upper_het_prob3 = &tally_tail_probabilities($depth,$agree,$heterozygous_probability3, 'u');
		
		## Expand
		
		my $lower_hom_prob1e = &expand($lower_hom_prob1);
		my $lower_het_prob1e = &expand($lower_het_prob1);
		my $upper_het_prob1e = &expand($upper_het_prob1);
		my $lower_het_prob2e = &expand($lower_het_prob2);
		my $upper_het_prob2e = &expand($upper_het_prob2);
		my $lower_het_prob3e = &expand($lower_het_prob3);
		my $upper_het_prob3e = &expand($upper_het_prob3);
		
		## Get the smallest probability for lower or upper tail
		
		my $prob_of_het1 = $upper_het_prob1e; 
		my $prob_of_het2 = $upper_het_prob2e; 
		my $prob_of_het3 = $upper_het_prob3e; 
		
		if($lower_het_prob1e < $upper_het_prob1e) { $prob_of_het1 = $lower_het_prob1e; }
		if($lower_het_prob2e < $upper_het_prob2e) { $prob_of_het2 = $lower_het_prob2e; }
		if($lower_het_prob3e < $upper_het_prob3e) { $prob_of_het3 = $lower_het_prob3e; }
		
		# Round to 10dp otherwise look-up tables get very large!
		
		my $rounded_homozygous = sprintf("%.10f", $lower_hom_prob1e);
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
		
		print OUT2 "$depth\t$agree\t$rounded_homozygous\t$rounded_heterozygous1\t$rounded_heterozygous2\t$rounded_heterozygous3\n";
	}
}
close OUT2;

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
	my ($homozygous_probability, $heterozygous_probability, $heterozygous_probability2, $heterozygous_probability3, $OUTFILE, $min, $max, $header) = @_;
	open OUT1, ">>$OUTFILE";
	if($header eq 'T') {
		warn "Making new file: $OUTFILE\n";
		print OUT1 "Depth\tAgree\tProb.Hom.\tProb.Het.Dip\tProb.Het.Trip1\tProb.Het.Trip2\n";
	}
	
	for(my $agree=0; $agree<=$max; $agree++) {
		for(my $depth=$min; $depth<=$max; $depth++) {
			if($agree > $depth) { $depth=$agree; }
			my $prob_hom = $$homozygous_probability{$depth}{$agree};
			my $prob_het = $$heterozygous_probability{$depth}{$agree};
			my $prob_het2 = $$heterozygous_probability2{$depth}{$agree};
			my $prob_het3 = $$heterozygous_probability3{$depth}{$agree};
			
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
	my (%homozygous_probability, %heterozygous_probability, %heterozygous_probability2, %heterozygous_probability3);
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
			
			$homozygous_probability{$depth}{$agree} = $rounded_homozygous;
			$heterozygous_probability{$depth}{$agree} = $rounded_heterozygous;
			$heterozygous_probability2{$depth}{$agree} = $rounded_heterozygous2;
			$heterozygous_probability3{$depth}{$agree} = $rounded_heterozygous3;
		}
		warn "$agree / $max\n";
	}
	$R->stopR();	
	return(\%homozygous_probability, \%heterozygous_probability, \%heterozygous_probability2, \%heterozygous_probability3);
}

sub expand {
        my $n = shift;
        $n =~ s/\[1\]|\s//g;
        return $n unless $n =~ /^(.*)e([-+]?)(.*)$/;
        my ($num, $sign, $exp) = ($1, $2, $3);
        my $sig = $sign eq '-' ? "." . ($exp - 1 + length $num) : '';
        return sprintf "%${sig}f", $n;
}
