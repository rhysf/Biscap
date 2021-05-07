#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules/";
use ascii;
my %ascii = ascii::ascii(); # a hash of ASCII symbols to values.

### r.farrer09@imperial.ac.uk

### Opening commands 

my $usage = "
Program:  Binomial SNP Caller from Pileup (BiSCaP)
Version:  0.11
Contact:  Rhys Farrer <r.farrer09\@imperial.ac.uk>
Usage:    $0 <commands>
Commands: -p\tAlignment in Pileup format 
          -r\tReference sequence in FASTA format
          
Optional: -m\tMinimum read depth to be considered a mutation [4]
          -e\tProbability of error (e.g. 0.1 or 0.01) [0.1]
          -g\tIf depth > max depth in look-up table, analyse up to max depth (y/n) [y]
          -a\tPrint out pileup homozygous Agree lines (y/n) [n]
          -i\tPrint out seperated pileup lines specifying polymorphisms in addition to default out (y/n) [n]
          -q\tRead Quality minimum cut-off (e.g. 10, 20...) [0]
          -o\tOutput folder location. default: folder of pileup
          -n\tSample name for VCF [name of pileup]
          -d\tPrint tallied read depths across genome (y/n) [n]
";
our($opt_p, $opt_r, $opt_m, $opt_e, $opt_a, $opt_q, $opt_i, $opt_o, $opt_n, $opt_g, $opt_d);
getopt('prmelsaqiog');
die $usage unless ($opt_p && $opt_r);   
if(!defined $opt_m) { $opt_m = 4;      }
if(!defined $opt_e) { $opt_e = 0.1;    }
if(!defined $opt_i) { $opt_i = 'n';    }
if(!defined $opt_q) { $opt_q = 0;      }
if(!defined $opt_a) { $opt_a = 'n';    }
if(!defined $opt_g) { $opt_g = 'y';    }
if(!defined $opt_n) { $opt_n = $opt_p; }
if(!defined $opt_d) { $opt_d = 'n';    }
if(defined $opt_o) { 
	if($opt_o !~ m/\/$/) { $opt_o .= '/'; }
	die "Output folder does not exist: $opt_o\n" unless (-d $opt_o); 
}
die "Cannot open $opt_p\n" unless (-e $opt_p);
die "Cannot open $opt_r\n" unless (-e $opt_r);
die "-m -e -q need to be numerical\n" if ((($opt_m) !~ /\d+/) || (($opt_e) !~ /\d+/) || (($opt_q) !~ /\d+/)); 
die "-e need to be numerical between 0 and 1\n" if (($opt_e > 1) || ($opt_e < 0));
die "-a needs to be y or n ($opt_a)\n" if (($opt_a !~ m/y/) && ($opt_a !~ m/n/));
die "-i needs to be y or n ($opt_i)\n" if (($opt_i !~ m/y/) && ($opt_i !~ m/n/));
my $Cumulative_Binoms = "$FindBin::Bin/Cumulative_Binomial_Probabilities_for_$opt_e.tab";
die "Cannot find $Cumulative_Binoms in directory of script. Remake probabilities using " unless (-e $Cumulative_Binoms);
my $settings = ('m-' . $opt_m . '-e-' . $opt_e . '-q-' . $opt_q);
&warn_settings($settings, $opt_p, $opt_r);
if(defined $opt_o) { warn "Output folder:   $opt_o\n"; }

### Save the length of reference to calculate uncovered bases and print VCF header
### and the sequence to pull out additional bases over insertions and deletions

my (%depth_of_coverage, %reference_length, %FASTA_ref);
$depth_of_coverage{0} = 0;
my $inseq = Bio::SeqIO->new('-file' => "<$opt_r",'-format' => 'fasta' ) ;
while (my $seq_obj = $inseq->next_seq ) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	$depth_of_coverage{0} += length($seq);
	$reference_length{$id} = length($seq);
	$FASTA_ref{$id} = $seq;
}
warn "Length of genome = $depth_of_coverage{0}\n";

### Prepare output files

my $OUTFIL = $opt_p;
if(defined $opt_o) {
	my @location_parts = split /\//, $OUTFIL;
	my $OUTFIL_name = $location_parts[(scalar(@location_parts) - 1)];
	$OUTFIL = ($opt_o . $OUTFIL_name);
}
if($opt_i eq 'y') { 
	my $HOMSNP = ($OUTFIL . '-Homozygous-SNPs-' 		. $settings . '.tab');
	my $HOMDEL = ($OUTFIL . '-Homozygous-Deletions-' 	. $settings . '.tab');
	my $HOMINS = ($OUTFIL . '-Homozygous-Insertions-' 	. $settings . '.tab');
	my $HETSNP = ($OUTFIL . '-Heterozygous-SNPs-' 		. $settings . '.tab');
	my $HETDEL = ($OUTFIL . '-Heterozygous-Deletions-' 	. $settings . '.tab');
	my $HETINS = ($OUTFIL . '-Heterozygous-Insertions-' 	. $settings . '.tab');
	open HOMSNP, ">$HOMSNP" or die "Couldn't open: $!";
	open HOMDEL, ">$HOMDEL" or die "Couldn't open: $!";
	open HOMINS, ">$HOMINS" or die "Couldn't open: $!";
	open HETSNP, ">$HETSNP" or die "Couldn't open: $!";
	open HETDEL, ">$HETDEL" or die "Couldn't open: $!";
	open HETINS, ">$HETINS" or die "Couldn't open: $!";
}
if($opt_a eq 'Y') { 
	my $HOMAGR = ($OUTFIL . '-Homozygous-agree-' 		. $settings . '.tab');
	open HOMAGR, ">$HOMAGR" or die "Couldn't open: $!"; 
}
if($opt_g ne 'y') { 
	my $OUTDIS = ($OUTFIL . '-Outside-Distribution-' 	. $settings . '.tab');
	open OUTDIS, ">$OUTDIS" or die "Couldn't open: $!"; 
}
if($opt_d eq 'y') { 
	my $DEPTHS = ($OUTFIL . '-Depth-per-base.tab');
	open DEPTHS, ">$DEPTHS" or die "Couldn't open: $!"; 
}

# mandatory output files

my $SUMMRY = ($OUTFIL . '-Summary-of-BiSCaP-' 		. $settings . '.tab');
my $ALLMUT = ($OUTFIL . '-Polymorphisms-' 		. $settings . '.VCF');
open ALLMUT, ">$ALLMUT" or die "Couldn't open: $!";
my $print_header = &print_VCF_header($0, $opt_r, \%reference_length, $settings, $opt_n);
open SUMMRY, ">$SUMMRY" or die "Couldn't open: $!";

### Save binomial distributions

my ($homozygous_probability, $heterozygous_probability, $heterozygous_probability2, $heterozygous_probability3);

warn "Reading from file: $Cumulative_Binoms\n";
open IN1, "<$Cumulative_Binoms";
my $max_Depth = 0;
while(my $line=<IN1>) {
	chomp $line;
	my @line_parts = split /\t/, $line;
	my ($depth, $agree, $lower_tail_ProbHom, $lower_or_upper_tail_ProbHet, $lower_or_upper_tail_ProbHet2, $lower_or_upper_tail_ProbHet3) = @line_parts;
	next if($line =~ m/^Depth/);
	$$homozygous_probability{$depth}{$agree} = $lower_tail_ProbHom;
	$$heterozygous_probability{$depth}{$agree} = $lower_or_upper_tail_ProbHet;
	$$heterozygous_probability2{$depth}{$agree} = $lower_or_upper_tail_ProbHet2;
	$$heterozygous_probability3{$depth}{$agree} = $lower_or_upper_tail_ProbHet3;
	if($depth > $max_Depth) { $max_Depth = $depth; }
}
close IN1;
warn "Max depth to be considered is $max_Depth. Run GBiD.pl for greater depth. Using reads up to this depth to determine genotype = $opt_g\n";

### Read from Pileup and determine the state at each base

my ($tally_HOMSNP, $tally_HOMINS, $tally_HOMDEL, $tally_HETSNP, $tally_HETINS, $tally_HETDEL) = (0, 0, 0, 0, 0, 0);
my ($tally_OUTDIS, $tally_HOMAGR, $tally_MINDEP, $tally_UNCHAR, $removed_bases) = (0, 0, 0, 0, 0);

warn "Reading from file: $opt_p...\n";
open IN2, "<$opt_p" or die "Can't open $opt_p : $!";
my @mapping_qualities;
PILEUP: while(my $line=<IN2>) {
	chomp $line;
	next PILEUP if($line =~ m/^\n$"/);
	
	my @bits = split /\t/, $line;
	my ($contig, $pos, $ref_base, $depth, $aligned, $read_qual) = @bits;
	die "Pileup line not separated by tab or missing column: $line\n" if (@bits != 6);
	my ($saved_depth, $saved_ref) = ($depth, $ref_base);
	$depth_of_coverage{$depth}++;
	$depth_of_coverage{0}--;
	
	### Record the mapping qualities
	
	my ($aligned2, $end_of_reads, $map_qualities) = &find_mapping_quality($aligned);
	push (@mapping_qualities, @{$map_qualities});
	
	if(($ref_base ne 'N') && ($depth >= $opt_m)) {
		
		### Remove lines with read-depths outside binomial distributions 
		
		if(($opt_g eq 'n') && (!exists $$homozygous_probability{$depth})) {
			print OUTDIS "$line\n"; 
			$tally_OUTDIS++; 
			next PILEUP;
		}
		
		### Count indels and compare probabilities
		
		my ($aligned3, $depth2, $indel_min_depth_pass, %indel_variants) = &count_indels($aligned2, $opt_m, $max_Depth, $depth, $opt_g);
		my ($consensus_indel, $indel_genotype, $indel_info) = &compare_probability(\%indel_variants, $depth2, $opt_e, $indel_min_depth_pass, $saved_ref);
		
		### Count bases and compare probabilities
		
		if($opt_q ne 0) { ($aligned3, $removed_bases, $read_qual) = &find_read_quality($aligned3, $read_qual, $opt_q); }
		my ($nt_min_depth_pass, $depth3, %nt_variants) = &count_bases($saved_ref, $depth2, $aligned3, $opt_m);
		my ($consensus_nt, $nt_genotype, $nt_info) = &compare_probability(\%nt_variants, $depth3, $opt_e, $nt_min_depth_pass, $saved_ref);
		$consensus_nt =~ s/!//g;

		### Deletions in VCF
		
		if($consensus_indel =~ m/-/) {
			($consensus_indel, $indel_genotype, $ref_base) = &deletions_in_VCF($consensus_indel, $ref_base, $contig, $pos, $indel_genotype, $nt_genotype, $consensus_nt, \%FASTA_ref);
		}

		### Insertions in VCF
		
		if($consensus_indel =~ m/\+/) {
			($consensus_indel, $indel_genotype, $indel_info) = &insertions_in_VCF($saved_ref, $consensus_indel, $consensus_nt, $indel_genotype, $nt_genotype, $indel_info);
		}
		
		### Pass min depth
		
		if(($indel_min_depth_pass eq 'N') && ($nt_min_depth_pass eq 'N')) {
			$tally_MINDEP++;
			next PILEUP;
		}
		
		### Reference base
		
		if(($indel_genotype eq 0) && ($nt_genotype eq 0)) { 
			$tally_HOMAGR++; 
			if($opt_a eq 'Y') { print HOMAGR "$line\n"; }
			next PILEUP; 
		}
		
		### Calculate mean mapping quality, base quality and max mapping quality 
		
		my ($mean_mapping_quality, $mean_base_quality) = (0,0);
		if($opt_q) {
			my $max_mapping_quality = (sort { $b <=> $a } @mapping_qualities)[0];
			my ($sum_base_qual, $sum_map_qual) = (0, 0);
			foreach(@mapping_qualities) { $sum_map_qual += $_; }
			foreach(@{$read_qual}) { $sum_base_qual += $_; }
			$mean_mapping_quality = int($sum_map_qual / scalar(@mapping_qualities));
			$mean_base_quality = int($sum_base_qual / scalar(@{$read_qual}));
		}
		
		### Record the polymorphism
		# To start with either indel or base-change
		# base-change
		if(($indel_genotype eq 0) && ($nt_genotype ne 0)) {
			
			### Ambigious base: Last check in case of 66:33 variant, only ref base within probability.
			if($ref_base eq $consensus_nt) { next PILEUP; }
			
			my $VCF_line = "$contig\t$pos\t.\t$ref_base\t$consensus_nt\t$mean_base_quality\t";
			$VCF_line .= ".\t$nt_info\tGT:DP:MMQ\t$nt_genotype:$saved_depth:$mean_mapping_quality";
			print ALLMUT "$VCF_line\n";
		} 
		# indel
		else {
			my $VCF_line = "$contig\t$pos\t.\t$ref_base\t$consensus_indel\t$mean_base_quality\t";
			$VCF_line .= ".\t$indel_info\tGT:DP:MMQ\t$indel_genotype:$saved_depth:$mean_mapping_quality";
			print ALLMUT "$VCF_line\n";
		}
		
		# Tally polymorphisms
		# Insertions
		if(($indel_info =~ m/INDEL/) && ($consensus_indel =~ m/[ACTG][ACTG]/)) {
			if($indel_genotype =~ m/\//) {
				if($opt_i eq 'y') { print HETINS "$line\n"; }
				$tally_HETINS++; 
			} else {
				if($opt_i eq 'y') { print HOMINS "$line\n"; }
				$tally_HOMINS++; 
			}
		}
		# Deletions
		if(($indel_info =~ m/INDEL/) && ($ref_base =~ m/[ACTG][ACTG]/)) {
			if($indel_genotype =~ m/\//) {
				if($opt_i eq 'y') { print HETDEL "$line\n"; }
				$tally_HETDEL++; 
			} else {
				if($opt_i eq 'y') { print HOMDEL "$line\n"; }
				$tally_HOMDEL++; 
			}
		}
		# SNP 
		if(($indel_info !~ m/INDEL/) && ($consensus_nt !~ m/[ACTG][ACTG]/) && ($ref_base !~ m/[ACTG][ACTG]/)) {
			if($nt_genotype =~ m/\//) {
				if($opt_i eq 'y') { print HETSNP "$line\n"; }
				$tally_HETSNP++; 
			} else {
				if($opt_i eq 'y') { print HOMSNP "$line\n"; }
				$tally_HOMSNP++; 
			}
		}
	}
	elsif($depth < $opt_m) { $tally_MINDEP++; }
	
	### Remove end of read mapping scores
	
	foreach(sort { $b <=> $a } @{$end_of_reads}) { 
		if(defined $mapping_qualities[$_]) {
			splice (@mapping_qualities, $_, 1); 
		}
	}
}
close IN2;

### Print the Read Depth of coverage

if($opt_d eq 'y') { 
	print DEPTHS "\"Read Depth\"\t\"Tally\"\n";
	foreach my $depths(sort {$a<=>$b} keys %depth_of_coverage) {
		my $tally = $depth_of_coverage{$depths};	
		print DEPTHS "$depths\t$tally\n";
	}
	close DEPTHS;
}

### Summary of SNP calling

my $warn_summary = &warn_summary($tally_HOMSNP, $tally_HOMINS, $tally_HOMDEL, 
	$tally_HOMAGR, $tally_HETSNP, $tally_HETINS, $tally_HETDEL, $tally_UNCHAR, 
	$tally_OUTDIS, $tally_MINDEP, $max_Depth, $removed_bases, $opt_q, $opt_m);
print SUMMRY "$warn_summary\n";

### Close output files

close ALLMUT; close HOMSNP; close HOMDEL; close HOMINS; close SUMMRY; close HETSNP; close HETDEL; close HETINS; 
if($opt_a eq 'Y') { close HOMAGR; }

### Sub routines

sub find_mapping_quality {
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

	return ($aligned, \@positions_of_ends, \@mapping_qualities);
}

sub count_bases {
	my ($ref_base, $depth, $aligned_in_find_base, $mindepth) = @_;
	
	my $uncharacterised = 'N';
	my $full_consensus = '';
	my %variants;
	my $type;
	
	### Count each aligned base
	
	my $match = ($aligned_in_find_base =~ tr/\.|\,//);
	my $mismatchA = ($aligned_in_find_base =~ tr/A|a//);
	my $mismatchC = ($aligned_in_find_base =~ tr/C|c//);
	my $mismatchT = ($aligned_in_find_base =~ tr/T|t//);
	my $mismatchG = ($aligned_in_find_base =~ tr/G|g//);
	my $deletion = ($aligned_in_find_base =~ tr/\*//);
	$depth -= $deletion;
	
	### Identify the three most common bases
	
	my (@nbase, @type_of_base);
	my $found = '';
	for(my $i=0; $i<3; $i++) {
		my $result = &MAX($match, $mismatchA, $mismatchC, $mismatchT, $mismatchG);
		if(($result eq 0) && ($found !~ m/$ref_base/)) { 
			$nbase[$i] = $match;
			$type_of_base[$i] = $ref_base;
			$match = 0;
			$found .= $ref_base;
		}
		elsif(($result eq 1) && ($found !~ m/A/)) { 
			$nbase[$i] = $mismatchA;
			$type_of_base[$i] = 'A';
			$mismatchA = 0;
			$found .= 'A';
		}
		elsif(($result eq 2) && ($found !~ m/C/)) { 
			$nbase[$i] = $mismatchC;
			$type_of_base[$i] = 'C';
			$mismatchC = 0;
			$found .= 'C';
		}
		elsif(($result eq 3) && ($found !~ m/T/)) { 
			$nbase[$i] = $mismatchT;
			$type_of_base[$i] = 'T';
			$mismatchT = 0;
			$found .= 'T';
		}
		elsif(($result eq 4) && ($found !~ m/G/)) { 
			$nbase[$i] = $mismatchG;
			$type_of_base[$i] = 'G';
			$mismatchG = 0;
			$found .= 'G';
		}
		else {
			$nbase[$i] = 0;
			$type_of_base[$i] = 'N';
		}
	}
	
	# Min depth filters on most common base
	if($nbase[0] < $mindepth) {
		$type = 'N';
		return ($type, %variants);
	}	
	# Otherwise use all available information to determine genotype
	$type = '!';
	$variants{$type}{'1'}{$type_of_base[0]} = $nbase[0];
	$variants{$type}{'2'}{$type_of_base[1]} = $nbase[1];
	$variants{$type}{'3'}{$type_of_base[2]} = $nbase[2];
	return ($type, $depth, %variants);
}

sub count_indels {
	my ($aligned, $mindepth, $maxdepth, $depth, $optg) = @_;
	my (%indel_type, %indel_count, %variants);
	my @newaligned = split //, $aligned;
	my $number_of_removed_bases = 0;
	my ($type1, $type2, $type3) = ('.','.','.');

	### Count and remove indels

	$indel_count{'+'} = 0; $indel_count{'-'} = 0;
	COUNTANDREMOVEINDELS: while($aligned =~ m/([\+|\-])(\d+)([ACGTNacgtn]+)/ig) {
		
		# if opt_g, use up to max depth as representative
		last COUNTANDREMOVEINDELS if(($optg eq 'y') && ((pos($aligned) - $number_of_removed_bases) > $maxdepth));
		
		# Count 
		my $length_of_indel = $2;
		my $bases = $3;
		my $indel = substr $bases, 0, $2;
		$indel =~ tr/a-z/A-Z/;
		my $type;
		if($1 =~ m/\-/) { $type = '-'; }
		else            { $type = '+'; }
		$indel_type{$type}{$indel}++; 
		$indel_count{$type}++;
		
		# Remove 
		my $string_length = (1 + length($length_of_indel) + $length_of_indel);
		my $pos_aligned = pos($aligned);
		
		if(length($bases) != $length_of_indel) { 
			$pos_aligned -= (length($bases) - $length_of_indel);
			
		}
		
		$number_of_removed_bases += $string_length;
		splice(@newaligned, ($pos_aligned - $number_of_removed_bases), $string_length);
		
	}
	if(length($aligned) != scalar(@newaligned)) { $aligned = join('', @newaligned); }

	### if opt_g, use up to max depth as representative
	
	if(($opt_g eq 'y') && ($depth > $maxdepth)) {
		my $new_aligned = (substr $aligned, 0, $maxdepth);
		$new_aligned =~ s/[\+|\-|\d]+//g;
		$aligned = $new_aligned;
		$depth = length($aligned);
	}
	
	### Identify the three most common indels
	
	if($indel_count{'+'} < $indel_count{'-'}) { $type1 = '-'; }
	else { $type1 = '+'; }
	my ($number_of_indel, $number_of_indel2, $number_of_indel3) = (0, 0, 0);
	my ($type_of_indel, $type_of_indel2, $type_of_indel3) = ('', '', '', '');
	foreach my $indel_nts(sort keys %{$indel_type{$type1}}) {
		my $frequency_of_indel = $indel_type{$type1}{$indel_nts};
		if($frequency_of_indel > $number_of_indel) { 
			$number_of_indel = $frequency_of_indel;
			$type_of_indel = $indel_nts;
			$indel_count{$type1} -= $frequency_of_indel;
		}
	}
	
	# Check if different type for 2nd most common
	
	if($indel_count{'+'} < $indel_count{'-'}) { $type2 = '-'; }
	else { $type2 = '+'; }
	IND2: foreach my $indel_nts(sort keys %{$indel_type{$type2}}) {
		my $frequency_of_indel = $indel_type{$type2}{$indel_nts};
		next IND2 if ($indel_nts eq $type_of_indel);
		if($frequency_of_indel > $number_of_indel2) { 
			$number_of_indel2 = $frequency_of_indel;
			$type_of_indel2 = $indel_nts;
			$indel_count{$type2} -= $frequency_of_indel;
		}
	}
	
	# Check if different type for 3rd most common
	
	if($indel_count{'+'} < $indel_count{'-'}) { $type3 = '-'; }
	else { $type3 = '+'; }
	IND3: foreach my $indel_nts(sort keys %{$indel_type{$type3}}) {
		my $frequency_of_indel = $indel_type{$type3}{$indel_nts};
		next IND3 if (($indel_nts eq $type_of_indel) || ($indel_nts eq $type_of_indel2));
		if($frequency_of_indel > $number_of_indel3) { 
			$number_of_indel3 = $frequency_of_indel;
			$type_of_indel3 = $indel_nts;
		}
	}
	
	# Min depth filter on most common indel
	if($number_of_indel < $mindepth) {
		return ($aligned, $depth, 'N', %variants);
	} 
	
	# Otherwise use all available information to determine genotype
	$variants{$type1}{'1'}{$type_of_indel} = $number_of_indel;
	$variants{$type2}{'2'}{$type_of_indel2} = $number_of_indel2;
	$variants{$type3}{'3'}{$type_of_indel3} = $number_of_indel3;
	return ($aligned, $depth, 'Y', %variants);
	
}

sub compare_probability {
	my ($variant_counts, $depth, $error_prob, $min_depth_pass, $refbase) = @_;
	my ($consensus, $genotype, $info) = ('',0,'.');	
	
	# Pass min depth
	if($min_depth_pass eq 'N') { return ($consensus, $genotype, $info); }
	
	### Pull info from hash
	
	my ($type_variant1, $type_variant2, $type_variant3) = ('', '', '');
	my ($type1, $type2, $type3) = ('', '', '');
	my ($num1, $num2, $num3) = (0, 0, 0);
	foreach my $variant_type(keys %{$variant_counts}) {
		foreach my $prevelance(keys %{$$variant_counts{$variant_type}}) {
			foreach my $type(keys %{$$variant_counts{$variant_type}{$prevelance}}) {
				my $count = $$variant_counts{$variant_type}{$prevelance}{$type};
				
				if($prevelance eq 1) {
					$type1 = $type;
					$num1 = $count;
					$type_variant1 = $variant_type;
				}
				if($prevelance eq 2) {
					$type2 = $type;
					$num2 = $count;
					$type_variant2 = $variant_type;
				}
				if($prevelance eq 3) {
					$type3 = $type;
					$num3 = $count;
					$type_variant3 = $variant_type;
				}
			}
		}
	}
	
	### Get probabilities 
	
	# hom
	my $bin_hom_agree = $$homozygous_probability{$depth}{$num1};
	# dip. het
	my $bin_het_agree = $$heterozygous_probability{$depth}{$num1};
	my $bin_het_agree2 = $$heterozygous_probability{$depth}{$num2};
	# trip. het (33:33:33)
	my $bin_het_agree3 = $$heterozygous_probability2{$depth}{$num1};
	my $bin_het_agree4 = $$heterozygous_probability2{$depth}{$num2};
	my $bin_het_agree5 = $$heterozygous_probability2{$depth}{$num3};
	# trip. het (66:33)
	my $bin_het_agree6 = $$heterozygous_probability3{$depth}{$num1};
	my $result = &MAX($bin_hom_agree, $bin_het_agree, $bin_het_agree2, $bin_het_agree3, $bin_het_agree4, $bin_het_agree5, $bin_het_agree6);
	
	### Compare probabilities
	
	# Homozygous variant
	if($result eq 0) { 
		if($bin_hom_agree >= $error_prob) {
			#warn "homozygous variant or refbase\n";
			$consensus = ($type_variant1 . $type1);
			if(($type_variant1 !~ m/-|\+/) && ($type1 eq $refbase)) { $genotype = '0'; }
			else { $genotype = '1'; }
		}
	}
	# Heterozygous variant (biallelic)
	elsif($result eq 1) { 
		if($bin_het_agree >= $error_prob) {
			
			# biallelic with other allele specifying no variant
			if($bin_het_agree2 < $error_prob) { 
				if(($type_variant1 !~ m/-|\+/) && ($type1 eq $refbase)) { $genotype = '0'; }
				else {
					$consensus = ($type_variant1 . $type1);
					$genotype = '0/1'; 
				}
			}
			# biallelic with other allele specifying variant
			elsif($bin_het_agree2 >= $error_prob) { 
				if(($type_variant1 !~ m/-|\+/) && ($type1 eq $refbase)) { 
					$consensus = ($type_variant2 . $type2);
					$genotype = '0/1'; 
				}
				elsif(($type_variant1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type_variant1 . $type1);
					$genotype = '0/1'; 
				}
				else {
					$consensus = ($type_variant1 . $type1 . ',' . $type_variant2 . $type2);
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
				$consensus = ($type_variant1 . $type1);
				$info = 'TRIPROB:0/1/?;';
				$genotype = '0/1'; 
			}
			# 3rd allele specifies no variant
			elsif(($bin_het_agree4 >= $error_prob) && ($bin_het_agree5 < $error_prob)) {
				$consensus = ($type_variant1 . $type1 . ',' . $type_variant2 . $type2);
				$info = 'TRIALLELIC:0/1/2;';
				$genotype = '0/1/2'; 
			}
			# 3rd allele specifies variant
			elsif(($bin_het_agree4 >= $error_prob) && ($bin_het_agree5 >= $error_prob)) {
				if(($type_variant1 !~ m/-|\+/) && ($type1 eq $refbase)) {
					$consensus = ($type_variant2 . $type2 . ',' . $type_variant3 . $type3);
					$info = 'TRIALLELIC:0/1/2;';
					$genotype = '0/1/2';
				}
				elsif(($type_variant1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type_variant1 . $type1 . ',' . $type_variant3 . $type3);
					$info = 'TRIALLELIC:0/1/2;';
					$genotype = '0/1/2';
				}
				else {
					$consensus = ($type_variant1 . $type1 . ',' . $type_variant2 . $type2 . ',' . $type_variant3 . $type3);
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
				$consensus = ($type_variant1 . $type1);
				$info = 'TRIPROB:0/1/1;';
				$genotype = '0/1'; 	
			}
			# 2nd allele specifies variant
			elsif($bin_het_agree4 >= $error_prob) {
				if(($type_variant1 !~ m/-|\+/) && ($type1 eq $refbase)) {
					$consensus = ($type_variant2 . $type2);
					$info = 'TRIPROB:0/0/1;';
					$genotype = '0/1';
				}
				elsif(($type_variant1 !~ m/-|\+/) && ($type2 eq $refbase)) { 
					$consensus = ($type_variant1 . $type1);
					$info = 'TRIPROB:0/1/1;';
					$genotype = '0/1';
				}
				else {
					if($num2 != 0) {
						$consensus = ($type_variant1 . $type1 . ',' . $type_variant2 . $type2);
						$info = 'TRIPROB:1/1/2;';
						$genotype = '1/2';
					} else {
						$consensus = ($type_variant1 . $type1);
						$info = 'TRIPROB:1/1/2;';
						$genotype = '0/1';
					}
				}
			}
		}
	}
	if((defined $genotype) && ($type_variant1 =~ m/-|\+/)) { 
		if($info eq '.') { $info = 'INDEL'; }
		else { $info .= 'INDEL'; }
	}
	if((!defined $consensus) || ($consensus eq '')) { $consensus = $refbase; }
	return ($consensus, $genotype, $info);
}

sub find_read_quality {
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

sub insertions_in_VCF {
	my ($nt_ref, $indel_consensus, $nt_consensus, $indel_geno, $nt_geno, $indel_info) = @_;
	my ($consensus_both, $both_genotype) = ('','');	
	$indel_consensus =~ s/\+//g;
	
	# Insertion with ref or SNP
	if(($nt_geno eq 0) || ($nt_geno eq 1)) {
		# Homozygous variant
		if($indel_geno eq 1) { $consensus_both = ($nt_consensus . $indel_consensus); }
		# Heterozygous variant
		else { 
			my @parts = split /,/, $indel_consensus;
			foreach(@parts) { $consensus_both .= ($nt_consensus . $_ . ','); }
			$consensus_both =~ s/,$//;
		}
		$both_genotype = $indel_geno;
	}
	# Insertion with heterozygous base
	elsif(($nt_geno eq '0/1') && ($indel_geno eq '0/1')) {
		my @parts = split /,/, $indel_consensus;
		$consensus_both .= ($nt_consensus . ',' . $nt_ref . $indel_consensus);
		$both_genotype = '1/2';
		
		# Found a 1/3 prob for indel and a het snp
		if($indel_info =~ m/\?/) {
			$indel_info = 'TRIALLELIC:1/1/2;INDEL';
		}
		
	}
	return($consensus_both, $both_genotype, $indel_info);
}

sub deletions_in_VCF {
	my ($consensus_indel, $ref_base, $refcontig, $pos, $indel_geno, $nt_geno, $consen_nt, $FASTA_seq_ref) = @_;
	
	# Find longest deletion for reference base (ignore refs)
			
	my $max_length = 0;
	my @dels = split /,/, $consensus_indel;
	$consensus_indel = '';
	DELS: foreach(@dels) { 
		next DELS if ($_ =~ m/\+/);
		
		# -1 for the '-'
		if($_ =~ m/-/) {
			if((length($_) - 1) > $max_length) { $max_length = (length($_) - 1); }
		}
	}
	
	# Query ref
	my $contig = $$FASTA_seq_ref{$refcontig};
	my $new_ref_base =  substr $contig, ($pos - 1), ($max_length + 1);
	my @full_seq_parts = split //, $new_ref_base;
	
	# Deletions with a heterozygous base
	if(($consen_nt ne $ref_base) && (($indel_geno =~ m/0\/1/) && ($indel_geno =~ m/0\/1/))) {
		push @dels, $consen_nt;
		$indel_geno = '1/2';
	}
	
	
	# If there are more than 1 length for deletion, get what remains of the longest deletion
	
	DELS: foreach(@dels) {
		# for het del/ins
		if($_ =~ m/\+/) {
			$_ =~ s/\+//;
			$consensus_indel .= ($new_ref_base . $_ . ',');
			next DELS;	
		}
		# for dels
		elsif($_ =~ m/-/) {
			$_ =~ s/-//;
			if(length($_) eq $max_length) { $consensus_indel .= $ref_base . ','; }
			else { 
				$consensus_indel .= $ref_base;
				my $additional_sequen = '';
				for(my $i=$max_length; $i>length($_); $i--) { 
					$additional_sequen = ($full_seq_parts[($i)] . $additional_sequen); 
				}
				$consensus_indel .= ($additional_sequen . ',');
			}
		}
		# For het positions
		else {
			$consensus_indel .= $consen_nt;
			for(my $i = 1; $i<length($new_ref_base); $i++) {
				$consensus_indel .= $full_seq_parts[($i)];
			}
		}
	}
	$consensus_indel =~ s/,$//;
	
	### New long (non-del) ref base
	
	$ref_base = $new_ref_base;
	
	return ($consensus_indel, $indel_geno, $ref_base);
}

sub MAX {
	my (@values) = @_;
	my ($position_of_number, $big_number) = (0, 0);
	my @positions_of_numbers;
	for(my $i=0; $i<(scalar(@values) - 1); $i++) {
		if(($values[$i] > $values[$i+1]) && ($values[$i] > $big_number)) {  $position_of_number = $i; $big_number = $values[$i]; }
		elsif(($values[$i+1] > $values[$i]) && ($values[$i+1] > $big_number)) { $position_of_number = ($i+1); $big_number = $values[($i+1)]; }
	}
	return $position_of_number;
}

sub print_VCF_header {
	my ($source, $reference, $FASTA_length_ref, $settings, $sample_name) = @_;
	my @reference_name_parts = split /\//, $reference;
	my $reference_file_name = $reference_name_parts[(scalar(@reference_name_parts) - 1)];
	my @source_name_parts = split /\//, $source;
	my $source_file_name = $source_name_parts[(scalar(@source_name_parts) - 1)];
	
	print ALLMUT ('##fileformat=VCFv4.1' . "\n");
	print ALLMUT ('##source=' . "$source_file_name\n");
	print ALLMUT ('##' . "$settings\n");
	print ALLMUT ('##reference=' . "$reference_file_name\n");
	foreach my $contig_name(sort keys(%{$FASTA_length_ref})) {
		my $length = $$FASTA_length_ref{$contig_name};
		print ALLMUT ('##contig=<ID=' . $contig_name . ',length=' . $length . '>' . "\n");
	}
	print ALLMUT ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' . "\n");
	print ALLMUT ('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' . "\n");
	print ALLMUT ('##FORMAT=<ID=MMQ,Number=1,Type=Integer,Description="Mean mapping quality score">' . "\n");
	print ALLMUT ('##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">' . "\n");
	print ALLMUT ('#' . "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n");	
	
	return 1;
}

sub warn_settings {
	my ($param_full, $pileup_name, $ref) = @_;
	my @par = split /-/, $param_full;
	my ($min_depth, $error_prob, $qual) = ($par[1], $par[3], $par[5]);
	my $settings = "Settings:\n";
	$settings .= "Pileup:          $pileup_name\n";
	$settings .= "Reference FASTA: $ref\n";
	$settings .= "Min Read Depth:  $min_depth\n";
	$settings .= "Predicted Error: $error_prob\n";
	$settings .= "Quality Filter:  $qual";
	warn "$settings\n";
	return 1;	
}

sub warn_summary { 
	my ($HOMSNP, $HOMINS, $HOMDEL, $HOMAGR, $HETSNP, $HETINS, $HETDEL, $UNCHAR, $OUTDIS, $MINDEP, $MAXDEP, $REMBAS, $optq, $optm) = @_;
	my $summary = "Finished\n\n-----------------------------------------\n";
	$summary .= "Homozygous SNPs identified: $HOMSNP\n";
	$summary .= "Homozygous INSs identified: $HOMINS\n";
	$summary .= "Homozygous DELs identified: $HOMDEL\n";
	$summary .= "Homozygous agree with ref : $HOMAGR\n";
	$summary .= "Heterozygous SNPs identified: $HETSNP\n"; 
	$summary .= "Heterozygous INSs identified: $HETINS\n"; 
	$summary .= "Heterozygous DELs identified: $HETDEL\n"; 
	$summary .= "-----------------------------------------\n";
	$summary .= "Unable to resolve (number of bases): $UNCHAR\n";
	
	if($tally_OUTDIS >= 1) { 
		$summary .= "Bases unresolved due to being greater than lookup table ($MAXDEP): $OUTDIS (increase by using GBiD.pl)\n"; 
	}
	$summary .= "\nIgnored bases lower than Quality $optq: $REMBAS\n";
	$summary .= "Ignored bases lower than min depth $optm: $MINDEP\n";
	warn "$summary\n";
	return 1;
}
