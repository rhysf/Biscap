#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::Tools::GFF;

# r.farrer09@imperial.ac.uk

### Opening commands 

my $usage = "
Usage:    $0 <commands>

Commands: -m\tSNP-calls from test.pileup in Variant Call Format (VCF)
";
our($opt_m);
getopt('m');
die $usage unless ($opt_m);  
die "Cannot open $opt_m : $!" unless (-e $opt_m);


### Go through VCF and tally known sites

my (@FP_lines, @FN_lines, @FN_lines_no_cov);
my $location_of_genotype_in_format = 0;
my @saved_results;

open IN1, "<$opt_m";
warn "Cataloging SNPs in from $opt_m...\n\n";
FOUNDPOLY: while (my $line = <IN1>) {
   	chomp $line;

   	### Get info from VCF
	
	next FOUNDPOLY if ($line =~ m/^\#/);
	my @bits = split /\t/, $line;
   	next FOUNDPOLY if (@bits < 9);
   	die "$0 Currently not made to check more than a VCF for single isolate\n" if (@bits > 10);
   	my ($supercontig, $position, $id, $reference, $consensus, $cons_qual, $filter, $info, $format, $sample_info) = @bits;
	my @sample_info_parts = split /:/, $sample_info;
	my $alleles = $sample_info_parts[$location_of_genotype_in_format];
	
	### Check calls
	
	# Homozygous SNPs
	if($position eq 6) {
		if(($reference eq 'T') && ($consensus eq 'A') && ($alleles eq 1)) { push @saved_results, "Homozygous SNP (T -> A) Checked"; }
		else { push @saved_results, "SNP (A) ERROR!\n"; }
	}
	elsif($position eq 7) {
		if(($reference eq 'G') && ($consensus eq 'C') && ($alleles eq 1)) { push @saved_results, "Homozygous SNP (G -> C) Checked"; }
		else { push @saved_results, "SNP (C) ERROR!\n"; }
	}
	elsif($position eq 8) {
		if(($reference eq 'C') && ($consensus eq 'G') && ($alleles eq 1)) { push @saved_results, "Homozygous SNP (C -> G) Checked"; }
		else { push @saved_results, "SNP (G) ERROR!\n"; }
	}
	elsif($position eq 9) {
		if(($reference eq 'A') && ($consensus eq 'T') && ($alleles eq 1)) { push @saved_results, "Homozygous SNP (A -> T) Checked\n"; }
		else { push @saved_results, "SNP (T) ERROR!\n"; }
	}
	
	# DELS
	elsif($position eq 11) {
		if(($reference eq 'CA') && ($consensus eq 'C') && ($alleles eq 1)) { push @saved_results, "Homozygous deletion (CA -> A) Checked"; }
		else { push @saved_results, "Homozygous deletion (CA -> A) ERROR!\n"; }
	}
	elsif($position eq 12) {
		if(($reference eq 'AA') && ($consensus eq 'A') && ($alleles eq '0/1')) { push @saved_results, "Bi-allelic heterozygous deletion with ref. (AA -> A) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous deletion with ref. (AA -> A) ERROR!\n"; }
	}
	elsif($position eq 13) {
		if(($reference eq 'ATG') && ($consensus eq 'AG,A') && ($alleles eq '1/2')) { push @saved_results, "Bi-allelic heterozygous deletion (ATG -> A, ATG -> A-G) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous deletion (ATG -> A, ATG -> A-G) ERROR!\n"; }
	}
	elsif($position eq 14) {
		if(($reference eq 'TG') && ($consensus eq 'T') && ($alleles eq '0/1') && ($info =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion with ref. (TG -> T, TG -> T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion with ref. (TG -> T, TG -> T) ERROR!\n"; }
	}
	elsif($position eq 15) {
		if(($reference eq 'GTA') && ($consensus eq 'GA,G') && ($alleles eq '1/2') && ($info =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion (GTA -> G, GTA -> G-A, GTA -> G-A) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion (GTA -> G, GTA -> G-A, GTA -> G-A) ERROR!\n"; }
	}
	elsif($position eq 16) {
		if(($reference eq 'TAA') && ($consensus eq 'TA,T') && ($alleles eq '0/1/2') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous deletion with ref. (TAA -> T, TAA -> T-A) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous deletion with ref. (TAA -> T, TAA -> T-A) ERROR!\n"; }
	}
	elsif($position eq 17) {
		if(($reference eq 'AAGG') && ($consensus eq 'AGG,AG,A') && ($alleles eq '1/2/3') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous deletion (AAGG -> A, AAGG -> A--G, AAGG -> A-GG) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous deletion (AAGG -> A, AAGG -> A--G, AAGG -> A-GG) ERROR!\n"; }
	}
	
	# Heterozygous positions
	elsif($position eq 20) {
		if(($reference eq 'G') && ($consensus eq 'C') && ($alleles eq '0/1')) { push @saved_results, "Heterozygous position with ref. (G -> G/C) Checked"; }
		else { push @saved_results, "Heterozygous position with ref. (G -> G/C) ERROR!\n"; }
	}
	elsif($position eq 21) {
		if(($reference eq 'T') && (($consensus eq 'C,G') || ($consensus eq 'G,C')) && ($alleles eq '1/2')) { push @saved_results, "Heterozygous position (T -> G/C) Checked"; }
		else { push @saved_results, "Heterozygous position (T -> G/C) ERROR!\n"; }
	}
	elsif($position eq 22) {
		if(($reference eq 'C') && ($consensus eq 'G') && ($alleles eq '0/1') && ($info =~ m/TRIPROB:0\/1\/1/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position with ref. (C -> C/G/G) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position with ref. (C -> C/G/G) ERROR!\n"; }
	}
	elsif($position eq 23) {
		if(($reference eq 'C') && ($consensus eq 'T') && ($alleles eq '0/1') && ($info =~ m/TRIPROB:0\/0\/1/)) { push @saved_results, "Bi-allelic (tri. prob.) position with ref. (C -> C/C/T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) position with ref. (C -> C/C/T) ERROR!\n"; }
	}
	elsif($position eq 24) {
		if(($reference eq 'G') && ($consensus eq 'C,T') && ($alleles eq '1/2') && ($info =~ m/TRIPROB:1\/1\/2/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position (G -> C/C/T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position (G -> C/C/T) ERROR!\n"; }
	}
	elsif($position eq 25) {
		if(($reference eq 'C') && (($consensus eq 'T,G') || ($consensus eq 'G,T')) && ($alleles eq '0/1/2')) { push @saved_results, "Tri-allelic heterozygous position with ref. (C -> C/T/G) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous position with ref. (C -> C/T/G) ERROR!\n"; }
	}
	elsif($position eq 26) {
		if(($reference eq 'G') && (($consensus eq 'A,C,T') || ($consensus eq 'A,T,C') || ($consensus eq 'C,A,T') || ($consensus eq 'T,A,C')) && ($alleles eq '1/2/3')) { push @saved_results, "Tri-allelic heterozygous position (G -> A/T/C) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous position (G -> A/T/C) ERROR!\n"; }
	}

	# INS
	elsif($position eq 30) {
		if(($reference eq 'A') && ($consensus eq 'AA') && ($alleles eq 1)) { push @saved_results, "Homozygous insertion (A -> AA) Checked"; }
		else { push @saved_results, "Homozygous insertion (A -> AA) ERROR!\n"; }
	}
	elsif($position eq 31) {
		if(($reference eq 'T') && ($consensus eq 'TA') && ($alleles eq '0/1')) { push @saved_results, "Bi-allelic heterozygous insertion with ref. (T -> TA) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous insertion with ref. (T -> TA) ERROR!\n"; }
	}
	elsif($position eq 32) {
		if(($reference eq 'C') && ($consensus eq 'CT,CTG') && ($alleles eq '1/2')) { push @saved_results, "Bi-allelic heterozygous insertion (C -> CT, C -> CTG) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous insertion (C -> CT, C -> CTG) ERROR!\n"; }
	}
	elsif($position eq 33) {
		if(($reference eq 'T') && ($consensus eq 'TG') && ($alleles eq '0/1') && ($info =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion with ref. (T -> TG, T -> TG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion with ref. (T -> TG, T -> TG) ERROR!\n"; }
	}
	elsif($position eq 34) {
		if(($reference eq 'G') && ($consensus eq 'GT,GTA') && ($alleles eq '1/2') && ($info =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion (G -> GT, G -> GTA) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion (G -> GT, G -> GTA) ERROR!\n"; }
	}
	elsif($position eq 35) {
		if(($reference eq 'T') && ($consensus eq 'TA,TAA') && ($alleles eq '0/1/2') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion with ref. (T -> TA, T -> TAA) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion with ref. (T -> TA, T -> TAA) ERROR!\n"; }
	}
	elsif($position eq 36) {
		if(($reference eq 'T') && ($consensus eq 'TA,TAG,TAGG') && ($alleles eq '1/2/3') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion (T -> TA, T -> TAG, T -> TAGG) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion (T -> TA, T -> TAG, T -> TAGG) ERROR!\n"; }
	}
	
	###  MIXES
	# INDELS
	elsif($position eq 41) {
		if(($reference eq 'TG') && ($consensus eq 'TGTG,T') && ($alleles eq '1/2')) { push @saved_results, "Heterozygous insertion and deletion (TG -> T, TG -> TGTG) Checked"; }
		else { push @saved_results, "Heterozygous insertion and deletion (TG -> T, TG -> TGTG) ERROR!\n"; }
	}
	elsif($position eq 42) {
		if(($reference eq 'GT') && ($consensus eq 'GTTG,G') && ($alleles eq '1/2') && ($info =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion and deletion (GT -> GTTG, GT -> GTTG, GT -> G) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion and deletion (GT -> GTTG, GT -> GTTG, GT -> G) ERROR!\n"; }
	}
	elsif($position eq 43) {
		if(($reference eq 'TT') && ($consensus eq 'TTTG,T') && ($alleles eq '0/1/2') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion, deletion and ref. (TT -> TTTG, TT -> T) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion, deletion and ref. (TT -> TTTG, TT -> T) ERROR!\n"; }
	}
	elsif($position eq 44) {
		if(($reference eq 'TG') && ($consensus eq 'TGT,TGTG,T') && ($alleles eq '1/2/3') && ($info =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertions and a deletion (TG -> TGT, TG -> TGTG, TG -> T) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertions and a deletion (TG -> TGT, TG -> TGTG, TG -> T) ERROR!\n"; }
	}
	# DEL/HET or INS/HET
	elsif($position eq 45) {
		if(($reference eq 'G') && ($consensus eq 'C,GTG') && ($alleles eq '1/2')) { push @saved_results, "Heterozygous position and insertion (G -> C, G -> GTG) Checked"; }
		else { push @saved_results, "Heterozygous position and insertion (G -> C, G -> GTG) ERROR!\n"; }
	}
	elsif($position eq 46) {
		if(($reference eq 'AG') && ($consensus eq 'A,CG') && ($alleles eq '1/2')) { push @saved_results, "Heterozygous position and deletion (AG -> A, A -> C) Checked"; }
		else { push @saved_results, "Heterozygous position and deletion (AG -> A, A -> C) ERROR!\n"; }
	}
	
	elsif($position eq 47) {
		if(($reference eq 'G') && ($consensus eq 'C,GTG') && ($alleles eq '1/2')) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (G -> C, G -> GTG, G -> GTG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (G -> C, G -> GTG, G -> GTG) ERROR!\n"; }
	}
	elsif($position eq 48) {
		if(($reference eq 'T') && ($consensus eq 'C,TTG') && ($alleles eq '1/2')) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (T -> C, T -> C, T -> TTG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (T -> C, T -> C, T -> TTG) ERROR!\n"; }
	}
	elsif($position eq 49) {
		if(($reference eq 'GG') && ($consensus eq 'GGA,G,CG') && ($alleles eq '1/2')) { push @saved_results, "Tri-allelic heterozygous position, insertion and deletion (GG -> GGA, GG -> G, GG -> CG) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous position, insertion and deletion (GG -> GGA, GG -> G, GG -> CG) ERROR!\n"; }
	}
	
	else {
		push @saved_results, "False positive: $line ERROR!\n";
	}
}
close IN1;

# report test results
foreach(@saved_results) {
	print "$_\n";
}
