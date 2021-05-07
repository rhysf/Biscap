#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF_lines;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -m <VCF outputfrom test.pileup>\n";
our($opt_m);
getopt('m');
die $usage unless ($opt_m);  
die "Cannot open $opt_m : $!" unless (-e $opt_m);

my (@FP_lines, @FN_lines, @FN_lines_no_cov);
my $location_of_genotype_in_format = 0;
my @saved_results;

# Go through VCF and tally known sites
open my $fh, '<', $opt_m or die "Cannot open $opt_m : $!\n";
warn "Cataloging variants from $opt_m...\n";
VCF: while (my $line = <$fh>) {
   	chomp $line;
	my $VCF_line = vcflines::read_VCF_lines($line);
	next VCF if($$VCF_line{'next'} eq 1);
   	
	### Check calls
	
	# Homozygous SNPs
	if($$VCF_line{'position'} eq 6) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'A') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous SNP (T -> A) Checked"; }
		else { push @saved_results, "SNP (A) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 7) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'C') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous SNP (G -> C) Checked"; }
		else { push @saved_results, "SNP (C) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 8) {
		if(($$VCF_line{'reference_VCF_format'} eq 'C') && ($$VCF_line{'consensus_VCF_format'} eq 'G') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous SNP (C -> G) Checked"; }
		else { push @saved_results, "SNP (G) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 9) {
		if(($$VCF_line{'reference_VCF_format'} eq 'A') && ($$VCF_line{'consensus_VCF_format'} eq 'T') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous SNP (A -> T) Checked\n"; }
		else { push @saved_results, "SNP (T) ERROR!\n"; }
	}
	
	# DELS
	elsif($$VCF_line{'position'} eq 11) {
		if(($$VCF_line{'reference_VCF_format'} eq 'CA') && ($$VCF_line{'consensus_VCF_format'} eq 'C') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous deletion (CA -> A) Checked"; }
		else { push @saved_results, "Homozygous deletion (CA -> A) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 12) {
		if(($$VCF_line{'reference_VCF_format'} eq 'AA') && ($$VCF_line{'consensus_VCF_format'} eq 'A') && ($$VCF_line{'GT0'} eq '0/1')) { push @saved_results, "Bi-allelic heterozygous deletion with ref. (AA -> A) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous deletion with ref. (AA -> A) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 13) {
		if(($$VCF_line{'reference_VCF_format'} eq 'ATG') && ($$VCF_line{'consensus_VCF_format'} eq 'AG,A') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Bi-allelic heterozygous deletion (ATG -> A, ATG -> A-G) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous deletion (ATG -> A, ATG -> A-G) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 14) {
		if(($$VCF_line{'reference_VCF_format'} eq 'TG') && ($$VCF_line{'consensus_VCF_format'} eq 'T') && ($$VCF_line{'GT0'} eq '0/1') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion with ref. (TG -> T, TG -> T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion with ref. (TG -> T, TG -> T) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 15) {
		if(($$VCF_line{'reference_VCF_format'} eq 'GTA') && ($$VCF_line{'consensus_VCF_format'} eq 'GA,G') && ($$VCF_line{'GT0'} eq '1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion (GTA -> G, GTA -> G-A, GTA -> G-A) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous deletion (GTA -> G, GTA -> G-A, GTA -> G-A) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 16) {
		if(($$VCF_line{'reference_VCF_format'} eq 'TAA') && ($$VCF_line{'consensus_VCF_format'} eq 'TA,T') && ($$VCF_line{'GT0'} eq '0/1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous deletion with ref. (TAA -> T, TAA -> T-A) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous deletion with ref. (TAA -> T, TAA -> T-A) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 17) {
		if(($$VCF_line{'reference_VCF_format'} eq 'AAGG') && ($$VCF_line{'consensus_VCF_format'} eq 'AGG,AG,A') && ($$VCF_line{'GT0'} eq '1/2/3') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous deletion (AAGG -> A, AAGG -> A--G, AAGG -> A-GG) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous deletion (AAGG -> A, AAGG -> A--G, AAGG -> A-GG) ERROR!\n"; }
	}
	
	# Heterozygous positions
	elsif($$VCF_line{'position'} eq 20) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'C') && ($$VCF_line{'GT0'} eq '0/1')) { push @saved_results, "Heterozygous position with ref. (G -> G/C) Checked"; }
		else { push @saved_results, "Heterozygous position with ref. (G -> G/C) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 21) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && (($$VCF_line{'consensus_VCF_format'} eq 'C,G') || ($$VCF_line{'consensus_VCF_format'} eq 'G,C')) && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Heterozygous position (T -> G/C) Checked"; }
		else { push @saved_results, "Heterozygous position (T -> G/C) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 22) {
		if(($$VCF_line{'reference_VCF_format'} eq 'C') && ($$VCF_line{'consensus_VCF_format'} eq 'G') && ($$VCF_line{'GT0'} eq '0/1') && ($$VCF_line{'info'} =~ m/TRIPROB:0\/1\/1/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position with ref. (C -> C/G/G) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position with ref. (C -> C/G/G) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 23) {
		if(($$VCF_line{'reference_VCF_format'} eq 'C') && ($$VCF_line{'consensus_VCF_format'} eq 'T') && ($$VCF_line{'GT0'} eq '0/1') && ($$VCF_line{'info'} =~ m/TRIPROB:0\/0\/1/)) { push @saved_results, "Bi-allelic (tri. prob.) position with ref. (C -> C/C/T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) position with ref. (C -> C/C/T) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 24) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'C,T') && ($$VCF_line{'GT0'} eq '1/2') && ($$VCF_line{'info'} =~ m/TRIPROB:1\/1\/2/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position (G -> C/C/T) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position (G -> C/C/T) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 25) {
		if(($$VCF_line{'reference_VCF_format'} eq 'C') && (($$VCF_line{'consensus_VCF_format'} eq 'T,G') || ($$VCF_line{'consensus_VCF_format'} eq 'G,T')) && ($$VCF_line{'GT0'} eq '0/1/2')) { push @saved_results, "Tri-allelic heterozygous position with ref. (C -> C/T/G) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous position with ref. (C -> C/T/G) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 26) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && (($$VCF_line{'consensus_VCF_format'} eq 'A,C,T') || ($$VCF_line{'consensus_VCF_format'} eq 'A,T,C') || ($$VCF_line{'consensus_VCF_format'} eq 'C,A,T') || ($$VCF_line{'consensus_VCF_format'} eq 'T,A,C')) && ($$VCF_line{'GT0'} eq '1/2/3')) { push @saved_results, "Tri-allelic heterozygous position (G -> A/T/C) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous position (G -> A/T/C) ERROR!\n"; }
	}

	# INS
	elsif($$VCF_line{'position'} eq 30) {
		if(($$VCF_line{'reference_VCF_format'} eq 'A') && ($$VCF_line{'consensus_VCF_format'} eq 'AA') && ($$VCF_line{'GT0'} eq 1)) { push @saved_results, "Homozygous insertion (A -> AA) Checked"; }
		else { push @saved_results, "Homozygous insertion (A -> AA) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 31) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'TA') && ($$VCF_line{'GT0'} eq '0/1')) { push @saved_results, "Bi-allelic heterozygous insertion with ref. (T -> TA) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous insertion with ref. (T -> TA) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 32) {
		if(($$VCF_line{'reference_VCF_format'} eq 'C') && ($$VCF_line{'consensus_VCF_format'} eq 'CT,CTG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Bi-allelic heterozygous insertion (C -> CT, C -> CTG) Checked"; }
		else { push @saved_results, "Bi-allelic heterozygous insertion (C -> CT, C -> CTG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 33) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'TG') && ($$VCF_line{'GT0'} eq '0/1') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion with ref. (T -> TG, T -> TG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion with ref. (T -> TG, T -> TG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 34) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'GT,GTA') && ($$VCF_line{'GT0'} eq '1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion (G -> GT, G -> GTA) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion (G -> GT, G -> GTA) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 35) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'TA,TAA') && ($$VCF_line{'GT0'} eq '0/1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion with ref. (T -> TA, T -> TAA) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion with ref. (T -> TA, T -> TAA) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 36) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'TA,TAG,TAGG') && ($$VCF_line{'GT0'} eq '1/2/3') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion (T -> TA, T -> TAG, T -> TAGG) Checked\n"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion (T -> TA, T -> TAG, T -> TAGG) ERROR!\n"; }
	}
	
	###  MIXES
	# INDELS
	elsif($$VCF_line{'position'} eq 41) {
		if(($$VCF_line{'reference_VCF_format'} eq 'TG') && ($$VCF_line{'consensus_VCF_format'} eq 'TGTG,T') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Heterozygous insertion and deletion (TG -> T, TG -> TGTG) Checked"; }
		else { push @saved_results, "Heterozygous insertion and deletion (TG -> T, TG -> TGTG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 42) {
		if(($$VCF_line{'reference_VCF_format'} eq 'GT') && ($$VCF_line{'consensus_VCF_format'} eq 'GTTG,G') && ($$VCF_line{'GT0'} eq '1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion and deletion (GT -> GTTG, GT -> GTTG, GT -> G) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous insertion and deletion (GT -> GTTG, GT -> GTTG, GT -> G) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 43) {
		if(($$VCF_line{'reference_VCF_format'} eq 'TT') && ($$VCF_line{'consensus_VCF_format'} eq 'TTTG,T') && ($$VCF_line{'GT0'} eq '0/1/2') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertion, deletion and ref. (TT -> TTTG, TT -> T) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertion, deletion and ref. (TT -> TTTG, TT -> T) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 44) {
		if(($$VCF_line{'reference_VCF_format'} eq 'TG') && ($$VCF_line{'consensus_VCF_format'} eq 'TGT,TGTG,T') && ($$VCF_line{'GT0'} eq '1/2/3') && ($$VCF_line{'info'} =~ m/TRI/)) { push @saved_results, "Tri-allelic heterozygous insertions and a deletion (TG -> TGT, TG -> TGTG, TG -> T) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous insertions and a deletion (TG -> TGT, TG -> TGTG, TG -> T) ERROR!\n"; }
	}
	# DEL/HET or INS/HET
	elsif($$VCF_line{'position'} eq 45) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'C,GTG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Heterozygous position and insertion (G -> C, G -> GTG) Checked"; }
		else { push @saved_results, "Heterozygous position and insertion (G -> C, G -> GTG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 46) {
		if(($$VCF_line{'reference_VCF_format'} eq 'AG') && ($$VCF_line{'consensus_VCF_format'} eq 'A,CG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Heterozygous position and deletion (AG -> A, A -> C) Checked"; }
		else { push @saved_results, "Heterozygous position and deletion (AG -> A, A -> C) ERROR!\n"; }
	}
	
	elsif($$VCF_line{'position'} eq 47) {
		if(($$VCF_line{'reference_VCF_format'} eq 'G') && ($$VCF_line{'consensus_VCF_format'} eq 'C,GTG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (G -> C, G -> GTG, G -> GTG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (G -> C, G -> GTG, G -> GTG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 48) {
		if(($$VCF_line{'reference_VCF_format'} eq 'T') && ($$VCF_line{'consensus_VCF_format'} eq 'C,TTG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (T -> C, T -> C, T -> TTG) Checked"; }
		else { push @saved_results, "Bi-allelic (tri. prob.) heterozygous position and insertion (T -> C, T -> C, T -> TTG) ERROR!\n"; }
	}
	elsif($$VCF_line{'position'} eq 49) {
		if(($$VCF_line{'reference_VCF_format'} eq 'GG') && ($$VCF_line{'consensus_VCF_format'} eq 'GGA,G,CG') && ($$VCF_line{'GT0'} eq '1/2')) { push @saved_results, "Tri-allelic heterozygous position, insertion and deletion (GG -> GGA, GG -> G, GG -> CG) Checked"; }
		else { push @saved_results, "Tri-allelic heterozygous position, insertion and deletion (GG -> GGA, GG -> G, GG -> CG) ERROR!\n"; }
	}
	
	else {
		push @saved_results, "False positive: $line ERROR!\n";
	}
}
close $fh;

# report test results
foreach(@saved_results) {
	print "$_\n";
}
