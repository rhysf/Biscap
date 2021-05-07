package ascii;
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

sub ascii {
	my %symbol_to_ascii;

	$symbol_to_ascii{NUL} = 0;
	$symbol_to_ascii{SOH} = 1;
	$symbol_to_ascii{STX} = 2;
	$symbol_to_ascii{ETX} = 3;
	$symbol_to_ascii{EOT} = 4;
	$symbol_to_ascii{ENQ} = 5;
	$symbol_to_ascii{ACK} = 6;
	$symbol_to_ascii{BEL} = 7;
	$symbol_to_ascii{BS} = 8;
	$symbol_to_ascii{TAB} = 9;
	$symbol_to_ascii{LF} = 10;
	$symbol_to_ascii{VT} = 11;
	$symbol_to_ascii{FF} = 12;
	$symbol_to_ascii{CR} = 13;
	$symbol_to_ascii{SO} = 14;
	$symbol_to_ascii{SI} = 15;

	$symbol_to_ascii{DLE} = 16;
	$symbol_to_ascii{DC1} = 17;
	$symbol_to_ascii{DC2} = 18;
	$symbol_to_ascii{DC3} = 19;
	$symbol_to_ascii{DC4} = 20;
	$symbol_to_ascii{NAK} = 21;
	$symbol_to_ascii{SYN} = 22;
	$symbol_to_ascii{ETB} = 23;
	$symbol_to_ascii{CAN} = 24;
	$symbol_to_ascii{EM} = 25;
	$symbol_to_ascii{SUB} = 26;
	$symbol_to_ascii{ESC} = 27;
	$symbol_to_ascii{FS} = 28;
	$symbol_to_ascii{GS} = 29;
	$symbol_to_ascii{RS} = 30;
	$symbol_to_ascii{US} = 31;

	$symbol_to_ascii{'\s'} = 32;
	$symbol_to_ascii{'!'} = 33;
	$symbol_to_ascii{'"'} = 34;
	$symbol_to_ascii{'#'} = 35;
	$symbol_to_ascii{'$'} = 36;
	$symbol_to_ascii{'%'} = 37;
	$symbol_to_ascii{'&'} = 38;
	$symbol_to_ascii{"'"} = 39;
	$symbol_to_ascii{'('} = 40;
	$symbol_to_ascii{')'} = 41;
	$symbol_to_ascii{'*'} = 42;
	$symbol_to_ascii{'+'} = 43;
	$symbol_to_ascii{','} = 44;
	$symbol_to_ascii{'-'} = 45;
	$symbol_to_ascii{'.'} = 46;
	$symbol_to_ascii{'/'} = 47;

	$symbol_to_ascii{'0'} = 48;
	$symbol_to_ascii{'1'} = 49;
	$symbol_to_ascii{'2'} = 50;
	$symbol_to_ascii{'3'} = 51;
	$symbol_to_ascii{'4'} = 52;
	$symbol_to_ascii{'5'} = 53;
	$symbol_to_ascii{'6'} = 54;
	$symbol_to_ascii{'7'} = 55;
	$symbol_to_ascii{'8'} = 56;
	$symbol_to_ascii{'9'} = 57;
	$symbol_to_ascii{':'} = 58;
	$symbol_to_ascii{';'} = 59;
	$symbol_to_ascii{'<'} = 60;
	$symbol_to_ascii{'='} = 61;
	$symbol_to_ascii{'>'} = 62;
	$symbol_to_ascii{'?'} = 63;

	$symbol_to_ascii{'@'} = 64;
	$symbol_to_ascii{A} = 65;
	$symbol_to_ascii{B} = 66;
	$symbol_to_ascii{C} = 67;
	$symbol_to_ascii{D} = 68;
	$symbol_to_ascii{E} = 69;
	$symbol_to_ascii{F} = 70;
	$symbol_to_ascii{G} = 71;
	$symbol_to_ascii{H} = 72;
	$symbol_to_ascii{I} = 73;
	$symbol_to_ascii{J} = 74;
	$symbol_to_ascii{K} = 75;
	$symbol_to_ascii{L} = 76;
	$symbol_to_ascii{M} = 77;
	$symbol_to_ascii{N} = 78;
	$symbol_to_ascii{O} = 79;

	$symbol_to_ascii{P} = 80;
	$symbol_to_ascii{Q} = 81;
	$symbol_to_ascii{R} = 82;
	$symbol_to_ascii{S} = 83;
	$symbol_to_ascii{T} = 84;
	$symbol_to_ascii{U} = 85;
	$symbol_to_ascii{V} = 86;
	$symbol_to_ascii{W} = 87;
	$symbol_to_ascii{X} = 88;
	$symbol_to_ascii{Y} = 89;
	$symbol_to_ascii{Z} = 90;
	$symbol_to_ascii{"["} = 91;
	$symbol_to_ascii{"\\"} = 92;
	$symbol_to_ascii{"]"} = 93;
	$symbol_to_ascii{'^'} = 94;
	$symbol_to_ascii{'_'} = 95;

	$symbol_to_ascii{'`'} = 96;
	$symbol_to_ascii{a} = 97;
	$symbol_to_ascii{b} = 98;
	$symbol_to_ascii{c} = 99;
	$symbol_to_ascii{d} = 100;
	$symbol_to_ascii{e} = 101;
	$symbol_to_ascii{f} = 102;
	$symbol_to_ascii{g} = 103;
	$symbol_to_ascii{h} = 104;
	$symbol_to_ascii{i} = 105;
	$symbol_to_ascii{j} = 106;
	$symbol_to_ascii{k} = 107;
	$symbol_to_ascii{l} = 108;
	$symbol_to_ascii{m} = 109;
	$symbol_to_ascii{n} = 110;
	$symbol_to_ascii{o} = 111;

	$symbol_to_ascii{p} = 112;
	$symbol_to_ascii{q} = 113;
	$symbol_to_ascii{r} = 114;
	$symbol_to_ascii{s} = 115;
	$symbol_to_ascii{t} = 116;
	$symbol_to_ascii{u} = 117;
	$symbol_to_ascii{v} = 118;
	$symbol_to_ascii{w} = 119;
	$symbol_to_ascii{x} = 120;
	$symbol_to_ascii{y} = 121;
	$symbol_to_ascii{z} = 122;
	$symbol_to_ascii{'{'} = 123;
	$symbol_to_ascii{'|'} = 124;
	$symbol_to_ascii{'}'} = 125;
	$symbol_to_ascii{'~'} = 126;
	$symbol_to_ascii{DEL} = 127;           

	return %symbol_to_ascii;
}

1;
