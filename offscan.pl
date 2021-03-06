#!/usr/bin/perl -w

=head1 ABOUT

  offscan.pl - A perl script used for large scale profiling possible
               target sites of CRISPR/CAS9 system by using seqmap

=head1 VERSION

  Version: 1.5.1
  Jul 24, 2017

=head1 SYNOPSIS

  $ perl ./offscan.pl <option>
  For usage, type './offscan.pl -h'
  To test a demo, type './offscan.pl -t'

=head1 DESCRIPTION

  This script takes a FASTA file of the target sites and queries the
  reference database by seqmap tolerant designated mismatches, and gives
  the human readable and computer mineable target/off target sites and
  flanking sequencing, also the summary file.
  
  Seqmap is describled:
  Jiang, H. and Wong, WH. (2008) SeqMap: mapping aassive amount of
  oligonucleotides to the genome, Bioinformatics, 24(20):2395-2396.
  And available at http://www-personal.umich.edu/~jianghui/seqmap/
  
  This script is firstly used in the following research:
  Guo, XG., Zhang, TJ., Hu, Z., Zhang, YQ., Shi ZY., Wang, QH., Cui, Y.,
  Wang, FQ., Zhao, H. and Chen, YL. (2014) Efficient RNA/Cas9-mediated
  genome editing in Xenopus tropicalis. Development, 141(3):707-714.

=head1 NOTE

  1. Max number of mismatch (option -n) is 5.
  2. The type memory size (option -m) is a number. Unit: megabyte.
  3. User can download the latest seqmap and replace the one included in
     directory 'bin'.
 
=head1 AUTHOR

  Wang, Qinhu
  Northwest A&F University
  E-mail: wangqinhu@nwafu.edu.cn
  
=head1 LICENSE
 
  For non-commercial use only.
  
=cut


use strict;
use Getopt::Std;


#-------------------------------------------------------------------------------
# Get Options
#-------------------------------------------------------------------------------
my %options;
getopts("g:q:r:d:l:n:o:s:m:f:a:b:th" , \%options);

my $grd = $options{g};                      # Query RNA in sgRNA Designer format
my $rna = $options{q};                      # Query RNA in FASTA format
my $rsz = $options{r} || '23';              # RNA size
my $ref = $options{d};                      # Database for searching, reference file
my $len = $options{l} || '500';             # Length of flanking sequence to extrat (for primer design)
my $mis = $options{n} || '5';               # Number of mismatch allowed
my $rpt = $options{o} || 'cas9off.xls';     # Output report file, tab-delimited file
my $smf = $options{s} || 'sum.xls';         # Summary file shows the statics result of the off-targeted sites
my $mem = $options{m} || '2048';            # Memory allocated, required by seqmap
my $fmt = $options{f} || "true";            # Format reference fasta
my $app = $options{a} || "";                # Path of seqmap
my $map = $options{b} || "seqmap.out";      # Map file, seqmap output


# Usage
&usage() if $options{h};

#-------------------------------------------------------------------------------
# Run demo
#-------------------------------------------------------------------------------
if ($options{t}) {
	$rna = 'demo/target.fa';
	$ref = 'demo/reference.fa';
	$rpt = 'demo.cas9off.xls';
	$smf = 'demo.sum.xls';
	$fmt = 'false';
}

# Check inputs
if (defined($grd) && defined($rna)) {
	warn "Aborted: Both [sgRNA Designer format, -g] and [fasta format, -q] were found.\nPlease use one of them!\n";
	&usage();
}

if (defined($grd)) {
	$rna = sgRNA2fasta($grd);
}

unless (defined($rna) && defined($ref)) {
	&usage();
}

# Judge OS
# Support Mac, Linux (32 and 64 bit) and Windows
my @osp = ("./bin/seqmap", "./bin/seqmap-1.0.12-linux", "./bin/seqmap-1.0.12-linux-64", "./bin/seqmap-1.0.12-windows.exe");
if ($app eq '') {
	$app = $osp[&ios()];
}

# Format DB
&format_seqid($options{d}) if $fmt eq "true";

# Check NGG
open (NGG, $rna) or die "Cannot open file $rna: $!\n";
while (my $line = <NGG>) {
	next if $line =~ /^\>/;
	$line =~ s/\s+//g;
	$line =~ s/\d+//g;
	unless ($line =~ /GG$/i) {
		print "$line\n";
		print "In RNA file $rna, line $.: No NGG found!\n";
		exit;
	}
	if (length($line) != $rsz) {
		print "In RNA file $rna, line $.: RNA length is not $rsz!\n" ;
		exit;
	}
}
close NGG;

#-------------------------------------------------------------------------------
#  Main
#-------------------------------------------------------------------------------

print "Mapping ...\n\n";

# Mapping
my $qlen = $rsz - 3;
system("$app $mis $rna $ref $map /cut:1,$qlen /output_all_matches /available_memory:$mem /silent");

print "\nParse map ...\n\n";

# Read genome
my %ref = ();
my $seqid = undef;
open (SEQ, $ref);
while (my $line = <SEQ>) {
	if ($line =~ /^>(\S+)/) {
		$seqid = $1;
	} else {
		$line =~ s/\s+//g;
		$line =~ s/\d+//g;
		$ref{$seqid} .= $line;
	}
}
close SEQ;

# Read map
my @map = ();
my $i = 0;
my @tab = ();
open (MAP, "$map") or die "Cannot open file $map: $!\n";
while (<MAP>) {
	next if /^trans_id/;
	next if /^\s*$/;
	my ($refid, $refcor,$refaln,$rnaid,$rnaseq,$numis,$str) = split /\s+/, $_;
	# Modify alignment and extract flanking sequences
	my $maln = undef;
	my $flk = undef;
	if ($str eq "+") {
		# If mstar <= 0
		my $mstar = $refcor - 1;
		$mstar = 0 if $mstar <= 0;
		# End if
	
		# get maln
		$maln = substr($ref{$refid}, $mstar, $rsz);
		
		# If rstar <= 0
		my $rstar = $refcor-$len - 1;
		$rstar = 0 if $rstar <= 0;
		# End if
		
		# get flk
		$flk = substr($ref{$refid}, $rstar, 2*$len+length($rnaseq) + 3);
	}
	if ($str eq "-") {
		# Addtional 3 bases inclued since we used 1-20 (but not 21-23) for search.
		$refcor = $refcor - 3;
		
		# If mstar <= 0
		my $mstar = $refcor - 1;
		$mstar = 0 if $mstar <= 0;
		# End if
		
		# get maln
		$maln = substr($ref{$refid}, $mstar, $rsz);
		$maln = &revcom_seq($maln);
		
		# If rstart <= 0
		my $rstar = $refcor-$len - 1;
		$rstar = 0 if $rstar <= 0;
		# End if
		
		#get flk
		$flk = substr($ref{$refid}, $rstar, 2*$len+length($rnaseq) + 3);	
		$flk = &revcom_seq($flk);
	}
	# Filter sequence without NGG
	next unless $maln =~ /GG$/i;
	$map[$i++] = {	"refid"  => $refid,
			"refcor" => $refcor,
			"refaln" => uc($maln),
			"rnaid"  => $rnaid,
			"rnaseq" => uc($rnaseq . "NGG"),
			"numis"  => $numis,
			"str"    => $str,
			"flk"    => uc($flk)
	};
	push @tab, $rnaid;
}
close MAP;

# Sort map
@map = sort by_map @map;

# Output map
open (OUT, ">$rpt") or die "Cannot open file $rpt: $!\n";
print OUT "Target ID\tTarget Site\tReference ID\tReference Alignment\tNo. of Mismatch\tStrand\tMatch Start\tFlanking Sequence\n";
foreach my $map (@map) {
	my $diff_aln = mark_diff_aln($map->{"refaln"}, $map->{"rnaseq"});
	print OUT $map->{"rnaid"}, "\t", $map->{"rnaseq"}, "\t", $map->{"refid"}, "\t",
		  $diff_aln, "\t", $map->{"numis"}, "\t", $map->{"str"}, "\t",
		  $map->{"refcor"}, "\t", $map->{"flk"}, "\n";
}
close OUT;

my %tab=();
foreach (@tab) {
	$tab{$_} = 1;
}

# Summarize
my %sum = ();

foreach my $rid (keys %tab) {
	foreach my $i (0..$mis) {
		$sum{$rid}{$i} = 0;
		foreach my $map (@map) {
			if ($map->{"rnaid"} eq $rid && $map->{"numis"} eq $i) {
				$sum{$rid}{$i}++;
			}
		}
	}
}

# Output summary
open (SUM, ">$smf");
for my $i (0..$mis) {
	print SUM "\t$i";
}
print SUM "\n";
foreach my $rid (sort keys %tab) {
	print SUM $rid, "\t";
	foreach my $i (0..$mis) {
		print SUM $sum{$rid}{$i}, "\t";
	}
	print SUM "\n";
}
close SUM;

unlink("$map");

print "Done!\n\n";

#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
sub by_map {
	$a->{"rnaid"} cmp $b->{"rnaid"}
	or $a->{"refid"} cmp $b->{"refid"}
	or $a->{"refcor"} cmp $b->{"refcor"}
	or $a->{"numis"} cmp $b->{"numis"}
	or $a->{"str"} cmp $b->{"str"}
}

sub revcom_seq {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ATCGNXatcgnx/TAGCNXtagcnx/;
	return $seq;
}

sub format_seqid {
	die "Database file required!\n" unless $options{d};
	my ($in) = @_;
	my $out = 'refseq.fa';
	print "Formatting DB ...\n";
	open (IN, $in) or die "Cannot open file $in: $!\n";
	open (OUT,">$out" ) or die "Cannot open file $out: $!\n";
	while (<IN>) {
		if (/^\>(\S+)/) {
			my $sid = $1;
			print OUT ">$sid\n";
		} else {
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	$ref = 'refseq.fa';
}

sub ios {
	my $os = $^O;
	my $osid = undef;
	if($os  eq  "darwin") {
		$osid = 0;
	} elsif($os  eq  "linux") {
		my $bit=`getconf LONG_BIT`;
		chomp $bit;
		if ($bit eq "32") {
			$osid = 1;
		} elsif ($bit eq "64")  {
			$osid = 2;		
		} else {
			die "Unknown CPU!\n";
		}
	} elsif($os  eq  "MSWin32") {
		$osid = 3;
	} else {
		die "Unknown OS detected by cas9off!\n";
	}
	return $osid;
}

sub sgRNA2fasta {
	my $in_file = shift;
	my $out_file = "sgRNA.fa";
	my $head = 4;
	my $tail = 3;
	my $gRNA = 23;
	my %seq = ();
	my $id = 'seq00000';
	open (IN, $in_file) or die "Cannot open $in_file : $!\n";
	open (OUT, ">$out_file") or die "Cannot open $out_file : $!\n";
	while (<IN>) {
		chomp;
		next if /^\s*$/;
		next if /^Spacer/;
		if (/\S{$head}(\S{$gRNA})\S{$tail}/) {
			$id++;
			print OUT ">$id\n$1\n";
		} else {
			print STDERR "Unknown gRNA found in $in_file!\n";
		}
	}
	close OUT;
	close IN;
	return $out_file;
}

sub mark_diff_aln {
	my ($aln, $rna) = @_;
	my $len = $rsz - 3;
	my @aln = split //, $aln;
	my @rna = split //, $rna;
	for (my $i = 0; $i < $len; $i++) {
		if ($aln[$i] ne $rna[$i]) {
			$aln[$i] = lc($aln[$i]);
		}
	}
	my $diff_aln = join('', @aln);
	return $diff_aln;
}

sub usage
{
	print <<USAGE;

cas9off version 1.5.1

Usage:

    offscan.pl <option>
    
    -h  Help
    -t  Test demo
    -g  Query RNA file, in sgRNA desinger output format
    -q  Query RNA file, in FASTA format
    -r  Query RNA length, default is 23 bp, including NGG
    -d  Database file, in FASTA format, contains reference sequences
    -l  Length of flanking sequences to extract, default is 500 bp
    -n  Number of mismatch (0-5), does not include NGG
    -o  Output filename, default is cas9off.xls
    -s  Summary filename, default is sum.xls
    -m  Memory size, default is 2048 [Mb (2Gb)]
    -f  Format reference sequence or not, can be 'true' [default] or 'false',
        critical if your FASTA file has annotation
    -a  Path of seqmap
    -b  Map file, seqmap output

    For citation:
	
    * Jiang, H. and Wong, WH. (2008) SeqMap: mapping aassive amount of 
      oligonucleotides to the genome, Bioinformatics, 24(20):2395-2396.
	
    * Guo, XG., Zhang, TJ., Hu, Z., Zhang, YQ., Shi ZY., Wang, QH., Cui, Y.,
      Wang, FQ., Zhao, H. and Chen, YL. (2014) Efficient RNA/Cas9-mediated
      genome editing in Xenopus tropicalis. Development, 141(3):707-714.

USAGE
	exit;
}


__END__
