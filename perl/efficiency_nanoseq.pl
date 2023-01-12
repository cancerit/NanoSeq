#!/usr/bin/env perl

########## LICENCE ##########
# Copyright (c) 2022 Genome Research Ltd
# 
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
# 
# This file is part of NanoSeq.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
###########################

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Which;
use Capture::Tiny qw(capture);

# 1. Calculate num reads in the deduplicated bam
# 2. Calculate num reads merged bam
# 3. Calculate duplicate rate and number of bases sequenced
# 4. Get RB conformations - and save
# 5. Calculate efficiency with my own binomial thingy
my %opts = (
    't' => 1,
);

GetOptions('d|dedup=s'   => \my $deduplicated_bam,
           'x|duplex=s'  => \my $merged_bam,
           'o|out=s'     => \my $output_prefix,
           'r|ref=s'     => \my $ref_genome,
           'p|panel=s'   => \my $panel,
           't|threads=i' => \$opts{'t'},
           'h|help'      => \$opts{'h'}
) or pod2usage(2);
pod2usage(-verbose => 1, -exitval => 0) if(defined $opts{'h'});
pod2usage(2) if ( not defined $output_prefix and not defined $deduplicated_bam and not defined $merged_bam and not defined $ref_genome and not defined $panel);

die ("\nOutput prefix not defined\n") unless( $output_prefix );
die ("\nMust define the reference\n") unless ( $ref_genome);
die ("\nReference $ref_genome not found\n") unless ( -e $ref_genome );
die ("\nReference $ref_genome index not found\n") unless ( -e $ref_genome . ".fai" );
die ("\nMust define a deduplicated BAM\n") unless ( $deduplicated_bam );
die ("\nFile $deduplicated_bam not found\n") unless ( -e $deduplicated_bam );
my ($ext) = $deduplicated_bam  =~ /(\.[^.]+)$/;
$ext =~ s/.$/i/;
die ("\nIndex for $deduplicated_bam not found\n") unless ( -e "$deduplicated_bam".$ext );
die ("\nMust define a tumour duplex BAM\n") unless ( $merged_bam );
die ("\nFile $merged_bam not found\n") unless ( -e $merged_bam );
($ext) = $merged_bam  =~ /(\.[^.]+)$/;
$ext =~ s/.$/i/;
die ("\nIndex for $merged_bam $ext not found!\n") unless ( -e "$merged_bam".$ext );
die ("\nsamtools not found in path\n") unless ( which 'samtools' );
die ("\nRscript not found in path\n") unless ( which 'Rscript' );
die ("\nefficiency_nanoseq.R must be in path\n") unless ( which 'efficiency_nanoseq.R' );
my $do_panel;
if ($panel ) {
  die ("\nTargeted panel $panel not found\n") unless( -e $panel );
  $do_panel = 1;
} else {
  $do_panel = 0;
  $panel = "";
}

my $threads = $opts{'t'};
my $rb_output     = "$output_prefix.RBs";
my $main_output   = "$output_prefix.tsv";
# Get the first contig/chr:
open(REFI,"<$ref_genome.fai") or die("Couldn't open reference index\n");
my $region        = ( split(/\t/,<REFI>) )[0];
close(REFI);

##########################################################################################
# Calculating number of reads and duplicate rates
my($num_unique_reads,$num_sequenced_reads,$dup_rate,$on_near_frac);
my($num_sequenced_reads_on_near, $num_unique_reads_on_near );

sub read_count {
    my ($threads, $bam, $panel, $extra_opts) = @_;
    $extra_opts = q{} unless(defined $extra_opts);
    my $tgt_opt = $panel eq "" ? "" : "-L $panel";
    $tgt_opt   .= " -f 2 -F 2828";
    print STDOUT "Calculating number of reads in $bam...\n";
    my $cmd = sprintf "samtools view $tgt_opt $extra_opts -@ %d -c %s", $threads, $bam;
    my ($stdout, $stderr, $exit) = capture {
        system($cmd);
    };
    die "Error calling $cmd, $stderr\n" if ( $exit != 0 );
    chomp( $stdout);
    return $stdout; # reads counted
}

$num_sequenced_reads = &read_count($threads, $merged_bam, "");
print STDOUT "  Num sequenced reads=$num_sequenced_reads\n";
if ( $do_panel ) {
  $num_sequenced_reads_on_near = &read_count($threads, $merged_bam, $panel);
  $num_unique_reads_on_near = &read_count($threads, $deduplicated_bam, $panel);
  $dup_rate = ($num_sequenced_reads_on_near-$num_unique_reads_on_near)/$num_sequenced_reads_on_near;
  $on_near_frac = $num_sequenced_reads_on_near/$num_sequenced_reads;
  print STDOUT "  Num sequenced reads (on+near)=$num_sequenced_reads_on_near\n";
  print STDOUT "  Num unique reads(on+near)=$num_unique_reads_on_near\n";
  print STDOUT "  Duplicate rate=$dup_rate\n";
  print STDOUT "  On+near fraction=$on_near_frac\n";
} else {
  $num_unique_reads = &read_count($threads, $deduplicated_bam, "");
  $dup_rate = ($num_sequenced_reads-$num_unique_reads)/$num_sequenced_reads;

  print STDOUT "  Num unique reads=$num_unique_reads\n";
  print STDOUT "  Duplicate rate=$dup_rate\n";
}


##########################################################################################
# Get read bundle comformations
my $bam = $merged_bam;
sub count_RB {
  my %rbs;
  my ($threads, $flag, $bam, $region, $panel) = @_;
  my $tgt_opt = $panel eq "" ? "" : "-L $panel";
  $tgt_opt   .= " -F 2828";
  my $cmd = "samtools view $tgt_opt -@ $threads -f $flag $bam $region";
  open(IN, "$cmd |") || die "Error launching $cmd\n"; 
  while(<IN>) {
    chomp;
    my @tmp = (split(/\t/,$_));
    my $rb;
    foreach my $t ( @tmp ) {
      if ( $t =~ /^RB:Z/ ) {
        $rb = $t;
        last;
      }
    }
    next unless(defined($rb));
    $rbs{$rb}++;
  }
  close(IN) or die ("error when calling $cmd : $?, $!\n");
  return( %rbs );
}
# Get first reads in reverse:
print STDOUT "RB comformation: 1st reads in reverse...\n";
my %rbf2r1 = &count_RB($threads, 82, $bam, $region, $panel);

print STDOUT "RB comformation: 2nd reads in reverse...\n";
my %rbf1r2 = &count_RB($threads, 146, $bam, $region, $panel);

foreach my $ikey (keys %rbf2r1) {
  $rbf1r2{$ikey} = 0 if ( ! exists( $rbf1r2{$ikey}));
}

foreach my $ikey (keys %rbf1r2) {
  $rbf2r1{$ikey} = 0 if ( ! exists( $rbf2r1{$ikey}));
}

open(OUT,">$rb_output") || die "Error openning output $rb_output\n";
foreach my $rb ( keys %rbf1r2 ) {
  print OUT "$rb\t",$rbf1r2{$rb},"\t",$rbf2r1{$rb},"\n";
}
close(OUT);

##########################################################################################
# Call R to get the two values that inform on strand misses:
my($reads_per_rb,$f_eff,$zib_eff,$ok_rbs,$total_rbs,$gc_both,$gc_single,$total_reads);
my $cmd = "efficiency_nanoseq.R $rb_output $ref_genome ";
print STDOUT "Running: $cmd\n";
my ($stdout, $stderr, $exit) = capture {
  system($cmd);
};
die "Error calling $cmd, $stderr\n" if ( $exit != 0 );
foreach (split('\n',$stdout)) {
  print STDOUT "  Routput: ",$_,"\n";
  chomp;
  if(/READS_PER_RB/) {
    $reads_per_rb = (split)[1];
  } elsif(/F-EFF/) {
    $f_eff = (split)[1];
  } elsif(/ZIB-EFF/) {
    $zib_eff = (split)[1];
  } elsif(/OK_RBS/) {
    $ok_rbs = (split)[1];
  } elsif(/TOTAL_RBS/) {
    $total_rbs = (split)[1];
  } elsif(/GC_BOTH/) {
    $gc_both = (split)[1];
  } elsif(/GC_SINGLE/) {
    $gc_single = (split)[1];
  } elsif(/TOTAL_READS/) {
    $total_reads = (split)[1];
  }
}
# Check it run properly:
open(OUT, ">$main_output") || die "Error writing to $main_output\n";
print OUT "# Whole-genome metrics:\n";
print OUT "NUM_READS_SEQUENCED\t$num_sequenced_reads\n";
if ($do_panel) {
  print OUT "ON+NEAR_FRACTION\t$on_near_frac\n";
  print OUT "ON+NEAR_UNIQUE_READS\t$num_unique_reads_on_near\n";
  print OUT "ON+NEAR_SEQUENCED_READS\t$num_sequenced_reads_on_near\n";
} else {
  print OUT "NUM_UNIQUE_READS\t$num_unique_reads\n";
}
print OUT "DUPLICATE_RATE\t$dup_rate\n";

print OUT "# RB metrics are reported for chr/contig $region only:\n";

my $bases_sequenced = $total_reads * 150;
my $bases_ok_rbs    = $ok_rbs * (300); # assuming mates don't overlap. Removing 50 bps for varios trimmings (rough estimate)
print OUT "TOTAL_RBS\t$total_rbs\n";
print OUT "TOTAL_READS_IN_RBS\t$total_reads\n";
print OUT "OK_RBS(2+2)\t$ok_rbs\n";
print OUT "READS_PER_RB\t$reads_per_rb\n";
print OUT "F-EFF\t$f_eff\n";
print OUT "EFFICIENCY\t",$bases_ok_rbs / $bases_sequenced,"\n";
print OUT "GC_BOTH\t$gc_both\n";
print OUT "GC_SINGLE\t$gc_single\n";
close(OUT);

__END__
=head1 NAME

efficiency_nanoseq.pl

=head1 SYNOPSIS

efficiency_nanoseq.pl [-h] [-t n] [-p panel ] -dedup BAM -duplex BAM -o prefix -r reference

    -dedup             -d   Deduplicated BAM
    -duplex            -x   Duplex BAM
    -ref               -r   reference file
    -out               -o   Ouptup prefix
    -panel             -p   Target panel (bed)

  Optional parameters:
    -threads           -t   Threads (1)
    -help              -h

=head1 OPTIONS

=over 8

=item B<-dedup>

Deduplicated BAM (neat) and index. The BAM after being processed by randomreadinbundle.

=item B<-duplex>

Duplex BAM and index. BAM prior to being processed by randomreadinbundle.

=item B<-out>

Output prefix

=item B<-ref>

Reference file and index

=item B<-panel>

Target panel (bed format)

=item B<-threads>

Use samtools with these many threads

=back

=cut

