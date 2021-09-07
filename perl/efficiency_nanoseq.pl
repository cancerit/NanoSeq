#!/usr/bin/perl -w

########## LICENCE ##########
# Copyright (c) 2020 Genome Research Ltd.
#
# Author: Cancer Genome Project <cgphelp@sanger.ac.uk>
#
# This file is part of NanoSeq.
#
# NanoSeq is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Which;
use Capture::Tiny qw(capture);
my $VERSION="1.0.0";

# 1. Calculate num reads neat bam
# 2. Calculate num reads merged bam
# 3. Calculate duplicate rate and number of bases sequenced
# 4. Get RB conformations - and save
# 5. Calculate efficiency with my own binomial thingy
my %opts;

GetOptions('n|normal=s'  => \my $neat_bam,
           'd|tumour=s'  => \my $merged_bam,
           'o|out=s'     => \my $output_prefix,
           'r|ref=s'     => \my $ref_genome,
           't|threads=i' => \$opts{'t'},
           'h|help'      => \$opts{'h'},
           'v|version'   => \$opts{'v'},
) or pod2usage(2);
pod2usage(-verbose => 1) if(defined $opts{'h'});
if(defined $opts{'v'}) {
  print sprintf "VERSION: %s\n", $VERSION;
  exit 0;
  }
pod2usage(2) if ( not defined $output_prefix and not defined $neat_bam and not defined $merged_bam and not defined $ref_genome );

die ("\nOutput prefix not defined\n") unless( $output_prefix );
die ("\nMust define the reference\n") unless ( $ref_genome);
die ("\nReference $ref_genome not found\n") unless ( -e $ref_genome );
die ("\nMust define a neat normal BAM\n") unless ( $neat_bam );
die ("\nFile $neat_bam not found\n") unless ( -e $neat_bam );
die ("\nIndex for $neat_bam not found\n") unless ( -e "$neat_bam".".bai" );
die ("\nMust define a tumour duplex BAM\n") unless ( $merged_bam );
die ("\nFile $merged_bam not found\n") unless ( -e $merged_bam );
die ("\nIndex for $merged_bam not found\n") unless ( -e "$merged_bam".".bai" );
die ("\nRscript not found in path\n") unless ( which 'Rscript' );
die ("\nefficiency_nanoseq.R must be in path\n") unless ( which 'efficiency_nanoseq.R' );

my $threads = 1;
$threads = $opts{'t'} if ( defined $opts{'t'} );
my $rb_output     = "$output_prefix.RBs";
my $main_output   = "$output_prefix.tsv";
# Get the first contig/chr:
my $region        = `samtools view -H $neat_bam | grep "^\@SQ" | head -1 | cut -f2`;
$region           =~ s/SN://;
$region           = chomp($region);

##########################################################################################
# Calculating number of reads and duplicate rates
my($num_unique_reads,$num_sequenced_reads,$dup_rate);

print STDOUT "Calculating number of reads in $neat_bam...\n";
my ($stdout, $stderr, $exit) = capture {
    system("samtools view -@ $threads -c $neat_bam");
};
die "Error calling samtools view -@ $threads -c $neat_bam, $stderr\n" if ( $exit != 0 );
chomp( $stdout);
$num_unique_reads = $stdout;

print STDOUT "  Num unique reads=$num_unique_reads\n";

print STDOUT "Calculating number of reads in $merged_bam...\n";
($stdout, $stderr, $exit) = capture {
    system("samtools view -@ $threads -c $merged_bam");
};
die "Error calling samtools view -@ $threads -c $merged_bam, $stderr\n" if ( $exit != 0 );
chomp( $stdout);
$num_sequenced_reads = $stdout;

print STDOUT "  Num sequenced reads=$num_sequenced_reads\n";

$dup_rate = ($num_sequenced_reads-$num_unique_reads)/$num_sequenced_reads;
print STDOUT "  Duplicate rate=$dup_rate\n";

##########################################################################################
# Get read bundle comformations
my %rbs;
my $bam = $merged_bam;
# Get first reads in reverse:
print STDOUT "RB comformation: 1st reads in reverse...\n";
open(IN, "samtools view -@ $threads -f 82 $bam $region |") || die "Error launching samtools view -@ $threads -f 82 $bam $region\n"; 
while(<IN>) {
	chomp;
	my @tmp = (split(/\t/,$_));
	my $rb;
	foreach my $t ( @tmp ) {
		$rb = $t if($t =~ /^RB/);
	}
	$rbs{$rb}->{"f2r1"}++;
}
close(IN) or die ("error when calling samtools: $?, $!\n");

print STDOUT "RB comformation: 2nd reads in reverse...\n";
open(IN, "samtools -@ $threads view -f 146 $bam $region |") || die "samtools view -@ $threads -f 146 $bam $region\n"; 
while(<IN>) {
	chomp;
	my @tmp = (split(/\t/,$_));
        my $rb;
        foreach my $t ( @tmp ) {
                $rb = $t if($t =~ /^RB/);
        }
	$rbs{$rb}->{"f1r2"}++;
}
close(IN) or die ("error when calling samtools: $?, $!\n");

open(OUT,">$rb_output") || die "Error openning output $rb_output\n";
foreach my $rb ( keys %rbs ) {
	$rbs{$rb}->{"f1r2"} = 0 if(!exists($rbs{$rb}->{"f1r2"}));
	$rbs{$rb}->{"f2r1"} = 0 if(!exists($rbs{$rb}->{"f2r1"}));
	print OUT "$rb\t",$rbs{$rb}->{"f1r2"},"\t",$rbs{$rb}->{"f2r1"},"\n";
}
close(OUT);

##########################################################################################
# Call R to get the two values that inform on strand misses:
my($reads_per_rb,$f_eff,$zib_eff,$ok_rbs,$total_rbs,$gc_both,$gc_single,$total_reads);
print STDOUT "Running: efficiency_nanoseq.R $rb_output $ref_genome\n";
open(IN, "efficiency_nanoseq.R $rb_output $ref_genome |") || die "Error running efficiency_nanoseq.R $rb_output $ref_genome\n";
while(<IN>) {
	print STDOUT "  Routput: ",$_;
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
close(IN)  or die ("error when calling samtools: $?, $!\n");
# Check it run properly:
open(OUT, ">$main_output") || die "Error writing to $main_output\n";
print OUT "# Whole-genome metrics:\n";
print OUT "NUM_UNIQUE_READS\t$num_unique_reads\n";
print OUT "NUM_SEQUENCED_READS\t$num_sequenced_reads\n";
print OUT "DUPLICATE_RATE\t$dup_rate\n";

print OUT "# RB metrics are reported for chr/contig $region only:\n";
if(!defined($zib_eff) || $zib_eff eq "") {
	print STDOUT "ERROR: R script didn't worke correctly\n";
} 
my $bases_sequenced = $total_reads * 150;
my $bases_ok_rbs    = $ok_rbs * (300-50); # assuming mates don't overlap. Removing 50 bps for varios trimmings (rough estimate)
print OUT "TOTAL_RBS\t$total_rbs\n";
print OUT "TOTAL_READS_IN_RBS\t$total_reads\n";
print OUT "OK_RBS(2+2)\t$ok_rbs\n";
print OUT "READS_PER_RB\t$reads_per_rb\n";
print OUT "F-EFF\t$f_eff\n";
#print OUT "ZIB-EFF\t$zib_eff\n";
print OUT "EFFICIENCY\t",$bases_ok_rbs / $bases_sequenced,"\n";
#print OUT "EFFICIENCY2\t",$bases_ok_rbs / ($num_sequenced_reads*150),"\n";
print OUT "GC_BOTH\t$gc_both\n";
print OUT "GC_SINGLE\t$gc_single\n";
close(OUT);

__END__
=head1 NAME

efficiency_nanoseq.pl

=head1 SYNOPSIS

efficiency_nanoseq.pl [-h -v] -normal BAM -tumour BAM -o prefix -r reference

=head1 OPTIONS

=over 8

=item B<-normal>

Normal BAM with index

=item B<-normal>

Tumour (Duplex) BAM with index

=item B<-out>

Output prefix

=item B<-ref>

Reference file

=item B<-threads>

Use samtools with these many threads

=back

=cut

