#!/usr/bin/perl -w

########## LICENCE ##########
# Copyright (c) 2020-2021 Genome Research Ltd
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
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;

my $VERSION="1.0.0";

my %opts;
my $MIN_SIZE_SUBFAM = 2; #Minimum number of family size 2 for 2+2, 3 for 3+3, etc
my $FILTER_5_PRIME  = 10; #Number of bases to be trimmed from 5'. Set as 0 if no filter is wanted. Example: 10 for the first 10 bases
my $FILTER_3_PRIME  = 135; #Number of bases to start trimming from. Set as 0 if no filter is wanted. For example: 140 for the 10 last bases of 150-bp reads
my $BULK_MIN_COV    = 16;

GetOptions('rb|reads-bundle=i'  => \$MIN_SIZE_SUBFAM,
           't3|trim3=i'  => \$FILTER_3_PRIME,
           't5|trim5=i'  => \$FILTER_5_PRIME, 
           'mc|min-coverage=i' => \$BULK_MIN_COV,
           'o|out=s' => \$opts{'o'},
           'h|help' => \$opts{'h'},
           'v|version' => \$opts{'v'},
) or pod2usage(2);
pod2usage(-verbose => 1) if(defined $opts{'h'});
if(defined $opts{'v'}) {
  print sprintf "VERSION: %s\n", $VERSION;
  exit 0;
  }
pod2usage(2) if( @ARGV != 1 );
die ("\nOutput file not defined\n") unless(defined $opts{'o'});
open(OUT,"| gzip >$opts{'o'}") or die ("\nCouldn't stream to output file $!\n");

my $FILE = $ARGV[0];
die ("\nInput file $FILE not found\n") unless ( -e $FILE );
open(IN, "zcat $FILE |") or die( "\nProblem with gunzip $FILE\n" );

my %complement = ( "C" => "G", "G" => "C", "A" => "T", "T" => "A", "x" => "x", "?" => "?", "-" => "-", "N" => "N" ); 

while(<IN>) {
	next if(/^#/); #in case it is receiving multiple csv files!
	chomp;

  my($chrom,$chromBeg,$chromEnd,$context,$commonSNP,$shearwater,$bulkASXS,$bulkNM,
       $bulkForwardA,$bulkForwardC,$bulkForwardG,$bulkForwardT,$bulkForwardIndel,
       $bulkReverseA,$bulkReverseC,$bulkReverseG,$bulkReverseT,$bulkReverseIndel,
       $dplxBreakpointBeg,$dplxBreakpointEnd,$dplxBarcode,$dplxOri,$dplxASXS,$dplxCLIP,$dplxNM,
       $dplxForwardA,$dplxForwardC,$dplxForwardG,$dplxForwardT,$dplxForwardIndel,
       $dplxReverseA,$dplxReverseC,$dplxReverseG,$dplxReverseT,$dplxReverseIndel,
       $dplxCQForwardA,$dplxCQForwardC,$dplxCQForwardG,$dplxCQForwardT,
       $dplxCQReverseA,$dplxCQReverseC,$dplxCQReverseG,$dplxCQReverseT,
       $bulkProperPair,$dplxProperPair) =split(/\t/,$_);
	   
	my $dplxForwardTotal = $dplxForwardA+$dplxForwardC+$dplxForwardG+$dplxForwardT+$dplxForwardIndel;
	my $dplxReverseTotal = $dplxReverseA+$dplxReverseC+$dplxReverseG+$dplxReverseT+$dplxReverseIndel;
	
	my $bulkForwardTotal = $bulkForwardA+$bulkForwardC+$bulkForwardG+$bulkForwardT+$bulkForwardIndel;
	my $bulkReverseTotal = $bulkReverseA+$bulkReverseC+$bulkReverseG+$bulkReverseT+$bulkReverseIndel;

	# Translate variable names (for compatibility with previous version):
  my $rbase  = (split(//,$context))[1];
  my $site   = $chromEnd;
  my $coord1 = $dplxBreakpointBeg;
  my $coord2 = $dplxBreakpointEnd;
  my $A      = $dplxForwardA;
  my $a      = $dplxReverseA;
  my $C      = $dplxForwardC;
  my $c      = $dplxReverseC;
  my $G      = $dplxForwardG;
  my $g      = $dplxReverseG;
  my $T      = $dplxForwardT;
  my $t      = $dplxReverseT;
	#my $r1     = $A+$C+$G+$T; # $fwdDP
	#my $r2     = $a+$c+$g+$t; # $revDP;
	my $r1     = $dplxForwardTotal;
	my $r2     = $dplxReverseTotal;
	my $dp     = ($r1+$r2)."[$r1,$r2]";
	my $chr    = $chrom;
	my $position_zerobased = $chromBeg;
	my $position_onebased  = $chromEnd;
	my $beg_breakpoint_coordinate = $dplxBreakpointBeg;
	my $end_breakpoint_coordinate = $dplxBreakpointEnd;

	my $bulk_fwd = $bulkForwardTotal;
	my $bulk_rev = $bulkReverseTotal;
	next if($bulk_fwd + $bulk_rev < $BULK_MIN_COV); # Bulk minimum coverage
	next if($dplxCLIP>0);
	next if($dplxNM > 20); # made very liberal to allow long indels. Check the impact!
	next if($dplxASXS < 50);

	if($r1 >= $MIN_SIZE_SUBFAM && $r2 >= $MIN_SIZE_SUBFAM) {
		my $bulktotal = $bulkForwardTotal+$bulkReverseTotal;
		if($dplxForwardIndel+$dplxReverseIndel < 0.9*($dplxForwardTotal+$dplxReverseTotal)) {
			next; # only interested in indels
		}
		my $qpos;
		my $orientation_type;
		my $qposF = $chromBeg-$dplxBreakpointBeg + 1;
		my $qposR = $dplxBreakpointEnd-$chromEnd;
		if($qposR < $qposF) {
			$qpos = $qposR;
			$orientation_type = "Rev";
		} else {
			$qpos = $qposF;
			$orientation_type = "Fwd";
		}
		if($FILTER_5_PRIME > 0) {
			if($qpos < $FILTER_5_PRIME) {
				next;
			}
		}
		if($FILTER_3_PRIME > 0) {
			if($qpos > $FILTER_3_PRIME) {
				next;
			}
		}
		my $dp = $r1+$r2;
		my $site_tags = "";
		$site_tags .= "$chr:$coord1-$coord2:$dplxBarcode;DP=$dp;QPOS=$qpos;";
		
		# Annotate the sequence context:
		my $signature;
		my $signature_trinuc;
		$signature_trinuc = $context    ;
		$signature_trinuc =~ s/\./$rbase/;
		$signature_trinuc .= ">indel";
		if($rbase =~ /[AGag]/) {
			$signature_trinuc = &reverse_signature($signature_trinuc);
		}
		$site_tags .= "$signature_trinuc;SW=$shearwater;cSNP=$commonSNP";

		# If seen in the bulk, flag it:
		if($bulkForwardIndel+$bulkReverseIndel > 1) {
			$site_tags .= ";BULK_SEEN($dplxForwardIndel+$dplxReverseIndel/$bulktotal)";
		} else {
      $site_tags .= ";";
    }
		print OUT "$chr\t",$site-1,"\t$site\t$site_tags\n"; #BED Format including site tags!
	}
}
OUT->flush;
close(OUT);

sub reverse_signature {
	my $signature = uc($_[0]);
	my @tmp = split(//,$signature);
	$signature = $complement{$tmp[2]}.$complement{$tmp[1]}.$complement{$tmp[0]}.">"."indel";
	return $signature;
}


__END__

=head1 NAME

indelCaller_step1.pl - Filter sites of interest to look for indels.

=head1 SYNOPSIS

indelCaller_step1.pl  [options] -o out.bed.gz  input.bed.gz

    -out               -o   Output file

  Optional parameters:
    -reads-bundle      -rb  Minimum reads in a bundle. (2)
    -trim3             -t3  Excess bases above this value are trimmed from 3' (135)
    -trim5             -t5  bases to trim from 5' reads (10)
    -min-coverage      -mc  minimum bulk coverage (16)
    -help              -h
    -version           -v

=head1 OPTIONS

=over 8

=item B<-out>

Output file. 
Final ouptput is a fileterd bed file.

=item B<-read-bundles>

Minimum number of reads per strand to consider a read bundle for indel processing.

=item B<-trim5>

Number of bases to trim from the 5' end of reads. Set to 0 for no filter

=item B<-trim3>

Total lenght of the read to consider, bases beyond this value are trimmed from the 3' end. Eg. trim 10 last bases if set to 150. Set to 0 for no filter

=item B<-min-coverage>

Minimum bulk coverage

=back

=cut
