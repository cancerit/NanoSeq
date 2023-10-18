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

my %opts;
my $min_size_subfam = 2; #Minimum number of family size 2 for 2+2, 3 for 3+3, etc
my $filter_5_prime  = 10; #Number of bases to be trimmed from 5'. Set as 0 if no filter is wanted. Example: 10 for the first 10 bases
my $filter_3_prime  = 135; #Number of bases to start trimming from. Set as 0 if no filter is wanted. For example: 140 for the 10 last bases of 150-bp reads
my $bulk_min_cov    = 20;
my $min_asxs        = 50;
my $max_clip        = 0.02;
my $max_vaf         = 0.2;

GetOptions('rb|reads-bundle=i'  => \$min_size_subfam,
           't3|trim3=i'  => \$filter_3_prime,
           't5|trim5=i'  => \$filter_5_prime, 
           'mc|min-coverage=i' => \$bulk_min_cov,
           'vaf|max-vaf=f' => \$max_vaf,
           'a|min-as-xs=i' => \$min_asxs,
           'c|max-clip=f' => \$max_clip,
           'o|out=s' => \$opts{'o'},
           'h|help' => \$opts{'h'}
) or pod2usage(2);
pod2usage(-verbose => 1, -exitval => 0) if(defined $opts{'h'});
pod2usage(2) if( @ARGV != 1 );
die ("\nOutput file not defined\n") unless(defined $opts{'o'});
open(OUT,"| gzip >$opts{'o'}") or die ("\nCouldn't stream to output file $!\n");

my $FILE = $ARGV[0];
die ("\nInput file $FILE not found\n") unless ( -e $FILE );
open(IN, "zcat $FILE |") or die( "\nProblem with gunzip $FILE\n" );

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
  next if($bulk_fwd + $bulk_rev < $bulk_min_cov); # Bulk minimum coverage
  next if($dplxCLIP > $max_clip);
  next if($dplxNM > 20); # made very liberal to allow long indels. Check the impact! --> it seems is working, good
  next if($dplxASXS < $min_asxs || $bulkASXS < $min_asxs ); #fa8: fixed bug, we needed to check AS-XS for the bulk too
  															#     Fixed but: from bitwise OR to logical OR (ainsss)

  if($r1 >= $min_size_subfam && $r2 >= $min_size_subfam) {
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
    if($filter_5_prime > 0) {
      if($qpos < $filter_5_prime) {
        next;
      }
    }
    if($filter_3_prime > 0) {
      if($qpos > $filter_3_prime) {
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
    $site_tags .= ";BBEG=$dplxBreakpointBeg";
    $site_tags .= ";BEND=$dplxBreakpointEnd";
    $site_tags .= ";DEPTH_FWD=$r1";
    $site_tags .= ";DEPTH_REV=$r2";
    $site_tags .= ";DEPTH_NORM_FWD=$bulkForwardTotal";
    $site_tags .= ";DEPTH_NORM_REV=$bulkReverseTotal";
    $site_tags .= ";DPLX_ASXS=$dplxASXS";
    $site_tags .= ";DPLX_CLIP=$dplxCLIP";
    $site_tags .= ";DPLX_NM=$dplxNM";
    $site_tags .= ";BULK_ASXS=$bulkASXS";
    $site_tags .= ";BULK_NM=$bulkNM";
    # If seen in the bulk, flag it:
    if($bulkForwardIndel+$bulkReverseIndel > $max_vaf * $bulktotal) {
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
  $signature = reverse(substr($signature,0,3) =~ tr/ACGT/TGCA/r).q{>indel};
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
    -max-vaf           -vaf maximum VAF (0.2)
    -min-as-xs         -a   minimum AS-XS (50)
    -max-clip          -c   maximum clip fraction (0.02)
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

=item B<-max-vaf>

Tag sites with this VAF (or greater) in the bulk sample

=item B<-min-asxs>

Miniumum alignment score for the bulk reads


=item B<-max-clip>

Maximum clipping fraction for the bulk reads

=back

=cut
