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

use strict;
use List::Util qw(min);
$|=1;

# fa8, Nov 2023
# This script reads the discardedvariants.csv file and creates a summary VCF
# Output is written to standard output
# E.g. perl post_process_discarded_variants.pl tmpNanoSeq/post/discardedvariants.csv reference_genome.fa > tmpNanoSeq/post/discardedvariants.vcf



my $input_file = $ARGV[0];
my $ref_genome = $ARGV[1];
if(!defined($ref_genome) || $ref_genome eq "") {
	die "Reference genome fasta file needs to be provided. Exiting...\n".
	"Usage:  perl post_process_discarded_variants.pl discardedvariants.csv reference_genome.fa > discardedvariants.vcf\n";
}

# Prepare header:
my $header = 
"##fileformat=VCFv4.2
##source=NanoSeq pipeline (discarded variants from post_process_discarded_variants.pl)
##reference=file://$ref_genome 
##INFO=<ID=PYR_SUB,Number=1,Type=String,Description=\"Pyrimidine-based trinucleotide substitution\">
##INFO=<ID=N,Number=1,Type=Integer,Description=\"Number of times this variant has been discarded\">
##INFO=<ID=dplx_clip_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the dplx_clip_filter\">
##INFO=<ID=alignment_score_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the alignment_score_filter (AS-XS)\">
##INFO=<ID=mismatch_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the mismatch_filter\">
##INFO=<ID=matched_normal_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the matched_normal_filter\">
##INFO=<ID=duplex_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the duplex_filter\">
##INFO=<ID=consensus_base_quality_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the consensus_base_quality_filter\">
##INFO=<ID=indel_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the indel_filter\">
##INFO=<ID=five_prime_trim_filter,Number=1,Type=Integer,Description=\"Number of times this variant has failed the five_prime_trim_filter\">
##INFO=<ID=three_prime_trim_filter,Number=1,Type=String,Description=\"Number of times this variant has failed the three_prime_trim_filter\">
##INFO=<ID=proper_pair_filter,Number=1,Type=String,Description=\"Number of times this variant has failed the proper_pair_filter\">
##INFO=<ID=vaf_filter,Description=\"VAF in matched normal higher than threshold\">
##INFO=<ID=QPOS,Number=.,Type=Integer,Description=\"Read position closest to 5-prime end. Up to 10 QPOS are reported\">
##INFO=<ID=NORM_VAF,Number=1,Type=Float,Description=\"VAF in matched normal\">
##INFO=<ID=MEAN_DX_ASXS,Number=1,Type=Float,Description=\"mean AS-XS for duplex\">
##INFO=<ID=MEAN_NORM_ASXS,Number=1,Type=Float,Description=\"mean AS-XS for normal\">
##INFO=<ID=MEAN_DX_NM,Number=1,Type=Float,Description=\"mean NM for duplex\">
##INFO=<ID=MEAN_NORM_NM,Number=1,Type=Float,Description=\"mean NM for normal\">
##INFO=<ID=NORM_COV,Number=1,Type=Integer,Description=\"Coverage in the matched normal\">
##FILTER=<ID=commonSNP,Description=\"Common SNP site\">
##FILTER=<ID=shearwater,Description=\"Noisy site\">
##FILTER=<ID=not_in_masks,Description=\"Not in the commonSNP and noise masks\">\n";
open(I, "$ref_genome.fai") || die "Error: cannot find file $ref_genome.fai\n";
while(<I>) {
	chomp;
	my($contig_name,$length) = (split)[0,1];
	$header .= "##contig=<ID=$contig_name,length=$length>\n";
}
close(I);
$header .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
print $header;

my %ds;
my %complement;
$complement{"A"} = "T";
$complement{"C"} = "G";
$complement{"G"} = "C";
$complement{"T"} = "A";
	print STDERR "Reading $input_file...\n";
	open(I, $input_file) || die "Error: cannot open $input_file for reading\n";
	while(<I>) {
		next if(/^#/);     # ignore header
		next if(/^chrom/); # ignore header
		chomp;
		my($chrom,$chromStart,$context,$commonSNP,$shearwater,$bulkASXS,$bulkNM,$bulkForwardA,
		   $bulkForwardC,$bulkForwardG,$bulkForwardT,$bulkForwardIndel,$bulkReverseA,$bulkReverseC,
		   $bulkReverseG,$bulkReverseT,$bulkReverseIndel,$dplxBreakpointBeg,$dplxBreakpointEnd,
		   $bundleType,$dplxASXS,$dplxCLIP,$dplxNM,
		   $dplxfwdA,$dplxfwdC,$dplxfwdG,$dplxfwdT,$dplxfwdIndel,
		   $dplxrevA,$dplxrevC,$dplxrevG,$dplxrevT,$dplxrevIndel,
		   $dplxCQfwdA,$dplxCQfwdC,$dplxCQfwdG,$dplxCQfwdT,$dplxCQrevA,$dplxCQrevC,$dplxCQrevG,$dplxCQrevT,
		   $bulkForwardTotal,$bulkReverseTotal,$dplxfwdTotal,$dplxrevTotal,$left,$right,
		   $qpos,$call,$isvariant,$pyrcontext,$stdcontext,$pyrsub,$stdsub,$ismasked,$dplxBarcode,
		   $dplx_clip_filter,$alignment_score_filter,$mismatch_filter,$matched_normal_filter,
		   $duplex_filter,$consensus_base_quality_filter,$indel_filter,$five_prime_trim_filter,
		   $three_prime_trim_filter,$proper_pair_filter,$vaf_filter) = split(/,/,$_);
		# Find which base has highest support (mut base):
		# my $mut = $call;
		# Debugging, debugging...
		# print STDERR "(qpos=$qpos,$call,$isvariant,std=$stdcontext,pyr=$pyrcontext,$pyrsub,$stdsub) jajaja\n";
		# print STDERR "chrom=$chrom,$chromStart,$context,$commonSNP,$shearwater,$bulkASXS,$bulkNM,$bulkForwardA,
		#    bulkForwardC=$bulkForwardC,$bulkForwardG,$bulkForwardT,$bulkForwardIndel,$bulkReverseA,$bulkReverseC,
		#    $bulkReverseG,$bulkReverseT,$bulkReverseIndel,$dplxBreakpointBeg,$dplxBreakpointEnd,
		#    $bundleType,$dplxASXS,$dplxCLIP,$dplxNM,
		#    $dplxfwdA,$dplxfwdC,$dplxfwdG,$dplxfwdT,$dplxfwdIndel,
		#    dplxrevA=$dplxrevA,$dplxrevC,$dplxrevG,$dplxrevT,$dplxrevIndel,
		#    $dplxCQfwdA,$dplxCQfwdC,$dplxCQfwdG,$dplxCQfwdT,$dplxCQrevA,$dplxCQrevC,$dplxCQrevG,$dplxCQrevT,
		#    $bulkForwardTotal,$bulkReverseTotal,$dplxfwdTotal,$dplxrevTotal,$left,$right,
		#    qpos=$qpos,$call,$isvariant,$pyrcontext,$stdcontext,$pyrsub,$stdsub,$ismasked,$dplxBarcode,
		#    $dplx_clip_filter,$alignment_score_filter,$mismatch_filter,$matched_normal_filter,
		#    $duplex_filter,$consensus_base_quality_filter,$indel_filter,$five_prime_trim_filter,
		#    etc=$three_prime_trim_filter,$proper_pair_filter\n";
		# print STDERR "commonSNP=$commonSNP,shearwater=$shearwater,bulkASXS=$bulkASXS\n";
		# print STDERR "$bulkForwardA,$bulkForwardC,$bulkForwardG,$bulkForwardT,$bulkForwardIndel,$bulkReverseA,$bulkReverseC,$bulkReverseG,$bulkReverseT,$bulkReverseIndel\n";
		# Create id:
		my $ref = (split(//,$context))[1];
		my($ref_prelim,$mut) = (split(//,$pyrsub))[1,4];
		if($ref_prelim ne $ref) {
			$mut = $complement{$mut};
		}
		my $pos = $chromStart + 1;
		my $id = $chrom.":".$pos.":".$ref.">".$mut;  
		# chromStart is 0-based
		if(!exists($ds{"$chrom:$pos:$ref:$mut"})) {
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_ASXS"                 } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_ASXS"               } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"counter"                      } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ"                  } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_NM"                   } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_NM"                 } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"dplx_clip_filter"             } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"alignment_score_filter"       } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"mismatch_filter"              } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"matched_normal_filter"        } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"duplex_filter"                } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"consensus_base_quality_filter"} = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"indel_filter"                 } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"five_prime_trim_filter"       } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"three_prime_trim_filter"      } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"proper_pair_filter"           } = 0;
			$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_filter"                   } = 0;
		}
		$ds{"$chrom:$pos:$ref:$mut"}->{"counter"}++;
		$ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_ASXS"   } = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_ASXS"  } + $dplxASXS;
		# print STDERR $ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_ASXS"}, "+ $bulkASXS\n";
		$ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_ASXS" } = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_ASXS"} + $bulkASXS;
		$ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_NM"     } = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_dx_NM"    } + $dplxNM;
		$ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_NM"   } = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_norm_NM"  } + $bulkNM;
		my $normal_coverage = $bulkForwardA+$bulkForwardC+$bulkForwardG+$bulkForwardT+$bulkForwardIndel+
		                      $bulkReverseA+$bulkReverseC+$bulkReverseG+$bulkReverseT+$bulkReverseIndel;
		# print STDERR "mut = $mut; normal_coverage = $normal_coverage\n";
		# print STDERR "$bulkForwardA+$bulkForwardC+$bulkForwardG+$bulkForwardT+$bulkForwardIndel+$bulkReverseA+$bulkReverseC+$bulkReverseG+$bulkReverseT+$bulkReverseIndel\n";
		$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_normal" } = "NA";
		if($call eq "A") {
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ"} = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ" } + min($dplxCQfwdA,$dplxCQrevA);
			$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_normal" } = ($bulkForwardA + $bulkReverseA)/$normal_coverage if($normal_coverage > 0);
		} elsif($call eq "C") {  
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ"} = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ" } + min($dplxCQfwdC,$dplxCQrevC);
			$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_normal" } = ($bulkForwardC + $bulkReverseC)/$normal_coverage if($normal_coverage > 0);
		} elsif($call eq "G") {  
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ"} = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ" } + min($dplxCQfwdG,$dplxCQrevG);
			$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_normal" } = ($bulkForwardG + $bulkReverseG)/$normal_coverage if($normal_coverage > 0);
		} elsif($call eq "T") {  
			$ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ"} = $ds{"$chrom:$pos:$ref:$mut"}->{"mean_min_BQ" } + min($dplxCQfwdT,$dplxCQrevT);
			$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_normal" } = ($bulkForwardT + $bulkReverseT)/$normal_coverage if($normal_coverage > 0);
		} else {
			die "Error: $call doesn't match A, C, G, T. Exiting...\n";
		
		}
		push(@{$ds{"$chrom:$pos:$ref:$mut"}->{"QPOS"}},$qpos);
		$ds{"$chrom:$pos:$ref:$mut"}->{"dplx_clip_filter"             } = $ds{"$chrom:$pos:$ref:$mut"}->{"dplx_clip_filter"             } + $dplx_clip_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"alignment_score_filter"       } = $ds{"$chrom:$pos:$ref:$mut"}->{"alignment_score_filter"       } + $alignment_score_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"mismatch_filter"              } = $ds{"$chrom:$pos:$ref:$mut"}->{"mismatch_filter"              } + $mismatch_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"matched_normal_filter"        } = $ds{"$chrom:$pos:$ref:$mut"}->{"matched_normal_filter"        } + $matched_normal_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"duplex_filter"                } = $ds{"$chrom:$pos:$ref:$mut"}->{"duplex_filter"                } + $duplex_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"consensus_base_quality_filter"} = $ds{"$chrom:$pos:$ref:$mut"}->{"consensus_base_quality_filter"} + $consensus_base_quality_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"indel_filter"                 } = $ds{"$chrom:$pos:$ref:$mut"}->{"indel_filter"                 } + $indel_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"five_prime_trim_filter"       } = $ds{"$chrom:$pos:$ref:$mut"}->{"five_prime_trim_filter"       } + $five_prime_trim_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"three_prime_trim_filter"      } = $ds{"$chrom:$pos:$ref:$mut"}->{"three_prime_trim_filter"      } + $three_prime_trim_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"proper_pair_filter"           } = $ds{"$chrom:$pos:$ref:$mut"}->{"proper_pair_filter"           } + $proper_pair_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"vaf_filter"                   } = $ds{"$chrom:$pos:$ref:$mut"}->{"vaf_filter"                   } + $vaf_filter;
		$ds{"$chrom:$pos:$ref:$mut"}->{"cov_normal"                   } = $normal_coverage;
		$ds{"$chrom:$pos:$ref:$mut"}->{"commonSNP"                    } = $commonSNP;
		$ds{"$chrom:$pos:$ref:$mut"}->{"shearwater"                   } = $shearwater;
		$ds{"$chrom:$pos:$ref:$mut"}->{"pyrsub"                       } = $pyrsub;
	}
	close(I);


print STDERR scalar(keys(%ds)), " different variants\n";
foreach my $id ( keys %ds ) {
	my($chr,$pos,$ref,$mut) = (split(/:/,$id))[0,1,2,3];
	print "$chr\t$pos\t.\t$ref\t$mut\t.";
	my $filter = "";
	if($ds{$id}->{"commonSNP" } == 1) {
		$filter = "commonSNP;";
	}
	if($ds{$id}->{"shearwater" } == 1) {
		$filter .= "shearwater;";
	}
	if($filter eq "") {
		$filter = "not_in_masks;";
	}
	chop($filter);
	print "\t$filter";
	print "\tPYR_SUB=",$ds{$id}->{"pyrsub"},";";
	print "N=",$ds{$id}->{"counter"},";";
	print "dplx_clip_filter="              ,$ds{$id}->{"counter"}-$ds{$id}->{"dplx_clip_filter"},";"              if($ds{$id}->{"dplx_clip_filter"}              != $ds{$id}->{"counter"});
	print "alignment_score_filter="        ,$ds{$id}->{"counter"}-$ds{$id}->{"alignment_score_filter"},";"        if($ds{$id}->{"alignment_score_filter"}        != $ds{$id}->{"counter"});
	print "mismatch_filter="               ,$ds{$id}->{"counter"}-$ds{$id}->{"mismatch_filter"},";"               if($ds{$id}->{"mismatch_filter"}               != $ds{$id}->{"counter"});
	print "matched_normal_filter="         ,$ds{$id}->{"counter"}-$ds{$id}->{"matched_normal_filter"},";"         if($ds{$id}->{"matched_normal_filter"}         != $ds{$id}->{"counter"});
	print "duplex_filter="                 ,$ds{$id}->{"counter"}-$ds{$id}->{"duplex_filter"},";"                 if($ds{$id}->{"duplex_filter"}                 != $ds{$id}->{"counter"});
	print "consensus_base_quality_filter=" ,$ds{$id}->{"counter"}-$ds{$id}->{"consensus_base_quality_filter"},";" if($ds{$id}->{"consensus_base_quality_filter"} != $ds{$id}->{"counter"});
	print "indel_filter="                  ,$ds{$id}->{"counter"}-$ds{$id}->{"indel_filter"},";"                  if($ds{$id}->{"indel_filter"}                  != $ds{$id}->{"counter"});
	print "five_prime_trim_filter="        ,$ds{$id}->{"counter"}-$ds{$id}->{"five_prime_trim_filter"},";"        if($ds{$id}->{"five_prime_trim_filter"}        != $ds{$id}->{"counter"});
	print "three_prime_trim_filter="       ,$ds{$id}->{"counter"}-$ds{$id}->{"three_prime_trim_filter"},";"       if($ds{$id}->{"three_prime_trim_filter"}       != $ds{$id}->{"counter"});
	print "proper_pair_filter="            ,$ds{$id}->{"counter"}-$ds{$id}->{"proper_pair_filter"},";"            if($ds{$id}->{"proper_pair_filter"}            != $ds{$id}->{"counter"});
	print "vaf_filter="                    ,$ds{$id}->{"counter"}-$ds{$id}->{"vaf_filter"},";"                    if($ds{$id}->{"vaf_filter"}                    != $ds{$id}->{"counter"});
	my @qpos = @{$ds{$id}->{"QPOS"}};
	my $max_to_print = 10;
	my $qpos_str = "";
	for(my $i=0; $i<=$#qpos; $i++) {
		$qpos_str .= $qpos[$i].",";
		last if($i >= $max_to_print);
	}
	chop($qpos_str);
	print "QPOS=$qpos_str;";
	$ds{$id}->{"mean_dx_ASXS"  } = $ds{$id}->{"mean_dx_ASXS"  } / $ds{$id}->{"counter"};
	$ds{$id}->{"mean_norm_ASXS"} = $ds{$id}->{"mean_norm_ASXS"} / $ds{$id}->{"counter"};
	$ds{$id}->{"mean_dx_NM"    } = $ds{$id}->{"mean_dx_NM"    } / $ds{$id}->{"counter"};
	$ds{$id}->{"mean_norm_NM"  } = $ds{$id}->{"mean_norm_NM"  } / $ds{$id}->{"counter"};
	if($ds{$id}->{"vaf_normal"} eq "NA") {
		print "NORM_VAF=",   $ds{$id}->{"vaf_normal"},";";
	} else {
		print "NORM_VAF=",   sprintf("%.2f", $ds{$id}->{"vaf_normal"    }),";";
	}
	print "NORM_COV=",       sprintf("%.2f", $ds{$id}->{"cov_normal"    }),";";
	print "MEAN_DX_ASXS=",   sprintf("%.2f", $ds{$id}->{"mean_dx_ASXS"  }),";";
	print "MEAN_NORM_ASXS=", sprintf("%.2f", $ds{$id}->{"mean_norm_ASXS"}),";";
	print "MEAN_DX_NM=",     sprintf("%.2f", $ds{$id}->{"mean_dx_NM"    }),";";
	print "MEAN_NORM_NM=",   sprintf("%.2f", $ds{$id}->{"mean_norm_NM"  }),";";
	print "\n";
}




__END__

Header of discarded_var files:
 1 chrom
 2 chromStart
 3 context
 4 commonSNP
 5 shearwater
 6 bulkASXS
 7 bulkNM
 8 bulkForwardA
 9 bulkForwardC
10 bulkForwardG
11 bulkForwardT
12 bulkForwardIndel
13 bulkReverseA
14 bulkReverseC
15 bulkReverseG
16 bulkReverseT
17 bulkReverseIndel
18 dplxBreakpointBeg
19 dplxBreakpointEnd
20 bundleType
21 dplxASXS
22 dplxCLIP
23 dplxNM
24 dplxfwdA
25 dplxfwdC
26 dplxfwdG
27 dplxfwdT
28 dplxfwdIndel
29 dplxrevA
30 dplxrevC
31 dplxrevG
32 dplxrevT
33 dplxrevIndel
34 dplxCQfwdA
35 dplxCQfwdC
36 dplxCQfwdG
37 dplxCQfwdT
38 dplxCQrevA
39 dplxCQrevC
40 dplxCQrevG
41 dplxCQrevT
42 bulkForwardTotal
43 bulkReverseTotal
44 dplxfwdTotal
45 dplxrevTotal
46 left
47 right
48 qpos
49 call
50 isvariant
51 pyrcontext
52 stdcontext
53 pyrsub
54 stdsub
55 ismasked
56 dplxBarcode
57 dplx_clip_filter
58 alignment_score_filter
59 mismatch_filter
60 matched_normal_filter
61 duplex_filter
62 consensus_base_quality_filter
63 indel_filter
64 five_prime_trim_filter
65 three_prime_trim_filter
66 proper_pair_filter
