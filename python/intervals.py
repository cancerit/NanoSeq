#!/usr/bin/python

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
###########################

"""Build a BED file of intervals
"""

rnames = { 1  : 249250621,
           2  : 243199373,
           3  : 198022430,
           4  : 191154276,
           5  : 180915260,
           6  : 171115067,
           7  : 159138663,
           8  : 146364022,
           9  : 141213431,
           10 : 135534747,
           11 : 135006516,
           12 : 133851895,
           13 : 115169878,
           14 : 107349540,
           15 : 102531392,
           16 : 90354753,
           17 : 81195210,
           18 : 78077248,
           19 : 59128983,
           20 : 63025520,
           21 : 48129895,
           22 : 51304566,
           "X": 155270560,
           "Y": 59373566 }


# interval size
chunk = 1000000

# write out new intervals in bed format
for rname in sorted(rnames.keys()):
  prev = 0
  length = rnames[rname]
  for curr in xrange(chunk, length, chunk):
    print '%s\t%s\t%s' % (rname, prev, curr - 1)
    prev = curr
  if prev < length:
    curr = prev + chunk
    print '%s\t%s\t%s' % (rname, prev, curr - 1)

