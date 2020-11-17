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

import argparse
import os


"""Utility to check that batches have been successfully processed by LSF.
"""


def get_intervals():
  """Get expected intervals from chunks.bed.
  Returns:
    list of intervals
  """
  intervals = []
  for line in open('../bed/1.megabase.intervals.bed', 'rU'):
    rname, beg, end = line.strip().split('\t')
    intervals.append('%s.%s.%s' % (rname, beg, end))
  return intervals


def get_results(dirname):
  """Returns a dict of LSF results from one directory
  Args:
    dirname: name of directory containing LSF .out files
  Returns:
    dict with key = interval, value = LSF result
  """
  results = {}
  for fn in os.listdir(dirname):
    if fn.endswith('.out'):
      out = open('%s/%s' % (dirname, fn), 'rU')
      for line in out:
        if line.startswith('# LSBATCH: User input'):
          [out.next() for _ in xrange(3)]
          result = out.next().strip()
          results[fn.replace('.out','')] = result
    if fn.endswith('.log'):
      if os.path.getsize('%s/%s' % (dirname, fn)) != 0:
        print fn.replace('.log',''), 'Log file has size > 0 bytes.'
  return results


def check_results(intervals, results):
  """Check LSF results for all intervals.
  Args:
    intervals: expected intervals
    results: dict with key = interval, value = LSF result
  """
  for interval in intervals:
    result = results.get(interval, None)
    if result == None:
      result = 'Missing.'
    if result != "Successfully completed.":
      print interval, result
  print 'Complete'


def check_batch(dirname):
  """Check for errors on batch LSF run.
  """
  intervals = get_intervals()
  results   = get_results(dirname)
  check_results(intervals, results)


if __name__ == "__main__":
  msg = ('Check that batches have been successfully processed by LSF')
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument('-d', '--dirname', required=True)
  args = parser.parse_args()
  check_batch(args.dirname)
