/**   LICENCE
* Copyright (c) 2020 Genome Research Ltd.
* 
* Author: Cancer Genome Project <cgphelp@sanger.ac.uk>
* 
* This file is part of NanoSeq.
* 
* NanoSeq is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
* 
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
* 
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "writeout.h"


int WriteOut::Round(double d) {
  int n = std::round(d);
  return (n == -0)? 0 : n;
}


std::string WriteOut::Counts(bundle *bin) {
  int alph_len = 5;
  char alph[5] = {'A', 'C', 'G', 'T', 'i'};
  std::stringstream ss;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < alph_len; j++) {
      ss << bin->counts[i][alph[j]];
      ss << "\t";
    }
  }
  return ss.str();
}


std::string WriteOut::BaseQuals(bundle *bin) {
  int alph_len = 4;
  std::stringstream ss;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < alph_len; j++) {
      ss << WriteOut::Round(bin->consensus[i][j]);
      ss << "\t";
    }
  }
  return ss.str();
}


float WriteOut::MeanOfOneVector(std::vector<int> v1) {
  float t = std::accumulate(v1.begin(), v1.end(), 0);
  float n = static_cast<float>(v1.size());
  if (n == 0) {
    return 0;
  } else {
    return t/n;
  }
}


float WriteOut::MeanOfTwoVectors(std::vector<int> v1, std::vector<int> v2) {
  for (int i = 0; i < v1.size(); i ++) {
    v2.push_back(v1[i]);
  }
  return WriteOut::MeanOfOneVector(v2);
}


std::string WriteOut::DplxCovariates(bundle *bin) {
  float asxs0 = WriteOut::MeanOfOneVector(bin->asxs[0]);
  float asxs1 = WriteOut::MeanOfOneVector(bin->asxs[1]);
  int asxs    = WriteOut::Round(std::min(asxs0, asxs1));
  float clip1 = WriteOut::MeanOfTwoVectors(bin->clip[0], bin->clip[1]);
  int clip    = WriteOut::Round(clip1);
  float nmms0 = WriteOut::MeanOfOneVector(bin->nmms[0]);
  float nmms1 = WriteOut::MeanOfOneVector(bin->nmms[1]);
  float nmms  = 0.1 * WriteOut::Round(std::max(nmms0, nmms1) * 10.0);
  std::stringstream ss;
  ss << asxs;
  ss << "\t";
  ss << clip;
  ss << "\t";
  ss << nmms;
  ss << "\t";
  return ss.str();
}


std::string WriteOut::BulkCovariates(bundle *bin) {
  float asxs0 = WriteOut::MeanOfOneVector(bin->asxs[0]);
  float asxs1 = WriteOut::MeanOfOneVector(bin->asxs[1]);
  // fa8 modified:
  // int asxs    = WriteOut::Round(std::min(asxs0, asxs1));
  int asxs;
  if(bin->asxs[0].size() == 0 && bin->asxs[1].size() != 0) {
  	asxs    = WriteOut::Round(asxs1);
  } else if(bin->asxs[0].size() != 0 && bin->asxs[1].size() == 0) {
  	asxs    = WriteOut::Round(asxs0);
  } else {
    asxs    = WriteOut::Round(std::min(asxs0, asxs1));
  }
  //end
  float nmms0 = WriteOut::MeanOfOneVector(bin->nmms[0]);
  float nmms1 = WriteOut::MeanOfOneVector(bin->nmms[1]);
  // fa8 modified:
  // float nmms  = 0.1 * WriteOut::Round(std::max(nmms0, nmms1) * 10.0);
  float nmms;
  if(bin->nmms[0].size() == 0 && bin->nmms[1].size() != 0) {
  	nmms    = WriteOut::Round(nmms1);
  } else if(bin->nmms[0].size() != 0 && bin->nmms[1].size() == 0) {
  	nmms    = WriteOut::Round(nmms0);
  } else {
    nmms    = WriteOut::Round(std::min(nmms0, nmms0));
  }  
  //end
  std::stringstream ss;
  ss << asxs;
  ss << "\t";
  ss << nmms;
  ss << "\t";
  return ss.str();
}


std::string WriteOut::ProperPairs(bundle *bin) {
  float ppair0 = MeanOfOneVector(bin->ppair[0]);
  float ppair1 = MeanOfOneVector(bin->ppair[1]);
  //fa8
  // float ppair  = 0.1 * WriteOut::Round(std::min(ppair0, ppair1) * 10.0);
  float ppair;
  if(bin->ppair[0].size() == 0 && bin->ppair[1].size() != 0) {
  	ppair    = WriteOut::Round(ppair1);
  } else if(bin->ppair[0].size() != 0 && bin->ppair[1].size() == 0) {
  	ppair    = WriteOut::Round(ppair0);
  } else {
    ppair    = WriteOut::Round(std::min(ppair0, ppair1));
  }
  // end
  std::stringstream ss;
  ss << ppair;
  ss << "\t";
  return ss.str();
}


std::string WriteOut::IdentifierString(bundle *bin) {
  std::stringstream ss;
  ss << bin->idf.beg;
  ss << "\t";
  ss << bin->idf.end;
  ss << "\t";
  ss << bin->idf.fwd_bc;
  ss << "|";
  ss << bin->idf.rev_bc;
  ss << "\t";
  ss << bin->bundle_type;
  ss << "\t";
  return ss.str();
}


void WriteOut::WriteRows(bundle bulk, bundles dplx, std::string posn, ogzstream & gzout ) {
  std::stringstream ss;
  std::string strb;
  ss << posn;
  ss << WriteOut::BulkCovariates(&bulk);
  ss << WriteOut::Counts(&bulk);
  std::string bulkrow = ss.str();
  for (auto bndls : dplx) {
    std::stringstream row;
    row << bulkrow;
    row << WriteOut::IdentifierString(&bndls.second);
    row << WriteOut::DplxCovariates(&bndls.second);
    row << WriteOut::Counts(&bndls.second);
    row << WriteOut::BaseQuals(&bndls.second);
    row << WriteOut::ProperPairs(&bulk);
    row << WriteOut::ProperPairs(&bndls.second);
    strb = row.str();
    strb.pop_back();//remove tab at the eol
    gzout << strb << std::endl ;
  }
}
