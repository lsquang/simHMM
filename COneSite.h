#ifndef __CONESITE__
#define __CONESITE__

#include <string>
#include <vector>
#include "const.h"


struct onesite
{
  std::string al1, al2, hap, rsid;
  std::string inputhap;
  vecdouble lk, pplk;  
  vecdouble prob0; //probability of 0
  //lk: n*3 probability of data given genotype 00 01 and 11
  //posterior probability of haplotype 00, 10, 10, 11
  int pos;
  double pos_cm;
};

bool operator == (const onesite &x, const onesite &y);
bool operator < (const onesite &x, const onesite &y);

#endif

