#ifndef __CREC__
#define  __CREC__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>


class CRec
{
 private:
  std::vector< std::pair <int,double> > recs;
 public:
  CRec();
  CRec(const char* fn);
  void Init(const char* fn);
  double GetPos(int pos) const;
};
#endif
