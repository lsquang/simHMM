#ifndef __COPTION__
#define  __COPTION__

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "COneSite.h"
#include "const.h"
#include "CRec.h"

class COption
{
 private:
  double err;
  vecdouble mu;
  CRec rec;

  double N;
  double Ne;

  int nmcmc;


 public:
  COption();
  std::string fn0, fn1, fnoutput;
  std::vector<onesite> snps, rmsnps;
  std::vector<int> positions, rmpositions;
  std::vector<std::string> includedsamples;
  std::vector<std::string> excludedsamples;

  void GetData(void);
  int minpos, maxpos, intv;

  void GetIncludedSamples(const char* fn);
  void GetExcludedSamples(const char* fn);

  int  checkpos(const onesite &e) const;
  bool RMSNPs(const onesite &e) const;
  bool RMSample(const std::string &s) const;

  double GetGenotypeError(void) const;
  const vecdouble* GetMus(void) const;
  const char* GetfnOutput(void) const;
  double GetNe(void) const;
  double GetN(void) const;
  int GetNMCMC(void) const;
};

#endif
