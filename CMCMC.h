#ifndef __CMCMC__
#define  __CMCMC__

#include "CData.h"
#include "COption.h"
#include "CHMMOneHap.h"
#include "CHMM.h"
#include <vector>
#include "const.h"
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */


class CMCMC
{
 private:
  std::vector<CHMM*> hmms; //hmm for each haplotype
  
  const CData *hdata;
  CData *dlk;
  const COption *op;
  vecdouble pos_cm; //position
  const CRec *geneticmap;
  double N,  Ne;


 public:
  CMCMC(const CData *idata, const COption *ioption, const CRec *irec, CData *ilk);
  ~CMCMC();
  
  void OneSite(int s);
};

#endif
