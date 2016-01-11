#ifndef __CHMMOneHap__
#define  __CHMMOneHap__


#include "const.h"
#include "COption.h"
#include "CHMM.h"
#include <string>
#include <string>
#include "CData.h"


class CHMMOneHap{

 private:
  CHMM hmm;
 public:  
  CHMMOneHap(const COption *op, const CData *idata, int iid);
};

#endif

