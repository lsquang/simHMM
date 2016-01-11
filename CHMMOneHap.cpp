#include "CHMMOneHap.h"


CHMMOneHap::CHMMOneHap(const COption *op, const CData *idata, int iid)
{
  vecstring panel;
  idata->GetPanel(iid, panel);
  std::string v = idata->GetHap(iid);
  double loglk = hmm.working(panel, op->GetMus(), op->GetError(), v);
}

