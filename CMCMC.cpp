#include "CMCMC.h"


CMCMC ::   CMCMC(const CData *ihdata,
		 const COption *iop,
		 const CRec *irec,
		 CData *idlk)
{
  geneticmap = irec;
  hdata = ihdata;
  op = iop;
  dlk = idlk;

  int nhap = hdata->GetNHap();
  int nsite = hdata->GetNSite();
  hmms.clear();
  //  printf("starting prepare HMM\n");

  vecdouble rec;
  N = hdata->GetNState();
  Ne = op->GetNe();

  //  std::cout << "N, Ne = " << N << " " << Ne << "\n";

  double t = hdata->GetPos_cm(0);
  for (int i = 1; i != nsite; ++i)
    {
      double t1 = hdata->GetPos_cm(i);
      double r = recombination(fabs(t1 - t));
      rec.push_back(r);
    }
      
  for (int i = 0; i != nhap; ++i)
    {
      //      if (i != 100) continue;
      CHMM *hmm = new CHMM;
      //      const vecint *hid = hdata->GetPanelID(i);

      hmm->Init(hdata,
		hdata->GetPanelID(i),
		hdata->GetHap(i), 		
		op->GetGenotypeError(), 
		rec);      
      hmms.push_back(hmm);
      
      if (i % 100 == 0) printf("%d/%d\n", i, nhap);
    }

  time_t begin, end;
  time(&begin);
  
  for (int s0 = 0; s0 != idlk->GetNSite(); ++s0)        
    {
      OneSite(s0);
      if (s0 % 10 == 0) 
	{
	  time(&end);
	  std::cout << "s/nsite = " << s0 << "/" << idlk->GetNSite() <<  " - " <<  end - begin << "\n";
	}
    }
}

CMCMC :: ~CMCMC()
{
  for (int i = 0; i != (int) hmms.size(); ++i) delete hmms[i];
  hmms.clear();
}



void CMCMC :: OneSite(int s0)
{
  int pos = dlk->GetPos(s0);
  int s   = hdata->GetLocation(pos);
  if (s == -3) return;//found exact location 
 
  double geneticpos = dlk->GetPos_cm(s0);
  double recleft = 0, recright = 0; 
  int nsite = hdata->GetNSite();

  if (s == -1)  recleft  = recombination(fabs(geneticpos - hdata->GetPos_cm( 0 )));
  if (s == -2)  recright = recombination(fabs(geneticpos - hdata->GetPos_cm( nsite-1 )));

  if (s >= 0 && s + 1 < nsite)
    {
      recleft  =  recombination(fabs(geneticpos - hdata->GetPos_cm(s)));
      recright = recombination(fabs(geneticpos - hdata->GetPos_cm(s+1)));      
    }

  int nmcmc = op->GetNMCMC();
  const std::vector<vecint> *panelid = hdata->GetPanelIDs();

  vecdouble prob0A, prob0B;
  int nc = 1, k = 0;
  int step = 0;

  int steps[]={50, 5, 2};
  double minp[] = {.25, 1e-3, 0};
   
  std::cout << "s = " << s<< "\n";
  for (int ii = 0; ii != 3; ++ii)
    {
      vecint sids;
      dlk->GetPhasingSamples(s0, sids, minp[ii]);//get samples at site s0 that p(het) <= minp[ii] --> sids

      step = 0;
      k = 0;
      while (k++ < nmcmc && step < steps[ii])
	{
	  nc = 0;
	  for (int ii = 0; ii != (int)sids.size(); ++ii)
	    {
	      int sid = sids[ii];
	      int hid = 2*sid;
	      //	      std::cout << "hid = " << hid << "\n";
	      dlk->GetProb0(s0, panelid->at(hid), prob0A);
	      double a = hmms[hid]->GetProb(s, recleft, recright, prob0A); //p
	      
	      //	      std::cout << "hid = " << hid+1 << "\n";
	      dlk->GetProb0(s0, panelid->at(hid + 1), prob0B);
	      double b = hmms[hid + 1]->GetProb(s, recleft, recright, prob0B); //p

	      //	      std::cout << "sid " << sid << "\n";
	      
	      if (dlk->Update(s0, sid, a, b))	    
		{
		  nc ++;
		}
	    }
	  if (nc == 0) step ++;
	  else step = 0;
	}
    }
}


