#include "CHMM.h"

void CHMM :: Init(const CData  *ihdata,
		  const vecint *ihid,
		  const std::string iv,		  
		  double igenotype_err,
		  const vecdouble &irec )
{
  genotype_err = igenotype_err;
  hid = ihid;

  hdata = ihdata;
  v = iv;
  nsite = (int) hdata->GetNSite();
  nhap  = (int) hid->size();
  for (int i = 0; i != (int) irec.size(); ++i) rec.push_back(irec[i]);

  N = nhap;
  view = false;

  working();
}


//p0: likelihood at site s. need to calculate
//p1: likelihood at site s+1. used to calculate p0
//pp: panel probability at s+1
//rec0: prob of no recombination (1-exp(-4*r*Ne/N))/N   + exp(-4*r*Ne/N)
//rec1: (1-exp(-4*r*Ne/N))/N

double CHMM::Backward(const vecdouble&pp1, 
		      double rec, 
		      vecdouble &p0, const vecdouble &p1)
{

  double sumlk;
  sumprod(p1, pp1, sumlk);
  p0.resize(nhap);
  double rec0 = TranNoRec(rec);
  for (int i = 0; i != nhap; ++i)
    p0[i] = (sumlk-p1[i]*pp1[i])*rec + rec0*p1[i]*pp1[i];
  return 1;
}


double CHMM :: Backward(void)
{
  vecdouble lk, backward, pp1;
  backwards.resize(nsite, backward);
  backwardnorm.resize(nsite, 1);

  backward.resize(nhap, 1);

  backwards[nsite-1] = backward;
  backwardnorm[nsite-1] = Norm(backwards[nsite-1]);
  backward = backwards[nsite-1];


  vecdouble *p0;
  const vecdouble *p1;

  for (int s = nsite - 1; s -1 >= 0; --s)
    {
      if ((nsite-1-s) % 2 == 0) { p1 = &backward; p0 = &lk;}
      else { p1 = &lk; p0 = &backward;}
      GetMatchProb(s, pp1);
      Backward(pp1, rec[s-1], *p0, *p1);
      backwardnorm[s-1] = Norm((*p0));
      if ((s-1) % scale == 0 || s-1 == 0) backwards[s-1] = (*p0);           
    }


  //  std::cout << "Finish forward\n";

  return 1;
 
}

double CHMM :: Backward(int st, int ed, vecdouble &bw)
{
  vecdouble lk, pp1;
  bw = backwards[ed];
  int p01 = 0;
  for (int s = ed; s - 1 >= st; --s)
    {
      vecdouble *p0 = &bw;
      const vecdouble *p1 = &lk;
      if (p01 == 0) { p0 = &lk; p1 = &bw;}      
      p01 = 1 - p01;
      GetMatchProb(s, pp1);
      Backward(pp1, rec[s-1], *p0, *p1);
      Norm(*p0);
    }
  if (p01 == 1) bw = lk;
  return 1;
}


double CHMM :: Forward(int st, int ed, vecdouble &fw)
{
  vecdouble lk, pp1;
  fw = forwards[st];  
  int p01 = 0;
  for (int s = st; s + 1 <= ed; ++s)
    {
      const vecdouble *p0 = &fw;
      vecdouble *p1 = &lk;
      if (p01 == 1) { p0 = &lk; p1 = &fw;}      
      p01 = 1 - p01;
      GetMatchProb(s+1, pp1);
      Forward(pp1, rec[s], *p0, *p1);
      Norm(*p1);
    }
  if (p01 == 1) fw = lk;
  return 1;
}

double CHMM :: Forward(void)
{
  vecdouble forward;
  forwards.resize(nsite, forward);
  forwardnorm.resize(nsite, 1) ;

  vecdouble pp1, lk;  
  GetMatchProb(0, forward);
  forwards[0] = forward;
  forwardnorm[0] = Norm(forwards[0]);

  forward = forwards[0];


  const vecdouble  *p0;
  vecdouble  *p1;
  
  for (int s = 0; s+1 < nsite; ++s) 
    {   
      if (s % 2 == 0) { p0 = &forward; p1 = &lk;}
      else { p0 = &lk; p1 = &forward;}

      GetMatchProb(s+1, pp1);
      Forward(pp1, rec[s], *p0, *p1);
  
      forwardnorm[s+1] = Norm(*p1);
      if ((s+1) % scale == 0 || s+1 == nsite-1)  forwards[s+1] = (*p1);
    }

  //  std::cout << "Finish forward\n";

  return 0;
  
}


double CHMM:: working(void)
{
  Forward();
  Backward();

  if (CHECK)
    {
      double sumlk = LogLK(nsite-1);
      //      std::cout << "sumlk = " << sumlk << "\n";
      int nc = 0;
      for (int s = 0; s != nsite-1; ++s)
	if (forwards[s].size() > 0 && backwards[s].size() > 0)
	  {
	    if (fabs(sumlk - LogLK(s)) > 1e-6)
	      {
		std::cout << "error S = " << sumlk << " " << LogLK(s) << "\n";
		exit(1);
	      }
	    nc ++;
	  }
      if (nc == 0)
	{
	  std::cout << "nc == " << nc << "\n";
	  exit(1);
	}

      vecdouble fw;
      Forward(0,10,fw);
    }
  return 0;
}

double CHMM :: Forward(int s)
{
  vecdouble pp1;
  GetMatchProb(s+1, pp1);
  Forward(pp1, rec[s], forwards[s], forwards[s+1]);
  forwardnorm[s+1] = Norm(forwards[s+1]); 
  return 0;
}


double CHMM :: Backward(int s)
{
  vecdouble pp1;
  GetMatchProb(s,pp1);
  Backward(pp1, rec[s-1], backwards[s-1], backwards[s]);
  backwardnorm[s-1] = Norm(backwards[s-1]);
  return 0;
}


double CHMM :: Norm(vecdouble &lk)
{
  double sumlk ;
  sumvector(lk, sumlk);
  for (int i = 0; i != (int) lk.size(); ++i) lk[i] /= sumlk;
  return sumlk;
}


double CHMM :: Norm(vecdouble &lk, const vecdouble &weight)
{
  double sumlk;
  sumprod(lk,weight, sumlk);
  for (int i = 0; i != (int) lk.size(); ++i) lk[i] /= sumlk;
  return sumlk;
}



//  ---- #forward ----
//p0: lk at s; availble
//p1: lk at s+1; need to estimate
//pp1: panel match at site s+1: prob(hap[s][i] == v[s+1]
//rec: recombination between s, s+1
double CHMM:: Forward(const vecdouble&pp1, double rec, const vecdouble &p0, vecdouble &p1)
{ 
  double sumlk;
  sumvector(p0, sumlk);
  double rec0 = TranNoRec(rec); //calculate the probability of no recombination
  p1.resize(nhap);
  for (int i = 0; i != nhap; ++i) 
    p1[i] = ((sumlk-p0[i])*rec + p0[i]*rec0)*pp1[i];
  //  for (int i = 0; i != nhap; ++i)    if (p1[i] < 1e-8) p1[i] = 1e-8;
  return 1;
}

double CHMM :: Prepare(int s, double recleft, double recright)
{
  if (s  > -1 && s < nsite-1)
    {
      int left = s, right = s + 1;
      while (forwards[left].size() == 0) left--;
      while (backwards[right].size() == 0) right++;
      Forward(left, s, FW);

      vecdouble bw;
      Backward(s+1, right, bw);

      vecdouble pp1;
      GetMatchProb(s+1, pp1); //match prob at s+1
      Backward(pp1, recright,  BW, bw); 
      S = s;
    }
  
  return 0;
}

//[s - recdist - k - ...  -s+1]
//prob0: prob(hi = '0' at k
//probleft: prob at site s+1: prob(hap[s+1][i] == v[s+1])
double CHMM :: GetProb(int s, 
		       double recleft,
		       double recright,
		       const vecdouble &prob0) 
{

  //  std::cout << "s,recleft,recright,prob0: " << s << " " << recleft << " " << recright << " " << prob0.size() << "\n";


  vecdouble prob1; //prob for v = 1
  for (int i = 0; i != (int) prob0.size(); ++i) prob1.push_back(1-prob0[i]);
  vecdouble pp1;

  double p0=1, p1=1;
  if (s == -1)  //left end
    { 
      vecdouble lk;
      GetMatchProb(0,pp1); 
      Backward(pp1, recleft, lk, backwards[0]);
      sumprod(lk, prob0, p0);
      sumprod(lk, prob1, p1);
      
    }

  if (s == -2) // right end
    {
      vecdouble lk;
      Forward(prob0, recright, forwards[nsite-1], lk);
      sumvector(lk, p0);
      Forward(prob1, recright, forwards[nsite-1], lk);
      sumvector(lk, p1);
    }

  
  if (s > -1 && s < nsite-1)
    {
      if (S != s) Prepare(s, recleft, recright);

      vecdouble lk0, lk1;
      //moving forward from s to [pos] with 0 and 1 value
      Forward (prob0,     recleft,   FW,  lk0);
      //      std::cout << "ok 0\n";
      Forward (prob1,     recleft,   FW,  lk1);

      sumprod(lk0, BW, p0);
      sumprod(lk1, BW, p1);
    }

  return (p0*prior0/(p0*prior0+p1*prior1));
}


CHMM:: CHMM()
{
  genotype_err = 1e-3;
  prior0 = .5;
  prior1 = 1 - prior0;  
  CHECK = true;
  scale = 10;
  S = -1;
};


double CHMM :: LogLK(int site)
{
  if (forwards[site].size() == 0 || backwards[site].size() ==0)
    {
      std::cout << "call the site = " << site << "\n";
      exit(1);
    }

  double sumlk = 0;
  for (int s = 0; s <= site; ++s) sumlk += log(forwardnorm[s]);
  for (int s = site; s <= nsite; ++s) sumlk += log(backwardnorm[s]);
  double tt = 0;
  for (int i = 0; i != nhap; ++i) tt += forwards[site][i] * backwards[site][i];
  sumlk += log(tt);
  return sumlk;
}
