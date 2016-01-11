#ifndef __CHMM__
#define  __CHMM__

#include "const.h"

#include <string>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "CRec.h"
#include "CData.h"

#define GenotypeError(h1,h2,err)((((h1)) == ((h2)))? ((1-err)):((err)))


#define GetMatchProb(s,prob) { prob.clear(); std::string hhh = hdata->GetHap(s, (*hid));   for (int ii = 0; ii != nhap; ++ii) prob.push_back(GenotypeError(hhh[ii], v[s], genotype_err));}

class CHMM
{
  bool CHECK;
  int scale;

  double prior0, prior1;
  int nhap, nsite;
  double genotype_err;
  vecdouble rec; //rec[s]: recombination between site s and s+1

  //std::vector<std::string> panel; //panel
  const CData *hdata;
  const vecint *hid; //haplotype id using for Panel

  double N;
  std::string v; //haplotype in consideration

  std::vector<vecdouble> forwards, backwards;
  vecdouble forwardnorm, backwardnorm;


  bool view;


  vecdouble FW, BW;
  int S;

 public:
  CHMM();
  double working(void);
  double GetProb(int s, 
		 double recleft,
		 double recright,
		 const vecdouble &prob0);

  double Norm(vecdouble &lk);
  double Norm(vecdouble &lk, const vecdouble &weight);
  
  void Init(const CData *hdata,
	    const vecint *ihid,
	    const std::string iv,//haplotype in considered
	    double igenotypeerr,  //genotype error; one value across all site
	    const vecdouble &irec); //recombination



  //p0: likelihood at site s. need to calculate
  //p1: likelihood at site s+1. used to calculate p0
  //pp: panel probability at s+1
  //rec: recombination between s - (s + 1)  
  double Backward(const vecdouble&pp1, double rec, vecdouble &p0, const vecdouble &p1);
  double Backward(int st, int ed, vecdouble &bw);
  double Backward(int s);
  double Backward(void);    
    
  //p0: lk at s; availble
  //p1: lk at s+1; need to estimate
  //pp1: panel match at site s+1: prob(hap[s][i] == v[s+1]
  //rec: recombination between s, s+1
  double Forward(const vecdouble&pp1, double rec, const vecdouble &p0, vecdouble &p1);

  double Forward(int st, int ed, vecdouble &fw);
  double Forward(int s);
  double Forward(void);
  double LogLK(int s);

  double Prepare(int s, double recleft, double recright);

};

#endif




