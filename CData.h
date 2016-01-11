
#ifndef __CDATA__
#define __CDATA__

//read data

#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::binary_search, std::sort
#include <vector> 
#include <gzstream/gzstream.h>

#include "COneSite.h" 
#include "COption.h"
#include <math.h>

#include <iostream>

#define split2vector(bf,a,det){ a.clear(); int l = strlen(bf); for (int ii = 0; ii < l; ++ii) if (bf[ii] > 32 && bf[ii] != det &&  (ii == 0 || bf[ii-1] == det)) a.push_back( &(bf[ii])); for (int ii = 0; ii < l; ++ii) if (bf[ii] == det) bf[ii] = '\0';}



class CData
{
  const COption *op;
  std::vector<std::string> samplenames;
  std::vector<onesite> data;  
  int lktype;

  int nstate;

  char *bf;
  std::vector<char*>a, b, c;

  int nstat;
  std::vector<vecint> panel_id;
  std::vector<vecint> dist;

  std::vector<std::string> genes;

  std::string chromosome;


 public:

  CData (const COption *inop);
  ~CData();
  void GetHaplotype(const char* fn_vcf);
  void GetLK(const char* fn_vcf);
  void MatchData(std::vector<std::string> samples);
  std::vector<std::string> GetSamples(void);
  const std::vector<onesite> *GetHap(void) const;  
  void GetHaplotypeID(int nstat);


  void Estimate_Panel_ID(int nstat);

  void GetPanel(int i, vecstring &panel) const;

  int GetNHap(void) const;
  int GetNSite(void) const;
  int GetPos(int s) const;
  int GetLocation(int pos) const;
  const vecdouble * GetLK(int s) const;

  std::string GetHap(int s, const vecint &pid) const ; //get haplotypes for selected haps of haplotype i
  std::string GetHap(int i) const;

  double GetPos_cm(int s) const;

  const std::vector<vecint> * GetPanelIDs(void) const;
  bool  Update(int s, int sid, double a, double b);//update site i. return true if changed

  void Output(const char* fn);
  void GetPhasingSamples(int s, vecint &sids);
  bool CheckRun(int s) const;

  void RandomHap(int s);
  void GetProb0(int s, const vecint &pid, vecdouble &prob0) const;

  void ViewFinal(int s) const;


  void GetProbLeft(int s, int i, const vecint &pid,  vecdouble &probleft);

  void CM_Location(const CRec &rec);
  int GetNState(void) const;
  void GetPhasingSamples(int s, vecint &sids, double pvalue) const;

  const vecint * GetPanelID(int i) const; //get panel id for haplotype i


};

#endif

