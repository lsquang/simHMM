#include "COption.h"


COption::COption()
{
  intv = 500000;
  minpos = -1;
  maxpos = 3000000000;
  fn0 = "";
  fn1 = "";
  err = 1e-3;
  mu.resize(100000, 1e-8);

  fn0 = "";
  fn1 = "";

  Ne = 17500;
  N = 1000;
  nmcmc = 100;

}

double COption :: GetGenotypeError(void) const { return err;}
const vecdouble *COption :: GetMus(void) const { return &mu;}


void COption::GetData(void)
{
  #define  M100 100000000
  char* al1  = (char*)calloc(M100, sizeof(char));
  char* al2  = (char*)calloc(M100, sizeof(char));
   
  onesite e;

  if (fn0 != "")
    {
      FILE *f = fopen(fn0.c_str(), "rb");
      while (fscanf(f, "%d %s %s", &e.pos, al1, al2) == 3)
	{
	  e.al1 = al1;
	  e.al2 = al2;
	  snps.push_back(e);
	  positions.push_back(e.pos);
	}
      fclose(f);
      printf("include variants: %d\n", (int) snps.size());
    }

  if (fn1 != "")
    {
      FILE *f = fopen(fn1.c_str(), "rb");
      while (fscanf(f, "%d %s %s", &e.pos, al1, al2) == 3)
	{
	  e.al1 = al1;
	  e.al2 = al2;
	  rmsnps.push_back(e);
	  rmpositions.push_back(e.pos);
	}
      fclose(f);
      printf("exclude variants: %d\n", (int) rmsnps.size());
    }

  free(al1);
  free(al2);
  std::sort(snps.begin(), snps.end());
  std::sort(rmsnps.begin(), rmsnps.end());
  std::sort(positions.begin(), positions.end());
  std::sort(rmpositions.begin(), rmpositions.end());
}



bool COption:: RMSNPs(const onesite &e) const
{
  if (std::binary_search(rmpositions.begin(), rmpositions.end(), e.pos)) return true;
  if (positions.size() > 0 && !std::binary_search(positions.begin(), positions.end(), e.pos)) return true;
  return false;
}

/*
bool COption:: RMSNPs(const onesite &e) const
{
  if (std::binary_search(rmsnps.begin(), rmsnps.end(), e)) return true;
  if (snps.size() > 0 && !std::binary_search(snps.begin(), snps.end(), e)) return true;
  return false;
}
*/
int COption:: checkpos(const onesite &e) const
{
  //[-2: maxpos -int , -1: maxpos ---- maxpos: 1, maxpos + intv : 2]
  int v = 0;

  if (e.pos < minpos-intv)  v= 2;
  if (e.pos < minpos)  v = -1;
  if (e.pos > maxpos)  v = 1;
  if (e.pos > maxpos + intv) v = 2;
  return v;  
}



void COption :: GetIncludedSamples(const char* fn)
{
  #define M100 100000000
  char *tempc = (char*)calloc(M100, sizeof(char));
  FILE *f = fopen(fn, "rb");
  while (fscanf(f, "%s", tempc) == 1) includedsamples.push_back(tempc);
  fclose(f);
  free(tempc);
  std::sort(includedsamples.begin(), includedsamples.end());
  //  for (int i = 0; i != (int) includedsamples.size(); ++i)    printf("%s\n", includedsamples[i].c_str());

}

void COption :: GetExcludedSamples(const char* fn)
{
  #define M100 100000000
  char *tempc = (char*)calloc(M100, sizeof(char));
  FILE *f = fopen(fn, "rb");
  while (fscanf(f, "%s", tempc) == 1) excludedsamples.push_back(tempc);
  fclose(f);
  free(tempc);

  std::sort(excludedsamples.begin(), excludedsamples.end());
}



bool COption :: RMSample(const std::string &s) const
{
  if (includedsamples.size() > 0 && !std::binary_search(includedsamples.begin(), includedsamples.end(), s)) return true;
  
  if (std::binary_search(excludedsamples.begin(), excludedsamples.end(), s)) return true;
  
  return false;

}




const char* COption :: GetfnOutput(void) const { return fnoutput.c_str();}


double COption :: GetN(void) const { return N;}
double COption :: GetNe(void) const { return Ne;}

int COption :: GetNMCMC(void) const {return nmcmc;}
