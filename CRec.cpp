#include "CRec.h"


CRec::CRec()
{
  recs.clear();
}

CRec::CRec(const char* fn)
{
  
  Init(fn);
}

void CRec::Init(const char* fn)
{ 
  #define M10 10000000
  char *bf = (char*)calloc(M10, sizeof(char));
  FILE *f = fopen(fn, "rb");
  int pos;
  double rate, posCM;
  recs.clear();
  if (fgets(bf, M10, f) == NULL) return;
  while (!feof(f) && fscanf(f, "%d %lf %lf", &pos, &rate, &posCM) == 3) 
    recs.push_back( std::make_pair(pos, posCM) );
  free(bf);
  std::sort(recs.begin(), recs.end());

  printf("Reading genetic maps: %d\n", (int) recs.size());

}


double CRec :: GetPos(int pos) const
{
  int st = 0, l = 1;
  while (l > 0)
    if (st+l >= (int) recs.size() || recs[st+l].first > pos) l /= 2;
    else 
      {
	st +=  l;
	l *= 2;
      }
  if (st + 1 >= (int) recs.size()) st = recs.size() - 2;
  double out = (recs[st].second +    ((double)(pos - recs[st].first) / (double)(recs[st+1].first - recs[st].first))  	  * (recs[st+1].second - recs[st].second));
  return out;
 }

