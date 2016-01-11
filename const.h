#ifndef __CONST__
#define __CONST__

#include <string>
#include <vector>


#define split2vector(bf,a,det){ a.clear(); int l = strlen(bf); for (int ii = 0; ii < l; ++ii) if (bf[ii] > 32 && bf[ii] != det &&  (ii == 0 || bf[ii-1] == det)) a.push_back( &(bf[ii])); for (int ii = 0; ii < l; ++ii) if (bf[ii] == det) bf[ii] = '\0';}




#define sumvector(a,sumlk){ sumlk = 0; for (int iii = 0; iii != (int) a.size(); ++iii) sumlk += a[iii];}
#define sumprod(a,b,sumlk){ sumlk = 0; for (int ii = 0; ii != (int) a.size(); ++ii) sumlk += (a[ii]*b[ii]);}

typedef std::vector<double> vecdouble;
typedef std::vector<std::string> vecstring;
typedef std::vector<int> vecint;

#define recombination(l) ((1-exp( -4*(l)*Ne/N)) / N)
#define TranNoRec(rec) (1-(N-1)*rec) //tranmission of no recombination: hap[i] --> hap[i]



#endif 
