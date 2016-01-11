#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>


#include "CData.h"
#include "COption.h"
#include "CMCMC.h"
#include <math.h>


int main(int ArgI, char* ArgC[])
{
  COption op;
  
  CRec geneticmap;

  
  for (int i = 0; i != ArgI-1; ++i)
    {
      if (strcmp(ArgC[i], "-interval") == 0 && i + 2 < ArgI) 
	{
	  op.minpos = atoi(ArgC[i+1]);
	  op.maxpos = atoi(ArgC[i+2]);
	}

      if (strcmp(ArgC[i], "-include-variant") == 0) 	  op.fn0 = ArgC[i+1];
      if (strcmp(ArgC[i], "-exclude-vatiant") == 0)	  op.fn1 = ArgC[i+1];
      if (strcmp(ArgC[i], "-genetic-map") == 0) geneticmap.Init(ArgC[i+1]);
      if (strcmp(ArgC[i], "-buffer") == 0) op.intv = atoi(ArgC[i+1]);
      if (strcmp(ArgC[i], "-include-sample") == 0) op.GetIncludedSamples(ArgC[i+1]);
      if (strcmp(ArgC[i], "-exclude-sample") == 0) op.GetExcludedSamples(ArgC[i+1]);
      if (strcmp(ArgC[i], "-output") == 0) op.fnoutput = ArgC[i+1];   
      
    }
	
  op.GetData();

  CData dhap (&op);
  CData dlk (&op);

  for (int i = 0; i != ArgI-1; ++i)
    {
      if (strcmp(ArgC[i], "-hap") == 0) dhap.GetHaplotype(ArgC[i+1]);
      if (strcmp(ArgC[i], "-lk") == 0) dlk.GetLK(ArgC[i+1]);
    }
  
  dhap.MatchData(dlk.GetSamples()); 
  dhap.CM_Location(geneticmap);
  dlk.CM_Location(geneticmap);

  dhap.Estimate_Panel_ID(200); 
  CMCMC mcmc(&dhap, &op, &geneticmap, &dlk);
  
  dlk.Output(op.GetfnOutput());
  return 0;
}
