#include "CData.h"
#include "COption.h"

CData::CData (const COption *inop)
{
  op = inop;
  bf = (char*)calloc(100000000, sizeof(char));
  genes.clear();
  genes.push_back("00");
  genes.push_back("10");
  genes.push_back("01");
  genes.push_back("11");
}


CData :: ~CData()
{
  free(bf);
}


const std::vector<onesite>* CData :: GetHap(void) const
{
  return (&data);
}

  

void CData :: GetHaplotype(const char* fn_vcf)
{
  fprintf(stderr, "loading %s\n", fn_vcf);

  bool first = false;
  if (data.size() == 0) first = true;

   igzstream in(fn_vcf);
   std::string line;
   while (getline(in,line) && (strstr(line.c_str(), "#CHROM") != line.c_str()));
   strcpy(bf, line.c_str());
   split2vector(bf,a,'\t');
   for (int i = 9; i != (int) (a.size()); ++i) samplenames.push_back(a[i]);

   fprintf(stderr, "load %d\n", (int) samplenames.size());

   while (getline(in,line))
     {
       strcpy(bf, line.c_str());
       split2vector(bf,a,'\t');
	  
       onesite e;
       e.pos = atoi(a[1]);
       e.al1 = a[3];
       e.al2 = a[4];

       if (op->checkpos(e) == -2) continue;
       if (op->checkpos(e) == 2) break;
       if (op->RMSNPs(e)) continue;

       e.hap = "";
       for (int i = 9; i != (int) a.size(); ++i)
	 {
	   e.hap += a[i][0];
	   e.hap += a[i][2];
	 }

       if (first) data.push_back(e);
       else
	 {
	   std::vector<onesite>::iterator low = std::lower_bound(data.begin(), data.end(), e);
	   if (*low == e) (*low).hap += e.hap;	   
	 }
    }
	 

  if (first) sort(data.begin(), data.end());
  fprintf(stderr, "Haplotype data load: %d samples, %d variants\n", (int) samplenames.size(), (int) data.size());

}


const vecdouble * CData :: GetLK(int s) const { return &data[s].lk;}

//-1: left
//-2: right
//-3: match
//other [s, pos, s+1]
int CData ::GetLocation(int pos) const
{
  
  if (pos < data[0].pos) return -1;
  if (pos > data[ data.size()-1].pos) return -2; //
  
  int st = 0, step = 1;
  while (step > 0)
    {
      int  i = st + step;
      if (i < (int) data.size() && data[i].pos < pos) { step *= 2; st = i;}
      else
	if (i < (int) data.size() && data[i].pos == pos) return -3; 
	else step /= 2;
    }
  return st;
}




void CData :: GetLK(const char* fn_vcf)
{
  
  vecdouble phred;
  for (int i = 0; i != 10000; ++i) phred.push_back(exp(log(10.0)* (((double) -i)/10)));

  #define MAXLEN 100000000

  igzstream in(fn_vcf);
  std::string line;
  char* bf = (char*)calloc(MAXLEN, sizeof(char));
  
  while (getline(in,line) && (strstr(line.c_str(), "#CHROM") != line.c_str()));
  
  strcpy(bf, line.c_str());
  split2vector(bf, a, '\t');
  for (int i = 9; i != (int) a.size(); ++i)        samplenames.push_back(a[i]);
  fprintf(stderr, "loading %d samples\n", (int) samplenames.size());
  

  chromosome = "";
  
  while (getline(in,line))
    {
      strcpy(bf, line.c_str());
      split2vector(bf, a, '\t');
      if (chromosome == "") chromosome = a[0];
      onesite e;
      e.pos = atoi(a[1]);
      e.rsid = a[2];
      e.al1 = a[3];
      e.al2 = a[4];
     
     
      if (strstr(e.al1.c_str(), ",") != NULL) continue; //multiple allelel variants
      if (strstr(e.al2.c_str(), ",") != NULL) continue; //multiple allelel variants
      if (op->checkpos(e) == -1) continue;
      if (op->checkpos(e) == 1) break;
      
      
      split2vector(a[8],b,':');
      int id = -1;
      lktype = -1;
      for (int i = 0; i != (int) b.size(); ++i)
	if (strcmp(b[i], "PL") == 0)
	  {
	    id = i;
	    lktype = 1;
	    break;
	  }
      
      for (int i = 0; i != (int) b.size(); ++i)
	if (strcmp(b[i], "GL") == 0)
	  {
	     id = i;
	     lktype = 0;
	  }
      
      if (lktype == -1)
	{
	  fprintf(stderr, "no GL or PL at %s", line.c_str());
	  continue;
	}

      int gtid = -1;
      for (int i = 0; i != (int) b.size(); ++i)
	if (strcmp(b[i], "GT") == 0)  gtid = i;


      e.lk.clear();
      e.hap = "";
      for (int i = 9; i != (int) a.size(); ++i)
	{
	  split2vector(a[i],b,':');
	  double p00 = 1;
	  double p01 = 1;
	  double p11 = 1;
	  if (id < (int) b.size())
	    {
	      split2vector(b[id], c, ',');
	      vecint v;
	      for (int kk = 0; kk != 3; ++kk) v.push_back(atoi(c[kk]));
	      for (int kk = 0; kk != 3; ++kk)
		if (v[kk] > 500) v[kk] = 500;

	      p00 =  phred[ v[0]];
	      p01 =  phred[ v[1]];
	      p11 =  phred[ v[2]];
	    }

	  double sum = p00+p01+p11;
	  e.lk.push_back(p00/sum);
	  e.lk.push_back(p01/sum);
	  e.lk.push_back(p11/sum);

	  if (gtid != -1 && gtid < (int) b.size()) 	      
	    {
	      e.hap += b[gtid][0];
	      e.hap += b[gtid][2]; //init haplotype
	    }
	  else
	    {
	      if (p00 + 0.01 >= p01 && p00 + 0.01 >= p11) e.hap += "00";
	      else
		if (p01 >= p00 && p01 >= p11) e.hap += "01";
		else e.hap += "11";
	    }
	}
      data.push_back(e);
    }
 
 
   std::vector<int> kid; //samples kept
   for (int i = 0; i != (int) samplenames.size(); ++i)
     if (!op->RMSample(samplenames[i])) kid.push_back(i);
   

   for (int s = 0; s != (int)data.size(); ++s)
     for (int i = 0; i != (int) samplenames.size(); ++i)
       if (data[s].hap[i*2] == '1' && data[s].hap[i*2+1] == '1' && data[s].lk[i*3] > .75)
	 {
	   std::cout << "error here " << data[s].pos << " " << s << " " << samplenames[i] << " " << " i = " << i << "\n";	   
	   for (int k = 0; k != 3; ++k) std::cout << " " << data[s].lk[i*3+k];
	   std::cout << "\n";
	   exit(1);
	 }


   if (kid.size() != samplenames.size())
     {   
       for (int i = 0; i != (int) kid.size(); ++i) samplenames[i] = samplenames[kid[i]];
       samplenames.erase(samplenames.begin() + kid.size(), samplenames.end());

       int nsample = (int)samplenames.size();   
       for (int s = 0; s != (int) data.size(); ++s)
	 {
	   for (int j = 0; j != (int) kid.size(); ++j)
	     {
	       for (int k = 0; k != 3;  ++k) data[s].lk[3*j+k] = data[s].lk[ kid[j]*3+k];
	       for (int k = 0; k != 2;  ++k) data[s].hap[j*2 + k] = data[s].hap[ kid[j]*2 + k];
	     }
	   data[s].lk.erase(data[s].lk.begin() + nsample*3, data[s].lk.end());
	   data[s].hap.erase(nsample*2, data[s].hap.size());
	 }
     }

   for (int s = 0; s != (int) data.size(); ++s) data[s].pplk.resize(samplenames.size()*4, 0);
   fprintf(stderr, "Likelihood data after filtering : %d samples, %d variants, lk size = %d\n", (int) data[0].hap.size(), (int) data.size(), (int) data[0].lk.size());


   
   for (int s = 0; s != (int)data.size(); ++s)
     for (int i = 0; i != (int) samplenames.size(); ++i)
       if (data[s].hap[i*2] == '1' && data[s].hap[i*2+1] == '1' && data[s].lk[i*3]  > .75)
	 {
	   std::cout << "error here: filter " << data[s].pos << " " << s << " " << samplenames[i] << " " << " i = " << i << "\n";	   
	   for (int k = 0; k != 3; ++k) std::cout << " " << data[s].lk[i*3+k];
	   std::cout << "\n";
	   exit(1);   
	 }
   

   for (int s = 0; s!= (int)data.size(); ++s)
     {
       data[s].prob0.resize(data[s].hap.size(), 1.0/3.0);
       for (int i = 0; i != (int) data[s].hap.size()/2; ++i)
	 {
	   data[s].prob0[i*2]   = data[s].lk[i*3] + 0.5* data[s].lk[i*3+1];
	   data[s].prob0[i*2+1] = data[s].lk[i*3] + 0.5* data[s].lk[i*3+1];
	 }
     }
     

   for (int s = 0; s != (int) data.size(); ++s)     data[s].inputhap = data[s].hap;
   for (int s = 0; s != (int) data.size(); ++s)     
     for (int i = 0; i != (int) data[s].hap.size(); ++i)
       if (data[s].hap[i] == '.') data[s].hap[i] = '0';


   //   std::cout << "s == 0" << data[0].prob0.size() << "\n";
   //   std::cout << "s == 0" << data[0].hap << "\n";


   //check

}


void CData :: ViewFinal(int s) const
{
  printf("site = %d\n", s);
  for (int i = 0; i + 1 < (int) data[s].hap.size(); i += 2)
    //    if (data[s].hap[i] != data[s].inputhap[i] || data[s].hap[i+1] != data[s].inputhap[i+1])
    if (data[s].inputhap[i] != data[s].inputhap[i+1])      
      {
	std:: cout << data[s].inputhap[i] << data[s].inputhap[i+1] << " -- ";
	std:: cout << data[s].hap[i] << data[s].hap[i+1] << " -- ";
	for (int k = 0; k != 3; ++k) std::cout << " " << data[s].lk[(i/2)*3 + k];
	std::cout << "; pplk = " ;
	for (int k = 0; k != 4; ++k) std::cout << " " << data[s].pplk[(i/2)*4 + k];
	std::cout << "\n";
      }
  

}

std::vector<std::string> CData :: GetSamples(void)
{
  return samplenames;
}


void CData :: MatchData(std::vector<std::string> samples)
{
  int l = 0;
  for (int i = 0; i != (int) data.size(); ++i)
    if (data[i].hap.size() == 2*samplenames.size()) data[l++] = data[i];
 
  fprintf(stderr, "remove %d sites\n", (int)data.size() - l);

  data.erase(data.begin() + l, data.end());

  std::vector<int> kid;
  for (int i = 0; i != (int) samples.size(); ++i)
    {
      int ok = 0;
      for (int j = 0; j != (int) samplenames.size(); ++j)
	if (samplenames[j] == samples[i]) 
	  {
	    kid.push_back(j);
	    ok = 1;
	    break;
	  }
      if (ok == 0)
	{
	  printf("no data for sample %s\n", samples[i].c_str());
	  exit(1);
	}
    }
  samplenames = samples;
  for (int i = 0; i != (int) data.size(); ++i)
    {
      std::string hap = "";
      
      for (int j = 0; j != (int) kid.size(); ++j)
	{
	  hap += data[i].hap[kid[j]*2];
	  hap += data[i].hap[kid[j]*2 + 1];
	}
      data[i].hap = hap;
    }

  fprintf(stderr, "The haplotype data match perfectly: %d %d\n", (int) data.size(), (int) data[0].hap.size());

}


const vecint *CData :: GetPanelID(int i) const { return &panel_id[i];}



std::string CData :: GetHap(int i) const
{ 
  std::string v = "";
  for (int s = 0; s != (int) data.size(); ++s) v = v + data[s].hap[i];
  return v;
}


void CData :: Estimate_Panel_ID(int instate)
{
  nstate = instate;
  int nhap = data[0].hap.size();
  int nsite = (int) data.size();
  printf("Estimate panels: %d haplotypes %d variants\n", nhap, nsite);

  vecint e;
  e.resize(nhap, 0);
  dist.resize(nhap, e);

  for (int s = 0; s != nsite; ++s)
    {
      vecint id0, id1;
      for (int i = 0; i != nhap; ++i)
	if (data[s].hap[i] == '0') id0.push_back(i);
	else id1.push_back(i);

      for (int i = 0; i != (int) id0.size(); ++i)
	for (int j  = 0; j != (int) id1.size(); ++j)
	  {
	    dist[id0[i]] [id1[j]] ++;
	    dist[id1[j]] [id0[i]] ++;
	  }
      
      if (s % 1000 == 0) printf("%d / %d\n", s, nsite);
    }
  e.clear();
  panel_id.resize(nhap,e);

  for (int i = 0; i != nhap; ++i)
    {
      std::vector< std::pair<int, int > > q;
      q.clear();
      for (int j  = 0; j != nhap; ++j)
	if ((i/2) != (j / 2))	  q.push_back( std::make_pair(dist[i][j], j));

      sort(q.begin(), q.end());
      for (int j = 0; j != nstate; ++j)	panel_id[i].push_back( q[j]. second);
      sort(panel_id[i].begin(), panel_id[i].end());
    }
  


}



int CData :: GetNSite(void) const { return (int) data.size();}
int CData :: GetNHap(void) const { return (int) data[0].hap.size();}
int CData :: GetPos(int s) const { return data[s].pos;}
double CData :: GetPos_cm(int s) const { return data[s].pos_cm;}



void CData :: GetProb0(int s, const vecint &pid, vecdouble &prob0) const
{
  //  std::cout << "pid = " << pid.size() << "\n";
  //  std::cout << "data[s].prob0 = " << data[s].prob0.size() << "\n";
  //  for (int i = 0; i != (int) pid.size(); ++i)       std::cout << pid[i] << " ";
  //  std::cout << "\n";
  

  prob0.resize(pid.size(),0);
  for (int i = 0; i != (int) pid.size(); ++i)       prob0[i] =  data[s].prob0 [ pid[i] ];
}

std::string CData :: GetHap(int s, const vecint &pid) const
{
  std::string h = "";
  for (int j = 0; j != (int) pid.size(); ++j) h += data[s].hap[ pid[j] ];
  return h;
}


const std::vector<vecint> * CData :: GetPanelIDs(void) const
{
  return &panel_id;
}

bool CData :: CheckRun(int s) const
{
  for (int i = 0; i != (int) data[s].hap.size(); ++i)
    if (data[s].hap[i] != data[s].hap[0]) return true;
  return false;
}


bool CData :: Update(int s, int sid, double a, double b)
{
  int ppid =  sid*4; //location for posterior propability
  int did  =  sid*3; //location for data from vcf: P(d|00), p(d|01), p(d|11)

  
  data[s].pplk[ppid  ] =  data[s].lk[did]      *a     *b; //00
  data[s].pplk[ppid+1] =  data[s].lk[did+1]    *(1-a) *b; //10
  data[s].pplk[ppid+2] =  data[s].lk[did+1]    *1     *(1-b); //01
  data[s].pplk[ppid+3] =  data[s].lk[did+2]    *(1-a) *(1-b); //11

  double sumlk = 0;
  for (int i = 0; i != 4; ++i) sumlk += data[s].pplk[ppid+i];
  for (int i = 0; i != 4; ++i) data[s].pplk[ppid+i] /= sumlk; 
   

  data[s].prob0[sid*2]   = data[s].pplk[ppid] + data[s].pplk[ppid+2]; // 00 + 01
  data[s].prob0[sid*2+1] = data[s].pplk[ppid] + data[s].pplk[ppid+1]; // 00 + 10
  
  int ml = 0;
  for (int l = 0; l != 4; ++l)
    if (data[s].pplk[ppid+l] > data[s].pplk[ppid+ml]) ml = l;
  

  bool change = false;
  if (data[s].hap[sid*2] != genes[ml][0] || data[s].hap[sid*2+1] != genes[ml][1])        change = true;

  if (1==2)
    //  if (data[s].hap[sid*2] != data[s].hap[sid*2+1] || change)
    if (change)
    {
      printf("\nsid = %d %c%c --> %s\n", sid, data[s].hap[sid*2],data[s].hap[sid*2+1], genes[ml].c_str());
      std::cout << "lk = " << data[s].lk[did] << " " <<  data[s].lk[did+1]  << " " <<  data[s].lk[did+2] << "\n";
      for (int i = 0; i != 4; ++i) std::cout <<  " " <<  data[s].pplk[ppid+i];
      printf("\n");
    }

  data[s].hap[sid*2] = genes[ml][0];
  data[s].hap[sid*2+1] = genes[ml][1];
  return change;  
}



void CData :: Output(const char* fn)
{
  ogzstream op (fn);
  
  op << "##fileformat=VCFv4.2\n";
  op << "##INFO=<ID=GambiaCP3Variants,Number=0,Type=Flag,Description=\"Novel variants found in Gambia and CP3\">\n";  op << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  op << "##FORMAT=<ID=PQ,Number=4,Type=Integer,Description=\"posterior probability: 00, 10, 01, and 11 \"\n";
  op << "##FORMAT=<FID=PQ,Number=4,Type=Integer,Description=\"posterior probability: 00, 10, 01, and 11 \"\n";
  
  op << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTqER\tINFO\tFORMAT";
  for (int i = 0; i != (int) samplenames.size(); ++i) op << "\t" + samplenames[i];
  op << "\n";
  
  for (int s = 0; s != (int) data.size(); ++s)
    {
      op << chromosome << "\t" << data[s].pos << "\t" << data[s].rsid;
      op << "\t" << data[s].al1 << "\t" << data[s].al2;
      op << "\t.\t.\tGambiaCP3Variants";
      op << "\tGT:PQ";

      for (int i = 0; i != (int) samplenames.size(); ++i) 
	{
	  op << "\t";
	  op << data[s].hap[i*2]<<"|"<<data[s].hap[i*2+1] << ":";
	  op << data[s].pplk[i*4];
	  for (int j = 1; j != 4;  ++j) op<< "," << data[s].pplk[i*4+j];
	}
      op << "\n";     
    }
}


void CData :: RandomHap(int s)
{
  for (int i = 0; i != (int) data[s].hap.size(); ++i)
    if (rand() % 2 == 0) data[s].hap[i] = '0';
    else data[s].hap[i] = '1';
}

void CData :: GetPhasingSamples(int s, vecint &sids, double pvalue) const
{
  sids.clear();
  for (int i = 0; i != (int) samplenames.size(); ++i)
    if (data[s].lk[i*3+1]  > pvalue) sids.push_back(i);

}


int CData :: GetNState(void) const { return nstate;}
void CData :: CM_Location(const CRec &rec)
{
  for (int s = 0; s != (int) data.size(); ++s)  data[s].pos_cm = rec.GetPos(data[s].pos);
  //  for (int s = 0;s  != (int) 10; ++s)      std::cout << data[s].pos << " " << data[s].pos_cm << "\n";
}

