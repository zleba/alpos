#include <iostream>
#include <cmath>
#include <vector>
#include "APFEL/APFEL.h"

using namespace std;

extern "C" void externalsetapfel_(double* x, double* q, double* xf);

double funGamama(double a0, double a1, double a2, double x) {
   double eps = 1e-9;
   if(x <= eps || x >= 1-eps) return 0;
   return a0 * pow(x, a1) * pow(1-x, a2) * exp(-0.01/(1-x));
}

// function to be called by APFEL
void externalsetapfel_(double* x, double* q, double* xf) {
   double g0  =    1.0;
  // double g0  =    0.14591;
   double g1  =    0;
   double g2  =    0;
   //double g2  =   -0.94705;

   double s0  =    1.0;
   double s1  =    0.0;
   double s2  =    0.0;
   //double s0  =    1.0587;
   //double s1  =    2.2964;
   //double s2  =    0.56894;

   for(int i = 0; i < 13; ++i) {
      if(i == 6)
         xf[i] = funGamama(g0, g1, g2, *x); //gluon
      else if(abs(i-6) <= 3) 
         xf[i] = funGamama(s0, s1, s2, *x)/6; //quark
      else
         xf[i] = 0;
   }
}


double getF2(string scheme, double mu, double q, double x)
{
   int iOrd = 0;
   APFEL::SetTheory("QCD");
   int nf = 0;
   if ( nf==0 ) APFEL::SetVFNS();
   else APFEL::SetFFNS(nf);

   /// hard coded limits
   APFEL::SetQLimits(0.1,300*300); // 0.4Gev to 100TeV
   /// set grids with mostly hard-coded parameters
   int nGridPts0=150, nGridPts1 = 100, nGridPts2=60;
   if ( (nGridPts2) > 0 )  APFEL::SetNumberOfGrids(3);
   else APFEL::SetNumberOfGrids(2);
   APFEL::SetGridParameters(1,(nGridPts0),3,1.e-3);
   APFEL::SetGridParameters(2,(nGridPts1),5,1.e-1);
   if ( (nGridPts2) > 0 )  APFEL::SetGridParameters(3,(nGridPts2),5,8.e-1);

   //APFEL::SetMassScheme("ZM-VFNS");
   APFEL::SetMassScheme(scheme);
   //double mu = 0.5;
   APFEL::SetRenQRatio(mu);
   APFEL::SetFacQRatio(mu);
   APFEL::SetRenFacRatio(1.0);

   APFEL::SetAlphaQCDRef(0.118, 91.1876);
   APFEL::SetPDFEvolution("exactalpha");
   APFEL::SetFastEvolution(1);
   //APFEL::EnableLeptonEvolution(PAR(EnableLeptonEvolution));
   APFEL::SetPDFSet("external");
   APFEL::SetAlphaEvolution("exact");
   APFEL::SetMaxFlavourAlpha(5);
   APFEL::SetMaxFlavourPDFs(5);
   //APFEL::EnableMassRunning(PAR(EnableMassRunning));
   APFEL::InitializeAPFEL_DIS();
   APFEL::SetPerturbativeOrder(iOrd);

   double q0  = sqrt(1.75);

   //double x   = 0.1;
   //double q = 20;
   APFEL::EvolveAPFEL(q0, q);
   vector<double> ret(13);
   APFEL::xPDFall(x, &ret[0]);
   for(auto v : ret) cout << v <<" ";
   cout << endl;


   APFEL::ComputeStructureFunctionsAPFEL(q0,q);

   //cout << APFEL::F2total(x) << endl;
   return APFEL::F2total(x);
}








int main()
{
   vector<double> fonl(3), zmns(3);
   vector<double> mu{0.5};
   double q = sqrt(200);
   double z = 0.5;
   for(int i = 0; i < mu.size(); ++i) {
      fonl[i] = getF2("FONLL-C", mu[i], q, z);
      zmns[i] = getF2("ZM-VFNS", mu[i], q, z);
      cout <<"RADEK " <<  mu[i] <<" : "<< fonl[i] << " "<< zmns[i] << endl;
   }

   for(int i = 0; i < 3; ++i)
      cout << mu[i] <<" : "<< fonl[i] << " "<< zmns[i] << endl;

   return 0;
}
