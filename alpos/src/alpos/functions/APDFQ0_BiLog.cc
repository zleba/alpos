#include "alpos/functions/APDFQ0_BiLog.h"

#include <iostream>
//#include "TMathBase.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_BiLog::fRequirements = {"xp","iPDF", 
							      "fs", 
							      "G1","G2","G3","G4","G5",
							      "uv1","uv2","uv3","uv4","uv5",
							      "dv1","dv2","dv3","dv4","dv5",
							      "ub1","ub2","ub3","ub4","ub5",
							      "db1","db2","db3","db4","db5"}; //!< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_BiLog::fStopFurtherNotification = {""}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_BiLog::fFunctionName = "PDFQ0_BiLog"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_BiLog::APDFQ0_BiLog(const std::string& name) : AParmFuncBase<double>(name) { 

}


// ______________________________________________________________________________________ //
APDFQ0_BiLog::~APDFQ0_BiLog() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_BiLog::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   
   fPar_g.resize(5);
   fPar_uv.resize(5);
   fPar_dv.resize(5);
   fPar_ub.resize(5);
   fPar_db.resize(5);

   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_BiLog::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_BiLog::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_BiLog::GetQuick(const vector<double>& ipdf_xp) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Advantage: Do not recalculate sum-rules
   //!   ::GetQuick(vector<double> ipdf_xp);

   int ipdf = int(ipdf_xp[0]);
   double xp = ipdf_xp[1];

   vector<double> ret(1);

   if(ipdf== 0) 
      ret[0] = splogn(xp,fPar_g); //ok
   else if(ipdf== 1) 
      ret[0] =  splogn(xp,fPar_dv); 
   else if(ipdf== 2) 
      ret[0] =  splogn(xp,fPar_uv);
   else if(ipdf== 3) {
      ret[0] =  2*ffs*splogn(xp,fPar_db);
   }
   else if(ipdf== 4) 
      ret[0] = splogn(xp, fPar_ub);
   else if(ipdf== 5) 
      ret[0] = splogn(xp,fPar_db);
   else ret[0] = 0;

   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_BiLog::Update() {

   int ipdf = PAR(iPDF);
   debug["Update"]<<"ipdf="<<ipdf<<"\tx="<<PAR(xp)<<"\tCHECK(ipdf)="<<CHECK(iPDF)<<"\tCHECK(x)="<<CHECK(xp)<<endl;
   
   if ( ipdf == -1 ) {
      // return QCDNUM vector 'def'
      vector<double> def = {
	 //tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t 
	 //-6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
	 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., // dval
	 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., // uval
	 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., // s+sbar
	 // 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., // Ubar 
	 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., // Ubar 
	 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., // Dbar 
	 0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,  //s-sbar 
      };
      fValue = def;
      fValue.resize(13*12);
      fError.resize(fValue.size());
    }
   else if ( ipdf == -2 ) {
      // return vector 'def' for all flavors
      vector<double> def = {
         //tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t 
         //-6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
         0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., // gluon
         0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., // dval
         0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., // uval
         0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., // s+sbar
         // 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., // Ubar 
         0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., // Ubar 
         0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., // Dbar 
         0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,  //s-sbar 
         // --- remaining linear combination
         0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c+ 
         0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c- 
         0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b+ 
         0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b- 
         0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t+ 
         -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t-
      };
      fValue = def;
      //fValue.resize(13*13);
      fError.resize(fValue.size());
   }
   else { 
      fValue.resize(1);
      fError.resize(1);
      
      // quite slow code, but prepare for 'Quick'-access
      double xp = PAR(xp);
      ffs = PAR(fs);
      //double q0 = PAR(Q0);

      fPar_uv.resize(5);
      fPar_uv[1] = PAR(uv2);
      fPar_uv[2] = PAR(uv3);
      fPar_uv[3] = PAR(uv4);
      fPar_uv[4] = PAR(uv5);
      double uvA = PAR(uv1);
      fPar_uv[0] = uvA==0 ? Get_uvA() : uvA;
      // Par_dv
      fPar_dv.resize(5);
      fPar_dv[1] = PAR(dv2);
      fPar_dv[2] = PAR(dv3);
      fPar_dv[3] = PAR(dv4);
      fPar_dv[4] = PAR(dv5);
      double dvA = PAR(dv1);
      fPar_dv[0] = dvA==0 ? Get_dvA() : dvA;
      // Par_ub
      fPar_ub.resize(5);
      fPar_ub[3] = PAR(ub4);
      fPar_ub[4] = PAR(ub5);
      double ubA = PAR(ub1);
      fPar_ub[0] = ubA==0 ? Get_UbarA() : ubA;
      fPar_ub[1] = ubA==0 ? PAR(db2) : PAR(ub2);
      fPar_ub[2] = ubA==0 ? PAR(db3) : PAR(ub3);

      // Par_db
      fPar_db.resize(5);
      fPar_db[0] = PAR(db1);
      fPar_db[1] = PAR(db2);
      fPar_db[2] = PAR(db3);
      fPar_db[3] = PAR(db4);
      fPar_db[4] = PAR(db5);
      // Par_g
      fPar_g.resize(5);
      fPar_g[1] = PAR(G2);
      fPar_g[2] = PAR(G3);
      fPar_g[3] = PAR(G4);
      fPar_g[4] = PAR(G5);
      double g1 = PAR(G1);
      fPar_g[0] = g1==0 ? CalcGluonASumRule() : g1;

      // debug["Update"]<<"fgA="<<fPar_g[0]<<"\tfdvA="<<fPar_dv[0]<<"\tfuvA="<<fPar_uv[0]<<"\tfUbarA="<<fPar_ub[0]<<\t"dbarA="<<fPar_db[0]<<endl;
      debug["Update"]<<"uv:  "<<fPar_uv[0]<<"  \t"<<fPar_uv[1]<<"  \t"<<fPar_uv[2]<<"  \t"<<fPar_uv[3]<<"  \t"<<fPar_uv[4]<<endl;
      debug["Update"]<<"dv:  "<<fPar_dv[0]<<"  \t"<<fPar_dv[1]<<"  \t"<<fPar_dv[2]<<"  \t"<<fPar_dv[3]<<"  \t"<<fPar_dv[4]<<endl;
      debug["Update"]<<"Ub:  "<<fPar_ub[0]<<"  \t"<<fPar_ub[1]<<"  \t"<<fPar_ub[2]<<"  \t"<<fPar_ub[3]<<"  \t"<<fPar_ub[4]<<endl;
      debug["Update"]<<"Db:  "<<fPar_db[0]<<"  \t"<<fPar_db[1]<<"  \t"<<fPar_db[2]<<"  \t"<<fPar_db[3]<<"  \t"<<fPar_db[4]<<endl;
      debug["Update"]<<"GL:  "<<fPar_g[0]<<"  \t"<<fPar_g[1]<<"  \t"<<fPar_g[2]<<"  \t"<<fPar_g[3]<<"  \t"<<fPar_g[4]<<endl;

       // set fValue
      fValue = GetQuick({double(ipdf),xp}); 
      
   }
   return true;
}



// __________________________________________________________________________________________ // 
double APDFQ0_BiLog::Get_UbarA(){
  //! return recent UbarA
  return PAR(db1)*(1.-PAR(fs)) /(1.-0);

}


// __________________________________________________________________________________________ // 
double APDFQ0_BiLog::Get_dvA(){
   //! return recent dvA
  return 1./SumRuleASpar(-1, PAR(dv1), PAR(dv2), PAR(dv3), PAR(dv4), PAR(dv5));
}


// __________________________________________________________________________________________ // 
double APDFQ0_BiLog::Get_uvA(){
   //! return recent uvA
  return 2./SumRuleASpar(-1, PAR(uv1), PAR(uv2), PAR(uv3), PAR(uv4), PAR(uv5));
}




// __________________________________________________________________________________________ //
double APDFQ0_BiLog::CalcGluonASumRule() {
  double sumMom = 2.*fPar_ub[0]*SumRuleASpar(0, fPar_ub) +
                  2.*fPar_db[0]*SumRuleASpar(0, fPar_db) +
                     fPar_uv[0]*SumRuleASpar(0, fPar_uv) +
                     fPar_dv[0]*SumRuleASpar(0, fPar_dv);
  double sumGlue = SumRuleASpar(0, PAR(G1),  PAR(G2),  PAR(G3),  PAR(G4),  PAR(G5));
  return (1. - sumMom)/sumGlue;
}

// __________________________________________________________________________________________ //
double APDFQ0_BiLog::SumRuleASpar(int n, double A, double B, double C, double D, double E) {


  A = 1.;
  B += n;

  //C---------------------------------------------------------------
  //*     Numerical Integration of the 
  //*     Special lognormal function 
  //*     using the Simpson method
  //*
  //*
  //*     A.Schoening, University Heidelberg, Physikalisches Institut
  //*     Creation: 12.6.2011
  //*  
  //*  ASPDF = A1*x**(A2-A3*log(x))*(1-x)**(A4-A5*log(1-x))
  //C---------------------------------------------------------------


  bool logflag = true;
  bool falling = false;
 
  double eps = 0.001;

  double f1 = splogn(eps, A, B, C, D, E);
  double f2 = splogn(0.5, A, B, C, D, E);
  double f3 = splogn(1.-eps, A, B, C, D, E);

  if(f1 > f2+f3) 
    {
      logflag = true;
      falling = true;
    }
  else if(f3 > f2+f1)
    {
      logflag = true;
      falling = false;
    }
  else 
    {
      logflag = false;
    }

  
  double h;                                    
  double sum;
  double xas, xnas;
  double splognni = 0;
  double xlmin = -10;
  int np = 1000;
  if(!logflag) // linear integration 
    {
      h = 1./np;
      sum = 0.;
      for(int i=1; i<=np-1; i+=2)
        {
          xas = i*h;
          sum += splogn(xas, A, B, C, D, E);
        }
      sum *= 2.;
      for(int i=2; i<=np-2; i+=2)
        {
          xas = i*h;
          sum += splogn(xas, A, B, C, D, E);
        }
      sum *= 2.;
      sum += splogn(0., A, B, C, D, E) + splogn(1., A, B, C, D, E);
      splognni = h/3. * sum;
    }
  else if(falling) // steeply falling distribution
    {
      h = -xlmin/np;
      sum = 0.;
      for(int i=1; i<=np-1; i+=2)
        {
          xas = pow(10, xlmin+i*h);
          sum += xas*splogn(xas, A, B, C, D, E);
        }
      sum *= 2.;
      for(int i=0; i<=np-2; i+=2)
        {
          xas = pow(10, xlmin+i*h);
          sum += xas*splogn(xas, A, B, C, D, E);
        }
      sum *= 2.;
      xas = 1.;
      sum += splogn(xas, A, B, C, D, E);
      splognni = h/3. * sum * log(10.);
    }
  else // steeply rising distribution
    {
      h = -xlmin/np;
      sum = 0.;
      for(int i=1; i<=np-1; i+=2)
      {
        xas = pow(10, xlmin+i*h);
        xnas = 1. - xas;
        sum += xas*splogn(xnas, A, B, C, D, E);
          }
      sum *= 2.;
      for(int i=0; i<=np-2; i+=2)
        {
          xas = pow(10, xlmin+i*h);
              xnas = 1. - xas;
              sum += xas*splogn(xnas, A, B, C, D, E);
        }
      sum *= 2.;
      xas = 0.;
      sum += splogn(xas, A, B, C, D, E);
      splognni = h/3. * sum * log(10.);
    }
  if(splognni == 0)
    {
      splognni = 0.1;
    }
  return splognni;
}



// __________________________________________________________________________________________ // 

// __________________________________________________________________________________________ // 

double APDFQ0_BiLog::splogn(double x, double A, double B, double C, double D, double E) {
   // C----------------------------------------------------
   // c     Special lognormal function 
   // c
   // c
   // c     A.Schoening, University Heidelberg, Physikalisches Institut
   // c     Creation: 12.6.2011
   // c     A1*x**(A2-A3*log(x))*(1-x)**(A4-A5*log(1-x))
   // c     translated to c++, 25.8.15
   // C-----------------------------------------------------

   double splgn = 0;
   if ( x>0 && x < 1) 
      //splgn = A * pow(x,(B-C*log(x))) * pow((1.-x),D-E*log(1.-x));
      splgn = A * pow(x,(B-C*log(x))) * pow((1.-x),D-E*log1p(-x));
 
   // constrain range
   if ( fabs(splgn)<1e30 && fabs(splgn)>1e-30 )
      return splgn;
   else
      return 0;
}

// ______________________________________________________________________________________ //
