#include "alpos/functions/APDFQ0_HERA.h"

#include <iostream>
//#include "TMathBase.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_HERA::fRequirements = {"xp","iPDF", 
							     "fs", 
							     "gA", "gB","gC", "gD","gE", "gF","gAP", "gBP","gCP", 
							     "uvA", "uvB", "uvC", "uvD", "uvE", "uvF", 
							     "dvA", "dvB", "dvC", "dvD", "dvE", "dvF", 
							     "UbarA", "UbarB", "UbarC", "UbarD", "UbarE", 
							     "DbarA", "DbarB", "DbarC","DbarD", "DbarE", 
							     "seaA","seaB","seaC","seaD","seaE"}; //< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_HERA::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_HERA::fFunctionName = "PDFQ0_HERA"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_HERA::APDFQ0_HERA(const std::string& name) : AParmFuncBase<double>(name) { 

}


// ______________________________________________________________________________________ //
APDFQ0_HERA::~APDFQ0_HERA() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_HERA::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   
   fPar_g.resize(10);
   fPar_uv.resize(10);
   fPar_dv.resize(10);
   fPar_ub.resize(10);
   fPar_db.resize(10);
   fPar_sea.resize(10);

   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_HERA::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_HERA::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_HERA::GetQuick(const vector<double>& ipdf_xp) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Advantage: Do not recalculate sum-rules
   //!   ::GetQuick(vector<double> ipdf_xp);

   int ipdf = int(ipdf_xp[0]);
   double xp = ipdf_xp[1];

   vector<double> ret(1);

   if(ipdf== 0) 
      ret[0] = DefaultHERAParam(xp,fPar_g);
   else if(ipdf== 1) 
      ret[0] =  DefaultHERAParam(xp,fPar_dv);
   else if(ipdf== 2) 
      ret[0] =  DefaultHERAParam(xp,fPar_uv);
   else if(ipdf== 3) {
      double fs = PAR(fs);
      ret[0] =  2*fs*DefaultHERAParam(xp,fPar_db);
   }
   else if(ipdf== 4) 
      ret[0] =  DefaultHERAParam(xp, fPar_ub);
   else if(ipdf== 5) 
      ret[0] = DefaultHERAParam(xp,fPar_db);
   else if(ipdf== 6) 
      ret[0] = DefaultHERAParam(xp,fPar_sea);
   else ret[0] = 0;

   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_HERA::Update() {

   int ipdf = PAR(iPDF);
   debug["Update"]<<"ipdf="<<ipdf<<"\tx="<<PAR(xp)<<"\tCHECK(ipdf)="<<CHECK(iPDF)<<"\tCHECK(x)="<<CHECK(xp)<<endl;
   
   if ( ipdf == -1 ) {
      // return QCDNUM vector 'def' (gluon is expected to be [0])
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
      //double q0 = PAR(Q0);

      // Par_uv
      fPar_uv.resize(10);
      fPar_uv[1] = PAR(uvB);
      fPar_uv[2] = PAR(uvC);
      fPar_uv[3] = PAR(uvD);
      fPar_uv[4] = PAR(uvE);
      fPar_uv[5] = PAR(uvF);
      double uvA = PAR(uvA);
      fPar_uv[0] = uvA==0 ? Get_uvA() : uvA;
      // Par_dv
      fPar_dv.resize(10);
      fPar_dv[1] = PAR(dvB);
      fPar_dv[2] = PAR(dvC);
      fPar_dv[3] = PAR(dvD);
      fPar_dv[4] = PAR(dvE);
      fPar_dv[5] = PAR(dvF);
      double dvA = PAR(dvA);
      fPar_dv[0] = dvA==0 ? Get_dvA() : dvA;
      // Par_ub
      fPar_ub.resize(10);
      fPar_ub[1] = PAR(UbarB);
      fPar_ub[2] = PAR(UbarC);
      fPar_ub[3] = PAR(UbarD);
      fPar_ub[4] = PAR(UbarE);
      double ubA = PAR(UbarA);
      fPar_ub[0] = ubA==0 ? Get_UbarA() : ubA;
      // Par_db
      fPar_db.resize(10);
      fPar_db[0] = PAR(DbarA);
      fPar_db[1] = PAR(DbarB);
      fPar_db[2] = PAR(DbarC);
      fPar_db[3] = PAR(DbarD);
      fPar_db[4] = PAR(DbarE);
      // Par_sea
      fPar_sea.resize(10);
      fPar_sea[0] = PAR(seaA);
      fPar_sea[1] = PAR(seaB);
      fPar_sea[2] = PAR(seaC);
      fPar_sea[3] = PAR(seaD);
      fPar_sea[4] = PAR(seaE);
      // Par_g
      fPar_g.resize(10);
      fPar_g[1] = PAR(gB);
      fPar_g[2] = PAR(gC);
      fPar_g[3] = PAR(gD);
      fPar_g[4] = PAR(gE);
      fPar_g[5] = PAR(gF);
      fPar_g[6] = PAR(gAP);
      fPar_g[7] = PAR(gBP);
      fPar_g[8] = PAR(gCP);
      double gA = PAR(gA);
      fPar_g[0] = gA==0 ? CalcGluonASumRule() : gA;

      // debug["Update"]<<"fgA="<<fPar_g[0]<<"\tfdvA="<<fPar_dv[0]<<"\tfuvA="<<fPar_uv[0]<<"\tfUbarA="<<fPar_ub[0]<<\t"dbarA="<<fPar_db[0]<<endl;
      debug["Update"]<<"uv:  "<<fPar_uv[0]<<"  \t"<<fPar_uv[1]<<"  \t"<<fPar_uv[2]<<"  \t"<<fPar_uv[3]<<"  \t"<<fPar_uv[4]
		     <<"  \t" <<fPar_uv[5]<<"  \t"<<fPar_uv[6]<<"  \t"<<fPar_uv[7]<<"  \t"<<fPar_uv[8]<<endl;
      debug["Update"]<<"dv:  "<<fPar_dv[0]<<"  \t"<<fPar_dv[1]<<"  \t"<<fPar_dv[2]<<"  \t"<<fPar_dv[3]<<"  \t"<<fPar_dv[4]
		     <<"  \t" <<fPar_dv[5]<<"  \t"<<fPar_dv[6]<<"  \t"<<fPar_dv[7]<<"  \t"<<fPar_dv[8]<<endl;
      debug["Update"]<<"Ub:  "<<fPar_ub[0]<<"  \t"<<fPar_ub[1]<<"  \t"<<fPar_ub[2]<<"  \t"<<fPar_ub[3]<<"  \t"<<fPar_ub[4]
		     <<"  \t" <<fPar_ub[5]<<"  \t"<<fPar_ub[6]<<"  \t"<<fPar_ub[7]<<"  \t"<<fPar_ub[8]<<endl;
      debug["Update"]<<"Db:  "<<fPar_db[0]<<"  \t"<<fPar_db[1]<<"  \t"<<fPar_db[2]<<"  \t"<<fPar_db[3]<<"  \t"<<fPar_db[4]
		     <<"  \t" <<fPar_db[5]<<"  \t"<<fPar_db[6]<<"  \t"<<fPar_db[7]<<"  \t"<<fPar_db[8]<<endl;
      debug["Update"]<<"GL:  "<<fPar_g[0]<<"  \t"<<fPar_g[1]<<"  \t"<<fPar_g[2]<<"  \t"<<fPar_g[3]<<"  \t"<<fPar_g[4]
		     <<"  \t" <<fPar_g[5]<<"  \t"<<fPar_g[6]<<"  \t"<<fPar_g[7]<<"  \t"<<fPar_g[8]<<endl;
      debug["Update"]<<"ST:  "<<fPar_sea[0]<<"  \t"<<fPar_sea[1]<<"  \t"<<fPar_sea[2]<<"  \t"<<fPar_sea[3]<<"  \t"<<fPar_sea[4]
		     <<"  \t" <<fPar_sea[5]<<"  \t"<<fPar_sea[6]<<"  \t"<<fPar_sea[7]<<"  \t"<<fPar_sea[8]<<endl;
      if(ipdf== 0) 
	 fValue[0] = DefaultHERAParam(xp,fPar_g);
      else if(ipdf== 1) 
	 fValue[0] =  DefaultHERAParam(xp,fPar_dv);
      else if(ipdf== 2) 
	 fValue[0] =  DefaultHERAParam(xp,fPar_uv);
      else if(ipdf== 3) {
	 double fs = PAR(fs);
	 fValue[0] =  2*fs*DefaultHERAParam(xp,fPar_db);
      }
      else if(ipdf== 4) 
	 fValue[0] =  DefaultHERAParam(xp, fPar_ub);
      else if(ipdf== 5) 
	 fValue[0] = DefaultHERAParam(xp,fPar_db);
      else if(ipdf== 6) 
	 fValue[0] = DefaultHERAParam(xp,fPar_sea);
      else fValue[0] = 0;
      
   }
   return true;
}



// __________________________________________________________________________________________ // 
double APDFQ0_HERA::Get_UbarA(){
   //! return recent UbarA
   return PAR(DbarA)*(1.-PAR(fs));// /(1.-0);
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::Get_dvA(){
   //! return recent dvA
   //! sum rule: D - Dbar = 1 -> dvA
   return 1./GetIntegralDefaultHERAParam(0,PAR(dvB),PAR(dvC));
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::Get_uvA(){
   //! return recent uvA
   //! sum rule: U - Ubar = 2  -> uvA
   return 2./GetIntegralDefaultHERAParam(0,PAR(uvB),PAR(uvC),0,PAR(uvE));
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::GetIntegralXHERA(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   // calculate integral of x*pdf for default parameterisation
   debug["GetIntegralXHERA"]<<"A="<<A<<"\tB="<<B<<"\tC="<<C<<"\tD="<<D<<"\tE="<<E<<"\tF="<<F<<endl;

   double sum = CalcIntegral(B,C);
   if ( D!= 0 )  sum += D*CalcIntegral(B+1,C);
   if ( E!= 0 )  sum += E*CalcIntegral(B+2,C);
   if ( F!= 0 )  sum += F*CalcIntegral(B+3,C);
   if ( AP!= 0 ) sum += AP*CalcIntegral(B+4,C);
   if ( BP!= 0 ) sum += BP*CalcIntegral(B+5,C);
   if ( CP!= 0 ) sum += CP*CalcIntegral(B+6,C);
   debug["GetIntegralDefaultHERAParam"]<<"sum: "<<sum<<"\t1/sum="<<1./sum<<"\t2/sum="<<2./sum<<endl;
   return sum;
   
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::GetIntegralDefaultHERAParam(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   //! Calculate integral of 'DefaultHERAParam'
   debug["GetIntegralDefaultHERAParam"]<<"A="<<A<<"\tB="<<B<<"\tC="<<C<<"\tD="<<D<<"\tE="<<E<<"\tF="<<F<<endl;

   double sum = CalcIntegral(B-1,C);
   if ( D!= 0 )  sum += D*CalcIntegral(B+1-1,C);
   if ( E!= 0 )  sum += E*CalcIntegral(B+2-1,C);
   if ( F!= 0 )  sum += F*CalcIntegral(B+3-1,C);
   if ( AP!= 0 ) sum += AP*CalcIntegral(B+4-1,C);
   if ( BP!= 0 ) sum += BP*CalcIntegral(B+5-1,C);
   if ( CP!= 0 ) sum += CP*CalcIntegral(B+6-1,C);
   debug["GetIntegralDefaultHERAParam"]<<"sum: "<<sum<<"\t1/sum="<<1./sum<<"\t2/sum="<<2./sum<<endl;
   return sum;

}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::CalcIntegral(double alpha, double beta){
   //! Calculates int_0^1 dx x^(alpha) (1-x)^(beta) 
   //! Requires alpha > -1 and beta > -1
   //! Adapted from HERAFitter
   const double eps = 1e-5;
   double aa = alpha + 1;
   double bb = beta + 1;
   if (aa<=0.) aa = eps;
   if (bb<=0.) bb = eps;
   double u = TMath::Gamma(aa);
   double v = TMath::Gamma(bb);
   double uv = TMath::Gamma(aa+bb);
   return u*v / uv;
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::CalcGluonASumRule(){
   //! calculate gluon normalization from sum rule
   // HERAFitter does:
   //    gA = (tNoGlue+parglue(7)*tgMRST)/tg
   // with
   //    if (FlexibleGluon) 
   //         tgMRST=CalcIntegral(parglue(8),parglue(9))
   //         tg = CalcIntXpdfFixN(parglue,7) 
   //    else
   //         tg = CalcIntXpdf(parglue)
   //         tgMRST = 0
   //    tNoGlue = 1.D0 - ( tUv + tDv + tSea + tPho)
   //    tUv  = paruval(1)*CalcIntXpdf(paruval)
   //    tDv  = pardval(1)*CalcIntXpdf(pardval) 
   //    tUb  = parubar(1)*CalcIntXpdf(parubar)
   //    tDb  = pardbar(1)*CalcIntXpdf(pardbar)
   //    tSea = 2.0d0 * (tUb + tDb + tStr + tPho)
   //    tStr = 0
   //    tPho = 0
   //    


   const double tStr   = (fPar_sea[1]!=0) ? fPar_sea[0] * GetIntegralXHERA(fPar_sea) : 0;
   const double tPho   = 0;
   const double tUbar  = fPar_ub[0]*GetIntegralXHERA(fPar_ub);
   const double tDbar  = fPar_db[0]*GetIntegralXHERA(fPar_db);
   const double tSea   = 2.*(tUbar+tDbar+tStr+tPho);
   const double tUv = fPar_uv[0] * GetIntegralXHERA(fPar_uv);
   const double tDv = fPar_dv[0] * GetIntegralXHERA(fPar_dv);

   double tNoGlue = 1. - (tUv + tDv + tSea + tPho); // 1.-(val+sea) [momentum]

   double tg = 0, tgMRST=0;
   const bool FlexGlue = true;
   if ( FlexGlue ) {
      tg = GetIntegralXHERA(fPar_g[0],fPar_g[1],fPar_g[2],fPar_g[3],fPar_g[4],fPar_g[5],fPar_g[6]);
      tgMRST = CalcIntegral(fPar_g[7],fPar_g[8]);
   }
   else {
      tg = GetIntegralXHERA(fPar_g);
   }
   
   // cout<<"sumrules......\t"<<tUv<<"\t"<<tDv<<endl;
   // cout<<"sumrules.bar..\t"<<tUbar<<"\t"<<tDbar<<endl;
   // cout<<"sumrules.nglu.\t"<<tNoGlue<<endl;
   // cout<<"sumrules..tg..\t"<<tg<<"  \t"<<tgMRST<<endl;
   if( tg<=0)tg=0.001;
   double gA = (tNoGlue+fPar_g[6]*tgMRST)/tg;
   return gA;
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERA::DefaultHERAParam(double x, double A,double B,double C,double D,double E,double F, double AP,double BP,double CP) {
   //! standard-like parameterization as used in HERAFitter: 
   //! (following is taken from HERAFitter) 
   //! C  AF = (a*x**b)*(1 - x)**c*(1 + d*x + e*x**2+f*x**3)- 
   //! C     - (ap*x**bp)*(1-x)**cp
   double xf = A*pow(x,B) * pow(1-x,C) * (1 + D*x + E*pow(x,2)+ F*pow(x,3)) - AP*pow(x,BP)*pow(1-x,CP);
   // cout<<"StandardHera, x="<<x<<", A="<<A<<", B="<<B<<", C="<<C<<", E="<<E<<", F="<<F<<"\txf: "<<xf<<endl;                                                                                                
   return xf;
}

// ______________________________________________________________________________________ //
