#include "alpos/functions/APDFQ0_HERAStyle.h"

#include <iostream>
//#include "TMathBase.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_HERAStyle::fRequirements = {"xp","iPDF", "fs", "gB","gC", "uvB", "uvC", "uvE", "dvB", "dvC", "UbarB", "UbarC", "DbarA", "DbarB", "DbarC"}; //< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_HERAStyle::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_HERAStyle::fFunctionName = "PDFQ0_HERA_10pts"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_HERAStyle::APDFQ0_HERAStyle(const std::string& name) : AParmFuncBase<double>(name) { 

}


// ______________________________________________________________________________________ //
APDFQ0_HERAStyle::~APDFQ0_HERAStyle() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_HERAStyle::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_HERAStyle::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_HERAStyle::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_HERAStyle::GetQuick(const vector<double>& ipdf_xp_q0) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Advantage: Do not recalculate sum-rules
   //!   ::GetQuick(vector<double> ipdf_xp);
   
   int ipdf = int(ipdf_xp_q0[0]);
   double xp = ipdf_xp_q0[1];
   //double q0 = ipdf_xp_q0[2]; // q0 is ignored

   vector<double> ret(1);
   if(ipdf== 0) 
      ret[0] = DefaultHERAParam(xp,fgA,PAR(gB),PAR(gC));
   else if(ipdf== 1) 
      ret[0] =  DefaultHERAParam(xp,fdvA,PAR(dvB),PAR(dvC));
   else if(ipdf== 2) 
      ret[0] =  DefaultHERAParam(xp,fuvA,PAR(uvB),PAR(uvC),0,PAR(uvE));
   else if(ipdf== 3) 
      ret[0] =  2*PAR(fs)*DefaultHERAParam(xp,PAR(DbarA),PAR(DbarB),PAR(DbarC));
   else if(ipdf== 4) 
      ret[0] =  DefaultHERAParam(xp, fUbarA ,PAR(UbarB),PAR(UbarC));
   else if(ipdf== 5) 
      ret[0] = DefaultHERAParam(xp,PAR(DbarA),PAR(DbarB),PAR(DbarC));
   else if(ipdf== 6) ret[0] = 0;
   else ret[0] = 0;

   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_HERAStyle::Update() {

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
	 // 0., 0., -1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., // 
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
      
      double xp = PAR(xp);
      //double q0 = PAR(Q0);

      fgA = CalcGluonASumRule();
      fdvA = Get_dvA();
      fuvA = Get_uvA();
      fUbarA = Get_UbarA();

      debug["Update"]<<"fgA="<<fgA<<"\tfdvA="<<fdvA<<"\tfuvA="<<fuvA<<"\tfUbarA="<<fUbarA<<endl;

      if(ipdf== 0) {
	 double gA = fgA;
	 fValue[0] = DefaultHERAParam(xp,gA,PAR(gB),PAR(gC));
      }
      else if(ipdf== 1) {
	 double dvA = fdvA;
	 fValue[0] =  DefaultHERAParam(xp,dvA,PAR(dvB),PAR(dvC));
      }
      else if(ipdf== 2) {
	 double uvA = fuvA;
	 fValue[0] =  DefaultHERAParam(xp,uvA,PAR(uvB),PAR(uvC),0,PAR(uvE));
      }
      else if(ipdf== 3) {
	 double fs = PAR(fs);
	 fValue[0] =  2*fs*DefaultHERAParam(xp,PAR(DbarA),PAR(DbarB),PAR(DbarC));
      }
      else if(ipdf== 4) {
	 double UbarA = fUbarA;
	 fValue[0] =  DefaultHERAParam(xp, UbarA ,PAR(UbarB),PAR(UbarC));
      }
      else if(ipdf== 5) fValue[0] = DefaultHERAParam(xp,PAR(DbarA),PAR(DbarB),PAR(DbarC));
      else if(ipdf== 6) fValue[0] = 0;
      else fValue[0] = 0;
      
   }
   return true;
}



// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::Get_UbarA(){
   //! return recent UbarA
   return PAR(DbarA)*(1.-PAR(fs));// /(1.-0);
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::Get_dvA(){
   //! return recent dvA
   //! sum rule: D - Dbar = 1 -> dvA
   return 1./GetIntegralDefaultHERAParam(0,PAR(dvB),PAR(dvC));
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::Get_uvA(){
   //! return recent uvA
   //! sum rule: U - Ubar = 2  -> uvA
   return 2./GetIntegralDefaultHERAParam(0,PAR(uvB),PAR(uvC),0,PAR(uvE));
}


// __________________________________________________________________________________________ // 
TF1& APDFQ0_HERAStyle::GetTF1(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   //! return suitable TF1
   if ( AP==0 &&  BP==0 &&  CP== 0 ){
      if ( D==0 && E==0 && F==0 )
	 return fTF1Def3;
      else 
	 return fTF1Def6;
   }
   else return fTF1Def9;
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::GetIntegralXHERA(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   // calculate integral of x*pdf for default parameterisation
   debug["GetIntegralXHERA"]<<"A="<<A<<"\tB="<<B<<"\tC="<<C<<"\tD="<<D<<"\tE="<<E<<"\tF="<<F<<endl;

   double sum = CalcIntegral(B,C);
   if ( D!= 0 )  sum += D*CalcIntegral(B+1,C);
   if ( E!= 0 )  sum += E*CalcIntegral(B+2,C);
   if ( F!= 0 )  sum += F*CalcIntegral(B+3,C);
   if ( AP!= 0 ) sum += AP*CalcIntegral(B+4,C);
   if ( BP!= 0 ) sum += BP*CalcIntegral(B+5,C);
   if ( CP!= 0 ) sum += CP*CalcIntegral(B+6,C);
   debug["GetIntegralXHERA"]<<"sum "<<sum<<"\t1/sum="<<1./sum<<"\t2/sum="<<2./sum<<endl;
   return sum;
   
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::GetIntegralDefaultHERAParam(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   //! Calculate integral of 'DefaultHERAParam'
   debug["GetIntegralDefaultHERAParam"]<<"A="<<A<<"\tB="<<B<<"\tC="<<C<<"\tD="<<D<<"\tE="<<E<<"\tF="<<F<<endl;

   double sum = CalcIntegral(B-1,C);
   if ( D!= 0 )  sum += D*CalcIntegral(B+1-1,C);
   if ( E!= 0 )  sum += E*CalcIntegral(B+2-1,C);
   if ( F!= 0 )  sum += F*CalcIntegral(B+3-1,C);
   if ( AP!= 0 ) sum += AP*CalcIntegral(B+4-1,C);
   if ( BP!= 0 ) sum += BP*CalcIntegral(B+5-1,C);
   if ( CP!= 0 ) sum += CP*CalcIntegral(B+6-1,C);
   debug["GetIntegralDefaultHERAParam"]<<"sum:  "<<sum<<"\t1/sum="<<1./sum<<"\t2/sum="<<2./sum<<endl;
   return sum;

   /*
   TF1& f = GetTF1(A,B,C,D,E,F,AP,BP,CP);
   // https://root.cern.ch/drupal/content/function-integration 
   f.SetParameter(0,A);
   f.SetParameter(1,B);
   f.SetParameter(2,C);
   if ( D!= 0 && E!=0 && F!=0 )  {
      f.SetParameter(3,D);
      f.SetParameter(4,E);
      f.SetParameter(5,F);
   }
   if ( AP!= 0 && BP!=0 && CP!=0 )  {
      f.SetParameter(6,AP);
      f.SetParameter(7,BP);
      f.SetParameter(8,CP);
   }

   ROOT::Math::WrappedTF1 wf1(f);
 
   // Create the Integrator
   ROOT::Math::GaussLegendreIntegrator ig;
   
   const double IntPrec = 1.e-5;
   const double nPts = 40;
   const double xmin = 1.e-7;

   // Set parameters of the integration
   ig.SetFunction(wf1);
   ig.SetRelTolerance(IntPrec);
   ig.SetNumberPoints(nPts);
 
   double integral = ig.Integral(xmin, 1);
   cout << "Integral: "<< integral << endl;
   return integral;
   */
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::CalcIntegral(double alpha, double beta){
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
double APDFQ0_HERAStyle::CalcGluonASumRule(){
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

   const double tStr   = 0;
   const double tPho   = 0;
   const double UbarA  = Get_UbarA();
   const double tUbar  = UbarA*GetIntegralXHERA(UbarA ,PAR(UbarB),PAR(UbarC));
   const double tDbar  = PAR(DbarA)*GetIntegralXHERA(PAR(DbarA),PAR(DbarB),PAR(DbarC));
   const double tSea   = 2.*(tUbar+tDbar+tStr+tPho);
      
   const double uvA = Get_uvA();
   const double tUv = uvA * GetIntegralXHERA(uvA,PAR(uvB),PAR(uvC),0,PAR(uvE));
   const double dvA = Get_dvA();
   const double tDv = dvA * GetIntegralXHERA(dvA,PAR(dvB),PAR(dvC));
   // cout<<"sumrules......\t"<<tUv<<"\t"<<tDv<<endl;
   // cout<<"sumrules.bar..\t"<<tUbar<<"\t"<<tDbar<<endl;

   double tNoGlue = 1. - (tUv + tDv + tSea + tPho);
   double tg = GetIntegralXHERA(0,PAR(gB),PAR(gC));
   // cout<<"sumrules.nglu.\t"<<tNoGlue<<endl;
   // cout<<"sumrules..tg..\t"<<tg<<endl;
   if( tg<=0)tg=0.001;
   double gA = (tNoGlue)/tg;
   return gA;
}


// __________________________________________________________________________________________ // 
double APDFQ0_HERAStyle::DefaultHERAParam(double x, double A,double B,double C,double D,double E,double F, double AP,double BP,double CP) {
   //! standard-like parameterization as used in HERAFitter: 
   //! (following is taken from HERAFitter) 
   //! C  AF = (a*x**b)*(1 - x)**c*(1 + d*x + e*x**2+f*x**3)- 
   //! C     - (ap*x**bp)*(1-x)**cp
   double xf = A*pow(x,B) * pow(1-x,C) * (1 + D*x + E*pow(x,2)+ F*pow(x,3)) - AP*pow(x,BP)*pow(1-x,CP);
   // cout<<"StandardHera, x="<<x<<", A="<<A<<", B="<<B<<", C="<<C<<", E="<<E<<", F="<<F<<"\txf: "<<xf<<endl;                                                                                                
   return xf;
}

// ______________________________________________________________________________________ //
