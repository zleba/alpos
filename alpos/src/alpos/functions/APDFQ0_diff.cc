#include "alpos/functions/APDFQ0_diff.h"

#include <iostream>
//#include "TMathBase.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_diff::fRequirements = {"xp","iPDF", "Ag", "Bg", "Cg", "Aq", "Bq", "Cq"};//< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_diff::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_diff::fFunctionName = "PDFQ0_diff"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_diff::APDFQ0_diff(const std::string& name) : AParmFuncBase<double>(name) { 

}


// ______________________________________________________________________________________ //
APDFQ0_diff::~APDFQ0_diff() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_diff::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_diff::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_diff::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_diff::GetQuick(const vector<double>& ipdf_xp_q0) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Advantage: Do not recalculate sum-rules
   //!   ::GetQuick(vector<double> ipdf_xp);
   
   int ipdf = int(ipdf_xp_q0[0]);
   double xp = ipdf_xp_q0[1];
   //double q0 = ipdf_xp_q0[2]; // q0 is ignored

   double Ag = PAR(Ag);
   double Bg = PAR(Bg);
   double Cg = PAR(Cg);
   double Aq = PAR(Aq);
   double Bq = PAR(Bq);
   double Cq = PAR(Cq);

   vector<double> ret(1);
   if(ipdf== 0) 
      ret[0] =   DefaultDiffParam(xp, Ag, Bg, Cg); //gluon
   else if(ipdf== 1) 
      ret[0] =  0; // singlet
   else if(ipdf== 2) 
      ret[0] =  0; //valence
   else if(ipdf>= 3) 
      ret[0] = 0;
   else ret[0] = 0;

   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_diff::Update() {

   int ipdf = PAR(iPDF);
   debug["Update"]<<"ipdf="<<ipdf<<"\tx="<<PAR(xp)<<"\tCHECK(ipdf)="<<CHECK(iPDF)<<"\tCHECK(x)="<<CHECK(xp)<<endl;
   
   if ( ipdf == -1 ) {
      error["Update"]<<"Error. ipdf=-1 not checked, supported and implemented, please use -3"<<endl;
      exit(3);
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
      error["Update"]<<"Error. ipdf=-2 not checked, supported and implemented, please use -3"<<endl;
      exit(3);
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
   else if ( ipdf == -3 ) {
      // return QCDNUM vector 'def'
      vector<double> def = {
	 //tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t 
	 //-6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
	 0., 0., 0., 1., 1., 1., 0., 1., 1., 1., 0., 0., 0., // light singlet
	 0., 0., 0., -1., -1., -1., 0., 1., 1., 1., 0., 0., 0., // valence
         0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0, // T3  = u^+ - d^+
         0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0, // V3  = u^- - d^-
         0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0, // T8  = u^+ + d^+ - 2 s^+
         0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0, // V8  = u^- + d^- - 2 s^-
         0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c+ 
         0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c- 
         0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b+ 
	 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b- 
	 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t+ 
	 -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t-
	 // 0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0, // T3  = u^+ - d^+
	 // 0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0, // V3  = u^- - d^-
         // 0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0, // T8  = u^+ + d^+ - 2 s^+
         // 0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0, // V8  = u^- + d^- - 2 s^-
         // 0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0, // T15 = u^+ + d^+ + s^+ - 3 c^+
         // 0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0, // V15 = u^- + d^- + s^- - 3 c^-
         // 0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0, // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+
         // 0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0, // V24 = u^- + d^- + s^- + c^- - 4 b^-
         // 5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5, // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
         // 5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5, // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-
      };
      fValue = def;
      fValue.resize(13*12,0);
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

      vector<double> ipdf_xp_q0{double(ipdf),xp,0};
      
      fValue = GetQuick(ipdf_xp_q0);
      fError.resize(fValue.size());
   }
   return true;
}



// __________________________________________________________________________________________ // 
double APDFQ0_diff::Get_UbarA(){
   //! return recent UbarA
   return 0;//PAR(DbarA)*(1.-PAR(fs));// /(1.-0);
}


// __________________________________________________________________________________________ // 
double APDFQ0_diff::Get_dvA(){
   //! return recent dvA
   //! sum rule: D - Dbar = 1 -> dvA
   return 0;//1./GetIntegralDefaultHERAParam(0,PAR(dvB),PAR(dvC));
}


// __________________________________________________________________________________________ // 
double APDFQ0_diff::Get_uvA(){
   //! return recent uvA
   //! sum rule: U - Ubar = 2  -> uvA
   return 0;//2./GetIntegralDefaultHERAParam(0,PAR(uvB),PAR(uvC),0,PAR(uvE));
}



// __________________________________________________________________________________________ // 
TF1& APDFQ0_diff::GetTF1(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   //! return suitable TF1
   if ( AP==0 &&  BP==0 &&  CP== 0 ){
      if ( D==0 && E==0 && F==0 )
	 return fTF1Def3;
      else 
	 return fTF1Def6;
   }
   else return fTF1Def9;
}


//not needed
// __________________________________________________________________________________________ // 
double APDFQ0_diff::GetIntegralXHERA(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
   return 0;
}


//not needed
// __________________________________________________________________________________________ // 
double APDFQ0_diff::GetIntegralDefaultHERAParam(double A, double B, double C,double D,double E,double F,double AP,double BP,double CP){
    return 0;
}

// __________________________________________________________________________________________ // 
//Not needed
double APDFQ0_diff::CalcIntegral(double alpha, double beta){
   return 0;
}

// __________________________________________________________________________________________ // 

//Not needed
double APDFQ0_diff::CalcGluonASumRule(){
   return 0;
}


// __________________________________________________________________________________________ // 

double APDFQ0_diff::DefaultDiffParam(double x, double A, double B, double C) {
   //! standard-like parameterization as used in HERAFitter: 
   //! (following is taken from HERAFitter) 
   //! C  AF = (a*x**b)*(1 - x)**c*(1 + d*x + e*x**2+f*x**3)- 
   //! C     - (ap*x**bp)*(1-x)**cp
   double xf = A*pow(x,B) * pow(1-x,C);
   return xf;
}



//Not needed
double APDFQ0_diff::DefaultHERAParam(double x, double A,double B,double C,double D,double E,double F, double AP,double BP,double CP) {
   return 0;
}

// ______________________________________________________________________________________ //
