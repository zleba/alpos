#include "alpos/functions/APDFQ0_diff.h"

#include <iostream>
//#include "TMathBase.h"
#include <TMath.h>
#include <TF1.h>
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_diff::fRequirements = {"xp","iPDF", 
							     "gTF1","sTF1","vTF1",
							     "g0", "g1", "g2", "g3", "g4",
							     "s0", "s1", "s2", "s3", "s4",
							     "v0", "v1", "v2", "v3", "v4",
};//< List of all AParm's which this function depends on
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

   CONST(gTF1);
   CONST(sTF1);
   CONST(vTF1);

   { // init and check validity of TF1's

      if ( TString(PAR_S(gTF1))!="default" ) {
	 fgTF1 = TF1("g",TString(PAR_S(gTF1)),0,1);
	 fgTF1.SetParameters(PAR(g0),PAR(g1),PAR(g2),PAR(g3),PAR(g4));
	 if ( std::isnan(fgTF1.Eval(0.01)) || std::isnan(fgTF1.Eval(PAR(xp) ))  ) {
	    error["Init"]<<"Funtion gTF1 is not a valid formula for a TF1: "<< PAR_S(gTF1) <<endl;
	    exit(3);
	 }
      }

      if ( TString(PAR_S(sTF1))!="default" ) {
	 fsTF1 = TF1("s",TString(PAR_S(sTF1)),0,1);
	 fsTF1.SetParameters(PAR(s0),PAR(s1),PAR(s2),PAR(s3),PAR(s4));
	 if ( std::isnan(fsTF1.Eval(0.01)) || std::isnan(fsTF1.Eval(PAR(xp) ))  ) {
	    error["Init"]<<"Funtion sTF1 is not a valid formula for a TF1: "<< PAR_S(sTF1) <<endl;
	    exit(3);
	 }
      }

      if ( TString(PAR_S(vTF1))!="default" ) {
	 fvTF1 = TF1("v",TString(PAR_S(vTF1)),0,1);
	 fvTF1.SetParameters(PAR(v0),PAR(v1),PAR(v2),PAR(v3),PAR(v4));
	 if ( std::isnan(fvTF1.Eval(0.01)) || std::isnan(fvTF1.Eval(PAR(xp) ))  ) {
	    error["Init"]<<"Funtion vTF1 is not a valid formula for a TF1: "<< PAR_S(vTF1) <<endl;
	    exit(3);
	 }
      }
   }

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

   vector<double> ret(1);
   if(ipdf== 0) 
      ret[0] = fgTF1.IsValid() ? fgTF1.Eval(xp) : DefaultDiffParam(xp,PAR(g0),PAR(g1),PAR(g2)); // gluon
   else if(ipdf== 1) 
      ret[0] = fsTF1.IsValid() ? fsTF1.Eval(xp) : DefaultDiffParam(xp,PAR(s0),PAR(s1),PAR(s2)); // singlet
   else if(ipdf== 2) 
      ret[0] = fvTF1.IsValid() ? fvTF1.Eval(xp) : DefaultDiffParam(xp,PAR(v0),PAR(v1),PAR(v2)); //valence
   else if(ipdf>= 3) // all other cases
      ret[0] = 0;
   else ret[0] = 0;

   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_diff::Update() {

   int ipdf = PAR(iPDF);
   debug["Update"]<<"ipdf="<<ipdf<<"\tx="<<PAR(xp)<<"\tCHECK(ipdf)="<<CHECK(iPDF)<<"\tCHECK(x)="<<CHECK(xp)<<endl;
   
   if ( ipdf == -3 ) {
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
   else if ( ipdf>= 0 ) { 
      fValue.resize(1);
      fError.resize(1);
      
      double xp = PAR(xp);
      //double q0 = PAR(Q0);

      if ( fgTF1.IsValid() )
	 fgTF1.SetParameters(PAR(g0),PAR(g1),PAR(g2),PAR(g3),PAR(g4));
      if ( fsTF1.IsValid() )
	 fsTF1.SetParameters(PAR(s0),PAR(s1),PAR(s2),PAR(s3),PAR(s4));
      if ( fvTF1.IsValid() )
	 fvTF1.SetParameters(PAR(v0),PAR(v1),PAR(v2),PAR(v3),PAR(v4));

      vector<double> ipdf_xp_q0{double(ipdf),xp,0};
      
      fValue = GetQuick(ipdf_xp_q0);
      fError.resize(fValue.size());
   }
   else {
      error["Update"]<<"Wrong ipdf: "<<ipdf<<endl;
   }
   return true;
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

// ______________________________________________________________________________________ //
