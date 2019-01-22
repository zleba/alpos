#include "alpos/functions/APDFQ0_LHAPDF.h"

#include <iostream>

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_LHAPDF::fRequirements = {"PDF","xp","Q0","iPDF","PDFLiCo"}; //< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_LHAPDF::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_LHAPDF::fFunctionName = "PDFQ0_LHAPDF"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_LHAPDF::APDFQ0_LHAPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("APDFQ0_LHAPDF");
}


// ______________________________________________________________________________________ //
APDFQ0_LHAPDF::~APDFQ0_LHAPDF() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_LHAPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   fLico = PAR_S(PDFLiCo);
   if ( fLico == "HERAStyle" ) {
      fdefQcdnum = {
	  //tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
          //-6  -5  -4  -3  -2  -1  0   1   2   3   4   5   6
	    // 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., //dval
	    // 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., //uval
	    // 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s+sbar
	    // 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., //Ubar
	    // 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., // Dbar
	    // 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s-sbar
	    // 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // cval
	    // 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., //cbar

	    0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., //dval
	    0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., //uval
	    0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., //Ubar
	    0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., // Dbar
	    // 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s+sbar
	    // 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s-sbar
	    0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s+
	    0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0., //s-
	    0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c+
	    0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., // c-

	    0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b+
	    0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., // b-
	    0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t+
	   -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., // t-
      };
      fdefQcdnum.resize(13*12);
      fdef = { 0,0,0,0,0,0,1,0,0,0,0,0,0};
      AlposTools::operator+=(fdef, fdefQcdnum);
   }
   else if ( fLico == "LHAPDFdefinition" ) {
      fdef = {
	  //tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t
          //-6  -5  -4  -3  -2  -1  0   1   2   3   4   5   6
	    1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	    0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	    0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	    0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	    0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
	    0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
	    0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
	    0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
	    0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
      };
   }
   else if ( fLico == "ApfelxxDefinition" ) {
      fdef = {
	   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, // gluon
	   1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, // SIGMA = \sum_q q^+      
	  -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1, // VALENCE = \sum_q q^-    
	   0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0, // T3  = u^+ - d^+         
	   0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0, // V3  = u^- - d^-         
	   0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0, // T8  = u^+ + d^+ - 2 s^+ 
	   0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0, // V8  = u^- + d^- - 2 s^- 
	   0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0, // T15 = u^+ + d^+ + s^+ - 3 c^+            
	   0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0, // V15 = u^- + d^- + s^- - 3 c^-            
	   0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0, // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+      
	   0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0, // V24 = u^- + d^- + s^- + c^- - 4 b^-      
	  -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5, // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
	   5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5, // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-

      };
   }
   else {
      error["Init"]<<"PDFLiCo '"<<fLico<<"' is unknown."<<endl;
      exit(1);
      return false;
   }
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_LHAPDF::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_LHAPDF::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_LHAPDF::GetQuick(const vector<double>& ipdf_xp_muf) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);


   int ipdf = int(ipdf_xp_muf[0]);
   double xp = ipdf_xp_muf[1];
   double q0 = ipdf_xp_muf[2];
   vector<double> ret = {0};
   vector<double> values = QUICK(PDF, ({xp,q0}) );
   // if ( ipdf == 0 ) ret[0] = values[6]; // 0 is always gluon
   // else {
   //    //AlposTools::CalcLicoFrom13partons(LiCo, fdef, values );
   //    for ( int i = 0; i<13 ; i++ ) {
   // 	 ret[0] += fdef[ipdf*13+i]*values[i];
   //    }
   // }
   for ( int i = 0; i<13 ; i++ ) {
      ret[0] += fdef[ipdf*13+i]*values[i];
   }
   return ret;
}


// ______________________________________________________________________________________ //
bool APDFQ0_LHAPDF::Update() {

   int ipdf = PAR(iPDF);
   
   if ( ipdf == -1 ) {
      if ( fLico == "HERAStyle") {
	 // return QCDNUM vector 'def'
	 fValue = fdefQcdnum;
	 fValue.resize(13*12);
	 fError.resize(fValue.size());
      } 
      else {
	 error["Update"]<<"Linear combination '"<<fLico<<"' is requested, but this the requested program requires the gluon as zeroth component."<<endl;
	 error["Update"]<<"Please use PDF linearcombination 'HERAStyle'"<<endl;
	 return false;
      }
   }
   else if ( ipdf == -2 ) {
      // return 'default' vector 'def'
      fValue = fdef;
      fValue.resize(13*13);
      fError.resize(fValue.size());
   }
   else { 
      fValue.resize(1);
      fError.resize(1);
      PAR(PDF); // update PDF
      double xp = PAR(xp);
      double q0 = PAR(Q0);
      vector<double> values = QUICK(PDF, ({xp,q0}) );
      fValue[0] = 0;
      // if ( ipdf == 0 ) fValue[0] = values[6]; // 0 is always gluon
      // else {
      // 	 //AlposTools::CalcLicoFrom13partons(LiCo, fdef, values );
      // 	 for ( int i = 0; i<13 ; i++ ) {
      // 	    fValue[0] += fdef[ipdf*13+i]*values[i];
      // 	 }
      // }
      for ( int i = 0; i<13 ; i++ ) {
	 fValue[0] += fdef[ipdf*13+i]*values[i];
      }
      // if     (ipdf== 0) fValue[0] = values[6];
      // else if(ipdf== 1) fValue[0] = values[7]-values[5];
      // else if(ipdf== 2) fValue[0] = values[8]-values[4];
      // else if(ipdf== 3) fValue[0] = values[9]+values[3];
      // else if(ipdf== 4) fValue[0] = values[4]+values[2];
      // else if(ipdf== 5) fValue[0] = values[3]+values[5];
      // else if(ipdf== 6) fValue[0] = values[9]-values[3];
      // else fValue[0] = 0.;
   }
   return true;
}

// ______________________________________________________________________________________ //
