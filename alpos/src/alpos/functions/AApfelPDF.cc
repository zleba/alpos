#include "alpos/functions/AApfel.h"
#include "APFEL/APFEL.h"
#include <iostream>


using namespace std;

const std::vector<std::string> AApfelPDF::fRequirements = {"ApfelInit","xp","muf"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AApfelPDF::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelPDF::fFunctionName = "ApfelPDF"; //< The function's name


// __________________________________________________________________________________________ //
AApfelPDF::AApfelPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AApfelPDF");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
AApfelPDF::~AApfelPDF() {
}


// ___________________________________________________________________________________________ //
bool AApfelPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.

   fq0 = PAR(ApfelInit.Q0);
   flastmuf = 0;

   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelPDF::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(2,xp,Q);
   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   double xp  = va_arg(ap, double);
   double muf = va_arg(ap, double);
   va_end(ap);

   return GetQuick({xp,muf});
    
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelPDF::GetQuick(const vector<double>& xp_muf) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> ret(fValue.size());
   if ( xp_muf.size() != 2) {
      cout<<"Error in AApfelPDF::GetQuick(vector). Quick acces is implemented for two parameter which are 'xp' and 'muf'."<<endl;
      return ret;
   }
   
   //--- slow but definitely correctly working
   // APFEL::EvolveAPFEL(fq0,xp_muf[1]);
   // for (int i = 0; i < 13; i++)  ret[i] = APFEL::xPDF(i-6,xp_muf[0]);
   // return ret;

   // //--- fast, but danger of possible bias
   // static double lastQ0=0;
   // double Q0 = fq0;//PAR(ApfelInit.Q0);
   // double muf = xp_muf[1];
   // if ( muf != lastQ0 ) 
   //    APFEL::EvolveAPFEL(Q0,muf);
   // lastQ0=muf;
   // for (int i = 0; i < 13; i++) ret[i] = APFEL::xPDF(i-6,xp_muf[0]);
   // return ret;

   //--- using new APFEL functionality
   // double muf = xp_muf[1];
   // if ( muf != APFEL::GetMuF0() ) 
   //    APFEL::EvolveAPFEL(fq0,muf);
   // for (int i = 0; i < 13; i++) ret[i] = APFEL::xPDF(i-6,xp_muf[0]);
   // return ret;


   //--- 02/19
   double muf = xp_muf[1];
   if ( muf != flastmuf ) {
      APFEL::EvolveAPFEL(fq0,muf);
      //info["Quick"]<<"Evolve! fq0="<<fq0<<", muf="<<muf<<"\tlastmuf="<<flastmuf<<endl;
      flastmuf=muf;
   }
   APFEL::xPDFall(xp_muf[0], &ret[0]);
   return ret;
   
}


// __________________________________________________________________________________________ //
bool AApfelPDF::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());
   
   UPDATE(ApfelInit); // ensure that apfel instance is up-to-date
   
   // Load evolution
   //double Q0 = PAR(ApfelInit.Q0);
   double Q  = PAR(muf);
   if ( Q != flastmuf ) {
      APFEL::EvolveAPFEL(fq0,Q);
      flastmuf = Q;
   }

   double xp = PAR(xp);

   APFEL::xPDFall(xp, &fValue[0]);
   // for (int i = 0; i < 13; i++) {
   //    fValue[i] = APFEL::xPDF(i-6,xp);
   // }

   // photon
   //APFEL::xgamma(xlha[i]) 

   return true;
}

