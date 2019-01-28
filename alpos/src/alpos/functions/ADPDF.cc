
#include "alpos/functions/ADPDF.h"
#include "apfel/dglapbuilder.h"

#include <iostream>


//! 
//! A small function, to put together
//! different
//! 

using namespace std;

const std::vector<std::string> ADPDF::fRequirements = {  
   "xpom", // xpom (dummy)
   "zpom", // zpom (dummy)
   "muf", // mu_f (dummy)
   "pom1", // pomeron 1, this must be a 'PDF'
   "pom2", // pomeron 2, this must be a 'PDF'
   "reg1", // reggeon 1, this must be a 'PDF'
}; //< List of all AParm's which this function depends on
const std::vector<std::string> ADPDF::fStopFurtherNotification = {"xpom","zpom","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ADPDF::fFunctionName = "ApfelxxPDF"; //< The function's name


// __________________________________________________________________________________________ //
ADPDF::ADPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("ADPDF");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
ADPDF::~ADPDF() {
}


// ___________________________________________________________________________________________ //
bool ADPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;

   CONST(pom1);
   CONST(pom2);
   CONST(reg1);
   
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> ADPDF::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(3,xp,Q);

   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   double xpom  = va_arg(ap, double);
   double zpom  = va_arg(ap, double);
   double muf = va_arg(ap, double);
   va_end(ap);

   return GetQuick({xpom,zpom,muf});
    
}


// ______________________________________________________________________________________ //
std::vector<double> ADPDF::GetQuick(const vector<double>& xpom_zpom_muf) {

   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q


   std::vector<double> pdf(fValue.size(),0);
   if ( xpom_zpom_muf.size() != 3) {
      error["GetQuick"]<<"Quick acces is implemented for two parameter which are 'xpom','zpom' and 'muf'."<<endl;
      return pdf;
   }

   // if calculation failed, return 'fudged' values
   if ( fValue[6]==1 ) { // PDF is nan (see below)
      vector<double> ret{0, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0};
      return ret;
   }    

   cout<<"Daniel ADPDF::GetQuick2 Todo!"<<endl;
   // quick (pdf1, pdf2, reg1)
   
   return pdf;
  
}


// __________________________________________________________________________________________ //
bool ADPDF::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   cout<<"Daniel ADPDF::Update Todo!"<<endl;

   vector<double> pom1 = VALUES(pom1);
   vector<double> pom2 = VALUES(pom1);
   vector<double> reg1 = VALUES(reg1);

   if ( pom2.size() == 1 && pom2[0]==0 ) {
      cout<<"test: no pom2 specified."<<endl;
   }

   fValue = pom1;
   

   return true;
}

