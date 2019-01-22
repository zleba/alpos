
#include "alpos/functions/AQcdnumPDF.h"

#include <iostream>


using namespace std;

extern "C" void fpdfxq_( int* iset, double* x, double* qmu2, double *pdfs, int *ichk );

const std::vector<std::string> AQcdnumPDF::fRequirements = {"QcdnumInit","xp","muf"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumPDF::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumPDF::fFunctionName = "QcdnumPDF"; //< The function's name


// __________________________________________________________________________________________ //
AQcdnumPDF::AQcdnumPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AQcdnumPDF");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
AQcdnumPDF::~AQcdnumPDF() {
}


// ___________________________________________________________________________________________ //
bool AQcdnumPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AQcdnumPDF::GetQuick(int n, ...) {
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
std::vector<double> AQcdnumPDF::GetQuick(const vector<double>& xp_muf) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> ret(fValue.size());
   if ( xp_muf.size() != 2) {
      error["GetQuick"]<<"Quick acces is implemented for two parameter which are 'xp' and 'muf'."<<endl;
      return ret;
   }

   int iset=1;//5; // iset: 5 is external PDF
   double muf2 = xp_muf[1]*xp_muf[1];
   double xp  = xp_muf[0];
   int ichk;
   fpdfxq_( &iset, &xp, &muf2, &ret[0], &ichk );

   return ret;
  
}


// __________________________________________________________________________________________ //
bool AQcdnumPDF::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());
   UPDATE(QcdnumInit);

   int iset= 1; // iset: 5 is external PDF
   double muf2 = PAR(muf)*PAR(muf);
   double xp = PAR(xp);
   int ichk;
   fpdfxq_( &iset, &xp, &muf2, &fValue[0], &ichk );

   return true;
}

