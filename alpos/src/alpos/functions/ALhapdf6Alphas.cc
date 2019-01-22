#include "alpos/functions/ALhapdf6Alphas.h"

#include <iostream>

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> ALhapdf6Alphas::fRequirements = {"AlphasMz","LHAPDFFile","PDFSet","mur"}; //< List of all AParm's which this function depends on
const std::vector<std::string> ALhapdf6Alphas::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ALhapdf6Alphas::fFunctionName = "LHAPDF6Alphas"; //< The function's name


// ______________________________________________________________________________________ //
ALhapdf6Alphas::ALhapdf6Alphas(const std::string& name) : AParmFuncBase<double>(name) { 
   //cout<<"ALhapdf6Alphas Constructor."<<endl;
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("ALhapdf6Alphas");
   fValue.resize(1);
   fError.resize(1);
}


// ______________________________________________________________________________________ //
ALhapdf6Alphas::~ALhapdf6Alphas() {
   if ( fPDFSet ) delete fPDFSet;
   if ( fPDF ) delete fPDF;
}


// ___________________________________________________________________________________________ //
bool ALhapdf6Alphas::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> ALhapdf6Alphas::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);
   std::vector<double> ret(fValue.size());
   if ( n != 1 ) {
      cout<<"Error in ALhapdf6Alphas::GetQuick(...). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   // for a simple example see: http://en.wikipedia.org/wiki/Variadic_function
   va_list ap;
   va_start(ap, n);
   double mur = va_arg(ap, double);
   va_end (ap);
   ret[0] = fPDF->alphasQ(mur);
   return ret;
}


// ______________________________________________________________________________________ //
std::vector<double> ALhapdf6Alphas::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);
   std::vector<double> ret(fValue.size());
   if ( mur.size() != 1 ) {
      cout<<"Error in ALhapdf6Alphas::GetQuick(vector<double>). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   ret[0] = fPDF->alphasQ(mur[0]);
   return ret;
  
}


// ______________________________________________________________________________________ //
bool ALhapdf6Alphas::Update() {
   if ( fPDFSet )delete fPDFSet;
   if ( fPDF) delete fPDF;
   fPDFSet = new LHAPDF::PDFSet(PAR_S(LHAPDFFile));
   fPDF = fPDFSet->mkPDF(PAR(PDFSet));

   // set "output" parameter AlphasMz
   SET(AlphasMz, fPDF->alphasQ(91.1876), 0)

   fValue[0] = fPDF->alphasQ(PAR(mur));
   fError[0] = 0.;

   return true;
}

// ______________________________________________________________________________________ //
