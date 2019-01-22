
#include "alpos/functions/ALhapdf6.h"

#include <iostream>


using namespace std;

const std::vector<std::string> ALhapdf6::fRequirements = {"LHAPDFFile",  // PDF name. See LHAPDF for possible sets.
							  "PDFSet", // PDFset. Typially 0 is central value
							  "xp", "muf", // x and muf values: Usually integration constants
                                                         }; //< List of all AParm's which this function depends on
const std::vector<std::string> ALhapdf6::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ALhapdf6::fFunctionName = "LHAPDF6"; //< The function's name


// __________________________________________________________________________________________ //
ALhapdf6::ALhapdf6(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("ALhapdf6");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
ALhapdf6::~ALhapdf6() {
   if ( fPDFSet ) delete fPDFSet;
   if ( fPDF ) delete fPDF;
}


// ___________________________________________________________________________________________ //
bool ALhapdf6::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> ALhapdf6::GetQuick(int n, ...) {
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
std::vector<double> ALhapdf6::GetQuick(const vector<double>& xp_muf) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> ret(fValue.size());
   if ( xp_muf.size() != 2) {
      cout<<"Error in ALhapdf6::GetQuick(vector). Quick acces is implemented for two parameter which are 'xp' and 'muf'."<<endl;
      return ret;
   }
   if ( !fPDFSet || !fPDF ) {
      error["GetQuick"]<<"Error in ALhapdf6::GetQuick(vector). PDF is not initialized. "<<endl;
      error["GetQuick"]<<"Supposely, this function is access from 'init' of another task."<<endl;
      error["GetQuick"]<<"Make further sure, that you have updated 'ALhapdf6, before access GetQuick()."<<endl;
      error["GetQuick"]<<"Exiting"<<endl;
      exit(1);
      return ret;
   }

   fPDF->xfxQ(xp_muf[0],xp_muf[1],ret);
   return ret;
  
}


// __________________________________________________________________________________________ //
bool ALhapdf6::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());
   
   //xfxQ2(double x, double q2, std::vector<double>& rtn)
   if ( fPDFSet ) delete fPDFSet;
   if ( fPDF) delete fPDF;
   fPDFSet = new LHAPDF::PDFSet(PAR_S(LHAPDFFile));
   fPDF = fPDFSet->mkPDF(PAR(PDFSet));
   fPDF->xfxQ(PAR(xp),PAR(muf),fValue);
   // fnPDFs = PDFSet->size();

   return true;
}

