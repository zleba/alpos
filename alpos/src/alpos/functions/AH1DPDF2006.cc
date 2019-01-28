
#include "alpos/functions/AH1DPDF2006.h"
#include "apfel/dglapbuilder.h"

#include <iostream>


//! 
//! A wrapper to h1pdf2006 FitA, and FitB
//! different
//! 

using namespace std;

extern "C" {
   // tint is the 'maximum t'
   void diffpdferr_(double* xpom, double*  zpom, double*  muf, int* ifit, int* ierr,  int* ireg, double *pdfs);
}


const std::vector<std::string> AH1DPDF2006::fRequirements = {  
   "xpom", // xpom (dummy)
   "zpom", // zpom (dummy)
   "muf", // mu_f (dummy)
   "iFit", //FitA=1 or FitB=2
   "ierr", //error set
   "ireg", // 0: both, 1: only IP, 2: only IR
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AH1DPDF2006::fStopFurtherNotification = {"xpom","zpom","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AH1DPDF2006::fFunctionName = "H1DPDF2006"; //< The function's name


// __________________________________________________________________________________________ //
AH1DPDF2006::AH1DPDF2006(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AH1DPDF2006");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
AH1DPDF2006::~AH1DPDF2006() {
}


// ___________________________________________________________________________________________ //
bool AH1DPDF2006::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;
   fValue.resize(13);
   fError.resize(13);
   
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AH1DPDF2006::GetQuick(int n, ...) {
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
std::vector<double> AH1DPDF2006::GetQuick(const vector<double>& xpom_zpom_muf) {

   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> pdf(fValue.size(),0);
   double xpom = xpom_zpom_muf[0];
   double zpom = xpom_zpom_muf[1];
   double muf  = xpom_zpom_muf[2];
   diffpdferr_(&xpom, &zpom, &muf, &fifit, &fierr, &fireg,  &fValue[0]);

   return pdf;
  
}


// __________________________________________________________________________________________ //
bool AH1DPDF2006::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   double xpom = PAR(xpom);
   double zpom = PAR(zpom);
   double muf  = PAR(muf);
   fifit = int(PAR(iFit));
   fierr = int(PAR(ierr));
   fireg = int(PAR(ireg));
   diffpdferr_(&xpom, &zpom, &muf, &fifit, &fierr, &fireg,  &fValue[0]);
   
   return true;
}
