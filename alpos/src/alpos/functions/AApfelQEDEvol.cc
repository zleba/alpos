#include "alpos/functions/AApfel.h"
#include "APFEL/APFEL.h"
#include <iostream>


using namespace std;

const std::vector<std::string> AApfelQEDEvol::fRequirements = {"ApfelInit","mur"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AApfelQEDEvol::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelQEDEvol::fFunctionName = "ApfelQEDEvol"; //< The function's name


// __________________________________________________________________________________________ //
AApfelQEDEvol::AApfelQEDEvol(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AApfelQEDEvol");
   fValue.resize(1);
   fError.resize(1);
}


// __________________________________________________________________________________________ //
AApfelQEDEvol::~AApfelQEDEvol() {
}


// ___________________________________________________________________________________________ //
bool AApfelQEDEvol::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.

   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelQEDEvol::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(2,xp,Q);
   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   vector<double> mur = {va_arg(ap, double)};
   va_end(ap);

   return GetQuick(mur);
    
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelQEDEvol::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   if ( mur.size() != 1) {
      error["GetQuick"]<<"Error in AApfelQEDEvol::GetQuick(vector). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      exit(1);
      return fValue;
   }

   std::vector<double> ret = {APFEL::AlphaQED(mur[0])};
   return ret;
  
}


// __________________________________________________________________________________________ //
bool AApfelQEDEvol::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());
   PAR(ApfelInit);
   fValue[0] = APFEL::AlphaQED(PAR(mur));
   return true;
}

