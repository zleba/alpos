#include "alpos/functions/ACRunDecFunction.h"

#include <iostream>

using namespace std;
const std::vector<std::string> ACRunDecFunction::fRequirements = {"AlphasMz","Mz","nLoop","nFlavor","mur"}; //< List of all AParm's which this function depends on
const std::vector<std::string> ACRunDecFunction::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ACRunDecFunction::fFunctionName = "CRunDec"; //< The function's name


// Todo (possibly):
// ACRunDecFunction constructor should made 'private'
// Only 'TheoryHandler' is allowed to 'construct' AParmFunc's
// and then also calls 'ARegisterRequirements()' and 'Init()'.
// Or this is anyhow done for ALL AParmBase's by the theory handler, but AParmConst has Init()={} and GetReq()={};


// ______________________________________________________________________________________ //
ACRunDecFunction::ACRunDecFunction(const std::string& name) : AParmFuncBase<double>(name) { 
   //cout<<"ACRunDecFunction Constructor."<<endl;
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("ACRunDecFunction");
   fValue.resize(1);
   fError.resize(1);
}


// ______________________________________________________________________________________ //
ACRunDecFunction::~ACRunDecFunction() {
   if (!fCRunDec) delete fCRunDec;
   fCRunDec=NULL;
}


// ___________________________________________________________________________________________ //
bool ACRunDecFunction::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   if (!fCRunDec) fCRunDec = new CRunDec();
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> ACRunDecFunction::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);
   std::vector<double> ret(fValue.size());
   if ( n != 1 ) {
      cout<<"Error in ACRunDecFunction::GetQuick(...). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }

   // for a simple example see: http://en.wikipedia.org/wiki/Variadic_function
   va_list ap;
   va_start(ap, n);
   double mur = va_arg(ap, double);
   va_end (ap);
   ret[0] = fCRunDec->AlphasExact(PAR(AlphasMz),PAR(Mz),mur,PAR(nFlavor),PAR(nLoop));
   return ret;
}


// ______________________________________________________________________________________ //
std::vector<double> ACRunDecFunction::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);
   std::vector<double> ret(fValue.size());
   if ( mur.size() != 1 ) {
      cout<<"Error in ACRunDecFunction::GetQuick(vector). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   ret[0] = fCRunDec->AlphasExact(PAR(AlphasMz),PAR(Mz),mur[0],PAR(nFlavor),PAR(nLoop));
   return ret;
}


// ______________________________________________________________________________________ //
bool ACRunDecFunction::Update() {
   fValue[0] = fCRunDec->AlphasExact(PAR(AlphasMz),PAR(Mz),PAR(mur),PAR(nFlavor),PAR(nLoop));
   debug["Update"]<<"CRunDec update Alpha_s value: as(asmz="<<PAR(AlphasMz)<<", mz="<<PAR(Mz)<<", mur="<<PAR(mur)<<") = "<<fValue[0]<<endl;
   fError[0] = 0.;
   return true;
}

// ______________________________________________________________________________________ //
