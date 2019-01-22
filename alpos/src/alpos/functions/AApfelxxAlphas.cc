// DB 07.2017

#include "alpos/functions/AApfelxxAlphas.h"
#include <iostream>

using namespace std;

// __________________________________________________________________________________________ //
const std::vector<std::string> AApfelxxAlphas::fRequirements = {
   "mur", // scale (dummy)
   "AlphaQCDRef", // coupling strength at reference scale
   "MuAlphaQCDRef", // reference scale
   "iOrd", // perturbative order
   "mc","mb","mt" // quark masses 
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AApfelxxAlphas::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelxxAlphas::fFunctionName = "ApfelxxAlphas"; //< The function's name


// __________________________________________________________________________________________ //
AApfelxxAlphas::AApfelxxAlphas(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   fValue.resize(1);
   fError.resize(1);
}


// __________________________________________________________________________________________ //
AApfelxxAlphas::~AApfelxxAlphas() {
}


// ___________________________________________________________________________________________ //
bool AApfelxxAlphas::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// __________________________________________________________________________________________ //
bool AApfelxxAlphas::Update() {

   debug["Update"]<<endl;

   const double AlphaQCDRef = PAR(AlphaQCDRef);
   const double MuAlphaQCDRef = PAR(MuAlphaQCDRef);
   const int PerturbativeOrder = PAR(iOrd);
   const vector<double> Masses = {0, 0, 0, PAR(mc), PAR(mb), PAR(mt)}; 

   fAs = std::unique_ptr<apfel::AlphaQCD>(new apfel::AlphaQCD{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder});
   //AlphaQCD as{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};

   fValue[0] = fAs->Evaluate(PAR(mur));
   return true;
}

// ______________________________________________________________________________________ //
std::vector<double> AApfelxxAlphas::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);

   if ( mur.size() != 1 ) {
      cout<<"Error in AApfelxxAlphas::GetQuick(vector<double>). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return std::vector<double>{0};
   }
   
   return std::vector<double>{fAs->Evaluate(mur[0])};
}

// ______________________________________________________________________________________ //
std::vector<double> AApfelxxAlphas::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);
   if ( n != 1 ) {
      cout<<"Error in ACRunDecFunction::GetQuick(...). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return std::vector<double>(1);
   }

   // for a simple example see: http://en.wikipedia.org/wiki/Variadic_function
   va_list ap;
   va_start(ap, n);
   // double mur = va_arg(ap, double);
   // va_end (ap);
   // ret[0] = fCRunDec->AlphasExact(PAR(AlphasMz),PAR(Mz),mur,PAR(nFlavor),PAR(nLoop));
   return std::vector<double>{fAs->Evaluate(va_arg(ap, double))};
}
