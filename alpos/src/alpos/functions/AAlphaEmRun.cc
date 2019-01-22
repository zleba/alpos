// DB 22.02.2015

#include "alpos/functions/AAlphaEmRun.h"
#include "alpos/ATheory.h"
#include <cmath>
#include <iostream>

using namespace std;

extern "C" {
   // void readewpar_();
   void eprc_init_(int* doprint);
   double aemrun_(double* q2);
}

// __________________________________________________________________________________________ //
const std::vector<std::string> AAlphaEmRun::fRequirements = {"mur"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AAlphaEmRun::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AAlphaEmRun::fFunctionName = "AemRun"; //< The function's name
int AAlphaEmRun::fNinstances = 0;

// __________________________________________________________________________________________ //
AAlphaEmRun::AAlphaEmRun(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("AAlphaEmRun");
   if ( ++fNinstances > 1 ) {
      cout<<"Error. AAlphaEmRun::AAlphaEmRun(). Only one AlphaEmRun function is allowed."<<endl;
      exit(1);
   }
   fValue.resize(1);
   fError.resize(1);
}


// __________________________________________________________________________________________ //
AAlphaEmRun::~AAlphaEmRun() {
}


// ___________________________________________________________________________________________ //
bool AAlphaEmRun::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   // 'real' QCDNUM init is called in constructor

   // readewpar_();
   int doprint = 1;
   eprc_init_(&doprint);
   return true;
}


// __________________________________________________________________________________________ //
bool AAlphaEmRun::Update() {

   double q2 = PAR(mur);
   q2*=q2;
   //double arun = aemrun_(&q2);
   fValue[0] = aemrun_(&q2);
   fError[0] = 0;
 
   return true;
}

// ______________________________________________________________________________________ // 
std::vector<double> AAlphaEmRun::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);
   std::vector<double> ret(fValue.size());
   if ( n != 1 ) {
      cout<<"Error in AAlphaEmRun::GetQuick(...). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   // for a simple example see: http://en.wikipedia.org/wiki/Variadic_function
   va_list ap;
   va_start(ap, n);
   double mur = va_arg(ap, double);
   va_end (ap);
   double q2 = mur*mur;
   ret[0] = aemrun_(&q2);
   return ret;
}


// ______________________________________________________________________________________ // 
std::vector<double> AAlphaEmRun::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using: 
   //!   ::GetQuick(vector<double> mur);   
   std::vector<double> ret(fValue.size());
   if ( mur.size() != 1 ) {
      cout<<"Error in AAlphaEmRun::GetQuick(vector<double>). Quick acces is implemented for one parameter which is 'mur'."<<endl;
      return ret;
   }
   double q2 = mur[0]*mur[0];
   ret[0] = aemrun_(&q2);
   return ret;
}
