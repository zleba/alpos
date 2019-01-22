// DB 02.2016

#include "alpos/functions/AFunctionTF1.h"
#include <iostream>
#include "TF1.h"
#include <cfloat>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AFunctionTF1::fRequirements = { "Formula",
//							       "InitialValue",
							       "x",// parameters
							       "PAR0", // parameters
							       "PAR1", // parameters
							       "PAR2", // parameters
							       "PAR3", // parameters
							       "PAR4",
							       "PAR5" }; //< List of all AParm's which this function depends on
const std::vector<std::string> AFunctionTF1::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AFunctionTF1::fFunctionName = "FunctionTF1"; //< The function's name


// __________________________________________________________________________________________ //
AFunctionTF1::AFunctionTF1(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
}


// __________________________________________________________________________________________ //
AFunctionTF1::~AFunctionTF1() {
   if ( fTF1 ) delete fTF1;
}


// ___________________________________________________________________________________________ //
bool AFunctionTF1::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.


   if ( fTF1 ) delete fTF1;
   info["Init"]<<"New TF1-function (if code crashes, add brackets around your formula)."<<endl;
   string formula = PAR_S(Formula);
   //fTF1 = new TF1(GetAlposName(),formula.c_str(),DBL_MIN,DBL_MAX);
   //cout<<"New TF1 with formula '"<<formula.c_str()<<"' and name: "<<GetAlposName()<<endl;
   fTF1 = new TF1(GetAlposName().c_str(),formula.c_str());
   
   CONST(Formula);
   
   fValue.resize(1);
   fError.resize(1);
   return true;
}


// __________________________________________________________________________________________ //
bool AFunctionTF1::Update() {
   //
   if ( !fTF1 ) {
      // A 'FunctionTF1' may be called accidentally before initalisation is done!
      // -> Alpos rule! Never access function values before final initialization.
      Init();
      fValue.resize(1);
      fValue[0] = 0;//PAR(InitialValue);
      return false;
   }
   fTF1->SetParameter(0,PAR(PAR0));
   fTF1->SetParameter(1,PAR(PAR1));
   fTF1->SetParameter(2,PAR(PAR2));
   fTF1->SetParameter(3,PAR(PAR3));
   fTF1->SetParameter(4,PAR(PAR4));
   fTF1->SetParameter(5,PAR(PAR5));
   fValue[0] = fTF1->Eval(PAR(x));

   return true;
}

