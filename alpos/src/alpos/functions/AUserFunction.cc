
#include "alpos/functions/AUserFunction.h"

#include <iostream>

using namespace std;

const std::vector<std::string> AUserFunction::fRequirements = {"par1","par2"}; //< List of all AParm's which this function depends on. These must be specified in the steering
const std::vector<std::string> AUserFunction::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AUserFunction::fFunctionName = "UserFunction"; //< The function's name


// Todo (possibly):
// AUserFunction constructor should made 'private'
// Only 'TheoryHandler' is allowed to 'construct' AParmFunc's
// and then also calls 'ARegisterRequirements()' and 'Init()'.
// Or this is anyhow done for ALL AParmBase's by the theory handler, but AParmConst has Init()={} and GetReq()={};

// ___________________________________________________________________________________________ //
AUserFunction::AUserFunction(const std::string& name) : AParmFuncBase<double>(name) { 
}


// ___________________________________________________________________________________________ //
AUserFunction::~AUserFunction() {
}


// ___________________________________________________________________________________________ //
bool AUserFunction::Init() {
   //! Init is once called for each function
   //! Initialize all needed member variables (typically nothing is needed here)
   //!
   //! return true if initialization was successful.
   //! 
   debug["Init"]<<"GetAlposName:" <<GetAlposName()<<endl;
  return true;
}


// ___________________________________________________________________________________________ //
bool AUserFunction::Update() {
   //!
   //! Provide calculations for return values of this function.
   //! Return values must be stored in member variable fValue[].
   //! 
   debug["Update"]<<"AlposName:" <<GetAlposName()<<endl;
   int nValues = 2;
   fValue.resize(nValues);
   fError.resize(nValues);

   // an example. return alpos parameters as function return-values
   double par1 = PAR(par1);
   double par2 = PAR(par2);

   fValue[0] = par1;
   fValue[1] = par2;

   return true;
}


// ___________________________________________________________________________________________ //
