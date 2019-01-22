
#include "alpos/functions/AExampleFunction.h"

#include <iostream>

using namespace std;

const std::vector<std::string> AExampleFunction::fRequirements = {"bla","blubb","par1","par2","par3","par4","par5"}; //< List of all AParm's which this function depends on
const std::vector<std::string> AExampleFunction::fStopFurtherNotification = {"blubb", "par4"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AExampleFunction::fFunctionName = "ExampleFunction"; //< The function's name


// Todo (possibly):
// AExampleFunction constructor should made 'private'
// Only 'TheoryHandler' is allowed to 'construct' AParmFunc's
// and then also calls 'ARegisterRequirements()' and 'Init()'.
// Or this is anyhow done for ALL AParmBase's by the theory handler, but AParmConst has Init()={} and GetReq()={};

// ___________________________________________________________________________________________ //
AExampleFunction::AExampleFunction(const std::string& name) : AParmFuncBase<double>(name) { 
   cout<<"AExampleFunction Constructor."<<endl;
   //ARegisterRequirements(this); // needed in every constructor
}


// ___________________________________________________________________________________________ //
AExampleFunction::~AExampleFunction() {
}


// ___________________________________________________________________________________________ //
bool AExampleFunction::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ___________________________________________________________________________________________ //
bool AExampleFunction::Update() {
   cout<<" AExampleFunction::Update(). GetAlposName:" <<GetAlposName()<<endl;
   fValue.resize(GetRequirements().size());
   fError.resize(GetRequirements().size());

   fValue[0] = PAR(bla);
   fValue[1] = PAR(par1);
   fValue[2] = PAR(par2);
   //fValue[3] = ((AParmD*)TheoryHandler::Handler()->GetParameter(GetFunctionName()+string(".")+string("par3") ))->GetValue();
   fValue[3] = PAR(par3);
   fValue[4] = PAR(par4);
   fValue[5] = PAR(par5);
   fValue[6] = PAR(blubb);

   return true;
}


// ___________________________________________________________________________________________ //
