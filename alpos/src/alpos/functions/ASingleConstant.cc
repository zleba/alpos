// DB 15.01.2015

#include "alpos/functions/ASingleConstant.h"
#include <iostream>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> ASingleConstant::fRequirements = {"theConst"}; //< List of all AParm's which this function depends on
const std::vector<std::string> ASingleConstant::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ASingleConstant::fFunctionName = "SingleConstant"; //< The function's name


// __________________________________________________________________________________________ //
ASingleConstant::ASingleConstant(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
}


// __________________________________________________________________________________________ //
ASingleConstant::~ASingleConstant() {
}


// ___________________________________________________________________________________________ //
bool ASingleConstant::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
 return true;
}


// __________________________________________________________________________________________ //
bool ASingleConstant::Update() {
   fValue = VALUES(theConst);
   fError.resize(fValue.size());

   return true;
}

