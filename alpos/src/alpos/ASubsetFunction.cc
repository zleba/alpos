// DB 12.02.2015

#include "alpos/ASubsetFunction.h"
#include <iostream>

using namespace std;
/**
 ASubsetFunction
*/

// __________________________________________________________________________________________ //
//const std::vector<std::string> ASubsetFunction::fRequirements = {}; //< List of all AParm's which this function depends on
const std::vector<std::string> ASubsetFunction::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ASubsetFunction::fFunctionName = "SubsetFunction"; //< The function's name


// __________________________________________________________________________________________ //
ASubsetFunction::ASubsetFunction(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("ASubsetFunction");
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
}


// __________________________________________________________________________________________ //
ASubsetFunction::~ASubsetFunction() {
}


// ___________________________________________________________________________________________ //
bool ASubsetFunction::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   
   //! For a SubsetFunction, at least once 'SetRequirementValidPoints' has to called.
   
   return true;
}


// __________________________________________________________________________________________ //
bool ASubsetFunction::Update() {
   //! update
   if ( fValue.size() == 0 ) return true; // it's all excluded! See: SetRequirementValidPoints()

   if ( GetRequirements().size() != 1) {
      error["Update"]<<"This function 'requires' exactly one other function."<<endl;
      exit(1);
   }

   string req = GetRequirements()[0];
   const vector<double>& vals = VALUES_ANY(GetAlposName()+std::string(".")+req);
   fValue = AlposTools::VectorSubset(vals,fPointValid);
   //const vector<double>& errs = VALUES_ANY(GetAlposName()+std::string(".")+req);
   fError.resize(fValue.size()); // not yet implemented

   if ( TheoryHandler::Handler()->GetFuncD(req)->N() != fPointValid.size() ) {
      error["Update"]<<"Number of data points in function ("<<TheoryHandler::Handler()->GetFuncD(req)->N()<<") is not identical to array of accepted/valid points ("<<fPointValid.size()<<")."<<endl;
      exit(1);
   }

   InitErrors();  // [BUGFIX]
   return true;
}


// __________________________________________________________________________________________ //
void ASubsetFunction::SetRequirementValidPoints(const std::string& req, const std::vector<bool>& valid) {
   //! Set required function
   //! and set array of ignored/accepted points
   if ( !fRequirements.empty() ){
      error["SetRequirementValidPoints"]<<"Already initialized."<<endl;
      exit(1);
   }

   fRequirements.resize(1);
   fRequirements[0] = req;
   // register requirement
   string alias= GetAlposName()+"."+req;
   TheoryHandler::Handler()->NewAlias(alias,req);
   ARegisterRequirements(this);


   fPointValid = valid;
   fValue = AlposTools::VectorSubset(vector<double>(valid.size()),fPointValid); // resize fValue, and keep
   fError = fValue;
   SetIsOutdated();
}

// __________________________________________________________________________________________ //
void ASubsetFunction::InitErrors(){
   //! init errors
   debug["InitErrors"]<<endl;
   AParmFuncBase* th = (AParmFuncBase*)TheoryHandler::Handler()->GetFuncD(GetAlposName() + std::string(".") + ASubsetFunction::fRequirements[0]);
   fAllErrors = th->GetAllErrors() ; // make a dirty copy of all errors
   for ( auto& i : fAllErrors ) i.second.ApplyPointValidMap(fPointValid);

   fSumErrors.clear(); // clear error cache
   fSumErrMats.clear(); // clear error cache
   // not needed under new interface -> calculation happens on first call to GetSumErrorMatrix(...)
   //CalculateMatrices();
}

 // __________________________________________________________________________________________ //
 std::string ASubsetFunction::GetFunctionName() const {
    string name = GetRealFunctionName();
    name +="->";
    const string& fn = TheoryHandler::Handler()->GetFuncD(fRequirements[0])->GetFunctionName();
    name+=fn;
    name += " [";
    
    int n=0;
    for ( const auto& i: fPointValid) if (i) n++;
    name += to_string(n);
    name += "/";
    name += to_string(fPointValid.size());
    name += "]";

    return name;

    return GetRealFunctionName();

    //! return the underlying function's name plus "(subset [n/N])"
//    //! 
//    if ( fRequirements.empty() ) 
//       return "ASubsetFunction not yet initialized. Please call SetRequirementValidPoints()";

//    const string& fn = TheoryHandler::Handler()->GetFuncD(fRequirements[0])->GetFunctionName();
//    string ret = fn + " (Subset [";
//    int n = 0;
//    for ( const auto& i: fPointValid) if (i) n++;
//    ret += to_string(n);
//    ret += "/";
//    ret += to_string(fPointValid.size());
//    ret += "])";
//    return ret;

}
