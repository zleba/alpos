// DB 12.02.2015

#include "alpos/ASubsetData.h"
#include <iostream>

using namespace std;

/**
 ASubsetData
*/

// __________________________________________________________________________________________ //
//const std::vector<std::string> ASubsetData::fRequirements = {}; //< List of all AParm's which this function depends on
const std::vector<std::string> ASubsetData::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ASubsetData::fFunctionName = "SubsetData"; //< The function's name


// __________________________________________________________________________________________ //
ASubsetData::ASubsetData(const std::string& name) : AData(name) { 
   SetClassName("ASubsetData");
   // Remember: no access to parameters possible in constructor!
}


// __________________________________________________________________________________________ //
ASubsetData::~ASubsetData() {
}


// ___________________________________________________________________________________________ //
bool ASubsetData::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //! For a SubsetFunction, at least once 'SetRequirementValidPoints' has to called.
   return true;
}


// __________________________________________________________________________________________ //
bool ASubsetData::Update() {
   debug["Update"]<<endl;
   fValue.clear();
   fError.clear();
   if ( GetRequirements().size() != 1) {
      error["Update"]<<"This function 'requires' exactly one other function."<<endl;
      exit(1);
   }
   string req = GetRequirements()[0];
   AData* dat = (AData*)TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+req);
   //const vector<double>& errs = VALUES_ANY(GetAlposName()+std::string(".")+req);
   fValue = AlposTools::VectorSubset(dat->GetValues(),fPointValid);
   fError.resize(fValue.size()); // not yet implemented

   // store reduced data table (all columns)
   fReducedDataTable = dat->GetDataTable();
   for (auto& col : fReducedDataTable) {
      fReducedDataTable[col.first] = AlposTools::VectorSubset(col.second, fPointValid);
   }

   InitErrors(); 
   return true;
}


// __________________________________________________________________________________________ //
void ASubsetData::SetRequirementValidPoints(const std::string& input, const std::vector<bool>& valid) {
   //! Set required function
   //! and set array of ignored/accepted points
   //! Adds "_Data" to req automatically
   
   if ( !ASubsetData::fRequirements.empty() ) {
      error["SetRequirementValidPoints"]<<"Already initialized."<<endl;
      exit(1);
   }
   string req = input+"_Data";

   ASubsetData::fRequirements.resize(1);
   ASubsetData::fRequirements[0] = req;
   // register requirement
   string alias = GetAlposName()+"."+req;
   debug["InitSuperFunctions"]<<" + New ALIAS '"<<alias<<"' for '"<<req<<"'."<<endl;
   TheoryHandler::Handler()->NewAlias(alias,req);
   ARegisterRequirements(this);

   if ( TheoryHandler::Handler()->GetFuncD(req)->N() != valid.size() ) {
      error["SetRequirementValidPoints"]<<"Number of data points in function is not identical to array of accepted/valid points. N="<<TheoryHandler::Handler()->GetFuncD(req)->N()<<", size="<<valid.size()<<endl;
      exit(1);
   }
   fPointValid = valid;

   // --- init errors (not necessarily needed) here
   //Update();
   vector<double> vals = TheoryHandler::Handler()->GetFuncD(req)->GetValues();
   fValue = AlposTools::VectorSubset(vals,fPointValid);
   InitErrors(); // initerrors needs fValue

   SetIsOutdated();
}


// __________________________________________________________________________________________ //
void ASubsetData::InitErrors(){
   //! init errors
   debug["InitErrors"]<<endl;
   AData* dat = (AData*)TheoryHandler::Handler()->GetFuncD(GetAlposName()+std::string(".")+ASubsetData::fRequirements[0]);
   fAllErrors = dat->GetAllErrors() ; // make a dirty copy of all errors
   for ( auto& i : fAllErrors ) i.second.ApplyPointValidMap(fPointValid);

   this->SetHasMultErrors(dat->HasMultErrors());

   fSumErrors.clear(); // clear error cache
   fSumErrMats.clear(); // clear error cache
   // not needed under new interface -> calculation happens on first call to GetSumErrorMatrix(...)
   //CalculateMatrices();
}


// __________________________________________________________________________________________ //
std::string ASubsetData::GetFunctionName() const {
   //! Overloaded function with more specific name:
   //! return the underlying function's name plus "(subset [n/N])"
   //! 
   return GetRealFunctionName();

}
