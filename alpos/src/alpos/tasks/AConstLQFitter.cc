#include "alpos/tasks/AConstLQFitter.h"
#include "alpos/AFactory.h"
#include <iostream>
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include <cmath>
#include <ctime>

//
/* 
 ATask

 */

using namespace std;

const string AConstLQFitter::fTaskType = "AConstLQFitter";

//____________________________________________________________________________________ //
//AConstLQFitter::AConstLQFitter(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AConstLQFitter::AConstLQFitter(const string& aname ) : ATask(aname) {
   //! constructor

   //! Important: create always new result-object here!
   fResult = new AConstLQFitterResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AConstLQFitter::~AConstLQFitter(){
   //! destructor.
   //! do not delete the AResult object
   if ( fFitter ) delete fFitter;
   if ( fConstraint ) delete fConstraint;
}


//____________________________________________________________________________________ //
bool AConstLQFitter::Init(){
  const auto& super = TheoryHandler::Handler()->GetSuperPair();
  const vector<string>& FitPar = STRING_ARR_NS(FitParameters, NS());  
  const string ConstrType = STRING_NS(ConstraintType, NS()); 
  
  // getting constraint definition
  if(ConstrType == "Normal") fConstraint = new ANormal(FitPar, super.first, super.second);
  else if(ConstrType == "LogNormal") fConstraint = new ALogNormal(FitPar, super.first, super.second);
  else if(ConstrType == "LogNormalNuisance") fConstraint = new ALogNormalNuisance(FitPar, super.first, super.second);
  else if(ConstrType == "NormalNuisance") fConstraint = new ANormalNuisance(FitPar, super.first, super.second);
  else if(ConstrType == "HeraFitter") fConstraint = new AHeraFitterConstr(FitPar, super.first, super.second);
  else {error["Init"]<<"Constraint type not defined!"<<endl; exit(1);}

  // setup minimizer
  const Eigen::VectorXd& X =  fConstraint->InitialX();
  fFitter = new Apccpp(X);                    // create fitter
  fFitter->SetPrintFlag(INT_NS(PrintFlag, NS()));   
    
  //cout << fFitter->evGetX() << endl;
  //cout << fConstraint->Constraint(fFitter->evGetX()) << endl;
  return true;
}


//____________________________________________________________________________________ //
bool AConstLQFitter::Execute(){
  Eigen::VectorXd F = fConstraint->Constraint(fFitter->evGetX());
  Eigen::MatrixXd Cov = fConstraint->Covariance(fFitter->evGetX());
  // loop for Apccpp
  do {
    do {
      F = fConstraint->Constraint(fFitter->evGetX());
    } while(fFitter->bApcpp(F, Cov));  // while iteration is not finished
    // recalculate covariance
    Cov = fConstraint->Covariance(fFitter->evGetX());
  } while(!fFitter->bIsFinished()); // while fit is not finished

  // print result
  if(BOOL_NS(PrintResults, NS())) fFitter->PrintResult();

  // done succesfully
  return true;
}


//____________________________________________________________________________________ //
