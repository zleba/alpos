#include "alpos/tasks/AApcFitter.h"
#include <iostream>
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AConstraint.h"

/* 
 AApcFitter

 add a docu here.

 */

using namespace std;

extern "C" {

  void apcero_(int* NX, double* X, double* VX, int* NF, double* F); 
  void apc_(int* NX, double* X, double* VX, int* NF, double* F, double* STATUS, int* ISP, int* IST, double* WORK);
  void apcres_(int* NX, double* X, double* VX, double* WORK);
  void apcpul_(int* NX, double* X, double* VX, double* PULL, double* WORK);
  void apcova_(int* NX, double* VX, double* WORK);
  void aprxvx_(int* NX, double* X, double* VX);
  void apcorr_(int* NX, double* VX);
}

const string AApcFitter::fTaskType = "ApcFit";

//____________________________________________________________________________________ //
//AApcFitter::AApcFitter(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AApcFitter::AApcFitter(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AApcFitter");
   //! Important: create always new result-object here!
   fResult = new AApcFitterResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AApcFitter::~AApcFitter(){
   //! destructor.
   //! Do not delete the AResult object!
  if ( fConstraint ) delete fConstraint;
}


//____________________________________________________________________________________ //
bool AApcFitter::Init(){

  const auto& super = TheoryHandler::Handler()->GetSuperPair();
  vector<string> FitPar = STRING_ARR_NS(FitParameters,NS());
  const string ConstrType = STRING_NS(ConstraintType, NS()); 

  // get constraint formulation
  if(ConstrType == "Normal") fConstraint = new ANormal(FitPar, super.first, super.second);
  else if(ConstrType == "LogNormal") fConstraint = new ALogNormal(FitPar, super.first, super.second);
  else if(ConstrType == "LogNormalNuisance") fConstraint = new ALogNormalNuisance(FitPar, super.first, super.second);
  else if(ConstrType == "NormalNuisance") fConstraint = new ANormalNuisance(FitPar, super.first, super.second);
  else if(ConstrType == "HeraFitter") fConstraint = new AHeraFitterConstr(FitPar, super.first, super.second);
  else {error["Init"]<<"Constraint type not defined!"<<endl; exit(1);}

  
  return true;
}


//____________________________________________________________________________________ //
bool AApcFitter::Execute(){
  Eigen::VectorXd X = fConstraint->InitialX();
  Eigen::VectorXd F = fConstraint->Constraint(X);
  
  // sizes have to be passed to apcalc
  int NX = X.size();
  int NF = F.size();
  int N = NX+NF;

    
  // arrays fo apcalc + X.data() and F.data()
  double VX[(NX*NX+NX)/2];
  double STATUS[6];
  double WORK[(N*N+N)/2 + 5*N + 5*NF];
  double PULL[NX];
  
  int ISP = INT_NS(Steering,NS());   // steering parameter
  int IST;                           // parameter for fit status
  
  // init VX to 0, not done by apcero as X is already filled
  for(int j=0; j<(NX*NX+NX)/2; ++j) {
    VX[j] = 0;
  }
    
  // obtain covariance matrix
  Eigen::MatrixXd Covar = fConstraint->Covariance(X);
  int i=0;      
  for(int j=0; j<Covar.cols(); ++j) {
    for(int k=0; k<=j; ++k) {
      VX[i] = Covar(j, k);
      ++i;
    }
  }
    


  do {
    do  {
      F = fConstraint->Constraint(X);
      apc_(&NX, &X.data()[0], &VX[0], &NF, &F.data()[0], &STATUS[0], &ISP, &IST, &WORK[0]);
    } while(IST < 0); // while iteration is not finished
    // recalculate covariance matrix
    Covar = fConstraint->Covariance(X);
    int i=0;      
    for(int j=0; j<Covar.cols(); ++j) {
      for(int k=0; k<=j; ++k) {
        VX[i] = Covar(j, k);
        ++i;
      }
    }
  } while(IST > 0);
  
      
  // print result
  if(BOOL_NS(PrintResult,NS())) {
     apcres_(&NX, &X.data()[0], &VX[0], &WORK[0]);            
  }
  // print correlations, not really usefull yet
  if(BOOL_NS(PrintCorrelation,NS())) {
    apcpul_(&NX, &X.data()[0], &VX[0], &PULL[0], &WORK[0]);
    apcova_(&NX, &VX[0], &WORK[0]);
    for(int i=0; i<NX; ++i){
      cout << VX[i] << endl;
    }
    apcorr_(&NX, &VX[0]);
  }
      

  return true;
}


//____________________________________________________________________________________ //
