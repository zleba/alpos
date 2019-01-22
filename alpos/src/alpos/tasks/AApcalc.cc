#include "alpos/tasks/AApcalc.h"
#include "alpos/AFactory.h"
#include <iostream>
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
/* 
 ATask

 */

using namespace std;

extern "C" {

  void apcero_(int* NX, double* X, double* VX, int* NF, double* F); 
  void apc_(int* NX, double* X, double* VX, int* NF, double* F, double* STATUS, int* ISP, int* IST, double* WORK);
  void apcres_(int* NX, double* X, double* VX, double* WORK);
  void apcpul_(int* NX, double* X, double* VX, double* PULL, double* WORK);
  void apcova_(int* NX, double* VX, double* WORK);
}


const string AApcalc::fTaskType = "AApcalc";

//____________________________________________________________________________________ //
//AApcalc::AApcalc(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AApcalc::AApcalc(const string& aname ) : ATask(aname) {
   SetClassName("AApcalc");
   //! constructor

   //! Important: create always new result-object here!
   fResult = new AApcalcResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AApcalc::~AApcalc(){
   //! destructor.
   //! do not delete the AResult object
}


//____________________________________________________________________________________ //
bool AApcalc::Init(){
   return true;
}


//____________________________________________________________________________________ //
bool AApcalc::Execute(){
  
  bool lognormal = false;
  if(STRING_NS(Chisq,NS()) == "LogNormal") {
     info["Execute"] << "AApcalc: Running in LogNormal mode." << endl;
     lognormal = true;
  }
  else if (STRING_NS(Chisq,NS()) == "Covariance") {
     info["Execute"] << "AApcalc: Running in 'Normal' mode." << endl;
  }
  else {
     error["Execute"]<<"Chisq must be 'LogNormal' or 'Covariance' for AApcalc. Exiting."<<endl;
     exit(1);
  }

  // get data and theory
  const auto& super = TheoryHandler::Handler()->GetSuperPair();
  vector<double> X = super.first->GetValues();
  vector<double> th = super.second->GetValues();
  int nData = X.size();
  if(lognormal) {
     info["Execute"] << "AApcalc: Taking log of data." << endl;
     for(int j=0; j<nData; ++j) {
	X[j] = log(X[j]);
     }
  } 
  
  // get errors
  const map<string,AError>& errors = super.first->GetAllErrors();
  vector<vector<double>> systErrors;
  vector<vector<double>> systErrorsUp;
  vector<vector<double>> systErrorsDown;
  for( auto ie : errors) {
    if(ie.second.GetCorrelatedFraction() > 0) {
      X.push_back(0);
      //systErrors.push_back(ie.second.GetErrorCorrRelAvg());  // old interface
      systErrors.push_back(ie.second.GetError("RelAvCor"));

      //systErrorsUp.push_back(ie.second.GetErrorCorrRelUp());  // old interface
      systErrorsUp.push_back(ie.second.GetError("RelUpCor"));

      //systErrorsDown.push_back(ie.second.GetErrorCorrRelDn());  // old interface
      systErrorsDown.push_back(ie.second.GetError("RelDnCor"));
      info["Execute"]<< "AApcalc: Include parameter for correlated part of uncertainty: " << ie.first << endl;
    }
    else {
      info["Execute"] << "AApcalc: Do not include parameter for uncorrelated uncertainty: " << ie.first << endl;
    }
  }
  int nSysts = systErrors.size();
  
  vector<string> parnames = STRING_ARR_NS(FitParameters,NS());
  int nFitPar = parnames.size();
  for(int j=0; j<parnames.size(); ++j) {
    X.push_back(PAR_ANY(parnames[j]));
    info["Execute"] << "AApcalc: Include free parameter " << parnames[j] << " with starting value " << PAR_ANY(parnames[j]) << endl;
  }


  // size of vectors
  int NX = nData+nSysts+nFitPar;
  int Ncov = nData+nSysts; // number of parameters w/ non-zero entries in cov matrix
  int NF = nData;
  int N = NX+NF;
  if(NX != X.size()) {
     info["Execute"] << "AApcalc: ERROR: Size of X does not fit number of included Parameters. Size of X: " << X.size() << " Incl. Param.: " << NX << endl;
     return false;
  }
  
  // arrays fo apcalc
  double VX[(NX*NX+NX)/2];
  double F[NF];
  double STATUS[6];
  double WORK[(N*N+N)/2 + 5*N + 5*NF];
  double PULL[NX];

  int ISP = 7;   // steering parameter
  int IST;  // parameter for fit status

    
  for(int j=0; j<(NX*NX+NX)/2; ++j)  VX[j] = 0;


  //const TMatrixD& Cov = super.first->GetCovariance();  // old interface
  const TMatrixD& Cov = super.first->GetSumErrorMatrix("AA", "AbsAvTot");
  //const TMatrixD& CovRel = super.first->GetCovarianceStatUncorRel();  // old interface
  TMatrixDSym CovRel(super.first->GetSumErrorMatrix("AS", "RelAv"));
  CovRel += super.first->GetSumErrorMatrix("AY", "RelAvUnc");
  int i=0;
  if(lognormal) {
    info["Execute"] << "AApcalc: Using relative covariance." << endl;
    for(int j=0; j<nData; ++j) {
      for(int k=0; k<=j; ++k) {
        VX[i] = CovRel[j][k];
        // cout << "Setting non-zero matrix element " << i << endl;
        ++i;
      }
    }
    // cout << "Matrix elements for data are set." << endl;
    int fla = nData;
    for(int k=i+fla; k<(Ncov*Ncov+Ncov)/2; k += fla+1) {
      VX[k] = 1;
      // cout << "Setting non-zero matrix element " << k << endl;
      ++fla;
    }
  }
  else {
    info["Execute"] << "AApcalc: Using absolute covariance." << endl;
    for(int j=0; j<nData; ++j) {
      for(int k=0; k<=j; ++k) {
        VX[i] = Cov[j][k];
        ++i;
      }
    }
  }

  info["Execute"] << "AApcalc: Starting Loop" << endl;
  do {
    do  {
      for(int j=0; j<parnames.size(); ++j) {
        SET_ANY(parnames[j], X[nData+nSysts+j], 0);
      }
      th = super.second->GetValues();
      for(int j=0; j<NF; ++j) {
        if(lognormal) {
          F[j] = X[j] - log(th[j]);
          for(int k=0; k<nSysts; ++k) {
            F[j] -= X[nData+k]*(((systErrorsUp[k])[j]-(systErrorsDown[k])[j])/2);
          }
        }
        else {
          F[j] = X[j] - th[j];
        }
      }
      apc_(&NX, &X[0], &VX[0], &NF, &F[0], &STATUS[0], &ISP, &IST, &WORK[0]);
    } while(IST < 0);
    //cout << "Reweight" << endl;
  } while(IST > 0);
  
  // set theory parameters to finally fitted values
  for(int j=0; j<parnames.size(); ++j) {
    SET_ANY(parnames[j], X[nData+nSysts+j], 0);
  }
  
  if(BOOL_NS(PrintResults,NS())) {
    apcres_(&NX, &X[0], &VX[0], &WORK[0]);            
  }
  apcpul_(&NX, &X[0], &VX[0], &PULL[0], &WORK[0]);
  apcova_(&NX, &VX[0], &WORK[0]);

  if(BOOL_NS(PrintCovariance,NS())) {
     warn["Execute"] << "PrintCovariance not yet implemented" << endl;
  }
  

  /*
  // Example: traingle with 3 sides and one angle measured
  int NX = 5; 
  int NF = 2;
  
  double X[NX];
  double VX[(NX*NX+NX)/2];
  double F[NF];
  double pfera;
  double sarea;
  int N = NX+NF;
  double STATUS[6];
  double WORK[(N*N+N)/2 + 5*N + 5*NF];
  int ISP = 7;
  int IST = -2;

  //cout << "AApcalc: Initializing Fit" << endl;
  //apcero_(&NX, &X[0], &VX[0], &NF, &F[0]);

  // triangle, 3 sides and one angle
  X[0] = 10.0; // a
  X[1] =  7.0; // b
  X[2] =  9.0; // c
  X[3] =  1.0; // gamma
  X[4] = 30.0; // Area, starting value
  // covariance matrix
  for(int i=0; i < (NX*NX+NX)/2; ++i) {
    VX[i] = 0;
  }
  VX[0] = 0.05*0.05;
  VX[2] = 0.2*0.2;
  VX[5] = 0.2*0.2;
  VX[9] = 0.02*0.02;

  cout << "AApcalc: Starting Loop" << endl;
  while(IST != 0) {
    pfera = 0.5*(X[0]+X[1]+X[2]);
    sarea = sqrt(pfera*(pfera-X[0])*(pfera-X[1])*(pfera-X[2]));
    F[0] = tan(0.5*X[3])-sarea/(pfera*(pfera-X[2]));
    F[1] = sarea-X[4];
    apc_(&NX, &X[0], &VX[0], &NF, &F[0], &STATUS[0], &ISP, &IST, &WORK[0]);
  }
  apcres_(&NX, &X[0], &VX[0], &WORK[0]); // print results
  */
   // --- job done successfully
   return true;
   
}


//____________________________________________________________________________________ //
