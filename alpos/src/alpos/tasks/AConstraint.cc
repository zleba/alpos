#include "alpos/tasks/AConstraint.h"
#include "alpos/tasks/APDFQ0_HERAStyleNoSum.h"



/* Treating uncertainties as noraml distributed
   uses full absolute covariance matrix
   initial X contains data vector da
   constaint is F_i = da_i - th_i

 */
ANormal::ANormal(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  const int Nda = Data->GetValues().size();
  //const TMatrixD& Cov = Data->GetCovariance();  // old-interface
  const TMatrixD& Cov = Data->GetSumErrorMatrix("AA", "AbsAvTot");
  fCovar = Eigen::MatrixXd::Map(&Cov.GetMatrixArray()[0], Nda, Nda);
}

Eigen::VectorXd ANormal::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  Eigen::VectorXd X(Da.size() + fFitPar.size());
  X.head(Da.size()) =  Eigen::VectorXd::Map(&Da[0], Da.size());
  for( unsigned int i = 0; i < fFitPar.size(); ++i) {
    X(Da.size()+i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd ANormal::Constraint(const Eigen::VectorXd& X) const {
  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X[X.size()-fFitPar.size()+i],0);
  }
  
  const std::vector<double>& Th = fTheo->GetValues();
  Eigen::VectorXd F(X.size()-fFitPar.size());
  for( unsigned int i = 0; i < F.size(); ++i) {
    F(i) = X(i) - Th[i];
  }
  return F;
}

Eigen::MatrixXd ANormal::Covariance(const Eigen::VectorXd& X) const {
  return fCovar;
}


// ___________________________________________________________________________//


/* Treating uncertainties as log normal distributed
   uses full relative covariance matrix
   initial X contains log of data vector da' = log(da)
   constaint is F_i = da_i' - log(th_i)

 */
ALogNormal::ALogNormal(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  const int Nda = Data->GetValues().size();
  //const TMatrixD& Cov = Data->GetCovarianceRel();  // old-interface
  const TMatrixD& Cov = Data->GetSumErrorMatrix("AA", "RelAvTot");
  fCovar = Eigen::MatrixXd::Map(&Cov.GetMatrixArray()[0], Nda, Nda);
}

Eigen::VectorXd ALogNormal::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  Eigen::VectorXd X(Da.size() + fFitPar.size());
  X.head(Da.size()) =  Eigen::ArrayXd::Map(&Da[0], Da.size()).log();
  for( unsigned int i = 0; i < fFitPar.size(); ++i) {
    X(Da.size()+i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd ALogNormal::Constraint(const Eigen::VectorXd& X) const {
  const int Nda = X.size()-fFitPar.size();

  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X(Nda+i),0);
  }
  
  const Eigen::VectorXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda).log();
  return X.head(Nda) - Th;
}

Eigen::MatrixXd ALogNormal::Covariance(const Eigen::VectorXd& X) const {
  return fCovar;
}

//____________________________________________________________________________//

/* Treating uncertainties as log normal distributed with nuisance parameters
   correlated systematic uncertainties
   uses relative covariance matrix of statistic and uncorrelated uncertainties
   initial X contains log of data vector da' = log(da)
   theory is modified to th' = (th + sum b_j MA_j)*exp(sum b_j M_j) + sum b_j A_j
   where MA, M  and A is the nature according to the datafile
   constaint is F_i = da_i' - log(th_i')

 */
ALogNormalNuisance::ALogNormalNuisance(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  
  const int Nda = Data->GetValues().size();

  const map<string, AError>& Errors = Data->GetAllErrors();
  int Nnuis = 0;
  for(auto ie: Errors)
    {
      if(ie.second.GetCorrelatedFraction() > 0)
        ++Nnuis;
    }
  fAddErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  fMultAddErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  fMultErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  int Count = 0;
  for(auto ie : Errors) {
    if(ie.second.GetCorrelatedFraction() > 0) {
      if(ie.second.GetNature() == "A") {
        fAddErr.col(Count) = 
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      }
      else if(ie.second.GetNature() == "MA") {
        fMultAddErr.col(Count) = 
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      }
      else {
        fMultErr.col(Count) =
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      } 
    }
  }

  //const TMatrixD& Cov = Data->GetUncertaintyUncorrRel();  // old-interface
  const TMatrixDSym& Cov = Data->GetSumErrorMatrix("AY", "RelAvUnc");
  fCovar = Eigen::MatrixXd::Identity(Nda + Nnuis, Nda + Nnuis);
  fCovar.topLeftCorner(Nda, Nda) = 
    Eigen::MatrixXd::Map(&Cov.GetMatrixArray()[0], Nda, Nda); 
}

Eigen::VectorXd ALogNormalNuisance::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  const int Nda = Da.size();
  const int Nnuis = fMultErr.cols();
  const int Nfit = fFitPar.size();
  Eigen::VectorXd X = 
    Eigen::VectorXd::Zero(Nda + Nfit + Nnuis); 
  X.head(Nda) =  Eigen::ArrayXd::Map(&Da[0], Nda).log();
  for( unsigned int i = 0; i < Nfit; ++i) {
    X(Nda + Nnuis + i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd ALogNormalNuisance::Constraint(const Eigen::VectorXd& X) const {
  const int Nnuis = fMultErr.cols();
  const int Nda = X.size()-fFitPar.size() - Nnuis;


  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X(Nda+Nnuis+i),0);
  }
  
  Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  Th += (fMultAddErr*X.segment(Nda, Nnuis)).array();
  Th *= (fMultErr*X.segment(Nda, Nnuis)).array().exp();
  Th += (fAddErr*X.segment(Nda, Nnuis)).array();    
  return X.head(Nda) - Th.log().matrix();
}

Eigen::MatrixXd ALogNormalNuisance::Covariance(const Eigen::VectorXd& X) const {
  return fCovar;
}

//____________________________________________________________________________//


/* Treating uncertainties as normal distributed with nuisance parameters
   correlated systematic uncertainties
   uses relative covariance matrix of statistic and uncorrelated uncertainties
   initial X contains data vector da
   theory is modified to th' = (th + sum b_j MA_j)*exp(sum b_j M_j) + sum b_j A_j
   where MA, M  and A is the nature according to the datafile
   constaint is F_i = da_i - th_i'
 */
ANormalNuisance::ANormalNuisance(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  
  const int Nda = Data->GetValues().size();

  const map<string, AError>& Errors = Data->GetAllErrors();
  int Nnuis = 0;
  for(auto ie: Errors)
    {
      if(ie.second.GetCorrelatedFraction() > 0)
        ++Nnuis;
    }
  fAddErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  fMultAddErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  fMultErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  int Count = 0;
  for(auto ie : Errors) {
    if(ie.second.GetCorrelatedFraction() > 0) {
      if(ie.second.GetNature() == "A") {
        fAddErr.col(Count) = 
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      }
      else if(ie.second.GetNature() == "MA") {
        fMultAddErr.col(Count) = 
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      }
      else {
        fMultErr.col(Count) =
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
      } 
    }
  }

  //const TMatrixD& Cov = Data->GetCovarianceStatUncorr();  // old-interface
  TMatrixDSym Cov(Data->GetSumErrorMatrix("AS", "AbsAv"));
  Cov += Data->GetSumErrorMatrix("AY", "AbsAvUnc");
  fCovar = Eigen::MatrixXd::Identity(Nda + Nnuis, Nda + Nnuis);
  fCovar.topLeftCorner(Nda, Nda) = 
    Eigen::MatrixXd::Map(&Cov.GetMatrixArray()[0], Nda, Nda); 
}

Eigen::VectorXd ANormalNuisance::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  const int Nda = Da.size();
  const int Nnuis = fMultErr.cols();
  const int Nfit = fFitPar.size();
  Eigen::VectorXd X = 
    Eigen::VectorXd::Zero(Nda + Nfit + Nnuis); 
  X.head(Nda) =  Eigen::ArrayXd::Map(&Da[0], Nda);
  for( unsigned int i = 0; i < Nfit; ++i) {
    X(Nda + Nnuis + i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd ANormalNuisance::Constraint(const Eigen::VectorXd& X) const {
  const int Nnuis = fMultErr.cols();
  const int Nda = X.size()-fFitPar.size() - Nnuis;


  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X(Nda+Nnuis+i),0);
  }
  
  Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  Th += (fMultAddErr*X.segment(Nda, Nnuis)).array();
  Th *= (fMultErr*X.segment(Nda, Nnuis)).array().exp();
  Th += (fAddErr*X.segment(Nda, Nnuis)).array();    
  return X.head(Nda) - Th.matrix();
}

Eigen::MatrixXd ANormalNuisance::Covariance(const Eigen::VectorXd& X) const {
  return fCovar;
}

//____________________________________________________________________________//


/* Treating uncertainties as normal distributed with nuisance parameters
   for correlated systematic uncertainties according to HERAFitter
   covarince matrix is diagonal with entries d_stat da*th + d_uncorr th*th
   where d_stat and d_uncorr are relaitve and uncorrelated systematics
   initial X contains data vector da
   theory is modified to th' = th(1 + sum b_j G_j) where G_j are the systematic
   uncertainties
   constaint is F_i = da_i - th_i'
 */
AHeraFitterConstr::AHeraFitterConstr(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  
  const int Nda = Data->GetValues().size();

  const map<string, AError>& Errors = Data->GetAllErrors();
  int Nnuis = 0;
  for(auto ie: Errors)
    {
      if(ie.second.GetCorrelatedFraction() > 0)
        ++Nnuis;
    }
  fSystErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  int Count = 0;
  for(auto ie : Errors) {
    if(ie.second.GetCorrelatedFraction() > 0) {
        fSystErr.col(Count) =
          //Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);  // old-interface
          Eigen::VectorXd::Map(&ie.second.GetError("RelAvCor")[0], Nda);
        ++Count;
    }
  }
}

Eigen::VectorXd AHeraFitterConstr::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  const int Nda = Da.size();
  const int Nnuis = fSystErr.cols();
  const int Nfit = fFitPar.size();
  Eigen::VectorXd X = 
    Eigen::VectorXd::Zero(Nda + Nfit + Nnuis); 
  X.head(Nda) =  Eigen::ArrayXd::Map(&Da[0], Nda);
  for( unsigned int i = 0; i < Nfit; ++i) {
    X(Nda + Nnuis + i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd  AHeraFitterConstr::Constraint(const Eigen::VectorXd& X) const {
  const int Nnuis = fSystErr.cols();
  const int Nda = X.size()-fFitPar.size() - Nnuis;


  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X(Nda+Nnuis+i), 0);
  }
  
  Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  Th *= 1 - (fSystErr*X.segment(Nda, Nnuis)).array();
  return X.head(Nda) - Th.matrix();
}

Eigen::MatrixXd AHeraFitterConstr::Covariance(const Eigen::VectorXd& X) const {
  const int Nda = fData->GetValues().size();
  const int Nnuis = fSystErr.cols();
  const Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  const Eigen::ArrayXd Da = Eigen::ArrayXd::Map(&fData->GetValues()[0], Nda);
  const Eigen::ArrayXd Uncorr = 
    //Eigen::ArrayXd::Map(&fData->GetUncertaintyUncorrRel()[0], Nda);  // old-interface
    Eigen::ArrayXd::Map(&fData->GetSumError("AY", "RelAvUnc")[0], Nda);
  //const Eigen::ArrayXd Stat = Eigen::ArrayXd::Map(&fData->GetUncertaintyStatRel()[0], Nda);   // old-interface
  const Eigen::ArrayXd Stat = Eigen::ArrayXd::Map(&fData->GetSumError("AS", "RelAvTot")[0], Nda);
  const Eigen::ArrayXd DiagCovar = Uncorr*Uncorr*Th*Th + Stat*Stat*Th*Da;
  Eigen::MatrixXd Covar = Eigen::MatrixXd::Identity(Nda+Nnuis, Nda+Nnuis);
  Covar.topLeftCorner(Nda, Nda) = DiagCovar.matrix().asDiagonal();
  return Covar;
}

//____________________________________________________________________________//

/* Implementing PDF fit with sum rules, not working yet

 */

/*
 ASumRulesHeraFitter::ASumRulesHeraFitter(const std::vector<string>& FitPar, AData* Data, AFuncD* Theo) 
  : AConstraint(FitPar, Data, Theo) {
  
  const int Nda = Data->GetValues().size();

  const map<string, AError>& Errors = Data->GetAllErrors();
  int Nnuis = 0;
  for(auto ie: Errors)
    {
      if(ie.second.GetCorrelatedFraction() > 0)
        ++Nnuis;
    }
  fSystErr = Eigen::MatrixXd::Zero(Nda, Nnuis);
  int Count = 0;
  for(auto ie : Errors) {
    if(ie.second.GetCorrelatedFraction() > 0) {
        fSystErr.col(Count) =
          Eigen::VectorXd::Map(&ie.second.GetErrorCorrRelAvg()[0], Nda);
        ++Count;
    }
  }
  
  for( unsigned int i = 0; i < fFitPar.size(); ++i) {
    if(fFitPar[i] == "PDFQ0_HERA_NoSum.gA") {
      fgAind = i;
    }
    else if(fFitPar[i] == "PDFQ0_HERA_NoSum.uvA") {
      fuvAind = i;
    }
    else if(fFitPar[i] == "PDFQ0_HERA_NoSum.dvA") {
      fdvAind = i;
    }
    else if(fFitPar[i] == "PDFQ0_HERA_NoSum.UbarA") {
      fUbarAind = i;
    }
  }
}

Eigen::VectorXd  ASumRulesHeraFitter::InitialX() const {
  const std::vector<double>& Da = fData->GetValues();
  const int Nda = Da.size();
  const int Nnuis = fSystErr.cols();
  const int Nfit = fFitPar.size();
  Eigen::VectorXd X = 
    Eigen::VectorXd::Zero(Nda + Nfit + Nnuis); 
  X.head(Nda) =  Eigen::ArrayXd::Map(&Da[0], Nda);
  for( unsigned int i = 0; i < Nfit; ++i) {
    X(Nda + Nnuis + i) = PAR_ANY(fFitPar[i]);
  }
  return X;
}

Eigen::VectorXd   ASumRulesHeraFitter::Constraint(const Eigen::VectorXd& X) const {
  const int Nnuis = fSystErr.cols();
  const int Nda = X.size()-fFitPar.size() - Nnuis;


  for ( unsigned int i = 0; i < fFitPar.size() ; ++i ) {
    SET_ANY(fFitPar[i],X(Nda+Nnuis+i), 0);
  }
  
  Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  Th *= 1 - (fSystErr*X.segment(Nda, Nnuis)).array();
  Eigen::VectorXd F(Nda+4); 
  F.head(Nda) =  X.head(Nda) - Th.matrix();
  
  
  F(Nda) = X(Nda+fgAind) -
    F(Nda+1) = X(Nda+fdvAind) -
    F(Nda+2) = X(Nda+fuvAind) -
    F(Nda+3) = X(Nda+fUbarAind) -
  
  return F;
}

Eigen::MatrixXd  ASumRulesHeraFitter::Covariance(const Eigen::VectorXd& X) const {
  const int Nda = fData->GetValues().size();
  const int Nnuis = fSystErr.cols();
  const Eigen::ArrayXd Th = Eigen::ArrayXd::Map(&fTheo->GetValues()[0], Nda);
  const Eigen::ArrayXd Da = Eigen::ArrayXd::Map(&fData->GetValues()[0], Nda);
  const Eigen::ArrayXd Uncorr = 
    Eigen::ArrayXd::Map(&fData->GetUncertaintyUncorrRel()[0], Nda);
  const Eigen::ArrayXd Stat = Eigen::ArrayXd::Map(&fData->GetUncertaintyStatRel()[0], Nda); 
  const Eigen::ArrayXd DiagCovar = Uncorr*Uncorr*Th*Th + Stat*Stat*Th*Da;
  Eigen::MatrixXd Covar = Eigen::MatrixXd::Identity(Nda+Nnuis, Nda+Nnuis);
  Covar.topLeftCorner(Nda, Nda) = DiagCovar.matrix().asDiagonal();
  return Covar;
}
*/

//____________________________________________________________________________//
