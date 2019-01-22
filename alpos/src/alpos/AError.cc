#include "alpos/AError.h"
#include "alpos/ATheory.h"
#include "alpos/AlposTools.h"
#include <algorithm>

/** 
 AError

 Base implementation of an error source

 */

using namespace std;

//____________________________________________________________________________________ //
AError::AError() : AlposObject() {
   SetClassName("AError");
}


//____________________________________________________________________________________ //
AError::AError(const string& ErrorName, const string& ErrorSet) : AlposObject(){
   SetClassName("AError");
   SetErrorName(ErrorName,ErrorSet);
   //cout<<"New error. ErrorName="<<fErrorName<<", errorSet="<<fErrorSet<<endl;
}


//____________________________________________________________________________________ //
AError::~AError(){

}


//____________________________________________________________________________________ // 
void AError::SetErrorName(const std::string& ErrorName, const std::string& ErrorSet) { 
   //! Set error name. ErrorName an ErrorSet is required  
   fErrorSet = ErrorSet; 
   fErrorName = (ErrorSet=="") ? ErrorName : ErrorSet+"_"+ErrorName; 
   //fIsStat = ( ErrorName.find("stat") != string::npos || ErrorName.find("Stat") != string::npos ); 
   //fIsStat = (ErrorName=="stat" || ErrorName=="Stat" );//
   //if(fIsStat) fType[1]='S';
} 


//____________________________________________________________________________________ //
double AError::GetAvg(double up, double dn) const {
   //! Calculate average of an up and down error following the averaging prescription
   //! See AError::UpDnAveraging for details
   double sign = GetSign(up,dn);
   // sgn   kLinear  kAbsLinear  kAbsSignLinear  kImprovedLinear  kSignImprovedLinear  kQuadratic  kImprovedQuadratic  kSignImprovedQuadratic kMax kSignMax kMin kSignMin
   switch ( fUpDnAvg ) {
   case kLinear: 
      return (up-dn)/2.;
   case kAbsLinear:
      return (fabs(up)+fabs(dn))/2.; 
   case kAbsSignLinear:
      return (fabs(up)+fabs(dn))/2.*sign;
   case kImprovedLinear:
      return (fmax(fmax(up,dn),0.)-fmin(fmin(up,dn),0.))/2.;
   case kSignImprovedLinear:
      return (fmax(fmax(up,dn),0.)-fmin(fmin(up,dn),0.))/2. *sign;
   case kQuadratic:
      return sqrt((up*up+dn*dn)/2.);
   case kImprovedQuadratic:
      return sqrt((pow(fmax(fmax(up,dn),0.),2)+pow(fmin(fmin(up,dn),0.),2))/2.);
   case kSignImprovedQuadratic:
      return sqrt((pow(fmax(fmax(up,dn),0.),2)+pow(fmin(fmin(up,dn),0.),2))/2.) *sign;      
   case kMax:
      return fmax(fabs(up),fabs(dn));
   case kSignMax:
      return fmax(fabs(up),fabs(dn))*sign;
   case kMin:
      return fmin(fabs(up),fabs(dn));
   case kSignMin:
      return fmin(fabs(up),fabs(dn))*sign;
   default: return 0;
   }
}


//____________________________________________________________________________________ //
std::vector<double> AError::GetAvg(const std::vector<double>& up, const std::vector<double> dn) const {
   //!< Get average of asymmetric uncrtainty following the avereraging prescription
   vector<double> ret;
   CalcAvg(ret,up,dn);
   return ret;
}


//____________________________________________________________________________________ //
void AError::CalcAvg(std::vector<double>& avg, const std::vector<double>& up, const std::vector<double> dn) const {
   //!< Get average of asymmetric uncrtainty following the avereraging prescription
   if ( avg.size() != up.size() ) avg.resize(up.size());
   for ( unsigned int i = 0 ; i<up.size() ; i++ ) {
      avg[i] = GetAvg(up[i],dn[i]);
   }
}


//____________________________________________________________________________________ //
void AError::SetSymError(const std::vector<double>& Err, const std::vector<double>& Sigma, bool RelativeVals, double CorrFrac, std::string nature){
   using namespace AlposTools;
   //! Set values of an symmetric error
   std::vector<double> ErrDn = Err;
   //for ( auto& i : ErrDn ) i*=-1
   ErrDn *=-1.;
   SetAsymError(Err,ErrDn,Sigma,RelativeVals,CorrFrac,nature);
}


//____________________________________________________________________________________ //
void AError::SetAsymError(const std::vector<double>& ErrUp, const std::vector<double>& ErrDn, const std::vector<double>& Sigma, bool RelativeVals, double CorrFrac, std::string nature, UpDnAveraging avg) {
   //! Set values of an asymmetric error
   using namespace AlposTools;
   debug["SetAsymError"]<<endl;

   fErrUp = ErrUp;
   fErrDn = ErrDn;
   if ( !RelativeVals ) {
      fErrUp /= Sigma;
      fErrDn /= Sigma;
   }
   fSigma = Sigma;
   fRelVals = RelativeVals; //! keep how the initial information was passed
   // check for invalid CorrFrac
   if ( (CorrFrac < 0)||(CorrFrac > 1) ) {
      error["SetAsymError"]<<"Cannot construct error '"<<fErrorName<<"': requested correlation fraction lies outside the range [0, 1]."<<endl;
      exit(1);  // makes no sense to continue
   }
   else {
      fCorrFrac = CorrFrac;
   }
   fNature = nature;
   fUpDnAvg = avg;

   fErrors["RelUpTotAll"] = fErrUp;
   fErrors["RelDnTotAll"] = fErrDn;

   // calculate "derived" covariance matrix + fCorr here
   TMatrixDSym mat(fSigma.size());
   //mat.Zero();

   // calculate and cache the 'averaged' error
   const std::vector<double> errRATR = GetError("RelAvTotAll");
   
   // fill uncorrelated part into 'derived' matrix
   for ( unsigned int x = 0 ; x<fSigma.size() ; x++ )  {
      mat[x][x] += pow(errRATR[x],2) * (1-fCorrFrac);
   }
   // fill correlated part into 'derived' matrix
   for ( unsigned int x = 0 ; x<fSigma.size() ; x++ )  {
      for ( unsigned int y = 0 ; y<fSigma.size() ; y++ )  {
         mat[x][y] += errRATR[x] * errRATR[y] * fCorrFrac;
      }
   }


   // --- 'derived' rel covariance matrix
   fErrorMats["RelAvTotAll"].ResizeTo(mat);
   fErrorMats["RelAvTotAll"]=mat;

   // --- calculate correlations (fCorr)
   fCorr.ResizeTo(mat);
   fCorr=mat;
   AlposTools::CovRelToCov(fCorr,fSigma);
   AlposTools::CovToCorr(fCorr);

   // TODO: remnants of old interface -> remove once these become accessible through new interface
   // --- set fMatErr, fMatErrAbs to zero
   fMatErr.clear();
   fMatErrAbs.clear();
   fMatErr    = vector<double>(fErrUp.size());
   fMatErrAbs = vector<double>(fErrUp.size());
}


//____________________________________________________________________________________ //
void AError::SetMatrixError(const std::vector<double>& Err, const std::vector<std::vector<double> >& mat, const std::vector<double>& Sigma, bool RelativeVals, bool IsCorrelationMatrix, std::string nature) {
   //! Set an error with correlation or covariance matrix

   TMatrixDSym Tmat(mat.size());
   for ( unsigned int x = 0 ; x<mat.size() ;x++ ){
      //for ( unsigned int y = 0 ; y<x+1 ;y++ ){
      for ( unsigned int y = 0 ; y<mat.size() ;y++ ){
	      Tmat[x][y] = mat[x][y];
      }
   }
   SetMatrixError(Err, Tmat, Sigma, RelativeVals, IsCorrelationMatrix, nature);
}


//____________________________________________________________________________________ //
void AError::SetMatrixError(const std::vector<double>& Err, const TMatrixDSym& mat, const std::vector<double>& Sigma, bool RelativeVals, bool IsCorrelationMatrix, std::string nature) {
   //! Set an error with correlation or covariance matrix
   //!
   //! A matrix error cannot have asymmetric uncertainties.
   //!  
   //! Err:      vector of errors (RelativeVals: relative or absolute values)
   //! mat:      Covariance or correlation matrix
   //! Sigma:    Cross sections for normalisation
   //! RelativeVals:        if 'Err' are in relative or aboslute units
   //! IsCorrelationMatrix: if false: it is a 'Covariance matrix' (in case of covariance matrix it is assumed that sqrt(cov[i][i]) ~= err*sigma [if err is relative])
   //! nature:   information how to use this error source
   //!

   using namespace AlposTools;
   debug["SetMatrixError"]<<endl;

   fSigma = Sigma;
   fRelVals = RelativeVals; // keep how the initial information was passed
   fCorrFrac = -1; // -1 denotes to use the matrix
   fNature = nature;

   // check if the input covariance matrix 'mat' is a square matrix
   if (mat.GetNrows() != mat.GetNcols()) {
      error["SetMatrixError"] << "Provided matrix is not a square matrix! Exiting..." << endl;
      exit(1);
   }

   fErrUp = Err;
   fErrDn = Err;
   if ( !RelativeVals ) {
      fErrUp /= Sigma;
      fErrDn /= Sigma;
   }

   fErrDn *= -1;  // convention: 'down' errors have negative sign

   // compute absolute errors
   vector<double> errAbs = fErrUp;
   errAbs *= fSigma;

   // check if input 'Err' is the same length as the covariance matrix 'mat'
   if (Err.size() != mat.GetNrows()) {
      error["SetMatrixError"] << "Provided error vector is not the same size as the provided matrix! Exiting..." << endl;
      exit(1);
   }

   // allocate matrices
   fErrorMats["RelAvTotAll"].ResizeTo(mat);
   fCorr.ResizeTo(mat);
   fCorr=mat; 
   
   if ( IsCorrelationMatrix ) {
      // set the covariance matrix
      fErrorMats["AbsAvTotAll"].ResizeTo(mat);
      fErrorMats["AbsAvTotAll"]=mat;
      AlposTools::CorrToCov(fErrorMats["AbsAvTotAll"], errAbs);
      fErrorMats["RelAvTotAll"]=fErrorMats["AbsAvTotAll"];
      AlposTools::CovToCovRel(fErrorMats["RelAvTotAll"], fSigma);
   }
   else { // input is covariance matrix
      // check if input 'Err' is compatible/identical to input covariance matrix 'mat'
      debug["SetMatrixError"]<<"Checking if input 'Err' is compatible/identical to input covariance matrix 'mat'."<<endl;
      for (unsigned int iObs = 0; iObs < errAbs.size(); iObs++) {
         if (!IdenticalToNSignificantDigits(errAbs[iObs], sqrt(mat[iObs][iObs]), 3)) {
            warn["SetMatrixError"] << "Error value in data table is not compatible with the covariance matrix! "
				   << " absErr[" << iObs << "]=" << errAbs[iObs] 
				   << "\tsqrt(mat[" << iObs << "][" << iObs << "])=" << sqrt(mat[iObs][iObs]) << "\tratio: "<<sqrt(mat[iObs][iObs])/errAbs[iObs]<<endl;
         }
      }
      // set the correlations
      fErrorMats["AbsAvTotAll"]=mat;
      fErrorMats["RelAvTotAll"]=mat;
      AlposTools::CovToCovRel(fErrorMats["RelAvTotAll"], fSigma);
      AlposTools::CovToCorr(fCorr);
   }

   // fErrorMats["AbsAvTotAll"].ResizeTo(fCov);
   // fErrorMats["AbsAvTotAll"] = fCov;
   fErrors["RelUpTotAll"] = fErrUp;
   fErrors["RelDnTotAll"] = fErrDn;
   fType[2] = 'C';

   // CovRel
   // fCovRel.ResizeTo(mat);
   // fCovRel=fCov;
   // AlposTools::CovToCovRel(fCovRel,fSigma);


   // not needed under new interface -> calculation happens on first call to GetErrorMatrix(...)
   //ClearDerivedErrorsCache();

   // --- set fMatErr, fMatErrAbs
   fMatErr = fErrUp;
   fMatErrAbs = fErrUp;
   fMatErrAbs *= fSigma;

}


//____________________________________________________________________________________ //
void AError::SetAveragingPrescription(UpDnAveraging Avg) {
   if ( GetIsMatType() ) {
      error["SetAveragingPrescription"]<<"Attempting to set the error averaging prescription for 'matrix'-type error '"<<fErrorName<<"'. Ignoring call."<<endl;
      return;  // don't change fUpDnAvg for 'matrix'-type errors -> is this sensible?
   }
   else {
      fUpDnAvg = Avg;
      ClearDerivedErrorsCache();  // invalidate cached 'Tot' error vectors especially
   }
}


//____________________________________________________________________________________ //
void AError::SetCorrelatedFraction(double corrfrac) {
   if ( GetIsMatType() ) {
      error["SetCorrelatedFraction"]<<"Attempting to set the correlation fraction for 'matrix'-type error '"<<fErrorName<<"'. Ignoring call."<<endl;
      return;  // don't change fCorrFrac for 'matrix'-type errors
   }
   else {
      if ( (corrfrac < 0)||(corrfrac > 1) ) {
         error["SetCorrelatedFraction"]<<"Attempt to set correlation fraction for error '"<<fErrorName<<"' to value outside the range [0, 1]. Ignoring call."<<endl;
      }
      else {
         fCorrFrac = corrfrac;
         ClearDerivedErrorsCache();  // invalidate cached 'Cor' and 'Unc' error vectors especially
      }
   }
}


//____________________________________________________________________________________ //
void AError::ClearDerivedErrorsCache() {
   //! clear cache for all 'derived' error vectors (does not affect cached matrices)
   //! only "RelUpTotAll" and "RelDnTotAll" are 'non-derived'

   std::vector<double> tmpUpErrors = fErrors["RelUpTotAll"];
   std::vector<double> tmpDnErrors = fErrors["RelDnTotAll"];
   fErrors.clear();
   fErrors["RelUpTotAll"] = tmpUpErrors;
   fErrors["RelDnTotAll"] = tmpDnErrors;

   // clear derived matrix cache
   // RelAvTotAll must always be present
   set<string> todel;
   for ( auto i : fErrorMats ) 
      if ( !GetIsMatType() && i.first!="RelAvTotAll" ) todel.insert(i.first);
   for ( auto i : todel ) fErrorMats.erase(i);
   
}


//____________________________________________________________________________________ //
void AError::ApplyPointValidMap(const std::vector<bool>& PointValid) {
   //!< Apply an array to remove unneeded points. 
   //! Should be called only once (if at all!)!
   //! The size of the arrays is then reduced to the number of valid points.

   if ( PointValid.size() != fSigma.size() ) {
      error["ApplyPointValidMap"]<<"Size of PointValid is not compatible with this AError."<<endl;
      error["ApplyPointValidMap"]<<"PointValid.size: "<<PointValid.size()<<",\tErr.size: "<<fSigma.size()<<",\t ErrorName: "<<GetErrorName()<<endl;
      exit(1);
   }
   fSigma = AlposTools::VectorSubset(fSigma,PointValid);
   fErrUp = AlposTools::VectorSubset(fErrUp,PointValid);
   fErrDn = AlposTools::VectorSubset(fErrDn,PointValid);
   if ( GetIsMatType() ) {
      // "reduce" the matrix by removing rows and columns corresponding to "invalid" points
      TMatrixDSym redmat(fSigma.size());
      //TMatrixDSym currentMat(fCorr);//fErrorMats["RelAvTotAll"]);fErrorMats["AbsAvTotAll"]);
      int nx=0;

      for ( int x = 0 ; x<fCorr.GetNcols() ; x++ ) {
	      int ny=0;
	      if ( PointValid[x] ) {
	         for ( int y = 0 ; y<fCorr.GetNrows() ; y++ ) {
	            if ( PointValid[y] ) {
		            redmat[nx][ny] = fCorr[x][y];
		            ny++;
	            }
	         }
	         nx++;
	      }
      }
      //      SetMatrixError(fErrUp, redmat, fSigma, true, false, fNature);
      // delete all, and then re-initialize
      fErrorMats.clear();
      fErrors.clear(); 
      SetMatrixError(fErrUp, redmat, fSigma, true, true, fNature);
   }
   else {
      // "reduce" the error vectors by removing entries corresponding to "invalid" points
      vector<double> ErrUp = AlposTools::VectorSubset(fErrors["RelUpTotAll"],PointValid);
      vector<double> ErrDn = AlposTools::VectorSubset(fErrors["RelDnTotAll"],PointValid);
      // delete everything and re-initialize
      fErrorMats.clear(); 
      fErrors.clear(); 
      SetAsymError(ErrUp, ErrDn, fSigma, true, fCorrFrac, fNature );
      // fErrors["RelUpTotAll"] = AlposTools::VectorSubset(fErrors["RelUpTotAll"],PointValid);
      // fErrors["RelDnTotAll"] = AlposTools::VectorSubset(fErrors["RelDnTotAll"],PointValid);
   }
   //ClearDerivedErrorsCache();

}


//____________________________________________________________________________________ //
const std::vector<double>& AError::GetError(std::string access, bool recalculate) const {
   //! acces to error
   //! access = 'AaaBb[Ccc[Ddd]]' with:
   //!    Aaa:  'Abs', 'Rel'        Obtain relative or absolute uncertinty 
   //!    Bb:   'Up', 'Dn', 'Av'    Obtain up-,down- or averaged uncertainty
   //!    Ccc:  'Tot', 'Cor', 'Unc' Obtain total or only correlated or uncorrelated part
   //!    (deprecated) Ddd:  'Reg', 'Inv'        Obtain 'inverse' or non-inverted (regular) uncertainty     
   //!    Ddd:  'Mul', 'Add', 'All'    Obtain only additive or multiplicative errors (or all [default])
   //!
   //! "RelUpTotAll" and "RelDnTotAll" are present as default
   //!
   //!
   //! \param recalculate If true, ignore the cache and force recalculation of the error.
   //!

   using namespace AlposTools;

   CheckAccess(access);

   // if not already calculated (first request) or if forced, do the calculation
   if ( (!fErrors.count(access)) || recalculate ) {
      // --- 'Up', 'Dn', 'Av'
      // use 'RelUpTotReg' and 'RelDnTotReg' as a starting point -> assume they always exist
      // TODO: maybe check explicitly if they exist
      if ( access.substr(3,2)=="Up" )	 fErrors[access] = fErrors["RelUpTotAll"];
      else if ( access.substr(3,2)=="Dn" )	 fErrors[access] = fErrors["RelDnTotAll"];
      else {
	 CalcAvg(fErrors["RelAvTotAll"],fErrors["RelUpTotAll"],fErrors["RelDnTotAll"]);
	 fErrors[access] = fErrors["RelAvTotAll"]; // keep it
      }

      // --- 'Abs', 'Rel'
      if ( access.substr(0,3)=="Abs") {
	 fErrors[access] *= fSigma;
	 // keep it
         // for absolute errors, multiply by reference value
	 string TmpAcc = "Abs" + access.substr(3,2) + "Tot";
	 fErrors[TmpAcc] = fErrors[access];
      }

      // --- 'Cor', 'Unc', 'Tot'
      if ( access.substr(5,3)=="Cor")       fErrors[access] *= sqrt(fCorrFrac);
      else if ( access.substr(5,3)=="Unc" ) fErrors[access] *= sqrt(1.-fCorrFrac);

      if ( access.substr(8,3)=="Mul" && !fIsMult) {
	 for ( auto& e : fErrors[access] ) e = 0;
      }
      else if ( access.substr(8,3)=="Add" && fIsMult) {
	 for ( auto& e : fErrors[access] ) e = 0;
      }
      
      // // --- keep 'Reg'
      // string RegAcc = access.substr(0,3+2+3) + "All";
      // fErrors[RegAcc] = fErrors[access];

      // --- 'Reg', 'Inv'
      // if ( access.substr(8,3) == "Inv" ) {
      // 	 fErrors[access].clear();
      // 	 fErrors[access].resize(fSigma.size(),1.);
      // 	 fErrors[access] /= fErrors[RegAcc]; // 1./sigma [or better return 1/sigma^2] ??
      // }
   }

   return fErrors[access];
}


//____________________________________________________________________________________ //
const TMatrixDSym& AError::GetErrorMatrix(std::string access, bool recalculate) const {
   //! acces to error
   //!
   //! access = 'AaaBb[Ccc[Ddd]]' with:
   //!           0..3..5...8...
   //!
   //!    Aaa:  'Abs', 'Rel'        Obtain relative or absolute uncertinty 
   //!    Bb:   'Up', 'Dn', 'Av'    Obtain up-,down- or averaged uncertainty
   //!    Ccc:  'Tot', 'Cor', 'Unc' Obtain total or only correlated or uncorrelated part
   //!    (deprecated) Ddd:  'Reg', 'Inv'        Obtain 'inverse' or non-inverted (regular) uncertainty     
   //!    Ddd:  'Mul', 'Add', 'All'    Obtain only additive or multiplicative errors (or all [default])
   //!
   //! "RelAvTotAll" is always present as default for matrix type uncertainties
   //!
   //! \param recalculate If true, ignore the cache and force recalculation of the error.
   //!


   CheckAccess(access);

   // if not already calculated (first request) or if forced, do the calculation
   if ((!fErrorMats.count(access)) || recalculate) {
      debug["GetErrorMatrix"] << "(Re)Calculating error matrix for error '" << fErrorName << "'. actype='" << access << "'" << endl;
      // allocate space
      fErrorMats[access].ResizeTo(fSigma.size(), fSigma.size());
      // check mult/add flag
      if (access.substr(8, 3) == "Mul" && !fIsMult) fErrorMats[access].Zero();
      else if (access.substr(8, 3) == "Add" && fIsMult) fErrorMats[access].Zero();
      else {
         if ( GetIsMatType() ) {
            // if matrix type:
            //   +  Up=Dn=Av,
            //   +  Correlations are specified
            //       =>  if 'Unc'       => Delete off-diagonals
            //       =>  if 'Cor'       => return empty matrix
            const TMatrixDSym& ratr = fErrorMats["RelAvTotAll"];
            if (access.substr(5, 3) == "Unc") {
               fErrorMats[access].Zero();
               for (unsigned int ix = 0; ix < fSigma.size(); ix++) {
                  // only fill the main diagonal
                  fErrorMats[access](ix, ix) = ratr(ix, ix);
               }
            }
            else if (access.substr(5, 3) == "Tot") {
               fErrorMats[access] = ratr; // get the whole matrix
	    }
            else if (access.substr(5, 3) == "Cor") {
               //TODO: warn about having requested 'Cor' for a 'matrix'-type error (?)
	       cout<<"temporary test-> mail to DB."<<endl;
               fErrorMats[access].Zero();
            }
	    if (access.substr(0, 3) == "Abs")
	       AlposTools::CovRelToCov(fErrorMats[access], fSigma);
            // if (access.substr(0, 3) == "Rel")
            //    AlposTools::CovToCovRel(fErrorMats[access], fSigma);
         }
         else {
            // if error is not of matrix type
            // built matrix from vector of errors
            // built and reset:

            // get total error
            const vector<double>& Err = GetError(access.substr(0, 5) + "TotAll");

            // construct matrix
            for (unsigned int x = 0; x < fSigma.size(); x++) {
               // 'uncorrelated' contribution to diagonal
               if ((access.substr(5, 3) != "Cor"))
                  fErrorMats[access](x, x) += pow(Err[x], 2) * (1. - fCorrFrac);

               // 'correlated' contribution
               if ((access.substr(5, 3) != "Unc")) {
                  // 'correlated' contribution to diagonal
                  fErrorMats[access](x, x) += pow(Err[x], 2) * (fCorrFrac);

                  // 'correlated' contribution to off-diagonals
                  for (unsigned int y = 0; y < x; y++) {
                     fErrorMats[access](x, y) = Err[x] * Err[y] * fCorrFrac;
                     fErrorMats[access](y, x) = Err[x] * Err[y] * fCorrFrac;
                  }
               }
            }
         }
      }
      // --- 'Reg', 'Inv'
      // if ( access.substr(8,3) == "Inv" ) {
      // 	 fErrorMats[access].ResizeTo(fErrorMats[RegAcc]);
      // 	 fErrorMats[access] = AlposTools::InvertChol(fErrorMats[RegAcc]);
      // }
   }
   else {
      debug["GetErrorMatrix"] << "Take error matrix from cache for error '" << fErrorName << "'. actype='" << access << "'"<<endl;
   }
   return fErrorMats[access];
}


//____________________________________________________________________________________ //
bool AError::CheckAccess(std::string& access) {
   //! Check 'access' to error
   //! add [Tot] and [All] if not present

   // throw error if 'old way' of accessing inverse matrices is used
   // FIXME: should phase out the "Inv" substring entirely...
   if ( access.size() == 11 && access.substr(8,3)=="Inv"  ) {
      say::error["AError::CheckAccess"]<<"Accessing inverse error matrices using GetError/GetSumError is deprecated. Use AInvMatrices::Instance()->GetInvMatrix() instead."<<endl;
      exit(1);
   }

   if ( access.size() != 5 && access.size() != 8 && access.size() !=11 ){
      say::warn["AError::CheckAccess"]<<"'access' has wrong size, must be 5,8 or 10, but is: "<<access<<endl;
      return false;
   }

   if ( access.size() == 5 )  access+="Tot";
   if ( access.size() == 8 )  access+="All";
   
   if ( access.substr(0,3)!="Abs" && access.substr(0,3)!="Rel" ) {
      //      if (access[0]!='A' && access[0]!='R' ) {// faster
      say::warn["AError::CheckAccess"]<<"First tag [first letter] must be 'Rel' or 'Abs' but is: "<<access.substr(0,3)<<endl;
      return false;
   }

   if ( access.substr(3,2)!="Up" && access.substr(3,2)!="Dn" && access.substr(3,2)!="Av" ) {
      say::warn["AError::CheckAccess"]<<"Second tag [4th letter] must be 'Up', 'Dn' or 'Av' but is: "<<access.substr(3,2)<<endl;
      return false;
   }

   if ( access.substr(5,3)!="Tot" && access.substr(5,3)!="Cor" && access.substr(5,3)!="Unc"  ) {
      say::warn["AError::CheckAccess"]<<"Third tag [6th letter] must be 'Tot', 'Cor' or 'Unc' but is: "<<access.substr(5,3)<<endl;
      return false;
   }

   if ( access.substr(8,3)!="Mul" && access.substr(8,3)!="Add" && access.substr(8,3)!="All"  ) {
      say::warn["AError::CheckAccess"]<<"Fourth tag [9th letter] must be 'Add', 'Mul' or 'All' but is: "<<access.substr(8,3)<<endl;
      return false;
   }

   // if ( access.substr(8,3)!="Reg" && access.substr(8,3)!="Inv"  ) {
   //    //      if (access[8]!='R' ...
   //    say::warn["AError::CheckAccess"]<<"Fourth tag [9th letter] must be 'Reg' or 'Inv' but is: "<<access.substr(8,3)<<endl;
   //    return false;
   // }

   return true;
   
}

//____________________________________________________________________________________ //
bool AError::RescaleError(const std::vector<double>& sigma_rescale) {
   //! Rescale this error to new reference cross sections.
   /*!
       \param sigma_rescale new cross sections to use for calculating
                            the absolute errors/error matrices

       This method overwrites the reference values (fSigma) and
       rescales all absolute errors/error matrices accordingly.
       Relative errors are not changed.

       \note Running this method destroys all cached errors/matrices!
    */

   using namespace AlposTools;

   // check input validity
   if (sigma_rescale.size() != fSigma.size()) {
      error["RescaleError"]<<"The reference cross section vector (size "<<fSigma.size()
			   <<") and the rescaled cross section vector (size "<<sigma_rescale.size()
			   <<") are incompatible for error "<< fErrorName<<endl;
      exit(1);  // remove?
      return false;
   }

   debug["RescaleError"]<<"Rescaling error reference values for error '"<<fErrorName
                        <<"' from {"<<fSigma[0]<<", ...} to {"<<sigma_rescale[0]<<", ...}"<<endl;

   // all uncertainties are stored as 'relative' uncertianty
   // the next access will provide a 'rescaled' uncertainty.
   ClearDerivedErrorsCache();

   /*
   auto errUp = fErrUp;
   auto errDn = fErrDn;

   if ( GetIsMatType() ) {
      TMatrixDSym rescaled_mat = GetErrorMatrix("AbsAvTotAll");//fErrorMats["AbsAvTotAll"];  // salvage reference error matrix 
      // remove all errors
      fErrorMats.clear();
      fErrors.clear();
      // rescale error matrix
      AlposTools::CovRescale(rescaled_mat, fSigma, sigma_rescale);
      if (!fRelVals) {
         // make errors absolute
         errUp *= sigma_rescale;
      }
      // re-initialize with rescaled cross sections as reference
      SetMatrixError(errUp, rescaled_mat, sigma_rescale, fRelVals, false, fNature);
   }
   else {
      // remove all errors
      fErrorMats.clear();
      fErrors.clear(); // TODO: for non-matrix-type errors, only clear cached absolute errors (?)
      if (!fRelVals) {
         // make errors absolute
         errUp *= sigma_rescale;
         errDn *= sigma_rescale;
      }
      // re-initialize with rescaled cross sections as reference
      SetAsymError(errUp, errDn, sigma_rescale, fRelVals, fCorrFrac, fNature);
   }
   */
   return true;

}

