// DB. 01/2015
#ifndef Alpos_AError
#define Alpos_AError

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

#include "alpos/AlposTools.h"
#include "alpos/AlposObject.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixTSym.h"


//typedef TMatrixTSym<double> TMatrixTSymD;

/**
 * AError.
 * A class to hold the uncertainties of one source
 *
 */

class AError : public AlposObject {
public:
   //! Specify how asymetric errors are symmetrized
   enum UpDnAveraging {  
      // Example           sgn   kLinear  kAbsLinear  kAbsSignLinear  kImprovedLinear  kSignImprovedLinear  kQuadratic  kImprovedQuadratic  kSignImprovedQuadratic kMax kSignMax kMin kSignMin
      // up:  5, dn -3  ->  +       4         4             4               4                  4              4.123            4.123                 4.123          5      5      3      3
      // up: -5, dn  3  ->  -      -4         4            -4               4                 -4              4.123            4.123                -4.123          5     -5      3     -3
      // --- samed sign cases:															          
      // up:  5, dn  3  ->  +       1         4             4               2.5                2.5            4.123            3.536                 3.536          5      5      3      3
      // up: -5, dn -3  ->  -      -1         4            -4               2.5               -2.5            4.123            3.536                -3.536          5     -5      3     -3
      // up   3, dn  5  ->  -      -1         4            -4               2.5               -2.5            4.123            3.536                -3.536          5     -5      3     -3
      // up  -3, dn -5  ->  +       1         4             4               2.5                2.5            4.123            3.536                 3.536          5      5      3      3
      // sgn = if (|up|>|dn|): 'sgn(up)' if (|dn|>|up) '-sgn(dn)'
      kLinear,  //!< linear, i.e. (up-dn)/2. (keep sign)
      kAbsLinear, //!< linear using abs values, i.e. (|up|+|dn|)/2. (always positive)
      kAbsSignLinear, //!< linear using abs values, i.e. (|up|+|dn|)/2.*sgn (keep sign)
      kImprovedLinear, //!< improved linear sum, i.e. (max(up,dn,0)-min(up,dn,0)/2. (always positive)
      kSignImprovedLinear, //!< sign improved linear sum, i.e. (max(up,dn,0)-min(up,dn,0)/2.*sgn (keep sign)
      kQuadratic, //!< quadratic average, i.e. sqrt((up*up+dn*dn)/2.) (always positive)
      kImprovedQuadratic, //!< improved quadratic average, i.e. sqrt((max(up,dn,0)**2 + min(up,dn,0)**2)/2) (always positive)
      kSignImprovedQuadratic, //!< improved quadratic average, i.e. sqrt((max(up,dn,0)**2 + min(up,dn,0)**2)/2)*sgn (keep sign)
      kMax, //!< take the larger value, i.e. max(|up|,|dn|) (always positive)
      kSignMax, //!< take the larger value and keep a sign, i.e. max(|up|,|dn|)*sign (keep sign)
      kMin,//!< take the larger value, i.e. max(|up|,|dn|) always positive (always positive)
      kSignMin //!< take the smaller value and keep a sign, i.e. max(|up|,|dn|)*sign (keep sign)
   };
   
public:
   AError();
   AError(const std::string& ErrorName, const std::string& ErrorSet="");
   ~AError();

   void SetErrorName(const std::string& ErrorName, const std::string& ErrorSet); //!< Set error name. ErrorName an ErrorSet is required

   //! initialize members
   //! An AError class is ONLY initialized by one of these functions!
   void SetSymError(const std::vector<double>& Err, const std::vector<double>& fSigma, bool RelVals=true, double CorrFrac=1., std::string nature = ""); //!< Set values of an symmetric error
   void SetAsymError(const std::vector<double>& ErrUp, const std::vector<double>& ErrDn, const std::vector<double>& fSigma, bool RelVals=true, double CorrFrac=1., std::string nature="", UpDnAveraging avg=kSignImprovedQuadratic); //!< Set values of an asymmetric error
   void SetMatrixError(const std::vector<double>& Err, const std::vector<std::vector<double> >& mat, const std::vector<double>& Sigma, bool RelativeVals=true, bool IsCorrelationMatrix=true, std::string nature="") ;
   void SetMatrixError(const std::vector<double>& Err, const TMatrixDSym& mat, const std::vector<double>& Sigma, bool RelativeVals=true, bool IsCorrelationMatrix=true, std::string nature="") ;
   
   // getters/setters
   std::string GetErrorName() const { return fErrorName;} //!< Get error name
   std::string GetErrorSet() const { return fErrorSet;} //!< Get error name
   std::string GetColUp() const { return fColUp;} //!< Get name of column in 'Data'-table of up uncertainties
   std::string GetColDn() const { return fColDn;} //!< Get name of column in 'Data'-table of down uncertainties
   bool GetIsStat() const { return fIsStat;} //!< Get if this is a statistical uncertainty 
   void SetIsStat(bool IsStat=true) { fIsStat=IsStat; fType[1] = fIsStat ? 'S' : 'Y'; } //!< Set this error to be a statistical error   
   bool GetIsTheo() const { return fIsTheo;} //!< Get if this is a statistical uncertainty 
   void SetIsTheo(bool IsTheo=true) { fIsTheo=IsTheo; fType[0] = fIsTheo ? 'T' : 'E'; } //!< Set this error to be a statistical error   
   void SetIsMult(bool IsMult=true) { fIsMult=IsMult; fType[3] = fIsMult ? 'M' : 'A'; } //!< Set this error to multiplicative or additive
   std::string GetType() const { return fType; }//!< return 'type' of error
   void SetIsAdditive(bool IsAdd=true) { SetIsMult(!IsAdd); } //!< Set this error to multiplicative or additive
   bool GetIsMult() const { return fIsMult;} //!< Get if this is a multiplicative or additive error
   void SetColUpDn(const std::string& ColUp, const std::string& ColDn) {fColUp = ColUp; fColDn = ColDn;} //!< Set ColUp and ColDn
   UpDnAveraging GetAveragingPrescription() const { return fUpDnAvg;} //!< Get averaging prescription for asymetric uncertainties
   void SetAveragingPrescription(UpDnAveraging Avg); //!< Set averaging prescription for asymetric uncertainties
   void SetCorrelatedFraction(double corrfrac); //!< Set if error is correalted (1) or uncorrelated (0);
   double GetCorrelatedFraction() const {return fCorrFrac;} //!< Get if error is correalted (1) or uncorrelated (0);
   inline double GetSign(double up, double dn) const { 
      if ( fabs(up)*fabs(dn) == 0 ) return 1;
      return (fabs(up)>fabs(dn)) ? up/fabs(up) : -1*dn/fabs(dn); }
   double GetAvg(double up, double dn) const; //!< Get average of asymmetric uncrtainty following the avereraging prescription
   std::vector<double> GetAvg(const std::vector<double>& up, const std::vector<double> dn) const; //!< Get average of asymmetric uncrtainty following the avereraging prescription
   void CalcAvg(std::vector<double>& avg, const std::vector<double>& up, const std::vector<double> dn) const; //!< Get average of asymmetric uncrtainty following the avereraging prescription
   const std::string& GetNature() const {return fNature;} //!< Get if error is correalted (1) or uncorrelated (0);
   void SetNature(const std::string& nat) { fNature = nat;} //!< Set 'nature' of uncertainty

   bool GetIsMatType() const { return (fCorrFrac<0 || fCorrFrac>1);}//!< Get if error is a 'matrix' type error, i.e. no correlated and uncorrelated errors are reasonably defined
   // Getters for errors and error matrices
   // TODO: remnants of old interface -> remove once these become accessible through new interface
   const std::vector<double>& GetMatError() const { return fMatErrAbs;} //!< Get matrix type uncertainites as uncorrelated
   const std::vector<double>& GetMatErrorRel() const { return fMatErr;} //!< Get matrix type uncertainites as uncorrelated
   const TMatrixDSym& GetCorrelations() const { return fCorr;} //!< Get correlation matrix

   void ApplyPointValidMap(const std::vector<bool>& PointValid); //!< Apply an array to remove unneeded points. Should be called only once (if at all!)!
   const std::vector<double>& GetError(std::string access, bool recalculate=false) const ; //!< Get error [new interface]
   const TMatrixDSym& GetErrorMatrix(std::string access, bool recalculate=false) const ; //!< Get error as matrix [new interface]
   static bool CheckAccess(std::string& access); //!< check 'access' tag for error

   bool RescaleError(const std::vector<double>& sigma_rescale);
   void ClearDerivedErrorsCache(); //!< clear all 'derived' errors from the cache

private:
   bool fIsStat = false; 
   bool fIsTheo = false; 
   bool fIsMult = true; //!< Error is of multiplicative or additive nature 
   bool fRelVals; //!< if the input values have been in percent
   std::string fNature; //!< absolute/multiplicative/poisson/... 
   std::string fErrorSet; //!< the error set
   std::string fErrorName; //!< Name of the error
   std::string fColUp; //!< label of the column in 'Data' for the up-error
   std::string fColDn; //!< label of the column in 'Data' for the down-error
   
   // error values
   std::vector<double> fErrUp; //!< The up error (relative)
   std::vector<double> fErrDn; //!< The up error (relative)
   std::vector<double> fSigma; //!< cross section
   UpDnAveraging fUpDnAvg = kSignImprovedQuadratic; //!< Average prescription
   double fCorrFrac = 1; //!< correlated fracton (0-1), use -1, if correlation matrix is used

   // TODO: remnants of old interface -> remove once these become accessible through new interface
   std::vector<double> fMatErr; //!< The (uncorrelated) error of 'matrix' type uncertainties (relative)
   std::vector<double> fMatErrAbs; //!< The (uncorrelated) error of 'matrix' type uncertainties (absolute)

   TMatrixDSym fCorr; //!< correlation matrix
   
   // new error interface
   mutable std::map<std::string,  TMatrixDSym > fErrorMats; //!< Member to store matrix error's
   mutable std::map<std::string, std::vector<double> > fErrors; //! new member store vector error's 
   std::string fType = "EYNM"; //!< Error type: {[E,T],[S,Y],[C,N],[A,M]}: [E,T]=Exp|Theo, [S,Y]=Stat|Sys, [C,N]=MatrixType|NoMatriType, [A,M]=Additive|Multiplicative

protected:

   void AbsToRel(std::vector<double>& inp) {using namespace AlposTools; inp /= fSigma;};
   void RelToAbs(std::vector<double>& inp) {using namespace AlposTools; inp *= fSigma;};
   

};

#endif
