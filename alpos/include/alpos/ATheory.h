// DB. 12/2014, 01/2015
#ifndef Alpos_ATheory
#define Alpos_ATheory

/** 
 * ATheory header
 *
 * contains class declarations and implementations for
 * objects representing Alpos theory parameters.
 * 
 */

#include <map>
#include <set>
#include <string>
#include <stdarg.h>
#include <sstream>
#include <vector>

#include "alpos/Alpos.h"
#include "alpos/ATheoryHandler.h"
#include "alpos/AlposObject.h"
#include "alpos/AlposTools.h"
#include "alpos/AError.h"
#include "alpos/AParmNamed.h"

//! getters and setters for Parameters/Constants/Functions values
// ---- Getters (Get the value of parameter or function, where this this function depends on):
//! return (double) value of a parameter
//! usage: e.g. PAR(Alpha_s) (i.e. no quotes around 'Alpha_s')
//! can only be used within a function
//! PAR(X) is identical to the usage of: VALUES(X)[0]
#define PAR(X)           TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string(#X))->GetValue()
#define PAR_S(X)         TheoryHandler::Handler()->GetParmS(this->GetAlposName()+std::string(".")+std::string(#X))->GetValue().c_str()
//! return (vector<double>) values of a function (can only be used within a function)
#define VALUES(X)        TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetValues()
#define UPDATE(X)        PAR(X)

//! Check if an update of this requested from this 'member' parameter (can only be used within a function)
//! Return true, if update is needed
#define CHECK(X)         CheckUpdateRequested(std::string(#X))

//! Set paramter to constant. E.g. if it is not considered in 'Update' and hence should not be changed.
#define CONST(X)         TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string(#X))->SetIsConstant()

//! returns 'quick' values of a function. It is recommended that before 'QUICK' is called frequently (i.e. within an integration) outside of the integration first call VALUES(X) to ensure that the values are up-to-date.
//! MIND: QUICK(Alpha_s,Q) will not work, but use: QUICK(Alpha_s,{Q}) or QUICK_VAR(Alpha_s,1,Q);
#define QUICK_VAR(X, N, ...)    TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetQuick(int(N),__VA_ARGS__)
// //#define QUICK(X,Y)         TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetQuick(Y)
#define BRACED_INIT_LIST(...) {__VA_ARGS__}
#define QUICK(X, A)    TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetQuick(BRACED_INIT_LIST A)
#define QUICK_VEC(X, A)    TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".")+std::string(#X))->GetQuick(vector<double>({A}))

// ---- Getters for 'any' parameters (should not be used in principle)
// usage: e.g. PAR("Alpha_s") (i.e. WITH quotes around 'Alpha_s')
#define PAR_ANY(X)      TheoryHandler::Handler()->GetParmD(X)->GetValue()
#define PAR_ANY_S(X)     TheoryHandler::Handler()->GetParmS(X)->GetValue().c_str()
#define VALUES_ANY(X)    TheoryHandler::Handler()->GetFuncD(X)->GetValues()
#define QUICK_ANY(X, Y)     TheoryHandler::Handler()->GetFuncD(X)->GetQuick(Y)

// ---- Setters for parameters where this function depends on
//! set (double) parameter value
//! usage: e.g. SET(Alpha_s,0.1185,0.0006),  (i.e. no quotes around 'Alpha_s')
#define SET(X, VAL, ERR)   TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string(#X))->SetValue(VAL,ERR,false);
//! set (string) parameter value
#define SET_S(X, VAL, ERR) TheoryHandler::Handler()->GetParmS(this->GetAlposName()+std::string(".")+std::string(#X))->SetValue(VAL,ERR,false);

// ---- Setters for 'any' parameters (should not be used in principle) 
//! usage: e.g. PAR("fnlo1.Alpha_s.Mz",91.16,0.01),  (i.e. WITH quotes around 'fnlo1.Alpha_s.Mz')
#define SET_ANY(X, VAL, ERR)   TheoryHandler::Handler()->GetParmD(X)->SetValue(VAL,ERR,false);
#define SET_ANY_S(X, VAL, ERR) TheoryHandler::Handler()->GetParmS(X)->SetValue(VAL,ERR,false);


/**********************************************************

  virtual template class AParmBase
  
  basic class to hold a theory parameter of type T.
  The specializations are available for int, double and string
  
  The class contains a value, error and a flag fIsConst if the
  value is changeable or not.

  The AParmBase holds a list of dependent AParmFunction 
  which are notified about a change if the parameter
  changes.

***********************************************************/

// ___________________________________________________________________________________________________
template<typename T>
class AParmBase : public AParmNamed {

public:
   typedef std::vector<T> vecT;
   AParmBase(const std::string& name);

   virtual ~AParmBase() { ; };

   // -> to AParmConstD void SetValue(T value, T error, bool force=false); //! Check if IsConst
   const vecT& GetValues() {
      if (!GetIsUpToDate()) {
         bool success = Update();
         if (!success)
            error["GetValues"] << "Update() of '" << GetAlposName() << "' was not successfull!" << std::endl;
         else {
            SetIsOutdated(nullptr, false);
            ApplyThFactors();
         }
      }
      return fValue;
   }

   const vecT& GetErrors() {
      if (!GetIsUpToDate()) {
         Update();
         SetIsOutdated(nullptr, false);
      }
      return fError;
   }

   const T& GetValue() { return GetValues()[0]; }

   const T& GetError() { return GetErrors()[0]; }

   void SetIsConstant(bool IsConst = true) { fIsConst = IsConst; }; //!< Set parameter to constant
   bool GetIsConstant() const { return fIsConst; } //!< Get if parameter is constant
   unsigned int N() const { return fValue.size(); } //!< Get size of values (i.e. number of bins)

   void ApplyThFactors(); //!< Apply the k-factors
   void AddThFactor(const std::string& kname, const T& thfac);

   //!< Add a k-factor
   void AddThFactor(const std::string& kname, const vecT& thfac);//!< Add a k-factor

   virtual void GetContent(std::vector<std::string>& val, std::vector<std::string>& err,
                           bool& konst); //!< Get the values, errors and const-flag as a strings
   virtual void GetContent(AThCont& cont) { GetContent(cont.Values, cont.Errors, cont.Const); }

   virtual void SetContent(const std::vector<std::string>& val, const std::vector<std::string>& err,
                           bool konst); //!< Get the values, errors and const-flag from strings
   virtual void SetContent(const AThCont& cont) { SetContent(cont.Values, cont.Errors, cont.Const); }

private:
   bool fIsConst = false;

protected:
   virtual bool Update() = 0;

   void SetValue(const T& val, const T& err);

   void SetValues(const vecT& val, const vecT& err);

   vecT fValue;
   // = {0};
   vecT fError;
   // = {0};
   std::map<std::string, vecT> fThFactors;
};


// ___________________________________________________________________________________________________
template<typename T>
AParmBase<T>::AParmBase(const std::string& name) : AParmNamed(name) {
   SetClassName("AParmBase");
   fValue.resize(1);
   fError.resize(1);
   // Call: TheoryHandler and Register, (vecT value,vecT error, bool IsConst=false);
   //std::cout<<"AParmBase<T>::AParmBase. Todo: Call: TheoryHandler and Register, (vecT value,vecT error, bool IsConst=false)"<<std::endl;
};


// ___________________________________________________________________________________________________
template<typename T>
void AParmBase<T>::SetValue(const T& val, const T& err) {
   if (fValue[0] == val && fError[0] == err) {
      // nothing to do
   }
   else {
      if (this->GetIsConstant()) {
         if (fValue[0] != val || fError[0] != err) {
            error["SetContent"] << "Attempt to change constant parameter (par='" << this->GetAlposName() <<
            "'). Exiting" << std::endl;
            exit(702);
         }
      }

      fValue[0] = val;
      fError[0] = err;
      SetIsOutdated(this); // notify dependencies.
      SetIsOutdated(nullptr, false); // a parameter is never outdated.
   }
};


// ___________________________________________________________________________________________________
template<typename T>
void AParmBase<T>::SetValues(const std::vector<T>& val, const std::vector<T>& err) {
   //! simple setter for all 'values'
   if (val.size() != err.size()) {
      error["SetValues"] << "val and err must be of same size." << std::endl;
   }

   for (unsigned int i = 0; i < val.size(); i++) {
      if (this->GetIsConstant()) {
         if (fValue[i] != val[i]) {
            error["SetContent"] << "Attempt to change constant parameter (par='" << this->GetAlposName() <<
            "'). Exiting" << std::endl;
            exit(701);
         }
      }
      fValue[i] = val[i];
      fError[i] = err[i];
      SetIsOutdated(this); // notify dependencies.
      SetIsOutdated(nullptr, false); // a parameter is never outdated.
   }
};


// ___________________________________________________________________________________________________
template<typename T>
void AParmBase<T>::GetContent(std::vector<std::string>& val, std::vector<std::string>& err, bool& konst) {
   // a getter using 'strings' as values
   const std::vector<T>& gval = this->GetValues();
   const std::vector<T>& gerr = this->GetErrors();
   // val.resize(gval.size());
   // err.resize(gerr.size());
   val.clear();
   err.clear();
   konst = this->GetIsConstant();
   for (auto i : gval) {
      std::stringstream ss;
      ss << i;
      //val[i] = ss.str();
      val.push_back(ss.str());
   }
   for (auto i : gerr) {
      std::stringstream ss;
      ss << i;
      //err[i] = ss.str();
      err.push_back(ss.str());
   }
}


// ___________________________________________________________________________________________________
template<typename T>
void AParmBase<T>::SetContent(const std::vector<std::string>& val, const std::vector<std::string>& err, bool konst) {
   // a setter using 'strings' as values
   if (val.size() != err.size()) {
      error["SetContent"] << "val and err must be of same size." << std::endl;
   }
   std::vector<T> vV(val.size()), vE(val.size());
   for (unsigned int i = 0; i < val.size(); i++) {
      std::stringstream sse, ssv;
      ssv << val[i];
      sse << err[i];
      ssv >> vV[i];
      ssv >> vE[i];
   }
   this->SetIsConstant(konst);
   this->SetValues(vV, vE);
}


// ___________________________________________________________________________________________________
//! A theory parameter (i.e. a constant or variable to a function) for Alpos
template<typename T>
//class AParm : public virtual AParmBase<T> {
class AParm : public AParmBase<T> {
public:
   AParm(const std::string& name) : AParmBase<T>(name) { PrimalScream::SetClassName("AParm"); };

   AParm(const std::string& name, const T& value, const T& error, bool IsConst);

   virtual ~AParm() { };

   virtual bool Update() {
      // dummy implementation. Update() should never be called on a parameter.
      PrimalScream::debug["Update"] <<
      "'Update()' should not be called on a constant/variable since it does not depend on anything, but was called on " <<
      this->GetAlposName() << " (or steering parameter is '0')" << std::endl;
      return true;
   };

   bool SetValue(const T& value, const T& error, bool force = false);

};

// ___________________________________________________________________________________________________
template<typename T>
AParm<T>::AParm(const std::string& name, const T& value, const T& error, bool IsConst) : AParm(name) {
   SetValue(value, error, true);
   AParmBase<T>::SetIsConstant(IsConst);
};


// ___________________________________________________________________________________________________
template<typename T>
bool AParm<T>::SetValue(const T& value, const T& error, bool force) {
   if ((AParmBase<T>::GetIsConstant() && !force) || !AParmBase<T>::GetIsConstant()) {
      AParmBase<T>::SetValue(value, error);
      return true;
   }
   else {
      PrimalScream::warn["SetValue"] << "Parm '" << AParmBase<T>::GetAlposName() <<
      "' is constant. Not allowed to change it. Ignoring call!" << std::endl;
      exit(703);
      return false;
   }
};

// ___________________________________________________________________________________________________


/**********************************************************
  virtual template class AParmFuncBase
  
  basic virtual class to hold a theory function which uses
  AParm or AParmFuncBases as input.
  
  The function calculates the values 'fValue' and may
  also provide errors.

  If it is not allowed to change the input parameters, set
  it to const.

  The AParmFuncBase holds a list of dependent AParmFunc 
  which are notified about a change if the parameter
  changes.

***********************************************************/

// ___________________________________________________________________________________________________
template<typename T>
class AParmFuncBase : public AParmBase<T> {

public:
   AParmFuncBase(const std::string& name) : AParmBase<T>(name) { PrimalScream::SetClassName("AParmFuncBase"); };

   virtual ~AParmFuncBase() { ; };

   void SetNeedsUpdate(AParmNamed* from = nullptr /*std::string& from*/) {
      AParmBase<T>::SetIsOutdated(from);
      //AParmBase<T>::NotifyDependencies();
   }

   virtual bool Init() = 0; //!< Initialize the function. Init is once called for each function.
   //void NotifyDependencies() {/*todo*/};
   virtual bool Update() = 0;
   // todo:   virtual void RegisterDependency(std::string name,AParmFuncFunction<T>* pFunc) = 0;
   //virtual void RegisterDependency(const std::string& name,AParmFuncBase<T>* pFunc) = 0;

   virtual std::vector<std::string> GetRequirements() const = 0; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const = 0; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const = 0; //!< The function's name, how it can be initialized
   virtual std::vector<T> GetQuick(int,
                                   ...) = 0; //!< The possibilty to implement a quick access without changing of any parameters or recent (member-)values
   virtual std::vector<T> GetQuick(
         const std::vector<double>&) = 0; //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   // -- error source adders/setters

   bool AddSymError(const std::string name, const std::vector<double>& Err, const std::vector<double>& fSigma,
                    bool RelVals = true, double CorrFrac = 1., std::string nature = "") {
      this->debug["AddSymError"] << "Adding symmetric non-matrix-type error source '" << name << "' to '" << this->GetAlposName() << "'..." << std::endl;
      if (this->HasErrorNamed(name)) this->warn["AddSymError"] << "Not adding error source." << std::endl;
      // add error to list of errors
      if ( Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
	 // ignore theory error
	 this->info["AddSymError"]<<"Adding matrix-type error source '" << name << "' skipped, because option 'IgnoreTheoryErrors' set to true."<<std::endl;
	 return true;
      }
      AError tmp(name);  // TODO: how to handle ErrorSet?
      fAllErrors[tmp.GetErrorName()] = tmp;
      // initialize error
      fAllErrors[tmp.GetErrorName()].SetSymError(Err, fSigma, RelVals, CorrFrac, nature);
      return true;
   }; //!< Add a symmetric error to the function/dataset

   bool AddAsymError(const std::string name, const std::vector<double>& ErrUp, const std::vector<double>& ErrDn,
                     const std::vector<double>& fSigma, bool RelVals = true, double CorrFrac = 1.,
                     std::string nature = "", AError::UpDnAveraging avg = AError::kSignImprovedQuadratic) {
      this->debug["AddAsymError"] << "Adding asymmetric non-matrix-type error source '" << name << "' to '" << this->GetAlposName() << "'..." << std::endl;
      if (this->HasErrorNamed(name)) this->warn["AddAsymError"] << "Not adding error source." << std::endl;
      if ( Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
	 // ignore theory error
	 this->info["AddAsymError"]<<"Adding matrix-type error source '" << name << "' skipped, because option 'IgnoreTheoryErrors' set to true."<<std::endl;
	 return true;
      }
      // add error to list of errors
      AError tmp(name);  // TODO: how to handle ErrorSet?
      fAllErrors[tmp.GetErrorName()] = tmp;
      // initialize error
      fAllErrors[tmp.GetErrorName()].SetAsymError(ErrUp, ErrDn, fSigma, RelVals, CorrFrac, nature, avg);
      return true;
   }; //!< Add an asymmetric error to the function/dataset

   bool AddMatrixError(const std::string name, const std::vector<double>& Err,
                       const std::vector<std::vector<double>>& mat, const std::vector<double>& Sigma,
                       bool RelativeVals = true, bool IsCorrelationMatrix = true, std::string nature = "") {
      this->debug["AddMatrixError"] << "Adding matrix-type error source '" << name << "' to '" << this->GetAlposName() << "'..." << std::endl;
      if (this->HasErrorNamed(name)) this->warn["AddMatrixError"] << "Not adding error source." << std::endl;
      if ( Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
	 // ignore theory error
	 this->info["AddMatrixError"]<<"Adding matrix-type error source '" << name << "' skipped, because option 'IgnoreTheoryErrors' set to true."<<std::endl;
	 return true;
      }
      // add error to list of errors
      AError tmp(name);  // TODO: how to handle ErrorSet?
      fAllErrors[tmp.GetErrorName()] = tmp;
      // initialize error
      fAllErrors[tmp.GetErrorName()].SetMatrixError(Err, mat, Sigma, RelativeVals, IsCorrelationMatrix, nature);
      return true;
   }; //!< Add a matrix-type error to the function/dataset (vector<vector<double>>)

   bool AddMatrixError(const std::string name, const std::vector<double>& Err, const TMatrixDSym& mat,
                       const std::vector<double>& Sigma, bool RelativeVals = true, bool IsCorrelationMatrix = true,
                       std::string nature = "") {
      this->debug["AddMatrixError"] << "Adding matrix-type error source '" << name << "' to '" << this->GetAlposName() << "'..." << std::endl;
      if (this->HasErrorNamed(name)) this->warn["AddMatrixError"] << "Not adding error source." << std::endl;
      if ( Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
	 // ignore theory error
	 this->info["AddMatrixError"]<<"Adding matrix-type error source '" << name << "' skipped, because option 'IgnoreTheoryErrors' set to true."<<std::endl;
	 return true;
      }
      // add error to list of errors
      AError tmp(name);  // TODO: how to handle ErrorSet?
      fAllErrors[tmp.GetErrorName()] = tmp;
      // initialize error
      fAllErrors[tmp.GetErrorName()].SetMatrixError(Err, mat, Sigma, RelativeVals, IsCorrelationMatrix, nature);
      return true;
   }; //!< Add a matrix-type error to the function/dataset (ROOT matrix)

   bool HasErrorNamed(const std::string name) {
      if (fAllErrors.count(name) > 0) {
         this->warn["HasErrorNamed"] << "An error source with the same name already exists: '" << name << "'!" << std::endl;
         return true;
      }
      return false;
   }

   // -- uncertainty-related methods (moved up from AData)
   const std::map<std::string, AError>& GetAllErrors() const { return fAllErrors; }

   // TODO: remnants of old interface -> remove once these become accessible through new interface
   const std::vector<double>& GetUncertaintyMat() const { return fErrMat; }; //!< Return (quadratic) sum of all matrix-type uncertainties (except stat)
   const std::vector<double>& GetUncertaintyMatRel() const { return fErrMatRel; }; //!< Return (quadratic) sum of all matrix-type uncertainties (except stat) as relative uncertainty    //const vector<double> GetUncertaintyCorr(); //!< this does not make sense
   //const vector<double> GetUncertaintyCorrRel(); //!< this does not make sense

   bool HasErrors() { return !fAllErrors.empty(); };

   void ClearErrorCache() {
      this->fSumErrors.clear();
      this->fSumErrMats.clear();
      for (auto& err : this->fAllErrors) 
	 err.second.ClearDerivedErrorsCache();
   }

   void RescaleErrors() {
      //! Rescale uncertainties to current function values
      /*!
       *  This rescales all uncertainties registered for this object
       *  to the current value vector. The absolute errors/error matrices
       *  are recalculated accordingly.
       *
       *  \warning This assumes that the relative uncertainties do not change
       *           significantly. This method should not be used if this is
       *           not the case.
       *
       *  \note Rescaling uncertainties invalidates the 'sum' matrix cache.
       *
       *  \see AError::RescaleError()
       */

      for (auto& err : this->fAllErrors) {
         this->debug["RescaleUncertainties"] << "Rescaling uncertainty '" << err.first << "'" << std::endl;
         err.second.RescaleError(this->fValue);
      }

      // clear 'sum' error cache
      this->fSumErrors.clear();
      this->fSumErrMats.clear();
   };

   //____________________________________________________________________________________ //
   static bool CheckType(std::string& type) {
      //! Check 'type' to error
      //! add [Tot] and [Reg] if not present

      if (type.size() != 2) {
         say::warn["AParmFuncBase::CheckType"] << "'type' has wrong size, must be 2, but is: " << type.size() <<
         std::endl;
         return false;
      }

      if (type[0] != 'E' && type[0] != 'T' && type[0] != 'A') {// faster
         say::warn["AParmFuncBase::CheckType"] << "First first letter must be 'E', 'T' or 'A', but is: " << type[0] <<
         std::endl;
         return false;
      }

      if (type[1] != 'S' && type[1] != 'Y' && type[1] != 'A') {// faster
         say::warn["AParmFuncBase::CheckType"] << "First first letter must be 'S', 'Y' or 'A', but is: " << type[1] <<
         std::endl;
         return false;
      }
      return true;
   }


// __________________________________________________________________________________________ //
   const std::vector<double>& GetSumError(std::string type, std::string access) const {
      //!< Access to errors
      CheckType(type);
      AError::CheckAccess(access);
      if (!fSumErrors.count(type + access)) CalculateSumError(type + access);
      return fSumErrors[type + access];
   }


// __________________________________________________________________________________________ //
   const TMatrixDSym& GetSumErrorMatrix(std::string type, std::string access) const {
      //!< Access to errors
      CheckType(type);
      AError::CheckAccess(access);
      if (!fSumErrMats.count(type + access)) CalculateSumErrMat(type + access);
      return fSumErrMats[type + access];
   }

   bool HasMultErrors() { return fHasMultErrors; };  //!< True if at least one error source is multiplicative

   void SetHasMultErrors(
         bool has_mult = true) { fHasMultErrors = has_mult; };  //!< Set whether at least one error source is multiplicative

protected:
   bool CheckUpdateRequested(const std::string& aname) const {
      using namespace std;
      //! check whether this particluar parameter has changed
      //! returns 'true' if updated is needed!
      const std::set<AParmNamed*> nots = this->GetCurrentNotifiers();
      std::string lname = this->GetAlposName() + "." + aname;
      if (nots.find(nullptr) != nots.end()) return true; // always update, if a nullptr is present.
      if (nots.find(TheoryHandler::Handler()->GetParameter(lname)) != nots.end()) {
         return true;
      }
      return false;
   }

   void PrintCurrentNotifiers() const {
      PrimalScream::debug["PrintCurrentNotifiers"] << "The function '" << this->GetAlposName() <<
      "' is currently notified for an update by:" << std::endl;
      for (auto i :  this->GetCurrentNotifiers())
         if (i != nullptr)
            PrimalScream::debug["PrintCurrentNotifiers"] << "  + " << i->GetAlposName() << std::endl;
         else
            PrimalScream::debug["PrintCurrentNotifiers"] << "  + 'nullptr'" << std::endl;
   }

   // -- uncertainty-related data members (moved up from AData)

   bool fHasMultErrors = false;  //!< True if at least one error source is multiplicative

   std::map<std::string, AError> fAllErrors; //!< all individual error sources
   //mutable std::map<std::string,AError> fSumErrors; //!< sum of individual errors
   //void CalculateSumError(std::string type ) const;
   //void CalculateSumErrMat(std::string type ) const;

   mutable std::map<std::string, std::vector<double> > fSumErrors; //!< sum of individual errors
   mutable std::map<std::string, TMatrixDSym> fSumErrMats; //!< sum of individual errors

   // TODO: remnants of old interface -> remove once these become accessible through new interface
   std::vector<double> fErrMat;
   std::vector<double> fErrMatRel;

   // __________________________________________________________________________________________ //
   void CalculateSumError(std::string actype) const {
      //! calculate summed errors
      //! Options:
      //!
      //!           1-2    3-5     6-7    8-10    11-13
      //! actype = 'XY' + 'Aaa' + 'Bb' + 'Ccc' + 'Ddd'
      //! with:
      //!       X: 'E'/'T': experimental/theoretical uncertainties only. 'A': all
      //!       Y: 'S'/'Y': statistical/systematic uncertainties only. 'A': all
      //!     Aaa: 'Abs'/'Rel': absolute or relative uncertainties
      //!      Bb: 'Av'/'Up'/'Dn': symmetrized, asymmetric up or down uncertainties
      //!     Ccc: 'Tot'/'Cor'/'Unc': total uncertainty, correlated fraction of uncertainty, uncorrelated fraction
      //!     Ddd: 'Reg'/'Inv': regular or inverse uncertainties
      //!
      //! Warning: These partially are not entirely reaonsable
      //! For instance, if up or dn errors are calculated
      //! the 'sign' is lost
      //! Or e.g. if the correlated stat. uncertainty is requested.

      using namespace AlposTools;

      //CheckType(type);
      if (fSumErrors.count(actype)) {
         this->warn["CalculateSumError"] << "Function should be called only once for each uncertainty. actype=" <<
         actype << std::endl;
      }

      fSumErrors[actype].resize(this->fValue.size());
      // fSumErrors[type] = AError(type,"");
      // fSumErrors[type].SetIsStat(type[1]=='S');
      // //fSumErrors[type].SetIsTheo(type[1]=='T'); // todo
      // if ( type[1]=='T'  ){
      //    cout<<"Error. Theoreical errors not yet implements 34238947"<<endl;
      //    exit(1);
      // }

      //std::vector<double> e2(this->fValue.size());
      for (const auto& ierr : fAllErrors) {
         bool UseSYA = false;
         if (actype[1] == 'S' && ierr.second.GetIsStat()) UseSYA = true; // stat
         else if (actype[1] == 'Y' && !ierr.second.GetIsStat()) UseSYA = true;// sys
         else if (actype[1] == 'A') UseSYA = true;

         bool UseETA = false;
         if (actype[0] == 'E' && !ierr.second.GetIsTheo()) UseETA = true; // stat
         else if (actype[0] == 'T' && ierr.second.GetIsTheo()) UseETA = true;// sys
         else if (actype[0] == 'A') UseETA = true;

         if (UseETA && UseSYA) {
            //this->debug["CalculateSumErrors"] << "Adding error source '" << ierr.first <<
            //"' to sum fSumErrors[actype='" << actype << "']" << std::endl;
            std::vector<double> e2 = ierr.second.GetError(actype.substr(2, 11));
            e2 *= e2;
            for (unsigned int i = 0; i < this->fValue.size(); i++) fSumErrors[actype][i] += e2[i]; // total
            //this->debug["CalculateSumErrors"] << "fSumErrors[actype='" << actype << "'] = {" <<
            //sqrt(fSumErrors[actype][0]) << ", " << sqrt(fSumErrors[actype][1]) << ", ...}" << std::endl;
         }
      }
      for (auto& ie : fSumErrors[actype]) ie = sqrt(ie);

      // if ( actype.substr(10,3)=="Inv") {
      //    string accReg = actype;
      //    accReg.replace(10,3,"Reg");
      //    fSumErrors[accReg] = fSumErrors[actype];
      //    for ( auto& ie : fSumErrors[actype] ) 1./ie;
      // }

   }


   // __________________________________________________________________________________________ //
   void CalculateSumErrMat(std::string actype) const {
      //! calculate summed error matrix
      //! Options:
      //!
      //!           1-2    3-5     6-7    8-10    11-13
      //! actype = 'XY' + 'Aaa' + 'Bb' + 'Ccc' + 'Ddd'
      //! with:
      //!       X: 'E'/'T': experimental/theoretical uncertainties only. 'A': all
      //!       Y: 'S'/'Y': statistical/systematic uncertainties only. 'A': all
      //!     Aaa: 'Abs'/'Rel': absolute or relative uncertainties
      //!      Bb: 'Av'/'Up'/'Dn': [symmetrized, asymmetric up or down uncertainties] -*- not used for matrices -*-
      //!     Ccc: 'Tot'/'Unc': whole covariance matrix or just diagonal  -*- for Ccc == 'Cor', see Warning2 below -*-
      //!     Ddd: 'Reg'/'Inv': regular or inverse covariance matrices
      //!
      //! Warning: These partially are not entirely reaonsable
      //! For instance, if up or dn errors are calculated
      //! the 'sign' is lost
      //! Or e.g. if the correlated stat. uncertainty is requested.
      //!
      //! Warning2: For Ccc == 'Cor', the error sum will exclude any
      //! 'matrix'-type error source, as well as any non-matrix error
      //! source whose correlation coefficient is 0!

      using namespace AlposTools;

      //CheckType(type);
      if (fSumErrMats.count(actype)) {
         this->warn["CalculateSumErrMat"] << "Function should be called only once for each uncertainty. actype=" <<
         actype << std::endl;
      }
      using namespace std;
      fSumErrMats[actype].ResizeTo(this->N(), this->N());
      fSumErrMats[actype].Zero();

      for (const auto& ierr : fAllErrors) {
         bool UseSYA = false;
         if (actype[1] == 'S' && ierr.second.GetIsStat()) UseSYA = true; // stat
         else if (actype[1] == 'Y' && !ierr.second.GetIsStat()) UseSYA = true; // sys
         else if (actype[1] == 'A') UseSYA = true;

         bool UseETA = false;
         if (actype[0] == 'E' && !ierr.second.GetIsTheo()) UseETA = true; // exp
         else if (actype[0] == 'T' && ierr.second.GetIsTheo()) UseETA = true; // theo
         else if (actype[0] == 'A') UseETA = true;

         if (UseETA && UseSYA) {
            this->debug["CalculateSumErrMat"] << "Adding error source '" << ierr.first << "' to sum fSumErrMats[actype='" << actype << "']" << std::endl;
            const TMatrixDSym& mat = ierr.second.GetErrorMatrix(actype.substr(2, 11)); // total
            if (fSumErrMats[actype].GetNcols() != mat.GetNcols() || fSumErrMats[actype].GetNrows() != mat.GetNrows() ) {
               this->error["CalculateSumErrMat"] << "fSumErrMats[actype='" << actype << "'] += err.GetErrorMatrix(" <<
		  actype.substr(2, 11) << "): Ncols/Nrows mismatch  nCols " << fSumErrMats[actype].GetNcols() << " vs " <<
		  mat.GetNcols() << "\tnRows " << fSumErrMats[actype].GetNrows() << " vs " <<
		  mat.GetNrows() <<std::endl;
	       cout<<"This is a bug! Fix it."<<endl;
	       //const TMatrixDSym& mat2 = ierr.second.GetErrorMatrix("AbsDnTot"); // total
	       //cout<<"Problem in function: "<<this->GetAlposName()<<"\t nReq= "<<this->GetRequirements().size()<<"\t errorsigmasize: "<<ierr.second.fSigma.size()<<endl;
	       //cout<<"fErrorssize: " <<ierr.second.fErrors.size()<<"\terrmat="<<ierr.second.fErrorMats.size()<<"\tfErrUp="<<ierr.second.fErrUp.size()<<"\tfErrDn.size()="<<ierr.second.fErrDn.size()<<"\t ferrstd= "<<ierr.second.fErrors["RelDnTotAll"].size()<<"\tmat2="<<mat2.GetNrows()<<endl;
	       //const TMatrixDSym& mat3 = ierr.second.GetErrorMatrix("RelAvTotAll"); // total  
	       //cout<<"   mat3: "<<mat.GetNrows()<<endl;
            }
	    fSumErrMats[actype] += ierr.second.GetErrorMatrix(actype.substr(2, 11));   // total
	    this->debug["CalculateSumErrMat"] << "fSumErrMats[actype='" << actype << "'] = {{" <<
	       fSumErrMats[actype](0, 0) << ", " << (this->N()>1?fSumErrMats[actype](0, 1):0) << ", ...}...}" << std::endl;
	 }
         else {
            this->debug["CalculateSumErrMat"] << "NOT adding error source '" << ierr.first <<
            "' to sum fSumErrMats[actype='" << actype << "'" << std::endl;
         }
      }
      // if ( actype.substr(10,3)=="Inv") {
      //    string accReg = actype;
      //    accReg.replace(10,3,"Reg");
      //    fSumErrMats[accReg].ResizeTo(fSumErrMats[actype]);
      //    fSumErrMats[accReg] = fSumErrMats[actype];
      //    fSumErrMats[actype] = AlposTools::InvertChol(fSumErrMats[accReg]);
      // }
      
   }

};


// ___________________________________________________________________________________________________
// A global template function.
// this could also go as a member funtion to ATheoryHandler.
template<typename T>
bool ARegisterRequirements(T* obj) {
   using namespace say;
   for (auto ir : obj->GetRequirements()) {
      debug["ARegisterRequirements"] << "Register dependency(" << obj->GetAlposName() << "," << ir << ")" << std::endl;
      std::string parname = obj->GetAlposName() + "." + ir;
      std::string alname = obj->GetFunctionName() + "." + ir;
      // if ( !TheoryHandler::Handler()->CheckParameter(alname) ) {
      // 	 error["ARegisterRequirements"]<<"Default parameter with name '"<<alname<<"' does not exist! Exiting."<<std::endl;
      // 	 exit(1);
      // 	 return false;
      // }
      debug["ARegisterRequirements"] << " + ARegisterRequirements(). Register dependency of " << obj->GetAlposName() <<
      " at " << parname << "." << std::endl;
      // if parameter does exist:       register dependency
      // if parameter does NOT exist:   register new alias with name
      // NO: (if parameter does NOT exist:   register at 'functype.par')
      if (TheoryHandler::Handler()->CheckParameter(parname)) {
         debug["ARegisterRequirements"] << "    Ok! Parameter '" << parname << "' does exist." << std::endl;
         TheoryHandler::Handler()->GetParameter(parname)->RegisterDependency(obj->GetAlposName(), obj);
      }
      else {
         if (!TheoryHandler::Handler()->CheckParameter(alname)) {
            error["ARegisterRequirements"] << "Default parameter with name '" << alname <<
            "' does not exist! Exiting." << std::endl;
            exit(1);
            return false;
         }
         debug["ARegisterRequirements"] << "    Ok! Parameter '" << parname <<
         "' does not exist. Registering new alias." << std::endl;
         bool success = TheoryHandler::Handler()->NewAlias(parname, alname);
         //TheoryHandler::Handler()->GetParameter(parname)->RegisterDependency(obj->GetAlposName(),obj);
         TheoryHandler::Handler()->GetParameter(alname)->RegisterDependency(obj->GetAlposName(), obj);
         if (!success) return false;
      }
   }
   return true;
};


//____________________________________________________________________________________ // 
template<typename T>
void AParmBase<T>::ApplyThFactors() {
   //!< Apply the k-factors  
   using namespace AlposTools;
   //if (fThFactors.size() != 0)
   //   this->debug["ApplyThFactors"] << "Applying " << fThFactors.size()<<" sets of theory factors:" << std::endl;
   for (const auto& i : fThFactors) {
      //this->debug["ApplyThFactors"] << "Applying factor set {" << i.second[0]<<", ...} "
      //                              << "to values {" << fValue[0] << ", ...} ..." << std::endl;
      fValue *= i.second;
      fError *= i.second;
   }
}


//____________________________________________________________________________________ //
template<typename T>
void AParmBase<T>::AddThFactor(const std::string& kName, const T& thfac) {
   //!< Add a k-factor   
   Update(); // to initialize the size of fValue    
   std::vector<T> k(fValue.size(), thfac);
   //for ( auto& i: k) i = thfac;
   AddThFactor(kName, k);
}


//____________________________________________________________________________________ //
template<typename T>
void AParmBase<T>::AddThFactor(const std::string& kName, const vecT& thfac) {
   //!< Add a k-factor
   fThFactors[kName] = thfac;
}

//____________________________________________________________________________________ //                                      



#endif
