// DB 15.01.2015
#ifndef Alpos_ASuperData
#define Alpos_ASuperData

/* 
 *  
 * ASuperData
 *  
 * The class, which constructs one large array of all data points
 * and calculates the 'common' errors
 * 
 */


#include "alpos/ATheory.h"
#include <vector>
#include "alpos/AData.h"
#include <TMatrixDfwd.h>
#include <TMatrixT.h>

class ASuperData : public AData {

public:
   ASuperData(const std::string& name);
   virtual ~ASuperData();

   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fDaRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name
   
   // --- SuperData functions:
   void AddRequirement(const std::string& req) { fDaRequirements.push_back(req);} //!< Add one requirement (Requirements still have to be registered) 
   void SetRequirements(const std::vector<std::string>& reqs) {fDaRequirements=reqs;} //!< Set the requirements (Requirements still have to be registered) 

   // access to data tables
   const std::vector<AData*> GetChildren() const {return fChildren;}
   //const std::map<std::string, std::vector<double> >& GetDataTable(int iDataset) const {return *fOrigDataTables[iDataset];} //FIXME: is this needed?

   // // --- convenient covariance matrices
   // const TMatrixDSym& GetCovariance() const { return fCov; } //!< Get covariance matrix of all uncertainties
   // const TMatrixDSym& GetCovarianceRel() const { return fCovRel; } //!< Get covariance matrix of all relative uncertainties
   // const TMatrixDSym& GetCovarianceStat() const { return fStat; } //!< Get covariance matrix of statistical uncertainties
   // const TMatrixDSym& GetCovarianceStatRel() const { return fStatRel; } //!< Get covariance matrix of relative statistical uncertainties
   // const TMatrixDSym& GetCovarianceStatUncorr() const { return fStatUncor; } //!< Get covariance matrix of uncorrelated and statistical uncertainities
   // const TMatrixDSym& GetCovarianceStatUncorRel() const { return fStatUncorRel; } //!< Get covariance matrix of relative uncorrelated and stat. uncertainties

   // const TMatrixD& GetInverseCovariance() const {return CheckReturnInverse(fInvCov,fCov); } //!< Get invers of covariance matrix
   // const TMatrixD& GetInverseCovarianceRel() const {return CheckReturnInverse(fInvCov,fCov); }; //!< Get inverse matrix
   // const TMatrixD& GetInverseCovarianceStat() const {return CheckReturnInverse(fInvCov,fCov); }; //!< Get inverse matrix
   // const TMatrixD& GetInverseCovarianceStatRel() const {return CheckReturnInverse(fInvCov,fCov); }; //!< Get inverse matrix
   // const TMatrixD& GetInverseCovarianceStatUncorr() const {return CheckReturnInverse(fInvCov,fCov); }; //!< Get inverse matrix
   // const TMatrixD& GetInverseCovarianceStatUncorRel() const {return CheckReturnInverse(fInvCov,fCov); }; //!< Get inverse matrix

private:
   std::vector<std::string> fDaRequirements; //!< List of all AParm's which this function depends on
   std::vector<AData*> fChildren;  //! pointers to AData objects of datasets
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies
   
   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);};
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);};

private:
   void CalculateSuperErrors(); //!< calcuate errors of the 'super'-array
   // const TMatrixD& CheckReturnInverse(TMatrixD& Inv, const TMatrixDSym& Mat) const;

   // // convenient matrices
   // TMatrixDSym fCov;
   // TMatrixDSym fCovRel;
   // TMatrixDSym fStat;
   // TMatrixDSym fStatRel;
   // TMatrixDSym fStatUncor;
   // TMatrixDSym fStatUncorRel;
   // // inverse matrices
   // mutable TMatrixD fInvCov;
   // mutable TMatrixD fInvCovRel;
   // mutable TMatrixD fInvStat;
   // mutable TMatrixD fInvStatRel;
   // mutable TMatrixD fInvStatUncor;
   // mutable TMatrixD fInvStatUncorRel;

};


#endif
