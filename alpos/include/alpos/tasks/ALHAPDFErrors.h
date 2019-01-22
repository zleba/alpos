//
// Created by DS on 4/18/16.
//

#ifndef ALPOS_ALHAPDFERRORS_H
#define ALPOS_ALHAPDFERRORS_H


/**
 * class ALHAPDFErrors
 * class ALHAPDFErrorsResult
 *
 * A first implementation of a fitter class
 * Perform a fit of the theory to data.
 *
 *
 */

#include <set>

#include "alpos/AlposObject.h"
//#include "alpos/functions/ALhapdf6.h"
#include "alpos/ATask.h"


// ____________________________________________________________________________ //
class ALHAPDFErrorsResult : public ATaskResult {
public:
   ALHAPDFErrorsResult(const std::string& aname, const std::string& taskName) : ATaskResult(aname, taskName) { };

   ~ALHAPDFErrorsResult() { };
protected:

};


// ____________________________________________________________________________ //
class ALHAPDFErrors : public ATask {
public:
   static const std::string fTaskType;

public:
   // ALHAPDFErrors();
   ALHAPDFErrors(
         const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //ALHAPDFErrors(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~ALHAPDFErrors();

   virtual bool Init();

   virtual bool Execute();

   virtual ATaskResult* GetResult() { return fResult; }; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType; }; //!< Get task type, often similar to the class name

protected:

//   //std::vector<double> UncertaintiesForEigenvector(int eigenvector_symm_idx);
//   std::vector<double> Uncertainties(double cl = 100 * boost::math::erf(1 / sqrt(2)), bool alternative = false);
//   std::vector<double> CorrelationMatrix();
//   std::vector<double> UncertaintiesForEigenvector(int eigenvector_up_idx, int eigenvector_down_idx, double cl = 100 * boost::math::erf(1 / sqrt(2)), bool alternative = false);
//   void CovarianceMatrix(const AFuncD& theory_function, double cl = 100 * boost::math::erf(1 / sqrt(2)), bool alternative = false);

   std::string fErrorPrefix;
   std::string fLhapdfFunctionName;
   std::string fCalculateAs;
   std::string fOverridePDFSet;

   ALHAPDFErrorsResult* fResult;

private:

   bool CalcAndAddAsCovariance(std::vector<AFuncD*> th_funcs);

   bool CalcHessianAndAddAsEigenvectors(std::vector<AFuncD*> th_funcs);
   bool CalcSymHessianAndAddAsEigenvectors(std::vector<AFuncD*> th_funcs);

   std::vector<std::vector<double>> CalcObservableArray(AFuncD* th_func);
   std::vector<std::vector<std::vector<double>>> CalcObservableArrays(std::vector<AFuncD*> th_funcs);

};

#endif //ALPOS_ALHAPDFERRORS_H
