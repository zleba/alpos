// DB. 01/2015
#ifndef Alpos_AChi2InterpolPDFas
#define Alpos_AChi2InterpolPDFas

/** 
 * class AChi2InterpolPDFas
 * class AChi2InterpolPDFasResult
 * 
 * Permorm a chi2 scan for all data(-subsets).
 * 
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

// ROOT
#include "Math/IParamFunction.h"

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
// #include "alpos/AChi2InterpolPDFasResult.h"

// ____________________________________________________________________________ //
class AChi2InterpolPDFasResult : public ATaskResult {
public:
   AChi2InterpolPDFasResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AChi2InterpolPDFasResult() {};

   double fChi2Minimum;
   double fParAtChi2Minimum;
   double fParErrorSymmetric;
   double fParErrorUp;
   double fParErrorDn;
protected:

};


// ____________________________________________________________________________ //
class AChi2InterpolPDFas : public ATask {
public:
   static const std::string fTaskType;

public:
   // AChi2InterpolPDFas();
   AChi2InterpolPDFas(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AChi2InterpolPDFas(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AChi2InterpolPDFas();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

   double Chi2Interpolation(double alpha_s) const;
   double Chi2Extrapolation(double alpha_s, int reference_pdf_member) const;

   std::string fChisqdef;
   std::string fPdfFunction;

   std::vector<double> fAsmzValues;
   std::vector<std::pair<double,double>> fAsmzChi2Pairs;
   std::vector<std::string> fPdfSets;
   std::vector<int> fPdfMembers;

   int GetNPoints() const {return fAsmzValues.size();};

protected:

   //std::string fTaskType;
   AChi2InterpolPDFasResult* fResult;


private:

   //void __minuitFCN(int& nPar, double* gradient, double& f, double par[], int& flags);
   //void minuitFCN(int& nPar, double* gradient, double& f, double par[], int& flags);
   //double __minuitFunction(const double *param);

   AFuncD* fAsRun;

};


class AChi2InterpolFCN : public ROOT::Math::IMultiGenFunction {
public:
   AChi2InterpolFCN(const AChi2InterpolPDFas* myChi2InterpolTask ) : fChi2InterpolTask(myChi2InterpolTask) {};
   virtual ~AChi2InterpolFCN() {};
   virtual unsigned int NDim() const { return 1; }; // one-dimensional fit
   virtual double DoEval(const double *params) const;
   ROOT::Math::IMultiGenFunction* Clone() const;

protected:
   const AChi2InterpolPDFas* fChi2InterpolTask;
};


#endif
