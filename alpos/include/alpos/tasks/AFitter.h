// DB. 01/2015
#ifndef Alpos_AFitter
#define Alpos_AFitter

/** 
 * class AFitter
 * class AFitterResult
 * 
 * A first implementation of a fitter class
 * Perform a fit of the theory to data.
 * 
 * 
 */

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/AChisq.h"
// #include "alpos/AFitterResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class AFitterResult : public ATaskResult {
public:
   AFitterResult(const std::string& aname, const std::string& taskName) : ATaskResult(aname, taskName) { };

   ~AFitterResult() { };
   double chisq;
   /// more parameters todo
protected:

};


// ____________________________________________________________________________ //
class AFitter : public ATask {
public:
   static const std::string fTaskType;

public:
   // AFitter();
   AFitter(const std::string& aname); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AFitter();

   void PrintNuisanceParameters() const;

   void FitDataTheoryPairs() const;

   virtual bool Init();

   virtual bool Execute();

   virtual ATaskResult* GetResult() { return fResult; }; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType; }; //!< Get task type, often similar to the class name
   ROOT::Fit::Fitter* GetFitter() { return fFitter; }

   //!< return Fit::Fitter object
   double GetChisq() const { return fFinalChisq; } //!< get final chisq
protected:

   // Alpos requirements
   AFitterResult* fResult;

   // other members
   ROOT::Fit::Fitter* fFitter = NULL;
   AChisqBase* fChisq = NULL;
   std::vector<double> fFitPar;

   double fInitialChisq;
   double fFinalChisq;

   // parameter "metadata"
   std::vector<std::string> fParNames;
   std::vector<bool> fParIsLimited;
   std::vector<bool> fParIsFixed;
   std::vector<double> fParStepSizes;
   std::vector<double> fParLowerLimits;
   std::vector<double> fParUpperLimits;
};


#endif
