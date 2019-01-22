// DB. 01/2015
#ifndef Alpos_AConstLQFitter
#define Alpos_AConstLQFitter

/** 
 * class AConstLQFitter
 * class AConstLQFitterResult
 * 
 * An implementation of the constrained least square fit
 * Perform a fit of the theory to data using Apccpp.
 * 
 */

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include <Apccpp.h>

#include "alpos/tasks/AConstraint.h"

// ____________________________________________________________________________ //
class AConstLQFitterResult : public ATaskResult {
public:
   AConstLQFitterResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AConstLQFitterResult() {};
protected:

};


// ____________________________________________________________________________ //
class AConstLQFitter : public ATask {
public:
   static const std::string fTaskType;

public:
   // AConstLQFitter();
   AConstLQFitter(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AConstLQFitter(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AConstLQFitter();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   // Alpos requirements
   AConstLQFitterResult* fResult;

   // other members
   Apccpp* fFitter = nullptr;
   AConstraint* fConstraint = nullptr;
   std::vector<double> fFitPar;

   double fInitialChisq;
   double fFinalChisq;
};


#endif
