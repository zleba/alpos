// DB. 01/2015
#ifndef Alpos_AApcalc
#define Alpos_AApcalc

/** 
 * class AApcalc
 * class AApcalcResult
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
// #include "alpos/AApcalcResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class AApcalcResult : public ATaskResult {
public:
   AApcalcResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AApcalcResult() {};
protected:

};


// ____________________________________________________________________________ //
class AApcalc : public ATask {
public:
   static const std::string fTaskType;

public:
   // AApcalc();
   AApcalc(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AApcalc(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AApcalc();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   // Alpos requirements
   AApcalcResult* fResult;

   // other members

};


#endif
