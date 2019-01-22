// DB. 01/2015
#ifndef Alpos_AApcFitter
#define Alpos_AApcFitter

/** 
 * class AApcFitter
 * class AApcFitterResult
 * 
 * An implementation of the constrained least square fit
 * Perform a fit of the theory to data using Apcalc.
 * 
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/tasks/AConstraint.h"



// ____________________________________________________________________________ //
class AApcFitterResult : public ATaskResult {
public:
   AApcFitterResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AApcFitterResult() {};
protected:

};


// ____________________________________________________________________________ //
class AApcFitter : public ATask {
public:
   static const std::string fTaskType;

public:
   // AApcFitter();
   AApcFitter(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AApcFitter(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AApcFitter();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   AApcFitterResult* fResult;
   AConstraint* fConstraint;
};


#endif
