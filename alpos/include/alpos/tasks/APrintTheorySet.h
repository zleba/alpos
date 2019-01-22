// DB. 01/2015
#ifndef Alpos_APrintTheorySet
#define Alpos_APrintTheorySet

/** 
 * class APrintTheorySet
 * 
 * A first implementation of a fitter class
 * Perform a fit of the theory to data.
 * 
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "alpos/ATask.h"
#include "alpos/ATheory.h"


// ____________________________________________________________________________ //
//! TaskResult not needed here
// class APrintTheorySetResult : public ATaskResult {
// public:
//    APrintTheorySetResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
//    ~APrintTheorySetResult() {};
// protected:
// };


// ____________________________________________________________________________ //
class APrintTheorySet : public ATask {
public:
   static const std::string& TaskType() { static const std::string TT = "PrintTheorySet"; return TT;}

public:
   // APrintTheorySet();
   APrintTheorySet(const std::string& aname) : ATask(aname) { fResult = new ATaskResult(aname,TaskType());  };
   virtual ~APrintTheorySet(){};

   virtual bool Init(){return true;};
   virtual bool Execute() { TheoryHandler::Handler()->PrintCurrentTheorySet();return true;};
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return TaskType();}; //!< Get task type, often similar to the class name

protected:
   ATaskResult* fResult;
};


#endif
