// DB. 01/2019
#ifndef Alpos_ASaveTheorySet
#define Alpos_ASaveTheorySet

/** 
 * class ASaveTheorySet
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
// class ASaveTheorySetResult : public ATaskResult {
// public:
//    ASaveTheorySetResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
//    ~ASaveTheorySetResult() {};
// protected:
// };


// ____________________________________________________________________________ //
class ASaveTheorySet : public ATask {
public:
   static const std::string& TaskType() { static const std::string TT = "SaveTheorySet"; return TT;}

public:
   // ASaveTheorySet();
   ASaveTheorySet(const std::string& aname) : ATask(aname) { fResult = new ATaskResult(aname,TaskType());  };
   virtual ~ASaveTheorySet(){};

   virtual bool Init(){return true;};
   virtual bool Execute() { TheoryHandler::Handler()->SaveCurrentTheorySet(GetTaskName());return true;};
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return TaskType();}; //!< Get task type, often similar to the class name

protected:
   ATaskResult* fResult;
};


#endif
