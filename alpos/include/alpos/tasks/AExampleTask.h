// DB. 01/2015
#ifndef Alpos_AExampleTask
#define Alpos_AExampleTask

/** 
 * class AExampleTask
 * class AExampleTaskResult
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
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
// #include "alpos/AExampleTaskResult.h"


// ____________________________________________________________________________ //
class AExampleTaskResult : public ATaskResult {
public:
   AExampleTaskResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AExampleTaskResult() {};
protected:

};


// ____________________________________________________________________________ //
class AExampleTask : public ATask {
public:
   static const std::string fTaskType;

public:
   // AExampleTask();
   AExampleTask(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AExampleTask(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AExampleTask();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   AExampleTaskResult* fResult;

};


#endif
