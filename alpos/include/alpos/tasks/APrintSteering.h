// DB. 01/2015
#ifndef Alpos_APrintSteering
#define Alpos_APrintSteering

/** 
 * class APrintSteering
 * 
 * A first implementation of a fitter class
 * Perform a fit of the theory to data.
 * 
 * 
 */

#include "alpos/ATask.h"
#include "fastnlotk/read_steer.h"

// ____________________________________________________________________________ //
class APrintSteering : public ATask {
public:
   static const std::string& TaskType() { static const std::string TT = "PrintSteering"; return TT;}

public:
   // APrintSteering();
   APrintSteering(const std::string& aname) : ATask(aname) { fResult = new ATaskResult(aname,TaskType()); };
   virtual ~APrintSteering(){};

   virtual bool Init(){return true;};
   virtual bool Execute() { PRINTALL();return true;};
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return TaskType();}; //!< Get task type, often similar to the class name

protected:
   ATaskResult* fResult;
};


#endif
