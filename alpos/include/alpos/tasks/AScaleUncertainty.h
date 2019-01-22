// DB. 01/2015
#ifndef Alpos_AScaleUncertainty
#define Alpos_AScaleUncertainty

/** 
 * class AScaleUncertainty
 * class AScaleUncertaintyResult
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
// #include "alpos/AScaleUncertaintyResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class AScaleUncertaintyResult : public ATaskResult {
public:
   AScaleUncertaintyResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AScaleUncertaintyResult() {};
   double chisq;
   /// more parameters todo
protected:

};


// ____________________________________________________________________________ //
class AScaleUncertainty : public ATask {
public:
   static const std::string fTaskType;

public:
   AScaleUncertainty(const std::string& aname ); //!< constructor
   virtual ~AScaleUncertainty();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

   // Alpos requirements
   AScaleUncertaintyResult* fResult;

};


#endif
