// DB. 08/2015
#ifndef Alpos_AReplaceDataWithTheoryValues
#define Alpos_AReplaceDataWithTheoryValues

/** 
 * class AReplaceDataWithTheoryValues
 * class AReplaceDataWithTheoryValuesResult
 * 
 * A first implementation of a statistical analysis
 * 
 */

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/AChisq.h"
// #include "alpos/AReplaceDataWithTheoryValuesResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class AReplaceDataWithTheoryValuesResult : public ATaskResult {
public:
   AReplaceDataWithTheoryValuesResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AReplaceDataWithTheoryValuesResult() {};
protected:

};


// ____________________________________________________________________________ //
class AReplaceDataWithTheoryValues : public ATask {
public:
   static const std::string fTaskType;

public:
   // AReplaceDataWithTheoryValues();
   AReplaceDataWithTheoryValues(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AReplaceDataWithTheoryValues();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   // Alpos requirements
   AReplaceDataWithTheoryValuesResult* fResult;

   // other members
};


#endif
