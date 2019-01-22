// DB. 08/2015
#ifndef Alpos_APrintDataTheory
#define Alpos_APrintDataTheory

/** 
 * class APrintDataTheory
 * class APrintDataTheoryResult
 * 
 * A first implementation of a statistical analysis
 * 
 */

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/AChisq.h"
// #include "alpos/APrintDataTheoryResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class APrintDataTheoryResult : public ATaskResult {
public:
   APrintDataTheoryResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~APrintDataTheoryResult() {};
protected:

};


// ____________________________________________________________________________ //
class APrintDataTheory : public ATask {
public:
   static const std::string fTaskType;

public:
   // APrintDataTheory();
   APrintDataTheory(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //APrintDataTheory(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~APrintDataTheory();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   // Alpos requirements
   APrintDataTheoryResult* fResult;

   // other members
   std::vector<std::string> fColumnNames;
   int fColumnWidth;

};


#endif
