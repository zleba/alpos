// DB. 08/2015
#ifndef Alpos_APrintErrorSummary
#define Alpos_APrintErrorSummary

/** 
 * class APrintErrorSummary
 * class APrintErrorSummaryResult
 * 
 * A first implementation of a statistical analysis
 * 
 */

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/AChisq.h"
// #include "alpos/APrintErrorSummaryResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class APrintErrorSummaryResult : public ATaskResult {
public:
   APrintErrorSummaryResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~APrintErrorSummaryResult() {};
protected:

};


// ____________________________________________________________________________ //
class APrintErrorSummary : public ATask {
public:
   static const std::string fTaskType;

public:
   // APrintErrorSummary();
   APrintErrorSummary(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //APrintErrorSummary(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~APrintErrorSummary();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name
   void FillSummary(AData* dat, AFuncD* th);
   void FillSummary(const std::string& dat, const std::string& th);

protected:

   // Alpos requirements
   APrintErrorSummaryResult* fResult;

   // other members
   std::map<std::string,std::map<std::string,std::string> > fsummary;
   std::vector<std::string> fdatsets;
};


#endif
