// DB. 01/2015
#ifndef Alpos_AStatAnalysis
#define Alpos_AStatAnalysis

/** 
 * class AStatAnalysis
 * class AStatAnalysisResult
 * 
 * A first implementation of a statistical analysis
 * 
 */

#include <string>
#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
#include "alpos/AChisq.h"
// #include "alpos/AStatAnalysisResult.h"

#include "Fit/Fitter.h"


// ____________________________________________________________________________ //
class AStatAnalysisResult : public ATaskResult {
public:
   AStatAnalysisResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AStatAnalysisResult() {};
protected:

};


// ____________________________________________________________________________ //
class AStatAnalysis : public ATask {
public:
   static const std::string fTaskType;

public:
   // AStatAnalysis();
   AStatAnalysis(const std::string& aname /*, const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AStatAnalysis(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AStatAnalysis();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name
   void Print(AData* dat, AFuncD* th);
   void Print(const std::string& dat, const std::string& th);

protected:

   // Alpos requirements
   AStatAnalysisResult* fResult;

   // other members
   AChisqBase* fChisq = NULL;
   //std::vector<double> fFitPar;

};


#endif
