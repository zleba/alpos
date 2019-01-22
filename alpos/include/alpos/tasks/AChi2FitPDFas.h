// DB. 01/2015
#ifndef Alpos_AChi2FitPDFas
#define Alpos_AChi2FitPDFas

/** 
 * class AChi2FitPDFas
 * class AChi2FitPDFasResult
 * 
 * Permorm a chi2 scan for all data(-subsets).
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
// #include "alpos/AChi2FitPDFasResult.h"


// ____________________________________________________________________________ //
class AChi2FitPDFasResult : public ATaskResult {
public:
   AChi2FitPDFasResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AChi2FitPDFasResult() {};

   double fChi2Minimum;
   double fParAtChi2Minimum;
   double fParErrorSymmetric;
   double fParErrorUp;
   double fParErrorDn;
protected:

};


// ____________________________________________________________________________ //
class AChi2FitPDFas : public ATask {
public:
   static const std::string fTaskType;

public:
   // AChi2FitPDFas();
   AChi2FitPDFas(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AChi2FitPDFas(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AChi2FitPDFas();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   AChi2FitPDFasResult* fResult;

};


#endif
