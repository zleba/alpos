// DB. 01/2015
#ifndef Alpos_AChi2Scan
#define Alpos_AChi2Scan

/** 
 * class AChi2Scan
 * class AChi2ScanResult
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
// #include "alpos/AChi2ScanResult.h"


// ____________________________________________________________________________ //
class AChi2ScanResult : public ATaskResult {
public:
   AChi2ScanResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AChi2ScanResult() {};
protected:

};


// ____________________________________________________________________________ //
class AChi2Scan : public ATask {
public:
   static const std::string fTaskType;

public:
   // AChi2Scan();
   AChi2Scan(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //AChi2Scan(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AChi2Scan();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   AChi2ScanResult* fResult;

};


#endif
