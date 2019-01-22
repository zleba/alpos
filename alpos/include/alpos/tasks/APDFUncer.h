// DB, K. Bjoerke, 08/2015
#ifndef Alpos_APDFUncer
#define Alpos_APDFUncer

/** 
 * class APDFUncer
 * class APDFUncerResult
 * 
 * Obtain PDF uncertainty (e.g. from an LHAPDF-type PDF)
 * for fit parameters.
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
// #include "alpos/APDFUncerResult.h"


// ____________________________________________________________________________ //
class APDFUncerResult : public ATaskResult {
public:
   APDFUncerResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~APDFUncerResult() {};
protected:

};


// ____________________________________________________________________________ //
class APDFUncer : public ATask {
public:
   static const std::string fTaskType;

public:
   APDFUncer(const std::string& aname);
   virtual ~APDFUncer();
   virtual bool Init();
   virtual bool Execute();
   virtual const std::string& GetTaskType() const {return fTaskType;};
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results

protected:
   APDFUncerResult* fResult;   

};


#endif
