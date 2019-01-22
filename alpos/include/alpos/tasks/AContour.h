// DB. 01/2015
#ifndef Alpos_AContour
#define Alpos_AContour

/** 
 * class AContour
 * class AContourResult
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
// #include "alpos/AContourResult.h"


// ____________________________________________________________________________ //
class AContourResult : public ATaskResult {
public:
   AContourResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~AContourResult() {};
protected:

};


// ____________________________________________________________________________ //
class AContour : public ATask {
public:
   static const std::string fTaskType;

public:
   // AContour();
   AContour(const std::string& aname); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~AContour();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   AContourResult* fResult;

};


#endif
