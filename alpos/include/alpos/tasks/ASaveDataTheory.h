// DB. 01/2015
#ifndef Alpos_SaveDataTheory
#define Alpos_SaveDataTheory

/** 
 * class ASaveDataTheory
 * class ASaveDataTheoryResult
 * 
 * A first implementation of a fitter class
 * Perform a fit of the theory to data.
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
// #include "alpos/ASaveDataTheoryResult.h"


// ____________________________________________________________________________ //
class ASaveDataTheoryResult : public ATaskResult {
public:
   ASaveDataTheoryResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~ASaveDataTheoryResult() {};
protected:

};


// ____________________________________________________________________________ //
class ASaveDataTheory : public ATask {
public:
   static const std::string fTaskType;

public:
   // ASaveDataTheory();
   ASaveDataTheory(const std::string& aname/* , const std::string& rsnmsp, const std::map<std::string,ATaskResult> const *previousResults*/); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //ASaveDataTheory(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~ASaveDataTheory();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   ASaveDataTheoryResult* fResult;

};


#endif
