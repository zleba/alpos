// DB. 01/2015
#ifndef Alpos_ATaskResult
#define Alpos_ATaskResult

/** 
 * purely virtual base class ATaskResult
 * 
 * Base class for all results from ATasks
 * Each task must provide a result.
 * Some tasks may not, and take the default.
 * ATaskResult.
 *
 * ATaskResult can be identified uniquely 
 * by its name, and then typcasted to the
 * right class, in order to access members.
 * 
 * Each result holds the theory set (ATheorySet)
 * with all theory parameters from before and after
 * the execution of the task (initial and final
 * theory set).
 *
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "alpos/ATheory.h"
#include "alpos/Alpos.h"

// class ATaskResult;
// struct AResultCompare {
//    bool operator()(const ATaskResult& lhs, const ATaskResult& rhs) {
//       return lhs->GetResultName()<rhs->GetResultName();
//    }
// }
class Alpos;


class ATaskResult : public AlposObject {
   friend Alpos;

public:
   // ATaskResult();
   ATaskResult(const std::string& aname,const std::string& taskName) : AlposObject(), fResultName(aname), fTaskName(taskName) {};
   ~ATaskResult() {};

   const std::string& GetResultName() const { return fResultName;}; //!< get task name, i.e the name if this instance
   const std::string& GetTaskName() const { return fTaskName;}; //!< get task name, i.e the name if this instance

   const ATheorySet& GetInitialTheorySet() const { return fInitialTS;} //!< Get theory set, before the task was executed
   const ATheorySet& GetFinalTheorySet() const { return fFinalTS;} //!< Get theory set, after the task was executed
   
private:
   std::string fResultName;
   std::string fTaskName;
   ATheorySet fInitialTS;
   ATheorySet fFinalTS;

private:
   void SetInitialTheorySet(const ATheorySet&  TS ) {fInitialTS=TS;}
   void SetFinalTheorySet(const ATheorySet&  TS ) {fFinalTS=TS;}

};

#endif
