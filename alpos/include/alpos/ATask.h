// DB. 01/2015
#ifndef Alpos_ATask
#define Alpos_ATask

/** 
 * purely virtual base class ATask
 * 
 * Base class for all tasks that Alpos can perform
 * Inherited classes must implement:
 *    Init(), Execute(), GetResult(), GetTaskName()
 * All tasks must further be implemented in Alpos.cc
 *  in Alpos::FindTask().
 * All tasks must further implement a method GetResult.
 * It is often convenient to inherit an suitable result
 * class from ATaskResult for a specific task.
 * A task is then exectued by calling Init() first,
 * then once Execute().
 * The 'result' object should be created in the heap
 * (i.e. with new) and is then owned by Alpos. Don't 
 * destroy it in the destrutor!
 * Tasks may be designed such, that they take results
 * from other tasks, which have been executed earlier.
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATaskResult.h"

class ATask : public AlposObject {
public:
   ATask(const std::string& aname ); //!< aname defines task name and also read_steer-namespace, i.e. where to look up the steering values
   virtual ~ATask();
   virtual bool Init() = 0;
   virtual bool Execute() = 0;
   virtual ATaskResult* GetResult() = 0; //!< Get the results
   virtual const std::string& GetTaskType() const = 0; //!< Get task type, often similar to the class name
   const std::string& GetTaskName() const { return fTaskName;}; //!< get task name, i.e the name if this instance
   const std::string& NS() const {return fTaskName;} //!< Get read_steer namespace
protected:
   std::string fTaskName;
};

#endif
