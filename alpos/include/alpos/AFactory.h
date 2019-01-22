// DB 15.01.2015
#ifndef Alpos_AFactory
#define Alpos_AFactory

#include "alpos/ATheory.h"
#include "alpos/AData.h"
#include "alpos/ATask.h"
#include "alpos/AChisq.h"
#include <vector>
// ROOT
#include "Math/IParamFunction.h"

namespace AFactory {
   
   //! The task and theoryfunction factory
   AFuncD* FunctionFactory(const std::string& functype,const std::string& funcname); //! init a new function
   ATask* TaskFactory(const std::string& tname, const std::string& ttype); //!< instantiat a new task
   AChisqBase* ChisqFactory(const std::string& chisq, const std::vector<std::string>& par, AData* data, AFuncD* theo);

}

#endif
