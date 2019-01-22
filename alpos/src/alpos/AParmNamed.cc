//
// Created by DS on 4/22/16.
//

#include "alpos/AParmNamed.h"
#include "alpos/ATheoryHandler.h"

//AParmNamed() : AParmNamed("UntitledParameter") {  };
AParmNamed::AParmNamed(const std::string& name) : AlposObject("AParmNamed"), fName(name) {
   TheoryHandler::Handler()->RegisterParameter(this);
}

AParmNamed::~AParmNamed() {
   //TheoryHandler::Handler()->EraseParameter(GetAlposName());
}

void AParmNamed::SetIsOutdated(AParmNamed* notifier, bool IsOutdated) {
   //! Set the AParmNamed to be 'outdated' (i.e. not up-to-date)
   //! When accessing the next time the 'values', then parameter will
   //! call an Update().
   //! 'notifier' specifies the pointer to the parameter, which has called SetIsOutdated
   //! a NULL ptr will result in a 'full' update, while others may
   //! be used for a more detailed Update().
   fIsUpToDate = !IsOutdated;
   if (IsOutdated) {
      NotifyDependencies(this);
      if (notifier != this) {
         NotifyDependencies(notifier);
         fNotifiers.insert(notifier);
      }
   }
   else ClearNotifiers();
}

void AParmNamed::NotifyDependencies(AParmNamed* notifier) {
   for (auto ir : fDependencies) {
      TheoryHandler::Handler()->GetParameter(ir.first)->SetIsOutdated(notifier);
   }
   // for ( auto ir : fDependencies2 ) {
   // 	 TheoryHandler::Handler()->GetParameter(ir)->SetIsOutdated(notifier);
   // }
}

// todo:   virtual void RegisterDependency(std::string name,AParmFunction<T>* pFunc) = 0;
//virtual void RegisterDependency(const std::string& name,AParmBase<T>* pFunc) = 0;
//virtual void AddDependency(const std::string& name) {fDependencies.insert(name);};
void AParmNamed::RegisterDependency(const std::string& name, AParmNamed* pFunc) {
   debug["RegisterDependency"] << "name=" << name << " at " << GetAlposName() << std::endl;
   //fDependencies2.insert(name);
   fDependencies[name] = pFunc;
}