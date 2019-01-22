//
// Created by DS on 4/22/16.
//

#ifndef ALPOS_APARMNAMED_H
#define ALPOS_APARMNAMED_H

#include "alpos/AlposObject.h"
#include <set>
#include <vector>

/**********************************************************

  AParmNamed: base class for alpos parameters

  all Alpos parameters are derived from this
  base class.

  the parameter observer pattern is implemented on this
  level, meaning that this class contains all methods necessary
  for updating dependent parameters when they change.

***********************************************************/

class AThCont;

class AParmNamed : public AlposObject {
public:
   AParmNamed(const std::string& name);
   ~AParmNamed();

   const std::string& GetAlposName() const { return fName; };

   virtual bool Update() { return true; };

   void SetIsOutdated(AParmNamed* notifier = nullptr, bool IsOutdated = true);

   bool GetIsUpToDate() const { return fIsUpToDate; }

   void NotifyDependencies(AParmNamed* notifier);

   // todo:   virtual void RegisterDependency(std::string name,AParmFunction<T>* pFunc) = 0;
   //virtual void RegisterDependency(const std::string& name,AParmBase<T>* pFunc) = 0;
   //virtual void AddDependency(const std::string& name) {fDependencies.insert(name);};
   void RegisterDependency(const std::string& name, AParmNamed* pFunc);

   std::map<std::string, AParmNamed*> GetAllDependencies() const { return fDependencies; }

   virtual void GetContent(std::vector<std::string>& val, std::vector<std::string>& err, bool& konst) {
      warn["GetContent"] << "GetContent should not be called on AParmNamed, but was called for (" << GetAlposName() <<
      ")." << std::endl;
   }  //!< Get the values, errors and const-flag as a strings

   virtual void GetContent(AThCont& cont) {
      warn["GetContent"] << "GetContent should not be called on AParmNamed, but was called for (" << GetAlposName() <<
      ")." << std::endl;
   }

   virtual void SetContent(const std::vector<std::string>& val, const std::vector<std::string>& err, bool konst) {
      warn["SetContent"] << "SetContent should not be called on AParmNamed, but was called for (" << GetAlposName() <<
      ")." << std::endl;
   }

   virtual void SetContent(const AThCont& cont) {
      warn["SetContent"] << "SetContent should not be called on AParmNamed, but was called for (" << GetAlposName() <<
      ")." << std::endl;
   }

private:
   std::string fName; //! or in alposObject
   bool fIsUpToDate = false;
   //std::set<std::string> fDependencies; // would have to call TheoryHandler to resolve name
   //std::set<AParmFunctionD*> fDependencies; // could directly call SetIsOutdate().
   //std::set<std::pair<std::string,AParmFunction<T>*> > fDependencies; // we can decide later which is the best call.
   //std::map<std::string,AParmFunction<T>*> fDependencies
   //std::set<std::string> fDependencies;

   std::map<std::string, AParmNamed*> fDependencies; //!< List of parm's which are going to be notified, if parm is being updated.
   //std::set<std::string> fDependencies2;
   std::set<AParmNamed*> fNotifiers = std::set<AParmNamed*>(
         {nullptr}); //!< pointers to parm which notified about an updated (only used reasonably in AParmFunc)

protected:
   void SetName(std::string name) { fName = name; }

   void ClearNotifiers() { fNotifiers.clear(); }

   const std::set<AParmNamed*>& GetCurrentNotifiers() const { return fNotifiers; }
};

#endif //ALPOS_APARMNAMED_H
