// DB 15.01.2015
#ifndef Alpos_ASuperTheory
#define Alpos_ASuperTheory

/* 
 An implementation of the most trivial function:
 It consists only of one value 'theConst'.
 This function may be used for data, which directly
 measure a constant.

 */


#include "alpos/ATheory.h"
#include <vector>

class ASuperTheory : public AFuncD /*AParmFuncBase<double>*/ {

public:
   ASuperTheory(const std::string& name);
   virtual ~ASuperTheory();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   
   // --- SuperTheory functions:
   void AddRequirement(const std::string& req) { fRequirements.push_back(req);} //!< Add one requirement (Requirements still have to be registered) 
   void SetRequirements(const std::vector<std::string>& reqs) {fRequirements=reqs;} //!< Set the requirements (Requirements still have to be registered) 

   // access to child theory objects
   const std::vector<AParmFuncBase<double>*> GetChildren() const { return fChildren; }

private:
   std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   std::vector<AParmFuncBase<double>*> fChildren;  //! pointers to child theory objects
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   void CalculateSuperErrors(); //!< calcuate errors of the 'super'-array

protected:
   

};


#endif
