// DB 15.01.2015
#ifndef Alpos_ASubsetFunction
#define Alpos_ASubsetFunction

/* 
 An implementation of the most trivial function:
 It consists only of one value 'theConst'.
 This function may be used for data, which directly
 measure a constant.

 */

#include "alpos/ATheory.h"

class ASubsetFunction : public AParmFuncBase<double> {

public:
   ASubsetFunction(const std::string& name);
   virtual ~ASubsetFunction();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const; //< Return the function name of the underlying function (type)
   std::string GetRealFunctionName() const { return fFunctionName;}; //< Return 'SubsetFunction'
   static const std::string fFunctionName; //< The function's name
   void InitErrors(); //!< Initialize fAllErrors
   void SetRequirementValidPoints(const std::string& req, const std::vector<bool>& valid);
   
private:
   std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   std::vector<bool> fPointValid;

protected:
   

};


#endif
