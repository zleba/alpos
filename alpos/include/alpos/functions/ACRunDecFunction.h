//
#ifndef Alpos_ACRunDecFunction
#define Alpos_ACRunDecFunction

/* 
 An simple example for a AParmFunction

 */


#include "alpos/ATheory.h"
#include "fastnlotk/CRunDec.h"

//class ACRunDecFunction : public AParmFuncBase<double> {
class ACRunDecFunction : public AFuncD {

public:
   ACRunDecFunction(const std::string& name);
   virtual ~ACRunDecFunction();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual std::vector<double> GetQuick(int n,...); //< The possibilty to implement a quick access without changing of any parameters or recent (member-)values
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   CRunDec* fCRunDec = NULL;

protected:
   

};


#endif
