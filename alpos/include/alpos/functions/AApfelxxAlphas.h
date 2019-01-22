// DB 07.2017
#ifndef Alpos_AApfelxxAlphas
#define Alpos_AApfelxxAlphas

/** 
 An implementation of the Apfel++ alpha_s evolution
 */

#include "alpos/ATheory.h"
#include "apfel/alphaqcd.h"
#include <memory>

class AApfelxxAlphas : public AParmFuncBase<double> {

public:
   AApfelxxAlphas(const std::string& name);
   virtual ~AApfelxxAlphas();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   virtual std::vector<double> GetQuick(int n,...); //< The possibilty to implement a quick access without changing of any parameters

protected:

   std::unique_ptr<apfel::AlphaQCD> fAs;
   AFuncD* fAsAlpos  = NULL;
};


#endif
