// DB 15.01.2015
#ifndef Alpos_AQcdnumAlphas
#define Alpos_AQcdnumAlphas

/* 
 An implementation of fastNLO as an Alpos function.
 */

#include "alpos/ATheory.h"
#include "alpos/functions/fastNLOAlpos.h"

class AQcdnumAlphas : public AParmFuncBase<double> {

public:
   AQcdnumAlphas(const std::string& name);
   virtual ~AQcdnumAlphas();

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

   virtual std::vector<double> GetQuick(int n,...) { std::cout<<"Notimplemented.8139"<<std::endl;exit(3);return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters

protected:
   

};


#endif
