// DB 15.01.2015
#ifndef Alpos_AQcdnumDISCS
#define Alpos_AQcdnumDISCS

/* 
 An alternative implementation of fastNLO

 */

#include "alpos/ATheory.h"
#include "fastnlotk/fastNLOReader.h"

class AQcdnumDISCS : public AParmFuncBase<double> {

public:
   AQcdnumDISCS(const std::string& name);
   virtual ~AQcdnumDISCS();

   // --- virtual functions from AParmFuncBase
   virtual bool Update();
   virtual bool Init();       //!< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name

   
private:
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   //v irtual functions from AParmFuncBase. Not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

protected:
   

};


#endif
