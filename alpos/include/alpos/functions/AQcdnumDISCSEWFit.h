// DB 15.02.2016
#ifndef Alpos_AQcdnumDISCSEWFit
#define Alpos_AQcdnumDISCSEWFit

/* 
 An alternative implementation of fastNLO

 */

#include "alpos/ATheory.h"
#include "fastnlotk/fastNLOReader.h"

class AQcdnumDISCSEWFit : public AParmFuncBase<double> {

public:
   AQcdnumDISCSEWFit(const std::string& name);
   virtual ~AQcdnumDISCSEWFit();

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
   double Gets2wOS_from_s2wMSbar(double q2, double s2wMSbar);


};


#endif
