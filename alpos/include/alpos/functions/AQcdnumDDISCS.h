//
#ifndef Alpos_AQcdnumDDISCS
#define Alpos_AQcdnumDDISCS

/** 
 *
 * An simple example for a AParmFunction
 *
 */

#include "alpos/ATheory.h"

//class AQcdnumDDISCS : public AParmFuncBase<double> {
class AQcdnumDDISCS : public AFuncD {

public:
   AQcdnumDDISCS(const std::string& name);
   virtual ~AQcdnumDDISCS();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   //double rfluxRaw(double x_pom, double a0, double ap, double b0, double tcut) const;
   //double rflux(double x_pom, double a0, double ap, double b0) const;
   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

protected:
   

};


#endif
