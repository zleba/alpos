//
#ifndef Alpos_ADPDF
#define Alpos_ADPDF

/** 
 Interface to the Apfel++ PDF evolution.
 */


#include "alpos/ATheory.h"
#include <memory>

class ADPDF : public AParmFuncBase<double> {

public:
   ADPDF(const std::string& name);
   virtual ~ADPDF();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual std::vector<double> GetQuick(int n,...); //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

   //double rfluxRaw(double x_pom, double a0, double ap, double b0, double tcut);
   double rflux(double xpom, double tcut, double a0, double ap, double b0, double xPomNorm);
   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies
   double fPomNorm, fRegNorm;

protected:
};


#endif
