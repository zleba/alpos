// DB 15.01.2015
#ifndef Alpos_AfastNLOInterpolPDFasNormDIS
#define Alpos_AfastNLOInterpolPDFasNormDIS

/* 
 An implementation of fastNLO as an Alpos function.
 */

#include "alpos/ATheory.h"
#include "alpos/functions/AfastNLOInterpolPDFas.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include "TGraph.h"

class AfastNLOInterpolPDFasNormDIS : public AfastNLOInterpolPDFas {

public:
   AfastNLOInterpolPDFasNormDIS(const std::string& name);
   virtual ~AfastNLOInterpolPDFasNormDIS();

   //virtual bool Update();
   //virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name

   bool ReInit();       //< Initialize the function

   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   // store NC DIS cross sections as a map from <q2min,q2max> pairs
   std::map < std::pair<double,double> , double > q2LoUpCSmap;

protected:

   double CalcNCDISCrossSection(double q2min, double q2max);

};


#endif
