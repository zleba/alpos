// DB 15.01.2015
#ifndef Alpos_AfastNLOInterpolPDFas
#define Alpos_AfastNLOInterpolPDFas

/* 
 An implementation of fastNLO as an Alpos function.
 */

#include "alpos/ATheory.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include "TGraph.h"

class AfastNLOInterpolPDFas : public AParmFuncBase<double> {

public:
   AfastNLOInterpolPDFas(const std::string& name);
   virtual ~AfastNLOInterpolPDFas();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual bool ReInit();       //< Initialize the function

   void PrintInterpolatedPoints(unsigned int nPoints = 200);  //< Sample interpolation function at 'nPoints' points and print cross section values
   void PrintSupportPoints();  //< Print cross sections for each of the requested alpha_s values/PDFs
   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

protected:
   void CalcPDFUncertaintiesAsCovariance(double cl = 100 * 0.6826894 /*boost::math::erf(1 / sqrt(2))*/, bool alternative = false);
   void CalcHessianPDFUncertaintiesAsEigenvectors(double cl = 100 * 0.6826894 /*boost::math::erf(1 / sqrt(2)) */ );

   void SetOrder(fastNLOLHAPDF* fnlo);
   bool PrintValues();

   // members
   std::vector<TGraph> fGraphs;
   std::vector<TSpline*> fSplines;
   bool fHasSplines = false;

   std::string fNeedPDFUncertaintiesAs = "none"; //< whether to calculate PDF uncertainties or not
   bool fCalculatedPDFUncertainties = false; //< whether PDF uncertainties have already been calculated or not
   bool fRescalePDFUncertainties = true; //< whether to rescale PDF uncertainties on Update()

};


#endif
