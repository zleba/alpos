// DB 15.01.2015
#ifndef Alpos_AfastNLORatio
#define Alpos_AfastNLORatio

/*!
 AfastNLORatio

 This function computes the ratio of two fastNLO tables.
 */

#include "alpos/ATheory.h"
#include "alpos/functions/fastNLOAlpos.h"
#include "fastnlotk/fastNLOLHAPDF.h"

class AfastNLORatio : public AParmFuncBase<double> {

public:
   AfastNLORatio(const std::string& name);

   virtual ~AfastNLORatio();

   virtual bool Update();

   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements; }; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification; }; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName; }; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   const std::vector<bool>& GetBinmap() const { return fBinmap; }; //!< get bin map of valid points obtained from datacard

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   bool CalcCrossSections();

   // members
   std::vector<fastNLOAlpos*> fnlosNumerator;
   std::vector<fastNLOAlpos*> fnlosDenominator;

   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n, ...) {
      return std::vector<double>(1);
   }; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) {
      return std::vector<double>(1);
   }; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

protected:

   void SetOrder();

   std::vector<bool> fBinmap;

};


#endif
