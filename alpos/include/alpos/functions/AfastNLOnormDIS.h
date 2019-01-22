// lc DS 11.12.2015
#ifndef Alpos_AfastNLOnormDIS
#define Alpos_AfastNLOnormDIS

/* 
 * An Alpos implementation of fastNLO for calculating
 * jet cross sections normalized to deep inelastic
 * scattering cross sections.
 */
#include <vector>
#include <map>
#include <utility>

#include "alpos/ATheory.h"
#include "alpos/functions/fastNLOAlpos.h"

class AfastNLOnormDIS : public AParmFuncBase<double> {

public:
   AfastNLOnormDIS(const std::string& name);
   virtual ~AfastNLOnormDIS();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   const std::vector<bool>& GetBinmap() const { return fBinmap;}; //!< get bin map of valid points obtained from datacard
   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   double EvolveAlphas(double Q) const ;
   virtual bool InitPDF();
   std::vector<double> GetXFX(double xp, double muf) const;

   // members
   //fastNLOAlpos* fnlo = NULL;
   std::vector<fastNLOAlpos*> fnlos;  // allow more than one instance of fastNLO
   //int nTables;  // number of fastNLO instances
   
   // not used here (therefore made private)
   virtual std::vector<double> GetQuick(int n,...) { return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&) { return std::vector<double>(1);}; //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values
   
   // store NC DIS cross sections as a map from <q2min,q2max> pairs
   std::map < std::pair<double,double> , double > q2LoUpCSmap;

protected:
   void SetOrder();
   double CalcNCDISCrossSection(double q2min, double q2max);
   std::vector<bool> fBinmap;
   std::vector <std::pair<double,double> > fQ2BinsLoUps;

};


#endif
