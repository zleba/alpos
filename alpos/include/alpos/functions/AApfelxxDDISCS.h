// DB 01.2019
#ifndef Alpos_AApfelxxDDISCS
#define Alpos_AApfelxxDDISCS

/* 
    Calculation of NC and CC cross sections
    using APFEL

 */

#include "alpos/ATheory.h"
#include <apfel/structurefunctionbuilder.h>
#include <vector>

class AApfelxxDDISCS : public AParmFuncBase<double> {

public:
   AApfelxxDDISCS(const std::string& name);
   virtual ~AApfelxxDDISCS();

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
   std::vector<double> q2,y,x;
   double charge, polty;
   bool IsRedCS, IsNC;

   // caches for faster access
   AFuncD* fPDF = NULL;
   AFuncD* fAs  = NULL;

   std::unique_ptr<const apfel::Grid> fGrid;
   std::function<apfel::StructureFunctionObjects(const double&, const std::vector<double>&)> fF2Obj;
   std::function<apfel::StructureFunctionObjects(const double&, const std::vector<double>&)> fFLObj;
   std::function<apfel::StructureFunctionObjects(const double&, const std::vector<double>&)> fF3Obj;

};


#endif
