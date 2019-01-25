//
#ifndef Alpos_AApfelxxPDF
#define Alpos_AApfelxxPDF

/** 
 Interface to the Apfel++ PDF evolution.
 */


#include "alpos/ATheory.h"
#include <apfel/grid.h>
#include <apfel/dglap.h>
#include <apfel/tabulateobject.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tools.h>
#include <memory>
#include <TMatrixD.h>

class AApfelxxPDF : public AParmFuncBase<double> {

public:
   AApfelxxPDF(const std::string& name);
   virtual ~AApfelxxPDF();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual std::vector<double> GetQuick(int n,...); //< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values
   
private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

protected:
   std::unique_ptr<const apfel::Grid> fGrid;
   std::map<int, apfel::DglapObjects> fDglapObj;
   std::unique_ptr<apfel::Dglap<apfel::Distribution>> fDglap;
   std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> fTabulatedPDFs;
   AFuncD* fAs  = NULL;
   AFuncD* fPDF0  = NULL;
   TMatrixD fPdf0ToApfl;
};


#endif
