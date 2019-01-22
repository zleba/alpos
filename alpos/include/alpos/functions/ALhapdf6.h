//
#ifndef Alpos_ALhapdf6
#define Alpos_ALhapdf6

/* 
 An simple example for a AParmFunction

 */


#include "alpos/ATheory.h"
#include <LHAPDF/LHAPDF.h>

class ALhapdf6 : public AParmFuncBase<double> {

public:
   ALhapdf6(const std::string& name);
   virtual ~ALhapdf6();
   virtual bool Init();       //!< Initialize the function
   virtual bool Update();     //!< Calculate fValues

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name
   virtual std::vector<double> GetQuick(int n,...); //!< The possibilty to implement a quick access without changing of any parameters
   virtual std::vector<double> GetQuick(const std::vector<double>&); //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values
   LHAPDF::PDFSet* GetPDFSet() {return fPDFSet;} //!< return LHAPDF::PDFSet
   LHAPDF::PDF* GetPDF() {return fPDF;} //!< return LHAPDF::PDF
   
private:
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

   // LHAPDF members
   LHAPDF::PDFSet* fPDFSet = NULL;
   LHAPDF::PDF* fPDF = NULL;

protected:
   

};


#endif
