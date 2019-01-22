//
#ifndef Alpos_APDFQ0_QcdnumExample
#define Alpos_APDFQ0_QcdnumExample

/** 
 Provide the PDFs at the starting scale from
 an 'external' PDF: typically LHAPDF

 Inputs are:
   + PDF        PDF function, e.g. LHAPDF6. Must provide 'GetQuick(x,muf(),ipdf);
   + PDFLiCo    String to specify the licaer combination of the 'ipdf' flag
 Internal variables:
   + x          x
   + Q0         the starting scale
   + iPDF       flag for the pdf linear combination. -1 returns array 'def' for QCDNUM

*/


#include "alpos/ATheory.h"

//class APDFQ0_QcdnumExample : public AParmFuncBase<double> {
class APDFQ0_QcdnumExample : public AFuncD {

public:
   APDFQ0_QcdnumExample(const std::string& name);
   virtual ~APDFQ0_QcdnumExample();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   virtual std::vector<double> GetQuick(int n,...); //< The possibilty to implement a quick access without changing of any parameters or recent (member-)values
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

protected:

};


#endif
