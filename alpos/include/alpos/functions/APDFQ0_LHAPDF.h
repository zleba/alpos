//
#ifndef Alpos_APDFQ0_LHAPDF
#define Alpos_APDFQ0_LHAPDF

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

//class APDFQ0_LHAPDF : public AParmFuncBase<double> {
class APDFQ0_LHAPDF : public AFuncD {

public:
   APDFQ0_LHAPDF(const std::string& name);
   virtual ~APDFQ0_LHAPDF();
   virtual bool Init();       //< Initialize the function
   virtual bool Update();

   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //!< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //!< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //!< Return function name (type)
   static const std::string fFunctionName; //!< The function's name
   virtual std::vector<double> GetQuick(int n,...); //!< The possibilty to implement a quick access without changing of any parameters or recent (member-)values
   virtual std::vector<double> GetQuick(const std::vector<double>&); //!< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values
   const std::vector<double>& GetDef() const { return fdef;} //!< get PDF linear combination
   const std::vector<double>& GetDefQcdnum() const { return fdefQcdnum;} //!< get PDF linear combination

private:
   static const std::vector<std::string> fRequirements; //!< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //!< List of Parm's which have changed, but this function does not notify further dependencies

protected:
   std::string fLico;
   std::vector<double> fdef;
   std::vector<double> fdefQcdnum;
};


#endif