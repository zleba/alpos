//
#ifndef Alpos_APDFQ0_HERAStyleNoSum
#define Alpos_APDFQ0_HERAStyleNoSum

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
#include "TF1.h"
#include "TString.h"


//class APDFQ0_HERAStyleNoSum : public AParmFuncBase<double> {
class APDFQ0_HERAStyleNoSum : public AFuncD {

public:
   APDFQ0_HERAStyleNoSum(const std::string& name);
   virtual ~APDFQ0_HERAStyleNoSum();
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

   double DefaultHERAParam(double x, double A,double B,double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //!< Default parametrisation as used at HERA
   double GetIntegralDefaultHERAParam(double A, double B, double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //! Calculate integral of 'DefaultHERAParam'
   double GetIntegralXHERA(double A, double B, double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //! Calculate integral of x*'DefaultHERAParam'

   // Sum rules...
   double CalcGluonASumRule(); // calc 'gA' parameter from sumrule
   double Get_UbarA();
   double Get_dvA();
   double Get_uvA();
   double CalcIntegral(double alpha, double beta);

   TF1& GetTF1(double A, double B, double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0);
   TF1 fTF1Def3 = TF1("HERA 3 param", "[0]*pow(x,[1]) * pow(1-x,[2]) ", 1e-7, 1);
   TF1 fTF1Def6 = TF1("HERA 6 param", "[0]*pow(x,[1]) * pow(1-x,[2]) * (1 + [3]*x + [4]*pow(x,2)+ [5]*pow(x,3))", 1e-7, 1);
   TF1 fTF1Def9 = TF1("HERA 9 param", "[0]*pow(x,[1]) * pow(1-x,[2]) * (1 + [3]*x + [4]*pow(x,2)+ [5]*pow(x,3)) - [6]*pow(x,[7])*pow(1-x,[8])", 1e-7, 1);
};


#endif
