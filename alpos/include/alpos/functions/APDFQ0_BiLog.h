//
#ifndef Alpos_APDFQ0_BiLog
#define Alpos_APDFQ0_BiLog


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


//class APDFQ0_BiLog : public AParmFuncBase<double> {
class APDFQ0_BiLog : public AFuncD {

public:
   APDFQ0_BiLog(const std::string& name);
   virtual ~APDFQ0_BiLog();
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


   double BiLog(double x, double A,double B,double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //!< Default parametrisation as used at HERA



   double SumRuleASpar(int n, double A, double B, double C, double D, double E);
   double SumRuleASpar(int n, const std::vector<double>& param) {return SumRuleASpar(n, param[0], param[1], param[2], param[3], param[4]);}
   double splogn(double x, double A, double B, double C, double D, double E);
   double splogn(double x,const std::vector<double>& param) {return splogn(x, param[0], param[1], param[2], param[3], param[4]);};



   /*
   double DefaultHERAParam(double x, const std::vector<double>& par){
      return DefaultHERAParam(x,par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8]);}; //!< Default parametrisation as used at HERA
   double GetIntegralDefaultHERAParam(double A, double B, double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //! Calculate integral of 'DefaultHERAParam'
   double GetIntegralDefaultHERAParam(const std::vector<double>& par){
      return GetIntegralDefaultHERAParam(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8]);}; //! Calculate integral of 'DefaultHERAParam'
   double GetIntegralXHERA(double A, double B, double C,double D=0,double E=0,double F=0,double AP=0,double BP=0,double CP=0); //! Calculate integral of x*'DefaultHERAParam'
   double GetIntegralXHERA(const std::vector<double>& par) {
      return GetIntegralXHERA(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8]);}//! Calculate integral of x*'DefaultHERAParam'

   double BiLog(double x, double A,double B,double C,double D,double E); //!< BiLog parameterization
   double BiLog(double x, const std::vector<double>& V) {return BiLog(x,V[0],V[1],V[2],V[3],V[4]);} //!< BiLog parameterization

   */
   // Sum rules...
   double CalcGluonASumRule(); // calc 'gA' parameter from sumrule
   double Get_UbarA();
   double Get_dvA();
   double Get_uvA();



   // arrays filed in 'update' for 'quick' access
   double ffs = 0;
   std::vector<double> fPar_g = std::vector<double>(5,0);
   std::vector<double> fPar_dv = std::vector<double>(5,0);
   std::vector<double> fPar_uv = std::vector<double>(5,0);
   std::vector<double> fPar_ub = std::vector<double>(5,0);
   std::vector<double> fPar_db = std::vector<double>(5,0);

};


#endif
