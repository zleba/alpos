#ifndef Alpos_ABelleTDCPV
#define Alpos_ABelleTDCPV

/* 
 An implementation for the Belle TDCPV fit
 */

#include "alpos/ATheory.h"
#include <TFile.h>
#include <TH2D.h>

class ABelleTDCPV : public AParmFuncBase<double> {

public:
   ABelleTDCPV(const std::string& name);
   virtual ~ABelleTDCPV();

   virtual bool Update();
   virtual bool Init();       //< Initialize the function
   virtual std::vector<std::string> GetRequirements() const { return fRequirements;}; //< List of all AParm's which this function depends on
   virtual std::vector<std::string> GetStopFurtherNotification() const { return fStopFurtherNotification;} ; //< List of Parm's which have changed, but this function does not notify further dependencies
   virtual std::string GetFunctionName() const { return fFunctionName;}; //< Return function name (type)
   static const std::string fFunctionName; //< The function's name
   
   virtual std::vector<double> GetQuick(const std::vector<double>&); //< Another possibilty to implement a quick access without changing of any parameters or recent (member-)values

private:
   static const std::vector<std::string> fRequirements; //< List of all AParm's which this function depends on
   static const std::vector<std::string> fStopFurtherNotification; //< List of Parm's which have changed, but this function does not notify further dependencies

   virtual std::vector<double> GetQuick(int n,...) { std::cout<<"Notimplemented.8139"<<std::endl;exit(3);return std::vector<double>(1);}; //< The possibilty to implement a quick access without changing of any parameters

protected:
   std::vector<TFile*> fMcFiles;
   //std::vector<TFile*> fDataFiles;
   double Pdt(double dt, double q, double tau, double dm, double A, double S);
   TH2D   MakedtddtTH2D(std::string name);
   int    GetClassId(double dtpar, int rqi, int ftag, int PDGtag);
   
   std::map<int,TH2D> fH2Dt_dtddt;
   std::map<int,TH2D> fH2MC_dtddt_p;
   std::map<int,TH2D> fH2MC_dtddt_m;
   

};


#endif
