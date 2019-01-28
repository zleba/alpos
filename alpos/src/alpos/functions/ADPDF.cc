
#include "alpos/functions/ADPDF.h"
#include "apfel/dglapbuilder.h"

#include <iostream>


//! 
//! A small function, to put together
//! different
//! 

using namespace std;

const std::vector<std::string> ADPDF::fRequirements = {  
   "xpom", // xpom (dummy)
   "zpom", // zpom (dummy)
   "muf", // mu_f (dummy)
   "pom1", // pomeron 1, this must be a 'PDF'
   "pom2", // pomeron 2, this must be a 'PDF'
   "reg1", // reggeon 1, this must be a 'PDF'
   "xPomFluxNorm", "tcut" ,// flux
   "Flux_pom1_a0","Flux_pom1_ap","Flux_pom1_b0",
   "Flux_pom2_a0","Flux_pom2_ap","Flux_pom2_b0",
   "Flux_reg1_a0","Flux_reg1_ap","Flux_reg1_b0","reg1_n",
}; //< List of all AParm's which this function depends on
const std::vector<std::string> ADPDF::fStopFurtherNotification = {"xpom","zpom","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ADPDF::fFunctionName = "DPDF"; //< The function's name


// __________________________________________________________________________________________ //
ADPDF::ADPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("ADPDF");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
ADPDF::~ADPDF() {
}


// ___________________________________________________________________________________________ //
bool ADPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;

   // CONST(pom1);
   // CONST(pom2);
   // CONST(reg1);
   
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> ADPDF::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(3,xp,Q);

   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   double xpom  = va_arg(ap, double);
   double zpom  = va_arg(ap, double);
   double muf = va_arg(ap, double);
   va_end(ap);
   return GetQuick({xpom,zpom,muf});
    
}


// ______________________________________________________________________________________ //
std::vector<double> ADPDF::GetQuick(const vector<double>& xpom_zpom_muf) {

   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> pdf(fValue.size(),0);
   if ( xpom_zpom_muf.size() != 3) {
      error["GetQuick"]<<"Quick acces is implemented for two parameter which are 'xpom','zpom' and 'muf'."<<endl;
      return pdf;
   }
   double xpom = xpom_zpom_muf[0];
   if ( xpom == 0 ) xpom = PAR(xpom);
   if ( xpom < 1.e-5 || xpom > 1 ) {
      error["Quick"]<<"Unreasonale xpom value:" << xpom<<endl; 
      cout<<"xpom in:  "<<xpom_zpom_muf[0]<<endl;
      cout<<"xpom par: "<<PAR(xpom)<<endl;
      exit(1);
   }
   double zpom = xpom_zpom_muf[1];
   double muf  = xpom_zpom_muf[2];

   // flux parameters
   double tcut = PAR(tcut);
   double xPomNorm= PAR(xPomFluxNorm);

   // pomeron
   pdf = QUICK(pom1,({zpom,muf}));
   //Pomeron flux
   double a0_IP = PAR(Flux_pom1_a0);
   double ap_IP = PAR(Flux_pom1_ap);
   double b0_IP = PAR(Flux_pom1_b0);
   double flxIP = rflux(xpom, tcut, a0_IP, ap_IP, b0_IP, xPomNorm);

   // cout<<"ADPDF! flxIP = " << flxIP <<endl;
   // cout<<"ADPDF! gluon   " <<pdf[6]<<endl;
   // cout<<"ADPDF! up      " <<pdf[7]<<endl;
   // cout<<"ADPDF! dn      " <<pdf[8]<<endl;

   for ( auto& p : pdf ) p*=flxIP;



   // pom2
   double a0_P2 = PAR(Flux_pom2_a0);
   if ( a0_P2!=0 ) {
      vector<double> pom2 = QUICK(pom2,({zpom,muf}));
      if ( pom2.size() != 13 ) {
	 error["Quick"]<<"pom2 is requested (Flux_pom2_a0!=0), but no reasonable PDF function is provided."<<endl;
	 exit(1);
      }
      double ap_P2 = PAR(Flux_pom2_ap);
      double b0_P2 = PAR(Flux_pom2_b0);
      double flxP2 = rflux(xpom, tcut, a0_P2, ap_P2, b0_P2, xPomNorm);
      for ( int i = 0 ; i < 13 ; i++ ) 
	 pdf[i] += pom2[i] * flxP2;
   }
   
   // reggeon
   vector<double> reg1 = QUICK(reg1,({zpom,muf}));
   //Reggeon flux 
   double n_IR = PAR(reg1_n);
   if ( n_IR!=0 ) {
      double a0_IR = PAR(Flux_reg1_a0);
      double ap_IR = PAR(Flux_reg1_ap);
      double b0_IR = PAR(Flux_reg1_b0);
      double flxIR = rflux(xpom, tcut, a0_IR, ap_IR, b0_IR, xPomNorm);
      if ( reg1.size()==13) {
	 for ( int i = 0 ; i < 13 ; i++ ) 
	    pdf[i] += reg1[i] * flxIR * n_IR;
      }
   }

   //cout<<"muf="<<muf<<"\txpom="<<xpom<<"\tzpom="<<zpom<<"\tpdf0="<<pdf[0]<<"\tpdf6="<<pdf[6]<<endl;
   
   return pdf;
  
}


// __________________________________________________________________________________________ //
bool ADPDF::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   // update
   SET(pom1.xp,PAR(zpom),0);
   SET(pom1.muf,PAR(muf),0);
   UPDATE(pom1);
   if ( PAR(Flux_pom2_a0) ) {
      SET(pom2.xp,PAR(zpom),0);
      SET(pom2.muf,PAR(muf),0);
      UPDATE(pom2);
   }
   if ( PAR(reg1_n) != 0 ) {
      SET(reg1.xp,PAR(zpom),0);
      SET(reg1.muf,PAR(muf),0);
      UPDATE(reg1);
   }

   fValue = GetQuick(vector<double>{PAR(xpom),PAR(zpom),PAR(muf)});

   cout<<this->GetAlposName();
   for ( auto p : fValue ) cout<<"\t"<<p;
   cout<<endl;

   return true;
}

//tcut is negative : tcut = -1
double ADPDF::rfluxRaw(double xpom, double a0, double ap, double b0, double tcut)
{
   const double mp = 0.93827231;

   //     calc min. kinematically  allowed t
   double tmin= -pow(mp*xpom,2)/(1.-xpom);

   //     c*xpom**(-(2apom-1))
   double fl =  exp((2.0*a0-1.)*log(1.0/xpom));
   double b=(b0+2.0*ap*log(1.0/xpom));

   //   at fixed t:  exp(Bt)
   //  fl = fl * exp(b*tcut);

   //   t-integrated: (1/B)*[exp(-B*tmax)-exp(-B*tmin)]
   fl = fl * (exp(tmin*b)-exp(tcut*b))/b;

   return fl;
}


double ADPDF::rflux(double xpom, double tcut, double a0, double ap, double b0, double xPomNorm)
{
   // double tcut = -1;
   // double xPomNorm = 0.003;
   const double dm =  rfluxRaw(xPomNorm, a0, ap, b0, tcut);
   double  norm=(1./(xPomNorm*dm)); //xpom * flux normalized to 1 at xpom = 0.003

   return  norm * rfluxRaw(xpom, a0, ap, b0, tcut);
}



