#if _CMAKE_FOUND_APFEL
#include "APFEL/APFEL.h"
#include "alpos/functions/AApfelDDISCS.h"
#include "alpos/functions/AApfel.h"
//#include "alpos/functions/AApfelInit.h"
#include <iostream>
#include "fastnlotk/read_steer.h"
#include <set>

using namespace std;

extern "C" {
   void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
   void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}
// double rfluxInt(double a0, double ap, double b0, double x_pom, double tAbsMin, double tAbsMax);
// double rflux(double a0, double ap, double b0, double x_pom, double tAbs);


// __________________________________________________________________________________________ //
const std::vector<std::string> AApfelDDISCS::fRequirements = {
   "ApfelInit",
   "e-charge",
   "e-polarity",
   "a0_IP", "ap_IP", "b0_IP",     //Pomeron flux       
   "a0_IR", "ap_IR", "b0_IR",     //Reggeon flux       
   "n_IR",                        //Reggeon suppression
   "ReggeonCS"                    // function to calculate reggeon cross section (use H1FitA for hard-coded ones)
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AApfelDDISCS::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelDDISCS::fFunctionName = "ApfelDDISCS"; //!< The function's name

// __________________________________________________________________________________________ //
AApfelDDISCS::AApfelDDISCS(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AApfelDDISCS");
}


// __________________________________________________________________________________________ //
AApfelDDISCS::~AApfelDDISCS() {
}


// ___________________________________________________________________________________________ //
bool AApfelDDISCS::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   AlposObject::debug["Init"]<<endl;

   // get points
   q2 = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   y = DOUBLE_COL_NS(Data,y,GetAlposName());
   xpom = DOUBLE_COL_NS(Data,xp,GetAlposName());
   beta = DOUBLE_COL_NS(Data,beta,GetAlposName());
   sigmaData = DOUBLE_COL_NS(Data,Sigma,GetAlposName());
   fValue.resize(q2.size());
   fError.resize(q2.size());


   static const double mp2 = pow(0.92, 2);
   if ( y.empty() ) {
      double sqs=DOUBLE_NS(sqrt-s,GetAlposName());
      double ss = sqs*sqs;
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) {
         double x = beta[i]*xpom[i];
         y[i] = q2[i]/(ss-mp2)/x;
         //y[i] = q2[i] / (x[i]*ss); 
      }
   }

   // conceptually, these should be taken as alpos parameters like: PAR(charge)
   // charge = DOUBLE_NS(e-charge,GetAlposName()); 
   // polty  = DOUBLE_NS(e-polarity,GetAlposName());
   charge = PAR(e-charge);
   CONST(e-charge);

   IsRedCS  = BOOL_NS(IsReducedCS,GetAlposName()); // access directly from steering
   IsNC     = BOOL_NS(IsNC,GetAlposName());

   // fDPDF = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".DPDF"));
   // fAs   = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".Alpha_s"));

   return true;
}


// __________________________________________________________________________________________ //
bool AApfelDDISCS::Update() {  //alpos
   debug["Update"]<<"AlposName: "<<GetAlposName()<<endl;

   // diff stuff from steering.
   const double sqrts = DOUBLE_NS(sqrt-s,GetAlposName());
   const double s = sqrts*sqrts;
   const double mp2 = pow(0.92, 2);
   const bool Is4D = EXIST_NS(Is4D,GetAlposName()) && BOOL_NS(Is4D,GetAlposName());
   vector<double> tAbsVal;
   if ( Is4D ) {
      tAbsVal = DOUBLE_COL_NS(Data,tAbs,GetAlposName());
   }

   //flux coefficients
   const double a0_IP = PAR(a0_IP);
   const double ap_IP = PAR(ap_IP);
   const double b0_IP = PAR(b0_IP);
   const double a0_IR = PAR(a0_IR);
   const double ap_IR = PAR(ap_IR);
   const double b0_IR = PAR(b0_IR);
   const double n_IR  = PAR(n_IR);
   const double tAbsMin = EXIST_NS(tAbsMin,GetAlposName()) ? DOUBLE_NS(tAbsMin,GetAlposName()) : 0;
   const double tAbsMax = EXIST_NS(tAbsMax,GetAlposName()) ? DOUBLE_NS(tAbsMax,GetAlposName()) : 1.;



   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(ApfelInit); 
   polty  = PAR(e-polarity);
   const double convfac= 0.389379338e9; //0.389379323e9;
   const double Mz = PAR(ApfelInit.mZ);//APFEL::GetZMass();
   const double Mw = PAR(ApfelInit.mW);//APFEL::GetWMass();
   const double Gf = PAR(ApfelInit.Gf);//APFEL::GetGFermi();
   const double q0 = PAR(ApfelInit.Q0);
   if ( charge!=-1 &&  charge!=1) {
      error["Update"]<<"Could not get charge of lepton."<<endl;
      exit(1);
   }

   if ( IsNC ) APFEL::SetProcessDIS("NC"); // 'EM', 'NC', 'CC'
   else APFEL::SetProcessDIS("CC");

   if ( charge == 1) APFEL::SetProjectileDIS("positron");
   else if ( charge == -1) APFEL::SetProjectileDIS("electron");
   else { error["Update"]<<"Could not get charge of lepton. ch="<<charge<<endl; exit(1);}
   APFEL::SetTargetDIS("proton");
   //APFEL::SelectCharge(string selch)://selects one particular charge in the NC structure functions ('selch' = 'down', 'up', 'strange', 'charm', 'bottom', 'top', 'all', default 'selch' = 'all')

   // --- reggeon cross section
   AParmNamed* AParmReggeonCS = TheoryHandler::Handler()->GetParameter(this->GetAlposName()+std::string(".ReggeonCS"));//->SetValue(VAL,ERR,false);                                        
   const bool AParmReggeonCSIsFunction = AParmReggeonCS->GetAlposName().find(".ReggeonCS") == string::npos;
   debug["Update"]<<"Parameter ReggeonCS is a function: "<<(AParmReggeonCSIsFunction?"true":"false")<<endl;
   const bool H1FitAReg  = !AParmReggeonCSIsFunction && (string(PAR_S(ReggeonCS))=="H1FitA");
   info["Init"]<<"Use ReggeonCS=='H1FitA': "<< (H1FitAReg ? "true" : "false" )<<endl;
   vector<double> sigma_reg; 
   if ( H1FitAReg ) {
      info["Update"]<<"Using pre-calculated reggeon cross sections from H1FitA"<<endl;
   }
   else {
      info["Update"]<<"Setting DataAlposName to reggeon cross section function and get values."<<endl;
      SET_S(ReggeonCS.DataAlposName,GetAlposName(),"");
      sigma_reg = VALUES(ReggeonCS);
   }

   // ------ calc structure functions
   set<double> q2val;
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) q2val.insert(q2[i]);
   const double qin = q0;
   double qfi;
   APFEL::SetPDFSet("external");
   for ( auto qq : q2val ) { // order all data points in Q2
      qfi = sqrt(qq);

      APFEL::ComputeStructureFunctionsAPFEL(qin,qfi);
      //cout << "RADEK_qin_qfi: " << qin <<" "<<  qfi << endl;
      APFEL::SetPDFSet("external");

      /*
      if( qfi > qin ) {
         qin = qfi;
         APFEL::SetPDFSet("apfel");
      }
      else {
         qin = q0;
         APFEL::SetPDFSet("external");
      }
      */
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
         if ( q2[i] != qq ) continue;
         double F2  = APFEL::F2total(beta[i]);
         double FL  = APFEL::FLtotal(beta[i]);
         double xF3 = APFEL::F3total(beta[i]);
         //cout<<"F2: "<<F2<<"\t charge: "<<charge<<endl;
         double x = beta[i]*xpom[i];
         double yy = y[i];//q2[i]/(s-mp2)/x; // y[i]
         double yplus  = 1+(1-yy)*(1-yy);
         double yminus = 1-(1-yy)*(1-yy);

         //cout<<"Apfel   Q2="<<q2[i]<<"\tf2="<<F2<<"\tfl="<<FL<<"\txpom="<<xpom[i]<<"\tbeta="<<beta[i]<<endl;

         // NC
         if ( IsNC ) {

            // --- Reduced x-section for pomeron 
            double flxIP = Is4D ?
               AlposTools::rflux   (a0_IP, ap_IP, b0_IP, xpom[i], tAbsVal[i]):
               AlposTools::rfluxInt(a0_IP, ap_IP, b0_IP, xpom[i], tAbsMin, tAbsMax);
            double xpSigRed_IP =  flxIP*xpom[i] * (F2  - yy*yy/yplus*FL);


            //--- Get the Reggeon structure function from the H12006
            double sigmaReg = 0;
            if ( H1FitAReg ) {
               static bool isFirst = true; //hack for faster calculation
               int ifit = 0;
               if(isFirst) { ifit = 1; isFirst = false; } //Reggeon should be the same for both FitA and FitB
               double xPq[13];
               double f2FitA[2], flFitA[2]; //0 = pomeron, 1 = reggeon
               double c2FitA[2], clFitA[2];
               qcd_2006_(&beta[i], &q2[i],  &ifit, xPq, f2FitA, flFitA, c2FitA, clFitA);
               double F2r = f2FitA[1];
               double FLr = flFitA[1];
               
               sigmaReg = F2r  - yy*yy/yplus*FLr;
               //if ( !H1FitAReg && (i<4 || isnan(sigmaReg) || isnan(sigma_reg[i]) ) ) cout<<"q2="<<q2[i]<<"\txpom="<<xpom[i]<<"\tbbeta="<<beta[i]<<"\tH1FitA="<<sigmaReg<<"\tRegCS="<<sigma_reg[i]<<endl;
            }
            else {
               sigmaReg = sigma_reg[i];
            }
            
            //  reggeon flux
            double flxIR = Is4D ?
               AlposTools::rflux   (a0_IR, ap_IR, b0_IR, xpom[i], tAbsVal[i]):
               AlposTools::rfluxInt(a0_IR, ap_IR, b0_IR, xpom[i], tAbsMin, tAbsMax);
            
            //Reduced x-section for reggeon
            double xpSigRed_IR =  flxIR*xpom[i] * (sigmaReg);

            // --- reduced cross section
            fValue[i] = xpSigRed_IP + n_IR*xpSigRed_IR;
            if ( isnan(fValue[i]) ) {
               cout<<"ISNAN: i="<<i<<endl;
               cout<<"flxIP: "<<flxIP<<endl;
               cout<<"xpom[i]: "<<xpom[i]<<endl;
               cout<<"(F2  - yy*yy/yplus*FL): "<<(F2  - yy*yy/yplus*FL)<<endl;
               cout<<"F2="<<F2<<endl;
               cout<<"FL="<<FL<<endl;
               cout<<"flxIR="<<flxIR<<endl;
               cout<<"xpSigRed_IP: "<<xpSigRed_IP<<endl;
               cout<<"n_IR="<<n_IR<<endl;
               cout<<"xpSigRed_IR="<<xpSigRed_IR<<endl;
               cout<<"FlxIP: "<<a0_IP<<", "<<ap_IP<<", "<<b0_IP<<", "<<xpom[i]<<", "<<tAbsMin<<", "<<tAbsMax<<endl;
               cout<<"Apfel   Q2="<<q2[i]<<"\tf2="<<F2<<"\tfl="<<FL<<"\txpom="<<xpom[i]<<"\tbeta="<<beta[i]<<endl;
               cout<<"q2="<<q2[i]
                   <<"\tyy="<<yy
                   <<"\ty="<<y[i]
                   <<"\txpom="<<xpom[i]<<"\tbeta="<<beta[i]<<endl;

            }

            // non-diff. DIS
	    //xF3 *= -1.*charge;
	    //fValue[i] = F2 + yminus/yplus*xF3 - y[i]*y[i]/yplus*FL;
	 }
	 else if ( !IsNC ) {  // CC
            error["Update"]<<"CC not implemented."<<endl;
            exit(1);
	    // F2 *= 0.5;
	    // FL *= 0.5;
	    // xF3 *= 0.5;
	    // if ( charge == 1 ) {
	    //    fValue[i] = 0.5*(yplus*F2 - yminus*xF3 - y[i]*y[i]*FL);
	    //    fValue[i] *= (1+polty);
	    // }
	    // else if ( charge==-1 ) {
	    //    fValue[i] = 0.5*(yplus*F2 + yminus*xF3 - y[i]*y[i]*FL);
	    //    fValue[i] *=(1-polty);
	    // }
	    // else { cout<<"Error. Wrong charge."<<endl;exit(1); }
	 }

	 // // ------ calc non-reduced CS if needed
	 // if ( !IsRedCS ) {
	 //    if ( IsNC ) {
	 //       double aem = 7.29735e-3; // 1/137.035999074 // 7.29927d-3;//
	 //       fValue[i] *= 2*M_PI*yplus/(beta[i]*q2[i]*q2[i])*convfac*aem*aem;
	 //    }
	 //    else { //CC
	 //       fValue[i] *= pow(Mw,4)/pow(Mw*Mw+q2[i],2)*Gf*Gf/(2*M_PI*beta[i])*convfac;
	 //    }
	 // }
      }
   }
   APFEL::SetPDFSet("external");
   // ------ done
   fError.resize(fValue.size());
   if ( std::isnan(fValue[0])) { // this statement fixes some odd compiler optimizations which may yield to nan coming from APFEL::Fxyz()
      error["Update"]<<endl;
      error["Update"]<<"Cross section[0] is isnan: "<<fValue[0]<<"\t dataset: "<<this->GetAlposName()<<endl; 
      error["Update"]<<endl;
      //TheoryHandler::Handler()->PrintCurrentTheorySet();
      fValue.clear();
      fValue.resize(fError.size());// set all elements to zero ... and continue
      //exit(1);

      for( auto dd : sigma_reg ) cout<<dd<<endl;
      exit(1);

   }
   return true;

}

//______________________________________________________________________________

#endif
