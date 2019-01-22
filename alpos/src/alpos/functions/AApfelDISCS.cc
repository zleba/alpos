// DB 16.01.2015
// DB Updates: 10/15

#include "APFEL/APFEL.h"
#include "alpos/functions/AApfelDISCS.h"
#include "alpos/functions/AApfel.h"
#include <iostream>
#include "fastnlotk/read_steer.h"
#include <set>

using namespace std;

// __________________________________________________________________________________________ //
const std::vector<std::string> AApfelDISCS::fRequirements = {
   "ApfelInit",
   "e-charge",
   "e-polarity",
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AApfelDISCS::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelDISCS::fFunctionName = "ApfelDISCS"; //!< The function's name

// __________________________________________________________________________________________ //
AApfelDISCS::AApfelDISCS(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AApfelDISCS");
}


// __________________________________________________________________________________________ //
AApfelDISCS::~AApfelDISCS() {
}


// ___________________________________________________________________________________________ //
bool AApfelDISCS::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   AlposObject::debug["Init"]<<endl;

   // get points
   q2 = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   x = DOUBLE_COL_NS(Data,x,GetAlposName());
   y = DOUBLE_COL_NS(Data,y,GetAlposName());
   fValue.resize(q2.size());
   fError.resize(q2.size());

   if ( y.empty() ) {
      double sqs=DOUBLE_NS(sqrt-s,GetAlposName());
      double ss = sqs*sqs;
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) { 
	 y[i] = q2[i] / (x[i]*ss);
      }
   }
   // conceptually, these should be taken as alpos parameters like: PAR(charge)
   // charge = DOUBLE_NS(e-charge,GetAlposName()); 
   // polty  = DOUBLE_NS(e-polarity,GetAlposName());
   charge = PAR(e-charge);
   
   IsRedCS  = BOOL_NS(IsReducedCS,GetAlposName()); // access directly from steering
   IsNC     = BOOL_NS(IsNC,GetAlposName());

   return true;
}


// __________________________________________________________________________________________ //
bool AApfelDISCS::Update() {  //alpos
   debug["Update"]<<"AlposName: "<<GetAlposName()<<endl;

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
   // APFEL::SetTargetDIS("proton");
   //APFEL::SelectCharge(string selch)://selects one particular charge in the NC structure functions ('selch' = 'down', 'up', 'strange', 'charm', 'bottom', 'top', 'all', default 'selch' = 'all')

   // ------ calc structure functions
   set<double> q2val;
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) q2val.insert(q2[i]);
   double qin = q0;
   double qfi;
   APFEL::SetPDFSet("external");
   for ( auto qq : q2val ) {
      qfi = sqrt(qq);
      APFEL::ComputeStructureFunctionsAPFEL(qin,qfi);
      if( qfi > qin ) {
	 qin = qfi;
	 APFEL::SetPDFSet("apfel");
      }
      else {
	 qin = q0;
	 APFEL::SetPDFSet("external");
      }
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 if ( q2[i] != qq ) continue;
	 double F2  = APFEL::F2total(x[i]);
	 double FL  = APFEL::FLtotal(x[i]);
	 double xF3 = APFEL::F3total(x[i]);
	 //cout<<"F2: "<<F2<<"\t charge: "<<charge<<endl;
	 double yplus  = 1+(1-y[i])*(1-y[i]);
	 double yminus = 1-(1-y[i])*(1-y[i]);

	 if ( IsNC ) {
	    xF3 *= -1.*charge;
	    fValue[i] = F2 + yminus/yplus*xF3 - y[i]*y[i]/yplus*FL;
	 }
	 else if ( !IsNC ) {  // CC
	    F2 *= 0.5;
	    FL *= 0.5;
	    xF3 *= 0.5;
	    if ( charge == 1 ) {
	       fValue[i] = 0.5*(yplus*F2 - yminus*xF3 - y[i]*y[i]*FL);
	       fValue[i] *= (1+polty);
	    }
	    else if ( charge==-1 ) {
	       fValue[i] = 0.5*(yplus*F2 + yminus*xF3 - y[i]*y[i]*FL);
	       fValue[i] *=(1-polty);
	    }
	    else { cout<<"Error. Wrong charge."<<endl;exit(1); }
	 }

	 // ------ calc non-reduced CS if needed
	 if ( !IsRedCS ) {
	    if ( IsNC ) {
	       double aem = 7.29735e-3; // 1/137.035999074 // 7.29927d-3;//
	       fValue[i] *= 2*M_PI*yplus/(x[i]*q2[i]*q2[i])*convfac*aem*aem;
	    }
	    else { //CC
	       fValue[i] *= pow(Mw,4)/pow(Mw*Mw+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	    }
	 }
      }
   }
   APFEL::SetPDFSet("external");
   // ------ done
   fError.resize(fValue.size());
   if ( std::isnan(fValue[0])) { // this statement fixes some odd compiler optimizations which may yield to nan coming from APFEL::Fxyz()
      error["Update"]<<endl;
      error["Update"]<<"Cross section is isnan: "<<fValue[0]<<"\t dataset: "<<this->GetAlposName()<<endl; 
      error["Update"]<<endl;
      TheoryHandler::Handler()->PrintCurrentTheorySet();
      fValue.clear();
      fValue.resize(fError.size());// set all elements to zero ... and continue
      //exit(1);
   }
   return true;

}

//______________________________________________________________________________

