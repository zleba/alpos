// DB 01.2019

#include <iostream>
#include <set>

#include "alpos/functions/AApfelxxReggeonCS.h"
#include "fastnlotk/read_steer.h"

#include <apfel/dglapbuilder.h>
#include <apfel/structurefunctionbuilder.h>
#include <apfel/grid.h>
#include <apfel/timer.h>
#include <apfel/tools.h>
#include <apfel/alphaqcd.h>
#include <apfel/tabulateobject.h>
#include <apfel/constants.h>

using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AApfelxxReggeonCS::fRequirements = {
   "RegPDF", // evolved PDF function
   "Alpha_s", // alpha_s evolution function
   "e-charge",
   "iOrd",
   "mc","mb","mt",
   "nGridFac",
   "DataAlposName"
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AApfelxxReggeonCS::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelxxReggeonCS::fFunctionName = "ApfelxxReggeonCS"; //!< The function's name

// __________________________________________________________________________________________ //
AApfelxxReggeonCS::AApfelxxReggeonCS(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AApfelxxReggeonCS");
}


// __________________________________________________________________________________________ //
AApfelxxReggeonCS::~AApfelxxReggeonCS() {
}


// ___________________________________________________________________________________________ //
bool AApfelxxReggeonCS::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;

   fPDF = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".RegPDF"));
   fAs   = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".Alpha_s"));
   int nGridFac = PAR(nGridFac);
   CONST(nGridFac);
   //const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};
   fGrid = unique_ptr<const apfel::Grid>(new apfel::Grid({
	    //apfel::SubGrid{10*nGridFac,1e-3,3}, 
	    apfel::SubGrid{5*nGridFac,1e-3,3}, 
	    apfel::SubGrid{5*nGridFac,1e-1,3}, 
	    apfel::SubGrid{5*nGridFac,6e-1,3}, 
	    apfel::SubGrid{5*nGridFac,8e-1,3}
	 }));
   const vector<double> Thresholds = {0, 0, 0, PAR(mc), PAR(mb), PAR(mt)};
   fF2Obj = apfel::InitializeF2NCObjectsZM(*fGrid, Thresholds);
   fFLObj = apfel::InitializeFLNCObjectsZM(*fGrid, Thresholds);
   fF3Obj = apfel::InitializeF3NCObjectsZM(*fGrid, Thresholds);

   return true;
}


// __________________________________________________________________________________________ //
bool AApfelxxReggeonCS::Update() {  //alpos
   //debug["Update"]<<"AlposName: "<<GetAlposName()<<", with DataAlposName: "<<PAR_S(DataAlposName)<<endl;

   cout<<"update"<<endl;
   AParmNamed* AParmDataAlposName = TheoryHandler::Handler()->GetParameter(this->GetAlposName()+std::string(".DataAlposName"));//->SetValue(VAL,ERR,false);
   cout<<"ParmName="<<AParmDataAlposName<<endl;
   cout<<"DataAlposName->GetAlposName()="<<AParmDataAlposName->GetAlposName()<<endl;
   AParmS* AParmSDataAlposName = (AParmS*)AParmDataAlposName;
   cout<<"AParmSDataAlposName. GetValue: "<<AParmSDataAlposName->GetValue()<<endl;

   string DataAlposName = PAR_S(DataAlposName);
   cout<<"DataAlposName = "<<DataAlposName<<endl;
   const bool DataPresent = EXIST_NS(Data,DataAlposName);
   info["Update"]<<"DataAlposName: "<<DataAlposName<<", Found 'Data': "<<(DataPresent?"true":"false")<<endl;
   if ( !DataPresent) {
      warn["Update"]<<"Could not read 'Data' table from steering namespace '"<<DataAlposName<<"'. Skipping calculation."<<endl;
      fValue.resize(1);
      fError.resize(1);
      return true;
   }

   // get points
   vector<double> q2 = DOUBLE_COL_NS(Data,Q2,DataAlposName);
   vector<double> y = DOUBLE_COL_NS(Data,y,DataAlposName);
   vector<double> xpom = DOUBLE_COL_NS(Data,xp,DataAlposName);
   vector<double> beta = DOUBLE_COL_NS(Data,beta,DataAlposName);
   vector<double> sigmaData = DOUBLE_COL_NS(Data,Sigma,DataAlposName);
   fValue.resize(q2.size() * 1);
   fError.resize(q2.size() * 1);

   static const double mp2 = pow(0.92, 2);
   if ( y.empty() ) {
      double sqs=DOUBLE_NS(sqrt-s,DataAlposName);
      double ss = sqs*sqs;
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) { 
	 double x = beta[i]*xpom[i];
	 y[i] = q2[i]/(ss-mp2)/x;
	    //y[i] = q2[i] / (x[i]*ss);
      }
   }
   // conceptually, these should be taken as alpos parameters like: PAR(charge)
   // charge = DOUBLE_NS(e-charge,DataAlposName); 
   // polty  = DOUBLE_NS(e-polarity,DataAlposName);
   charge = PAR(e-charge);


   UPDATE(Alpha_s);
   UPDATE(RegPDF);

   // --- Apfel++ following structurefunction_test.cc
   // Initializers
   //const double mu0 = sqrt(2);
   const int PerturbativeOrder = PAR(iOrd);
   const vector<double> Masses = {0, 0, 0, PAR(mc), PAR(mb), PAR(mt)};
   const vector<double> Thresholds = Masses;

   // --- Running coupling
   const auto asalpos = [&] (double const& mu) -> double{ 
      return fAs->GetQuick(vector<double>{mu})[0]; 
   };

   // --- PDFs
   const auto PDFsalpos = [&] (double const& x, double const& mu) -> map<int,double>{
      return AlposTools::LicoLhaToApfelxxMap(fPDF->GetQuick({x,mu}));
   };


   // Charges
   function<vector<double>(double const&)> fBq = [Thresholds] (double const& Q) -> vector<double> {
      vector<double> Bq(Thresholds.size());
      for (auto i = 0; i < (int) Thresholds.size(); i++)
   	 Bq[i] = (Q > Thresholds[i] ? apfel::QCh2[i] : 0);
      return Bq;
   };
   // function<vector<double>(double const&)> fDq = [Thresholds] (double const&) -> vector<double>{ return {0, 0, 0, 0, 0, 0}; };
   // --- Charges w/ EW contributions
   // Relevant constants for the computation of the EW charges.

   // const double MZ         = 91.1876;
   // const double MZ2        = MZ * MZ;
   // const double Sin2ThetaW = 0.23126;
   // const double VD         = - 0.5 + 2. * Sin2ThetaW / 3.;
   // const double VU         = + 0.5 - 4. * Sin2ThetaW / 3.;
   // const vector<double> Vq = {VD, VU, VD, VU, VD, VU};
   // const double AD         = - 0.5;
   // const double AU         = + 0.5;
   // const vector<double> Aq = {AD, AU, AD, AU, AD, AU};
   // // Unpolarized electron target.
   // const int ie     = - 1; // Electron
   // const double Ve  = - 0.5 + 2. * Sin2ThetaW;
   // const double Ae  = - 0.5;
   // const double pol = 0;   // No polarization
   // // Effective charges.
   // function<vector<double>(double const&)> fBq = [=] (double const& Q) -> vector<double>  {
   //    const double Q2  = Q * Q;
   //    const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
   //    const double PZ2 = PZ * PZ;
   //    vector<double> Bq;
   //    for (auto i = 0; i < (int) Thresholds.size(); i++)
   //    {
   //       const double b = apfel::QCh2[i] 
   //          - 2 * apfel::QCh[i] * Vq[i] * ( Ve + ie * pol * Ae ) * PZ
   //          + ( Ve * Ve + Ae * Ae + ie * pol * 2 * Ve * Ae )
   //          * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
   //       Bq.push_back((Q > Thresholds[i] ? b : 0));
   //    }
   //    return Bq;
   // };
   // function<vector<double>(double const&)> fDq = [=] (double const& Q) -> vector<double> {
   //    const double Q2  = Q * Q;
   //    const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
   //    const double PZ2 = PZ * PZ;
   //    vector<double> Dq;
   //    for (auto i = 0; i < (int) Thresholds.size(); i++)
   //    {
   //       const double d = - 2 * apfel::QCh[i] * Aq[i] * ( Ae + ie * pol * Ve ) * PZ
   //          + 2 * Vq[i] * Aq[i] * ( 2 * Ve * Ae 
   //      			    + ie * pol * ( Ve * Ve + Ae * Ae ) ) * PZ2;
   //       Dq.push_back((Q > Thresholds[i] ? d : 0));
   //    }
   //    return Dq;
   // };
   // --- --- --- 


   // --- Apfel++ following structurefunction_test.cc
   //  Initialize structure functions
   const auto F2 = apfel::BuildStructureFunctions(fF2Obj, PDFsalpos, PerturbativeOrder, asalpos, fBq);
   const auto FL = apfel::BuildStructureFunctions(fFLObj, PDFsalpos, PerturbativeOrder, asalpos, fBq);
   //const auto F3 = apfel::BuildStructureFunctions(fF3Obj, PDFsalpos, PerturbativeOrder, asalpos, fDq);

   
   // calculate q2 evolution only once.
   // take care for xpom !!
   map<pair<double,double> ,double > all_q2_beta_F2;
   map<pair<double,double> ,double > all_q2_beta_FL;

   // ------ calc structure functions
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
      double Q = sqrt(q2[i]);
      double yplus  = 1+(1-y[i])*(1-y[i]);
      double yminus = 1-(1-y[i])*(1-y[i]);

      if ( all_q2_beta_F2.count({q2[i],beta[i]}) == 0 ) {
         all_q2_beta_F2[{q2[i],beta[i]}] = F2.at(0).Evaluate(Q).Evaluate(beta[i]);
         all_q2_beta_FL[{q2[i],beta[i]}] = FL.at(0).Evaluate(Q).Evaluate(beta[i]);
      }
      
      double f2 = all_q2_beta_F2[{q2[i],beta[i]}];
      double fl = all_q2_beta_FL[{q2[i],beta[i]}];
      double f3 = 0;//q2_f3.at({q2[i],xpom[i]}).Evaluate(beta[i]);

      //cout<<"Apfel++Reggeon Q2="<<Q*Q<<"\tf2="<<f2<<"\tfl="<<fl<<"\tf3="<<f3<<"\tq="<<charge<<endl;
   
      //if ( IsNC ) {
      f3 *= -1.*charge;
      fValue[i] = f2 + yminus/yplus*f3 - y[i]*y[i]/yplus*fl;
      //fValue[i] *= xpom[i]; // ?

      // // ------ calc non-reduced CS if needed
      // if ( !IsRedCS ) {
      // 	 if ( IsNC ) {
      // 	    double aem = 7.29735e-3; // 1/137.035999074 // 7.29927d-3;//
      // 	    fValue[i] *= 2*M_PI*yplus/(x[i]*q2[i]*q2[i])*convfac*aem*aem;
      // 	 }
      // 	 else { //CC
      // 	    fValue[i] *= pow(Mw,4)/pow(Mw*Mw+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
      // 	 }
      // }
   }

   // ------ done
   fError.resize(fValue.size());
   if ( std::isnan(fValue[0])) { // this statement fixes some odd compiler optimizations which may yield to nan coming from APFEL::Fxyz()
      error["Update"]<<endl;
      error["Update"]<<"Cross section is isnan: "<<fValue[0]<<"\t dataset: "<<this->GetAlposName()<<endl; 
      error["Update"]<<endl;
      TheoryHandler::Handler()->PrintCurrentTheorySet();
      fValue.clear();
      fValue.resize(fError.size());// set all elements to zero ... and continue
      exit(1);
   }
   return true;

}

//______________________________________________________________________________

