// DB 01.2019

#include <iostream>
#include <set>

#include "alpos/functions/AApfelxxDDISCS.h"
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

// LH Toy PDFs
/*
code duplication with AApfelxxDISCS
double xupv(double const& x)  { return 5.107200 * pow(x,0.8) * pow((1-x),3); }
double xdnv(double const& x)  { return 3.064320 * pow(x,0.8) * pow((1-x),4); }
double xglu(double const& x)  { return 1.7 * pow(x,-0.1) * pow((1-x),5); }
double xdbar(double const& x) { return 0.1939875 * pow(x,-0.1) * pow((1-x),6); }
double xubar(double const& x) { return xdbar(x) * (1-x); }
double xsbar(double const& x) { return 0.2 * ( xdbar(x) + xubar(x) ); }
double LHToyPDFs(int const& i, double const& x, double const&)
{
   // Gluon
   if      (i == 0)
      return xglu(x);
   // Singlet, T15, T24, T35
   else if (i == 1 || i == 7 || i == 9 || i == 11 )
      return xdnv(x) + 2 * xdbar(x) + xupv(x) + 2 * xubar(x) + 2 * xsbar(x);
   // T3
   else if (i == 3)
      return xupv(x) + 2 * xubar(x) - xdnv(x) - 2 * xdbar(x);
   // T8
   else if (i == 5)
      return xupv(x) + 2 * xubar(x) + xdnv(x) + 2 * xdbar(x) - 4 * xsbar(x);
   // Valence, V8, V15, V24, V35
   else if (i == 2 || i == 6 || i == 8 || i == 10 || i == 12)
      return xupv(x) + xdnv(x);
   // V3
   else if (i == 4)
      return xupv(x) - xdnv(x);
   else
      return 0;
}
*/


// __________________________________________________________________________________________ //
const std::vector<std::string> AApfelxxDDISCS::fRequirements = {
   "DPDF", // evolved PDF function
   "Alpha_s", // alpha_s evolution function
   "e-charge",
//   "e-polarity",
   "iOrd",
   "mc","mb","mt",
   "nGridFac"
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AApfelxxDDISCS::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelxxDDISCS::fFunctionName = "ApfelxxDDISCS"; //!< The function's name

// __________________________________________________________________________________________ //
AApfelxxDDISCS::AApfelxxDDISCS(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AApfelxxDDISCS");
}


// __________________________________________________________________________________________ //
AApfelxxDDISCS::~AApfelxxDDISCS() {
}


// ___________________________________________________________________________________________ //
bool AApfelxxDDISCS::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;

   // CONST(DPDF);
   // CONST(Alpha_s);

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

   fDPDF = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".DPDF"));
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
bool AApfelxxDDISCS::Update() {  //alpos
   debug["Update"]<<"AlposName: "<<GetAlposName()<<endl;

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   polty  = PAR(e-polarity);

   UPDATE(Alpha_s);

   // --- Apfel++ following structurefunction_test.cc
   // Initializers
   //const double mu0 = sqrt(2);
   const int PerturbativeOrder = PAR(iOrd);
   const vector<double> Masses = {0, 0, 0, PAR(mc), PAR(mb), PAR(mt)};
   const vector<double> Thresholds = Masses;

   // --- Running coupling
   //const double AlphaQCDRef = 0.35;
   //const double MuAlphaQCDRef = sqrt(2);
   //apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Masses, PerturbativeOrder};
   //const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
   //const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
   // --- Running coupling alpos
   const auto asalpos = [&] (double const& mu) -> double{ return fAs->GetQuick(vector<double>{mu})[0]; };

   // --- PDFs
   // Initialize DGLAP evolution
   // const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};
   // auto EvolvedPDFs = apfel::DglapBuildQCD(g, LHToyPDFs, mu0, Masses, Thresholds, PerturbativeOrder, as);
   // // Tabulate PDFs 
   // const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

   // Evolved PDFs
   //const auto PDFs = [&] (int const& i, double const& x, double const& Q) -> double{ return TabulatedPDFs.EvaluatexQ(i,x,Q); };
   // --- evolved alpos
   //const auto PDFsalpos = [&] (int const& i, double const& x, double const& mu) -> double{
   //   return AlposTools::LicoLhaToApfelxx(fPDF->GetQuick({x,mu}))[i];
   //   // vector<double> xfx = fPDF->GetQuick({x,mu});
   //   // vector<double> xfxApfl = ;
   //   // return xfxApfl[i];
   //   //return fPDF->GetQuick({x,mu})[i]; 
   //};
   // const auto PDFsalpos = [&] (int const& i, double const& x, double const& mu) -> double{ 
   //    return AlposTools::LicoLhaToApfelxx(fPDF->GetQuick({x,mu}))[i];
   const auto PDFsalpos = [&] (double const& x, double const& mu) -> map<int,double>{
      return AlposTools::LicoLhaToApfelxxMap(fDPDF->GetQuick({0,x,mu}));
      // vector<double> xfx = fPDF->GetQuick({x,mu});
      // vector<double> xfxApfl = ;
      // return xfxApfl[i];
      //return fPDF->GetQuick({x,mu})[i]; 
   };


   // Charges
   // function<vector<double>(double const&)> fBq = [Thresholds] (double const& Q) -> vector<double> {
   //    vector<double> Bq(Thresholds.size());
   //    for (auto i = 0; i < (int) Thresholds.size(); i++)
   // 	 Bq[i] = (Q > Thresholds[i] ? apfel::QCh2[i] : 0);
   //    return Bq;
   // };
   // function<vector<double>(double const&)> fDq = [Thresholds] (double const&) -> vector<double>{ return {0, 0, 0, 0, 0, 0}; };
   // --- Charges w/ EW contributions
   // Relevant constants for the computation of the EW charges.
   const double MZ         = 91.1876;
   const double MZ2        = MZ * MZ;
   const double Sin2ThetaW = 0.23126;
   const double VD         = - 0.5 + 2. * Sin2ThetaW / 3.;
   const double VU         = + 0.5 - 4. * Sin2ThetaW / 3.;
   const vector<double> Vq = {VD, VU, VD, VU, VD, VU};
   const double AD         = - 0.5;
   const double AU         = + 0.5;
   const vector<double> Aq = {AD, AU, AD, AU, AD, AU};
   // Unpolarized electron target.
   const int ie     = - 1; // Electron
   const double Ve  = - 0.5 + 2. * Sin2ThetaW;
   const double Ae  = - 0.5;
   const double pol = 0;   // No polarization
   // Effective charges.
   function<vector<double>(double const&)> fBq = [=] (double const& Q) -> vector<double>  {
      const double Q2  = Q * Q;
      const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
      const double PZ2 = PZ * PZ;
      vector<double> Bq;
      for (auto i = 0; i < (int) Thresholds.size(); i++)
      {
	 const double b = apfel::QCh2[i] 
	    - 2 * apfel::QCh[i] * Vq[i] * ( Ve + ie * pol * Ae ) * PZ
	    + ( Ve * Ve + Ae * Ae + ie * pol * 2 * Ve * Ae )
	    * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
	 Bq.push_back((Q > Thresholds[i] ? b : 0));
      }
      return Bq;
   };
   function<vector<double>(double const&)> fDq = [=] (double const& Q) -> vector<double> {
      const double Q2  = Q * Q;
      const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
      const double PZ2 = PZ * PZ;
      vector<double> Dq;
      for (auto i = 0; i < (int) Thresholds.size(); i++)
      {
	 const double d = - 2 * apfel::QCh[i] * Aq[i] * ( Ae + ie * pol * Ve ) * PZ
	    + 2 * Vq[i] * Aq[i] * ( 2 * Ve * Ae 
				    + ie * pol * ( Ve * Ve + Ae * Ae ) ) * PZ2;
	 Dq.push_back((Q > Thresholds[i] ? d : 0));
      }
      return Dq;
   };
   // --- --- --- 


   // --- Apfel++ following structurefunction_test.cc
   //  Initialize structure functions
   // const auto F2 = apfel::StructureFunctionBuildNC(fF2Obj, PDFsalpos, Thresholds, PerturbativeOrder, asalpos, fBq);
   // const auto FL = apfel::StructureFunctionBuildNC(fFLObj, PDFsalpos, Thresholds, PerturbativeOrder, asalpos, fBq);
   // const auto F3 = apfel::StructureFunctionBuildNC(fF3Obj, PDFsalpos, Thresholds, PerturbativeOrder, asalpos, fDq);
   const auto F2 = apfel::BuildStructureFunctions(fF2Obj, PDFsalpos, PerturbativeOrder, asalpos, fBq);
   const auto FL = apfel::BuildStructureFunctions(fFLObj, PDFsalpos, PerturbativeOrder, asalpos, fBq);
   //const auto F3 = apfel::BuildStructureFunctions(fF3Obj, PDFsalpos, PerturbativeOrder, asalpos, fDq);


   // calculate q2 evolution only once.
   // take care for xpom !!
   map<pair<double,double>,apfel::Distribution> q2_f2, q2_fl, q2_f3;
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
      double Q   = sqrt(q2[i]);
      SET(DPDF.xpom,xpom[i],0); 
      UPDATE(DPDF);
      if ( q2_f2.count({q2[i],xpom[i]}) == 0 ) {
         q2_f2.insert({{q2[i],xpom[i]},F2.at(0).Evaluate(Q)});
         q2_fl.insert({{q2[i],xpom[i]},FL.at(0).Evaluate(Q)});
         //q2_f3.insert({{q2[i],xpom[i]},F3.at(0).Evaluate(Q)});
      }
   }

   // ------ calc structure functions
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
      double Q = sqrt(q2[i]);
      double yplus  = 1+(1-y[i])*(1-y[i]);
      double yminus = 1-(1-y[i])*(1-y[i]);

      //    SET(DPDF.xpom,xpom[i],0); 
      //    UPDATE(DPDF);
      //    double f2 = F2.at(0).Evaluate(Q).Evaluate(beta[i]);
      //    double fl = FL.at(0).Evaluate(Q).Evaluate(beta[i]);
      //    double f3 = 0;//F3.at(0).Evaluate(Q).Evaluate(beta[i]);
      //    cout<<"q2="<<q2[i]<<"\tf2="<<f2<<"\tfl="<<fl<<endl;

      double f2 = q2_f2.at({q2[i],xpom[i]}).Evaluate(beta[i]);
      double fl = q2_fl.at({q2[i],xpom[i]}).Evaluate(beta[i]);
      double f3 = 0;//q2_f3.at({q2[i],xpom[i]}).Evaluate(beta[i]);

      //cout<<"Apfel++ Q2="<<Q*Q<<"\tf2="<<f2<<"\tfl="<<fl<<"\tf3="<<f3<<"\tq="<<charge<<"\tpol="<<polty<<endl;
   
      if ( IsNC ) {
	 f3 *= -1.*charge;
	 fValue[i] = f2 + yminus/yplus*f3 - y[i]*y[i]/yplus*fl;
	 fValue[i] *= xpom[i];
      }
      else if ( !IsNC ) {  // CC
	 error["Update"]<<"ERROR CC not not validated!"<<endl;
	 exit(3);
	 f2 *= 0.5;
	 fl *= 0.5;
	 f3 *= 0.5;
	 if ( charge == 1 ) {
	    fValue[i] = 0.5*(yplus*f2 - yminus*f3 - y[i]*y[i]*fl);
	    fValue[i] *= (1+polty);
	 }
	 else if ( charge==-1 ) {
	    fValue[i] = 0.5*(yplus*f2 + yminus*f3 - y[i]*y[i]*fl);
	    fValue[i] *=(1-polty);
	 }
	 else { cout<<"Error. Wrong charge."<<endl;exit(1); }
      }

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

