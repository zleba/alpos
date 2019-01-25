
#include "alpos/functions/AApfelxxPDF.h"
#include "apfel/dglapbuilder.h"

#include <iostream>


//! Interface to the Apfel++ PDF evolution

using namespace std;

const std::vector<std::string> AApfelxxPDF::fRequirements = {  
   "xp", // xp (dummy)
   "muf", // mu_f (dummy)
   "Alpha_s", // alpha_s function
   "PDFmu0",  // PDF (function) at the starting scale
   "nx1","xmin1","intdegree1", // 1st grid
   "nx2","xmin2","intdegree2", // 2. grid
   "nx3","xmin3","intdegree3", // 3. grid
   "nx4","xmin4","intdegree4", // 4. grid
   "mu0", // reference scale (evoultion starting scale)
   "iOrd", // order
   "mc","mb","mt", // heavy quark masses
   
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AApfelxxPDF::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelxxPDF::fFunctionName = "ApfelxxPDF"; //< The function's name


// __________________________________________________________________________________________ //
AApfelxxPDF::AApfelxxPDF(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AApfelxxPDF");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
AApfelxxPDF::~AApfelxxPDF() {
}


// ___________________________________________________________________________________________ //
bool AApfelxxPDF::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<endl;
   fAs    = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".Alpha_s"));
   fPDF0  = TheoryHandler::Handler()->GetFuncD(this->GetAlposName()+std::string(".PDFmu0"));


   CONST(nx1);
   CONST(xmin1);
   CONST(intdegree1);
   CONST(nx2);
   CONST(xmin2);
   CONST(intdegree2);
   CONST(nx3);
   CONST(xmin3);
   CONST(intdegree3);
   CONST(nx4);
   CONST(xmin4);
   CONST(intdegree4);

   fGrid = unique_ptr<const apfel::Grid>(new apfel::Grid({
	    apfel::SubGrid{int(PAR(nx1)),PAR(xmin1),int(PAR(intdegree1))}, 
	       apfel::SubGrid{int(PAR(nx2)),PAR(xmin2),int(PAR(intdegree2))}, 
		  apfel::SubGrid{int(PAR(nx3)),PAR(xmin3),int(PAR(intdegree3))}, 
		     apfel::SubGrid{int(PAR(nx4)),PAR(xmin4),int(PAR(intdegree4))}, 
			}));
   //fDglap = unique_ptr<apfel::Dglap>(&new apfel::DglapBuildQCD());
   //const vector<double> Masses = {0, 0, 0, sqrt(2), 4.5, 175};
   const vector<double> Masses = {0, 0, 0, PAR(mc),PAR(mb),PAR(mt)};
   const vector<double> Thresholds = Masses;
   fDglapObj = apfel::InitializeDglapObjectsQCD(*fGrid, Masses, Thresholds);

   CONST(mu0);
   CONST(iOrd);
   
   //fDglap = unique_ptr<apfel::Dglap>(&(apfel::DglapBuildQCD(g, LHToyPDFs, mu0, Masses, Thresholds, PerturbativeOrder, asalpos)));
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelxxPDF::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(2,xp,Q);
   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   double xp  = va_arg(ap, double);
   double muf = va_arg(ap, double);
   va_end(ap);

   return GetQuick({xp,muf});
    
}


// ______________________________________________________________________________________ //
std::vector<double> AApfelxxPDF::GetQuick(const vector<double>& xp_muf) {

   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   std::vector<double> pdf(fValue.size());
   if ( xp_muf.size() != 2) {
      error["GetQuick"]<<"Quick acces is implemented for two parameter which are 'xp' and 'muf'."<<endl;
      return pdf;
   }

   // if calculation failed, return 'fudged' values
   if ( fValue[6]==1 ) { // PDF is nan (see below)
      vector<double> ret{0, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0};
      return ret;
   }    

   for ( int i=0; i<13; i++ )
      pdf[i] = fTabulatedPDFs->EvaluatexQ(i,xp_muf[0],xp_muf[1]);
   
   return AlposTools::LicoApfelxxToLha(pdf);
  
}


// __________________________________________________________________________________________ //
bool AApfelxxPDF::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   // init linear-rotation
   if ( fPdf0ToApfl.GetNcols() == 0 ) {
      // --- get linear combination of inital PDF
      //vector<double> def = {0,0,0,0,0,0,1,0,0,0,0,0,0}; // gluon is by definition 0th element
      //SET_ANY(this->GetAlposName()+".PDFmu0.iPDF",-2,0); // -2 ! this returns full linear-combination including gluon

      SET(PDFmu0.iPDF,-2,0); // set PDFQ0Param to 'def' mode. 
      vector<double> def = VALUES(PDFmu0);

      // vector<double> def = VALUES_ANY(this->GetAlposName()+".PDFmu0");
      
      SET(PDFmu0.iPDF,6,0); //reset
      //SET_ANY(this->GetAlposName()+".PDFmu0.iPDF",6,0); //reset     

      TMatrixD LiCo(13,13,&def[0]);
      // cout<<"LiCo: "<<endl; //debug
      // LiCo.Print(); //debug
      TMatrixD InvLiCo(13,13);
      InvLiCo = AlposTools::InvertLU(LiCo);
      // cout<<"InvLiCo: "<<endl;
      // InvLiCo.Print();
      
      static const std::vector<std::vector<double> > LhToApfl {
	 {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}, // gluon
	 {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}, // SIGMA = \sum_q q^+     
	 { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}, // VALENCE = \sum_q q^-   
	 {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}, // T3  = u^+ - d^+        
	 {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}, // V3  = u^- - d^-        
	 {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}, // T8  = u^+ + d^+ - 2 s^+
	 {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}, // V8  = u^- + d^- - 2 s^-
	 {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}, // T15 = u^+ + d^+ + s^+ - 3 c^+            
	 {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}, // V15 = u^- + d^- + s^- - 3 c^-            
	 {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}, // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+      
	 {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}, // V24 = u^- + d^- + s^- + c^- - 4 b^-      
	 { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}, // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
	 {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,}  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-
      };
      fPdf0ToApfl.Clear();
      fPdf0ToApfl.ResizeTo(13,13);
      for ( int i = 0 ; i<13; i++ ) {
      	 for ( int j = 0 ; j<13; j++ ) {
      	    for ( int k = 0 ; k<13; k++ ) {
      	       fPdf0ToApfl[i][j] += LhToApfl[i][k] * InvLiCo[k][j];
      	    }
      	 }
      }
      for ( int i = 0 ; i<13; i++ ) {
	 for ( int j = 0 ; j<13; j++ ) {
	    if ( fabs(fPdf0ToApfl[i][j])<1.e-14 ) fPdf0ToApfl[i][j]=0;
	 }
      }
      // cout<<"Pdf0ToApfl: "<<endl;
      // fPdf0ToApfl.Print();
   }

   UPDATE(Alpha_s);
   UPDATE(PDFmu0);

   const double mu0 = PAR(mu0);
   const int PerturbativeOrder = PAR(iOrd);
   const vector<double> Masses = {0, 0, 0, PAR(mc),PAR(mb),PAR(mt)};
   const vector<double> Thresholds = Masses;

   // Initialize DGLAP evolution
   //const auto asalpos = [&] (double const& mu) -> double{ return fAs->GetQuick(1,mu)[0]; };
   const auto asalpos = [&] (double const& mu) -> double{ return fAs->GetQuick(vector<double>{mu})[0]; };
   //const auto pdfmu0  = [&] (int const& i, double const& x, double const& ) -> double{ return fPDFmu0->GetQuick(1,mu)[0]; };
   //fDglap = apfel::DglapBuildQCD(*fGrid, LHToyPDFs, mu0, Masses, Thresholds, PerturbativeOrder, asalpos);
   const auto pdf0alpos = [&] (double const& x, double const& mu) -> map<int,double>{
      static const vector<bool> found{true,true,true,true,true,true,true,true,true,true,true,true,true};
      vector<double> xfpdf0(13);
      for ( int ii = 0 ; ii<13 ;ii++ ) {
	 xfpdf0[ii] = fPDF0->GetQuick(vector<double>{double(ii),x,mu})[0];
	 //cout<<"ii="<<ii<<"\tx="<<x<<"\tmu="<<mu<<"\txpdf0[0]="<<xfpdf0[ii]<<endl;
      }
      vector<double> xf13(13);
      AlposTools::Calc13partonsFromLico(xf13,fPdf0ToApfl,xfpdf0,found);
      
      map<int,double> xf13map;
      for ( int ii = 0 ; ii<13 ;ii++ )
	 xf13map.insert({ii,xf13[ii]});
      return xf13map;
   };
   //fDglap = apfel::DglapBuildQCD(*fGrid, pdf0alpos, mu0, Masses, Thresholds, PerturbativeOrder, asalpos);
   //apfel::DglapObjects fDglapObj = InitializeDglapObjectsQCD(*fGrid);
   //fDglap = apfel::DglapBuild(fDglapObj, pdf0alpos, mu0, Masses, Thresholds, PerturbativeOrder, asalpos);
   fDglap = apfel::BuildDglap(fDglapObj, pdf0alpos, mu0, PerturbativeOrder, asalpos);

   // Tabulate PDFs
   //const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*fDglap, 50, 1, 1000, 3};
   fTabulatedPDFs = unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>>(new apfel::TabulateObject<apfel::Set<apfel::Distribution>>{*fDglap, 50, 1, 1000, 3});


   fValue.resize(13);
   double xp = PAR(xp);
   double muf= PAR(muf);
   for ( int i=0; i<13; i++ ) {
      fValue[i] = fTabulatedPDFs->EvaluatexQ(i,xp,muf);
   }
   
   fValue = AlposTools::LicoApfelxxToLha(fValue);

   if (  !isfinite(fValue[6]) ) {
      error["Update"]<<"non-finite value detected."<<endl;
      error["Update"]<<"muf="<<muf<<"\txp="<<xp<<endl;
      error["Update"]<<"Return values from fTabulatedPDFs:"<<endl;
      for ( int i=0; i<13; i++ ) {
	 cout<<"\t"<<fTabulatedPDFs->EvaluatexQ(i,xp,muf)<<endl;
      }
      error["Update"]<<"Alpos PDFmu0:"<<endl;
      for ( int i=0; i<13; i++ ) {
	 cout<<"\t"<<pdf0alpos(xp,muf)[i]<<endl;
      }
      error["Update"]<<"Continueing with fudged values."<<endl;
      fValue = vector<double>{0, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0};
   }

   return true;
}

