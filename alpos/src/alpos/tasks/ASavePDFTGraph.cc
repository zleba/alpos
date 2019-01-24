#include "alpos/tasks/ASavePDFTGraph.h"
#include <iostream>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "alpos/AFactory.h"
#include "alpos/AChisq.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/functions/ALhapdf6.h"

/** 
 ASavePDFTGraph

 * Save TGraphs of PDFs

 */

using namespace std;

const string ASavePDFTGraph::fTaskType = "SavePDFTGraph";

//____________________________________________________________________________________ //
ASavePDFTGraph::ASavePDFTGraph(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("ASavePDFTGraph");
   //! Important: create always new result-object here!
   fResult = new ASavePDFTGraphResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
ASavePDFTGraph::~ASavePDFTGraph(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool ASavePDFTGraph::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   debug["Init"]<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool ASavePDFTGraph::Execute(){
   debug["Execute"]<<endl;

   // init output
   //TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   TDirectory* file = Alpos::Current()->Settings()->rootoutput;
   file->mkdir(GetTaskName().c_str())->cd();
   string pdffunc = STRING_NS(PDF,NS());
   PAR_ANY(pdffunc); // update
   info["Execute"]<<"Evaluating PDF from theory parameter: "<<pdffunc<<" which is of type '"<<
      TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName()<<"'."<<endl;

   // x and q2 spacing
   vector<double> q2val = DOUBLE_ARR_NS(Q2Values,NS());
   int nx = INT_NS(nx,NS());
   double xmin = DOUBLE_NS(xmin,NS());
   double lxstep = (log(1)-log(xmin))/nx;
   vector<double> xval = {xmin};
   for ( int ix = 0 ; ix < nx-1 ; ix++ ) {
      xval.push_back(exp((log(xval.back())+lxstep)));
   }
   xval.push_back(1);

   int nPDF=1;
   bool LHAType=false;
   if ( TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName().find(ALhapdf6::fFunctionName) !=std::string::npos ) {
      ALhapdf6* lha = (ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc);
      nPDF = lha->GetPDFSet()->size();
      info["Execute"]<<"Found "<<nPDF<<" PDF members."<<endl;
      LHAType=true;
   }

   int Mem0 = -1; // PDFSet in the beginning
   if ( LHAType ) Mem0 = PAR_ANY(pdffunc+".PDFSet");


   map<string,vector<double> > pdfdef = GetPDFdef();


   if ( nPDF > 1 && !LHAType ) warn["Execute"]<<"only loop over LHAPDF members impemented."<<endl;
   map<string,map<double,map<double, vector<double> > > > AllValues; // <parname,q2,xp,val>
   for ( int iMem = 0 ;iMem<nPDF ; iMem++ ) { // iMem-loop is slowest!
      if ( LHAType ) { SET_ANY(pdffunc+".PDFSet",iMem,0); }
      PAR_ANY(pdffunc); // update PDF

      for ( auto q2 : q2val ) {
	 TString dirname = Form("Q2_%.1f",q2);
	 if ( !file->GetDirectory(dirname) ) file->mkdir(dirname);
	 file->GetDirectory(dirname)->mkdir(Form("PDF_%d",iMem))->cd();

	 for ( auto ipdf : pdfdef ) {

	    TGraph gPDF;
	    gPDF.SetName(ipdf.first.c_str());
	    for ( auto xp : xval ) {
	       
	       //SET_ANY("LHAPDF.xp",xp,0);
	       //vector<double> xfx = VALUES_ANY("LHAPDF");
	       vector<double> xp_muf={xp,sqrt(q2)};
	       vector<double> xfx = QUICK_ANY(pdffunc,xp_muf);
	       double xf = 0;
	       for ( unsigned int jf=0 ; jf< ipdf.second.size() ; jf++ ) {
		  xf += ipdf.second[jf]*xfx[jf];
	       }
	       gPDF.SetPoint(gPDF.GetN(),xp,xf);
	       AllValues[ipdf.first][q2][xp].push_back(xf);
	    }
	    gPDF.Write();
	 }
      }
   }

   // errors
   if ( LHAType ) {
      ALhapdf6* lha = (ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc);
      for ( auto q2 : q2val ) {
	 TString dirname = Form("Q2_%.1f",q2);
	 file->GetDirectory(dirname)->mkdir("PDF_Errors")->cd();
	 for ( auto ipdf : pdfdef ) {
	    TGraphAsymmErrors gPDF;
	    gPDF.SetName(ipdf.first.c_str());
	    for ( auto xp : xval ) {
	       const vector<double>& values = AllValues[ipdf.first][q2][xp];
	       //cout<<"ipdf.first="<<ipdf.first<<"\tiMem="<<iMem<<"\tq2="<<q2<<"\txp="<<xp<<"\tsize="<<values.size()<<endl;
	       LHAPDF::PDFUncertainty PDFUnc = lha->GetPDFSet()->uncertainty(values);
	       int nP = gPDF.GetN();
	       //cout<<"nP="<<nP<<"\txp="<<xp<<"\tcentral="<<PDFUnc.central<<"\tup="<<PDFUnc.errminus<<"\tdn="<<PDFUnc.errplus<<endl;
	       gPDF.SetPoint(nP,xp,PDFUnc.central);
	       gPDF.SetPointError(nP,0,0,PDFUnc.errminus,PDFUnc.errplus);
	    }
	    gPDF.Write();
	 }
      }
   }
   else {
      info["Execute"]<<"No error TGraphs are drawn."<<endl;
   }

   if ( Mem0>=0 ) SET_ANY(pdffunc+".PDFSet",Mem0,0);

   // close rootfile
   file->Write();
   //if ( file != Alpos::Current()->Settings()->rootoutput) file->Close();
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk.\n"<<endl;
   return true;
}


//____________________________________________________________________________________ //
map<string,vector<double> > ASavePDFTGraph::GetPDFdef() const {
   //! return definition for all PDF linear combinations
   map<string,vector<double> > pdfdef;
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   pdfdef[   "tb"] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "bb"] = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "cb"] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "sb"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "ub"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "db"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "g"] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "d"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ,0, 0 ,0 };
   pdfdef[    "u"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[    "s"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[    "c"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   pdfdef[    "b"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 1 ,0 };
   pdfdef[    "t"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,1 };
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   pdfdef["g*0.05"] = { 0, 0, 0, 0, 0, 0,0.05, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "dv"] = { 0, 0, 0, 0, 0,-1, 0, 1, 0, 0 ,0, 0 ,0 };
   pdfdef[   "uv"] = { 0, 0, 0, 0,-1, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[   "s+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "s-"] = { 0, 0, 0,-1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c-"] = { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   pdfdef[ "Ubar"] = { 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[ "Dbar"] = { 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "U"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,1, 0 ,1 };
   pdfdef[    "D"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 ,0, 1 ,0 };
   pdfdef[   "U+"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,1, 0 ,1 };
   pdfdef[   "D+"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 ,0, 1 ,0 };
   //
   pdfdef["gluon"] =   {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}; // gluon
   pdfdef["SIGMA"] =   {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}; // SIGMA = \sum_q q^+
   pdfdef["VALENCE"]=  { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}; // VALENCE = \sum_q q^-
   pdfdef[   "T3"] =   {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}; // T3  = u^+ - d^+
   pdfdef[   "V3"] =   {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}; // V3  = u^- - d^-
   pdfdef[   "T8"] =   {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}; // T8  = u^+ + d^+ - 2 s^+
   pdfdef[   "V8"] =   {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}; // V8  = u^- + d^- - 2 s^-
   pdfdef[  "T15"] =   {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}; // T15 = u^+ + d^+ + s^+ - 3 c^+
   pdfdef[  "V15"] =   {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}; // V15 = u^- + d^- + s^- - 3 c^-
   pdfdef[  "T24"] =   {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}; // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+
   pdfdef[  "V24"] =   {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}; // V24 = u^- + d^- + s^- + c^- - 4 b^-
   pdfdef[  "T35"] =   { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}; // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
   pdfdef[  "V35"] =   {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,};  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-

   return pdfdef;

}
