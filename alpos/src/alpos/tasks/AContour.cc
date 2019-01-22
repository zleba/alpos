#include "alpos/tasks/AContour.h"
#include <iostream>
#include <TDirectory.h>
#include <TGraph.h>
#include <TAxis.h>
#include "alpos/AFactory.h"
#include "alpos/AChisq.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AFitter.h"

#include "Math/Minimizer.h"
#include "TMinuitMinimizer.h"
/* 
 AContour

 Add a docu here

 */

using namespace std;

const string AContour::fTaskType = "Contour";

//____________________________________________________________________________________ //
//AContour::AContour(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AContour::AContour(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AContour");
   //! Important: create always new result-object here!
   fResult = new AContourResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AContour::~AContour(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool AContour::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   info<<"Hello. AContour::Init()."<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool AContour::Execute(){
   debug["Execute"]<<endl;
   //info["Execute"]<<"RootFilename   "<<STRING_NS(RootFilename,NS())<<endl;
   int nPoints  = INT_NS(nPoints,NS());
   string par1 = STRING_NS(par1,NS());
   string par2 = STRING_NS(par2,NS());
   string afitter = STRING_NS(afitter,NS());
   info["Execute"]<<"par1       "<<par1<<endl;
   info["Execute"]<<"par2       "<<par2<<endl;
   info["Execute"]<<"nPoints:   "<<nPoints<<""<<endl;
   string Minimizer = STRING_NS(Minimizer,afitter);
   if ( Minimizer != "TMinuit" ) {
      error["Execute"]<<"Contour works only with TMinuit as minimizer, but 'Minimizer' in '"<<afitter<<"' is '"<<Minimizer<<"'. Exiting."<<endl;
      exit(1);
   }
   vector<string> fitpar = STRING_ARR_NS(FitParameters,afitter);
   int ipar1=-1, ipar2=-1;
   int ii=0;
   for ( auto ipar : fitpar ) {
      if ( ipar == par1 ) ipar1=ii;
      if ( ipar == par2 ) ipar2=ii;
      ii++;
   }
   if ( ipar1==-1 ) {
      error["Execute"]<<"Could not find requester parameter 1 '"<<par1<<"' in list of fitparameters."<<endl;
      return false;
   }
   if ( ipar2==-1 ) {
      error["Execute"]<<"Could not find requester parameter 2 '"<<par2<<"' in list of fitparameters."<<endl;
      return false;
   }


   // get chisq min
   vector<string> emptyparnames;
   const auto& super = TheoryHandler::Handler()->GetSuperPair();
   AFitter* fitter = (AFitter*)AFactory::TaskFactory(afitter,"AFitter");
   fitter->Init();
   fitter->Execute();

   ROOT::Fit::Fitter* rfitter  = fitter->GetFitter();
   ROOT::Math::Minimizer* minim =rfitter->GetMinimizer();
   // todo: Check if minimizer is 'minuit'
   TMinuitMinimizer* minuit = (TMinuitMinimizer*)minim;
   unsigned int pts = nPoints;
   vector<double> xi(pts);
   vector<double> xj(pts);
   double errdef = 2.30;
   if ( EXIST_NS(sig,NS()) ) { 
      double sig=DOUBLE_NS(sig,NS());
      error["Execute"]<<"The steering parameter 'sig' is outdated. Use parameter errdef instead."<<endl;
      error["Execute"]<<"For 2 fit-parameters"<<endl;
      error["Execute"]<<"Use   errdef    CL (inside hypercontour)"<<endl;
      error["Execute"]<<"Use     1              39 %  "<<endl;
      error["Execute"]<<"Use     2.30           68.3 %  "<<endl;
      error["Execute"]<<"Use     2.41           70 %  "<<endl;
      error["Execute"]<<"Use     4.61           90 %  "<<endl;
      error["Execute"]<<"Use     5.99           95 %  "<<endl;
      error["Execute"]<<"Use     9.21           99 %  "<<endl;
      error["Execute"]<<"Details found in minuit manual or book of G. Cowan."<<endl;
      // df:  degrees of freedom for χ2 curve 
      // P:     area under the χ2 curve with df degrees of freedom to the right 
      //  Tail probabilities P
      // df / P: 0.25    0.20    0.15    0.10    0.05    0.025    0.02    0.01    0.005    0.0025    0.001    0.0005 
      //  1      1.32   1.64    2.07    2.71    3.84     5.02    5.41    6.63     7.88      9.14    10.83    12.12
      //  2      2.77   3.22    3.79    4.61    5.99     7.38    7.82    9.21    10.60     11.98    13.82    15.20
      exit(10);
   }
   else if ( EXIST_NS(errdef,NS()) ) {
      errdef = DOUBLE_NS(errdef,NS());
   }
   else {
      warn["Execute"]<<"No parameter 'errdef' specified, using default 2.3"<<endl;
   }
   info["Execute"]<<"errdef     "<<errdef<<endl;
   //minuit->SetErrorDef(pow(sig,2)); // set significance level
   minuit->SetErrorDef(errdef); // set significance level
   minuit->Contour(ipar1,ipar2,pts,&xi[0],&xj[0]);
   minuit->SetErrorDef(1.); // reset

   info["Execute"]<<"contour done!"<<endl;
   for ( int i = 0 ; i<pts ; i++ ) {
      cout<<"Contour: i="<<i<<"\txi="<<xi[i]<<"\txj="<<xj[i]<<endl;
   }
   cout<<"init tgraph."<<endl;
   xi.push_back(xi[0]);
   xj.push_back(xj[0]);
   TGraph g(nPoints+1,&xi[0],&xj[0]);
   g.GetXaxis()->SetTitle(par1.c_str());
   g.GetYaxis()->SetTitle(par2.c_str());
   g.SetName(Form("contour-%s-%s;%s;%s",par1.c_str(),par2.c_str(),par1.c_str(),par2.c_str()));

   // --- output
   TDirectory* file = Alpos::Current()->Settings()->rootoutput;
   //TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   file->mkdir(GetTaskName().c_str())->cd();
   g.Write();
   info["Execute"]<<"TGraph of countour written to TDirectory: "<<GetTaskName()<<endl;
   
   // output as ascii file
   string gDir = Alpos::Current()->Settings()->outputdir;
   //gSystem->mkdir(gDir+"/"+GetTaskName(),true);
   string outfile = gDir+"/"+GetTaskName();
   ofstream out(outfile.c_str());
   out<<"n\t"<<par1.c_str()<<"\t"<<par2.c_str()<<endl;
   for ( int i = 0 ; i<pts ; i++ ) {
      out<<i<<"\t"<<xi[i]<<"\t"<<xj[i]<<endl;
   }
   info["Execute"]<<"Ascii output file written: "<<outfile<<endl;
   
   
   // // close rootfile
   // file->Close();
   file->Write();
   
   //delete fitter->GetResult();
   delete fitter;
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk.\n"<<endl;

   // reset theory to inital values
   TheoryHandler::Handler()->SetTheorySet(this->GetResult()->GetInitialTheorySet());

   return true;

}


//____________________________________________________________________________________ //
