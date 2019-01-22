#include "alpos/tasks/AChi2Scan.h"
#include <iostream>
#include <TFile.h>
#include <TGraph.h>
#include <TAxis.h>
#include "alpos/AFactory.h"
#include "alpos/AChisq.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AFitter.h"

/* 
 AChi2Scan

 Add a docu here

 */

using namespace std;

const string AChi2Scan::fTaskType = "Chi2Scan";

//____________________________________________________________________________________ //
//AChi2Scan::AChi2Scan(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AChi2Scan::AChi2Scan(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AChi2Scan");
   //! Important: create always new result-object here!
   fResult = new AChi2ScanResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AChi2Scan::~AChi2Scan(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool AChi2Scan::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   info<<"Hello. AChi2Scan::Init()."<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool AChi2Scan::Execute(){
   debug["Execute"]<<endl;
   info["Execute"]<<"RootFilename   "<<STRING_NS(RootFilename,NS())<<endl;
   double range = DOUBLE_NS(SigmaRange,NS());
   int nPoints  = INT_NS(nPoints,NS());
   info["Execute"]<<"Range:     +/- "<<range<<" sigma"<<endl;
   info["Execute"]<<"nPoints:   "<<nPoints<<""<<endl;
   string chisqdef = STRING_NS(Chisq,NS());
   cout<<endl;
   bool DoFit = BOOL_NS(DoFit,NS());
   string fitterNS = STRING_NS(Fitter,NS());
   vector<string> FitPars = STRING_ARR_NS(FitParameters,fitterNS);
   vector<string> ScanPars = STRING_ARR_NS(FitParameters,NS());

   // init output
   TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   
   // get chisq min
   vector<string> emptyparnames;
   const auto& super = TheoryHandler::Handler()->GetSuperPair();
   AChisqBase* chi = AFactory::ChisqFactory(chisqdef,emptyparnames,super.first,super.second);
   double chisq0 = chi->DoEval(NULL);
   info["Execute"]<<"Chisq0: "<<chisq0<<endl;

   map<string,double> val0;
   map<string,double> err0;
   for ( auto ipar : FitPars) {
      val0[ipar] = TheoryHandler::Handler()->GetParmD(ipar)->GetValue();
      err0[ipar] = TheoryHandler::Handler()->GetParmD(ipar)->GetError();
   }

   // --- loop over parameters
   for ( auto ipar : ScanPars) {
      info["Execute"]<<"Scanning parameter: "<<ipar<<endl;
      const double val = val0[ipar];//TheoryHandler::Handler()->GetParmD(ipar)->GetValue();
      const double err = err0[ipar];//TheoryHandler::Handler()->GetParmD(ipar)->GetError();
      if ( err == 0 ) {
	 error["Execute"]<<"Cannot determine 'range' because error is zero."<<endl;
	 exit(13);
      }
      range*=2 ;// symmetric on both sides
      double step = 2*err/nPoints*range;
      double OneSig = err/val;
      debug["Execute"]<<"val="<<val<<"\terr="<<err<<"\tstep="<<step<<endl;
      info["Execute"]<<"val="<<val<<"\terr="<<err<<"\tstep="<<step<<endl;
      TGraph gChi2Val;      // gChi2Val.SetName(ipar.c_str());//Form("%s_Chi2",ipar.c_str()));
      TGraph gChi2Sig;
      TGraph gChi2MinVal;
      TGraph gChi2MinSig;

      for ( int ipt = 0 ; ipt < nPoints+1 ; ipt++ ) {
	 double vt = val-range*err+ ipt*step;
	 double sig = ipt*range/nPoints-range/2.;
	 //cout<<"vt="<<vt<<" \tsig="<<sig<<endl;
	 //reset all for next parameter scan (if fit failed)
	 for ( auto xpar : FitPars)  {
	    info["Execute"]<<"Reseting parameter: "<<xpar<<"\to "<<val0[xpar]<<" +/- "<<0 /*err0[xpar]*/<<endl;
	    SET_ANY(xpar,val0[xpar],0);
	 }
	 SET_ANY(ipar,vt,0);
	 
	 // calculate chisq
	 double chisq = -1;
	 if ( !DoFit ) {
	    vector<string> emptyparnames;
	    const auto& super = TheoryHandler::Handler()->GetSuperPair();
	    AChisqBase* chi = AFactory::ChisqFactory(chisqdef,emptyparnames,super.first,super.second);
	    chisq = chi->DoEval(NULL);
	    int ndf = chi->Data()->N();
	    cout<<"chisq fix   = "<<chisq<<endl;
	 }
	 else {
	    AFitter fitter(fitterNS);
	    fitter.Init();
	    // --- fix parameter
	    vector<ROOT::Fit::ParameterSettings>& parsettings = fitter.GetFitter()->Config().ParamsSettings();
	    for ( auto& iset : parsettings ) {
	       if ( iset.Name() == ipar )  {
		  info["Execute"]<<"Fixing parameter '"<<ipar<<"'."<<endl;
		  iset.Fix();
	       }
	    }
	    fitter.Execute();
	    chisq = fitter.GetChisq();
	    cout<<"chisq refit = "<<chisq<<endl;
	    //delete fitter.fResults;
	 }
	 gChi2Val.SetPoint(ipt,vt,chisq);
	 gChi2Sig.SetPoint(ipt,sig,chisq);
	 gChi2MinVal.SetPoint(ipt,vt,chisq-chisq0);
	 gChi2MinSig.SetPoint(ipt,sig,chisq-chisq0);
      }

      //reset all for next parameter scan (if fit failed)
      for ( auto xpar : FitPars) {
	 SET_ANY(xpar,val0[xpar],err0[xpar]);
      }


      // write TGraph to file
      if ( !file->GetDirectory("Chi2") ) file->mkdir("Chi2")->cd();
      else file->cd("Chi2");
      gChi2Val.SetName(ipar.c_str());
      gChi2Val.GetXaxis()->SetTitle(ipar.c_str());
      gChi2Val.GetYaxis()->SetTitle("#chi^{2}");
      //gChi2Val.SetTitle(Form("%s;%s;%s",ipar.c_str(),ipar.c_str(),"#chi^{2}"));
      gChi2Val.Write();


      if ( !file->GetDirectory("Chi2_vs_Sigma") ) file->mkdir("Chi2_vs_Sigma")->cd();
      else file->cd("Chi2_vs_Sigma");
      gChi2Sig.SetName(ipar.c_str());//Form("%s_Chi2_Sigma",ipar.c_str()));
      gChi2Sig.GetXaxis()->SetTitle(Form("#sigma(%s)",ipar.c_str()));
      gChi2Sig.GetYaxis()->SetTitle("#chi^{2}");
      gChi2Sig.Write();

      if ( !file->GetDirectory("Chi2Min") ) file->mkdir("Chi2Min")->cd();
      else file->cd("Chi2Min");
      gChi2MinVal.SetName(ipar.c_str());//Form("%s_Chi2Min",ipar.c_str()));
      gChi2MinVal.GetXaxis()->SetTitle(ipar.c_str());
      gChi2MinVal.GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{min}");
      gChi2MinVal.Write();

      if ( !file->GetDirectory("Chi2Min_vs_Sigma") ) file->mkdir("Chi2Min_vs_Sigma")->cd();
      else file->cd("Chi2Min_vs_Sigma");
      gChi2MinSig.SetName(ipar.c_str());//Form("%s_Chi2Min_Sigma",ipar.c_str()));
      gChi2MinSig.GetXaxis()->SetTitle(Form("#sigma(%s)",ipar.c_str()));
      gChi2MinSig.GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{min}");
      gChi2MinSig.Write();


      file->Write();
   }


   // close rootfile
   file->Close();
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk.\n"<<endl;

   return true;

}


//____________________________________________________________________________________ //
