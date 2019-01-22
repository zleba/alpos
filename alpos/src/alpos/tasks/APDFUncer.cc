#include "alpos/tasks/APDFUncer.h"
#include <iostream>
#include <algorithm>
#include <TFile.h>
#include <TH1F.h>
#include "alpos/AFactory.h"
#include "alpos/functions/ALhapdf6.h"
/** 
 APDFUncer

 Add a docu here

 */

using namespace std;

const string APDFUncer::fTaskType = "PDFUncertainty";

//____________________________________________________________________________________ //
APDFUncer::APDFUncer(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("PDFUncer");
   //! Important: create always new result-object here!
   fResult = new APDFUncerResult(aname,fTaskType);
}


//____________________________________________________________________________________ //
APDFUncer::~APDFUncer(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool APDFUncer::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   

   return true;
}



//____________________________________________________________________________________ //
bool APDFUncer::Execute(){
   debug["Execute"]<<"Now getting 'WelcomeString' from steering and printing it:"<<endl;

   const string fittername = STRING_NS(fittername,NS());
   const string fittertype = STRING_NS(fittertype,NS());
   const vector<string> fitpar = STRING_ARR_NS(FitParameters,fittername);
   const string pdffunc = STRING_NS(PDF,NS());
   info["Execute"]<<"Evaluating PDF from theory parameter: "<<pdffunc<<" which is of type '"<<
	   TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName()<<"'."<<endl;

   bool LHAType=TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName().find(ALhapdf6::fFunctionName) !=std::string::npos;
   if ( LHAType ) {
	   ALhapdf6* lha = (ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc);
	   int nPDF = lha->GetPDFSet()->size();
	   info["Execute"]<<"Found "<<nPDF<<" PDF members."<<endl;
	   // loop over all PDF members
	   map<string,vector<double> > vE; //<parname,errors>
	   for ( auto ip : fitpar ) vE[ip].resize(nPDF);
	   for ( int iMem = 0 ;iMem<nPDF ; iMem++ ) {
		   SET_ANY(pdffunc+".PDFSet",iMem,0);
		   ATask* fitter = AFactory::TaskFactory(fittername,fittertype);
		   fitter->Init();
		   fitter->Execute();
		   for ( auto ip : fitpar ) vE[ip][iMem] = PAR_ANY(ip);
		   //delete fitter->GetResult();
		   delete fitter;
	   }
	   // get errors and print them
	   cout<<endl;
	   info["Execute"]<<"Error type of PDF set '"<<lha->GetPDFSet()->name()<<"' is '"<<lha->GetPDFSet()->errorType()<<"' at a confidence level of "<<lha->GetPDFSet()->errorConfLevel()<<"%."<<endl;
	   for ( auto ip : fitpar ) {
		   const vector<double>& values = vE[ip];
		   LHAPDF::PDFUncertainty PDFUnc = lha->GetPDFSet()->uncertainty(values);
		   info["Execute"]<<"PDF uncertainty. Parameter ["<<ip<<"]:"<<endl;
		   info["Execute"]<<"                 central: "<<PDFUnc.central<<" ("<<vE[ip][0]<<"),\terrplus: " <<PDFUnc.errplus<<"\terrminus: "<<PDFUnc.errminus<<endl;
	   }
   }
   else {
	   warn["Execute"]<<"Only LHAPDF style PDFs are currently implemneted."<<endl;
   }
   cout<<endl;
   info["Execute"]<<"Task done. Resetting theory to initial values."<<endl;

   // reset theory to inital values
   TheoryHandler::Handler()->SetTheorySet(this->GetResult()->GetInitialTheorySet());

   return true;
   /*
   double fpdfstart = DOUBLE_NS(PDFSetStart,NS());
   double fpdfend = DOUBLE_NS(PDFSetEnd,NS());
   string fparname = STRING_NS(ParName,NS());
   double fparstart = DOUBLE_NS(ParStart,NS());
   double fparend = DOUBLE_NS(ParEnd,NS());
   string fchisqdef = STRING_NS(Chisq,NS());
   int fnbins = INT_NS(NumberOfBins,NS());

   cout<<endl;
   cout<<"  PDFSetStart:     "<<fpdfstart<<endl;
   cout<<"  PDFSetEnd:       "<<fpdfend<<endl;
   cout<<endl;

   // init output
   TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   
   const auto& super = TheoryHandler::Handler()->GetSuperPair();


   vector<string> emptyfparnames;
   string chisqdef = STRING_NS(Chisq,NS());
   AChisqBase* chi = AFactory::ChisqFactory(chisqdef,emptyfparnames,super.first,super.second);
   const int fval0 = PAR_ANY("LHAPDF6.PDFSet");

   TH1F *histogram = new TH1F("histogram", "PDF Uncertainty histogram", fnbins, fparstart, fparend);
   
   for (int i=fpdfstart; i<fpdfend+1; i++){
      // Used to not include certain PDF sets. Used when one gets error from a certain PDF set.
      //if ( i == 45 || i == 60 || i == 83 || i == 93 ) continue;

      SET_ANY("LHAPDF6.PDFSet",i,0);

      ATask* task = AFactory::TaskFactory("MyFit","AFitter");
      task->Init();
      task->Execute();

      histogram->Fill(PAR_ANY(fparname));
      delete task;
   }

   SET_ANY("LHAPDF6.PDFSet",fval0,0);

   histogram->Write();   

   // close rootfile
   file->Close();
   info["Execute"]<<"File "<<file->GetName()<<" written to disk."<<endl;


   // --- job done successfully
   return true;
   */
}


//____________________________________________________________________________________ //
