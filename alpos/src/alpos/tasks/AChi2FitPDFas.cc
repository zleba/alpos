#include "alpos/tasks/AChi2FitPDFas.h"
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <alpos/ASubsetFunction.h>
#include <alpos/functions/ACRunDecFunction.h>
#include <alpos/functions/ALhapdf6Alphas.h>
#include "alpos/AFactory.h"
#include "alpos/AChisq.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AFitter.h"

/*!
 AChi2FitPDFas

 This task computes $\chi^2$ using a number of different PDFs with different
 assumptions for $\alpha_s(M_Z)$ and fits a user-specified function to
 the resulting points. It then determines the minimum of this function
 and estimates the uncertainty of the minimum via the
 $\chi^2_\mathrm{min} + 1$ prescription.

 Each of the PDFs used is characterized by the value of $\alpha_s(M_Z)$
 which was assumed when fitting the PDF. In the following, this is
 referred to as $\alpha_s^\mathrm{(PDF)}(M_Z)$. This task is then useful for
 studying the influence the choice of PDF would have on the $\chi^2$ value,
 and hence on the agreement between data and theory.

 A Table containing the names of the PDF set and the respective member of
 the PDF set to be used for each $\chi^2$ point must be given in the steering
 file. The task output is a ROOT file containing a TGraphErrors of the $\chi^2$
 at the requested values of $\alpha_s^\mathrm{(PDF)}(M_Z)$.

 Optionally, a user-specified function (e.g. 'pol2') may be fitted to the resulting
 $\chi^2$ points. This fit is then typically used to estimate the parameter
 error by determining the values of $\alpha_s(M_Z)$ for which $\chi^2$ is increased
 by one with respect to the minimum.
 This is the procedure applied in [arXiv:1410.6765].

 NOTE: This task should only be used with function AfastNLO // TODO: check this

 */

using namespace std;

const string AChi2FitPDFas::fTaskType = "Chi2FitPDFas";


//____________________________________________________________________________________ //
//AChi2FitPDFas::AChi2FitPDFas(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AChi2FitPDFas::AChi2FitPDFas(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AChi2FitPDFas");
   //! Important: create always new result-object here!
   fResult = new AChi2FitPDFasResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AChi2FitPDFas::~AChi2FitPDFas(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool AChi2FitPDFas::Init(){
   return true;
}


//____________________________________________________________________________________ //
bool AChi2FitPDFas::Execute(){
   debug["Execute"]<<"Executing Task '"<<fTaskName<<"' of type '"<<fTaskType<<"'."<<endl;

   // specify chi2 definition
   string chisqdef = STRING_NS(Chisq,NS());
   info["Execute"]<<"Using chisqdef  '"<<chisqdef<<"'"<<endl;

   // specify PDF function
   string pdfFunction = STRING_NS(PDF,NS());
   info["Execute"]<<"Using pdfFunction  '"<<pdfFunction<<"'"<<endl;
   debug["Execute"]<<"Updating pdfFunction..."<<endl;
   PAR_ANY_S(pdfFunction);  // update LHAPDF

   // get parameter points, pdf filesets and members
   vector<double> asmzValues = DOUBLE_COL_NS(ParameterPDFTable,AlphasMz,NS());
   vector<string> pdfSets = STRING_COL_NS(ParameterPDFTable,PDFSet,NS());
   vector<int> pdfMembers = INT_COL_NS(ParameterPDFTable,PDFMember,NS());

   int nPoints  = asmzValues.size();

   debug["Execute"]<<"No. of alpha_s(M_Z) points: "<<nPoints<<endl;
   for ( int iParPoint = 0 ; iParPoint < nPoints ; iParPoint++ ) {
      debug["Execute"]<<"    Point #"<<iParPoint<<", val="<<asmzValues[iParPoint]<<", pdfset="<<pdfSets[iParPoint]<<", pdfmember="<<pdfMembers[iParPoint]<<endl;
   }

   // get interpolation function
   string fitFuncName = STRING_NS(FitFunc,NS());

   // (re)create ROOT file
   debug["Execute"]<<"Creating ROOT file '"<<STRING_NS(RootFilename,NS())<<"'"<<endl;
   TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");

   // initialize TGraph with (ScanParameter, chi2) points
   TGraph gChi2Val;
   gChi2Val.SetName(string("alpha_s(M_Z)").c_str());//Form("%s_Chi2",ScanPar.c_str()));
   gChi2Val.GetXaxis()->SetTitle(string("alpha_s(M_Z)").c_str());
   gChi2Val.GetYaxis()->SetTitle("#chi^{2}");

   // get PDF function type
   std::string pdfFunctionType = TheoryHandler::Handler()->GetFuncD(pdfFunction)->GetFunctionName();

   // retain initial values of PDF function parameters
   const unsigned int initialPDFMember = PAR_ANY(pdfFunction+std::string(".PDFSet"));
   const std::string initialPDFFile = PAR_ANY_S(pdfFunction+std::string(".LHAPDFFile"));

   // loop through parameter values
   info["Execute"]<<"Scanning chi2."<<endl;
   for ( int ipt = 0 ; ipt < nPoints ; ipt++ ) {
      PAR_ANY_S(pdfFunction);  // update LHAPDF

      // store parameter value in local variable
      double parValue = asmzValues[ipt];
      // set parameter value
      debug["Execute"]<<"Setting PDF scan parameter to "<<parValue<<endl;
      debug["Execute"]<<"PDF is from Alpos parameter '"<<pdfFunction<<"', which is of type '"<<pdfFunctionType<<"'."<<endl;

      // choose different PDF in LHAPDF before fit
      SET_ANY_S(pdfFunction+std::string(".LHAPDFFile"), pdfSets[ipt], std::string(""));
      SET_ANY(pdfFunction+std::string(".PDFSet"), pdfMembers[ipt], 0);

      debug["Execute"]<<"Updating function "<<pdfFunction<<", which is of type "<<pdfFunctionType<<endl;
      PAR_ANY_S(pdfFunction);  // update LHAPDF

      debug["Execute"]<<"Calculating chi2 for PDF scan parameter "<<parValue<<endl;

      // Output information about:  (1) Alpha_s evolution code used
      //                            (2) Value of Alpha_s(M_Z) used

      // get DataTheoryPairs
      const auto &dtps = TheoryHandler::Handler()->GetDataTheoryPairs();
      if ( dtps.size() > 1 ) {
         warn["Update"]<<"Multiple datasets not fully supported by Chi2FitPDFas. Results may be incorrect."<<endl;
      }

      // get the first DataTheoryPairs
      const auto &dtp = TheoryHandler::Handler()->GetDataTheoryPairs().begin();

      //AFuncD* theo = (*(dtps.begin())).second.second; // get theoryfunction from first Dataset
      AFuncD* theo = dtp->second.second; // get theoryfunction from first Dataset
      std::string theoFuncName = theo->GetFunctionName();

      // check if function is actually a subset
      if (theoFuncName.find(ASubsetFunction::fFunctionName) != std::string::npos) {
         // use the first (and only!) requirement name (this is the non-subset theory function)
         theoFuncName = theo->GetRequirements()[0];
      }
      AFuncD* asRun = TheoryHandler::Handler()->GetFuncD(theoFuncName+std::string(".")+std::string("Alpha_s"));
      string asRunFunctionName = asRun->GetFunctionName();
      debug["Execute"]<<"Using alpha_s evolution from: "<<asRunFunctionName<<endl;

      debug["Execute"]<<"Updating function "<<asRun->GetAlposName()<<" which is of type "<<asRun->GetFunctionName()<<endl;
      PAR_ANY(asRun->GetAlposName());  // important: Update Alphas evolution function!

      if ( asRunFunctionName == ACRunDecFunction::fFunctionName ) {
         // Notify CRunDec (if used) that AlphasMz has changed
         debug["Execute"]<<"Updating alpha_s(M_Z) in CRunDec."<<endl;
         SET_ANY(asRun->GetAlposName()+std::string(".")+std::string("AlphasMz"),asmzValues[ipt],0);

         double asmz = PAR_ANY(asRun->GetAlposName()+std::string(".")+std::string("AlphasMz"));
         double mz = PAR_ANY(asRun->GetAlposName()+std::string(".")+std::string("Mz"));
         debug["Execute"]<<"      alpha_s(M_Z)        is: "<<asmz<<endl;
         debug["Execute"]<<"              M_Z         is: "<<mz<<endl;
      }
      else if ( asRunFunctionName == ALhapdf6Alphas::fFunctionName ) {
         // LHAPDF6Alphas does not need to be notified about changes in AlphasMz,
         // since it doesn't take the Alpos Parameter AlphasMz as an input
         // but only uses it to communicate the AlphasMz value used by LHAPDF

         double asmz = PAR_ANY(asRun->GetAlposName()+std::string(".")+std::string("AlphasMz"));
         //double mz = PAR_ANY(asRun->GetAlposName()+std::string(".")+std::string("Mz"));
         //debug["Execute"]<<"Can't get Alpha_s(M_Z) from '"<<asRunFunctionName<<"': no Get...() methods available."<<endl;
         debug["Execute"]<<"Can't get     M_Z  from '"<<asRunFunctionName<<"': no Get...() method available."<<endl;
         debug["Execute"]<<"    alpha_s(?M_Z?)        is: "<<asmz<<endl;
      }
      else {
         error["Execute"]<<"Alpha_s from function '"<<asRunFunctionName<<"' not supported: must be either 'CRunDec' or 'LHAPDF6Alphas."<<endl;
      }

	   // calculate chisq
	   double chisq = -1;

      // initialize a ChisqBase to calculate chi2
      vector<string> emptyparnames;
      const auto& super = TheoryHandler::Handler()->GetSuperPair();

      AChisqBase* chi = AFactory::ChisqFactory(chisqdef,emptyparnames,super.first,super.second);
      chisq = chi->DoEval(NULL);  // just evaluate chi2 w/o fitting

      info["Execute"]<<"chi2 value for alpha_s(M_Z) = "<<parValue<<" is: "<<chisq<<endl;

	   gChi2Val.SetPoint(ipt,parValue,chisq);  // fill chi2 point into TGraph
   } // end loop through pdf parameter points

   info["Execute"] << "Fitting chi2 points. Fit function is: "<<fitFuncName<<endl;

   // do fit and get result
   double parValMin = asmzValues[0];
   double parValMax = asmzValues[nPoints-1];
   TFitResultPtr fitResult = gChi2Val.Fit(fitFuncName.c_str(),"Q", "", parValMin, parValMax);
   TF1 *myfunc = gChi2Val.GetFunction(fitFuncName.c_str());

   info["Execute"] << "Fit result is: "<<endl;
   myfunc->Print();

   double chi2FitFuncMinimum = myfunc->GetMinimum();
   double xAtMinimum = myfunc->GetMinimumX();
   double firstXAtChi2PlusOne = myfunc->GetX(chi2FitFuncMinimum + 1, parValMin, xAtMinimum);
   double secondXAtChi2PlusOne = myfunc->GetX(chi2FitFuncMinimum + 1, xAtMinimum, parValMax);

   double sigmaPlus = secondXAtChi2PlusOne - xAtMinimum;
   double sigmaMinus = firstXAtChi2PlusOne - xAtMinimum;

   double sigmaSymm = (sigmaPlus - sigmaMinus)/2;
   double sigmaSymmParabolic = 1.0 / sqrt(abs(myfunc->GetParameter(2)));

   double sigmaDiff = (sigmaPlus + sigmaMinus)/2;
   double sigmaAsymmetry = sigmaDiff/sigmaSymm;

   // some debug output
   info["Execute"] << "Minimum is at x = "<<xAtMinimum<<endl;
   if (xAtMinimum <= parValMin) {
      warn["Execute"] << "x_min <= " << parValMin << "! (outside point range)" << endl;
   }
   else if (xAtMinimum >= parValMax) {
      warn["Execute"] << "x_min >= " << parValMax << "! (outside point range)" << endl;
   }

   info["Execute"] << "Minimum chi2 = "<<chi2FitFuncMinimum<<endl;
   info["Execute"] << "Searching for 'Delta_chi2 = 1' numerically gives: "<<sigmaPlus<<" (up), "<<sigmaMinus<<" (dn)"<<endl;
   info["Execute"] << "Parabolic error is   +/- "<<sigmaSymmParabolic<<endl;

   if (abs(sigmaAsymmetry) > 1e-4) {
      warn["Execute"] << "Error is asymmetric! "<<sigmaAsymmetry<<endl;
      warn["Execute"] << "Solving 'Delta_chi2 = 1' numerically and symmetrizing gives: +/- "<<sigmaSymm<<endl;
      warn["Execute"] << "Error difference is "<<sigmaDiff<<endl;
   }

   // write TGraph to file
   if ( !file->GetDirectory("Chi2") ) file->mkdir("Chi2")->cd();
   else file->cd("Chi2");
   gChi2Val.Write();

   // close rootfile
   file->Close();
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk."<<endl;


   // --- set values in TaskResults
   fResult->fChi2Minimum = chi2FitFuncMinimum;
   fResult->fParAtChi2Minimum = xAtMinimum;
   fResult->fParErrorSymmetric = sigmaSymm;
   fResult->fParErrorUp = sigmaPlus;
   fResult->fParErrorDn = sigmaMinus;

   // --- restore initial PDF function parameters
   SET_ANY_S(pdfFunction+std::string(".LHAPDFFile"), initialPDFFile, std::string(""));
   SET_ANY(pdfFunction+std::string(".PDFSet"), initialPDFMember, 0);

   // --- set Alpos parameter "AlphasMz" to best-fit value
   const auto &dtp = TheoryHandler::Handler()->GetDataTheoryPairs().begin();
   AFuncD* theo = dtp->second.second;  // get theoryfunction from first Dataset
   string theoFuncName = theo->GetFunctionName();

   // check if function is actually a subset
   if (theoFuncName.find(ASubsetFunction::fFunctionName) != std::string::npos) {
      // use the first (and only!) requirement name (this is the non-subset theory function)
      theoFuncName = theo->GetRequirements()[0];
   }

   AFuncD* asRun = TheoryHandler::Handler()->GetFuncD(theoFuncName+std::string(".Alpha_s"));

   info["Execute"] << "Setting Alpos parameter '" << asRun->GetAlposName()+std::string(".AlphasMz") << "' to best fit value (" << xAtMinimum << ")..." << endl;
   SET_ANY(asRun->GetAlposName()+std::string(".AlphasMz"), xAtMinimum, 0);

   return true;
}


//____________________________________________________________________________________ //
