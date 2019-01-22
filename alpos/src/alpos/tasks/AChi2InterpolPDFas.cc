#include "alpos/tasks/AChi2InterpolPDFas.h"
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <alpos/ASubsetFunction.h>
#include "alpos/AFactory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/tasks/AFitter.h"

/*!
 AChi2InterpolPDFas

 This task computes $\chi^2$ using a number of different PDFs with different
 assumptions for $\alpha_s(M_Z)$ and interpolates between the points using
 the refined interpolation described in [1]. It then performs a minimization
 of this function to provide a best estimate for $\alpha_s(M_Z)$ and its
 the uncertainty.

 Each of the PDFs used is characterized by the value of $\alpha_s(M_Z)$
 which was assumed when fitting the PDF. In the following, this is
 referred to as $\alpha_s^\mathrm{(PDF)}(M_Z)$. This task is then useful for
 studying the influence the choice of PDF would have on the $\chi^2$ value,
 and hence on the agreement between data and theory.

 A Table containing the names of the PDF set and the respective member of
 the PDF set to be used for each $\chi^2$ point must be given in the steering
 file. The task output is a ROOT file containing a TGraphErrors of the $\chi^2$
 at the requested values of $\alpha_s^\mathrm{(PDF)}(M_Z)$.

 The $\chi^2$ at a value of $\alpha_s(M_Z)$ that lies between two support points
 is calculated by evolving $\alpha_s$. A separate evolution is done from each of
 the two support points with each evolution assuming the $\alpha_s(M_Z)$ value
 of the corresponding point, thus yielding two estimates for $\chi^2$.
 The final $\chi^2$ value is calculated by a simple linear interpolation
 between these estimates.

 For calculating the $\chi^2$ at a value of $\alpha_s(M_Z)$ that lies outside
 the point range, the $\alpha_s$ evolution is done only once, assuming a
 value of $\alpha_s(M_Z)$ that corresponds to the nearest outermost point.
 For more information on the interpolation method, see [1].

 NOTE: This task should only be used with function AfastNLO // TODO: check this

 //TODO: is there a public link to the DIS15 proceedings available?
 [1] https://ekptrac.physik.uni-karlsruhe.de/trac/fastNLO/browser/trunk/doc/DIS15-proceedings/dis2015_stober_arxiv.pdf

 */

using namespace std;

const string AChi2InterpolPDFas::fTaskType = "Chi2InterpolPDFas";

//____________________________________________________________________________________ //
//AChi2InterpolPDFas::AChi2InterpolPDFas(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AChi2InterpolPDFas::AChi2InterpolPDFas(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AChi2InterpolPDFas");
   //! Important: create always new result-object here!
   fResult = new AChi2InterpolPDFasResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AChi2InterpolPDFas::~AChi2InterpolPDFas(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool AChi2InterpolPDFas::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   debug["Init"]<<"Initializing Task '"<<fTaskName<<"' of type '"<<fTaskType<<"'."<<endl;

   // specify chi2 definition
   fChisqdef = STRING_NS(Chisq,NS());
   debug["Init"]<<"Using chisqdef  '"<<fChisqdef<<"'"<<endl;

   // specify PDF function
   fPdfFunction = STRING_NS(PDF,NS());
   debug["Init"]<<"Using pdfFunction  '"<<fPdfFunction<<"'"<<endl;

   // get parameter points, pdf filesets and members
   fAsmzValues = DOUBLE_COL_NS(ParameterPDFTable,AlphasMz,NS());
   fPdfSets = STRING_COL_NS(ParameterPDFTable,PDFSet,NS());
   fPdfMembers = INT_COL_NS(ParameterPDFTable,PDFMember,NS());

   fAsmzChi2Pairs.resize(fAsmzValues.size());

   return true;
}


//____________________________________________________________________________________ //
bool AChi2InterpolPDFas::Execute(){
   debug["Execute"]<<"Executing Task '"<<fTaskName<<"' of type '"<<fTaskType<<"'."<<endl;

   debug["Execute"]<<"Updating pdfFunction..."<<endl;
   PAR_ANY_S(fPdfFunction);  // update LHAPDF

   int nPoints = GetNPoints();

   debug["Execute"]<<"No. of alpha_s(M_Z) points: "<<nPoints<<endl;
   for ( int iParPoint = 0 ; iParPoint < nPoints ; iParPoint++ ) {
      debug["Execute"]<<"    Point #"<<iParPoint<<", val="<<fAsmzValues[iParPoint]<<", pdfset="<<fPdfSets[iParPoint]<<", pdfmember="<<fPdfMembers[iParPoint]<<endl;
   }

   // (re)create ROOT file
   debug["Execute"]<<"Creating ROOT file '"<<STRING_NS(RootFilename,NS())<<"'"<<endl;
   TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");

   // initialize TGraph with (alphas_mz, chi2) points
   TGraph gChi2Val;
   gChi2Val.SetName(string("alpha_s(M_Z)").c_str());//Form("%s_Chi2",ScanPar.c_str()));
   gChi2Val.GetXaxis()->SetTitle(string("alpha_s(M_Z)").c_str());
   gChi2Val.GetYaxis()->SetTitle("#chi^{2}");

   // loop through parameter values
   info["Execute"]<<"Scanning chi2."<<endl;
   for ( int ipt = 0 ; ipt < nPoints ; ipt++ ) {
      PAR_ANY_S(fPdfFunction);  // update LHAPDF

      // store parameter value in local variable
      double parValue = fAsmzValues[ipt];
      // set parameter value
      debug["Execute"]<<"Setting PDF scan parameter to "<<parValue<<endl;

      std::string pdfFunctionType = TheoryHandler::Handler()->GetFuncD(fPdfFunction)->GetFunctionName();
      debug["Execute"]<<"PDF is from theory parameter: "<<fPdfFunction<<" which is of type '"<<pdfFunctionType<<"'."<<endl;

      // choose different PDF in LHAPDF before fit
      SET_ANY_S(pdfFunctionType+std::string(".LHAPDFFile"), fPdfSets[ipt], std::string(""));
      SET_ANY(pdfFunctionType+std::string(".PDFSet"), fPdfMembers[ipt], 0);

      debug["Execute"]<<"Updating function "<<fPdfFunction<<" which is of type "<<pdfFunctionType<<endl;
      PAR_ANY_S(fPdfFunction);  // update LHAPDF

      debug["Execute"]<<"Calculating chi2 for PDF scan parameter "<<parValue<<endl;

      // Output information about:  (1) Alpha_s evolution code used
      //                            (2) Value of Alpha_s(M_Z) used

      // get DataTheoryPairs
      const auto &dtps = TheoryHandler::Handler()->GetDataTheoryPairs();
      if ( dtps.size() > 1 ) {
         warn["Update"]<<"Multiple datasets not fully supported by Chi2FitpDFas. Results may be incorrect."<<endl;
      }

      // get the first DataTheoryPairs
      const auto &dtp = TheoryHandler::Handler()->GetDataTheoryPairs().begin(); //TODO: handle more than one dataset (?)


      // determine which alpha_s evolution code is used
      AFuncD* theo = dtp->second.second; // get theoryfunction from first Dataset
      std::string theoFuncName = theo->GetFunctionName();

      // check if function is actually a subset
      if (theoFuncName.find(ASubsetFunction::fFunctionName) != std::string::npos) {
         // use the first (and only!) requirement name (this is the non-subset theory function)
         theoFuncName = theo->GetRequirements()[0];
      }
      fAsRun = TheoryHandler::Handler()->GetFuncD(theoFuncName+std::string(".")+std::string("Alpha_s"));
      string asRunFunctionName = fAsRun->GetFunctionName();
      debug["Execute"]<<"Using alpha_s evolution from: "<<asRunFunctionName<<endl;

      debug["Execute"]<<"Updating function "<<fAsRun->GetAlposName()<<" which is of type "<<fAsRun->GetFunctionName()<<endl;
      PAR_ANY(fAsRun->GetAlposName());  // important: Update Alphas evolution function!

      if ( asRunFunctionName == "CRunDec" ) {
         // Notify CRunDec (if used) that AlphasMz has changed
         debug["Execute"]<<"Updating alpha_s(M_Z) in CRunDec."<<endl;

         SET_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("AlphasMz"),fAsmzValues[ipt],0);

         double asmz = PAR_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("AlphasMz"));
         double mz = PAR_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("Mz"));
         debug["Execute"]<<"      alpha_s(M_Z)        is: "<<asmz<<endl;
         debug["Execute"]<<"              M_Z         is: "<<mz<<endl;
      }
      else if ( asRunFunctionName == "LHAPDF6Alphas" ) {
         // LHAPDF6Alphas does not need to be notified about changes in AlphasMz,
         // since it doesn't take the Alpos Parameter AlphasMz as an input
         // but only uses it to communicate the AlphasMz value used by LHAPDF

         double asmz = PAR_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("AlphasMz"));
         //double mz = PAR_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("Mz"));
         //debug["Execute"]<<"Can't get Alpha_s(M_Z) from '"<<asRunFunctionName<<"': no Get...() methods available."<<endl;
         debug["Execute"]<<"Can't get     M_Z  from '"<<asRunFunctionName<<"': no Get...() method available."<<endl;
         debug["Execute"]<<"  alpha_s(91.1876)        is: "<<asmz<<endl;
      }

	   // calculate chisq
	   double chisq = -1;

      // // don't do fit; just evaluate chi2 by definition
      // debug["Execute"]<<"(DoAlphasFit="<<DoAlphasFit<<") Not doing fit for alpha_s(M_Z) = "<<parValue<<endl;

      // initialize a ChisqBase to calculate chi2
      vector<string> emptyparnames;
      const auto& super = TheoryHandler::Handler()->GetSuperPair();

      AChisqBase* chi = AFactory::ChisqFactory(fChisqdef,emptyparnames,super.first,super.second);
      chisq = chi->DoEval(NULL);  // just evaluate chi2 w/o fitting
      fAsmzChi2Pairs[ipt] = make_pair(parValue,chisq);
      delete chi;  // cleanup

      info["Execute"]<<"chi2 value for alpha_s(M_Z) = "<<parValue<<" is: "<<chisq<<endl;
      gChi2Val.SetPoint(ipt,parValue,chisq);  // fill chi2 point into TGraph

   } // end loop through pdf parameter points

   info["Execute"] << "Interpolating chi2 points."<<endl;

   // plot TGraph with resulting interpolation
   TGraph gChi2ValInterpol;
   bool PlotInterpolation = EXIST_NS(PlotInterpolation,NS()) ? BOOL_NS(PlotInterpolation,NS()) : false;
   if (PlotInterpolation) {
      gChi2ValInterpol.SetName(string("alpha_s(M_Z)_Interpolation").c_str());//Form("%s_Chi2",ScanPar.c_str()));
      gChi2ValInterpol.GetXaxis()->SetTitle(string("alpha_s(M_Z)").c_str());
      gChi2ValInterpol.GetYaxis()->SetTitle("#chi^{2}");
      int nPlotPoints = EXIST_NS(PlotPoints, NS()) ? INT_NS(PlotPoints, NS()) : 50;
      //TODO: remove hard coded padding
      double asmzRangeMin = fAsmzValues[0] - 0.002;
      double asmzRangeMax = fAsmzValues[fAsmzValues.size() - 1] + 0.002;
      for (int iPlot = 0; iPlot < nPlotPoints; iPlot++) {
         double alphas = asmzRangeMin + (asmzRangeMax - asmzRangeMin) * iPlot / nPlotPoints;
         double chi2 = Chi2Interpolation(alphas);
         gChi2ValInterpol.SetPoint(iPlot, alphas, chi2);  // fill chi2 point into TGraph
      }
   }

   ROOT::Fit::Fitter* chi2Fitter = new ROOT::Fit::Fitter();
   chi2Fitter->Config().SetMinimizer("TMinuit");

   // --- Configure fitter
   chi2Fitter->Config().MinimizerOptions().SetPrintLevel(1);
   chi2Fitter->Config().MinimizerOptions().SetTolerance(0.01);
   chi2Fitter->Config().MinimizerOptions().SetStrategy(1);
   chi2Fitter->Config().MinimizerOptions().SetMaxFunctionCalls(300000);
   chi2Fitter->Config().MinimizerOptions().SetMaxIterations(300000);

   // --- set fit parameters
   vector<ROOT::Fit::ParameterSettings> paramSettings(1);
   paramSettings[0].Set("AlphasMz", 0.118, 1.e-6, fAsmzValues[0], fAsmzValues[fAsmzValues.size()-1]); /*name,val,step,lower,upper*/
   //paramSettings[0].SetLimits(fAsmzValues[0], fAsmzValues[fAsmzValues.size()-1]);

   chi2Fitter->Config().SetParamsSettings(paramSettings);

   const AChi2InterpolFCN* chi2FCN = new AChi2InterpolFCN(this);
   chi2Fitter->SetFCN(*chi2FCN);
   bool fitOK = chi2Fitter->FitFCN();

   ROOT::Fit::FitResult chi2FitResult = chi2Fitter->Result();

   // --- set errors (and also set again central values)

   double asmzMin = chi2FitResult.GetParams()[0];
   double asmzMinErr = chi2FitResult.GetErrors()[0];

   info["Execute"] << "Interpolation result is: "<<endl;
   //minimizer->Print();
   chi2FitResult.Print(cout);

   // some debug output
   debug["Execute"] << "Minimum is at alpha_s(m_Z) = "<<asmzMin<<endl;
   debug["Execute"] << "Symmetric error is  "<<asmzMinErr<<endl;

   // write TGraph to file
   if ( !file->GetDirectory("Chi2") ) file->mkdir("Chi2")->cd();
   else file->cd("Chi2");
   gChi2Val.Write();

   if (PlotInterpolation) {
      if (!file->GetDirectory("Chi2")) file->mkdir("Chi2")->cd();
      else file->cd("Chi2");
      gChi2ValInterpol.Write();
   }

   // close rootfile
   file->Close();
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk."<<endl;

   // --- set values in TaskResults
   fResult->fParAtChi2Minimum = asmzMin;
   fResult->fParErrorSymmetric = asmzMinErr;

   return true;
}

//____________________________________________________________________________________ //
double AChi2InterpolPDFas::Chi2Interpolation(double alphas_mz) const {
   // interpolates chi2(alphas_mz) between specified (alphas_mz_PDF, chi2) points
   // NOTE: this uses Chi2Extrapolation to "evolve" chi2 from the two nearest
   //       (alphas_mz_PDF, chi2) points, and then interpolates between these two values

   //debug["Chi2Interpolation"]<<"Interpolating for alphas_mz = "<<alphas_mz<<endl;
   // find nearest two chi2 points

   // TODO: make sure list is ordered? check if ordered by default
   //debug["Chi2Interpol"]<<"Looking for alphas_mz = "<<alphas_mz<<" in range ["<<fAsmzChi2Pairs.begin()->first<<", "<<fAsmzChi2Pairs.rbegin()->first<<")"<<endl;
   // find (alphas_mz_PDF, chi2) point immediately after alphas_mz
   vector<pair<double,double>>::const_iterator iterAsmzChi2Pairs;
   iterAsmzChi2Pairs = lower_bound(fAsmzChi2Pairs.begin(), fAsmzChi2Pairs.end(), make_pair(alphas_mz, double(-INFINITY)));

   //debug["Chi2Interpol"]<<"Found alphas_mz = "<<alphas_mz<<" in range ["<<prev(iterAsmzChi2Pairs)->first<<", "<<iterAsmzChi2Pairs->first<<")"<<endl;

   if ( iterAsmzChi2Pairs == fAsmzChi2Pairs.begin() ) {
      // if alphas_mz below PDF asmz range -> "downward" extrapolation
      debug["Chi2Interpolation"]<<"alphas_mz ("<<alphas_mz<<") below PDF asmz range -> \"downward\" extrapolation."<<endl;
      int iPdf = 0;
      return Chi2Extrapolation(alphas_mz, iPdf);
   }
   else if ( iterAsmzChi2Pairs == fAsmzChi2Pairs.end() ) {
      // if alphas_mz above PDF asmz range -> "upward" extrapolation
      debug["Chi2Interpolation"]<<"alphas_mz ("<<alphas_mz<<") above PDF asmz range -> \"upward\" extrapolation."<<endl;
      int iPdf = fAsmzValues.size()-1;
      return Chi2Extrapolation(alphas_mz, iPdf);
   }
   else {
      // interpolation
      debug["Chi2Interpolation"]<<"alphas_mz ("<<alphas_mz<<") in PDF asmz range -> interpolation."<<endl;
      // find nearest lower and higher asmz values
      double alphas_mz_A = prev(iterAsmzChi2Pairs)->first;
      double alphas_mz_B = (iterAsmzChi2Pairs)->first;
      //double chisq_A = prev(iterAsmzChi2Pairs)->second;
      //double chisq_B = (iterAsmzChi2Pairs)->second;

      // get indices from iterator
      int iPdfB = distance(fAsmzChi2Pairs.begin(), iterAsmzChi2Pairs);
      int iPdfA = distance(fAsmzChi2Pairs.begin(), prev(iterAsmzChi2Pairs)); // equivalently: iPdfA = iPdfB - 1
      debug["Chi2Interpolation"]<<alphas_mz_A<<" < alphas_mz < "<<alphas_mz_B<<endl;
      debug["Chi2Interpolation"]<<iPdfA<<" = iPdfA, iPdfB = "<<iPdfB<<endl;
      // get chisq values to interpolate between
      double chisq_A = Chi2Extrapolation(alphas_mz, iPdfA);
      double chisq_B = Chi2Extrapolation(alphas_mz, iPdfB);

      //debug["Chi2Interpolation"]<<"Interpolating points ("<<alphas_mz_A<<", "<<chisq_A<<") and ("<<alphas_mz_B<<", "<<chisq_B<<")"<<endl;
      //debug["Chi2Interpolation"]<<"Interpolated point is ("<<alphas_mz<<", "<<(chisq_A + (chisq_B - chisq_A) * (alphas_mz - alphas_mz_A)/(alphas_mz_B - alphas_mz_A))<<")"<<endl;

      // simple linear interpolation
      return (chisq_A + (chisq_B - chisq_A) * (alphas_mz - alphas_mz_A)/(alphas_mz_B - alphas_mz_A));
   }
}

//____________________________________________________________________________________ //
double AChi2InterpolPDFas::Chi2Extrapolation(double alphas_mz, int reference_pdf_member) const {
   // extrapolates chi2(alphas_mz) from a specified (alphas_mz, chi2) point
   // NOTE: this sets AsRun.AlphasMz and calculates chi2 directly

   debug["Chi2Extrapolation"]<<"Extrapolating chi2 from point #"<<reference_pdf_member<<endl;
   PAR_ANY(fAsRun->GetAlposName());  // important: Update Alphas evolution function!

   // choose different PDF in LHAPDF before fit
   // TODO: works only for CRunDec -> adapt to work with LHAPDF6Alphas (?)
   string pdfFunctionType = TheoryHandler::Handler()->GetFuncD(fPdfFunction)->GetFunctionName();
   SET_ANY_S(pdfFunctionType+std::string(".LHAPDFFile"), fPdfSets[reference_pdf_member], std::string(""));
   SET_ANY(pdfFunctionType+std::string(".PDFSet"), fPdfMembers[reference_pdf_member], 0);

   // set AsRun.AlphasMz to requested value
   SET_ANY(fAsRun->GetAlposName()+std::string(".")+std::string("AlphasMz"), alphas_mz,0);

   vector<string> emptyparnames;
   const auto& super = TheoryHandler::Handler()->GetSuperPair();

   AChisqBase* chi = AFactory::ChisqFactory(fChisqdef,emptyparnames,super.first,super.second);
   double chisq = chi->DoEval(NULL);  // just evaluate chi2 w/o fitting

   // calculate correction factor, so that interpolation function goes through the support points
   // TODO: need chi2 with "native" LHAPDF alpha_s evolution instead of AsRun for correction
   // NOTE: currently, as only AsRun is used by fastNLOAlpos, no correction is necessary
   double corr_factor = 1; // correction factor cannot be calculated yet, assume 1

   return corr_factor*chisq;

}


//// AChi2InterpolFCN - FCN to minimize for best alpha_s(M_Z) estimate

//____________________________________________________________________________________ //
double AChi2InterpolFCN::DoEval(const double *param) const {
   double ret = fChi2InterpolTask->Chi2Interpolation(param[0]);
   cout<<"DoEval: FCN(alphas_mz="<<param[0]<<") = "<<ret<<endl;
   return ret;
}

//____________________________________________________________________________________ //
ROOT::Math::IMultiGenFunction* AChi2InterpolFCN::Clone() const {
   return new AChi2InterpolFCN(*this);
}