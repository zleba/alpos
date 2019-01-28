#include "alpos/tasks/AFitter.h"
#include "alpos/AFactory.h"
#include <iostream>
#include <alpos/functions/AfastNLOInterpolPDFas.h>
#include <alpos/ASubsetFunction.h>
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "TMinuitMinimizer.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"

/* 
 ATask

 */

using namespace std;

const string AFitter::fTaskType = "AFitter";

//____________________________________________________________________________________ //
//AFitter::AFitter(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AFitter::AFitter(const string& aname) : ATask(aname) {
   //! constructor
   SetClassName("AFitter");

   //! Important: create always new result-object here!
   fResult = new AFitterResult(aname, GetTaskType());
}


//____________________________________________________________________________________ //
AFitter::~AFitter() {
   //! destructor.
   //! do not delete the AResult object
   if (fFitter) delete fFitter;
   if (fChisq) delete fChisq;
}


//____________________________________________________________________________________ //
bool AFitter::Init() {
   // --- setting up minimzer
   debug["Init"] << "Init minimizer for fitter '" << GetTaskName() << "'" << endl;

   fFitter = new ROOT::Fit::Fitter();
   string minimizer = STRING_NS(Minimizer, NS());
   fFitter->Config().SetMinimizer(minimizer.c_str());

   // --- Configure fitter
   int PrintLevel = INT_NS(PrintLevel, NS());
   double Tolerance = DOUBLE_NS(Tolerance, NS());
   int Strategy = INT_NS(Strategy, NS());
   fFitter->Config().MinimizerOptions().SetPrintLevel(PrintLevel);
   fFitter->Config().MinimizerOptions().SetTolerance(Tolerance);
   fFitter->Config().MinimizerOptions().SetStrategy(Strategy);
   fFitter->Config().MinimizerOptions().SetMaxFunctionCalls(300000);
   fFitter->Config().MinimizerOptions().SetMaxIterations(300000);


   // --- get fit parameters from steering

   // exit if no fit parameters
   if (!EXIST_NS(FitParameters, NS())) {
      error["Init"] << "No fit parameters given for fitter '" << GetTaskName() << "'. Exiting." << endl;
      exit(1);
   }

   // look for array named "FitParameters"
   if (EXISTARRAY_NS(FitParameters, NS())) {
      // fit parameters given as simple list of parameter names
      fParNames = STRING_ARR_NS(FitParameters, NS());

      // parameter "metadata"
      fParIsLimited = vector<bool>(fParNames.size(), false);
      fParIsFixed = vector<bool>(fParNames.size(), false);
      fParStepSizes = vector<double>(fParNames.size(), 1e-6);
      fParLowerLimits = vector<double>(fParNames.size(), 0.);
      fParUpperLimits = vector<double>(fParNames.size(), 0.);

      // NOTE: parameter range will be unlimited, initial step size will be 1E-5 * par value or 1E-6 if zero
      for (unsigned int iPar = 0; iPar < fParNames.size(); iPar++) {
         if (PAR_ANY(fParNames[iPar]) != 0)
            fParStepSizes[iPar] = fabs(PAR_ANY(fParNames[iPar])) * 1E-5;
      }
   }
   else {
      // otherwise, assume fit parameters are given as a table
      // NOTE: there is no EXISTTABLE_NS !
      vector<string> columnHeaders = TABLEHEADER_NS(FitParameters, NS());

      // column 'Name' is mandatory
      fParNames = STRING_COL_NS(FitParameters, Name, NS());

      // parameter "metadata"
      fParIsLimited = vector<bool>(fParNames.size(), false);
      fParStepSizes = vector<double>(fParNames.size(), 1e-5);
      fParLowerLimits = vector<double>(fParNames.size(), 0.);
      fParUpperLimits = vector<double>(fParNames.size(), 0.);

      // optional columns
      for (string colName : columnHeaders) {
         if (colName == "Step") {
            fParStepSizes = DOUBLE_COL_NS(FitParameters, Step, NS());
         }
         else if (colName == "LowerLimit") {
            fParLowerLimits = DOUBLE_COL_NS(FitParameters, LowerLimit, NS());
         }
         else if (colName == "UpperLimit") {
            fParUpperLimits = DOUBLE_COL_NS(FitParameters, UpperLimit, NS());
         }
         else {
            warn["Init"] << "Unknown column '" << colName << "' found in 'FitParameters' table. " <<
            "Should be one of 'Name', 'Step', 'LowerLimit', 'UpperLimit'. Ignoring..." << endl;
         }
      }

      for (unsigned int iPar = 0; iPar < fParNames.size(); iPar++) {
         if ((fParLowerLimits[iPar] == 0) && (fParUpperLimits[iPar] == 0))
            fParIsLimited[iPar] = false;
         else if (fParLowerLimits[iPar] >= fParUpperLimits[iPar]) {
            error["Init"] << "Lower limit '" << fParLowerLimits[iPar] <<
            "' is greater than or equal to the upper limit '" << fParUpperLimits[iPar] <<
            "' for parameter #" << iPar << " ('" << fParNames[iPar] << "')" << endl;
         }
         else {
            fParIsLimited[iPar] = true;
         }
      }
   }

   // --- set chisq function
   string chisqdef = STRING_NS(Chisq, NS());
   const auto& super = TheoryHandler::Handler()->GetSuperPair();
   fChisq = AFactory::ChisqFactory(chisqdef, fParNames, super.first, super.second);
   // (AData*)TheoryHandler::Handler()->GetFuncD("SuperData"),
   // TheoryHandler::Handler()->GetFuncD("SuperTheory"));
   if (!fChisq) {
      error["Init"] << "AFitter:Init(). Failed to initialize chisq with name: " << chisqdef << endl;
      exit(1);
   }
   vector<string> allparnames = fChisq->GetFitPar();
   vector<bool> allparfixed = vector<bool>(allparnames.size(), false);

   // handle fixing of parameters (only applies to parameters of interest as declared in steering)
   if (EXISTARRAY_NS(FixedParameters, NS())) {
      // fixed parameters given as simple list of parameter names -> fix to current value
      vector<string> fixedParNames = STRING_ARR_NS(FixedParameters, NS());

      for (unsigned int iFixedPar = 0; iFixedPar < fixedParNames.size(); iFixedPar++) {
         // check if fit parameter exists
         std::vector<string>::iterator fixParAt = std::find(allparnames.begin(), allparnames.end(), fixedParNames[iFixedPar]);
         if (fixParAt != allparnames.end()) {
            // match parameter number and set flag
            unsigned int iPar = fixParAt - allparnames.begin();
            allparfixed[iPar] = true;
            // NOTE: fParIsFixed isn't used, but keep in case we want to use it in the future
            if (iPar < fParIsFixed.size()) {
               fParIsFixed[iPar] = true;
            }
         }
         else {
            // parameter not found -> warn
            warn["Init"] << "Requested fixing of parameter '" << fixedParNames[iFixedPar] << "', but no such parameter found in FitParameters or among nuisance parameters! Ignoring..." << endl;
         }
      }
   }
   // TODO: add possibility to override current value before fix

   // fFitPar.resize(fChisq->NDim());
   // // --- init fitpar with current values
   // // all other parameters are initialized with zero
   // for ( unsigned int i = 0 ; i<fParNames.size() ; i++ )
   //    fFitPar[i] = PAR_ANY(fParNames[i]);
   // fFitter->Config().SetParamsSettings(fFitPar.size(),&fFitPar[0]);

   // --- set fit parameters
   vector<ROOT::Fit::ParameterSettings> parsettings(fChisq->NDim());
   unsigned int i = 0;
   for (auto& iset : parsettings) {
      if (i < fParNames.size()) {
         // main parameters (parameters of interest) may be limited
         if (fParIsLimited[i])
            iset.Set(fParNames[i], PAR_ANY(fParNames[i]), fParStepSizes[i], fParLowerLimits[i],
                     fParUpperLimits[i]); /*name,val,step,lower,upper*/
         else
            iset.Set(fParNames[i], PAR_ANY(fParNames[i]), fParStepSizes[i]); /*name,val,step*/
      }
      else {
         // additional parameters (e.g. nuisance parameters) are not limited
         iset.Set(allparnames[i], 0, 1.e-3); /*name,val,step,lower,upper*/
      }

      // do actual fixing of parameters, if so requested
      if (allparfixed[i]) {
         iset.Fix();
      }

      i++;
   }

   fFitter->Config().SetParamsSettings(parsettings);

   // --- Set chisq function
   fFitter->SetFCN(*fChisq);

   return true;
}


//____________________________________________________________________________________ //
bool AFitter::Execute() {


   // --- Get initial chisq
   fFitter->EvalFCN();
   fInitialChisq = fFitter->Result().MinFcnValue();


   // --- Do fit
   bool fitOK = fFitter->FitFCN();
   if (!fitOK) {
      warn["Execute"] << "Warning! AFitter::Execute(). Fit failed." << endl;
   }

   // --- Get final chisq
   ROOT::Fit::FitResult result = fFitter->Result();
   fFinalChisq = result.MinFcnValue();
   info["Execute"] << "Info. AFitter::Execute(). Initial chisq: " << fInitialChisq << endl;
   info["Execute"] << "Info. AFitter::Execute(). Final Chisq:   " << fFinalChisq << endl;

   // --- set errors (and also set again central values)
   const double* parsptr = fFitter->Result().GetParams();
   const double* errsptr = fFitter->Result().GetErrors();
   
   vector<string> allparnames = fChisq->GetFitPar();
   vector <double> pars(allparnames.size()); 
   vector <double> errs(allparnames.size()); 
   for ( unsigned int ip =0; ip<allparnames.size() ; ip++) {
      pars[ip] = parsptr[ip];
      errs[ip] = errsptr[ip];
   }

   if (EXIST_NS(CalculateMinosErrors, NS()) && BOOL_NS(CalculateMinosErrors, NS())) {
      fFitter->CalculateMinosErrors();
   }

   // --- More printing options
   if (BOOL_NS(PrintResults, NS())) {
      fFitter->Result().Print(cout);
      cout << endl;
      cout << "****************************************" << endl;
      cout << "Printing nuisance parameters:" << endl;
      fChisq->PrintNuisanceParameters(true);
   }

   if (EXIST_NS(CalculateMinosErrors, NS()) && BOOL_NS(CalculateMinosErrors, NS())) {
      std::cout << endl;
      std::cout << "Printing Minos errors:" << endl;
      for (unsigned int iPar = 0; iPar < allparnames.size(); iPar++) {
         if (fFitter->Result().HasMinosError(iPar)) {
            double eup = fFitter->Result().UpperError(iPar);
            double edn = fFitter->Result().LowerError(iPar);
            std::cout << allparnames[iPar] << " = " << pars[iPar] << "  " << std::showpos
                      <<  eup << " (up)  "
                      <<  edn << " (dn)" << std::noshowpos << std::endl;
         }
      }
   }

   if (BOOL_NS(PrintCovariance, NS()))
      fFitter->Result().PrintCovMatrix(cout);

   // --- set values to TaskResults
   fResult->chisq = fFinalChisq;




   // --- if requested, do fits for individual DataTheory pairs
   if (EXIST_NS(FitDataTheoryPairs, NS()) && BOOL_NS(FitDataTheoryPairs, NS())) {
      FitDataTheoryPairs();
   }
   // --- set parameters to best-fit values
   debug["execute"] << "Setting Alpos parameters to best fit values..." << endl;
   for ( unsigned int ip =0; ip<allparnames.size() ; ip++) {
      if (TheoryHandler::Handler()->FindParameter(allparnames[ip]) != "") {
         info["execute"] << "Setting fitted value '" << allparnames[ip] << "' to: " << pars[ip] << " +/- " <<
         errs[ip] << endl;
         SET_ANY(allparnames[ip], pars[ip], errs[ip]);
      }
      else {
         //FIXME: warning masks real cause of error -> implement way of checking if nuisance parameter
         warn["execute"] << "Alpos parameter '" << allparnames[ip] << "' not found: Assuming nuisance parameter and skipping..." << endl;
      }
   }

   //--- store results in the root output file
   TDirectory* rootfile        =  Alpos::Current()->Settings()->rootoutput;
   if ( rootfile ) {
      rootfile->cd();
      TH1D *histFitParameters=new TH1D
	 ("fitparameters",";parameter",allparnames.size(),-0.5,allparnames.size()-0.5);
      TH2D *histFitCorrelations=new TH2D
	 ("fitcorrelations",";parameter",
	  allparnames.size(),-0.5,allparnames.size()-0.5,
	  allparnames.size(),-0.5,allparnames.size()-0.5);
      TFitResult fitResult(fFitter->Result());
      fitResult.SetName("fitresult");
      fitResult.Write();

      /* if(fFitter->GetMinimizer()->CovMatrixStatus()<2) {
	 cout<<"status before Hesse: "<<fFitter->GetMinimizer()->CovMatrixStatus()<<"\n"
	 fFitter->GetMinimizer()->Hesse();
	 cout<<"status after Hesse: "<<fFitter->GetMinimizer()->CovMatrixStatus()<<"\n"
	 }
	 if(fFitter->GetMinimizer()->CovMatrixStatus()<3) {
	 cout<<"status before GetCovMatrix: "<<fFitter->GetMinimizer()->CovMatrixStatus()<<"\n"
	 fFitter->GetMinimizer()->GetCovMatrix(rhoij.GetMatrixArray());
	 }
      */

      for ( unsigned int ip = 0; ip<allparnames.size(); ip++) {
         if (fParIsFixed[ip]) {
            cout << " " << allparnames[ip] << " " << pars[ip] << " (fixed)" << std::endl;
            continue;
         }
         cout << " " << allparnames[ip] << " " << pars[ip] << " " << errs[ip] << " " << fitResult.Correlation(ip, ip) << "\n";
      }

      for ( unsigned int ip = 0; ip<allparnames.size(); ip++) {
	 histFitParameters->SetBinContent(ip+1, pars[ip]);
	 histFitParameters->SetBinError(ip+1, errs[ip]);
	 histFitParameters->GetXaxis()->SetBinLabel(ip+1,allparnames[ip].c_str());
	 for ( unsigned int jp =0; jp<allparnames.size() ; jp++) {
            histFitCorrelations->SetBinContent(ip+1, jp+1, 0);
            if ((!fParIsFixed[ip]) && (!fParIsFixed[jp])) {
               histFitCorrelations->SetBinContent(ip+1, jp+1, fitResult.Correlation(ip, jp));
            }
	 }
      }

      histFitParameters->Write();
      histFitCorrelations->Write();
   }

   // --- job done successfully
   return fitOK;
}


//____________________________________________________________________________________ //
void AFitter::PrintNuisanceParameters() const {
   //! Print nuisance parameters if these have been calculated
   map<string, double> nui = fChisq->GetNuisanceParameters();
   cout << endl;
   if (nui.empty()) {
      info["PrintNuisanceParameters"] << "No Nuisance parameters have been calculated by chisq-definition '" <<
      STRING_NS(Chisq, NS()) << "'." << endl;
   }
   else {
      info["PrintNuisanceParameters"] << "Printing nuisance parameters calculated by '" << STRING_NS(Chisq, NS()) <<
      "'." << endl;
      for (auto ie : nui) {
         printf("   %-22s\t\t% 5.3f\n", ie.first.c_str(), ie.second);
      }
   }

}

//____________________________________________________________________________________ //


void AFitter::FitDataTheoryPairs() const {
   //! Execute fit for subsets and/or DataTheoryPairs
   //! This can be called in addition to the normal SuperData|SuperTheory fitter

   info["FitDataTheoryPairs"] << "Doing subfits of DataTheoryPairs as requested." << endl;

   // --- setting up minimzer
   debug["FitDataTheoryPairs"] << "Init sub-minimizer." << endl;

   ROOT::Fit::Fitter* subFitter = new ROOT::Fit::Fitter();
   string minimizer = STRING_NS(Minimizer, NS());
   subFitter->Config().SetMinimizer(minimizer.c_str());


   // --- Configure fitter
   int PrintLevel = INT_NS(PrintLevel, NS());
   double Tolerance = DOUBLE_NS(Tolerance, NS());
   int Strategy = INT_NS(Strategy, NS());
   subFitter->Config().MinimizerOptions().SetPrintLevel(PrintLevel);
   subFitter->Config().MinimizerOptions().SetTolerance(Tolerance);
   subFitter->Config().MinimizerOptions().SetStrategy(Strategy);
   subFitter->Config().MinimizerOptions().SetMaxFunctionCalls(300000);
   subFitter->Config().MinimizerOptions().SetMaxIterations(300000);

   // --- set chisq function
   string chisqdef = STRING_NS(Chisq, NS());

   AChisqBase* subChisq;

   for (const auto& dtpId : TheoryHandler::Handler()->GetDataTheoryPairs()) {
      // initialize a sub-fit for each DataTheoryPair
      //////////////////////////////////////////////////

      debug["FitDataTheoryPairs"] << "Fitting DataTheoryPair '" << dtpId.first << "'." << endl;
      subChisq = AFactory::ChisqFactory(chisqdef, fParNames, dtpId.second.first, dtpId.second.second);


      if (!subChisq) {
         error["FitDataTheoryPairs"] << "Failed to initialize chisq with name: " << chisqdef << endl;
         exit(1);
      }
      vector<string> allsubparnames = subChisq->GetFitPar();

      // --- set fit parameters
      vector<ROOT::Fit::ParameterSettings> parsettings(subChisq->NDim());
      unsigned int i = 0;
      for (auto& iset : parsettings) {
         if (i < fParNames.size()) {
            // main parameters (parameters of interest) may be limited
            if (fParIsLimited[i])
               iset.Set(fParNames[i], PAR_ANY(fParNames[i]), fParStepSizes[i], fParLowerLimits[i],
                        fParUpperLimits[i]); /*name,val,step,lower,upper*/
            else
               iset.Set(fParNames[i], PAR_ANY(fParNames[i]), fParStepSizes[i]); /*name,val,step*/
         }
         else {
            // additional parameters (e.g. nuisance parameters) are not limited
            iset.Set(allsubparnames[i], 0, 1.e-3); /*name,val,step*/
         }
         i++;
      }

      subFitter->Config().SetParamsSettings(parsettings);

      // --- Set chisq function
      subFitter->SetFCN(*subChisq);

      // execute the sub-fit for each DataTheoryPair
      //////////////////////////////////////////////////
      // --- Get initial chisq
      subFitter->EvalFCN();
      double subInitialChisq = subFitter->Result().MinFcnValue();

      // --- Do fit
      bool fitOK = subFitter->FitFCN();
      if (!fitOK) {
         warn["FitDataTheoryPairs"] << "Subfit failed." << endl;
      }

      // --- Get final chisq
      ROOT::Fit::FitResult result = subFitter->Result();
      double subFinalChisq = result.MinFcnValue();
      info["FitDataTheoryPairs"] << "Subfit initial chisq: " << subInitialChisq << endl;
      info["FitDataTheoryPairs"] << "Subfit final chisq:   " << subFinalChisq << endl;

      // --- set errors (and also set again central values)
      const double* pars = subFitter->Result().GetParams();
      const double* errs = subFitter->Result().GetErrors();
      //vector<string> allsubparnames = fChisq->GetFitPar();
      // for ( unsigned int ip ; ip<allsubparnames.size() ; ip++) {
      //    info["execute"]<<"Setting fitted value '"<<allsubparnames[ip]<<"' to: "<<pars[ip]<<" +/- "<<errs[ip]<<endl;
      //    SET_ANY(allsubparnames[ip],pars[ip],errs[ip]);
      // }

      // --- More printing options
      //if ( BOOL_NS(PrintResults,NS()) ) {
      subFitter->Result().Print(cout);
      PrintNuisanceParameters();
      //}

      //if ( BOOL_NS(PrintCovariance,NS()) )
      subFitter->Result().PrintCovMatrix(cout);

   }

   //TODO: same for Subsets:
   // for ( const auto& is: TheoryHandler::Handler()->GetAllSubsetPairs() )

   delete subFitter;
   delete subChisq;
}
