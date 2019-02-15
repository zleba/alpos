/**
 * 
 *  AFactory
 * 
 *  Instantiate a task or a function
 * 
 **/

#include "alpos/AFactory.h"
#include "fastnlotk/speaker.h"
#include <string>
#include <cmath>

// functions
#include "alpos/functions/AExampleFunction.h"
#include "alpos/functions/ACRunDecFunction.h"
#if _CMAKE_FOUND_APFEL
#include "alpos/functions/AApfel.h"
#include "alpos/functions/AApfelDISCS.h"
#include "alpos/functions/AApfelDDISCS.h"
#include "alpos/functions/AApfelDISCSEWFit.h"
#endif //_CMAKE_FOUND_APFEL
#include "alpos/functions/ALhapdf6.h"
#include "alpos/functions/ALhapdf6Alphas.h"
#include "alpos/functions/AStrowp1.h"
#include "alpos/functions/AH1DPDF2006.h"
//#include "alpos/AAlphasDependentPDF.h"
#include "alpos/functions/AfastNLO.h"
#include "alpos/functions/AfastNLOalt.h"
#include "alpos/functions/AfastNLODiffDIS.h"
#include "alpos/functions/AfastNLOnormDIS.h"
#include "alpos/functions/AfastNLOnormDISalt.h"
//#include "alpos/AfastNLODiffDIS.h"
#include "alpos/functions/AfastNLOInterpolPDFas.h"
#include "alpos/functions/AfastNLOInterpolPDFasNormDIS.h"
#include "alpos/functions/AfastNLORatio.h"
#if _CMAKE_FOUND_APPLGRID
    #include "alpos/functions/AApplgrid.h"
#endif //_CMAKE_FOUND_APPLGRID

#if _CMAKE_FOUND_APFELxx
    #include "alpos/functions/AApfelxxPDF.h"
    #include "alpos/functions/AApfelxxAlphas.h"
    #include "alpos/functions/AApfelxxDISCS.h"
    #include "alpos/functions/AApfelxxDDISCS.h"
#endif //_CMAKE_FOUND_Apfelxx
#include "alpos/functions/ADPDF.h"

#include "alpos/AData.h"
#include "alpos/functions/ASingleConstant.h"
#include "alpos/functions/AFunctionTF1.h"
#include "alpos/ASuperTheory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASubsetFunction.h"
#include "alpos/ASubsetData.h"
#include "alpos/functions/AQcdnumInit.h"
#include "alpos/functions/AQcdnumAlphas.h"
#include "alpos/functions/AQcdnumPDF.h"
#include "alpos/functions/AQcdnumDISCS.h"
#include "alpos/functions/AQcdnumDDISCS.h"
#include "alpos/functions/APDFQ0_diff.h"
#include "alpos/functions/AQcdnumDISCSEWFit.h"
#include "alpos/functions/AAlphaEmRun.h"
#include "alpos/functions/AEprc.h"
#include "alpos/functions/APDFQ0_LHAPDF.h"
#include "alpos/functions/APDFQ0_QcdnumExample.h"
#include "alpos/functions/APDFQ0_HERAStyle.h"
#include "alpos/functions/APDFQ0_HERA.h"
#include "alpos/functions/APDFQ0_BiLog.h"
#include "alpos/functions/ATwoPomeronModel.h"
#include "alpos/functions/ATwoPomeronParams.h"
#include "alpos/functions/AUserFunction.h"
#include "alpos/functions/AFixedValues.h"

// tasks
#include "alpos/tasks/AExampleTask.h"
#include "alpos/tasks/AFitter.h"
#include "alpos/tasks/AApcalc.h"
#include "alpos/tasks/AApcFitter.h"
#include "alpos/tasks/AConstLQFitter.h"
#include "alpos/tasks/APrintTheorySet.h"
#include "alpos/tasks/ASaveTheorySet.h"
#include "alpos/tasks/APrintSteering.h"
#include "alpos/tasks/AStatAnalysis.h"
#include "alpos/tasks/APrintErrorSummary.h"
#include "alpos/tasks/APrintDataTheory.h"
#include "alpos/tasks/AChi2FitPDFas.h"
#include "alpos/tasks/AChi2InterpolPDFas.h"
#include "alpos/tasks/AChi2Scan.h"
#include "alpos/tasks/AContour.h"
#include "alpos/tasks/AReplaceDataWithTheoryValues.h"
#include "alpos/tasks/AApcFitter.h"
#include "alpos/tasks/ASavePDFTGraph.h"
#include "alpos/tasks/ASaveDPDFTGraph.h"
#include "alpos/tasks/ASaveDataTheory.h"

#include "alpos/tasks/AScaleUncertainty.h"
#include "alpos/tasks/APDFUncer.h"
#include "alpos/tasks/ALHAPDFErrors.h"

// chisq functions
#include "alpos/AChisq.h"


using namespace std;

//____________________________________________________________________________________ //
AFuncD* AFactory::FunctionFactory(const std::string& functype,const std::string& funcname) {
   //! init a new function
   //! Initialize a function with a certain name
   //! All available function must be included here and as '#include' above

   say::info["AFactory::FunctionFactory"]<<"Initializing function of type '"<<functype<<"', as name '"<<funcname <<"'."<<endl;
   
   AFuncD* ptr = NULL;
   if      ( functype == AExampleFunction::fFunctionName )  ptr = new AExampleFunction(funcname);
   else if ( functype == ACRunDecFunction::fFunctionName )  ptr = new ACRunDecFunction(funcname);
#if _CMAKE_FOUND_APFEL
   else if ( functype == AApfelInit::fFunctionName )        ptr = new AApfelInit(funcname);
   else if ( functype == AApfelPDF::fFunctionName )         ptr = new AApfelPDF(funcname);
   else if ( functype == AApfelAs::fFunctionName )          ptr = new AApfelAs(funcname);
   else if ( functype == AApfelQEDEvol::fFunctionName )     ptr = new AApfelQEDEvol(funcname);
   else if ( functype == AApfelDISCS::fFunctionName )       ptr = new AApfelDISCS(funcname);
   else if ( functype == AApfelDISCSEWFit::fFunctionName )  ptr = new AApfelDISCSEWFit(funcname);
   else if ( functype == AApfelDDISCS::fFunctionName )      ptr = new AApfelDDISCS(funcname);
#endif //_CMAKE_FOUND_APFEL
#if _CMAKE_FOUND_APFELxx
   else if ( functype == AApfelxxDISCS::fFunctionName )       ptr = new AApfelxxDISCS(funcname);
   else if ( functype == AApfelxxDDISCS::fFunctionName )      ptr = new AApfelxxDDISCS(funcname);
   else if ( functype == AApfelxxAlphas::fFunctionName )      ptr = new AApfelxxAlphas(funcname);
   else if ( functype == AApfelxxPDF::fFunctionName )         ptr = new AApfelxxPDF(funcname);
#endif //_CMAKE_FOUND_APFELxx
   else if ( functype == ADPDF::fFunctionName )             ptr = new ADPDF(funcname);
   else if ( functype == ALhapdf6::fFunctionName )          ptr = new ALhapdf6(funcname);
   else if ( functype == ALhapdf6Alphas::fFunctionName )    ptr = new ALhapdf6Alphas(funcname);
   else if ( functype == AStrowp1::fFunctionName )          ptr = new AStrowp1(funcname);
   else if ( functype == AH1DPDF2006::fFunctionName )       ptr = new AH1DPDF2006(funcname);
   //else if ( functype == AAlphasDependentPDF::fFunctionName)ptr = new AAlphasDependentPDF(funcname);
   else if ( functype == AfastNLO::fFunctionName )          ptr = new AfastNLO(funcname);
   else if ( functype == AfastNLOalt::fFunctionName )       ptr = new AfastNLOalt(funcname);
   else if ( functype == AfastNLODiffDIS::fFunctionName )   ptr = new AfastNLODiffDIS(funcname);
#if _CMAKE_FOUND_QCDNUM
   else if ( functype == AfastNLOnormDIS::fFunctionName )   ptr = new AfastNLOnormDIS(funcname);
   else if ( functype == AfastNLOnormDISalt::fFunctionName )   ptr = new AfastNLOnormDISalt(funcname);
   //else if ( functype == AfastNLODiffDIS::fFunctionName )   ptr = new AfastNLODiffDIS(funcname);
#endif //_CMAKE_FOUND_QCDNUM
   else if ( functype == AfastNLOInterpolPDFas::fFunctionName )       ptr = new AfastNLOInterpolPDFas(funcname);
#if _CMAKE_FOUND_QCDNUM
   else if ( functype == AfastNLOInterpolPDFasNormDIS::fFunctionName )       ptr = new AfastNLOInterpolPDFasNormDIS(funcname);
#endif //_CMAKE_FOUND_QCDNUM
   else if ( functype == AfastNLORatio::fFunctionName )     ptr = new AfastNLORatio(funcname);
#if _CMAKE_FOUND_APPLGRID
   else if ( functype == AApplgrid::fFunctionName )         ptr = new AApplgrid(funcname);
#endif //_CMAKE_FOUND_APPLGRID
   else if ( functype == ASuperTheory::fFunctionName )      ptr = new ASuperTheory(funcname); // should not be used in steering
   else if ( functype == ASuperData::fFunctionName )        ptr = new ASuperData(funcname); // should not be used in steering
   else if ( functype == ASubsetFunction::fFunctionName )   ptr = new ASubsetFunction(funcname); // should not be used in steering
   else if ( functype == ASubsetData::fFunctionName )       ptr = new ASubsetData(funcname); // should not be used in steering
   else if ( functype == AData::fFunctionName )             ptr = new AData(funcname); //< A data set
   else if ( functype == ASingleConstant::fFunctionName )   ptr = new ASingleConstant(funcname);
   else if ( functype == AFunctionTF1::fFunctionName )      ptr = new AFunctionTF1(funcname);
#if _CMAKE_FOUND_QCDNUM
   else if ( functype == AQcdnumInit::fFunctionName )       ptr = new AQcdnumInit(funcname);
   else if ( functype == AQcdnumAlphas::fFunctionName )     ptr = new AQcdnumAlphas(funcname);
   else if ( functype == AQcdnumPDF::fFunctionName )        ptr = new AQcdnumPDF(funcname);
   else if ( functype == AQcdnumDISCS::fFunctionName )      ptr = new AQcdnumDISCS(funcname);
   else if ( functype == AQcdnumDDISCS::fFunctionName )     ptr = new AQcdnumDDISCS(funcname);

   else if ( functype == AQcdnumDISCSEWFit::fFunctionName ) ptr = new AQcdnumDISCSEWFit(funcname);
#endif //_CMAKE_FOUND_QCDNUM
   else if ( functype == AAlphaEmRun::fFunctionName )       ptr = new AAlphaEmRun(funcname);
   else if ( functype == AEprc::fFunctionName )	            ptr = new AEprc(funcname);
   else if ( functype == APDFQ0_LHAPDF::fFunctionName )     ptr = new APDFQ0_LHAPDF(funcname);
   else if ( functype == APDFQ0_QcdnumExample::fFunctionName )     ptr = new APDFQ0_QcdnumExample(funcname);
   else if ( functype == APDFQ0_HERAStyle::fFunctionName )  ptr = new APDFQ0_HERAStyle(funcname);
   else if ( functype == APDFQ0_diff::fFunctionName )       ptr = new APDFQ0_diff(funcname);
   else if ( functype == APDFQ0_HERA::fFunctionName )       ptr = new APDFQ0_HERA(funcname);
   else if ( functype == APDFQ0_BiLog::fFunctionName )      ptr = new APDFQ0_BiLog(funcname);
   else if ( functype == ATwoPomeronModel::fFunctionName )      ptr = new ATwoPomeronModel(funcname);
   else if ( functype == ATwoPomeronParams::fFunctionName )     ptr = new ATwoPomeronParams(funcname);
   else if ( functype == AUserFunction::fFunctionName )         ptr = new AUserFunction(funcname);
   else if ( functype == AFixedValues::fFunctionName )          ptr = new AFixedValues(funcname);
   // ...
   // ...
   // ... else if ( functype == "YourFunction" ) new ANewFunction(funcname);
   // ...
   // ...
   else {
      say::error["AFactory::FunctionFactory"]<<"Failed to identify function '"<<functype<<"'. Exiting."<<endl;
      exit(1);
      return NULL;
   }
   say::info["AFactory::FunctionFactory"]<<"Initializing function of type '"<<functype<<"', as name '"<<funcname <<"' DONE."<<endl;
   
   return ptr;
}



//____________________________________________________________________________________ //
ATask* AFactory::TaskFactory(const std::string& tname, const std::string& ttype) {
   //! instantiate a new task with name 'tname' of type 'ttype'
   say::info["AFactory::TaskFactory"]<<"Initializing task of type '"<<ttype<<"', as name '"<<tname <<"'."<<endl;
   if      ( ttype == AExampleTask::fTaskType )       return new AExampleTask(tname);
   else if ( ttype == AFitter::fTaskType )            return new AFitter(tname);//AFitter(tname,fSteerfile,&fResults);
   else if ( ttype == AApcalc::fTaskType )            return new AApcalc(tname);//AFitter(tname,fSteerfile,&fResults);
   else if ( ttype == AApcFitter::fTaskType )         return new AApcFitter(tname);//AFitter(tname,fSteerfile,&fResults);
   else if ( ttype == AConstLQFitter::fTaskType )     return new AConstLQFitter(tname);//AFitter(tname,fSteerfile,&fResults);
   else if ( ttype == APrintTheorySet::TaskType() )   return new APrintTheorySet(tname);
   else if ( ttype == ASaveTheorySet::TaskType() )    return new ASaveTheorySet(tname);
   else if ( ttype == APrintSteering::TaskType() )    return new APrintSteering(tname);
   else if ( ttype == AStatAnalysis::fTaskType )      return new AStatAnalysis(tname);
   else if ( ttype == APrintErrorSummary::fTaskType ) return new APrintErrorSummary(tname);
   else if ( ttype == APrintDataTheory::fTaskType )   return new APrintDataTheory(tname);
   else if ( ttype == AChi2FitPDFas::fTaskType )      return new AChi2FitPDFas(tname);
   else if ( ttype == AChi2InterpolPDFas::fTaskType ) return new AChi2InterpolPDFas(tname);
   else if ( ttype == AChi2Scan::fTaskType )          return new AChi2Scan(tname);
   else if ( ttype == AContour::fTaskType )           return new AContour(tname);
   else if ( ttype == AReplaceDataWithTheoryValues::fTaskType )           return new AReplaceDataWithTheoryValues(tname);
   else if ( ttype == AApcFitter::fTaskType )         return new AApcFitter(tname);
   else if ( ttype == ASavePDFTGraph::fTaskType )     return new ASavePDFTGraph(tname);
   else if ( ttype == ASaveDPDFTGraph::fTaskType )    return new ASaveDPDFTGraph(tname);
   else if ( ttype == ASaveDataTheory::fTaskType )    return new ASaveDataTheory(tname);
   else if ( ttype == AScaleUncertainty::fTaskType )  return new AScaleUncertainty(tname);
   else if ( ttype == APDFUncer::fTaskType )          return new APDFUncer(tname);
   else if ( ttype == ALHAPDFErrors::fTaskType )      return new ALHAPDFErrors(tname);
   // ...
   // ...
   // ...
   else {
      say::warn["AFactory::TaskFactory"]<<"Failed to identify task '"<<ttype<<"'."<<endl;
      return nullptr;
   }
}


//____________________________________________________________________________________ //
//ROOT::Math::IMultiGenFunction*
AChisqBase* AFactory::ChisqFactory(const std::string& chisq, const std::vector<std::string>& par, AData* data, AFuncD* theo ){
   //! instantiate a chi^2 class
   say::debug["AFactory::ChisqFactory"]<<"Initializing chisq of type '"<<chisq<<"'."<<endl;
   // --- chisq using full 'lognormal' prescription
   if      ( chisq == AChisqLogNormal::GetChisqName() ) return new AChisqLogNormal(par,data,theo);
   else if ( chisq == AChisqLogNormalNuisance::GetChisqName() ) return new AChisqLogNormalNuisance(par,data,theo);
   else if ( chisq == AChisqLogNormalNuisanceFit::GetChisqName() ) return new AChisqLogNormalNuisanceFit(par,data,theo);
   else if ( chisq == AChisqLogNormalStatUncorr::GetChisqName() ) return new AChisqLogNormalStatUncorr(par,data,theo);
   else if ( chisq == AChisqNormalLogNormal::GetChisqName() ) return new AChisqNormalLogNormal(par,data,theo);
   // --- chisq without covariance matrices
   else if ( chisq == AChisqSimple::GetChisqName() ) return new AChisqSimple(par,data,theo);
   else if ( chisq == AChisqSimpleNuisanceMultFit::GetChisqName() ) return new AChisqSimpleNuisanceMultFit(par,data,theo);
   else if ( chisq == AChisqSimpleNuisanceAddFit::GetChisqName() ) return new AChisqSimpleNuisanceAddFit(par,data,theo);
   // else if ( chisq == AChisqSimpleNuisanceMult::GetChisqName() ) return new AChisqSimpleNuisanceMult(par,data,theo);
   // else if ( chisq == AChisqSimpleNuisanceAdd::GetChisqName() ) return new AChisqSimpleNuisanceAdd(par,data,theo);
   // --- chisq with covariance matrices and nuisance parameters
   else if ( chisq == AChisqCov::GetChisqName() ) return new AChisqCov(par,data,theo);
   else if ( chisq == AChisqCMS::GetChisqName() ) return new AChisqCMS(par,data,theo);
   else if ( chisq == AChisqCovMult::GetChisqName() ) return new AChisqCovMult(par,data,theo);
   else if ( chisq == AChisqCovStatUncorr::GetChisqName() ) return new AChisqCovStatUncorr(par,data,theo);
   else if ( chisq == AChisqNuisanceAdd::GetChisqName() ) return new AChisqNuisanceAdd(par,data,theo);
   else if ( chisq == AChisqNuisanceMult::GetChisqName() )  return new AChisqNuisanceMult(par,data,theo);
   // else if ( chisq == AChisqNuisanceMultFit::GetChisqName() ) return new AChisqNuisanceMultFit(par,data,theo);
   // else if ( chisq == AChisqNuisanceAddFit::GetChisqName() ) return new AChisqNuisanceAddFit(par,data,theo);
   // --- Attempt for HERAFitter style chisq's
   else if ( chisq == AChisqHERAFitterDefault::GetChisqName() ) return new AChisqHERAFitterDefault(par,data,theo); 
   else if ( chisq == AChisqHERAFitterLogDefault::GetChisqName() ) return new AChisqHERAFitterLogDefault(par,data,theo); 
   else if ( chisq == AChisqHERAFitterDefaultMatrix::GetChisqName() ) return new AChisqHERAFitterDefaultMatrix(par,data,theo);
   else if ( chisq == AChisqHERAFitterDefaultFit::GetChisqName() ) return new AChisqHERAFitterDefaultFit(par,data,theo);
   else if ( chisq == AChisqHERAFitterFull::GetChisqName() ) return new AChisqHERAFitterFull(par,data,theo);
   else if ( chisq == AChisqHERAFitterFullImproved::GetChisqName() ) return new AChisqHERAFitterFullImproved(par,data,theo);
   else if ( chisq == AChisqD0Fit::GetChisqName() ) return new AChisqD0Fit(par,data,theo);
   else if ( chisq == AChisqD0StatCorrFit::GetChisqName() ) return new AChisqD0StatCorrFit(par,data,theo);
   else if ( chisq == APull::GetChisqName() ) return new APull(data,theo);
   else  {
      say::error["AFactory::ChisqFactory"]<<"Failed to identify chisq class '"<<chisq<<"'. Exiting."<<endl;
      exit(1);
      return nullptr;
   }
   return nullptr;
}

