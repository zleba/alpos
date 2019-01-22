// DB 15.01.2015

#include "alpos/functions/AfastNLOInterpolPDFasNormDIS.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include <iostream>
#include <fastnlotk/read_steer.h>
#include <alpos/functions/ALhapdf6.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <TSpline.h>
#include <alpos/functions/AfastNLOInterpolPDFas.h>
#include "TF1.h"
#include "TList.h"

// need structure function calculation from QCDNUM
extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);


using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLOInterpolPDFasNormDIS::fRequirements = {//"Filename",                    // fastNLO table name
                        "LHAPDFasSeries",              // LHAPDF set with alpha_s series
                        "AlphasMz",                    // alpha_s(mz) of current cross sections
                        "ScaleFacMuR","ScaleFacMuF",   // scale factors
                        "Units",                       // units consistent wiht input data
                        "iOrd",                        // order of calculation
                        "iThr",                        // order of threshold corrections
                        "MuRFuncForm","MuFFuncForm",   // mur and muf functional form for fastNLO flexible scale tables
                        "FitFunc",                     // ROOT::TF1 fit function for interpolation ('pol2')
                        "nQ2",                         // # points for Q2 grid
                        "nx",                          // # points for x  grid
                        "AemRun",						    // alpha_em running coupling
                        "QcdnumInit"						 // QcdnumInit singleton
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOInterpolPDFasNormDIS::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLOInterpolPDFasNormDIS::fFunctionName = "fastNLOInterpolPDFasNormDIS"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLOInterpolPDFasNormDIS::AfastNLOInterpolPDFasNormDIS(const std::string& name) : AfastNLOInterpolPDFas(name) {
   // Remember: no access to parameters possible in constructor!
   SetClassName("AfastNLOInterpolPDFasNormDIS");
}


// __________________________________________________________________________________________ //
AfastNLOInterpolPDFasNormDIS::~AfastNLOInterpolPDFasNormDIS() {
   // nothing to do
}



// ___________________________________________________________________________________________ //
bool AfastNLOInterpolPDFasNormDIS::ReInit() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   info["ReInit"]<<"Re-initializing AfastNLOInterpolPDFas"<<endl;

   using namespace AlposTools;

   //string filename = PAR_S(Filename);
   string lhapdfstring = PAR_S(LHAPDFasSeries);

   vector<string> filenames;
   if (EXIST_NS(Filename, GetAlposName()))
      filenames.push_back(STRING_NS(Filename, GetAlposName()));
   else if (EXIST_NS(Filenames, GetAlposName())) {
      filenames = STRING_ARR_NS(Filenames, GetAlposName());
   }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit);
   UPDATE(AemRun);

   fGraphs.clear();
   // split lhapdfstring into vector<string>
   vector<string> lhapdfnames;
   //boost::split(lhapdfnames, lhapdfstring, boost::is_any_of(", "), boost::token_compress_on);

   //RADEK hack instead of boost
   char *token;
   char *str = new char[10000];
   strcpy(str, lhapdfstring.c_str());
   assert(str != NULL);
   while ((token = strsep(&str, ",")) != NULL)
       lhapdfnames.push_back(token);
   delete [] str;





   // some PDFs are published as separate filesets for each value of alpha_s
   // e.g. for CT10..., which uses CT10..._as_0112 etc.
   // others have PDFs for all values of alpha_s in a single fileset
   // e.g. for MSTW2008nlo_asmzrange
   //     -> need to handle cases differently

   int nTotalAlphasPoints = 0;       // # of total alpha_s points for interpolation (from all PDF filesets)
   bool qSinglePDFFileset;

   // check # of pdf filesets and # of members
   if (lhapdfnames.size() == 1) {
      // one PDF fileset with several alpha_s values
      info["ReInit"]<<"Using single multi-member PDF fileset '"<<lhapdfnames[0]<<"' as an alpha_s(Mz)-series PDF."<<endl;

      qSinglePDFFileset = true;
   }
   else {
      // several PDF filesets, each with one alpha_s value
      info["ReInit"]<<"Using several single-member PDF filesets, one for each alpha_s(Mz)."<<endl;

      qSinglePDFFileset = false;
   }

   // preliminary run through all lhapdfnames
   //vector<fastNLOLHAPDF*> fnlos;
   set<double> processedAsmzValues;  // set to register alpha_s(M_Z) values
   for (string lhapdfname : lhapdfnames) {
      debug["ReInit"]<<"Reading LHAPDF Info for PDF fileset '"<<lhapdfname<<"'..."<<endl;
      LHAPDF::PDFSet* lhapdfset = new LHAPDF::PDFSet(lhapdfname);

      int nMembers = lhapdfset->size();

      // check # of pdf filesets and # of members
      if (qSinglePDFFileset) {
         // single-file, multi-member mode
         if (nMembers == 1) {
            error["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has a single member."
            <<"No interpolation possible. Exiting."<<endl;
            exit(1);
         }
         else {
            nTotalAlphasPoints = nMembers;
            debug["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has "<<nMembers<<" members. "<<endl;

         }
      }
      else {
         // multi-file, single-member mode
         if (nMembers != 1) {
            warn["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has more than one member, but only"
            <<"one PDF is needed per alpha_s(M_Z) point. Using member #0..."<<endl;
            nMembers = 1;
         }
         else {
            debug["ReInit"]<<"PDF fileset '"<<lhapdfname<<"' has 1 member."<<endl;
         }
         nTotalAlphasPoints += 1;
      }

      // check for alpha_s(M_Z) duplicates
      for ( int iMember = 0 ; iMember<nMembers ; iMember++ ) {
         LHAPDF::PDF *lhapdfmember = lhapdfset->mkPDF(iMember);
         double asmz = lhapdfmember->alphasQ(91.1876);

         if ( processedAsmzValues.find(asmz) != processedAsmzValues.end() ) {
            // already processed a PDF member with the same alpha_s(M_Z) -> exclude from TGraph
            error["ReInit"]<<"Already encountered this value of alpha_s(M_Z)! "
            <<"Multiple chi2 points for the same value of alpha_s(M_Z) not supported. Exiting."<<endl;
            exit(1);
         }
         processedAsmzValues.insert(asmz);   // store processed alpha_s(M_Z)
      }
   }  // end preliminary loop through lhapdfnames
   debug["ReInit"]<<"Preliminary loop done."<<endl;
   info["ReInit"]<<"Interpolating "<<nTotalAlphasPoints<<" alpha_s(M_Z) points."<<endl;

   //debug["ReInit"]<<"Loading fastNLOLHAPDF for PDF fileset '"<<lhapdfname<<"'..."<<endl;
   debug["ReInit"]<<"Loading fastNLOLHAPDF."<<endl;


   std::vector<fastNLOLHAPDF*> fnlos;

   int nBins = 0;
   for (auto fn : filenames) {
      info["Init"] << "Reading table file: " << fn << endl;
      fastNLOLHAPDF* f = new fastNLOLHAPDF(fn, lhapdfnames[0], 0);
      f->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
      if (f->GetIsFlexibleScaleTable()) {
         f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm, GetAlposName())));
         f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm, GetAlposName())));
      }
      SetOrder(f);
      f->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));
      nBins += f->GetNObsBin();
      fnlos.push_back(f);
   }

   fGraphs.resize(nBins);

   // -------- Calculate NC DIS cross section for normalization
   // WIP: get rid of some hard-coded values


   // erase previous NC DIS cross sections
   q2LoUpCSmap.clear();

   vector <pair<double,double>> binQ2LoUps;

   for ( auto f : fnlos ) {
      // go through fastNLO bins
      for ( auto binsLoUp : f->GetObsBinsBounds(0) ) {
         binQ2LoUps.push_back(binsLoUp);
         if (q2LoUpCSmap.find(binsLoUp) == q2LoUpCSmap.end()) {
            // if we haven't already calculated this CS
            debug["Update"]<<"NC DIS cross section for "<<binsLoUp.first<<" < Q2 < "<<binsLoUp.second<<" not cached. Calculating..."<<endl;
            q2LoUpCSmap[binsLoUp] = CalcNCDISCrossSection(binsLoUp.first, binsLoUp.second);
         }
         else {
            debug["Update"]<<"NC DIS cross section for "<<binsLoUp.first<<" < Q2 < "<<binsLoUp.second<<" cached. Skipping calculation..."<<endl;
         }
      }
   }

   // --------- End NC DIS cross section calculation

   // run through all lhapdfnames and
   for ( int iLhapdfname = 0 ; iLhapdfname < lhapdfnames.size() ; iLhapdfname++ ) {
      info["ReInit"]<<"Processing PDF fileset '"<<lhapdfnames[iLhapdfname]<<"'..."<<endl;

      for (auto fnlo : fnlos) {
         fnlo->SetLHAPDFFilename(lhapdfnames[iLhapdfname]);
      }

      info["ReInit"]<<"Running over all " << fnlos[0]->GetNPDFMembers() <<" PDF members."<<endl;
      for ( int iFnloMember = 0 ; iFnloMember<fnlos[0]->GetNPDFMembers() ; iFnloMember++ ) {
         info["ReInit"]<<"Processing PDF set '"<<lhapdfnames[iLhapdfname]<<"', member #"<<iFnloMember<<endl;

         vector<double> cs;

         for (auto fnlo : fnlos) {
            fnlo->SetLHAPDFMember(iFnloMember);
            fnlo->CalcCrossSection();
            cs += fnlo->GetCrossSection();
         }

         double asmz = fnlos[0]->GetAlphasMz();
         info["ReInit"]<<"Member #"<<iFnloMember<<" has alpha_s(M_Z) = "<<asmz<<endl;

         // decide how to fill CS Points into fGraphs
         int iGraphPoint;
         if (qSinglePDFFileset) {
            iGraphPoint = iFnloMember;
         }
         else {
            iGraphPoint = iLhapdfname;
         }

         // Fill points into TGraph
         debug["ReInit"]<<"Adding CS points to TGraph"<<endl;
         for ( int iCSBin = 0 ; iCSBin<cs.size() ; iCSBin++ ) {
            fGraphs[iCSBin].SetPoint(iGraphPoint, asmz, cs[iCSBin]/q2LoUpCSmap[binQ2LoUps[iCSBin]]);
            debug["ReInit"]<<"fGraphs["<<iCSBin<<"].SetPoint("<<iGraphPoint<<", "<<asmz<<", cs["<<iCSBin<<"]="<<cs[iCSBin]/q2LoUpCSmap[binQ2LoUps[iCSBin]]<<");"<<endl;
         }

         // only process first member, if not in single PDF set mode
         if (!qSinglePDFFileset) {
            break;
         }

      } // end loop over PDFMembers

   } // end loop over fnlos

   // sort the graphs by asmz in ascending order
   for (unsigned int iCSBin = 0; iCSBin < fGraphs.size(); iCSBin++) {
      fGraphs[iCSBin].Sort();
   }

   // do interpolation fits
   string FitFunc = PAR_S(FitFunc);
//   for ( auto& fg : fGraphs ) {

   // should not be necessary, but include for safety
   for (unsigned int iSpline = 0; iSpline < fSplines.size(); iSpline++) {
      if (fSplines[iSpline]) delete fSplines[iSpline];
   }

   fSplines.resize(nBins);
   for (unsigned int iCSBin = 0; iCSBin < nBins; iCSBin++) {
      auto& fg = fGraphs[iCSBin];
      if ((FitFunc == "TSpline3") || (FitFunc == "Spline3") || (FitFunc == "spline3")) {
         //if (fSplines.size() != nBins) fSplines.resize(nBins);
         fSplines[iCSBin] = new TSpline3(fg.GetName(), &fg);
         this->fHasSplines = true;
      }
      else if ((FitFunc == "TSpline5") || (FitFunc == "Spline5") || (FitFunc == "spline5")) {
         //if (fSplines.size() != nBins) fSplines.resize(nBins);
         fSplines[iCSBin] = new TSpline5(fg.GetName(), &fg);
         this->fHasSplines = true;
      }
      else {
         fg.Fit(FitFunc.c_str(), "Q");
         this->fHasSplines = false;
      }
   }

   // cleanup
   for (auto fnlo : fnlos) {
      delete fnlo;
   }

   fValue.resize(nBins);
   fError.resize(nBins);

   return true;
}


//___________________________________________________________________//
double AfastNLOInterpolPDFasNormDIS::CalcNCDISCrossSection(double q2min, double q2max) {

   //vector<double> q2bins = {150, 200, 270, 400, 700, 5000, 15000};
   //vector<double> CSNCDIS(q2bins.size()-1);

   double CSNCDIS = 0;
   double ymin = 0.2;
   double ymax = 0.7;
   double S = 318.11947*318.11947;

   static const double convfac = 0.389379338e9; //0.389379323e9;
   const int nq2split=PAR(nQ2); // integration constants
   const int nxsplit=PAR(nx);
   int npts = nxsplit*nq2split;
   vector<double> q2(npts), dq2(npts), x(npts), dx(npts), y(npts);
   //vector<double> aem(nq2split*(q2bins.size()-1));
   map<double,double> aem; // <q2,aem>
   //for ( unsigned int ib = 0 ; ib<q2bins.size()-1 ;ib++ ) {
   //double q2min = q2bins[ib];
   //double q2max = q2bins[ib+1];
   int j=0;
   double q2_2 = q2min;
   for ( int iq2=1 ;iq2<=nq2split ; iq2++ ) {
      double q2_1 = q2_2;
      q2_2 = exp( log(q2min)+(log(q2max)-log(q2min))/nq2split*iq2 );
      double q2j = exp( log(q2_1) + 0.5*(log(q2_2) - log(q2_1)) ) ;
      aem[q2j] = QUICK(AemRun,({sqrt(q2j)}))[0];
      for ( int ix=0 ; ix<nxsplit ; ix++ ) {
         dq2[j] = q2_2 - q2_1;
         q2[j]  = q2j;//exp( log(q2_1) + 0.5*(log(q2_2) - log(q2_1)) ) ;
         double xmax = q2[j] / (S * ymin);
         double xmin = q2[j] / (S * ymax);
         x[j]  = xmin + (xmax-xmin)/nxsplit*(ix+0.5);
         dx[j] = (xmax-xmin) / nxsplit;
         y[j]  = q2[j] / (S * x[j]);
         j++;
      }
   }

   int iFL=1, iF2=2;//, iF3=3;
   double A_u = 4./9.;
   double A_d = 1./9.;
   vector<double> part = {0.,A_d,A_u,A_d,A_u,A_d,0.,A_d,A_u,A_d,A_u,A_d,0. };
   vector<double> FL(q2.size()), F2(q2.size());
   int nchk=0;
   zmstfun_(&iFL,&part[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
   zmstfun_(&iF2,&part[0],&x[0],&q2[0],&F2[0],&npts,&nchk);

   for ( int ip = 0 ; ip<npts ;ip++ ) {
      double yplus = 1+(1-y[ip])*(1-y[ip]);
      double CSi = F2[ip]*yplus - y[ip]*y[ip]*FL[ip];
      //double alem = QUICK(AemRun,({sqrt(q2[ip])}))[0];
      double alem = aem[q2[ip]];
      // double alem = aemrun(q2(j))
      //cout<<"Q="<<sqrt(q2[ip])<<"\talem="<<alem<<endl;
      //double alem = 1./137.;
      double factor=2.*M_PI*alem*alem/(x[ip]*q2[ip]*q2[ip])*convfac;
      CSi *= factor* dq2[ip]* dx[ip];
      CSNCDIS += CSi;
   }
   //cout<<endl;
   cout<<"CS NCDIS: "<<CSNCDIS<<endl;
   //printf("CS NCDIS: %14.12f\n",CSNCDIS[ib]);

   return CSNCDIS;
}
