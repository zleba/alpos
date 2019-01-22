// DB 16.01.2015

#include "alpos/functions/AfastNLOnormDISalt.h"
#include "alpos/functions/AQcdnumInit.h"
#include <iostream>

using namespace std;

extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);

// __________________________________________________________________________________________ //
// const std::vector<std::string> AfastNLOnormDISalt::fRequirements = {"QcdnumInit","Filename","fastNLO-PDF","fastNLO-Alpha_s","fastNLO-ScaleFacMuR","fastNLO-ScaleFacMuF","fastNLO-iOrd","fastNLO-MuRFuncForm","fastNLO-MuFFuncForm"}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOnormDISalt::fRequirements = {
   "QcdnumInit",                // instance of QcdnumInit
   "QcdnumAs",                  // alpha_s evolution from qcdnum ('QcdnumAlphas')
   "PDF",                       // PDF. Could be LHAPDF, or QcdnumPDF if PDF fit is used
   "AemRun",                    // Alpha_em running
   "Filename",                  // fastNLO table
   "fastNLO-ScaleFacMuR",       // scale factor ren. scale
   "fastNLO-ScaleFacMuF",       // scale factor fact. scale
   "fastNLO-iOrd",              // order of calculation of fastNLO table (NCDIS calculation taken from QCDNUM)
   "fastNLO-MuRFuncForm",       // mur functional form
   "fastNLO-MuFFuncForm",       // muf functional form
   "nx","nQ2",                  // x,Q2 integration points
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOnormDISalt::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLOnormDISalt::fFunctionName = "fastNLOnormDISalt"; //!< The function's name

// __________________________________________________________________________________________ //
AfastNLOnormDISalt::AfastNLOnormDISalt(const std::string& name) : AParmFuncBase<double>(name), fastNLOReader() {
   AlposObject::SetClassName("AfastNLOnormDISalt");
   // fastNLO initializations
   fUnits               = fastNLO::kAbsoluteUnits; // valid for H1-jets
   fMuRFunc             = fastNLO::kQuadraticMean;
   fMuFFunc             = fastNLO::kScale1;
   // fastNLO inherited
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   fUseHoppet           = false;
}


// __________________________________________________________________________________________ //
AfastNLOnormDISalt::~AfastNLOnormDISalt() {
}


// ___________________________________________________________________________________________ //
bool AfastNLOnormDISalt::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   AlposObject::debug["Init"]<<endl;
   string filename = PAR_S(Filename);
   // redo fastNLO initialisation
   //   fastNLOBase::SetFilename(filename); // just 'set' it
   ReadTable();
   fastNLOReader::SetFilename(filename); // do some remaining initialization
   // Things which do  not make sense to change later
   // SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
   // if (GetIsFlexibleScaleTable()) {
   //    SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
   //    SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
   // }

   return true;
}


// __________________________________________________________________________________________ //
bool AfastNLOnormDISalt::Update() {  //alpos

   if ( CHECK(Filename) )  fastNLOReader::SetFilename(PAR_S(Filename));
   if ( CHECK(fastNLO-ScaleFacMuR) || CHECK(fastNLO-ScaleFacMuF) )
      SetScaleFactorsMuRMuF(PAR(fastNLO-ScaleFacMuR),PAR(fastNLO-ScaleFacMuF));
   if ( CHECK(fastNLO-iOrd) ) {
      int iOrd = PAR(fastNLO-iOrd);
      //SetContributionON(fastNLO::kFixedOrder , 0 , true);  // LO is always ON
      if ( iOrd==0 ) {
         SetContributionON(fastNLO::kFixedOrder , 1, false);
         SetContributionON(fastNLO::kFixedOrder , 2, false);
      }
      else if ( iOrd==1 ){
         SetContributionON(fastNLO::kFixedOrder , 1,  true);
         SetContributionON(fastNLO::kFixedOrder , 2, false);
      }
      else if ( iOrd==2 ){
         SetContributionON(fastNLO::kFixedOrder , 1,  true);
         SetContributionON(fastNLO::kFixedOrder , 2,  true);
      }
   }

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(PDF);
   UPDATE(QcdnumAs);
   UPDATE(AemRun);

   // get fastNLO cross sections
   CalcCrossSection(); // fastNLOReader::
   fValue = GetCrossSection(); // fastNLOReader::
   fError.resize(fValue.size());

   // -------- NC DIS cross section
   vector<double> q2bins = {150, 200, 270, 400, 700, 5000, 15000 };
   vector<double> CSNCDIS(q2bins.size()-1);
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
   for ( unsigned int ib = 0 ; ib<q2bins.size()-1 ;ib++ ) {
      double q2min = q2bins[ib];
      double q2max = q2bins[ib+1];
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
         CSNCDIS[ib] += CSi;
      }
      //cout<<"CS NCDIS["<<ib<<"]: "<<CSNCDIS[ib]<<endl;
      printf("CS NCDIS: %14.12f\n",CSNCDIS[ib]);
   }


   int nPt = fValue.size()/CSNCDIS.size();
   int j=0;
   for ( unsigned int ib = 0 ; ib<CSNCDIS.size() ;ib++ ) {
      for ( int ip = 0 ; ip<nPt ; ip++ ) {
         fValue[j++]/=CSNCDIS[ib];
      }
   }

   return true;
}

//______________________________________________________________________________


double AfastNLOnormDISalt::EvolveAlphas(double Q) const { //fastNLO
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return QUICK(QcdnumAs, ({Q}) )[0];
}


//______________________________________________________________________________


bool AfastNLOnormDISalt::InitPDF() { //fastNLO
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // nothing todo
   return true;
}


//______________________________________________________________________________



vector<double> AfastNLOnormDISalt::GetXFX(double x, double Q) const { //fastNLO
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return QUICK(PDF, ({x,Q}) );
}


//______________________________________________________________________________
