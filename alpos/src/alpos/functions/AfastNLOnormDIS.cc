// lc DS 11.12.2015

#include "alpos/functions/AfastNLOnormDIS.h"
#include "alpos/functions/AQcdnumInit.h"  // need QCDNUM for NC DIS cross section
#include "fastnlotk/read_steer.h"
#include <iostream>

#include "alpos/functions/ALhapdf6.h"  // need QCDNUM for NC DIS cross section

using namespace std;


// alpha_em running klasen, kramer, et al.
double alphem_KK(double Q2) {
   //include 'common.f'
   //real*8 alphem,alko,aempi,rpigg,scal
   const double alko = 1./137.0359895;       //  ! Electron charge
   double aempi = alko/(3.*M_PI);
   double scal = Q2;
   double rpigg = 0;
   if ( scal < 2.e-6) 
      rpigg=0.;
   else if (scal < 0.09) 
      rpigg= aempi*(13.4916+log(scal))+0.00835*log(1.+scal);
   else if ( scal < 9.) 
      rpigg=aempi*(16.3200+2.*log(scal))+0.00238*log(1.+3.927*scal);
   else if( scal < 1.e4) 
      rpigg=aempi*(13.4955+3.*log(scal))+0.00165+0.00299*log(1.+scal);
   else
      rpigg=aempi*(13.4955+3.*log(scal))+0.00221+0.00293*log(1.+scal);

   
   double alphem=alko/(1.-rpigg);
   return alphem;
}
//c=============================================== 

// need structure function calculation from QCDNUM
extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);


// __________________________________________________________________________________________ //
const std::vector<std::string> AfastNLOnormDIS::fRequirements = {//"Filename",                      // fastNLO table filenmae
							  "PDF",                            // a 'PDF' function with Quick-function. Moslty LHAPDF6
							  "Alpha_s",                        // alpha_s(mu_r) function (from QCDNUM).
							  "ScaleFacMuR","ScaleFacMuF",      // Scale factors for ren. and fact. scale
							  "iOrd",                           // order
							  "iThr",                           // Use threshold corrections if available in fastNLO table
							  "nQ2",                            // # points for Q2 grid
							  "nx",                             // # points for x  grid
							  "AemRun",						         // alpha_em running coupling
							  "QcdnumInit",						   // QcdnumInit singleton
							  "ymin", "ymax", "sqrts",          // inelasticity and sqrt(s)
							  "MuRFuncForm","MuFFuncForm",      // mu_r and mu_f functional form for fastNLO flexible scale tables
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AfastNLOnormDIS::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AfastNLOnormDIS::fFunctionName = "fastNLOnormDIS"; //< The function's name


// __________________________________________________________________________________________ //
AfastNLOnormDIS::AfastNLOnormDIS(const std::string& name) : AParmFuncBase<double>(name) {
   SetClassName("AfastNLOnormDIS");
}


// __________________________________________________________________________________________ //
AfastNLOnormDIS::~AfastNLOnormDIS() {
   //if ( fnlo ) delete fnlo;
   for ( auto i : fnlos ) delete i;
}


// ___________________________________________________________________________________________ //
bool AfastNLOnormDIS::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Requires in datafile defintions  of
   //!  'Filename' (one entry), 'Filenames' (array) or 'Tables' (table)
   //! and further 'Units' and if applicable: 'MuRFuncForm' and 'MuFFuncForm'
   //!

   using namespace AlposTools;

   double Units = DOUBLE_NS(Units,GetAlposName());

   vector<string> filenames;
   vector<int> firstbins, lastbins;
   if ( EXIST_NS(Filename,GetAlposName() ))
      filenames.push_back(STRING_NS(Filename,GetAlposName()));
   else if ( EXIST_NS(Filenames,GetAlposName() ) )
      filenames = STRING_ARR_NS(Filenames,GetAlposName()); // direct access to array
   else if ( EXIST_NS(Tables,GetAlposName()) ) {
      filenames = STRING_COL_NS(FnloTables,Filenames,GetAlposName());
      firstbins = INT_COL_NS(FnloTables,FirstBin,GetAlposName());
      lastbins  = INT_COL_NS(FnloTables,LastBin,GetAlposName());
   }

   // init binmap if specified
   if ( !firstbins.empty() ) {
      for ( unsigned int ig = 0; ig<filenames.size() ; ig++ ){
         vector<bool> bmap(fnlos[ig]->GetNObsBin());
         for ( unsigned int ib=0 ; ib<bmap.size() ; ib++ ) {
            bmap[ib] = (ib>=firstbins[ig] && ib <= lastbins[ig] );
         }
         fBinmap += bmap;
      }
   }


   //string filename = PAR_S(Filename);
   //cout<<"fastnlo: Filename: "<<filename<<endl;
   for ( string& fn : filenames ) {
      AlposTools::CheckFileExit(fn);
      info["Init"]<<"Reading table file: "<<fn<<endl;
      fastNLOAlpos* f = new fastNLOAlpos(fn);
      fnlos.push_back(f);
      f->SetAlposName(this->GetAlposName());
      f->SetUnits(static_cast<fastNLO::EUnits>(Units));
      if ( f->GetIsFlexibleScaleTable()) {
	 f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
	 f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
         // debug["Init"]<<"read DOUBLE_NS(MuRFuncForm,"<<GetAlposName()<<") = "<<DOUBLE_NS(MuRFuncForm,GetAlposName())<<endl;
         // f->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuRFuncForm,GetAlposName())));
         // debug["Init"]<<"     PAR(MuRFuncForm,...) would be "<<PAR(MuRFuncForm)<<endl;
         
         // debug["Init"]<<"read DOUBLE_NS(MuFFuncForm,"<<GetAlposName()<<") = "<<DOUBLE_NS(MuFFuncForm,GetAlposName())<<endl;
         // f->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(DOUBLE_NS(MuFFuncForm,GetAlposName())));
         // debug["Init"]<<"     PAR(MuFFuncForm,...) would be "<<PAR(MuFFuncForm)<<endl;
      }
      info["Init"]<<"Read "<<f->GetNObsBin()<<" rows"<<endl;
   }
   
   
   
   

   /*
   //if ( fnlo ) delete fnlo;
   fnlo = new fastNLOAlpos(filenames[0]);
   fnlo->SetAlposName(this->GetAlposName());
   //Things which do not make sense to change later
   //fnlo->SetFilename(PAR_S(Filename)); // filename should not be changed
   fnlo->SetUnits(static_cast<fastNLO::EUnits>(PAR(Units)));
   if (fnlo->GetIsFlexibleScaleTable()) {
      fnlo->SetMuRFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuRFuncForm)));
      fnlo->SetMuFFunctionalForm(static_cast<fastNLO::EScaleFunctionalForm>(PAR(MuFFuncForm)));
   }
   CONST(MuRFuncForm);
   CONST(MuFFuncForm);
   CONST(Units);
   CONST(Filename);
   */
   CONST(iOrd);
   CONST(iThr);


   //fQ2BinsLoUps
   // --- get q2 binning
   vector<double> q2min = DOUBLE_COL_NS(Data,q2min,GetAlposName());
   vector<double> q2max = DOUBLE_COL_NS(Data,q2max,GetAlposName());
   for ( unsigned int iq = 0 ; iq<q2min.size() ; iq++ ) {
      fQ2BinsLoUps.push_back(make_pair(q2min[iq],q2max[iq]));
   }

   return true;
}


// __________________________________________________________________________________________ //
bool AfastNLOnormDIS::Update() {

   if ( fValue.empty() || (fValue.size()==1 && fValue[0]==0 ) ) SetOrder();

   using namespace AlposTools;

   if ( CHECK(ScaleFacMuR) || CHECK(ScaleFacMuF) )
      for ( auto fnlo : fnlos )
	     fnlo->SetScaleFactorsMuRMuF(PAR(ScaleFacMuR),PAR(ScaleFacMuF));

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit);
   UPDATE(PDF);
   UPDATE(Alpha_s);
   UPDATE(AemRun);

   // get cross sections
   fValue.clear();
   for (int i=0; i<fnlos.size(); i++) {
      debug["Update"]<<"Re-reading fastNLO table #"<<i<<endl;
      fnlos[i]->CalcCrossSection();
      fValue += fnlos[i]->GetCrossSection();
   }
   
   // apply binmap if applicable
   if ( !fBinmap.empty() ) {
      int ii=0;
      for ( int ib = 0 ; ib<fValue.size() ; ib++ ){
         if ( fBinmap[ib] ) fValue[ii++]=fValue[ib];
      }
      fValue.resize(ii);
   }

   // -------- Calculate NC DIS cross section for normalization
   // WIP: get rid of some hard-coded values

   
   // erase previous NC DIS cross sections
   q2LoUpCSmap.clear();
   
   vector <pair<double,double>> binQ2LoUps;
   for ( const pair<double,double>& binsLoUp : fQ2BinsLoUps ) {
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
   /*
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
   */

   // do the normalization
   for (int ip=0; ip<fValue.size(); ip++) {
      debug["Update"]<<"Normalizing CS: fValue["<<ip<<"] = "<<fValue[ip]<<" / "<<q2LoUpCSmap[binQ2LoUps[ip]]<<" ("<<binQ2LoUps[ip].first<<" < Q2 < "<<binQ2LoUps[ip].second<<")"<<endl;
      fValue[ip]/=q2LoUpCSmap[binQ2LoUps[ip]];
   }

   fError.resize(fValue.size());


   /*
   {
      ALhapdf6* lhapdfInstance = (ALhapdf6*) TheoryHandler::Handler()->GetFuncD("LHAPDF6");
      int initialPdfMember = PAR_ANY("LHAPDF6.PDFSet");

      // PDF set meta-information (e.g. Confidence Level, error type)  
      std::string pdfSetName = lhapdfInstance->GetPDFSet()->name();
      unsigned int nMembers = lhapdfInstance->GetPDFSet()->size();
      double setCL = lhapdfInstance->GetPDFSet()->errorConfLevel() / 100.;
      double reqCL = std::erf(1 / sqrt(2));  // required confidence level: 68.2689... % 
      double errScale = 1.0;  // scaling due to set CL != required CL

      
      double binning[9] = {5.5,8,11,16,22,30,42,60,80};
      cout<<"Calcualte PDF uncertainty for NC DIS:"<<endl;
      vector<vector<double> > cs0(8); // [bin][PDFset]
      for ( int ib = 0 ; ib<8 ; ib++ ) cs0[ib].resize(101);
      
      for ( int iPDFset = 0 ; iPDFset<=100 ; iPDFset++ ) {
	 //cs0.resize(cs0.size()+1);
	 cs0.resize(cs0.size()+1);
	 q2LoUpCSmap.clear();
	 SET(PDF.PDFSet,iPDFset,0);
	 UPDATE(PDF);
	 UPDATE(QcdnumInit);
	 //for ( auto binsLoUp : fnlos[0]->GetObsBinsBounds(0) ) {
	 for ( int ib = 0 ; ib<8 ; ib++ ){
	    cs0[ib][iPDFset] = CalcNCDISCrossSection(binning[ib], binning[ib+1]);
	    cout<<"CS (PDFSset "<<iPDFset<<") Q2bin("<<binning[ib]<<"--"<<binning[ib+1]<<"): "<<cs0[ib][iPDFset]<<endl;
	    //cs0.back().push_back( CalcNCDISCrossSection(binsLoUp.first, binsLoUp.second) );
	    //cout<<"CS (PDFSset "<<iPDFset<<") Q2bin("<<binsLoUp.first<<"--"<<binsLoUp.second<<"): "<<cs0.back().back()<<endl;
	    //q2LoUpCSmap[binsLoUp] = cs0.back().back();
	 }
      }
      

      
      int ib=0;
      for ( int ib = 0 ; ib<8 ; ib++ ){
	 LHAPDF::PDFUncertainty err_struct = lhapdfInstance->GetPDFSet()->uncertainty(cs0[ib], 100*reqCL,false);
	 // // TODO: use Alpos averaging, or get symmetrized from LHAPDF? 
	 // errVals[iObs] = err_struct.errsymm;
	 cout<<"bin "<<ib<<"\tcs0="<<cs0[ib][0]<<"\terrsym="<<err_struct.errsymm<<"\te-"<<err_struct.errminus<<"\te+"<<err_struct.errplus<<endl;
      }


   }   
   exit(2);
   */

   return true;
}


// __________________________________________________________________________________________ //
void AfastNLOnormDIS::SetOrder() {
   //! Set correct order of fastNLO calculation
   for ( auto fnlo : fnlos ) {
   // if ( CHECK(iOrd) ) {
      //! Check on existence of various pQCD contributions in table (Id = -1 if not existing)
      //! Check on existence of LO (Id = -1 if not existing)
      int ilo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
      if (ilo < 0) {
         error["SetOrder"] << "LO not found, nothing to be done!" << endl;
         exit(1);
      } else {
         info["SetOrder"] << "The LO contribution has Id: " << ilo << endl;
      }
      //! Check on existence of NLO (Id = -1 if not existing)
      int inlo  = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
      if (inlo < 0) {
         info["SetOrder"] << "No NLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NLO contribution has Id: " << inlo << endl;
      }
      //! Check on existence of NNLO (Id = -1 if not existing)
      int innlo = fnlo->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToNextToLeading);
      if (innlo < 0) {
         info["SetOrder"] << "No NNLO contribution found!" << endl;
      } else {
         info["SetOrder"] << "The NNLO contribution has Id: " << innlo << endl;
      }
      //! Check on existence of threshold corrections
      int ithc1 = fnlo->ContrId(fastNLO::kThresholdCorrection, fastNLO::kLeading);
      int ithc2 = fnlo->ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
      if (ithc1 < 0) {
         info["SetOrder"] << "1-loop threshold corrections not found!" << endl;
      } else {
         info["SetOrder"] << "1-loop threshold corrections have Id: " << ithc1 << endl;
      }
      if (ithc2 < 0) {
         info["SetOrder"] << "2-loop threshold corrections not found!" << endl;
      } else {
         info["SetOrder"] << "2-loop threshold corrections have Id: " << ithc2 << endl;
      }

      //! Switch selected contributions ON, if possible
      //! LO & NLO are ON by default
      //! Fixed-order
      int iOrd = PAR(iOrd);
      if ( iOrd==0 ) {
         if ( !(inlo<0) )  {fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, false);}
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==1 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) ) {fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, false);}
      } else if ( iOrd==2 ) {
         if ( !(inlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, inlo, true);
         } else {
            error["SetOrder"] << "NLO requested, but not found. Nothing to be done!" << endl;
            exit(1);
         }
         if ( !(innlo<0) )  {
            fnlo->SetContributionON(fastNLO::kFixedOrder, innlo, true);
         } else {
            warn["SetOrder"] << "NNLO requested, but not found. Nothing to be done!" << endl;
            //exit(1);
         }
      }
      // Threshold corrections
//      if ( CHECK(iThr) ) {
         int iThr = PAR(iThr);
         if ( !(iThr<0) ) {
            if ( iThr==0 && iOrd==0 ) {
               if ( !(ithc1<0) )  {
                  fnlo->SetContributionON(fastNLO::kThresholdCorrection, ithc1, true);
               } else {
                  error["SetOrder"] << "LO threshold corrections requested, but not found. Nothing to be done!" << endl;
                  exit(1);
               }
            } else if ( iThr==1 && iOrd==1 ) {
               if ( !(ithc2<0) )  {
                  fnlo->SetContributionON(fastNLO::kThresholdCorrection, ithc2, true);
               } else {
                  error["SetOrder"] << "NLO threshold corrections requested, but not found. Nothing to be done!" << endl;
                  exit(1);
               }
            } else {
               error["SetOrder"] << "Inconsistent request for threshold corrections. Nothing to be done!" << endl;
               exit(1);
            }
         }
//      }
//   }
   }
}

//___________________________________________________________________//
double AfastNLOnormDIS::CalcNCDISCrossSection(double q2min, double q2max) {
   
   //vector<double> q2bins = {150, 200, 270, 400, 700, 5000, 15000};
   //vector<double> CSNCDIS(q2bins.size()-1);
   
   double CSNCDIS = 0;

   // //double S = 318.11947*318.11947;
   // double S = 27.6*920*4;
   // double ymin = 0.2;
   // double ymax = 0.6;

   double ymin = PAR(ymin);
   double ymax = PAR(ymax);
   double S    = pow(PAR(sqrts),2);

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
	 //double alem = 1./137.035999074;
         //double alem = aemrun(q2(j))
         //double alem = 1./137.;
	 //double alem = alphem_KK(q2[ip]);;
	 //cout<<"Q="<<sqrt(q2[ip])<<"\talem="<<alem<<endl;

         double factor=2.*M_PI*alem*alem/(x[ip]*q2[ip]*q2[ip])*convfac;
         CSi *= factor* dq2[ip]* dx[ip];
	      CSNCDIS += CSi;
      }
      //cout<<endl;
      cout<<"CS NCDIS: "<<CSNCDIS<<endl;
      //printf("CS NCDIS: %14.12f\n",CSNCDIS[ib]);
      
      return CSNCDIS;
}


//______________________________________________________________________________


double AfastNLOnormDIS::EvolveAlphas(double Q) const { //fastNLO
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return QUICK(Alpha_s, ({Q}) )[0];
}


//______________________________________________________________________________


bool AfastNLOnormDIS::InitPDF() { //fastNLO
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // nothing todo
   return true;
}


//______________________________________________________________________________



vector<double> AfastNLOnormDIS::GetXFX(double x, double Q) const { //fastNLO
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return QUICK(PDF, ({x,Q}) );
}


//______________________________________________________________________________
