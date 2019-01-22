#include "alpos/functions/AApfel.h"
#include "APFEL/APFEL.h"
#include <iostream>


using namespace std;
//extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf);

extern "C" void externalsetapfel_(double* x, double* q, double* xf);

// function to be called by APFEL
void externalsetapfel_(double* x, double* q, double* xf) {

   if (*x >= 1 || *x < 0) {
      //say::warn["externalsetapfel"]<<"x greater than 1. Unphysical. Returning zero. x="<<*x<<endl;
      for ( int i = 0 ; i<14 ; i++ )  xf[i] = 0;
      return;
   }

   double q0 = *q;//xf[0]; 
   static vector<double> def = {0,0,0,0,0,0,1,0,0,0,0,0,0};
   static vector<bool> ignore(13); //lico
   static TMatrixD InvLiCo(13,13);
   static vector<bool> found(13); // parton
   static vector<bool> occupied(13); //lico
   // --- get LiCo and InvLiCo
   if ( def.size()==13 ) { // get def in first initalization
      SET_ANY(AApfelInit::fFunctionName+".PDFQ0.iPDF",-1,0);
      AlposTools::operator+=(def, VALUES_ANY(AApfelInit::fFunctionName+".PDFQ0"));
      SET_ANY(AApfelInit::fFunctionName+".PDFQ0.iPDF",6,0); //reset
      int idef = 0;
      for ( int ip = 0 ; ip<13 ;ip++ )  found[ip]=false;
      for ( int ip = 0 ; ip<13 ;ip++ )  occupied[ip]=false;
      for ( int i = 0 ; i<13 ; i++ )  {
	 ignore[i] = true;
	 for ( int ip = 0 ; ip<13 ; ip++ )  {
	    if ( def[idef] != 0 )  found[ip] = true;
	    if ( def[idef] !=0 )   occupied[i] = true;
	    if ( def[idef++] != 0 ) ignore[i]=false;
	 }
      }
      for ( int ip = 0 ; ip<13 ; ip++ ) {
	 if ( found[ip] ) continue;
	 else {
	    //cout<<"Parton "<<ip<<" not found."<<endl;
	    if (occupied[13]){
	       say::error["ApfelInit::externalsetapfel"]<<"Inconsistent number of linear combinations."<<endl;
	       exit(1);
	    }
	    for ( int i1 = 0 ; i1<13 ; i1++ )  {
	       if ( !occupied[i1] ) {
		  def[i1*13+ip] = 1;
		  occupied[i1] = true;
		  break;
	       }
	    }
	 }
      }
      if (!occupied[12]){
	 say::error["ApfelInit::externalsetapfel"]<<"Inconsistent number of linear combinations. Too little linear combinations for specified partons."<<endl;
	 exit(1);
      }
      
      // --- InvLiCo
      TMatrixD LiCo(13,13,&def[0]);
      InvLiCo = AlposTools::InvertLU(LiCo);

      // ---- test LiCo
      vector<double> test13 = {1,2,3,4,5,6,7,8,9,10,11,12,13};
      vector<double> test6(13);
      AlposTools::CalcLicoFrom13partons(test6,def,test13);
      vector<double> testxf(13);
      AlposTools::Calc13partonsFromLico(testxf,InvLiCo,test6,found);
      say::info["externalsetapfel_"]<<"Test LiCo (should be 1-13, and '.' if parton is not considered at starting scale): \n             ";
      for ( int ip = 0 ; ip<13 ; ip++ ) {
	 if ( testxf[ip] )cout<<testxf[ip]<<"  ";
	 else cout<<".  ";
      }
      cout<<endl;
      for ( int ip = 0 ; ip<13 ; ip++ ) {
	 if ( testxf[ip] && testxf[ip]!=test13[ip] ) {
	    say::error["externalsetapfel_"]<<"PDF linear combination seems to be unreasonable."<<endl;
	    exit(3001);
	 }
      }
      // ----
   }


   // --- get pdf LiCo
   vector<double> xf6(13);
   for ( int i = 0 ; i<13 ; i++ )  {
      if ( occupied[i] ) { // skip unneeded calls
	 vector<double> ixq = {double(i),*x,q0};
	 xf6[i] = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDFQ0")->GetQuick(ixq)[0];
	 // SET_ANY(AApfelInit::fFunctionName+".PDFQ0.xp",*x,0);
	 // SET_ANY(AApfelInit::fFunctionName+".PDFQ0.iPDF",i,0);
	 // xf6[i] = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDFQ0")->GetValues()[0];
      }
   }

   // --- get xfx13
   vector<double>xf13(13);
   AlposTools::Calc13partonsFromLico(xf13,InvLiCo,xf6,found);
   
   // copy array
   for ( int i = 0 ; i<13 ; i++ )  xf[i] = xf13[i];
   xf[14]=0;

   //if (*x==0.1) cout<<"x="<<*x<<"\tg: "<<xf[6]<<"\tgB: "<<PAR_ANY("PDFQ0_HERA.gC")<<"\tgC: "<<PAR_ANY("PDFQ0_HERA.gC")<<endl;
   return ;

   // --- test with LHAPDF
   vector<double> xforig = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDFQ0.PDF")->GetQuick(2,*x,q0);
   for ( int i = 0 ; i<13 ; i++ )  {
      vector<double> ixq = {double(i),*x,q0};
      xf6[i] = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDFQ0")->GetQuick(ixq)[0];
      cout<<"xf13: "<<xf13[i]<<"\torig: "<<xforig[i]<<endl;
   }
   
   // --- simple usage 
   /*
   vector<double> xxx = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDF")->GetQuick(2,*x,q0); // never 'specialize' the ApfelInit instance

   for ( int i = 0 ; i<13 ; i++ )  xf[i] = xxx[i];
   xf[14]=0;
   return;
   */

   // --- Apfel example:
   /*
   // -------------------- 
   cout<<"q0="<<q0<<"\tx="<<*x<<"\tq="<<*q<<"\tg="<<xf[6]<<"  \tu="<<xf[5]<<","<<xf[7]<<"  \td="<<xf[4]<<","<<xf[6]<<"  \ts="<<xf[3]<<","<<xf[7]<<endl;

   // APFEL example TabulationExtern.f
   double N_uv = 5.107200;
   double    auv  = 0.8;
   double   buv  = 3;
   double   N_dv = 3.064320;
   double   adv  = 0.8;
   double   bdv  = 4;
   double   N_g  = 1.7;
   double    ag   = -0.1;
   double    bg   = 5;
   double    N_db = 0.1939875;
   double    adb  = -0.1;
   double    bdb  = 6;
   double    fs   = 0.2;
// *
// *     User defined PDFs
// *
   double    xuv   = N_uv * pow(*x,auv) * pow(( 1 - *x ),buv);
   double    xdv   = N_dv * pow(*x,adv) * pow(( 1 - *x ),bdv);
   double    xg    = N_g  * pow(*x,ag)  * pow(( 1 - *x ),bg) ;
   double    xdbar = N_db * pow(*x,adb) * pow(( 1 - *x ),bdb);
   double    xubar = xdbar * ( 1 - *x );
   double    xs    = fs * ( xdbar + xubar );
   double    xsbar = xs;

   for ( int i = 0 ; i<14 ; i++ )  xf[i] = 0;
   xf[3+6]  = xs;
   xf[2+6]  = xuv + xubar;
   xf[1+6]  = xdv + xdbar;
   xf[0+6]  = xg;
   xf[-1+6] = xdbar;
   xf[-2+6] = xubar;
   xf[-3+6] = xsbar;
   cout<<"x = " <<*x<<"\t\tgluon: " <<xf[0+6]<<endl;
   */
}


const std::vector<std::string> AApfelInit::fRequirements = {
   "iOrd", // perturbative order of the evolution
   "Q0", //
   "nf", //  Fixed-Flavour Number Scheme with 'nf' active flavours.
   "PDFQ0", // the PDF at the starting scale
   "mZ", "mW", //boson masses
   "Gf", // Fermi constant
   "MassScheme", // mass scheme to be used to compute the structure functions ('ms' = 'ZM-VFNS', 'FFNS', 'FONLL-A', 'FONLL-B', 'FONLL-C')
   // "mp" // proton mass
   // CKM elements
   "SinThetaW" ,//the value of sin(theta_W) [0.23126]
   "TargetMassCorrections", //enables or disables the target mass corrections to the DIS structure functions due to the finite mass of the proton
   "theory", // theory to be used in the evolution. (theory = 'QCD','QED','QUniD', default 'QCD')
   "QCDRefAs","QCDRefQRef", // reference values of alphas at the scale 'Qref' [GeV] to 'alpharef' (default 'alpharef' = 0.35, 'Qref' = sqrt(2) GeV)
   "QEDRefAem","QEDRefQRef", // reference values of alpha at the scale 'Qref' [GeV] to 'alpharef' (default 'alpharef' = 7.496252d-3, 'Qref' = 1.777 GeV).
   "LambdaQCDRef","LambdaQCDnref", // LambdaQCD [GeV] with 'nref' flavours to 'lambdaref' (default 'lambdaref' = 0.220, 'nref' = 5)
   "mc","mb","mt", // heavy quark thresholds  in GeV in the Pole-mass scheme
   "mcMSbar","mbMSbar","mtMSbar", // heavy quark thresholds [GeV] in the MSbar scheme.
   "mtau", // mtau
   "MaxFlavourAlpha", // maximum number of active flavours in the couplings evolution (including the masses)
   "MaxFlavourPDFs", // maximum number of active flavours in the PDF evolution (default 'nf' = 6)
   //   "TimeLikeEvolution", // time-like evolution (frag. functions) 
   "RenFacRatio", //         1
   "SmallxResummation", //  'LL', 'NLL' or 0 to disable Small x resummation
   "AlphaEvolution", //solution of the beta function equations for the running couplings ('evol' = 'exact','expanded','lambda')
   "PDFEvolution", // sets the solution of the DGLAP equations for PDFs ('evolp' = 'exactmu', 'expandalpha', 'expandalpha') 
   "QLimitMin", "QLimitMax",// range where it is possible to perform the evolution
   "EnableMassRunning", //enables/disables the running of the MSbar masses
   "FastEvolution", // enable fast evolution
   "nGridPts0","nGridPts1","nGridPts2", // number of grid points for the three sub-grids (other grid paramters are hardcoded)
   "EnableLeptonEvolution", // enables the evolution of the lepton PDFs when the fast QUniD is used (default false)
}; //< List of all AParm's which this function depends on

const std::vector<std::string> AApfelInit::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AApfelInit::fFunctionName = "ApfelInit"; //< The function's name


// __________________________________________________________________________________________ //
AApfelInit::AApfelInit(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AApfelInit");
}


// __________________________________________________________________________________________ //
AApfelInit::~AApfelInit() {
}


// ___________________________________________________________________________________________ //
bool AApfelInit::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.

/**
   Get help for APFEL:
      $ apfel-config --list-funcs
*/   

   APFEL::SetPerturbativeOrder(PAR(iOrd));
   APFEL::SetTheory(PAR_S(theory));
   int nf = PAR(nf);
   if ( nf==0 ) APFEL::SetVFNS();
   else APFEL::SetFFNS(nf);

   /// hard coded limits
   APFEL::SetQLimits(0.1,1e6); // 0.4Gev to 100TeV
   /// set grids with mostly hard-coded parameters
   if ( PAR(nGridPts2) > 0 )  APFEL::SetNumberOfGrids(3);
   else APFEL::SetNumberOfGrids(2);
   APFEL::SetGridParameters(1,PAR(nGridPts0),3,5.e-7);
   APFEL::SetGridParameters(2,PAR(nGridPts1),5,1.e-1);
   if ( PAR(nGridPts2) > 0 )  APFEL::SetGridParameters(3,PAR(nGridPts2),5,8.e-1);

   CONST(iOrd);
   CONST(theory);
   CONST(nf);
   CONST(nGridPts0);
   CONST(nGridPts1);
   CONST(nGridPts2);

   APFEL::SetMassScheme(PAR_S(MassScheme));
   CONST(MassScheme);

   if ( PAR_S(SmallxResummation) == "LL" || PAR_S(SmallxResummation)=="NLL" ) {
      APFEL::SetSmallxResummation(true,PAR_S(SmallxResummation));
      cout<<"SmallxResummation not yet fully working in APFEL."<<endl;
      cout<<"Don't use this option."<<endl;
      exit(1);
   }

   // Initializes integrals on the grids
   //APFEL::InitializeAPFEL();
   // APFEL::InitializeAPFEL_DIS();

   APFEL::SetPDFEvolution(PAR_S(PDFEvolution));
   APFEL::SetFastEvolution(PAR(FastEvolution));
   //APFEL::EnableLeptonEvolution(PAR(EnableLeptonEvolution));
   APFEL::SetPDFSet("external");
   APFEL::SetAlphaEvolution(PAR_S(AlphaEvolution));
   APFEL::SetMaxFlavourAlpha(PAR(MaxFlavourAlpha));
   APFEL::SetMaxFlavourPDFs(PAR(MaxFlavourPDFs));
   APFEL::EnableMassRunning(PAR(EnableMassRunning));

   CONST(PDFEvolution);
   CONST(FastEvolution);
   CONST(AlphaEvolution);
   CONST(MaxFlavourAlpha);
   CONST(MaxFlavourPDFs);
   CONST(EnableMassRunning);
   // // Set q0 to evolution evolution
   // double eps = 1e-10;
   // double Q0 = PAR(Q0) - eps;
   // double val = -1;
   // externalsetapfel_(&val,&val,&Q0);

   return true;
}

// __________________________________________________________________________________________ //
bool AApfelInit::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   if ( CHECK(RenFacRatio) || CHECK(mc) || CHECK(mb) || CHECK(mt)
   	|| CHECK(mcMSbar) || CHECK(mbMSbar) || CHECK(mtMSbar)   ){ // called always at least once
      APFEL::SetPoleMasses(PAR(mc),PAR(mb),PAR(mt));
      //APFEL::SetMSbarMasses(PAR(mcMSbar),PAR(mbMSbar),PAR(mtMSbar));
      APFEL::SetRenFacRatio(PAR(RenFacRatio));
      APFEL::InitializeAPFEL_DIS();
      APFEL::EnableWelcomeMessage(false);
      APFEL::SetPerturbativeOrder(PAR(iOrd));
   }

   APFEL::SetAlphaQCDRef(PAR(QCDRefAs), PAR(QCDRefQRef));
   APFEL::SetAlphaQEDRef(PAR(QEDRefAem),PAR(QEDRefQRef));
   APFEL::SetLambdaQCDRef(PAR(LambdaQCDRef),int(PAR(LambdaQCDnref)));
   APFEL::SetTauMass(PAR(mtau));

   APFEL::SetGFermi(PAR(Gf));
   APFEL::SetSin2ThetaW(PAR(SinThetaW));
   APFEL::SetWMass(PAR(mW));
   APFEL::SetZMass(PAR(mZ));
   
   
   UPDATE(PDFQ0); // important for 'quick access' in externalsetapfel



   // double eps = 1e-10;
   // double Q0 = sqrt(2) - eps;//PAR(Q0) - eps;
   // double Q  = 100;//PAR(muf);
   //APFEL::EvolveAPFEL(Q0,Q);

   /*
   // Tabulate PDFs for the LHA x values
   double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
		    1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
   cout << "Q0   = " <<Q0<<endl;
   cout << "Q    = " <<Q<<endl;
   cout << "alpha_QCD(mu2F) = " << APFEL::AlphaQCD(Q) << endl;
   cout << "alpha_QED(mu2F) = " << APFEL::AlphaQED(Q) << endl;
   cout << endl;
   cout << "   x   "
	<< setw(11) << "   u-ubar    "
	<< setw(11) << "   d-dbar    "
	<< setw(11) << "  2(ubr+dbr) "
	<< setw(11) << "   c+cbar    "
	<< setw(11) << "   gluon     "
	<< setw(11) << "   photon    " << endl;
   cout << scientific;
   for (int i = 2; i < 11; i++)
      cout << setprecision(1)
	   << xlha[i] << "\t" << setprecision(4)
	   << setw(11) <<  APFEL::xPDF(2,xlha[i]) - APFEL::xPDF(-2,xlha[i]) << "  "
	   << setw(11) <<  APFEL::xPDF(1,xlha[i]) - APFEL::xPDF(-1,xlha[i]) << "  "
	   << setw(11) << 2*(APFEL::xPDF(-1,xlha[i]) + APFEL::xPDF(-2,xlha[i])) << "  "
	   << setw(11) <<  APFEL::xPDF(4,xlha[i]) + APFEL::xPDF(-4,xlha[i]) << "  "
	   << setw(11) <<  APFEL::xPDF(0,xlha[i]) << "  "
	   << setw(11) <<  APFEL::xgamma(xlha[i]) << "  "
	   << endl;
   */

   // double xf=1.e-2;
   // double qtst = 100;
   // vector<double> xxx = TheoryHandler::Handler()->GetFuncD(AApfelInit::fFunctionName+".PDF")->GetQuick(2,xf,qtst); // never 'specialize' the ApfelInit instance
   // APFEL::EvolveAPFEL(Q0-eps,qtst);
   // vector<double> ret(13);
   // for (int i = 0; i < 13; i++) {
   //    cout<<i<<"\tLHAPDF: "<<xxx[i]<<"\tAPFEL: "<< APFEL::xPDF(i-6,xf) <<endl;
   // }

   fValue = {1};
   fError = {0};

   return true;
}

