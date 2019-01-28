// DB 15.01.2015

#include "QCDNUM/QCDNUM.h"
#include "alpos/functions/AQcdnumInit.h"
#include "alpos/ATheory.h"
#include <cmath>
#include <iostream>
#include "TMatrixD.h"
#include "TSystem.h"

using namespace std;


extern "C" {
   // call qcinit(6,' ');//        !initialize 
   // call setord(I_FIT_ORDER);//         !LO, NLO, NNLO       
   //void qcdnum_ini_();
   // void getintegrateddisxsection_(double*);

   void qcinit_(int*,char*);
   void setalf_(double*,double*);
   void getalf_(double*,double*);
   void setcbt_(int* nfin,int* mc,int* mb,int* mt);
   double asfunc_(double* q2, int* nf, int* ierr);
   double setord_(int* iord);
   void gxmake_(double* xmin,int *iwt,int* nx,int* nin,int* nxout,int* ispline);//        !x-grid 
   //void gqmake_(double* qarr,int* wgt,int* n,int* nqin,int* nqout);//          !mu2-grid  
   void gqmake_(double* qarr,double* wgt,int* n,int* nqin,int* nqout);//          !mu2-grid  
   int iqfrmq_(double* imc);// returns closes grid indices
   void fillwt_(int* itype,int* id1,int* id2,int* nw);
   void zmfillw_(int* nw);
   void zswitch_(int* IPDFSET);
   void hswitch_(int* IPDFSET);

   void eprc_init_(int* boool);
   void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);
   void pdfinp_( void (*) (double*,double*,double*), int*,double*,double*,int*);
   //void pdfext_( void (*) (double*,double*,double*), int*,int*,double*,double*); // deprecated
   void extpdf_(double (*)(int*, double* ,double* ,bool*), int*, int*, double*, double*);
   void zmdefq2_(double* a, double* b);   
   void setabr_(double* a, double* b);
   void evolfg_( int* itype, double (*) (int* ipdf, double* x), double* def, int* iq0, double* epsi );
   void myevolfg_( int* itype, double (*) (int* ipdf, double* x), int* iq0, double* epsi );
   void printdef_(double* def);
   
   void fpdfxq_(int*,double* x, double* q, double* pdf, int*);
}

namespace AWrap {
   void GetXFX(double* x, double* mu2, double* xfx) {
      double q = sqrt(*mu2);
      //cout<<"GetXFXL x,mu: "<<*x<<"\t"<<q<<endl;
      vector<double> xxx = TheoryHandler::Handler()->GetFuncD(AQcdnumInit::fFunctionName+".PDF")->GetQuick(2,*x,q); // never 'specialize' the QcdnumInit instance
      for ( int i = 0 ; i<13 ; i++ ) xfx[i] = xxx[i];
   }

   // adopted from qcdnum/testjobsCxx/pdfsetsCxx.cc
   /*----------------------------------------------------------------------
    * Return q,qbar,g for ipdf = [-6:6], proton singlet [7], nonsinglet [8]
    *----------------------------------------------------------------------
    */
   double fun(int* ii, double* xx, double* qq, bool* fst) {
      // The index ipdf runs from -6 to 6 , indexed according to (5.2) in the PDG convention. 
      // For the extra pdfs, if any, the index is 7 ≤ pdf ≤ +n.  
      // The first time func is called by extpdf the flag first is set to true
      // and you can initialise the function, if needed.
      static double proton[13];
      double q = sqrt(*qq);
      vector<double> xxx = TheoryHandler::Handler()->GetFuncD(AQcdnumInit::fFunctionName+".PDF")->GetQuick(2,*xx,q); // never 'specialize' the QcdnumInit instance
      int ipdf = *ii;
      
      // bool first = *fst;
      // if(first) {
      // 	 for(int i=1; i<=13; i++) { QCDNUM::qstore("Read",i,proton[i-1]); }
      // }
      double f = 0;
      //if(ipdf >= -6 && ipdf <= 6) f = QCDNUM::fvalxq(1,ipdf,*xx,q,1);
      // if(ipdf == 7)               f = QCDNUM::sumfxq(1,proton,2,*xx,q,1);
      // if(ipdf == 8)               f = QCDNUM::sumfxq(1,proton,3,*xx,q,1);
      if(ipdf >= -6 && ipdf <= 6) f = xxx[ipdf+6];//QCDNUM::fvalxq(1,ipdf,*xx,q,1);
      else {
	 cout<<" ERROR ! AQcdnumInit.cc"<<endl;
	 exit(3);
      }
      return f;
   }

   double pdfinputQCDNUMExample(int* ipdf, double* xptr){
      using namespace std;
      double x = *xptr;
      if(*ipdf== 0) return 1.7 * pow(x,-0.1) * pow(1.-x,5);
      if(*ipdf== 1) return 3.064320*pow(x,0.8) * pow(1.-x,4.);
      if(*ipdf== 2) return 5.107200*pow(x,0.8)* pow(1.-x,3.);
      if(*ipdf== 3) return 0.;
      if(*ipdf== 4) return  0.1939875 * pow(x,-0.1) * pow(1.-x,6.);
      if(*ipdf== 5) return (0.1939875 * pow(x,-0.1) * pow(1.-x,6.))*(1-x);
      if(*ipdf== 6) return 0.2 * ((0.1939875 * pow(x,-0.1) * pow(1.-x,6.))+((0.1939875 * pow(x,-0.1) * pow(1.-x,6.))*(1-x)));
      if(*ipdf== 7) return 0.;
      if(*ipdf== 8) return 0.;
      if(*ipdf== 9) return 0.;
      if(*ipdf==10) return 0.;
      if(*ipdf==11) return 0.;
      if(*ipdf==12) return 0.;
      return 0;
   }

   double pdfinput(int* ipdf, double* xptr){
      //! pointer to AQcdnumInput object has passed once with 
      //!    AWrap::pdfinput(&(-1),(double*)this); // pass this pointer to pdfinut wrapper
      //! Then, always AQcdnumInput is called for defintion of PDFs.
      //! For more flexibility, it also possible to specify
      //! one specific member function of AQcdnumInit. Call:
      //!    AWrap::pdfinput(&(-1),(double*)this); // pass this pointer to pdfinut wrapper
      //!    AWrap::pdfinput(&(-2),(double*)<pointer-to-memberfunction>(int,double));
      using namespace std;
      static AQcdnumInit* AQI = NULL;
      typedef double (AQcdnumInit::*PDFInp) (int,double) ;
      static PDFInp* PDFI = NULL;
      // pass pointer to AQcdnumInit with ipdf=-1
      if(*ipdf==-1) {
	 AQI = (AQcdnumInit*)xptr;
	 return 0;
      }
      else if ( *ipdf==-2 ) {
	 PDFI = (PDFInp*)xptr;
	 return 0;
      }
      // check if AQI pointer was initialized
      if ( !AQI && !PDFI ) {
	 cout<<"ERROR. AWrap::pdfinput. AQcdnumInit not set!"<<endl;
	 exit(1);
      }
      if ( AQI && PDFI==NULL ) return AQI->pdfinput(*ipdf,*xptr);
      if ( AQI && PDFI ) return (AQI->*(PDFInp) *PDFI )(*ipdf,*xptr);
      return 0;
   }
}

// __________________________________________________________________________________________ //
const std::vector<std::string> AQcdnumInit::fRequirements = {"AlphasMz",
							     "Mz", // boson mass
							     "iOrd", // order, LO, NLO, NNLO
							     "PDFQ0Param", // function for PDF parameterization
							     "nxGrid", // number of x-grid points
							     "nqGrid", // number of q2-grid points
							     "Q0", // starting scale for evolution
							     "InitEvolution", // use QCDNUM for PDF evolution (1: yes, 0: no)
							     "PDF", // if InitEvolution=0: specify external PDF-function
							     "mcharm","mbottom","mtop", // quark masses
							     "nfFix", // number of flavors (nfFix=0: ZMVFNS)
							     "ScaleFacMuR", // renormalization scale factor
							     "ScaleFacMuF", // factorization scale factor
}; //< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumInit::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumInit::fFunctionName = "QcdnumInit"; //< The function's name
int AQcdnumInit::fNinstances = 0;

// __________________________________________________________________________________________ //
AQcdnumInit::AQcdnumInit(const std::string& name) : AParmFuncBase<double>(name) { 
   // Remember: no access to parameters possible in constructor!
   //ARegisterRequirements(this); // needed in every constructor
   SetClassName("AQcdnumInit");
   if ( ++fNinstances > 1 ) {
      cout<<"Error. AQcdnumInit::AQcdnumInit(). Only on QcdnumInit function is allowed."<<endl;
      exit(1);
   }
   fValue.resize(1);
   fError.resize(1);
   // qcdnum itself
   int lun = 6;// qcdnum output on standard screen
   char nix[] = " ";
   qcinit_(&lun,nix);
}


// __________________________________________________________________________________________ //
AQcdnumInit::~AQcdnumInit() {
}


// __________________________________________________________________________________________ //
int AQcdnumInit::SetXGrid(vector<double> xgrid, vector<int> wgt, int nNodes){
   //! Set x grid in qcdnum
   //! return number of x nodes
   int nxout = -1;
   int ispline = 2 ;
   int nx = int(xgrid.size());
   gxmake_(&xgrid[0],&wgt[0],&nx,&nNodes,&nxout,&ispline);
   info["SetXGrid"]<<"Made x grid with "<<nxout<<" grid points."<<endl;
   return nxout;
}


// __________________________________________________________________________________________ //
int AQcdnumInit::SetQGrid(vector<double> qarr, vector<double> wgt, int nNodes){
   //! Set Q grid in qcdnum
   //! return number of Q nodes
   int nq = int(qarr.size());
   int nxout = -1;
   gqmake_(&qarr[0],&wgt[0],&nq,&nNodes,&nxout);//  muf grid
   info["SetQGrid"]<<"Made muf grid with "<<nxout<<" grid points."<<endl;
   return nxout;
}


// ___________________________________________________________________________________________ //
bool AQcdnumInit::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   // 'real' QCDNUM init is called in constructor

   debug["Init"]<<"AQcdnumInit::Init()"<<endl;
   info["Init"]<<"The usage of QCDNUM always requires, that the QCDNUM internal alpha_s evolution is used for consistency reason."<<endl;
  
   // --- set order
   int iOrd = PAR(iOrd)+1;
   setord_(&iOrd);

   // --- read some input values
   int finitevol = PAR(InitEvolution);
   fInitEvol = finitevol == 1 ? true : false;
   fQ0 = PAR(Q0); // Q0 should be below m_c

   // --- pass this pointer to pdfinut wrapper
   int mone=-1;
   AWrap::pdfinput(&mone,(double*)this);


   // --- init x and muf grids to set nf in alpha_s
   vector<double> xmin = {9.9e-7, 0.01e0, 0.10e0, 0.40e0, 0.70e0};
   vector<int>    iwgt = {1,2,4,8,16};
   int nin = PAR(nxGrid); // number of x-grid points
   SetXGrid(xmin, iwgt, nin);


   // --- check Q0 and mc input:
   if ( PAR(Q0) > PAR(mcharm) ) {
      error["Init"]<<"Staring scale 'Q0' must be below charm mass thrshold 'mcharm'. Evolution not possible. Exiting."<<endl;
      exit(1);
   }

   // --- init muf grid
   // // HERAFitter:
   // vector<double> qarr = {1.0  ,      1.8999999761581421   ,    1.959999918937683    ,   22.56250 ,   30276.0, 205000000.0};
   // vector<double> dwgt = {1,4,2,1.5,1,1};
   map<double,double> q2map = {{1,1},{205000000.0,1}}; // ordered map for nodes and weights
   q2map[PAR(Q0)*PAR(Q0)] = 4;
   if ( abs(PAR(mcharm)-PAR(Q0)) > 1.e-3 )
      q2map[PAR(mcharm)*PAR(mcharm)] = 2;
   q2map[PAR(mbottom)*PAR(mbottom)] = 1.5;
   q2map[PAR(mtop)*PAR(mtop)] = 1;

   vector<double> qarr, dwgt;
   for ( auto iq : q2map ) {
      qarr.push_back(iq.first);
      dwgt.push_back(iq.second);
   }
   int nqin = PAR(nqGrid);
   SetQGrid(qarr,dwgt,nqin);
      
   // --- calculate and fill weight tables
   int itype=0,id1,id2,nw;
   fillwt_(&itype,&id1,&id2,&nw); // calculate weight table
   zmfillw_(&nw); // fill weight table

   if( PAR(nfFix) >= 3) { //needed only for FFNS
       //Fill Structure functions for Charm & Bottm
       double aq2 = 1., bq2 = 0.; //q2 = a*muF^2 + b, or muF = 1/a*q2 - b/a
       int nwords; //output
       double hqmass[] = { PAR(mcharm),  PAR(mbottom), 0. };
       int iF2FLsets = 3; //fill F2 & FL

       //Read weights from file to boost the execution
       // Try to read the weight file and create one if that fails
       string dirName = Alpos::Current()->Settings()->Alpos_dir+"/temp";
       gSystem->mkdir(dirName.c_str(), true);
       string hqname = dirName + "/hqstf.wgt";
       int ilun = 22;
       int ierr;
       //Try to read the existing grid
       QCDNUM::hqreadw(ilun, hqname,  nwords, ierr);//nwords,ierr are outputs
       if(ierr == 0) {  //if reading OK, do some additional consistency checks
           double a, b;
           double qmas[3];
           QCDNUM::hqparms(qmas, a, b);
           if(qmas[0] != hqmass[0]) ierr = 1;
           if(qmas[1] != hqmass[1]) ierr = 1;
           if(qmas[2] != hqmass[2]) ierr = 1;
           if(a != aq2)             ierr = 1;
           if(b != bq2)             ierr = 1;
       }
       //if if grid reading failed, fill and write grid again
       if(ierr != 0) {

           //  define relation between muf2 and Q2
           //    q2 = a*muf2 + b
           //  default: a=1, b=0
           //  can only be varied if muf==mur
           bq2 = 0;
           aq2 = 1.0/pow(PAR(ScaleFacMuF),2);
           QCDNUM::hqfillw(iF2FLsets, hqmass, aq2, bq2, nwords);//nwords is output
           QCDNUM::hqdumpw(22, hqname);
       }

       cout << "In total nwords "<< nwords << " read" << endl;
   }

   // --- set alpha_s MZ (we must use the internal alpha_s code)
   //     The internal alpha_s routine is available through 'AQcdnumAlphas'
   Setcbt();

   CONST(InitEvolution);
   CONST(nxGrid);
   CONST(nqGrid);

   ios::sync_with_stdio(); // to not mix FORTRAN and C++ output

   return true;
}


// __________________________________________________________________________________________ //
bool AQcdnumInit::Update() {
   static int cc=0;
   debug["Update"]<<"Call: "<<cc++<<endl;

   // --- set cbt, asmz, ord
   Setcbt();
   SetAsMz();
   int iOrd = PAR(iOrd)+1;
   setord_(&iOrd);

   // factorization scale
   double bf = 0;
   double af = 1.0/pow(PAR(ScaleFacMuF),2);
   //  define relation between muf2 and Q2
   //    q2 = a*muf2 + b
   //  default: a=1, b=0
   //  can only be varied if muf==mur
   zmdefq2_(&af, &bf);   

   // renormalization scale
   double br = 0;
   double ar = pow(PAR(ScaleFacMuR)/PAR(ScaleFacMuF),2);
   // define the relation between the factorisation scale muf and q: 
   //    mur2=a*muf2+b
   // default: ar=1, br-0
   setabr_(&ar, &br);


   if ( !fInitEvol ) {
      // --- update without evolution
      debug["Update"]<<"Update should not be called repeatedly, since no predictions are provided in 'InitEvolution=0' mode."<<endl;
      // --- update pdf
      PAR(PDF); // 'init' or update PDF

      // set external PDF (this requires access to the PDF, which is not yet initialized during the 'init' step)
      double epsi; // Maximum deviation of the quadratic spline interpolation from linear interpotion
      int IPDFSet = 5 ;//external PDF (5-24)
      double offset = 0.001; // to catch discontinuities at the thresholds
      int nExtraPDF = 0; //Number of extra pdfs to be imported, beyond the 13 quark and gluon pdfs
      int nwds;
      //pdfinp_(AWrap::GetXFX, &IPDFSet, &offset,&epsi ,&nwds);
      //pdfext_(AWrap::GetXFX, &IPDFSet, &nExtraPDF, &offset,&epsi); // QCDNUM earlier versions
      extpdf_(AWrap::fun, &IPDFSet, &nExtraPDF, &offset,&epsi);
      info["Update"]<<"PDF imported. Max deviation of linear w.r.t. quadratic spline: "<<epsi<<endl;
      zswitch_(&IPDFSet);
      hswitch_(&IPDFSet);
   }
   else {
      // Update for InitEvolution
      int IPDFSet = 1 ;//external PDF: 5
      zswitch_(&IPDFSet);
      hswitch_(&IPDFSet);
      if(PAR(nfFix) >= 3)
          QCDNUM::hswitch(IPDFSet);
      // double offset = 0.001;
      // int nwds;
      // pdfinp_(AWrap::GetXFX, &IPDFSet, &offset,&epsi ,&nwds);

      SET(PDFQ0Param.iPDF,-3,0); // set PDFQ0Param to 'def' mode.
      //SET(PDFQ0Param.Q0,fQ0,0); 
      vector<double> def = VALUES(PDFQ0Param); // get 'def'

      /*
      TMatrixD mat(13,13);
      assert(def.size() == 13*13);
      cout << "Control output " << endl;
      for(int i = 0; i < def.size(); ++i) {
          int ii = i / 13;
          int jj = i % 13;
          mat(ii,jj) = def[i];
      }
      cout << mat.Determinant() << endl;
      */
      //for(auto  d : def)
          //cout << d <<" "<< endl;

      // cout<<"PDF param def:";
      // for ( auto i : def ) cout<<"\t"<<i;
      // cout<<endl;

      double epsi;
      double q02 = fQ0*fQ0;
      int iq0 = iqfrmq_(&q02);
      int iset = 1; //1: un-polarised, 2: polarised, 3: custom, ...
      
      SET(PDFQ0Param.iPDF,0,0); // set PDFQ0Param to 'gluon' for update.
      PAR(PDFQ0Param); // update PDFQ0

      debug["Update"]<<"Calling evolfg_"<<endl;
      evolfg_( &iset, AWrap::pdfinput, &def[0], &iq0, &epsi );
	 

      fValue[0] = 1;
      fError[0]=1;
   }

   return true;
}


// __________________________________________________________________________________________ //
double AQcdnumInit::pdfinput(int ipdf, double xp){
   //! PDFs at starting scale and pass to QCDNUM

   // SET(PDFQ0Param.iPDF,ipdf,0);
   // SET(PDFQ0Param.xp,xp,0);
   // return PAR(PDFQ0Param); 
   return QUICK(PDFQ0Param, ({double(ipdf),xp,fQ0}) )[0];
}


// __________________________________________________________________________________________ //
void AQcdnumInit::SetAsMz(){
   //! set alpha_s at mz
   double Mz2 = PAR(Mz);
   Mz2*=Mz2;
   double AsMz = PAR(AlphasMz);
   setalf_(&AsMz,&Mz2);
}


// __________________________________________________________________________________________ //
void AQcdnumInit::Setcbt(){
   // call setcbt with recent parameters.
   double mc2 = pow(PAR(mcharm),  2);
   double mb2 = pow(PAR(mbottom), 2);
   double mt2 = pow(PAR(mtop),    2);
   int imc = iqfrmq_(&mc2);
   int imb = iqfrmq_(&mb2);
   int imt = iqfrmq_(&mt2);
   int nfin  = PAR(nfFix);
   setcbt_(&nfin,&imc,&imb,&imt);
}

// __________________________________________________________________________________________ //

