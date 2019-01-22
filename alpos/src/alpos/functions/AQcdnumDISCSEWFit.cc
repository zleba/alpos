// DB 16.02.2016

#include "alpos/functions/AQcdnumDISCSEWFit.h"
#include "alpos/functions/AQcdnumInit.h"
#include <iostream>
#include "fastnlotk/read_steer.h"
#include "TGraph.h"
#include <map>

using namespace std;

extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);
extern "C" void eprc_get_mw_(double* mw0);
extern "C" void eprc_set_mw_(double* mw0);
extern "C" void eprc_get_mt_(double* mt0);
extern "C" void eprc_set_mt_(double* mt0);
extern "C" double eprc_sw2_cm_(double* q2);


// __________________________________________________________________________________________ //
const std::vector<std::string> AQcdnumDISCSEWFit::fRequirements = {"QcdnumInit",
								   "Mz","Mw", // boson masses
								   "MwProp", // boson masses
								   "Gf", // gf
								   "sin2thw", // sin2th_w effective
								   "EPRC",
								   "e-polarity", // polarity
								   //"e-charge", // charge
								   "I3u","I3d","I3Ru","I3Rd",
								   "au","ad","vu0","vd0", // boson-quark couplings
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumDISCSEWFit::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumDISCSEWFit::fFunctionName = "QcdnumDISCSEWFit"; //!< The function's name

// __________________________________________________________________________________________ //
AQcdnumDISCSEWFit::AQcdnumDISCSEWFit(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AQcdnumDISCSEWFit");
}


// __________________________________________________________________________________________ //
AQcdnumDISCSEWFit::~AQcdnumDISCSEWFit() {
}


// ___________________________________________________________________________________________ //
double AQcdnumDISCSEWFit::Gets2wOS_from_s2wMSbar(double q2, double sw2MSbar) { 
   //const double q2v[] = {3.5, 5, 6.5, 8.5, 12, 15, 20, 25, 35, 45, 60, 90, 100, 120, 150, 200, 250, 300, 400, 500, 650, 800, 1000, 1200, 1500, 2000, 3000, 5000, 8000, 12000, 15000, 20000, 30000, 50000,  };
   static map<double,TGraph*> lumap;
   static double mz = PAR(EPRC.mZ);//PAR(Mz);
   //if ( PAR(Mz) != mz ) {
   if ( PAR(EPRC.mZ) != mz ) {
      warn["Gets2wOS_from_s2wMSbar"]<<"mZ has changed! Deleteing lumap."<<endl;
      lumap.clear();
      mz = PAR(Mz);
   }
   if ( !lumap.count(q2) ) {
      info["Gets2wOS_from_s2wMSbar"]<<"Init new lookup map to map sw2MSbar to sw2OS at q2="<<q2<<", and mz="<<mz<<endl;
      // get TGraph
      lumap[q2] = new TGraph();
      TGraph* g = lumap[q2];
      const double mwlo = 77;
      const double mwup = 85;
      const double step = 0.02;
      //const double mz = PAR(Mz);
      double mw0 = 0;
      eprc_get_mw_(&mw0);
      double mt0 = 0;
      eprc_get_mt_(&mt0);
      double mtmsbar = 164;
      eprc_set_mt_(&mtmsbar);
      //info["Gets2wOS_from_s2wMSbar"]<<"mw0 = "<<mw0<<endl;
      int ip=0;
      for ( double mw = mwlo; mw<=mwup+1.e-7 ; mw+=step ) {
	 eprc_set_mw_(&mw);
	 double sw2OS = 1. - mw*mw/mz/mz;
	 double sw2MS = eprc_sw2_cm_(&q2);
	 //cout<<"q2="<<q2<<"\tq="<<sqrt(q2)<<"\tmw="<<mw<<"\tsw2OS="<<sw2OS<<"\tsw2MS="<<sw2MS<<"\tsw2MS/sw2OS="<<sw2MS/sw2OS<<endl;
	 lumap[q2]->SetPoint(ip++,sw2MS,sw2OS);
      }
      eprc_set_mw_(&mw0);
      eprc_set_mt_(&mt0);
   }
   return lumap[q2]->Eval(sw2MSbar);
}


// ___________________________________________________________________________________________ //
bool AQcdnumDISCSEWFit::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   AlposObject::debug["Init"]<<endl;
   //CONST(e-charge);
   return true;
}


// __________________________________________________________________________________________ //
bool AQcdnumDISCSEWFit::Update() {  //alpos
   debug["Update"]<<"AlposName: "<<GetAlposName()<<endl;

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit); 

   // calculation of Gf for mW-Gf plot
   // {   
   //    for ( double mZ = 91.1 ; mZ<=91.34 ; mZ+=0.1 ) {
   // 	 //const double mZ = 91.1;
   // 	 const double mz2 = mZ*mZ;
   // 	 SET(EPRC.mZ,mZ,0);
   // 	 cout<<"\n >> mZ = "<<mZ<<endl;
   // 	 for ( double mW = 78; mW<=86.01 ; mW+=0.2 ) {
   // 	    double mw2=mW*mW;
   // 	    SET(EPRC.mW,mW,0);
   // 	    const vector<double>& ewpar = VALUES(EPRC);
   // 	    const double dr = ewpar[1];
   // 	    const double aem = PAR(EPRC.alpha);
   // 	    double s2thw0 = 1. - mw2/mz2;
   // 	    double Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
   // 	    printf("mW = %5.2f     Gf = %10.4e\n",mW,Gf);
   // 	 }
   //    }
   // 	 cout<<"ende gF calculation!"<<endl;
   //    exit(1);
   // }


   // --------------------------------------------
   // --- get q2, x, y grid points
   vector<double> q2 = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   vector<double> x = DOUBLE_COL_NS(Data,x,GetAlposName());
   vector<double> y = DOUBLE_COL_NS(Data,y,GetAlposName());
   fValue.resize(q2.size());
   fError.resize(q2.size());   
   
   // --------------------------------------------
   // --- calculate y, if y was not provided in data table
   if ( y.empty() ) {
      double sqs=DOUBLE_NS(sqrt-s,GetAlposName());
      double ss = sqs*sqs;
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) { 
	 y[i] = q2[i] / (x[i]*ss);
      }
   }

   // --------------------------------------------
   // --- get charge, polarity, and other flags from data table
   double charge = DOUBLE_NS(e-charge,GetAlposName()); //PAR(e-charge);//
   double polty  = PAR(e-polarity) ;//DOUBLE_NS(e-polarity,GetAlposName());
   bool IsRedCS  = BOOL_NS(IsReducedCS,GetAlposName());
   bool IsNC     = BOOL_NS(IsNC,GetAlposName());

   if ( charge==0 ) {
      error["Update"]<<"Could not get charge of lepton."<<endl;
      exit(1);
   }

   // --------------------------------------------
   // ------ calculate structure functions
   // --------------------------------------------
   int iFL=1, iF2=2, ixF3=3,nchk=0;
   vector<double> FL(q2.size()), F2(q2.size()), xF3(q2.size());
   vector<double> FLm(q2.size()), F2m(q2.size()), xF3m(q2.size());
   vector<double> CEP2F,CEM2F,CEP3F,CEM3F;

   vector<double> CEP2F_q,CEM2F_q,CEP3F_q,CEM3F_q;
   vector<double> CEP2F_b,CEM2F_b,CEP3F_b,CEM3F_b;
   vector<double> FL_q(q2.size()), F2_q(q2.size()), xF3_q(q2.size());
   vector<double> FLm_q(q2.size()), F2m_q(q2.size()), xF3m_q(q2.size());
   vector<double> FL_b(q2.size()), F2_b(q2.size()), xF3_b(q2.size());
   vector<double> FLm_b(q2.size()), F2m_b(q2.size()), xF3m_b(q2.size());
   //      vector<double> F2Gamma(q2.size()),FLGamma(q2.size());

   // --------------------------------------------
   // --- define linear combinations
   if ( IsNC ) {
      // --- general linear combinations
      CEP2F = {0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM2F = {0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0. };
      CEP3F = {0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM3F = {0.,-1.,0.,-1.,0.,-1.,0.,1.,0.,1.,0.,1.,0. };
      // --- linear combinations for EW-corrections
      // q
      CEP2F_q = {0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM2F_q = {0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,1.,0. };
      CEP3F_q = {0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM3F_q = {0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,1.,0. };
      // qbar
      CEP2F_b = {0.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0. };
      CEM2F_b = {0.,1.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0. };
      CEP3F_b = {0.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0. };
      CEM3F_b = {0.,-1.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.,0. };
   }
   else {
      // --- general linear combinations
      CEP2F = {0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM2F = {0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0. };
      CEP3F = {0.,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM3F = {0.,0. ,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0. };
      // --- linear combinations for EW-corrections
      // q
      CEP2F_q = {0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM2F_q = {0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEP3F_q = {0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM3F_q = {0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0. };
      // qbar
      CEP2F_b = {0.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0. };
      CEM2F_b = {0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0. };
      CEP3F_b = {0.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,0. };
      CEM3F_b = {0.,0.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.,0. };
   }

   // --------------------------------------------
   // ---- get structure functions
   int npts=int(q2.size());
   if ( IsNC ) {
      // --- structure functions
      zmstfun_(&iFL, &CEP2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
      zmstfun_(&iF2, &CEP2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEP3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
      zmstfun_(&iFL, &CEM2F[0],&x[0],&q2[0],&FLm[0],&npts,&nchk);
      zmstfun_(&iF2, &CEM2F[0],&x[0],&q2[0],&F2m[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEM3F[0],&x[0],&q2[0],&xF3m[0],&npts,&nchk);
      // --- structure functions for quark/anti-quarks separately
      /*
      // q
      zmstfun_(&iFL, &CEP2F_q[0],&x[0],&q2[0],&FL_q[0],&npts,&nchk);
      zmstfun_(&iF2, &CEP2F_q[0],&x[0],&q2[0],&F2_q[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEP3F_q[0],&x[0],&q2[0],&xF3_q[0],&npts,&nchk);
      zmstfun_(&iFL, &CEM2F_q[0],&x[0],&q2[0],&FLm_q[0],&npts,&nchk);
      zmstfun_(&iF2, &CEM2F_q[0],&x[0],&q2[0],&F2m_q[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEM3F_q[0],&x[0],&q2[0],&xF3m_q[0],&npts,&nchk);
      // qbar
      zmstfun_(&iFL, &CEP2F_b[0],&x[0],&q2[0],&FL_b[0],&npts,&nchk);
      zmstfun_(&iF2, &CEP2F_b[0],&x[0],&q2[0],&F2_b[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEP3F_b[0],&x[0],&q2[0],&xF3_b[0],&npts,&nchk);
      zmstfun_(&iFL, &CEM2F_b[0],&x[0],&q2[0],&FLm_b[0],&npts,&nchk);
      zmstfun_(&iF2, &CEM2F_b[0],&x[0],&q2[0],&F2m_b[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEM3F_b[0],&x[0],&q2[0],&xF3m_b[0],&npts,&nchk);
      */
   }
   else {
      if ( charge > 0 ) { // CC e+
	 // --- structure functions
	 zmstfun_(&iFL, &CEP2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEP2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEP3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
	 /*
	 // --- structure functions for quark/anti-quarks separately
	 // q
	 zmstfun_(&iFL, &CEP2F_q[0],&x[0],&q2[0],&FL_q[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEP2F_q[0],&x[0],&q2[0],&F2_q[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEP3F_q[0],&x[0],&q2[0],&xF3_q[0],&npts,&nchk);
	 // qbar
	 zmstfun_(&iFL, &CEP2F_b[0],&x[0],&q2[0],&FL_b[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEP2F_b[0],&x[0],&q2[0],&F2_b[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEP3F_b[0],&x[0],&q2[0],&xF3_b[0],&npts,&nchk);
	 */
      }
      else { // CC e-
	 // --- structure functions
	 zmstfun_(&iFL, &CEM2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEM2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEM3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
	 /*
	 // --- structure functions for quark/anti-quarks separately
	 // q
	 zmstfun_(&iFL, &CEM2F_q[0],&x[0],&q2[0],&FL_q[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEM2F_q[0],&x[0],&q2[0],&F2_q[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEM3F_q[0],&x[0],&q2[0],&xF3_q[0],&npts,&nchk);
	 // qbar
	 zmstfun_(&iFL, &CEM2F_b[0],&x[0],&q2[0],&FL_b[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEM2F_b[0],&x[0],&q2[0],&F2_b[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEM3F_b[0],&x[0],&q2[0],&xF3_b[0],&npts,&nchk);
	 */
      }
   }

   // cout<<"\n ----- " <<GetAlposName()<<" ------------"<<endl;
   // for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
   //    cout<<i<<", q2="<<q2[i]<<", x="<<x[i]<<", y="<<y[i]<<", f2="<<F2[i]<<endl;
   // }   
   // // exit(1);

   
   // --------------------------------------------
   // --- parameters for cross section calculation
   double s2thweff = PAR(sin2thw);
   double slopeFitQ20 = s2thweff; 
   double slopeFitS0 = PAR(I3Ru); // should be initialized to sw2
   double slopeFitS1 = PAR(I3Rd); 
   //const double cos2thw = 1-s2thweff;
   double Mz = PAR(Mz);
   double mz2 = Mz*Mz;
   const double aem = PAR(EPRC.alpha);
   double Mw = PAR(Mw);
   double mw2 = Mw*Mw;
   double Gf = PAR(Gf);
   if ( Gf != 0 && Mz != 0 ) { // MOMS
      for ( int i=0; i<3 ; i++ ) { // 3 iteration because of mW dependence of dR
	 const vector<double>& ewpar = VALUES(EPRC);
	 double dr = ewpar[1];
	 double A = aem * M_PI/ (sqrt(2)*Gf);
	 mw2 = mz2/2.* (1 + sqrt(1-4*A/(mz2*(1-dr))));
	 Mw = sqrt(mw2);
	 SET(EPRC.mW,Mw,0);
	 SET(Mw,Mw,0);
      }
      //cout<<"Gf="<<Gf<<"\tMw="<<Mw<<endl;
   }
   else if ( Gf != 0 && Mz == 0 ) {
      // OMS(!) with (aem,mz,mw,dr) but fit gf(mz,...)
      for ( int i=0; i<3 ; i++ ) { // 3 iteration because of mW dependence of dR   
         const vector<double>& ewpar = VALUES(EPRC);
         double dr = ewpar[1];
         double A = aem * M_PI/ (sqrt(2)*Gf);
	 mz2 = mw2/(1 - A/mw2/(1-dr));
	 Mz = sqrt(mz2);
         SET(EPRC.mZ,Mz,0);
         //SET(Mz,Mz,0); //DO NOT SET IT!
      }

   }
   
   double s2thw0 = 1. - mw2/mz2;
   static const double convfac= 0.389379338e9; //0.389379323e9;


   // --------------------------------------------
   // --- cross section calculation
   // --------------------------------------------
   map<double,double> mq2r, ms2weff, kq;
   if ( IsNC ) {
      static const double euq =  2./3.;
      static const double edq = -1./3.;
      static const double e2u = euq*euq;
      static const double e2d = edq*edq;
      double au, ad, vu, vd, ve;
      
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) { // loop over all data points
	 // --- access values from EPRC or used already calculated ones
         if ( mq2r.count(q2[i])==0 ) { 
            SET(EPRC.mur,sqrt(q2[i]),0);
            const vector<double>& ewpar = VALUES(EPRC);
            mq2r[q2[i]]    = ewpar[1];
            ms2weff[q2[i]] = ewpar[2];
            kq[q2[i]] = ewpar[2]/s2thw0;
         }

	 /*
	 static set<double> allq2;
	 int s0 = allq2.size();
	 allq2.insert(q2[i]);
	 if ( s0 != allq2.size() ) {
	    for ( auto qq : allq2 ) cout<<"  "<<qq;
	    cout<<endl;
	    for ( auto qq : allq2 )  cout<<"  "<<sqrt(qq);
	    cout<<endl;
	 }
	 */

	 // --------------------------------------------
	 // --- couplings 
	 // --- old code with incorrect EW corrections
	 // --- and some functionality depending of 
	 // --- fit methodlogy
	 // --------------------------------------------

	 // --- ve 
         //ve = -0.5 + 2.*s2thweff; // HERE s2thweff must NOT be taken from EPRC: kq (right)
         //ve = -0.5 + 2.*ms2weff[q2[i]]; // HERE s2thweff must NOT be taken from EPRC: kq (wrong)
	 ve = -0.5 + 2.*s2thw0;

	 // --- quark couplings 
	 // --- various versions depending on the fit-type
	 if ( slopeFitQ20>=1 ) {
	    double sw2MSbar = slopeFitS0 + slopeFitS1*(log(q2[i])-log(slopeFitQ20));
	    s2thw0 = Gets2wOS_from_s2wMSbar(q2[i],sw2MSbar);// get OS sw2 from lu-table
	    ve = -0.5 + 2.*s2thw0;
	    au =  0.5;
	    ad = -0.5;
	    vu =  0.5 - 2*euq*s2thw0; // all on-shell for the moment 1
	    vd = -0.5 - 2*edq*s2thw0; // all on-shell for the moment 1
	    mw2 = (1. - s2thw0) * mz2 ;
	    Mw = sqrt(mw2);
	 }
	 else if ( slopeFitQ20 < 0 ) { // fit three points
	    if ( q2[i]<=1100 ) Mw = PAR(I3u);
	    else  if ( q2[i]>1100  && q2[i]<=4000 ) Mw = PAR(I3d);
	    else if ( q2[i]>4000 ) Mw = PAR(I3Ru);
	    //else if ( q2[i]>13000 ) Mw = PAR(I3Rd);


	    // // 7pts
	    // // if ( q2[i] <200  ) Mw = PAR(EPRC.au);
	    // // else if  ( q2[i]>=200 && q2[i]<=800 ) Mw = PAR(EPRC.ad);
	    // if ( q2[i] <=800  ) Mw = PAR(EPRC.au);
	    // else if  ( q2[i]==1000 || q2[i]==1200 || q2[i]==1500 ) Mw = PAR(EPRC.vu);
	    // else if  ( q2[i]==2000 ) Mw = PAR(EPRC.vd);
	    // else if ( q2[i]==3000 ) Mw = PAR(I3u);
	    // else if ( q2[i]==5000 ) Mw = PAR(I3d);
	    // else if ( q2[i]==8000 ) Mw = PAR(I3Ru);
	    // else if ( q2[i]>8000 ) Mw = PAR(I3Rd);
	    

	    SET(EPRC.mW,Mw,0);
            const vector<double>& ewpar = VALUES(EPRC);
            mq2r[q2[i]]    = ewpar[1];

	    mw2 = Mw*Mw;
	    s2thw0 = 1.-mw2/mz2;
	    ve = -0.5 + 2.*s2thw0;
	    au =  0.5;
	    ad = -0.5;
	    vu =  0.5 - 2*euq*s2thw0; // all on-shell for the moment 1
	    vd = -0.5 - 2*edq*s2thw0; // all on-shell for the moment 1
	 }
	 else if ( PAR(I3Ru) == 0 ) { 
	    au = PAR(au) == 0 ?   0.5 /*PAR(I3u)*/ : PAR(au); // 0.5; I3u
	    ad = PAR(ad) == 0 ?  -0.5 /*PAR(I3d)*/ : PAR(ad); //-0.5; I3d

	    // --- kq expressed
	    // vu = PAR(vu0)!= 0 ? PAR(vu0) - (1-kq[q2[i]])*2*euq*s2thw0 : PAR(I3u) - 2*euq*s2thweff;//PAR(I3u) - 2*euq*ms2weff[q2[i]]
	    // vd = PAR(vd0)!= 0 ? PAR(vd0) - (1-kq[q2[i]])*2*edq*s2thw0 : PAR(I3d) - 2*edq*s2thweff;//PAR(I3d) - 2*edq*ms2weff[q2[i]]
	    vu = PAR(vu0) == 0 ?  0.5 - 2*euq*s2thw0 : PAR(vu0);//PAR(I3u) - 2*euq*ms2weff[q2[i]]
	    vd = PAR(vd0) == 0 ? -0.5 - 2*edq*s2thw0 : PAR(vd0);//PAR(I3d) - 2*edq*ms2weff[q2[i]]

	 }
	 else {
	    // --- attempt to fit (assumed) right-handed couplings
	    if ( PAR(I3Rd)==0 ) {
	       // fit I3R ad I3L
	       au =     PAR(I3u) - PAR(I3Ru);
	       vu =     PAR(I3u) + PAR(I3Ru) - 2*euq*s2thweff;
	       ad = -1.*PAR(I3u) + PAR(I3Ru);
	       vd = -1.*PAR(I3u) - PAR(I3Ru) - 2*edq*s2thweff;
	    }
	    else {
	       // fit I3Ru ad I3Rd
	       au =   0.5 - PAR(I3Ru);
	       ad =  -0.5 - PAR(I3Rd);
	       vu =   0.5 + PAR(I3Ru) - 2*euq*s2thweff;
	       vd =  -0.5 + PAR(I3Rd) - 2*edq*s2thweff;
	    }
	 }
	 // --------------------------------------------

	 // au = 0.5;
	 // ad =-0.5;
	 // vu = 0.2028037;
	 // vd = -0.35140185;

	 // au = 0.5;
	 // ad =-0.5;

	 // vu = 0.5 - eps*2*euq*s2thw0;
	 // vd = -0.5 - eps*2*edq*s2thw0;
	  
	 // --------------------------------------------
	 // --- NC Propagator in OMS or MOMS scheme
	 // --- if Gf is defined or not...
	 double kzs;
         if ( Gf == 0  || PAR(Mz)==0) {
            // OMS (aem,Mz,Mw,dR)
            double dr = mq2r[q2[i]];
            kzs = mz2/(4*mw2*s2thw0)/(1.-dr);
         }
         else {
            // MOMS (aem,Gf,Mz,Mw)
            kzs = Gf*mz2/2/sqrt(2)/M_PI/aem;
         }
	 double kz= q2[i]/(q2[i]+mz2)*kzs;

	 //cout<<"Q2: "<<q2[i]<<"\tkz="<<kz<<"\tvu="<<vu<<"\tvd="<<vd<<"\tvu0="<<0.5<<"\ts2e="<<ms2weff[q2[i]]<<"\teps="<<eps<<"\tdr="<<mq2r[q2[i]]<<endl;


	 // --------------------------------------------
	 // --- structure functions
         const double ae = ( charge>0 ) ? 0.5 : -0.5; // ae including sign for e+/e- str.fun.
	 double PagZ=  (-ve + polty*ae)*kz;
	 double PaZ2 = (ve*ve + ae*ae - 2*polty*ve*ae)*kz*kz;
	 double A_u = e2u + 2.*euq*vu*PagZ + (vu*vu+au*au)*PaZ2;
	 double A_d = e2d + 2.*edq*vd*PagZ + (vd*vd+ad*ad)*PaZ2;
	 F2[i]   = A_u*F2[i]   + A_d*F2m[i];
	 FL[i]   = A_u*FL[i]   + A_d*FLm[i];	    

	 double PbgZ = 2.*(-ae + polty*ve)*kz;
	 double PbZ2  = (2.*ve*ae - polty*(ve*ve+ae*ae) )*kz*kz*2.;
	 double B_u = euq*au*PbgZ + vu*au*PbZ2;
	 double B_d = edq*ad*PbgZ + vd*ad*PbZ2;
	 xF3[i]  = B_u*xF3[i]  + B_d*xF3m[i];

	 /*
	 // --------------------------------------------
	 // --- structure functions with full EW
	 // --- corrections (full flavor breakdown)
	 // --------------------------------------------

	 // --------------------------------------------
	 // --- form factors
	 // --- todo: get values from EPRC
	 double r_euq = 1;
	 double r_eub = 1;
	 double r_edq = 1;
	 double r_edb = 1;
	 double k_e = 1;
	 double k_uq = 1;
	 double k_ub = 1;
	 double k_dq = 1;
	 double k_db = 1;
	 double k_euq = 1;
	 double k_eub = 1;
	 double k_edq = 1;
	 double k_edb = 1;


	 // --------------------------------------------
	 // couplings with form factors
	 // following EPRC manual, eq. 14, 15 and 16
	 ve = -0.5 + 2.*s2thw0*k_e; // fabs(eq) ??
	 double vuq =  0.5 - 2*fabs(euq)*s2thw0*k_uq; 
	 double vub =  0.5 - 2*fabs(euq)*s2thw0*k_ub; 
	 double vdq = -0.5 + 2*fabs(edq)*s2thw0*k_dq; 
	 double vdb = -0.5 + 2*fabs(edq)*s2thw0*k_db; 

	 double I3e = -0.5; // I3e  ?!
	 double I3u =  0.5;
	 double I3d = -0.5;
	 double vveuq = I3e*vuq + ve*I3u - I3e*I3u*(1-16*fabs(euq*charge)*s2thw0*s2thw0*k_euq);
	 double vveub = I3e*vub + ve*I3u - I3e*I3u*(1-16*fabs(euq*charge)*s2thw0*s2thw0*k_eub);
	 double vvedq = I3e*vdq + ve*I3d - I3e*I3d*(1-16*fabs(edq*charge)*s2thw0*s2thw0*k_edq);
	 double vvedb = I3e*vdb + ve*I3d - I3e*I3d*(1-16*fabs(edq*charge)*s2thw0*s2thw0*k_edb);

	 // --- For tests to have identical values with 
	 // --- initial code. Differences:
	 // ---   + vveq = ve*vq
	 // ---   + vq was using sin2thweff (!)
	 // vuq =  0.5 - 2*fabs(euq)*s2thweff*k_uq; // this uses sin2thw-eff!
	 // vub =  0.5 - 2*fabs(euq)*s2thweff*k_ub; // this uses sin2thw-eff!
	 // vdq = -0.5 + 2*fabs(edq)*s2thweff*k_dq; // this uses sin2thw-eff!
	 // vdb = -0.5 + 2*fabs(edq)*s2thweff*k_db; // this uses sin2thw-eff!
	 // vveuq = ve*vuq;
	 // vveub = ve*vub;
	 // vvedq = ve*vdq;
	 // vvedb = ve*vdb;

	 
	 // --------------------------------------------
	 // --- couplings of quarks and leptons to the currents
	 double A_uq = e2u + 2.*euq*(-vveuq + polty*ae*vuq)*kz*r_euq + (vuq*vuq+au*au)*(ve*ve + ae*ae - 2*polty*ve*ae)*kz*kz*r_euq*r_euq;
	 double A_ub = e2u + 2.*euq*(-vveub + polty*ae*vub)*kz*r_eub + (vub*vub+au*au)*(ve*ve + ae*ae - 2*polty*ve*ae)*kz*kz*r_eub*r_eub;
	 double A_dq = e2d + 2.*edq*(-vvedq + polty*ae*vdq)*kz*r_edq + (vdq*vdq+ad*ad)*(ve*ve + ae*ae - 2*polty*ve*ae)*kz*kz*r_edq*r_edq;
	 double A_db = e2d + 2.*edq*(-vvedb + polty*ae*vdb)*kz*r_edb + (vdb*vdb+ad*ad)*(ve*ve + ae*ae - 2*polty*ve*ae)*kz*kz*r_edb*r_edb;

	 double B_uq =       2.*euq*au*(-ae + polty*ve)*kz*r_euq + (2.*ve*ae*vuq*au - polty*(ve*ve+ae*ae)*vuq*au )*kz*kz*r_euq*r_euq*2.;
	 double B_ub =       2.*euq*au*(-ae + polty*ve)*kz*r_eub + (2.*ve*ae*vub*au - polty*(ve*ve+ae*ae)*vub*au )*kz*kz*r_eub*r_eub*2.;
	 double B_dq =       2.*edq*ad*(-ae + polty*ve)*kz*r_edq + (2.*ve*ae*vdq*ad - polty*(ve*ve+ae*ae)*vdq*ad )*kz*kz*r_edq*r_edq*2.;
	 double B_db =       2.*edq*ad*(-ae + polty*ve)*kz*r_edb + (2.*ve*ae*vdb*ad - polty*(ve*ve+ae*ae)*vdb*ad )*kz*kz*r_edb*r_edb*2.;

	 // --- structure functions with form factors
	 // --- todo. These are not passed further
	 double F2EW = A_uq*F2_q[i]  + A_ub*F2_b[i]  + A_dq*F2m_q[i]  + A_db*F2m_b[i];
	 double FLEW = A_uq*FL_q[i]  + A_ub*FL_b[i]  + A_dq*FLm_q[i]  + A_db*FLm_b[i];
	 double F3EW = B_uq*xF3_q[i] + B_ub*xF3_b[i] + B_dq*xF3m_q[i] + B_db*xF3m_b[i];

	 if ( fabs(F2EW/F2[i]-1) > 1.e-6)  cout<<"Difference in old/new code: F2: F2EW="<<F2EW<<"\tF2="<<F2[i]<<endl;
	 if ( fabs(FLEW/FL[i]-1) > 1.e-6)  cout<<"Difference in old/new code: FL: FLEW="<<FLEW<<"\tF2="<<FL[i]<<endl;
	 if ( fabs(F3EW/xF3[i]-1) > 1.e-6) cout<<"Difference in old/new code: F3: FEEW="<<F3EW<<"\tF2="<<xF3[i]<<endl;
	 */

	 // F2Gamma[i] = 4./9. * F2[i] + F2m[i]/9.;
	 // FLGamma[i] = 4./9. * FL[i] + FLm[i]/9.;
      }
   }
   else { //CC
      //for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 // F2Gamma[i] = 4./9*F2[i] + F2m[i]/9.;
	 // FLGamma[i] = 4./9*FL[i] + FLm[i]/9.;
      //}
   }
  

   // --------------------------------------------
   // ------ calc CS
   if ( IsNC ) { //NC
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 double yplus  = 1+(1-y[i])*(1-y[i]);
	 double yminus = 1-(1-y[i])*(1-y[i]);
	 fValue[i] = F2[i] + yminus/yplus*xF3[i] - y[i]*y[i]/yplus*FL[i];
	 if ( !IsRedCS ) 
	    fValue[i] *= 2*M_PI/(x[i]*q2[i]*q2[i])*convfac*aem*aem;
      }
   }
   else { //CC
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 double yplus  = 1+(1-y[i])*(1-y[i]);
	 double yminus = 1-(1-y[i])*(1-y[i]);
	 if ( charge>0 )  // CCe+ 
	    fValue[i] = 0.5*(yplus*F2[i] - yminus*xF3[i] - y[i]*y[i]*FL[i]);
	 else
	    fValue[i] = 0.5*(yplus*F2[i] + yminus*xF3[i] - y[i]*y[i]*FL[i]);
	 fValue[i] *=  ( 1 + polty*charge);

	 // --- EW HO form factors from EPRC
	 //     following EPRC manual eq. 27
	 //     todo: These are not passed further
	 double rW_euq2 = 1;
	 double rW_eub2 = 1;
	 double rW_edq2 = 1;
	 double rW_edb2 = 1;
	 double sEW = 0; // sigma with EW corrections from EPRC
	 if ( charge>0 )  // CCe+
	    sEW = 0.5*(yplus*(F2_q[i]*rW_edq2+F2_b[i]*rW_eub2) - yminus*(xF3_q[i]*rW_edq2+xF3_b[i]*rW_eub2) - y[i]*y[i]*(FL_q[i]*rW_edq2+FL_b[i]*rW_eub2));
	 else  //CC e-
	    sEW = 0.5*(yplus*(F2_q[i]*rW_euq2+F2_b[i]*rW_edb2) + yminus*(xF3_q[i]*rW_edq2+xF3_b[i]*rW_eub2) - y[i]*y[i]*(FL_q[i]*rW_euq2+FL_b[i]*rW_edb2));
	 sEW *=  ( 1 + polty*charge);
      }
      
      if ( !IsRedCS ) {
	    if ( slopeFitQ20>=1 ) {
	       //SET(EPRC.mur,Mz,0); // any scale
	       const vector<double>& ewpar = VALUES(EPRC);
	       double dr  = ewpar[1];
	       for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
		  double sw2MSbar = slopeFitS0 + slopeFitS1*(log(q2[i])-log(slopeFitQ20));
		  s2thw0 = Gets2wOS_from_s2wMSbar(q2[i], sw2MSbar);
		  mw2 = (1. - s2thw0) * mz2;
		  double mwprop = sqrt( mw2 );
		  Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
		  fValue[i] *= pow(mwprop,4)/pow(mwprop*mwprop+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	       }
	    }
	    if ( slopeFitQ20<0 ) { // fit three points
	       for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	    if ( q2[i]<=1100 ) Mw = PAR(I3u);
	    else  if ( q2[i]>1100  && q2[i]<=4000 ) Mw = PAR(I3d);
	    else if ( q2[i]>4000 ) Mw = PAR(I3Ru);
	    //else if ( q2[i]>13000 ) Mw = PAR(I3Rd);

	    // 7pts
	    // if ( q2[i] <=800  ) Mw = PAR(EPRC.au);
	    // else if  ( q2[i]==1000 || q2[i]==1200 || q2[i]==1500 ) Mw = PAR(EPRC.vu);
	    // else if  ( q2[i]==2000 ) Mw = PAR(EPRC.vd);
	    // else if ( q2[i]==3000 ) Mw = PAR(I3u);
	    // else if ( q2[i]==5000 ) Mw = PAR(I3d);
	    // else if ( q2[i]==8000 ) Mw = PAR(I3Ru);
	    // else if ( q2[i]>8000 ) Mw = PAR(I3Rd);

		  mw2 = Mw*Mw;
		  s2thw0 = 1.-mw2/mz2;
		  double mwprop = Mw;

		  SET(EPRC.mW,Mw,0);
		  const vector<double>& ewpar = VALUES(EPRC);
		  double dr  = ewpar[1];

		  Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
		  fValue[i] *= pow(mwprop,4)/pow(mwprop*mwprop+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	       }
	    }
	    else {
	       double mwprop = PAR(MwProp);
	       if ( Gf==0 || PAR(Mz)==0 ) {
		  // OMS (aem,Mz,Mw,dR)
		  SET(EPRC.mur,Mz,0); // any scale
		  const vector<double>& ewpar = VALUES(EPRC);
		  double dr  = ewpar[1];
		  Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
		  //Gf = M_PI*aem/sqrt(2)/s2thw0/(mwprop*mwprop)*(1/(1-dr));
	       }
	       for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
		  fValue[i] *= pow(mwprop,4)/pow(mwprop*mwprop+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	       }	       
	       
	 //       double mwprop = PAR(MwProp);
	 // if ( Gf==0 ) {
	 //    // OMS (aem,Mz,Mw,dR)
	 //    SET(EPRC.mur,Mz,0); // any scale
	 //    const vector<double>& ewpar = VALUES(EPRC);
	 //    double dr  = ewpar[1];
	 //    Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
	 // }
	 // for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 //    if ( slopeFitQ20>=10 ) {
	 //       SET(EPRC.mur,Mz,0); // any scale
	 //       const vector<double>& ewpar = VALUES(EPRC);
	 //       double dr  = ewpar[1];
	 //       double sw2MSbar = slopeFitS0 + slopeFitS1*(log(q2[i])-log(slopeFitQ20));
	 //       s2thw0 = Gets2wOS_from_s2wMSbar(q2[i], sw2MSbar);
	 //       mwprop = sqrt( (1. - s2thw0) * mz2 );
	 //       mw2 = mwprop*mwprop;
	 //       Gf = M_PI*aem/sqrt(2)/s2thw0/mw2*(1/(1-dr));
	 //       fValue[i] *= pow(mwprop,4)/pow(mwprop*mwprop+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	 //    }
	 //    else {
	 //       fValue[i] *= pow(mwprop,4)/pow(mwprop*mwprop+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
	 //    }
	 }
      }
   }


   // --------------------------------------------
   // ------ calc non-reduced CS if requested

   // ------ done
   fError.resize(fValue.size());

   
    // cout<<"\n ----- " <<GetAlposName()<<" ------------"<<endl;
    // cout<<"ch: "<<charge<<", pol="<<polty<<", red: "<<IsRedCS<<", IsNC: "<<IsNC<<endl;
    // vector<double> data = DOUBLE_COL_NS(Data,Sigma,GetAlposName());
    // for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
    //    //printf(" idx:%12d     %18.16f\n",i,fValue[i]);
    //    cout<<i<<", q2="<<q2[i]<<", x="<<x[i]<<", y="<<y[i]<<", f2="<<F2[i]<<", cs="<<fValue[i]<<" ("<<data[i]<<" [da])"<<endl;   
    // }
    // exit(1);


   return true;
}

//______________________________________________________________________________

