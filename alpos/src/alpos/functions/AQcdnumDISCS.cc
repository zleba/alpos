// DB 16.01.2015

#include "alpos/functions/AQcdnumDISCS.h"
#include "alpos/functions/AQcdnumInit.h"
#include <iostream>
#include "fastnlotk/read_steer.h"

using namespace std;

extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);

// __________________________________________________________________________________________ //
const std::vector<std::string> AQcdnumDISCS::fRequirements = {"QcdnumInit",
							      "Mz","Mw", // boson masses
							      "Gf", // gf
							      "sin2thw", // sin2th_w effective
							      "sin2thweffFix", // sin2th_w_effective fixed for EW fit. Set to 0 if ignore it
							      "EPRC",
							      "e-polarity", // polarity
							      //"e-charge", // charge
							      "au","ad","vu0","vd0", // boson-quark couplings
}; //!< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumDISCS::fStopFurtherNotification = {}; //!< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumDISCS::fFunctionName = "QcdnumDISCS"; //!< The function's name

// __________________________________________________________________________________________ //
AQcdnumDISCS::AQcdnumDISCS(const std::string& name) : AParmFuncBase<double>(name) { 
   AlposObject::SetClassName("AQcdnumDISCS");
}


// __________________________________________________________________________________________ //
AQcdnumDISCS::~AQcdnumDISCS() {
}


// ___________________________________________________________________________________________ //
bool AQcdnumDISCS::Init() { //alpos
   //! Init is once called for each function
   //! return true if initialization was successful.
   AlposObject::debug["Init"]<<endl;
   //CONST(e-charge);
   return true;
}


// __________________________________________________________________________________________ //
bool AQcdnumDISCS::Update() {  //alpos
   debug["Update"]<<"AlposName: "<<GetAlposName()<<endl;

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit); 

   // get 
   vector<double> q2 = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   vector<double> x = DOUBLE_COL_NS(Data,x,GetAlposName());
   vector<double> y = DOUBLE_COL_NS(Data,y,GetAlposName());
   fValue.resize(q2.size());
   fError.resize(q2.size());

   if ( y.empty() ) {
      double sqs=DOUBLE_NS(sqrt-s,GetAlposName());
      double ss = sqs*sqs;
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) { 
	 y[i] = q2[i] / (x[i]*ss);
      }
   }

   double charge = DOUBLE_NS(e-charge,GetAlposName()); //PAR(e-charge);//
   double polty  = PAR(e-polarity) ;//DOUBLE_NS(e-polarity,GetAlposName());
   bool IsRedCS  = BOOL_NS(IsReducedCS,GetAlposName());
   bool IsNC     = BOOL_NS(IsNC,GetAlposName());

   if ( charge==0 ) {
      error["Update"]<<"Could not get charge of lepton."<<endl;
      exit(1);
   }

   // ------ calc structure functions
   int iFL=1, iF2=2, ixF3=3,nchk=0;
   vector<double> FL(q2.size()), F2(q2.size()), xF3(q2.size());
   vector<double> FLm(q2.size()), F2m(q2.size()), xF3m(q2.size());
   vector<double> CEP2F,CEM2F,CEP3F,CEM3F;
   if ( IsNC ) {
      CEP2F = {0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM2F = {0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0. };
      CEP3F = {0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,1.,0.,0. };
      CEM3F = {0.,-1.,0.,-1.,0.,-1.,0.,1.,0.,1.,0.,1.,0. };
   }
   else {
      CEP2F = {0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM2F = {0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0. };
      CEP3F = {0.,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.,0. };
      CEM3F = {0.,0. ,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0. };
   }
   //      vector<double> F2Gamma(q2.size()),FLGamma(q2.size());
   int npts=int(q2.size());
   if ( IsNC ) {
      zmstfun_(&iFL, &CEP2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
      zmstfun_(&iF2, &CEP2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEP3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
      zmstfun_(&iFL, &CEM2F[0],&x[0],&q2[0],&FLm[0],&npts,&nchk);
      zmstfun_(&iF2, &CEM2F[0],&x[0],&q2[0],&F2m[0],&npts,&nchk);
      zmstfun_(&ixF3,&CEM3F[0],&x[0],&q2[0],&xF3m[0],&npts,&nchk);
   }
   else {
      if ( charge > 0 ) { // CC e+
	 zmstfun_(&iFL, &CEP2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEP2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEP3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
      }
      else { // CC e+
	 zmstfun_(&iFL, &CEM2F[0],&x[0],&q2[0],&FL[0],&npts,&nchk);
	 zmstfun_(&iF2, &CEM2F[0],&x[0],&q2[0],&F2[0],&npts,&nchk);
	 zmstfun_(&ixF3,&CEM3F[0],&x[0],&q2[0],&xF3[0],&npts,&nchk);
      }
   }

   // cout<<"\n ----- " <<GetAlposName()<<" ------------"<<endl;
   // for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
   //    cout<<i<<", q2="<<q2[i]<<", x="<<x[i]<<", y="<<y[i]<<", f2="<<F2[i]<<endl;
   // }   
   // // exit(1);
      
   // set<double> allq2;
   // for ( auto qq : q2 ) allq2.insert(qq);
   // cout<<"\n ----- " <<GetAlposName()<<" ------------"<<endl;
   // for ( auto qq : allq2 ) cout<<"\t"<<qq;
   // cout<<endl;

   double s2thweff = PAR(sin2thw);
   const double cos2thw = 1-s2thweff;
   const double s2thweffslope = PAR(sin2thweffFix); // not used (only to indicate 'EWFit' option)
   const double Mz = PAR(Mz);
   const double Mw = PAR(Mw);
   const double s2thw0 = 1. - Mw*Mw/(Mz*Mz);
   const double Gf = PAR(Gf);
   static const double convfac= 0.389379338e9; //0.389379323e9;

   if ( IsNC ) {
      static const double euq =  2./3.;
      static const double edq = -1./3.;
      static const double e2u = euq*euq;
      static const double e2d = edq*edq;
      double dR = 0;
      double au, ad, vu, vd, ve;
      double PZ=0;
      map<double,double> mq2r, ms2weff;
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 if ( s2thweffslope == 0 ) {
	    //! Use SM parameters of couplings
	    ve = -0.5 + 2.*s2thweff;
	    au =  0.5; // I_3 //PAR(au);// 0.5;
	    ad = -0.5; // I_3 //PAR(ad);//-0.5;
	    vu =  au - 2*euq*s2thweff;
	    vd =  ad - 2*edq*s2thweff;
	    // vu = PAR(vu0) - 2*euq*s2thweff;
	    // vd = PAR(vd0) - 2*edq*s2thweff;
	    PZ = 4. * s2thweff * cos2thw * (1.+Mz*Mz/q2[i]) * (1.-dR); // propagator factor for gamma/Z
	    PZ=1./PZ;
	 }
	 else {
	       if ( mq2r.count(q2[i])==0 ) {
		  SET(EPRC.mur,sqrt(q2[i]),0);
		  const vector<double>& ewpar = VALUES(EPRC);
		  mq2r[q2[i]]    = ewpar[1];
		  ms2weff[q2[i]] = ewpar[2];
		  //cout<<"q2="<<q2[i]<<"\tdR="<<dR<<endl;
	       }
	       //dR = VALUES(EPRC)[1];
	       dR = mq2r[q2[i]];
	       
	       if ( s2thweffslope==1 ) {
		  // take sin2thw from EPRC
		  s2thweff = ms2weff[q2[i]];
	       }
	       else {
		  // fit sin2thw
		  const double q2ref = 5000;//25;
		  //s2thweff = PAR(sin2thw) + s2thweffslope*log(q2[i]);
		  s2thweff = PAR(sin2thw) + s2thweffslope*(log(q2[i])-log(q2ref));
	       }

	       //! EW-fit style: 
	       au = PAR(au);// 0.5;
	       ad = PAR(ad);//-0.5;
	       ve = -0.5 + 2*s2thweff;
	       double eps = s2thweff/s2thw0 - 1.;
	       vu = PAR(vu0) - eps*2*euq*s2thw0;
	       vd = PAR(vd0) - eps*2*edq*s2thw0;


	       PZ = 4. * s2thw0 * (1.-s2thw0) * (1.+Mz*Mz/q2[i]) * (1.-dR); // propagator factor for gamma/Z
	       PZ=1./PZ;
	       //cout<<"Q2: "<<q2[i]<<"\tPZ="<<PZ<<"\tvu="<<vu<<"\tvd="<<vd<<"\tvu0="<<PAR(vu0)<<"\ts2e="<<s2thweff<<"\teps="<<eps<<"\tdr="<<dR<<endl;
	 }
	 

	 const double ae = ( charge>0 ) ? 0.5 : -0.5;
	 
	 double PagZ= 2.*(-ve + polty*ae)*PZ;
	 double PaZ2 = (ve*ve + ae*ae - 2*polty*ve*ae)*PZ*PZ;
	 double A_u = e2u + euq*vu*PagZ + (vu*vu+au*au)*PaZ2;
	 double A_d = e2d + edq*vd*PagZ + (vd*vd+ad*ad)*PaZ2;

	 // -- no gammZ and Z exchange as in case of Apfel++
	 //cout<<"QCDNUM  Q2="<<q2[i]<<"\tf2p="<<F2[i]<<"\tf2m="<<F2m[i]<<"\tflp="<<FL[i]<<"\tflm="<<FLm[i]<<endl;
	 // double PagZ= 2.*(-ve + polty*ae)*PZ;
	 // double PaZ2 = (ve*ve + ae*ae - 2*polty*ve*ae)*PZ*PZ;
	 // double A_u = e2u;
	 // double A_d = e2d;

	 F2[i]   = A_u*F2[i]   + A_d*F2m[i];
	 FL[i]   = A_u*FL[i]   + A_d*FLm[i];	    

	 double PbgZ = 2.*(-ae + polty*ve)*PZ;
	 double PbZ2  = (2.*ve*ae - polty*(ve*ve+ae*ae) )*PZ*PZ*2.;
	 double B_u = euq*au*PbgZ + vu*au*PbZ2;
	 double B_d = edq*ad*PbgZ + vd*ad*PbZ2;
	 xF3[i]  = B_u*xF3[i]  + B_d*xF3m[i];

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

   // ------ calc reduced CS
   for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
      double yplus  = 1+(1-y[i])*(1-y[i]);
      double yminus = 1-(1-y[i])*(1-y[i]);
      cout<<"QCDNUM  Q2="<<q2[i]<<"\tf2="<<F2[i]<<"\tfl="<<FL[i]<<"\tf3="<<xF3[i]<<"\tq="<<charge<<"\tpol="<<polty<<endl;

      if ( IsNC ) { //NC
	 fValue[i] = F2[i] + yminus/yplus*xF3[i] - y[i]*y[i]/yplus*FL[i];
      }
      else { //CC
	 if ( charge>0 ) { // CCe+ 
	    fValue[i] = 0.5*(yplus*F2[i] - yminus*xF3[i] - y[i]*y[i]*FL[i]);
	    fValue[i] *= (1+polty);
	 }
	 else { //CC e-
	    fValue[i] = 0.5*(yplus*F2[i] + yminus*xF3[i] - y[i]*y[i]*FL[i]);
	    fValue[i] *=(1-polty);
	 }
      }
   }

   // ------ calc non-reduced CS if needed
   if ( !IsRedCS ) {
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
 	 if ( IsNC ) {
	    double aem = 7.29735e-3; // 1/137.035999074 // 7.29927d-3;//
 	    fValue[i] *= 2*M_PI/(x[i]*q2[i]*q2[i])*convfac*aem*aem;
 	 }
 	 else { //CC
 	    fValue[i] *= pow(Mw,4)/pow(Mw*Mw+q2[i],2)*Gf*Gf/(2*M_PI*x[i])*convfac;
 	 }
      }
   }
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

