
#include "alpos/functions/AQcdnumDDISCS.h"
#include "alpos/functions/AQcdnumInit.h"


extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);
extern "C" void hqstfun_(int* istf, int *iflavour, double* def, double* x, double* Q2, double* f, int* n, int* nchk);
extern "C" void allfxq_( int* ityp, double* x, double* q, double* pdf, int *n, int* ichk); 

extern "C" {
    void qcd_2006_(double *z,double *q2, int *ifit, double *xPq, double *f2, double *fl, double *c2, double *cl);
    void h12006flux_(double *xpom, double *t, int *Int, int *ifit, int *ipom, double *flux);
}


  //hqstfun_(iF2, &icharm,&CEP2F[0],&beta[0],&q2[0],&F2c[0],&npts,&ichk);

double rflux(double x_pom, double a0, double ap, double b0);


#include <iostream>

using namespace std;

const std::vector<std::string> AQcdnumDDISCS::fRequirements = {"QcdnumInit",
							       "Aq","Bq","Cq","Ag","Bg","Cg", //pdf parameters
                                                               "a0_IP", "ap_IP", "b0_IP",     //Pomeron flux
                                                               "a0_IR", "ap_IR", "b0_IR",     //Reggeon flux
                                                               "n_IR",                        //Reggeon suppression
                                                                }; //< List of all AParm's which this function depends on
const std::vector<std::string> AQcdnumDDISCS::fStopFurtherNotification = {"blubb", "par4"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AQcdnumDDISCS::fFunctionName = "QcdnumDDISCS"; //< The function's name


// Todo (possibly):
// AQcdnumDDISCS constructor should made 'private'
// Only 'TheoryHandler' is allowed to 'construct' AParmFunc's
// and then also calls 'ARegisterRequirements()' and 'Init()'.
// Or this is anyhow done for ALL AParmBase's by the theory handler, but AParmConst has Init()={} and GetReq()={};

// ___________________________________________________________________________________________ //
AQcdnumDDISCS::AQcdnumDDISCS(const std::string& name) : AParmFuncBase<double>(name) { 
   //ARegisterRequirements(this); // needed in every constructor
}


// ___________________________________________________________________________________________ //
AQcdnumDDISCS::~AQcdnumDDISCS() {
}


// ___________________________________________________________________________________________ //
bool AQcdnumDDISCS::Init() {
   debug["Update"] << "Init AQcdnumDDISCS. Nothing to do." << endl;
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ___________________________________________________________________________________________ //
bool AQcdnumDDISCS::Update() {
   debug["Update"]<<" AQcdnumDDISCS::Update(). GetAlposName:" <<GetAlposName()<<endl;

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit);
   
   vector<double> xpom     = DOUBLE_COL_NS(Data,xp,GetAlposName());
   vector<double> q2       = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   vector<double> beta     = DOUBLE_COL_NS(Data,beta,GetAlposName());
   vector<double> sigmaVec = DOUBLE_COL_NS(Data,Sigma,GetAlposName());


   // cout << "q2Vector size " << q2.size() << endl;
   // for(int i = 0; i < xpom.size(); ++i)
   //     cout << xpom[i] <<" "<<  q2[i] <<" "<< beta[i] <<" "<< sigmaVec[i]<<  endl;

   fValue.resize(q2.size());
   fError.resize(q2.size());

   //fValue[0] = PAR(Aq);
   //fValue[1] = PAR(Bq);
   //fValue[2] = PAR(Cq);
   //fValue[3] = PAR(Ag);
   //fValue[4] = PAR(Bg);
   //fValue[5] = PAR(Cg);


   const double a0_IP = PAR(a0_IP);
   const double ap_IP = PAR(ap_IP);
   const double b0_IP = PAR(b0_IP);

   const double a0_IR = PAR(a0_IR);
   const double ap_IR = PAR(ap_IR);
   const double b0_IR = PAR(b0_IR);

   const double n_IR  = PAR(n_IR);


   int iFL=1, iF2=2, nchk=0;
   vector<double> FLl(q2.size(),0.), F2l(q2.size(),0.);
   vector<double> FLc(q2.size(),0.), F2c(q2.size(),0.);
   vector<double> FLb(q2.size(),0.), F2b(q2.size(),0.);

   int npts=int(q2.size());

   //Quark's charges^2
   vector<double> CEP2F = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
   for(int i=0; i<13; i++) CEP2F[i] /= 9;

   int ichk = 0;
   //Light
   zmstfun_(&iF2,&CEP2F[0],&beta[0],&q2[0],&F2l[0], &npts,&ichk);
   zmstfun_(&iFL,&CEP2F[0],&beta[0],&q2[0],&FLl[0], &npts,&ichk);

   if(PAR(QcdnumInit.nfFix) >= 3) {
       int icharm = 1, ibottom = -2;
       //charm
       hqstfun_(&iF2, &icharm,&CEP2F[0],&beta[0],&q2[0],&F2c[0],&npts,&ichk);
       hqstfun_(&iFL, &icharm,&CEP2F[0],&beta[0],&q2[0],&FLc[0],&npts,&ichk);

       //bottom
       hqstfun_(&iF2, &ibottom,&CEP2F[0],&beta[0],&q2[0],&F2b[0],&npts,&ichk);//no check on nf = 4
       hqstfun_(&iFL, &ibottom,&CEP2F[0],&beta[0],&q2[0],&FLb[0],&npts,&ichk);//no check on nf = 4
   }

   // { // ugly check:a
   //    vector<double> qqq{1.75,8.5,20,800};
   //    int iset=1;//5; // iset: 5 is external PDF                                                //    double xp  = 0.1;
   //    int ichk;
   //    int nul = 0;
   //    vector<double> xfx(13);
   //    for ( auto qq : qqq ) {
   // 	 double muf2 = qq;
   // 	 allfxq_( &iset, &xp, &muf2, &xfx[0], &nul, &ichk );
   // 	 cout<<"q2= "<<muf2<<": ";
   // 	 for ( auto p : xfx ) cout<<"\t"<<p;
   // 	 cout<<endl;
   //    }
   // }
   // //exit(3);




   const double mp2 = pow(0.92, 2);

   // ------ calc reduced CS
   for (unsigned int i =0; i<q2.size(); i++) {

       double Ep = (q2[i] < 120) ? 820 : 920;
       double Ee = 27.5;
       const double s = 4*Ep * Ee;

       double x = beta[i]*xpom[i];
       double y = q2[i]/(s-mp2)/x;

       double yplus  = 1+pow(1-y,2);
       //cout<<"QCDNUM  Q2="<<q2[i]<<"\tf2="<<F2[i]<<"\tfl="<<FL[i]<<"\tf3="<<xF3[i]<<"\tq="<<charge<<"\tpol="<<polty<<endl;

       double F2 = F2l[i] + F2c[i] + F2b[i];
       double FL = FLl[i] + FLc[i] + FLb[i];

       //Pomeron flux
       double flxIP = rflux(xpom[i], a0_IP, ap_IP, b0_IP);

       //Reduced x-section for pomeron
       double xpSigRed_IP =  flxIP*xpom[i] * (F2  - y*y/yplus*FL);

       //cout<<"QCDNUM  Q2="<<q2[i]<<"\tf2*flxIP="<<F2*flxIP<<"\tfl*flx="<<FL*flxIP<<endl;


       //Reggeon flux
       double flxIR = rflux(xpom[i], a0_IR, ap_IR, b0_IR);

       //Get the Reggeon structure function from the H12006
       static bool isFirst = true; //hack for faster calculation
       int ifit = 0;
       if(isFirst) { ifit = 1; isFirst = false; } //Reggeon should be the same for both FitA and FitB

       double xPq[13];
       double f2FitA[2], flFitA[2]; //0 = pomeron, 1 = reggeon
       double c2FitA[2], clFitA[2];
       qcd_2006_(&beta[i], &q2[i],  &ifit, xPq, f2FitA, flFitA, c2FitA, clFitA);
       double F2r = f2FitA[1];
       double FLr = flFitA[1];

       //Reduced x-section for reggeon
       double xpSigRed_IR =  flxIR*xpom[i] * (F2r  - y*y/yplus*FLr);


       fValue[i] = xpSigRed_IP + n_IR*xpSigRed_IR;


       // double sRedP = (f2[0]) - y*y/(1 + pow(1-y,2)) * (fl[0]);
       // double sRedR = (f2[1]) - y*y/(1 + pow(1-y,2)) * (fl[1]);
       // double sRedC = (c2[0]) - y*y/(1 + pow(1-y,2)) * (cl[0]); //should not be included
       // fValue[i]= xpom*(sRedP+sRedR);


       //double redfac = (2*M_PI*(1./137)*(1./137) / (x*q2[i]*q2[i]) );

       // cout<<"fValue="<<fValue[i]<<"\tsigma="<<sigmaVec[i]<<"\tratio="<<fValue[i]/sigmaVec[i]
       // 	   <<"\tq2="<<q2[i]
       // 	   <<"\tflxIP="<<flxIP
       // 	   <<"\tx="<<x
       // 	   <<"\ty="<<y
       // 	   <<"\tb="<<beta[i]
       // 	   <<"\txp="<<xpom[i]
       // 	  //<<"\tredfac = "<<redfac<<"\t inv="<<1./redfac
       // 	   <<endl;


   }

   return true;
}


// ___________________________________________________________________________________________ //
//
//
//


//tcut is negative : tcut = -1
static double rfluxRaw(double x_pom, double a0, double ap, double b0, double tcut)
{
    const double mp = 0.93827231;

    //     calc min. kinematically  allowed t
    double tmin= -pow(mp*x_pom,2)/(1.-x_pom);

    //     c*xpom**(-(2apom-1))
    double fl =  exp((2.0*a0-1.)*log(1.0/x_pom));
    double b=(b0+2.0*ap*log(1.0/x_pom));

    //   at fixed t:  exp(Bt)
    //  fl = fl * exp(b*tcut);

    //   t-integrated: (1/B)*[exp(-B*tmax)-exp(-B*tmin)]
    fl = fl * (exp(tmin*b)-exp(tcut*b))/b;

    return fl;
}



double rflux(double x_pom, double a0, double ap, double b0)
{
    double tcut = -1;
    double xPomNorm = 0.003;
    const double dm =  rfluxRaw(xPomNorm, a0, ap, b0, tcut);
    double  norm=(1./(xPomNorm*dm)); //xpom * flux normalized to 1 at xpom = 0.003

    return  norm * rfluxRaw(x_pom, a0, ap, b0, tcut);
}



