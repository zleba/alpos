
#include "alpos/functions/AQcdnumDDISCS.h"
#include "alpos/functions/AQcdnumInit.h"


extern "C" void zmstfun_(int* istf, double* def, double* x, double* Q2, double* f, int* n, int* nchk);
extern "C" void hqstfun_(int* istf, int *iflavour, double* def, double* x, double* Q2, double* f, int* n, int* nchk);

  //hqstfun_(iF2, &icharm,&CEP2F[0],&beta[0],&q2[0],&F2c[0],&npts,&ichk);

double rflux(double x_pom, double a0, double ap, double b0);


#include <iostream>

using namespace std;

const std::vector<std::string> AQcdnumDDISCS::fRequirements = {"Aq","Bq","Cq","Ag","Bg","Cg", //pdf parameters
                                                               "a0_IP", "ap_IP", "b0_IP",     //Pomeron flux
                                                               "a0_IR", "ap_IR", "b0_IR",     //Reggeon flux
                                                               "n_IR"                         //Reggeon suppression
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
   cout<<"AQcdnumDDISCS Constructor."<<endl;
   //ARegisterRequirements(this); // needed in every constructor
}


// ___________________________________________________________________________________________ //
AQcdnumDDISCS::~AQcdnumDDISCS() {
}


// ___________________________________________________________________________________________ //
bool AQcdnumDDISCS::Init() {
    cout << "Init AQcdnumDDISCS" << endl;
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ___________________________________________________________________________________________ //
bool AQcdnumDDISCS::Update() {
   cout<<" AQcdnumDDISCS::Update(). GetAlposName:" <<GetAlposName()<<endl;

   // 'Update' PDF and Alpha_s values to ensure that 'Quick'-access are correct.
   UPDATE(QcdnumInit); 

   fValue.resize(GetRequirements().size());
   fError.resize(GetRequirements().size());

   vector<double> xpom   = DOUBLE_COL_NS(Data,xp,GetAlposName());
   vector<double> q2   = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   vector<double> beta = DOUBLE_COL_NS(Data,beta,GetAlposName());
   vector<double> sigmaVec = DOUBLE_COL_NS(Data,Sigma,GetAlposName());
   cout << "q2Vector size " << q2.size() << endl;

   for(int i = 0; i < xpom.size(); ++i)
       cout << xpom[i] <<" "<<  q2[i] <<" "<< beta[i] <<" "<< sigmaVec[i]<<  endl;


   //fValue[0] = PAR(Aq);
   //fValue[1] = PAR(Bq);
   //fValue[2] = PAR(Cq);
   //fValue[3] = PAR(Ag);
   //fValue[4] = PAR(Bg);
   //fValue[5] = PAR(Cg);


   const double a0_IP = PAR(a0_IP);
   const double ap_IP = PAR(ap_IP);
   const double b0_IP = PAR(b0_IP);



   int iFL=1, iF2=2, nchk=0;
   vector<double> FLl(q2.size()), F2l(q2.size());
   vector<double> FLc(q2.size()), F2c(q2.size());
   vector<double> FLb(q2.size()), F2b(q2.size());

   int npts=int(q2.size());

   //Quark's charges^2
   vector<double> CEP2F = { 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. };
   for(int i=0; i<13; i++) CEP2F[i] /= 9;

   int ichk = 0;
   //Light
   zmstfun_(&iF2,&CEP2F[0],&beta[0],&q2[0],&F2l[0], &npts,&ichk);
   zmstfun_(&iFL,&CEP2F[0],&beta[0],&q2[0],&FLl[0], &npts,&ichk);

   int icharm = 1, ibottom = -2;

   //charm
   hqstfun_(&iF2, &icharm,&CEP2F[0],&beta[0],&q2[0],&F2c[0],&npts,&ichk);
   hqstfun_(&iFL, &icharm,&CEP2F[0],&beta[0],&q2[0],&FLc[0],&npts,&ichk);

   //bottom
   hqstfun_(&iF2, &ibottom,&CEP2F[0],&beta[0],&q2[0],&F2b[0],&npts,&ichk);//no check on nf = 4
   hqstfun_(&iFL, &ibottom,&CEP2F[0],&beta[0],&q2[0],&FLb[0],&npts,&ichk);//no check on nf = 4



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

       fValue[i] = flxIP * (F2  - y*y/yplus*FL);
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



