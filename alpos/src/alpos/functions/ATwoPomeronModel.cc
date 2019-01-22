
#include "alpos/functions/ATwoPomeronModel.h"
#include "fastnlotk/read_steer.h"

#include <iostream>

using namespace std;

const std::vector<std::string> ATwoPomeronModel::fRequirements = {"beta0","alphaPrime0","beta1","alphaPrime1","mp","TwoPomeronPars"}; //< List of all AParm's which this function depends on
const std::vector<std::string> ATwoPomeronModel::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ATwoPomeronModel::fFunctionName = "TwoPomeronModel"; //< The function's name


// ___________________________________________________________________________________________ //
ATwoPomeronModel::ATwoPomeronModel(const std::string& name) : AParmFuncBase<double>(name) { 
   //ARegisterRequirements(this); // needed in every constructor
}


// ___________________________________________________________________________________________ //
ATwoPomeronModel::~ATwoPomeronModel() {
}


// ___________________________________________________________________________________________ //
bool ATwoPomeronModel::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Do not access values of other input 'functions' here.
   //! 
   return true;
}


// ___________________________________________________________________________________________ //
bool ATwoPomeronModel::Update() {
   //! The 'Update()' function must fill the member variable fValue
   //! which should fit the cross sections
   //!

   // --- Access q2, x and y _directly_ from the data-files
   vector<double> q2 = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   vector<double> x  = DOUBLE_COL_NS(Data,x,GetAlposName());
   vector<double> y  = DOUBLE_COL_NS(Data,y,GetAlposName());
   double sqs        = DOUBLE_NS(sqrt-s,GetAlposName());
   double ss         = sqs*sqs;

   // --- calculate y (if needed)
   if ( y.empty() ) {
      y.resize(q2.size());
      for ( unsigned int i = 0 ; i<q2.size() ; i++ ) 
         y[i] = q2[i] / (x[i]*ss);
   }
   // --- get other parameters from steering (or from alpos)
   double charge = DOUBLE_NS(e-charge,GetAlposName()); //PAR(e-charge);//
   double polty  = DOUBLE_NS(e-polarity,GetAlposName()); //PAR(e-polarity) ;
   bool IsRedCS  = BOOL_NS(IsReducedCS,GetAlposName());
   bool IsNC     = BOOL_NS(IsNC,GetAlposName());
   // ---- sanity check
   if ( charge==0 ) {
      error["Update"]<<"Could not get charge of lepton."<<endl;
      exit(1);
   }

   // --- resize return vectors
   fValue.resize(q2.size()); 
   fError.resize(q2.size());

   // --- get parameters for calculation from alpos
   double mp     = PAR(mp);
   double beta0  = PAR(beta0);  // soft IP
   double aP0    = PAR(alphaPrime0); // soft IP
   double beta1  = PAR(beta1); // hard IP
   double aP1    = PAR(alphaPrime1); // hard IP

   // --- perform calculation of 
   if ( IsNC ) {
      for ( unsigned int i =0 ; i<q2.size() ; i++ ) {
	 // --- get pomeron parameters
	 SET(TwoPomeronPars.Q2,q2[i],0);
	 SET(TwoPomeronPars.x,x[i],0);
	 SET(TwoPomeronPars.y,y[i],0);
	 SET(TwoPomeronPars.s,ss,0);
	 const vector<double>& IPparm = VALUES(TwoPomeronPars);

	 double aHat0 = IPparm[0];
	 double bHat0 = IPparm[1];
	 double eps0  = IPparm[2];
	 double aHat1 = IPparm[3];
	 double bHat1 = IPparm[4];
	 double eps1  = IPparm[5];
	 
	 // calculate individual factors
	 const double W  = q2[i]*(1./x[i]-1);
	 const double B0 = 3*beta0*pow(W*W*aP0,eps0)*cos(M_PI/2*eps0);
	 const double B1 = 3*beta1*pow(W*W*aP1,eps1)*cos(M_PI/2*eps1);

	 // --- calculate W1 & W2
	 const double mp2 = mp*mp;
	 const double pq  = (W*W+q2[i]-mp*mp)/2; // not squared but, (pq)*2
	 const double pq2 = pq*2; // not squared but, (pq)*2


	 const double W1_0 = B0*( (bHat0-2*q2[i]*aHat0)*pow(pq2,2) + (bHat0-q2[i]*aHat0)*2*q2[i]*mp*mp);
	 const double W1_1 = B1*( (bHat1-2*q2[i]*aHat1)*pow(pq2,2) + (bHat1-q2[i]*aHat1)*2*q2[i]*mp*mp );
	 const double W1 = 1/(2*M_PI*mp*W*W)*(W1_0+W1_1);
	 const double W2 = mp/(2*M_PI*W*W)*( B0*4*q2[i]*bHat0 + B1*4*q2[i]*bHat1 );

 
	 // --- calculate F1 & F2
	 // const double F2_IP0 = mp*B0*4*bHat0*x[i];	
	 // const double F2_IP1 = mp*B1*4*bHat1*x[i];
	 // double F2 = mp*mp/M_PI/W/W*(F2_IP0+F2_IP1);

	 // const double FL_IP0 = mp*B0*( 
	 //    (bHat0-2*q2[i]*aHat0)*q2[i]*q2[i]/pow(mp,4)/x[i] + 
	 //    (bHat0-q2[i]*q2[i]*aHat0)*q2[i]*q2[i]*x[i]/mp/mp -
	 //    4*bHat0*x[i]);
	 // const double FL_IP1 = mp*B1*( 
	 //    (bHat1-2*q2[i]*aHat1)*q2[i]*q2[i]/pow(mp,4)/x[i] + 
	 //    (bHat1-q2[i]*q2[i]*aHat1)*q2[i]*q2[i]*x[i]/mp/mp -
	 //    4*bHat1*x[i]);

	 // double FL = mp*mp/M_PI/W/W*(FL_IP0+FL_IP1);
	 
	 // --- check F2 and FL from W or directly
	 double F2_fromW = pq/mp2*W2;
	 //double FL_fromW = pq/mp2*W2 - q2[i]/pq*W1;
	 double FL_fromW = F2_fromW - 2*x[i]*W1;


	 // --- calculate sigma-reduced
	 double sigred = F2_fromW - y[i]*y[i]/(1+(1-y[i])*(1-y[i]))*FL_fromW;

	 // double sigred from eq14
	 double YY = y[i]*y[i]/(1+(1-y[i])*(1-y[i]));
	 double sr_0 = mp*B0*(q2[i]*bHat0/mp2 - 2*YY*q2[i]*q2[i]/mp2*aHat0);
	 double sr_1 = mp*B1*(q2[i]*bHat1/mp2 - 2*YY*q2[i]*q2[i]/mp2*aHat1);
	 double sigred14 = (W*W+q2[i])/(M_PI*W*W)*( sr_0+sr_1);

	 //cout<<"sigred="<<sigred<<"\t2: "<<sigred14<<endl;
	 // --- pass sigma reduced back to alpos
	 //fValue[i] = sigred;
	 fValue[i] = sigred14;
	 if ( !IsRedCS ) {
	    error["Update"]<<"todo. uncalculated sig-reduced."<<endl;
	    // double aem = 7.29735d-3;
	    // double convfac = ...;
	    // fValue[i] *= 2*M_PI/(x[i]*q2[i]*q2[i])*convfac*aem*aem;
	    exit(1);
	 }
      }
   }
   else {
      error["Update"]<<"CC not implemented."<<endl;
      exit(1);
   }
  
   // --- calculation succeeded?
   return true;
}


// ___________________________________________________________________________________________ //
