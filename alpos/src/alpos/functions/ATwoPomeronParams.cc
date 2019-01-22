
#include "alpos/functions/ATwoPomeronParams.h"

#include <iostream>

using namespace std;

const std::vector<std::string> ATwoPomeronParams::fRequirements = {
   "Q2","y","x","s","mp",
   "a0","b0","e0","a1","b1","e1","eta0a","eta1a","eta0b","eta1b", // set these parameters unequal to 0 if parameterisation should be used
   "a0_q2_0p045","b0_q2_0p045","e0_q2_0p045","a1_q2_0p045","b1_q2_0p045","e1_q2_0p045", // Q2=0.045
   "a0_q2_0p065","b0_q2_0p065","e0_q2_0p065","a1_q2_0p065","b1_q2_0p065","e1_q2_0p065", // Q2=0.065
   "a0_q2_0p085","b0_q2_0p085","e0_q2_0p085","a1_q2_0p085","b1_q2_0p085","e1_q2_0p085", // Q2=0.085
   "a0_q2_0p11","b0_q2_0p11","e0_q2_0p11","a1_q2_0p11","b1_q2_0p11","e1_q2_0p11", // Q2=0.11
   "a0_q2_0p15","b0_q2_0p15","e0_q2_0p15","a1_q2_0p15","b1_q2_0p15","e1_q2_0p15", // Q2=0.15
   "a0_q2_0p20","b0_q2_0p20","e0_q2_0p20","a1_q2_0p20","b1_q2_0p20","e1_q2_0p20", // Q2=0.2
   "a0_q2_0p25","b0_q2_0p25","e0_q2_0p25","a1_q2_0p25","b1_q2_0p25","e1_q2_0p25", // Q2=0.25
   "a0_q2_0p35","b0_q2_0p35","e0_q2_0p35","a1_q2_0p35","b1_q2_0p35","e1_q2_0p35", // 
   "a0_q2_0p40","b0_q2_0p40","e0_q2_0p40","a1_q2_0p40","b1_q2_0p40","e1_q2_0p40", // 
   "a0_q2_0p50","b0_q2_0p50","e0_q2_0p50","a1_q2_0p50","b1_q2_0p50","e1_q2_0p50", // 
   "a0_q2_0p65","b0_q2_0p65","e0_q2_0p65","a1_q2_0p65","b1_q2_0p65","e1_q2_0p65", // 
   "a0_q2_0p85","b0_q2_0p85","e0_q2_0p85","a1_q2_0p85","b1_q2_0p85","e1_q2_0p85", // 
   "a0_q2_1p2","b0_q2_1p2","e0_q2_1p2","a1_q2_1p2","b1_q2_1p2","e1_q2_1p2", // 
   "a0_q2_1p5","b0_q2_1p5","e0_q2_1p5","a1_q2_1p5","b1_q2_1p5","e1_q2_1p5", // 
   "a0_q2_2p0","b0_q2_2p0","e0_q2_2p0","a1_q2_2p0","b1_q2_2p0","e1_q2_2p0", // 
   "a0_q2_2p5","b0_q2_2p5","e0_q2_2p5","a1_q2_2p5","b1_q2_2p5","e1_q2_2p5", // 
   "a0_q2_2p7","b0_q2_2p7","e0_q2_2p7","a1_q2_2p7","b1_q2_2p7","e1_q2_2p7", // 
   "a0_q2_3p5","b0_q2_3p5","e0_q2_3p5","a1_q2_3p5","b1_q2_3p5","e1_q2_3p5", // 
   "a0_q2_4p5","b0_q2_4p5","e0_q2_4p5","a1_q2_4p5","b1_q2_4p5","e1_q2_4p5", // 
   "a0_q2_6p5","b0_q2_6p5","e0_q2_6p5","a1_q2_6p5","b1_q2_6p5","e1_q2_6p5", // 
   "a0_q2_8p5","b0_q2_8p5","e0_q2_8p5","a1_q2_8p5","b1_q2_8p5","e1_q2_8p5", // 
   "a0_q2_10","b0_q2_10","e0_q2_10","a1_q2_10","b1_q2_10","e1_q2_10", // 
   "a0_q2_12","b0_q2_12","e0_q2_12","a1_q2_12","b1_q2_12","e1_q2_12", // 
}; //< List of all AParm's which this function depends on
const std::vector<std::string> ATwoPomeronParams::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string ATwoPomeronParams::fFunctionName = "TwoPomeronParams"; //< The function's name


// ___________________________________________________________________________________________ //
ATwoPomeronParams::ATwoPomeronParams(const std::string& name) : AParmFuncBase<double>(name) { 
   //ARegisterRequirements(this); // needed in every constructor
}


// ___________________________________________________________________________________________ //
ATwoPomeronParams::~ATwoPomeronParams() {
}


// ___________________________________________________________________________________________ //
bool ATwoPomeronParams::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   //!
   //! Do not access values of other input 'functions' here.
   //! 
   return true;
}


// ___________________________________________________________________________________________ //
bool ATwoPomeronParams::Update() {
   //! The 'Update()' function must fill the member variable fValue
   //! which should fit the cross sections
   //!


   // ---
   double aHat0 = 0;
   double bHat0 = 0;
   double eps0  = 0;
   double aHat1 = 0;
   double bHat1 = 0;
   double eps1  = 0;

   // --- input values
   double q2 = PAR(Q2);
   double x  = PAR(x);
   double y  = PAR(y);
   double ss = PAR(s);
   double mp = PAR(mp);
   
   double mp2= mp*mp;
   
   aHat0 = ( PAR(a0)!=0 ) ? 
      //PAR(a0)/q2*pow(2*mp*mp/q2,1+PAR(eta0a)) :
      PAR(a0)*pow(2*mp2/q2,PAR(eta0a))/(q2*q2/mp2) :
      GetParValue("a0",q2);

   aHat1 = ( PAR(a1)!=0 ) ?
      //PAR(a1)/q2*pow(2*mp*mp/q2,1+PAR(eta1a)) :
      PAR(a1)*pow(2*mp2/q2,PAR(eta1a))/(q2*q2/mp2) :
      GetParValue("a1",q2);

   bHat0 = ( PAR(b0)!=0 ) ? 
      //PAR(b0)*pow(2*mp*mp/q2,1+PAR(eta0b)) :
      PAR(b0)*pow(2*mp2/q2,PAR(eta0b))/(q2/mp2) :
      GetParValue("b0",q2);

   bHat1 = ( PAR(b1)!=0 ) ? 
      PAR(b1)*pow(2*mp2/q2,PAR(eta1b))/(q2/mp2) :
      GetParValue("b1",q2);

   eps0 = PAR(e0);//GetParValue("e0",q2);

   eps1 = PAR(e1);//GetParValue("e1",q2);
   
   // --- return values
   fValue.resize(6);
   fError.resize(6);
   fValue[0] =  aHat0;
   fValue[1] =  bHat0;
   fValue[2] =  eps0 ;
   fValue[3] =  aHat1;
   fValue[4] =  bHat1;
   fValue[5] =  eps1 ;

  
   // --- calculation succeeded?
   return true;
}

// ___________________________________________________________________________________________ //
double ATwoPomeronParams::GetParValue(const string& tag, double q2) {
   // return parameter

   string q2s ="";

   if      (  AreSame(q2, 0.045 ) ) q2s = "0p045";
   else if (  AreSame(q2, 0.065 ) ) q2s = "0p065";
   else if (  AreSame(q2, 0.085 ) ) q2s = "0p085";
   else if (  AreSame(q2, 0.11  ) ) q2s = "0p11"; 
   else if (  AreSame(q2, 0.15  ) ) q2s = "0p15"; 
   else if (  AreSame(q2, 0.20  ) ) q2s = "0p20"; 
   else if (  AreSame(q2, 0.25  ) ) q2s = "0p25"; 
   else if (  AreSame(q2, 0.35  ) ) q2s = "0p35"; 
   else if (  AreSame(q2, 0.40  ) ) q2s = "0p40"; 
   else if (  AreSame(q2, 0.50  ) ) q2s = "0p50"; 
   else if (  AreSame(q2, 0.65  ) ) q2s = "0p65"; 
   else if (  AreSame(q2, 0.85  ) ) q2s = "0p85"; 
   else if (  AreSame(q2, 1.2   ) ) q2s = "1p2";  
   else if (  AreSame(q2, 1.5   ) ) q2s = "1p5";  
   else if (  AreSame(q2, 2.0   ) ) q2s = "2p0";  
   else if (  AreSame(q2, 2.5   ) ) q2s = "2p5";  
   else if (  AreSame(q2, 2.7   ) ) q2s = "2p7";  
   else if (  AreSame(q2, 3.5   ) ) q2s = "3p5";  
   else if (  AreSame(q2, 4.5   ) ) q2s = "4p5";  
   else if (  AreSame(q2, 6.5   ) ) q2s = "6p5";  
   else if (  AreSame(q2, 8.5   ) ) q2s = "8p5";  
   else if (  AreSame(q2,10.    ) ) q2s = "10";   
   else if (  AreSame(q2,12.    ) ) q2s = "12";   

   if ( q2s=="" ) {
      // warn["GetValue"]<<"Could not identify q2 value. Q2="<<q2<<"."<<endl;
      // warn["GetValue"]<<"Returning value for Q2=12GeV2"<<endl;
      q2s="12";
   }
   
   //return PAR(tag+"_q2_"+q2s);
   return TheoryHandler::Handler()->GetParmD(this->GetAlposName()+std::string(".")+std::string(tag+"_q2_"+q2s))->GetValue(); // PAR(xy) not working
}


// ___________________________________________________________________________________________ //
