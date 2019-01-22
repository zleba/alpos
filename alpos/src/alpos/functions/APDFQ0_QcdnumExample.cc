#include "alpos/functions/APDFQ0_QcdnumExample.h"

#include <iostream>
#include <cmath>

using namespace std;


// ______________________________________________________________________________________ //
const std::vector<std::string> APDFQ0_QcdnumExample::fRequirements = {"xp","Q0","iPDF"}; //< List of all AParm's which this function depends on
const std::vector<std::string> APDFQ0_QcdnumExample::fStopFurtherNotification = {"mur"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string APDFQ0_QcdnumExample::fFunctionName = "PDFQ0_QcdnumExample"; //< The function's name


// ______________________________________________________________________________________ //
APDFQ0_QcdnumExample::APDFQ0_QcdnumExample(const std::string& name) : AParmFuncBase<double>(name) { 

}


// ______________________________________________________________________________________ //
APDFQ0_QcdnumExample::~APDFQ0_QcdnumExample() {
}


// ___________________________________________________________________________________________ //
bool APDFQ0_QcdnumExample::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   fValue.resize(1);
   fError.resize(1);
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_QcdnumExample::GetQuick(int n,...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(1, double mur);

   cout<<"Error in APDFQ0_QcdnumExample::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
std::vector<double> APDFQ0_QcdnumExample::GetQuick(const vector<double>& mur) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> mur);

   cout<<"Error in APDFQ0_QcdnumExample::GetQuick(vector<double>). Quick acces is not implemented. Exiting."<<endl;
   exit(1);
   return vector<double>();
}


// ______________________________________________________________________________________ //
bool APDFQ0_QcdnumExample::Update() {

   int ipdf = PAR(iPDF);
   
   if ( ipdf == -1 ) {
      // return QCDNUM vector 'def'
      vector<double> def = {
	 //tb  bb  cb  sb  ub  db  g   d   u   s   c   b   t                                          
	 //-6  -5  -4  -3  -2  -1  0   1   2   3   4   5   6                                          
	 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,                                        
	 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,                                        
	 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,                                        
	 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,                                        
	 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0.,                                        
	 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.                                         
      };
      fValue = def;
      fValue.resize(13*12);
      fError.resize(fValue.size());
   }
   else { 
      fValue.resize(1);
      fError.resize(1);
      
      double xp = PAR(xp);
      double q0 = PAR(Q0);
      if(ipdf== 0) fValue[0] = 1.7 * pow(xp,-0.1) * pow(1.-xp,5);
      else if(ipdf== 1) fValue[0] = 3.064320*pow(xp,0.8) * pow(1.-xp,4.);
      else if(ipdf== 2) fValue[0] = 5.107200*pow(xp,0.8)* pow(1.-xp,3.);
      else if(ipdf== 3) fValue[0] = 0.;
      else if(ipdf== 4) fValue[0] =  0.1939875 * pow(xp,-0.1) * pow(1.-xp,6.);
      else if(ipdf== 5) fValue[0] = (0.1939875 * pow(xp,-0.1) * pow(1.-xp,6.))*(1-xp);
      else if(ipdf== 6) fValue[0] = 0.2 * ((0.1939875 * pow(xp,-0.1) * pow(1.-xp,6.))+((0.1939875 * pow(xp,-0.1) * pow(1.-xp,6.))*(1-xp)));
      else fValue[0] = 0;

   }
   return true;
}

// ______________________________________________________________________________________ //
