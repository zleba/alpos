
#include "alpos/functions/AStrowp1.h"

#include <iostream>

extern "C" void strowp1_(double* X, double* SCALE, double* UPV, double* DNV, double* SEA, double* STR, double* CHM, double* GL);

using namespace std;

const std::vector<std::string> AStrowp1::fRequirements = {"xp", "muf", // x and muf values: Usually integration constants
                                                         }; //< List of all AParm's which this function depends on
const std::vector<std::string> AStrowp1::fStopFurtherNotification = {"xp","muf"}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AStrowp1::fFunctionName = "STROWP1"; //< The function's name


// __________________________________________________________________________________________ //
AStrowp1::AStrowp1(const std::string& name) : AParmFuncBase<double>(name) { 
   SetClassName("AStrowp1");
   fValue.resize(13);
   fError.resize(13);
}


// __________________________________________________________________________________________ //
AStrowp1::~AStrowp1() {
}


// ___________________________________________________________________________________________ //
bool AStrowp1::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   return true;
}


// ______________________________________________________________________________________ //
std::vector<double> AStrowp1::GetQuick(int n, ...) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(2,xp,Q);
   std::vector<double> ret(fValue.size());
   va_list ap;
   va_start(ap, n); /* Requires the last fixed parameter (to get the address) */
   double xp  = va_arg(ap, double);
   double muf = va_arg(ap, double);
   va_end(ap);

   return GetQuick({xp,muf});
    
}


// ______________________________________________________________________________________ //
std::vector<double> AStrowp1::GetQuick(const vector<double>& xp_muf) {
   //! The possibilty to implement a quick access without changing of any parameters
   //! Use the quick access to calculate alpha_s(mur) using:
   //!   ::GetQuick(vector<double> xp_mur;
   //!   Input parameters must be:
   //!   xp_mur[0] = xp
   //!   xp_muf[0] = Q

   double x   = xp_muf[0];
   double muf = xp_muf[1];
   // --- cache
   vector<double>& rr = fValCache[{x,muf}];
   if ( rr.size() ) return rr;
   // ---

   if ( xp_muf.size() != 2) {
      cout<<"Error in AStrowp1::GetQuick(vector). Quick acces is implemented for two parameter which are 'xp' and 'muf'."<<endl;
      return vector<double>(fValue.size());
   }

   double UPV, DNV, SEA, STR, CHM, GL;
   strowp1_(&x, &muf, &UPV, &DNV, &SEA, &STR, &CHM, &GL);
   // C...Put into output arrays.
   //    UPV = XQ(1)
   //    DNV = XQ(1)
   //    SEA = XQ(3)/SIXD
   //    STR = XQ(3)/SIXD
   //    CHM = XQ(4)
   //         BOT = ZEROD
   //         TOP = ZEROD
   //    GL  = XQ(2)

   //tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
   std::vector<double> ret(fValue.size());
   ret[0] = 0;
   ret[1] = 0;
   ret[2] = CHM;
   ret[3] = STR;
   ret[4] = SEA;
   ret[5] = SEA;

   ret[6] = GL;
 
   ret[7] = DNV+SEA; 
   ret[8] = UPV+SEA; 
   ret[9] = STR;
   ret[10] = CHM;
   ret[11] = 0;
   ret[12] = 0;

   fValCache[{x,muf}] = ret;// cache

   return ret;
  
}


// __________________________________________________________________________________________ //
bool AStrowp1::Update() {
   debug["Update"]<<"GetAlposName:" <<GetAlposName()<<endl;
   //fValue.resize(GetRequirements().size());
   //fError.resize(GetRequirements().size());

   vector<double> xp_muf{PAR(xp),PAR(muf)} ;
   fValue = GetQuick(xp_muf);

   return true;
}

