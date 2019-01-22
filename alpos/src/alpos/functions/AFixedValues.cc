#include "alpos/functions/AFixedValues.h"
#include <iostream>
#include <fstream>

using namespace std;

const std::vector<std::string> AFixedValues::fRequirements = {"filename","format","tag"}; //< List of all AParm's which this function depends on. These must be specified in the steering
const std::vector<std::string> AFixedValues::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AFixedValues::fFunctionName = "FixedValues"; //< The function's name


// ___________________________________________________________________________________________ //
AFixedValues::AFixedValues(const std::string& name) : AParmFuncBase<double>(name) { 
}


// ___________________________________________________________________________________________ //
AFixedValues::~AFixedValues() {
}


// ___________________________________________________________________________________________ //
bool AFixedValues::Init() {
   //! Init is once called for each function
   //! Initialize all needed member variables (typically nothing is needed here)
   //!
   //! return true if initialization was successful.
   //! 
   debug["Init"]<<"GetAlposName:" <<GetAlposName()<<endl;
   
   string filename = PAR_S(filename);
   string format   = PAR_S(format);
   string tag = PAR_S(tag);

   CONST(filename);
   CONST(format);
   CONST(tag);

   info["Init"]<<"Reading fixed values from file: "<<filename<<endl;
   info["Init"]<<"Reading file in format: "<<format<<endl;
   
   if ( format=="plain" ) {
      double q2in, xin,yin,csin,thin;
      double tmp;
      ifstream ds(filename.c_str());
      if(!ds.is_open()) {
	 error["Init"]<<"failed to open file\n";
	 exit(2);
      }
      //char weg[1000];
      //ds.getline(weg,1000);
      fValue.clear();
      while(!ds.eof()) {
	 ds>>tmp;
	 if(ds.fail()) {
	    info["Init"]<<"Read "<<fValue.size()<<" values."<<endl;
	    break;
	 }
	 fValue.push_back(tmp);
      }
      fError.resize(fValue.size());
   }
   else if (format=="colum") {
      
   }
   else {
      error["Init"]<<"Did not recognize format of filename."<<endl;
      exit(1);
   }

   return true;
}


// ___________________________________________________________________________________________ //
bool AFixedValues::Update() {
   //!
   //! Nothing to be done here. Since fValues are expected to be constant
   //! 
   debug["Update"]<<"AlposName:" <<GetAlposName()<<endl;
   return true;
}


// ___________________________________________________________________________________________ //
