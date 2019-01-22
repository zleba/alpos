
#include "alpos/functions/AQcdnumDDISCS.h"

#include <iostream>

using namespace std;

const std::vector<std::string> AQcdnumDDISCS::fRequirements = {"Aq","Bq","Cq","Ag","Bg","Cg"}; //< List of all AParm's which this function depends on
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
   fValue.resize(GetRequirements().size());
   fError.resize(GetRequirements().size());

   vector<double> q2Vec = DOUBLE_COL_NS(Data,Q2,GetAlposName());
   cout << "q2Vector size " << q2Vec.size() << endl;

   for(auto q2 : q2Vec)
       cout << q2 <<" "<< endl;

   //vector<double> x = DOUBLE_COL_NS(Data,x,GetAlposName());
   //vector<double> y = DOUBLE_COL_NS(Data,y,GetAlposName());



   fValue[0] = PAR(Aq);
   fValue[1] = PAR(Bq);
   fValue[2] = PAR(Cq);
   //fValue[3] = ((AParmD*)TheoryHandler::Handler()->GetParameter(GetFunctionName()+string(".")+string("par3") ))->GetValue();
   fValue[3] = PAR(Ag);
   fValue[4] = PAR(Bg);
   fValue[5] = PAR(Cg);

   return true;
}


// ___________________________________________________________________________________________ //
