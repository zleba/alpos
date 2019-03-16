#include "alpos/ATheoryHandler.h"
#include "fastnlotk/read_steer.h"

#include "alpos/AFactory.h"
#include "alpos/ASuperTheory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include "alpos/AlposTools.h"


/* 
   A temporary implementation of the theory handler

*/

using namespace std;

//____________________________________________________________________________________ //
TheoryHandler* TheoryHandler::instance = NULL;

// // // //  _______ ___     _______ _______ _______  
// // // // /   _   Y   )   Y   _   Y   _   Y   _   Y
// // // // |   l   |   |   |   l   |   |   |   l___l
// // // // |   _   |   |___|   ____|   |   |____   |
// // // // |   |   |   l   |   |   |   l   |   l   |
// // // // l___|   |_______l___|   l_______l_______/
// // // //     `---'                                  


//              ,/k.
//             /  ih,             
//          ,-' ,  `:7b           
//      _.-/   '  /b.`.4p,
//   --"  ,    ,-' ^6x, `."^=._


//____________________________________________________________________________________ //
TheoryHandler::TheoryHandler() : AlposObject("TheoryHandler") {
   cout<<endl;
   cout<<"          ,"<<endl;
   cout<<"      \\   |    /                                                                                      "<<endl;
   cout<<" `.    \\  '   /  .'         _______ ___     _______ _______ _______             , __;_,__(            "<<endl;
   cout<<"   `. .-*\"\"*-. .'          /   _   Y   )   /   _   Y   _   Y   _   \\         _/./  .    /__.-         "<<endl;
   cout<<"*-._ /        \\ _.-*\"      |   l   |   |   |   l   |   |   |   l___l        /'   .  <  _ (_           "<<endl;
   cout<<"     (        )            |   _   |   |___|   ____|   |   |____   |        L.__/__) _./              "<<endl; 
   cout<<".-*\" \\        / \"*-._      |   |   |   l   |   |   |   l   |   l   |        (-(- D                    "<<endl;  
   cout<<"   .' `-.__.-' `.          l___|   |_______l___|   l_______l_______/        /    |  ..-.              "<<endl;
   cout<<" .'   /   .  \\   `.            `---'                                         \\_.''       \\            "<<endl;  
   cout<<"     /    |   \\            An object-oriented data to theory               __|            \\           "<<endl;
   cout<<"          '                comparison and fitting framework              ,'   \\  -._.\\     \\          "<<endl;
   cout<<"                                                                        /      \\      \\     |         "<<endl;
   cout<<"                                                                       /    )   \\     .`.   )         "<<endl; 
   cout<<"                                                               .--.   /    /       --'   |  `         "<<endl;           
   cout<<"                                                               Oc_ `-'  _.' \\.   /       |  |         "<<endl;           
   cout<<"                                                                  `l _)/     `-..       .|  :         "<<endl;             
   cout<<"                                                                   `._j          |   `./ L  J>        "<<endl;           
   cout<<"   @                                               ,sdPBbs.                      L   _/  `--'         "<<endl;   
   cout<<"    @      .-'^--.                               ,d$$$$$$$$b.                   (    (                "<<endl;  
   cout<<"          (     , )               ,/'b.        ,d$P'`Y'`Y'`?$b,       .         \\   J               ,d"<<endl; 
   cout<<"   ___I_ (  ,      )            _/    7h,    ,/    / \\      |  \\,   ,d*\\         )  \\    ,/'\\.    ,d  "<<endl;
   cout<<"  /\\-_--\\ ( \\'./  .'   .d^*b, .'    , _ #b--'     ^  |      \\    |./  ' b,_     / ` )) _/   .h,_ /   \\"<<endl;
   cout<<" /  \\_-__\\ '-| |-'    /  , _\\/     _/'.   /'  ,      /   \\   ^    `-._    \\  ,-'/,  /   , _/ - '/' ^  "<<endl;
   cout<<" |[]| [] |   | |  __/'  /   ,, __/'    ,-'   ,   /    |    |     \\     `--/' C/-..,'   _/'.  ,-'   '  "<<endl;
   cout<<"-------------------------------------------------------------------------------------------------------"<<endl;
   cout<<endl;
   cout<<endl;
}


//____________________________________________________________________________________ //
TheoryHandler::~TheoryHandler(){

}


//____________________________________________________________________________________ //
TheoryHandler* TheoryHandler::Handler(){
   if (instance==NULL) instance = new TheoryHandler();
   return instance;

}


//____________________________________________________________________________________ //
bool TheoryHandler::InitTheory(const std::string& steerfile) {
   //! Init all theory functions and parameters.
   debug["InitTheory"]<<"("<<steerfile<<")"<<endl;

   if ( steerfile != "" ) {
      string actualsteerfile = steerfile;
      AlposTools::CheckFileExit(actualsteerfile);
      READ(actualsteerfile);
   }
   else {
      debug["InitTheory"]<<"No steering file given, but this is alright."<<endl;
   }


   if ( !fAllParmNames.empty() ) {
      error["InitTheory"]<<"InitTheory should only be called once! Exiting."<<endl;
      exit(1);
   }
   set<string> AllFuncNames;

   //!  fAllParNames contains
   //!   1)  All primitive values from steering
   //!   2)  All function names
   //!   3)  All default function parameters (like: CRunDec.AlphasMz)
   //!   4)  All function instance parameters (like: As1.AlphasMz, PDF1.Filename)

   //1)
   debug["InitTheory"]<<"  *** First fill fAllParmNames *** \n"<<endl;
   fAllParmNames = read_steer::Steering()->GetAvailableLabels();

   debug["InitTheory"]<<"  *** First initializing functions *** \n"<<endl;
   vector<vector<string> > lfunc = read_steer::getstringtable("InitFunctions");
   vector<AFuncD*> funcs;
   for ( const auto& line : lfunc ) {
      info["InitTheory"]<<"Init a function '"<<line[1]<<"', with name: "<<line[0]<<endl;
      // 2,3,4)
      funcs.push_back(AFactory::FunctionFactory(line[1],line[0]));
      if (!CheckFunction(funcs.back())) {
         error["InitTheory"]<<"Failed to initialize function '"<<line[0]<<"' of type '"<<line[1]<<"'."<<endl;
      }
      AllFuncNames.insert(funcs.back()->GetAlposName());
      // init DataTheoryPairs (quite dirty)
      if ( line[1]=="Data" ) 
         fDataTheoryPairs[line[0]] = make_pair((AData*)funcs.back(),GetFuncD(line[0]));
   }
   fAllParmNames.insert("");
   fAllParmNames.erase("");

   // a short crosscheck
   if (say::debug.GetSpeak()) PrintListOfAllParms();


   debug["InitTheory"]<<"  *** Registering constants *** \n"<<endl;
   set<string> pars = read_steer::Steering()->GetAvailableLabels();
   for ( const auto& lab : pars ) {
      if ( lab=="")  continue;
      string val = read_steer::Steering()->gets(lab);
      if ( fAllParmNames.count(val) != 0 ) {
         debug["InitTheory"]<<" + Skipping initialization of " <<lab<<endl;
         // register alias later
      }
      else {
         if (read_steer::CheckNumber(read_steer::Steering()->gets(lab)) ) {
            // numeric constant
            debug["InitTheory"]<<" + New (double) parameter with name: "<<lab<<", and value: "<<read_steer::Steering()->getd(lab)<<endl;
            new AParmD(lab,read_steer::Steering()->getd(lab),0,false);
         }
         else { // string constant
            debug["InitTheory"]<<"  + New (string) parameter with name '"<<lab<<"' and value '"<<val<<"'."<<endl;
            string err = "";
            new AParmS(lab,val,err,false);	    
         }
      }
   }


   //! Now init aliases
   //! Aliases are allowed to point to:
   //!    + functions (i.e. As1)
   //!    + Parameters with a value (i.e. asMz = 91.16)
   //!    + Aliase which are already initialized
   debug["InitTheory"]<<"  *** Now initializing aliases of constants *** \n"<<endl;
   for ( const auto& lab : pars ) {
      if ( lab=="") continue;
      string val = read_steer::Steering()->gets(lab);
      if ( fAllParmNames.count(val) != 0 ) {
         if (read_steer::CheckNumber(read_steer::Steering()->gets(lab)) ) {
            // nothing to do. It is a numeric value
         }
         else {
            if ( TheoryHandler::Handler()->CheckParameter(val) ) { //read_steer::Steering()->exist(val) ) {
               //if ( read_steer::Steering()->exist(val) || TheoryHandler::Handler()->CheckParameter(val) ) { //! allow forward declaration
               debug["InitTheory"]<<" + New ALIAS '"<<lab<<"' for '"<<val<<"'."<<endl;
               TheoryHandler::Handler()->NewAlias(lab,val);
            }
            else { 
               warn["InitTheory"]<<"Could not initialize parameter (neither with value, or as alias). ParName: "<<lab<<", value="<<val<<endl;
               warn["InitTheory"]<<"Aliase are only allowed to point to 1) parameters with a value, 2) initializations of functions, 3) to already otherwise defined aliase or parameters."<<endl;
            }
         }
      }
   }

   //! finalize initialization of functions
   //!   + Register dependencies
   //!   + Call 'Init'
   debug["InitTheory"]<<"  *** Now finalize initialization of functions *** \n"<<endl;
   //for ( auto line : lfunc ) 
   for ( const auto& line : AllFuncNames ) {
      debug["InitTheory"]<<">> Registering dependencies for function: "<<line<<endl;
      bool success = ARegisterRequirements(GetFuncD(line));
      debug["InitTheory"]<<"Call Init() function: "<<line<<endl;
      success += GetFuncD(line)->Init();
      if ( success ) 
         debug["InitTheory"]<<"Function "<<line <<" initialized successfully."<<endl;
      else 
         error["InitTheory"]<<"Function "<<line<<" could not be initialized successfully."<<endl;

   }

   info["InitTheory"]<<"Theory initialized."<<endl;

   return true;
}


//____________________________________________________________________________________ //
//bool TheoryHandler::InitSuperFunctions(const vector<std::string>& reqs) {
bool TheoryHandler::InitSuperFunctions() {
   //! Init SuperTheory and SuperData, which are functions that
   //! calculate one array/vector which holds all values of the specified functions

   // make a check here, if the 'DataTheoryPair'-map is reasonable
   for ( auto id : GetDataTheoryPairs()) {
      string dan = id.second.first->GetAlposName();
      dan.resize(dan.size()-5);
      if ( id.second.second->GetAlposName() != id.first ||  dan != id.first  ) {
         warn["InitSuperFunctions"]<<"The map 'DataTheoryPairs' has different map-key than the function name would suggest! Exiting."<<endl;
         warn["InitSuperFunctions"]<<"This may cause problems, since map.first is not the same as the data or theory-name."<<endl;
         cout<<"map-key: " <<id.first<<endl;
         cout<<"     da: "<<id.second.first->GetAlposName()<<endl;
         cout<<"     th: "<<id.second.second->GetAlposName()<<endl;
         exit(1);
      }
   }

   vector<string> reqs;
   for ( auto id : GetDataTheoryPairs()) {
      //reqs.push_back(id.first); // supposely the same
      reqs.push_back(id.second.second->GetAlposName()); // get the name of the theory function
   }
   debug["InitSuperFunctions"]<<"reqs.size()="<<reqs.size()<<endl;
   for ( auto i : reqs )
      debug["InitSuperFunctions"]<<"       "<<i<<endl;


   // Init SuperTheory
   ASuperTheory* SupTh =  (ASuperTheory*)AFactory::FunctionFactory("SuperTheory","SuperTheory");
   CheckFunction(SupTh);
   SupTh->SetRequirements(reqs);
   for ( const auto& req : reqs ) {
      string alias = "SuperTheory."+req;
      debug["InitSuperFunctions"]<<" + New ALIAS '"<<alias<<"' for '"<<req<<"'."<<endl;
      TheoryHandler::Handler()->NewAlias(alias,req);
   }
   ARegisterRequirements(SupTh);
   //SupTh->Init(); // call Init() (or not)!?

   // Init SuperData
   ASuperData* SupDa =  (ASuperData*)AFactory::FunctionFactory("SuperData","SuperData");
   CheckFunction(SupDa);
   for ( const auto& req : reqs ) {
      // todo: decide whether 'SuperData' or 'SuperData_Data' (or maybe even SuperTheory_Data') shall be reigstered
      // todo: decide whtere requirements of SuperData are named <DataSet> or <DataSet>_Data (while both have to point/alias ot the same).
      //((ASuperTheory*)SupDa)->AddRequirement(req+"_Data"); 
      //string alias = "SuperData."+req+"_Data";
      SupDa->AddRequirement(req); // use same names as in SuperTheory     
      string alias = "SuperData."+req;
      debug["InitSuperFunctions"]<<" + New ALIAS '"<<alias<<"' for '"<<(req+"_Data")<<"'."<<endl;
      TheoryHandler::Handler()->NewAlias(alias,req+"_Data");
   }
   ARegisterRequirements(SupDa);
   // SupDa->Init(); // call Init() (or not)!?

   fSuperPair = make_pair(SupDa,SupTh);

   return true;
}


//____________________________________________________________________________________ //
bool TheoryHandler::InitSubsetFunctions(const vector<std::string>& reqs, bool DoCuts, bool DoSubsets) {
   //! initialize the subset functions
   //! Subsets are specified in the steering cards
   //! A special subset is the subset, where the cuts are applied: DoCuts should always be 'true'
   for ( const auto& ireq : reqs ) {
      AData* da = (AData*)GetFuncD(ireq+"_Data");
      for ( const auto& isub : da->GetSubsets() ) {
         bool UseSubset = false;
         bool IsCut = ( isub.first.find("cuts") != string::npos );
         if ( DoCuts && IsCut ) UseSubset=true;
         if ( DoSubsets && !IsCut ) UseSubset=true;
         if ( UseSubset ) {
            // subsetdata
            string subsetname = ireq+"_"+isub.first;
            ASubsetData* SubDa     =  (ASubsetData* )AFactory::FunctionFactory("SubsetData",subsetname);
            SubDa->SetRequirementValidPoints(ireq,isub.second);  
            // subset theory
            ASubsetFunction* SubTh =  (ASubsetFunction* )AFactory::FunctionFactory("SubsetFunction",subsetname);
            SubTh->SetRequirementValidPoints(ireq,isub.second);
            // register subset
            if ( IsCut ) {
               // put no-cut data/theory pair into subsetpair
               // make a new subset without cuts (called <xyz>_NoCuts) 
               {
                  vector<bool> nocuts(isub.second.size(),true);
                  string subnamenocut = ireq+"_NoCuts";
                  ASubsetData*     SubDaNoCut =  (ASubsetData* )    AFactory::FunctionFactory("SubsetData",subnamenocut);
                  ASubsetFunction* SubThNoCut =  (ASubsetFunction* )AFactory::FunctionFactory("SubsetFunction",subnamenocut);
                  SubDaNoCut->SetRequirementValidPoints(ireq,nocuts);  
                  SubThNoCut->SetRequirementValidPoints(ireq,nocuts);
                  fSubsetPairs[ireq][subnamenocut] = make_pair(SubDaNoCut,SubThNoCut);
               }
               // replace 'DataTheoryPair' with subsets including cuts
               fDataTheoryPairs.erase(ireq);
               //fDataTheoryPairs[ireq] = make_pair(SubDa,SubTh);
               fDataTheoryPairs[subsetname] = make_pair(SubDa,SubTh); // give it the new proper name
            }
            else { // a simple subset 
               fSubsetPairs[ireq][subsetname] = make_pair(SubDa,SubTh);
            }


         }
      }
   }
   return true;
}


//____________________________________________________________________________________ //
void TheoryHandler::PrintListOfAllParms() const {
   //! Print the names of all available parameters
   cout<<"\n---------------------------------------\n"<<endl;
   cout<<"List of all parameters in this setup"<<endl;
   for ( const auto& in : fAllParmNames ) 
      cout<<"   "<<in<<endl;
   cout<<"---------------------------------------\n"<<endl;
}


//____________________________________________________________________________________ //
void TheoryHandler::PrintCurrentTheorySet(ATheorySet* setptr, std::ostream& strm) const{
   //! Print the current theory set (i.e. all parameters and its values)
   //! or if 'set' is specified, print this set


   static const string ten = "++++++++++";
   strm<<"\n  +"<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<endl;

   ATheorySet set = GetTheorySet();
   if (setptr==nullptr) setptr = &set;

   //for ( int zwo = 0 ; zwo<2 ;zwo++ ){
   for ( const auto& i : setptr->GetSet() ) {
      string name = i.first;
      string alias = TheoryHandler::Handler()->GetParameter(name)->GetAlposName();
      if ( alias == name ) 
         alias = "";
      else {
         alias = "--> "+alias;
         name += "  "+alias;
      }
      string cnst = i.second.Const ? "(const)" : "";
      if ( i.second.Const ) name += "  (const)";

      string vals = "[empty]";
      if ( !i.second.Values.empty() ) { // if no datasets are instantiated
         vals = i.second.Values[0]; // single values
         if ( !(i.second.Errors[0] == "0" ||  i.second.Errors[0] == "" )) vals += " +- "+i.second.Errors[0]; // errors
      }

      // array (functions)
      if ( i.second.Values.size() > 1 ) {
         vals += " [0],  "+i.second.Values[1];
         if ( !(i.second.Errors[1] == "0" ||  i.second.Errors[1] == "" )) vals += " +- "+i.second.Errors[1];
      }
      if ( i.second.Values.size() > 2 ) {
         vals += " [1],  "+i.second.Values[2];
         if ( !(i.second.Errors[2] == "0" ||  i.second.Errors[2] == "" )) vals += " +- "+i.second.Errors[2];
      }
      else if ( i.second.Values.size() > 1 ) vals += " [1]";

      if ( i.second.Values.size() == 4 ) {
         vals += " [2],  "+i.second.Values[2];
         if ( !(i.second.Errors[2] == "0" ||  i.second.Errors[2] == "" )) vals += " +- "+i.second.Errors[2];
         vals += " [3];";
      }
      else if (i.second.Values.size() > 4 ){
         vals += " [2],  ... , ";
         vals += i.second.Values.back();
         if ( !(i.second.Errors.back() == "0" ||  i.second.Errors.back() == "" )) vals += " +- "+i.second.Errors.back();
         vals += " ["+to_string(i.second.Values.size()-1)+"];";
      }

      //strm<<"  +  "<<i.first<<"  \t\t"<<i.second.Values[0]<<"\t const="<<i.second.Const<<endl;
      // 'ordered' or not? (alpha_s has for instance only 'one' value)
      // if ( zwo==0 && i.second.Values.size() == 1 )
      //    printf("  +  %-50s  %-7s  %s\n",
      // 	   name.c_str(), cnst.c_str(),vals.c_str());
      // else if (zwo==1 && i.second.Values.size() > 1)
      //    printf("  +  %-50s  %-7s  %s\n",
      // 	   name.c_str(), cnst.c_str(),vals.c_str());
      //}
      char buf[1000];
      sprintf(buf,"  +  %-55s  %-7s  %s\n",
            name.c_str(), cnst.c_str(),vals.c_str());
      strm<<buf;
}
strm<<"  +"<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<ten<<endl;
strm<<endl;
}


//____________________________________________________________________________________ //
void TheoryHandler::SaveCurrentTheorySet(const string& outname, ATheorySet* setptr) const{
   //! Print the current theory set (i.e. all parameters and its values)
   //! or if 'set' is specified, print this set

   // init
   string outputdirectory =  Alpos::Current()->Settings()->outputdir;
   ATheorySet set = GetTheorySet();
   if (setptr==nullptr) setptr = &set;

   // -------- write to ascii file in 'print' format
   {
      string printfile  = outputdirectory+outname+"_cout.txt";
      info["SaveCurrentTheorySet"]<<"Writing Theory set in 'cout'-format to ascii file : "<<printfile<<endl;
      ofstream pfile(printfile.c_str());
      PrintCurrentTheorySet(setptr,pfile);
   }

   // -------- write to ascii file
   {
      string outputfile  = outputdirectory+outname+".txt";
      info["SaveCurrentTheorySet"]<<"Writing Theory set to ascii file : "<<outputfile<<endl;
      //ofstream file((outputdirectory+GetTaskName()+"_"+dataChildren[iChild]->GetAlposName()+".txt").c_str());
      ofstream file(outputfile.c_str());
      for ( const auto& i : setptr->GetSet() ) {
         string name = i.first;
         // string alias = TheoryHandler::Handler()->GetParameter(name)->GetAlposName();
         // if ( alias == name ) 
         //    alias = "";
         // else {
         //    alias = "--> "+alias;
         //    name += "  "+alias;
         // }
         // string cnst = i.second.Const ? "(const)" : "";
         // if ( i.second.Const ) name += "  (const)";

         char buf [55];
         sprintf(buf,"%-60s",
               name.c_str(), "");

         file << buf;
         for ( auto val : i.second.Values ) file<<"\t"<<val;
         file<<endl;
      }
   }

   // -------- write to root file
   {
      TDirectory* rootfile  = Alpos::Current()->Settings()->rootoutput; 
      TDirectory* outtdir   = rootfile->mkdir(outname.c_str());
      info["SaveCurrentTheorySet"]<<"Writing Theory set to TDirectory "<<outtdir->GetName()<<" in root file :"<< rootfile->GetName() <<endl;
      outtdir->cd();


      // reset
      rootfile->cd();
   }

}


//____________________________________________________________________________________ //
ATheorySet TheoryHandler::GetTheorySet() const {
   //! Get all theory parameters and values
   ATheorySet set;
   for ( const auto& i : fParms ) {
      AThCont cont; 
      i.second->GetContent(cont);
      set.AddContent(i.first,cont);
   } 
   return set;
}


//____________________________________________________________________________________ //
void TheoryHandler::SetTheorySet(const ATheorySet& set) {
   //! Set a full set of theory parameters
   for ( const auto& i : set.GetSet() ) {
      GetParameter(i.first)->SetContent(i.second);
   }
}


//____________________________________________________________________________________ //
bool TheoryHandler::CheckFunction(AFuncD* func) {
   //! Check if function is resonable
   //! Add parameters to list of parameters

   if ( !func ) return false;
   fAllParmNames.insert(func->GetAlposName());
   for( const auto& ir : func->GetRequirements() ) {
      fAllParmNames.insert(func->GetAlposName()+"."+ir);
      fAllParmNames.insert(func->GetFunctionName()+"."+ir);
   }
   debug["CheckFunction"]<<"Todo. Check if 'TheoryHandler::CheckFunction() GetRequirements' provides same 'requirements' as func->GetRequirements()!"<<endl;

   return true;
}


//____________________________________________________________________________________ //
set<string> TheoryHandler::GetRequirements(const string& functype) {
   // get requirements:
   //   dirty: get 'default' values from 'default' namespace
   //   (since 'GetRequirements()' can only be called on an instance)
   //   it would be sufficient to access <AlposSet>::fRequirements                                        

   set<string> keys = read_steer::Steering()->GetAvailableLabels(); // look up in default NS
   debug["GetRequirements"]<<"functype: "<<functype<<endl;
   debug["GetRequirements"]<<"keys.size(): "<<keys.size()<<endl;

   set<string> ret;
   string prefix = functype + ".";
   for ( const auto& lbl : keys ) {
      // if (!lbl.compare(0, prefix.size(), prefix)) {// lbl begins with "'functype+'.'"
      // 	 string pname = atoi(lbl.substr(prefix.size()).c_str());
      // 	 cout<<"found requriement "<<pname <<" for functype: "<<functype<<endl;
      // 	 ret.insert(pname);
      // }
      debug["GetRequirements"]<<"Checking lbl: "<<lbl<<" for prefix: "<<prefix<<endl;
      if(lbl.substr(0, prefix.size()) == prefix) {
         std::string pname = lbl.substr(prefix.size());
         debug["GetRequirements"]<<"Found requriement "<<pname <<" for functype: "<<functype<<endl;
         ret.insert(pname);
      }
   }

   return ret;
}


//____________________________________________________________________________________ //
bool TheoryHandler::RegisterParameter(AParmNamed* par) {
   debug["RegisterParameter"]<<"AlposName="<<par->GetAlposName()<<std::endl;
   //! Register a new parameter
   if (CheckParameter(par->GetAlposName()) ) 
      error["RegisterParameter"]<<"A parameter with the same name ('"<<par->GetAlposName()<<"') is already registered."<<std::endl;
   else {
      debug["RegisterParameter"]<<"Registering new parameter with name: " <<par->GetAlposName()<<std::endl;
      fParms[par->GetAlposName()]=par;
      //par->Init();
   }
   return true;
}


//____________________________________________________________________________________ //
bool TheoryHandler::NewAlias(const std::string& alias,const std::string& aliasfor) {
   //! Register a new alias with name 'alias' for the parameter 'aliasfor'
   //! The parameter 'aliasfor' should already exist.
   //! Functions returns false if an error occurs.
   if ( !CheckParameter(aliasfor) ) {
      error["NewAlias"]<<"Parameter with name '"<<aliasfor<<"' does not exist."<<endl;
      return false;
   }
   debug["NewAlias"]<<" > New Alias: "<<alias<<" -> "<<aliasfor<<endl;
   fParms[alias] = fParms[aliasfor];
   return true;
}


/*
//____________________________________________________________________________________ //
set<string> TheoryHandler::GetListOfAllParm() const {
//! get a set of all (future) existing parameter names 
//!  These are:
//!   1)  All primitive values from steering
//!   2)  All function names
//!   3)  All default function parameters (like: CRunDec.AlphasMz)
//!   4)  All function instance parameters (like: As1.AlphasMz, PDF1.Filename)

// 1) All primitive values from steering
set<string> ret = read_steer::Steering()->GetAvailableLabels();

// 2) All function names
vector<vector<string> > fcn = read_steer::getstringtable("InitFunctions");
for ( const auto& line : fcn ) 
ret.insert(line[0]);

// 3) All default function parameters (like: CRunDec.AlphasMz)
// for ( const auto& line : fcn ) 
//    line[1]
// 	 hier

// 4) All function instance parameters (like: As1.AlphasMz, PDF1.Filename)

// return;
return ret;
}
*/

//____________________________________________________________________________________ //
bool TheoryHandler::CheckParameter(const std::string& name) const {
   if (fParms.count(name)>0)
      return true;
   else return false;
}



//____________________________________________________________________________________ //
bool TheoryHandler::EraseParameter(const std::string& name) {
   //! Delete a parameter
   if (CheckParameter(name)){
      fParms.erase(name);
      return true;
   }
   else {
      warn["EraseParameter"]<<"TheoryHandler::EraseParameter. Parameter '"<<name<<"' does not exist."<<std::endl;
      return false;
   }
}


//____________________________________________________________________________________ //
AParmNamed* TheoryHandler::GetParameter(const std::string& name) {
   //std::cout<<"TH:GetParameter(name="<<name<<")"<<std::endl;
   const auto& itf = fParms.find(name);
   if ( itf != fParms.end() ) {
      //debug["GetParameter"]<<"Returning fparms["<<name<<"]"<<std::endl;
      return fParms.find(name)->second;//[name]; 
   }
   else {
      string rname = FindParameter(name); 
      //cout<<"rname = "<<rname<<"\t (looking for: "<<name<<")"<<endl;
      if ( rname != "" ) {
         return fParms.find(rname)->second;  
      }
      else { 
         error["GetParameter"]<<"Parameter '"<<name<<"' does not exist. Exiting"<<std::endl;
         exit(1);// it turned out that in most cases it is not reasonable to continue the program.
         // 
         new AParmNamed(name);
         error["GetParameter"]<<"Adding new AParmNamed with name: "<<name<<std::endl;
         error["GetParameter"]<<"It will not be possible to call 'GetValue()' or 'SetValue()' on this parameter."<<std::endl;
         return fParms[name];
         // GetParameter(name)->SetIsOutdated();
      }
   }

   // if (CheckParameter(name)) {
   //    //std::cout<<"TH:GetParameter() returning fparms["<<name<<"]"<<std::endl;
   //    return fParms.find(name)->second;//[name];
   // }
   // else {
   //    string rname = FindParameter(name);
   //    if ( rname!="" ) {
   // 	 return fParms.find(rname)->second;
   //    }
   //    else {
   // 	 std::cout<<"TheoryHandler::GetParameter. Parameter '"<<name<<"' does not exist."<<std::endl;
   // 	 std::cout<<" +++ Adding new AParmNamed with name: "<<name<<std::endl;
   // 	 new AParmNamed(name);
   // 	 return fParms[name];
   // 	 // GetParameter(name)->SetIsOutdated();
   //    }
   // }
   return nullptr;
}


//____________________________________________________________________________________ //
string TheoryHandler::FindParameter(const std::string& name) const {
   //! Find a parameter's base name
   //! i.e. for name "fastNLO.Alphas.xp" may point to "As1.xp"
   char* myString = strdup(name.c_str());
   char* p = strtok(myString, ".");
   vector<string> toks;
   while (p) {
      //printf ("Token: %s\n", p);
      toks.push_back(p);
      p = strtok(NULL, ".");
   }
   free(myString);
   debug["FindParameter"]<<"("<<name<<")"<<endl;
   if ( toks.size()==1 ) { // nothing todo
      if (CheckParameter(name) ) 
         return fParms.find(name)->second->GetAlposName();
      else return "";
   }
   else if (toks.size() == 2 ) {
      string temp = toks[0]+"."+toks[1];
      if (CheckParameter(temp) ) 
         return fParms.find(name)->second->GetAlposName();
      else {
         warn["FindParameter"]<<"No parameter with name '"<<temp<<"' was found when looking for '"<<name<<"'."<<endl;
         if ( fParms.find(toks[0]) != fParms.end() ) {
            string baseFunction = TheoryHandler::Handler()->GetFuncD(toks[0])->GetFunctionName();
            warn["FindParameter"]<<"    Most likely '"<<toks[1]<<"' is not a 'member' of function '"<<baseFunction<<"'."<<endl;
         }
         else {
            warn["FindParameter"]<<"    No parameter "<<toks[0]<<" found at all."<<endl;
         }
         return "";
      }
   }
   else { // recursively
      string temp = toks[0]+"."+toks[1];
      string klar = FindParameter(temp);
      //cout<<" ** Recursively. temp="<<temp<<", klar="<<klar<<endl;
      if ( klar=="" ) return ""; // nothing found
      if ( klar==temp) {
         warn["FindParameter"]<<"Parameter not defined or you are looking for an ill-defined name (name="<<name<<"). '"<<temp<<"' may probably be a parameter but not a function."<<endl;
         return "";
      }
      if ( strstr(klar.c_str(), ".") != NULL ) {
         info["FindParameter"]<<"Parameter may not be further resolvable (name="<<name<<"). Found: '"<<temp<<" -> "<<klar<<"', which may probably be a parameter but not a function."<<endl;
         //return "";
      }
      klar = fParms.find(klar)->second->GetAlposName()+"."+toks[2];
      return FindParameter(klar);
   }

   //return "bla";
}


//____________________________________________________________________________________ //
bool TheoryHandler::ResolveName(std::string& name) const {
   //! Resolve the parameter's name into its base name
   //! Return false, if parameter was not found
   name = FindParameter(name);
   return (name != "");
}


//____________________________________________________________________________________ //
double TheoryHandler::GetParmDValue(const std::string& name) {
   //! Get the value of a parameter of type (double) 
   return ((AParmD*)GetParameter(name))->GetValue();
}


//____________________________________________________________________________________ //
double TheoryHandler::GetParmDError(const std::string& name){
   //! Get the error of a parameter of type (double) 
   return ((AParmD*)GetParameter(name))->GetError();
}


//____________________________________________________________________________________ //
const std::string& TheoryHandler::GetParmSValue(const std::string& name){
   //! Get the value of a parameter of type (string) 
   return ((AParmS*)GetParameter(name))->GetValue();
}


//____________________________________________________________________________________ //
const std::string& TheoryHandler::GetParmSError(const std::string& name){
   //! Get the error of a parameter of type (string) 
   return ((AParmS*)GetParameter(name))->GetError();
}

//____________________________________________________________________________________ //







// template <typename T>
// bool AParmFuncBase<T>::CheckUpdateRequested(const std::string& aname) const{
//    //! check whether this particluar parameter has changed
//    //! returns 'true' if updated is needed!

//    const set<AParmNamed*> nots = GetNotifiers();
//    string lname = GetAlposName() + "." + aname;
//    cout<<"TH::CheckUpdate: seeking:" <<aname<<", which is "<<lname<<"\t todo: not yet done."<<endl;
//    if ( nots.find(lname) != nots.end() ) {
//       cout<<"                 return true."<<endl;
//       return true;
//    }
//       cout<<"                 return false ."<<endl;
//    return false;
// }


//____________________________________________________________________________________ //


