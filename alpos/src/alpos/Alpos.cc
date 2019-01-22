#include <ctime>
#include <cstdlib>
#include "alpos/Alpos.h"
#include "alpos/ATheory.h"
#include "alpos/AFactory.h"
#include "alpos/AlposTools.h"
#include "TSystem.h"

/* 
 Alpos
 DB. Jan 2015.

 */

using namespace std;

Alpos* Alpos::fCurrent = NULL;

//____________________________________________________________________________________ //
Alpos::Alpos() : AlposObject() {

}


//____________________________________________________________________________________ //
Alpos::Alpos(const string& steerfile, int argc, char** argv) : AlposObject("Alpos"){
   fCurrent = this;
   read_steer::SetDefaultNamespace("AlposTheory");
   PARSE(argc,argv);
   ReadSteering(steerfile);
   AddTheorySteeringFlagsToDefault();
   TheoryHandler::Handler()->InitTheory();
   AddThFactorsToTheoryFunctions();
   TheoryHandler::Handler()->InitSubsetFunctions(GetDataTheorySets(),true,BOOL_NS(InitSubsets,steerfile));
   //TheoryHandler::Handler()->InitSuperFunctions(GetDataTheorySets());
   TheoryHandler::Handler()->InitSuperFunctions();
}


//____________________________________________________________________________________ //
Alpos::~Alpos(){
}



//____________________________________________________________________________________ //
bool Alpos::ReadSteering(const string& steerfile){
   fSteerfile = steerfile;
   READ_NS(steerfile,steerfile); 
   InitSettings();
   SetVerbosity(STRING_NS(GlobalVerbosity,fSteerfile));
   fDaThSets             = STRING_COL_NS(DataTheorySets,AlposName,fSteerfile);
   vector<string> files  = STRING_COL_NS(DataTheorySets,SteerFile,fSteerfile);
   vector<string> thfunc = STRING_COL_NS(DataTheorySets,TheoryFunction,fSteerfile);
   vector<vector<string>> dathtab = STRING_TAB_NS(DataTheorySets,fSteerfile);
   // --- ignore data sets with 'ignore' statements
   {
      vector<string> tmp_DaThSets, tmp_files, tmp_thfunc;
      vector<vector<string>> tmp_dathtab;
      set<int> ignore;
      for ( unsigned int i=0 ; i<fDaThSets.size() ;i++ ) {
	 for ( string value : dathtab[i] ) {
	    if ( value.find("gnore")!= string::npos ) {
	       info["ReadSteering"]<<"Found 'ignore' statement ("<<value<<") for data set "<<fDaThSets[i]<<". Ignoring it."<<endl;
	       ignore.insert(i);
	       break;
	    }
	 }
	 if ( ignore.count(i) == 0 ) {
	    tmp_DaThSets.push_back(fDaThSets[i]);
	    tmp_files.push_back(files[i]);
	    tmp_thfunc.push_back(thfunc[i]);
	    tmp_dathtab.push_back(dathtab[i]);
	 }
      }
      fDaThSets = tmp_DaThSets;
      files = tmp_files;
      thfunc = tmp_thfunc;
      dathtab = tmp_dathtab;
   }
   // ---
   if ( thfunc.size()!=files.size() ){
      error["ReadSteering"]<<"Column 'TheoryFunction' has different number of entries than colummn 'SteerFile' in table 'DataTheorySets'. Exiting."<<endl;
      exit(1);
   }
   // ---
   for ( unsigned int i=0 ; i<fDaThSets.size() ;i++ ) {
      if ( fDaThSets[i].find(".") != string::npos ) {
	 error["ReadSteering"]<<"Dataset names are not allowed to contain character '.'; Please change: "<<fDaThSets[i]<<endl;
	 exit(277);
      }
      ADD_NS("TheoryFunction",thfunc[i],fDaThSets[i]); // add 'TheoryFunction'
      if ( dathtab[i].size() > 3 ) { // add cuts and other steering parameters
	 vector<string> vrest(dathtab[i].begin() + 3, dathtab[i].end());
	 vector<string> cuts;
	 vector<string> param;
	 for ( auto ispec : vrest ) {
	    if ( ispec.find("==") != string::npos || 
		 ispec.find(">") != string::npos ||
		 ispec.find(">=") != string::npos ||
		 ispec.find("<") != string::npos ||
		 ispec.find("<=") != string::npos ||
		 ispec.find("!=") != string::npos ||
		 ispec.find(".gt.") != string::npos ||
		 ispec.find(".eq.") != string::npos ||
		 ispec.find(".lt.") != string::npos ||
		 ispec.find(".ge.") != string::npos ||
		 ispec.find(".le.") != string::npos ) { // cut
	       cuts.push_back(ispec);
	    }
	    else if (ispec.find("=") != string::npos ) { // steering parameter
	       param.push_back(fDaThSets[i]+"::"+ispec);
	    }
	    else {
	       warn["ReadSteering"]<<"Found additional parameter '"<<ispec<<"' for dataset '"<<fDaThSets[i]<<"', but it's not a cut or a steering parameter! Ignoring it."<<endl;
	    }
	 }
	 // steering parameters
	 PARSEV(param);
	 // cuts
	 ADDARRAY_NS("CutsMainSteering",cuts,fDaThSets[i]);
      }
   }
   for ( unsigned int i=0 ; i<fDaThSets.size() ;i++ ) {
      info["ReadSteering"]<<"Reading file '"<<files[i]<<"'."<<endl;
      debug["ReadSteering"]<<"Reading file into namespace:" <<fDaThSets[i]<<endl;
      AlposTools::CheckFileExit(files[i]);
      READ_NS(files[i],fDaThSets[i]);
   }
   return true;
}


//____________________________________________________________________________________ //
void Alpos::InitSettings(){
   //! Initialize fSettings from input steering.
   //! (ReadSteering() must have been called earlier).

   // --- error averaging
   if ( EXIST_NS(ErrorSymmetrization,fSteerfile) ) {
      static const std::map<string,AError::UpDnAveraging > UpDnMap {
	 { "Linear", AError::kLinear },
	 { "AbsLinear", AError::kAbsLinear },
	 { "AbsSignLinear", AError::kAbsSignLinear },
	 { "ImprovedLinear", AError::kImprovedLinear },
	 { "SignImprovedLinear", AError::kSignImprovedLinear },
	 { "Quadratic", AError::kQuadratic },
	 { "ImprovedQuadratic", AError::kImprovedQuadratic },
	 { "SignImprovedQuadratic", AError::kSignImprovedQuadratic },
	 { "Max", AError::kMax },
	 { "SignMax", AError::kSignMax },
	 { "Min", AError::kMin },
	 { "SignMin", AError::kSignMin },
	    };

      const string& avg = STRING_NS(ErrorSymmetrization,fSteerfile);
      if ( UpDnMap.count(avg)==0 ) {
	 info["InitSettings"]<<"Unrecognized value for 'ErrorSymmetrization' of '"<<avg<<"'. Using default: 'SignImprovedQuadratic'."<<endl;
         fSettings.ErrorSymmetrization = AError::kSignImprovedQuadratic;
      }
      else
	 fSettings.ErrorSymmetrization = UpDnMap.at(avg);
   }
   else {
      info["InitSettings"]<<"Steering flag 'ErrorSymmetrization' not found. Using default: 'SignImprovedQuadratic'."<<endl;
      fSettings.ErrorSymmetrization = AError::kSignImprovedQuadratic;
   }

   // --- numeric parameters
   // matrix inversion tolerance
   fSettings.matInvTol = EXIST_NS(MatrixInversionTolerance, fSteerfile) ? 
      DOUBLE_NS(MatrixInversionTolerance, fSteerfile) : 1e-6;
      
   // --- outputs
   fSettings.outputdir = "";
   fSettings.rootoutput = NULL;
   if ( EXIST_NS(Output,fSteerfile) && STRING_NS(Output,fSteerfile)!= "" ) {
      string out = STRING_NS(Output,fSteerfile);
      //if (out=="") out="./";
      string rootout;
      //string out2 = out.substr(0,out.rfind("/"));
      string last  = out.substr(out.rfind("/")+1,out.size());
      if ( last.find(".")!=string::npos && out!="." && out !="./" && out !="./.") {
	 fSettings.outputdir = out.substr(0,out.rfind("/"));
	 fSettings.outputdir += "/";
	 rootout = out.substr(0,out.rfind("."))+".root";
      }
      else {
	 if ( out.substr(out.size()-1,1) != "/" ) out+="/";
	 fSettings.outputdir = out;
	 rootout = fSettings.outputdir + "/alpos_results.root";
      }
      gSystem->mkdir(fSettings.outputdir.c_str(),true);
      fSettings.rootoutput = TFile::Open(rootout.c_str(),"RECREATE");
      if ( !fSettings.rootoutput ) {
	 error["InitSettings"]<<"root output file could not be created at: "<<rootout<<endl;
	 exit(16);
      }
      info["InitSettings"]<<"Output directory set to: "<<fSettings.outputdir<<endl;
      info["InitSettings"]<<"Output root-file set to: "<<fSettings.rootoutput->GetName()<<endl;
      // --- copy steering file
      if ( fSettings.outputdir != "./" ) {
	 string ex = "/bin/cp "+fSteerfile+" " +fSettings.outputdir+"steering.str";
	 gSystem->Exec(ex.c_str()); // or use glob instead
	 ex = "/bin/cp "+fSteerfile+" " +fSettings.outputdir+".";
	 gSystem->Exec(ex.c_str());
	 info["InitSettings"]<<"Input steering file copied to: "<<fSettings.outputdir<<endl;
      }
   }
   if ( fSettings.rootoutput==NULL ) {
      info["InitSettings"]<<"No output file will be written."<<endl;
      // it may avaoid crashes if a 'dummy' TDirectory is generated
      //fSettings.rootoutput = new TDirectory("AlposNoOutputfile","AlposNoOutputfile");
   }

   // --- ignore theory uncertainties
   fSettings.IgnoreTheoryErrors = false;
   if ( EXIST_NS(IgnoreTheoryErrors,fSteerfile) ) {
      fSettings.IgnoreTheoryErrors=BOOL_NS(IgnoreTheoryErrors,fSteerfile);
      info["InitSettings"]<<"Ignore theory uncertainties (flag found in steering): "<<(fSettings.IgnoreTheoryErrors ? "true" :  "false")<<endl;
   }

   // --- ignore theory uncertainties
   fSettings.Alpos_dir = getenv("PWD");
   //if ( fSettings.Alpos_dir=="" ) fSettings.Alpos_dir = system("pwd");
   // boost::filesystem::path full_path( boost::filesystem::current_path() );
   char* env_ad = getenv("ALPOS_DIR");
   if ( env_ad ) fSettings.Alpos_dir = env_ad;
   if ( EXIST_NS(ALPOS_DIR,fSteerfile) ) fSettings.Alpos_dir = STRING_NS(ALPOS_DIR,fSteerfile);
   if ( fSettings.Alpos_dir.back() != '/' ) fSettings.Alpos_dir+="/";
   info["InitSettings"]<<"ALPOS_DIR will be: "<<fSettings.Alpos_dir<<endl;


}


//____________________________________________________________________________________ //
void Alpos::SetVerbosity(const string& sverb){
   if (sverb=="DEBUG" || sverb=="Debug" || sverb=="debug")
      speaker::SetGlobalVerbosity(say::DEBUG);
   else if (sverb=="MANUAL" || sverb=="Manual" || sverb=="manual")
      speaker::SetGlobalVerbosity(say::MANUAL);
   else if (sverb=="INFO" || sverb=="Info" || sverb=="info")
      speaker::SetGlobalVerbosity(say::INFO);
   else if (sverb=="WARNING" || sverb=="Warning" || sverb=="warning")
      speaker::SetGlobalVerbosity(say::WARNING);
   else if (sverb=="ERROR" || sverb=="Error" || sverb=="error")
      speaker::SetGlobalVerbosity(say::ERROR);
   else if (sverb=="SILENT" || sverb=="Silent" || sverb=="silent")
      speaker::SetGlobalVerbosity(say::SILENT);
   else
      speaker::SetGlobalVerbosity(say::INFO);
}


//____________________________________________________________________________________ //
bool Alpos::AddTheorySteeringFlagsToDefault(){
   //!  + init data objects
   //!  + init theory functions
   //! as specified in the steering file(s)

   // --- add theory functions to list
   for ( const auto& is : fDaThSets ) {
      string functype = STRING_NS(TheoryFunction,is);
      debug["AddTheorySteeringFlagsToDefault"]<<"Init theory of type '"<<functype<<"' with name "<<is<<endl;
      ExtractTheoryFuncParms(is);
      //todo: AddFunctionTo read_steer InitFunctions { is functype }: 
      vector<string> entry = {is,functype};
      read_steer::appendtotable("InitFunctions",entry);
   }

   // --- add data functions to list
   for ( const auto& is : fDaThSets ) {
      debug["AddTheorySteeringFlagsToDefault"]<<"Init data with name "<<is<<endl;
      //vector<string> entry = {is+"_Data","Data"};
      vector<string> entry = {is,"Data"};
      read_steer::appendtotable("InitFunctions",entry);
   }
  
   return true;
}


//____________________________________________________________________________________ //
void Alpos::ExtractTheoryFuncParms(const string& AlposSet) {
   //! provide function parameters, which are defined in a DataTheory file
   //! in the format such that these are correctly  recognized by TheoryHandler::Handler()->InitFunction()
   // 1. get 'requirements' first
   // 2. look which parameters are defined
   // .. and provide these in the common namespace
   
   //string functype = STRING_NS(TheoryFunction,AlposSet);
   string functype = STRING_NS(TheoryFunction,AlposSet);
   // get requirements:
   set<string> reqs = TheoryHandler::Handler()->GetRequirements(functype);
   for ( const auto& req : reqs ) {
      // EXIST_NS(req,AlposSet) not working since called by name EXIST_NS(req,AlposSet) 
      if ( read_steer::getexist(req,AlposSet) ) { //if parm exist, add 'correctly' to default namespace
	 debug["ExtractTheoryFuncParms"]<<"Adding parameter to 'default' steering namespace! req="<<req<<endl;
	 ADD(AlposSet+"."+req,read_steer::getstring(req,AlposSet)); // ADD_NS, STRING_NS
      }
   }
}


//____________________________________________________________________________________ //
void Alpos::AddThFactorsToTheoryFunctions() {
   //! Add the list of k-factors as specified
   //! in the data-theory steering files
   //! to the theory functions

   for ( const auto& is : fDaThSets ) {
      string functype = STRING_NS(TheoryFunction,is);
      debug["AddThFactorsToTheoryFunctions"]<<"Init theory of type "<<functype<<" with name "<<is<<endl;
      vector<string> kNames = STRING_ARR_NS(TheoryFactors,is);
      for ( const auto& kn : kNames ){
	 if ( read_steer::CheckNumber(kn) ) {  //! ThFactor is specified directly  
	    debug["AddThFactorsToTheoryFunctions"]<<"Applying theory factor for '"<<is<<"'to all datapoints of: "<<stod(kn)<<endl;
	    TheoryHandler::Handler()->GetFuncD(is)->AddThFactor(kn,stod(kn));
	 }
	 else { // ThFactor is specified in column
	    vector<double>  k = read_steer::getdoublecolumn("Data",kn,is);
	    if ( k.empty() ){
	       warn["AddThFactorsToTheoryFunctions"]<<"Could not identify 'Data' column with name "<<kn<<" to get multiplicative theory-factors. Exiting."<<endl;
	       exit(1);
	    }
	    TheoryHandler::Handler()->GetFuncD(is)->AddThFactor(kn,k);
	 }
      }
   }
}


//____________________________________________________________________________________ //
void Alpos::DoTasks() {
   //! Run all tasks as specified in the steering
   fCurrent = this;
   
   // Add current theory set to the result
   ATaskResult* initalTS = new ATaskResult("VeryInitialAndFinalTheorySets","Alpos");
   initalTS->SetInitialTheorySet(  TheoryHandler::Handler()->GetTheorySet() );
   fResults["Alpos"] = initalTS;

   // run tasks
   vector<vector<string> > rstasks = STRING_TAB_NS(Tasks,fSteerfile);
   vector<string> tns  = STRING_COL_NS(Tasks,TaskName,fSteerfile);
   vector<string> tts  = STRING_COL_NS(Tasks,TaskType,fSteerfile);
   //   for ( const auto& ts : rstasks ) {
   int nsuctasks = 0;
   for ( unsigned int it = 0 ; it< rstasks.size() ; it++ ) {
      const string& taskname  = tns[it];
      const string& tasktype  = tts[it];

      info>>"\n"<<endl;
      info>>" +-----------------------------------------------------------------------------------------------------------+"<<endl;
      info>>" + ---------------------------- Alpos::DoTasks(). Instantiating new task. "<<endl;
      info>>" + ----------------------------        Task type: "<< tasktype <<endl;
      info>>" + ----------------------------        Task name: "<< taskname <<endl;
      info>>" +-----------------------------------------------------------------------------------------------------------+"<<endl;
      info>>"\n"<<endl;
      
      const time_t start = time(0);
      ATask* task = AFactory::TaskFactory(taskname,tasktype);
      if ( !task ) {
	 // this is not a task.
	 bool parexists = TheoryHandler::Handler()->CheckParameter(tasktype);
	 if ( parexists ) {
	    info["DoTasks"]<<"Task type was identified to be a parameter."<<endl;
	    info["DoTasks"]<<"Now setting parameter '"<<tasktype<<"' to new value: '"<<taskname<<"' (be careful with strings and numerical values.)"<<endl;
	    if ( read_steer::CheckNumber(taskname) ) 
	       TheoryHandler::Handler()->GetParmD(tasktype)->SetValue(stod(taskname),0,false);
	    else 
	       TheoryHandler::Handler()->GetParmS(tasktype)->SetValue(taskname,0,false);
	    info["DoTasks"]<<"Parameter set successfully. Now proceeding to next task."<<endl;
	    continue;
	 }
	 else {
	    warn["DoTasks"]<<"Cannot identify task of type (or parameter named) '"<<tasktype<<"' (name/value='"<<taskname<<"'). Continueing."<<endl;
	    continue;
	 }
      }

      info["DoTasks"]<<"Task '"<<taskname<<"' instantiated."<<endl;
      info["DoTasks"]<<"Running the task '"<<taskname<<"' of type '"<<tasktype<<"'."<<endl;
      //
      fResults[taskname] = task->GetResult();
      debug["DoTasks"]<<"Task '"<<taskname<<"' returned results: "<<fResults[taskname]<<endl;
      bool success = true;
      task->GetResult()->SetInitialTheorySet(  TheoryHandler::Handler()->GetTheorySet() );
      debug["DoTasks"]<<"TaskResult '"<<taskname<<"' returned results: "<<fResults[taskname]<<endl;
      success &= task->Init();
      debug["DoTasks"]<<"task->Init() called for task '"<<taskname<<"'. Return value: "<<success<<endl;
      success &= task->Execute();
      debug["DoTasks"]<<"task->Execute() called for task '"<<taskname<<"'. Return value: "<<success<<endl;
      task->GetResult()->SetFinalTheorySet( TheoryHandler::Handler()->GetTheorySet() );

      // -- timing
      int ss = difftime(time(0),start);
      int mm = ss/60;
      int hh = ss/60/60;
      mm-=hh*60;
      ss-=mm*60;
      char str[20];
      sprintf(str,"%02d:%02d:%02d",hh,mm,ss);

      // -- info message
      cout<<endl;
      if ( success ) info["DoTasks"]<<"Task '"<<taskname<< "' done successfully. Runtime "<<str<<endl;
      else  warn["DoTasks"]<<"Task '"<<taskname<< "'  FAILED. Runtime "<<str<<endl;
      if ( success ) nsuctasks++;

      debug["DoTasks"]<<"Now deleteing '"<<tasktype<<"' object."<<endl; // bad destructors often crash
      delete task;
   }
   // add also the very final theory set 
   initalTS->SetFinalTheorySet(  TheoryHandler::Handler()->GetTheorySet() );
   info>>"\n-----------------------------------------------------------------------------------------------------"<<endl;
   info["DoTasks"]<<"All tasks done ["<<nsuctasks<<"/"<<rstasks.size()<<"].\n"<<endl;
}


//____________________________________________________________________________________ //
