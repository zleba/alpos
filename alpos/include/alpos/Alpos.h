// DB. 01/2015
#ifndef Alpos_Alpos
#define Alpos_Alpos

/** 
 * 
 * Alpos
 * The guy, who does the work!
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "alpos/AError.h"
#include "fastnlotk/read_steer.h"
#include "TFile.h"

class ATask;
class ATaskResult;


struct ASettings {
   AError::UpDnAveraging ErrorSymmetrization = AError::kSignImprovedQuadratic; //!< Error averaging
   std::string outputdir; //!< Output director for ascii output
   TDirectory* rootoutput = NULL; //!< Output root-file
   bool IgnoreTheoryErrors = false; //!< Ignore all theory errors
   std::string Alpos_dir = "";
   double matInvTol = 1e-6;
};

class Alpos : public AlposObject {

public:
   Alpos();
   Alpos(const std::string& steerfile, int argc=0, char** argv=NULL);
   ~Alpos();
   static Alpos* Current() { return fCurrent; } //!< Get 'current' Alpos instance

   bool ReadSteering(const std::string& steerfile);
   void InitSettings();
   bool AddTheorySteeringFlagsToDefault();
   void AddThFactorsToTheoryFunctions();

   const std::string& GetNamespace() const {return fSteerfile;} //!< Get namespace (here identical to the filename) where the steering parameters are stored
   const std::vector<std::string>& GetDataTheorySets() const { return fDaThSets;} //!< Get DataTheorySets
   void DoTasks(); //!< Run the tasks as specified in the steering
   void SetVerbosity(const std::string& verb);
   ASettings* Settings() { return &fSettings;} //!< Get Settings

private:
   static Alpos* fCurrent;

   void ExtractTheoryFuncParms(const std::string& AlposSet);

   std::string fSteerfile;
   std::vector<std::string> fDaThSets;
   std::map<std::string,ATaskResult*> fResults;
   ASettings fSettings;
};

#endif
