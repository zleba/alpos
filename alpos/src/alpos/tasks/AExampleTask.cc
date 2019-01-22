#include "alpos/tasks/AExampleTask.h"
#include <iostream>

/* 
 ATask

 */

using namespace std;

const string AExampleTask::fTaskType = "AExampleTask";

//____________________________________________________________________________________ //
//AExampleTask::AExampleTask(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
AExampleTask::AExampleTask(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("AExampleTask");
   //! Important: create always new result-object here!
   fResult = new AExampleTaskResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
AExampleTask::~AExampleTask(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool AExampleTask::Init(){
   info<<"Hello. AExampleTask::Init()."<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool AExampleTask::Execute(){
   debug["Execute"]<<"Now getting 'WelcomeString' from steering and printing it:"<<endl;
   cout<<endl;
   cout<<"  "<<STRING_NS(WelcomeString,NS())<<endl;
   cout<<endl;
   return true;
}


//____________________________________________________________________________________ //
