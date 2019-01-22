//! DB. Jan 2015
#include <iostream>
#include "alpos/Alpos.h"
#include <ctime>

int main(int argc, char** argv) {
   using namespace std;
   
   // ---  Parse commmand line                                                                            
   if (argc <= 1) {
      cout<<" Program: alpos"<<endl;
      cout<<"   usage:"<<endl;
      cout<<"   alpos <Steering.str> [steering parameters[...]]"<<endl;
      return 0;
   }
   /**
      //     Parsing steering values through the command line:
      //     Specify a value when executing alpos over the command line like:
      //            >$ alpos <steering.str> WikiLink::WelcomeString="www.google.com" steerfile=afilecontainingfurthervalues.str->WikiLink
      //                                                                                                                                                    
      //     This example will overwrite the steering parameter
      //     WikiLink {{{
      //            WelcomeString "www.google.com"
      //     }}}
      //     Further it will read in another steering file, with further values
      //     which may be used by the WikiLink task
    */
   

   

   // --- initialize and execute Alpos
   ios::sync_with_stdio(); // to not mix FORTRAN and C++ output
   time_t start = time(0);
   string steerfile = argv[1];
   Alpos alpos(steerfile,argc-1,&argv[1]);
   //TheoryHandler::Handler()->PrintCurrentTheorySet();
   alpos.DoTasks();

   // double sss = difftime(time(0),start);
   // int mmm = sss/60;
   // int hhh = sss/60/60;
   // mmm-=hhh*60;
   // sss-=mmm*60;
   // printf("Alpos total runtime %02d:%02d:%02.2f\n",hhh,mmm,sss);
   
   // print elapsed time
   clock_t tc = clock();
   float ss = ((float)tc)/CLOCKS_PER_SEC;
   int mm = ss/60;
   int hh = ss/60/60;
   mm-=hh*60;
   ss-=mm*60;
   printf("Alpos total runtime %02d:%02d:%02.2f\n",hh,mm,ss);

   return 0;
}
