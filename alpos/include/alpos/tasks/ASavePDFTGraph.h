// DB. 08/2015
#ifndef Alpos_ASavePDFTGraph
#define Alpos_ASavePDFTGraph

/** 
 * class ASavePDFTGraph
 * class ASavePDFTGraphResult
 * 
 * Save TGraphs of PDFs
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
// #include "alpos/ASavePDFTGraphResult.h"


// ____________________________________________________________________________ //
class ASavePDFTGraphResult : public ATaskResult {
public:
   ASavePDFTGraphResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~ASavePDFTGraphResult() {};
protected:

};


// ____________________________________________________________________________ //
class ASavePDFTGraph : public ATask {
public:
   static const std::string fTaskType;

public:
   // ASavePDFTGraph();
   ASavePDFTGraph(const std::string& aname); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //ASavePDFTGraph(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~ASavePDFTGraph();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   //std::string fTaskType;
   std::map<std::string,std::vector<double> > GetPDFdef() const ;
   ASavePDFTGraphResult* fResult;

};


#endif
