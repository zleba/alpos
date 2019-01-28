// DB. 08/2015
#ifndef Alpos_ASaveDPDFTGraph
#define Alpos_ASaveDPDFTGraph

/** 
 * class ASaveDPDFTGraph
 * class ASaveDPDFTGraphResult
 * 
 * Save TGraphs of DPDFs
 * 
 */

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <fstream>

#include "alpos/AlposObject.h"
#include "fastnlotk/read_steer.h"

#include "alpos/ATask.h"
// #include "alpos/ASaveDPDFTGraphResult.h"


// ____________________________________________________________________________ //
class ASaveDPDFTGraphResult : public ATaskResult {
public:
   ASaveDPDFTGraphResult(const std::string& aname,const std::string& taskName) : ATaskResult(aname,taskName){};
   ~ASaveDPDFTGraphResult() {};
protected:

};


// ____________________________________________________________________________ //
class ASaveDPDFTGraph : public ATask {
public:
   static const std::string fTaskType;

public:
   // ASaveDPDFTGraph();
   ASaveDPDFTGraph(const std::string& aname); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   //ASaveDPDFTGraph(const std::string& aname, const std::string& rsnmsp); //!< rsnmsp: read_steer-namespace, i.e. where to look up the steering values
   virtual ~ASaveDPDFTGraph();

   virtual bool Init();
   virtual bool Execute();
   virtual ATaskResult* GetResult() {return fResult;}; //!< Get the results
   virtual const std::string& GetTaskType() const { return fTaskType;}; //!< Get task type, often similar to the class name

protected:

   std::vector<double> GetLogNodes(double min, double max, int npts) const;
   void WriteDPDFGrid(std::ostream& strm, const std::vector<int>& DPDFid, const std::vector<double>& xpt, const std::vector<double>& qpt, std::vector<double>* Grid0=NULL);
   void WriteAlphasGrid(std::ostream& strm, const std::vector<double>& qpt, std::vector<double>* Grid0=NULL);
   //std::string fTaskType;
   std::map<std::string,std::vector<double> > GetDPDFdef() const ;
   ASaveDPDFTGraphResult* fResult;

};


#endif
