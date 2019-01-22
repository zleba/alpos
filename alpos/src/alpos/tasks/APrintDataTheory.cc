#include "alpos/tasks/APrintDataTheory.h"
#include "alpos/AFactory.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/AData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include "TH1D.h"

/** 
 APrintDataTheory

 Print data, total uncertainties and theory predictions in a tabular format.
 This is provided mainly to be of use when plotting results.

 */

using namespace std;

const string APrintDataTheory::fTaskType = "PrintDataTheory";

//____________________________________________________________________________________ //
//APrintDataTheory::APrintDataTheory(const string& aname, const string& rsnmsp/*, const std::map<std::string,ATaskResult> const *previousResults*/) : ATask(aname,rsnmsp/*,previousResults*/) {
APrintDataTheory::APrintDataTheory(const string& aname ) : ATask(aname) {
   //! constructor
   //! Important: create always new result-object here!
   fResult = new APrintDataTheoryResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
APrintDataTheory::~APrintDataTheory() {
   //! destructor.
   //! do not delete the AResult object
}


//____________________________________________________________________________________ //
bool APrintDataTheory::Init() {
   // --- set up task

   // read parameters from steering

   if ( EXIST_NS(Columns, NS()) ) {
      fColumnNames = STRING_ARR_NS(Columns,NS());
      info["Init"]<<"Requested colums: "<<endl;
      for (std::string s : fColumnNames ) {
         info["Init"]<<"    "<<s<<endl;
      }
   }
   else {
      info["Init"]<<"No columns specified, choosing default: 'iData', 'Data', 'Theory', 'ErrTot'"<<endl;;
      fColumnNames.push_back("iData");
      fColumnNames.push_back("Data");
      fColumnNames.push_back("Theory");
      fColumnNames.push_back("ErrTot");
   }

   if ( EXIST_NS(ColumnWidth, NS()) ) {
      fColumnWidth = INT_NS(ColumnWidth,NS());
      info["Init"]<<"Read column width: "<<fColumnWidth<<endl;
   }
   else {
      info["Init"]<<"No column width specified, choosing default: 12"<<endl;;
      fColumnWidth = 12;
   }


   return true;
}

//____________________________________________________________________________________ //
bool APrintDataTheory::Execute() {

   cout<<endl;
   cout<<" +--------------------  DataTheory  ---------------------- "<<endl;

   // pointers to superdata and supertheory
   const auto& supdata = TheoryHandler::Handler()->GetSuperPair().first;
   const auto& suptheo = TheoryHandler::Handler()->GetSuperPair().second;

   // get pointers to DataTables
   std::vector<AData*> dataChildren = supdata->GetChildren();           //TODO: handle subsets correctly
   std::vector<AParmFuncBase<double>*> theoChildren = suptheo->GetChildren();   //TODO: handle subsets correctly

   // write to file
   TDirectory* rootfile  = Alpos::Current()->Settings()->rootoutput;
   TDirectoryFile *dir   = NULL;
   if ( rootfile ) {
      rootfile->cd();
       dir = new TDirectoryFile(GetTaskName().c_str(),GetTaskName().c_str());
   }

   // go through each DataTable in SuperData (one per AData object)
   int rowsAddedSoFar = 0;
   for (int iChild = 0 ; iChild < dataChildren.size() ; iChild++ ) {
      // current data table
      const std::map<std::string, vector<double> > dataTable = dataChildren[iChild]->GetDataTable();

      // an irange {0..nDatapoints} (index of point in current dataset)
      std::vector<double> ptsIndices(dataChildren[iChild]->N());
      std::iota(ptsIndices.begin(), ptsIndices.end(), 0);

      // an irange {0..nDatapoints} (index of point in SuperData)
      std::vector<double> ptsSuperIndices(dataChildren[iChild]->N());
      std::iota(ptsSuperIndices.begin(), ptsSuperIndices.end(), rowsAddedSoFar);

      // pointers to Data and Theory object 'sum' error arrays
      // NOTE: Columns 'ErrData' and 'ErrTheoFunc' refer to where the errors are stored, not their type ('E' vs 'T')
      //        This is a little confusing -> change column nomenclature to e.g. 'ErrExp' and 'ErrTheo'?
      const std::vector<double>* errsData = &dataChildren[iChild]->GetSumError("AA", "AbsAvTot");
      const std::vector<double>* errsTheo = &theoChildren[iChild]->GetSumError("AA", "AbsAvTot");

      // list of total errors (AData and theory object errors added in quadrature)
      std::vector<double> errsTotal;
      errsTotal.resize(errsData->size());
      // add errors in quadrature
      for (unsigned int iRow = 0; iRow < errsTotal.size(); iRow++) {
         errsTotal[iRow] = std::sqrt(errsData->at(iRow)*errsData->at(iRow) + errsTheo->at(iRow)*errsTheo->at(iRow));
      }

      // containers for 'derived' columns (populate later, if requested)
      std::vector<double> dataTheoryRatios;
      std::vector<double> dataTheoryPulls;

      // store column headers, printf formats and pointers to column contents
      std::vector<std::string> columnHeads;
      std::vector<std::string> columnHeadFormats;
      std::vector<std::string> columnFormats;
      std::vector<const std::vector<double> *> columns;

      //display name of AData object
      info["Execute"] << "Printing data table columns for data table '" << dataChildren[iChild]->GetAlposName() << "'" << endl;

      // populate columns for current DataTable
      for (int iCol = 0; iCol < fColumnNames.size(); iCol++) {
         debug["Execute"] << "Looking up column #" << iCol << ": '" << fColumnNames[iCol] << "'" << endl;
         if (dataTable.find(fColumnNames[iCol]) != dataTable.end()) {
            debug["Execute"] << "Column '" << fColumnNames[iCol] << "' found in data table." << endl;
            columnHeads.push_back(fColumnNames[iCol]);
            columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
            columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
            columns.push_back(&(dataTable.at(fColumnNames[iCol])));
         }
         else {
            // FIXME: Add 'ErrExp' and 'ErrTheo' to separate errors by type? (see above)
            debug["Execute"] << "Column '" << fColumnNames[iCol] <<
            "' NOT found in data table. Exploring derived columns" << endl;
            // handle derived columns (not in data table)
            if (fColumnNames[iCol] == "iSuperData") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%12s ");
               columnFormats.push_back("%12d ");
               columns.push_back(&ptsSuperIndices);
            }
            else if (fColumnNames[iCol] == "iData") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%7s ");
               columnFormats.push_back("%7d ");
               columns.push_back(&ptsIndices);
            }
            else if (fColumnNames[iCol] == "ErrTot") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               // calculated these before loop, since they are needed for "Pulls", too
               columns.push_back(&errsTotal);
            }
            else if (fColumnNames[iCol].substr(0, 7) == "ErrData") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               if (fColumnNames[iCol].size() == 19) {
                  // requested data object errors with access string (like "ErrData(ASAbsAvTot)")
                  const std::vector<double>* errsDataAccStr = &dataChildren[iChild]->GetSumError(fColumnNames[iCol].substr(8, 2), fColumnNames[iCol].substr(10, 8));
                  columns.push_back(errsDataAccStr);
               }
               else if (fColumnNames[iCol].size() == 7) {
                  // requested data object errors without access string, use "AAAbsAvTot"
                  columns.push_back(errsData);
               }
               else {
                  // malformed error request
                  debug["Execute"] << "Malformed error request '" << fColumnNames[iCol] << ". Ignoring" <<
                  endl;
               }
            }
            else if (fColumnNames[iCol].substr(0, 11) == "ErrTheoFunc") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
	       // const std::vector<double>* errsTheo = &theoChildren[iChild]->GetSumError("AA", "AbsAvTot");//SS
	       // columns.push_back(errsTheo); //SS
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               if (fColumnNames[iCol].size() == 23) {
                  // requested data object errors with access string (like "ErrTheoFunc(ASAbsAvTot)")
                  const std::vector<double>* errsTheoAccStr = &theoChildren[iChild]->GetSumError(fColumnNames[iCol].substr(12, 2), fColumnNames[iCol].substr(14, 8));
                  columns.push_back(errsTheoAccStr);
               }
               else if (fColumnNames[iCol].size() == 11) {
                  // requested data object errors without access string, use "AAAbsAvTot"
                  columns.push_back(errsTheo);
               }
               else {
                  // malformed error request
                  debug["Execute"] << "Malformed error request '" << fColumnNames[iCol] << ". Ignoring" <<
                  endl;
               }
            }
            else if (fColumnNames[iCol] == "Theory") {
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               const std::vector<double>* theoPts = &theoChildren[iChild]->GetValues();
               columns.push_back(theoPts);
            }
            else if (fColumnNames[iCol] == "Data") {
               // add 'Data' as an alias for 'Sigma'
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               const std::vector<double>* dataPts = &dataChildren[iChild]->GetValues();
               columns.push_back(dataPts);
            }
            else if (fColumnNames[iCol] == "Ratio") {
               // add a 'Ratio' colum, meaning Data/Theory
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               const std::vector<double>* dataPts = &dataChildren[iChild]->GetValues();
               const std::vector<double>* theoPts = &theoChildren[iChild]->GetValues();
               dataTheoryRatios.resize(dataPts->size());
               // calculate ratios: ratio = data/theory
               for (unsigned int iRow = 0; iRow < dataTheoryRatios.size(); iRow++) {
                  dataTheoryRatios[iRow] = dataPts->at(iRow) / theoPts->at(iRow);
               }
               columns.push_back(&dataTheoryRatios);
            }
            else if (fColumnNames[iCol] == "Pull") {
               // add a 'Ratio' colum, meaning Data/Theory
               columnHeads.push_back(fColumnNames[iCol]);
               columnHeadFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "s ");
               columnFormats.push_back("%" + std::to_string(max(fColumnWidth, (int)fColumnNames[iCol].size())) + "g ");
               const std::vector<double>* dataPts = &dataChildren[iChild]->GetValues();
               const std::vector<double>* theoPts = &theoChildren[iChild]->GetValues();
               dataTheoryPulls.resize(dataPts->size());
               // calculate pulls: pull = (data-theory)/error
               for (unsigned int iRow = 0; iRow < dataTheoryRatios.size(); iRow++) {
                  dataTheoryPulls[iRow] = (dataPts->at(iRow) - theoPts->at(iRow)) / errsTotal[iRow];
               }
               columns.push_back(&dataTheoryPulls);
            }
            else {
               debug["Execute"] << "Column '" << fColumnNames[iCol] << "' NOT found among derived columns. Ignoring" <<
               endl;
            }
         }
      }// end column loop

      // begin printing current DataTable
      cout << endl;

      // buffer
      char buf[20000];
      int n=0;

      // open file
      string outputdirectory =  Alpos::Current()->Settings()->outputdir;
      ofstream file((outputdirectory+GetTaskName()+"_"+dataChildren[iChild]->GetAlposName()+".txt").c_str());


      // print column headers
      for (int iCol = 0; iCol < columnHeads.size(); iCol++) {
         n+=sprintf(&buf[n],columnHeadFormats[iCol].c_str(), columnHeads[iCol].c_str());
      }
      cout<<buf<<endl;
      file<<buf<<endl;
      n=0;

      // write data to histograms
      if ( dir ) {
	 dir->cd();
	 TDirectoryFile *subdir=new TDirectoryFile(dataChildren[iChild]->GetAlposName().c_str(),dataChildren[iChild]->GetAlposName().c_str());
	 subdir->cd();
	 int iErrTotCol=-1;
	 for (int iCol = 0; iCol < columnHeads.size(); iCol++) {
	    string name(columnHeads[iCol]);
	    if(name=="ErrTot") iErrTotCol=iCol;
	 }
	 for (int iCol = 0; iCol < columnHeads.size(); iCol++) {
	    string name(columnHeads[iCol]);
	    if(name=="iData") continue; // skip bin number
	    TH1D *h=new TH1D(name.c_str(),";bin",columns[0]->size(),-0.5,-0.5+columns[0]->size());
	    for (int iRow = 0; iRow < columns[0]->size(); iRow++) {
	       h->SetBinContent(iRow+1,columns[iCol]->at(iRow));
	       if((name=="Data")&&(iErrTotCol>=0)) {
		  h->SetBinError(iRow+1,columns[iErrTotCol]->at(iRow));
	       }
	    }
	 }
      }

      // print table contents
      for (int iRow = 0; iRow < columns[0]->size(); iRow++) {
         for (int iCol = 0; iCol < columnHeads.size(); iCol++) {
            if ((columnHeads[iCol] == "iData") || (columnHeads[iCol] == "iSuperData")) {
               // print 'iData' and 'iSuperData' indices as integers
               n+=sprintf(&buf[n],columnFormats[iCol].c_str(), (int) (*columns[iCol])[iRow] );
            }
            else {
	       n+=sprintf(&buf[n],columnFormats[iCol].c_str(), (*columns[iCol])[iRow] );

            }
         }
	 cout<<buf<<endl;
	 file<<buf<<endl;
	 n=0;
      }

      rowsAddedSoFar += dataChildren[iChild]->N();  // keep track of # of rows added for each dataset
      cout<<endl<<endl;
   }

   if ( dir ) dir->Write();
   cout<<" +---------------------------------------------------------- "<<endl;
   cout<<endl;

   // --- job done successfully
   return true;

}


//____________________________________________________________________________________ //
