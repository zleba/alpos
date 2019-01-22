// DB 15.01.2015

#include "alpos/AData.h"
#include "alpos/ASubsetData.h"
#include "alpos/ASubsetFunction.h"
#include "alpos/AFactory.h"
#include "alpos/Alpos.h"
//#include "alpos/AlposTools.h"
#include "fastnlotk/read_steer.h"
#include <iostream>
#include <algorithm>


using namespace std;


// __________________________________________________________________________________________ //
const std::vector<std::string> AData::fRequirements = {}; //< List of all AParm's which this function depends on
const std::vector<std::string> AData::fStopFurtherNotification = {}; //< List of Parm's which have changed, but this function does not notify further dependencies
const std::string AData::fFunctionName = "Data"; //< The function's name


// __________________________________________________________________________________________ //
AData::AData(const std::string& name) : AParmFuncBase<double>(name+"_Data") { 
   SetClassName("AData");
   fDataName = name;
}


// __________________________________________________________________________________________ //
AData::AData(const std::string& name, bool klar) : AParmFuncBase<double>(name) { 
   SetClassName("AData");
   fDataName = name;
}


// __________________________________________________________________________________________ //
AData::~AData() {
}


// ___________________________________________________________________________________________ //
bool AData::Init() {
   //! Init is once called for each function
   //! return true if initialization was successful.
   debug["Init"]<<"Init AData: "<<fDataName<<endl;

   vector<double> Sigma = DOUBLE_COL_NS(Data,Sigma,fDataName);

   fValue = Sigma;
   fError.resize(fValue.size());
   
   if ( fValue.empty() ) {
      warn["Init"]<<"The column 'Sigma' seems not to exist in table Data in dataset '"<<fDataName<<"', or has no values. Exiting"<<endl;
      exit(1);
   }
   vector<vector<double > > data = DOUBLE_TAB_NS(Data,fDataName);
   vector<string> header = TABLEHEADER_NS(Data,fDataName);
   // check briefly size of data table
   if ( data.front().size() != data.back().size() ) {
      warn["Init"]<<"The input data table 'Data' in "<<fDataName<<" does not contain the same number of entries for each column or row."<<endl;
   }

   // --- add 'Data' to fOrigData
   AddColumnToDataTable("Sigma",Sigma); // first 'Sigma'
   for ( unsigned int i=0 ; i<header.size() ; i++ ){
      if ( header[i]!="Sigma")
	 AddColumnToDataTable(header[i],read_steer::getdoublecolumn("Data",header[i],fDataName));
   }

   // --- init errors
   ReadErrors(fDataName); // init fAllErrors

   // --- init subsets;
   InitSubsets();

   // --- init 'subset' with applied cuts
   InitSubsetWithCuts();

   debug["Init"]<<"fValue.size() = "<<fValue.size()<<endl;
   fError.resize(fValue.size());

   // data should not change:
   this->SetIsConstant();

   return true;
}


// __________________________________________________________________________________________ //
bool AData::Update() {
   // nothing todo
   return true;
}


// __________________________________________________________________________________________ //
//std::map<std::string,AError> 
void AData::ReadErrors(const std::string& DataName) {
   //! Read errors from the steering file of dataset 'DataName'
   //! returns a vector of all AError's

   using namespace AlposTools;
   debug["ReadErrors"]<<"DataName: "<<DataName<<endl;
   //std::map<std::string,AError> AllErrors; // return value;

   // get values from steering
   // check header:
   vector<double> Sigma    = DOUBLE_COL_NS(Data,Sigma,DataName);
   vector<string> errname  = STRING_COL_NS(Errors,ErrorName,DataName);
   vector<string> col      = STRING_COL_NS(Errors,Column,DataName);
   //vector<string> covma    = STRING_COL_NS(Errors,Correlation,DataName);
   vector<string> nature   = STRING_COL_NS(Errors,Nature,DataName); // not yet implemented
   if ( nature.empty() ) {
      warn["ReadErrors"]<<"Column 'Nature' in table 'Errors' in steering file '"<<DataName<<" not found."<<endl;
      nature.resize(errname.size(),"not specified");
   }
   if ( nature.size() != errname.size() ) {
      warn["ReadErrors"]<<"Column 'Nature' in table 'Errors' in steering file '"<<DataName<<" has different number of entries than ErrorName."<<endl;
      nature.resize(Sigma.size(),"not specified");
   }
   string ErrorSet         = STRING_NS(ErrorSet,DataName);
   string units            = STRING_NS(ErrorUnit,DataName);
   bool relval = true;
   
   vector<string> header = TABLEHEADER_NS(Errors,DataName);
   bool OldFormat = std::find(header.begin(),header.end(),"Correlation") != header.end();
   //vector<string> cors = STRING_COL_NS(Errors,Correlation,DataName);
   vector<string> type = OldFormat ? STRING_COL_NS(Errors,Correlation,DataName) : STRING_COL_NS(Errors,Type,DataName);


   // loop over all specified error
   for ( unsigned int i=0; i<errname.size() ;i++ ){
      // Build AError object, and read values
      vector<double> ErrUp, ErrDn;
      string eup = col[i];
      string edn = col[i];;
      size_t pos = col[i].find(":");
      if ( pos != string::npos ) {
	 eup = col[i].substr(0,pos);
	 edn = col[i].substr(pos+1);
	 debug["ReadErrors"]<<"Reading columns '"<<eup<<"' and '"<<edn<<"' as up and down errors."<<endl;
      }
      if ( read_steer::CheckNumber(eup) ) { //! Error is specified directly
	 ErrUp.clear(); ErrDn.clear();
	 double evalup = stod(eup);
	 ErrUp.resize(fValue.size(),evalup);
	 if ( pos != string::npos ) {
	    double evaldn = stod(edn);
	    ErrDn.resize(fValue.size(),evaldn);
	 }
	 else
	    ErrDn=ErrUp;
	 if ( ErrDn[0]*ErrUp[0]>0 ) ErrDn *= -1.;
      }
      else { // error is specified in column
	 ErrUp = read_steer::getdoublecolumn("Data",eup,DataName);
	 ErrDn = read_steer::getdoublecolumn("Data",edn,DataName);
	 if ( ErrUp.empty() || ErrDn.empty() ) {
	    error["ReadErrors"]<<"Error '"<<col[i]<<"' could not be read properly. Exiting."<<endl;
	    exit(1);
	 }
	 if ( pos == string::npos ) {
	    int i0 = 0;
	    while (ErrDn[i0]==0 && ErrUp[i0]==0 && i0<ErrDn.size())i0++;
	    if ( ErrDn[i0]*ErrUp[i0]>0 ) ErrDn *= -1.;
	    //if ( ErrDn[0]*ErrUp[0]>0 ) ErrDn *= -1.;
	 }
      }

      // adjust units
      if ( units=="Relative" || units=="relative") relval=true;
      else if (units=="Absolute" || units=="absolute") relval=false;
      else if (units=="Percent" || units=="percent") {
	 relval=true;
	 ErrUp /= 100.;
	 ErrDn /= 100.;
      }
      else {
	 error["ReadErrors"]<<"The key 'ErrorUnit' in dataset '"<<DataName<<" must be 'Absolute', 'Relative' or 'Percent', but is: "<<units<<". Exiting."<<endl; exit(1);
      }

      // --- instantiate new error
      AError tmp(errname[i],ErrorSet);

      // --- specify 'type'
      //if ( OldFormat ) {
      if (type[i].find("Matrix") !=string::npos) {
	 type[i].replace(type[i].find("Matrix"),6,"C");
      }
      if (type[i].find("matrix") !=string::npos ) {
	 type[i].replace(type[i].find("matrix"),6,"C");
      }

      // --- sanity checks
      if ( type[i].find("S") != string::npos && type[i].find("Y") !=string::npos ) {
	 error["ReadErrors"]<<"Error type may be either of statistical ('S') or systematic ('Y') nature but is both in: "<<type[i]<<" for error source '"<<tmp.GetErrorName()<<"'."<<endl;
	 exit(3);
      }
      if ( type[i].find("E") != string::npos && type[i].find("T") !=string::npos ) {
	 error["ReadErrors"]<<"Error type may be either of experimental ('E') or theoretical ('T') nature but is both in: "<<type[i]<<" for error source '"<<tmp.GetErrorName()<<"'."<<endl;
	 exit(3);
      }
      if ( type[i].find("T") != string::npos && Alpos::Current()->Settings()->IgnoreTheoryErrors ) {
	 info["ReadErrors"]<<"Global flag 'IgnoreTheoryErrors' set to true. Ignoring error '"<<errname[i]<<"' because it is of type 'theory'."<<endl;
	 continue;
      }
      if ( type[i].find("N") != string::npos && type[i].find("C") !=string::npos ) {
	 error["ReadErrors"]<<"Error type may be either specified as matrix ('C') or column ('N') but is both in: "<<type[i]<<" for error source '"<<tmp.GetErrorName()<<"'."<<endl;
	 exit(3);
      }
      if ( type[i].find("A") != string::npos && type[i].find("M") !=string::npos ) {
	 error["ReadErrors"]<<"Error type may be either specified as multiplicative ('M') or additive ('A') but is both in: "<<type[i]<<" for error source '"<<tmp.GetErrorName()<<"'."<<endl;
	 exit(3);
      }

      // --- add error to list of errors
      fAllErrors[tmp.GetErrorName()] = tmp;
      AError& err = fAllErrors[tmp.GetErrorName()];

      
      // --- stat or systematic
      if ( type[i].find("S") != string::npos || type[i].find("s") !=string::npos ) err.SetIsStat();
      if ( /*OldFormat &&*/ i==0 && ( err.GetErrorName().find("stat") != string::npos || err.GetErrorName().find("Stat") != string::npos ) ) err.SetIsStat();
      if ( type[i].find("Y") != string::npos || type[i].find("Y") !=string::npos ) err.SetIsStat(false); // default
      // --- experimental or theoretical nature
      if ( type[i].find("E") != string::npos || type[i].find("e") !=string::npos ) err.SetIsTheo(false); //default
      if ( type[i].find("T") != string::npos || type[i].find("t") !=string::npos ) err.SetIsTheo();
      // --- additive or multiplicative error
      if ( type[i].find("M") != string::npos || type[i].find("m") !=string::npos ) err.SetIsMult(); // default
      if ( type[i].find("A") != string::npos || type[i].find("a") !=string::npos ) err.SetIsMult(false); // default

      if (err.GetIsMult())
         this->SetHasMultErrors();

      // --- specify correlations (corrfrac or matrix)
      // Matrix type uncertainty
      if ( type[i].find("C") != string::npos || type[i].find("c") !=string::npos ){
	 string format  = read_steer::getstring(errname[i]+"_Matrix_Format",DataName);
	 string mattype = read_steer::getstring(errname[i]+"_Matrix_Type",DataName);
	 vector<vector<double> > mvals = read_steer::getdoubletable(errname[i]+"_Matrix",DataName);
	 vector<string> mhead = read_steer::gettableheader(errname[i]+"_Matrix",DataName);
	 vector<string> dhead = TABLEHEADER_NS(Data,fDataName); // 'Data'
	 vector<vector<double> > dat = DOUBLE_TAB_NS(Data,fDataName); // 'Data'
	 //TMatrixDSym mat = ReadMatrix(mvals,format); // call for "Matrix"
	 TMatrixDSym Mat = ReadMatrix(mvals,format,mhead,dhead,dat); 	 
	 info["ReadErrors"]<<"mat: ncol="<<Mat.GetNcols() <<"\tnrow="<< Mat.GetNrows()<<endl;
	 // recalculate 'percent' notation
	 if ( mattype == "CorrelationPercent" || mattype=="correlationpercent" ) {
	    for ( unsigned int x = 0 ; x<dat.size() ; x++ ) {
	       for ( unsigned int y = 0 ; y<dat.size() ; y++ ) 
	       	  Mat[x][y] = Mat[x][y]/100.;
	    }
	    mattype = "Correlation";
	 }

	 //
	 bool bMIsCorr = true;
	 if ( mattype == "Covariance" ) {
	    debug["ReadErrors"]<<"A check if the diagonal elements are about the value in the 'Data'-table could be helpful."<<endl;
	    bMIsCorr=false;	    // i.e. mat[i][i] ~= Data_Col<errorname>[i]**2
	 }
	 else if ( mattype == "Correlation" || mattype == "correlation" ) {
	    // set diagonal elements.
	    bMIsCorr = true;
	    for ( unsigned int x = 0 ; x<dat.size() ; x++ )
	       Mat[x][x] = 1.;
	 }

	 debug["ReadErrors"]<<"mattype="<<mattype<<"\tbMIsCorr="<<bMIsCorr<<", relval="<<relval<<"\te[0]="<<ErrUp[0]<<"\te[3]="<<ErrUp[3]<<endl; //
	 // set error
	 err.SetMatrixError(ErrUp, Mat,Sigma,relval,bMIsCorr,nature[i]);
	 err.SetColUpDn(eup,edn);

      }
      // 'normal' uncertainty
      else {//if ( read_steer::CheckNumber(cors[i]) ) {
	 double corrfrac = 1;
	 if ( type[i].find("0") != string::npos || type[i].find(".") !=string::npos || type[i].find("1") !=string::npos  ) {
	    if ( type[i].find("0") != string::npos && type[i].find(".") ==string::npos )
	       corrfrac = 0;
	    else if ( type[i].find("1") !=string::npos && type[i].find(".") ==string::npos) 
	       corrfrac = 1;
	    else if ( type[i].find("1.") !=string::npos ) 
	       corrfrac = 1;
	    else {
	       size_t send ;
	       corrfrac = stod(type[i].substr(type[i].find(".")),&send);
	    }
	 }
	 else {
	    debug["ReadErrors"]<<"Take correlation coefficient '1' for error source '"<< err.GetErrorName()<<"'."<<endl;
	    err.SetCorrelatedFraction(1);
	 }
	 // set error
	 //err.SetAsymError(ErrUp,ErrDn,Sigma,relval,corrfrac,nature[i],AError::kSignImprovedQuadratic);
	 err.SetAsymError(ErrUp,ErrDn,Sigma,relval,corrfrac,nature[i],Alpos::Current()->Settings()->ErrorSymmetrization);
	 err.SetColUpDn(eup,edn);
	 err.SetCorrelatedFraction(corrfrac);
      }
   } // errors done.
   

   fSumErrMats.clear();
   fSumErrors.clear();

   // not needed
   //CalculateMatrices();

   // TODO: remnants of old interface -> remove once these become accessible through new interface
   // fErrMat, fErrMatRel (matrix-type uncertainties) excluding stat
   fErrMat.clear();
   fErrMatRel.clear();
   fErrMat.resize(fValue.size());
   fErrMatRel.resize(fValue.size());
   for ( const auto& ierr : fAllErrors ) {
      if ( !ierr.second.GetIsStat() && ierr.second.GetIsMatType()) {
         for ( unsigned int iea = 0 ; iea<fValue.size() ; iea++ ) {
            fErrMat[iea]    += pow(ierr.second.GetMatError()[iea],2) ;
            fErrMatRel[iea] += pow(ierr.second.GetMatErrorRel()[iea],2) ;
         }
      }
   }
   for ( unsigned int iea = 0 ; iea<fValue.size() ; iea++ ) {
      fErrMat[iea]    = sqrt(fErrMat[iea]);
      fErrMatRel[iea] = sqrt(fErrMatRel[iea]);
   }

   //return AllErrors;
}


// __________________________________________________________________________________________ //
TMatrixDSym AData::ReadMatrix(const std::vector<vector<double> >& values, const std::string& format, const std::vector<std::string>& matHeader, const vector<std::string>& DataHeader, const vector<vector<double> >& DataVals ) {
   //!  Read a (symmetrix) matrix in the given 'format'
   //!  Supported formats are: "Matrix" and "SingleValues"
   //!  'SingleValues' require in addition the DataTable (i.e. the DataHeader and DataVals) from
   //!  the steering file.
   //!
   //!  Format: "Matrix"
   //!     Specify each element of the (half-)matrix in the lower-left half-matrix
   //!     or the full (symmetric) matrix, e.g.:
   //!     Stat_Matrix {{
   //!        # empty line for 'header'
   //!        1.0
   //!        0.3  1.0
   //!        0.1  0.4  1.0
   //!     }}
   //!       
   //!  Format: "SingleValues"
   //!     Specify each value of the (symmetric) matrix in a single line and 
   //!     give details in each row by using the 'names' of the 'Data' table columns.
   //!     Use a set of column-names to uniquely identify a row (the first occurence is used).
   //!     Use the same identifier(s) twice, e.g.:
   //!     Stat_Matrix {{
   //!         q2min  ptmin   q2min   ptmin   value
   //!          150.0  7.0    150.0   11.0    0.16
   //!          400.0  7.0    400.0   18.0    0.06
   //!         5000.0  7.0   5000.0   11.0    0.16
   //!     }}
   //!     or (often simpler):
   //!     Stat_Matrix {{
   //!         ID      ID   value
   //!          1      2     0.3
   //!          1      3     0.1
   //!          2      3     0.4
   //!     }}
   //!    
   //!   

   TMatrixDSym mat;
   if ( format=="Matrix" || format=="matrix" ) {
      //check correct size
      if ( values.size() != values.back().size() ){
	 error["ReadMatrix"]<<"The input table 'values', has the wrong format: n-rows="<<values.size()<<", n-columns="<<values.back().size()<<endl;
	 return mat;
      }
      mat.ResizeTo(values.size(),values.size());
      //loop over values and fill into matrix:
      for ( unsigned int y = 0 ; y<values.size() ; y++ ) {
	 // check again correct size:
	 //if ( values[y].size() != y+1 && values[y].size() != values.size() ) 
	 if ( values[y].size() < y+1 || values[y].size() > values.size() ) 
	    error["ReadMatrix"]<<"Wrong format of input table in column "<<y<<endl;
	 for ( unsigned int x = 0 ; x<y+1 ; x++ ) {
	    mat[y][x]  = values[y][x];
	    mat[x][y]  = values[y][x];
	 }
      }
   }
   else if ( format=="SingleValues" || format=="singlevalues" ) {
      // get number of identifiers:
      if ( (values[0].size()-1)%2 != 0 )
	 error["ReadMatrix"]<<"Wrong size of columns."<<endl;
      mat.ResizeTo(DataVals.size(),DataVals.size());
      int nid = (values[0].size()-1)/2;
      vector<int> pos;
      for ( int i=0 ; i<nid ;i++ ) {
	 int ipos=-1;
	 debug["ReadMatrix"]<<"Looking for column with name '"<<matHeader[i]<<"' in 'Data'-table."<<endl;
	 if ( matHeader[i] != matHeader[i+nid] )
	    error["ReadMatrix"]<<"Identifiers must be present twice, but the second one ws not found for "<<matHeader[i]<<endl;
	 for ( unsigned int h = 0 ; h<DataHeader.size() ; h++ ) {
	    if ( matHeader[i]==DataHeader[h] ) {
	       if ( ipos == -1 ) {
		  ipos=h;
		  pos.push_back(h);
		  //cout<<"found col. matHeader[i]="<<matHeader[i]<<", pos="<<h<<endl;
	       }
	       else 
		  warn["ReadMatrix"]<<"Found column '"<<matHeader[i]<<"' multiple times in 'Data'-table."<<endl;
	    }
	 }
	 if ( ipos==-1) {
	    error["ReadMatrix"]<<"Could not identify the requested column '"<<matHeader[i]<<"'."<<endl;
	    exit(1);
	 }
      }
      if ( int(pos.size()) != nid ) {
	 error["ReadMatrix"]<<"Could not identify the requested columns."<<endl;
      }

      // fill matrix
      int x=-1,y=-1;
      for ( unsigned int i = 0 ; i<values.size() ;i++ ) { // all specified values
	 for ( unsigned int d = 0 ; d<DataVals.size() ;d++ ) { // find the two rows
	    if ( nid == 1 ) {
	       if ( DataVals[d][pos[0]] == values[i][0])      x=d;
	       if ( DataVals[d][pos[0]] == values[i][0+nid] )  y=d;
	    }
	    else if ( nid == 2 ) {
	       if ( DataVals[d][pos[0]] == values[i][0]     && DataVals[d][pos[1]] == values[i][1] )     x=d;
	       if ( DataVals[d][pos[0]] == values[i][0+nid] && DataVals[d][pos[1]] == values[i][1+nid] ) y=d;
	    }
	    else if ( nid == 3 ) {
	       if ( DataVals[d][pos[0]] == values[i][0]     && DataVals[d][pos[1]] == values[i][1]     && DataVals[d][pos[2]] == values[i][2] )     x=d;
	       if ( DataVals[d][pos[0]] == values[i][0+nid] && DataVals[d][pos[1]] == values[i][1+nid] && DataVals[d][pos[2]] == values[i][2+nid] ) y=d;
	    }
	    else if ( nid == 4 ) {
	       if ( DataVals[d][pos[0]] == values[i][0] && DataVals[d][pos[1]] == values[i][1] && DataVals[d][pos[2]] == values[i][2] && DataVals[d][pos[3]] == values[i][3] ) x=d;
	       if ( DataVals[d][pos[0]] == values[i][0+nid] && DataVals[d][pos[1]] == values[i][1+nid] && DataVals[d][pos[2]] == values[i][2+nid] && DataVals[d][pos[3]] == values[i][3+nid] ) y=d;
	    }
	 }
	 if ( x==-1 ) { error["ReadMatrix"]<<"Could not find x-index in row "<<i<<endl; continue; }
	 if ( y==-1 ) { error["ReadMatrix"]<<"Could not find y-index in row "<<i<<endl; continue; }
	 //cout<<"Found values. x="<<x<<"\ty="<<y<<"\tval="<<values[i].back()<<endl;
	 mat[x][y] = values[i].back();
	 mat[y][x] = values[i].back();
      }
   }
   else {
      error["ReadMatrix"]<<"Unknown format: '"<<format<<"'. Don't know what to do."<<endl;
   }
   return mat;
}



// __________________________________________________________________________________________ //
bool AData::AddColumnToDataTable(const std::string& colname, const std::vector<double>& values) {
   //! Add one additional column to the data table
   //! return true, if no conflict or problem occured
   //! return false, if a columns with same name already exists
   //! return false, if size of array does not match the size of the array 'Sigma'
   
   if (fOrigData.count(colname) > 0 ) {
      warn["AddColumnToDataTable"]<<"A column with the same name already exists: '"<<colname<<"'. Ignoring call."<<endl;
      return false;
   }
   if ( fOrigData.count("Sigma") == 0 && colname!= "Sigma" ) {
      info["AddColumnToDataTable"]<<"Cannot check size of input array since no column with name 'Sigma' exists yet in data table."<<endl;
   }
   if ( fOrigData.count("Sigma") > 0 && fOrigData["Sigma"].size() != values.size() ) {
      warn["AddColumnToDataTable"]<<"Size of values ("<<values.size()<<") of column '"<<colname<<"' is different than numbers of entries in 'Sigma' ("<<fOrigData["Sigma"].size()<<")."<<endl;
      warn["AddColumnToDataTable"]<<"This may yield sever problems later."<<endl;
   }
   fOrigData[colname] = values;
   return true;
}


// __________________________________________________________________________________________ //
const TMatrixD& AData::CheckReturnInverse(TMatrixD& Inv, const TMatrixDSym& Mat) const {
   //!< return the inverse matrix
   //!< if matrix does not exist, invert it
   if ( Inv.GetNrows() != Mat.GetNrows() ) {
      // --- invert matrix
      info["CheckReturnInverse"]<<"Inverting matrix."<<endl;
      //TMatrixD inv = AlposTools::InvertLU(Mat);
      TMatrixD inv = AlposTools::InvertChol(Mat);
      Inv.ResizeTo(inv);
      Inv = TMatrixD(inv);
   }
   return Inv;
}


// __________________________________________________________________________________________ //
void AData::InitSubsets(){
   //! Init map of subset 'ValidPoints'
   const vector<string>& subs = STRING_ARR_NS(Subsets,fDataName) ;
   for ( const auto& sub : subs ) {
      const vector<double>& val = read_steer::getdoublecolumn("Data",sub,fDataName);
      set<double> subval;
      for ( const auto& v: val ) subval.insert(v);
      for ( const auto& subv : subval ) {
	 int n=0;
	 vector<bool> valpts(fValue.size());
	 for ( unsigned int iv = 0 ; iv<val.size() ;iv++ ) {
	    if ( val[iv] == subv ) { valpts[iv] = true; n++; }
	    else valpts[iv] = false;
	 }
	 double ssub = subv;
	 string subname = sub+"=="+to_string(ssub);
	 fSubsets[subname] = valpts;
	 debug["InitSubsets"]<<"Set subset '"<<subname<<"' with "<< n << " valid points."<<endl;
      }
   }

   /*
   //! Init Subset functions
   for ( const auto& isub : fSubsets ) {
      // subsetdata
      string subsetname = fDataName+"_"+isub.first;
      ASubsetData* SupDa =  (ASubsetData* )AFactory::FunctionFactory("SubsetData",subsetname);
      SupDa->SetRequirementValidPoints(fDataName,isub.second);
      // subset theory
      ASubsetFunction* SupTh =  (ASubsetFunction* )AFactory::FunctionFactory("SubsetFunction",subsetname);
      SupTh->SetRequirementValidPoints(fDataName,isub.second);
   }
   */
}


struct sCut {
   std::string valname;
   std::string cmp;
   double vcut;
};

// __________________________________________________________________________________________ //
void AData::InitSubsetWithCuts(){
   using namespace AlposTools;
   
   //! Init vector of 'ValidPoints'
   vector<string> cutsIn;
   if (EXISTARRAY_NS(Cuts,fDataName))
      cutsIn += STRING_ARR_NS(Cuts,fDataName) ;
   if (EXISTARRAY_NS(CutsMainSteering,fDataName))
      cutsIn += STRING_ARR_NS(CutsMainSteering,fDataName) ;

   if (cutsIn.empty()) return;// nothing todo

   // building vector<sCut>
   string subname="cuts:";
   vector<sCut> cuts;
   for ( auto cut : cutsIn ) {
      subname+=cut;
      subname+="&&";
      string valname=cut, cmpop, scut;
      const vector<string> cmpops{"==",">=","<=",">","<","!=",".gt.",".lt.",".eq.",".ge.",".le."}; // mind the order!
      for ( auto op : cmpops ) {
	 if ( cut.find(op)!=string::npos ) {
	    read_steer::separatetag(valname, scut, op);
	    cmpop=op;
	    break; // respect the order!
	 }
      }
      sCut acut;
      acut.valname=valname;
      acut.cmp=cmpop;
      acut.vcut=stod(scut);
      info["InitSubsetWithCuts"]<<"Found cut with properties: "<<valname<<"\t"<<cmpop<<"\t"<<acut.vcut<<endl;
      cuts.push_back(acut);
   }   
   subname = subname.substr(0, subname.size()-2);

   // build vector of valid points
   vector<bool> valpts(fValue.size(),true);
   for ( const auto& icut : cuts ) { // cut1 && cut2 && cut3 ...
      const vector<double>& val = read_steer::getdoublecolumn("Data",icut.valname,fDataName);
      for ( unsigned int iv = 0 ; iv<val.size() ;iv++ ) {
	 if ( !valpts[iv] ) continue;
	 double cval = icut.vcut;
	 if ( icut.cmp=="==" && !(val[iv]==cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp==".eq." && !(val[iv]==cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp=="!=" && !(val[iv]!=cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp==">=" && !(val[iv]>=cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp=="<=" && !(val[iv]<=cval) ) 
	    valpts[iv]=false; 
	 if ( icut.cmp==">" && !(val[iv]>cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp=="<" && !(val[iv]<cval) ) 
	    valpts[iv]=false; 
	 if ( icut.cmp==".ge." && !(val[iv]>=cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp==".le." && !(val[iv]<=cval) ) 
	    valpts[iv]=false; 
	 if ( icut.cmp==".gt." && !(val[iv]>cval) )  
	    valpts[iv]=false; 
	 if ( icut.cmp==".lt." && !(val[iv]<cval) ) 
	    valpts[iv]=false; 
      }      
   }
   int n=0;
   for ( auto iv : valpts ) if (iv) n++; // just count
   info["InitSubsets"]<<"Set subset including cuts with "<< n << " valid points."<<endl;
   fSubsets[subname] = valpts;

   /*
   for ( const auto& sub : subs ) {
      const vector<double>& val = read_steer::getdoublecolumn("Data",sub,fDataName);
      set<double> subval;
      for ( const auto& v: val ) subval.insert(v);
      for ( const auto& subv : subval ) {
	 int n=0;
	 vector<bool> valpts(fValue.size());
	 for ( unsigned int iv = 0 ; iv<val.size() ;iv++ ) {
	    if ( val[iv] == subv ) { valpts[iv] = true; n++; }
	    else valpts[iv] = false;
	 }
	 double ssub = subv;
	 string subname = sub+"="+to_string(ssub);
	 fSubsets[subname] = valpts;
	 debug["InitSubsets"]<<"Set subset '"<<subname<<"' with "<< n << " valid points."<<endl;
      }
   }
   */

   /*
   //! Init Subset functions
   for ( const auto& isub : fSubsets ) {
      // subsetdata
      string subsetname = fDataName+"_"+isub.first;
      ASubsetData* SupDa =  (ASubsetData* )AFactory::FunctionFactory("SubsetData",subsetname);
      SupDa->SetRequirementValidPoints(fDataName,isub.second);
      // subset theory
      ASubsetFunction* SupTh =  (ASubsetFunction* )AFactory::FunctionFactory("SubsetFunction",subsetname);
      SupTh->SetRequirementValidPoints(fDataName,isub.second);
   }
   */
}


// __________________________________________________________________________________________ //
