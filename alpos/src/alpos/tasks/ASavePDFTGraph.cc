#include "alpos/tasks/ASavePDFTGraph.h"
#include <iostream>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMatrixDSymEigen.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TSystem.h>
#include "alpos/AFactory.h"
#include "alpos/AChisq.h"
#include "alpos/ASuperData.h"
#include "alpos/ASuperTheory.h"
#include "alpos/functions/ALhapdf6.h"

/** 
 ASavePDFTGraph

 * Save TGraphs of PDFs

 */

using namespace std;

const string ASavePDFTGraph::fTaskType = "SavePDFTGraph";

//____________________________________________________________________________________ //
ASavePDFTGraph::ASavePDFTGraph(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("ASavePDFTGraph");
   //! Important: create always new result-object here!
   fResult = new ASavePDFTGraphResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
ASavePDFTGraph::~ASavePDFTGraph(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool ASavePDFTGraph::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   debug["Init"]<<endl;
   
   return true;
}


//____________________________________________________________________________________ //
bool ASavePDFTGraph::Execute(){
   debug["Execute"]<<endl;
   
   // init output
   //TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   TDirectory* file = Alpos::Current()->Settings()->rootoutput;
   if ( !file ) {
      info["Execute"]<<"No output directory or file specified in steering. Skipping task."<<endl;
      return true;
   }
   TDirectory* taskdir = file->mkdir(GetTaskName().c_str());
   taskdir->cd();
   // the PDF-function
   string pdffunc = STRING_NS(PDF,NS());
   PAR_ANY(pdffunc); // update
   info["Execute"]<<"Evaluating PDF from theory parameter: "<<pdffunc<<" which is of type '"<<
      TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName()<<"'."<<endl;
   // alpha_s evolution (as it is used for PDF evolution)
   string asfunc = STRING_NS(Alphas,NS());
   if ( asfunc=="" ) {
      warn["Execute"]<<"Alphas evolution not specified. Please provde Alphas-evolution function with key 'Alphas'"<<endl;
   }
   else {
      PAR_ANY(asfunc); // update
      info["Execute"]<<"Evaluating Alphas evolution from theory parameter: "<<asfunc<<" which is of type '"<<
	 TheoryHandler::Handler()->GetFuncD(asfunc)->GetFunctionName()<<"'."<<endl;
   }
   string WriteLHAPDF = STRING_NS(LHAPDFoutput,NS());
   if ( WriteLHAPDF!="" && asfunc=="" ) {
      error["Execute"]<<"LHAPDF output specified, but no Alphas function given. Cannot write LHAPDF grids."<<endl;
      WriteLHAPDF="";
   }
   else if ( WriteLHAPDF!="" ) info["Execute"]<<"LHAPDF compatible files will be written into directory: "<<WriteLHAPDF<<endl;

   
   // x and q2 spacing
   vector<double> q2val = DOUBLE_ARR_NS(Q2Values,NS());
   int nx = INT_NS(nx,NS());
   double xmin = DOUBLE_NS(xmin,NS());
   double lxstep = (log(1)-log(xmin))/nx;
   vector<double> xval = {xmin};
   for ( int ix = 0 ; ix < nx-1 ; ix++ ) {
      xval.push_back(exp((log(xval.back())+lxstep)));
   }
   xval.push_back(1);

   int nPDF=1;
   bool LHAType =false;
   bool LHAError=false;
   bool FitType=false;
   bool WeightNNPDF = false;
   vector<double> vWgtNNPDF;
   vector<double> vlWgtNNPDF;
   vector<double> vChi2NNPDF;
   if ( EXIST_NS(WeightNNPDF,NS()) ) WeightNNPDF = BOOL_NS(WeightNNPDF,NS());
   if ( WeightNNPDF ) info["execute"]<<"Calculate weights for any NNPDF replica."<<endl;

   int Mem0 = -1; // PDFSet in the beginning
   // --- LHAType
   if ( TheoryHandler::Handler()->GetFuncD(pdffunc)->GetFunctionName().find(ALhapdf6::fFunctionName) !=std::string::npos ) {
      LHAType=true;
      ALhapdf6* lha = (ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc);
      nPDF = lha->GetPDFSet()->size();
      info["Execute"]<<"Found "<<nPDF<<" PDF members."<<endl;
      Mem0 = PAR_ANY(pdffunc+".PDFSet");
      LHAError = !( lha->GetPDFSet()->errorType() == "unknown" || lha->GetPDFSet()->errorType() == "UNKNOWN" );
   }

   // --- WeightNNPDF setings
   // loop over all NNPDF replicas, and calculate weight.
   // later: use these weights for tgraphs
   if ( WeightNNPDF ) {
      LHAPDF::PDFSet* lhapdfset = LHAType ?
	 ((ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc))->GetPDFSet() : 
	 NULL;
      if ( !lhapdfset ) {
	 error["execute"]<<"This is not a LHAPDF PDF set. Ignoring 'WeigthNNPDF' flag."<<endl;
	 WeightNNPDF=false;
      }
      if ( lhapdfset->errorType() == "replicas" ) 
	 info["execute"]<<"This is a (NNPDF) PDF replica set. Alright."<<endl;
      else {
	 warn["execute"]<<"This is not a replica PDF!"<<endl;
	 //error["execute"]<<"This is not a replica PDF! Ignoring 'WeigthNNPDF' flag."<<endl;
	 //WeightNNPDF=false;
      }
      WeightNNPDF &= nPDF>1 && LHAError;
      if ( !WeightNNPDF ) 
	 error["execute"]<<"Incompatible settings (nPDF<=0, or !LHAError). Ignoring 'WeigthNNPDF' flag."<<endl;
   }


   // --- FitType
   TMatrixD V,VT;//EVmat;
   TMatrixDSym cov;
   TVectorD eigV;
   // --- n-1 Eig (for alpha_s absorbed in a single shift)
   TMatrixD Vm1,VTm1;//EVmat;
   TMatrixDSym covm1;
   TVectorD eigVm1,V1;
   
   TFitResult* fitresult = (TFitResult*)file->Get("fitresult");
   vector<double> initialparams;
   vector<string> fitparamnames;
   const bool DoAlphasSingleShift = false;
   if ( fitresult ) {
      FitType = true;
      int nxn = fitresult->GetCovarianceMatrix().GetNcols();
      nPDF = 1+nxn*2;
      cov.ResizeTo(nxn,nxn);
      cov = fitresult->GetCovarianceMatrix();
      // Eigenvalues and eigenvectors of a real symmetric matrix:
      // If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
      // diagonal and the eigenvector matrix V is orthogonal. 
      // That is, the diagonal values of D are the eigenvalues, and V*V' = I, 
      // where I is the identity matrix. 
      // The columns of V represent the eigenvectors in the sense that A*V = V*D. 

      TMatrixDSymEigen eigen(cov);
      eigV.ResizeTo(nxn);
      eigV = eigen.GetEigenValues();
      info["execute"]<<"Printing Eigenvalues: "<<endl;
      eigV.Print();

      V.ResizeTo(nxn,nxn);
      V = eigen.GetEigenVectors();
      // VT.ResizeTo(nxn,nxn);
      // VT = V;
      // VT = VT.Transpose(VT);

      if ( DoAlphasSingleShift ) { // a dedicated 'shift' for first parameter
	 // vector with alpha_s shift
	 V1.ResizeTo(nxn);
	 for ( int ix = 0 ; ix<nxn ; ix++ ) 
	    V1[ix]=cov[0][ix]/sqrt(cov[0][0]);

	 // covariance matrix minus cov of alpha_s shift
	 covm1.ResizeTo(nxn-1,nxn-1);
	 for ( int ix = 0 ; ix<nxn-1 ; ix++ ) {
	    for ( int iy = 0 ; iy<nxn-1 ; iy++ ) {
	       covm1[ix][iy] = cov[ix+1][iy+1] - V1[ix+1]*V1[iy+1];
	    }
	 }
	 TMatrixDSymEigen eigenm1(covm1);
	 eigVm1.ResizeTo(nxn-1);
	 eigVm1 = eigenm1.GetEigenValues();
	 info["execute"]<<"Printing Eigenvalues (n-1): "<<endl;
	 eigVm1.Print();
	 
	 Vm1.ResizeTo(nxn-1,nxn-1);
	 Vm1 = eigenm1.GetEigenVectors();
      }

      // --- debug
      // cout<<"EigenVectors:"<<endl;
      // V.Print();
      // EVmat.ResizeTo(nxn,nxn);
      // for ( int x = 0 ; x<nxn ; x++ ) EVmat(x,x) = eigen.GetEigenValues()[x];
      // Sym = VT*cov*V;
      // for ( int x = 0 ; x<nxn ; x++ ) {
      // 	 for ( int y = 0 ; y<nxn ; y++ ) {
      // 	    // remove numerical inacuracies. all off-diagonal are by construction zero.
      // 	    if ( x!=y ) Sym(x,y)=0;
      // 	 }
      // }
      // TMatrixD orig = V*Sym*VT;
      // TMatrixD Vinv = AlposTools::InvertLU(V);
      // TMatrixD Vinv2(V);
      // Vinv2.Invert();
      // cout<<"Symmetric Matrix:"<<endl;
      // Sym.Print();
      // cout<<"Orig:"<<endl;
      // orig.Print();
      // ---

      for ( int i = 0 ; i<(nPDF-1)/2 ; i++ ) {
	 fitparamnames.push_back(fitresult->ParName(i));
	 initialparams.push_back(fitresult->GetParams()[i]);
	 if ( PAR_ANY(fitparamnames.back()) != initialparams.back() ) {
	    error["execute"]<<"It is assumed, that the Alpos::TheorySet as input to this task is unchanged w.r.t. the fit results of the previous fit."<<endl;
	    error["execute"]<<"ParName: "<<fitparamnames.back()<<"\tAlposPar="<<PAR_ANY(fitparamnames.back())<<"\tFitResult="<<initialparams.back()<<endl;
	    exit(3);
	 }
      }
   }

   // -------------------------------------------------------
   // ---  write LHAPDF info file
   double XMin = 1.e-6;
   double XMax = 1;
   double QMin = 1.3;//4.472140;
   double QMax = 14142.1;
   double MZ   = 91.1876;
   double MC   = 1.3;
   double MB   = 4.5;
   double MT   = 173.;
   double ORD  = 2.;
   int nX = 120;
   int nQ = 120;
   if ( EXIST_NS(LHAPDF.XMin,NS())) XMin = DOUBLE_NS(LHAPDF.XMin,NS());
   //if ( EXIST_NS(LHAPDF.XMax,NS())) XMax = DOUBLE_NS(LHAPDF.XMax,NS());
   if ( EXIST_NS(LHAPDF.QMin,NS())) QMin = DOUBLE_NS(LHAPDF.QMin,NS());
   if ( EXIST_NS(LHAPDF.QMax,NS())) QMax = DOUBLE_NS(LHAPDF.QMax,NS());
   if ( EXIST_NS(LHAPDF.MZ,NS()))   MZ   = DOUBLE_NS(LHAPDF.MZ,NS());
   if ( EXIST_NS(LHAPDF.nX,NS()))   nX   = DOUBLE_NS(LHAPDF.nX,NS());
   if ( EXIST_NS(LHAPDF.nQ,NS()))   nQ   = DOUBLE_NS(LHAPDF.nQ,NS());
   if ( EXIST_NS(LHAPDF.MC,NS()))   MC   = DOUBLE_NS(LHAPDF.MC,NS());
   if ( EXIST_NS(LHAPDF.MB,NS()))   MB   = DOUBLE_NS(LHAPDF.MB,NS());
   if ( EXIST_NS(LHAPDF.MT,NS()))   MT   = DOUBLE_NS(LHAPDF.MT,NS());
   if ( EXIST_NS(LHAPDF.ORD,NS()))  ORD  = DOUBLE_NS(LHAPDF.ORD,NS());
   if ( WriteLHAPDF!="" ) {
      info["Execute"]<<"Writing LHAPDF grid."<<endl;
      info["Execute"]<<"LHAPDF grid: x-nodes: "<<nX<<" ["<<XMin<<","<<XMax<<"]."<<endl;
      info["Execute"]<<"LHAPDF grid: q-nodes: "<<nQ<<" ["<<QMin<<","<<QMax<<"]."<<endl;
      info["Execute"]<<"LHAPDF MZ: "<<MZ<<endl;
      info["Execute"]<<"LHAPDF HQ masses: "<<MC<<", "<<MB<<", "<<MT<<endl;
      info["Execute"]<<"LHAPDF order: "<<ORD<<endl;
      string outputdir = Alpos::Current()->Settings()->outputdir;
      gSystem->mkdir((outputdir+"/"+WriteLHAPDF).c_str(),true);
      ofstream infofile((outputdir+"/"+WriteLHAPDF+"/"+WriteLHAPDF+".info").c_str(),std::ofstream::out);
      infofile<<"SetDesc: \"Add your description here\""<<endl;
      infofile<<"SetIndex: 00000"<<endl;
      infofile<<"Authors: Alpos auto-generated file"<<endl;
      infofile<<"Reference: [Add reference here]"<<endl;
      infofile<<"Format: lhagrid1"<<endl;
      infofile<<"DataVersion: 1"<<endl;
      infofile<<"NumMembers: "<<nPDF <<endl;
      infofile<<"Particle: 2212"<<endl;
      infofile<<"Flavors: [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]"<<endl;
      infofile<<"OrderQCD: "<<ORD<<endl;
      infofile<<"FlavorScheme: variable"<<endl;
      infofile<<"NumFlavors: 6"<<endl;
      infofile<<"ErrorType: hessian"<<endl;
      infofile<<"XMin: "<<XMin<<endl;
      infofile<<"XMax: "<<XMax<<endl;
      infofile<<"QMin: "<<QMin<<endl;
      infofile<<"QMax: "<<QMax<<endl;
      infofile<<"MZ: "<<MZ<<endl;
      infofile.close();
   }


   // ----
   map<string,vector<double> > pdfdef = GetPDFdef();

   //if ( nPDF > 1 && !( LHAType || FitType )  ) warn["Execute"]<<"only loop over LHAPDF  members or results of a PDF fit impemented."<<endl;
   map<string,map<double,map<double, vector<double> > > > AllValues; // <parname,q2,xp,val>
   vector<TGraph*> gFitParams,gShiftParams;
   vector<TH1D*>   hFitParams,hShiftParams;
   vector<double> TmpLHAPDFset0;
   vector<double> TmpLHAPDFAlphaS_Val0;
   double asmz0=0;
   for ( int iMem = 0 ;iMem<nPDF ; iMem++ ) { // iMem-loop is slowest!
      info["execute"]<<"Evaluating PDFset "<<iMem<<"/"<<nPDF-1<<endl;
      if ( LHAType ) { 
	 SET_ANY(pdffunc+".PDFSet",iMem,0); 
      }
      else if ( FitType ) {
	 int nEig = iMem!=0 ? (iMem-1)/2 : -1;
	 if ( nEig >= 0 ) {
	    //SET_ANY(pdffunc+".PDFSet",iMem,0); 
	    int nxn = (nPDF-1)/2;
	    bool IsUp = (iMem%2 == 1);
	    /*
	    TMatrixD CovEV(nxn,nxn);
	    CovEV(nEig,nEig) = eigV[nEig];
	    // cout<<"Matrix 'CovEV': "<<endl;
	    // CovEV.Print();
	    CovEV = V*CovEV*VT;
	    cout<<"Matrix 'CovEV' in FitParam space: "<<endl;
	    CovEV.Print();
	    //cout<<"Up/Dn; IsUp="<<IsUp<<endl;
	    for ( int i = 0 ; i <nxn ;i++ ) {
	       double ishift = sqrt(CovEV(i,i));
	       if ( !IsUp ) ishift*=-1;
	       cout<<"Parameter:  "<<fitparamnames[i]<<"\tshift="<<ishift<<endl;
	       SET_ANY(fitparamnames[i],initialparams[i]+ishift,0); 
	    }
	    */
	    //TVectorD shiftV1(nxn);
	    gFitParams.  push_back(new TGraph());
	    gShiftParams.push_back(new TGraph());
	    gFitParams.  back()->SetName("GraphFitParametersShifted");
	    gShiftParams.back()->SetName("GraphFitParametersShift");
	    //hFitParams.push_back(new TH1D(Form("HistFitParametersShifted_%d",iMem),";parameter",nxn,-0.5,nxn-0.5));
	    hFitParams.  push_back(new TH1D("HistFitParametersShifted",";parameter",nxn,-0.5,nxn-0.5));
	    hShiftParams.push_back(new TH1D("HistFitParametersShift",";parameter",nxn,-0.5,nxn-0.5));
	    for ( int i = 0 ; i <nxn ;i++ ) {
	       //shiftV1[i] = sqrt(eigV[nEig])*V(i,nEig);
	       //double ishift = sqrt(eigV[nEig])*V(i,nEig); // scale eigenvector matrix with shifts
	       //ishift*=0.1;// scale by factor of 10. Undo this when plotting results!

	       double ishift= 0;
	       if ( DoAlphasSingleShift ) {
		  if ( nEig==0 ) ishift = V1[i]; // scale eigenvector matrix with shifts
		  else if ( i==0 ) ishift = 0;
	       	  else ishift = sqrt(eigVm1[nEig-1])*Vm1(i-1,nEig-1);
	       }
	       else {
		  ishift = sqrt(eigV[nEig])*V(i,nEig); // scale eigenvector matrix with shifts
	       }

	       if ( !IsUp ) ishift*=-1;
	       //cout<<"Parameter:  "<<fitparamnames[i]<<"\tshift="<<ishift<<endl;
	       SET_ANY(fitparamnames[i],initialparams[i]+ishift,0); 
	       gFitParams.  back()->SetPoint(i,i,initialparams[i]+ishift);
	       gShiftParams.back()->SetPoint(i,i,ishift);
	       hFitParams.  back()->SetBinContent(i+1,initialparams[i]+ishift);
	       hShiftParams.back()->SetBinContent(i+1,ishift);
	       hFitParams.  back()->GetXaxis()->SetBinLabel(i+1,fitparamnames[i].c_str());
	       hShiftParams.back()->GetXaxis()->SetBinLabel(i+1,fitparamnames[i].c_str());
	       //cout<<"vShift1: "<<Vshifts1(i,nEig)<<", "<<Vshifts1(nEig,i)<<"\tVShifts2="<<Vshifts2(i,nEig)<<", "<<Vshifts2(nEig,i)<<endl;
	    }	    
	    // cout<<"shiftV1"<<endl;
	    // shiftV1.Print();
	 }
      }
      PAR_ANY(pdffunc); // update PDF
      if ( asfunc != "" ) PAR_ANY(asfunc); // update Alphas

      // --- Q2/mu_f loop
      for ( auto q2 : q2val ) {
	 TString dirname = Form("Q2_%.1f",q2);
	 if ( !taskdir->GetDirectory(dirname) ) taskdir->mkdir(dirname);
	 TDirectory* q2dir = taskdir->GetDirectory(dirname);
	 //taskdir->GetDirectory(dirname)->mkdir(Form("PDF_%d",iMem))->cd();
	 if ( !q2dir->GetDirectory(Form("PDF_%d",iMem))) q2dir->mkdir(Form("PDF_%d",iMem));
	 q2dir->cd(Form("PDF_%d",iMem));

	 // --- store shift of parameters in each directory
	 if ( gFitParams.size() ) {
	    gFitParams.  back()->Write();
	    gShiftParams.back()->Write();
	    hFitParams.  back()->Write();
	    hShiftParams.back()->Write();
	 }

	 // --- store all PDFs
	 for ( auto ipdf : pdfdef ) {
	    TGraph gPDF;
	    gPDF.SetName(ipdf.first.c_str());
	    for ( auto xp : xval ) {
	       //SET_ANY("LHAPDF.xp",xp,0);
	       //vector<double> xfx = VALUES_ANY("LHAPDF");
	       vector<double> xp_muf={xp,sqrt(q2)};
	       vector<double> xfx = QUICK_ANY(pdffunc,xp_muf);
	       double xf = 0;
	       for ( unsigned int jf=0 ; jf< ipdf.second.size() ; jf++ ) {
		  xf += ipdf.second[jf]*xfx[jf];
	       }
	       gPDF.SetPoint(gPDF.GetN(),xp,xf);
	       AllValues[ipdf.first][q2][xp].push_back(xf);
	    }
	    gPDF.Write();
	 }

	 // --- store alpha_s(mZ) and alpha_s(mu)
	 if ( asfunc!= "" ) {
	    TGraph gAsMz,gAsMu;
	    gAsMz.SetName("Alpha_s(Mz)");
	    gAsMu.SetName("Alpha_s(Q)");
	    double asmu = QUICK_ANY(asfunc,vector<double>{sqrt(q2)})[0];
	    double asmz = QUICK_ANY(asfunc,vector<double>{91.1876})[0];
	    gAsMz.SetPoint(0,1,asmz);
	    gAsMu.SetPoint(0,1,asmu);
	    gAsMz.Write();
	    gAsMu.Write();
	 }
	 gDirectory->Write();
      }

      // ---WeightNNPDF loop
      if (WeightNNPDF) {
	 info["Execute"]<<"Calculate chi2 for PDF set "<<iMem<<"/"<<nPDF-1<<endl;
	 if ( !EXIST_NS(Chisq,NS()) ) { 
	       error["exexute"]<<"Chisq definition needed. Please specify (e.g. 'chisq LogNormal')"<<endl;
	       exit(1);
	    }
	 vector<string> parnames;
	 string chisqdef = STRING_NS(Chisq,NS());
	 AChisqBase* chi = AFactory::ChisqFactory(chisqdef,parnames,
						  (AData*)TheoryHandler::Handler()->GetFuncD("SuperData"),
						  TheoryHandler::Handler()->GetFuncD("SuperTheory"));
	 double chisq = chi->DoEval(NULL);
	 //std::map<std::string, double> nui = chi->GetNuisanceParameters();

	 // --- calculate weight
	 // https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Maximum-likelihood_estimation_for_the_multivariate_normal_distribution
	 // http://nnpdf.mi.infn.it/research/reweighting/
	 // arxiv.org/abs/1108.1758
	 //
	 // https://indico.desy.de/indico/event/9388/session/2/contribution/33/material/slides/0.pdf 
	 // Bayesian Reweighting
	 // First implementation suggested by Giele and Keller, who thought of it as 
	 // method for producing new PDF fits
	 // 			[W. Giele & S. Keller, hep-ph/9803393]
 	 //
	 // Reformulated in the context of the NNPDF fits
	 // (based on Monte Carlo methodology for uncertainties estimation) and applied for the first time to 
	 // include data in a global PDF fit (NNPDF2.2)
	 // 			[R.D. Ball et al., arXiv:1012.0836]
	 // 		        [R.D. Ball et al., arXiv:1108.1758]
	 // !
	 // Recently extended by G. Watt and R. Thorne to PDF based on the Hessian method
	 //  for estimation of uncertainties
	 // 		   [G. Watt & R.S. Thorne, arXiv:1205.4024]
	 // 		   [LHCb, De Lorenzi et. al, arXiv:1011.4260]
 	 // 
	 double constfactor = 1.;// constant factor (2pi)^{-np/2} PROD{sqrt(1/det(V))}
	 //double wgt = exp(-0.5*chisq) * constfactor;
	 double N = nPDF;
	 double n = TheoryHandler::Handler()->GetFuncD("SuperTheory")->GetValues().size();
	 //double wgt = pow(chisq,0.5*(n-1)) * exp(-0.5*chisq); // the first totally 'explodes', while the second one vanishes!
	 double wgt = exp ( (n-1)/2. * log(chisq) - chisq/2.);
	 double lwgt = (n-1)/2. * log(chisq) - chisq/2.;
	 info["execute"]<<"Calculating WeightNNPDF. iMem: "<<iMem<<"\t\tchi2="<<chisq<<"\texp(-0.5chi2)="<<exp(-0.5*chisq)<<"\tPS="<<pow(chisq,0.5*(n-1))<<"\twgt="<<wgt<<"\tlog(wgt)="<<lwgt<<endl;
	 vWgtNNPDF.push_back(wgt);
	 vlWgtNNPDF.push_back(lwgt);
	 vChi2NNPDF.push_back(chisq);
      }

      // LHAPDF output file
      if ( WriteLHAPDF!="" ) {
	 string outputdir = Alpos::Current()->Settings()->outputdir;
	 ofstream gridfile((outputdir+"/"+WriteLHAPDF+"/"+WriteLHAPDF+"_"+Form("%04d",iMem)+".dat").c_str(),std::ofstream::out);
	 //gridfile

	 if ( iMem==0 ) gridfile<<"PdfType: central"<<endl;
	 else           gridfile<<"PdfType: error"<<endl;
	 gridfile<<"Format: lhagrid1"<<endl;
	 gridfile<<"XMin: "<<XMin<<endl;
	 gridfile<<"XMax: "<<XMax<<endl;
	 gridfile<<"QMin: "<<QMin<<endl;
	 gridfile<<"QMax: "<<QMax<<endl;
	 gridfile<<"MZ: "<<MZ<<endl;
	 //
	 gridfile<<"MUp: 0"<<endl;
	 gridfile<<"MDown: 0"<<endl;
	 gridfile<<"MStrange: 0"<<endl;
	 gridfile<<"MCharm: "<<MC<<endl;
	 gridfile<<"MBottom: "<<MB<<endl;
	 gridfile<<"MTop: "<<MT<<endl;
	 double asmz = QUICK_ANY(asfunc,vector<double>{MZ})[0];
	 double ErrorDef =  EXIST_NS(LHAPDF.ErrorDef,NS()) ?  DOUBLE_NS(LHAPDF.ErrorDef,NS()) : 1.;
	 if (asmz0 == 0 ) asmz0=asmz;
	 else if ( ErrorDef!=1 ) {
	    double ScaleErr = 1./sqrt(ErrorDef);
	    asmz = asmz0+(asmz-asmz0)*ScaleErr;
	 }
	 double asorg=QUICK_ANY(asfunc,vector<double>{MZ})[0];
	 cout<<"asmz0="<<asmz0
	     <<"\tasmz="<<asmz
	     <<"\tasmzorg="<<asorg
	     <<"\t(as^-as0)="<<asmz-asmz0
	     <<"\t(asor-as0)="<<asorg-asmz0
	     <<"\t(as^/as0-1)="<<asmz/asmz0-1
	     <<"\t(asor/as0-1)="<<asorg/asmz0-1
	     <<endl;

	 gridfile<<"AlphaS_MZ: "<<asmz<<endl;
	 gridfile<<"AlphaS_OrderQCD: "<<ORD<<endl;
	 //
	 gridfile<<"AlphaS_Type: ipol"<<endl;
	
	 vector<double> xpt1 =  this->GetLogNodes(XMin,0.1,nX/2);
	 vector<double> xpt2 =  this->GetLogNodes(0.1,XMax,nX/2);
	 xpt1.resize(xpt1.size()-1);// remove last node
	 for ( auto xx : xpt2 ) xpt1.push_back(xx);
	 vector<double> qpt =  this->GetLogNodes(QMin,QMax,nQ);
	 vector<int> PDFid{-6,-5, -4, -3, -2, -1,21, 1, 2, 3, 4, 5, 6};
	 // write alpha_s and PDF grid:
	 this->WriteAlphasGrid(gridfile,qpt,&TmpLHAPDFAlphaS_Val0);
	 this->WritePDFGrid(gridfile,PDFid,xpt1,qpt,&TmpLHAPDFset0);
	 gridfile.close();
      }

   }

   if (WeightNNPDF) {
      info["Execute"]<<"Calculate weights for all PDF sets."<<endl;
      
      double N = nPDF;
      double n = TheoryHandler::Handler()->GetFuncD("SuperTheory")->GetValues().size();
      double sum = 0;
      for ( const auto& wgt : vWgtNNPDF ) sum+=wgt;
      cout<<"sum="<<sum<<endl;
      for ( auto& wgt : vWgtNNPDF ) {
       	 wgt = wgt/sum*N;
       	 cout<<"wgt="<<wgt<<endl;
       }

      // eq 38,39 in 1012.0836
      double ekavg = 1;
      for ( const auto& lwgt : vlWgtNNPDF ) ekavg = lwgt/N;
      double sumek = 0;
      for ( const auto& lwgt : vlWgtNNPDF )  sumek += exp(lwgt-ekavg);
      double norm = N/sumek;
      cout<<"norm="<<norm<<"\tsum(e_k)="<<sumek<<"\t<e_k>="<<ekavg<<endl;
      for ( auto& lwgt : vlWgtNNPDF ) {
      	 cout<<"lwgt="<<lwgt<<"\twgt="<<exp(lwgt)<<endl;
      	 double newwgt = norm*exp(lwgt - ekavg);
      	 cout<<"wgt="<<newwgt<<endl;
	 //lwgt=newwgt;
      }


      // eq 38,39 in 1012.0836
      // and exp(norm)...
      for ( const auto lwgt : vlWgtNNPDF )  {
	 double sumExpA = 0;
	 for ( const auto ai : vlWgtNNPDF )  {
	    sumExpA += exp(ai-lwgt);
	 }
	 cout<<"wgt = "<<1./sumExpA*N<<endl;
      }

      for ( unsigned int i = 0 ; i<vWgtNNPDF.size() ; i++ ) {
	 double sumExpA = 0;
	 for ( const auto ai : vlWgtNNPDF )   sumExpA += exp(ai-vlWgtNNPDF[i]);
	 //cout<<"wgt = "<<1./sumExpA*N<<endl;
	 vWgtNNPDF[i]=1./sumExpA*N;
      }


      // shannon entropy
      double Neff= 0;
      for ( auto& wgt : vWgtNNPDF ) Neff += wgt*log(N/wgt);
      Neff/=N;
      Neff = exp(Neff);
      info["Execute"]<<"Shannon entropy is: "<<Neff<<", for N="<<N<<endl;
      
   }


   // for ( auto g : gFitParams )   delete g;
   // for ( auto g : gShiftParams ) delete g;
   // for ( auto h : hFitParams )   delete h;
   // for ( auto h : hShiftParams ) delete h;

   
   // TGraphs with error bands
   if ( nPDF>1 && ( LHAError || FitType) ) {
      LHAPDF::PDFSet* lhapdfset = LHAType ?
	 ((ALhapdf6*)TheoryHandler::Handler()->GetFuncD(pdffunc))->GetPDFSet() : 
	 NULL;
      for ( auto q2 : q2val ) {
	 TString dirname = Form("Q2_%.1f",q2);
	 TDirectory* q2dir = taskdir->GetDirectory(dirname);
	 // if ( !q2dir->GetDirectory("PDF_Errors") ) q2dir->mkdir("PDF_Errors")->cd();
	 // q2dir->cd("PDF_Errors");
	 q2dir->mkdir("PDF_ErrorsSymm");
	 q2dir->mkdir("PDF_ErrorsAsym");

	 for ( auto ipdf : pdfdef ) {
	    TGraphAsymmErrors gPDFSymm,gPDFAsym;
	    gPDFSymm.SetName(ipdf.first.c_str());
	    gPDFAsym.SetName(ipdf.first.c_str());
	    const map<double,map<double, vector<double> > >& iPdfValues = AllValues.at(ipdf.first);
	    for ( auto xp : xval ) {
	       int nP = gPDFSymm.GetN();
	       //const vector<double>& values = AllValues[ipdf.first][q2][xp]; // this is a bit slow
	       const vector<double>& values = iPdfValues.at(q2).at(xp);
	       double errdn=0,errup=0,errsym=0;
	       if (LHAError) {
		  //cout<<"ipdf.first="<<ipdf.first<<"\tiMem="<<iMem<<"\tq2="<<q2<<"\txp="<<xp<<"\tsize="<<values.size()<<endl;
		  LHAPDF::PDFUncertainty PDFUnc = lhapdfset->uncertainty(values);
		  //cout<<"nP="<<nP<<"\txp="<<xp<<"\tcentral="<<PDFUnc.central<<"\tup="<<PDFUnc.errminus<<"\tdn="<<PDFUnc.errplus<<endl;
		  //cent = PDFUnc.central;
		  errdn = PDFUnc.errminus;
		  errup = PDFUnc.errplus;
	       }
	       else { 
		  // calculate asymmetric 'hessian' uncertainties
		  // see e.g.:
		  // hep-ph/0101032
		  // arXiv:1101.0536
		  // code adapted from: LHAPDF::src/PDFSet.cc 
		  int nxn = (nPDF-1)/2;
		  for (int ieig = 0; ieig < nxn; ieig++) {
		     errup  += AlposTools::sq(max(max(values[2*ieig+1]-values[0],values[2*ieig+2]-values[0]), 0.));
		     errdn  += AlposTools::sq(max(max(values[0]-values[2*ieig+1],values[0]-values[2*ieig+2]), 0.));
		     errsym += AlposTools::sq(values[2*ieig+1]-values[2*ieig+2]);
		  }
		  errsym = 0.5*sqrt(errsym);
		  errup  = sqrt(errup);
		  errdn  = sqrt(errdn);
	       }
	       double cent = values[0];
	       gPDFSymm.SetPoint(nP,xp,cent);
	       gPDFAsym.SetPoint(nP,xp,cent);
	       gPDFSymm.SetPointError(nP,0,0,errsym,errsym);
	       gPDFAsym.SetPointError(nP,0,0,errdn,errup);
	    }
	    // save to file
	    q2dir->cd("PDF_ErrorsAsym");
	    gPDFAsym.Write();
	    if ( FitType ) {
	       q2dir->cd("PDF_ErrorsSymm");
	       gPDFSymm.Write();
	    }
	 }
	 taskdir->Write();
      }
      info["Execute"]<<"TGraphs with errors calculated."<<endl;
   }
   else {
      info["Execute"]<<"No error TGraphs are drawn."<<endl;
   }

   // reset
   if ( Mem0>=0 ) SET_ANY(pdffunc+".PDFSet",Mem0,0);
   if ( FitType ) {
      int nxn = (nPDF-1)/2;
      for ( int i = 0 ; i <nxn ;i++ ) 
	 SET_ANY(fitparamnames[i],initialparams[i],0); 
   }

   // close rootfile
   file->Write();
   //if ( file != Alpos::Current()->Settings()->rootoutput) file->Close();
   cout<<endl;
   info["Execute"]<<"File "<<file->GetName()<<" written to disk.\n"<<endl;
   // TheoryHandler::Handler()->SetTheorySet(this->GetResult()->GetInitialTheorySet());
   // info["Execute"]<<"Initial parameters resetted."<<endl;

   return true;
}


//____________________________________________________________________________________ //
map<string,vector<double> > ASavePDFTGraph::GetPDFdef() const {
   //! return definition for all PDF linear combinations
   map<string,vector<double> > pdfdef;
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   pdfdef[   "tb"] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "bb"] = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "cb"] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "sb"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "ub"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "db"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "g"] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "d"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ,0, 0 ,0 };
   pdfdef[    "u"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[    "s"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[    "c"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   pdfdef[    "b"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 1 ,0 };
   pdfdef[    "t"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,1 };
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   pdfdef["g*0.05"] = { 0, 0, 0, 0, 0, 0,0.05, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "dv"] = { 0, 0, 0, 0, 0,-1, 0, 1, 0, 0 ,0, 0 ,0 };
   pdfdef[   "uv"] = { 0, 0, 0, 0,-1, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[   "s+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "s-"] = { 0, 0, 0,-1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c-"] = { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   pdfdef[ "Ubar"] = { 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[ "Dbar"] = { 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "U"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,1, 0 ,1 };
   pdfdef[    "D"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 ,0, 1 ,0 };
   pdfdef[   "U+"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,1, 0 ,1 };
   pdfdef[   "D+"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 ,0, 1 ,0 };
   //
   pdfdef["gluon"] =   {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}; // gluon
   pdfdef["SIGMA"] =   {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}; // SIGMA = \sum_q q^+     
   pdfdef["VALENCE"]=  { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}; // VALENCE = \sum_q q^-   
   pdfdef[   "T3"] =   {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}; // T3  = u^+ - d^+        
   pdfdef[   "V3"] =   {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}; // V3  = u^- - d^-        
   pdfdef[   "T8"] =   {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}; // T8  = u^+ + d^+ - 2 s^+
   pdfdef[   "V8"] =   {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}; // V8  = u^- + d^- - 2 s^-
   pdfdef[  "T15"] =   {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}; // T15 = u^+ + d^+ + s^+ - 3 c^+            
   pdfdef[  "V15"] =   {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}; // V15 = u^- + d^- + s^- - 3 c^-            
   pdfdef[  "T24"] =   {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}; // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+      
   pdfdef[  "V24"] =   {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}; // V24 = u^- + d^- + s^- + c^- - 4 b^-      
   pdfdef[  "T35"] =   { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}; // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
   pdfdef[  "V35"] =   {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,};  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-

   return pdfdef;

}


//____________________________________________________________________________________ //
std::vector<double> ASavePDFTGraph::GetLogNodes(double min, double max, int npts) const {
   std::vector<double> ret(npts+1);
   double lmin = log(min);
   double lmax = log(max);
   for ( int i=0 ; i<npts+1 ; i++ ) {
      ret[i] = exp(lmin+(lmax-lmin)/npts*i);
   }
   return ret;
}
	
//____________________________________________________________________________________ //
void ASavePDFTGraph::WriteAlphasGrid(ostream& strm, const vector<double>& qpt, vector<double>* grid0){
   //! write alpha_s grid to stream
   if ( !EXIST_NS(Alphas,NS()) ) {
      error["WriteAlphasGrid"]<<"Alpha_s function needs to be given to ASavePDFTGraph, when writing LHAPDF file."<<endl;
      exit(5);
   }

   string asfunc = STRING_NS(Alphas,NS());

   bool StoreGrid0 = grid0!=NULL && grid0->size()==0;
   bool DoScaleErr = grid0!=NULL && grid0->size()!=0 && EXIST_NS(LHAPDF.ErrorDef,NS());
   double ScaleErr = 1.;
   if ( DoScaleErr ) {
      double ErrorDef = DOUBLE_NS(LHAPDF.ErrorDef,NS());
      info["WritePDFGrid"]<<"Scaling alpha_s grid according to "<<NS()<<"::LHAPDF.ErrorDef="<<ErrorDef<<endl;
      ScaleErr=1./sqrt(ErrorDef);
   }

   //AlphaS_Qs: [ .. , .. , .. ]
   strm<<"AlphaS_Qs: [ ";
   for ( auto iq : qpt )  {
      strm<<iq;
      if ( iq!=qpt.back() ) strm<<", ";
      else strm<<" ]"<<endl;
   }
   //AlphaS_Vals: [ .., .., .. ]
   strm<<"AlphaS_Vals: [";	
   int cc=0;
   grid0->reserve(qpt.size());
   for ( auto iq : qpt )  {
      double asmu = QUICK_ANY(asfunc,vector<double>{iq})[0];
      if ( StoreGrid0 ) 
	 grid0->push_back(asmu);
      else if ( DoScaleErr ) {
	 asmu = grid0->at(cc)+(asmu-grid0->at(cc))*ScaleErr;
	 cc++;
      }
      strm<<asmu;
      if ( iq!=qpt.back() ) strm<<", ";
      else strm<<" ]"<<endl;
   }
   strm<<"---"<<endl; /// delimiter
}

//____________________________________________________________________________________ //
void ASavePDFTGraph::WritePDFGrid(ostream& strm, const vector<int>& PDFid, const vector<double>& xpt, const vector<double>& qpt, vector<double>* grid0){
   //! write PDF grid to stream
   // x0 x1 x2 ... xN   /// e.g. 1.E-7 ... 1.0
   // q0 q1 q2 ...  // as q0=QMin qN=QMax as above 
   // -5 -4 -3 -2 -1 1 2 3 4 5 21
   // x0,q0 ....
   // x0,q1 ....

   info["WritePDFGrid"]<<"Writing PDF grid with "<<xpt.size()<<" x nodes and "<<qpt.size()<<" q nodes."<<endl;
   // write x nodes
   for ( auto ix : xpt ) strm<<" "<<ix;
   strm<<endl;
   // write Q nodes
   for ( auto iq : qpt ) strm<<" "<<iq;
   strm<<endl;
   // write PDF ids
   for ( auto ip : PDFid ) strm<<" "<<ip;
   strm<<endl;
   
   string pdffunc = STRING_NS(PDF,NS());

   bool StoreGrid0 = grid0!=NULL && grid0->size()==0;
   bool DoScaleErr = grid0!=NULL && grid0->size()!=0 && EXIST_NS(LHAPDF.ErrorDef,NS());
   double  ScaleErr = 1;
   if ( DoScaleErr ) {
      double ErrorDef = DOUBLE_NS(LHAPDF.ErrorDef,NS());
      info["WritePDFGrid"]<<"Scaling grid according to "<<NS()<<"::LHAPDF.ErrorDef="<<ErrorDef<<endl;
      ScaleErr=1./sqrt(ErrorDef);
   }

   int cc=0;
   grid0->reserve(xpt.size()*qpt.size());
   for ( auto ix : xpt ) {
      for ( auto iq : qpt ) {
	 vector<double> xp_muf{ix,iq};
	 vector<double> xfx = QUICK_ANY(pdffunc,xp_muf);
	 for ( auto ip : PDFid ) {
	    int ipdf = ip;
	    if ( ipdf==21 ) ipdf=0; // adjust convention
	    double pdf = xfx[ipdf+6];
	    if ( StoreGrid0 ) 
	       grid0->push_back(pdf);
	    else if ( DoScaleErr ) {
	       pdf = grid0->at(cc)+(pdf-grid0->at(cc))*ScaleErr;
	       cc++;
	    }
	    //if ( pdf<0 ) pdf=0; // force positiveness
	    if ( fabs(pdf) < 1.e-12 ) pdf=0;// drop negligible values
	    strm<<" "<<pdf;
	 }
	 strm<<endl;
      }
   }
   strm<<"---"<<endl;
}
