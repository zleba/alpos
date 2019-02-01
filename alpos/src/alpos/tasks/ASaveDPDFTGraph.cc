#include "alpos/tasks/ASaveDPDFTGraph.h"
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
#include "alpos/ASubsetFunction.h"
#include <numeric>

/** 
    ASaveDPDFTGraph

    * Save TGraphs of PDFs

    */

using namespace std;

const string ASaveDPDFTGraph::fTaskType = "SaveDPDFTGraph";

//____________________________________________________________________________________ //
ASaveDPDFTGraph::ASaveDPDFTGraph(const string& aname ) : ATask(aname) {
   //! constructor
   // You may use the 'speaker' for variuos verbosity levels:
   SetClassName("ASaveDPDFTGraph");
   //! Important: create always new result-object here!
   fResult = new ASaveDPDFTGraphResult(aname,GetTaskType());
}


//____________________________________________________________________________________ //
ASaveDPDFTGraph::~ASaveDPDFTGraph(){
   //! destructor.
   //! Do not delete the AResult object!
}


//____________________________________________________________________________________ //
bool ASaveDPDFTGraph::Init(){
   // Init() and Execute() are directly executed after each other. So to say: there is no difference!
   debug["Init"]<<endl;
   return true;
}


//____________________________________________________________________________________ //
bool ASaveDPDFTGraph::Execute(){
   debug["Execute"]<<endl;
   

   // -------------------------------------------------------------------------- //
   // init output
   // -------------------------------------------------------------------------- //
   //TFile* file = TFile::Open(CHAR_NS(RootFilename,NS()),"RECREATE");
   TDirectory* file = Alpos::Current()->Settings()->rootoutput;
   if ( !file ) {
      info["Execute"]<<"No output directory or file specified in steering. Skipping task."<<endl;
      return true;
   }
   TDirectory* taskdir = file->mkdir(GetTaskName().c_str());
   taskdir->cd();
   // the PDF-function
   string dpdffunc = STRING_NS(DPDF,NS());
   PAR_ANY(dpdffunc); // update
   info["Execute"]<<"Evaluating DPDF from theory parameter: "<<dpdffunc<<" which is of type '"<<
      TheoryHandler::Handler()->GetFuncD(dpdffunc)->GetFunctionName()<<"'."<<endl;
   // alpha_s evolution (as it is used for PDF evolution)

   // ascii/LHAPDF
   /*
     string asfunc = STRING_NS(Alphas,NS());
     if ( asfunc=="" ) {
     warn["Execute"]<<"Alphas evolution not specified. Please provde Alphas-evolution function with key 'Alphas'"<<endl;
     }
     else {
     PAR_ANY(asfunc); // update
     info["Execute"]<<"Evaluating Alphas evolution from theory parameter: "<<asfunc<<" which is of type '"<<
     TheoryHandler::Handler()->GetFuncD(asfunc)->GetFunctionName()<<"'."<<endl;
     }
     string WriteAscii = STRING_NS(LHAPDFoutput,NS());
     if ( WriteAscii!="" && asfunc=="" ) {
     error["Execute"]<<"Ascii/LHAPDF output specified, but no Alphas function given. Cannot write LHAPDF grids."<<endl;
     WriteLHAPDF="";
     }
     else if ( WriteLHAPDF!="" ) info["Execute"]<<"LHAPDF compatible files will be written into directory: "<<WriteLHAPDF<<endl;
   */

   
   // -------------------------------------------------------------------------- //
   // x and q2 spacing
   vector<double> q2val   = DOUBLE_ARR_NS(Q2Values,NS());
   vector<double> xpomval = DOUBLE_ARR_NS(xpomValues,NS());
   int nzp = INT_NS(nzpom,NS());
   double zmin = DOUBLE_NS(zmin,NS());
   double lxstep = (log(1)-log(zmin))/nzp;
   double lxlinstep = (1-zmin)/nzp;
   
   vector<double> zval = {zmin};
   for ( int ix = 0 ; ix < nzp-1 ; ix++ ) {
      double zNext = exp(log(zval.back())+lxstep);
      if ( zNext - zval.back() > lxlinstep ) break;
      zval.push_back(zNext);
   }
   while ( zval.back() < 1 ) { // lin spacing at high z
      zval.push_back(zval.back()+lxlinstep);
   }
   zval.back() = 1;
   // -------------------------------------------------------------------------- //

   // -------------------------------------------------------------------------- //
   // --- FitType
   int nPDF=1;
   bool FitType=false;
   int Mem0 = -1; // PDFSet in the beginning

   // linear algebra
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
      info["Execute"]<<"Found fit result. Now calculating eigenvectors."<<endl;
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
      //   for ( int y = 0 ; y<nxn ; y++ ) {
      //      // remove numerical inacuracies. all off-diagonal are by construction zero.
      //      if ( x!=y ) Sym(x,y)=0;
      //   }
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
   // -------------------------------------------------------------------------- //


   // -------------------------------------------------------------------------- //
   // --- cross sections
     

   // for ( const auto& id: TheoryHandler::Handler()->GetAllSubsetPairs() )
   //    crosssections[id.second.second->GetAlposName()][0] = id.second.second->GetValues();
    
   //   Print(id.second.first,id.second.second);
   // const auto& vv =  super.first->GetValues();
   // for ( double v : vv ) cout<<"\t"<<v;
   // cout<<endl;
   // // single datasets
   // for ( const auto& id : TheoryHandler::Handler()->GetDataTheoryPairs() )
   //   Print(id.second.first,id.second.second);
   // // subsets 
   // for ( const auto& is: TheoryHandler::Handler()->GetAllSubsetPairs() )
   //    Print(is.second.first,is.second.second);
   // -------------------------------------------------------------------------- //


   // -------------------------------------------------------------------------- //
   // ----
   map<string,vector<double> > pdfdef = GetDPDFdef();

   //if ( nPDF > 1 && !( LHAType || FitType )  ) warn["Execute"]<<"only loop over LHAPDF  members or results of a PDF fit impemented."<<endl;
   map<string,map<double,map<double, vector<double> > > > AllValPom; // <parname,q2,xp,val>
   map<string,map<double,map<double, vector<double> > > > AllValReg; // <parname,q2,xp,val>
   map<string,map<double,map<double,map<double, vector<double> > > > > AllValues; // <parname,xpom,q2,xp,val>
   map<string,map<int,vector<double> > > crosssections; // theoname, iMem, values<>
   vector<TGraph*> gFitParams,gShiftParams;
   vector<TH1D*>   hFitParams,hShiftParams;
   vector<double> TmpLHAPDFset0;
   vector<double> TmpLHAPDFAlphaS_Val0;
   double asmz0=0;
   info["Execute"]<<"Calculating DPDFs for "<<xpomval.size()<<" xpom values, and "<<q2val.size()<<" mu_f^2 values."<<endl;
   for ( int iMem = 0 ;iMem<nPDF ; iMem++ ) { // iMem-loop is slowest!
      info["execute"]<<"Evaluating PDFset "<<iMem<<"/"<<nPDF-1<<endl;
      if ( FitType ) {
         int nEig = iMem!=0 ? (iMem-1)/2 : -1;
         if ( nEig >= 0 ) {
            //SET_ANY(dpdffunc+".PDFSet",iMem,0); 
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
      PAR_ANY(dpdffunc); // update PDF
      //if ( asfunc != "" ) PAR_ANY(asfunc); // update Alphas

      // store cross section predictions 
      const auto& super = TheoryHandler::Handler()->GetSuperPair();
      crosssections[super.second->GetAlposName()][iMem] = super.second->GetValues();
      for ( const auto& id : TheoryHandler::Handler()->GetDataTheoryPairs() )
         crosssections[id.second.second->GetAlposName()][iMem] = id.second.second->GetValues();
      for ( const auto& id: TheoryHandler::Handler()->GetAllSubsetPairs() )
         crosssections[id.second.second->GetAlposName()][iMem] = id.second.second->GetValues();


      // --- Q2/mu_f loop
      for ( auto xp : xpomval ) {
         for ( auto q2 : q2val ) {
            TString dirname = Form("Q2_%g__xpom_%0.4f",q2,xp);
            if ( !taskdir->GetDirectory(dirname) ) taskdir->mkdir(dirname);
            TDirectory* q2dir = taskdir->GetDirectory(dirname);
            //taskdir->GetDirectory(dirname)->mkdir(Form("PDF_%d",iMem))->cd();
            if ( !q2dir->GetDirectory(Form("DPDF_%d",iMem))) q2dir->mkdir(Form("DPDF_%d",iMem));
            q2dir->cd(Form("DPDF_%d",iMem));

            // --- store shift of parameters in each directory
            if ( gFitParams.size() ) {
               gFitParams.  back()->Write();
               gShiftParams.back()->Write();
               hFitParams.  back()->Write();
               hShiftParams.back()->Write();
            }

            // --- store all PDFs
            for ( auto ipdf : pdfdef ) {
               TGraph gPDF, gReg, gPom;
               gPDF.SetName("DPDF_"+TString(ipdf.first.c_str()));
               for ( auto zp : zval ) {
                  //SET_ANY("LHAPDF.xp",xp,0);
                  //vector<double> xfx = VALUES_ANY("LHAPDF");
                  vector<double> xpom_zp_muf={xp,zp,sqrt(q2)};
                  vector<double> xfx = QUICK_ANY(dpdffunc,xpom_zp_muf);
                  double xf = 0;
                  for ( unsigned int jf=0 ; jf< ipdf.second.size() ; jf++ ) {
                     xf += ipdf.second[jf]*xfx[jf];
                  }
                  gPDF.SetPoint(gPDF.GetN(),zp,xf);
                  AllValues[ipdf.first][xp][q2][zp].push_back(xf);
               }
               gPDF.Write();
            }
            //gDirectory->Write();
         }
      }
      //taskdir->Write();

      // --- Q2/mu_f loop for pom and reg
      for ( auto q2 : q2val ) {
         TString dirname = Form("Q2_%g",q2);
         if ( !taskdir->GetDirectory(dirname) ) taskdir->mkdir(dirname);
         TDirectory* q2dir = taskdir->GetDirectory(dirname);
         if ( !q2dir->GetDirectory(Form("DPDF_%d",iMem))) q2dir->mkdir(Form("DPDF_%d",iMem));
         q2dir->cd(Form("DPDF_%d",iMem));

         // --- store shift of parameters in each directory
         if ( gFitParams.size() ) {
            gFitParams.  back()->Write();
            gShiftParams.back()->Write();
            hFitParams.  back()->Write();
            hShiftParams.back()->Write();
         }

         // --- store all PDFs
         for ( auto ipdf : pdfdef ) {
            TGraph gReg, gPom;
            gPom.SetName("Pom_"+TString(ipdf.first.c_str()));
            gReg.SetName("Reg_"+TString(ipdf.first.c_str()));
            for ( auto zp : zval ) {
               //SET_ANY("LHAPDF.xp",xp,0);
               //vector<double> xfx = VALUES_ANY("LHAPDF");
               vector<double> zp_muf={zp,sqrt(q2)};
               vector<double> pom = QUICK_ANY(dpdffunc+".pom1",zp_muf);
               vector<double> reg = QUICK_ANY(dpdffunc+".reg1",zp_muf);
               double xfp = 0;
               double xfr = 0;
               for ( unsigned int jf=0 ; jf< ipdf.second.size() ; jf++ ) {
                  xfp += ipdf.second[jf]*pom[jf];
                  xfr += ipdf.second[jf]*reg[jf];
               }
               gPom.SetPoint(gPom.GetN(),zp,xfp);
               gReg.SetPoint(gReg.GetN(),zp,xfr);
               AllValPom[ipdf.first][q2][zp].push_back(xfp);
               AllValReg[ipdf.first][q2][zp].push_back(xfr);
            }
            gPom.Write();
            gReg.Write();
         }
         gDirectory->Write();
      }
      //taskdir->Write();
   }


   // for ( auto g : gFitParams )   delete g;
   // for ( auto g : gShiftParams ) delete g;
   // for ( auto h : hFitParams )   delete h;
   // for ( auto h : hShiftParams ) delete h;

   
   // TGraphs with error bands for DPDF
   if ( nPDF>1 && ( FitType) ) {
      // --- DPDF
      info["Execute"]<<"Calculating error bands for "<<xpomval.size()<<" xpom values, and "<<q2val.size()<<" mu_f^2 values."<<endl;
      for ( auto xp : xpomval ) {
         for ( auto q2 : q2val ) {
            TString dirname = Form("Q2_%g__xpom_%0.4f",q2,xp);
            TDirectory* q2dir = taskdir->GetDirectory(dirname);
            //TString dirname = Form("Q2_%.1f",q2);
            // if ( !q2dir->GetDirectory("PDF_Errors") ) q2dir->mkdir("PDF_Errors")->cd();
            // q2dir->cd("PDF_Errors");
            q2dir->mkdir("DPDF_ErrorsSymm");
            q2dir->mkdir("DPDF_ErrorsAsym");

            for ( auto ipdf : pdfdef ) {
               TGraphAsymmErrors gPDFSymm,gPDFAsym;
               gPDFSymm.SetName(ipdf.first.c_str());
               gPDFAsym.SetName(ipdf.first.c_str());
               const map<double,map<double,map<double, vector<double> > > >& iPdfValues = AllValues.at(ipdf.first);
               for ( auto zp : zval ) {
                  int nP = gPDFSymm.GetN();
                  //const vector<double>& values = AllValues[ipdf.first][q2][xp]; // this is a bit slow
                  const vector<double>& values = iPdfValues.at(xp).at(q2).at(zp);
                  double errdn=0,errup=0,errsym=0;
                  {
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
                  gPDFSymm.SetPoint(nP,zp,cent);
                  gPDFAsym.SetPoint(nP,zp,cent);
                  gPDFSymm.SetPointError(nP,0,0,errsym,errsym);
                  gPDFAsym.SetPointError(nP,0,0,errdn,errup);
               }
               // save to file
               q2dir->cd("DPDF_ErrorsAsym");
               gPDFAsym.Write();
               if ( FitType ) {
                  q2dir->cd("DPDF_ErrorsSymm");
                  gPDFSymm.Write();
               }
            }
         }
      }
      //taskdir->Write();

      // --- DPDF
      info["Execute"]<<"Calculating error bands for pom and reg for "<<q2val.size()<<" mu_f^2 values."<<endl;
      for ( auto q2 : q2val ) {
         TString dirname = Form("Q2_%g",q2);
         //TString dirname = Form("Q2_%.1f",q2);
         TDirectory* q2dir = taskdir->GetDirectory(dirname);
         // if ( !q2dir->GetDirectory("PDF_Errors") ) q2dir->mkdir("PDF_Errors")->cd();
         // q2dir->cd("PDF_Errors");
         q2dir->mkdir("Pom_ErrorsSymm");
         q2dir->mkdir("Pom_ErrorsAsym");
         q2dir->mkdir("Reg_ErrorsSymm");
         q2dir->mkdir("Reg_ErrorsAsym");

         for ( auto ipdf : pdfdef ) {
            TGraphAsymmErrors gPomSymm,gPomAsym;
            TGraphAsymmErrors gRegSymm,gRegAsym;
            gPomSymm.SetName("Pom_"+TString(ipdf.first.c_str()));
            gPomAsym.SetName("Pom_"+TString(ipdf.first.c_str()));
            gRegSymm.SetName("Reg_"+TString(ipdf.first.c_str()));
            gRegAsym.SetName("Reg_"+TString(ipdf.first.c_str()));
            const map<double,map<double, vector<double> > >& iPomValues = AllValPom.at(ipdf.first);
            const map<double,map<double, vector<double> > >& iRegValues = AllValReg.at(ipdf.first);
            for ( auto zp : zval ) {
               int nP = gPomSymm.GetN();
               //const vector<double>& values = AllValues[ipdf.first][q2][xp]; // this is a bit slow
               const vector<double>& pomvalues = iPomValues.at(q2).at(zp);
               const vector<double>& regvalues = iRegValues.at(q2).at(zp);
               double Perrdn=0,Perrup=0,Perrsym=0;
               double Rerrdn=0,Rerrup=0,Rerrsym=0;
               {
                  // calculate asymmetric 'hessian' uncertainties
                  // see e.g.:
                  // hep-ph/0101032
                  // arXiv:1101.0536
                  // code adapted from: LHAPDF::src/PDFSet.cc 
                  int nxn = (nPDF-1)/2;
                  for (int ieig = 0; ieig < nxn; ieig++) {
                     Perrup  += AlposTools::sq(max(max(pomvalues[2*ieig+1]-pomvalues[0],pomvalues[2*ieig+2]-pomvalues[0]), 0.));
                     Perrdn  += AlposTools::sq(max(max(pomvalues[0]-pomvalues[2*ieig+1],pomvalues[0]-pomvalues[2*ieig+2]), 0.));
                     Perrsym += AlposTools::sq(pomvalues[2*ieig+1]-pomvalues[2*ieig+2]);
                     Rerrup  += AlposTools::sq(max(max(regvalues[2*ieig+1]-regvalues[0],regvalues[2*ieig+2]-regvalues[0]), 0.));
                     Rerrdn  += AlposTools::sq(max(max(regvalues[0]-regvalues[2*ieig+1],regvalues[0]-regvalues[2*ieig+2]), 0.));
                     Rerrsym += AlposTools::sq(regvalues[2*ieig+1]-regvalues[2*ieig+2]);
                  }
                  Perrsym = 0.5*sqrt(Perrsym);
                  Perrup  = sqrt(Perrup);
                  Perrdn  = sqrt(Perrdn);
                  Rerrsym = 0.5*sqrt(Rerrsym);
                  Rerrup  = sqrt(Rerrup);
                  Rerrdn  = sqrt(Rerrdn);
               }
               {
                  double cent = pomvalues[0];
                  gPomSymm.SetPoint(nP,zp,cent);
                  gPomAsym.SetPoint(nP,zp,cent);
                  gPomSymm.SetPointError(nP,0,0,Perrsym,Perrsym);
                  gPomAsym.SetPointError(nP,0,0,Perrdn,Perrup);
               }
               {
                  double cent = pomvalues[0];
                  gRegSymm.SetPoint(nP,zp,cent);
                  gRegAsym.SetPoint(nP,zp,cent);
                  gRegSymm.SetPointError(nP,0,0,Rerrsym,Rerrsym);
                  gRegAsym.SetPointError(nP,0,0,Rerrdn,Rerrup);
               }
            }
            // save to file
            q2dir->cd("Pom_ErrorsAsym");
            gPomAsym.Write();
            q2dir->cd("Reg_ErrorsAsym");
            gRegAsym.Write();
            if ( FitType ) {
               q2dir->cd("Pom_ErrorsSymm");
               gPomSymm.Write();
               q2dir->cd("Reg_ErrorsSymm");
               gRegSymm.Write();
            }
         }
      }
      //taskdir->Write();
      info["Execute"]<<"TGraphs with errors calculated."<<endl;
   }
   else {
      info["Execute"]<<"No error TGraphs are drawn."<<endl;
   }

   // --------------------------------------------------
   // store cross sections
   //map<string,map<int,vector<double> > > crosssections; // theoname, iMem, values<>
   for ( const auto& iMemCS : crosssections ) {
      string csname = iMemCS.first;
      for ( const auto& csm : iMemCS.second ) {
         int iMem = csm.first;
         const vector<double>& cs = csm.second;
         TString dirname = Form("Theo_Eig_%d",iMem);
         if ( !taskdir->FindObjectAny(dirname) ) 
            taskdir->mkdir(dirname);
         taskdir->cd(dirname);
         std::vector<double> x(cs.size()) ;
         std::iota (std::begin(x), std::end(x), 0); // Fill with 0, 1, ..., cs.size()
         TGraph g(cs.size(), &x[0] ,&cs[0] );
         g.SetName(csname.c_str());
         g.Write();
      }
   }   


   // --------------------------------------------------
   info["Execute"]<<"Writing to disk."<<endl;
   taskdir->Write();
   
   // reset
   if ( Mem0>=0 ) SET_ANY(dpdffunc+".PDFSet",Mem0,0);
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
map<string,vector<double> > ASaveDPDFTGraph::GetDPDFdef() const {
   //! return definition for all PDF linear combinations
   map<string,vector<double> > pdfdef;
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   // pdfdef[   "tb"] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   // pdfdef[   "bb"] = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "cb"] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "sb"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "ub"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[   "db"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "g"] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ,0, 0 ,0 };
   pdfdef[    "d"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ,0, 0 ,0 };
   pdfdef[    "u"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[    "s"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[    "c"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   // pdfdef[    "b"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 1 ,0 };
   // pdfdef[    "t"] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0 ,1 };
   //                     tb bb cb sb ub db  g  d  u  s  c  b  t      
   pdfdef["g*0.05"] = { 0, 0, 0, 0, 0, 0,0.05, 0, 0, 0 ,0, 0 ,0 };
   // pdfdef[   "dv"] = { 0, 0, 0, 0, 0,-1, 0, 1, 0, 0 ,0, 0 ,0 };
   // pdfdef[   "uv"] = { 0, 0, 0, 0,-1, 0, 0, 0, 1, 0 ,0, 0 ,0 };
   pdfdef[   "s+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "s-"] = { 0, 0, 0,-1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c+"] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 ,0, 0 ,0 };
   pdfdef[   "c-"] = { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0 ,1, 0 ,0 };
   // pdfdef[ "Ubar"] = { 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 ,0, 0 ,0 };
   // pdfdef[ "Dbar"] = { 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 ,0, 0 ,0 };
   // pdfdef[    "U"] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,1, 0 ,1 };
   // pdfdef[    "D"] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 ,0, 1 ,0 };
   // pdfdef[   "U+"] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ,1, 0 ,1 };
   // pdfdef[   "D+"] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 ,0, 1 ,0 };
   //
   pdfdef["gluon"] =   {  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,}; // gluon
   pdfdef["SIGMA"] =   {  1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,}; // SIGMA = \sum_q q^+     
   pdfdef["VALENCE"]=  { -1,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1, 1,}; // VALENCE = \sum_q q^-   
   pdfdef[   "T3"] =   {  0, 0, 0, 0, 1,-1, 0,-1, 1, 0, 0, 0, 0,}; // T3  = u^+ - d^+        
   pdfdef[   "V3"] =   {  0, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0, 0, 0,}; // V3  = u^- - d^-        
   // pdfdef[   "T8"] =   {  0, 0, 0,-2, 1, 1, 0, 1, 1,-2, 0, 0, 0,}; // T8  = u^+ + d^+ - 2 s^+
   // pdfdef[   "V8"] =   {  0, 0, 0, 2,-1,-1, 0, 1, 1,-2, 0, 0, 0,}; // V8  = u^- + d^- - 2 s^-
   // pdfdef[  "T15"] =   {  0, 0,-3, 1, 1, 1, 0, 1, 1, 1,-3, 0, 0,}; // T15 = u^+ + d^+ + s^+ - 3 c^+            
   // pdfdef[  "V15"] =   {  0, 0, 3,-1,-1,-1, 0, 1, 1, 1,-3, 0, 0,}; // V15 = u^- + d^- + s^- - 3 c^-            
   // pdfdef[  "T24"] =   {  0,-4, 1, 1, 1, 1, 0, 1, 1, 1, 1,-4, 0,}; // T24 = u^+ + d^+ + s^+ + c^+ - 4 b^+      
   // pdfdef[  "V24"] =   {  0, 4,-1,-1,-1,-1, 0, 1, 1, 1, 1,-4, 0,}; // V24 = u^- + d^- + s^- + c^- - 4 b^-      
   // pdfdef[  "T35"] =   { -5, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,-5,}; // T35 = u^+ + d^+ + s^+ + c^+ + b^+ - 5 t^+
   // pdfdef[  "V35"] =   {  5,-1,-1,-1,-1,-1, 0, 1, 1, 1, 1, 1,-5,};  // V35 = u^- + d^- + s^- + c^- + b^- - 5 t^-

   return pdfdef;

}


//____________________________________________________________________________________ //
std::vector<double> ASaveDPDFTGraph::GetLogNodes(double min, double max, int npts) const {
   std::vector<double> ret(npts+1);
   double lmin = log(min);
   double lmax = log(max);
   for ( int i=0 ; i<npts+1 ; i++ ) {
      ret[i] = exp(lmin+(lmax-lmin)/npts*i);
   }
   return ret;
}
 
//____________________________________________________________________________________ //
void ASaveDPDFTGraph::WriteAlphasGrid(ostream& strm, const vector<double>& qpt, vector<double>* grid0){
   //! write alpha_s grid to stream
   if ( !EXIST_NS(Alphas,NS()) ) {
      error["WriteAlphasGrid"]<<"Alpha_s function needs to be given to ASaveDPDFTGraph, when writing LHAPDF file."<<endl;
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
void ASaveDPDFTGraph::WriteDPDFGrid(ostream& strm, const vector<int>& PDFid, const vector<double>& xpt, const vector<double>& qpt, vector<double>* grid0){
   return ;
   /*
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
   
   string dpdffunc = STRING_NS(DPDF,NS());

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
   vector<double> xfx = QUICK_ANY(dpdffunc,xp_muf);
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
   */
}
