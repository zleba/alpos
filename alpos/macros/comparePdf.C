R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

R__LOAD_LIBRARY($PROJECT_DIR/alposBuild/libaem.so)

R__ADD_INCLUDE_PATH($PROJECT_DIR/dpdfLib)
R__LOAD_LIBRARY($PROJECT_DIR/dpdfLib/DPDFset.so)
R__LOAD_LIBRARY(/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-c9645/x86_64-slc6-gcc7-opt/lib/libLHAPDF.so)

#include "DPDFset.h"
#include "plottingHelper.h"
   using namespace PlottingHelper;//pollute the namespace!


   TString rn() { return Form("%d", rand()); }


   struct point {
      double xp,  q2,  beta, xpSig;
      double th, thErr;
      double thOrgA, thOrgB;
      double errStat, errSys, errTot, errUnc;
      std::vector<double> errs;//10 items
   };

TGraph *corrMuShift(TGraph *gr, TGraph *grRef);
pair<TGraph*,TGraph*>  GetUpDnCurve(TGraphAsymmErrors* gr);

struct sysShift {
   std::vector<point> data;
   TH1D *hPars, *hFixedPars;
   TH2D *hCorrs;
   map<double, vector<TGraph*>> singletQ2;
   map<double, vector<TGraph*>> gluonQ2;
   double chi2;
   int ndf;
   bool readData(TString inFile);
};

bool sysShift::readData(TString inFile)
{
   TFile *file = TFile::Open(inFile);

   hPars  = dynamic_cast<TH1D*>(file->Get("fitparameters"));
   hFixedPars  = dynamic_cast<TH1D*>(file->Get("ASaveDataTheory/hFixedPars"));
   hCorrs = dynamic_cast<TH2D*>(file->Get("fitcorrelations"));

   if(!hPars || !hCorrs) {
      cout << "Not loaded " << __LINE__ <<" : "<< inFile << endl;
      return false;
   }
   TGraph *grXsec = dynamic_cast<TGraph*>(file->Get("SaveDPDFTGraph/Theo_Eig_0/SuperTheory"));
   if(!grXsec) {
      cout <<"ndf issue " <<  inFile << endl;
      ndf = 0;
      //return false;
   }
   else {
      ndf  = grXsec->GetN() - hPars->GetNbinsX();
   }

   chi2 = (dynamic_cast<TFitResult*>(file->Get("fitresult")))->MinFcnValue();


   int nShifts = 2*hPars->GetNbinsX() + 1;
   vector<double> q2Vals = {1.75, 8.5, 20, 90, 800};
   //vector<double> q2Vals = {2.5, 8.5, 20, 90, 800};
   for(double q2 : q2Vals) {
      singletQ2[q2].resize(nShifts, nullptr);
      gluonQ2[q2].resize(nShifts, nullptr);

      //cout << "Start loop " << q2 << endl;
      for(int i = 0; i < nShifts; ++i) {
         //cout << Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_gluon", q2, i ) << endl;
         //cout << Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_SIGMA", q2, i ) << endl;
         gluonQ2[q2][i]   = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_gluon", q2, i )));
         singletQ2[q2][i] = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_SIGMA", q2, i )));
         if(!gluonQ2[q2][i] || !singletQ2[q2][i]) {
            cout << "Not loaded " << __LINE__  <<" "<< q2 <<" "<< i <<" : "<< inFile << endl;
            return false;
            //gluonQ2[q2][i]   = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_gluon", q2, i )));
            //singletQ2[q2][i] = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_SIGMA", q2, i )));
            //exit(1);
         }
      }
      //cout << "End loop " << q2 << endl;
   }
   //cout << "RADEK end" << endl;
   return true;
}


struct PDF {
   vector<sysShift> shifts;
   TString name;

   PDF(TString inDir) { readData(inDir); }
   PDF() {}

   //static PDF loadLHAPDF(TString name);

   pair<TGraphAsymmErrors*,TGraphAsymmErrors*> getBandModel(double q2);
   void readData(TString inDir, int nErr = 100);
   /*
      void init() {
      gr = new TGraph(20);
      int i = 0;
      for(double a = -5; a <= 5; a += 0.1) {
      gr->SetPoint(i, a, sin(a));
      }
      }
      */

};

PDF loadLHAPDF(TString name, TString name2 = "")
{
   PDF pdf;
   DPDFset  dpdf(name.Data(), name2.Data());
   pdf.name = name;

   pdf.shifts.resize(1);// dpdf.size());

   vector<double> q2s = {1.75, 8.5, 20, 90, 800};

   const int npt     = 300;
   const double zmin = 0.001;
   for(auto q2 : q2s) {
      for(int i = 0; i < dpdf.size(); ++i) {
         //double q2 = q2s[i];
         pdf.shifts[0].singletQ2[q2].resize(dpdf.size());
         pdf.shifts[0].singletQ2[q2][i] = new TGraphAsymmErrors(npt);
         pdf.shifts[0].gluonQ2[q2].resize(dpdf.size());
         pdf.shifts[0].gluonQ2[q2][i]   = new TGraphAsymmErrors(npt);

         double z = zmin;
         for(int k = 0; k < npt; ++k) {
            double z = zmin * pow(1/zmin, k/(npt-1.));

            double xp = 0.003;
            vector<double> pdfsNow(13);
            dpdf.zfzQ2xp(i, z, q2, xp, 0, 1, pdfsNow);
            double gluon    = xp*pdfsNow[6];
            double singlet  = 0;
            for(int l = 0; l < pdfsNow.size(); ++l) if(l != 6) singlet += xp*pdfsNow[l];

            if(z < 0.999) {
            pdf.shifts[0].singletQ2[q2][i]->SetPoint(k, z, singlet);
            pdf.shifts[0].gluonQ2[q2][i]->SetPoint(k, z, gluon);
            }
         }
      }
   }

   return pdf;


}


void PDF::readData(TString inDir, int nErr)
{
   name = inDir(inDir.Last('/')+1, inDir.Last('.') - inDir.Last('/') - 1 );

   for(int i = 0; i <= nErr; ++i) {
      TString fName = inDir+Form("/steering.str%d_dir/out.root",i);

      ifstream testFile(fName.Data());
      if(testFile.fail()) {
         cout <<"FAILED: " <<  fName << endl;
         break;
      }
      sysShift shift;
      if(shift.readData(fName)) {

         /*
         //Correction for scale variation
         if(i == 15 || i == 16) {
            for(const auto & s : shift.singletQ2) {
               double q2 = s.first;
               for(int k = 0; k < shift.singletQ2[q2].size(); ++k) {
                  shift.singletQ2[q2][k] = corrMuShift(shift.singletQ2[q2][k], shifts[0].singletQ2[q2][k]);
                  shift.gluonQ2[q2][k] = corrMuShift(shift.gluonQ2[q2][k], shifts[0].gluonQ2[q2][k]);
               }
            }
         }
         */

         shifts.push_back(shift);
      }
      else {
         break;
      }
      //shifts.resize(shifts.size()+1);
      //shifts.back().readData(fName);
   }
   cout << "After reading " << endl;

   if(shifts.size() == 0) {
      cout << "Nothing read" << endl;
      exit(1);
   }

}






struct Canvas {

   TCanvas *can;
   bool doRatio, inLog, compErrors;
   int nPdfs = 0;
   map<double,TGraph *> grSRef, grGRef;
   //vector<TString> pdfsNames;
   TLegend *leg;

   Canvas() {
      can = new TCanvas(rn(), "", 600, 600);

      TH1D *h = new TH1D(rn(), "", 1, 1e-5, 1);
      h->Draw();
      h->GetYaxis()->SetRangeUser(-3, 3);
   }
};

Canvas can() {
   Canvas c;
   return c;
}

struct saver {
   TString fName;
};
saver save(TString fName_) {
   saver s;
   s.fName = fName_;
   return s;
}

/*
   Canvas operator<<(Canvas can, PDF pdf)
   {
   can.can->cd();
   pdf.gr->Draw("l");
   return can;
//return canvas();
}
*/

Canvas operator<<(Canvas can, saver s)
{
   //TString n1 = can.pdfsNames[0];
   //TString n2 = can.pdfsNames[1];

   can.leg->Draw();

   cout << "Saving " << endl;
   can.can->SaveAs(s.fName);
   return can;
}

Canvas pdfPlot(TString config);
//Canvas &operator<<(Canvas &canvas, PDF &pdf);



//Make template
Canvas pdfPlot(TString config) 
{
   Canvas canvas;
   canvas.inLog =  config.Contains("Log"); 
   canvas.doRatio = config.Contains("Ratio");
   canvas.compErrors = config.Contains("CompErrors");
   if(canvas.compErrors) canvas.doRatio = true;

   vector<double> q2s = {1.75, 8.5, 20, 90, 800};

   canvas.can = new TCanvas(rn(),"", 600, 600);
   TCanvas *can = canvas.can;
   gStyle->SetOptStat(0);
   SetLeftRight(0.1, 0.14);
   SetTopBottom(0.1, 0.1);

   double zMin = 4e-3;
   DivideTransparent(group(1, 0.5, 2), group(1, 0, q2s.size()));

   for(int i = 0; i < q2s.size(); ++i) {
      //Singlet
      can->cd(2*i + 1);
      TH1D *hFrS = new TH1D(rn(), "", 1, zMin, 1);
      hFrS->Draw("axis");

      if(canvas.doRatio) GetYaxis()->SetRangeUser(0.4, 1.6);
      else        GetYaxis()->SetRangeUser(0, 0.27);
      //GetYaxis()->SetRangeUser(0.9, 1.1);
      GetYaxis()->SetNdivisions(503);
      GetXaxis()->SetNdivisions(404);
      SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});

      if(canvas.inLog) gPad->SetLogx();

      if(i == 0) {
         DrawLatexUp(-1, "Singlet");
         if(canvas.doRatio) GetYaxis()->SetTitle("#Sigma(z,Q^{2}) / #Sigma_{0}(z,Q^{2})");
         else        GetYaxis()->SetTitle("z #Sigma(z,Q^{2})");
      }
      if(i == q2s.size()-1) {
         GetXaxis()->SetTitle("z");
      }
      else {
         GetXaxis()->SetLabelOffset(50000);
      }

      can->cd(2*i + 2);
      TH1D *hFrG = new TH1D(rn(), "", 1, zMin, 1);
      hFrG->Draw("axis");

      if(canvas.doRatio) GetYaxis()->SetRangeUser(0.4, 1.6);
      else        GetYaxis()->SetRangeUser(0, 2.25);

      //GetYaxis()->SetRangeUser(0.9, 1.1);
      GetYaxis()->SetNdivisions(503);
      GetXaxis()->SetNdivisions(404);
      SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});
      if(canvas.inLog) gPad->SetLogx();

      if(i == 0) {
         DrawLatexUp(-1, "Gluon");
         if(canvas.doRatio) GetYaxis()->SetTitle("g(z,Q^{2}) / g_{0}(z,Q^{2})");
         else        GetYaxis()->SetTitle("z g(z,Q^{2})");
      }
      if(i == q2s.size()-1) {
         GetXaxis()->SetTitle("z");
      }
      else {
         GetXaxis()->SetLabelOffset(50000);
      }

      DrawLatexRight(1, Form("Q^{2}=%g",q2s[i]), -1, "l");
   }

   canvas.leg = new TLegend(0.3, 0.92, 0.95, 0.99);
   canvas.leg->SetBorderSize(0);
   canvas.leg->SetFillStyle(0);

   return canvas;
}



TGraphAsymmErrors *getBand(vector<TGraph*> graphs);
TGraphAsymmErrors *addBands(vector<TGraphAsymmErrors*> graphs);
TGraphAsymmErrors *GetFraction(TGraphAsymmErrors *gr, TGraph *grRef);


struct Style {
   int lineCol;
   int innerCol;
   int outerCol;
   int fillStyle;
};

vector<Style> pdfStyles = {   {kBlue,kBlue,kRed, 1001}, {kGreen,kGreen,kGreen+3, 3345 /*3244*/}   };


Canvas operator <<(Canvas canvas, PDF pdf)
{
   TCanvas *can = canvas.can;

   vector<double> q2s = {1.75, 8.5, 20, 90, 800};

   double zMin = 4e-3;

   auto pdfStyle = pdfStyles[canvas.nPdfs];

   for(int i = 0; i < q2s.size(); ++i) {
      //Fill Graph
      cout << "Helenka " <<  pdf.shifts[0].singletQ2.at(q2s[i]).size() << " "<<pdf.shifts[0].singletQ2.at(q2s[i])[0]<< endl;
      TGraphAsymmErrors *grS = getBand(pdf.shifts[0].singletQ2.at(q2s[i]));
      TGraphAsymmErrors *grG = getBand(pdf.shifts[0].gluonQ2.at(q2s[i]));

      TGraphAsymmErrors *grSm, *grGm;
      tie(grGm, grSm) = pdf.getBandModel(q2s[i]);

      TGraphAsymmErrors *grStot = addBands({grS, grSm});
      TGraphAsymmErrors *grGtot = addBands({grG, grGm});

      grS->SetLineColor(pdfStyle.lineCol);
      grG->SetLineColor(pdfStyle.lineCol);

      grG->SetLineWidth(3);
      grS->SetLineWidth(3);
      grGtot->SetLineWidth(3);
      grStot->SetLineWidth(3);

      if(canvas.doRatio) {
         TGraph *grSRefNow, *grGRefNow;
         if(canvas.compErrors) {
            grSRefNow = (TGraph*) grS->Clone();
            grGRefNow =(TGraph*) grG->Clone();
         }
         else {
            if(canvas.nPdfs == 0) {
               canvas.grSRef[q2s[i]] = (TGraph*) grS->Clone();
               canvas.grGRef[q2s[i]] = (TGraph*) grG->Clone();
            }
            grSRefNow = canvas.grSRef[q2s[i]];
            grGRefNow = canvas.grGRef[q2s[i]];
         }
         grS    = GetFraction(grS, grSRefNow);
         grStot = GetFraction(grStot, grSRefNow);

         grG    = GetFraction(grG, grGRefNow);
         grGtot = GetFraction(grGtot, grGRefNow);
      }


      can->cd(2*i + 1);
      double tr =  pdfStyle.fillStyle == 1001 ? 0.5 : 0.5;

      grStot->SetFillColorAlpha(pdfStyle.outerCol, tr);
      grStot->SetFillStyle(pdfStyle.fillStyle);
      if(pdfStyle.fillStyle == 1001) grStot->Draw("le03 same");
      else                           grStot->Draw("le02 same");

      auto plotBorders = [&](TGraphAsymmErrors *grNow) {
         if(pdfStyle.fillStyle != 1001 || true) {
            TGraph *grStotup, *grStotdn;
            tie(grStotup, grStotdn) =  GetUpDnCurve(grNow);
            grStotup->Draw("l same");
            grStotdn->Draw("l same");
         }
      };
      plotBorders(grStot);

      grS->SetFillColorAlpha(pdfStyle.innerCol, tr);
      grS->SetFillStyle(pdfStyle.fillStyle);
      if(pdfStyle.fillStyle == 1001) grS->Draw("le03 same ][");
      else                           grS->Draw("le02 same ][");
      plotBorders(grS);

      can->cd(2*i + 2);

      grGtot->SetFillColorAlpha(pdfStyle.outerCol, tr);
      grGtot->SetFillStyle(pdfStyle.fillStyle);
      if(pdfStyle.fillStyle == 1001) grGtot->Draw("le03 same ][");
      else                           grGtot->Draw("le02 same ][");
      plotBorders(grGtot);

      grG->SetFillColorAlpha(pdfStyle.innerCol, tr);
      grG->SetFillStyle(pdfStyle.fillStyle);
      if(pdfStyle.fillStyle == 1001) grG->Draw("le03 same ][");
      else                           grG->Draw("le02 same ][");
      plotBorders(grG);

      if(i == q2s.size() - 1)
         //canvas.leg->AddEntry(grG, Form("%s (#chi^{2}/ndf=%.2f, n_{sh}=%zu)", pdf.name.Data(), pdf.shifts[0].chi2 / pdf.shifts[0].ndf, pdf.shifts.size()), "lf");
         canvas.leg->AddEntry(grG, Form("%s (#chi^{2}/ndf=%.2f)", pdf.name.Data(), pdf.shifts[0].chi2 / pdf.shifts[0].ndf), "lf");


   }
   cout << "Marketka" << endl;

   ++canvas.nPdfs;

   return canvas;

}

TGraphAsymmErrors *getBand(vector<TGraph*> graphs)
{
   TGraphAsymmErrors *gr = new TGraphAsymmErrors(graphs[0]->GetN());
   //cout << "Radek " << graphs.size() << endl;
   for(int i = 0; i < graphs[0]->GetN(); ++i) {
      double x, vCnt;
      graphs[0]->GetPoint(i, x, vCnt);

      double errP = 0, errM = 0;
      for(auto g : graphs) {
         double xTmp, vNow;
         g->GetPoint(i, xTmp, vNow) ;
         double v = vNow - vCnt;
         errP = hypot(errP, max(0.0, v));
         errM = hypot(errM, max(0.0,-v));
      }

      double xLeft = x, xRight = x;
      if(i > 0) {
         double vTmp;
         graphs[0]->GetPoint(i-1, xLeft, vTmp);
      }
      if(i < graphs[0]->GetN() - 1) {
         double vTmp;
         graphs[0]->GetPoint(i+1, xRight, vTmp);
      }

      gr->SetPoint(i, x, vCnt);
      gr->SetPointError(i, (x-xLeft)/2, (xRight-x)/2, errM, errP);

   }
   return gr;

}

TGraph *corrMuShift(TGraph *gr, TGraph *grRef)
{
   if(gr->GetN() != grRef->GetN()) {
      cout << "Problem " <<__LINE__ << endl;
      exit(0);
   }

   double fact = 1;//log(2)/log(1.1);

   TGraph *grNew = (TGraph*) gr->Clone(rn());
   for(int i = 0; i < grNew->GetN(); ++i) {
      double x, y, x0, y0;
      gr->GetPoint(i, x, y);
      grRef->GetPoint(i, x0, y0);

      double d = y - y0;

      double yn = y0 + d*fact;
      grNew->SetPoint(i, x, yn);
   }
   return grNew;
}

TGraphAsymmErrors *addBands(vector<TGraphAsymmErrors*> graphs)
{
   TGraphAsymmErrors *gr = (TGraphAsymmErrors*) graphs[0]->Clone(rn());
   for(int i = 0; i < graphs[0]->GetN(); ++i) {
      double x, vCnt;
      double errP=0, errM = 0;
      double el, er, eM, eP;
      for(int j = 0; j < graphs.size(); ++j) {
         errP = hypot(errP, graphs[j]->GetErrorYhigh(i));
         errM = hypot(errM, graphs[j]->GetErrorYlow(i));
      }
      gr->SetPointEYhigh(i, errP);
      gr->SetPointEYlow(i, errM);
   }
   return gr;
}

pair<TGraph*,TGraph*>  GetUpDnCurve(TGraphAsymmErrors* gr)
{
   TGraph *grUp = new TGraph(gr->GetN());
   TGraph *grDn = new TGraph(gr->GetN());

   grUp->SetLineColor(gr->GetFillColor());
   grDn->SetLineColor(gr->GetFillColor());
   grUp->SetLineStyle(1);
   grDn->SetLineStyle(1);
   grUp->SetLineWidth(2);
   grDn->SetLineWidth(2);

   for(int i = 0; i < gr->GetN(); ++i) {
      double x, y;
      gr->GetPoint(i, x, y);
      double errP = gr->GetErrorYhigh(i);
      double errM = gr->GetErrorYlow(i);
      grUp->SetPoint(i, x, y+errP);
      grDn->SetPoint(i, x, y-errM);
   }
   return {grUp, grDn};
}




pair<TGraphAsymmErrors*,TGraphAsymmErrors*> PDF::getBandModel(double q2)
{
   vector<TGraph*> gluons, singlets;
   for(int i = 0; i < shifts.size(); ++i) {
      singlets.push_back(shifts[i].singletQ2.at(q2)[0]);
      gluons.push_back(shifts[i].gluonQ2.at(q2)[0]);
   }

   cout << " radek " << singlets.size() <<" "<< gluons.size() <<endl;
   return {getBand(gluons), getBand(singlets)};
}

TGraphAsymmErrors *GetFraction(TGraphAsymmErrors *gr, TGraph *grRef)
{
   TGraphAsymmErrors *grRat = (TGraphAsymmErrors*) gr->Clone(rn());
   for(int i = 0; i < gr->GetN(); ++i) {
      double x, v;
      double xR, vRef;
      gr->GetPoint(i, x,  v);
      grRef->GetPoint(i, xR, vRef);
      if(x != xR) exit(1);

      double refInv = (vRef == 0) ? 0 : 1./vRef;

      double vRat = v * refInv;
      double eH = gr->GetErrorYhigh(i) * refInv;
      double eL = gr->GetErrorYlow(i) * refInv;

      grRat->SetPoint(i, x, vRat);
      grRat->SetPointEYhigh(i, eH);
      grRat->SetPointEYlow(i, eL);
   }
   return grRat;
}


void plotAllTypes(PDF &pdf1, PDF &pdf2, TString name)
{
      pdfPlot("Log")           << pdf1  << pdf2 << save(name+"_Log.pdf");
      cout << "Log Rat plotting" << endl;
      pdfPlot("Log Ratio")     << pdf1  << pdf2 << save(name+"_LogRat.pdf");
      cout << "Log CompErrors" << endl;
      pdfPlot("Log CompErrors")<< pdf1  << pdf2 << save(name+"_LogErrs.pdf");
      cout << "Lin" << endl;
      pdfPlot("Lin")           << pdf1  << pdf2 << save(name+"_Lin.pdf");
      pdfPlot("Lin Ratio")     << pdf1  << pdf2 << save(name+"_LinRat.pdf");
      pdfPlot("Lin CompErrors")<< pdf1  << pdf2 << save(name+"_LinErrs.pdf");
}







void comparePdf()
{
   cout << "Start" << endl;
   //canvas();

   vector<vector<TString>> compNames = {
      /*
      //Theor comp.
      {"Ext_nloF_heraI", "Ext_nlo_heraI" }, 
      {"Ext_nlo_heraI", "Ext_nnlo_heraI" },
      {"Ext_nloF_heraC", "Ext_nlo_heraC" },
      {"Ext_nlo_heraC", "Ext_nnlo_heraC" },
      {"Ext_nloF_heraCjets", "Ext_nlo_heraCjets" },
      {"Ext_nlo_heraCjets", "Ext_nnlo_heraCjets" },

      //Data comp.
      {"Ext_nloF_heraI", "Ext_nloF_heraC" },
      {"Ext_nlo_heraC", "Ext_nlo_heraCjets" },
      {"Ext_nlo_heraI", "Ext_nlo_heraC" },
      {"Ext_nnlo_heraI", "Ext_nnlo_heraC" },
      {"Ext_nnlo_heraC", "Ext_nnlo_heraCjets" },
      {"Ext_nnlo_heraI", "Ext_nnlo_heraIjets" }
      */

      //Apfel vs QCDNUM comp.
      /*
         {"Ext_nlo_heraI", "AExt_nlo_heraI" }, 
         {"Ext_nnlo_heraI", "AExt_nnlo_heraI" },
         {"Ext_nlo_heraC", "AExt_nlo_heraC" },
         {"Ext_nnlo_heraC", "AExt_nnlo_heraC" },

         {"Ext_nlo_heraIjets", "AExt_nlo_heraIjets" }, 
         {"Ext_nnlo_heraIjets", "AExt_nnlo_heraIjets" },
         {"Ext_nlo_heraCjets", "AExt_nlo_heraCjets" },
         {"Ext_nnlo_heraCjets", "AExt_nnlo_heraCjets" },
         */

      { "AExt_nnlo_heraC", "BExt_nnlo_heraC"}, //Theory impact
      { "AExt_nlo_heraC", "BExt_nlo_heraC"}, //Theory impact

      { "AExt_nnlo_heraCjets", "BExt_nnlo_heraCjets"}, //Theory impact
      { "AExt_nlo_heraCjets", "BExt_nlo_heraCjets"}, //Theory impact


      //{ "AExt_nnlo_heraCjets", "AExt_nlo_heraCjets"}, //Theory impact
      //{ "AExt_nnlo_heraC", "AExt_nlo_heraI"} //Theory impact

      //Theor comp.
//      { "AExt_nnlo_heraC",  "AExt_nnlo_heraII"} //old fit vs new fit

      //{ "AExt_nnlo_fps3D",  "AExt_nnlo_fps4D"},//old fit vs new fit
      //{ "AExt_nnlo_heraC",  "AExt_nnlo_fps4D"},//old fit vs new fit

      //{ "AExt_nnlo_heraC",  "AExt_nnlo_heraCfps3D"},//old fit vs new fit
      //{ "AExt_nnlo_heraC",  "AExt_nnlo_heraCfps4D"}//old fit vs new fit


      /*
      { "AExt_nnlo_heraCjets",  "AExt_nlo_heraI"}, //old fit vs new fit
      { "AExt_nnlo_heraCjets", "AExt_nlo_heraCjets"}, //Theory impact
      { "AExt_nnlo_heraC", "AExt_nlo_heraC"}, //Theory impact
      { "AExt_nlo_heraCjets", "AExt_nlo_heraI" }, //Data impact NLO
      { "AExt_nnlo_heraCjets", "AExt_nnlo_heraI" }, //Data impact NNLO
      */



    /* {"Ext_nloF_heraI", "AExt_nlo_heraI" }, 
      {"AExt_nlo_heraI", "AExt_nnlo_heraI" },
      {"Ext_nloF_heraC", "AExt_nlo_heraC" },
      {"AExt_nlo_heraC", "AExt_nnlo_heraC" },
      {"Ext_nloF_heraCjets", "AExt_nlo_heraCjets" },
      {"AExt_nlo_heraCjets", "AExt_nnlo_heraCjets" },

      //Data comp.
      {"Ext_nloF_heraI", "Ext_nloF_heraC" },
      {"AExt_nlo_heraC", "AExt_nlo_heraCjets" },
      {"AExt_nlo_heraI", "AExt_nlo_heraC" },
      {"AExt_nnlo_heraI", "AExt_nnlo_heraC" }, */

      /*
      { "AExt_nnlo_heraCjets", "AExt_nnlo_heraC" }, //jets impact NNLO
      { "AExt_nlo_heraCjets",  "AExt_nlo_heraC"  }  //jets impact NLO
      */
      //{ "AExt_nnlo_heraIjets", "AExt_nnlo_heraI" }, //jets impact NNLO
      //{ "AExt_nlo_heraIjets",  "AExt_nlo_heraI"  }  //jets impact NLO

      //{ "AExt_nnlo_heraC", "AExt_nnlo_heraI" },
      //{ "AExt_nlo_heraC",  "AExt_nlo_heraI"  }




   };

   TString pathIn = "../farm/variants/";

   /*

   PDF pdfFitBreg =  loadLHAPDF("lhaTest_reg");
   PDF pdfFitBpom =  loadLHAPDF("lhaTest_pom");
   PDF pdfABCDE   =  loadLHAPDF("ABCDE_pom");
   //PDF pdfZEUS =  loadLHAPDF("ZEUS_DPDF_2009SJ_NLO_pom");
   PDF pdfNew =  loadLHAPDF("lhaTest_reg");
   TString tag1 = "AExt_nlo_heraCjets";
   TString tag2 = "AExt_nnlo_heraC";
   PDF pdfH1(pathIn + tag1 + ".str_dir");
   PDF pdfH2(pathIn + tag2 + ".str_dir");
   //pdfPlot("Log CompErrors")       <<  pdfH1 << pdfZEUS   << save("zeus.pdf");
   pdfPlot("Lin")       <<  pdfH2 << pdfABCDE   << save("xcheck_abcde.pdf");
   //TString name = "plots/" +  tag +"__"+ comp[1];
   //plotAllTypes(pdfH1, pdfH2, "plotsZEUS/h1vszeus");

   return;
   */

   for(auto comp: compNames) {
      PDF pdf1(pathIn + comp[0] + ".str_dir");
      PDF pdf2(pathIn + comp[1] + ".str_dir");

      TString name = "plots/" +  comp[0] +"__"+ comp[1];
      //pdf.init();
      pdfPlot("Log")           << pdf1  << pdf2 << save(name+"_Log.pdf");
      pdfPlot("Log Ratio")     << pdf1  << pdf2 << save(name+"_LogRat.pdf");
      pdfPlot("Log CompErrors")<< pdf1  << pdf2 << save(name+"_LogErrs.pdf");
      pdfPlot("Lin")           << pdf1  << pdf2 << save(name+"_Lin.pdf");
      pdfPlot("Lin Ratio")     << pdf1  << pdf2 << save(name+"_LinRat.pdf");
      pdfPlot("Lin CompErrors")<< pdf1  << pdf2 << save(name+"_LinErrs.pdf");
   }
}
