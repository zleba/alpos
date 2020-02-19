R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)
#include "plottingHelper.h"
using namespace PlottingHelper;


//Plot the results of the fit
//i.e. the fit quality of the data on the inclusive or jet data points



#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TNtuple.h"

#include <vector>
#include <map>
#include <algorithm>


using namespace PlottingHelper;//pollute the namespace!
using namespace std;

TString rn() {return Form("%d",rand());}


#include <vector>


struct pointInc {
    double xp,  q2,  beta, tAbs, xpSig;
    double th, thErr;
    double thOrgA, thOrgB;
    double errStat, errSys, errTot, errUnc;
    bool isInside;
    std::vector<double> errs;//10 items
};


struct dPlotter {
    std::vector<pointInc> data;
    std::map<pair<double,double>, TH1D*> jetsData, jetsTh;

    TH1D *hPars;
    TH2D *hCorrs;
    TString outDir;

    void readDataInc(TString inFile, vector<TString> samples);
    void readDataJets(TString inFile, vector<TString> samples);
    void plotBeta(TString fTag, double xpom);
    void plotBetaRat(TString fTag, double xpom);
    void plotBetaMore(TString fTag, vector<double> xpomVec);

    void plotQ2(TString fTag, double xpom);
    void plotXpom();
    void plotXpomT();
    void plotBslope();
    void plotPDFs(bool inLog);
    void plotParameters();
    void plotCorrelations();


    void plotJets();


};

void plotJets(dPlotter &nlo, dPlotter &nnlo);


void dplotter(TString inFile = "../testA/alpos.out.root")
{
    dPlotter dplt;

    inFile = "../farm/variants/AExt_nnlo_heraCjets.str_dir/steering.str0_dir/out.root";
    //inFile = "../farm/variants/Ext_nloF_heraI.str_dir/steering.str0_dir/out.root";

    inFile = "../farm/variants/AExt_nnlo_heraCfps4D.str_dir/steering.str0_dir/out.root";


    //inFile = "../steering.str0_FactorizationNNLO/out.root";

    //dplt.readDataInc(inFile, {"H1incDDIS_HERA_I_LAr_cuts", "H1incDDIS_HERA_I_SpacMB_cuts", "H1incDDIS_HERA_I_SpacTrg_cuts"});
    //dplt.readDataInc(inFile, {});

    dplt.readDataInc(inFile, {"H1incDDIS_Comb"});

    dplt.plotXpom();

    dplt.plotBeta("comb", 0.0003);
    dplt.plotBeta("comb", 0.001);
    dplt.plotBeta("comb", 0.003);
    dplt.plotBeta("comb", 0.01);
    dplt.plotBeta("comb", 0.03);


    /*
    dplt.plotBetaRat("comb", 0.0003);
    dplt.plotBetaRat("comb", 0.001);
    dplt.plotBetaRat("comb", 0.003);
    dplt.plotBetaRat("comb", 0.01);
    dplt.plotBetaRat("comb", 0.03);
    */


    dplt.plotQ2("comb", 0.0003);
    dplt.plotQ2("comb", 0.001);
    dplt.plotQ2("comb", 0.003);
    dplt.plotQ2("comb", 0.01);
    dplt.plotQ2("comb", 0.03);


    dPlotter dpltFPS;
    dpltFPS.readDataInc(inFile, {"H1FPS_4D"});
    dpltFPS.plotXpomT();
    dpltFPS.plotBslope();
    //return;

    dPlotter dplt252;
    dplt252.readDataInc(inFile, {"H1incDDIS_LowE252"});
    //dplt252.plotBeta("e252", 0.0005);
    //dplt252.plotBeta("e252", 0.003);
    dplt252.plotBetaMore("e252", {0.003, 0.0005});


    dPlotter dplt225;
    dplt225.readDataInc(inFile, {"H1incDDIS_LowE225"});
    dplt225.plotBetaMore("e225", {0.003, 0.0005});

    //dplt225.plotBeta("e225", 0.0005);
    //dplt225.plotBeta("e225", 0.003);

    //dplt225.plotBetaRat("e225", 0.0005);
    //dplt225.plotBetaRat("e225", 0.003);

    /*
    dPlotter dpltNNLO;
    dpltNNLO.readDataJets("../farm/variants/AExt_nnlo_heraCjets.str_dir/steering.str0_dir/out.root", {"H1_LRG_DiffDijets"});
    dPlotter dpltNLO;
    dpltNLO.readDataJets("../farm/variants/AExt_nlo_heraCjets.str_dir/steering.str0_dir/out.root", {"H1_LRG_DiffDijets"});
    plotJets(dpltNLO, dpltNNLO);
    return;
    */

    //dplt.readDataInc(inFile, {"H1incDDIS_LowE252_cuts"});
    //dplt.plotBeta(0.003);
    //dplt.plotBeta(0.0005);


    //dplt.plotParameters();
    //dplt.plotCorrelations();

    //dplt.plotBeta(0.03f);

}

void dPlotter::readDataInc(TString inFile, vector<TString> samples)
{
    TFile *file = TFile::Open(inFile);


    hPars  = dynamic_cast<TH1D*>(file->Get("fitparameters"));
    hCorrs = dynamic_cast<TH2D*>(file->Get("fitcorrelations"));
    assert(hPars);
    assert(hCorrs);



    //Loop over all
    //vector<TString> dataNames = {"sample_H1incDDIS_HERA_I_LAr_cuts", "sample_H1incDDIS_HERA_I_SpacMB_cuts", "sample_H1incDDIS_HERA_I_SpacTrg_cuts"};
    for(auto n : samples) {
       TTree *tuple = dynamic_cast<TTree*>(file->Get("ASaveDataTheory/sample_" + n));
       //TNtuple *tuple = (TNtuple*) (file->Get("ASaveDataTheory/sample_" + n));


       pointInc pt;
       Char_t isIn;
       tuple->SetBranchAddress("xp",&pt.xp);
       tuple->SetBranchAddress("isInside",&isIn);
       tuple->SetBranchAddress("Q2",&pt.q2);
       tuple->SetBranchAddress("beta",&pt.beta);
       if(n.Contains("FPS_4D"))
           tuple->SetBranchAddress("tAbs",&pt.tAbs);
       tuple->SetBranchAddress("xpSigData",&pt.xpSig);
       tuple->SetBranchAddress("xpSigDataErr",&pt.errTot);
       tuple->SetBranchAddress("xpSigTh",&pt.th);
       tuple->SetBranchAddress("xpSigThErr",&pt.thErr);

       tuple->SetBranchAddress("xpSigThOrgA",&pt.thOrgA);
       tuple->SetBranchAddress("xpSigThOrgB",&pt.thOrgB);


       Int_t nentries = (Int_t)tuple->GetEntries();
       for (Int_t i=0;i<nentries;i++) {
           tuple->GetEntry(i);
            pt.isInside = isIn;
           //cout << "Bla " << pt.q2 <<" "<< pt.isInside << endl;
           //cout << "Hel " << pt.xp <<" "<< pt.xpSig  <<" "<< pt.errTot <<   endl;
           data.push_back(pt);
       }
       delete tuple;
    }

    outDir =  inFile(0, inFile.Last('/'));
    outDir += "/dPlots";

    gSystem->mkdir(outDir, true);
}


void dPlotter::readDataJets(TString inFile, vector<TString> samples)
{
    TFile *file = TFile::Open(inFile);

    hPars  = dynamic_cast<TH1D*>(file->Get("fitparameters"));
    hCorrs = dynamic_cast<TH2D*>(file->Get("fitcorrelations"));
    assert(hPars);
    assert(hCorrs);

    //Loop over all
    //vector<TString> dataNames = {"sample_H1incDDIS_HERA_I_LAr_cuts", "sample_H1incDDIS_HERA_I_SpacMB_cuts", "sample_H1incDDIS_HERA_I_SpacTrg_cuts"};
    for(auto n : samples) {
       TNtuple *tuple = (TNtuple*) file->Get("ASaveDataTheory/sample_" + n);

       double q2_max, q2_min, pt_min, pt_max;

       tuple->SetBranchAddress("q2_min",&q2_min);
       tuple->SetBranchAddress("q2_max",&q2_max);

       tuple->SetBranchAddress("pt_min",&pt_min);
       tuple->SetBranchAddress("pt_max",&pt_max);


       double Sig, errTot, th;
       tuple->SetBranchAddress("SigData",&Sig);
       tuple->SetBranchAddress("SigDataErr",&errTot);
       tuple->SetBranchAddress("SigTh",&th);
       //tuple->SetBranchAddress("SigThErr",&pt.thErr);


       const vector<double> bins  = {5.5, 7, 9, 15};


       Int_t nentries = (Int_t)tuple->GetEntries();
       for (Int_t i=0;i<nentries;i++) {
           tuple->GetEntry(i);
           //cout << "Hel " << pt.xp <<" "<< pt.xpSig  <<" "<< pt.errTot <<   endl;

          if(jetsData.count({q2_min,q2_max}) == 0) {
             jetsData[{q2_min,q2_max}] = new TH1D(rn(), "", bins.size()-1,bins.data());
          }
          if(jetsTh.count({q2_min,q2_max}) == 0) {
             jetsTh[{q2_min,q2_max}] = new TH1D(rn(), "", bins.size()-1,bins.data());
          }

          int binId = jetsData[{q2_min,q2_max}]->FindBin((pt_min+pt_max)/2);

          errTot *= Sig; //To absolute
          jetsData[{q2_min,q2_max}]->SetBinContent(binId, Sig);
          jetsData[{q2_min,q2_max}]->SetBinError(binId, errTot);
          jetsTh[{q2_min,q2_max}]->SetBinContent(binId, th);

       }
       delete tuple;
    }

    outDir =  inFile(0, inFile.Last('/'));
    outDir += "/dPlots";

    gSystem->mkdir(outDir, true);
}

void dPlotter::plotBetaMore(TString fTag, vector<double> xpomVec)
{
    vector<map<double, vector<pointInc>>> dataMap(xpomVec.size());

    for(int ixp = 0; ixp < xpomVec.size(); ++ixp) {
       for(pointInc &p : data)
           if(p.xp == xpomVec[ixp]) {
               dataMap[ixp][p.q2].push_back(p);
           }
    }

    int nRows = xpomVec.size();
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", (1+4)*200, (1+nRows)*200);
    SetLeftRight(0.5/(1+4) +0.05, 0.5/(1+4) -0.05);
    SetTopBottom(0.5/(1+nRows), 0.5/(1+nRows));

    DivideTransparent(group(1, 0, 4), group(1, 0, nRows));

    cout << "next " << dataMap.size() <<  endl;
    for(int ixp = 0; ixp < dataMap.size(); ++ixp) {
       int i = ixp*4;
       for(auto item : dataMap[ixp]) {
           double q2    = item.first;
           auto &points = item.second;
           can->cd(i+1);

           TGraphErrors *gData = new TGraphErrors(points.size());

           for(int j = 0; j < points.size(); ++j) {
               gData->SetPoint(j, points[j].beta, points[j].xpSig);

               double er = points[j].xpSig * points[j].errTot;
               gData->SetPointError(j, 0, er);
               //cout <<"ahoj  " << j<<" "<< points[j].beta <<" "<< points[j].xpSig << " " << endl;
           }

           map<double,double> zMin = { {0.0003, 0.02},
                                       {0.0005, 0.02},
                                       {0.001, 0.015},
                                       {0.003, 0.0015},
                                       {0.01, 0.0015},
                                       {0.03, 0.0011}
                                      };

           TH1D *hFr = new TH1D(rn(), "", 1, zMin.at(xpomVec[ixp]), 1);
           hFr->Draw("axis");
           

           //TGraph *gTh    = new TGraph(points.size());

           /*
           TGraph *gThOrgA = new TGraph(points.size());
           TGraph *gThOrgB = new TGraph(points.size());
           for(int j = 0; j < points.size(); ++j) {
               gThOrgA->SetPoint(j, points[j].beta, points[j].thOrgA);
               gThOrgB->SetPoint(j, points[j].beta, points[j].thOrgB);
               gTh->SetPoint(j, points[j].beta, points[j].th);
           }
           */


           TGraph *gTh = new TGraph();
           for(int j = 0; j < points.size(); ++j) {
              if(!points[j].isInside) continue;
              //cout << "RADEK " << points[j].q2 <<" "<< points[j].isInside << endl;
              gTh->SetPoint(gTh->GetN(), points[j].beta, points[j].th);
           }


           TGraph *gThOut = new TGraph();
           for(int j = 0; j < points.size(); ++j) {
              if(!points[j].isInside ||  (j<points.size()-1 && !points[j+1].isInside) ) 
                 gThOut->SetPoint(gThOut->GetN(), points[j].beta, points[j].th);
           }


           TString st = (gTh->GetN()+gThOut->GetN()) > 1 ? "l same" : "l* same";


           //MyTheory
           gTh->SetLineColor(kBlue);
           gTh->SetMarkerColor(kBlue);
           gTh->SetLineWidth(3);
           if(gTh->GetN() > 0)
              gTh->Draw(st);


           gThOut->SetLineColor(kBlue);
           gThOut->SetMarkerColor(kBlue);
           gThOut->SetLineWidth(3);
           gThOut->SetLineStyle(2);
           if(gThOut->GetN() > 0)
              gThOut->Draw(st);



           gData->SetMarkerColor(kRed);
           gData->SetLineColor(kRed);
           gData->SetMarkerStyle(20);
           gData->SetMarkerSize(1.5);
           gData->Draw("pe same");

           gPad->SetLogx();

           GetYaxis()->SetNdivisions(503);

           GetYaxis()->SetRangeUser(0, 0.085);

           SetFTO({30}, {14}, {1.55, 2.2, 0.4, 3.4});
           GetXaxis()->SetTitleSize(1.4*GetXaxis()->GetTitleSize());
           GetYaxis()->SetTitleSize(1.4*GetYaxis()->GetTitleSize());

           DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2), 28);
           //DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2), GetXaxis()->GetLabelSize());



           if(i < dataMap.size() - 4  && dataMap.size() > 4) {
               GetXaxis()->SetTickSize(0);
               GetXaxis()->SetLabelOffset(5000);
           }

           if(i % 4 != 0) {
               GetYaxis()->SetLabelOffset(5000);
               GetYaxis()->SetTickSize(0);
           }

           if(ixp != dataMap.size()-1 && dataMap[ixp+1].size() >= i+1) {
               GetXaxis()->SetTickSize(0);
               GetXaxis()->SetLabelOffset(5000);
           }



           if(i == 0)
               GetYaxis()->SetTitle("x_{IP} #sigma_{r}^{D(3)}");
           if(i == dataMap.size() -1)
               GetXaxis()->SetTitle("#beta");



           if(i == 0) {
               can->cd();
               TLegend *leg = new TLegend(0.355-0.205, 1-can->GetTopMargin()+0.01, 0.95-0.205, 1);
               leg->SetBorderSize(0);
               leg->SetTextSize(PxFontToRel(32));
               //leg->SetTextSize(0.035);
               leg->SetNColumns(2);
               leg->SetMargin(0.2);
               leg->AddEntry(gData, "H1 Data #sqrt{s} = "+ fTag(1,1000) + " GeV", "pe");
               //leg->AddEntry(gThOrgA,   "Fit A", "lp");
               //leg->AddEntry(gTh,       "Our Fit A", "lp");
               leg->AddEntry((TObject*)nullptr,       "", "");
               leg->AddEntry(gTh,       "H1 Fit 2019 NNLO prelim.", "lp");
               leg->AddEntry(gThOut,     "extrap.", "lp");


               //leg->AddEntry(gThOrgB,   "Fit B", "lp");
               leg->Draw();
           }



           ++i;
       }
    }

    for(int ixp = 0; ixp < xpomVec.size(); ++ixp) {
       int pos = dataMap[ixp].size();
       //DrawLatexRight(can->GetPad(ixp*4 + pos), 0.9, Form("x_{IP} = %g", xpomVec[ixp]), -1, "l");
       DrawLatex(can->GetPad(ixp*4 + pos), 1.1+(3-pos), 0.5, Form("x_{IP} = %g", xpomVec[ixp]), -1, "l");
    }


    can->cd();
    TLegend *leg = new TLegend(0.4, 1-can->GetTopMargin(), 0.9, 1);


    can->SaveAs(Form(outDir + "/%s_betaMore_xpom.pdf", fTag.Data()));
    //can->SaveAs(Form(outDir + "/%s_beta_xpom%g.pdf", fTag.Data(), xpom));

}








void dPlotter::plotBeta(TString fTag, double xpom)
{
    map<double, vector<pointInc>> dataMap;

    for(pointInc &p : data)
        if(p.xp == xpom) {
            dataMap[p.q2].push_back(p);
        }

    int nRows = ceil(dataMap.size() / 4.);
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", (1+4)*200, (1+nRows)*200);
    SetLeftRight(0.5/(1+4) +0.05, 0.5/(1+4) -0.05);
    SetTopBottom(0.5/(1+nRows), 0.5/(1+nRows));

    DivideTransparent(group(1, 0, 4), group(1, 0, nRows));

    cout << "next " << dataMap.size() <<  endl;
    int i = 0;
    for(auto item : dataMap) {
        double q2    = item.first;
        auto &points = item.second;
        can->cd(i+1);

        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].beta, points[j].xpSig);

            double er = points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
            //cout <<"ahoj  " << j<<" "<< points[j].beta <<" "<< points[j].xpSig << " " << endl;
        }

        map<double,double> zMin = { {0.0003, 0.02},
                                    {0.0005, 0.02},
                                    {0.001, 0.015},
                                    {0.003, 0.0015},
                                    {0.01, 0.0015},
                                    {0.03, 0.0011}
                                   };

        TH1D *hFr = new TH1D(rn(), "", 1, zMin[xpom], 1);
        hFr->Draw("axis");
        

        //TGraph *gTh    = new TGraph(points.size());

        /*
        TGraph *gThOrgA = new TGraph(points.size());
        TGraph *gThOrgB = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gThOrgA->SetPoint(j, points[j].beta, points[j].thOrgA);
            gThOrgB->SetPoint(j, points[j].beta, points[j].thOrgB);
            gTh->SetPoint(j, points[j].beta, points[j].th);
        }
        */


        TGraph *gTh = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           if(!points[j].isInside) continue;
           //cout << "RADEK " << points[j].q2 <<" "<< points[j].isInside << endl;
           gTh->SetPoint(gTh->GetN(), points[j].beta, points[j].th);
        }


        TGraph *gThOut = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           if(!points[j].isInside ||  (j<points.size()-1 && !points[j+1].isInside) ) 
              gThOut->SetPoint(gThOut->GetN(), points[j].beta, points[j].th);
        }


        TString st = (gTh->GetN()+gThOut->GetN()) > 1 ? "l same" : "l* same";


        //MyTheory
        gTh->SetLineColor(kBlue);
        gTh->SetMarkerColor(kBlue);
        gTh->SetLineWidth(3);
        if(gTh->GetN() > 0)
           gTh->Draw(st);


        gThOut->SetLineColor(kBlue);
        gThOut->SetMarkerColor(kBlue);
        gThOut->SetLineWidth(3);
        gThOut->SetLineStyle(2);
        if(gThOut->GetN() > 0)
           gThOut->Draw(st);



        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->SetMarkerSize(1.5);
        gData->Draw("pe same");

        gPad->SetLogx();

        GetYaxis()->SetNdivisions(503);

        GetYaxis()->SetRangeUser(0, 0.085);

        SetFTO({30}, {14}, {1.55, 2.2, 0.4, 3.4});
        GetXaxis()->SetTitleSize(1.4*GetXaxis()->GetTitleSize());
        GetYaxis()->SetTitleSize(1.4*GetYaxis()->GetTitleSize());

        DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2), 28);
        //DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2), GetXaxis()->GetLabelSize());



        if(i < dataMap.size() - 4  && dataMap.size() > 4) {
            GetXaxis()->SetTickSize(0);
            GetXaxis()->SetLabelOffset(5000);
        }

        if(i % 4 != 0) {
            GetYaxis()->SetLabelOffset(5000);
            GetYaxis()->SetTickSize(0);
        }

        if(i == 0)
            GetYaxis()->SetTitle("x_{IP} #sigma_{r}^{D(3)}");
        if(i == dataMap.size() -1)
            GetXaxis()->SetTitle("#beta");

        if(i == 0) {
            can->cd();
            TLegend *leg = new TLegend(0.365, 1-can->GetTopMargin()+0.01, 0.95, 1);
            leg->SetBorderSize(0);
            leg->SetMargin(0.2);
            leg->SetTextSize(PxFontToRel(32));
            //leg->SetTextSize(0.035);
            leg->SetNColumns(2);
            leg->AddEntry(gData, "H1 Data #sqrt{s} = 319 GeV", "pe");
            leg->AddEntry((TObject*)nullptr,       "", "");

            leg->AddEntry(gTh,       "H1 Fit 2019 NNLO prelim.", "lp");

            leg->AddEntry(gThOut,     "extrap.", "lp");


            //leg->AddEntry(gThOrgB,   "Fit B", "lp");
            leg->Draw();
        }


        ++i;
    }

    DrawLatexUp(can->GetPad(1), 0.9, Form("x_{IP} = %g", xpom), -1, "l");

    can->cd();
    TLegend *leg = new TLegend(0.4, 1-can->GetTopMargin(), 0.9, 1);


    can->SaveAs(Form(outDir + "/%s_beta_xpom%g.pdf", fTag.Data(), xpom));

}

void dPlotter::plotBetaRat(TString fTag, double xpom)
{
    map<double, vector<pointInc>> dataMap;

    for(pointInc &p : data)
        if(p.xp == xpom) {
            dataMap[p.q2].push_back(p);
        }

    int nRows = ceil(dataMap.size() / 4.);
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", (1+4)*200, (1+nRows)*200);
    SetLeftRight(0.5/(1+4) +0.05, 0.5/(1+4) -0.05);
    SetTopBottom(0.5/(1+nRows), 0.5/(1+nRows));

    cout << "Middle" << endl;
    DivideTransparent(group(1, 0, 4), group(1, 0, nRows));

    cout << "next " << dataMap.size() <<  endl;
    int i = 0;
    for(auto item : dataMap) {
        double q2    = item.first;
        auto &points = item.second;
        can->cd(i+1);

        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            //gData->SetPoint(j, points[j].beta, points[j].xpSig);
            gData->SetPoint(j, points[j].beta, 1);

            //double er = points[j].xpSig * points[j].errTot;
            //gData->SetPointError(j, 0, er);
            gData->SetPointError(j, 0, points[j].errTot);
            //cout <<"ahoj  " << j<<" "<< points[j].beta <<" "<< points[j].xpSig << " " << endl;
        }

        map<double,double> zMin = { {0.0003, 0.02},
                                    {0.0005, 0.02},
                                    {0.001, 0.015},
                                    {0.003, 0.0015},
                                    {0.01, 0.0015},
                                    {0.03, 0.0011}
                                   };




        TH1D *hFr = new TH1D(rn(), "", 1, zMin[xpom], 1);
        hFr->Draw("axis");
        

        //TGraph *gTh    = new TGraph(points.size());

        /*
        TGraph *gThOrgA = new TGraph(points.size());
        TGraph *gThOrgB = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gThOrgA->SetPoint(j, points[j].beta, points[j].thOrgA);
            gThOrgB->SetPoint(j, points[j].beta, points[j].thOrgB);
            gTh->SetPoint(j, points[j].beta, points[j].th);
        }
        */


        TGraph *gTh = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           if(!points[j].isInside) continue;
           //cout << "RADEK " << points[j].q2 <<" "<< points[j].isInside << endl;
           gTh->SetPoint(gTh->GetN(), points[j].beta, points[j].th / points[j].xpSig);
        }


        TGraph *gThOut = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           if(!points[j].isInside ||  (j<points.size()-1 && !points[j+1].isInside) ) 
              gThOut->SetPoint(gThOut->GetN(), points[j].beta, points[j].th / points[j].xpSig);
        }


        TString st = (gTh->GetN()+gThOut->GetN()) > 1 ? "l same" : "l* same";


        //MyTheory
        gTh->SetLineColor(kBlue);
        gTh->SetMarkerColor(kBlue);
        gTh->SetLineWidth(2);
        if(gTh->GetN() > 0)
           gTh->Draw(st);


        gThOut->SetLineColor(kBlue);
        gThOut->SetMarkerColor(kBlue);
        gThOut->SetLineWidth(2);
        gThOut->SetLineStyle(2);
        if(gThOut->GetN() > 0)
           gThOut->Draw(st);



        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");

        gPad->SetLogx();

        GetYaxis()->SetNdivisions(303);

        double eps = 0.001;
        GetYaxis()->SetRangeUser(0.43+eps, 1.57 -eps);

        SetFTO({30}, {14}, {1.55, 1.9, 0.4, 2.7});
        GetXaxis()->SetTitleSize(1.4*GetXaxis()->GetTitleSize());
        GetYaxis()->SetTitleSize(1.4*GetYaxis()->GetTitleSize());

        DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2), 28);



        if(i < dataMap.size() - 4  && dataMap.size() > 4) {
            GetXaxis()->SetTickSize(0);
            GetXaxis()->SetLabelOffset(5000);
        }

        if(i % 4 != 0) {
            GetYaxis()->SetLabelOffset(5000);
            GetYaxis()->SetTickSize(0);
        }

        if(i == 0)
            GetYaxis()->SetTitle("#sigma/#sigma^{data}");
        if(i == dataMap.size() -1)
            GetXaxis()->SetTitle("#beta");

        if(i == 0) {
            can->cd();
            TLegend *leg = new TLegend(0.4, 1-can->GetTopMargin()+0.01, 0.9, 1);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);
            leg->AddEntry(gData, "H1 Data", "pe");
            //leg->AddEntry(gThOrgA,   "Fit A", "lp");
            //leg->AddEntry(gTh,       "Our Fit A", "lp");
            leg->AddEntry(gTh,       "H1 Fit 2019 NNLO", "lp");
            leg->AddEntry((TObject*)nullptr,       "", "");
            leg->AddEntry(gThOut,     "(extrapol. fit)", "lp");
            //leg->AddEntry(gThOrgB,   "Fit B", "lp");
            leg->Draw();
        }


        ++i;
    }

    DrawLatexUp(can->GetPad(1), 1.1, Form("x_{IP} = %g", xpom), -1, "l");

    can->cd();
    TLegend *leg = new TLegend(0.4, 1-can->GetTopMargin(), 0.9, 1);


    can->SaveAs(Form(outDir + "/%s_betaRat_xpom%g.pdf", fTag.Data(), xpom));

}





//To print the curve tag
pair<double,double> getRight(TGraphErrors *gr)
{
    pair<double,double> pRight(-100,0);
    for(int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        if(pRight.first < x) {
            pRight = {x,y};
        }
    }
    return pRight;
}


void dPlotter::plotQ2(TString fTag, double xpom)
{
    map<double, vector<pointInc>> dataMap;

    for(pointInc &p : data)
        if(p.xp == xpom)
            dataMap[p.beta].push_back(p);

    int nPoints = dataMap.size();

    gStyle->SetOptStat(0);
    int ySize = xpom > 0.002 ? 800 : 500;
    TCanvas *can = new TCanvas(rn(),"", 500, ySize);
    SetLeftRight(0.15, 0.07);
    SetTopBottom(0.08, 0.12);

    gPad->SetLogx();
    gPad->SetLogy();


    /*
       map<double,double> zMin = { {0.0003, 0.02},
       {0.001, 0.015},
       {0.003, 0.0015},
       {0.01, 0.0015},
       {0.03, 0.0015}
       };
       */


    //Plotting style
    TH1D *hFr = new TH1D(rn(), "", 1, 2, 1.5e4);
    hFr->Draw("axis");

    SetFTO({20}, {14}, {1.4, 2.2, 0.4, 2.8});
    GetXaxis()->SetTitleSize(1.3*GetXaxis()->GetTitleSize());
    GetYaxis()->SetTitleSize(1.3*GetYaxis()->GetTitleSize());

    map<double,vector<double>> yMinMax = {
             {0.0003, {0.04, 1}},
             {0.001,  {2e-2, 30}},
             {0.003,  {2e-2, 2e2}},
             {0.01,   {1e-2, 1e4}},
             {0.03,   {1e-2, 2e5}},
         };




    GetYaxis()->SetRangeUser(yMinMax.at(xpom)[0], yMinMax.at(xpom)[1]);
    GetYaxis()->SetTitle("3^{i} * x_{IP} #sigma_{r}^{D(3)}");
    GetXaxis()->SetTitle("Q^{2} [GeV]");



    int i = 0;
    for(auto item : dataMap) {
        double beta  = item.first;
        auto &points = item.second;

        int idx = dataMap.size() - 1 - i;
        double fac = pow(3, idx);
        TGraphErrors *gData = new TGraphErrors(points.size());


        sort(points.begin(), points.end(), [](const pointInc &a, const pointInc &b) {return a.q2 < b.q2;});

        //Check that it's rising 

        for(int k = 0; k < points.size() - 1; ++k) {
            if(points[k].q2 == points[k+1].q2)
               cout << "Big problem " << endl;
        }


        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].q2, fac*points[j].xpSig);

            double er = fac*points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
        }
        

        TGraph *gTh = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           if(!points[j].isInside) continue;
           //cout << "RADEK " << points[j].q2 <<" "<< points[j].isInside << endl;
           gTh->SetPoint(gTh->GetN(), points[j].q2, fac*points[j].th);
        }

        TGraph *gThOrgA = new TGraph();
        TGraph *gThOrgB = new TGraph();
        for(int j = 0; j < points.size(); ++j) {
           gThOrgA->SetPoint(gThOrgA->GetN(), points[j].q2, fac*points[j].thOrgA);
           gThOrgB->SetPoint(gThOrgB->GetN(), points[j].q2, fac*points[j].thOrgB);
        }



        TGraph *gThOut = new TGraph();
        bool lastIn = false;
        for(int j = 0; j < points.size(); ++j) {
           if(points[j].isInside && (lastIn || gThOut->GetN() == 0)) continue;
           gThOut->SetPoint(gThOut->GetN(), points[j].q2, fac*points[j].th);
           if(points[j].isInside) lastIn = true;
        }


        /*
        gThOrgA->SetLineColor(kBlue);
        gThOrgA->SetMarkerColor(kBlue);
        gThOrgA->Draw("l* same");
        */

        TString st = (gTh->GetN()+gThOut->GetN()) > 1 ? "l same" : "l* same";

        /*
        gThOrgB->SetLineColor(kGreen);
        gThOrgB->SetMarkerColor(kGreen);
        gThOrgB->SetLineStyle(1);
        if(gThOrgB->GetN() > 0)
           gThOrgB->Draw(st);
        */


        gTh->SetLineColor(kBlue);
        gTh->SetMarkerColor(kBlue);
        gTh->SetLineWidth(2);
        if(gTh->GetN() > 0)
           gTh->Draw(st);

        gThOut->SetLineColor(kBlue);
        gThOut->SetMarkerColor(kBlue);
        gThOut->SetLineStyle(2);
        gThOut->SetLineWidth(2);
        if(gThOut->GetN() > 0)
           gThOut->Draw(st);


        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");



        double xT, yT;
        tie(xT,yT) =  getRight(gData);
        TLatex *lat = new TLatex;
        lat->SetTextSize(PxFontToRel(16));
        lat->SetTextAlign(12);
        if(idx == 13) yT *= 1.3; //Hack for the highest point
        lat->DrawLatex(xT*1.3, yT, Form("#beta=%g (i=%lu)", beta, idx));


        if(i == 0) {
            can->cd();

            double legSize = xpom < 0.002 ? 0.10 : 0;

            TLegend *leg = new TLegend(0.50, 1-can->GetTopMargin()-0.15 - legSize, 0.9, 1-can->GetTopMargin()-0.01);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(1.0*GetYaxis()->GetLabelSize());
            leg->SetNColumns(1);
            leg->AddEntry(gData, "H1 Data", "pe");
            leg->AddEntry((TObject*)nullptr, "#sqrt{s} = 319 GeV", "");
            //leg->AddEntry(gThOrgA,   "Fit A", "lp");
            //leg->AddEntry(gTh,       "Our Fit A", "lp");
            leg->AddEntry(gTh,       "H1 Fit 2019 NNLO", "l");
            leg->AddEntry((TObject*)nullptr,   "      preliminary", "");
            leg->AddEntry(gThOut,    "extrapolation", "l");

            //leg->AddEntry(gThOrgB,   "Fit B", "lp");
            leg->Draw();
        }


        ++i;
    }

    DrawLatexUp(0.9, Form("x_{IP} = %g",xpom), -1, "l");


    can->SaveAs(Form(outDir +  "/%s_q2_xpom%g.pdf", fTag.Data(), xpom));

}

void dPlotter::plotJets()
{
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 800);
    SetLeftRight(0.1, 0.1);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0, 2), group(1, 0, 3));

    int i = 0;
    for(auto data : jetsData) {
       double q2Min, q2Max;
       tie(q2Min,q2Max) = data.first;
       TH1D *hData = data.second;
       TH1D *hTh   = jetsTh[data.first];

       can->cd(i+1);
       gPad->SetLogy();
       hData->Draw("e");
       hTh->Draw("same");
       ++i;
    }
    can->SaveAs(outDir + "/jets.pdf");
}

void plotJets(dPlotter &nlo, dPlotter &nnlo)
{
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 800);
    SetLeftRight(0.15, 0.05);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0, 2), group(1, 0, 3));

    int i = 0;
    for(auto data : nlo.jetsData) {
       double q2Min, q2Max;
       tie(q2Min,q2Max) = data.first;
       TH1D *hData = data.second;
       TH1D *hNlo    = nlo.jetsTh[data.first];
       TH1D *hNNlo   = nnlo.jetsTh[data.first];

       can->cd(i+1);
       gPad->SetLogy();
       hData->SetMarkerColor(kBlack);
       hData->SetMarkerSize(0.8);
       hData->SetMarkerStyle(20);
       hData->Draw("P X0 e");
       hData->SetLineColor(kBlack);

       hNlo->Draw("same");
       hNNlo->Draw("same");
       hNNlo->SetLineColor(kRed);

       if(i == 4)
          GetXaxis()->SetTitle("p^{*}_{T,1} [GeV]");
       if(i == 0)
          GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp^{*}_{T,1}dQ^{2}} [pb/GeV^{3}]");
       GetYaxis()->SetNdivisions(303);
       GetXaxis()->SetNdivisions(303);
       SetFTO({20}, {10}, {1.4, 2.2, 0.4, 3.3});

        if(i < nlo.jetsData.size() - 2  && nlo.jetsData.size() > 2) {
            GetXaxis()->SetTickSize(0);
            GetXaxis()->SetLabelOffset(5000);
        }
        if(i % 2 != 0) {
            GetYaxis()->SetTickSize(0);
            GetYaxis()->SetLabelOffset(5000);
        }


        /*
        if(i == 1) {
           TLegend *leg = newLegend(kPos9);
           leg->AddEntry(hData, "H1 HERAII jet data", "le");
           leg->AddEntry(hNlo, "NLO prediction");
           leg->AddEntry(hNNlo, "NNLO prediction");
           DrawLegends({leg});
        }
        */

       DrawLatexUp(-1, Form("%g < Q^{2} < %g GeV^{2}",q2Min, q2Max), -1, "c");

       if(i == 4) {
          can->cd();
          TLegend *leg = new TLegend(0.6, 0.2, 0.95, 0.3);
          leg->SetTextSize(0.03);
          leg->SetBorderSize(0);
          leg->AddEntry(hData, "H1 HERAII jet data", "pe");
          leg->AddEntry(hNlo, "NLO prediction", "l");
          leg->AddEntry(hNNlo, "NNLO prediction", "l");
          leg->Draw();
       }




       ++i;
    }




    can->SaveAs("jetsNLOvsNNLO.pdf");
}


//The xpom plot in the H1 style for FPS 4D data
void dPlotter::plotXpomT()
{
    gStyle->SetOptStat(0);

    //for(auto &p : pars) p = smear(p);

    map<pair<double,double>, vector<pointInc>> pMap;

    for(auto &p : data) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    map<pair<double,double>,map<double,TGraphErrors*>> grDataMap;
    map<pair<double,double>,map<double,TGraph*>> grThMap;


    double chi2 = 0;
    int nDf = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        double q2, beta;
        tie(beta,q2) = pm.first;
        
        if(ps.size() < 2) continue;

        map<double,TGraphErrors*> grData;
        map<double,TGraph*> grTh;

        //Loop over points with equal beta,q2 for data
        for(auto &p : ps) {
            if(grData.count(p.tAbs) == 0) grData[p.tAbs] = new TGraphErrors();
            grData.at(p.tAbs)->SetPoint(grData.at(p.tAbs)->GetN(), p.xp, p.xpSig);
            grData.at(p.tAbs)->SetPointError(grData.at(p.tAbs)->GetN()-1, 0,  p.xpSig*p.errTot);

            if(grTh.count(p.tAbs) == 0) grTh[p.tAbs] = new TGraph();
            grTh.at(p.tAbs)->SetPoint(grTh.at(p.tAbs)->GetN(), p.xp, p.th);

        }
        assert(grDataMap.count(pm.first)==0);
        grDataMap[pm.first] = grData;
        grThMap[pm.first] = grTh;
    }


    vector<double> q2Set, betaSet;
    for(auto pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        if(ps.size() < 2) continue;

        betaSet.push_back(pm.first.first);
        q2Set.push_back(pm.first.second);
    }

    auto clean = [](vector<double> &v) {
        std::sort(v.begin(), v.end()); // 1 1 2 2 3 3 3 4 4 5 5 6 7 
        auto last = std::unique(v.begin(), v.end());
        // v now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
        v.erase(last, v.end());
    };
    clean(betaSet);
    clean(q2Set);

    cout << "Sizes " << q2Set.size() <<" "<< betaSet.size() << endl;
    auto can = new TCanvas("can","", 600, 670);
    SetLeftRight(0.15, 0.05);
    double f = 0.28358;
    SetTopBottom(0.12/0.20*f, 0.08/0.20*f);

    //DividePad(vector<double>(q2Set.size(),1), vector<double>(betaSet.size(),1));
    DivideTransparent(group(1,0,q2Set.size()), group(1,0,betaSet.size()));

    int iq2 = 0, ibeta = 0;
    double chi2tot = 0;
    int ndf = 0;

    for(int iq2 = 0; iq2 < q2Set.size(); ++iq2)
    for(int ibeta = 0; ibeta < betaSet.size(); ++ibeta) {
        can->cd(iq2*q2Set.size() + ibeta + 1);

        TH1D *hAx = new TH1D(rn(), "", 1, 0.0013, 0.15);
        hAx->Draw("axis");

        gPad->SetLogx();
        GetYaxis()->SetRangeUser(-0.01, 0.12);

        int ifit = 2;
        double beta = betaSet[ibeta];
        double q2 = q2Set[iq2];


        if(grDataMap.count({beta, q2}) > 0) {
            auto gr   = grDataMap.at({beta, q2});
            auto grTh = grThMap.at({beta, q2});

            /*
            vector<double> aVec, data, regg, err;
            //Fill data and theory points
            for(auto p : pMap.at({beta,q2})) { //points with equal beta, q2



                //Subtract Reggeon from fitted data
                double xPq[13], f2[2], fl[2], c2[2], cl[2];
                double z = p.beta;
                double q2 = p.q2;
                qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
                ifit = 0;
                //const double bNorm = fluxI(0.003, a0_R, ap_R/b0_R);

                const double a0_R = 0.5;
                const double ap_R = 0.3;
                const double b0_R = 1.6;

                const double a0_P = 1.11101;
                const double ap_P = 0.06;
                const double b0_P = 5.5;

                double bFluxFitB = flux(p.xp,p.tAbs, a0_R, ap_R, b0_R);
                double aFluxFitB = flux(p.xp,p.tAbs, a0_P, ap_P, b0_P);
                //p.sigma -= f2[1]*bFlux*pars[3];

                const double s = 319*319;
                double y = q2 / (s*p.xp*p.beta);
                //double FL = y*y / (1+pow(1-y,2))*fl[0];
                double flC = y*y / (1+pow(1-y,2));

                double corr = (f2[1]-flC*fl[1])*bFluxFitB*pars[3] - flC*fl[0]*aFluxFitB;

                aVec.push_back(aFlux);
                regg.push_back(corr);
                //cout << "Radek " << p.sigma <<" "<< p.err/p.sigma << endl;
                double errNow = smear(p.err/p.sigma) * p.sigma;
                data.push_back(smear(p.sigma) - corr);

                data.push_back(p.xpSig);
                err.push_back(p.errTot);
            }
            double norm;
            double chiNow =  chi2Raw(data, err,  aVec, &norm);
            chi2tot += chiNow;
            ndf += data.size() -1;

            map<double, TGraph*> grTh;
            //Fill theory plot
            int i = 0;
            for(auto p : pMap.at({beta,q2})) { //points with equal beta, q2
                //aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
                if(grTh.count(p.tAbs) == 0) grTh[p.tAbs] = new TGraph;
                int idx = grTh.at(p.tAbs)->GetN();
                grTh.at(p.tAbs)->SetPoint(idx, p.xp, aVec[i]*norm + regg[i]);
                ++i;
            }
            */    



            for(auto g : gr) {
                g.second->Draw("*e same");
                if(g.first == 0.2) g.second->SetLineColor(kBlue);
                else if(g.first == 0.4) g.second->SetLineColor(kRed);
                else if(g.first == 0.6) g.second->SetLineColor(kBlack);
                if(g.first == 0.2) g.second->SetMarkerColor(kBlue);
                else if(g.first == 0.4) g.second->SetMarkerColor(kRed);
                else if(g.first == 0.6) g.second->SetMarkerColor(kBlack);

                if(g.first == 0.2) g.second->SetMarkerStyle(20);
                else if(g.first == 0.4) g.second->SetMarkerStyle(24);
                else if(g.first == 0.6) g.second->SetMarkerStyle(22);

                g.second->SetMarkerSize(0.7);
                grTh.at(g.first)->SetLineColor(g.second->GetLineColor());
                grTh.at(g.first)->Draw("l same");
            }


            //Legend
            if(iq2 == q2Set.size()-1 && ibeta == betaSet.size()-1) {
                auto back = gPad;
                can->cd();
                auto leg = new TLegend(0.15, 0.835, 0.25, 0.935);        
                leg->SetTextSize(PxFontToRel(20));
                leg->SetBorderSize(0);
                leg->AddEntry(gr.at(0.2), "|t| = 0.2 GeV^{-2}", "p");
                leg->AddEntry(gr.at(0.4), "|t| = 0.4 GeV^{-2}", "p");
                leg->AddEntry(gr.at(0.6), "|t| = 0.6 GeV^{-2}", "p");
                leg->Draw();

                auto leg2 = new TLegend(0.55, 0.835, 0.65, 0.935);        
                leg2->SetTextSize(PxFontToRel(20));
                leg2->SetBorderSize(0);
                leg2->SetHeader("H1 FPS");
                leg2->AddEntry(grTh.at(0.6), "H1 Fit 2020 NNLO prelim.", "l");
                leg2->Draw();


                back->cd();
            }
        }


        if(iq2 == q2Set.size()-1 && ibeta == betaSet.size()-1) {
            GetXaxis()->SetTitle("x_{IP}");
        }
        if(iq2 == 0 && ibeta == 0) {
            GetYaxis()->SetTitle("x_{IP}#sigma_{r}^{D(4)} (GeV^{-2})");
        }


        double fSize = 16;
        SetFTO({fSize}, {7}, {1.4, 1.7, 0.3, 2.7});

        GetXaxis()->SetTitleSize(PxFontToRel(26));
        GetYaxis()->SetTitleSize(PxFontToRel(26));

        GetXaxis()->SetNdivisions(303);
        GetYaxis()->SetNdivisions(603, kTRUE);

        if(iq2 != q2Set.size()-1) GetXaxis()->SetLabelSize(0);

        if(iq2 == 0)   DrawLatexUp(-1, Form("#beta=%g", beta), fSize);
        if(ibeta == 0) {
            if(iq2 == 0) DrawLatexUp(-2, Form("%g GeV^{2}", q2), fSize);
            else DrawLatexUp(-1, Form("%g GeV^{2}", q2), fSize);
        }



    }
    cout << "H1 result: chi2= "<< chi2tot << " / " << ndf << endl;

    can->SaveAs("xpomRegge.pdf");

}


void dPlotter::plotBslope()
{
    gStyle->SetOptStat(0);

    map<pair<double,double>, vector<pointInc>> pMap;

    for(auto &p : data) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    map<pair<double,double>,TGraphErrors*> grSlopes, grSlopesTh;

    double chi2 = 0;
    int nDf = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        double q2, beta;
        tie(beta,q2) = pm.first;
        
        map<double, TGraphErrors*> sigmaData, sigmaTh;

        //Loop over points with equal beta,q2 for data
        for(auto &p : ps) {

            if(sigmaData.count(p.xp) == 0) sigmaData[p.xp] = new TGraphErrors();
            if(sigmaTh.count(p.xp) == 0)   sigmaTh[p.xp] = new TGraphErrors();

            sigmaData.at(p.xp)->SetPoint(sigmaData.at(p.xp)->GetN(), p.tAbs, p.xpSig);
            sigmaData.at(p.xp)->SetPointError(sigmaData.at(p.xp)->GetN()-1, 0, p.xpSig*p.errTot);

            sigmaTh.at(p.xp)->SetPoint(sigmaTh.at(p.xp)->GetN(), p.tAbs, p.th);
            sigmaTh.at(p.xp)->SetPointError(sigmaTh.at(p.xp)->GetN()-1, 0, p.th*p.errTot);

        }


        for(auto gr : sigmaData) {
            TFitResultPtr res = gr.second->Fit("expo", "ES");
            double slope = -res->Parameter(1); 
            double err   = res->ParError(1);

            if(grSlopes.count(pm.first) == 0)
                grSlopes[pm.first] = new TGraphErrors();

            grSlopes.at(pm.first)->SetPoint(grSlopes.at(pm.first)->GetN(), gr.first, slope);
            grSlopes.at(pm.first)->SetPointError(grSlopes.at(pm.first)->GetN()-1, 0, err);
        }

        for(auto gr : sigmaTh) {
            TFitResultPtr res = gr.second->Fit("expo", "ES");
            double slope = -res->Parameter(1); 
            double err   = res->ParError(1);

            if(grSlopesTh.count(pm.first) == 0)
                grSlopesTh[pm.first] = new TGraphErrors();

            grSlopesTh.at(pm.first)->SetPoint(grSlopesTh.at(pm.first)->GetN(), gr.first, slope);
            grSlopesTh.at(pm.first)->SetPointError(grSlopesTh.at(pm.first)->GetN()-1, 0, 0*err);
        }
    }


    auto can = new TCanvas("can","", 600, 670);
    SetLeftRight(0.15, 0.05);
    double f = 0.28358;
    SetTopBottom(0.12/0.20*f, 0.08/0.20*f);



    //set<double> q2Set, betaSet;
    //for(auto s : grSlopes) {
    //    betaSet.insert(s.first.first);
    //    q2Set.insert(s.first.second);
    //}

    vector<double> q2Set, betaSet;
    for(auto pm : grSlopes) {
        auto ps = pm.second; //points with equal beta, q2
        //if(ps.size() < 2) continue;

        betaSet.push_back(pm.first.first);
        q2Set.push_back(pm.first.second);
    }

    auto clean = [](vector<double> &v) {
        std::sort(v.begin(), v.end()); // 1 1 2 2 3 3 3 4 4 5 5 6 7 
        auto last = std::unique(v.begin(), v.end());
        // v now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
        v.erase(last, v.end());
    };
    clean(betaSet);
    clean(q2Set);





    cout << "Sizes " << q2Set.size() <<" "<< betaSet.size() << endl;

    //DividePad(vector<double>(q2Set.size(),1), vector<double>(betaSet.size(),1));
    DivideTransparent(group(1,0,q2Set.size()), group(1,0,betaSet.size()));

    //int iq2 = 0, ibeta = 0;
    //for(auto q2 : q2Set) {
    //    ibeta = 0;
    //    for(auto beta : betaSet) {
    for(int iq2 = 0; iq2 < q2Set.size(); ++iq2)
    for(int ibeta = 0; ibeta < betaSet.size(); ++ibeta) {

            can->cd(iq2*q2Set.size() + ibeta + 1);

            TH1D *hAx = new TH1D(rn(), "", 1, 0.0013, 0.15);
            hAx->Draw("axis");

            gPad->SetLogx();
            GetYaxis()->SetRangeUser(0.0001, 9.9999);

            int ifit = 2;
            double beta = betaSet[ibeta];
            double q2 = q2Set[iq2];


            if(grSlopes.count({beta,q2}) >= 1) {
                //can->cd(iq2*betaSet.size() + ibeta + 1);
                auto gr = grSlopes.at({beta,q2});
                auto grTh = grSlopesTh.at({beta,q2});

                cout << "RADEK " << ibeta <<" "<< iq2 <<" "<< endl;

                //TH1D *h = new TH1D(rn(), "", 1, 1e-3, 1e-1);
                //h->Draw("axis");

                gr->SetLineColor(kRed);
                gr->SetMarkerColor(kRed);
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.7);

                gr->Draw("pe same");
                grTh->SetLineColor(kBlack);
                grTh->SetMarkerColor(kBlack);
                grTh->SetLineWidth(2);
                if(grTh->GetN() > 1)
                    grTh->Draw("l same");
                else
                    grTh->Draw("* same");

                //gPad->SetLogx();


            }

            //hAx->SetMinimum(0);
            //hAx->SetMaximum(8);


            if(iq2 == q2Set.size()-1 && ibeta == betaSet.size()-1) {
                GetXaxis()->SetTitle("x_{IP}");
            }
            if(iq2 == 0 && ibeta == 0) {
                GetYaxis()->SetTitle("B (GeV^{-2})");
            }




            double fSize = 16;
            SetFTO({fSize}, {7}, {1.45, 1.7, 0.3, 2.0});

            GetXaxis()->SetTitleSize(PxFontToRel(26));
            GetYaxis()->SetTitleSize(PxFontToRel(26));

            GetXaxis()->SetNdivisions(303);
            GetYaxis()->SetNdivisions(505);

            if(iq2 != q2Set.size()-1) GetXaxis()->SetLabelSize(0);
            if(ibeta != 0) GetYaxis()->SetLabelOffset(11000);

            if(iq2 == 0)   DrawLatexUp(-1, Form("#beta=%g", beta), fSize);
            if(ibeta == 0) {
                if(iq2 == 0) DrawLatexUp(-2, Form("%g GeV^{2}", q2), fSize);
                else DrawLatexUp(-1, Form("%g GeV^{2}", q2), fSize);
            }

            //Legend
            if(iq2 == q2Set.size()-1 && ibeta == betaSet.size()-1) {

                auto gr = grSlopes.at({beta,q2});
                auto grTh = grSlopesTh.at({beta,q2});

                auto back = gPad;
                can->cd();
                auto leg = new TLegend(0.15, 0.805, 0.25, 0.905);        
                leg->SetTextSize(PxFontToRel(20));
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(gr, "H1 FPS HERA II", "p");
                leg->Draw();

                auto leg2 = new TLegend(0.55, 0.805, 0.65, 0.905);        
                leg2->SetTextSize(PxFontToRel(20));
                leg2->SetBorderSize(0);
                leg2->SetFillStyle(0);
                //leg2->SetHeader("H1 FPS");
                leg2->AddEntry(grTh, "H1 Fit 2020 NNLO prelim.", "l");
                leg2->Draw();

                back->cd();
            }





    }

    /*

    map<double,double> sumSlopes, sumWgt;
    //Construct weight average
    for(auto s : grSlopes) {
        auto gr = s.second;
        for(int i = 0; i < gr->GetN(); ++i) {
            double xpom, sl;
            gr->GetPoint(i, xpom, sl);
            double err2 = pow(gr->GetErrorY(i),2);
            sumSlopes[xpom] += sl/err2;
            sumWgt[xpom] += 1./err2;
        }
    }


    TGraphErrors *gr = new TGraphErrors();

    for(auto sl: sumSlopes) {
        double xp = sl.first; 
        double sumSl = sl.second;
        double sumW = sumWgt.at(xp);

        double err = 1./sqrt(sumW);

        double rat = sumSl / sumW;
        cout <<"Slopes " <<  xp <<" "<< rat <<" "<<  err <<    endl;

        gr->SetPoint(gr->GetN(), log(1./xp), rat);
        gr->SetPointError(gr->GetN()-1, 0,  err);
    }

    auto dan = new TCanvas("dan","", 600, 600);

    gr->SetMarkerStyle(20);
    gr->Draw("ae*");
    gr->Fit("pol1");
    */

    can->SaveAs("slopeRegge.pdf");


}











void dPlotter::plotXpom()
{
    map<pair<double,double>, vector<pointInc>> dataMap;

    //vector<double> betas = {0.01, 0.04, 0.1, 0.2, 0.4, 0.65, 0.9};
    vector<double> betas = { 0.011, 0.043, 0.11, 0.2, 0.43, 0.67, 0.8};

    vector<double> q2s   = {3.5, 5, 6.5, 8.5, 12, 15, 20, 25, 35, 45, 60, 90, 200, 400, 800, 1600};

    map<double, set<double>> q2Beta;
    for(pointInc &p : data) {
        dataMap[{p.q2, p.beta}].push_back(p);
        q2Beta[p.q2].insert(p.beta);
    }

    int maxBetaN=0;
    for(auto it : q2Beta) {
        maxBetaN = max(maxBetaN, int(it.second.size()));
        cout << "q2 " << it.first << " " << it.second.size() << endl;
        for(auto b : it.second)
            cout <<"beta " <<  b << " " << endl;
    }


    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 800);
    SetLeftRight(0.1, 0.1);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0,  maxBetaN), group(1, 0, q2Beta.size()));

    int iq2 = 0, ibeta = 0;
    for(auto qq : q2Beta) {
        double q2   = qq.first;
        for(auto beta  : qq.second) {
            can->cd(iq2*maxBetaN + ibeta + 1);

            TH1D *hFr = new TH1D(rn(), "", 1, 1e-4, 4e-2);
            hFr->Draw("axis");
            GetYaxis()->SetRangeUser(0, 0.085);

            gPad->SetLogx();
            GetYaxis()->SetNdivisions(303);
            SetFTO({14}, {6}, {1.4, 2.2, 0.4, 3.9});


            /*
            if(iq != q2s.size() -1)
                GetXaxis()->SetTickSize(0);
            if(ib != 0)
                GetYaxis()->SetTickSize(0);
            */


            if(!dataMap.count({q2,beta})) {
                cout << "Not existing " << q2 <<" "<< beta << endl;
                continue;
            }
            else {
                cout << "Yes existing " << q2 <<" "<< beta << endl;
            }
            auto &points = dataMap.at({q2,beta});


            TGraphErrors *gData = new TGraphErrors(points.size());

            for(int j = 0; j < points.size(); ++j) {
                gData->SetPoint(j, points[j].xp, points[j].xpSig);

                double er = points[j].xpSig * points[j].errTot;
                gData->SetPointError(j, 0, er);
            }



            gData->SetMarkerColor(kRed);
            gData->SetLineColor(kRed);
            gData->SetMarkerSize(0.3);
            gData->SetMarkerStyle(20);
            gData->Draw("pe same");

            TGraph *gTh = new TGraph(points.size());
            for(int j = 0; j < points.size(); ++j) {
                gTh->SetPoint(j, points[j].xp, points[j].th);
            }
            gTh->SetLineColor(kBlue);
            gTh->SetMarkerColor(kBlue);
            gTh->Draw("*l same");




            DrawLatex(0.5, 0.5, Form("%zu", points.size()));

            ++ibeta;

        }
        ++iq2;
    }

    //DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));


    can->SaveAs(outDir + "/xpomGrid.pdf");

}

void dPlotter::plotParameters()
{
    gStyle->SetOptStat(0);

    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.13, 0.16);

    map<TString, vector<double>> pars;
    //Parameters  for fitA
    pars["g0"]  =  {  0.14591    ,   0.33171E-01   };
    pars["g2"]  =   { -0.94705    ,   0.20309       };
    pars["s0"] = {  1.0587     ,   0.322116   };
    pars["s1"] = {   2.2964    ,   0.36439       };
    pars["s2"] = {  0.56894    ,   0.14969       };
    pars["n_IR"] =     {  0.16966E-02,   0.41732E-03   };
    pars["a0_IP"] =   {   1.1182    ,   0.81319E-02 };

    TH1D *hParsFitA = (TH1D*) hPars->Clone("FitA");
    for(int i = 1; i <= hPars->GetNbinsX(); ++i) {
        TString s = hPars->GetXaxis()->GetBinLabel(i);
        int nFound = 0;
        for(auto &p : pars) {
            if(s.Contains(p.first)) {
               hParsFitA->SetBinContent(i, p.second[0]);
               hParsFitA->SetBinError(i, p.second[1]);
               ++nFound;
            }
        }
        assert(nFound == 1);
    }
    hParsFitA->SetLineColor(kRed);
    hParsFitA->SetLineStyle(2);

    hPars->Draw();
    hParsFitA->Draw("same");
    GetXaxis()->SetTitle("");
    GetYaxis()->SetTitle("Value");

    TLegend *leg = newLegend(kPos9);
    leg->AddEntry(hPars, "Our fit");
    leg->AddEntry(hParsFitA, "H1 2006 FitA");
    DrawLegends({leg});

    can->SaveAs(outDir + "/pars.pdf");

}

void dPlotter::plotCorrelations()
{
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.24, 0.16);
    SetTopBottom(0.2, 0.2);

    map<TString,TString> nMap = { {"PDFQ0_diff.g0", "A_{g}"},
                             {"PDFQ0_diff.g1", "B_{g}"},
                             {"PDFQ0_diff.g2", "C_{g}"},
                             {"PDFQ0_diff.s0", "A_{s}"},
                             {"PDFQ0_diff.s1", "B_{s}"},
                             {"PDFQ0_diff.s2", "C_{s}"},
                             {"QcdnumDDISCS.n_IR", "n_{IR}"},
                             {"QcdnumDDISCS.a0_IP", "a_{IP}(0)"} };

    
    for(int i = 1; i <= hCorrs->GetNbinsX(); ++i) {
       TString s = hCorrs->GetXaxis()->GetBinLabel(i);
       if(nMap.count(s)) {
          hCorrs->GetXaxis()->SetBinLabel(i, nMap[s]);
          hCorrs->GetYaxis()->SetBinLabel(i, nMap[s]);
       }
    }



    hCorrs->Draw("colz text");
    hCorrs->GetXaxis()->SetTitle("");
    hCorrs->GetZaxis()->SetRangeUser(-1, 1);
    
    double orSize = GetXaxis()->GetLabelSize();
    GetXaxis()->SetLabelSize(1.5*orSize);
    GetYaxis()->SetLabelSize(1.5*orSize);

    can->SaveAs(outDir + "/corrs.pdf");

}
