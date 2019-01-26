R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)
#include "plottingHelper.h"
using namespace PlottingHelper;


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TNtuple.h"

#include <vector>
#include <map>


using namespace PlottingHelper;//pollute the namespace!
using namespace std;

TString rn() {return Form("%d",rand());}


#include <vector>


struct point {
    double xp,  q2,  beta, xpSig;
    double th, thErr;
    double thOrg;
    double errStat, errSys, errTot, errUnc;
    std::vector<double> errs;//10 items
};


struct dPlotter {
    std::vector<point> data;
    TString outDir;

    void readData(TString inFile);
    void plotBeta(double xpom);
    void plotQ2(double xpom);
    void plotXpom();
    void plotPDFs(bool inLog);
};



void dplotter(TString inFile = "test/alpos.out.root")
{
    dPlotter dplt;
    dplt.readData(inFile);
    dplt.plotBeta(0.0003);
    dplt.plotBeta(0.001);
    dplt.plotBeta(0.003);
    dplt.plotBeta(0.01);
    dplt.plotBeta(0.03);


    dplt.plotQ2(0.0003);
    dplt.plotQ2(0.001);
    dplt.plotQ2(0.003);
    dplt.plotQ2(0.01);
    dplt.plotQ2(0.03);



    //dplt.plotBeta(0.03f);

}

void dPlotter::readData(TString inFile)
{
    TFile *file = TFile::Open(inFile);
    TNtuple *tuple = (TNtuple*) file->Get("ASaveDataTheory/ThDataTab");


    point pt;
    tuple->SetBranchAddress("xp",&pt.xp);
    tuple->SetBranchAddress("Q2",&pt.q2);
    tuple->SetBranchAddress("beta",&pt.beta);
    tuple->SetBranchAddress("xpSigData",&pt.xpSig);
    tuple->SetBranchAddress("xpSigDataErr",&pt.errTot);
    tuple->SetBranchAddress("xpSigTh",&pt.th);
    tuple->SetBranchAddress("xpSigThErr",&pt.thErr);

    Int_t nentries = (Int_t)tuple->GetEntries();
    for (Int_t i=0;i<nentries;i++) {
        tuple->GetEntry(i);
        //cout << "Hel " << pt.xp <<" "<< pt.xpSig  <<" "<< pt.errTot <<   endl;
        data.push_back(pt);
    }

    outDir =  inFile(0, inFile.Last('/'));
    outDir += "/dPlots";

    gSystem->mkdir(outDir, true);

}


void dPlotter::plotBeta(double xpom)
{
    map<double, vector<point>> dataMap;

    for(point &p : data)
        if(p.xp == xpom) {
            dataMap[p.q2].push_back(p);
        }

    int nRows = ceil(dataMap.size() / 4.);
    //cout << dataMap.size() << endl;

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", (1+4)*200, (1+nRows)*200);
    SetLeftRight(0.5/(1+4), 0.5/(1+4));
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
            gData->SetPoint(j, points[j].beta, points[j].xpSig);

            double er = points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
            //cout <<"ahoj  " << j<<" "<< points[j].beta <<" "<< points[j].xpSig << " " << endl;
        }

        map<double,double> zMin = { {0.0003, 0.02},
                                    {0.001, 0.015},
                                    {0.003, 0.0015},
                                    {0.01, 0.0015},
                                    {0.03, 0.0015}
                                   };

        TH1D *hFr = new TH1D(rn(), "", 1, zMin[xpom], 1);
        hFr->Draw("axis");
        

        TGraph *gThOrg = new TGraph(points.size());
        TGraph *gTh    = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gThOrg->SetPoint(j, points[j].beta, points[j].thOrg);
            gTh->SetPoint(j, points[j].beta, points[j].th);
        }
        gThOrg->SetLineColor(kBlue);
        gThOrg->SetMarkerColor(kBlue);
        gThOrg->Draw("l* same");

        //MyTheory
        gTh->SetLineColor(kGreen);
        gTh->SetMarkerColor(kGreen);
        gTh->Draw("l* same");


        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");



        gPad->SetLogx();

        GetYaxis()->SetNdivisions(503);

        GetYaxis()->SetRangeUser(0, 0.085);

        SetFTO({24}, {14}, {1.4, 2.2, 0.4, 3.9});

        DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));


        if(i < dataMap.size() - 4) {
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

        ++i;
    }

    DrawLatexUp(can->GetPad(1), 1.1, Form("x_{IP} = %g", xpom), -1, "l");

    can->SaveAs(Form(outDir + "/xpom%g.pdf", xpom));

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


void dPlotter::plotQ2(double xpom)
{
    map<double, vector<point>> dataMap;

    for(point &p : data)
        if(p.xp == xpom)
            dataMap[p.beta].push_back(p);

    int nPoints = dataMap.size();

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 500, 800);
    SetLeftRight(0.15, 0.07);
    SetTopBottom(0.1, 0.1);

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
    TH1D *hFr = new TH1D(rn(), "", 1, 2, 1e4);
    hFr->Draw("axis");

    SetFTO({16}, {14}, {1.4, 2.2, 0.4, 3.9});
    GetYaxis()->SetRangeUser(1e-2, 1e5);
    GetYaxis()->SetTitle("3^{i} * x_{IP} #sigma_{r}^{D(3)}");
    GetXaxis()->SetTitle("Q^{2} [GeV]");



    int i = 0;
    for(auto item : dataMap) {
        double beta  = item.first;
        auto &points = item.second;

        double fac = pow(3, dataMap.size() - 1 - i);
        TGraphErrors *gData = new TGraphErrors(points.size());

        for(int j = 0; j < points.size(); ++j) {
            gData->SetPoint(j, points[j].q2, fac*points[j].xpSig);

            double er = fac*points[j].xpSig * points[j].errTot;
            gData->SetPointError(j, 0, er);
        }
        

        TGraph *gTh = new TGraph(points.size());
        TGraph *gThOrg = new TGraph(points.size());
        for(int j = 0; j < points.size(); ++j) {
            gTh->SetPoint(j, points[j].q2, fac*points[j].th);
            gThOrg->SetPoint(j, points[j].q2, fac*points[j].thOrg);
        }

        gThOrg->SetLineColor(kBlue);
        gThOrg->SetMarkerColor(kBlue);
        gThOrg->Draw("l* same");

        gTh->SetLineColor(kGreen);
        gTh->SetMarkerColor(kGreen);
        gTh->Draw("l* same");

        gData->SetMarkerColor(kRed);
        gData->SetLineColor(kRed);
        gData->SetMarkerStyle(20);
        gData->Draw("pe same");



        double xT, yT;
        tie(xT,yT) =  getRight(gData);
        TLatex *lat = new TLatex;
        lat->SetTextSize(PxFontToRel(16));
        lat->SetTextAlign(12);
        lat->DrawLatex(xT*1.3, yT, Form("#beta=%g (i=%lu)", beta, dataMap.size() - 1 - i));


        ++i;
    }

    DrawLatexUp(1, Form("x_{IP} = %g",xpom), -1, "l");


    can->SaveAs(Form(outDir +  "/q2%g.pdf", xpom));

}


void dPlotter::plotXpom()
{
    map<pair<double,double>, vector<point>> dataMap;

    //vector<double> betas = {0.01, 0.04, 0.1, 0.2, 0.4, 0.65, 0.9};
    vector<double> betas = { 0.011, 0.043, 0.11, 0.2, 0.43, 0.67, 0.8};

    vector<double> q2s   = {3.5, 5, 6.5, 8.5, 12, 15, 20, 25, 35, 45, 60, 90, 200, 400, 800, 1600};


    for(point &p : data)
        dataMap[{p.q2, p.beta}].push_back(p);


    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 800);
    SetLeftRight(0.1, 0.1);
    SetTopBottom(0.1, 0.1);

    DivideTransparent(group(1, 0, betas.size()), group(1, 0, q2s.size()));

    for(int iq = 0; iq < q2s.size();   ++iq)
    for(int ib = 0; ib < betas.size(); ++ib) {
        double q2   = q2s[iq];
        double beta = betas[ib];

        can->cd(iq*betas.size() + ib + 1);

        TH1D *hFr = new TH1D(rn(), "", 1, 1e-4, 4e-2);
        hFr->Draw("axis");
        GetYaxis()->SetRangeUser(0, 0.085);

        gPad->SetLogx();
        GetYaxis()->SetNdivisions(303);
        SetFTO({14}, {6}, {1.4, 2.2, 0.4, 3.9});


        if(iq != q2s.size() -1)
            GetXaxis()->SetTickSize(0);
        if(ib != 0)
            GetYaxis()->SetTickSize(0);



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

    }

    //DrawLatexUp(-1.3, Form("Q^{2} = %g GeV^{2}", q2));


    can->SaveAs(outDir + "/xpomGrid.pdf");

}
