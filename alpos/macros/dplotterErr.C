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
    double thOrgA, thOrgB;
    double errStat, errSys, errTot, errUnc;
    std::vector<double> errs;//10 items
};


struct sysShift {
    std::vector<point> data;
    TH1D *hPars;
    TH2D *hCorrs;
    void readData(TString inFile);
};



struct dPlotter {
    TString outDir = ".";
    vector<sysShift> shifts;
    void readData(TString inFile, int nErr);

    void plotBeta(double xpom);
    void plotQ2(double xpom);
    void plotXpom();
    void plotPDFs(bool inLog);
    void plotParameters(int sh);
    void plotCorrelations();


    void plotParametersShifts(int sh);
    TGraphAsymmErrors *getPamametersBand(int shMax = 999);


};



void dplotterErr(TString inFile = "../farm/testNewNLO/H1diff_templ.str")
{
    dPlotter dplt;
    dplt.readData(inFile, 9);

    for(int i = 1; i <= 9; ++i) {
        dplt.plotParameters(i);
        dplt.plotParametersShifts(i);
    }
    //dplt.plotParameters(2);
    //dplt.plotCorrelations();
    //dplt.plotBeta(0.03f);

}


void sysShift::readData(TString inFile)
{
    TFile *file = TFile::Open(inFile);
    cout << inFile << endl;

    hPars  = dynamic_cast<TH1D*>(file->Get("fitparameters"));
    hCorrs = dynamic_cast<TH2D*>(file->Get("fitcorrelations"));
    assert(hPars);
    assert(hCorrs);

    return;

    TNtuple *tuple = (TNtuple*) file->Get("ASaveDataTheory/ThDataTab");

    point pt;
    tuple->SetBranchAddress("xp",&pt.xp);
    tuple->SetBranchAddress("Q2",&pt.q2);
    tuple->SetBranchAddress("beta",&pt.beta);
    tuple->SetBranchAddress("xpSigData",&pt.xpSig);
    tuple->SetBranchAddress("xpSigDataErr",&pt.errTot);
    tuple->SetBranchAddress("xpSigTh",&pt.th);
    tuple->SetBranchAddress("xpSigThErr",&pt.thErr);

    tuple->SetBranchAddress("xpSigThOrgA",&pt.thOrgA);
    tuple->SetBranchAddress("xpSigThOrgB",&pt.thOrgB);


    Int_t nentries = (Int_t)tuple->GetEntries();
    for (Int_t i=0;i<nentries;i++) {
        tuple->GetEntry(i);
        //cout << "Hel " << pt.xp <<" "<< pt.xpSig  <<" "<< pt.errTot <<   endl;
        data.push_back(pt);
    }

}

void dPlotter::readData(TString inFile, int nErr)
{
    shifts.resize(2*nErr+1);
    shifts[0].readData(inFile+"0_dir/out.root");
    for(int i = 0; i < nErr; ++i) {
        shifts[2*i+1].readData(inFile+(i+1)+"u_dir/out.root");
        shifts[2*i+2].readData(inFile+(i+1)+"d_dir/out.root");
    }
}


TGraphAsymmErrors *dPlotter::getPamametersBand(int shMax)
{
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(shifts[0].hPars->GetNbinsX());
    for(int i = 0; i < shifts[0].hPars->GetNbinsX(); ++i) {

        double vCnt = shifts[0].hPars->GetBinContent(i+1);

        double errP = 0, errM = 0;
        for(int j = 1; j < shifts.size() && j <= 2*shMax; ++j) {
            double v = shifts[j].hPars->GetBinContent(i+1) - vCnt;
            errP = hypot(errP, max(0.0, v));
            errM = hypot(errM, max(0.0,-v));
        }

        double x = shifts[0].hPars->GetBinCenter(i+1);
        double w = shifts[0].hPars->GetBinWidth(i+1);
        gr->SetPoint(i, x, vCnt);
        gr->SetPointError(i, w/2, w/2, errM, errP);

    }
    return gr;

}


void dPlotter::plotParameters(int sh)
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

    TH1D *hParsFitA = (TH1D*) shifts[0].hPars->Clone("FitA"+rn());
    hParsFitA->Reset();
    for(int i = 1; i <= hParsFitA->GetNbinsX(); ++i) {
        TString s = hParsFitA->GetXaxis()->GetBinLabel(i);
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

    //Remove Errors
    for(int i = 1; i <shifts.size(); ++i) {
        if(!shifts[i].hPars)
            shifts[i].hPars =  (TH1D*) shifts[0].hPars->Clone();

        for(int j = 1; j <= shifts[0].hPars->GetNbinsX(); ++j)
            shifts[i].hPars->SetBinError(j, 1e-8);
    }


    shifts[0].hPars->Draw();
    shifts[0].hPars->SetLineColor(kBlue);
    shifts[0].hPars->SetLineWidth(2);

    shifts[2*(sh-1)+1].hPars->Draw("same");
    shifts[2*(sh-1)+2].hPars->Draw("same");
    shifts[2*(sh-1)+1].hPars->SetLineWidth(2);
    shifts[2*(sh-1)+2].hPars->SetLineWidth(2);
    shifts[2*(sh-1)+1].hPars->SetLineStyle(2);
    shifts[2*(sh-1)+2].hPars->SetLineStyle(2);

    TGraphAsymmErrors *grAll = getPamametersBand(9);
    grAll->SetFillColorAlpha(kBlue, 0.3);
    grAll->SetLineColor(kBlue);
    grAll->SetFillStyle(1001);
    grAll->Draw("same e2");

    TGraphAsymmErrors *gr = getPamametersBand(8);
    gr->SetFillColorAlpha(kBlue, 0.4);
    gr->SetLineColor(kBlue);
    gr->SetFillStyle(1001);
    gr->Draw("same e2");



    hParsFitA->Draw("same");
    hParsFitA->SetLineColor(kRed);
    hParsFitA->SetLineStyle(2);
    hParsFitA->SetLineWidth(2);

    GetXaxis()->SetTitle("");
    GetYaxis()->SetTitle("Value");
    GetYaxis()->SetRangeUser(-2.1, 3.1);

    TLegend *leg = newLegend(kPos9);
    leg->AddEntry(shifts[0].hPars, "Our fit");
    leg->AddEntry(hParsFitA, "H1 2006 FitA");
    DrawLegends({leg});

    can->SaveAs(outDir + Form("/pars%d.pdf",sh));
}


void dPlotter::plotParametersShifts(int sh)
{
    gStyle->SetOptStat(0);

    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.13, 0.16);

    TH1D *h   = (TH1D*) shifts[0].hPars->Clone();
    TH1D *hUp = (TH1D*) shifts[2*(sh-1)+1].hPars->Clone();
    TH1D *hDn = (TH1D*) shifts[2*(sh-1)+2].hPars->Clone();
    
    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double v = h->GetBinContent(i);

        double vUp = (hUp->GetBinContent(i) - h->GetBinContent(i)) / h->GetBinError(i);
        double vDn = (hDn->GetBinContent(i) - h->GetBinContent(i)) / h->GetBinError(i);

        hUp->SetBinContent(i, vUp);
        hDn->SetBinContent(i, vDn);
        
        hUp->SetBinError(i, 0);
        hDn->SetBinError(i, 0);

    }


    //shifts[0].hPars->Draw();


    hUp->Draw();
    hUp->SetLineColor(kRed);
    hDn->Draw("same");
    hDn->SetLineStyle(2);


    GetXaxis()->SetTitle("");
    GetYaxis()->SetTitle("Value");
    GetYaxis()->SetRangeUser(-6, 6);

    /*
    TLegend *leg = newLegend(kPos9);
    leg->AddEntry(shifts[0].hPars, "Our fit");
    leg->AddEntry(hParsFitA, "H1 2006 FitA");
    DrawLegends({leg});
    */

    can->SaveAs(outDir + Form("/parsErr%d.pdf",sh));
}








void dPlotter::plotCorrelations()
{
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.2f");
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.24, 0.16);
    SetTopBottom(0.2, 0.2);

    //hCorrs->Draw("colz text");
    //hCorrs->GetXaxis()->SetTitle("");
    //hCorrs->GetZaxis()->SetRangeUser(-1, 1);
    

    can->SaveAs(outDir + "/corrs.pdf");

}
