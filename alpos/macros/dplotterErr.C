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


};



void dplotterErr(TString inFile = "../farm/testNewNLOvfns/H1diff_templ.str")
{
    dPlotter dplt;
    dplt.readData(inFile, 10);

    for(int i = 1; i <= 10; ++i) {
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


void dPlotter::plotParameters(int sh)
{
    gStyle->SetOptStat(0);

    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.13, 0.16);

    map<TString, vector<double>> pars;
    //Parameters  for fitA
    pars["g0"]  =  {  0.14591    ,   0.33171E-01   };
    pars["g2"]  =   { -0.94705    ,   0.20309       };
    pars["q0"] = {  1.0587     ,   0.322116   };
    pars["q1"] = {   2.2964    ,   0.36439       };
    pars["q2"] = {  0.56894    ,   0.14969       };
    pars["n_IR"] =     {  0.16966E-02,   0.41732E-03   };
    pars["a0_IP"] =   {   1.1182    ,   0.81319E-02 };

    TH1D *hParsFitA = (TH1D*) shifts[0].hPars->Clone("FitA");
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
    hParsFitA->SetLineColor(kRed);
    hParsFitA->SetLineStyle(2);

    shifts[0].hPars->Draw();

    shifts[2*(sh-1)+1].hPars->Draw("same");
    shifts[2*(sh-1)+2].hPars->Draw("same");

    hParsFitA->Draw("same");

    GetXaxis()->SetTitle("");
    GetYaxis()->SetTitle("Value");

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
