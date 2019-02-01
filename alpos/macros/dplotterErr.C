R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

R__LOAD_LIBRARY($PlH_DIR/../alposBuild/libaem.so)

extern "C" void qcd_2006_(double *z,double *q2, int *ifit, double *xPq,       double *f2, double *fl, double *c2, double *cl);

#include "plottingHelper.h"


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TNtuple.h"

#include <vector>
#include <map>
#include <cassert>


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
    map<double, vector<TGraph*>> singletQ2;
    map<double, vector<TGraph*>> gluonQ2;
    double chi2;
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

    void plotPDF(double q2, int fl, bool inLog);


    pair<TGraphAsymmErrors*,TGraphAsymmErrors*> getBandModel(double q2);

    void plotParametersShifts(int sh);
    TGraphAsymmErrors *getPamametersBand(int shMax = 999);


};



void dplotterErr(TString inFile = "../farm/testPdfNLO/H1diffQcdnum_templ.str")
{
    dPlotter dplt;

    dplt.readData(inFile, 7);
    dplt.plotPDFs(false);
    dplt.plotPDF(1.8, 0, false);

    return;

    for(int i = 1; i <= 7; ++i) {
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
    if(!hPars || !hCorrs) {
        cout << "Not loaded " << __LINE__ << endl;
        exit(1);
    }

    chi2 = (dynamic_cast<TFitResult*>(file->Get("fitresult")))->MinFcnValue();


    int nShifts = 2*hPars->GetNbinsX() + 1;
    vector<double> q2Vals = {1.8, 8.5, 20, 90, 800};
    for(double q2 : q2Vals) {
        singletQ2[q2].resize(nShifts, nullptr);
        gluonQ2[q2].resize(nShifts, nullptr);

        for(int i = 0; i < nShifts; ++i) {
            cout << Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_gluon", q2, i ) << endl;
            cout << Form("SaveDPDFTGraph/Q2_%g/DPDF_%d/Pom_d", q2, i ) << endl;
            gluonQ2[q2][i]   = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%.1f/DPDF_%d/Pom_gluon", q2, i )));
            singletQ2[q2][i] = dynamic_cast<TGraph*>(file->Get(Form("SaveDPDFTGraph/Q2_%.1f/DPDF_%d/Pom_d", q2, i )));
            if(!gluonQ2[q2][i] || !singletQ2[q2][i]) {
                cout << "Not loaded " << __LINE__  << endl;
                exit(1);
            }

            //Rescale singlet by 6:
            for(int k = 0; k < singletQ2[q2][i]->GetN(); ++k) {
                double x, y;
                singletQ2[q2][i]->GetPoint(k, x, y);
                singletQ2[q2][i]->SetPoint(k, x, 6*y);
            }
        }
    }



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
        gr->SetPointEYlow(i, errP);
    }
    return gr;
}

double GetMaximum(TGraph *gr)
{
    double m = -1e40;
    for(int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        m = max(m, y);
    }
    return m;
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


TGraphAsymmErrors *getBand(vector<TGraph*> graphs)
{
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(graphs[0]->GetN());
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


pair<TGraphAsymmErrors*,TGraphAsymmErrors*> dPlotter::getBandModel(double q2)
{
    vector<TGraph*> gluons, singlets;
    for(int i = 0; i < shifts.size(); ++i) {
        singlets.push_back(shifts[i].singletQ2.at(q2)[0]);
        gluons.push_back(shifts[i].singletQ2.at(q2)[0]);
    }

    return {getBand(gluons), getBand(singlets)};
}




void dPlotter::plotParameters(int sh)
{
    gStyle->SetOptStat(0);
    double varLeft  = shifts[2*sh-1].chi2 - shifts[0].chi2;
    double varRight = shifts[2*sh].chi2 - shifts[0].chi2;
    double vFit = (varLeft - varRight) / (varLeft + varRight);
    cout << "Shift " <<sh <<" : "<<  varLeft <<" "<< varRight << " " << -vFit << endl;

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

    TGraphAsymmErrors *grAll = getPamametersBand(7);
    grAll->SetFillColorAlpha(kBlue, 0.3);
    grAll->SetLineColor(kBlue);
    grAll->SetFillStyle(1001);
    grAll->Draw("same e2");

    TGraphAsymmErrors *gr = getPamametersBand(6);
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



void dPlotter::plotPDFs(bool inLog)
{
    vector<double> q2s;
    for(auto v : shifts[0].gluonQ2) {
        q2s.push_back(v.first);
    }

    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.1, 0.14);
    SetTopBottom(0.1, 0.1);

    double zMin = 4e-3;
    DivideTransparent(group(1, 0.5, 2), group(1, 0, q2s.size()));

    for(int i = 0; i < q2s.size(); ++i) {
        //Fill Graph

        TGraphAsymmErrors *grS = getBand(shifts[0].singletQ2.at(q2s[i]));
        TGraphAsymmErrors *grG = getBand(shifts[0].gluonQ2.at(q2s[i]));

        TGraphAsymmErrors *grSm, *grGm;
        tie(grGm, grSm) = getBandModel(q2s[i]);

        TGraphAsymmErrors *grStot = addBands({grS, grSm});
        TGraphAsymmErrors *grGtot = addBands({grG, grGm});

        grS->SetLineColor(kBlue);
        grG->SetLineColor(kBlue);

        TGraphAsymmErrors *grSfA = new TGraphAsymmErrors(grS->GetN());
        TGraphAsymmErrors *grGfA = new TGraphAsymmErrors(grS->GetN());

        int ifit = 1;
        for(int j = 0; j < grS->GetN(); ++j) {
            double xPq[13], f2[2], fl[2], c2[2], cl[2];
            double z, vTemp;
            grS->GetPoint(j, z, vTemp);
            qcd_2006_(&z,&q2s[i], &ifit, xPq, f2, fl, c2, cl);
            grSfA->SetPoint(j, z, 6*xPq[7]);
            grGfA->SetPoint(j, z, xPq[6]);
        }


        can->cd(2*i + 1);
        TH1D *hFrS = new TH1D(rn(), "", 1, zMin, 1);
        hFrS->Draw("axis");


        grStot->SetFillColorAlpha(kRed, 0.5);
        grStot->SetFillStyle(1001);
        grStot->Draw("le3 same");


        grS->SetFillColorAlpha(kBlue, 0.5);
        grS->SetFillStyle(1001);
        grS->Draw("le3 same");

        grSfA->Draw("l same");


        GetYaxis()->SetRangeUser(0, 0.27);
        //GetYaxis()->SetRangeUser(0.9, 1.1);
        GetYaxis()->SetNdivisions(503);
        GetXaxis()->SetNdivisions(404);
        SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});

        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Singlet");
            GetYaxis()->SetTitle("z #Sigma(z,Q^{2})");
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

        grGtot->SetFillColorAlpha(kRed, 0.5);
        grGtot->SetFillStyle(1001);
        grGtot->Draw("le3 same");

        grG->SetFillColorAlpha(kBlue, 0.5);
        grG->SetFillStyle(1001);
        grG->Draw("le3 same");

        grGfA->Draw("l same");

        GetYaxis()->SetRangeUser(0, 2.25);
        //GetYaxis()->SetRangeUser(0.9, 1.1);
        GetYaxis()->SetNdivisions(503);
        GetXaxis()->SetNdivisions(404);
        SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});
        if(inLog) gPad->SetLogx();

        if(i == 0) {
            DrawLatexUp(-1, "Gluon");
            GetYaxis()->SetTitle("z g(z,Q^{2})");
        }
        if(i == q2s.size()-1) {
            GetXaxis()->SetTitle("z");
        }
        else {
            GetXaxis()->SetLabelOffset(50000);
        }

        DrawLatexRight(1, Form("Q^{2}=%g",q2s[i]), -1, "l");

        if(i == 3) {
            auto *leg = newLegend(kPos9);
            leg->AddEntry(grG,   "OurFit");
            DrawLegends({leg});
        }
    }

    if(inLog)
        can->SaveAs(outDir + "/pdfsLog.pdf");
    else
        can->SaveAs(outDir + "/pdfsLin.pdf");

}



void dPlotter::plotPDF(double q2, int flav, bool inLog)
{
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas(rn(),"", 600, 600);
    SetLeftRight(0.1, 0.14);
    SetTopBottom(0.1, 0.1);

    double zMin = 4e-3;
    DivideTransparent({1.}, {1,0,0.7});

    //Fill Graph
    TGraphAsymmErrors *gr, *grTot;
    if(flav == 0)
        gr = getBand(shifts[0].gluonQ2.at(q2));
    else
        gr = getBand(shifts[0].singletQ2.at(q2));

    TGraphAsymmErrors *grSm, *grGm;
    tie(grGm, grSm) = getBandModel(q2);

    if(flav == 0)
        grTot = addBands({gr, grGm});
    else
        grTot = addBands({gr, grSm});


    gr->SetLineColor(kBlue);

    TGraphAsymmErrors *grfA = new TGraphAsymmErrors(gr->GetN());

    int ifit = 1;
    for(int j = 0; j < gr->GetN(); ++j) {
        double xPq[13], f2[2], fl[2], c2[2], cl[2];
        double z, vTemp;
        gr->GetPoint(j, z, vTemp);
        qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
        if(flav == 0)
            grfA->SetPoint(j, z, xPq[6]);
        else
            grfA->SetPoint(j, z, 6*xPq[7]);
    }

    can->cd(1);

    TH1D *hFr = new TH1D(rn(), "", 1, zMin, 1);
    hFr->Draw("axis");


    grTot->SetFillColorAlpha(kRed, 0.5);
    grTot->SetFillStyle(1001);
    grTot->Draw("le3 same");


    gr->SetFillColorAlpha(kBlue, 0.5);
    gr->SetFillStyle(1001);
    gr->Draw("le3 same");

    grfA->Draw("l same");


    double Max = max(GetMaximum(grfA), GetMaximum(gr));
    //GetYaxis()->SetRangeUser(0.9, 1.1);
    GetYaxis()->SetNdivisions(503);
    GetXaxis()->SetNdivisions(404);
    SetFTO({15}, {6}, {1.4, 2.2, 0.4, 3.9});

    if(inLog) gPad->SetLogx();

        //DrawLatexUp(-1, "Singlet");
    GetYaxis()->SetTitle("z #Sigma(z,Q^{2})");

    auto *leg = newLegend(kPos9);
    leg->AddEntry(gr,   "OurFit");
    DrawLegends({leg}, true);

    GetYaxis()->SetRangeUser(0, 1.2*Max);

    can->cd(2);

    TH1D *hFrRat = new TH1D(rn(), "", 1, zMin, 1);
    hFrRat->Draw("axis");

    TGraphAsymmErrors *grRat = GetFraction(gr, gr);
    TGraphAsymmErrors *grTotRat = GetFraction(grTot, gr);
    TGraphAsymmErrors *grfARat = GetFraction(grfA, gr);

    grTotRat->SetFillColorAlpha(kRed, 0.5);
    grTotRat->SetFillStyle(1001);
    grTotRat->Draw("le3 same");


    grRat->SetFillColorAlpha(kBlue, 0.5);
    grRat->SetFillStyle(1001);
    grRat->Draw("le3 same");

    grfARat->Draw("l same");


    GetYaxis()->SetRangeUser(0.5, 1.5);
    GetXaxis()->SetTitle("z");

    if(inLog)
        can->SaveAs(outDir + "/pdfLog.pdf");
    else
        can->SaveAs(outDir + "/pdfLin.pdf");

}
