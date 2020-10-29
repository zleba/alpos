R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

//Plot the coverage of the phase space by the data points



#include "TString.h"
#include <fstream>
#include <iostream>
#include "plottingHelper.h"


using namespace std;
using namespace PlottingHelper;

struct point {
    point(double _xp, double _q2, double _beta) : xp(_xp), q2(_q2), beta(_beta) {}
    double xp, q2, beta;
};

vector<point> readInclusive(TString fName)
{
    ifstream file(fName);

    vector<point> points;

    std::string line;
    bool isIn = false;
    int iIn = 0;
    while (std::getline(file, line)) {
        TString rLine = line;
        if(rLine.Contains("Data {{"))
            isIn = true;
        if(rLine.Contains("}}") && isIn)
            isIn = false;

        if(!isIn) continue;
        if(rLine.Contains("#")) continue;


        //Reading the properties
        if(iIn <= 1) {
            cout << "Reading" << endl;

        }
        else {
            std::istringstream iss(line);
            double xp, q2, beta;
            if(!fName.Contains("FPS"))
                iss >> xp >> q2 >> beta;
            else
                iss >> q2 >> beta >> xp;
            
            cout << xp <<" "<< q2 <<" "<< beta << endl;
            points.push_back(point(xp, q2, beta));


        }

        ++iIn;

    }

    return points;

}

void plotBetaQ2(vector<point> &pointsComb, vector<point> &points225, vector<point> &points252, vector<point> &pointsFPS)
{
    gStyle->SetOptStat(0);

    auto can = new TCanvas("can","", 600, 600);

    TGraph *grComb = new TGraph();
    for(auto &p : pointsComb) grComb->SetPoint(grComb->GetN(), p.beta, p.q2);

    TGraph *gr225 = new TGraph();
    for(auto &p : points225) gr225->SetPoint(gr225->GetN(), p.beta, p.q2);

    TGraph *gr252 = new TGraph();
    for(auto &p : points252) gr252->SetPoint(gr252->GetN(), p.beta, p.q2);

    TGraph *grFPS = new TGraph();
    for(auto &p : pointsFPS) grFPS->SetPoint(grFPS->GetN(), p.beta, p.q2);


    TH1D *hFr = new TH1D("hFr", "",  1, 1e-3, 1);
    hFr->Draw("axis");

    grComb->SetMarkerStyle(20);
    grComb->Draw("p same");
    gr225->SetMarkerColor(kRed);
    gr225->SetMarkerStyle(21);
    gr225->Draw("p same");
    gr252->SetMarkerColor(kBlue);
    gr252->SetMarkerStyle(22);
    gr252->Draw("p same");

    grFPS->SetMarkerColor(kYellow);
    grFPS->SetMarkerStyle(23);
    grFPS->Draw("p same");


    gPad->SetLogx();
    gPad->SetLogy();

    GetYaxis()->SetRangeUser(2, 2000);

    GetXaxis()->SetTitle("#beta");
    GetYaxis()->SetTitle("Q^{2}");
}


void plotBetaQ2Xpom(vector<point> &pointsComb, vector<point> &points225, vector<point> &points252, vector<point> &pointsFPS)
{
    gStyle->SetOptStat(0);

    auto can = new TCanvas("can","", 600, 600);

    TGraph2D *grComb = new TGraph2D();
    for(auto &p : pointsComb) grComb->SetPoint(grComb->GetN(), p.beta, p.q2, p.xp);

    TGraph2D *gr225 = new TGraph2D();
    for(auto &p : points225) gr225->SetPoint(gr225->GetN(), p.beta, p.q2, p.xp);

    TGraph2D *gr252 = new TGraph2D();
    for(auto &p : points252) gr252->SetPoint(gr252->GetN(), p.beta, p.q2, p.xp);

    TGraph2D *grFPS = new TGraph2D();
    for(auto &p : pointsFPS) grFPS->SetPoint(grFPS->GetN(), p.beta, p.q2, p.xp);


    //TH1D *hFr = new TH1D("hFr", "",  1, 1e-3, 1);
    //hFr->Draw("axis");

    grComb->SetMarkerStyle(20);
    grComb->Draw("p");

    /*
    gr225->SetMarkerColor(kRed);
    gr225->SetMarkerStyle(21);
    gr225->Draw("p same");
    gr252->SetMarkerColor(kBlue);
    gr252->SetMarkerStyle(22);
    gr252->Draw("p same");
    */

    grFPS->SetMarkerColor(kRed);
    grFPS->SetMarkerStyle(23);
    grFPS->Draw("p same");

    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();

    //GetYaxis()->SetRangeUser(2, 2000);

    grComb->GetXaxis()->SetTitle("#beta");
    grComb->GetYaxis()->SetTitle("Q^{2}");
    grComb->GetZaxis()->SetTitle("x_{IP}");

    grComb->SetMaximum(0.08);
    grComb->SetMinimum(0.001);


}





void plotCoverage()
{

    auto h1Comb = readInclusive("../datafiles/h1/1203.4495/h1_ddis_Comb.dat");
    auto h1_225 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_225.dat");
    auto h1_252 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_252.dat");
    auto h1fps  = readInclusive("../datafiles/h1/1010.1476/H1-DDIS-HERAII-FPS.dat");

    //plotBetaQ2(h1Comb, h1_225, h1_252, h1fps);
    plotBetaQ2Xpom(h1Comb, h1_225, h1_252, h1fps);


}
