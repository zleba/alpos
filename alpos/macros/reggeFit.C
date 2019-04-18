R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

#include "TString.h"
#include <fstream>
#include <iostream>
#include <set>
#include "plottingHelper.h"


using namespace std;
using namespace PlottingHelper;

TString rn() {return Form("%d", rand()); }


struct point {
    point(double _xp, double _q2, double _beta) : xp(_xp), q2(_q2), beta(_beta) {}
    point()  {}
    double xp, q2, beta, tAbs;
    double sigma, err;
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
            point p;

            double mx, stat, unc;
            if(!fName.Contains("FPS")) {
                iss >> p.xp >> p.q2 >> p.beta >> mx >> p.sigma >> stat >> unc;
                p.tAbs = -1;
                p.err = hypot(stat,unc) * p.sigma/100.;
            }
            else if(fName.Contains("FPS-4D")) {

                double mx, mBin, exc, relUnc, relUncMC;
                iss >> p.q2 >> p.beta >> p.xp >> p.tAbs>> mx >> mBin >> p.sigma >> exc >> relUnc >> relUncMC;
                p.err = p.sigma*hypot(relUnc,0.7*relUncMC)/100.;

            }
            else {
                //iss >> p.xp >> p.q2 >> p.beta;

                //iss >> p.xp >> p.q2 >> p.beta >> mx >> p.sigma >> stat >> unc;
                double mx, mBin, exc, relUnc, relUncMC;
                iss >> p.q2 >> p.beta >> p.xp >> mx >> mBin >> p.sigma >> exc >> relUnc >> relUncMC;
                p.err = p.sigma*hypot(relUnc,0.7*relUncMC)/100.;
                p.tAbs = -1;
            }
            
            //if(p.xp < 0.04)
                points.push_back(p);


        }

        ++iIn;

    }

    return points;

}

double chi2Raw(vector<double> &data, vector<double> &err,  vector<double> &aVec, vector<double> &bVec)
{
    assert(data.size() == aVec.size());
    assert(data.size() == err.size());

    //At least 3 points
    if(data.size() < 3) return 0;

    double Ya = 0, Yb = 0;
    double Aa = 0, Ba = 0;
    double Ab = 0, Bb = 0;

    //cout << "DataSize " << data.size() << endl;
    for(int i = 0; i < data.size(); ++i) {
        double err4 = pow(err[i],4);
        Ya += data[i]*aVec[i]/err4;
        Yb += data[i]*bVec[i]/err4;

        Aa += aVec[i]*aVec[i]/err4;
        Ba += bVec[i]*aVec[i]/err4;

        Ab += aVec[i]*bVec[i]/err4;
        Bb += bVec[i]*bVec[i]/err4;
    }

    //cout << "Bla_a " << Aa <<" "<< Ba << endl;
    //cout << "Bla_b " << Ab <<" "<< Bb << endl;

    double D  = Aa*Bb - Ba*Ab;
    double Da = Ya*Bb - Ba*Yb;
    double Db = Aa*Yb - Ab*Ya;

    //cout << "Discriminant " << D << endl;
    double A = Da / D;
    double B = Db / D;

    double chi2 = 0;
    for(int i = 0; i < data.size(); ++i) {
        double err2 = pow(err[i],2);
        //cout << i <<" "<< data[i] << " "<< aVec[i] << " "<< bVec[i] << endl;
        chi2 += pow(data[i] - A*aVec[i] - B*bVec[i], 2) / err2;
    }

    //cout << chi2 << endl;
    return chi2;

}

double chi2Raw(vector<double> &data, vector<double> &err,  vector<double> &aVec)
{
    assert(data.size() == aVec.size());
    assert(data.size() == err.size());

    //At least 2 points
    if(data.size() < 2) return 0;

    double Ya = 0;
    double Aa = 0;

    //cout << "DataSize " << data.size() << endl;
    for(int i = 0; i < data.size(); ++i) {
        double err4 = pow(err[i],4);
        Ya += data[i]*aVec[i]/err4;
        Aa += aVec[i]*aVec[i]/err4;
    }

    double A = Ya / Aa;

    double chi2 = 0;
    for(int i = 0; i < data.size(); ++i) {
        double err2 = pow(err[i],2);
        //cout << i <<" "<< data[i] << " "<< aVec[i] << " "<< bVec[i] << endl;
        chi2 += pow(data[i] - A*aVec[i], 2) / err2;
    }

    //cout << chi2 << endl;
    return chi2;

}





double fluxI(double xp, double a0, double aRat)
{
    return xp * pow(1./xp, 2*a0-1) * 1. / (1. - 2*aRat*log(xp));
}


double flux(double xp, double tAbs, double a0, double ap, double b0)
{
    return xp * pow(1./xp, 2*(a0-ap*tAbs)- 1) * exp(-b0*tAbs);
}



double getChi2(const vector<point> &points,  vector<double> pars)
{
    map<pair<double,double>, vector<point>> pMap;

    for(auto &p : points) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    double chi2 = 0;
    int nDf = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        
        vector<double> data, err, aVec, bVec;
        for(auto &p : ps) {

            double aFlux=0, bFlux=0;

            if(pars.size() == 4) {
                aFlux = fluxI(p.xp, pars[0], pars[1]);
                bFlux = fluxI(p.xp, pars[2], pars[3]);
            }
            else if(pars.size() == 6) {
                aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
                bFlux = flux(p.xp,p.tAbs, pars[3], pars[4], pars[5]);
            }

            else if(pars.size() == 2) {
                aFlux = fluxI(p.xp, pars[0], pars[1]);
            }
            else if(pars.size() == 3) {
                aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
            }



            //cout << "bla " << p.xp << " "<< aFlux <<" "<< bFlux << " : "<< p.sigma << " "<< p.err << endl;

            data.push_back(p.sigma);
            err.push_back(p.err);
            aVec.push_back(aFlux);
            bVec.push_back(bFlux);
        }
        if(pars.size() == 4 || pars.size() == 6)
            nDf += max(int(data.size() - 2), 0);
        else
            nDf += max(int(data.size() - 1), 0);

        //cout << "Chi2Now " << chi2Raw(data, err,  aVec, bVec) << endl;
        if(pars.size() == 4 || pars.size() == 6)
            chi2 += chi2Raw(data, err,  aVec, bVec);
        else
            chi2 += chi2Raw(data, err,  aVec);
    }
    cout <<"chi2="<< chi2 << " " << nDf << endl;
    //exit(0);
    return chi2;
}




void fit(vector<point> &pointsComb)
{




}

struct GlobalChi2 {
    GlobalChi2()  {}
    GlobalChi2(int nFl, vector<point> points_) : nFluxes(nFl), points(points_)  {}
    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)

    int nFluxes;
    vector<point> points;

    double operator() (const double *par) const {
        //cout <<"Hela " <<  par[0] << " "<< par[1] << " "<< par[2] << endl;
        //return hypot(par[0],par[1]) + hypot(par[1],par[2]);

        if(nFluxes == 2) {
            if(points[0].tAbs < -0.5) return getChi2(points, {par[0], par[1], par[2], par[3]});
            else      return getChi2(points, {par[0], par[1], par[2], par[3], par[4], par[5]});
        }
        else  { //Single Flux
            if(points[0].tAbs < -0.5) return getChi2(points, {par[0], par[1]});
            else      return getChi2(points, {par[0], par[1], par[2]});
        }
    }
};


void plotBslope(vector<point> &points)
{
    gStyle->SetOptStat(0);

    map<pair<double,double>, vector<point>> pMap;

    for(auto &p : points) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    map<pair<double,double>,TGraphErrors*> grSlopes;

    double chi2 = 0;
    int nDf = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        double q2, beta;
        tie(beta,q2) = pm.first;
        
        map<double, TGraphErrors*> sigmaT;

        for(auto &p : ps) {
            double aFlux, bFlux;

            if(sigmaT.count(p.xp) == 0)
                sigmaT[p.xp] = new TGraphErrors();

            sigmaT.at(p.xp)->SetPoint(sigmaT.at(p.xp)->GetN(), p.tAbs, p.sigma);
            sigmaT.at(p.xp)->SetPointError(sigmaT.at(p.xp)->GetN()-1, 0, p.err);
        }

        for(auto gr : sigmaT) {
            TFitResultPtr res = gr.second->Fit("expo", "ES");
            double slope = -res->Parameter(1); 
            double err   = res->ParError(1);

            if(grSlopes.count(pm.first) == 0)
                grSlopes[pm.first] = new TGraphErrors();

            grSlopes.at(pm.first)->SetPoint(grSlopes.at(pm.first)->GetN(), gr.first, slope);
            grSlopes.at(pm.first)->SetPointError(grSlopes.at(pm.first)->GetN()-1, 0, err);
        }
    }


    auto can = new TCanvas("can","", 600, 600);




    set<double> q2Set, betaSet;
    for(auto s : grSlopes) {
        betaSet.insert(s.first.first);
        q2Set.insert(s.first.second);
    }

    cout << "Sizes " << q2Set.size() <<" "<< betaSet.size() << endl;

    DividePad(vector<double>(q2Set.size(),1), vector<double>(betaSet.size(),1));

    int iq2 = 0, ibeta = 0;
    for(auto q2 : q2Set) {
        ibeta = 0;
        for(auto beta : betaSet) {
            if(grSlopes.count({beta,q2}) >= 1) {
                can->cd(iq2*betaSet.size() + ibeta + 1);
                auto gr = grSlopes.at({beta,q2});

                cout << "RADEK " << ibeta <<" "<< iq2 <<" "<< endl;

                TH1D *h = new TH1D(rn(), "", 1, 1e-3, 1e-1);
                h->Draw("axis");
                gr->Draw("*e same");
                h->SetMinimum(0);
                h->SetMaximum(8);
                gPad->SetLogx();

            }
            else {
                cout << "Nothing for " << beta <<" "<< q2 << endl;
            }
            ++ibeta;
        }
        ++iq2;
    }


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





void reggeFit()
{

    //auto h1Comb = readInclusive("../datafiles/h1/1203.4495/h1_ddis_Comb.dat");
    //auto h1_225 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_225.dat");
    //auto h1_252 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_252.dat");

    //auto h1fps  = readInclusive("../datafiles/h1/1010.1476/H1-DDIS-HERAII-FPS.dat");
    auto h1fps  = readInclusive("../datafiles/h1/1010.1476/H1-DDIS-HERAII-FPS-4D.dat");

    auto h1Data = h1fps;

    //plotBslope(h1fps);
    //return;


    /*
    sort(h1fps.begin(), h1fps.end(), [](point p1, point p2) {
        return (p1.q2*100 + p1.beta) > (p2.q2*100 + p2.beta);
            });


    for(auto p : h1fps) {
        cout << p.q2 << " "<< p.beta <<" "<< p.xp << endl;
    }
    */


    //return;


    int nFluxes = 2;
    GlobalChi2 globalChi2(nFluxes, h1Data);
    ROOT::Fit::Fitter fitter;

    int Npar = (h1Data[0].tAbs < -0.5) ? 2 : 3;
    Npar *= nFluxes;
    

    vector<double> pars;
    if( Npar == 4) {
        pars = { 1.1,0.02, 
                 0.5, 0.3};
    }
    else if( Npar == 6) {
        pars = { 1.1,0.25, 6,
                 0.5, 0.5, 1.5};
    }
    else if(Npar == 3) {
        pars = { 1.05,0.04, 6 };
    }
    else if(Npar == 2) {
        pars = { 1.1, 0.02 };
    }


    fitter.Config().SetParamsSettings(Npar,pars.data());
    // create before the parameter settings in order to fix or set range on them
    // fix 5-th parameter
    if( Npar == 4) {

        fitter.Config().ParSettings(0).SetLimits(0.9,5.5);
        fitter.Config().ParSettings(1).SetLimits(-0.51,0.53);

        fitter.Config().ParSettings(2).SetLimits(0.05,0.9);
        fitter.Config().ParSettings(3).SetLimits(-8.3,8.9);
    }
    else if(Npar == 6) {
        fitter.Config().ParSettings(0).SetLimits(0.9,5.5);
        fitter.Config().ParSettings(1).SetLimits(-0.51,0.53);
        fitter.Config().ParSettings(2).SetLimits(2,9);

        fitter.Config().ParSettings(3).SetLimits(0.1,0.95);
        fitter.Config().ParSettings(4).SetLimits(-0.9,0.9);
        fitter.Config().ParSettings(5).SetLimits(0.1,4.6);

        //fitter.Config().ParSettings(1).Fix();
        fitter.Config().ParSettings(3).Fix();

    }
    else if(Npar == 3) {
        fitter.Config().ParSettings(0).SetLimits(0.2,5.5);
        fitter.Config().ParSettings(1).SetLimits(-0.71,0.73);
        fitter.Config().ParSettings(2).SetLimits(1,10);

        //fitter.Config().ParSettings(0).Fix();
    }

    //return;


    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(Npar,globalChi2);//,0,dataB.Size()+dataSB.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
    result.PrintCovMatrix(cout);

    for(int i = 0; i < Npar; ++i)
        cout << i <<" "<< result.GlobalCC(i) << endl;




    //plotBetaQ2(h1Comb, h1_225, h1_252, h1fps);
    //plotBetaQ2Xpom(h1Comb, h1_225, h1_252, h1fps);


}
