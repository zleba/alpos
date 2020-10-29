R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)


R__LOAD_LIBRARY($PROJECT_DIR/alposBuild/libaem.so)

extern "C" void qcd_2006_(double *z,double *q2, int *ifit, double *xPq,   double *f2, double *fl, double *c2, double *cl);


//Regge fit of H1 data
//The F2 component needs to be subtracted from some QCD fit


#include "TString.h"
#include <fstream>
#include <iostream>
#include <set>
#include "plottingHelper.h"


using namespace std;
using namespace PlottingHelper;

TString rn() {return Form("%d", rand()); }
double rflux(double a0, double ap, double b0, double x_pom, double tAbs);






double getRoundStep(double a) //double a = 3.0023;
{
    double b;
    double t = a / 1e2;

    while(1) {
        double x = modf (t, &b);
        if(x > 0.5) x = 1 - x;
        if(x < 1e-8) break;
        t *= 10;
    }

    //cout << t <<" "<<  a/t <<  endl;
    return a/t;
}

double smear(double a)
{
    //return a;
    double s = getRoundStep(a);
    //cout << "Helenka " << getRoundStep(0.0094) << endl;
    //exit(0);
    //cout <<"Check " << setprecision(15) <<  a <<" "<< s << endl;
    return a + gRandom->Uniform(-s/2, s/2);
}

//data point
struct point {
    point(double _xp, double _q2, double _beta) : xp(_xp), q2(_q2), beta(_beta) {}
    point()  {}
    double xp, q2, beta, tAbs;
    double sigma, err;
};

//Read data from text file to vector<point>
vector<point> readInclusive(TString fName)
{
    gRandom->SetSeed(0);
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
                iss >> p.q2 >> p.beta >> p.xp >> p.tAbs>> mx >> mBin >> p.sigma >>  relUnc >> relUncMC;
                double sys, tot;
                iss >> sys >> tot;
                //p.err = p.sigma*hypot(relUnc,0.7*relUncMC)/100.;
                //p.err = p.sigma*hypot(relUnc,0*relUncMC)/100.;

                p.sigma = smear(p.sigma);

                p.err = p.sigma*smear(tot) / 100;
                //p.err = p.sigma*relUnc / 100;
                //cout << p.q2 <<" " << p.beta <<" : "<< relUncMC<<" "<< sys <<" "<< tot << endl;

            }
            else {
                //iss >> p.xp >> p.q2 >> p.beta;

                //iss >> p.xp >> p.q2 >> p.beta >> mx >> p.sigma >> stat >> unc;
                double mx, mBin, exc, relUnc, relUncMC;
                iss >> p.q2 >> p.beta >> p.xp >> mx >> mBin >> p.sigma >> exc >> relUnc >> relUncMC;
                p.err = p.sigma*hypot(relUnc,0.7*relUncMC)/100.;
                p.tAbs = -1;
            }
            
            const double s = 4*920*27.6;
            const double x = p.xp*p.beta;
            const double y = p.q2/(s*x);
            //if(y < 0.3) //to remove Fl sensitive region
                points.push_back(p);


        }

        ++iIn;

    }

    return points;

}


//"Subtract" fL to get F2 from sigmaRed
vector<point> subtractFL(vector<point> points)
{
    int ifit = 2;
    for(auto &p : points) {
        double xPq[13], f2[2], fl[2], c2[2], cl[2];
        double z;
        z  = p.beta;
        double q2 = p.q2;
        qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
        ifit = 0;

        const double s = 319*319;
        double y = q2 / (s*p.xp*p.beta);
        double FL = y*y / (1+pow(1-y,2))*fl[0];
        double F2 = f2[0];


        p.sigma = F2/(F2 - FL) * p.sigma;
    }
    return points;
}


//fit data+err vector with A*aVec + B*bVec, where aVec,bVec is the theory and A,B are free parameters (i.e. single q2,beta point)
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

//fit data+err vector with A*aVec, where aVec is the theory and A is the free parameter (i.e. single q2,beta point)
double chi2Raw(vector<double> &data, vector<double> &err,  vector<double> &aVec, double *norm = nullptr)
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
    if(norm) *norm = A;

    double chi2 = 0;
    for(int i = 0; i < data.size(); ++i) {
        double err2 = pow(err[i],2);
        //cout << i <<" "<< data[i] << " "<< aVec[i] << " "<< bVec[i] << endl;
        chi2 += pow(data[i] - A*aVec[i], 2) / err2;
    }

    //cout << chi2 << endl;
    return chi2;
}




//t-integrated flux
double fluxI(double xp, double a0, double aRat)
{
    return xp * pow(1./xp, 2*a0-1) * 1. / (1. - 2*aRat*log(xp));
}


//flux
double flux(double xp, double tAbs, double a0, double ap, double b0)
{
    return xp * rflux(a0, ap, b0, xp, tAbs);
    //return xp * pow(1./xp, 2*(a0-ap*tAbs)- 1) * exp(-b0*tAbs);
}



double getChi2(TString type, const vector<point> &points,  vector<double> pars)
{
    //group points according to beta,q2
    map<pair<double,double>, vector<point>> pMap;

    for(auto &p : points) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    double chi2 = 0;
    int nDf = 0;
    int nPoints = 0;
    int nGroups = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        
        vector<double> data, err, aVec, bVec;
        for(auto &p : ps) { //loop over points with equal beta,q2

            double aFlux=0, bFlux=0;

            if(type == "tInt_PR") { //fit t-integrated data (IP + IR)
                aFlux = fluxI(p.xp, pars[0], pars[1]);
                bFlux = fluxI(p.xp, pars[2], pars[3]);
            }
            else if(type == "tDep_PR") { //fit of t-dependent data (IP + IR)
                aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
                bFlux = flux(p.xp,p.tAbs, pars[3], pars[4], pars[5]);
            }

            else if(type == "tInt_P") { //fit t-integrated data (IP)
                aFlux = fluxI(p.xp, pars[0], pars[1]);
            }
            else if(type == "tDep_P") {
                aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
            }

            //cout << "bla " << p.xp << " "<< aFlux <<" "<< bFlux << " : "<< p.sigma << " "<< p.err << endl;

            data.push_back(p.sigma);
            err.push_back(p.err);
            aVec.push_back(aFlux);
            bVec.push_back(bFlux);
        }
        if(type.Contains("tDep")) {
            nDf += max(int(data.size() - 2), 0);
            if(data.size() - 2 > 0) {
                ++nGroups;
                nPoints += data.size();
            }
        }
        else
            nDf += max(int(data.size() - 1), 0);

        //cout << "Chi2Now " << chi2Raw(data, err,  aVec, bVec) << endl;
        if(type.Contains("_PR"))
            chi2 += chi2Raw(data, err,  aVec, bVec);
        else
            chi2 += chi2Raw(data, err,  aVec);
    }
    cout <<"chi2= "<< chi2 << " " << nDf << endl;
    cout <<"nGroups=" << nGroups << " nPoints "<< nPoints << endl;
    //exit(0);
    return chi2;
}

double getChi2RegFixed(TString type, const vector<point> &points,  vector<double> pars)
{
    //group points according to beta,q2
    map<pair<double,double>, vector<point>> pMap;

    for(auto &p : points) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    //cout <<"#points " <<  points.size() << endl;
    //exit(0);


    double chi2 = 0;
    int nDf = 0;
    int ifit = 2;
    int nPoints = 0;
    int nGroups = 0;

    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        
        vector<double> data, err, aVec, bVec;
        for(auto &p : ps) { //loop over points with equal beta,q2

            double aFlux=0, bFlux=0, aFluxFitB = 0;

            const double a0_R = 0.5;
            const double ap_R = 0.3;
            const double b0_R = 1.6;

            const double a0_P = 1.11101;
            const double ap_P = 0.06;
            const double b0_P = 5.5;


            //All fluxes multiplied by xIP
            if(type == "tInt_Rfix") { //fit t-integrated data (IP + IR)
                aFlux = fluxI(p.xp, pars[0], pars[1]);
                bFlux = fluxI(p.xp, a0_R, ap_R/b0_R);
            }
            else if(type == "tDep_Rfix") { //fit of t-dependent data (IP + IR)
                aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);
                bFlux = flux(p.xp,p.tAbs, a0_R, ap_R, b0_R);
                aFluxFitB = flux(p.xp,p.tAbs, a0_P, ap_P, b0_P);
            }

            //cout << "bla " << p.xp << " "<< aFlux <<" "<< bFlux << " : "<< p.sigma << " "<< p.err << endl;

            //subtract reggeon and FlPom from data
            double xPq[13], f2[2], fl[2], c2[2], cl[2];
            double z = p.beta;
            double q2 = p.q2;
            qcd_2006_(&z,&q2, &ifit, xPq, f2, fl, c2, cl);
            ifit = 0;


            const double s = 319*319;
            double y = q2 / (s*p.xp*p.beta);
            //double FL = y*y / (1+pow(1-y,2))*fl[0];
            double flC = y*y / (1+pow(1-y,2));

            double bFluxFitB = bFlux;
            double corr = (f2[1]-flC*fl[1])*bFluxFitB*pars[3] - flC*fl[0]*aFluxFitB;


            //const double bNorm = fluxI(0.003, a0_R, ap_R/b0_R);
            if(type == "tInt_Rfix") p.sigma -= f2[1]*bFlux*pars[2];
            //else                    p.sigma -= f2[1]*bFlux*pars[3];
            else                    p.sigma -= corr;

            data.push_back(p.sigma);
            err.push_back(p.err);
            aVec.push_back(aFlux);
            bVec.push_back(bFlux);
        }
        if(type.Contains("_PR"))
            nDf += max(int(data.size() - 2), 0);
        else {
            nDf += max(int(data.size() - 1), 0);
            if(data.size() - 1 > 0) {
                ++nGroups;
                nPoints += data.size();
            }
        }
            //nDf += max(int(data.size() ), 0);

        //cout << "Chi2Now " << chi2Raw(data, err,  aVec, bVec) << endl;
        chi2 += chi2Raw(data, err,  aVec);
    }
    cout <<"chi2= "<< chi2 << " " << nDf << endl;
    cout <<"nGroups=" << nGroups << " nPoints "<< nPoints << endl;
    //exit(0);
    return chi2;
}



void fit(vector<point> &pointsComb)
{




}

struct GlobalChi2 {
    GlobalChi2()  {}
    GlobalChi2(TString Type, vector<point> points_) : type(Type), points(points_)
    {
        if(type == "tInt_PR" || type == "tDep_PR") {
            nFluxes = 2; 
            if(type.Contains("tInt")) nPars = 4;
            if(type.Contains("tDep")) nPars = 6;
        }
        else if(type == "tInt_P" || type == "tDep_P") {
            nFluxes = 1; 
            if(type.Contains("tInt")) nPars = 2;
            if(type.Contains("tDep")) nPars = 3;
        }
        else if(type == "tDep_Rfix") {
            nFluxes = 1; 
            nPars = 4;
        }
        else if(type == "tInt_Rfix") {
            nFluxes = 1;
            nPars = 3;
        }
        else {
            cout << "Not existing type : " << type << endl;
            exit(1);
        }

    }
    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)

    TString type;
    int nFluxes, nPars;
    vector<point> points;

    double operator() (const double *par) const {
        //cout <<"Hela " <<  par[0] << " "<< par[1] << " "<< par[2] << endl;
        //return hypot(par[0],par[1]) + hypot(par[1],par[2]);
        vector<double> pars(par, par + nPars);
        

        if(type.Contains("Rfix")) {
            return getChi2RegFixed(type, points,  pars);
        } 
        else {
            return getChi2(type, points, pars);
        }
    }
};


void plotBslope(vector<point> &points, vector<double> pars)
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
        
        map<double, TGraphErrors*> sigmaData, sigmaTh;

        //Loop over points with equal beta,q2 for data
        for(auto &p : ps) {

            if(sigmaData.count(p.xp) == 0) sigmaData[p.xp] = new TGraphErrors();
            if(sigmaTh.count(p.xp) == 0)   sigmaTh[p.xp] = new TGraphErrors();

            sigmaData.at(p.xp)->SetPoint(sigmaData.at(p.xp)->GetN(), p.tAbs, p.sigma);
            sigmaData.at(p.xp)->SetPointError(sigmaData.at(p.xp)->GetN()-1, 0, p.err);

            sigmaTh.at(p.xp)->SetPoint(sigmaTh.at(p.xp)->GetN(), p.tAbs, p.sigma);
            //sigmaData.at(p.xp)->SetPointError(sigmaData.at(p.xp)->GetN()-1, 0, p.err);

        }

        //Put the x-sections and errors into vector to get the normalization
        vector<double> sigVec, errVec, flxAvec, flxBvec;
        for(auto &p : ps) {
            double aFlux, bFlux;
            if(pars.size() <= 3) {

            }

            if(sigmaData.count(p.xp) == 0) sigmaData[p.xp] = new TGraphErrors();
            if(sigmaTh.count(p.xp) == 0)   sigmaTh[p.xp] = new TGraphErrors();

            sigmaData.at(p.xp)->SetPoint(sigmaData.at(p.xp)->GetN(), p.tAbs, p.sigma);
            sigmaData.at(p.xp)->SetPointError(sigmaData.at(p.xp)->GetN()-1, 0, p.err);

            sigmaTh.at(p.xp)->SetPoint(sigmaTh.at(p.xp)->GetN(), p.tAbs, p.sigma);
            //sigmaData.at(p.xp)->SetPointError(sigmaData.at(p.xp)->GetN()-1, 0, p.err);
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

//My xpom plot: Reproduction of the H1 result
void plotXpom(vector<point> &points, vector<double> pars = {})
{
    gRandom->SetSeed(0);
    gStyle->SetOptStat(0);

    //for(auto &p : pars) p = smear(p);

    map<pair<double,double>, vector<point>> pMap;

    for(auto &p : points) {
        pMap[{p.beta, p.q2}].push_back(p);
    }

    map<pair<double,double>,map<double,TGraphErrors*>> grDataMap;

    double chi2 = 0;
    int nDf = 0;
    for(auto &pm : pMap) {
        auto ps = pm.second; //points with equal beta, q2
        double q2, beta;
        tie(beta,q2) = pm.first;
        
        if(ps.size() < 2) continue;

        map<double,TGraphErrors*> grData;

        //Loop over points with equal beta,q2 for data
        for(auto &p : ps) {
            if(grData.count(p.tAbs) == 0) grData[p.tAbs] = new TGraphErrors();
            grData.at(p.tAbs)->SetPoint(grData.at(p.tAbs)->GetN(), p.xp, p.sigma);
            grData.at(p.tAbs)->SetPointError(grData.at(p.tAbs)->GetN()-1, 0, p.err);
        }
        assert(grDataMap.count(pm.first)==0);
        grDataMap[pm.first] = grData;
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
        cout << "Before" << endl;
        cout << "after" << endl;

        TH1D *hAx = new TH1D(rn(), "", 1, 0.0013, 0.15);
        hAx->Draw("axis");

        gPad->SetLogx();
        GetYaxis()->SetRangeUser(-0.01, 0.12);

        int ifit = 2;
        double beta = betaSet[ibeta];
        double q2 = q2Set[iq2];
        if(grDataMap.count({beta, q2}) > 0) {
            auto gr = grDataMap.at({beta, q2});

            vector<double> aVec, data, regg, err;
            //Fill data and theory points
            for(auto p : pMap.at({beta,q2})) { //points with equal beta, q2
                double aFlux = flux(p.xp,p.tAbs, pars[0], pars[1], pars[2]);

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
                err.push_back(errNow);
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

                auto leg2 = new TLegend(0.7, 0.835, 0.8, 0.935);        
                leg2->SetTextSize(PxFontToRel(20));
                leg2->SetBorderSize(0);
                leg2->SetHeader("H1 FPS");
                leg2->AddEntry(grTh.at(0.6), "Regge fit IP+IR", "l");
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












void reggeFit()
{

    //auto h1Comb = readInclusive("../datafiles/h1/1203.4495/h1_ddis_Comb.dat");
    //auto h1_225 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_225.dat");
    //auto h1_252 = readInclusive("../datafiles/h1/1107.3420/h1_ddis_lowE_252.dat");

    //auto h1fps  = readInclusive("../datafiles/h1/1010.1476/H1-DDIS-HERAII-FPS.dat");




    auto h1fps  = readInclusive("../datafiles/h1/1010.1476/H1-DDIS-HERAII-FPS-4D.dat");
    //return;
    //h1fps = subtractFL(h1fps);


    auto h1Data = h1fps;


    plotXpom(h1Data, {1.10, 0.04, 5.73, 0.87e-3});
    //plotXpom(h1Data, {1.09546 , -0.0126443 , 6.27418 , 0.000860536});
    //return;

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


    //GlobalChi2 globalChi2("tDep_Rfix", h1Data);
    GlobalChi2 globalChi2("tDep_P", h1Data);
    ROOT::Fit::Fitter fitter;

    vector<double> pars;
    if(globalChi2.type == "tInt_PR") {
        pars = { 1.1,0.02, 
                 0.5, 0.3};
    }
    else if(globalChi2.type == "tDep_PR") {
        pars = { 1.1,0.25, 6,
                 0.5, 0.5, 1.5};
    }
    else if(globalChi2.type == "tDep_P") {
        pars = { 1.05,0.04, 6 };
    }
    else if(globalChi2.type == "tInt_P") {
        pars = { 1.1, 0.02 };
    }
    else if(globalChi2.type == "tDep_Rfix") {
        pars = {1.10,0.04, 5.73, 0.001}; //a0, ap, b0, nReg
    }
    else if(globalChi2.type == "tInt_Rfix") {
        pars = {1.10,0.04/5.73, 0.001};
    }
    else {
        cout << "Issue " << __LINE__ << endl;
        exit(1);
    }

    fitter.Config().SetParamsSettings(pars.size(),pars.data());

    /*

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
        //fitter.Config().ParSettings(3).Fix();

    }
    else if(Npar == 3) {
        fitter.Config().ParSettings(0).SetLimits(0.2,5.5);
        fitter.Config().ParSettings(1).SetLimits(-0.71,0.73);
        fitter.Config().ParSettings(2).SetLimits(1,10);
        //fitter.Config().ParSettings(0).Fix();
    }

    */

    //return;


    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(pars.size(),globalChi2);//,0,dataB.Size()+dataSB.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
    result.PrintCovMatrix(cout);

    cout << "Zita: " << result.Value(0) << " "<< result.Value(1) <<" "<<  result.Value(2) <<" "<< result.Value(3) <<   " : ";
    cout << result.ParError(0) << " "<< result.ParError(1) <<" "<<  result.ParError(2) <<" " <<result.ParError(3) <<   endl;
    //cout << result.Chi2() << endl;
    for(int i = 0; i < pars.size(); ++i)
        cout << i <<" "<< result.GlobalCC(i) << endl;




    //plotBetaQ2(h1Comb, h1_225, h1_252, h1fps);
    //plotBetaQ2Xpom(h1Comb, h1_225, h1_252, h1fps);


}




//____________________________________________________________________________________ //
double rfluxRawInt(double a0, double ap, double b0,  double x_pom, double tAbsMin, double tAbsMax) {
   const double mp = 0.93827231;

   //     calc min. kinematically  allowed t
   double tAbsMinKin = pow(mp*x_pom,2)/(1.-x_pom);
   tAbsMin = std::max(tAbsMin, tAbsMinKin);
   assert(tAbsMin < tAbsMax);

   //     c*xpom**(-(2apom-1))
   double fl =  exp((2.0*a0-1.)*log(1.0/x_pom));
   double b=(b0+2.0*ap*log(1.0/x_pom));

   //   at fixed t:  exp(Bt)
   //  fl = fl * exp(b*tcut);

   //   t-integrated: (1/B)*[exp(-B*tmax)-exp(-B*tmin)]
   fl = fl * (exp(-tAbsMin*b)-exp(-tAbsMax*b))/b;
   if ( std::isinf(fl) || fl==0 || std::isnan(fl)) {
      std::cout<<"[rfluxRawInt] rflux is not a reasonable value: "<<fl<<"\t fl0 = "<<exp((2.0*a0-1.)*log(1.0/x_pom))<<", b="<<b<<", a0="<<a0<<std::endl;
   }
   return fl;
}

double rfluxRaw(double a0, double ap, double b0, double x_pom, double tAbs) {
   //     c*xpom**(-(2apom-1))
   double fl =  exp((2.0*a0-1.)*log(1.0/x_pom));
   double b=(b0+2.0*ap*log(1.0/x_pom));

   //   at fixed t:  exp(Bt)
   //  fl = fl * exp(b*tcut);
   fl = fl * exp(-b*tAbs);

   return fl;
}

double rfluxInt(double a0, double ap, double b0, double x_pom, double tAbsMin, double tAbsMax) {
   double tAbscutNorm = 1;
   double xPomNorm = 0.003;
   const double dm =  rfluxRawInt(a0, ap, b0, xPomNorm,  0, tAbscutNorm);
   double  norm=(1./(xPomNorm*dm)); //xpom * flux normalized to 1 at xpom = 0.003

   double rFlux = norm * rfluxRawInt(a0, ap, b0, x_pom, tAbsMin, tAbsMax);
   if ( std::isnan(rFlux) ) {
      std::cout<<"[rfluxInt] rFlux isnan: "<<rFlux<<". Input: a0="<<a0<<", ap="<<ap<<", b0="<<b0<<std::endl;
      std::cout<<"[rfluxInt] rFlux isnan.           dm="<<dm<<", xPomNorm="<<xPomNorm<<std::endl;
   }
   return rFlux;
}

double rflux(double a0, double ap, double b0, double x_pom, double tAbs) {
   double tAbscutNorm = 1;
   double xPomNorm = 0.003;
   const double dm =  rfluxRawInt(a0, ap, b0, xPomNorm,  0, tAbscutNorm);
   double  norm=(1./(xPomNorm*dm)); //xpom * flux normalized to 1 at xpom = 0.003

   double rFlux = norm * rfluxRaw(a0, ap, b0, x_pom, tAbs);
   if ( std::isnan(rFlux) ) {
      std::cout<<"[rfluxInt] rFlux isnan: "<<rFlux<<". Input: a0="<<a0<<", ap="<<ap<<", b0="<<b0<<std::endl;
      std::cout<<"[rfluxInt] rFlux isnan.           dm="<<dm<<", xPomNorm="<<xPomNorm<<std::endl;
   }
   return rFlux;
}

