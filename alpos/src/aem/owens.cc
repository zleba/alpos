#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;


vector<double>  xfxQownens(double X, double SCALE)
{


    static const double COW[4][3][5] = {
    {
//...Expansion coefficients for up and down valence quark distributions.
     {  4.0000e-01,  7.0000e-01,  0.0000e+00,  0.0000e+00,  0.0000e+00},
     { -6.2120e-02,  6.4780e-01,  0.0000e+00,  0.0000e+00,  0.0000e+00},
     { -7.1090e-03,  1.3350e-02,  0.0000e+00,  0.0000e+00,  0.0000e+00}
    },
      //DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
//...Expansion coefficients for gluon distribution.
    {
     {  8.8800e-01,  0.0000e+00,  3.1100e+00,  6.0000e+00,  0.0000e+00},
     { -1.8020e+00, -1.5760e+00, -1.3170e-01,  2.8010e+00, -1.7280e+01},
     {  1.8120e+00,  1.2000e+00,  5.0680e-01, -1.2160e+01,  2.0490e+01}
    },
      //DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
//...Expansion coefficients for (up+down+strange) quark sea distribution.
    {
     {  9.0000e-01,  0.0000e+00,  5.0000e+00,  0.0000e+00,  0.0000e+00},
     { -2.4280e-01, -2.1200e-01,  8.6730e-01,  1.2660e+00,  2.3820e+00},
     {  1.3860e-01,  3.6710e-03,  4.7470e-02, -2.2150e+00,  3.4820e-01}
    },
     // DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
//...Expansion coefficients for charm quark sea distribution.
    {
     {  0.0000e+00, -2.2120e-02,  2.8940e+00,  0.0000e+00,  0.0000e+00},
     {  7.9280e-02, -3.7850e-01,  9.4330e+00,  5.2480e+00,  8.3880e+00},
     { -6.1340e-02, -1.0880e-01, -1.0852e+01, -7.1870e+00, -1.1610e+01}
    }
    }; 

//     +       COW(3,5,4),TS(6),XQ(9)

        double TS[6], XQ[9];


       double ALAM=0.2, Q02=4., QMAX2=2.e3;

//...Pion structure functions from Owens.
//...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

//...Determine set, Lambda and s expansion variable.
        double Q2 = SCALE*SCALE;
        double Q2IN = min( QMAX2,max( Q02,Q2));
        double SD = log( log( Q2IN/pow(ALAM,2))/ log( Q02/pow(ALAM,2)));

//...Calculate structure functions.
        for(int KFL = 0; KFL < 4; ++KFL) {
            for(int IS = 0; IS < 5; ++IS)
                TS[IS]=COW[KFL][0][IS]+COW[KFL][1][IS]*SD+ COW[KFL][2][IS]*SD*SD;
            if(KFL == 0) {

                double DENOM = tgamma(TS[0])*tgamma(TS[1]+1.0)/ tgamma(TS[0]+TS[1]+1.0);
                XQ[KFL]=pow(X,TS[0]) * pow(1.-X,TS[1]) /DENOM;
            } else {
                XQ[KFL]=TS[0] * pow(X,TS[1]) * pow(1.-X,TS[2]) * (1.+TS[3]*X+TS[4]* pow(X,2));
            }
        }

//...Put into output arrays.
    double    UPV = XQ[0];
    double    DNV = XQ[0];
    double    SEA = XQ[2]/6.;
    double    STR = XQ[2]/6.;
    double    CHM = XQ[3];
    double    BOT = 0.0;
    double    TOP = 0.0;
    double    GL  = XQ[1];
    vector<double> xqx = {TOP, BOT, CHM, STR, SEA, SEA, GL, SEA+DNV, SEA+UPV, STR, CHM, BOT, TOP};
    return xqx;
}

/*
#include <iostream>
int main()
{

    double SCALE = 2, X = 0.1;
    double w = 1e-5;
    double f = 0;
    for(X = 1e-5; X < 1; X += w) {
        double UPV, DNV, SEA, STR, CHM, GL;
        STROWP1(X, SCALE, UPV, DNV, SEA, STR, CHM, GL);
        //cout << X << " "<< CHM << endl;
        f += (3*SEA + (UPV+SEA) + (DNV+SEA) + STR + 2*CHM + GL) * w;
        //f += CHM * w;
    }
    cout << f << endl;
    
    double UPV, DNV, SEA, STR, CHM, GL;
    X = 0.1;
    SCALE = 10;
    STROWP1(X, SCALE, UPV, DNV, SEA, STR, CHM, GL);
    cout << UPV <<" " << DNV <<" "<< SEA <<" " << STR <<" "<< CHM <<" "<< GL << endl;
}
*/
