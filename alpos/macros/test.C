R__LOAD_LIBRARY($PROJECT_DIR/alposBuild/libaem.so)

extern "C" void qcd_2006_(double *z,double *q2, int *ifit, double *xPq,       double *f2, double *fl, double *c2, double *cl);
extern "C" void h12006pdf_(double *z, double *q2, int *ifit, int *ipdf, double *xpq, double *f2, double *fl, double *c2, double *cl);

void test()
{

    double xPq[13], f2[2], fl[2], c2[2], cl[2];
    double z = 0.49995, vTemp;
    double q2 = 8.5;
    int ifit = 1;

    ifit = 1;

    for(int ipdf = 0; ipdf <= 32; ++ipdf) {
        h12006pdf_(&z,&q2, &ifit, &ipdf, xPq,f2,fl,c2,cl);
        cout << xPq[7] << endl;
    }

    qcd_2006_(&z,&q2, &ifit, xPq,    f2, fl, c2, cl);
    cout << xPq[7] << endl;

}
