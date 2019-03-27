from ROOT import gStyle, TCanvas, TH1D, Form, gSystem, TGraph, gPad, kBlue
import os
gSystem.Load(os.environ['PlH_DIR'] +'/plottingHelper_C.so')
from ROOT import PlottingHelper as ph


import subprocess
tab = subprocess.Popen("cd /afs/desy.de/user/z/zlebcr/h1/diff/alpos/farm/variants/AExt_nnlo_heraC_scan.str_dir/logs/ && grep ' | Summary: SuperTheory         ' * | sort -n | sed 's/\.out:/ /' | awk '{print $1, $5, $6}'", shell=True, stdout=subprocess.PIPE).stdout.read()

nums = []

q2 = [3.5, 4.0, 5.0, 6.5, 8.5, 11.5, 12.0, 15.0, 20.0, 25.0, 35.0, 44, 45.0, 60.0, 90.0, 200.0, 400.0, 800.0, 1600.0,]
#q2 = [3.5, 5.0, 6.5, 8.5, 12.0, 15.0, 20.0, 25.0, 35.0, 45.0, 60.0, 90.0, 200.0, 400.0, 800.0, 1600.0,]

gr    = TGraph()
grNdf = TGraph()
for line in tab.splitlines():
    l = line.split()
    Id = int(l[0])
    chi2 = float(l[1])
    ndf = float(l[2])

    print q2[Id], chi2, ndf
    gr.SetPoint(gr.GetN(), q2[Id], chi2/max(ndf-8, 1))
    grNdf.SetPoint(grNdf.GetN(), q2[Id], ndf-8)





gStyle.SetOptStat(0)

def vec(vv):
    from ROOT import std
    vvv = std.vector("double")()
    for v in vv:
        vvv.push_back(v)
    return vvv


can = TCanvas("can", "", 500, 600)
ph.SetLeftRight(0.23, 0.05)
ph.DivideTransparent(vec([1]), vec([1,0,0.8]))
#hh = []

can.cd(1)
h = TH1D("p1",";;", 1, 3, 100)
h.Draw("axis")
gr.Draw("l* same")
gr.SetLineWidth(2)
gr.SetLineColor(kBlue)
gPad.SetLogx()
ph.GetYaxis().SetRangeUser(0.01, 2)

ph.GetYaxis().SetTitle("#chi^{2}/ndf")
ph.GetYaxis().SetNdivisions(305)

ph.SetFTO(vec([26]), vec([15]), vec([1.2, 2, 0.3, 3.5]))
ph.GetXaxis().SetLabelOffset(1000)

can.cd(2)
h2 = TH1D("p2",";;", 1, 3, 100)
h2.Draw("axis")
grNdf.Draw("l* same")
grNdf.SetLineWidth(2)
grNdf.SetLineColor(kBlue)

ph.SetFTO(vec([26]), vec([15]), vec([1.2, 2, 0.3, 3.5]))

gPad.SetLogx()
gPad.SetLogy()
ph.GetYaxis().SetRangeUser(9, 400.)
ph.GetXaxis().SetNoExponent()
ph.GetXaxis().SetMoreLogLabels()
ph.GetYaxis().SetNoExponent()
ph.GetYaxis().SetMoreLogLabels()
ph.GetYaxis().SetTitle("ndf")
ph.GetXaxis().SetTitle("Q^{2}_{cut} [GeV]")


#ph.DrawLatexUp(can.GetPad(1), can.GetPad(5), 2, "This is a testing grid in Python");
can.SaveAs("chi2Ndf.pdf");

