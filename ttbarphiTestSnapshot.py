import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup
from ttbarphiclass import ttbarphiClass


start = time.time()

CompileCpp('ttbarphimodules.cc')

selection = ttbarphiClass('ttbarphi-all_18.txt',18,1,1)
selection.Preselection()
selection.Selection(0.8)
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.Snapshot()

print ('%s sec'%(time.time()-start))

histList = []

h1 = selection.a.DataFrame.Histo1D(('mttbar','',100,0,10000),'mttbar')
histList.append(h1)

c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}.pdf'.format(name))
