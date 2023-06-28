import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup
from ttbarphiclass import ttbarphiClass


start = time.time()

CompileCpp('ttbarphimodules.cc')

selection = ttbarphiClass('ttbar-semilep_17.root',17,1,1)
selection.Preselection()
selection.Selection(0.8)
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.Snapshot()

print ('%s sec'%(time.time()-start))

histList = []

h1 = selection.a.DataFrame.Histo1D(('mttbar','',50,0,3000),'mttbar')
histList.append(h1)

c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}.pdf'.format(name))
