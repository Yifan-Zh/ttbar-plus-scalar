import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup
from ttbarphiclass import ttbarphiClass


start = time.time()

CompileCpp('ttbarphimodules.cc')

filename = 'QCDHT2000_17'

selection = ttbarphiClass('{}.txt'.format(filename),17,1,1)
selection.Preselection()
selection.Selection()
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.Snapshot()

print ('%s sec'%(time.time()-start))

histgroup = HistGroup('{}'.format(filename))

histList = []

h1 = selection.a.DataFrame.Histo1D(('PhiInvMass','',100,0,100),'PhiInvMass')
histList.append(h1)

h1.GetValue()

histgroup.Add('PhiInvMass',h1)

outfile = ROOT.TFile.Open('rootfiles/kinDist_{}.root'.format(filename),'RECREATE')
outfile.cd()
histgroup.Do('Write')
outfile.Close()

h2 = selection.a.DataFrame.Histo1D(('WhichLepton','',100,0,10),'WhichLepton')
histList.append(h2)

c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}.pdf'.format(name))
