import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup, Correction
from ttbarphiclass import ttbarphiClass


start = time.time()

CompileCpp('ttbarphimodules.cc')

filename = 'ttbarphi-signal'

selection = ttbarphiClass('{}.root'.format(filename),18,1,1)
selection.Preselection()
selection.Selection()
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.ApplyStandardCorrections(snapshot = True)
print('corrections are {}'.format(selection.GetXsecScale()))
selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale())

selection.Snapshot()


print ('%s sec'%(time.time()-start))



histgroup = HistGroup('{}'.format(filename))

histList = []

hist = selection.a.MakeTemplateHistos(ROOT.TH1F("PhiInvMass","Invariant Mass;m_{SD}(GeV);N_{Events}",100,0,100),'PhiInvMass')
#h1 = selection.a.DataFrame.Histo1D(('PhiInvMass','',100,0,100),'PhiInvMass')






outfile = ROOT.TFile.Open('rootfiles/kinDist_{}.root'.format(filename),'RECREATE')
outfile.cd()
hist.Do('Write')
outfile.Close()

tempfile = ROOT.TFile.Open('rootfiles/kinDist_{}.root'.format(filename))
h1 = tempfile.Get("PhiInvMass__nominal")
histList.append(h1)

h2 = selection.a.DataFrame.Histo1D(('WhichLepton','',100,0,10),'WhichLepton')
histList.append(h2)

c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}.pdf'.format(name))
