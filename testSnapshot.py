import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup
from classfunction import ttbarClass


start = time.time()

CompileCpp('ttbarmodules.cc')

selection = ttbarClass('ttbar-semilep_17.root',17,1,1)
selection.Preselection()
selection.Selection(0.8)
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.Snapshot()

print ('%s sec'%(time.time()-start))

histList = []

h1 = selection.a.DataFrame.Histo1D(('mttbar','',50,0,3000),'mttbar')
histList.append(h1)

h2 = selection.a.DataFrame.Histo1D(('Top_msoftdrop','',50,0,1000),'Top_msoftdrop')
histList.append(h2)

h3 = selection.a.DataFrame.Histo1D(('LepCandidate_mass','',50,0,1200),'LepCandidate_mass')
histList.append(h3)

h4 = selection.a.DataFrame.Histo1D(('Top_pt','',50,0,1500),'Top_pt')
histList.append(h4)

h5 = selection.a.DataFrame.Histo1D(('LepCandidate_pt','',50,0,1500),'LepCandidate_pt')
histList.append(h5)

h6 = selection.a.DataFrame.Histo1D(('Lepton_id','',5,0,5),'MyLepton_id')
histList.append(h6)

h7 = selection.a.DataFrame.Histo1D(('Lepton_mass','',100,-0.3,0.3),'MyLepton_mass')
histList.append(h7)

h8 = selection.a.DataFrame.Histo1D(('Lepton_pt','',50,0,1000),'MyLepton_pt')
histList.append(h8)

h9 = selection.a.DataFrame.Histo1D(('Bot_mass','',100,0,200),'Bot_mass')
histList.append(h9)

h10 = selection.a.DataFrame.Histo1D(('Lepton_eta','',20,-5,5),'MyLepton_eta')
histList.append(h10)

h11 = selection.a.DataFrame.Histo1D(('Bot_eta','',20,-5,5),'Bot_eta')
histList.append(h11)

h12 = selection.a.DataFrame.Histo1D(('Bot_pt','',50,0,1000),'Bot_pt')
histList.append(h12)

h13 = selection.a.DataFrame.Histo1D(('Neutrino_pt','',50,0,1000),'Neutrino_pt')
histList.append(h13)

h14 = selection.a.DataFrame.Histo1D(('bJetFromJets','',50,0,10),'bJetFromJets')
histList.append(h14)

h15 = selection.a.DataFrame.Histo1D(('GenB_genPartIdxMother','',40,0,20),'GenB_genPartIdxMother')
histList.append(h15)

h16 = selection.a.DataFrame.Histo1D(('Total_pt','',50,0,1000),'Total_pt')
histList.append(h16)

h17 = selection.a.DataFrame.Histo1D(('Mass Difference in two top candidate','',50,0,1000),'deltaMass')
histList.append(h17)


c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}.pdf'.format(name))

