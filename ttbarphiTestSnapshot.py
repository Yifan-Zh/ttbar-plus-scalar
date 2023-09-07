import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup, Correction
from ttbarphiclass import ttbarphiClass


start = time.time()

CompileCpp('ttbarphimodules.cc')

filename = 'WJetsHT400_17'

selection = ttbarphiClass('{}.txt'.format(filename),17,1,1)
selection.Preselection()
selection.Selection()
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.ApplyStandardCorrections(snapshot = True)
print('corrections are {}'.format(selection.GetXsecScale()))
selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale())

selection.Snapshot(signal = True)


print ('%s sec'%(time.time()-start))



histgroup = HistGroup('{}'.format(filename))

histList = []

#hist = selection.a.MakeTemplateHistos(ROOT.TH1F("PhiInvMass","Invariant Mass;m_{SD}(GeV);N_{Events}",100,0,100),'PhiInvMass')

h1 = selection.a.DataFrame.Histo1D(('PhiInvMass','Invariant Mass;m_{SD}(GeV);N_{Events}',100,0,100),'PhiInvMass')
histList.append(h1)

i1 = h1.Integral(0,100)
i2 = h1.Integral(3,7)

h2 = selection.a.DataFrame.Histo1D(('PhiLepton1_eta','Eta;eta;N_{Events}',100,0,10),'PhiLepton1_eta')
histList.append(h2)

h3 = selection.a.DataFrame.Histo1D(('PhiLepton2_eta','Eta;eta;N_{Events}',100,0,10),'PhiLepton2_eta')
histList.append(h3)

h4 = selection.a.DataFrame.Histo1D(('PhiLepton1_MotherType','Mother;MotherId;N_{Events}',30,0,15),'PhiLepton1_MotherType')
histList.append(h4)

h5 = selection.a.DataFrame.Histo2D(('LepCandidate1_pt_vs_DeltaR','Pt vs DeltaR;Pt_{GeV};DeltaR',20,0,60,40,0,4),'PhiLepton1_pt','PhiDeltaR')
histList.append(h5)

h6 = selection.a.DataFrame.Histo1D(('PhiLepton1_pt','Pt;pt;N_{Events}',100,0,100),'PhiLepton1_pt')
histList.append(h6)

h7 = selection.a.DataFrame.Histo1D(('PhiDeltaR','DeltaR;PhiDeltaR;N_{Events}',100,0,5),'PhiDeltaR')
histList.append(h7)

h8 = selection.a.DataFrame.Histo1D(('PhiLepton1_PfIso04','Isolation value;PhiLepton1_PfIso04;N_{Events}',100,0,2),'PhiLepton1_PfIso04')
histList.append(h8)

h9 = selection.a.DataFrame.Histo1D(('PhiInvMass0To10','Invariant Mass (0-10 GeV);m_{SD}(GeV);N_{Events}',200,0,10),'PhiInvMass')
histList.append(h9)

print('the total number of event is {}'.format(i1))
print('the effective number of event is {}'.format(i2))

#h2 = selection.a.DataFrame.Histo1D(('weight__Pileup_up','',100,0,100),'weight__Pileup_up')
#histList.append(h2)

#h3 = selection.a.DataFrame.Histo1D(('weight__Pileup_down','',100,0,100),'weight__Pileup_down')
#histList.append(h3)

#h4 = selection.a.DataFrame.Histo1D(('weight__nominal','',100,0,100),'weight__nominal')
#histList.append(h4)

'''outfile = ROOT.TFile.Open('rootfiles/kinDist_{}.root'.format(filename),'RECREATE')
outfile.cd()
hist.Do('Write')
outfile.Close()

tempfile = ROOT.TFile.Open('rootfiles/kinDist_{}.root'.format(filename))
h1 = tempfile.Get("PhiInvMass__nominal")
histList.append(h1)'''

h10 = selection.a.DataFrame.Histo1D(('WhichLepton','',100,0,10),'WhichLepton')
histList.append(h10)

c = ROOT.TCanvas('c','c')
c.cd()

for h in histList:
    c.Clear()
    name = h.GetName()
    h.Draw()
    c.Print('{}_{}.pdf'.format(filename,name))
    
