from argparse import Namespace
from glob import glob

from TIMBER.Tools.Common import DictStructureCopy, CompileCpp, ExecuteCmd, OpenJSON, StitchQCD
from TIMBER.Tools.Plot import CompareShapes
from TIMBER.Analyzer import Correction
import multiprocessing, ROOT, time

all_hist = {'bkg':{},'sig':{},'data':None}

for name in ['ttbarphi-signal','ttbar-semi','ttbar-had']:
    
    tempFile = ROOT.TFile.Open('rootfiles/kinDist_{}_17.root'.format(name))
    h = tempFile.Get('PhiInvMass')
    if name == 'ttbarphi-signal':
        all_hist['sig'][name] = h
        print('creating dictionary for {}, created under sig'.format(name))
    else:
        all_hist['bkg'][name] = h
        print('creating dictionary for {}, created under bkg'.format(name))

hist = {'bkg':{},'sig':{},'data':None}

tempfile1 = ROOT.TFile.Open('rootfiles/kinDist_ttbarphi-signal_17.root')
tempfile2 = ROOT.TFile.Open('rootfiles/kinDist_ttbar-had_17.root')
tempfile3 = ROOT.TFile.Open('rootfiles/kinDist_ttbar-semi_17.root')
tempfile4 = ROOT.TFile.Open('rootfiles/kinDist_QCDHT700_17.root')
tempfile5 = ROOT.TFile.Open('rootfiles/kinDist_QCDHT1000_17.root')
tempfile6 = ROOT.TFile.Open('rootfiles/kinDist_QCDHT1500_17.root')
tempfile7 = ROOT.TFile.Open('rootfiles/kinDist_QCDHT2000_17.root')
h1 = tempfile1.Get('PhiInvMass')
h2 = tempfile2.Get('PhiInvMass')
h3 = tempfile3.Get('PhiInvMass')
h4 = tempfile4.Get('PhiInvMass')
h5 = tempfile5.Get('PhiInvMass')
h6 = tempfile6.Get('PhiInvMass')
h7 = tempfile7.Get('PhiInvMass')
hist['sig']['ttbarphi-signal'] = h1
hist['bkg']['ttbar-had'] = h2
hist['bkg']['ttbar-semi'] = h3
hist['bkg']['QCD700'] = h4
hist['bkg']['QCD1000'] = h5
hist['bkg']['QCD1500'] = h6
hist['bkg']['QCD2000'] = h7

CompareShapes ('plots_PhiInvMass.pdf',17,'PhiInvMass[GeV]',
               bkgs=hist['bkg'],
               signals=hist['sig'],
               colors={'QCD700':ROOT.kOrange,'QCD1000':ROOT.kOrange-1,'QCD1500':ROOT.kOrange+1,'QCD2000':ROOT.kOrange-4,'ttbar-semi':ROOT.kRed-4,'ttbar-had':ROOT.kRed,'ttbarphi-signal':ROOT.kBlack},
               names={},
               scale = False, stackBkg = True,
               forceForward = True)



