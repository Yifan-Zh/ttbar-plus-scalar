import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
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