import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from classfunction import ttbarClass

CompileCpp('ttbarphimodules.cc')

print('CompileCompleted. Check if error occurs')