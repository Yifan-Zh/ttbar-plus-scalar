import ROOT, time
ROOT.gROOT.SetBatch(True)
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from TIMBER.Analyzer import HistGroup, Correction
from ttbarphiclass import ttbarphiClass

parser = ArgumentParser()
parser.add_argument('-s', type=str, dest='setname',
                    action='store', required=True,
                    help='Setname to process.')
parser.add_argument('-y', type=str, dest='era',
                    action='store', required=True,
                    help='Year of set (16, 17, 18).')
parser.add_argument('-j', type=int, dest='ijob',
                    action='store', default=1,
                    help='Job number')
parser.add_argument('-n', type=int, dest='njobs',
                    action='store', default=1,
                    help='Number of jobs')
args = parser.parse_args()
start = time.time()

CompileCpp('ttbarphimodules.cc')


selection = ttbarphiClass('raw_nano/{}_{}.txt'.format(args.setname,args.era),args.era,args.ijob,args.njobs)
selection.Preselection()
selection.Selection()
selection.JetsCandidateKinematicinfo()
selection.MassReconstruction()
selection.ApplyStandardCorrections(snapshot = True)
print('corrections are {}'.format(selection.GetXsecScale()))
selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale())

selection.Snapshot()


print ('%s sec'%(time.time()-start))
    
