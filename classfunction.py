import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from helpers import SplitUp
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Dijet'

class TbWClass:
    def __init__(self,inputfile,year,ijob,njobs):#initializer
        if inputfile.endswith('.txt'): 
            infiles = SplitUp(inputfile, njobs)[ijob-1]
	    # there is an issue with empty events TTrees, so make sure they don't make it through to the analyzer (mainly seen in V+Jets, esp at low HT)
            invalidFiles = []
            for iFile in infiles:
		#print('Adding {} to Analyzer'.format(iFile))
                f = ROOT.TFile.Open(iFile)
                if not f.Get('Events'):
                    print('\tWARNING: {} has no Events TTree - will not be added to analyzer'.format(iFile))
                    invalidFiles.append(iFile)
                    continue
                if not f.Get('Events').GetEntry():
                    print('\tWARNING: {} has an empty Events TTree - will not be added to analyzer'.format(iFile))
                    invalidFiles.append(iFile)
                f.Close()
            inputFiles = [i for i in infiles if i not in invalidFiles]
            if len(inputFiles) == 0:
                print("\n\t WARNING: None of the files given contain an Events TTree.")
            self.a = analyzer(inputFiles)
        else:
            infiles = inputfile
            self.a = analyzer(infiles)
            #self.a will be the Timber analyzer class
        if inputfile.endswith('.txt'):
            self.setname = inputfile.split('/')[-1].split('_')[0]
        else:
            self.setname = inputfile.split('/')[-1].split('_')[1]
        self.year = str(year)	# most of the time this class will be instantiated from other scripts with CLI args, so just convert to string for internal use
        self.ijob = ijob
        self.njobs = njobs
        self.config = OpenJSON('THconfig.json')
        self.cuts = self.config['CUTS']
        self.newTrigs = self.config['TRIGS']	
        self.trigs = {
	    16:['HLT_PFHT800','HLT_PFHT900'],
            17:["HLT_PFHT1050","HLT_AK8PFJet500","HLT_AK8PFHT750_TrimMass50","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet400_TrimMass30"],
            19:['HLT_PFHT1050','HLT_AK8PFJet500'], # just use 19 for trigger script for 17b, 17all
            #18:["HLT_PFHT1050","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet500","HLT_AK8PFJet400_TrimMass30","HLT_AK8PFHT750_TrimMass50"]
            18:['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFHT850_TrimMass50','HLT_PFHT1050']
        }

        if 'Data' in inputfile:		# SingleMuonDataX_year and DataX_year are possible data inputfile names
            self.a.isData = True
        else:
            self.a.isData = False

    #this is the end of initializer. Now we can apply preselection:
    # preselection we apply are the following:
    # we first want to mark all the event with:
    def Preselection(self):
        self.a.Cut('nFatJet','nFatJet > 0')# at least 1 AK8 jet
        self.a.Cut('nJet','nJet > 0') # at least 1 AK4 jet
        self.a.Cut('nLepton','nElectron > 0 || nMuon > 0') #make sure at least one lepton exist. Save some effort in c++ code
        self.a.Define('DijetIds','PickDijets(FatJet_phi,Jet_phi,Electron_pt,Muon_pt,Jet_btagCMVA)') #Jet selection parameter creation{FatJetId,JetId,Leptonid,Leptonpt,ElectronId,MuonId}. We demand lepton pt>50GeV, at least one Jet is b-tagged.
        self.a.Cut('preselected','DijetIds[0]> -1 && DijetIds[1] > -1 && DijetIds[2] == 1') #Cut the data according to our standard.
        return self.a.GetActiveNode()
    
    #now we define the selection according to the following standard: top tagging AK8, and a 2D cut on lepton+b
    def Selection(self,Ttag):
        self.a.Cut('TopTagging','FatJet_deepTagMD_TvsQCD[DijetIds[0]] > {}'.format(Ttag))
        self.a.ObjectFromCollection('bJet','Jet','DijetIds[1]')#isolate the b jet for 2D cut analysis purposes
        self.a.ObjectFromCollection('Lep','')#how to make a column for lepton momentum?
        return self.a.GetActiveNode()

        
