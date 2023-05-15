import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from helpers import SplitUp
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Trijet'

class ttbarClass:
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

    #this is the end of initializer. Now we can apply preselection. Cpp modules will be compiled for each task, but we will not compile them in class function.

    def AddCutflowColumn(self, var, varName):
        print('Adding cutflow information....\n\t{}\t{}'.format(varName, var))
        self.a.Define('{}'.format(varName), str(var))

    # preselection we apply are the following:
    # we first want to mark all the event with:
    def Preselection(self):
        self.a.Cut('nFatJet','nFatJet > 0')# at least 1 AK8 jet
        self.a.Cut('nJet','nJet > 0') # at least 1 AK4 jet
        self.a.Cut('nLepton','nElectron > 0 || nMuon > 0') #make sure at least one lepton exist. Save some effort in c++ code
        self.a.Define('DijetIds','PickDijets(FatJet_phi,Jet_phi,Electron_pt,Muon_pt,Jet_btagCMVA)') #Output: Jet selection parameter of the form{FatJetId,JetId,Leptonid,Leptonpt,ElectronId,MuonId}. We demand lepton pt>50GeV, at least one AK4Jet(named Jet) is b-tagged.
        self.a.Cut('preselected','DijetIds[0]> -1 && DijetIds[1] > -1 && DijetIds[2] > -1') #Cut the data according to our standard (FatJet, Jet, Lepton condtion respectively)
        return self.a.GetActiveNode()
    
    #now we define the selection according to the following standard: top tagging AK8, and a 2D cut on lepton+b
    def Selection(self,Ttagparam):
        self.a.Cut('TopTagging','FatJet_particleNet_TvsQCD[DijetIds[0]] > {}'.format(Ttagparam))
        self.a.ObjectFromCollection('bJet','Jet','DijetIds[1]')#isolate the b jet for 2D cut analysis purposes
        self.a.Define('ElectronId','CreateIntColumn(4,DijetIds)')#create a column "electron id" based on the 4th element of DijetIds
        self.a.Define('MuonId','CreateIntColumn(5,DijetIds)')#create a column for muon as well.
        self.a.Define('Pick2DCut','TwoDCut(ElectronId,MuonId,Electron_jetPtRelv2,Muon_jetPtRelv2,bJet_phi,Electron_phi,Muon_phi)')#creat the boolian parameter for 2D cut using C++ script
        self.a.Cut('2DCut','Pick2DCut[0] == 1 || Pick2DCut[1] == 1')#if either condition is met, we keep the event.
        return self.a.GetActiveNode()
    
    #now we need to make the plot. For purpose of invariant mass reconstruction, we need to specify the lepton pt, eta, phi and mass manually. We will do this using a user defined C++ code.
    def JetsCandidateKinematicinfo(self):
        #first give relatvent information of lepton
        self.a.Define('LeptonId','CreateColumn(2,DijetIds)')#create the electron/muon indicator
        self.a.Define('Lepton_pt','GetFloatLeptonProperty(LeptonId,ElectronId,MuonId,Electron_pt,Muon_pt)')
        self.a.Define('Lepton_eta','GetFloatLeptonProperty(LeptonId,ElectronId,MuonId,Electron_eta,Muon_eta)')
        self.a.Define('Lepton_phi','GetFloatLeptonProperty(LeptonId,ElectronId,MuonId,Electron_phi,Muon_phi)')
        self.a.Define('Lepton_mass','GetFloatLeptonProperty(LeptonId,ElectronId,MuonId,Electron_mass,Muon_mass)')
        # we do the same for AK8/4 candidate
        self.a.ObjectFromCollection('Top','FatJet','DijetIds[0]')
        self.a.ObjectFromCollection('Bot','Jet','DijetIds[1]')

        #for Neutrino
        self.a.Define('Neutrino_pt','CreateFloatColumn(0,MET_pt)')
        self.a.Define('Neutrino_phi','CreateFloatColumn(0,MET_phi)')
        Neutrinomass = float(0)
        Neutrinoeta = float(0)
        self.a.Define('Neutrino_eta',str(Neutrinoeta))
        self.a.Define('Neutrino_mass',str(Neutrinomass))

        return self.a.GetActiveNode()
    
    
    #We are ready to make plots. current plots: leptonic candidate pt, hardonic candidate pt, inv mass of ttbar
    #hardronic candidate: easy. Just AK8 Jet mass. No need to do any combination whatsover
    #leptonic candidate: need information from btaggedAK4, Lepton, and MET. Assume 0 mass nutrino, moves in xy plane only so eta =0.

    def MassReconstruction(self):
        self.a.Define('Top_vect','hardware::TLvector(Top_pt, Top_eta, Top_phi, Top_msoftdrop)')
        self.a.Define('Bot_vect','hardware::TLvector(Bot_pt, Bot_eta, Bot_phi, Bot_msoftdrop)')
        self.a.Define('Lep_vect','hardware::TLvector(Lepton_pt, Lepton_eta, Lepton_phi, Lepton_mass)')
        self.a.Define('Neut_vect','hardware::TLvector(Neutrino_pt, Neutrino_eta, Neutrino_phi, Neutrino_mass)')
        self.a.Define('LepCandidate_pt','Bot_pt[0]+Lepton_pt[0]+Neutrino_pt[0]')
        self.a.Define('mttbar','hardware::InvariantMass({Top_vect, Bot_vect, Lep_vect, Neut_vect})')#invariant mass of the resonance particle


