import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from helpers import SplitUp
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Dijet'

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
        self.a.Define('DijetIds','PickDijetsV2(FatJet_phi,Jet_phi,Electron_pt,Muon_pt,Jet_btagCSVV2)') #Output: Jet selection parameter of the form{FatJetId,JetId,Leptonid,Leptonpt,ElectronId,MuonId}. We demand lepton pt>50GeV, at least one AK4Jet(named Jet) is b-tagged.
        self.a.Cut('preselected','DijetIds[0]> -1 && DijetIds[1] > -1 && DijetIds[2] > -1') #Cut the data according to our standard (FatJet, Jet, Lepton condtion respectively)
        self.a.Define('bJetFromJets','DijetIds[1]')#take a look at which jet is being selected as the bjet
        return self.a.GetActiveNode()
    
    #now we define the selection according to the following standard: top tagging AK8, and a 2D cut on lepton+b
    def Selection(self,Ttagparam):
        self.a.Cut('TopTagging','FatJet_particleNet_TvsQCD[DijetIds[0]] > {}'.format(Ttagparam))
        self.a.ObjectFromCollection('bJet','Jet','DijetIds[1]')#isolate the b jet for 2D cut analysis purposes
        self.a.Define('Pick2DCut','TwoDCutV2(DijetIds[2],DijetIds[4],DijetIds[5],Electron_pt,Muon_pt,bJet_pt,Electron_phi,Muon_phi,bJet_phi,Electron_eta,Muon_eta,bJet_eta)')
        self.a.Cut('2DCut','Pick2DCut[0] == 1 || Pick2DCut[1] == 1')#if either condition is met, we keep the event.
        return self.a.GetActiveNode()
    
    #now we need to make the plot. For purpose of invariant mass reconstruction, we need to specify the lepton pt, eta, phi and mass manually. We will do this using a user defined C++ code.
    def JetsCandidateKinematicinfo(self):
        #first give relatvent information of lepton; do not use Lepton_*, will cause a bug in snapshot
        self.a.Define('MyLepton_id','DijetIds[2]')
        self.a.Define('MyLepton_pt','GetFloatLeptonProperty(DijetIds[2],DijetIds[4],DijetIds[5],Electron_pt,Muon_pt)')
        self.a.Define('MyLepton_eta','GetFloatLeptonProperty(DijetIds[2],DijetIds[4],DijetIds[5],Electron_eta,Muon_eta)')
        self.a.Define('MyLepton_phi','GetFloatLeptonProperty(DijetIds[2],DijetIds[4],DijetIds[5],Electron_phi,Muon_phi)')
        self.a.Define('MyLepton_mass','GetFloatLeptonProperty(DijetIds[2],DijetIds[4],DijetIds[5],Electron_mass,Muon_mass)')
        # we do the same for AK8/4 candidate
        self.a.ObjectFromCollection('Top','FatJet','DijetIds[0]')
        self.a.ObjectFromCollection('Bot','Jet','DijetIds[1]')

        #debug purpose:
        #self.a.Cut('BotMassCut','Bot_mass < 30')
        #self.a.Cut('BotIdCut,','DijetIds[1] > 1')#this is based on observation that heavy tail jets in mass distribution are mostly 0th or 1st energetic jet in the system.
        #more debuging:
        # find all b quarks, store their mothers' indices and pt
        self.a.Define("GenB_genPartIdxMother","GenPart_genPartIdxMother[GenPart_pdgId == 5 or GenPart_pdgId == -5]")
        self.a.Define("GenB_pt","GenPart_pt[GenPart_pdgId ==5 or GenPart_pdgId == -5]")
        #find mother's PDG IDs
        self.a.Define("GenB_pdgIdMother","FindMothersPdgId(GenPart_pdgId,GenB_genPartIdxMother)")
        self.a.Define("GenBfromT_pt","GenB_pt[GenB_pdgIdMother==6]")


        #for Neutrino:note, the simple method, assuming eta=0 will not work (because it is not) Need to solve conservation of 3 component of 4-vector
        self.a.Define('Neutrino_pt','MET_pt')
        self.a.Define('Neutrino_phi','MET_phi')
        #self.AddCutflowColumn(float(0.0),'Neutrino_eta')
        self.a.Define('Neutrino_eta','NeutrinoEta(MyLepton_pt,MyLepton_phi,MyLepton_eta,MET_pt,MET_phi)')#if someone is reading this, ask Amitav for the paper.
        self.AddCutflowColumn(float(0.0),'Neutrino_mass')

        #it might help to check how many leptons are left in the event
        self.a.Define('nResultLepton','nElectron + nMuon')

        return self.a.GetActiveNode()
    
    
    #We are ready to make plots. current plots: leptonic candidate pt, hardonic candidate pt, inv mass of ttbar
    #hardronic candidate: easy. Just AK8 Jet mass. No need to do any combination whatsover
    #leptonic candidate: need information from btaggedAK4, Lepton, and MET. Assume 0 mass nutrino, moves in xy plane only so eta =0.

    def MassReconstruction(self):
        self.a.Define('Top_vect','hardware::TLvector(Top_pt, Top_eta, Top_phi, Top_msoftdrop)')
        self.a.Define('Bot_vect','hardware::TLvector(Bot_pt, Bot_eta, Bot_phi, Bot_mass)')
        self.a.Define('Lep_vect','hardware::TLvector(MyLepton_pt, MyLepton_eta, MyLepton_phi, MyLepton_mass)')
        self.a.Define('Neut_vect','hardware::TLvector(Neutrino_pt, Neutrino_eta, Neutrino_phi, Neutrino_mass)')
        self.a.Define('mttbar','hardware::InvariantMass({Top_vect, Bot_vect, Lep_vect, Neut_vect})')#invariant mass of the resonance particle
        #for calculate the leptonic candidate, will need phi value of b quark, lepton, and neutrino
        self.a.Define('LepCandidate_mass','hardware::InvariantMass({Bot_vect,Lep_vect,Neut_vect})')#ignore Neutrino for now
        self.a.Define('LepCandidate_pt','LeptonicCandidatePt(Bot_pt, MyLepton_pt, Bot_phi, MyLepton_phi, Neutrino_pt, Neutrino_phi)')# the total transverse momentum of pt candidate?
        self.a.Define('Total_pt','TotalPt(Top_pt,Top_phi,Bot_pt,Bot_phi,MyLepton_pt,MyLepton_phi, Neutrino_pt, Neutrino_phi)')#total pt of hadronic top and leptonic top
        self.a.Define('deltaMass','abs(Top_msoftdrop-LepCandidate_mass)')#check the difference in invariant mass of two top. They should not be very different.
        #Debug purpose: get rid of all events where mass difference in two tops are greater than 200GeV
        self.a.Cut('DeltaMassCut','deltaMass < 100')
        return self.a.GetActiveNode()
    
    def Snapshot(self,node=None,colNames=[]):
        startNode = self.a.GetActiveNode()
        if node == None: node = self.a.GetActiveNode()
        #colNames[str]:give what variales to keep at the snapshot

        columns = [
            'Top_pt','Top_msoftdrop','mttbar','LepCandidate_pt','LepCandidate_mass',
            'Bot_mass','Bot_pt','Bot_eta','nResultLepton','Neutrino_pt',
            'MyLepton_pt','MyLepton_mass','MyLepton_id','MyLepton_eta',
            'bJetFromJets','GenB_genPartIdxMother','Total_pt',
            'deltaMass'
        ]

        if (len(colNames) > 0):
            columns.extend(colNames)

        self.a.SetActiveNode(node)
        self.a.Snapshot(columns,'ttbarsnapshot_%s_%s_%sof%s.root'%(self.setname,self.year,self.ijob,self.njobs),'Events',openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

        


