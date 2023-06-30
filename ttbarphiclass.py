import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from helpers import SplitUp
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Dijet'

class ttbarphiClass:
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
        self.a.Define('DijetIds','PickDijetsV2(FatJet_phi,Jet_phi,Jet_btagCSVV2)') #Output: Jet selection parameter of the form{FatJetId,JetId}. We demand at least one AK4Jet is b-tagged.
        self.a.Cut('preselected','DijetIds[0]> -1 && DijetIds[1] > -1') #Cut the data according to our standard (FatJet, Jet condtion respectively)

        self.a.Define('nTotalLepton','nElectron + nMuon')#since we have a light scalar decay into two lepton, we should have at least 3 lepton
        self.a.Cut('LeptonNumberCut','nTotalLepton > 2')
        # we would find the leading leptons according to their pt. The output has the form {Electron/Muon(represented by 1/2), relative postion insde the corresponding vector}
        #for example, if the leading leptons are Electron[2],Muon[3],Electroon[4], then it would be {1,2,1,2,3,4}
        self.a.Define('LeadingThreeLepton','FindLeadLepton(Electron_pt,Muon_pt)')
        # make sure the least energetic one have at least 50 GeV
        self.a.Define('LeptonMinPtConstriant','MinPtConstraint(Electron_pt,Muon_pt,LeadingThreeLepton[2],LeadingThreeLepton[5])')
        self.a.Cut('PreselectionPtCut','LeptonMinPtConstriant == 1')
        self.a.DataFrame.Count().GetValue()
        print ("Pass Preselection stage")
        return self.a.GetActiveNode()
    
    #now we define the selection according to the following standard: top tagging AK8, identify which lepton come from t decay and a 2D cut on lepton+b
    #the 2 general cases are: 3 same lepton/1 electron 2 muon or vice versa
    #in the second case, the singled one must come from a leptonic top decay, and thus must satisfy the 2Dcut standard. The other two should move in the same direction and close to AK8/4 Jets with a relatively small invariant mass
    #in the first case, it's more complicated, we'd like to have a couple conditions satisfied simultaneously:
    # the particle-antiparticle pair cloest to each other will be reconstructed as phi, they need to satisfy the same condition as in the second case.
    # we want the output to carry the following information:
    # {pass selection or not, which one comes from top, which ones come from phi}
    #the "close to jet test" will be postponed to after reconstruction. Otherwise we will have to reconstruct the object twice.
    def Selection(self,Ttagparam):
        self.a.Cut('TopTagging','FatJet_particleNet_TvsQCD[DijetIds[0]] > {}'.format(Ttagparam))
        self.a.ObjectFromCollection('bJet','Jet','DijetIds[1]')#isolate the 2 jets for 2D cut analysis purposes
        self.a.ObjectFromCollection('Top','FatJet','DijetIds[0]')
        self.a.DataFrame.Count().GetValue()
        print ("Pass TopTagging")
        #now we start to handle the leptons. We'll handle this part in c++
        self.a.Define('LeptonTestAndReOrdering','LeptonCategorize(LeadingThreeLepton,Electron_pt,Muon_pt,bJet_pt,Electron_phi,Muon_phi,bJet_phi,Electron_eta,Muon_eta,bJet_eta,Electron_charge,Muon_charge)')
        self.a.Cut('PassAllExceptJetRelAngle','LeptonTestAndReOrdering[0] == 1')
        self.a.DataFrame.Count().GetValue()
        print ("Pass selection stage")
        return self.a.GetActiveNode()
    
    #now we need to make the plot. For purpose of invariant mass reconstruction, we need to specify the lepton pt, eta, phi and mass manually. We will do this using a user defined C++ code.
    def JetsCandidateKinematicinfo(self):
        #first give relatvent information of lepton; do not use Lepton_*, will cause a bug in snapshot
        #we now define the kinematics variables of the three lepton. The one from W followed by ones from phi
        self.a.Define('WLepton_pt','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + 3],Electron_pt,Muon_pt)')
        self.a.Define('WLepton_eta','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + 3],Electron_eta,Muon_eta)')
        self.a.Define('WLepton_phi','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + 3],Electron_phi,Muon_phi)')
        self.a.Define('WLepton_mass','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + 3] ,Electron_mass,Muon_mass)')
        #for the pairs:
        self.a.Define('PhiLepton1_pt','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + 3],Electron_pt,Muon_pt)')
        self.a.Define('PhiLepton1_eta','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + 3],Electron_eta,Muon_eta)')
        self.a.Define('PhiLepton1_phi','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + 3],Electron_phi,Muon_phi)')
        self.a.Define('PhiLepton1_mass','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + 3] ,Electron_mass,Muon_mass)')

        self.a.Define('PhiLepton2_pt','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[3]],LeadingThreeLepton[LeptonTestAndReOrdering[3] + 3],Electron_pt,Muon_pt)')
        self.a.Define('PhiLepton2_eta','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[3]],LeadingThreeLepton[LeptonTestAndReOrdering[3] + 3],Electron_eta,Muon_eta)')
        self.a.Define('PhiLepton2_phi','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[3]],LeadingThreeLepton[LeptonTestAndReOrdering[3] + 3],Electron_phi,Muon_phi)')
        self.a.Define('PhiLepton2_mass','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[3]],LeadingThreeLepton[LeptonTestAndReOrdering[3] + 3] ,Electron_mass,Muon_mass)')

        # we do the same for b quark jet candidate
        self.a.ObjectFromCollection('Bot','Jet','DijetIds[1]')


        #for Neutrino:note, the simple method, assuming eta=0 will not work (because it is not) Need to solve conservation of 3 component of 4-vector
        self.a.Define('Neutrino_pt','MET_pt')
        self.a.Define('Neutrino_phi','MET_phi')
        self.a.Define('Neutrino_eta','NeutrinoEta(WLepton_pt,WLepton_phi,WLepton_eta,MET_pt,MET_phi)')#if someone is reading this, ask Amitav for the paper.
        self.AddCutflowColumn(float(0.0),'Neutrino_mass')

        return self.a.GetActiveNode()
    
    
    #We are ready to make plots. current plots: leptonic candidate pt, hardonic candidate pt, inv mass of ttbar
    #hardronic candidate: easy. Just AK8 Jet mass. No need to do any combination whatsover
    #leptonic candidate: need information from btaggedAK4, Lepton, and MET. Assume 0 mass nutrino, moves in xy plane only so eta =0.

    def MassReconstruction(self):
        self.a.Define('HadronicTop_vect','hardware::TLvector(Top_pt, Top_eta, Top_phi, Top_msoftdrop)')
        self.a.Define('Bot_vect','hardware::TLvector(Bot_pt, Bot_eta, Bot_phi, Bot_mass)')
        self.a.Define('WLep_vect','hardware::TLvector(WLepton_pt, WLepton_eta, WLepton_phi, WLepton_mass)')
        self.a.Define('PhiLep1_vect','hardware::TLvector(PhiLepton1_pt,PhiLepton2_eta,PhiLepton1_phi,PhiLepton1_mass)')
        self.a.Define('PhiLep2_vect','hardware::TLvector(PhiLepton2_pt,PhiLepton2_eta,PhiLepton2_phi,PhiLepton2_mass)')
        self.a.Define('Neut_vect','hardware::TLvector(Neutrino_pt, Neutrino_eta, Neutrino_phi, Neutrino_mass)')

        #Cut if the reconstructed phi is far away from both jets
        self.a.Define('LeptonicTop_vect','Bot_vect + WLep_vect + Neut_vect')#Neutrino removed
        self.a.Define('Phi_vect','PhiLep1_vect + PhiLep2_vect')
        self.a.Cut('CloseToEitherOfTop','(abs(hardware::DeltaPhi(Phi_vect.Phi(),LeptonicTop_vect.Phi())) < 0.785) || (abs(hardware::DeltaPhi(Phi_vect.Phi(),HadronicTop_vect.Phi())) < 0.785)')

        self.a.Define('mttbar','hardware::InvariantMass({HadronicTop_vect, Bot_vect, WLep_vect, Neut_vect,PhiLep1_vect,PhiLep2_vect})')#invariant mass of the resonance particle
        return self.a.GetActiveNode()
    
    def Snapshot(self,node=None,colNames=[]):
        startNode = self.a.GetActiveNode()
        if node == None: node = self.a.GetActiveNode()
        #colNames[str]:give what variales to keep at the snapshot

        columns = [
            'mttbar'
        ]

        if (len(colNames) > 0):
            columns.extend(colNames)

        self.a.SetActiveNode(node)
        self.a.Snapshot(columns,'ttbarphisnapshot_%s_%s_%sof%s.root'%(self.setname,self.year,self.ijob,self.njobs),'Events',openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

        


