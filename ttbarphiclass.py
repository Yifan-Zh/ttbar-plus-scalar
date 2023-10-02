import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU,AutoPU
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
                print('Adding {} to Analyzer'.format(iFile))
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
            self.setname = inputfile.split('/')[-1].split('_')[0]
        self.year = str(year)	# most of the time this class will be instantiated from other scripts with CLI args, so just convert to string for internal use
        self.ijob = ijob
        self.njobs = njobs
        print ('setname is {}, year is {}'.format(self.setname,self.year))
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
        
        # the following is added due to a current bug in corrections on signal files
        if 'signal' in inputfile:
            self.a.isSignal = True
        else:
            self.a.isSignal = False

    #this is the end of initializer. Now we can apply preselection. Cpp modules will be compiled for each task, but we will not compile them in class function.

    def AddCutflowColumn(self, var, varName):
        print('Adding cutflow information....\n\t{}\t{}'.format(varName, var))
        self.a.Define('{}'.format(varName), str(var))

    def getNweighted(self):
        if not self.a.isData:
            return self.a.DataFrame.Sum("genWeight").GetValue()
        else:
            return self.a.DataFrame.Count().GetValue()
        
    #standard correction to MC/Data:

    def ApplyStandardCorrections(self,snapshot=False):
        if snapshot:
            if self.a.isData:
                lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16])
                self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"}))

            else:
                #self.a = AutoPU(self.a, self.year, ULflag=True)
                self.a.AddCorrection(
                    Correction('Pdfweight','TIMBER/Framework/include/PDFweight_uncert.h',[self.a.lhaid],corrtype='uncert')
                )
                if self.year == '16' or self.year == '17' or self.year == '16APV':
		        # The column names are: L1PreFiringWeight_Up, L1PreFiringWeight_Dn, L1PreFiringWeight_Nom
                    L1PreFiringWeight = Correction("L1PreFiringWeight","TIMBER/Framework/TopPhi_modules/BranchCorrection.cc",constructor=[],mainFunc='evalWeight',corrtype='weight',columnList=['L1PreFiringWeight_Nom','L1PreFiringWeight_Up','L1PreFiringWeight_Dn'])
                    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})	
                '''
                elif self.year == '18':
                    self.a.AddCorrection(
                        Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr')
                    )
                '''
        
        else:
            if not self.a.isData:
                self.a.AddCorrection(Correction('Pileup',corrtype='weight'))
                self.a.AddCorrection(Correction('Pdfweight',corrtype='uncert'))
                if self.year == '16' or self.year == '17' or self.year == '16APV':
                    #self.a.AddCorrection(Correction('Prefire',corrtype='weight'))
                    self.a.AddCorrection(Correction('L1PreFiringWeight',corrtype='weight'))
                #elif self.year == '18':
                    #self.a.AddCorrection(Correction('HEM_drop',corrtype='corr'))
                if 'ttbar' in self.setname:
                    self.a.AddCorrection(Correction('TptReweight',corrtype='weight'))
                
        return self.a.GetActiveNode()


    # this is the preselection for non-boosted case. In this case DO NOT use top tagger (it won't work)
    # we first want to mark all the event with:
    def Preselection(self):
        self.NPROC = self.getNweighted()
        self.AddCutflowColumn(self.NPROC,"NPROC")
        #we need either:at least 1 AK8 + 1 AK4, or at least 3 AK4 jet For a semileptonic decay.
        #IMPORTANT:in fact, it might be better to demand 4 or more AK4 jets; it's fairly unlikely the two jets from W decay will merge
        
        #we want the jets the have at least pt >35GeV, FatJet with pt > 70GeV
        self.a.SubCollection('EnergeticFatJet','FatJet','FatJet_pt > 70')
        self.a.SubCollection('EnergeticJet','Jet','Jet_pt > 35')
        self.a.Cut('nFatJet','(nEnergeticFatJet > 0 && nEnergeticJet > 0) || (nEnergeticJet > 2)')
        #in any case, we should be able to pickup at least one b-jet using b tagger. Don't do it, cause loss in signal efficiency. Some b jets will not be picked up
        #self.a.Define('nbjets','Findnbjets(EnergeticJet_btagCSVV2,0.46)')
        #self.a.Cut('nbjetCut','nbjets > 0')
        self.NJETS = self.getNweighted()
        self.AddCutflowColumn(self.NJETS,"NJETS")

        # we do not want electrons produced by photon pairs. They will mess up our data because they have 0 invariant mass.
        self.a.SubCollection('NonConvertedElectron','Electron','Electron_convVeto == 1')
        self.a.Cut('nLepton','nNonConvertedElectron > 0 || nMuon > 0') #make sure at least one lepton exist
        #self.a.SubCollection('NotHvyMuon','Muon','Muon_genPartFlav == 0 || Muon_genPartFlav == 1')#exclude the muons coming from b hadron decay as they will also have very low mass. This is for MC data only. For actually data, similar effect can be achieved by doing isolation.
        self.a.SubCollection('NotHvyMuon','Muon','Muon_pt > 0.01 && Muon_pfRelIso04_all <0.15')#tight isolation

        #we do not want to reconstruct ttbar in this case. It's very difficult to do without the boosted condition


        self.a.Define('nTotalLepton','nNonConvertedElectron + nNotHvyMuon')
        self.a.Cut('LeptonNumberCut','nTotalLepton > 2')
        self.NLeptons = self.getNweighted()
        self.AddCutflowColumn(self.NLeptons,"NLeptons")         

        # we would find the leading leptons according to their pt. The output has the form {Electron/Muon(represented by 1/2), relative postion insde the corresponding vector}
        # for example, if the leading leptons are Electron[2],Muon[3],Electroon[4], then it would be {1,2,1,2,3,4}
        
        #note:currently, modified to examine Muons only. This means that we should have at least 2 Muon per event.
	    #we want to focus on the Muon NOT from heavy flavor particles
        self.a.Cut('MuonNumberCut','nNotHvyMuon > 1')
        self.a.Define('LeadingThreeLepton','FindLeadLepton(NonConvertedElectron_pt,NotHvyMuon_pt)')#not the leading three anymore, it's actually just all muons
        self.a.Define('nLeadingLeptons','LeadingThreeLepton.size()/2')

        #self.NPreselection = self.getNweighted()
        #self.AddCutflowColumn(self.NPreselection,"NPreselection")
        print ("Pass Preselection stage")
        return self.a.GetActiveNode()
    

    def Selection(self):

        #now we start to handle the leptons. We'll handle this part in c++
        self.a.Define('LeptonTestAndReOrdering','FindPhiLepton(LeadingThreeLepton,NonConvertedElectron_pt,NotHvyMuon_pt,NonConvertedElectron_phi,NotHvyMuon_phi,NonConvertedElectron_eta,NotHvyMuon_eta,NonConvertedElectron_charge,NotHvyMuon_charge)')
        self.a.Cut('PassAllSelection','LeptonTestAndReOrdering[0] == 1')
        #self.NPassAllSelection = self.getNweighted()
        #self.AddCutflowColumn(self.NPassAllSelection,"NPassAllSelection")  
        print ("Pass selection stage")
        return self.a.GetActiveNode()
    
    #now we need to make the plot. For purpose of invariant mass reconstruction, we need to specify the lepton pt, eta, phi and mass manually. We will do this using a user defined C++ code.
    def JetsCandidateKinematicinfo(self):
        #first give relatvent information of lepton; do not use Lepton_*, will cause a bug in snapshot
        #we now define the kinematics variables of the three lepton. The one from W followed by ones from phi

        self.a.Define('PhiLepton1_pt','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons],NonConvertedElectron_pt,NotHvyMuon_pt)')
        self.a.Define('PhiLepton1_eta','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons],NonConvertedElectron_eta,NotHvyMuon_eta)')
        self.a.Define('PhiLepton1_phi','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons],NonConvertedElectron_phi,NotHvyMuon_phi)')
        self.a.Define('PhiLepton1_mass','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons] ,NonConvertedElectron_mass,NotHvyMuon_mass)')
        self.a.Define('PhiLepton1_PfIso04','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons],NonConvertedElectron_pfRelIso03_all,NotHvyMuon_pfRelIso04_all)')

        self.a.Define('PhiLepton2_pt','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + nLeadingLeptons],NonConvertedElectron_pt,NotHvyMuon_pt)')
        self.a.Define('PhiLepton2_eta','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + nLeadingLeptons],NonConvertedElectron_eta,NotHvyMuon_eta)')
        self.a.Define('PhiLepton2_phi','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + nLeadingLeptons],NonConvertedElectron_phi,NotHvyMuon_phi)')
        self.a.Define('PhiLepton2_mass','GetFloatLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[2]],LeadingThreeLepton[LeptonTestAndReOrdering[2] + nLeadingLeptons] ,NonConvertedElectron_mass,NotHvyMuon_mass)')

        self.a.Define('PhiLepton1_MotherType','GetIntLeptonProperty(LeadingThreeLepton[LeptonTestAndReOrdering[1]],LeadingThreeLepton[LeptonTestAndReOrdering[1] + nLeadingLeptons],NonConvertedElectron_genPartFlav,NotHvyMuon_genPartFlav)')


        return self.a.GetActiveNode()
    
    
    #We are ready to make plots. current plots: leptonic candidate pt, hardonic candidate pt, inv mass of ttbar
    #hardronic candidate: easy. Just AK8 Jet mass. No need to do any combination whatsover
    #leptonic candidate: need information from btaggedAK4, Lepton, and MET. Assume 0 mass nutrino, moves in xy plane only so eta =0.

    def MassReconstruction(self):

        self.a.Define('PhiLep1_vect','hardware::TLvector(PhiLepton1_pt,PhiLepton1_eta,PhiLepton1_phi,PhiLepton1_mass)')
        self.a.Define('PhiLep2_vect','hardware::TLvector(PhiLepton2_pt,PhiLepton2_eta,PhiLepton2_phi,PhiLepton2_mass)')
 
        self.a.Define('PhiInvMass','hardware::InvariantMass({PhiLep1_vect,PhiLep2_vect})')#invariant mass of the resonance particle
        self.a.Define('WhichLepton','LeadingThreeLepton[LeptonTestAndReOrdering[2]]')
        self.a.Define('PhiDeltaR','hardware::DeltaR(PhiLep1_vect,PhiLep2_vect)')
        #self.a.Define('PhiPt','PhiCandidatePt(PhiLepton1_pt,PhiLepton2_pt,PhiLepton1_phi,PhiLepton2_phi)')
        #self.a.Cut('WhateverDebugThisIs','PhiLepton1_pt > 30 && PhiLepton2_pt > 30')
        self.NFinalEvent = self.getNweighted()
        self.AddCutflowColumn(self.NFinalEvent,"NFinalEvent")
        return self.a.GetActiveNode()
    
    def Snapshot(self,node=None,colNames=[]):
        startNode = self.a.GetActiveNode()
        if node == None: node = self.a.GetActiveNode()
        #colNames[str]:give what variales to keep at the snapshot

        columns = [
            'PhiInvMass','WhichLepton','PhiDeltaR',
            'PhiLepton1_pt','PhiLepton1_eta','PhiLepton1_phi','PhiLepton1_mass','PhiLepton1_MotherType','PhiLepton1_PfIso04',
            'PhiLepton2_pt','PhiLepton2_eta','PhiLepton2_phi','PhiLepton2_mass'
        ]
        
       # columns = ['nMuon']

        
        if not self.a.isData:
            if self.a.isSignal == False:
                columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down'])
                columns.extend(['weight__Pileup_up','weight__Pileup_down','weight__nominal','weight__Pdfweight_down','weight__Pdfweight_up'])

                if self.year == '16' or self.year == '17' or self.year == '16APV':
                    columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])
                    columns.extend(['weight__L1PreFiringWeight_down','weight__L1PreFiringWeight_up'])
	
                '''
                elif self.year == '18':
                    columns.append('HEM_drop__nom')
                '''
            elif self.a.isSignal == True:
                columns.extend(['Pileup__nom','Pileup__up','Pileup__down'])
                columns.extend(['weight__Pileup_up','weight__Pileup_down','weight__nominal'])

        if (len(colNames) > 0):
            columns.extend(colNames)

        self.a.SetActiveNode(node)
        self.a.Snapshot(columns,'ttbarphisnapshot_%s_%s_%sof%s.root'%(self.setname,self.year,self.ijob,self.njobs),'Events',openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year if 'APV' not in self.year else 16)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        return lumi*xsec/self.a.genEventSumw

