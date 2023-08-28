#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi, and some trig/hyperbolic functions.
#include <map>
#include <algorithm>
#include <Math/Vector4D.h>
#include <utility>

using namespace ROOT::VecOps;



RVec<int> PickDijetsV2(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Jet_btagCMVA, RVec<float> FatJet_particleNet_TvsQCD){
    int FatJetidx = -1;
    int Jetidx = -1;
    bool exitFatJetloop = false;

    for (int iFatJet =0; iFatJet < FatJet_phi.size() && exitFatJetloop == false;iFatJet++){//find the back to back Fatjet and jets
        for (int iJet = 0; iJet < Jet_phi.size();iJet++){
            if ((abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2) && Jet_btagCMVA[iJet] > 0.8 && FatJet_particleNet_TvsQCD[iFatJet] > 0.0){
                FatJetidx = iFatJet;
                Jetidx = iJet;
                exitFatJetloop = true;
                break;
            }
        }
    }

    return {FatJetidx,Jetidx};

}

//Code assisted by Xianglong Wang
//this method can fail if multiple lepton have same velocity. This is supposed to be extremely unlikely
RVec<int> FindLeadLepton(RVec<float> Electron_pt, RVec<float> Muon_pt){
    //(Leptonpt,(i,j)),i indicate its position inside the Rvector, j indicate which Rvector it belongs to
    std::map<float, std::pair<int, int>> LeptonIndex;
    for (int i = 0; i < Electron_pt.size(); i++){
        LeptonIndex[Electron_pt[i]] = std::make_pair(i,1);//(i,1) stand for element of electrons
    }

    for (int i = 0; i < Muon_pt.size(); i++){
        LeptonIndex[Muon_pt[i]] = std::make_pair(i,2);//(i,2) stand for element of muons
    }

    std::vector<float> Lepton_pt (Muon_pt.size());

    for (int i = 0; i < Muon_pt.size(); i++){
        Lepton_pt[i] = Muon_pt[i];
    }

    std::sort(Lepton_pt.begin(),Lepton_pt.end(),std::greater<float>());

    RVec<float> LeadLeptonInfo (2*Lepton_pt.size());

    for (int i = 0; i < Lepton_pt.size(); i++){
        LeadLeptonInfo[i] = LeptonIndex[Lepton_pt[i]].second;
        LeadLeptonInfo[i+Lepton_pt.size()] = LeptonIndex[Lepton_pt[i]].first;
    }


    return LeadLeptonInfo;

}

int MinPtConstraint(RVec<float> Electron_pt, RVec<float> Muon_pt, int LeptonType, int LeptonId){
    //based on whether this is an electron or a muon, compute whether its pt > 25 or not. Once this is satisfied then the other two automatically have pt > 25.
    int PassMinPt = 0;
    if (LeptonType == 1 && Electron_pt[LeptonId] > 25){
        PassMinPt = 1;
    }
    else{
        if (LeptonType == 2 && Muon_pt[LeptonId] > 25){
            PassMinPt = 1;
        }
    }
    return PassMinPt;
}


float LepbJetPtRel(float Lepton_pt, float Lepton_phi, float Lepton_eta, float bJet_pt, float bJet_phi, float bJet_eta){//compute the relative lepton momentum prependicular to bJet
    float Lepton_px = Lepton_pt * cos(Lepton_phi);
    float Lepton_py = Lepton_pt * sin(Lepton_phi);
    float Lepton_pz = Lepton_pt * sinh(Lepton_eta);
    float Lepton_p = Lepton_pt * cosh(Lepton_eta);
    float bJet_px = bJet_pt * cos(bJet_phi);
    float bJet_py = bJet_pt * sin(bJet_phi);
    float bJet_pz = bJet_pt * sinh(bJet_eta);
    float bJet_p = bJet_pt * cosh(bJet_eta);
    float Rel_px = Lepton_px - (bJet_px * (Lepton_px*bJet_px + Lepton_py*bJet_py + Lepton_pz*bJet_pz))/(bJet_p*bJet_p);
    float Rel_py = Lepton_py - (bJet_py * (Lepton_px*bJet_px + Lepton_py*bJet_py + Lepton_pz*bJet_pz))/(bJet_p*bJet_p);
    float Rel_pz = Lepton_pz - (bJet_pz * (Lepton_px*bJet_px + Lepton_py*bJet_py + Lepton_pz*bJet_pz))/(bJet_p*bJet_p);
    float Rel_psquare = Rel_px*Rel_px + Rel_py*Rel_py + Rel_pz*Rel_pz;
    return Rel_psquare;
}

int TwoDCutV2(int LeptonType, int LeptonId, RVec<float> Electron_pt, RVec<float> Muon_pt, float bJet_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, float bJet_phi, RVec<float> Electron_eta, RVec<float> Muon_eta, float bJet_eta){
    //this is the program to do 2DCut in the following LeptonCategorize function. It's written as a separate function to keep the readability of code.
    int Crel_pt = -1;//relative pt criteria
    int Crel_phi = -1;//relative phi criteria
    int iElectronId = 0;
    int iMuonId = 0;
    int Pass2DCut = 0;

    if (LeptonType == 1){
        iElectronId = LeptonId;
        if(LepbJetPtRel(Electron_pt[iElectronId], Electron_phi[iElectronId], Electron_eta[iElectronId], bJet_pt, bJet_phi, bJet_eta) > (25*25)){
            Crel_pt = 1;
        }
        if(abs(hardware::DeltaPhi(Electron_phi[iElectronId],bJet_phi)) > 0.4){
            Crel_phi = 1;
        }
    }
    else{
        if(LeptonType == 2){
            iMuonId = LeptonId;
            if(LepbJetPtRel(Muon_pt[iMuonId], Muon_phi[iMuonId], Muon_eta[iMuonId], bJet_pt, bJet_phi, bJet_eta) > (25*25)){
                Crel_pt = 1;
            }
            if(abs(hardware::DeltaPhi(Muon_phi[iMuonId],bJet_phi)) > 0.4){
                Crel_phi = 1;
            }
        }
    }

    if (Crel_pt == 1 || Crel_phi == 1){
        Pass2DCut = 1;
    }

    return Pass2DCut;
}


std::vector<std::vector<int>> MakePairs (std::vector<int> Indices){//warning:make sure the input has at least 2 elements!

    std::vector<std::vector<int>> Pairs;

    for (int i = 0; i < Indices.size(); i++){
        for (int j = i + 1; j < Indices.size(); j++){
            Pairs.push_back ({Indices[i],Indices[j]});
        }
    }

    return Pairs;
}

std::vector<std::vector<int>> FindSameFlavorPairs (std::vector<int> LeptonType){
//this code will find out a collection of electron and muon, and run combinations over the collection if they have at least two element. we expect size of LeptonType >=3 due to cut in preselection
    std::vector<int> ElectronIndices;
    std::vector<int> MuonIndices;
    std::vector<std::vector<int>> ElectronPairs;
    std::vector<std::vector<int>> MuonPairs;
    std::vector<std::vector<int>> LeptonPairs;
    for (int i = 0; i < LeptonType.size(); i++){
        if (LeptonType[i] == 1){//this is electron
            ElectronIndices.push_back (i);
        }
        else{
            MuonIndices.push_back (i);
        }
    }

    if (ElectronIndices.size() < 2){
        MuonPairs = MakePairs (MuonIndices);
        LeptonPairs.insert(LeptonPairs.end(),MuonPairs.begin(),MuonPairs.end());
    }
    else{
        if (MuonIndices.size() < 2){
            ElectronPairs = MakePairs (ElectronIndices);
            LeptonPairs.insert(LeptonPairs.end(),ElectronPairs.begin(),ElectronPairs.end());
        }
        else{
            MuonPairs = MakePairs (MuonIndices);
            ElectronPairs = MakePairs (ElectronIndices);
            LeptonPairs.insert(LeptonPairs.end(),MuonPairs.begin(),MuonPairs.end());
            LeptonPairs.insert(LeptonPairs.end(),ElectronPairs.begin(),ElectronPairs.end());
        }
    }

    return LeptonPairs;
}

int GetIntLeptonProperty(int LeptonType,int LeptonId,RVec<int> ElectronProperty, RVec<int> MuonProperty){

    int LeptonFloat = -1;//by setting default value to -1 will give us a warning messeage if no lepton satisfying condition exist. Such event should be cut in preselections
    
    if(LeptonType == 1){//this is an electron
        LeptonFloat = ElectronProperty[LeptonId];//then the lepton's property used to construct 4 vector should be relavent electron property
    }
    else{
        if(LeptonType == 2){//this is a muon
            LeptonFloat = MuonProperty[LeptonId];
        }
    }

    return LeptonFloat;
}

float GetFloatLeptonProperty(int LeptonType,int LeptonId,RVec<float> ElectronProperty, RVec<float> MuonProperty){

    float LeptonFloat = -1.0;//by setting default value to -1 will give us a warning messeage if no lepton satisfying condition exist. Such event should be cut in preselections
    
    if(LeptonType == 1){//this is an electron
        LeptonFloat = ElectronProperty[LeptonId];//then the lepton's property used to construct 4 vector should be relavent electron property
    }
    else{
        if(LeptonType == 2){//this is a muon
            LeptonFloat = MuonProperty[LeptonId];
        }
    }

    return LeptonFloat;
}

std::vector<std::vector<int>> FindOppositeSignPairs (RVec<int> LeadLeptonInfo, std::vector<std::vector<int>> Pairs, RVec<int> Electron_charge, RVec<int> Muon_charge){
    //this will slim the Pairs down to pairs with opposite signs.
    std::vector<std::vector<int>> OppositeSignPairs;
    int Lepton1Charge;
    int Lepton2Charge;
    for (int i = 0; i < Pairs.size(); i++){
        Lepton1Charge = GetIntLeptonProperty(LeadLeptonInfo[Pairs[i][0]],LeadLeptonInfo[Pairs[i][0] + (LeadLeptonInfo.size()/2)],Electron_charge,Muon_charge);
        Lepton2Charge = GetIntLeptonProperty(LeadLeptonInfo[Pairs[i][1]],LeadLeptonInfo[Pairs[i][1] + (LeadLeptonInfo.size()/2)],Electron_charge,Muon_charge);
        if (Lepton1Charge + Lepton2Charge == 0){
            OppositeSignPairs.push_back (Pairs[i]);
        }
    }
    return OppositeSignPairs;
}

float FindPhiMass (std::vector<int> Pair, RVec<int> LeadLeptonInfo, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, RVec<float> Electron_eta, RVec<float> Muon_eta){

    float PhiLepton1Pt = GetFloatLeptonProperty (LeadLeptonInfo[Pair[0]],LeadLeptonInfo[Pair[0] + (LeadLeptonInfo.size()/2)],Electron_pt,Muon_pt);
    float PhiLepton2Pt = GetFloatLeptonProperty (LeadLeptonInfo[Pair[1]],LeadLeptonInfo[Pair[1] + (LeadLeptonInfo.size()/2)],Electron_pt,Muon_pt);
    float PhiLepton1Eta = GetFloatLeptonProperty (LeadLeptonInfo[Pair[0]],LeadLeptonInfo[Pair[0] + (LeadLeptonInfo.size()/2)],Electron_eta,Muon_eta);
    float PhiLepton2Eta = GetFloatLeptonProperty (LeadLeptonInfo[Pair[1]],LeadLeptonInfo[Pair[1] + (LeadLeptonInfo.size()/2)],Electron_eta,Muon_eta);
    float PhiLepton1Phi = GetFloatLeptonProperty (LeadLeptonInfo[Pair[0]],LeadLeptonInfo[Pair[0] + (LeadLeptonInfo.size()/2)],Electron_phi,Muon_phi);
    float PhiLepton2Phi = GetFloatLeptonProperty (LeadLeptonInfo[Pair[1]],LeadLeptonInfo[Pair[1] + (LeadLeptonInfo.size()/2)],Electron_phi,Muon_phi);
    float PhiLepton1M;
    float PhiLepton2M;
    if (LeadLeptonInfo[Pair[0]] == 1){
        PhiLepton1M = 0.000511;
        PhiLepton2M = 0.000511;
    }
    else{
        PhiLepton1M = 0.105658;
        PhiLepton2M = 0.105658;
    }
    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(PhiLepton1Pt,PhiLepton1Eta,PhiLepton1Phi,PhiLepton1M);
    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(PhiLepton2Pt,PhiLepton2Eta,PhiLepton2Phi,PhiLepton2M);

    float InvariantMass = hardware::InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector});
    return InvariantMass;
}

std::vector<int> FindLowerPhiMass (std::vector<std::vector<int>> Pairs, RVec<int> LeadLeptonInfo, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, RVec<float> Electron_eta, RVec<float> Muon_eta){
    std::vector<float> PhiMass;
    int LowestMassPair;
    float LowestMassValue;

    for (int i = 0; i < Pairs.size();i++){
        PhiMass.push_back(FindPhiMass(Pairs[i],LeadLeptonInfo, Electron_pt, Muon_pt, Electron_phi, Muon_phi, Electron_eta, Muon_eta));
    }

    LowestMassPair = 0;
    LowestMassValue = PhiMass[0];

    for (int i = 0; i < Pairs.size(); i++){
        if (PhiMass[i] < LowestMassValue){
            LowestMassValue = PhiMass[i];
            LowestMassPair = i;
        }
    }

    return Pairs[LowestMassPair];
}

//Don't bother with the lepton from top in non-boosted case. Find all the relavent pairs and just compute if 1. they have opposite charge 2. their invariant mass.
RVec<int> FindPhiLepton (RVec<int> LeadLeptonInfo, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, RVec<float> Electron_eta, RVec<float> Muon_eta, RVec<int> Electron_charge, RVec<int> Muon_charge){
    //first, find out all the pairs using LeadLeptonInfo: the first half of it will contain only 1 or 2 based on whether it's electron or muon
    //it's possible that an event will contain no qualified pairs because they are all of the same sign. We want to abort such data if happened.
    std::vector<int> LeptonType;
    std::vector<std::vector<int>> SameFlavorPairs;
    std::vector<std::vector<int>> SameFlavorOppositeSignPairs;
    std::vector<int> LowerPhiMassPair = {0,0};
    RVec<int> PhiLeptonPair = {0,0,0}; //The first digit is used to make sure all test are passed.
    for (int i = 0; i < (LeadLeptonInfo.size()/2); i++){
        LeptonType.push_back(LeadLeptonInfo[i]);
    }

    SameFlavorPairs = FindSameFlavorPairs (LeptonType);
    SameFlavorOppositeSignPairs = FindOppositeSignPairs (LeadLeptonInfo,SameFlavorPairs,Electron_charge,Muon_charge);
    //now, we compute the mass of all these pairs, and choose the pair with lowest invariant mass
    if (SameFlavorOppositeSignPairs.size() > 0){
        PhiLeptonPair[0] = 1;
        LowerPhiMassPair = FindLowerPhiMass (SameFlavorOppositeSignPairs, LeadLeptonInfo, Electron_pt, Muon_pt, Electron_phi, Muon_phi,Electron_eta, Muon_eta);
    }

    PhiLeptonPair[1] = LowerPhiMassPair[0];
    PhiLeptonPair[2] = LowerPhiMassPair[1];

    return PhiLeptonPair;
}

float NeutrinoEta(float Lepton_pt, float Lepton_phi, float Lepton_eta, float MET_pt, float MET_phi){//find eta using mass of Wboson
    float W_mass = 80.4;
    float Lepton_px = Lepton_pt * cos(Lepton_phi);
    float Lepton_py = Lepton_pt * sin(Lepton_phi);
    float Lepton_pz = Lepton_pt * sinh(Lepton_eta);
    float Lepton_p = Lepton_pt * cosh(Lepton_eta);
    float MET_px = MET_pt * cos(MET_phi);
    float MET_py = MET_pt * sin(MET_phi);
    float Lepton_Esquare = Lepton_p * Lepton_p;// leptons has rest energy on order of MeV, for GeV events lets' pretend they are 0
    float Lambda = (W_mass * W_mass)/2 + Lepton_px*MET_px + Lepton_py*MET_py;
    float Neutrino_pz = ((Lambda * Lepton_pz)/(Lepton_pt * Lepton_pt)) - sqrt(pow((Lambda * Lepton_pz)/(Lepton_pt * Lepton_pt),2) - (Lepton_Esquare * MET_pt * MET_pt - Lambda * Lambda)/(Lepton_pt * Lepton_pt));
    float Neutrino_eta = atanh(Neutrino_pz/sqrt(MET_pt * MET_pt + Neutrino_pz * Neutrino_pz));
    return Neutrino_eta;
}

float LeptonicCandidatePt(float Bot_pt, float Lepton_pt, float Bot_phi, float Lepton_phi, float Neutrino_pt, float Neutrino_phi){
    float px = 0;
    float py = 0;
    float pt_tot = 0;
    px = Bot_pt * cos(Bot_phi) + Lepton_pt * cos(Lepton_phi) + Neutrino_pt * cos(Neutrino_phi);
    py = Bot_pt * sin(Bot_phi) + Lepton_pt * sin(Lepton_phi) + Neutrino_pt * sin(Neutrino_phi);
    pt_tot = sqrt(px*px + py*py);
    return pt_tot;
}

float TotalPt(float Top_pt,float Top_phi,float Bot_pt,float Bot_phi,float Lepton_pt,float Lepton_phi,float Neutrino_pt,float Neutrino_phi){
    float px = 0;
    float py = 0;
    float pt_tot = 0;
    px = Top_pt * cos(Top_phi) + Bot_pt * cos(Bot_phi) + Lepton_pt * cos(Lepton_phi) + Neutrino_pt * cos(Neutrino_phi);
    py = Top_pt * sin(Top_phi) + Bot_pt * sin(Bot_phi) + Lepton_pt * sin(Lepton_phi) + Neutrino_pt * sin(Neutrino_phi);
    pt_tot = sqrt(px*px + py*py);
    return pt_tot;
}

const ROOT::RVec<int> FindMothersPdgId(const ROOT::RVec<int>& genpart_id, const ROOT::RVec<int>& selected_genpart_mother_indices){

    std::size_t Ni = selected_genpart_mother_indices.size();
    RVec<int> mother_pdgids(Ni);    
    for(std::size_t i=0; i<Ni; i++) {
        mother_pdgids[i] = genpart_id[selected_genpart_mother_indices[i]];
    }
    return mother_pdgids;

}