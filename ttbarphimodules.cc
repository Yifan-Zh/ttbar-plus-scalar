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
//we want to 1. examine if there's at least a 50GeV lepton 2. If so, find the back to back AK8 and AK4 jets (we know there would at most be one)
//Problem left to fix: modify the code to pick to most energetic lepton


RVec<int> PickDijetsV2(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Jet_btagCMVA){
    int FatJetidx = -1;
    int Jetidx = -1;
    bool exitLeptonloop = false;
    bool exitJetloop = false;

    for (int iJet = 0; iJet < Jet_phi.size() && exitJetloop ==false; iJet++){//find the back to back jets
        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
            if (abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2 && Jet_btagCMVA[iJet] > 0.8){
                FatJetidx = iFatJet;
                Jetidx = iJet;
                exitJetloop =true;
                break;
            }
        }
    }

    return {FatJetidx,Jetidx};

}

//Code assisted by Xianglong Wang
//this method can fail if multiple lepton have same velocity. This is supposed to be extremely unlikely
RVec<int> FindLeadLepton(RVec<float> Electron_pt, RVec<float> Muon_pt){
    //(Leptonpt,(i,j)),i indicate its position inside the vector, j indicate which vector it belongs to
    std::map<float, std::pair<int, int>> LeptonIndex;
    for (int i = 0; i < Electron_pt.size(); i++){
        LeptonIndex[Electron_pt[i]] = std::make_pair(i,1);//(i,1) stand for element of electrons
    }

    for (int i = 0; i < Muon_pt.size(); i++){
        LeptonIndex[Muon_pt[i]] = std::make_pair(i,2);//(i,2) stand for element of muons
    }

    RVec<float> Lepton_pt (Electron_pt.size() + Muon_pt.size());

    for (int i = 0; i < Electron_pt.size(); i++){
        Lepton_pt[i] = Electron_pt[i];
    }

    for (int i = 0; i < Muon_pt.size(); i++){
        Lepton_pt[i + Electron_pt.size()] = Muon_pt[i];
    }

    std::sort(Lepton_pt.begin(),Lepton_pt.end(),std::greater<float>());
    Lepton_pt.resize(3);
    return {LeptonIndex[Lepton_pt[0]].second,LeptonIndex[Lepton_pt[1]].second,LeptonIndex[Lepton_pt[2]].second,LeptonIndex[Lepton_pt[0]].first,LeptonIndex[Lepton_pt[1]].first,LeptonIndex[Lepton_pt[2]].first};

}

int MinPtConstraint(RVec<float> Electron_pt, RVec<float> Muon_pt, int LeptonType, int LeptonId){
    //based on whether this is an electron or a muon, compute whether its pt >50 or not. Once this is satisfied then the other two automatically have pt > 50.
    int PassMinPt = 0;
    if (LeptonType == 1 && Electron_pt[LeptonId] > 50){
        PassMinPt = 1;
    }
    else{
        if (LeptonType == 2 && Muon_pt[LeptonId] > 50){
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

RVec<int> RemainingLepton (int TopLepton){
    //this code takes in the id of identified lepton from top, and returns the remaining lepton. For example, if the 2nd energetic lepton is from top, then it will return {0,2}
    RVec<int> Remain = {-1,-1};

    if (TopLepton == 0){
        Remain = {1,2};
    }
    if (TopLepton == 1){
        Remain = {0,2};
    }
    if (TopLepton == 2){
        Remain = {0,1};
    }
    return Remain;
}

RVec<int> FindParticlePairs(int LeptonType, int Lepton0Id, int Lepton1Id, int Lepton2Id, RVec<int> Electron_charge, RVec<int> Muon_charge){
    //there will be 2 different combination. We shall find them out first. 1 = {0,1}, 2 = {0,2}, 3= {1,2}. 
    //If we get a vector {1,2} this means pair {0,1} and pair {0,2} are pair with different charge and therefore particle-antiparticle pairs
    RVec<int> Pairs = {0,0};
    if (LeptonType == 1){
        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton1Id] == 0) && (Electron_charge[Lepton0Id] + Electron_charge[Lepton2Id] ==0)){
            Pairs = {1,2};
        }
        
        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton1Id] == 0) && (Electron_charge[Lepton1Id] + Electron_charge[Lepton2Id] ==0)){
            Pairs = {1,3};
        }

        if ((Electron_charge[Lepton0Id] + Electron_charge[Lepton2Id] == 0) && (Electron_charge[Lepton1Id] + Electron_charge[Lepton2Id] ==0)){
            Pairs = {2,3};
        }
    }

    if (LeptonType == 2){
        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton1Id] == 0) && (Muon_charge[Lepton0Id] + Muon_charge[Lepton2Id] ==0)){
            Pairs = {1,2};
        }
        
        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton1Id] == 0) && (Muon_charge[Lepton1Id] + Muon_charge[Lepton2Id] ==0)){
            Pairs = {1,3};
        }

        if ((Muon_charge[Lepton0Id] + Muon_charge[Lepton2Id] == 0) && (Muon_charge[Lepton1Id] + Muon_charge[Lepton2Id] ==0)){
            Pairs = {2,3};
        }
    }

    return Pairs;
}

//this code find out the top lepton based on the particle-antiparitle pair. 
int LeptonFromTop(int Pair){
    int TopLepton = -1;
    if (Pair == 1){
        TopLepton = 2;
    }
    if (Pair == 2){
        TopLepton = 1;
    }
    if (Pair == 3){
        TopLepton = 0;
    }
    return TopLepton;
}

float FindDeltaPhi (int PairIndex, float Lepton1Phi, float Lepton2Phi, float Lepton3Phi){
    float PairDeltaPhi;
    if (PairIndex == 1){
        PairDeltaPhi = abs(hardware::DeltaPhi(Lepton1Phi,Lepton2Phi));
    }
    
    if (PairIndex == 2){
        PairDeltaPhi = abs(hardware::DeltaPhi(Lepton1Phi,Lepton3Phi));
    }

    if (PairIndex ==3){
        PairDeltaPhi = abs(hardware::DeltaPhi(Lepton2Phi,Lepton3Phi));
    }
    
    return PairDeltaPhi;
}

int FindLowerDeltaPhi(RVec<int> Pairs, float Lepton1Phi, float Lepton2Phi, float Lepton3Phi){
    //find the deltaphi for each pair separately.
    int Pair = -1;
    float DeltaPhi1 = FindDeltaPhi(Pairs[0],Lepton1Phi,Lepton2Phi,Lepton3Phi);
    float DeltaPhi2 = FindDeltaPhi(Pairs[1],Lepton1Phi,Lepton2Phi,Lepton3Phi);
    if (DeltaPhi1 < DeltaPhi2){
        Pair = Pairs[0];
    }
    else{
        Pair = Pairs[1];
    }

    return Pair;
}

RVec<int> LeptonCategorize(RVec<int> LeadLeptonInfo, RVec<float> Electron_pt, RVec<float> Muon_pt, float bJet_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, float bJet_phi, RVec<float> Electron_eta, RVec<float> Muon_eta, float bJet_eta, RVec<int> Electron_charge, RVec<int> Muon_charge){
    //this code first find out whether we have 3 lepton of same class or not, and then apply the relative selecton standard
    //initialize all the needed check conditions. By default, all of then are set to not passed.
    int PassAllTest = 0;
    int TopLepton = -1;//initialize the lepton from W
    int Pass2DCut = 0;
    int PassPairCheck = 0;
    int PassScalarMassCut = 0;
    int PassRelativePhi = 0;
    int PassJetRelativePhi = 0;
    float PhiLepton1Phi;
    float PhiLepton2Phi;
    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(0,0,0,0);
    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(0,0,0,0);
    bool ThreeSameLepton;
    RVec<int> LeptonPairId = {-1,-1};
    int Lepton1Charge;
    int Lepton2Charge;
    int Lepton3Charge;
    float Lepton1Phi;
    float Lepton2Phi;
    float Lepton3Phi;
    int LowerDeltaPhiPair;

    if (LeadLeptonInfo[0] == LeadLeptonInfo [1] && LeadLeptonInfo[1] == LeadLeptonInfo[2]){//three same lepton
        bool ThreeSameLepton = true;
    }
    else{
        bool ThreeSameLepton = false;
    }
    //let's work out the different case first.
    if (ThreeSameLepton == false){
        if (LeadLeptonInfo[0] == LeadLeptonInfo[1]) {
            // First and second elements are equal
            if (LeadLeptonInfo[0] != LeadLeptonInfo[2]) {
                // Third element is different
                TopLepton = 2; // Index 2 represents the third element
            }
        }
        else {
                // First and second elements are different
            if (LeadLeptonInfo[0] == LeadLeptonInfo[2]) {
                // Second element is different
                TopLepton = 1; // Index 1 represents the second element
            }
            else {
                // First element is different
                TopLepton = 0; // Index 0 represents the first element
            }
        }

        //the type is given by LeadLeptonInfo[0~2], and their id is given by LeadLeptonInfo[3~5]
        if (TwoDCutV2(LeadLeptonInfo[TopLepton],LeadLeptonInfo[TopLepton + 3],Electron_pt,Muon_pt,bJet_pt,Electron_phi,Muon_phi,bJet_phi,Electron_eta,Muon_eta,bJet_eta) == 1){
            Pass2DCut = 1;
        }
        else{
            Pass2DCut = 0;
        }
        //now, check if the remaining two can give a invariant mass < 20 GeV
        LeptonPairId = RemainingLepton(TopLepton);
        //extrac the type of lepton and their position inside either Electron or Muon set
        int PhiLeptonType1 = LeadLeptonInfo[LeptonPairId[0]];
        int PhiLeptonId1 = LeadLeptonInfo[LeptonPairId[0] + 3];
        int PhiLeptonType2 = LeadLeptonInfo[LeptonPairId[1]];
        int PhiLeptonId2 = LeadLeptonInfo[LeptonPairId[1] + 3];
        //check if they actually form a pair of particle-antiparticle by checking their charge.
        if (PhiLeptonType1 == 1){
            if (Electron_charge[PhiLeptonId1] + Electron_charge[PhiLeptonId2] == 0){
                PassPairCheck = 1;
            }
        }
        else{
            if (PhiLeptonType1 == 2){
                if(Muon_charge[PhiLeptonId1] + Muon_charge[PhiLeptonId2] == 0){
                    PassPairCheck = 1;
                }
            }
        }
    
        //construct the TLvector for these two leptons and then compute their invariant mass.
        if (PhiLeptonType1 == 1){
            ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
        }
        else{
            ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
        }
        if (PhiLeptonType2 == 1){
            ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
        }
        else{
            ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
        }

        if (hardware::InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
            PassScalarMassCut = 1;
        }


        //check if they are close to either of the jets and close to each other
        //first, collect their phi info
        if (PhiLeptonType1 == 1){
            PhiLepton1Phi = Electron_phi[PhiLeptonId1];
        }
        else{
            PhiLepton1Phi = Muon_phi[PhiLeptonId1];
        }

        if (PhiLeptonType2 == 1){
            PhiLepton2Phi = Electron_phi[PhiLeptonId2];
        }
        else{
            PhiLepton2Phi = Muon_phi[PhiLeptonId2];
        }
        //check their relative phi
        if (abs(hardware::DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
            PassRelativePhi = 1;
        }
        else{
            PassRelativePhi = 0;
        }

        if (PassPairCheck == 1 && Pass2DCut == 1 && PassScalarMassCut == 1 && PassRelativePhi == 1){
            //Pass all the test
            PassAllTest = 1;
        }
    }
    //now we handle the case where we have eee or mumumu
    if (ThreeSameLepton == true){
        //first, get rid of the case e-e-e- or mu-,mu-,mu-, cannot form particle-antiparticle pairs in this case
        //we will initialize the conditon outside of if. We do not want to run test over already failed pairs.
        int PassAntiparticleCheck = 0;
        RVec<int> Pairs = {0,0};



        if (LeadLeptonInfo[0] == 1){//all electrons
            Lepton1Charge = Electron_charge[LeadLeptonInfo[3]];
            Lepton2Charge = Electron_charge[LeadLeptonInfo[4]];
            Lepton3Charge = Electron_charge[LeadLeptonInfo[5]];
            Lepton1Phi = Electron_phi[LeadLeptonInfo[3]];
            Lepton2Phi = Electron_phi[LeadLeptonInfo[4]];
            Lepton3Phi = Electron_phi[LeadLeptonInfo[5]];
            if (Lepton1Charge == Lepton2Charge && Lepton2Charge == Lepton3Charge){
                PassAntiparticleCheck = 0;
            }
            else{
                PassAntiparticleCheck = 1;
                //there will be 2 different combination. We shall find them out first. 1 = {0,1}, 2 = {0,2}, 3= {1,2}. 
                //If we get a vector {1,2} this means pair {0,1} and pair {0,2} are pair with different charge and therefore particle-antiparticle pairs
                //we chose the pair that give the lowest pt lepton from top
                Pairs = FindParticlePairs(LeadLeptonInfo[0],LeadLeptonInfo[3],LeadLeptonInfo[4],LeadLeptonInfo[5],Electron_charge,Muon_charge);
                LowerDeltaPhiPair = FindLowerDeltaPhi(Pairs, Lepton1Phi,Lepton2Phi,Lepton3Phi);
                TopLepton = LeptonFromTop(LowerDeltaPhiPair);
                if (TwoDCutV2(LeadLeptonInfo[TopLepton],LeadLeptonInfo[TopLepton + 3],Electron_pt,Muon_pt,bJet_pt,Electron_phi,Muon_phi,bJet_phi,Electron_eta,Muon_eta,bJet_eta) == 1){
                    Pass2DCut = 1;
                }
                else{
                    Pass2DCut = 0;
                }
                //copy from first part with different leptons, should be optimized.
                LeptonPairId = RemainingLepton(TopLepton);
                //extrac the type of lepton and their position inside either Electron or Muon set
                int PhiLeptonType1 = LeadLeptonInfo[LeptonPairId[0]];
                int PhiLeptonId1 = LeadLeptonInfo[LeptonPairId[0] + 3];
                int PhiLeptonType2 = LeadLeptonInfo[LeptonPairId[1]];
                int PhiLeptonId2 = LeadLeptonInfo[LeptonPairId[1] + 3];

                PassPairCheck = 1;//this case automatically pass the particle-antiparticle pair check by its construction.
                //now, check if the remaining two can give a invariant mass < 20 GeV
                //construct the TLvector for these two leptons and then compute their invariant mass.
                if (PhiLeptonType1 == 1){
                    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
                }
                else{
                    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
                }
                if (PhiLeptonType2 == 1){
                    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
                }
                else{
                    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
                }

                if (hardware::InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
                    PassScalarMassCut = 1;
                }


                //check if they are close to each other
                //first, collect their phi info
                if (PhiLeptonType1 == 1){
                    PhiLepton1Phi = Electron_phi[PhiLeptonId1];
                }
                else{
                    PhiLepton1Phi = Muon_phi[PhiLeptonId1];
                }

                if (PhiLeptonType2 == 1){
                    PhiLepton2Phi = Electron_phi[PhiLeptonId2];
                }
                else{
                    PhiLepton2Phi = Muon_phi[PhiLeptonId2];
                }
                //check their relative phi
                if (abs(hardware::DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
                    PassRelativePhi = 1;
                }
                else{
                    PassRelativePhi = 0;
                }

                if (PassPairCheck == 1 && Pass2DCut == 1 && PassScalarMassCut == 1 && PassRelativePhi == 1){
                    //Pass all the test
                    PassAllTest = 1;
                }

            }
        }

        if (LeadLeptonInfo[0] == 2){//all Muons
            Lepton1Charge = Muon_charge[LeadLeptonInfo[3]];
            Lepton2Charge = Muon_charge[LeadLeptonInfo[4]];
            Lepton3Charge = Muon_charge[LeadLeptonInfo[5]];
            Lepton1Phi = Muon_phi[LeadLeptonInfo[3]];
            Lepton2Phi = Muon_phi[LeadLeptonInfo[4]];
            Lepton3Phi = Muon_phi[LeadLeptonInfo[5]];
            if (Lepton1Charge == Lepton2Charge && Lepton2Charge == Lepton3Charge){
                PassAntiparticleCheck = 0;
            }
            else{
                PassAntiparticleCheck = 1;
                //there will be 2 different combination. We shall find them out first. 1 = {0,1}, 2 = {0,2}, 3= {1,2}. 
                //If we get a vector {1,2} this means pair {0,1} and pair {0,2} are pair with different charge and therefore particle-antiparticle pairs
                //we chose the pair that give the lowest pt lepton from top
                Pairs = FindParticlePairs(LeadLeptonInfo[0],LeadLeptonInfo[3],LeadLeptonInfo[4],LeadLeptonInfo[5],Electron_charge,Muon_charge);
                LowerDeltaPhiPair = FindLowerDeltaPhi(Pairs, Lepton1Phi,Lepton2Phi,Lepton3Phi);
                TopLepton = LeptonFromTop(LowerDeltaPhiPair);
                if (TwoDCutV2(LeadLeptonInfo[TopLepton],LeadLeptonInfo[TopLepton + 3],Electron_pt,Muon_pt,bJet_pt,Electron_phi,Muon_phi,bJet_phi,Electron_eta,Muon_eta,bJet_eta) == 1){
                    Pass2DCut = 1;
                }
                else{
                    Pass2DCut = 0;
                }
                //copy from first part with different leptons, should be optimized
                LeptonPairId = RemainingLepton(TopLepton);
                //extrac the type of lepton and their position inside either Electron or Muon set
                int PhiLeptonType1 = LeadLeptonInfo[LeptonPairId[0]];
                int PhiLeptonId1 = LeadLeptonInfo[LeptonPairId[0] + 3];
                int PhiLeptonType2 = LeadLeptonInfo[LeptonPairId[1]];
                int PhiLeptonId2 = LeadLeptonInfo[LeptonPairId[1] + 3];

                PassPairCheck = 1;//this case automatically pass the particle-antiparticle pair check by its construction.
                //now, check if the remaining two can give a invariant mass < 20 GeV
                //construct the TLvector for these two leptons and then compute their invariant mass.
                if (PhiLeptonType1 == 1){
                    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
                }
                else{
                    ROOT::Math::PtEtaPhiMVector PhiLepton1TLvector(Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
                }
                if (PhiLeptonType2 == 1){
                    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
                }
                else{
                    ROOT::Math::PtEtaPhiMVector PhiLepton2TLvector(Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
                }

                if (hardware::InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
                    PassScalarMassCut = 1;
                }


                //check if they are close to each other
                //first, collect their phi info
                if (PhiLeptonType1 == 1){
                    PhiLepton1Phi = Electron_phi[PhiLeptonId1];
                }
                else{
                    PhiLepton1Phi = Muon_phi[PhiLeptonId1];
                }

                if (PhiLeptonType2 == 1){
                    PhiLepton2Phi = Electron_phi[PhiLeptonId2];
                }
                else{
                    PhiLepton2Phi = Muon_phi[PhiLeptonId2];
                }
                //check their relative phi
                if (abs(hardware::DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
                    PassRelativePhi = 1;
                }
                else{
                    PassRelativePhi = 0;
                }

                if (PassPairCheck == 1 && Pass2DCut == 1 && PassScalarMassCut == 1 && PassRelativePhi == 1 && PassJetRelativePhi == 1){
                    //Pass all the test
                    PassAllTest = 1;
                }

            }
        }
    }

    return {PassAllTest, TopLepton, LeptonPairId[0], LeptonPairId[1]};
}



//lepton pt, eta, phi and mass of the lepton. This function should take in all related value of electron/muon and choose based on the value of Electron/Muon id

float GetFloatLeptonProperty(int LeptonType,int LeptonId,RVec<float> ElectronProperty, RVec<float> MuonProperty){

    float LeptonFloat = -1.0;//by setting default value to -1 will give us a warning messeage if no lepton satisfying condition exist. Such event should be cut in preselections
    
    if(LeptonType == 1){//this is an electron
        LeptonFloat = ElectronProperty[LeptonId];//then the lepton's property used to construct 4 vector should be relavent electron property
    }
    else{
        if(LeptonId == 2){//this is a muon
            LeptonFloat = MuonProperty[LeptonId];
        }
    }

    return LeptonFloat;
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

};