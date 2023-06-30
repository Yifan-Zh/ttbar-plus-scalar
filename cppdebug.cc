#include "ROOT/RVec.hxx"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi, and some trig/hyperbolic functions.
#include <map>
#include <algorithm>
#include <Math/Vector4D.h>
#include <utility>
#include <TMath.h>
#include <iostream>
#include <cassert>

using namespace ROOT::VecOps;
// Include the function definition here

float DeltaPhi(float phi1,float phi2) {
    float result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
}

double InvariantMass(RVec<ROOT::Math::PtEtaPhiMVector> vects) {
    ROOT::Math::PtEtaPhiMVector sum;
    sum.SetCoordinates(0,0,0,0);
    for (size_t i = 0; i < vects.size(); i++) {
        sum = sum + vects[i];
    }
    return sum.M();
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
        if(abs(DeltaPhi(Electron_phi[iElectronId],bJet_phi)) > 0.4){
            Crel_phi = 1;
        }
    }
    else{
        if(LeptonType == 2){
            iMuonId = LeptonId;
            if(LepbJetPtRel(Muon_pt[iMuonId], Muon_phi[iMuonId], Muon_eta[iMuonId], bJet_pt, bJet_phi, bJet_eta) > (25*25)){
                Crel_pt = 1;
            }
            if(abs(DeltaPhi(Muon_phi[iMuonId],bJet_phi)) > 0.4){
                Crel_phi = 1;
            }
        }
    }

    if (Crel_pt == 1 || Crel_phi == 1){
        Pass2DCut = 1;
    }

    return Pass2DCut;
}

void testTwoDCutV2(){
    //test case 1, both condition pass.
    int LeptonType = 1;
    int LeptonId = 1;
    RVec<float> ElectronPt = {10.5, 70.21, 8.9, 12.1};
    RVec<float> MuonPt = { 7.6, 11.3, 9.8, 11.2 };
    float bJetPt = 166.8;
    RVec<float> ElectronPhi = {1.2, 3.1, 0.65, 0.1};
    RVec<float> MuonPhi = {1.0, 2.1, 0.85, 0.14};
    float bJetPhi = 0.08;
    RVec<float> ElectronEta = {0.1,0.2,0.3,0.4};
    RVec<float> MuonEta = {0.01,0.11,0.22,0.33};
    float bJetEta = 1;

    int expected = 1;
    int result = TwoDCutV2(LeptonType,LeptonId,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta);
    assert(result == expected);
    std::cout << "TwoDCut Test Case 1 Pass" <<std::endl;

    //test case 2, pass relative angle only
    LeptonType = 2;
    LeptonId = 0;
    ElectronPt = {10.5, 66.8, 8.9, 12.1};
    MuonPt = { 65.2, 11.3, 9.8, 11.2 };
    bJetPt = 66.8;
    ElectronPhi = {1.2, 3.1, 0.65, 0.1};
    MuonPhi = {1.0, 2.1, 0.85, 0.14};
    bJetPhi = 0.08;
    ElectronEta = {0.1,0.2,0.3,0.4};
    MuonEta = {0.01,0.11,0.22,0.33};
    bJetEta = 1;

    expected = 1;
    result = TwoDCutV2(LeptonType,LeptonId,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta);
    assert(result == expected);
    std::cout << "TwoDCut Test Case 2 Pass" <<std::endl;

    //test case 3, pass relativept only
    LeptonType = 2;
    LeptonId = 0;
    ElectronPt = {10.5, 66.8, 8.9, 12.1};
    MuonPt = { 8000.2, 11.3, 9.8, 11.2 };
    bJetPt = 66.8;
    ElectronPhi = {1.2, 3.1, 0.65, 0.1};
    MuonPhi = {0.4, 2.1, 0.85, 0.14};
    bJetPhi = 0.08;
    ElectronEta = {0.1,0.2,0.3,0.4};
    MuonEta = {0.01,0.11,0.22,0.33};
    bJetEta = 1;

    expected = 1;
    result = TwoDCutV2(LeptonType,LeptonId,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta);
    assert(result == expected);
    std::cout << "TwoDCut Test Case 3 Pass" <<std::endl;

    //test case 4, pass neither
    LeptonType = 2;
    LeptonId = 0;
    ElectronPt = {10.5, 66.8, 8.9, 12.1};
    MuonPt = { 65.9, 11.3, 9.8, 11.2 };
    bJetPt = 66.8;
    ElectronPhi = {1.2, 3.1, 0.65, 0.1};
    MuonPhi = {0.4, 2.1, 0.85, 0.14};
    bJetPhi = 0.08;
    ElectronEta = {0.1,0.2,0.3,0.4};
    MuonEta = {0.99,0.11,0.22,0.33};
    bJetEta = 1;

    expected = 0;
    result = TwoDCutV2(LeptonType,LeptonId,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta);
    assert(result == expected);
    std::cout << "TwoDCut Test Case 4 Pass" <<std::endl;

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
        PairDeltaPhi = abs(DeltaPhi(Lepton1Phi,Lepton2Phi));
    }
    
    if (PairIndex == 2){
        PairDeltaPhi = abs(DeltaPhi(Lepton1Phi,Lepton3Phi));
    }

    if (PairIndex ==3){
        PairDeltaPhi = abs(DeltaPhi(Lepton2Phi,Lepton3Phi));
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

void testFindLowerDeltaPhi(){
    
    //test case 1
    RVec<int> Pairs = {1,2};
    float Lepton1Phi = 0.5;
    float Lepton2Phi = 1.3;
    float Lepton3Phi = 1.1;
    int expected = 2;
    int result = FindLowerDeltaPhi(Pairs, Lepton1Phi, Lepton2Phi, Lepton3Phi);
    assert (result == expected);
    std::cout << "FindLowerDeltaPhi Test Case 1 Pass" <<std::endl;

    //test case 2
    Pairs = {1,3};
    Lepton1Phi = 0.5;
    Lepton2Phi = 1.3;
    Lepton3Phi = 0.22;
    expected = 1;
    result = FindLowerDeltaPhi(Pairs, Lepton1Phi, Lepton2Phi, Lepton3Phi);
    assert (result == expected);
    std::cout << "FindLowerDeltaPhi Test Case 2 Pass" <<std::endl;

    //test case 3
    Pairs = {2,3};
    Lepton1Phi = 0.5;
    Lepton2Phi = 1.3;
    Lepton3Phi = 1.1;
    expected = 3;
    result = FindLowerDeltaPhi(Pairs, Lepton1Phi, Lepton2Phi, Lepton3Phi);
    assert (result == expected);
    std::cout << "FindLowerDeltaPhi Test Case 3 Pass" <<std::endl;
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
        ThreeSameLepton = true;
    }
    else{
        ThreeSameLepton = false;
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
            PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
        }
        else{
            PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
        }
        if (PhiLeptonType2 == 1){
            PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
        }
        else{
            PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
        }

        if (InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
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
        if (abs(DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
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
                    PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
                }
                else{
                    PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
                }
                if (PhiLeptonType2 == 1){
                    PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
                }
                else{
                    PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
                }

                if (InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
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
                if (abs(DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
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
                    PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId1],Electron_eta[PhiLeptonId1],Electron_phi[PhiLeptonId1],0.000511);
                }
                else{
                    PhiLepton1TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId1],Muon_eta[PhiLeptonId1],Muon_phi[PhiLeptonId1],0.105658);
                }
                if (PhiLeptonType2 == 1){
                    PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Electron_pt[PhiLeptonId2],Electron_eta[PhiLeptonId2],Electron_phi[PhiLeptonId2],0.000511);
                }
                else{
                    PhiLepton2TLvector = ROOT::Math::PtEtaPhiMVector (Muon_pt[PhiLeptonId2],Muon_eta[PhiLeptonId2],Muon_phi[PhiLeptonId2],0.105658);
                }

                if (InvariantMass({PhiLepton1TLvector,PhiLepton2TLvector}) < 20){
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
                if (abs(DeltaPhi(PhiLepton1Phi,PhiLepton2Phi)) < M_PI/4){
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
    }
    std::cout<<PassPairCheck<<" "<<Pass2DCut<<" "<<PassScalarMassCut<<" "<<PassRelativePhi<<std::endl;
    return {PassAllTest, TopLepton, LeptonPairId[0], LeptonPairId[1]};
}

void testLeptonCategorize(){

    //test case 1
    RVec<int> LeadLeptonInfo = {1,1,2,0,1,0};//lepton 2 (the third one) is from W
    RVec<float> ElectronPt = {150.1,140.1,130.1};
    RVec<float> MuonPt = {122.4,83.2,58.7};
    float bJetPt = 177.6;
    RVec<float> ElectronPhi = {-1.8,0.62,1.73};
    RVec<float> MuonPhi = {-0.56,0.239,2.1};
    float bJetPhi = 0.8;
    RVec<float> ElectronEta = {0.5,1.1,2.1};
    RVec<float> MuonEta = {0.33,0.521,1.53};
    float bJetEta = 0.52;
    RVec<int> ElectronCharge = {1,-1,1};
    RVec<int> MuonCharge = {1,-1,-1};

    std::cout << "Test 1:No grammar/logistic error in LeptonCategorize" << std::endl;

    RVec<int> result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    //test case 2
    LeadLeptonInfo = {1,2,2,0,1,2};//lepton 0 is from W
    ElectronPt = {150.1,140.1,130.1};
    MuonPt = {122.4,53.2,52.9};
    bJetPt = 177.6;
    ElectronPhi = {-1.8,0.62,1.73};
    MuonPhi = {-0.56,2.39,2.1};
    bJetPhi = 0.8;
    ElectronEta = {0.5,1.1,2.1};
    MuonEta = {0.33,0.521,0.53};
    bJetEta = 0.52;
    ElectronCharge = {1,-1,1};
    MuonCharge = {1,-1,1};

    std::cout << "Test 2:No grammar/logistic error in LeptonCategorize" << std::endl;

    result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    //test case 3
    LeadLeptonInfo = {2,2,2,0,1,2};//lepton 0 is from W
    ElectronPt = {150.1,140.1,130.1};
    MuonPt = {122.4,53.2,52.9};
    bJetPt = 177.6;
    ElectronPhi = {-1.8,0.62,1.73};
    MuonPhi = {-0.56,2.39,2.1};
    bJetPhi = 0.8;
    ElectronEta = {0.5,1.1,2.1};
    MuonEta = {0.33,0.521,0.53};
    bJetEta = 0.52;
    ElectronCharge = {1,-1,1};
    MuonCharge = {1,-1,1};

    std::cout << "Test 3:No grammar/logistic error in LeptonCategorize" << std::endl;
    result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    //test case 4, expect to fail
    LeadLeptonInfo = {2,2,2,0,1,2};//identical leptons
    ElectronPt = {150.1,140.1,130.1};
    MuonPt = {122.4,53.2,52.9};
    bJetPt = 177.6;
    ElectronPhi = {-1.8,0.62,1.73};
    MuonPhi = {-0.56,2.39,2.1};
    bJetPhi = 0.8;
    ElectronEta = {0.5,1.1,2.1};
    MuonEta = {0.33,0.521,0.53};
    bJetEta = 0.52;
    ElectronCharge = {1,-1,1};
    MuonCharge = {-1,-1,-1};

    std::cout << "Test 4:No grammar/logistic error in LeptonCategorize" << std::endl;
    result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    //test case 5, expect to fail
    LeadLeptonInfo = {1,2,2,0,1,2};//does not form a pair
    ElectronPt = {150.1,140.1,130.1};
    MuonPt = {122.4,53.2,52.9};
    bJetPt = 177.6;
    ElectronPhi = {-1.8,0.62,1.73};
    MuonPhi = {-0.56,2.39,2.1};
    bJetPhi = 0.8;
    ElectronEta = {0.5,1.1,2.1};
    MuonEta = {0.33,0.521,0.53};
    bJetEta = 0.52;
    ElectronCharge = {1,-1,1};
    MuonCharge = {1,-1,-1};

    std::cout << "Test 5:No grammar/logistic error in LeptonCategorize" << std::endl;

    result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    //test case 6, expect to fail due to relative phi
    LeadLeptonInfo = {1,2,2,0,1,2};//lepton 0 is from W
    ElectronPt = {150.1,140.1,130.1};
    MuonPt = {122.4,53.2,52.9};
    bJetPt = 177.6;
    ElectronPhi = {-1.8,0.62,1.73};
    MuonPhi = {-0.56,2.39,0.1};
    bJetPhi = 0.8;
    ElectronEta = {0.5,1.1,2.1};
    MuonEta = {0.33,0.521,0.53};
    bJetEta = 0.52;
    ElectronCharge = {1,-1,1};
    MuonCharge = {1,-1,1};

    std::cout << "Test 6:No grammar/logistic error in LeptonCategorize" << std::endl;

    result = LeptonCategorize (LeadLeptonInfo,ElectronPt,MuonPt,bJetPt,ElectronPhi,MuonPhi,bJetPhi,ElectronEta,MuonEta,bJetEta,ElectronCharge,MuonCharge);
    for (int i = 0; i < result.size(); i++){
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

}

int main(){
    std::cout<< "Compile complete" <<std::endl;
    testTwoDCutV2();
    testFindLowerDeltaPhi();
    testLeptonCategorize();

    return 0;
}