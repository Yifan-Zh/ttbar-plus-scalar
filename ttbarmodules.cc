#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi

using namespace ROOT::VecOps;
//we want to 1. examine if there's at least a 50GeV lepton 2. If so, find the back to back AK8 and AK4 jets (we know there would at most be one)
//Problem left to fix: modify the code to pick to most energetic lepton

RVec<int> PickDijets(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Jet_btagCMVA){
    int FatJetidx = -1;
    int Jetidx = -1;
    int Leptonidx = -1;
    int Electronidx = -1;
    int Muonidx = -1;
    int Leptonpt = 0;
    int C_Lepton_pt = -1;
    if (Electron_pt.size() < 1){
        if (Muon_pt.size() < 1){Leptonidx = -1;}
        else {
            for (int iMuon = 0; iMuon < Muon_pt.size(); iMuon++){
                if (Muon_pt[iMuon]>50){
                    Muonidx = iMuon;//give the first Muon sastifying our condition
                    Leptonidx = 2;//represent Muon as 2
                    Leptonpt = Muon_pt[iMuon];//The momentum of lepton is given by this muon
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2 && Jet_btagCMVA[iJet] > 0.8){
                                FatJetidx = iFatJet;
                                Jetidx = iJet;
                                break;
                            }
                        }
                    }
                    break;
                }

            }
           
        }
    }
    else {
        for (int iElectron = 0; iElectron < Electron_pt.size(); iElectron++){
                if (Electron_pt[iElectron]>50){
                    Electronidx = iElectron;//give the first Electron sastifying our condition
                    Leptonidx = 1;//represent Electron as 1
                    Leptonpt = Electron_pt[iElectron];//momentum is given by electron
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2 && Jet_btagCMVA[iJet] > 0.8){
                                FatJetidx = iFatJet;
                                Jetidx = iJet;
                                break;
                            }
                        }
                    }
                    break;
                }

        }
    }

    if(Leptonpt>50){
        C_Lepton_pt = 1;
    }
    return {FatJetidx,Jetidx,Leptonidx,C_Lepton_pt,Electronidx,Muonidx};//should fix this later. The C_Lepton_pt standard is useless.

}


RVec<int> PickDijetsV2(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Jet_btagCMVA){
    int FatJetidx = -1;
    int Jetidx = -1;
    int Leptonidx = -1;
    int Electronidx = -1;
    int Muonidx = -1;
    int Leptonpt = 0;
    int C_Lepton_pt = -1;
    bool exitLeptonloop = false;
    bool exitJetloop = false;

    if (Electron_pt.size() < 1){
        if (Muon_pt.size() < 1){Leptonidx = -1;}
        else {
            for (int iMuon = 0; iMuon < Muon_pt.size(); iMuon++){//no electron, but have muon
                if (Muon_pt[iMuon]>50){
                    Muonidx = iMuon;//give the first Muon sastifying our condition
                    Leptonidx = 2;//represent Muon as 2
                    Leptonpt = Muon_pt[iMuon];//The momentum of lepton is given by this muon
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
                    break;
                }

            }
           
        }
    }
    else {
        if (Muon_pt.size() < 1){//no muon, but have electron
            for (int iElectron = 0; iElectron < Electron_pt.size(); iElectron++){
                if (Electron_pt[iElectron]>50){
                    Electronidx = iElectron;//give the first Electron sastifying our condition
                    Leptonidx = 1;//represent Electron as 1
                    Leptonpt = Electron_pt[iElectron];//momentum is given by electron
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
                    break;
                }

            }
        }

        else {//both exsit
            for(int iElectron = 0; iElectron < Electron_pt.size() && exitLeptonloop ==false; iElectron++){
                for (int iMuon = 0; iMuon < Muon_pt.size(); iMuon++){
                    if(Electron_pt[iElectron]>Muon_pt[iMuon] && Electron_pt[iElectron] > 50){
                        Electronidx = iElectron;
                        Muonidx = iMuon;
                        Leptonidx = 1;
                        for (int iJet = 0; iJet < Jet_phi.size() && exitJetloop == false; iJet++){//find the back to back jets
                            for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                                if (abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2 && Jet_btagCMVA[iJet] > 0.8){
                                    FatJetidx = iFatJet;
                                    Jetidx = iJet;
                                    exitJetloop = true;
                                    break;
                                }
                            }
                        }
                        exitLeptonloop = true;
                        break;
                    }
                    else{
                        if(Muon_pt[iMuon]>Electron_pt[iElectron] && Muon_pt[iMuon] > 50){
                            Electronidx = iElectron;
                            Muonidx = iMuon;
                            Leptonidx = 2;
                            for (int iJet = 0; iJet < Jet_phi.size() && exitJetloop == false; iJet++){//find the back to back jets
                                for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                                    if (abs(hardware::DeltaPhi(FatJet_phi[iFatJet],Jet_phi[iJet])) > M_PI/2 && Jet_btagCMVA[iJet] > 0.8){
                                        FatJetidx = iFatJet;
                                        Jetidx = iJet;
                                        exitJetloop = true;
                                        break;
                                    }
                                }
                            }
                            exitLeptonloop = true;
                            break;
                        }

                    }
                }
            }

        }
    }

    if(Leptonpt>50){
        C_Lepton_pt = 1;
    }
    return {FatJetidx,Jetidx,Leptonidx,C_Lepton_pt,Electronidx,Muonidx};//should fix this later. The C_Lepton_pt standard is useless.

}

//a short program to extract ith component of an column with integer type data. Use when ObjectFromCollection does not satisfy our purpose
RVec<int> CreateIntColumn(int column, RVec<int> CalledColumn){
    int value = 0;
    value = CalledColumn[column];
    return {value};
}

//similar program to extrac float type data
RVec<float> CreateFloatColumn(int column, RVec<float> CalledColumn){
    float value = 0;
    value = CalledColumn[column];
    return {value};
}

RVec<int> TwoDCut(int LeptonId, int ElectronId, int MuonId, RVec<float> Electron_jetPtRelv2, RVec<float> Muon_jetPtRelv2, float bJet_phi, RVec<float> Electron_phi, RVec<float> Muon_phi){// for some reason, ObjectFromCollection return normal number instead of RVector
    int Crel_pt = -1;//relative pt criteria
    int Crel_phi = -1;//relative phi criteria
    int iElectronId = 0;
    int iMuonId = 0;
    //This selection prioritize the electron over muon. This need to be fixed later.
    //ElectronId = -1 means no electron exist in our system. Do not access any information about electron in this case
    if (LeptonId == 1){
        iElectronId = ElectronId;
        if(Electron_jetPtRelv2[iElectronId] > 25){
            Crel_pt = 1;
        }
        if(abs(hardware::DeltaPhi(Electron_phi[iElectronId],bJet_phi)) > 0.4){
            Crel_phi = 1;
        }
    }
    else{
        if(LeptonId == 2){//uncessary since all {-1,-1} have been cut in preselection
            iMuonId = MuonId;
            if(Muon_jetPtRelv2[iMuonId] >25){
                Crel_pt = 1;
            }
            if(abs(hardware::DeltaPhi(Muon_phi[iMuonId],bJet_phi)) > 0.4){
                Crel_phi = 1;
            }
        }
    }

    return {Crel_pt,Crel_phi};
}

//lepton pt, eta, phi and mass of the lepton. This function should take in all related value of electron/muon and choose based on the value of Electron/Muon id

float GetFloatLeptonProperty(int LeptonId, int ElectronId, int MuonId, RVec<float> ElectronProperty, RVec<float> MuonProperty){
    int LeptonIndex = -1;//this parameter will be used to point to lepton properties;
    float LeptonFloat = -1;//by setting default value to -1 will give us a warning messeage if no lepton satisfying condition exist. Such event should be cut in preselections
    
    if(LeptonId == 1){//this is an electron
        LeptonIndex = ElectronId;
        LeptonFloat = ElectronProperty[LeptonIndex];//then the lepton's property used to construct 4 vector should be relavent electron property
    }
    else{
        if(LeptonId == 2){//this is a muon
            LeptonIndex = MuonId;
            LeptonFloat = MuonProperty[LeptonIndex];
        }
    }

    return LeptonFloat;
}

float LeptonicCandidatePt(float Bot_pt, float Lepton_pt, float Bot_phi, float Lepton_phi){
    float px = 0;
    float py = 0;
    float pt_tot = 0;
    px = Bot_pt * cos(Bot_phi) + Lepton_pt * cos(Lepton_phi);
    py = Bot_pt * sin(Bot_phi) + Lepton_pt * sin(Lepton_phi);
    pt_tot = sqrt(px*px + py*py);
    return pt_tot;
}
