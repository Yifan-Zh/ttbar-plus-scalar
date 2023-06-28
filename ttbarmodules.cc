#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi, and some trig/hyperbolic functions.

using namespace ROOT::VecOps;
//we want to 1. examine if there's at least a 50GeV lepton 2. If so, find the back to back AK8 and AK4 jets (we know there would at most be one)
//Problem left to fix: modify the code to pick to most energetic lepton


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
            for(int iElectron = 0; iElectron < Electron_pt.size() && exitLeptonloop == false; iElectron++){
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



RVec<int> TwoDCut(int LeptonId, int ElectronId, int MuonId, RVec<float> Electron_jetPtRelv2, RVec<float> Muon_jetPtRelv2, float bJet_phi, RVec<float> Electron_phi, RVec<float> Muon_phi){// for some reason, ObjectFromCollection return normal number instead of RVector
    int Crel_pt = -1;//relative pt criteria
    int Crel_phi = -1;//relative phi criteria
    int iElectronId = 0;
    int iMuonId = 0;
    //This selection has a problem, it only examine the property of the most energetic lepton in the system.
    //Might need to be changed later such that it will go through all the electrons and muons.
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
        if(LeptonId == 2){//no need to check LeptonId=-1 since all {-1,-1} have been cut in preselection
            iMuonId = MuonId;
            if(Muon_jetPtRelv2[iMuonId] >25){//problem might be caused by this
                Crel_pt = 1;
            }
            if(abs(hardware::DeltaPhi(Muon_phi[iMuonId],bJet_phi)) > 0.4){
                Crel_phi = 1;
            }
        }
    }

    return {Crel_pt,Crel_phi};
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

RVec<int> TwoDCutV2(int LeptonId, int ElectronId, int MuonId, RVec<float> Electron_pt, RVec<float> Muon_pt, float bJet_pt, RVec<float> Electron_phi, RVec<float> Muon_phi, float bJet_phi, RVec<float> Electron_eta, RVec<float> Muon_eta, float bJet_eta){// for some reason, ObjectFromCollection return normal number instead of RVector
    int Crel_pt = -1;//relative pt criteria
    int Crel_phi = -1;//relative phi criteria
    int iElectronId = 0;
    int iMuonId = 0;
    //This selection has a problem, it only examine the property of the most energetic lepton in the system.
    //Might need to be changed later such that it will go through all the electrons and muons.
    if (LeptonId == 1){
        iElectronId = ElectronId;
        if(LepbJetPtRel(Electron_pt[iElectronId], Electron_phi[iElectronId], Electron_eta[iElectronId], bJet_pt, bJet_phi, bJet_eta) > (25*25)){
            Crel_pt = 1;
        }
        if(abs(hardware::DeltaPhi(Electron_phi[iElectronId],bJet_phi)) > 0.4){
            Crel_phi = 1;
        }
    }
    else{
        if(LeptonId == 2){//no need to check LeptonId=-1 since all {-1,-1} have been cut in preselection
            iMuonId = MuonId;
            if(LepbJetPtRel(Muon_pt[iMuonId], Muon_phi[iMuonId], Muon_eta[iMuonId], bJet_pt, bJet_phi, bJet_eta) > (25*25)){
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

    float LeptonFloat = -1.0;//by setting default value to -1 will give us a warning messeage if no lepton satisfying condition exist. Such event should be cut in preselections
    
    if(LeptonId == 1){//this is an electron
        LeptonFloat = ElectronProperty[ElectronId];//then the lepton's property used to construct 4 vector should be relavent electron property
    }
    else{
        if(LeptonId == 2){//this is a muon
            LeptonFloat = MuonProperty[MuonId];
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

int FindbJets (RVec<float> Jet_TaggerScore){
    int n = 0;//number count initializer
    for (int i = 0; i < Jet_TaggerScore.size(); i++){
        if (Jet_TaggerScore[i] > 0.8){//loose standard tagger score
            n = n + 1;
        }
    }
    return n;
}


const ROOT::RVec<int> bJetsWithAK4Around(const ROOT::RVec<float>& bquark_pt,const ROOT::RVec<float>& bquark_eta, const ROOT::RVec<float>& bquark_phi, const ROOT::RVec<float>& jet_eta, const ROOT::RVec<float>& jet_phi){
    std::size_t quarksize = bquark_eta.size();
    std::size_t jetsize = jet_eta.size();
    RVec<int> bwithAK4(quarksize);
    for (std::size_t i = 0; i < quarksize; i++){
        for (std::size_t j = 0; j < jetsize; j++){
            if (abs(ROOT::VecOps::DeltaR(bquark_eta[i],jet_eta[j],bquark_phi[i],jet_phi[j])) > 0.4 && bquark_pt[i] > 20){
                bwithAK4[i] = -1;//if we have a bquark with >20GeV pt, and it does not have AK4 jets accompanying it, then sometihing went off.
            }
            else{
                bwithAK4[i] = 1;//either its pt is so low so we don't care or it actual pass the selection standard.
                break;//go to next quark
            }
        }
    }

    return bwithAK4;
}