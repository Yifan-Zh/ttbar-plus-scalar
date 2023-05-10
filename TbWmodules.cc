#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi

using namespace ROOT::VecOps;

RVec<int> PickDijets(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int jet0Idx = -1;
    int jet1Idx = -1;
    for (int ijet = 0; ijet < pt.size(); ijet++) {
        if (jet1Idx == -1) {
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
                if (jet0Idx == -1) {
                    jet0Idx = ijet;
                } else {
                    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
                        jet1Idx = ijet;
                        break;
                    }
                }
            }
        }       
    }
    return {jet0Idx,jet1Idx};
}


RVec<int> PickDijets(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Electron_pt, RVec<float> Muon_pt){
    int FatJetidx = -1;
    int Jetidx = -1;
    int Leptonidx = -1;
    int Electronidx = -1;
    int Muonidx = -1;
    if (Electron_pt.size() < 1){
        if (Muon_pt.size() < 1){break;}
        else {
            for (iMuon = 0; iMuon < Muon_pt.size(); iMuon++){
                if (Muon_pt[iMuon]>50){
                    Muonidx = iMuon//give the first Muon sastifying our condition
                    Leptonidx = 1;//represent Muon as 1
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(FatJet_phi[iFatJet]-Jet_phi[iJet] > M_PI/2)){
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
        for (iElectron = 0; iElectron < Electron_pt.size(); iElectron++){
                if (Electron_pt[iElectron]>50){
                    Electronidx = iElectron//give the first Electron sastifying our condition
                    Leptonidx = 1;//represent Electron as 1
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(FatJet_phi[iFatJet]-Jet_phi[iJet] > M_PI/2)){
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
    return{FatJetidx,Jetidx,Leptonidx,Electronidx,Muonidx}

}