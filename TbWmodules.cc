#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..
#include <cmath>//need to use the value of pi

using namespace ROOT::VecOps;
//we want to 1. examine if there's at least a 50GeV lepton 2. If so, find the back to back AK8 and AK4 jets (we know there would at most be one)

RVec<int> PickDijets(RVec<float> FatJet_phi, RVec<float> Jet_phi, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Jet_btagCMVA){
    int FatJetidx = -1;
    int Jetidx = -1;
    int Leptonidx = -1;
    int Electronidx = -1;
    int Muonidx = -1;
    int Leptonpt=0
    if (Electron_pt.size() < 1){
        if (Muon_pt.size() < 1){break;}
        else {
            for (iMuon = 0; iMuon < Muon_pt.size(); iMuon++){
                if (Muon_pt[iMuon]>50){
                    Muonidx = iMuon//give the first Muon sastifying our condition
                    Leptonidx = 1;//represent Muon as 1
                    Leptonpt = Muon_pt[iMuon]//The momentum of lepton is given by this muon
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(FatJet_phi[iFatJet]-Jet_phi[iJet] > M_PI/2) && Jet_btagCMVA[iJet] > 0.8){
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
                    Leptonpt = Electron_pt[iElectron]//momentum is given by electron
                    for (int iJet = 0; iJet < Jet_phi.size(); iJet++){//find the back to back jets
                        for (int iFatJet =0; iFatJet < FatJet_phi.size(); iFatJet++){
                            if (abs(FatJet_phi[iFatJet]-Jet_phi[iJet] > M_PI/2 && Jet_btagCMVA[iJet] > 0.8)){
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
    return{FatJetidx,Jetidx,Leptonidx,Leptonpt,Electronidx,Muonidx}

}
