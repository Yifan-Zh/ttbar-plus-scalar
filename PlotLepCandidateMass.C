#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"

void PlotLepCandidateMass() {
        TFile *f = TFile::Open("ttbarsnapshot_17.root_17_1of1.root","READ");
        TTree *tree = (TTree*)f->Get("Events");
        TH1F *myh = new TH1F("myh","LepCandidate_mass",50,0,1500);
        tree->Draw("LepCandidate_mass>>myh");
        TCanvas *c = new TCanvas("c","c");
        c->cd();
        myh->Draw();
        c->Print("LepCandidate_mass.pdf");
}
