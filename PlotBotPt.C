#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"

void PlotBotPt() {
        TFile *f = TFile::Open("ttbarsnapshot_17.root_17_1of1.root","READ");
        TTree *tree = (TTree*)f->Get("Events");
        TH1F *myh = new TH1F("myh","Bot_pt",50,0,1500);
        tree->Draw("Bot_pt>>myh");
        TCanvas *c = new TCanvas("c","c");
        c->cd();
        myh->Draw();
        c->Print("Bot_pt.pdf");
}