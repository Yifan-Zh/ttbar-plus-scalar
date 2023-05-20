#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"

void PlotTopPt() {
        TFile *f = TFile::Open("ttbarsnapshot_17.root_17_1of1.root","READ");
        TTree *tree = (TTree*)f->Get("Events");
        TH1F *myh = new TH1F("myh","Top_pt",50,0,1500);
        tree->Draw("Top_pt>>myh");
        TCanvas *c = new TCanvas("c","c");
        c->cd();
        myh->Draw();
        c->Print("Top_pt.pdf");
}
