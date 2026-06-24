#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"

#include <iostream>

static TH1F *loadNormalizedHist(const char *fname, const char *hname,
                                const char *newName, Color_t color)
{
    TFile *file = TFile::Open(fname);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << fname << std::endl;
        return nullptr;
    }

    TH1F *hist = (TH1F*)file->Get(hname);
    if (!hist) {
        std::cerr << "Cannot find " << hname << " in " << fname << std::endl;
        file->Close();
        return nullptr;
    }

    TH1F *clone = (TH1F*)hist->Clone(newName);
    clone->SetDirectory(nullptr);
    file->Close();

    const double max = clone->GetMaximum();
    if (max > 0.) clone->Scale(1. / max);

    clone->SetStats(0);
    clone->SetLineColor(color);
    clone->SetMarkerColor(color);
    clone->SetLineWidth(2);
    return clone;
}

static void drawPair(TPad *pad, TH1F *h1, TH1F *h2,
                     const char *title, const char *label1, const char *label2)
{
    pad->cd();
    h1->SetTitle(title);
    h1->GetYaxis()->SetTitle("Normalized counts");
    h1->GetYaxis()->SetRangeUser(0., 1.15);
    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.62, 0.74, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h1, label1, "l");
    leg->AddEntry(h2, label2, "l");
    leg->Draw();
}

void vertex_recon_plot(const char *dir = "./")
{
    TH1F *h24685 = loadNormalizedHist(Form("%s/vertex_recon_moller_24685.root", dir),
                                      "h_moller_vz", "h_moller_vz_24685", kBlue + 1);
    TH1F *h24941 = loadNormalizedHist(Form("%s/vertex_recon_moller_24941.root", dir),
                                      "h_moller_vz", "h_moller_vz_24941", kRed + 1);
    TH1F *h24683 = loadNormalizedHist(Form("%s/vertex_recon_moller_24683.root", dir),
                                      "h_moller_vz", "h_moller_vz_24683", kBlue + 1);
    TH1F *h24942 = loadNormalizedHist(Form("%s/vertex_recon_moller_24942.root", dir),
                                      "h_moller_vz", "h_moller_vz_24942", kRed + 1);

    if (!h24685 || !h24941 || !h24683 || !h24942) return;

    TCanvas *c = new TCanvas("c_moller_vz_compare", "Moller vertex Z comparison", 1200, 500);
    c->Divide(2, 1);
    drawPair((TPad*)c->cd(1), h24685, h24941,
             "Moller Vertex Z: 24685 vs 24941;z_{vertex} (mm);Normalized counts",
             "Run 24685", "Run 24941");
    drawPair((TPad*)c->cd(2), h24683, h24942,
             "Moller Vertex Z: 24683 vs 24942;z_{vertex} (mm);Normalized counts",
             "Run 24683", "Run 24942");

    c->SaveAs("vertex_recon_moller_compare.png");
    c->SaveAs("vertex_recon_moller_compare.pdf");

    TFile *fout = TFile::Open("vertex_recon_moller_compare.root", "RECREATE");
    h24685->Write();
    h24941->Write();
    h24683->Write();
    h24942->Write();
    c->Write();
    fout->Close();

    std::cout << "Saved vertex_recon_moller_compare.png, .pdf, and .root" << std::endl;
}
