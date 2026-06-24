
// Overlay two live-charge-normalized Mott yields.
//
// The input histograms were already normalized by live charge when they were
// produced in background.C, so no additional scaling is applied here.
//
// Usage:
//   root -l 'emptyTarget/compare_flowRate.C'

#include <algorithm>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

// Fill in the two input file names here.
static const char *file600Flow = "../../C_600cc.root";
static const char *file1300Flow = "../../C_1300cc.root";

static TH1 *loadNormalizedMottYield(const char *fileName, const char *tag)
{
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << fileName << std::endl;
        delete file;
        return nullptr;
    }

    TH1 *hist = dynamic_cast<TH1 *>(file->Get("mott_yield"));
    if (!hist) {
        std::cerr << "mott_yield not found in " << fileName << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    TH1 *clone = dynamic_cast<TH1 *>(
        hist->Clone(Form("mott_yield_%s_live_charge", tag)));
    clone->SetDirectory(nullptr);
    if (clone->GetSumw2N() == 0) clone->Sumw2();

    file->Close();
    delete file;

    std::cout << "Loaded " << fileName
              << " (integral = " << clone->Integral() << ")" << std::endl;
    return clone;
}

void compare_flowRate()
{
    gStyle->SetOptStat(0);

    TH1 *h600Flow = loadNormalizedMottYield(file600Flow, "600flow");
    TH1 *h800Flow = loadNormalizedMottYield(file1300Flow, "1300flow");
    if (!h600Flow || !h800Flow) return;

    const double thetaMin = 0.5;
    const double thetaMax = 3.8;

    h600Flow->SetTitle(
        "Type C Mott Yield: 600flow vs 1300flow;"
        "Reconstructed Scattering Angle [deg];"
        "Count / live charge");
    h600Flow->SetLineColor(kBlue + 1);
    h600Flow->SetMarkerColor(kBlue + 1);
    h600Flow->SetMarkerStyle(20);
    h600Flow->SetMarkerSize(0.8);
    h600Flow->SetLineWidth(2);

    h800Flow->SetLineColor(kRed + 1);
    h800Flow->SetMarkerColor(kRed + 1);
    h800Flow->SetMarkerStyle(21);
    h800Flow->SetMarkerSize(0.8);
    h800Flow->SetLineWidth(2);

    h600Flow->GetXaxis()->SetRangeUser(thetaMin, thetaMax);
    h600Flow->GetXaxis()->CenterTitle();
    h600Flow->GetYaxis()->CenterTitle();
    h600Flow->GetXaxis()->SetTitleSize(0.045);
    h600Flow->GetYaxis()->SetTitleSize(0.045);
    h600Flow->GetYaxis()->SetTitleOffset(1.25);

    double maximum = std::max(h600Flow->GetMaximum(), h800Flow->GetMaximum());
    h600Flow->SetMinimum(0.);
    h600Flow->SetMaximum(1.20 * maximum);

    TCanvas *canvas = new TCanvas(
        "c_compare_flowRate", "Mott Yield Comparison", 1000, 700);
    canvas->SetGrid();
    canvas->SetLeftMargin(0.13);
    canvas->SetBottomMargin(0.12);

    h600Flow->Draw("E1");
    h800Flow->Draw("E1 SAME");

    TLegend *legend = new TLegend(0.67, 0.74, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h600Flow, "600flow", "PE");
    legend->AddEntry(h800Flow, "1300flow", "PE");
    legend->Draw();

    //canvas->SaveAs("mott_yield_600flow_800flow_live_charge.png");
}
