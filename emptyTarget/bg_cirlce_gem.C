// Generate the circle A/B/C GEM Mott-background comparison directly from
// three recon ROOT files.
//
// Edit only the configuration block below, then run:
//   root -l -q emptyTarget/bg_cirlce_gem.C

#include "../EventData.h"
#include "../PhysicsTools.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

// ── Configuration: change the A/B/C file names and live charges here ─────
static const char *fileA = "../data/empty_target/circle_A.filtered.root";
static const char *fileB = "../data/empty_target/circle_B.filtered.root";
static const char *fileC = "../data/empty_target/circle_C.filtered.root";

static const double liveChargeA = 1.0;
static const double liveChargeB = 1.0;
static const double liveChargeC = 1.0;

static const char *outputRoot = "circle_gem_background.root";
static const char *outputPng  = "circle_gem_background.png";

static const double beamEnergy = 2239.51;  // MeV
static const bool applyVertexCut = false;

static TH1F *fillNormalizedMottYield(const char *fileName,
                                     double liveCharge,
                                     const char *tag)
{
    if (liveCharge <= 0.) {
        std::cerr << "Invalid live charge for " << tag << ": "
                  << liveCharge << std::endl;
        return nullptr;
    }

    TChain tree("recon");
    if (tree.Add(fileName) == 0 || tree.GetEntries() == 0) {
        std::cerr << "Cannot read recon entries from " << fileName << std::endl;
        return nullptr;
    }

    TH1F *yield = new TH1F(
        Form("mott_yield_%s", tag),
        Form("Circle %s Mott Yield;Reconstructed Scattering Angle [deg];"
             "Count / live charge", tag),
        120, 0., 6.);
    yield->SetDirectory(nullptr);
    yield->Sumw2();

    ReconEventData ev;
    setupReconBranches(&tree, ev);

    const float resolution = 0.037f;
    const float angleEdge[4] = {0.8f, 1.2f, 2.0f, 3.0f};
    const float vertexResolution[5] = {417.f, 252.f, 165.f, 101.f, 76.f};
    const Long64_t entries = tree.GetEntries();

    std::cout << "Processing " << tag << ": " << fileName
              << " (" << entries << " entries), live charge = "
              << liveCharge << std::endl;

    for (Long64_t i = 0; i < entries; ++i) {
        if (i % 10000 == 0)
            std::cout << "  " << i << " / " << entries << "\r" << std::flush;

        tree.GetEntry(i);
        for (int j = 0; j < ev.matchNum; ++j) {
            const float x1 = ev.mHit_gx[j][1];
            const float y1 = ev.mHit_gy[j][1];
            const float z1 = ev.mHit_gz[j][1];
            const float x2 = ev.mHit_gx[j][0];
            const float y2 = ev.mHit_gy[j][0];
            const float z2 = ev.mHit_gz[j][0];

            if (std::abs(z1) < 1.e-12f || std::abs(z2) < 1.e-12f)
                continue;

            const float x = x1 * 6270.f / z1;
            const float y = y1 * 6270.f / z1;
            const float energy = ev.mHit_E[j];
            const float theta =
                std::atan2(std::sqrt(x * x + y * y), 6270.f) *
                180.f / M_PI;

            if (std::abs(x) < 20.25f * 2.5f &&
                std::abs(y) < 20.25f * 2.5f)
                continue;
            if (std::abs(x) > 20.25f * 16.f ||
                std::abs(y) > 20.25f * 16.f)
                continue;

            const float dx = x2 - x1;
            const float dy = y2 - y1;
            const float dz = z2 - z1;
            const float denominator = dx * dx + dy * dy;
            if (denominator < 1.e-12f) continue;

            const float t = -(x1 * dx + y1 * dy) / denominator;
            const float vertexZ = z1 + t * dz;

            const float gem1x = 6270.f * x1 / z1;
            const float gem1y = 6270.f * y1 / z1;
            const float gem2x = 6270.f * x2 / z2;
            const float gem2y = 6270.f * y2 / z2;
            if (std::abs(gem1x - gem2x) > 50.f ||
                std::abs(gem1y - gem2y) > 50.f)
                continue;

            int resolutionBin = 4;
            if (theta < angleEdge[0]) resolutionBin = 0;
            else if (theta < angleEdge[1]) resolutionBin = 1;
            else if (theta < angleEdge[2]) resolutionBin = 2;
            else if (theta < angleEdge[3]) resolutionBin = 3;

            const bool vertexIn =
                !applyVertexCut ||
                std::abs(vertexZ) <= 3.f * vertexResolution[resolutionBin];

            if (vertexIn &&
                isMott(energy, beamEnergy, resolution))
                yield->Fill(theta);
        }
    }
    std::cout << "  " << entries << " / " << entries << std::endl;

    yield->Scale(1. / liveCharge);
    std::cout << "  normalized Mott integral = " << yield->Integral()
              << std::endl;
    return yield;
}

static TH1 *rebinAboveTheta(TH1 *hist, double theta, int group,
                            const char *newName)
{
    const TAxis *axis = hist->GetXaxis();
    const int nBins = axis->GetNbins();
    std::vector<double> edges;
    edges.push_back(axis->GetBinLowEdge(1));

    int firstGroupedBin = 1;
    while (firstGroupedBin <= nBins &&
           axis->GetBinUpEdge(firstGroupedBin) <= theta + 1.e-9) {
        edges.push_back(axis->GetBinUpEdge(firstGroupedBin));
        ++firstGroupedBin;
    }

    for (int bin = firstGroupedBin; bin <= nBins; bin += group) {
        const int lastBin = std::min(bin + group - 1, nBins);
        edges.push_back(axis->GetBinUpEdge(lastBin));
    }

    TH1 *rebinned = hist->Rebin(
        static_cast<int>(edges.size()) - 1, newName, edges.data());
    rebinned->SetDirectory(nullptr);
    return rebinned;
}

static void setBoverAmBErrors(TH1 *ratio, const TH1 *histA,
                             const TH1 *histB)
{
    for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
        const double A = histA->GetBinContent(bin);
        const double B = histB->GetBinContent(bin);
        const double sigmaA = histA->GetBinError(bin);
        const double sigmaB = histB->GetBinError(bin);
        const double denominator = A - B;

        if (denominator == 0.) {
            ratio->SetBinContent(bin, 0.);
            ratio->SetBinError(bin, 0.);
            continue;
        }

        ratio->SetBinError(
            bin,
            std::sqrt(B * B * sigmaA * sigmaA +
                      A * A * sigmaB * sigmaB) /
                (denominator * denominator));
    }
}

static void setBmCoverAmBErrors(TH1 *ratio, const TH1 *histA,
                               const TH1 *histB, const TH1 *histC)
{
    for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
        const double A = histA->GetBinContent(bin);
        const double B = histB->GetBinContent(bin);
        const double C = histC->GetBinContent(bin);
        const double sigmaA = histA->GetBinError(bin);
        const double sigmaB = histB->GetBinError(bin);
        const double sigmaC = histC->GetBinError(bin);
        const double denominator = A - B;

        if (denominator == 0.) {
            ratio->SetBinContent(bin, 0.);
            ratio->SetBinError(bin, 0.);
            continue;
        }

        ratio->SetBinError(
            bin,
            std::sqrt((B - C) * (B - C) * sigmaA * sigmaA +
                      (A - C) * (A - C) * sigmaB * sigmaB +
                      denominator * denominator * sigmaC * sigmaC) /
                (denominator * denominator));
    }
}

void bg_cirlce_gem()
{
    gStyle->SetOptStat(0);

    TH1F *histAOriginal =
        fillNormalizedMottYield(fileA, liveChargeA, "A");
    TH1F *histBOriginal =
        fillNormalizedMottYield(fileB, liveChargeB, "B");
    TH1F *histCOriginal =
        fillNormalizedMottYield(fileC, liveChargeC, "C");
    if (!histAOriginal || !histBOriginal || !histCOriginal) return;

    // Match compare_flowRate.C: retain the original bins through 2 degrees,
    // then combine every three bins above 2 degrees.
    TH1 *histA = rebinAboveTheta(
        histAOriginal, 2.0, 3, "mott_yield_A_rebin3_above2deg");
    TH1 *histB = rebinAboveTheta(
        histBOriginal, 2.0, 3, "mott_yield_B_rebin3_above2deg");
    TH1 *histC = rebinAboveTheta(
        histCOriginal, 2.0, 3, "mott_yield_C_rebin3_above2deg");

    TH1 *signal = dynamic_cast<TH1 *>(histA->Clone("mott_signal_AminusB"));
    signal->SetDirectory(nullptr);
    signal->Add(histB, -1.);

    TH1 *bMinusC = dynamic_cast<TH1 *>(histB->Clone("mott_background_BminusC"));
    bMinusC->SetDirectory(nullptr);
    bMinusC->Add(histC, -1.);

    TH1 *ratioB = dynamic_cast<TH1 *>(histB->Clone("mott_ratio_B_over_AmB"));
    TH1 *ratioC = dynamic_cast<TH1 *>(histC->Clone("mott_ratio_C_over_AmB"));
    TH1 *ratioBmC =
        dynamic_cast<TH1 *>(bMinusC->Clone("mott_ratio_BmC_over_AmB"));
    ratioB->SetDirectory(nullptr);
    ratioC->SetDirectory(nullptr);
    ratioBmC->SetDirectory(nullptr);

    ratioB->Divide(histB, signal);
    ratioC->Divide(histC, signal);
    ratioBmC->Divide(bMinusC, signal);
    setBoverAmBErrors(ratioB, histA, histB);
    setBmCoverAmBErrors(ratioBmC, histA, histB, histC);

    TH1 *ratios[3] = {ratioB, ratioC, ratioBmC};
    const int colors[3] = {kRed, kGreen + 2, kBlue + 1};
    const int markers[3] = {20, 21, 24};
    const char *labels[3] = {
        "B/(A-B)", "C/(A-B)", "(B-C)/(A-B)"
    };

    for (int i = 0; i < 3; ++i) {
        ratios[i]->Scale(100.);
        ratios[i]->SetLineColor(colors[i]);
        ratios[i]->SetMarkerColor(colors[i]);
        ratios[i]->SetMarkerStyle(markers[i]);
        ratios[i]->SetLineWidth(2);
        ratios[i]->SetStats(0);
    }

    ratioB->SetTitle(
        "Circle GEM Mott Background Ratio;"
        "Reconstructed Scattering Angle [deg];"
        "N_{background} / N_{signal} (%)");
    ratioB->GetXaxis()->SetRangeUser(0.5, 3.8);
    ratioB->GetYaxis()->SetRangeUser(-0.2, 2.0);
    ratioB->GetXaxis()->CenterTitle();
    ratioB->GetYaxis()->CenterTitle();
    ratioB->GetXaxis()->SetTitleSize(0.045);
    ratioB->GetYaxis()->SetTitleSize(0.045);
    ratioB->GetXaxis()->SetTitleOffset(0.9);
    ratioB->GetYaxis()->SetTitleOffset(0.9);
    ratioB->GetXaxis()->SetMoreLogLabels();
    ratioB->GetXaxis()->SetNoExponent();

    TCanvas *canvas = new TCanvas(
        "c_circle_gem_background", "Circle GEM Mott Background Ratio",
        1000, 700);
    canvas->SetLogx();
    canvas->SetGrid();
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.04);
    canvas->SetBottomMargin(0.12);

    ratioB->Draw("E");
    ratioC->Draw("E SAME");
    ratioBmC->Draw("E SAME");

    TLine *lineOne = new TLine(0.5, 1., 3.8, 1.);
    lineOne->SetLineStyle(2);
    lineOne->SetLineColor(kGray + 1);
    lineOne->Draw("SAME");

    TLine *lineZero = new TLine(0.5, 0., 3.8, 0.);
    lineZero->SetLineStyle(3);
    lineZero->SetLineColor(kGray + 1);
    lineZero->Draw("SAME");

    TF1 *fit = new TF1(
        "fit_mott_BmC_over_AmB", "[0]+[1]*log(x)", 0.5, 3.8);
    fit->SetLineColor(kBlack);
    fit->SetLineStyle(2);
    fit->SetLineWidth(3);
    ratioBmC->Fit(fit, "RQ0");
    fit->Draw("SAME");

    TLegend *legend = new TLegend(0.57, 0.67, 0.90, 0.90);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    for (int i = 0; i < 3; ++i)
        legend->AddEntry(ratios[i], labels[i], "PE");
    legend->AddEntry(
        fit,
        Form("fit: %.3f %+.3f ln#theta",
             fit->GetParameter(0), fit->GetParameter(1)),
        "l");
    legend->Draw();

    canvas->SaveAs(outputPng);

    TFile output(outputRoot, "RECREATE");
    histAOriginal->Write();
    histBOriginal->Write();
    histCOriginal->Write();
    histA->Write();
    histB->Write();
    histC->Write();
    signal->Write();
    bMinusC->Write();
    ratioB->Write();
    ratioC->Write();
    ratioBmC->Write();
    fit->Write();
    canvas->Write();
    TParameter<double>("live_charge_A", liveChargeA).Write();
    TParameter<double>("live_charge_B", liveChargeB).Write();
    TParameter<double>("live_charge_C", liveChargeC).Write();
    output.Close();

    std::cout << "Saved plot to " << outputPng << std::endl;
    std::cout << "Saved histograms and canvas to " << outputRoot << std::endl;
}
