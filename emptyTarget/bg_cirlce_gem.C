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

static const double Ebeam = 2239.51;  // MeV
static const bool applyVertexCut = false;

struct CircleYields {
    TH1F *mott = nullptr;
    TH1F *moller = nullptr;
};

static CircleYields fillNormalizedYields(const char *fileName,
                                         double liveCharge,
                                         const char *tag)
{
    CircleYields yields;
    if (liveCharge <= 0.) {
        std::cerr << "Invalid live charge for " << tag << ": "
                  << liveCharge << std::endl;
        return yields;
    }

    TChain tree("recon");
    if (tree.Add(fileName) == 0 || tree.GetEntries() == 0) {
        std::cerr << "Cannot read recon entries from " << fileName << std::endl;
        return yields;
    }

    TH1F *mottYield = new TH1F(
        Form("mott_yield_%s", tag),
        Form("Circle %s Mott Yield;Reconstructed Scattering Angle [deg];"
             "Count / live charge", tag),
        120, 0., 6.);
    mottYield->SetDirectory(nullptr);
    mottYield->Sumw2();

    TH1F *mollerYield = new TH1F(
        Form("moller_yield_%s", tag),
        Form("Circle %s Moller Yield;Reconstructed Scattering Angle [deg];"
             "Count / live charge", tag),
        120, 0., 6.);
    mollerYield->SetDirectory(nullptr);
    mollerYield->Sumw2();

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
        std::vector<HCHit> moller_Hits_candidate;
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
                isMott(energy, Ebeam, resolution))
                mottYield->Fill(theta);

            if (energy > 80. &&
                energy < Ebeam - 2. * resolution * Ebeam /
                             std::sqrt(Ebeam / 1000.))
                moller_Hits_candidate.emplace_back(x, y, 6270.f, energy);

        }
        //select Moller event and fill moller hist
        std::sort(moller_Hits_candidate.begin(),
                  moller_Hits_candidate.end(),
                  [](const HCHit &a, const HCHit &b) {
                      return a.energy > b.energy;
                  });

        MollerData mollerData_event;
        const int nCand = moller_Hits_candidate.size();
        for (int ii = 0; ii < nCand; ++ii) {
            const HCHit &hi = moller_Hits_candidate[ii];
            const float theta_i =
                std::atan2(std::sqrt(hi.x * hi.x + hi.y * hi.y), hi.z) *
                180.f / M_PI;
            for (int jj = nCand - 1; jj > ii; --jj) {
                const HCHit &hj = moller_Hits_candidate[jj];
                const float theta_j =
                    std::atan2(std::sqrt(hj.x * hj.x + hj.y * hj.y),
                               hj.z) *
                    180.f / M_PI;
                if (isMoller_kinematic(theta_i, hi.energy,
                                       theta_j, hj.energy,
                                       Ebeam, resolution)) {
                    MollerEvent candidate = {
                        DataPoint(hi.x, hi.y, hi.z, hi.energy),
                        DataPoint(hj.x, hj.y, hj.z, hj.energy)
                    };
                    if (std::fabs(GetMollerPhiDiff(candidate)) > 10.f)
                        continue;
                    mollerData_event.emplace_back(candidate);
                }
            }
        }

        if (mollerData_event.size() == 0) continue;

        if (mollerData_event.size() > 1) {
            auto getPt = [](const MollerEvent &mev) -> float {
                const float sin_t1 =
                    std::sqrt(mev.first.x * mev.first.x +
                              mev.first.y * mev.first.y) /
                    mev.first.z;
                const float sin_t2 =
                    std::sqrt(mev.second.x * mev.second.x +
                              mev.second.y * mev.second.y) /
                    mev.second.z;
                return std::fabs(mev.first.E * sin_t1 -
                                 mev.second.E * sin_t2);
            };
            auto best = std::min_element(
                mollerData_event.begin(), mollerData_event.end(),
                [&](const MollerEvent &a, const MollerEvent &b) {
                    return getPt(a) < getPt(b);
                });
            MollerEvent bestPair = *best;
            mollerData_event.clear();
            mollerData_event.push_back(bestPair);
        }

        const MollerEvent &mev = mollerData_event.front();
        const float t1 =
            std::atan2(std::sqrt(mev.first.x * mev.first.x +
                                 mev.first.y * mev.first.y),
                       mev.first.z) *
            180.f / M_PI;
        const float t2 =
            std::atan2(std::sqrt(mev.second.x * mev.second.x +
                                 mev.second.y * mev.second.y),
                       mev.second.z) *
            180.f / M_PI;
        mollerYield->Fill(t1);
        mollerYield->Fill(t2);
    
    }
    std::cout << "  " << entries << " / " << entries << std::endl;

    mottYield->Scale(1. / liveCharge);
    mollerYield->Scale(1. / liveCharge);
    std::cout << "  normalized Mott integral = " << mottYield->Integral()
              << std::endl;
    std::cout << "  normalized Moller integral = " << mollerYield->Integral()
              << std::endl;
    yields.mott = mottYield;
    yields.moller = mollerYield;
    return yields;
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

    CircleYields yieldsA = fillNormalizedYields(fileA, liveChargeA, "A");
    CircleYields yieldsB = fillNormalizedYields(fileB, liveChargeB, "B");
    CircleYields yieldsC = fillNormalizedYields(fileC, liveChargeC, "C");
    if (!yieldsA.mott || !yieldsB.mott || !yieldsC.mott ||
        !yieldsA.moller || !yieldsB.moller || !yieldsC.moller)
        return;

    TH1F *histAOriginal = yieldsA.mott;
    TH1F *histBOriginal = yieldsB.mott;
    TH1F *histCOriginal = yieldsC.mott;
    TH1F *mollerAOriginal = yieldsA.moller;
    TH1F *mollerBOriginal = yieldsB.moller;
    TH1F *mollerCOriginal = yieldsC.moller;

    // Match compare_flowRate.C: retain the original bins through 2 degrees,
    // then combine every three bins above 2 degrees.
    TH1 *histA = rebinAboveTheta(
        histAOriginal, 2.0, 3, "mott_yield_A_rebin3_above2deg");
    TH1 *histB = rebinAboveTheta(
        histBOriginal, 2.0, 3, "mott_yield_B_rebin3_above2deg");
    TH1 *histC = rebinAboveTheta(
        histCOriginal, 2.0, 3, "mott_yield_C_rebin3_above2deg");
    TH1 *mollerA = rebinAboveTheta(
        mollerAOriginal, 2.0, 3, "moller_yield_A_rebin3_above2deg");
    TH1 *mollerB = rebinAboveTheta(
        mollerBOriginal, 2.0, 3, "moller_yield_B_rebin3_above2deg");
    TH1 *mollerC = rebinAboveTheta(
        mollerCOriginal, 2.0, 3, "moller_yield_C_rebin3_above2deg");

    TH1 *signal = dynamic_cast<TH1 *>(histA->Clone("mott_signal_AminusB"));
    signal->SetDirectory(nullptr);
    signal->Add(histB, -1.);

    TH1 *mollerSignal =
        dynamic_cast<TH1 *>(mollerA->Clone("moller_signal_AminusB"));
    mollerSignal->SetDirectory(nullptr);
    mollerSignal->Add(mollerB, -1.);

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
    mollerAOriginal->Write();
    mollerBOriginal->Write();
    mollerCOriginal->Write();
    histA->Write();
    histB->Write();
    histC->Write();
    mollerA->Write();
    mollerB->Write();
    mollerC->Write();
    signal->Write();
    mollerSignal->Write();
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
