
// Overlay two live-charge-normalized Mott yields.
//
// The input histograms were already normalized by live charge when they were
// produced in background.C, so no additional scaling is applied here.
//
// Usage:
//   root -l 'emptyTarget/compare_flowRate.C'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

// Fill in the A/B/C input file names for each flow setting here.
static const char *fileA600Flow = "A_600cc_vertex.root";
static const char *fileB600Flow = "B_600cc_vertex.root";
static const char *fileC600Flow = "C_600cc_vertex.root";

static const char *fileA1300Flow = "A_1300cc_vertex.root";
static const char *fileB1300Flow = "B_1300cc_vertex.root";
static const char *fileC1300Flow = "C_1300cc_vertex.root";

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

static TH1 *loadHistogram(const char *fileName, const char *histName,
                          const char *cloneName)
{
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << fileName << std::endl;
        delete file;
        return nullptr;
    }

    TH1 *hist = dynamic_cast<TH1 *>(file->Get(histName));
    if (!hist) {
        std::cerr << histName << " not found in " << fileName << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    TH1 *clone = dynamic_cast<TH1 *>(hist->Clone(cloneName));
    clone->SetDirectory(nullptr);
    if (clone->GetSumw2N() == 0) clone->Sumw2();

    file->Close();
    delete file;
    std::cout << "Loaded " << histName << " from " << fileName << std::endl;
    return clone;
}

static void drawNormalizedMollerVertexAComparison()
{
    TH1 *vertex600 = loadHistogram(
        fileA600Flow, "moller_vertex_z", "moller_vertex_z_A_600flow");
    TH1 *vertex1300 = loadHistogram(
        fileA1300Flow, "moller_vertex_z", "moller_vertex_z_A_1300flow");
    if (!vertex600 || !vertex1300) return;

    const double maximum600 = vertex600->GetMaximum();
    const double maximum1300 = vertex1300->GetMaximum();
    if (maximum600 <= 0. || maximum1300 <= 0.) {
        std::cerr << "Cannot normalize moller_vertex_z: histogram maximum "
                     "is not positive"
                  << std::endl;
        return;
    }

    vertex600->Scale(1. / maximum600);
    vertex1300->Scale(1. / maximum1300);

    vertex600->SetTitle(
        "Type A Moller Vertex Z: 600flow vs 1300flow;"
        "Moller Vertex Z [mm];Normalized counts");
    vertex600->SetLineColor(kBlue + 1);
    vertex600->SetLineWidth(2);
    vertex600->SetStats(0);

    vertex1300->SetLineColor(kRed + 1);
    vertex1300->SetLineWidth(2);
    vertex1300->SetStats(0);

    vertex600->SetMinimum(0.);
    vertex600->SetMaximum(1.15);
    vertex600->GetXaxis()->CenterTitle();
    vertex600->GetYaxis()->CenterTitle();
    vertex600->GetXaxis()->SetTitleSize(0.045);
    vertex600->GetYaxis()->SetTitleSize(0.045);
    vertex600->GetYaxis()->SetTitleOffset(1.15);

    TCanvas *canvas = new TCanvas(
        "c_moller_vertex_z_A_flow_compare",
        "Type A Moller Vertex Z Comparison",
        1000, 700);
    canvas->SetGrid();
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.12);

    vertex600->Draw("HIST");
    vertex1300->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.66, 0.75, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(vertex600, "A, 600flow", "l");
    legend->AddEntry(vertex1300, "A, 1300flow", "l");
    legend->Draw();

    canvas->SaveAs("moller_vertex_z_A_600flow_1300flow.png");
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
           axis->GetBinUpEdge(firstGroupedBin) <= theta + 1e-9) {
        edges.push_back(axis->GetBinUpEdge(firstGroupedBin));
        firstGroupedBin++;
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

static void fixBoverAmB(TH1 *ratio, TH1 *histA, TH1 *histB)
{
    for (int bin = 1; bin <= ratio->GetNbinsX(); bin++) {
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

static void fixBmCoverAmB(TH1 *ratio, TH1 *histA, TH1 *histB, TH1 *histC)
{
    for (int bin = 1; bin <= ratio->GetNbinsX(); bin++) {
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

static bool drawMottBackgroundRatio(TPad *pad,
                                    const char *flowTag,
                                    const char *flowLabel,
                                    const char *fileA,
                                    const char *fileB,
                                    const char *fileC)
{
    TH1 *histA = loadNormalizedMottYield(fileA, Form("%s_A", flowTag));
    TH1 *histB = loadNormalizedMottYield(fileB, Form("%s_B", flowTag));
    TH1 *histC = loadNormalizedMottYield(fileC, Form("%s_C", flowTag));
    if (!histA || !histB || !histC) return false;

    // Keep the original binning through 2 degrees. Above 2 degrees, merge
    // every three adjacent bins before constructing the signal and ratios.
    histA = rebinAboveTheta(
        histA, 2.0, 3, Form("mott_%s_A_rebin3_above2deg", flowTag));
    histB = rebinAboveTheta(
        histB, 2.0, 3, Form("mott_%s_B_rebin3_above2deg", flowTag));
    histC = rebinAboveTheta(
        histC, 2.0, 3, Form("mott_%s_C_rebin3_above2deg", flowTag));

    TH1 *signal = dynamic_cast<TH1 *>(
        histA->Clone(Form("mott_signal_%s_AminusB", flowTag)));
    signal->SetDirectory(nullptr);
    signal->Add(histB, -1.);

    TH1 *ratioB = dynamic_cast<TH1 *>(
        histB->Clone(Form("mott_ratio_%s_B_over_AmB", flowTag)));
    TH1 *ratioC = dynamic_cast<TH1 *>(
        histC->Clone(Form("mott_ratio_%s_C_over_AmB", flowTag)));
    TH1 *bMinusC = dynamic_cast<TH1 *>(
        histB->Clone(Form("mott_%s_BminusC", flowTag)));
    bMinusC->Add(histC, -1.);
    TH1 *ratioBmC = dynamic_cast<TH1 *>(
        bMinusC->Clone(Form("mott_ratio_%s_BmC_over_AmB", flowTag)));

    ratioB->SetDirectory(nullptr);
    ratioC->SetDirectory(nullptr);
    ratioBmC->SetDirectory(nullptr);
    ratioB->Divide(histB, signal, 1., 1., "");
    ratioC->Divide(histC, signal, 1., 1., "");
    ratioBmC->Divide(bMinusC, signal, 1., 1., "");

    fixBoverAmB(ratioB, histA, histB);
    fixBmCoverAmB(ratioBmC, histA, histB, histC);

    TH1 *ratios[3] = {ratioB, ratioC, ratioBmC};
    const int colors[3] = {kRed, kGreen + 2, kBlue + 1};
    const int markers[3] = {20, 21, 24};
    const char *labels[3] = {
        "B/(A-B)", "C/(A-B)", "(B-C)/(A-B)"
    };

    for (int i = 0; i < 3; i++) {
        ratios[i]->Scale(100.);
        ratios[i]->SetLineColor(colors[i]);
        ratios[i]->SetMarkerColor(colors[i]);
        ratios[i]->SetMarkerStyle(markers[i]);
        ratios[i]->SetLineWidth(2);
        ratios[i]->SetStats(0);
    }

    ratioB->SetTitle(
        Form("%s Mott Background Ratio;"
             "Reconstructed Scattering Angle [deg];"
             "N_{background} / N_{signal} (%%)",
             flowLabel));
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

    pad->cd();
    pad->SetLogx();
    pad->SetGrid();
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.04);
    pad->SetBottomMargin(0.12);

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
        Form("fit_mott_BmC_over_AmB_%s", flowTag),
        "[0]+[1]*log(x)", 0.5, 3.8);
    fit->SetLineColor(kBlack);
    fit->SetLineStyle(2);
    fit->SetLineWidth(3);
    ratioBmC->Fit(fit, "RQ0");
    fit->Draw("SAME");

    TLegend *legend = new TLegend(0.57, 0.67, 0.90, 0.90);
    legend->SetBorderSize(0);
    for (int i = 0; i < 3; i++)
        legend->AddEntry(ratios[i], labels[i], "PE");
    legend->AddEntry(
        fit,
        Form("fit: %.3f %+.3f ln#theta",
             fit->GetParameter(0), fit->GetParameter(1)),
        "l");
    legend->Draw();

    return true;
}

void compare_flowRate()
{
    gStyle->SetOptStat(0);

    TH1 *h600Flow = loadNormalizedMottYield(fileC600Flow, "600flow_C");
    TH1 *h1300Flow = loadNormalizedMottYield(fileC1300Flow, "1300flow_C");
    if (!h600Flow || !h1300Flow) return;

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

    h1300Flow->SetLineColor(kRed + 1);
    h1300Flow->SetMarkerColor(kRed + 1);
    h1300Flow->SetMarkerStyle(21);
    h1300Flow->SetMarkerSize(0.8);
    h1300Flow->SetLineWidth(2);

    h600Flow->GetXaxis()->SetRangeUser(thetaMin, thetaMax);
    h600Flow->GetXaxis()->CenterTitle();
    h600Flow->GetYaxis()->CenterTitle();
    h600Flow->GetXaxis()->SetTitleSize(0.045);
    h600Flow->GetYaxis()->SetTitleSize(0.045);
    h600Flow->GetYaxis()->SetTitleOffset(1.25);

    double maximum = std::max(h600Flow->GetMaximum(), h1300Flow->GetMaximum());
    h600Flow->SetMinimum(0.);
    h600Flow->SetMaximum(1.20 * maximum);

    TCanvas *canvas = new TCanvas(
        "c_compare_flowRate", "Mott Yield Comparison", 1000, 700);
    canvas->SetGrid();
    canvas->SetLeftMargin(0.13);
    canvas->SetBottomMargin(0.12);

    h600Flow->Draw("E1");
    h1300Flow->Draw("E1 SAME");

    TLegend *legend = new TLegend(0.67, 0.74, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h600Flow, "600flow", "PE");
    legend->AddEntry(h1300Flow, "1300flow", "PE");
    legend->Draw();

    //canvas->SaveAs("mott_yield_600flow_800flow_live_charge.png");

    TCanvas *backgroundCanvas = new TCanvas(
        "c_mott_background_flow_compare",
        "Mott Background Ratio: 600flow vs 1300flow",
        1200, 600);
    backgroundCanvas->Divide(2, 1);

    const bool drew600 = drawMottBackgroundRatio(
        static_cast<TPad *>(backgroundCanvas->cd(1)),
        "600flow", "600flow",
        fileA600Flow, fileB600Flow, fileC600Flow);
    const bool drew1300 = drawMottBackgroundRatio(
        static_cast<TPad *>(backgroundCanvas->cd(2)),
        "1300flow", "1300flow",
        fileA1300Flow, fileB1300Flow, fileC1300Flow);

    if (drew600 && drew1300)
        backgroundCanvas->SaveAs(
            "mott_background_ratio_600flow_1300flow.png");

    drawNormalizedMollerVertexAComparison();
}
