
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TList.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"

static const char *circleLabels[] = {
    "25360", "25356", "25340", "25282", "25335", "25329", "25318"
};
static const char *referenceRun = "25340";

static double fitStatError(const TF1 *fit, double theta)
{
    // For f(theta) = p0 + p1 * log(theta). The stored TF1 contains the
    // fitted parameter errors, but not their covariance.
    const double logTheta = std::log(theta);
    const double error0 = fit->GetParError(0);
    const double error1 = fit->GetParError(1);
    return std::sqrt(error0 * error0 +
                     logTheta * logTheta * error1 * error1);
}

static double powerLawFitStatError(const TF1 *fit, double theta)
{
    // For f(theta) = A * theta^p, ignoring the unavailable covariance term.
    const double amplitude = fit->GetParameter(0);
    const double exponent = fit->GetParameter(1);
    const double errorAmplitude = fit->GetParError(0);
    const double errorExponent = fit->GetParError(1);
    const double thetaPower = std::pow(theta, exponent);
    const double derivativeAmplitude = thetaPower;
    const double derivativeExponent =
        amplitude * thetaPower * std::log(theta);
    return std::sqrt(
        std::pow(derivativeAmplitude * errorAmplitude, 2) +
        std::pow(derivativeExponent * errorExponent, 2));
}

static TGraphErrors *histToGraphErrors(const TH1 *hist,
                                       const char *name,
                                       bool keepYErrors)
{
    TGraphErrors *graph = new TGraphErrors(hist->GetNbinsX());
    graph->SetName(name);

    int point = 0;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        const double theta = hist->GetBinCenter(bin);
        if (theta < 0.5 || theta > 3.8) continue;

        graph->SetPoint(point, theta, hist->GetBinContent(bin));
        graph->SetPointError(point, 0.,
                             keepYErrors ? hist->GetBinError(bin) : 0.);
        ++point;
    }

    graph->Set(point);
    return graph;
}

void circle_compare(const char *folder = "./",
                    const char *output = "circle_compare.png",
                    const char *fitRatioOutput = "circle_fit_ratio.png",
                    const char *bFitRatioOutput = "circle_B_fit_ratio.png",
                    const char *mottBYieldFitRatioOutput =
                        "circle_mott_yield_B_fit_ratio.png",
                    const char *mottSignalRatioOutput =
                        "circle_mott_signal_AmB_ratio.png")
{
    // ROOT's normal X11 canvas can ignore alpha transparency for filled
    // histograms.  Prefer the GL painter when available, and use a
    // non-solid fill style so the uncertainty band never becomes an opaque
    // block even on older ROOT backends.
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetHatchesLineWidth(2);
    gStyle->SetHatchesSpacing(0.45);

    std::vector<TString> fileNames;

    TSystemDirectory directory("circle_root_directory", folder);
    TList *files = directory.GetListOfFiles();
    if (!files) {
        std::cerr << "Cannot open directory: " << folder << std::endl;
        return;
    }

    TIter next(files);
    while (TSystemFile *entry = dynamic_cast<TSystemFile *>(next())) {
        const TString name = entry->GetName();
        if (entry->IsDirectory()) continue;
        if (!name.BeginsWith("circle_") || !name.EndsWith(".root")) continue;
        fileNames.push_back(name);
    }

    std::sort(fileNames.begin(), fileNames.end());

    if (fileNames.empty()) {
        std::cerr << "No circle_*.root files found in " << folder << std::endl;
        return;
    }

    if (fileNames.size() > 8) {
        std::cout << "Found " << fileNames.size()
                  << " files; only the first 8 will be drawn." << std::endl;
        fileNames.resize(8);
    }

    // First load and clone every complete canvas into memory.  The clones no
    // longer depend on the input ROOT files after those files are closed.
    std::vector<TCanvas *> sourceCanvases;
    std::vector<TString> sourceNames;
    std::vector<TString> sourceLabels;
    std::vector<TH1 *> sourceMottBYields;
    std::vector<TH1 *> sourceMottSignals;

    for (const TString &name : fileNames) {
        const TString path = TString::Format("%s/%s", folder, name.Data());
        TFile input(path, "READ");
        if (input.IsZombie()) {
            std::cerr << "Cannot open " << path << std::endl;
            continue;
        }

        TCanvas *source =
            dynamic_cast<TCanvas *>(input.Get("c_circle_gem_background"));
        if (!source) {
            std::cerr << "c_circle_gem_background not found in "
                      << path << std::endl;
            continue;
        }

        TCanvas *clone = dynamic_cast<TCanvas *>(
            source->Clone(Form("c_circle_gem_background_%zu",
                               sourceCanvases.size() + 1)));
        if (!clone) {
            std::cerr << "Cannot clone canvas from " << path << std::endl;
            continue;
        }

        TH1 *mottBYieldIn =
            dynamic_cast<TH1 *>(input.Get("mott_yield_B"));
        TH1 *mottBYield = nullptr;
        if (mottBYieldIn) {
            mottBYield = dynamic_cast<TH1 *>(
                mottBYieldIn->Clone(
                    Form("mott_yield_B_source_%zu",
                         sourceCanvases.size() + 1)));
            mottBYield->SetDirectory(nullptr);
        } else {
            std::cerr << "mott_yield_B not found in " << path << std::endl;
        }

        TH1 *mottSignalIn =
            dynamic_cast<TH1 *>(input.Get("mott_signal_AminusB"));
        TH1 *mottSignal = nullptr;
        if (mottSignalIn) {
            mottSignal = dynamic_cast<TH1 *>(
                mottSignalIn->Clone(
                    Form("mott_signal_AminusB_source_%zu",
                         sourceCanvases.size() + 1)));
            mottSignal->SetDirectory(nullptr);
        } else {
            std::cerr << "mott_signal_AminusB not found in "
                      << path << std::endl;
        }

        sourceCanvases.push_back(clone);
        sourceNames.push_back(name);
        sourceMottBYields.push_back(mottBYield);
        sourceMottSignals.push_back(mottSignal);
        TString circleNumber = name;
        circleNumber.ReplaceAll("circle_", "");
        circleNumber.ReplaceAll(".root", "");
        const int circleIndex = circleNumber.Atoi() - 1;
        if (circleIndex >= 0 &&
            circleIndex < static_cast<int>(
                sizeof(circleLabels) / sizeof(circleLabels[0]))) {
            sourceLabels.emplace_back(circleLabels[circleIndex]);
        } else {
            sourceLabels.push_back(name);
        }
        std::cout << "Loaded complete canvas from " << path << std::endl;
    }

    if (sourceCanvases.empty()) {
        std::cerr << "No c_circle_gem_background canvases were loaded."
                  << std::endl;
        return;
    }

    // Arrange runs numerically from small to large.
    std::vector<std::size_t> order(sourceCanvases.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) {
                  return sourceLabels[a].Atoi() < sourceLabels[b].Atoi();
              });

    std::vector<TCanvas *> sortedCanvases;
    std::vector<TString> sortedNames;
    std::vector<TString> sortedLabels;
    std::vector<TH1 *> sortedMottBYields;
    std::vector<TH1 *> sortedMottSignals;
    for (std::size_t index : order) {
        sortedCanvases.push_back(sourceCanvases[index]);
        sortedNames.push_back(sourceNames[index]);
        sortedLabels.push_back(sourceLabels[index]);
        sortedMottBYields.push_back(sourceMottBYields[index]);
        sortedMottSignals.push_back(sourceMottSignals[index]);
    }
    sourceCanvases.swap(sortedCanvases);
    sourceNames.swap(sortedNames);
    sourceLabels.swap(sortedLabels);
    sourceMottBYields.swap(sortedMottBYields);
    sourceMottSignals.swap(sortedMottSignals);

    // Only after all source canvases have been loaded do we create and draw
    // the combined 2 x 4 canvas.
    TCanvas *comparison = new TCanvas(
        "c_circle_compare",
        "Circle GEM Background Comparison",
        1800, 900);
    comparison->Divide(4, 2, 0.002, 0.002);

    for (std::size_t i = 0; i < sourceCanvases.size(); ++i) {
        TPad *pad =
            dynamic_cast<TPad *>(comparison->cd(static_cast<int>(i) + 1));
        sourceCanvases[i]->DrawClonePad();
        pad->cd();

        TLatex label;
        label.SetNDC();
        label.SetTextFont(62);
        label.SetTextSize(0.04);
        label.DrawLatex(0.15, 0.87, sourceLabels[i]);

        pad->Modified();
        pad->Update();
    }

    comparison->Modified();
    comparison->Update();
    comparison->SaveAs(output);
    std::cout << "Saved comparison canvas to " << output << std::endl;

    // Evaluate each (B-C)/(A-B) fit at the original histogram bin centers
    // and divide it by the reference-run fit at the same angle.
    std::size_t referenceIndex = sourceCanvases.size();
    for (std::size_t i = 0; i < sourceLabels.size(); ++i) {
        if (sourceLabels[i] == referenceRun) {
            referenceIndex = i;
            break;
        }
    }
    if (referenceIndex == sourceCanvases.size()) {
        std::cerr << "Reference run " << referenceRun
                  << " was not loaded." << std::endl;
        return;
    }

    TF1 *referenceFit = dynamic_cast<TF1 *>(
        sourceCanvases[referenceIndex]->GetPrimitive(
            "fit_mott_BmC_over_AmB"));
    TH1 *referenceHist = dynamic_cast<TH1 *>(
        sourceCanvases[referenceIndex]->GetPrimitive(
            "mott_ratio_BmC_over_AmB"));
    if (!referenceFit || !referenceHist) {
        std::cerr << "Cannot find the B-C fit or histogram in "
                  << sourceNames[referenceIndex] << std::endl;
        return;
    }

    std::vector<TH1 *> fitRatios;
    std::vector<TString> fitRatioLabels;
    double ratioMinimum = 1.;
    double ratioMaximum = 1.;

    for (std::size_t i = 0; i < sourceCanvases.size(); ++i) {
        TF1 *fit = dynamic_cast<TF1 *>(
            sourceCanvases[i]->GetPrimitive("fit_mott_BmC_over_AmB"));
        TH1 *originalHist = dynamic_cast<TH1 *>(
            sourceCanvases[i]->GetPrimitive("mott_ratio_BmC_over_AmB"));
        if (!fit || !originalHist) {
            std::cerr << "Cannot find the B-C fit or histogram in "
                      << sourceNames[i] << std::endl;
            continue;
        }

        TH1 *ratio = dynamic_cast<TH1 *>(referenceHist->Clone(
            Form("fit_ratio_circle%zu_over_circle1", i + 1)));
        ratio->SetDirectory(nullptr);
        ratio->Reset("ICES");

        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            const double theta = ratio->GetBinCenter(bin);
            if (theta < 0.5 || theta > 3.8) continue;

            const double referenceValue = referenceFit->Eval(theta);
            if (referenceValue == 0.) continue;

            const double fitValue = fit->Eval(theta);
            const double value = fitValue / referenceValue;
            const double fitError = fitStatError(fit, theta);
            const double referenceError =
                fitStatError(referenceFit, theta);

            // Central values and statistical uncertainties both come from
            // the fitted functions.
            double valueError = 0.;
            if (i == referenceIndex) {
                // Show the relative uncertainty of the reference fit around
                // the reference/reference line.
                valueError = referenceError / std::abs(referenceValue);
            } else {
                valueError =
                    std::abs(value) *
                    std::sqrt(
                        std::pow(fitError / fitValue, 2) +
                        std::pow(referenceError / referenceValue, 2));
            }

            ratio->SetBinContent(bin, value);
            ratio->SetBinError(bin, valueError);
            ratioMinimum = std::min(ratioMinimum, value - valueError);
            ratioMaximum = std::max(ratioMaximum, value + valueError);
        }

        fitRatios.push_back(ratio);
        fitRatioLabels.push_back(sourceLabels[i]);
    }

    if (fitRatios.empty()) {
        std::cerr << "No B-C fit ratios could be constructed." << std::endl;
        return;
    }

    const int colors[] = {
        kRed + 1, kBlue + 1, kGreen + 2,
        kMagenta + 1, kOrange + 7, kCyan + 2, kBlack
    };
    const int markers[] = {20, 21, 22, 23, 24, 25, 26};
    const double fitBandAlpha = 0.8;
    const int fitBandFillStyle = 3004;

    double span = ratioMaximum - ratioMinimum;
    if (span < 0.02) span = 0.02;

    TCanvas *fitRatioCanvas = new TCanvas(
        "c_circle_fit_ratio",
        Form("Circle B-C Fit Ratio Relative to %s", referenceRun),
        1100, 750);
    fitRatioCanvas->SetGrid();
    fitRatioCanvas->SetLeftMargin(0.12);
    fitRatioCanvas->SetBottomMargin(0.12);

    TLegend *fitRatioLegend = new TLegend(0.67, 0.62, 0.88, 0.88);
    fitRatioLegend->SetBorderSize(0);
    fitRatioLegend->SetFillStyle(0);

    std::vector<TGraphErrors *> fitRatioLines;
    for (std::size_t i = 0; i < fitRatios.size(); ++i) {
        TH1 *ratio = fitRatios[i];
        TGraphErrors *band = histToGraphErrors(
            ratio, Form("%s_band", ratio->GetName()), true);
        TGraphErrors *line = histToGraphErrors(
            ratio, Form("%s_line", ratio->GetName()), false);

        band->SetFillColorAlpha(colors[i % 7], fitBandAlpha);
        band->SetFillStyle(fitBandFillStyle);
        band->SetLineColor(kWhite);
        line->SetLineColor(colors[i % 7]);
        line->SetLineWidth(2);
        fitRatioLines.push_back(line);

        if (i == 0) {
            band->SetTitle(Form(
                "Circle (B-C)/(A-B) Fit Relative to %s;"
                "Reconstructed Scattering Angle [deg];"
                "Fit_{run} / Fit_{%s}",
                referenceRun, referenceRun));
            band->SetMinimum(ratioMinimum - 0.15 * span);
            band->SetMaximum(ratioMaximum + 0.15 * span);
            band->GetXaxis()->SetLimits(0.5, 3.8);
            band->GetXaxis()->CenterTitle();
            band->GetYaxis()->CenterTitle();
            band->Draw("A3");
        } else {
            band->Draw("3 SAME");
        }

        fitRatioLegend->AddEntry(
            line,
            Form("%s / %s",
                 fitRatioLabels[i].Data(),
                 referenceRun),
            "L");
    }

    for (TGraphErrors *line : fitRatioLines) {
        line->Draw("L SAME");
    }

    TLine *unity = new TLine(0.5, 1., 3.8, 1.);
    unity->SetLineColor(kBlack);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->Draw("SAME");
    fitRatioLegend->Draw();

    fitRatioCanvas->Modified();
    fitRatioCanvas->Update();
    fitRatioCanvas->SaveAs(fitRatioOutput);
    std::cout << "Saved B-C fit ratios to " << fitRatioOutput << std::endl;

    // Make the same fit-ratio comparison for B/(A-B).
    std::vector<TF1 *> bFits(sourceCanvases.size(), nullptr);
    std::vector<TH1 *> bHists(sourceCanvases.size(), nullptr);
    for (std::size_t i = 0; i < sourceCanvases.size(); ++i) {
        bHists[i] = dynamic_cast<TH1 *>(
            sourceCanvases[i]->GetPrimitive("mott_ratio_B_over_AmB"));
        if (!bHists[i]) {
            std::cerr << "Cannot find mott_ratio_B_over_AmB in "
                      << sourceNames[i] << std::endl;
            continue;
        }

        bFits[i] = new TF1(
            Form("fit_mott_B_over_AmB_%zu", i),
            "[0]+[1]*log(x)", 0.5, 3.8);
        bFits[i]->SetLineColor(kBlack);
        bHists[i]->Fit(bFits[i], "RQ0");
    }

    TF1 *bReferenceFit = bFits[referenceIndex];
    TH1 *bReferenceHist = bHists[referenceIndex];
    if (!bReferenceFit || !bReferenceHist) {
        std::cerr << "Cannot construct the B/(A-B) reference fit for "
                  << referenceRun << std::endl;
        return;
    }

    std::vector<TH1 *> bFitRatios;
    std::vector<TString> bFitRatioLabels;
    double bRatioMinimum = 1.;
    double bRatioMaximum = 1.;

    for (std::size_t i = 0; i < sourceCanvases.size(); ++i) {
        if (!bFits[i] || !bHists[i]) continue;

        TH1 *ratio = dynamic_cast<TH1 *>(bReferenceHist->Clone(
            Form("B_fit_ratio_%s_over_%s",
                 sourceLabels[i].Data(), referenceRun)));
        ratio->SetDirectory(nullptr);
        ratio->Reset("ICES");

        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            const double theta = ratio->GetBinCenter(bin);
            if (theta < 0.5 || theta > 3.8) continue;

            const double referenceValue = bReferenceFit->Eval(theta);
            const double fitValue = bFits[i]->Eval(theta);
            if (referenceValue == 0. || fitValue == 0.) continue;

            const double value = fitValue / referenceValue;
            const double fitError = fitStatError(bFits[i], theta);
            const double referenceError =
                fitStatError(bReferenceFit, theta);

            double valueError = 0.;
            if (i == referenceIndex) {
                valueError = referenceError / std::abs(referenceValue);
            } else {
                valueError =
                    std::abs(value) *
                    std::sqrt(
                        std::pow(fitError / fitValue, 2) +
                        std::pow(referenceError / referenceValue, 2));
            }

            ratio->SetBinContent(bin, value);
            ratio->SetBinError(bin, valueError);
            bRatioMinimum = std::min(bRatioMinimum, value - valueError);
            bRatioMaximum = std::max(bRatioMaximum, value + valueError);
        }

        bFitRatios.push_back(ratio);
        bFitRatioLabels.push_back(sourceLabels[i]);
    }

    if (bFitRatios.empty()) {
        std::cerr << "No B/(A-B) fit ratios could be constructed."
                  << std::endl;
        return;
    }

    double bSpan = bRatioMaximum - bRatioMinimum;
    if (bSpan < 0.02) bSpan = 0.02;

    TCanvas *bFitRatioCanvas = new TCanvas(
        "c_circle_B_fit_ratio",
        Form("Circle B/(A-B) Fit Ratio Relative to %s", referenceRun),
        1100, 750);
    bFitRatioCanvas->SetGrid();
    bFitRatioCanvas->SetLeftMargin(0.12);
    bFitRatioCanvas->SetBottomMargin(0.12);

    TLegend *bFitRatioLegend = new TLegend(0.67, 0.62, 0.88, 0.88);
    bFitRatioLegend->SetBorderSize(0);
    bFitRatioLegend->SetFillStyle(0);

    std::vector<TGraphErrors *> bFitRatioLines;
    for (std::size_t i = 0; i < bFitRatios.size(); ++i) {
        TH1 *ratio = bFitRatios[i];
        TGraphErrors *band = histToGraphErrors(
            ratio, Form("%s_band", ratio->GetName()), true);
        TGraphErrors *line = histToGraphErrors(
            ratio, Form("%s_line", ratio->GetName()), false);

        band->SetFillColorAlpha(colors[i % 7], fitBandAlpha);
        band->SetFillStyle(fitBandFillStyle);
        band->SetLineColor(kWhite);
        line->SetLineColor(colors[i % 7]);
        line->SetLineWidth(2);
        bFitRatioLines.push_back(line);

        if (i == 0) {
            band->SetTitle(Form(
                "Circle B/(A-B) Fit Relative to %s;"
                "Reconstructed Scattering Angle [deg];"
                "Fit_{run} / Fit_{%s}",
                referenceRun, referenceRun));
            band->SetMinimum(bRatioMinimum - 0.15 * bSpan);
            band->SetMaximum(bRatioMaximum + 0.15 * bSpan);
            band->GetXaxis()->SetLimits(0.5, 3.8);
            band->GetXaxis()->CenterTitle();
            band->GetYaxis()->CenterTitle();
            band->Draw("A3");
        } else {
            band->Draw("3 SAME");
        }

        bFitRatioLegend->AddEntry(
            line,
            Form("%s / %s",
                 bFitRatioLabels[i].Data(), referenceRun),
            "L");
    }

    for (TGraphErrors *line : bFitRatioLines) {
        line->Draw("L SAME");
    }

    TLine *bUnity = new TLine(0.5, 1., 3.8, 1.);
    bUnity->SetLineColor(kBlack);
    bUnity->SetLineStyle(2);
    bUnity->SetLineWidth(2);
    bUnity->Draw("SAME");
    bFitRatioLegend->Draw();

    bFitRatioCanvas->Modified();
    bFitRatioCanvas->Update();
    bFitRatioCanvas->SaveAs(bFitRatioOutput);
    std::cout << "Saved B/(A-B) fit ratios to "
              << bFitRatioOutput << std::endl;

    // Compare the live-charge-normalized mott_yield_B fits by run.
    std::vector<TF1 *> mottBYieldFits(sourceMottBYields.size(), nullptr);
    for (std::size_t i = 0; i < sourceMottBYields.size(); ++i) {
        if (!sourceMottBYields[i]) continue;

        mottBYieldFits[i] = new TF1(
            Form("fit_mott_yield_B_%s", sourceLabels[i].Data()),
            "[0]*pow(x,[1])", 0.5, 3.8);
        const int binAtOne = sourceMottBYields[i]->FindBin(1.0);
        mottBYieldFits[i]->SetParameters(
            sourceMottBYields[i]->GetBinContent(binAtOne), -2.0);
        mottBYieldFits[i]->SetLineColor(kBlack);
        sourceMottBYields[i]->Fit(mottBYieldFits[i], "RQ0");
    }

    TF1 *mottBYieldReferenceFit = mottBYieldFits[referenceIndex];
    TH1 *mottBYieldReferenceHist = sourceMottBYields[referenceIndex];
    if (!mottBYieldReferenceFit || !mottBYieldReferenceHist) {
        std::cerr << "Cannot construct the mott_yield_B reference fit for "
                  << referenceRun << std::endl;
        return;
    }

    std::vector<TH1 *> mottBYieldFitRatios;
    std::vector<TString> mottBYieldFitRatioLabels;
    double mottBYieldRatioMinimum = 1.;
    double mottBYieldRatioMaximum = 1.;

    for (std::size_t i = 0; i < sourceMottBYields.size(); ++i) {
        if (!sourceMottBYields[i] || !mottBYieldFits[i]) continue;

        TH1 *ratio = dynamic_cast<TH1 *>(mottBYieldReferenceHist->Clone(
            Form("mott_yield_B_fit_ratio_%s_over_%s",
                 sourceLabels[i].Data(), referenceRun)));
        ratio->SetDirectory(nullptr);
        ratio->Reset("ICES");

        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            const double theta = ratio->GetBinCenter(bin);
            if (theta < 0.5 || theta > 3.8) continue;

            const double referenceValue =
                mottBYieldReferenceFit->Eval(theta);
            const double fitValue = mottBYieldFits[i]->Eval(theta);
            if (referenceValue == 0. || fitValue == 0.) continue;

            const double value = fitValue / referenceValue;
            const double fitError =
                powerLawFitStatError(mottBYieldFits[i], theta);
            const double referenceError =
                powerLawFitStatError(mottBYieldReferenceFit, theta);

            double valueError = 0.;
            if (i == referenceIndex) {
                valueError = referenceError / std::abs(referenceValue);
            } else {
                valueError =
                    std::abs(value) *
                    std::sqrt(
                        std::pow(fitError / fitValue, 2) +
                        std::pow(referenceError / referenceValue, 2));
            }

            ratio->SetBinContent(bin, value);
            ratio->SetBinError(bin, valueError);
            mottBYieldRatioMinimum =
                std::min(mottBYieldRatioMinimum, value - valueError);
            mottBYieldRatioMaximum =
                std::max(mottBYieldRatioMaximum, value + valueError);
        }

        mottBYieldFitRatios.push_back(ratio);
        mottBYieldFitRatioLabels.push_back(sourceLabels[i]);
    }

    if (mottBYieldFitRatios.empty()) {
        std::cerr << "No mott_yield_B fit ratios could be constructed."
                  << std::endl;
        return;
    }

    double mottBYieldSpan =
        mottBYieldRatioMaximum - mottBYieldRatioMinimum;
    if (mottBYieldSpan < 0.02) mottBYieldSpan = 0.02;

    TCanvas *mottBYieldFitRatioCanvas = new TCanvas(
        "c_circle_mott_yield_B_fit_ratio",
        Form("Circle Mott Yield B Fit Ratio Relative to %s", referenceRun),
        1100, 750);
    mottBYieldFitRatioCanvas->SetGrid();
    mottBYieldFitRatioCanvas->SetLeftMargin(0.12);
    mottBYieldFitRatioCanvas->SetBottomMargin(0.12);

    TLegend *mottBYieldFitRatioLegend =
        new TLegend(0.67, 0.62, 0.88, 0.88);
    mottBYieldFitRatioLegend->SetBorderSize(0);
    mottBYieldFitRatioLegend->SetFillStyle(0);

    std::vector<TGraphErrors *> mottBYieldFitRatioLines;
    for (std::size_t i = 0; i < mottBYieldFitRatios.size(); ++i) {
        TH1 *ratio = mottBYieldFitRatios[i];
        TGraphErrors *band = histToGraphErrors(
            ratio, Form("%s_band", ratio->GetName()), true);
        TGraphErrors *line = histToGraphErrors(
            ratio, Form("%s_line", ratio->GetName()), false);

        band->SetFillColorAlpha(colors[i % 7], fitBandAlpha);
        band->SetFillStyle(fitBandFillStyle);
        band->SetLineColor(kWhite);
        line->SetLineColor(colors[i % 7]);
        line->SetLineWidth(2);
        mottBYieldFitRatioLines.push_back(line);

        if (i == 0) {
            band->SetTitle(Form(
                "Circle Mott Yield B Fit Relative to %s;"
                "Reconstructed Scattering Angle [deg];"
                "Fit_{run} / Fit_{%s}",
                referenceRun, referenceRun));
            band->SetMinimum(
                mottBYieldRatioMinimum - 0.15 * mottBYieldSpan);
            band->SetMaximum(
                mottBYieldRatioMaximum + 0.15 * mottBYieldSpan);
            band->GetXaxis()->SetLimits(0.5, 3.8);
            band->GetXaxis()->CenterTitle();
            band->GetYaxis()->CenterTitle();
            band->Draw("A3");
        } else {
            band->Draw("3 SAME");
        }

        mottBYieldFitRatioLegend->AddEntry(
            line,
            Form("%s / %s",
                 mottBYieldFitRatioLabels[i].Data(), referenceRun),
            "L");
    }

    for (TGraphErrors *line : mottBYieldFitRatioLines) {
        line->Draw("L SAME");
    }

    TLine *mottBYieldUnity = new TLine(0.5, 1., 3.8, 1.);
    mottBYieldUnity->SetLineColor(kBlack);
    mottBYieldUnity->SetLineStyle(2);
    mottBYieldUnity->SetLineWidth(2);
    mottBYieldUnity->Draw("SAME");
    mottBYieldFitRatioLegend->Draw();

    mottBYieldFitRatioCanvas->Modified();
    mottBYieldFitRatioCanvas->Update();
    mottBYieldFitRatioCanvas->SaveAs(mottBYieldFitRatioOutput);
    std::cout << "Saved mott_yield_B fit ratios to "
              << mottBYieldFitRatioOutput << std::endl;

    // Direct per-bin comparison of the live-charge-normalized A-B Mott
    // signal. No fit is used.
    TH1 *mottSignalReference = sourceMottSignals[referenceIndex];
    if (!mottSignalReference) {
        std::cerr << "Cannot find the A-B Mott signal for reference run "
                  << referenceRun << std::endl;
        return;
    }

    std::vector<TH1 *> mottSignalRatios;
    std::vector<TString> mottSignalRatioLabels;
    double signalRatioMinimum = 1.;
    double signalRatioMaximum = 1.;

    for (std::size_t i = 0; i < sourceMottSignals.size(); ++i) {
        TH1 *signal = sourceMottSignals[i];
        if (!signal) continue;

        TH1 *ratio = dynamic_cast<TH1 *>(mottSignalReference->Clone(
            Form("mott_signal_AmB_ratio_%s_over_%s",
                 sourceLabels[i].Data(), referenceRun)));
        ratio->SetDirectory(nullptr);
        ratio->Reset("ICES");

        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            const double theta = ratio->GetBinCenter(bin);
            if (theta < 0.5 || theta > 3.8) continue;

            const double numerator = signal->GetBinContent(bin);
            const double denominator =
                mottSignalReference->GetBinContent(bin);
            if (denominator == 0.) continue;

            const double numeratorError = signal->GetBinError(bin);
            const double denominatorError =
                mottSignalReference->GetBinError(bin);
            const double value = numerator / denominator;

            double valueError = 0.;
            if (i == referenceIndex) {
                valueError =
                    denominatorError / std::abs(denominator);
            } else {
                valueError = std::sqrt(
                    std::pow(numeratorError / denominator, 2) +
                    std::pow(numerator * denominatorError /
                                 (denominator * denominator), 2));
            }

            ratio->SetBinContent(bin, value);
            ratio->SetBinError(bin, valueError);
            signalRatioMinimum =
                std::min(signalRatioMinimum, value - valueError);
            signalRatioMaximum =
                std::max(signalRatioMaximum, value + valueError);
        }

        mottSignalRatios.push_back(ratio);
        mottSignalRatioLabels.push_back(sourceLabels[i]);
    }

    if (mottSignalRatios.empty()) {
        std::cerr << "No A-B Mott signal ratios could be constructed."
                  << std::endl;
        return;
    }

    double signalRatioSpan = signalRatioMaximum - signalRatioMinimum;
    if (signalRatioSpan < 0.02) signalRatioSpan = 0.02;

    TCanvas *mottSignalRatioCanvas = new TCanvas(
        "c_circle_mott_signal_AmB_ratio",
        Form("Circle A-B Mott Yield Ratio Relative to %s", referenceRun),
        1100, 750);
    mottSignalRatioCanvas->SetGrid();
    mottSignalRatioCanvas->SetLeftMargin(0.12);
    mottSignalRatioCanvas->SetBottomMargin(0.12);

    TLegend *mottSignalRatioLegend =
        new TLegend(0.67, 0.62, 0.88, 0.88);
    mottSignalRatioLegend->SetBorderSize(0);
    mottSignalRatioLegend->SetFillStyle(0);

    for (std::size_t i = 0; i < mottSignalRatios.size(); ++i) {
        TH1 *ratio = mottSignalRatios[i];
        ratio->SetStats(0);
        ratio->SetLineColor(colors[i % 7]);
        ratio->SetMarkerColor(colors[i % 7]);
        ratio->SetMarkerStyle(markers[i % 7]);
        ratio->SetMarkerSize(0.75);
        ratio->SetLineWidth(2);
        ratio->GetXaxis()->SetRangeUser(0.5, 3.8);

        if (i == 0) {
            ratio->SetTitle(Form(
                "Circle A-B Mott Yield Relative to %s;"
                "Reconstructed Scattering Angle [deg];"
                "(A-B)_{run} / (A-B)_{%s}",
                referenceRun, referenceRun));
            ratio->SetMinimum(
                signalRatioMinimum - 0.15 * signalRatioSpan);
            ratio->SetMaximum(
                signalRatioMaximum + 0.15 * signalRatioSpan);
            ratio->GetXaxis()->CenterTitle();
            ratio->GetYaxis()->CenterTitle();
            ratio->Draw("E1 LP");
        } else {
            ratio->Draw("E1 LP SAME");
        }

        mottSignalRatioLegend->AddEntry(
            ratio,
            Form("%s / %s",
                 mottSignalRatioLabels[i].Data(), referenceRun),
            "LP");
    }

    TLine *mottSignalUnity = new TLine(0.5, 1., 3.8, 1.);
    mottSignalUnity->SetLineColor(kBlack);
    mottSignalUnity->SetLineStyle(2);
    mottSignalUnity->SetLineWidth(2);
    mottSignalUnity->Draw("SAME");
    mottSignalRatioLegend->Draw();

    mottSignalRatioCanvas->Modified();
    mottSignalRatioCanvas->Update();
    mottSignalRatioCanvas->SaveAs(mottSignalRatioOutput);
    std::cout << "Saved A-B Mott yield bin-content ratios to "
              << mottSignalRatioOutput << std::endl;
}
