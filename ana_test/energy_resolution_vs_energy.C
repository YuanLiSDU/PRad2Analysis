#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "../EventData.h"
#include "../PhysicsTools.h"

namespace EnergyResolutionVsEnergy {

// Input follows position_sigma.C: one beam energy for each input file.
const int kNFiles = 3;
const char *kInputFiles[kNFiles] = {
    "../data/recon/0.7GeV/prad_024655_recon.root",
    "../data/recon/2.2GeV/prad_025202_filtered.root",
    "../data/recon/3.5GeV/prad_024917_recon.root"
};
const double kBeamEnergy[kNFiles] = {728.9, 2239.5, 3484.0}; // MeV
const char *kBeamLabel[kNFiles] = {"0p7GeV", "2p2GeV", "3p5GeV"};

// Set to 1 to process every event.  The default matches position_sigma.C.
const int kEntryDivisor = 5;
const double kHyCalZ = 6270.0; // mm
const double kEnergyBinWidth = 10.0; // MeV
const double kEnergyMax = 4200.0; // MeV
const double kMollerSelectionResolution = 0.035;
std::mutex gLogMutex;

const int kNAngleBins = 33;
const double kAngleEdge[kNAngleBins + 1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

struct PeakResult {
    int fileIndex = -1;
    int angleBin = -1;
    std::string reaction;
    double theta = 0.0;
    double thetaError = 0.0;
    double expectedEnergy = 0.0;
    double mean = 0.0;
    double meanError = 0.0;
    double sigma = 0.0;
    double sigmaError = 0.0;
    double ratio = 0.0;
    double ratioError = 0.0;
    double resolution = 0.0;      // sigma / mean, percent
    double resolutionError = 0.0; // percent
};

bool inHyCalAcceptance(double x, double y)
{
    const double module = 20.75; // mm
    return (std::fabs(x) > 2.5 * module || std::fabs(y) > 2.5 * module)
        && std::fabs(x) < 16.0 * module
        && std::fabs(y) < 16.0 * module;
}

int findAngleBin(double theta)
{
    const double *upper = std::upper_bound(kAngleEdge, kAngleEdge + kNAngleBins + 1,
                                            theta);
    const int bin = static_cast<int>(upper - kAngleEdge) - 1;
    return (bin >= 0 && bin < kNAngleBins) ? bin : -1;
}

bool hasRequiredGemMatch(uint32_t matchFlag)
{
    const bool upstream = (matchFlag & (1u << 0)) || (matchFlag & (1u << 1));
    const bool downstream = (matchFlag & (1u << 2)) || (matchFlag & (1u << 3));
    return upstream && downstream;
}

bool projectDownstreamGem(const ReconEventData &ev, int cluster, double &x,
                          double &y, double &z)
{
    const int gem = (ev.matchFlag[cluster] & (1u << 2)) ? 2 : 3;
    const double gemZ = ev.matchGEMz[cluster][gem];
    if (!std::isfinite(gemZ) || std::fabs(gemZ) < 1.0e-6) return false;

    const double scale = kHyCalZ / gemZ;
    x = ev.matchGEMx[cluster][gem] * scale;
    y = ev.matchGEMy[cluster][gem] * scale;
    z = kHyCalZ;
    return std::isfinite(x) && std::isfinite(y);
}

void processInputFile(int fileIndex, TH1F **energyHistogram,
                      TH1F **mollerEnergyHistogram)
{
    TFile *inputFile = TFile::Open(kInputFiles[fileIndex]);
    if (!inputFile || inputFile->IsZombie()) {
        std::lock_guard<std::mutex> lock(gLogMutex);
        std::cerr << "Cannot open input file: " << kInputFiles[fileIndex] << std::endl;
        if (inputFile) inputFile->Close();
        return;
    }
    TTree *tree = dynamic_cast<TTree *>(inputFile->Get("recon"));
    if (!tree) {
        std::lock_guard<std::mutex> lock(gLogMutex);
        std::cerr << "Cannot find tree 'recon' in " << kInputFiles[fileIndex] << std::endl;
        inputFile->Close();
        return;
    }

    ReconEventData event;
    setupReconBranches(tree, event);
    const Long64_t entries = tree->GetEntries() / kEntryDivisor;
    {
        std::lock_guard<std::mutex> lock(gLogMutex);
        std::cout << "[worker " << fileIndex + 1 << "/" << kNFiles << "] Processing "
                  << kInputFiles[fileIndex] << " (" << entries
                  << " entries, Ebeam=" << kBeamEnergy[fileIndex] << " MeV)"
                  << std::endl;
    }

    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if (entry > 0 && entry % 1000000 == 0) {
            std::lock_guard<std::mutex> lock(gLogMutex);
            std::cout << "[worker " << fileIndex + 1 << "] " << entry << " / "
                      << entries << std::endl;
        }

        for (int cluster = 0; cluster < event.n_clusters; ++cluster) {
            if (event.cl_nblocks[cluster] < 5) continue;
            if (!hasRequiredGemMatch(event.matchFlag[cluster])) continue;

            double x = 0.0, y = 0.0, z = 0.0;
            if (!projectDownstreamGem(event, cluster, x, y, z)) continue;
            if (!inHyCalAcceptance(x, y)) continue;

            // Deliberately no module-local xd/yd cut here.
            const double theta = std::atan2(std::hypot(x, y), z)
                               * 180.0 / TMath::Pi();
            const int angleBin = findAngleBin(theta);
            if (angleBin < 0) continue;
            energyHistogram[angleBin]->Fill(event.cl_energy[cluster]);
        }

        // Double-arm Moller selection.  Fill each electron at its own angle,
        // so every angle bin has an independent ee energy spectrum.
        if (event.n_clusters != 2 || event.matchNum != 2) continue;

        double x[2] = {}, y[2] = {}, z[2] = {}, theta[2] = {};
        bool validArms = true;
        for (int arm = 0; arm < 2; ++arm) {
            if (event.cl_nblocks[arm] < 5
                || !hasRequiredGemMatch(event.matchFlag[arm])
                || !projectDownstreamGem(event, arm, x[arm], y[arm], z[arm])
                || !inHyCalAcceptance(x[arm], y[arm])) {
                validArms = false;
                break;
            }
            theta[arm] = std::atan2(std::hypot(x[arm], y[arm]), z[arm])
                       * 180.0 / TMath::Pi();
        }
        if (!validArms) continue;

        if (!isMoller_kinematic(theta[0], event.cl_energy[0],
                                theta[1], event.cl_energy[1],
                                kBeamEnergy[fileIndex],
                                kMollerSelectionResolution)) {
            continue;
        }

        MollerEvent moller = {
            DataPoint(x[0], y[0], z[0], event.cl_energy[0]),
            DataPoint(x[1], y[1], z[1], event.cl_energy[1])
        };
        if (std::fabs(GetMollerPhiDiff(moller)) > 10.0) continue;

        for (int arm = 0; arm < 2; ++arm) {
            const int angleBin = findAngleBin(theta[arm]);
            if (angleBin >= 0)
                mollerEnergyHistogram[angleBin]->Fill(event.cl_energy[arm]);
        }
    }

    inputFile->Close();
    {
        std::lock_guard<std::mutex> lock(gLogMutex);
        std::cout << "[worker " << fileIndex + 1 << "] Finished "
                  << kInputFiles[fileIndex] << std::endl;
    }
}

bool fitPeak(TH1F *hist, double expectedEnergy, int fileIndex, int angleBin,
             const char *reaction, PeakResult &result)
{
    if (!hist || hist->GetEntries() < 30 || expectedEnergy <= 0.0
        || expectedEnergy >= kEnergyMax) {
        return false;
    }

    // HyCal stochastic resolution estimate, used only to locate and initialize the peak.
    const double sigmaGuess = 0.03 * expectedEnergy
                            / std::sqrt(expectedEnergy / 1000.0);
    const double searchLow = std::max(0.0, expectedEnergy - 3.0 * sigmaGuess);
    const double searchHigh = std::min(kEnergyMax, expectedEnergy + 3.0 * sigmaGuess);
    int firstBin = std::max(1, hist->FindBin(searchLow));
    int lastBin = std::min(hist->GetNbinsX(), hist->FindBin(searchHigh));
    if (lastBin < firstBin || hist->Integral(firstBin, lastBin) < 30.0) return false;

    int peakBin = firstBin;
    for (int bin = firstBin + 1; bin <= lastBin; ++bin) {
        if (hist->GetBinContent(bin) > hist->GetBinContent(peakBin)) peakBin = bin;
    }
    if (hist->GetBinContent(peakBin) < 4.0) return false;

    const double peak = hist->GetXaxis()->GetBinCenter(peakBin);
    const double fitLow = std::max(0.0, peak - 2.0 * sigmaGuess);
    const double fitHigh = std::min(kEnergyMax, peak + 2.0 * sigmaGuess);
    TF1 *fit = new TF1(Form("fit_%s_f%d_b%d", reaction, fileIndex, angleBin),
                       "gaus", fitLow, fitHigh);
    fit->SetParameters(hist->GetBinContent(peakBin), peak, sigmaGuess);

    TFitResultPtr fitResult = hist->Fit(fit, "RQ0S+");
    const double mean = fit->GetParameter(1);
    const double sigma = std::fabs(fit->GetParameter(2));
    const double meanError = fit->GetParError(1);
    const double sigmaError = fit->GetParError(2);
    const bool valid = static_cast<int>(fitResult) == 0
        && fit->GetNDF() > 0
        && std::isfinite(mean) && std::isfinite(sigma)
        && std::isfinite(meanError) && std::isfinite(sigmaError)
        && mean > 0.0 && sigma > 0.0
        && meanError > 0.0 && sigmaError > 0.0
        && std::fabs(mean - expectedEnergy) < 3.0 * sigmaGuess
        && sigma < 3.0 * sigmaGuess;
    if (!valid) return false;

    result.fileIndex = fileIndex;
    result.angleBin = angleBin;
    result.reaction = reaction;
    result.theta = 0.5 * (kAngleEdge[angleBin] + kAngleEdge[angleBin + 1]);
    result.thetaError = 0.5 * (kAngleEdge[angleBin + 1] - kAngleEdge[angleBin]);
    result.expectedEnergy = expectedEnergy;
    result.mean = mean;
    result.meanError = meanError;
    result.sigma = sigma;
    result.sigmaError = sigmaError;
    result.ratio = mean / expectedEnergy;
    result.ratioError = meanError / expectedEnergy;
    result.resolution = 100.0 * sigma / mean;
    result.resolutionError = result.resolution
        * std::sqrt(std::pow(sigmaError / sigma, 2) + std::pow(meanError / mean, 2));
    return true;
}

TGraphErrors *makeThetaGraph(const std::vector<PeakResult> &points, bool ratio,
                             const char *name)
{
    TGraphErrors *graph = new TGraphErrors();
    graph->SetName(name);
    for (const PeakResult &point : points) {
        const double y = ratio ? point.ratio : point.resolution;
        const double yError = ratio ? point.ratioError : point.resolutionError;
        graph->SetPoint(graph->GetN(), point.theta, y);
        graph->SetPointError(graph->GetN() - 1, point.thetaError, yError);
    }
    return graph;
}

void styleGraph(TGraphErrors *graph, Color_t color, Style_t marker)
{
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.85);
    graph->SetLineWidth(1);
}

void drawThetaSummary(const std::vector<PeakResult> ep[kNFiles],
                      const std::vector<PeakResult> ee[kNFiles],
                      const char *outputDirectory, TFile *outputFile)
{
    const Color_t colors[kNFiles] = {kBlue + 1, kMagenta + 1, kRed + 1};
    const Style_t markers[kNFiles] = {20, 21, 22};
    TGraphErrors *ratioEp[kNFiles] = {};
    TGraphErrors *resolutionEp[kNFiles] = {};
    TGraphErrors *ratioEe[kNFiles] = {};
    TGraphErrors *resolutionEe[kNFiles] = {};

    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        ratioEp[fileIndex] = makeThetaGraph(
            ep[fileIndex], true, Form("ratio_ep_%s", kBeamLabel[fileIndex]));
        resolutionEp[fileIndex] = makeThetaGraph(
            ep[fileIndex], false, Form("resolution_ep_%s", kBeamLabel[fileIndex]));
        ratioEe[fileIndex] = makeThetaGraph(
            ee[fileIndex], true, Form("ratio_double_arm_ee_%s", kBeamLabel[fileIndex]));
        resolutionEe[fileIndex] = makeThetaGraph(
            ee[fileIndex], false,
            Form("resolution_double_arm_ee_%s", kBeamLabel[fileIndex]));
        styleGraph(ratioEp[fileIndex], colors[fileIndex], markers[fileIndex]);
        styleGraph(resolutionEp[fileIndex], colors[fileIndex], markers[fileIndex]);
        styleGraph(ratioEe[fileIndex], colors[fileIndex], markers[fileIndex] + 4);
        styleGraph(resolutionEe[fileIndex], colors[fileIndex], markers[fileIndex] + 4);
    }

    double maxResolution = 0.0;
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex)
        for (const PeakResult &point : ep[fileIndex])
            maxResolution = std::max(maxResolution, point.resolution + point.resolutionError);
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex)
        for (const PeakResult &point : ee[fileIndex])
            maxResolution = std::max(maxResolution,
                                     point.resolution + point.resolutionError);
    if (maxResolution <= 0.0) maxResolution = 10.0;

    TCanvas *canvas = new TCanvas("c_energy_summary_vs_theta",
                                  "Energy response and resolution vs theta", 900, 900);
    canvas->Divide(1, 2, 0.0, 0.0);
    canvas->cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.03);
    TH1F *ratioFrame = gPad->DrawFrame(kAngleEdge[0], 0.8, kAngleEdge[kNAngleBins], 1.2);
    ratioFrame->SetTitle(";Scattering angle #theta (deg);E_{recon mean}/E_{expect}");
    ratioFrame->GetXaxis()->SetLabelSize(0.0);
    TLine *unity = new TLine(kAngleEdge[0], 1.0, kAngleEdge[kNAngleBins], 1.0);
    unity->SetLineStyle(7);
    unity->Draw("SAME");

    TLegend *legend = new TLegend(0.14, 0.55, 0.55, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.033);
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        ratioEp[fileIndex]->Draw("P SAME");
        legend->AddEntry(ratioEp[fileIndex],
                         Form("ep, E_{beam} = %.3g GeV", kBeamEnergy[fileIndex] / 1000.0),
                         "p");
        ratioEe[fileIndex]->Draw("P SAME");
        legend->AddEntry(ratioEe[fileIndex],
                         Form("double-arm ee, E_{beam} = %.3g GeV",
                              kBeamEnergy[fileIndex] / 1000.0), "p");
    }
    legend->Draw();

    canvas->cd(2);
    gPad->SetGrid();
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.14);
    TH1F *resolutionFrame = gPad->DrawFrame(kAngleEdge[0], 0.0,
                                             kAngleEdge[kNAngleBins],
                                             1.25 * maxResolution);
    resolutionFrame->SetTitle(";Scattering angle #theta (deg);Energy resolution #sigma/E_{recon mean} (%)");
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        resolutionEp[fileIndex]->Draw("P SAME");
        resolutionEe[fileIndex]->Draw("P SAME");
    }

    canvas->SaveAs(Form("%s/summary_vs_theta.png", outputDirectory));
    outputFile->cd();
    canvas->Write();
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        ratioEp[fileIndex]->Write();
        resolutionEp[fileIndex]->Write();
        ratioEe[fileIndex]->Write();
        resolutionEe[fileIndex]->Write();
    }
}

void drawResolutionVsEnergy(const std::vector<PeakResult> ep[kNFiles],
                            const std::vector<PeakResult> ee[kNFiles],
                            const char *outputDirectory, TFile *outputFile)
{
    for (int angleBin = 0; angleBin < kNAngleBins; ++angleBin) {
        TGraphErrors *epGraph = new TGraphErrors();
        epGraph->SetName(Form("resolution_vs_energy_ep_bin%02d", angleBin));
        TGraphErrors *eeGraph = new TGraphErrors();
        eeGraph->SetName(Form("resolution_vs_energy_ee_bin%02d", angleBin));

        double minEnergy = std::numeric_limits<double>::max();
        double maxEnergy = 0.0;
        double maxResolution = 0.0;
        for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
            for (const PeakResult &point : ep[fileIndex]) {
                if (point.angleBin != angleBin) continue;
                const double energyGeV = point.expectedEnergy / 1000.0;
                epGraph->SetPoint(epGraph->GetN(), energyGeV, point.resolution);
                epGraph->SetPointError(epGraph->GetN() - 1, 0.0,
                                       point.resolutionError);
                minEnergy = std::min(minEnergy, energyGeV);
                maxEnergy = std::max(maxEnergy, energyGeV);
                maxResolution = std::max(maxResolution,
                                         point.resolution + point.resolutionError);
            }
        }
        for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex)
            for (const PeakResult &point : ee[fileIndex]) {
                if (point.angleBin != angleBin) continue;
                const double energyGeV = point.expectedEnergy / 1000.0;
                eeGraph->SetPoint(eeGraph->GetN(), energyGeV, point.resolution);
                eeGraph->SetPointError(eeGraph->GetN() - 1, 0.0,
                                       point.resolutionError);
                minEnergy = std::min(minEnergy, energyGeV);
                maxEnergy = std::max(maxEnergy, energyGeV);
                maxResolution = std::max(maxResolution,
                                         point.resolution + point.resolutionError);
            }
        if (epGraph->GetN() + eeGraph->GetN() == 0) continue;

        styleGraph(epGraph, kBlue + 1, 20);
        styleGraph(eeGraph, kGreen + 2, 24);
        epGraph->Sort();
        eeGraph->Sort();

        // Use the same energy-dependence fit as position_sigma.C.
        TGraphErrors *combinedGraph = new TGraphErrors();
        combinedGraph->SetName(Form("resolution_vs_energy_combined_bin%02d", angleBin));
        TGraphErrors *sourceGraphs[2] = {epGraph, eeGraph};
        for (TGraphErrors *source : sourceGraphs) {
            for (int pointIndex = 0; pointIndex < source->GetN(); ++pointIndex) {
                double x = 0.0, y = 0.0;
                source->GetPoint(pointIndex, x, y);
                combinedGraph->SetPoint(combinedGraph->GetN(), x, y);
                combinedGraph->SetPointError(combinedGraph->GetN() - 1,
                                             source->GetErrorX(pointIndex),
                                             source->GetErrorY(pointIndex));
            }
        }
        combinedGraph->Sort();

        TF1 *energyFit = nullptr;
        if (combinedGraph->GetN() >= 2) {
            energyFit = new TF1(Form("resolution_energy_fit_bin%02d", angleBin),
                                "[0]/sqrt(x) + [1]/x + [2]", 0.1, 10.0);
            energyFit->SetParNames("k1", "k2", "k3");
            energyFit->SetParameters(combinedGraph->GetY()[0]
                                     * std::sqrt(combinedGraph->GetX()[0]),
                                     0.0, 0.0);
            energyFit->SetLineColor(kRed + 1);
            energyFit->SetLineWidth(2);
            combinedGraph->Fit(energyFit, "Q");
        } else {
            std::cerr << "Angle bin " << angleBin << " has only "
                      << combinedGraph->GetN()
                      << " valid point; need at least 2 for the resolution fit."
                      << std::endl;
        }

        const double xPadding = std::max(0.05, 0.08 * (maxEnergy - minEnergy));
        const double xLow = std::max(0.0, minEnergy - xPadding);
        const double xHigh = maxEnergy + xPadding;
        const double yHigh = maxResolution > 0.0 ? 1.3 * maxResolution : 10.0;

        TCanvas *canvas = new TCanvas(Form("c_resolution_vs_energy_bin%02d", angleBin),
                                      "Resolution vs energy", 800, 650);
        canvas->SetGrid();
        TH1F *frame = canvas->DrawFrame(xLow, 0.0, xHigh, yHigh);
        frame->SetTitle(Form("#theta #in [%.3f, %.3f] deg;E_{expect} (GeV);Energy resolution #sigma/E_{recon mean} (%%)",
                             kAngleEdge[angleBin], kAngleEdge[angleBin + 1]));
        epGraph->Draw("P SAME");
        eeGraph->Draw("P SAME");
        if (energyFit) energyFit->Draw("SAME");

        TLegend *legend = new TLegend(0.62, 0.73, 0.89, 0.89);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        if (epGraph->GetN() > 0) legend->AddEntry(epGraph, "ep", "p");
        if (eeGraph->GetN() > 0)
            legend->AddEntry(eeGraph, "double-arm ee", "p");
        if (energyFit) legend->AddEntry(energyFit, "Combined fit", "l");
        legend->Draw();

        if (energyFit) {
            TLatex *fitLabel = new TLatex();
            fitLabel->SetNDC();
            fitLabel->SetTextSize(0.030);
            fitLabel->DrawLatex(
                0.15, 0.87,
                "R(E) = k_{1}/#sqrt{E [GeV]} + k_{2}/E [GeV] + k_{3}");
            fitLabel->DrawLatex(
                0.15, 0.82,
                Form("k_{1} = %.4f #pm %.4f %% #sqrt{GeV}",
                     energyFit->GetParameter(0), energyFit->GetParError(0)));
            fitLabel->DrawLatex(
                0.15, 0.77,
                Form("k_{2} = %.4f #pm %.4f %% GeV",
                     energyFit->GetParameter(1), energyFit->GetParError(1)));
            fitLabel->DrawLatex(
                0.15, 0.72,
                Form("k_{3} = %.4f #pm %.4f %%",
                     energyFit->GetParameter(2), energyFit->GetParError(2)));
        }

        canvas->SaveAs(Form("%s/resolution_vs_energy_bin%02d_%.3f_%.3fdeg.png",
                            outputDirectory, angleBin, kAngleEdge[angleBin],
                            kAngleEdge[angleBin + 1]));
        outputFile->cd();
        canvas->Write();
        epGraph->Write();
        eeGraph->Write();
        combinedGraph->Write();
        if (energyFit) energyFit->Write();
    }
}

} // namespace EnergyResolutionVsEnergy

void energy_resolution_vs_energy()
{
    using namespace EnergyResolutionVsEnergy;

    ROOT::EnableThreadSafety();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    const char *outputDirectory = "energy_resolution_vs_energy_plots";
    gSystem->mkdir(outputDirectory, true);

    TH1F *energyHistogram[kNFiles][kNAngleBins] = {};
    TH1F *mollerEnergyHistogram[kNFiles][kNAngleBins] = {};
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        for (int angleBin = 0; angleBin < kNAngleBins; ++angleBin) {
            energyHistogram[fileIndex][angleBin] = new TH1F(
                Form("energy_%s_bin%02d", kBeamLabel[fileIndex], angleBin),
                Form("E_{beam}=%.1f MeV, #theta #in [%.3f, %.3f] deg;Reconstructed energy (MeV);Counts / %.0f MeV",
                     kBeamEnergy[fileIndex], kAngleEdge[angleBin],
                     kAngleEdge[angleBin + 1], kEnergyBinWidth),
                static_cast<int>(kEnergyMax / kEnergyBinWidth), 0.0, kEnergyMax);
            energyHistogram[fileIndex][angleBin]->SetDirectory(nullptr);

            mollerEnergyHistogram[fileIndex][angleBin] = new TH1F(
                Form("moller_energy_%s_bin%02d", kBeamLabel[fileIndex], angleBin),
                Form("Double-arm Moller, E_{beam}=%.1f MeV, #theta #in [%.3f, %.3f] deg;Reconstructed energy (MeV);Counts / %.0f MeV",
                     kBeamEnergy[fileIndex], kAngleEdge[angleBin],
                     kAngleEdge[angleBin + 1], kEnergyBinWidth),
                static_cast<int>(kEnergyMax / kEnergyBinWidth), 0.0, kEnergyMax);
            mollerEnergyHistogram[fileIndex][angleBin]->SetDirectory(nullptr);
        }
    }

    // One independent input file per worker.  All ROOT drawing and fitting stays
    // on the main thread after the workers have joined.
    std::thread workers[kNFiles];
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        workers[fileIndex] = std::thread(processInputFile, fileIndex,
                                         energyHistogram[fileIndex],
                                         mollerEnergyHistogram[fileIndex]);
    }
    for (std::thread &worker : workers) worker.join();

    std::cout << "All three workers finished; starting fits and plotting."
              << std::endl;

    TFile *outputFile = TFile::Open(
        Form("%s/energy_resolution_vs_energy.root", outputDirectory), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Cannot create output ROOT file in " << outputDirectory << std::endl;
        return;
    }

    std::vector<PeakResult> epResult[kNFiles];
    std::vector<PeakResult> eeResult[kNFiles];
    for (int fileIndex = 0; fileIndex < kNFiles; ++fileIndex) {
        for (int angleBin = 0; angleBin < kNAngleBins; ++angleBin) {
            const double theta = 0.5 * (kAngleEdge[angleBin] + kAngleEdge[angleBin + 1]);
            TH1F *hist = energyHistogram[fileIndex][angleBin];

            PeakResult ep;
            const double epExpected = ExpectedEnergy(theta, kBeamEnergy[fileIndex], "ep");
            if (fitPeak(hist, epExpected, fileIndex, angleBin, "ep", ep))
                epResult[fileIndex].push_back(ep);

            PeakResult ee;
            const double eeExpected = ExpectedEnergy(theta, kBeamEnergy[fileIndex], "ee");
            if (fitPeak(mollerEnergyHistogram[fileIndex][angleBin], eeExpected,
                        fileIndex, angleBin, "ee", ee)) {
                eeResult[fileIndex].push_back(ee);
            }

            outputFile->cd();
            hist->Write();
            mollerEnergyHistogram[fileIndex][angleBin]->Write();
        }
        std::cout << "Valid ep fits at " << kBeamEnergy[fileIndex] / 1000.0
                  << " GeV: " << epResult[fileIndex].size() << std::endl;
        std::cout << "Valid double-arm ee fits at "
                  << kBeamEnergy[fileIndex] / 1000.0 << " GeV: "
                  << eeResult[fileIndex].size() << std::endl;
    }

    drawThetaSummary(epResult, eeResult, outputDirectory, outputFile);
    drawResolutionVsEnergy(epResult, eeResult, outputDirectory, outputFile);

    outputFile->Close();
    std::cout << "Saved plots and ROOT objects in " << outputDirectory << "/" << std::endl;
}
