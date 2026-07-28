#include <iostream>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"


double fit_resolution(TH1 *hist, TF1 *&fit, const TString &fit_name, int color)
{
    const double peak = hist->GetBinCenter(hist->GetMaximumBin());
    fit = new TF1(fit_name, "gaus", peak - 70.0, peak + 70.0);
    fit->SetParameters(hist->GetMaximum(), peak, 30.0);
    fit->SetLineColor(color);
    fit->SetLineWidth(2);
    fit->SetLineStyle(2);

    const int fit_status = hist->Fit(fit, "RQ0");
    if (fit_status != 0) {
        std::cerr << "Warning: Gaussian fit failed for " << hist->GetName()
                  << " (status " << fit_status << ")" << std::endl;
        delete fit;
        fit = nullptr;
        return -1.0;
    }

    const double mean = fit->GetParameter(1);       // MeV
    const double sigma = std::abs(fit->GetParameter(2)); // MeV
    if (mean <= 0.0) return -1.0;

    // Project convention: (sigma / E) * sqrt(E [GeV]), expressed in percent.
    return 100.0 * sigma / mean * std::sqrt(mean / 1000.0);
}


void resolution()
{
    const TString macro_dir = gSystem->DirName(__FILE__);

    //with radiation
    const TString file_25655 = macro_dir + "/../data/x17/prad_025810_quick_check.root";
    const TString file_25560 = macro_dir + "/../data/x17/prad_025560_quick_check.root";

    TFile *f25655 = TFile::Open(file_25655, "READ");
    TFile *f25560 = TFile::Open(file_25560, "READ");
    if (!f25655 || f25655->IsZombie() || !f25560 || f25560->IsZombie()) {
        std::cerr << "Error: cannot open one or both quick-check ROOT files." << std::endl;
        if (f25655) f25655->Close();
        if (f25560) f25560->Close();
        return;
    }

    const char *modules[] = {"W459", "W460", "W493", "W494"};
    TH1 *h25655[4] = {nullptr};
    TH1 *h25560[4] = {nullptr};
    TF1 *fit25655[4] = {nullptr};
    TF1 *fit25560[4] = {nullptr};
    double res25655[4] = {-1.0, -1.0, -1.0, -1.0};
    double res25560[4] = {-1.0, -1.0, -1.0, -1.0};

    for (int i = 0; i < 4; ++i) {
        const TString hist_path = TString::Format("module_energy/h_%s", modules[i]);
        TH1 *source_25655 = dynamic_cast<TH1 *>(f25655->Get(hist_path));
        TH1 *source_25560 = dynamic_cast<TH1 *>(f25560->Get(hist_path));
        if (!source_25655 || !source_25560) {
            std::cerr << "Error: missing histogram " << hist_path << std::endl;
            f25655->Close();
            f25560->Close();
            return;
        }

        h25655[i] = dynamic_cast<TH1 *>(source_25655->Clone(
            TString::Format("h_%s_25655", modules[i])));
        h25560[i] = dynamic_cast<TH1 *>(source_25560->Clone(
            TString::Format("h_%s_25560", modules[i])));
        h25655[i]->SetDirectory(nullptr);
        h25560[i]->SetDirectory(nullptr);

        const double max_25655 = h25655[i]->GetMaximum();
        const double max_25560 = h25560[i]->GetMaximum();
        if (max_25655 > 0.0) h25655[i]->Scale(1.0 / max_25655);
        if (max_25560 > 0.0) h25560[i]->Scale(1.0 / max_25560);

        res25655[i] = fit_resolution(
            h25655[i], fit25655[i], TString::Format("fit_%s_25655", modules[i]), kRed + 1);
        res25560[i] = fit_resolution(
            h25560[i], fit25560[i], TString::Format("fit_%s_25560", modules[i]), kBlue + 1);
    }

    f25655->Close();
    f25560->Close();

    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("c_resolution", "Module energy comparison", 1400, 1000);
    canvas->Divide(2, 2, 0.01, 0.01);

    for (int i = 0; i < 4; ++i) {
        canvas->cd(i + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        h25655[i]->SetTitle(TString::Format("%s;Energy (MeV);Normalized counts", modules[i]));
        h25655[i]->SetLineColor(kRed + 1);
        h25655[i]->SetLineWidth(3);
        h25655[i]->SetMaximum(1.15);
        h25655[i]->SetMinimum(0.0);
        h25655[i]->GetXaxis()->SetRangeUser(1700.0, 2500.0);

        h25560[i]->SetLineColor(kBlue + 1);
        h25560[i]->SetLineWidth(3);

        h25655[i]->Draw("HIST");
        h25560[i]->Draw("HIST SAME");
        if (fit25655[i]) fit25655[i]->Draw("SAME");
        if (fit25560[i]) fit25560[i]->Draw("SAME");

        TLegend *legend = new TLegend(0.48, 0.70, 0.88, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        const TString label25655 = res25655[i] >= 0.0
            ? TString::Format("Now 40%% drop: #sigma/#sqrt{E} = %.3f%%", res25655[i])
            : "Run 25655: fit failed";
        const TString label25560 = res25560[i] >= 0.0
            ? TString::Format("Before: #sigma/#sqrt{E} = %.3f%%", res25560[i])
            : "Run 25560: fit failed";
        legend->AddEntry(h25655[i], label25655, "l");
        legend->AddEntry(h25560[i], label25560, "l");
        legend->Draw();
    }

    canvas->SaveAs(macro_dir + "/resolution_module_energy.png");
}
