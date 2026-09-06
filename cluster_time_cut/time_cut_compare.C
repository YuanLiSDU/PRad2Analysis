#include <algorithm>
#include <numeric>
#include <vector>

const int module = 660;

bool fit_resolution(TH1F *hist, double &resolution)
{
    if (!hist || hist->GetEntries() <= 2000) {
        return false;
    }

    double mean = hist->GetBinCenter(hist->GetMaximumBin());
    if (mean <= 0.) {
        return false;
    }

    TFitResultPtr fit_result = hist->Fit("gaus", "QS", "", mean - 80., mean + 80.);
    if (fit_result.Get() == nullptr || fit_result->Status() != 0) {
        return false;
    }

    TF1 *fit = hist->GetFunction("gaus");
    if (!fit) {
        return false;
    }

    fit->SetLineStyle(2);
    fit->SetLineWidth(2);
    fit->SetLineColor(kBlack);

    const double fit_mean = fit->GetParameter(1);
    const double fit_sigma = fit->GetParameter(2);
    if (fit_mean <= 0. || fit_sigma <= 0.) {
        return false;
    }

    resolution = fit_sigma / fit_mean * sqrt(fit_mean / 1000.) * 100.;
    return true;
}

TH1F *build_resolution_hist(const char *name, const char *title, const std::vector<double> &values,
                            int bins, double min, double max, Color_t color)
{
    TH1F *hist = new TH1F(name, title, bins, min, max);
    hist->SetDirectory(nullptr);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);

    for (double value : values) {
        hist->Fill(value);
    }

    return hist;
}

void time_cut_compare()
{
    TFile *file_0ns = TFile::Open("0ns_quick.root");
    TFile *file_4ns = TFile::Open("4ns_quick.root");
    TFile *file_16ns = TFile::Open("16ns_quick.root");

    TH1F *h_0ns = (TH1F*)file_0ns->Get(Form("module_energy/h_W%d", module));
    TH1F *h_4ns = (TH1F*)file_4ns->Get(Form("module_energy/h_W%d", module));
    TH1F *h_16ns = (TH1F*)file_16ns->Get(Form("module_energy/h_W%d", module));

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->SetGrid();
    h_0ns->SetLineColor(kRed);
    h_4ns->SetLineColor(kBlue);
    h_16ns->SetLineColor(kRed);

    double resolution_0ns = 0.;
    double resolution_4ns = 0.;
    double resolution_16ns = 0.;

    fit_resolution(h_0ns, resolution_0ns);
    fit_resolution(h_4ns, resolution_4ns);
    fit_resolution(h_16ns, resolution_16ns);

    h_4ns->GetXaxis()->SetRangeUser(1800., 2460.);

    //h_0ns->Draw();
    h_4ns->Draw("");
    h_16ns->Draw("SAME");

    TF1 *fit_4ns = h_4ns->GetFunction("gaus");
    TF1 *fit_16ns = h_16ns->GetFunction("gaus");
    if (fit_4ns) {
        fit_4ns->SetLineStyle(2);
        fit_4ns->SetLineColor(h_4ns->GetLineColor());
        fit_4ns->Draw("SAME");
    }
    if (fit_16ns) {
        fit_16ns->SetLineStyle(2);
        fit_16ns->SetLineColor(h_16ns->GetLineColor());
        fit_16ns->Draw("SAME");
    }

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legend->AddEntry(h_0ns, Form("0 ns (%.3f%%)", resolution_0ns), "l");
    legend->AddEntry(h_4ns, Form("4 ns (%.3f%%)", resolution_4ns), "l");
    legend->AddEntry(h_16ns, Form("16 ns (%.3f%%)", resolution_16ns), "l");
    legend->Draw();

}

void compare_resolution_distribution()
{
    TFile *file_4ns = TFile::Open("4ns_quick.root");
    TFile *file_16ns = TFile::Open("16ns_quick.root");
    if (!file_4ns || !file_16ns) {
        return;
    }

    TDirectory *dir_4ns = file_4ns->GetDirectory("module_energy");
    TDirectory *dir_16ns = file_16ns->GetDirectory("module_energy");
    if (!dir_4ns || !dir_16ns) {
        return;
    }

    std::vector<double> resolutions_4ns;
    std::vector<double> resolutions_16ns;

    TIter next_key(dir_4ns->GetListOfKeys());
    TKey *key = nullptr;
    while ((key = (TKey*)next_key())) {
        TString hist_name = key->GetName();
        if (!hist_name.BeginsWith("h_W")) {
            continue;
        }

        TH1F *h_4ns = (TH1F*)dir_4ns->Get(hist_name);
        TH1F *h_16ns = (TH1F*)dir_16ns->Get(hist_name);

        double resolution_4ns = 0.;
        double resolution_16ns = 0.;

        const bool fit_4ns = fit_resolution(h_4ns, resolution_4ns);
        const bool fit_16ns = fit_resolution(h_16ns, resolution_16ns);
        if (fit_4ns && fit_16ns) {
            resolutions_4ns.push_back(resolution_4ns);
            resolutions_16ns.push_back(resolution_16ns);
        }
    }

    if (resolutions_4ns.empty() || resolutions_16ns.empty()) {
        return;
    }

    auto mean_below_4 = [](const std::vector<double> &values) {
        std::vector<double> filtered;
        filtered.reserve(values.size());
        for (double value : values) {
            if (value < 4.0) {
                filtered.push_back(value);
            }
        }
        if (filtered.empty()) {
            return 0.0;
        }
        return std::accumulate(filtered.begin(), filtered.end(), 0.0) / static_cast<double>(filtered.size());
    };

    const double mean_4ns_below_4 = mean_below_4(resolutions_4ns);
    const double mean_16ns_below_4 = mean_below_4(resolutions_16ns);

    const double min_resolution = std::min(*std::min_element(resolutions_4ns.begin(), resolutions_4ns.end()),
                                           *std::min_element(resolutions_16ns.begin(), resolutions_16ns.end()));
    const double max_resolution = std::max(*std::max_element(resolutions_4ns.begin(), resolutions_4ns.end()),
                                           *std::max_element(resolutions_16ns.begin(), resolutions_16ns.end()));
    const double hist_min = std::max(0., min_resolution - 0.2);
    const double hist_max = max_resolution + 0.2;

    TH1F *h_res_4ns = build_resolution_hist("h_res_4ns",
                                            "Module resolution comparison;Resolution (%);Modules",
                                            resolutions_4ns, 80*4, hist_min, hist_max, kBlue);
    TH1F *h_res_16ns = build_resolution_hist("h_res_16ns",
                                             "Module resolution comparison;Resolution (%);Modules",
                                             resolutions_16ns, 80*4, hist_min, hist_max, kGreen + 2);

    const double max_count = std::max(h_res_4ns->GetMaximum(), h_res_16ns->GetMaximum());
    h_res_4ns->SetStats(0);
    h_res_4ns->SetMaximum(max_count * 1.15);

    TCanvas *c = new TCanvas("c_resolution_compare", "c_resolution_compare", 900, 700);
    h_res_4ns->GetXaxis()->SetRangeUser(2.5, 4.0);
    h_res_4ns->Draw("HIST");
    h_res_16ns->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.62, 0.72, 0.88, 0.88);
    legend->AddEntry(h_res_4ns, Form("4 ns (%zu fits, mean = %.3f%%)", resolutions_4ns.size(), mean_4ns_below_4), "l");
    legend->AddEntry(h_res_16ns, Form("16 ns (%zu fits, mean = %.3f%%)", resolutions_16ns.size(), mean_16ns_below_4), "l");
    legend->Draw();
}