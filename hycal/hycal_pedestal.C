
#include "EventData.h"
#include <map>
#include <vector>
#include <algorithm>

void setupRawBranches(TTree *tree, RawEventData &ev);

void hycal_pedestal()
{
    TChain *chain = new TChain("events");
    chain->Add("data/raw/prad_023899.00000_raw.root");
    TTree *t = chain;

    RawEventData data;
    setupRawBranches(t, data);

    std::map<int, TH1F*> h_ped_mean, h_ped_rms;

    int nEntries = t->GetEntries();
    std::cout << "Total entries: " << nEntries << std::endl;

    for (int i = 0; i < nEntries; i++) {
        t->GetEntry(i);
        if (i % 10000 == 0)
            std::cout << "Processing " << i << "/" << nEntries << "\r" << std::flush;

        for (int j = 0; j < data.nch; j++) {
            int id = data.module_id[j];
            if (h_ped_mean.find(id) == h_ped_mean.end()) {
                h_ped_mean[id] = new TH1F(Form("h_pmean_%d", id),
                                          Form("Module %d  ped mean; ADC counts; Events", id),
                                          200, 0, 600);
                h_ped_rms[id]  = new TH1F(Form("h_prms_%d",  id),
                                          Form("Module %d  ped rms;  ADC counts; Events", id),
                                          200, 0, 10);
            }
            h_ped_mean[id]->Fill(data.ped_mean[j]);
            h_ped_rms[id]->Fill(data.ped_rms[j]);
        }
    }
    std::cout << "\nDone. Modules found: " << h_ped_mean.size() << std::endl;

    // ── helper lambda: Gaussian fit → {mean, sigma}; fallback to hist stats ──
    auto gaussFit = [](TH1F* h) -> std::pair<double,double> {
        if (h->GetEntries() < 10) return {h->GetMean(), h->GetRMS()};
        double mu  = h->GetBinCenter(h->GetMaximumBin());
        double sig = h->GetRMS();
        TF1 *f = new TF1("_gaus","gaus", mu - 3*sig, mu + 3*sig);
        f->SetParameters(h->GetMaximum(), mu, sig);
        h->Fit(f, "RQ0");   // R=range, Q=quiet, 0=no draw
        double fm = f->GetParameter(1);
        double fs = TMath::Abs(f->GetParameter(2));
        delete f;
        return {fm, fs};
    };

    // ── separate W and G module IDs, sorted ──────────────────────────────────
    std::vector<int> w_ids, g_ids;
    for (auto &kv : h_ped_mean) {
        int id = kv.first;
        if (id >= 1001 && id <= 2156)      w_ids.push_back(id);  // W modules
        else if (id >= 1 && id <= 1728)    g_ids.push_back(id);  // G modules
        // skip id >= 3000 (LMS / veto / other)
    }
    std::sort(w_ids.begin(), w_ids.end());
    std::sort(g_ids.begin(), g_ids.end());
    std::cout << "W modules: " << w_ids.size()
              << "  G modules: " << g_ids.size() << std::endl;

    // ── write JSON ────────────────────────────────────────────────────────────
    FILE *jf = fopen("data/pedestal_prad_023899.json", "w");
    fprintf(jf, "[\n");
    bool first_entry = true;

    auto writeEntry = [&](const std::string &name,
                          double pm_mean, double pm_sigma,
                          double rms_mean, double rms_sigma) {
        if (!first_entry) fprintf(jf, ",\n");
        first_entry = false;
        fprintf(jf,
            "  {\"name\":\"%s\","
            "\"ped_mean\":%.4f,\"ped_mean_sigma\":%.4f,"
            "\"rms_mean\":%.4f,\"rms_sigma\":%.4f}",
            name.c_str(), pm_mean, pm_sigma, rms_mean, rms_sigma);
    };

    // W modules: Gaussian fit both ped_mean and ped_rms
    for (int id : w_ids) {
        std::string name = Form("W%d", id - 1000);
        auto [pm, ps] = gaussFit(h_ped_mean[id]);
        auto [rm, rs] = gaussFit(h_ped_rms[id]);
        writeEntry(name, pm, ps, rm, rs);
    }
    // G modules: histogram mean / std dev (no fit)
    for (int id : g_ids) {
        std::string name = Form("G%d", id);
        double pm  = h_ped_mean[id]->GetMean();
        double ps  = h_ped_mean[id]->GetRMS();
        double rm  = h_ped_rms[id]->GetMean();
        double rs  = h_ped_rms[id]->GetRMS();
        writeEntry(name, pm, ps, rm, rs);
    }

    fprintf(jf, "\n]\n");
    fclose(jf);
    std::cout << "JSON saved: data/pedestal_prad_023899.json" << std::endl;

    // ── display W500-W600 on 10 canvases (5 ped_mean + 5 ped_rms) ────────────
    std::vector<int> show_ids;
    for (int id : w_ids)
        if (id >= 1500 && id <= 1600) show_ids.push_back(id);
    std::cout << "Displaying " << show_ids.size() << " modules (W500-W600)" << std::endl;

    const int kPadsX = 4, kPadsY = 5, kPerCanvas = kPadsX * kPadsY;
    const int kCanvases = 5;
    for (int ic = 0; ic < kCanvases; ic++) {
        int idx0 = ic * kPerCanvas;
        if (idx0 >= (int)show_ids.size()) break;
        int first_w = show_ids[idx0] - 1000;
        int last_idx = std::min(idx0 + kPerCanvas - 1, (int)show_ids.size() - 1);
        int last_w  = show_ids[last_idx] - 1000;

        TCanvas *cm = new TCanvas(Form("c_pmean_%d", ic),
                                  Form("Ped Mean – W%d to W%d", first_w, last_w),
                                  1600, 1000);
        cm->Divide(kPadsX, kPadsY);
        TCanvas *cr = new TCanvas(Form("c_prms_%d", ic),
                                  Form("Ped RMS – W%d to W%d", first_w, last_w),
                                  1600, 1000);
        cr->Divide(kPadsX, kPadsY);

        for (int ip = 0; ip < kPerCanvas; ip++) {
            int idx = idx0 + ip;
            if (idx >= (int)show_ids.size()) break;
            int id = show_ids[idx];
            cm->cd(ip + 1); h_ped_mean[id]->Draw();
            cr->cd(ip + 1); h_ped_rms[id]->Draw();
        }
    }
}

void setupRawBranches(TTree *tree, RawEventData &ev)
{
    // event header
    tree->SetBranchAddress("event_num",        &ev.event_num);
    tree->SetBranchAddress("trigger_type",     &ev.trigger_type);
    tree->SetBranchAddress("trigger_bits",     &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",        &ev.timestamp);

    // HyCal per-channel data  (branches stored under "hycal." prefix)
    tree->SetBranchAddress("hycal.nch",        &ev.nch);
    tree->SetBranchAddress("hycal.module_id",  ev.module_id);
    tree->SetBranchAddress("hycal.ped_mean",   ev.ped_mean);
    tree->SetBranchAddress("hycal.ped_rms",    ev.ped_rms);
};