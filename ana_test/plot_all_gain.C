#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <sys/stat.h>

// Extract the quoted string value from a JSON line like:  "key": "value",
static std::string jsonStrVal(const std::string &line) {
    size_t c = line.find(':');
    if (c == std::string::npos) return "";
    size_t q1 = line.find('"', c + 1);
    if (q1 == std::string::npos) return "";
    size_t q2 = line.find('"', q1 + 1);
    if (q2 == std::string::npos) return "";
    return line.substr(q1 + 1, q2 - q1 - 1);
}

// Extract the integer value from a JSON line like:  "key": 42,
static int jsonIntVal(const std::string &line) {
    size_t c = line.find(':');
    if (c == std::string::npos) return -1;
    std::string v = line.substr(c + 1);
    v.erase(0, v.find_first_not_of(" \t"));
    return std::stoi(v);
}

void plot_all_gain() {

    const char *dir     = "../data/gain_corr/24827";
    const char *outdir  = "gain_all_plots";
    const char *hv_json = "../hycal/hv_settings_20260509_220811.json";
    mkdir(outdir, 0755);

    // ---- Parse HV JSON: build (crate, slot) -> sorted list of 0-based W indices ----
    struct CrateSlot {
        std::string crate;
        int slot;
        bool operator<(const CrateSlot &o) const {
            return crate != o.crate ? crate < o.crate : slot < o.slot;
        }
    };
    std::map<CrateSlot, std::vector<int>> cs_mods;

    {
        std::ifstream fin(hv_json);
        if (!fin.is_open()) {
            std::cerr << "Cannot open " << hv_json << std::endl;
            return;
        }
        std::string line, cur_crate, cur_name;
        while (std::getline(fin, line)) {
            size_t p = line.find_first_not_of(" \t");
            if (p == std::string::npos) continue;
            std::string tl = line.substr(p);
            if      (tl.rfind("\"crate\"", 0) == 0) cur_crate = jsonStrVal(line);
            else if (tl.rfind("\"name\"",  0) == 0) cur_name  = jsonStrVal(line);
            else if (tl.rfind("\"slot\"",  0) == 0) {
                int cur_slot = jsonIntVal(line);
                // add only W modules (name = "W" followed by digits)
                if (cur_name.size() > 1 && cur_name[0] == 'W') {
                    bool alldig = true;
                    for (size_t i = 1; i < cur_name.size(); i++)
                        if (!isdigit(cur_name[i])) { alldig = false; break; }
                    if (alldig)
                        cs_mods[{cur_crate, cur_slot}].push_back(std::stoi(cur_name.substr(1)) - 1);
                }
                cur_name.clear();
            }
        }
    }

    if (cs_mods.empty()) { std::cerr << "No crate/slot data parsed." << std::endl; return; }
    std::cout << "Parsed " << cs_mods.size() << " crate+slot groups." << std::endl;

    // ---- Read all gain_corr data ----
    std::vector<TString> files;
    {
        TSystemDirectory sd(dir, dir);
        TList *fl = sd.GetListOfFiles();
        TIter nx(fl);
        TSystemFile *sf;
        while ((sf = (TSystemFile*)nx())) {
            TString fn = sf->GetName();
            if (fn.EndsWith("_corr.root"))
                files.push_back(TString(dir) + "/" + fn);
        }
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) { std::cerr << "No _corr.root files found." << std::endl; return; }

    const int NMOD = 1156;
    std::vector<std::vector<double>> avg(NMOD);
    float gain_corr_W[NMOD][3];

    for (const TString &fpath : files) {
        TFile *f = TFile::Open(fpath);
        if (!f || f->IsZombie()) { delete f; continue; }
        TTree *t = (TTree*)f->Get("gain_corr");
        if (!t) { delete f; continue; }
        t->SetBranchAddress("gain_corr_W", gain_corr_W);
        Long64_t nentries = t->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            for (int m = 0; m < NMOD; m++)
                avg[m].push_back((1./gain_corr_W[m][1] + 1./gain_corr_W[m][2]) / 2.0);
        }
        delete f;
    }

    int N = (int)avg[0].size();
    if (N == 0) { std::cerr << "No data entries." << std::endl; return; }
    std::cout << "Total time samples: " << N << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.09);

    const int COLS = 8;   // fixed 8 columns; rows vary with module count

    int ig = 0;
    for (auto &[cs, mods_unsorted] : cs_mods) {
        std::vector<int> mods = mods_unsorted;
        std::sort(mods.begin(), mods.end());
        int nmods = (int)mods.size();
        int rows  = (nmods + COLS - 1) / COLS;

        // extract crate number from "PRadHV_N"
        int crate_num = 0;
        size_t us = cs.crate.rfind('_');
        if (us != std::string::npos) crate_num = std::stoi(cs.crate.substr(us + 1));

        TString canv_title = Form("Crate %d  Slot %d  (%d modules)", crate_num, cs.slot, nmods);
        TString fname      = Form("%s/crate%d_slot%02d.png", outdir, crate_num, cs.slot);

        TCanvas *c = new TCanvas(Form("c_%d", ig++), canv_title.Data(), COLS * 200, rows * 200);
        c->Divide(COLS, rows, 0.002, 0.002);

        for (int ip = 0; ip < nmods; ip++) {
            int midx = mods[ip];
            if (midx < 0 || midx >= NMOD) continue;

            c->cd(ip + 1);
            gPad->SetMargin(0.20, 0.04, 0.18, 0.16);

            TGraph *gr = new TGraph(N);
            double mean = 0., rms = 0.;
            for (int i = 0; i < N; i++) {
                gr->SetPoint(i, i + 1, avg[midx][i]);
                mean += avg[midx][i];
            }
            mean /= N;
            for (int i = 0; i < N; i++)
                rms += (avg[midx][i] - mean) * (avg[midx][i] - mean);
            rms = (N > 0) ? std::sqrt(rms / N) : 0.;
            double pct = (mean != 0.) ? rms / std::fabs(mean) * 100. : 0.;

            gr->SetTitle(Form("W%d  %.2f%%;Idx;gain", midx + 1, pct));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.3);
            gr->SetMarkerColor(kBlue + 1);
            gr->SetLineColor(kBlue + 1);
            gr->SetLineWidth(1);

            double ymin = *std::min_element(avg[midx].begin(), avg[midx].end());
            double ymax = *std::max_element(avg[midx].begin(), avg[midx].end());
            double margin = (ymax > ymin) ? (ymax - ymin) * 0.15 : 0.01;
            gr->GetYaxis()->SetRangeUser(ymin - margin, ymax + margin);
            gr->GetYaxis()->SetNdivisions(504);
            gr->GetXaxis()->SetNdivisions(504);
            gr->Draw("AL");
        }

        c->SaveAs(fname);
        std::cout << "Saved: " << fname << "  (" << nmods << " modules)" << std::endl;
        delete c;
    }

    std::cout << "Done. " << ig << " canvases saved to " << outdir << "/" << std::endl;
}
