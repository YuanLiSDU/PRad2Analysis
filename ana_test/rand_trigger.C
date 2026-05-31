#include "../EventData.h"
#include "../PhysicsTools.h"
#include "ROOT/TThreadExecutor.hxx"

const int max_files = 100;
const int n_threads = 8;

// Load all "factor" values for W modules from a calibration JSON file into a
// vector indexed by W number (idx 0 -> W1, idx 1 -> W2, ...). The JSON layout
// is an array of objects, each spanning multiple lines, with "name": "W###"
// and "factor": <float>.
static const std::vector<float>& LoadWFactors(const char* json_file){
    static std::map<std::string, std::vector<float>> cache;
    auto it = cache.find(json_file);
    if(it != cache.end()) return it->second;

    std::vector<float>& vec = cache[json_file];
    std::ifstream fin(json_file);
    if(!fin.is_open()){
        Warning("LoadWFactors", "Cannot open %s", json_file);
        return vec;
    }

    std::string line;
    int  cur_w     = -1;
    float cur_fac  = 0.f;
    bool has_fac   = false;
    while(std::getline(fin, line)){
        size_t pn = line.find("\"name\"");
        if(pn != std::string::npos){
            size_t q1 = line.find('"', line.find(':', pn) + 1);
            size_t q2 = (q1 == std::string::npos) ? std::string::npos : line.find('"', q1 + 1);
            if(q1 != std::string::npos && q2 != std::string::npos){
                std::string nm = line.substr(q1 + 1, q2 - q1 - 1);
                if(nm.size() > 1 && nm[0] == 'W'){
                    cur_w = std::atoi(nm.c_str() + 1);
                } else {
                    cur_w = -1;
                }
            }
        }
        size_t pf = line.find("\"factor\"");
        if(pf != std::string::npos){
            size_t colon = line.find(':', pf);
            if(colon != std::string::npos){
                cur_fac  = std::stof(line.substr(colon + 1));
                has_fac  = true;
            }
        }
        if(cur_w > 0 && has_fac){
            if((int)vec.size() < cur_w) vec.resize(cur_w, 0.f);
            vec[cur_w - 1] = cur_fac;
            cur_w   = -1;
            has_fac = false;
        }
    }
    return vec;
}

// Return the soft (FPGA peak integral) calibration factor for W module `w_id`
// (1-based, i.e. raw module_id - 1000). Returns 0 if not found.
float GetSoftFactor(int w_id, const char* json_file = "../calibration/calibration_factor_3p5.json"){
    const auto& v = LoadWFactors(json_file);
    if(w_id < 1 || w_id > (int)v.size()) return 0.f;
    return v[w_id - 1];
}

// Return the firm (DAQ FADC integral) calibration factor for W module `w_id`.
float GetFirmFactor(int w_id, const char* json_file = "../calibration/calibration_factor_fadc_3p5.json"){
    const auto& v = LoadWFactors(json_file);
    if(w_id < 1 || w_id > (int)v.size()) return 0.f;
    return v[w_id - 1];
}

void rand_trigger(const char* input_dir = "../data/0.7GeV"){
    // prime calibration caches in main thread (avoid race on first access)
    GetSoftFactor(1);
    GetFirmFactor(1);

    // ---- collect input files ----
    std::vector<std::string> file_list;
    TSystemDirectory dir(input_dir, input_dir);
    TList *files = dir.GetListOfFiles();
    if(!files){
        std::cerr << "Cannot open directory " << input_dir << std::endl;
        return;
    }
    TIter next(files);
    while(TSystemFile *file = (TSystemFile*)next()){
        if(file->IsDirectory()) continue;
        TString name = file->GetName();
        if(!name.EndsWith(".root")) continue;
        file_list.emplace_back(Form("%s/%s", input_dir, name.Data()));
        std::cout << "Added: " << file_list.back() << std::endl;
        if((int)file_list.size() >= max_files) break;
    }
    if(file_list.empty()){
        std::cerr << "No root files found, abort." << std::endl;
        return;
    }

    // ---- per-file worker (runs on worker threads) ----
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(false);

    auto worker = [](const std::string& fname) -> std::array<TH1F*, 2> {
        TFile fin(fname.c_str());
        if(fin.IsZombie()) return {nullptr, nullptr};
        TTree *tt = (TTree*)fin.Get("events");
        if(!tt) return {nullptr, nullptr};

        RawEventData ev;
        setupRawBranches(tt, ev);

        auto *hs = new TH1F("", "", 5000, 0.f, 5000.f);
        auto *hf = new TH1F("", "", 5000, 0.f, 5000.f);
        auto *ht = new TH1F("", "", 100, 0.f, 400.f);

        Long64_t nE = tt->GetEntries();
        for(Long64_t i = 0; i < nE; i++){
            tt->GetEntry(i);

            bool isPulser = (ev.trigger_bits & 1u << 15) != 0;
            bool isRawsum = (ev.trigger_bits & 1u << 8) != 0;
            if(!isPulser && !isRawsum) continue;
            if(!isPulser) continue;

            float rawsum_soft = 0.f, rawsum_firm = 0.f;

            for(int j = 0; j < ev.nch; j++){
                int mod_id = ev.module_id[j];
                if(mod_id <= 1000 || mod_id >= 3000) continue;
                mod_id -= 1000;

                float gain = ev.gain_factor[j];
                float integral_soft = 0.f, integral_firm = 0.f;
                float time_soft = 0.f;

                const float t_min = 0.f;
                const float t_max = 1000.f;
                int   best_p      = -1;
                float best_height = -1.f;
                for(int p = 0; p < ev.npeaks[j]; p++){
                    float time = ev.peak_time[j][p];
                    if(time < t_min || time > t_max) continue;
                    float height = ev.peak_height[j][p];
                    if(height > best_height){
                        best_height = height;
                        best_p      = p;
                    }
                }
                if(best_p > 0){
                    integral_soft  = ev.peak_integral[j][best_p];
                    integral_soft *= gain;
                    time_soft      = ev.peak_time[j][best_p];
                }
                integral_firm = ev.daq_peak_integral[j][0];

                rawsum_soft += integral_soft * GetSoftFactor(mod_id);
                rawsum_firm += integral_firm * GetFirmFactor(mod_id);

                ht->Fill(time_soft);
            }

            hs->Fill(rawsum_soft);
            hf->Fill(rawsum_firm);
        }
        return {hs, hf, ht};
    };

    // ---- parallel map over files ----
    ROOT::TThreadExecutor pool(n_threads);
    auto results = pool.Map(worker, file_list);

    // ---- merge results in main thread ----
    TH1::AddDirectory(true);
    TH1F *h_rawsum_soft = new TH1F("h_rawsum_soft", "RawSum from software peaks;RawSum;Counts", 5000, 0.f, 5000.f);
    TH1F *h_rawsum_firm = new TH1F("h_rawsum_firm", "RawSum from firmware peaks;RawSum;Counts", 5000, 0.f, 5000.f);
    TH1F *h_time_soft = new TH1F("h_time_soft", "Peak time distribution (software);Time (ns);Counts", 100, 0.f, 400.f);
    for(auto& pr : results){
        if(pr[0]){ h_rawsum_soft->Add(pr[0]); delete pr[0]; }
        if(pr[1]){ h_rawsum_firm->Add(pr[1]); delete pr[1]; }
        if(pr[2]){ h_time_soft->Add(pr[2]); delete pr[2]; }
    }

    TCanvas *c = new TCanvas("c", "RawSum Comparison", 1200, 600);
    c->Divide(2, 1);
    c->cd(1);
    h_rawsum_soft->Draw("E");
    c->cd(2);
    h_rawsum_firm->Draw("E");
    c->SaveAs("rand_trigger_comparison.png");

    TCanvas *c2 = new TCanvas("c2", "Software Peak Time Distribution", 600, 600);
    h_time_soft->Draw("E");
    c2->SaveAs("rand_trigger_time.png");

}