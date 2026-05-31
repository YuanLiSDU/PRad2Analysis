#include "../EventData.h"
#include "../PhysicsTools.h"
#include <thread>
#include <mutex>
#include <atomic>

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

    const size_t nF = file_list.size();
    std::vector<TH1F*> hs_list(nF, nullptr), hf_list(nF, nullptr), ht_list(nF, nullptr);
    std::atomic<size_t> next_idx{0};
    std::atomic<size_t> done_cnt{0};
    std::mutex log_mtx;

    auto worker = [&](){
        while(true){
            size_t i = next_idx.fetch_add(1);
            if(i >= nF) break;
            const std::string& fname = file_list[i];

            TFile fin(fname.c_str());
            if(fin.IsZombie()) continue;
            TTree *tt = (TTree*)fin.Get("events");
            if(!tt) continue;

            RawEventData ev;
            setupRawBranches(tt, ev);

            TH1F *hs = new TH1F("", "", 5000, 0.f, 5000.f);
            TH1F *hf = new TH1F("", "", 5000, 0.f, 5000.f);
            TH1F *ht = new TH1F("", "", 100, 0.f, 400.f);

            Long64_t nE = tt->GetEntries();
            for(Long64_t k = 0; k < nE; k++){
                tt->GetEntry(k);

                bool isPulser = (ev.trigger_bits & 1u << 15) != 0;
                bool isRawsum = (ev.trigger_bits & 1u << 8) != 0;
                bool isAlpha  = (ev.trigger_bits & 1u << 25) != 0;
                if(!isPulser && !isRawsum) continue;
                if(!isPulser && !isAlpha) continue;

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

                    if(time_soft > 0.f) ht->Fill(time_soft);
                }

                if(rawsum_soft > 0.f) hs->Fill(rawsum_soft);
                if(rawsum_firm > 0.f) hf->Fill(rawsum_firm);
            }
            hs_list[i] = hs;
            hf_list[i] = hf;
            ht_list[i] = ht;

            size_t d = ++done_cnt;
            std::lock_guard<std::mutex> lk(log_mtx);
            std::cout << "[" << d << "/" << nF << "] done: " << fname << std::endl;
        }
    };

    // ---- launch threads ----
    int nth = std::min<int>(n_threads, (int)nF);
    std::vector<std::thread> threads;
    threads.reserve(nth);
    for(int t = 0; t < nth; t++) threads.emplace_back(worker);
    for(auto& th : threads) th.join();

    // ---- merge results in main thread ----
    TH1::AddDirectory(true);
    TH1F *h_rawsum_soft = new TH1F("h_rawsum_soft", "RawSum from software peaks;RawSum;Counts", 5000, 0.f, 5000.f);
    TH1F *h_rawsum_firm = new TH1F("h_rawsum_firm", "RawSum from firmware peaks;RawSum;Counts", 5000, 0.f, 5000.f);
    TH1F *h_time_soft = new TH1F("h_time_soft", "Peak time distribution (software);Time (ns);Counts", 100, 0.f, 400.f);
    for(size_t i = 0; i < nF; i++){
        if(hs_list[i]){ h_rawsum_soft->Add(hs_list[i]); delete hs_list[i]; }
        if(hf_list[i]){ h_rawsum_firm->Add(hf_list[i]); delete hf_list[i]; }
        if(ht_list[i]){ h_time_soft->Add(ht_list[i]); delete ht_list[i]; }
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