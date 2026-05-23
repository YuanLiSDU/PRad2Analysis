#include "../EventData.h"
#include "../PhysicsTools.h"

int example_mod = 633+34-1;  // W module number: change this to switch modules
float shift_x = 0.7f, shift_y = 1.83f;
const int blocks = 8; // number of blocks per side for correction map (e.g. 8 for 8x8 blocks)

// ── Module geometry helper ────────────────────────────────────────────────────
struct ModInfo {
    TString name;          // e.g. "W635"
    float   x  = 0, y  = 0;   // module center (mm)
    float   sx = 0, sy = 0;   // module half-widths' parent size (mm)
    bool    found = false;
};

ModInfo GetModInfo(int w_num, const char* geo_file = "../hycal_modules.json") {
    ModInfo info;
    info.name = Form("W%d", w_num);
    std::ifstream fin(geo_file);
    if(!fin.is_open()) {
        Warning("GetModInfo", "Cannot open %s", geo_file);
        return info;
    }
    auto get_val = [](const std::string& line, const std::string& key) -> float {
        std::string tok = "\"" + key + "\": ";
        size_t pos = line.find(tok);
        if(pos == std::string::npos) return 0.f;
        pos += tok.size();
        return std::stof(line.substr(pos, line.find_first_of(",}", pos) - pos));
    };
    std::string target = "\"n\": \"" + std::string(info.name.Data()) + "\"";
    std::string line;
    while(std::getline(fin, line)) {
        if(line.find(target) != std::string::npos) {
            info.x  = get_val(line, "x");
            info.y  = get_val(line, "y");
            info.sx = get_val(line, "sx");
            info.sy = get_val(line, "sy");
            info.found = true;
            break;
        }
    }
    if(!info.found) Warning("GetModInfo", "Module W%d not found in %s", w_num, geo_file);
    return info;
}
// ─────────────────────────────────────────────────────────────────────────────

void energy_corr_test(){

    ModInfo mod = GetModInfo(example_mod);
    if(!mod.found) return;
    float x_lo = mod.x - mod.sx/2., x_hi = mod.x + mod.sx/2.;
    float y_lo = mod.y - mod.sy/2., y_hi = mod.y + mod.sy/2.;
    printf("Module %s: center (%.3f, %.3f)  size %.2f x %.2f mm\n",
           mod.name.Data(), mod.x, mod.y, mod.sx, mod.sy);

    // ── File 1: prad_024512_recon.root ───────────────────────────────────────
    TFile *f = TFile::Open("../data/calib/prad_024512_recon_new.root");
    TTree *tree = (TTree*)f->Get("recon");

    ReconEventData *evp = new ReconEventData();
    ReconEventData &ev = *evp;
    setupReconBranches(tree, ev);

    TH1F *mod_E = new TH1F("h1",
        Form("Energy for %s;Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200);
    
    TH1F *mod_E_corr = new TH1F("h1_corr",
        Form("Energy for %s after block-wise correction;Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200); 

    TH2F *mod_E_xy_gem = new TH2F("E_xy_gem",
        Form("Gaus fit mean/E_{exp} per %d#times%d Block in %s (GEM);x_{d};y_{d};E_{fit}/E_{exp}", blocks, blocks, mod.name.Data()),
        blocks, -0.5, 0.5, blocks, -0.5, 0.5);
    TProfile2D *mod_Eexp_xy_gem = new TProfile2D("Eexp_xy_gem",
        "Average E_{exp} per block (GEM);x_{d};y_{d};E_{exp} (MeV)",
        blocks, -0.5, 0.5, blocks, -0.5, 0.5);
    TH1F *mod_E_bloclk_gem[blocks][blocks];
    for(int bx = 0; bx < blocks; bx++){
        for(int by = 0; by < blocks; by++){
            mod_E_bloclk_gem[bx][by] = new TH1F(Form("E_block_gem_%d_%d", bx, by),
                Form("Energy block GEM (%d,%d) %s;Energy (MeV);Counts", bx, by, mod.name.Data()),
                420, 0, 4200);
        }
    }

    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
    
        if(ev.matchNum == 1 && ev.n_clusters == 1){ // look at GEM-matched events in the module
            for(int m = 0; m < ev.matchNum; m++){
                float gem_x = ev.mHit_gx[m][0];
                float gem_y = ev.mHit_gy[m][0];
                float gem_z = ev.mHit_gz[m][0];
                int det_id = ev.mHit_gid[m][0];
                GEMHit gem_hit{gem_x, gem_y, gem_z, static_cast<uint8_t>(det_id)};
                projectToHyCal(gem_hit);
                gem_hit.x += shift_x;
                gem_hit.y += shift_y;

                if(gem_hit.x < x_lo || gem_hit.x > x_hi || gem_hit.y < y_lo || gem_hit.y > y_hi)
                    continue; // only consider hits that fall within the module area after projection

                mod_E->Fill(ev.mHit_E[m]);

                float xd = (gem_hit.x - mod.x) / mod.sx;
                float yd = (gem_hit.y - mod.y) / mod.sy;
                float theta_deg = atan(sqrt(mod.x*mod.x + mod.y*mod.y) / gem_hit.z) * 180. / TMath::Pi();
                float E_expected = ExpectedEnergy(theta_deg, 3488.43f, "ep");
                int bx_gem = mod_E_xy_gem->GetXaxis()->FindBin(xd) - 1;
                int by_gem = mod_E_xy_gem->GetYaxis()->FindBin(yd) - 1;
                if(E_expected > 0 && ev.mHit_E[m] > E_expected - 200.f && ev.mHit_E[m] < E_expected + 200.f){
                    mod_Eexp_xy_gem->Fill(xd, yd, E_expected);
                    if(bx_gem >= 0 && bx_gem < blocks && by_gem >= 0 && by_gem < blocks)
                        mod_E_bloclk_gem[bx_gem][by_gem]->Fill(ev.mHit_E[m]);
                }
            }
        }
    }

    float correction[blocks][blocks] = {};

    // ── Gaussian fit per block (GEM) → fill mod_E_xy_gem ────────────────────
    for(int bx = 0; bx < blocks; bx++){
        for(int by = 0; by < blocks; by++){
            if(mod_E_bloclk_gem[bx][by]->GetEntries() < 20) continue;
            double eexp = mod_Eexp_xy_gem->GetBinContent(bx+1, by+1);
            if(eexp <= 0) continue;
            double sigma_hint = 0.03 * eexp * sqrt(eexp / 1000.);
            if(sigma_hint < 20.) sigma_hint = 20.;
            TF1 *gfit_gem = new TF1(Form("gfit_gem_%d_%d", bx, by), "gaus",
                                     eexp - 1.*sigma_hint, eexp + 1.*sigma_hint);
            gfit_gem->SetParameters(mod_E_bloclk_gem[bx][by]->GetMaximum(), eexp, sigma_hint);
            mod_E_bloclk_gem[bx][by]->Fit(gfit_gem, "RQ0");
            if(gfit_gem->GetParameter(2) > 0){
                mod_E_xy_gem->SetBinContent(bx+1, by+1, gfit_gem->GetParameter(1) / eexp);
                correction[bx][by] = gfit_gem->GetParameter(1) / eexp;
            }
        }
    }

    // Apply block-wise correction to all hits in the module
    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
    
        if(ev.matchNum == 1 && ev.n_clusters == 1){ // look at GEM-matched events in the module
            for(int m = 0; m < ev.matchNum; m++){
                float gem_x = ev.mHit_gx[m][0];
                float gem_y = ev.mHit_gy[m][0];
                float gem_z = ev.mHit_gz[m][0];
                int det_id = ev.mHit_gid[m][0];
                GEMHit gem_hit{gem_x, gem_y, gem_z, static_cast<uint8_t>(det_id)};
                projectToHyCal(gem_hit);
                gem_hit.x += shift_x;
                gem_hit.y += shift_y;

                if(gem_hit.x < x_lo || gem_hit.x > x_hi || gem_hit.y < y_lo || gem_hit.y > y_hi)
                    continue; // only consider hits that fall within the module area after projection

                float xd = (gem_hit.x - mod.x) / mod.sx;
                float yd = (gem_hit.y - mod.y) / mod.sy;
                int bx_gem = mod_E_xy_gem->GetXaxis()->FindBin(xd) - 1;
                int by_gem = mod_E_xy_gem->GetYaxis()->FindBin(yd) - 1;
                double corr = (bx_gem >= 0 && bx_gem < blocks && by_gem >= 0 && by_gem < blocks) ? correction[bx_gem][by_gem] : 1.0;
                if(corr > 0)                    
                    mod_E_corr->Fill(ev.mHit_E[m] / corr);
                //else if(corr == 0)               
                 //   mod_E_corr->Fill(ev.mHit_E[m]);
                else                              
                    Warning("energy_corr_test", "Negative correction factor %.3f at block (%d, %d)", corr, bx_gem, by_gem);
            }
        }
    }

    TCanvas *c = new TCanvas("c", "Energy comparison", 1200, 800);
    mod_E->SetLineColor(kBlue);
    mod_E->SetLineWidth(2);
    mod_E->Draw("HIST");
    mod_E_corr->SetLineColor(kRed);
    mod_E_corr->SetLineWidth(2);
    mod_E_corr->Draw("HIST SAME");

    mod_E->Fit("gaus", "Q","r", 3400, 3530);
    mod_E_corr->Fit("gaus", "Q","r", 3400, 3530);

    float mean = mod_E->GetFunction("gaus")->GetParameter(1);
    float sigma = mod_E->GetFunction("gaus")->GetParameter(2);
    float res = sigma / mean * sqrt(mean / 1000.);

    float mean_corr = mod_E_corr->GetFunction("gaus")->GetParameter(1);
    float sigma_corr = mod_E_corr->GetFunction("gaus")->GetParameter(2);
    float res_corr = sigma_corr / mean_corr * sqrt(mean_corr / 1000.);

    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(mod_E, Form("Before correction: res = %.2f%%", res * 100), "L");
    leg->AddEntry(mod_E_corr, Form("After correction: res = %.2f%%", res_corr * 100), "L");
    leg->Draw();

    TCanvas *c2 = new TCanvas("c2", "Block-wise correction map", 1200, 800);
    gStyle->SetPaintTextFormat(".3f");
    mod_E_xy_gem->SetMarkerSize(2.0);
    mod_E_xy_gem->Draw("colz text");



}