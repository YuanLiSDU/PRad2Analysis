#include "../EventData.h"
#include "../PhysicsTools.h"
#include <fstream>
#include <string>

int example_mod = 633;  // W module number: change this to switch modules
float shift_x = 0.7f, shift_y = 1.83f;

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

void energy_plot(){

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

    TH2F *hit_all = new TH2F("hit_all", "Hit distribution for all events;X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *hit_all_gem = new TH2F("hit_all_gem", "Hit distribution for all events (GEM projection);X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *h2 = new TH2F("h2", "Energy vs. Angle;Energy (MeV);Scattering Angle (degrees)", 100, 0, 5, 4200/2, 0, 4200);
    TH1F *mod_E = new TH1F("h1",
        Form("Energy for %s;Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200);
    TH1F *mod_E_center = new TH1F("h2",
        Form("Energy for %s (center);Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200);
    TH1F *mod_E_ratio = new TH1F("h3",
        Form("E/E_{exp} for %s;E_{rec}/E_{exp};Counts", mod.name.Data()), 500, 0, 2);
    TH1F *mod_E_bloclk[5][5];
    for(int bx = 0; bx < 5; bx++){
        for(int by = 0; by < 5; by++){
            mod_E_bloclk[bx][by] = new TH1F(Form("E_block_%d_%d", bx, by),
                Form("Energy block (%d,%d) %s;Energy (MeV);Counts", bx, by, mod.name.Data()),
                420, 0, 4200);
        }
    }
    TProfile2D *mod_E_xy = new TProfile2D("E_xy",
        Form("E_{rec}/E_{exp} per 5#times5 Block in %s;x_{d};y_{d};E_{rec}/E_{exp}", mod.name.Data()),
        5, -0.5, 0.5, 5, -0.5, 0.5);
    TH2F *mod_hit_xy = new TH2F("hit_xy",
        Form("Hit distribution in %s;x_{d};y_{d}", mod.name.Data()),
        75, -0.5, 0.5, 75, -0.5, 0.5);

    TProfile2D *mod_E_xy_gem = new TProfile2D("E_xy_gem",
        Form("E_{rec}/E_{exp} per 5#times5 Block in %s (GEM projection);x_{d};y_{d};E_{rec}/E_{exp}", mod.name.Data()),
        5, -0.5, 0.5, 5, -0.5, 0.5);
    
    //delta X between Hycal X and GEM X
    TProfile2D *mod_delta_x = new TProfile2D("delta_x", "Delta X between HyCal and GEM;x/d;y/d;X_{HyCal} - X_{GEM} (mm)",
        50, -0.5, 0.5, 1, -0.5, 0.5);
    TProfile2D *mod_delta_y = new TProfile2D("delta_y", "Delta Y between HyCal and GEM;x/d;y/d;Y_{HyCal} - Y_{GEM} (mm)",
        1, -0.5, 0.5, 50, -0.5, 0.5);
    TH1F *mod_delta_x_proj = new TH1F("delta_x_proj", "Delta X between HyCal and GEM (projected);X_{HyCal} - X_{GEM} (mm);Counts",
        100, -5, 5);
    TH1F *mod_delta_y_proj = new TH1F("delta_y_proj", "Delta Y between HyCal and GEM (projected);Y_{HyCal} - Y_{GEM} (mm);Counts",
        100, -5, 5);
    TProfile *h1_delta_x = new TProfile("h1_delta_x", "Mean #DeltaX vs x_{d};x_{d};#LT X_{HyCal} - X_{GEM,proj} #GT (mm)",
        50, -0.5, 0.5);
    TProfile *h1_delta_y = new TProfile("h1_delta_y", "Mean #DeltaY vs y_{d};y_{d};#LT Y_{HyCal} - Y_{GEM,proj} #GT (mm)",
        50, -0.5, 0.5);

    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
        
        for(int j = 0; j < ev.n_clusters; j++){
            if(ev.cl_nblocks[j] < 3) continue; // only consider clusters with 3 or more blocks to reduce noise
            
            ev.cl_x[j] += shift_x; // apply position shifts to align with expected module center
            ev.cl_y[j] += shift_y;

            hit_all->Fill(ev.cl_x[j], ev.cl_y[j]);
            
            float theta_deg = atan(sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]) / ev.cl_z[j]) * 180. / TMath::Pi();
            h2->Fill(theta_deg, ev.cl_energy[j]);

            if(ev.cl_x[j] > x_lo && ev.cl_x[j] < x_hi &&
               ev.cl_y[j] > y_lo && ev.cl_y[j] < y_hi && ev.n_clusters == 1){ // focus on single-cluster events in module 635
                float xd = (ev.cl_x[j] - mod.x) / mod.sx;
                float yd = (ev.cl_y[j] - mod.y) / mod.sy;
                mod_E->Fill(ev.cl_energy[j]);
                if(fabs(xd) < 0.3 && fabs(yd) < 0.3) mod_E_center->Fill(ev.cl_energy[j]);
                float E_expected = ExpectedEnergy(theta_deg, 3488.43f, "ep");
                if(E_expected > 0){
                    mod_E_ratio->Fill(ev.cl_energy[j] / E_expected);
                    mod_E_xy->Fill(xd, yd, ev.cl_energy[j] / E_expected);
                }
                int bx = mod_E_xy->GetXaxis()->FindBin(xd) - 1;
                int by = mod_E_xy->GetYaxis()->FindBin(yd) - 1;
                if(bx >= 0 && bx < 5 && by >= 0 && by < 5)
                    mod_E_bloclk[bx][by]->Fill(ev.cl_energy[j]);
                mod_hit_xy->Fill(xd, yd);
            }
        }
        if(ev.matchNum == 1 && ev.n_clusters == 1){ // also look at GEM-matched events in the module
            for(int m = 0; m < ev.matchNum; m++){
                float gem_x = ev.mHit_gx[m][0];
                float gem_y = ev.mHit_gy[m][0];
                float gem_z = ev.mHit_gz[m][0];
                int det_id = ev.mHit_gid[m][0];
                GEMHit gem_hit{gem_x, gem_y, gem_z, static_cast<uint8_t>(det_id)};
                projectToHyCal(gem_hit);
                gem_hit.x += shift_x;
                gem_hit.y += shift_y;
                mod_delta_x_proj->Fill(ev.mHit_x[m] - ev.mHit_gx[m][0]);
                mod_delta_y_proj->Fill(ev.mHit_y[m] - ev.mHit_gy[m][0]);
                hit_all_gem->Fill(gem_hit.x, gem_hit.y);

                float xd = (gem_hit.x - mod.x) / mod.sx;
                float yd = (gem_hit.y - mod.y) / mod.sy;
                float theta_deg = atan(sqrt(mod.x*mod.x + mod.y*mod.y) / gem_hit.z) * 180. / TMath::Pi();
                float E_expected = ExpectedEnergy(theta_deg, 3488.43f, "ep");
                if(E_expected > 0){
                    mod_E_xy_gem->Fill(xd, yd, ev.mHit_E[m] / E_expected);
                }
                mod_delta_x->Fill(xd, yd, ev.mHit_x[m] + shift_x - gem_hit.x);
                mod_delta_y->Fill(xd, yd, ev.mHit_y[m] + shift_y - gem_hit.y);
                h1_delta_x->Fill(xd, ev.mHit_x[m] + shift_x - gem_hit.x);
                h1_delta_y->Fill(yd, ev.mHit_y[m] + shift_y - gem_hit.y);
            }
        }
    }

    // ── Drawing ───────────────────────────────────────────────────────────────
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    h2->Draw("colz");

    // Energy comparison
    TCanvas *c2 = new TCanvas("c2", Form("Energy comparison for %s", mod.name.Data()), 1200, 800);
    mod_E_center->Scale(mod_E->GetMaximum() / mod_E_center->GetMaximum()); // scale center histogram to match the max of all hits for better visual comparison
    mod_E->SetLineColor(kBlue);
    mod_E->SetLineWidth(2);
    mod_E->Draw("HIST");
    mod_E_center->SetLineColor(kRed);
    mod_E_center->SetLineWidth(2);
    mod_E_center->Draw("HIST SAME");
    mod_E_center->Fit("gaus", "Q","r", 3400, 3500);
    
    float mean = mod_E_center->GetFunction("gaus")->GetParameter(1);
    float sigma = mod_E_center->GetFunction("gaus")->GetParameter(2);
    float res = sigma / mean * sqrt(mean / 1000.);
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(mod_E, "All hits", "L");
    leg->AddEntry(mod_E_center, Form("center hits (res = %.2f%%)", res * 100), "L");
    leg->Draw();

    // E/E_exp comparison
    TCanvas *c3 = new TCanvas("c3", Form("E/E_exp comparison for %s", mod.name.Data()), 1200, 800);
    mod_E_ratio->SetLineColor(kBlue);
    mod_E_ratio->SetLineWidth(2);
    mod_E_ratio->Draw("HIST");

    TCanvas *c4 = new TCanvas("c4", Form("E_rec/E_exp map for %s", mod.name.Data()), 800, 700);
    gStyle->SetPaintTextFormat(".3f");
    mod_E_xy->SetMarkerSize(2.0);
    mod_E_xy->Draw("colz text");

    TCanvas *c5 = new TCanvas("c5", Form("Hit distribution %s", mod.name.Data()), 800, 700);
    mod_hit_xy->Draw("colz");

    TCanvas *c6 = new TCanvas("c6", "Hit distribution for all events", 800, 700);
    hit_all->Draw("colz");

    TCanvas *c8 = new TCanvas("c8", "Hit distribution for all events (GEM projection)", 800, 700);
    hit_all_gem->Draw("colz");

    TCanvas *c9 = new TCanvas("c9", "Delta X between HyCal and GEM", 800, 700);
    mod_delta_x->Draw("colz");

    TCanvas *c10 = new TCanvas("c10", "Delta Y between HyCal and GEM", 800, 700);
    mod_delta_y->Draw("colz");

    TCanvas *c7 = new TCanvas("c7", Form("E_rec/E_exp map for %s (GEM projection)", mod.name.Data()), 800, 700);
    gStyle->SetPaintTextFormat(".3f");
    mod_E_xy_gem->SetMarkerSize(2.0);
    mod_E_xy_gem->Draw("colz text");

    TCanvas *c11 = new TCanvas("c11", "Delta X between HyCal and GEM (projected)", 800, 700);
    mod_delta_x_proj->Draw("HIST");

    TCanvas *c12 = new TCanvas("c12", "Delta Y between HyCal and GEM (projected)", 800, 700);
    mod_delta_y_proj->Draw("HIST");

    TCanvas *c13 = new TCanvas("c13", "Mean Delta X vs x_d", 800, 700);
    h1_delta_x->Draw("E1");

    TCanvas *c14 = new TCanvas("c14", "Mean Delta Y vs y_d", 800, 700);
    h1_delta_y->Draw("E1");

}