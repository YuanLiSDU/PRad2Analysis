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
    TFile *f = TFile::Open("../data/calib/prad_024512_recon.root");
    TTree *tree = (TTree*)f->Get("recon");

    ReconEventData *evp = new ReconEventData();
    ReconEventData &ev = *evp;
    setupReconBranches(tree, ev);

    TH2F *hit_all = new TH2F("hit_all", "Hit distribution for all events;X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *h2 = new TH2F("h2", "Energy vs. Angle;Energy (MeV);Scattering Angle (degrees)", 100, 0, 5, 4200/2, 0, 4200);
    TH1F *mod_E = new TH1F("h1",
        Form("Energy for %s;Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200);
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
        Form("E_{rec}/E_{exp} per 5#times5 Block in %s;X (mm);Y (mm);E_{rec}/E_{exp}", mod.name.Data()),
        5, x_lo, x_hi, 5, y_lo, y_hi);
    TH2F *mod_hit_xy = new TH2F("hit_xy",
        Form("Hit distribution in %s;X (mm);Y (mm)", mod.name.Data()),
        75, x_lo, x_hi, 75, y_lo, y_hi);
    TH1F *mod_x_recon = new TH1F("x_recon", Form("Reconstructed X in %s;X (mm);Counts", mod.name.Data()), 100, x_lo-5, x_hi+5);

    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
        
        for(int j = 0; j < ev.n_clusters; j++){
            if(ev.cl_nblocks[j] < 3) continue; // only consider clusters with 3 or more blocks to reduce noise
            
            ev.cl_x[j] += shift_x; // apply position shifts to align with expected module center
            ev.cl_y[j] += shift_y;

            if(ev.cl_nblocks[j] > 5) hit_all->Fill(ev.cl_x[j], ev.cl_y[j]);
            
            float theta_deg = atan(sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]) / ev.cl_z[j]) * 180. / TMath::Pi();
            h2->Fill(theta_deg, ev.cl_energy[j]);

            if(ev.cl_x[j] > x_lo && ev.cl_x[j] < x_hi &&
               ev.cl_y[j] > y_lo && ev.cl_y[j] < y_hi && ev.n_clusters == 1){ // focus on single-cluster events in module 635
                mod_E->Fill(ev.cl_energy[j]);
                float E_expected = ExpectedEnergy(theta_deg, 3488.43f, "ep");
                if(E_expected > 0){
                    mod_E_ratio->Fill(ev.cl_energy[j] / E_expected);
                    mod_E_xy->Fill(ev.cl_x[j], ev.cl_y[j], ev.cl_energy[j] / E_expected);
                }
                int bx = mod_E_xy->GetXaxis()->FindBin(ev.cl_x[j]) - 1;
                int by = mod_E_xy->GetYaxis()->FindBin(ev.cl_y[j]) - 1;
                if(bx >= 0 && bx < 5 && by >= 0 && by < 5)
                    mod_E_bloclk[bx][by]->Fill(ev.cl_energy[j]);
                mod_hit_xy->Fill(ev.cl_x[j], ev.cl_y[j]);
                mod_x_recon->Fill(ev.cl_x[j]);
            }
        }   
    }

    // ── File 2: prad_024512_recon_new.root ───────────────────────────────────
    TFile *f2 = TFile::Open("../data/calib/prad_024512_recon_new.root");
    TTree *tree2 = (TTree*)f2->Get("recon");

    ReconEventData *evp2 = new ReconEventData();
    ReconEventData &ev2 = *evp2;
    setupReconBranches(tree2, ev2);

    TH1F *mod_E_new = new TH1F("h1_new",
        Form("Energy for %s (new);Energy (MeV);Counts", mod.name.Data()), 420, 0, 4200);
    TH1F *mod_E_ratio_new = new TH1F("h3_new",
        Form("E/E_{exp} for %s (new);E_{rec}/E_{exp};Counts", mod.name.Data()), 500, 0, 2);

    for (Long64_t i=0; i<tree2->GetEntries(); i++){
        tree2->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree2->GetEntries() << "\r" << flush;

        for(int j = 0; j < ev2.n_clusters; j++){
            if(ev2.cl_nblocks[j] < 3) continue;

            ev2.cl_x[j] += shift_x;
            ev2.cl_y[j] += shift_y;

            if(ev2.cl_x[j] > x_lo && ev2.cl_x[j] < x_hi &&
               ev2.cl_y[j] > y_lo && ev2.cl_y[j] < y_hi && ev2.n_clusters == 1){
                mod_E_new->Fill(ev2.cl_energy[j]);
                float theta_deg = atan(sqrt(ev2.cl_x[j]*ev2.cl_x[j] + ev2.cl_y[j]*ev2.cl_y[j]) / ev2.cl_z[j]) * 180. / TMath::Pi();
                float E_expected = ExpectedEnergy(theta_deg, 3488.43f, "ep");
                if(E_expected > 0)
                    mod_E_ratio_new->Fill(ev2.cl_energy[j] / E_expected);
            }
        }
    }

    // ── Drawing ───────────────────────────────────────────────────────────────
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    h2->Draw("colz");

    // Energy comparison
    TCanvas *c2 = new TCanvas("c2", Form("Energy comparison for %s", mod.name.Data()), 1200, 800);
    mod_E->SetLineColor(kBlue);
    mod_E->SetLineWidth(2);
    mod_E_new->SetLineColor(kRed);
    mod_E_new->SetLineWidth(2);
    mod_E->Draw("HIST");
    //mod_E_new->Scale(0.9);
    mod_E_new->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.88);
    leg2->AddEntry(mod_E,     "recon (old)", "l");
    leg2->AddEntry(mod_E_new, "recon (new)", "l");
    leg2->Draw();

    // E/E_exp comparison
    TCanvas *c3 = new TCanvas("c3", Form("E/E_exp comparison for %s", mod.name.Data()), 1200, 800);
    mod_E_ratio->SetLineColor(kBlue);
    mod_E_ratio->SetLineWidth(2);
    mod_E_ratio_new->SetLineColor(kRed);
    mod_E_ratio_new->SetLineWidth(2);
    mod_E_ratio->Draw("HIST");
    //mod_E_ratio_new->Scale(0.9);
    mod_E_ratio_new->Draw("HIST SAME");
    TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.88);
    leg3->AddEntry(mod_E_ratio,     "recon (old)", "l");
    leg3->AddEntry(mod_E_ratio_new, "recon (new)", "l");
    leg3->Draw();

    TCanvas *c4 = new TCanvas("c4", Form("E_rec/E_exp map for %s", mod.name.Data()), 800, 700);
    gStyle->SetPaintTextFormat(".3f");
    mod_E_xy->SetMarkerSize(2.0);
    mod_E_xy->Draw("colz text");

    TCanvas *c5 = new TCanvas("c5", Form("Hit distribution %s", mod.name.Data()), 800, 700);
    mod_hit_xy->Draw("colz");

    TCanvas *c6 = new TCanvas("c6", "Hit distribution for all events", 800, 700);
    hit_all->Draw("colz");

    TCanvas *c7 = new TCanvas("c7", Form("Reconstructed X in %s", mod.name.Data()), 800, 700);
    mod_x_recon->Draw("E");

}