#include "../EventData.h"
#include "../PhysicsTools.h"

float x_lo = 20.75 * 2.5, x_hi = 20.75 * 16.;
float y_lo = 20.75 * 2.5, y_hi = 20.75 * 16.;

void gem_efficency(){

    // ── File 1: prad_024512_recon.root ───────────────────────────────────────
    TFile *f = TFile::Open("../data/calib/prad_024713_recon.root");
    TTree *tree = (TTree*)f->Get("recon");

    ReconEventData *evp = new ReconEventData();
    ReconEventData &ev = *evp;
    setupReconBranches(tree, ev);

    TH2F *hit_all = new TH2F("hit_all", "Hit distribution for all events;X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *hit_all_cut = new TH2F("hit_all_cut", "Hit distribution(single cluster & 3+ blocks & E ~ 3488 MeV);X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *hit_all_cut_dead = new TH2F("hit_all_cut_dead", "Hit distribution(single cluster & 3+ blocks & E ~ 3488 MeV, exclud dead area);X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *hit_all_cut_gem = new TH2F("hit_all_cut_gem", "Hit distribution(single cluster & 3+ blocks & E ~ 3488 MeV, 2 GEM matching);X (mm);Y (mm)", 700, -350, 350, 700, -350, 350);
    TH2F *h2 = new TH2F("2h", "Energy vs. Angle;Energy (MeV);Scattering Angle (degrees)", 100, 0, 5, 4200/2, 0, 4200);
    TH2D *h2_cut = new TH2D("2h_cut", "Energy vs. Angle with cuts;Energy (MeV);Scattering Angle (degrees)", 100, 0, 5, 4200/2, 0, 4200);

    for (Long64_t i=0; i<tree->GetEntries()/2; i++){
        tree->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;

        hit_all->Fill(ev.cl_x[0], ev.cl_y[0]);

        float theta_deg = atan(sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180. / TMath::Pi();
        h2->Fill(theta_deg, ev.cl_energy[0]);
        
        if(ev.n_clusters != 1) continue; 
        if(ev.cl_nblocks[0] < 3) continue; // only consider clusters with 3 or more blocks to reduce noise
        if(fabs(ev.cl_energy[0] - 3488.43f) > 200.f) continue; // apply loose energy cut to focus on elastic peak

        if( (fabs(ev.cl_x[0]) < x_lo && fabs(ev.cl_y[0]) < y_lo) ||
            fabs(ev.cl_x[0]) > x_hi || fabs(ev.cl_y[0]) > y_hi) continue;
        
        hit_all_cut->Fill(ev.cl_x[0], ev.cl_y[0]);
        h2_cut->Fill(theta_deg, ev.cl_energy[0]);

        if(!(ev.cl_x[0] > -185. && ev.cl_x[0] < -162. &&ev.cl_y[0] < 0.))
            hit_all_cut_dead->Fill(ev.cl_x[0], ev.cl_y[0]);

        if(ev.matchNum == 1){
            float gem_x = ev.mHit_gx[0][0];
            float gem_y = ev.mHit_gy[0][0];
            float gem_z = ev.mHit_gz[0][0];
            int det_id = ev.mHit_gid[0][0];
            GEMHit gem_hit{gem_x, gem_y, gem_z, static_cast<uint8_t>(det_id)};
            projectToHyCal(gem_hit);

            hit_all_cut_gem->Fill(gem_hit.x, gem_hit.y);

        }
    }

    // ── Drawing 

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    hit_all->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    hit_all_cut->Draw("COLZ");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
    hit_all_cut_gem->Draw("COLZ");

    TCanvas *c4 = new TCanvas("c4", "c4", 1200, 800);
    h2->Draw("COLZ");

    TCanvas *c5 = new TCanvas("c5", "c5", 1200, 800);
    h2_cut->Draw("COLZ");

    TCanvas *c6 = new TCanvas("c6", "c6", 800, 800);
    hit_all_cut_dead->Draw("COLZ");

    // binomial error: sigma_p = sqrt(p*(1-p)/n)
    auto binom_err = [](double k, double n) -> double {
        double p = k / n;
        return sqrt(p * (1. - p) / n);
    };

    //calculate efficiency gem efficiency = N_cut_gem / N_cut
    double n_all      = hit_all_cut->GetEntries();
    double n_dead     = hit_all_cut_dead->GetEntries();
    double n_gem      = hit_all_cut_gem->GetEntries();

    double efficiency         = n_gem / n_all;
    double efficiency_cutDead = n_gem / n_dead;

    // Left GEM (x < -30) and Right GEM (x > 30) efficiency
    int nx = hit_all_cut_dead->GetNbinsX();
    int ny = hit_all_cut_dead->GetNbinsY();
    int xbin_left  = hit_all_cut_dead->GetXaxis()->FindBin(-30.0) - 1;
    int xbin_right = hit_all_cut_dead->GetXaxis()->FindBin( 30.0) + 1;

    double left_all  = hit_all_cut_dead->Integral(1, xbin_left, 1, ny);
    double left_gem  = hit_all_cut_gem->Integral(1, xbin_left, 1, ny);
    double right_all = hit_all_cut_dead->Integral(xbin_right, nx, 1, ny);
    double right_gem = hit_all_cut_gem->Integral(xbin_right, nx, 1, ny);

    double eff_left  = left_gem  / left_all;
    double eff_right = right_gem / right_all;

    cout << fixed << setprecision(4);
    cout << endl;
    cout << "GEM efficiency                     = "
         << efficiency * 100         << " +/- " << binom_err(n_gem, n_all)  * 100 << " %" << endl;
    cout << "GEM efficiency (excl. dead area)   = "
         << efficiency_cutDead * 100 << " +/- " << binom_err(n_gem, n_dead) * 100 << " %" << endl;
    cout << "Left  GEM efficiency (x < -30 mm)  = "
         << eff_left  * 100 << " +/- " << binom_err(left_gem,  left_all)  * 100 << " %"
         << "  (" << left_gem  << " / " << left_all  << ")" << endl;
    cout << "Right GEM efficiency (x >  30 mm)  = "
         << eff_right * 100 << " +/- " << binom_err(right_gem, right_all) * 100 << " %"
         << "  (" << right_gem << " / " << right_all << ")" << endl;

}