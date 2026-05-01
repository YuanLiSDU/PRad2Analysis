
#include "../EventData.h"
#include "../PhysicsTools.h"

float gx[4] = {0, -1.13, 3.2, 2.37};
float gy[4] = {0,  1.5, 1.7, 1.66};
float gz[4] = {5812.0,  5852.0,   5412.4,   5458.9};

void setupReconBranches(TTree *tree, ReconEventData &ev);
void TransformToGEMFrame(float &x, float &y, float &z, int gem_id);

void gem_eff(){
    TChain *chain = new TChain("recon");
    chain->Add("../data/prad_024004_recon.root");
    TTree *t = chain;

    ReconEventData data;
    setupReconBranches(t, data);

    MollerData H_mollers, G_mollers[4];
    MollerEvent h_m, g_m[4];

    TH2F *gem_should_hit[4];
    TH2F *gem_match_hit[4];
    TH2F *gem_efficiency[4];
    float gem_x_range_lo[4] = {-250., -0, -250., -0.};
    float gem_x_range_hi[4] = {0., 250, 0, 250.};
    float eff_x_lo[4] = {-250., 0., -250., 0.};
    float eff_x_hi[4] = {0., 250., 0., 250.};
    float gem_y_range_lo = -250., gem_y_range_hi = 250.; 
    for (int i = 0; i < 4; i++) {
        gem_should_hit[i] = new TH2F(Form("h2_gem%d_should_hit", i), Form("GEM%d Should Hit Distribution; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
        gem_match_hit[i] = new TH2F(Form("h2_gem%d_match_hit", i), Form("GEM%d Matched Hit Distribution; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
        gem_efficiency[i] = new TH2F(Form("h2_gem%d_efficiency", i), Form("GEM%d Efficiency; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
    }

    TH2F *gem_2match_should_hit[4];
    TH2F *gem_2match_match_hit[4];
    TH2F *gem_2match_efficiency[4];
    for (int i = 0; i < 4; i++) {
        gem_2match_should_hit[i] = new TH2F(Form("h2_gem%d_2match_should_hit", i), Form("GEM%d 2 Should Hit Distribution; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
        gem_2match_match_hit[i] = new TH2F(Form("h2_gem%d_2match_match_hit", i), Form("GEM%d 2 Matched Hit Distribution; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
        gem_2match_efficiency[i] = new TH2F(Form("h2_gem%d_2match_efficiency", i), Form("GEM%d 2 Hit Efficiency; x (mm); y (mm)", i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_range_lo, gem_y_range_hi);
    }

    int nEntries = t->GetEntries();
    std::cout << "Total entries in recon tree: " << nEntries << std::endl;
    for (int i = 0; i < nEntries; i++) {
        t->GetEntry(i);
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << "\r" << std::flush;
        }

        for (int j = 0; j < data.n_clusters; j++) {
            if( fabs(data.cl_energy[j] - 2108.f) > 200.) continue; 
            if( fabs(data.cl_x[j]) < 20.75*2.5 && fabs(data.cl_y[j]) < 20.75*2.5) continue; // exclude central region
            for (int k = 0; k < 4; k++) {
                TransformToGEMFrame(data.cl_x[j], data.cl_y[j], data.cl_z[j], k);
                gem_should_hit[k]->Fill(data.cl_x[j], data.cl_y[j]);
                bool gem_match = (data.matchFlag[j] & (1 << k)) != 0; // check if cluster j matches with GEM k
                if(gem_match) 
                    gem_match_hit[k]->Fill(data.matchGEMx[j][k], data.matchGEMy[j][k]);
            }
            // gem 0, require hycal and gem 2 to match
            if( data.matchFlag[j] & ( 1 << 2) ){
                TransformToGEMFrame(data.cl_x[j], data.cl_y[j], data.cl_z[j], 0);
                gem_2match_should_hit[0]->Fill(data.cl_x[j], data.cl_y[j]);
                if( data.matchFlag[j] & ( 1 << 0) ){
                    gem_2match_match_hit[0]->Fill(data.matchGEMx[j][0], data.matchGEMy[j][0]);
                }
            }
            // gem 2, require hycal and gem 0 to match
            if( data.matchFlag[j] & ( 1 << 0) ){
                TransformToGEMFrame(data.cl_x[j], data.cl_y[j], data.cl_z[j], 2);
                gem_2match_should_hit[2]->Fill(data.cl_x[j], data.cl_y[j]);
                if( data.matchFlag[j] & ( 1 << 2) ){
                    gem_2match_match_hit[2]->Fill(data.matchGEMx[j][2], data.matchGEMy[j][2]);
                }
            }
            // gem 1, require hycal and gem 3 to match
            if( data.matchFlag[j] & ( 1 << 3) ){
                TransformToGEMFrame(data.cl_x[j], data.cl_y[j], data.cl_z[j], 1);
                gem_2match_should_hit[1]->Fill(data.cl_x[j], data.cl_y[j]);
                if( data.matchFlag[j] & ( 1 << 1) ){
                    gem_2match_match_hit[1]->Fill(data.matchGEMx[j][1], data.matchGEMy[j][1]);
                }
            }
            // gem 3, require hycal and gem 1 to match
            if( data.matchFlag[j] & ( 1 << 1) ){
                TransformToGEMFrame(data.cl_x[j], data.cl_y[j], data.cl_z[j], 3);
                gem_2match_should_hit[3]->Fill(data.cl_x[j], data.cl_y[j]);
                if( data.matchFlag[j] & ( 1 << 3) ){
                    gem_2match_match_hit[3]->Fill(data.matchGEMx[j][3], data.matchGEMy[j][3]);
                }
            }
        }
    }

    float overall_eff[4];
    for (int i = 0; i < 4; i++) {
        int should_hit = gem_should_hit[i]->Integral();
        int match_hit = gem_match_hit[i]->Integral();
        overall_eff[i] = (should_hit > 0) ? static_cast<float>(match_hit) / should_hit : 0.0;
        std::cout << Form("GEM%d Overall Efficiency: %.2f%% (%d/%d)", i, overall_eff[i]*100, match_hit, should_hit) << std::endl;
    }

    float overall_2match_eff[4];
    for (int i = 0; i < 4; i++) {
        int should_hit = gem_2match_should_hit[i]->Integral();
        int match_hit = gem_2match_match_hit[i]->Integral();
        overall_2match_eff[i] = (should_hit > 0) ? static_cast<float>(match_hit) / should_hit : 0.0;
        std::cout << Form("GEM%d 2-Match Overall Efficiency: %.2f%% (%d/%d)", i, overall_2match_eff[i]*100, match_hit, should_hit) << std::endl;
    }

    TCanvas *c_gem_hit[4];
    for (int i = 0; i < 4; i++) {
        c_gem_hit[i] = new TCanvas(Form("c_gem%d_hit", i), Form("GEM%d Hit Position", i), 800, 600);
        c_gem_hit[i]->Divide(2, 1);
        c_gem_hit[i]->cd(1);
        gem_should_hit[i]->SetLogz();
        gem_should_hit[i]->Draw("COLZ");
        c_gem_hit[i]->cd(2);
        gem_match_hit[i]->SetLogz();
        gem_match_hit[i]->Draw("COLZ");
    }

    TCanvas *c_eff = new TCanvas("c_eff", "GEM Efficiency", 1200, 800);
    c_eff->Divide(2, 2);
    for (int i = 0; i < 4; i++) {
        gem_efficiency[i]->Divide(gem_match_hit[i], gem_should_hit[i], 1, 1, "B");
        c_eff->cd(i+1);
        gem_efficiency[i]->SetStats(0);
        gem_efficiency[i]->Draw("COLZ");
    }

    TCanvas *c_gem_2match_hit[4];
    for (int i = 0; i < 4; i++) {
        c_gem_2match_hit[i] = new TCanvas(Form("c_gem%d_2match_hit", i), Form("GEM%d 2-Match Hit Position", i), 800, 600);
        c_gem_2match_hit[i]->Divide(2, 1);
        c_gem_2match_hit[i]->cd(1);
        gem_2match_should_hit[i]->SetLogz();
        gem_2match_should_hit[i]->Draw("COLZ");
        c_gem_2match_hit[i]->cd(2);
        gem_2match_match_hit[i]->SetLogz();
        gem_2match_match_hit[i]->Draw("COLZ");
    }

    TCanvas *c_2match_eff = new TCanvas("c_2match_eff", "GEM 2-Match Efficiency", 1200, 800);
    c_2match_eff->Divide(2, 2);
    for (int i = 0; i < 4; i++) {
        gem_2match_efficiency[i]->Divide(gem_2match_match_hit[i], gem_2match_should_hit[i], 1, 1, "B");
        c_2match_eff->cd(i+1);
        gem_2match_efficiency[i]->SetStats(0);
        gem_2match_efficiency[i]->Draw("COLZ");
    }
}

void setupReconBranches(TTree *tree, ReconEventData &ev){
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("total_energy", &ev.total_energy);
    // HyCal cluster branches
    tree->SetBranchAddress("n_clusters",   &ev.n_clusters);
    tree->SetBranchAddress("cl_x",         ev.cl_x);
    tree->SetBranchAddress("cl_y",         ev.cl_y);
    tree->SetBranchAddress("cl_z",         ev.cl_z);
    tree->SetBranchAddress("cl_energy",    ev.cl_energy);
    tree->SetBranchAddress("cl_nblocks",   ev.cl_nblocks);
    tree->SetBranchAddress("cl_center",    ev.cl_center);
    tree->SetBranchAddress("cl_flag",      ev.cl_flag);
    // Matching results
    tree->SetBranchAddress("matchFlag",    ev.matchFlag);
    tree->SetBranchAddress("matchGEMx",    ev.matchGEMx);
    tree->SetBranchAddress("matchGEMy",    ev.matchGEMy);
    tree->SetBranchAddress("matchGEMz",    ev.matchGEMz);
    tree->SetBranchAddress("match_num",     &ev.matchNum);
    //quick and simple matching results for quick check
    tree->SetBranchAddress("mHit_E", ev.mHit_E);
    tree->SetBranchAddress("mHit_x", ev.mHit_x);
    tree->SetBranchAddress("mHit_y", ev.mHit_y);
    tree->SetBranchAddress("mHit_z", ev.mHit_z);
    tree->SetBranchAddress("mHit_gx", ev.mHit_gx);
    tree->SetBranchAddress("mHit_gy", ev.mHit_gy);
    tree->SetBranchAddress("mHit_gz", ev.mHit_gz);
    tree->SetBranchAddress("mHit_gid", ev.mHit_gid);
    // GEM part
    tree->SetBranchAddress("n_gem_hits",   &ev.n_gem_hits);
    tree->SetBranchAddress("det_id",       ev.det_id);
    tree->SetBranchAddress("gem_x",        ev.gem_x);
    tree->SetBranchAddress("gem_y",        ev.gem_y);
    tree->SetBranchAddress("gem_x_charge", ev.gem_x_charge);
    tree->SetBranchAddress("gem_y_charge", ev.gem_y_charge);
    tree->SetBranchAddress("gem_x_peak",   ev.gem_x_peak);
    tree->SetBranchAddress("gem_y_peak",   ev.gem_y_peak);
    tree->SetBranchAddress("gem_x_size",   ev.gem_x_size);
    tree->SetBranchAddress("gem_y_size",   ev.gem_y_size);
    tree->SetBranchAddress("gem_x_mTbin",   ev.gem_x_mTbin);
    tree->SetBranchAddress("gem_y_mTbin",   ev.gem_y_mTbin);
};

void TransformToGEMFrame(float &x, float &y, float &z, int gem_id) {
    float ratio = gz[gem_id] / z;
    x = x * ratio;
    y = y * ratio;
    z = gz[gem_id];
}