
#include "EventData.h"
#include "PhysicsTools.h"

void setupReconBranches(TTree *tree, ReconEventData &ev);

void lowCountMod(){
    TChain *chain = new TChain("recon");
    chain->Add("data/prad_024004_recon.root");
    TTree *t = chain;

    ReconEventData data;
    setupReconBranches(t, data);

    TH1D *mod_counts = new TH1D("h_mod_counts", "Event counts per module; Event counts; Counts", 500, 0, 2e4);

    double mod_event_counts[1156] = {0}; // Assuming module IDs range from 0 to 1155

    int nEntries = t->GetEntries();
    std::cout << "Total entries in recon tree: " << nEntries << std::endl;
    for (int i = 0; i < nEntries; i++) {
        t->GetEntry(i);
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << "\r" << std::flush;
        }
        if( data.n_clusters == 1) {
            int mod_id = data.cl_center[0] - 1000;
            if (mod_id >= 0 && mod_id < 1156) {
                mod_event_counts[mod_id]++;
            }
        }
    }

    for(int mod_id = 0; mod_id < 1156; mod_id++) {
        mod_counts->Fill(mod_event_counts[mod_id]);
    }

    TCanvas *c_mod_counts = new TCanvas("c_mod_counts", "Event Counts per Module", 800, 600);
    mod_counts->Draw();

}

void setupReconBranches(TTree *tree, ReconEventData &ev)
{
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