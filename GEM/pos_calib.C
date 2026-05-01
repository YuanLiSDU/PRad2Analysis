
#include "../EventData.h"
#include "../PhysicsTools.h"

void setupReconBranches(TTree *tree, ReconEventData &ev);
void fitAndDraw(TH1F* hist, const float fit_range);

void pos_calib()
{   
    TChain *chain = new TChain("recon");
    chain->Add("../data/prad_024004_recon.root");
    TTree *t = chain;

    ReconEventData data;
    setupReconBranches(t, data);

    float gx[4] = {1.5, -1.13, 3.2, 2.37};
    float gy[4] = {0,  1.5, 1.7, 1.66};
    float gz[4] = {0,  0,   0,   0};

    TH2D *ep_pos = new TH2D("h2_ep_pos", "single cluster event distribution (HyCal); x (mm); y (mm)", 100, -350, 350, 100, -350, 350);
    TH2D *all_hit = new TH2D("h2_all_hit", "All cluster distribution (HyCal); x (mm); y (mm)", 100, -350, 350, 100, -350, 350);
    TH2D *hit_Ecut = new TH2D("h2_hit_Ecut", "Cluster distribution with E cut (HyCal); x (mm); y (mm)", 100, -350, 350, 100, -350, 350);
    TH2D *ee_pos = new TH2D("h2_ee_pos", "e-e event distribution (HyCal); x (mm); y (mm)", 100, -350, 350, 100, -350, 350);

    TH2F *h_moller_pos = new TH2F("h2_moller", "Moller event distribution (HyCal); x (mm); y (mm)", 100, -350, 350, 100, -350, 350);
    TH1F *h_moller_phi_diff = new TH1F("h1_moller_phi_diff", "Moller event phi difference distribution; #Delta#phi (degrees); Counts", 40, -15, 15);
    TH1F *h_moller_z = new TH1F("h1_moller_z", "Moller vertex distribution(HyCal); z (mm); Counts", 400, 5700, 6700);
    TH1F *h_moller_x = new TH1F("h1_moller_x", "Moller center X distribution(HyCal); x (mm); Counts", 200, -30, 30);
    TH1F *h_moller_y = new TH1F("h1_moller_y", "Moller center Y distribution(HyCal); y (mm); Counts", 200, -30, 30);
    TH2F *h_moller_center = new TH2F("h2_moller_center", "Moller center distribution(HyCal); x (mm); y (mm)", 100, -200, 200, 100, -200, 200);

    TH2F *g_moller_center[4];
    TH1F *g_moller_z[4];
    TH1F *g_moller_x[4];
    TH1F *g_moller_y[4];
    for (int i = 0; i < 4; i++) {
        g_moller_center[i] = new TH2F(Form("h2_moller_gem%d", i), Form("Moller center distribution (GEM%d); x (mm); y (mm)", i), 100, -200, 200, 100, -200, 200);
        g_moller_z[i] = new TH1F(Form("h1_moller_z_gem%d", i), Form("Moller vertex distribution (GEM%d); z (mm); Counts", i), 250, 5000, 6500);
        g_moller_x[i] = new TH1F(Form("h1_moller_x_gem%d", i), Form("Moller center X distribution (GEM%d); x (mm); Counts", i), 120, -30, 30);
        g_moller_y[i] = new TH1F(Form("h1_moller_y_gem%d", i), Form("Moller center Y distribution (GEM%d); y (mm); Counts", i), 60, -30, 30);
    }

    TH2F *gem_hit_pos[4];
    for (int i = 0; i < 4; i++) {
        gem_hit_pos[i] = new TH2F(Form("h2_gem%d_hits", i), Form("GEM%d hit distribution (total E ~ 2100); x (mm); y (mm)", i), 100, -200, 200, 100, -200, 200);
    }

    TH1F *h_phi_align[4];
    for (int i = 0; i < 4; i++) {
        h_phi_align[i] = new TH1F(Form("h1_phi_align_gem%d", i), Form("Anzimuth alignment for GEM%d; #Delta#phi (degrees); Counts", i), 150, -15, 15);
    }

    MollerData H_mollers, G_mollers[4];
    MollerEvent h_m, g_m[4];
    int index = 0;

    int nEntries = t->GetEntries();
    std::cout << "Total entries in recon tree: " << nEntries << std::endl;
    for (int i = 0; i < nEntries; i++)
    {
        t->GetEntry(i);
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << "\r" << std::flush;
        }

        if (fabs(data.total_energy - 2108.f) < 200.) {
            for(int j = 0; j < data.n_gem_hits; j++) {
                int det_id = data.det_id[j];
                    gem_hit_pos[det_id]->Fill(data.gem_x[j], data.gem_y[j]);
            }
        }

        if(data.n_clusters == 1){
            ep_pos->Fill(data.cl_x[0], data.cl_y[0]);
        }

        for(int j = 0; j < data.n_clusters; j++){
            if(fabs(data.cl_energy[j] - 2108.f) < 200.)
                hit_Ecut->Fill(data.cl_x[j], data.cl_y[j]);
            all_hit->Fill(data.cl_x[j], data.cl_y[j]);
        }

        bool gem_match[2][4];
        if(data.n_clusters != 2) continue; // only look at events with 2 clusters

        h_m.first = DataPoint(data.cl_x[0], data.cl_y[0]+0.4, data.cl_z[0], data.cl_energy[0]);
        h_m.second = DataPoint(data.cl_x[1], data.cl_y[1]+0.4, data.cl_z[1], data.cl_energy[1]);
        
        float Epair = h_m.first.E + h_m.second.E;
        float phi_diff = GetMollerPhiDiff(h_m);

        if( fabs(Epair -2108.f) > 200. || fabs(phi_diff) > 10.0) continue;
        if( fabs(data.cl_x[0]) < 20.75*2.5 && fabs(data.cl_y[0]) < 20.75*2.5) continue; 
        if( fabs(data.cl_x[1]) < 20.75*2.5 && fabs(data.cl_y[1]) < 20.75*2.5) continue;

        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 4; k++){
                gem_match[j][k] = (data.matchFlag[j] & (1 << k)) != 0; // check if cluster j matches with GEM k
            }
        }

        H_mollers.push_back(h_m);
        for (int k = 0; k < 4; k++) {
            if(gem_match[0][k] && gem_match[1][k]) {
                g_m[k].first = DataPoint(data.matchGEMx[0][k]+gx[k], data.matchGEMy[0][k]+gy[k], data.matchGEMz[0][k]+gz[k], index);
                g_m[k].second = DataPoint(data.matchGEMx[1][k]+gx[k], data.matchGEMy[1][k]+gy[k], data.matchGEMz[1][k]+gz[k], index);
                if( fabs(g_m[k].first.x) < 37. && fabs(g_m[k].first.y) < 37.) continue;
                if( fabs(g_m[k].second.x) < 37. && fabs(g_m[k].second.y) < 37.) continue;
                G_mollers[k].push_back(g_m[k]);
            }
        }
        index++;
    }

    //moller analysis and filling histograms
    for (int i = 1; i < H_mollers.size(); i++) {
        auto event = H_mollers[i];
        float z = GetMollerZdistance(event, 2108.f);
        float phi_diff = GetMollerPhiDiff(event);
        h_moller_pos->Fill(event.first.x, event.first.y);
        h_moller_pos->Fill(event.second.x, event.second.y);
        h_moller_z->Fill(z);
        h_moller_phi_diff->Fill(phi_diff);

        auto center = GetMollerCenter(H_mollers[i-1], H_mollers[i]);
        if(center[0] == 0.f && center[1] == 0.f) continue;
        h_moller_x->Fill(center[0]);
        h_moller_y->Fill(center[1]);
        h_moller_center->Fill(center[0], center[1]);
    }
    for (int k = 0; k < 4; k++) {
        for (int i = 1; i < G_mollers[k].size(); i++) {
            auto event = G_mollers[k][i];
            float z = GetMollerZdistance(event, 2108.f);
            float phi_diff = GetMollerPhiDiff(event);
            g_moller_z[k]->Fill(z);

            auto center = GetMollerCenter(G_mollers[k][i-1], G_mollers[k][i]);
            if(center[0] == 0.f && center[1] == 0.f) continue;
            g_moller_center[k]->Fill(center[0], center[1]);
            g_moller_x[k]->Fill(center[0]);
            g_moller_y[k]->Fill(center[1]);
            
            auto h_event = H_mollers[G_mollers[k][i].first.E];
            float hm_x[2], hm_y[2];
            float phi_hycal = atan2(h_event.first.y, h_event.first.x) * 180.f / M_PI / 2. + 
                            atan2(-h_event.second.y, -h_event.second.x) * 180.f / M_PI / 2.;
            float phi_gem = atan2(event.first.y, event.first.x) * 180.f / M_PI / 2. + 
                            atan2(-event.second.y, -event.second.x) * 180.f / M_PI / 2.;
            float delta_phi = phi_gem - phi_hycal;
            // Fold into [-90, 90] to avoid ~180-degree ambiguity from branch cuts
            if (delta_phi > 90.f)       delta_phi -= 180.f;
            else if (delta_phi < -90.f) delta_phi += 180.f;
            h_phi_align[k]->Fill(delta_phi);
        }
    }

    //summarize results and save plots
    cout << "Total Moller events: " << H_mollers.size() << std::endl;
    for (int k = 0; k < 4; k++) {
        cout << "Moller events with GEM" << k+1 << " match: " << G_mollers[k].size() << std::endl;
    }
    /*
    TCanvas *c_moller_pos = new TCanvas("c_moller_pos", "Moller Position", 800, 600);
    h_moller_pos->Draw("COLZ"); c_moller_pos->SaveAs("../data/moller_pos.png");

    TCanvas *c_moller_z = new TCanvas("c_moller_z", "Moller Z", 800, 600);
    h_moller_z->Draw();
    fitAndDraw(h_moller_z, 200.f);
    c_moller_z->SaveAs("../data/moller_z.png");

    TCanvas *c_moller_phi = new TCanvas("c_moller_phi", "Moller Phi Diff", 800, 600);
    h_moller_phi_diff->Draw();
    fitAndDraw(h_moller_phi_diff, 3.f);
    c_moller_phi->SaveAs("../data/moller_phi_diff.png");

    TCanvas *c_moller_x = new TCanvas("c_moller_x", "Moller X", 800, 600);
    h_moller_x->Draw();
    fitAndDraw(h_moller_x, 3.f);
    c_moller_x->SaveAs("../data/moller_x.png");

    TCanvas *c_moller_y = new TCanvas("c_moller_y", "Moller Y", 800, 600);
    h_moller_y->Draw();
    fitAndDraw(h_moller_y, 3.f);
    c_moller_y->SaveAs("../data/moller_y.png");

    TCanvas *c_moller_center = new TCanvas("c_moller_center", "Moller Center", 800, 600);
    h_moller_center->Draw("COLZ"); c_moller_center->SaveAs("../data/moller_center.png");

    for (int k = 0; k < 4; k++) {
        TCanvas *cg1 = new TCanvas(Form("c_gem%d_center", k), Form("GEM%d Moller Center", k), 800, 600);
        g_moller_center[k]->Draw("COLZ"); cg1->SaveAs(Form("../data/gem%d_moller_center.png", k));

        TCanvas *cg2 = new TCanvas(Form("c_gem%d_z", k), Form("GEM%d Moller Z", k), 800, 600);
        g_moller_z[k]->Draw();
        if(k!=0) fitAndDraw(g_moller_z[k], 50.f);
        cg2->SaveAs(Form("../data/gem%d_moller_z.png", k));

        TCanvas *cg3 = new TCanvas(Form("c_gem%d_x", k), Form("GEM%d Moller X", k), 800, 600);
        g_moller_x[k]->Draw();
        if(k!=0) fitAndDraw(g_moller_x[k], 2.f);
        cg3->SaveAs(Form("../data/gem%d_moller_x.png", k));

        TCanvas *cg4 = new TCanvas(Form("c_gem%d_y", k), Form("GEM%d Moller Y", k), 800, 600);
        g_moller_y[k]->Draw();
        if(k!=0) fitAndDraw(g_moller_y[k], 4.f);
        cg4->SaveAs(Form("../data/gem%d_moller_y.png", k));
    }

    for (int i = 0; i < 4; i++) {
        TCanvas *ch = new TCanvas(Form("c_gem%d_hit", i), Form("GEM%d Hit Position", i), 800, 600);
        gem_hit_pos[i]->Draw("COLZ"); 
        ch->SetLogz();
        ch->SaveAs(Form("../data/gem%d_hit_pos.png", i));
    }

    for (int i = 0; i < 4; i++) {
        TCanvas *cp = new TCanvas(Form("c_gem%d_phi_align", i), Form("GEM%d Phi Alignment", i), 800, 600);
        h_phi_align[i]->Draw();
        if(i!=0) fitAndDraw(h_phi_align[i], 1.f);
        cp->SaveAs(Form("../data/gem%d_phi_align.png", i));
    }*/

    TCanvas *c_ep_pos = new TCanvas("c_ep_pos", "e'p Position", 800, 600);
    ep_pos->Draw("COLZ");

    TCanvas *c_all_hit = new TCanvas("c_all_hit", "All Hits", 800, 600);
    all_hit->Draw("COLZ");

    TCanvas *c_hit_Ecut = new TCanvas("c_hit_Ecut", "Hits with E cut", 800, 600);
    hit_Ecut->Draw("COLZ");
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

void fitAndDraw(TH1F* hist, const float fit_range){
    float mean = hist->GetBinCenter(hist->GetMaximumBin());
    hist->Fit("gaus", "rq", "", mean-fit_range, mean+fit_range);
    hist->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.15, 0.85, Form("%.2f mm +- %.2f mm", hist->GetFunction("gaus")->GetParameter(1), hist->GetFunction("gaus")->GetParError(1)));
}