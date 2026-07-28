#include "../EventData.h"
#include "../PhysicsTools.h"

void gem_hit(){
    
    TFile *f = TFile::Open("/home/liyuan/evviewer/x17_rate_bg/data/3cl/25766/prad_025766_recon_001.root");
    TTree *tree = (TTree*)f->Get("recon");

    ReconEventData ev;
    setupReconBranches(tree, ev);

    TH2F *hit_gem[4], *hit_gem_mott[4], *hit_gem_match[4];
    for(int i=0; i<4; i++){
        hit_gem[i] = new TH2F(Form("hit_gem_%d", i), Form("GEM%d Hit Distribution;X (mm);Y (mm)", i), 1600, -800, 800, 1600, -800, 800);
        hit_gem_mott[i] = new TH2F(Form("hit_gem_mott_%d", i), Form("GEM%d Mott Hit Distribution;X (mm);Y (mm)", i), 1600, -800, 800, 1600, -800, 800);
        hit_gem_match[i] = new TH2F(Form("hit_gem_match_%d", i), Form("GEM%d Hit Match Distribution;X (mm);Y (mm)", i), 1600, -800, 800, 1600, -800, 800);
    }
    TH2F *hit_x_y_charge[4], *hit_x_y_charge_mott[4], *hit_x_y_charge_match[4];
    for(int i=0; i<4; i++){
        hit_x_y_charge[i] = new TH2F(Form("hit_x_y_charge_%d", i), Form("GEM%d Hits Charge;X_charge (ADC);Y_charge (ADC)", i), 500, 0, 10000, 500, 0, 10000);
        hit_x_y_charge_mott[i] = new TH2F(Form("hit_x_y_charge_mott_%d", i), Form("GEM%d Mott Hits Charge;X_charge (ADC);Y_charge (ADC)", i), 500, 0, 10000, 500, 0, 10000);
        hit_x_y_charge_match[i] = new TH2F(Form("hit_x_y_charge_match_%d", i), Form("GEM%d Matched Hits Charge;X_charge (ADC);Y_charge (ADC)", i), 500, 0, 10000, 500, 0, 10000);
    }
    TH2F *hit_x_y_peak[4], *hit_x_y_peak_mott[4], *hit_x_y_peak_match[4];
    for(int i=0; i<4; i++){
        hit_x_y_peak[i] = new TH2F(Form("hit_x_y_peak_%d", i), Form("GEM%d Hits Peak;X_peak (ADC);Y_peak (ADC)", i), 200, 0, 2500, 200, 0, 2500);
        hit_x_y_peak_mott[i] = new TH2F(Form("hit_x_y_peak_mott_%d", i), Form("GEM%d Mott Hits Peak;X_peak (ADC);Y_peak (ADC)", i), 200, 0, 2500, 200, 0, 2500);
        hit_x_y_peak_match[i] = new TH2F(Form("hit_x_y_peak_match_%d", i), Form("GEM%d Matched Hits Peak;X_peak (ADC);Y_peak (ADC)", i), 200, 0, 2500, 200, 0, 2500);
    }
    TH2F *peak_vs_charge_x[4], *peak_vs_charge_y[4], *peak_vs_charge_x_mott[4], *peak_vs_charge_y_mott[4];
    for(int i=0; i<4; i++){
        peak_vs_charge_x[i] = new TH2F(Form("peak_vs_charge_x_%d", i), Form("GEM%d Peak vs Charge (X);X_peak (ADC);X_charge (ADC)", i), 200, 0, 2500, 500, 0, 10000);
        peak_vs_charge_y[i] = new TH2F(Form("peak_vs_charge_y_%d", i), Form("GEM%d Peak vs Charge (Y);Y_peak (ADC);Y_charge (ADC)", i), 200, 0, 2500, 500, 0, 10000);
        peak_vs_charge_x_mott[i] = new TH2F(Form("peak_vs_charge_x_mott_%d", i), Form("GEM%d Peak vs Charge (X) Mott;X_peak (ADC);X_charge (ADC)", i), 200, 0, 2500, 500, 0, 10000);
        peak_vs_charge_y_mott[i] = new TH2F(Form("peak_vs_charge_y_mott_%d", i), Form("GEM%d Peak vs Charge (Y) Mott;Y_peak (ADC);Y_charge (ADC)", i), 200, 0, 2500, 500, 0, 10000);
    }
    TH2F *hit_x_y_time[4], *hit_x_y_time_mott[4], *hit_x_y_time_match[4];
    for(int i=0; i<4; i++){
        hit_x_y_time[i] = new TH2F(Form("hit_x_y_time_%d", i), Form("GEM%d Hits Time;X_time (ns);Y_time (ns)", i), 6, 0, 6, 6, 0, 6);
        hit_x_y_time_mott[i] = new TH2F(Form("hit_x_y_time_mott_%d", i), Form("GEM%d Mott Hits Time;X_time (ns);Y_time (ns)", i), 6, 0, 6, 6, 0, 6);
        hit_x_y_time_match[i] = new TH2F(Form("hit_x_y_time_match_%d", i), Form("GEM%d Matched Hits Time;X_time (ns);Y_time (ns)", i), 6, 0, 6, 6, 0, 6);
    }

    // Build boundary lines from two points for each line.
    const double lower_x1 = 500.0, lower_y1 = 100.0;
    const double lower_x2 = 6000.0, lower_y2 = 6500.0;
    const double upper_x1 = 100.0, upper_y1 = 800.0;
    const double upper_x2 = 5000.0, upper_y2 = 10000.0;

    const double lower_k = (lower_y2 - lower_y1) / (lower_x2 - lower_x1);
    const double lower_b = lower_y1 - lower_k * lower_x1;
    const double upper_k = (upper_y2 - upper_y1) / (upper_x2 - upper_x1);
    const double upper_b = upper_y1 - upper_k * upper_x1;

    cout << "[gem_hit] lower boundary: y = " << lower_k << " * x + " << lower_b
         << " (through (" << lower_x1 << ", " << lower_y1 << ") and ("
         << lower_x2 << ", " << lower_y2 << "))" << endl;
    cout << "[gem_hit] upper boundary: y = " << upper_k << " * x + " << upper_b
         << " (through (" << upper_x1 << ", " << upper_y1 << ") and ("
         << upper_x2 << ", " << upper_y2 << "))" << endl;

    TF1 *charge_lower[4], *charge_upper[4];
    const double fit_x_min[4] = {100.0, 100.0, 100.0, 100.0};
    const double fit_x_max[4] = {6000.0, 6000.0, 6000.0, 6000.0};
    for(int i=0; i<4; i++){
        charge_lower[i] = new TF1(Form("charge_lower_%d", i), "[0]*x+[1]", fit_x_min[i], fit_x_max[i]);
        charge_upper[i] = new TF1(Form("charge_upper_%d", i), "[0]*x+[1]", fit_x_min[i], fit_x_max[i]);
        charge_lower[i]->SetParameters(lower_k, lower_b);
        charge_upper[i]->SetParameters(upper_k, upper_b);
        charge_lower[i]->SetLineColor(kRed+1);
        charge_upper[i]->SetLineColor(kRed+1);
        charge_lower[i]->SetLineWidth(3);
        charge_upper[i]->SetLineWidth(3);
        charge_lower[i]->SetLineStyle(2);
        charge_upper[i]->SetLineStyle(2);
    }

    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);

        if(i % 1000 == 0)
            cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;

        bool is_sum = (ev.trigger_bits & (1u << 8)) != 0;
        bool is_3cl = (ev.trigger_bits & (1u << 11)) != 0;

        bool is_1cl = ev.n_clusters == 1;
        bool is_mott = is_1cl && (ev.cl_energy[0] > 2000 && ev.cl_energy[0] < 2400);
        bool is_gemMatch[4] = {false, false, false, false};
        float matchX[4] = {0.f, 0.f, 0.f, 0.f};
        float matchY[4] = {0.f, 0.f, 0.f, 0.f};
        
        auto charge_match = [&](float gem_x_charge, float gem_y_charge){
            double y_lower = lower_k * gem_x_charge + lower_b;
            double y_upper = upper_k * gem_x_charge + upper_b;
            return (gem_y_charge > y_lower && gem_y_charge < y_upper);
        };

        if(ev.matchNum == 1){
            for(int j=0; j<4; j++){
                if(ev.matchFlag[0] & (1u << j)){
                    is_gemMatch[j] = true;
                    matchX[j] = ev.matchGEMx[0][j];
                    matchY[j] = ev.matchGEMy[0][j];
                }
            }
        }

        if(is_sum) {
            int nGEMHits = ev.n_gem_hits;
            for(int j=0; j<nGEMHits; j++){
                int det_id = ev.det_id[j];
                float gem_x = ev.gem_x[j];
                float gem_y = ev.gem_y[j];
                hit_gem[det_id]->Fill(gem_x, gem_y);
                float gem_x_charge = ev.gem_x_charge[j];
                float gem_y_charge = ev.gem_y_charge[j];
                hit_x_y_charge[det_id]->Fill(gem_x_charge, gem_y_charge);
                float gem_x_peak = ev.gem_x_peak[j];
                float gem_y_peak = ev.gem_y_peak[j];
                hit_x_y_peak[det_id]->Fill(gem_x_peak, gem_y_peak);
                float gem_x_time = ev.gem_x_mTbin[j];
                float gem_y_time = ev.gem_y_mTbin[j];
                hit_x_y_time[det_id]->Fill(gem_x_time, gem_y_time);
                peak_vs_charge_x[det_id]->Fill(gem_x_peak, gem_x_charge);
                peak_vs_charge_y[det_id]->Fill(gem_y_peak, gem_y_charge);
                if(is_mott){
                    peak_vs_charge_x_mott[det_id]->Fill(gem_x_peak, gem_x_charge);
                    peak_vs_charge_y_mott[det_id]->Fill(gem_y_peak, gem_y_charge);
                }
                bool ADC_cut = (gem_x_charge > 220 && gem_y_charge > 220 &&
                    gem_x_peak > 150 && gem_y_peak > 150);

                if(is_mott && ADC_cut){
                    if(is_gemMatch[det_id]){
                        if(matchX[det_id] != 0.f && matchY[det_id] != 0.f){
                            if(fabs(gem_x - matchX[det_id]) < 3 && fabs(gem_y - matchY[det_id]) < 3){
                                // Matched GEM hit
                                hit_gem_mott[det_id]->Fill(gem_x, gem_y);
                                hit_x_y_charge_mott[det_id]->Fill(gem_x_charge, gem_y_charge);
                                hit_x_y_peak_mott[det_id]->Fill(gem_x_peak, gem_y_peak);
                                hit_x_y_time_mott[det_id]->Fill(gem_x_time, gem_y_time);
                            }
                        }
                    }
                }
                if(charge_match(gem_x_charge, gem_y_charge) && is_mott && ADC_cut){
                    hit_x_y_charge_match[det_id]->Fill(gem_x_charge, gem_y_charge);
                    hit_x_y_peak_match[det_id]->Fill(gem_x_peak, gem_y_peak);
                    hit_x_y_time_match[det_id]->Fill(gem_x_time, gem_y_time);
                    hit_gem_match[det_id]->Fill(gem_x, gem_y);
                }
            }
        }

    }

    TCanvas *c = new TCanvas("c", "GEM Hits", 1200, 800);
    c->Divide(2,2);
    for(int i=0; i<4; i++){
        c->cd(i+1);
        hit_gem[i]->Draw("COLZ");
    }

    TCanvas *c2 = new TCanvas("c2", "GEM Hits Charge", 1200, 800);
    c2->Divide(2,2);
    for(int i=0; i<4; i++){
        c2->cd(i+1);
        hit_x_y_charge[i]->Draw("COLZ");
    }

    TCanvas *c3 = new TCanvas("c3", "GEM Hits Peak", 1200, 800);
    c3->Divide(2,2);
    for(int i=0; i<4; i++){
        c3->cd(i+1);
        hit_x_y_peak[i]->Draw("COLZ");
    }

    TCanvas *c4 = new TCanvas("c4", "GEM Mott Hits", 1200, 800);
    c4->Divide(2,2);
    for(int i=0; i<4; i++){
        c4->cd(i+1);
        hit_gem_mott[i]->Draw("COLZ");
    }

    TCanvas *c5 = new TCanvas("c5", "GEM Mott Hits Charge", 1200, 800);
    c5->Divide(2,2);
    for(int i=0; i<4; i++){
        c5->cd(i+1);
        hit_x_y_charge_mott[i]->Draw("COLZ");
        charge_lower[i]->Draw("SAME");
        charge_upper[i]->Draw("SAME");
    }

    TCanvas *c6 = new TCanvas("c6", "GEM Mott Hits Peak", 1200, 800);
    c6->Divide(2,2);
    for(int i=0; i<4; i++){
        c6->cd(i+1);
        hit_x_y_peak_mott[i]->Draw("COLZ");
    }

    TCanvas *c7 = new TCanvas("c7", "GEM Hits Time", 1200, 800);
    c7->Divide(2,2);
    for(int i=0; i<4; i++){
        c7->cd(i+1);
        hit_x_y_time[i]->Draw("COLZ");
    }

    TCanvas *c8 = new TCanvas("c8", "GEM Mott Hits Time", 1200, 800);
    c8->Divide(2,2);
    for(int i=0; i<4; i++){
        c8->cd(i+1);
        hit_x_y_time_mott[i]->Draw("COLZ");
    }

    TCanvas *c9 = new TCanvas("c9", "GEM Hits Charge Match", 1200, 800);
    c9->Divide(2,2);
    for(int i=0; i<4; i++){
        c9->cd(i+1);
        hit_x_y_charge_match[i]->Draw("COLZ");
    }

    TCanvas *c10 = new TCanvas("c10", "GEM Hits Peak Match", 1200, 800);
    c10->Divide(2,2);
    for(int i=0; i<4; i++){
        c10->cd(i+1);
        hit_x_y_peak_match[i]->Draw("COLZ");
    }

    TCanvas *c11 = new TCanvas("c11", "GEM Hits Time Match", 1200, 800);
    c11->Divide(2,2);
    for(int i=0; i<4; i++){
        c11->cd(i+1);
        hit_x_y_time_match[i]->Draw("COLZ");
    }

    TCanvas *c12 = new TCanvas("c12", "GEM Hits Match", 1200, 800);
    c12->Divide(2,2);
    for(int i=0; i<4; i++){
        c12->cd(i+1);
        hit_gem_match[i]->Draw("COLZ");
    }

    TCanvas *c13 = new TCanvas("c13", "GEM Peak vs Charge X", 1200, 800);
    c13->Divide(2,2);
    for(int i=0; i<4; i++){
        c13->cd(i+1);
        peak_vs_charge_x[i]->Draw("COLZ");
    }

    TCanvas *c14 = new TCanvas("c14", "GEM Peak vs Charge Y", 1200, 800);
    c14->Divide(2,2);
    for(int i=0; i<4; i++){
        c14->cd(i+1);
        peak_vs_charge_y[i]->Draw("COLZ");
    }

    TCanvas *c15 = new TCanvas("c15", "GEM Peak vs Charge X Mott", 1200, 800);
    c15->Divide(2,2);
    for(int i=0; i<4; i++){
        c15->cd(i+1);
        peak_vs_charge_x_mott[i]->Draw("COLZ");
    }

    TCanvas *c16 = new TCanvas("c16", "GEM Peak vs Charge Y Mott", 1200, 800);
    c16->Divide(2,2);
    for(int i=0; i<4; i++){
        c16->cd(i+1);
        peak_vs_charge_y_mott[i]->Draw("COLZ");
    }
}