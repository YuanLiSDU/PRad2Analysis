#include "../EventData.h"
#include "../PhysicsTools.h"

float Ebeam = 2200.51f; // MeV
float resolution = 0.037f;

// Returns true if (x, y) [mm] lies inside the HyCal active acceptance:
// outside the beam hole (2.5 module widths) and inside the outer edge (16 module widths).
bool inHyCal(double xmm, double ymm) {
    const double module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 17. && fabs(ymm) < module * 17.);
}

bool inVeto(double xmm, double ymm) {
    return (xmm > 200) || (xmm < -240) || (ymm > 200) || (ymm < -240);
}

struct VetoHit {
    int channel = -1;
    int npeaks = 0;
    std::vector<float> times;
    std::vector<float> adcs;
};

void veto_ana_2p2(){

    TFile* file = TFile::Open("../data/2.2GeV/prad_025202_filter.root");
    TTree* tree = (TTree*)file->Get("recon");

    ReconEventData ev;
    setupReconBranches(tree, ev);

    TH2F *h2_hit_xy_hycal = new TH2F("h2_hit_xy_hycal", "HyCal Hit XY;X [mm];Y [mm]", 700, -350, 350, 700, -350, 350);
    TH2F *h2_hit_xy_hycal_coin = new TH2F("h2_hit_xy_hycal_coin", "coincidence HyCal Hit XY;X [mm];Y [mm]", 700, -350, 350, 700, -350, 350);

    TH1F *h1_veto_time[4], *h1_veto_adc[4], *h1_veto_hycal_time_diff[4];
    TH1F *h1_veto_time_coinArea[4], *h1_veto_adc_coinArea[4], *h1_veto_hycal_time_diff_coinArea[4];
    for(int i = 0; i < 4; i++) {
        h1_veto_time[i] = new TH1F(Form("h1_veto_time_%d", i), Form("Veto Channel %d Time;Time [ns];Counts", i), 100, 0, 400);
        h1_veto_adc[i] = new TH1F(Form("h1_veto_adc_%d", i), Form("Veto Channel %d ADC;ADC;Counts", i), 100, 0, 20000);
        h1_veto_time_coinArea[i] = new TH1F(Form("h1_veto_time_coinArea_%d", i), Form("Veto Channel %d Time Coincidence Area;Time [ns];Counts", i), 100, 0, 400);
        h1_veto_adc_coinArea[i] = new TH1F(Form("h1_veto_adc_coinArea_%d", i), Form("Veto Channel %d ADC Coincidence Area;ADC;Counts", i), 100, 0, 20000);
        h1_veto_hycal_time_diff[i] = new TH1F(Form("h1_veto_hycal_time_diff_%d", i), Form("Veto Channel %d HyCal Time Difference;Time Difference [ns];Counts", i), 100, -200, 200);
        h1_veto_hycal_time_diff_coinArea[i] = new TH1F(Form("h1_veto_hycal_time_diff_coinArea_%d", i), Form("Veto Channel %d HyCal Time Difference Coincidence Area;Time Difference [ns];Counts", i), 100, -200, 200);
    }

    TH1F *h1_hycal_time = new TH1F("h1_hycal_time", "HyCal Cluster Time;Time [ns];Counts", 100, 0, 400);

    TH2F *h2_hit_xy_veto[4];
    for(int i = 0; i < 4; i++) {
        h2_hit_xy_veto[i] = new TH2F(Form("h2_hit_xy_veto%d", i+1), Form("Veto %d Hit XY;X [mm];Y [mm]", i+1), 700, -350, 350, 700, -350, 350);
    }
    int n_events = tree->GetEntries()/10;

    for(int i = 0; i < n_events; i++) {
        tree->GetEntry(i);
        // Add analysis code here

        if(i % 1000 == 0) {
            std::cout << "Processed " << i << " entries/" << n_events << "\r" << std::flush;
        }

        if(ev.n_clusters != 1) continue;

        float x = ev.cl_x[0];
        float y = ev.cl_y[0];
        if(!inHyCal(x, y)) continue;

        VetoHit veto_hits[4];

        for(int j = 0; j < ev.veto_nch; j++) {
            int channel = ev.veto_id[j]-1;
            int npeaks = ev.veto_npeaks[j];
            veto_hits[channel].channel = channel;
            veto_hits[channel].npeaks = npeaks;
            veto_hits[channel].times.clear();
            veto_hits[channel].adcs.clear();
            for(int k = 0; k < npeaks; k++) {
                veto_hits[channel].times.push_back(ev.veto_peak_time[j][k]);
                veto_hits[channel].adcs.push_back(ev.veto_peak_integral[j][k]);
                float diff = ev.cl_time[0] - ev.veto_peak_time[j][k];
                h1_veto_time[channel]->Fill(ev.veto_peak_time[j][k]);
                h1_veto_adc[channel]->Fill(ev.veto_peak_integral[j][k]);
                h1_veto_hycal_time_diff[channel]->Fill(diff);
            }
        }
        
        h2_hit_xy_hycal->Fill(x, y);
        h1_hycal_time->Fill(ev.cl_time[0]);

        if(inVeto(x, y)){
            for(int j = 0; j < 4; j++) {
                if(veto_hits[j].npeaks > 0) {
                    for(int p = 0; p < veto_hits[j].npeaks; p++) {
                        h1_veto_time_coinArea[j]->Fill(veto_hits[j].times[p]);
                        h1_veto_adc_coinArea[j]->Fill(veto_hits[j].adcs[p]);
                        float diff = ev.cl_time[0] - veto_hits[j].times[p];
                        if(veto_hits[j].adcs[p] > 2000) continue;
                        h1_veto_hycal_time_diff_coinArea[j]->Fill(diff);
                    }
                }
            }
        }

        for(int j = 0; j < 4; j++) {
            if(veto_hits[j].npeaks > 0) {
                h2_hit_xy_veto[j]->Fill(x, y);
            }
        }

    }

    TCanvas* c1 = new TCanvas("c1", "HyCal Hit XY", 800, 600);
    h2_hit_xy_hycal->Draw("COLZ");

    TCanvas* c2 = new TCanvas("c2", "Veto Time", 800, 600);
    c2->cd();
    c2->SetLogy();
    c2->Divide(2,2);
    for(int i = 0; i < 4; i++) {
        c2->cd(i+1);
        gPad->SetLogy();
        h1_veto_time[i]->Draw();
        h1_veto_time_coinArea[i]->SetLineColor(kRed);
        h1_veto_time_coinArea[i]->Draw("SAME");
        TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(h1_veto_time[i], "All", "l");
        leg->AddEntry(h1_veto_time_coinArea[i], "Veto area", "l");
        leg->Draw();
    }

    TCanvas* c3 = new TCanvas("c3", "Veto ADC", 800, 600);
    c3->Divide(2,2);
    c3->SetLogy();
    for(int i = 0; i < 4; i++) {
        c3->cd(i+1);
        gPad->SetLogy();
        h1_veto_adc[i]->Draw();
        h1_veto_adc_coinArea[i]->SetLineColor(kRed);
        h1_veto_adc_coinArea[i]->Draw("SAME");
        TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(h1_veto_adc[i], "All", "l");
        leg->AddEntry(h1_veto_adc_coinArea[i], "Veto area", "l");
        leg->Draw();
    }

    TCanvas* c4 = new TCanvas("c4", "Veto HyCal Time Difference", 800, 600);
    c4->Divide(2,2);
    for(int i = 0; i < 4; i++) {
        c4->cd(i+1);
        h1_veto_hycal_time_diff[i]->Draw();
        h1_veto_hycal_time_diff_coinArea[i]->SetLineColor(kRed);
        h1_veto_hycal_time_diff_coinArea[i]->Draw("SAME");
        TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(h1_veto_hycal_time_diff[i], "All", "l");
        leg->AddEntry(h1_veto_hycal_time_diff_coinArea[i], "Veto area", "l");
        leg->Draw();
    }

    TCanvas* c5 = new TCanvas("c5", "HyCal Time", 800, 600);
    c5->SetLogy();
    h1_hycal_time->Draw();

    TCanvas* c6 = new TCanvas("c6", "Veto Hit XY", 800, 600);
    c6->Divide(2,2);
    for(int i = 0; i < 4; i++) {
        c6->cd(i+1);
        h2_hit_xy_veto[i]->Draw("COLZ");
    }
}
