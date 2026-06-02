#include "../EventData.h"
#include "../PhysicsTools.h"

float Ebeam = 3488.43f; // MeV
float resolution = 0.037f;

static TF1 *fitVertexPeak(TH1F *hist, Color_t color)
{
    if (!hist || hist->GetEntries() < 10) return nullptr;

    double peak = hist->GetBinCenter(hist->GetMaximumBin());
    double sigma = hist->GetRMS();
    if (sigma <= 0.) sigma = 200.;

    TF1 *fit = new TF1(Form("fit_%s", hist->GetName()), "gaus", peak - 2. * sigma, peak + 2. * sigma);
    fit->SetLineColor(color);
    fit->SetLineWidth(2);
    fit->SetParameters(hist->GetMaximum(), peak, sigma);

    hist->Fit(fit, "RQ0");
    double mean_fit = fit->GetParameter(1);
    double sigma_fit = fabs(fit->GetParameter(2));
    fit->SetRange(mean_fit - 2. * sigma_fit, mean_fit + 2. * sigma_fit);
    hist->Fit(fit, "RQ0+");

    return fit;
}

void vertex_recon(){

    TFile *file = TFile::Open("../data/calib/prad_024713_recon.root");
    TTree *tree = (TTree*)file->Get("recon");

    ReconEventData ev;
    setupReconBranches(tree, ev);

    TH1F *h_vz  = new TH1F("h_vz",  "Reconstructed Vertex Z of e-p Events at 3.5 GeV;z_{vertex} (mm);Counts", 500, -3500, 6500);
    TH1F *h_dca = new TH1F("h_dca", "Track-to-beam DCA;DCA (mm);Counts", 600, 0, 60);

    TH1F *h_vz_1deg = new TH1F("h_vz_1deg",  "Vertex Z from closest approach (theta~1deg);z_{vtx} (mm);Counts", 500, -3500, 6500);
    TH1F *h_vz_2deg = new TH1F("h_vz_2deg",  "Vertex Z from closest approach (theta~2deg);z_{vtx} (mm);Counts", 500, -3500, 6500);
    TH1F *h_vz_3deg = new TH1F("h_vz_3deg",  "Vertex Z from closest approach (theta~3deg);z_{vtx} (mm);Counts", 500, -3500, 6500);

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if(i%10000 == 0)
            std::cout << "  " << i << " / " << tree->GetEntries() << "\r" << std::flush;

        if(ev.n_clusters != 1) continue;
        if(ev.matchNum != 1) continue;
        if(!isMott(ev.cl_energy[0], Ebeam, resolution)) continue;
        if( fabs(ev.cl_x[0]) < 20.25 * 2.5 && fabs(ev.cl_y[0]) < 20.25 * 2.5 ) continue;
        if( fabs(ev.cl_x[0]) > 20.25 * 16. || fabs(ev.cl_y[0]) > 20.25 * 16. ) continue;

        float x[4], y[4], z[4], E;

        x[0] = 0.f; y[0] = 0.f; z[0] = 0.f; // target center
        x[1] = ev.mHit_gx[0][1]; y[1] = ev.mHit_gy[0][1]; z[1] = ev.mHit_gz[0][1];
        x[2] = ev.mHit_gx[0][0]; y[2] = ev.mHit_gy[0][0]; z[2] = ev.mHit_gz[0][0];
        x[3] = ev.cl_x[0]; y[3] = ev.cl_y[0]; z[3] = ev.cl_z[0];
        E = ev.cl_energy[0];
        float theta = atan2(std::sqrt(x[3]*x[3] + y[3]*y[3]), z[3]) * 180.f / M_PI;

        float vx = 0.f, vy = 0.f, vz = 0.f;

        // Track line from two GEM hits: P(t) = P1 + t * (P2 - P1)
        const float x1 = x[1], y1 = y[1], z1 = z[1];
        const float dx = x[2] - x[1];
        const float dy = y[2] - y[1];
        const float dz = z[2] - z[1];

        // Closest approach to beam axis (x=0, y=0): minimize x(t)^2 + y(t)^2
        const float den = dx * dx + dy * dy;
        if (den < 1e-12f) continue; // nearly parallel to beam axis in x-y projection

        const float t = -(x1 * dx + y1 * dy) / den;
        vx = x1 + t * dx;
        vy = y1 + t * dy;
        vz = z1 + t * dz;

        const float dca = std::sqrt(vx * vx + vy * vy);
        h_vz->Fill(vz);
        h_dca->Fill(dca);

        if (fabs(theta - 1.f) < 0.2f) {
            h_vz_1deg->Fill(vz);
        } else if (fabs(theta - 2.f) < 0.3f) {
            h_vz_2deg->Fill(vz);
        } else if (fabs(theta - 3.f) < 0.5f) {
            h_vz_3deg->Fill(vz);
        }

    }

    TCanvas *c = new TCanvas("c_vertex", "Vertex reconstruction", 1200, 500);
    c->Divide(2, 1);
    h_vz->SetStats(0);
    h_vz->SetLineColor(kBlack);
    h_vz_1deg->SetLineColor(kRed);
    h_vz_2deg->SetLineColor(kGreen+2);
    h_vz_3deg->SetLineColor(kBlue);
    h_vz->SetLineWidth(2.5);
    h_vz_1deg->SetLineWidth(2.5);
    h_vz_2deg->SetLineWidth(2.5);
    h_vz_3deg->SetLineWidth(2.5);

    c->cd(1);
    h_vz->Draw();
    h_vz_1deg->Draw("SAME");
    h_vz_2deg->Draw("SAME");
    h_vz_3deg->Draw("SAME");

    TF1 *fit_all  = fitVertexPeak(h_vz, kBlack);
    TF1 *fit_1deg = fitVertexPeak(h_vz_1deg, kRed);
    TF1 *fit_2deg = fitVertexPeak(h_vz_2deg, kGreen+2);
    TF1 *fit_3deg = fitVertexPeak(h_vz_3deg, kBlue);

    TLegend *leg_vz = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg_vz->SetBorderSize(0);
    leg_vz->AddEntry(h_vz, Form("All (#sigma = %.1f mm)", fit_all ? fit_all->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_1deg, Form("theta ~ 1 deg (#sigma = %.1f mm)", fit_1deg ? fit_1deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_2deg, Form("theta ~ 2 deg (#sigma = %.1f mm)", fit_2deg ? fit_2deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_3deg, Form("theta ~ 3 deg (#sigma = %.1f mm)", fit_3deg ? fit_3deg->GetParameter(2) : 0.), "l");
    leg_vz->Draw();
    c->cd(2); h_dca->Draw();
    c->SaveAs("vertex_recon.png");

}