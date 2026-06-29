#include "../EventData.h"
#include "../PhysicsTools.h"

float Ebeam = 2239.5f; // MeV
float resolution = 0.035f;

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

    TFile *file = TFile::Open("../data/recon/2.2GeV/prad_025202_filtered.root");
    TTree *tree = (TTree*)file->Get("recon");

    ReconEventData ev;
    setupReconBranches(tree, ev);

    TH1F *h_vz  = new TH1F("h_vz",  "Reconstructed Vertex Z of e-p Events at 2.2 GeV;z_{vertex} (mm);Counts", 300, -5500, 6500);
    TH1F *h_dca = new TH1F("h_dca", "Track-to-beam DCA;DCA (mm);Counts", 1200, 0, 120);
    TH2F *h_dca_vs_vz = new TH2F("h_dca_vs_vz", "Track-to-beam DCA vs Vertex Z;z_{vertex} (mm);DCA (mm)", 300, -5500, 6500, 1200, 0, 120);

    TH2F *h_vx_vy = new TH2F("h_vx_vy", "Vertex X vs Vertex Y;x_{vertex} (mm);y_{vertex} (mm)", 500, -50, 50, 500, -50, 50);

    float angle_edge[4] = {0.8f, 1.2f, 2.0f, 3.0f};

    TH1F *h_vz_0p7deg = new TH1F("h_vz_0p7deg",  "Vertex Z from closest approach (theta~0.7deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_1p0deg = new TH1F("h_vz_1p0deg",  "Vertex Z from closest approach (theta~1.0deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_1p6deg = new TH1F("h_vz_1p6deg",  "Vertex Z from closest approach (theta~1.6deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_2p5deg = new TH1F("h_vz_2p5deg",  "Vertex Z from closest approach (theta~2.5deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_3p5deg = new TH1F("h_vz_3p5deg",  "Vertex Z from closest approach (theta~3.5deg);z_{vtx} (mm);Counts", 300, -5500, 6500);

    TH1F *h_vz_0p7deg_proj = new TH1F("h_vz_0p7deg_proj",  "Vertex Z from projected track (theta~0.7deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_1p0deg_proj = new TH1F("h_vz_1p0deg_proj",  "Vertex Z from projected track (theta~1.0deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_1p6deg_proj = new TH1F("h_vz_1p6deg_proj",  "Vertex Z from projected track (theta~1.6deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_2p5deg_proj = new TH1F("h_vz_2p5deg_proj",  "Vertex Z from projected track (theta~2.5deg);z_{vtx} (mm);Counts", 300, -5500, 6500);
    TH1F *h_vz_3p5deg_proj = new TH1F("h_vz_3p5deg_proj",  "Vertex Z from projected track (theta~3.5deg);z_{vtx} (mm);Counts", 300, -5500, 6500);

    TH1F *h_phi_diff = new TH1F("h_phi_diff", "Difference in phi between two GEM hits;#Delta#phi (deg);Counts", 600, -30, 30);

    TH1F *deltaX_gems = new TH1F("deltaX_gems", "Difference in X between two GEM hits;#DeltaX (mm);Counts", 200, -50, 50);
    TH1F *deltaY_gems = new TH1F("deltaY_gems", "Difference in Y between two GEM hits;#DeltaY (mm);Counts", 200, -50, 50);

    for (Long64_t i = 0; i < tree->GetEntries()/10; i++) {
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

        float gem1x = 6270.f / z[1] * x[1];
        float gem1y = 6270.f / z[1] * y[1];
        float gem2x = 6270.f / z[2] * x[2];
        float gem2y = 6270.f / z[2] * y[2];
        
        deltaX_gems->Fill(gem1x - gem2x);
        deltaY_gems->Fill(gem1y - gem2y);

        if (std::abs(gem1x - gem2x) > 5.f || std::abs(gem1y - gem2y) > 5.f) continue;

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

        float dca = std::sqrt(vx * vx + vy * vy);
        dca = 1.f;
        if (dca > 40.f) continue;
        h_vz->Fill(vz);
        h_dca->Fill(dca);
        h_dca_vs_vz->Fill(vz, dca);
        h_vx_vy->Fill(vx, vy);

        if (theta < angle_edge[0] && dca < 40) {
            h_vz_0p7deg->Fill(vz);
        } else if (theta >= angle_edge[0] && theta < angle_edge[1] && dca < 30) {
            h_vz_1p0deg->Fill(vz);
        } else if (theta >= angle_edge[1] && theta < angle_edge[2] && dca < 10) {
            h_vz_1p6deg->Fill(vz);
        } else if (theta >= angle_edge[2] && theta < angle_edge[3] && dca < 10) {
            h_vz_2p5deg->Fill(vz);
        } else if (theta >= angle_edge[3] && dca < 10) {
            h_vz_3p5deg->Fill(vz);
        }

        // Fill phi difference histogram
        float phi1 = atan2(y[1], x[1]) * 180.f / M_PI;
        float phi2 = atan2(y[2], x[2]) * 180.f / M_PI;
        float dphi = phi1 - phi2;
        h_phi_diff->Fill(dphi);

        // radial distance projection to vertex Z recon
        float r1 = std::sqrt(x[1] * x[1] + y[1] * y[1]);
        float r2 = std::sqrt(x[2] * x[2] + y[2] * y[2]);
        
        float z_proj = (r2 * z[1] - r1 * z[2]) / (r2 - r1);

        if(theta < angle_edge[0]) {
            h_vz_0p7deg_proj->Fill(z_proj);
        } else if (theta >= angle_edge[0] && theta < angle_edge[1]) {
            h_vz_1p0deg_proj->Fill(z_proj);
        } else if (theta >= angle_edge[1] && theta < angle_edge[2]) {
            h_vz_1p6deg_proj->Fill(z_proj);
        } else if (theta >= angle_edge[2] && theta < angle_edge[3]) {
            h_vz_2p5deg_proj->Fill(z_proj);
        } else if (theta >= angle_edge[3]) {
            h_vz_3p5deg_proj->Fill(z_proj);
        }
    }

    TCanvas *c = new TCanvas("c_vertex", "Vertex reconstruction", 1200, 500);
    c->Divide(2, 1);
    h_vz->SetStats(0);
    h_vz->SetLineColor(kBlack);
    h_vz_0p7deg->SetLineColor(kRed);
    h_vz_1p0deg->SetLineColor(kGreen+2);
    h_vz_1p6deg->SetLineColor(kBlue);
    h_vz_2p5deg->SetLineColor(kMagenta);
    h_vz_3p5deg->SetLineColor(kCyan);
    h_vz->SetLineWidth(2);
    h_vz_0p7deg->SetLineWidth(2);
    h_vz_1p0deg->SetLineWidth(2);
    h_vz_1p6deg->SetLineWidth(2);
    h_vz_2p5deg->SetLineWidth(2);
    h_vz_3p5deg->SetLineWidth(2);

    c->cd(1);
    h_vz->Draw();
    h_vz_0p7deg->Draw("SAME");
    h_vz_1p0deg->Draw("SAME");
    h_vz_1p6deg->Draw("SAME");
    h_vz_2p5deg->Draw("SAME");
    h_vz_3p5deg->Draw("SAME");

    TF1 *fit_all  = fitVertexPeak(h_vz, kBlack);
    TF1 *fit_0p7deg = fitVertexPeak(h_vz_0p7deg, kRed);
    TF1 *fit_1p0deg = fitVertexPeak(h_vz_1p0deg, kGreen+2);
    TF1 *fit_1p6deg = fitVertexPeak(h_vz_1p6deg, kBlue);
    TF1 *fit_2p5deg = fitVertexPeak(h_vz_2p5deg, kMagenta);
    TF1 *fit_3p5deg = fitVertexPeak(h_vz_3p5deg, kCyan);

    TLegend *leg_vz = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg_vz->SetBorderSize(0);
    leg_vz->AddEntry(h_vz, Form("All (#sigma = %.1f mm)", fit_all ? fit_all->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_0p7deg, Form("theta ~ 0.7 deg (#sigma = %.1f mm)", fit_0p7deg ? fit_0p7deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_1p0deg, Form("theta ~ 1.0 deg (#sigma = %.1f mm)", fit_1p0deg ? fit_1p0deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_1p6deg, Form("theta ~ 1.6 deg (#sigma = %.1f mm)", fit_1p6deg ? fit_1p6deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_2p5deg, Form("theta ~ 2.5 deg (#sigma = %.1f mm)", fit_2p5deg ? fit_2p5deg->GetParameter(2) : 0.), "l");
    leg_vz->AddEntry(h_vz_3p5deg, Form("theta ~ 3.5 deg (#sigma = %.1f mm)", fit_3p5deg ? fit_3p5deg->GetParameter(2) : 0.), "l");
    leg_vz->Draw();
    c->cd(2); h_dca->Draw();
    c->SaveAs("vertex_recon.png");

    TCanvas *c_gems = new TCanvas("c_gems", "GEM hits difference", 1200, 500);
    c_gems->Divide(2, 1);
    c_gems->cd(1); deltaX_gems->Draw();
    c_gems->cd(2); deltaY_gems->Draw();

    TCanvas *c_dca_vs_vz = new TCanvas("c_dca_vs_vz", "DCA vs Vertex Z", 1200, 500);
    h_dca_vs_vz->Draw();
    c_dca_vs_vz->SaveAs("dca_vs_vz.png");

    TCanvas *c_vx_vy = new TCanvas("c_vx_vy", "Vertex X vs Vertex Y", 1200, 500);
    h_vx_vy->Draw("COLZ");
    c_vx_vy->SaveAs("vx_vs_vy.png");

    // Compare the two vertex-Z reconstruction methods in each angle range.
    TH1F *h_vz_closest[5] = {
        h_vz_0p7deg, h_vz_1p0deg, h_vz_1p6deg,
        h_vz_2p5deg, h_vz_3p5deg
    };
    TH1F *h_vz_projected[5] = {
        h_vz_0p7deg_proj, h_vz_1p0deg_proj, h_vz_1p6deg_proj,
        h_vz_2p5deg_proj, h_vz_3p5deg_proj
    };
    const char *angle_label[5] = {
        "#theta < 0.8 deg", "0.8 #leq #theta < 1.2 deg",
        "1.2 #leq #theta < 2.0 deg", "2.0 #leq #theta < 3.0 deg",
        "#theta #geq 3.0 deg"
    };
    const char *angle_tag[5] = {"0p7", "1p0", "1p6", "2p5", "3p5"};

    for (int i = 0; i < 5; ++i) {
        TCanvas *c_compare = new TCanvas(
            Form("c_vz_compare_%s", angle_tag[i]),
            Form("Vertex Z comparison, %s", angle_label[i]), 900, 650);

        h_vz_closest[i]->SetTitle(
            Form("Vertex Z comparison (%s);z_{vtx} (mm);Counts", angle_label[i]));
        h_vz_closest[i]->SetStats(0);
        h_vz_projected[i]->SetStats(0);
        h_vz_closest[i]->SetLineColor(kBlue + 1);
        h_vz_projected[i]->SetLineColor(kRed + 1);
        h_vz_closest[i]->SetLineWidth(2);
        h_vz_projected[i]->SetLineWidth(2);
        const double max_count = h_vz_closest[i]->GetMaximum() >
                                 h_vz_projected[i]->GetMaximum()
                               ? h_vz_closest[i]->GetMaximum()
                               : h_vz_projected[i]->GetMaximum();
        if (max_count > 0.) h_vz_closest[i]->SetMaximum(1.15 * max_count);

        h_vz_closest[i]->Draw("HIST");
        h_vz_projected[i]->Draw("HIST SAME");

        TLegend *leg_compare = new TLegend(0.61, 0.75, 0.88, 0.88);
        leg_compare->SetBorderSize(0);
        leg_compare->AddEntry(h_vz_closest[i], "Closest approach", "l");
        leg_compare->AddEntry(h_vz_projected[i], "Projected track", "l");
        leg_compare->Draw();

        c_compare->SaveAs(Form("vertex_z_compare_%sdeg.png", angle_tag[i]));
    }

    TCanvas *c_phi_diff = new TCanvas("c_phi_diff", "Phi Difference", 1200, 500);
    h_phi_diff->Draw();
    c_phi_diff->SaveAs("phi_diff.png");
}
