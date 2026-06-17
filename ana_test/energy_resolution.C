#include "../EventData.h"
#include "../PhysicsTools.h"
#include <sys/stat.h>

//shift to hycal local coordinates
float shift_x = 0.87f;
float shift_y = 1.0f;

//anagles bin edges
const int Nbins = 33;
const Double_t binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

// Returns true if (x, y) [mm] lies inside the HyCal active acceptance:
// outside the beam hole (2.5 module widths) and inside the outer edge (16 module widths).
bool inHyCal(double xmm, double ymm) {
    const double module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
};

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

float mod_x[1157] = {};
float mod_y[1157] = {};

void load_modXY(){
    for(int w = 1; w <= 1156; w++){
        ModInfo mod = GetModInfo(w);
        if(mod.found) {
            mod_x[w] = mod.x;
            mod_y[w] = mod.y;
        }
    }
}


void energy_resolution(){
    load_modXY();

    TFile *file_3p5 = TFile::Open("../data/recon/3.5GeV/prad_024917_recon.root");
    TTree *tree_3p5 = (TTree*)file_3p5->Get("recon");

    ReconEventData *evp_3p5 = new ReconEventData();
    ReconEventData &ev_3p5 = *evp_3p5;
    setupReconBranches(tree_3p5, ev_3p5);

    // Histograms for energy distribution in each angle bin for 3.5 GeV data
    // bin width = 4200/420 = 10 MeV
    const float bw = 10.0;
    TH1F *h1_hist_3p5[Nbins];
    TH1F *h1_hist_3p5_center[Nbins];
    for(int i = 0; i < Nbins; i++){
        h1_hist_3p5[i] = new TH1F(Form("h1_hist_3p5_%d", i), 
            Form("Energy for 3.5 GeV events in angle bin %d;Energy (MeV);Counts / %.1f MeV", i, bw), 
            420, 0, 4200);
        h1_hist_3p5_center[i] = new TH1F(Form("h1_hist_3p5_center_%d", i),
            Form("Energy for 3.5 GeV events in angle bin %d (centered);Energy (MeV);Counts / %.1f MeV", i, bw),
            420, 0, 4200);
    }
    TH1F *angles = new TH1F("angles", "Scattering Angle Distribution for 3.5 GeV;Theta (deg);Counts", Nbins, binEdge);
    TH2F *E_vs_theta_3p5 = new TH2F("E_vs_theta_3p5", "Energy vs Scattering Angle for 3.5 GeV;Theta (deg);Energy (MeV)", Nbins, binEdge, 420, 0, 4200);

    TH1F *deltaX = new TH1F("deltaX", "X_{cluster} - X_{projected};#DeltaX (mm);Counts", 3000, -300, 300);
    TH1F *deltaY = new TH1F("deltaY", "Y_{cluster} - Y_{projected};#DeltaY (mm);Counts", 3000, -300, 300);

    for (Long64_t i=0; i<tree_3p5->GetEntries(); i++){
        tree_3p5->GetEntry(i);
        if(i % 10000 == 0)
            cout << "Processing event " << i << " / " << tree_3p5->GetEntries() << "\r" << flush;
        
        for(int j = 0; j < ev_3p5.n_clusters; j++){
            if(ev_3p5.cl_nblocks[j] < 5) continue;
            bool gem_match0 = ((ev_3p5.matchFlag[j] & (1u<<0)) != 0); // matched with GEM0
            bool gem_match1 = ((ev_3p5.matchFlag[j] & (1u<<1)) != 0); // matched with GEM1
            bool gem_match2 = ((ev_3p5.matchFlag[j] & (1u<<2)) != 0); // matched with GEM2
            bool gem_match3 = ((ev_3p5.matchFlag[j] & (1u<<3)) != 0); // matched with GEM3
            if( !((gem_match0 || gem_match1) && (gem_match2 || gem_match3)) ) continue;

            float E = ev_3p5.cl_energy[j];
            float x, y, z;
            if(gem_match2) {
                x = ev_3p5.matchGEMx[j][2];
                y = ev_3p5.matchGEMy[j][2];
                z = ev_3p5.matchGEMz[j][2];
            } else if(gem_match3) {
                x = ev_3p5.matchGEMx[j][3];
                y = ev_3p5.matchGEMy[j][3];
                z = ev_3p5.matchGEMz[j][3];
            }
            float scale = 6270.f / z;
            x *= scale; y *= scale; z = 6270.f;
            deltaX->Fill(x - ev_3p5.cl_x[j]);
            deltaY->Fill(y - ev_3p5.cl_y[j]);
            //x += shift_x; y += shift_y; // apply position shifts to align with expected module center
            if(!inHyCal(x, y)) continue;

            float theta_deg = atan(sqrt(x*x + y*y) / z) * 180. / TMath::Pi();
            int bin = angles->FindBin(theta_deg)-1;
            if(bin < 0 || bin >= Nbins) continue;

            E_vs_theta_3p5->Fill(theta_deg, E);
            h1_hist_3p5[bin]->Fill(E);

            // Also fill the "centered" histogram if the hit is within 0.3 module widths of the expected center
            //ModInfo mod = GetModInfo(ev_3p5.cl_center[j]-1000);
            const float module = 20.75; // mm
            // require hit to be in central 3x3 of a 5x5 grid (|xd|,|yd| < 0.3)
            float xd = (x - mod_x[ev_3p5.cl_center[j]-1000]) / module;
            float yd = (y - mod_y[ev_3p5.cl_center[j]-1000]) / module;
            if (std::abs(xd) >= 0.25f || std::abs(yd) >= 0.25f) continue;
            h1_hist_3p5_center[bin]->Fill(E);
        }
    }

    cout << "\nFinished filling histograms for 3.5 GeV data, total events processed: " << tree_3p5->GetEntries() << endl;
    
    // Create output directory
    const char *outdir = "energy_resolution_plots";
    mkdir(outdir, 0755);

    TCanvas *cETheta = new TCanvas("cETheta", "Energy vs Theta", 800, 600);
    cETheta->SetGrid();
    E_vs_theta_3p5->Draw("COLZ");
    cETheta->SaveAs(Form("%s/E_vs_theta_3p5.png", outdir));
    cout << "Energy vs Theta plot saved to " << outdir << "/E_vs_theta_3p5.png" << endl;

    TCanvas *cDelta = new TCanvas("cDelta", "Position Residuals", 1200, 500);
    cDelta->Divide(2,1);
    cDelta->cd(1);
    deltaX->SetLineColor(kBlue + 1);
    deltaX->SetLineWidth(2);
    deltaX->SetTitle("X_{cluster} - X_{projected};#DeltaX (mm);Counts");
    deltaX->Draw("HIST");
    cDelta->cd(2);
    deltaY->SetLineColor(kBlue + 1);
    deltaY->SetLineWidth(2);
    deltaY->SetTitle("Y_{cluster} - Y_{projected};#DeltaY (mm);Counts");
    deltaY->Draw("HIST");

    gStyle->SetOptFit(1111);

    // Result vectors for summary plot (filled in per-bin loop below)
    std::vector<double> vTheta_ep, vErrTheta_ep, vRatio_ep, vRatioErr_ep, vRes_ep, vResErr_ep;
    std::vector<double> vTheta_ee, vErrTheta_ee, vRatio_ee, vRatioErr_ee, vRes_ee, vResErr_ee;

    TCanvas *c[Nbins];
    for(int bin_idx = 0; bin_idx < Nbins; bin_idx++){
        TH1F *h = h1_hist_3p5_center[bin_idx];
        if(h->GetEntries() < 10) continue;

        // Bin center angle for expected peak calculation
        double theta_center = 0.5 * (binEdge[bin_idx] + binEdge[bin_idx+1]);
        double E_ep_expect = ExpectedEnergy(theta_center, 3485.f, "ep");
        double E_ee_expect = ExpectedEnergy(theta_center, 3485.f, "ee");

        // Fit ep peak: find max in [0.7, 1.3]*E_ep_expect, then fit ±2*0.03/sqrt(E[GeV])*E around it
        TF1 *fit_ep = nullptr;
        if (E_ep_expect > 0 && E_ep_expect < 4200) {
            double hw_ep    = 2.0 * 0.03 / sqrt(E_ep_expect / 1000.) * E_ep_expect;
            int bin_lo_ep = h->FindBin(E_ep_expect - 2.*hw_ep);
            int bin_hi_ep = h->FindBin(E_ep_expect + 2.*hw_ep);
            int peak_bin_ep = bin_lo_ep;
            for (int b = bin_lo_ep; b <= bin_hi_ep; b++)
                if (h->GetBinContent(b) > h->GetBinContent(peak_bin_ep)) peak_bin_ep = b;
            double peak_ep  = h->GetXaxis()->GetBinCenter(peak_bin_ep);
            fit_ep = new TF1(Form("fit_ep_%d", bin_idx), "gaus", peak_ep - hw_ep, peak_ep + hw_ep);
            fit_ep->SetParameters(h->GetBinContent(peak_bin_ep), peak_ep, peak_ep * 0.03);
            h->Fit(fit_ep, "RQ+");
        }

        // Fit ee peak: find max in [0.7, 1.3]*E_ee_expect, then fit ±2*0.03/sqrt(E[GeV])*E around it
        TF1 *fit_ee = nullptr;
        if (E_ee_expect > 0 && E_ee_expect < 4200) {
            double hw_ee    = 2.0 * 0.03 / sqrt(E_ee_expect / 1000.) * E_ee_expect;
            int bin_lo_ee = h->FindBin(E_ee_expect - 4.*hw_ee);
            int bin_hi_ee = h->FindBin(E_ee_expect + 4.*hw_ee);
            int peak_bin_ee = bin_lo_ee;
            for (int b = bin_lo_ee; b <= bin_hi_ee; b++)
                if (h->GetBinContent(b) > h->GetBinContent(peak_bin_ee)) peak_bin_ee = b;
            double peak_ee  = h->GetXaxis()->GetBinCenter(peak_bin_ee);
            fit_ee = new TF1(Form("fit_ee_%d", bin_idx), "gaus", peak_ee - hw_ee, peak_ee + hw_ee);
            fit_ee->SetParameters(h->GetBinContent(peak_bin_ee), peak_ee, peak_ee * 0.03);
            h->Fit(fit_ee, "RQ+");
        }

        c[bin_idx] = new TCanvas(Form("c_erecon_%d", bin_idx), "E_recon", 800, 600);
        c[bin_idx]->SetGrid();
        c[bin_idx]->SetLeftMargin(0.13);
        c[bin_idx]->SetBottomMargin(0.13);

        h->SetStats(0);
        h->SetLineColor(kAzure + 1);
        h->SetLineWidth(2);
        h->SetFillColor(kAzure - 9);
        h->SetFillStyle(1001);
        h->SetTitle(Form("Reconstructed Energy #in [%.3f, %.3f] deg;Energy (MeV);Counts / %.1f MeV",
                         binEdge[bin_idx], binEdge[bin_idx+1], bw));
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
        h->Draw("HIST");

        // Vertical dashed lines for expected energies
        double ymax = h->GetMaximum();
        TLine *line_ep_exp = nullptr, *line_ee_exp = nullptr;
        if (E_ep_expect > 0 && E_ep_expect < 4200) {
            line_ep_exp = new TLine(E_ep_expect, 0, E_ep_expect, ymax);
            line_ep_exp->SetLineColor(kRed);
            line_ep_exp->SetLineWidth(2);
            line_ep_exp->SetLineStyle(7);
            line_ep_exp->Draw("SAME");
        }
        if (E_ee_expect > 0 && E_ee_expect < 4200) {
            line_ee_exp = new TLine(E_ee_expect, 0, E_ee_expect, ymax);
            line_ee_exp->SetLineColor(kGreen + 2);
            line_ee_exp->SetLineWidth(2);
            line_ee_exp->SetLineStyle(7);
            line_ee_exp->Draw("SAME");
        }

        // Build legend
        TLegend *leg = new TLegend(0.14, 0.62, 0.62, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.036);

        if (fit_ep) {
            fit_ep->SetLineColor(kRed);
            fit_ep->SetLineWidth(2);
            fit_ep->Draw("SAME");
            double mu_ep    = fit_ep->GetParameter(1);
            double sigma_ep = fabs(fit_ep->GetParameter(2));
            double res_ep   = (mu_ep > 0) ? sigma_ep / mu_ep * 100. : 0.;
            if (line_ep_exp)
                leg->AddEntry(line_ep_exp, Form("ep expect: %.1f MeV", E_ep_expect), "l");
            leg->AddEntry(fit_ep, Form("ep fit: #mu=%.1f, #sigma=%.1f MeV, #sigma/#mu=%.2f%%/#sqrt{E}", mu_ep, sigma_ep, res_ep*sqrt(mu_ep/1000.)), "l");
            // Store for summary plot
            if (mu_ep > 0) {
                double mu_error  = fit_ep->GetParError(1);
                double sigma_error = fit_ep->GetParError(2);
                double res_error    = res_ep * sqrt(pow(sigma_error/sigma_ep,2) + pow(mu_error/mu_ep,2));
                vTheta_ep.push_back(theta_center);
                vErrTheta_ep.push_back(0.);
                vRatio_ep.push_back(mu_ep / E_ep_expect);
                vRatioErr_ep.push_back(mu_error / E_ep_expect);
                vRes_ep.push_back(res_ep * sqrt(mu_ep / 1000.));
                vResErr_ep.push_back(res_error * sqrt(mu_ep / 1000.));
            }
        }

        if (fit_ee) {
            fit_ee->SetLineColor(kGreen + 2);
            fit_ee->SetLineWidth(2);
            fit_ee->Draw("SAME");
            double mu_ee    = fit_ee->GetParameter(1);
            double sigma_ee = fabs(fit_ee->GetParameter(2));
            double res_ee   = (mu_ee > 0) ? sigma_ee / mu_ee * 100. : 0.;
            if (line_ee_exp)
                leg->AddEntry(line_ee_exp, Form("ee expect: %.1f MeV", E_ee_expect), "l");
            leg->AddEntry(fit_ee, Form("ee fit: #mu=%.1f, #sigma=%.1f MeV, #sigma/#mu=%.2f%%/#sqrt{E}", mu_ee, sigma_ee, res_ee*sqrt(mu_ee/1000.)), "l");
            // Store for summary plot
            if (mu_ee > 0) {
                double mu_error  = fit_ee->GetParError(1);
                double sigma_error = fit_ee->GetParError(2);
                double res_error    = res_ee * sqrt(pow(sigma_error/sigma_ee,2) + pow(mu_error/mu_ee,2));
                vTheta_ee.push_back(theta_center);
                vErrTheta_ee.push_back(0.);
                vRatio_ee.push_back(mu_ee / E_ee_expect);
                vRatioErr_ee.push_back(mu_error / E_ee_expect);
                vRes_ee.push_back(res_ee * sqrt(mu_ee / 1000.));
                vResErr_ee.push_back(res_error  * sqrt(mu_ee / 1000.));
            }
        }

        leg->Draw();

        c[bin_idx]->SaveAs(Form("%s/erecon_bin%02d_%.3f_%.3fdeg.png", outdir, bin_idx,
                       binEdge[bin_idx], binEdge[bin_idx+1]));
    }

    cout << "All plots saved to " << outdir << "/" << endl;

    // ---- Summary plot: E_rec/E_exp (top) and resolution/sqrt(E) (bottom) vs angle ----
    // (results already collected in vTheta_ep/ee etc. during the per-bin loop above)
    int N_ep = (int)vTheta_ep.size();
    int N_ee = (int)vTheta_ee.size();

    TGraphErrors *gr_ratio_ep = new TGraphErrors(N_ep, vTheta_ep.data(), vRatio_ep.data(),
                                                  vErrTheta_ep.data(), vRatioErr_ep.data());
    TGraphErrors *gr_ratio_ee = new TGraphErrors(N_ee, vTheta_ee.data(), vRatio_ee.data(),
                                                  vErrTheta_ee.data(), vRatioErr_ee.data());
    TGraphErrors *gr_res_ep   = new TGraphErrors(N_ep, vTheta_ep.data(), vRes_ep.data(),
                                                  vErrTheta_ep.data(), vResErr_ep.data());
    TGraphErrors *gr_res_ee   = new TGraphErrors(N_ee, vTheta_ee.data(), vRes_ee.data(),
                                                  vErrTheta_ee.data(), vResErr_ee.data());

    auto styleGr = [](TGraphErrors *g, Color_t col, int marker) {
        g->SetMarkerColor(col); g->SetLineColor(col);
        g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
        g->SetLineWidth(1);
    };
    styleGr(gr_ratio_ep, kRed,       20);
    styleGr(gr_ratio_ee, kGreen + 2, 21);
    styleGr(gr_res_ep,   kRed,       20);
    styleGr(gr_res_ee,   kGreen + 2, 21);

    TCanvas *cSum = new TCanvas("cSum", "Energy Resolution Summary", 800, 900);
    cSum->SetLeftMargin(0);

    // --- Top pad ---
    TPad *pad1 = new TPad("pad1", "", 0, 0.45, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.13);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.08);
    pad1->Draw();
    pad1->cd();

    TH1F *frame1 = pad1->DrawFrame(binEdge[0], 0.85, binEdge[Nbins], 1.15);
    frame1->SetTitle(";; E_{rec} / E_{exp}");
    frame1->GetYaxis()->SetTitleSize(0.07);
    frame1->GetYaxis()->SetTitleOffset(0.85);
    frame1->GetYaxis()->SetLabelSize(0.06);
    frame1->GetXaxis()->SetLabelSize(0);
    frame1->GetXaxis()->SetTickLength(0.06);

    TLine *line1 = new TLine(binEdge[0], 1.0, binEdge[Nbins], 1.0);
    line1->SetLineStyle(7); line1->SetLineWidth(2); line1->SetLineColor(kBlack);
    line1->Draw("SAME");

    gr_ratio_ep->Draw("P SAME");
    gr_ratio_ee->Draw("P SAME");

    TLegend *leg1 = new TLegend(0.14, 0.72, 0.45, 0.90);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.065);
    leg1->AddEntry(gr_ratio_ep, "e-p elastic", "p");
    leg1->AddEntry(gr_ratio_ee, "e-e (M#oslash ller)", "p");
    leg1->Draw();

    // --- Bottom pad ---
    cSum->cd();
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.45);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.18);
    pad2->SetLeftMargin(0.13);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    double res_max = 0;
    for (auto v : vRes_ep) if (v > res_max) res_max = v;
    for (auto v : vRes_ee) if (v > res_max) res_max = v;
    res_max = res_max * 1.4;
    if (res_max < 0.01) res_max = 5.0;

    TH1F *frame2 = pad2->DrawFrame(binEdge[0], 0, binEdge[Nbins], res_max);
    frame2->SetTitle(";Scattering Angle (deg);#sigma/E #times #sqrt{E} (% #sqrt{GeV})");
    frame2->GetXaxis()->SetTitleSize(0.08);
    frame2->GetXaxis()->SetLabelSize(0.07);
    frame2->GetYaxis()->SetTitleSize(0.07);
    frame2->GetYaxis()->SetTitleOffset(0.85);
    frame2->GetYaxis()->SetLabelSize(0.06);
    frame2->GetXaxis()->CenterTitle();
    frame2->GetYaxis()->CenterTitle();

    gr_res_ep->Draw("P SAME");
    gr_res_ee->Draw("P SAME");

    TLegend *leg2 = new TLegend(0.14, 0.72, 0.45, 0.92);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.065);
    leg2->AddEntry(gr_res_ep, "e-p elastic", "p");
    leg2->AddEntry(gr_res_ee, "e-e (M#oslash ller)", "p");
    leg2->Draw();

    cSum->SaveAs(Form("%s/energy_resolution_summary.png", outdir));
    cout << "Summary plot saved to " << outdir << "/energy_resolution_summary.png" << endl;

    TFile *f_out = new TFile(Form("%s/energy_resolution_results.root", outdir), "RECREATE");
    for(int i = 0; i < Nbins; i++) {
        auto h = h1_hist_3p5_center[i];
        if (h->GetEntries() > 0) {
            h->Write();
            if (c[i]) c[i]->Write();
        }
    }
    if (cSum) cSum->Write();
    f_out->Close();
}