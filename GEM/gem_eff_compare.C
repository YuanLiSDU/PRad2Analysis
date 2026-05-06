
// gem_eff_compare.C
// Compare two run configurations:
//   config1 = first .root argument, config2 = second .root argument
// Draws:
//   1. 2-match efficiency per-chamber comparison (TGraphErrors overlay)
//   2. 2-match efficiency 2D maps (config1 top row, config2 bottom row)
//   3. Inter-layer DeltaX 1D histograms per GEM, both configs overlaid (normalized)
//   4. Inter-layer DeltaY 1D histograms per GEM, both configs overlaid (normalized)
//
// Usage:
//   root -l -q 'gem_eff_compare.C' config1.root config2.root

#include "../EventData.h"
#include "../PhysicsTools.h"
#include <algorithm>

// ---- Geometry (same as gem_eff.C) ----------------------------------------
static const float gz_c[4] = {5817.0, 5853.8, 5414.0, 5458.0};

static const float cmp_x_lo[4]  = {-280., -10., -280., -10.};
static const float cmp_x_hi[4]  = {  10., 280.,   10., 280.};
static const float cmp_y_lo     = -280., cmp_y_hi = 280.;
static const int   cmp_x_bins   = 80,   cmp_y_bins = 160;
static const float cmp_dxy_lo   = -30., cmp_dxy_hi = 30.;
static const int   cmp_dxy_bins = 300;
static const float Ebeam_cmp    = 3488.43f;

static void TrToGEM(float &x, float &y, float &z, int id) {
    float r = gz_c[id] / z;
    x = x * r; y = y * r; z = gz_c[id];
}

// ---------------------------------------------------------------------------
// Fill per-file histograms needed for the comparison:
//   h_2sh[4], h_2mh[4]  – 2-match should-hit / matched-hit (2D)
//   h_dx[4],  h_dy[4]   – inter-layer DeltaX / DeltaY (1D)
// ---------------------------------------------------------------------------
static void fillCompareHists(const TString &fname,
                              TH2F *h_2sh[4], TH2F *h_2mh[4],
                              TH1F *h_dx[4],  TH1F *h_dy[4])
{
    TChain chain("recon");
    chain.Add(fname);
    ReconEventData data;
    setupReconBranches(&chain, data);

    auto isGEM0Dead = [](float x, float y) -> bool {
        return (x >= -170.f && x <= -150.f && y < 0.f);
    };

    int N = chain.GetEntries();
    std::cout << "Processing " << fname << " (" << N << " entries) ..." << std::endl;

    for (int i = 0; i < N; i++) {
        chain.GetEntry(i);
        if (i % 10000 == 0) std::cout << "  " << i << "/" << N << "\r" << std::flush;

        for (int j = 0; j < data.n_clusters; j++) {
            if (fabs(data.cl_energy[j] - Ebeam_cmp) > 250.) continue;
            if (fabs(data.cl_x[j]) < 20.75*2.5 && fabs(data.cl_y[j]) < 20.75*2.5) continue;

            float cx = data.cl_x[j], cy = data.cl_y[j], cz = data.cl_z[j];

            // ---- 2-match efficiency histograms --------------------------------
            // GEM0: require GEM2 or GEM3
            if (data.matchFlag[j] & (1<<2) || data.matchFlag[j] & (1<<3)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 0);
                if (!isGEM0Dead(tx, ty)) {
                    h_2sh[0]->Fill(tx, ty);
                    if (data.matchFlag[j] & (1<<0)) h_2mh[0]->Fill(data.matchGEMx[j][0], data.matchGEMy[j][0]);
                }
            }
            // GEM1: require GEM3 or GEM2
            if (data.matchFlag[j] & (1<<3) || data.matchFlag[j] & (1<<2)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 1);
                h_2sh[1]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<1)) h_2mh[1]->Fill(data.matchGEMx[j][1], data.matchGEMy[j][1]);
            }
            // GEM2: require GEM0 or GEM1
            if (data.matchFlag[j] & (1<<0) || data.matchFlag[j] & (1<<1)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 2);
                h_2sh[2]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<2)) h_2mh[2]->Fill(data.matchGEMx[j][2], data.matchGEMy[j][2]);
            }
            // GEM3: require GEM1 or GEM0
            if (data.matchFlag[j] & (1<<1) || data.matchFlag[j] & (1<<0)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 3);
                h_2sh[3]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<3)) h_2mh[3]->Fill(data.matchGEMx[j][3], data.matchGEMy[j][3]);
            }

            // ---- Inter-layer DeltaX / DeltaY ----------------------------------
            // GEM0 (layer A left): prefer GEM2, fallback GEM3
            if ((data.matchFlag[j] & (1<<0)) &&
                !isGEM0Dead(data.matchGEMx[j][0], data.matchGEMy[j][0])) {
                int p = (data.matchFlag[j] & (1<<2)) ? 2 : (data.matchFlag[j] & (1<<3)) ? 3 : -1;
                if (p >= 0) {
                    h_dx[0]->Fill(data.matchGEMx[j][0] - data.matchGEMx[j][p] * gz_c[0] / gz_c[p]);
                    h_dy[0]->Fill(data.matchGEMy[j][0] - data.matchGEMy[j][p] * gz_c[0] / gz_c[p]);
                }
            }
            // GEM1 (layer A right): prefer GEM3, fallback GEM2
            if (data.matchFlag[j] & (1<<1)) {
                int p = (data.matchFlag[j] & (1<<3)) ? 3 : (data.matchFlag[j] & (1<<2)) ? 2 : -1;
                if (p >= 0) {
                    h_dx[1]->Fill(data.matchGEMx[j][1] - data.matchGEMx[j][p] * gz_c[1] / gz_c[p]);
                    h_dy[1]->Fill(data.matchGEMy[j][1] - data.matchGEMy[j][p] * gz_c[1] / gz_c[p]);
                }
            }
            // GEM2 (layer B left): prefer GEM0, fallback GEM1
            if (data.matchFlag[j] & (1<<2)) {
                int p = (data.matchFlag[j] & (1<<0)) ? 0 : (data.matchFlag[j] & (1<<1)) ? 1 : -1;
                if (p >= 0) {
                    h_dx[2]->Fill(data.matchGEMx[j][2] - data.matchGEMx[j][p] * gz_c[2] / gz_c[p]);
                    h_dy[2]->Fill(data.matchGEMy[j][2] - data.matchGEMy[j][p] * gz_c[2] / gz_c[p]);
                }
            }
            // GEM3 (layer B right): prefer GEM1, fallback GEM0
            if (data.matchFlag[j] & (1<<3)) {
                int p = (data.matchFlag[j] & (1<<1)) ? 1 : (data.matchFlag[j] & (1<<0)) ? 0 : -1;
                if (p >= 0) {
                    h_dx[3]->Fill(data.matchGEMx[j][3] - data.matchGEMx[j][p] * gz_c[3] / gz_c[p]);
                    h_dy[3]->Fill(data.matchGEMy[j][3] - data.matchGEMy[j][p] * gz_c[3] / gz_c[p]);
                }
            }
        }
    }
    std::cout << std::endl;
}

// ===========================================================================
void gem_eff_compare()
{
    // ---- Parse command-line arguments: expect exactly two .root files ------
    TString cfg_file[2];
    int nCfg = 0;
    int argc = gApplication->Argc();
    char **argv = gApplication->Argv();
    for (int i = 1; i < argc && nCfg < 2; i++) {
        TString arg(argv[i]);
        if (arg.EndsWith(".root") && !arg.BeginsWith("-"))
            cfg_file[nCfg++] = arg;
    }
    if (nCfg < 2) {
        std::cerr << "Usage: root -l -q 'gem_eff_compare.C' config1.root config2.root\n";
        return;
    }
    std::cout << "Config1: " << cfg_file[0] << "\n"
              << "Config2: " << cfg_file[1] << "\n";

    const char *cfg_label[2] = {"Config1", "Config2"};
    const int   cfg_color[2] = {kBlue, kRed};
    const int   cfg_marker[2]= {20,    21};

    // ---- Allocate histograms ------------------------------------------------
    TH2F *h_2sh [2][4], *h_2mh [2][4];
    TH2F *h_2eff[2][4];
    TH1F *h_dx  [2][4], *h_dy  [2][4];

    for (int c = 0; c < 2; c++) {
        TString tag = (c == 0) ? "c1" : "c2";
        for (int i = 0; i < 4; i++) {
            h_2sh[c][i]  = new TH2F(Form("h2sh_%s_%d",  tag.Data(), i), "",
                                    cmp_x_bins, cmp_x_lo[i], cmp_x_hi[i],
                                    cmp_y_bins, cmp_y_lo, cmp_y_hi);
            h_2mh[c][i]  = new TH2F(Form("h2mh_%s_%d",  tag.Data(), i), "",
                                    cmp_x_bins, cmp_x_lo[i], cmp_x_hi[i],
                                    cmp_y_bins, cmp_y_lo, cmp_y_hi);
            h_2eff[c][i] = new TH2F(Form("h2eff_%s_%d", tag.Data(), i),
                                    Form("%s  GEM%d 2-Match Efficiency; x (mm); y (mm)", cfg_label[c], i),
                                    cmp_x_bins, cmp_x_lo[i], cmp_x_hi[i],
                                    cmp_y_bins, cmp_y_lo, cmp_y_hi);
            h_dx[c][i]   = new TH1F(Form("hdx_%s_%d",   tag.Data(), i),
                                    Form("%s  GEM%d Inter-layer #DeltaX; #DeltaX (mm); Counts", cfg_label[c], i),
                                    cmp_dxy_bins, cmp_dxy_lo, cmp_dxy_hi);
            h_dy[c][i]   = new TH1F(Form("hdy_%s_%d",   tag.Data(), i),
                                    Form("%s  GEM%d Inter-layer #DeltaY; #DeltaY (mm); Counts", cfg_label[c], i),
                                    cmp_dxy_bins, cmp_dxy_lo, cmp_dxy_hi);
        }
    }

    // ---- Fill ---------------------------------------------------------------
    for (int c = 0; c < 2; c++) {
        fillCompareHists(cfg_file[c], h_2sh[c], h_2mh[c], h_dx[c], h_dy[c]);
        for (int i = 0; i < 4; i++)
            h_2eff[c][i]->Divide(h_2mh[c][i], h_2sh[c][i], 1, 1, "B");
    }

    // ========================================================================
    // Canvas 1: 2-match efficiency per chamber (point comparison)
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_2match_cmp", "GEM 2-Match Efficiency Comparison", 700, 520);
        cv->SetGrid();
        double xpts[4] = {0, 1, 2, 3}, xerr0[4] = {};
        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("GEM 2-Match Efficiency Comparison; GEM chamber; 2-Match Efficiency");
        TLegend *leg = new TLegend(0.12, 0.15, 0.42, 0.35);
        for (int c = 0; c < 2; c++) {
            double ypts[4], yerr[4];
            for (int k = 0; k < 4; k++) {
                double n2s = h_2sh[c][k]->Integral(), n2m = h_2mh[c][k]->Integral();
                ypts[k] = (n2s > 0) ? n2m / n2s : 0.;
                yerr[k] = (n2s > 0) ? sqrt(ypts[k] * (1. - ypts[k]) / n2s) : 0.;
            }
            TGraphErrors *gr = new TGraphErrors(4, xpts, ypts, xerr0, yerr);
            gr->SetName(Form("gr2m_%d", c));
            gr->SetMarkerStyle(cfg_marker[c]); gr->SetMarkerSize(1.5);
            gr->SetMarkerColor(cfg_color[c]);  gr->SetLineColor(cfg_color[c]);
            mg->Add(gr, "PL");
            leg->AddEntry(gr, cfg_label[c], "PL");
        }
        mg->Draw("A");
        mg->GetXaxis()->SetLimits(-0.5, 3.5);
        mg->GetYaxis()->SetRangeUser(0., 1.05);
        // GEM id labels on x-axis
        for (int k = 0; k < 4; k++)
            mg->GetXaxis()->SetBinLabel(mg->GetXaxis()->FindBin(k), Form("GEM%d", k));
        leg->Draw();
        cv->Update();
    }

    // ========================================================================
    // Canvas 2: 2-match efficiency 2D maps
    //   Row 1 (pads 1-4): Config1 GEM0..3
    //   Row 2 (pads 5-8): Config2 GEM0..3
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_2eff_2d", "GEM 2-Match Efficiency 2D", 1600, 800);
        cv->Divide(4, 2);
        for (int c = 0; c < 2; c++) {
            for (int i = 0; i < 4; i++) {
                cv->cd(c * 4 + i + 1);
                h_2eff[c][i]->SetStats(0);
                h_2eff[c][i]->GetZaxis()->SetRangeUser(0., 1.);
                h_2eff[c][i]->Draw("COLZ");
            }
        }
    }

    // ========================================================================
    // Canvas 3: Inter-layer DeltaX per GEM, both configs overlaid (normalized)
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_dx_cmp", "GEM Inter-layer #DeltaX Comparison", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            // Normalize to unit area for shape comparison
            for (int c = 0; c < 2; c++) {
                double s = h_dx[c][i]->Integral();
                if (s > 0) h_dx[c][i]->Scale(1. / s);
                h_dx[c][i]->SetLineColor(cfg_color[c]);
                h_dx[c][i]->SetLineWidth(2);
            }
            h_dx[0][i]->SetTitle(Form("GEM%d  Inter-layer #DeltaX; #DeltaX (mm); Normalized", i));
            double ymax = std::max(h_dx[0][i]->GetMaximum(), h_dx[1][i]->GetMaximum());
            h_dx[0][i]->GetYaxis()->SetRangeUser(0., ymax * 1.2);
            h_dx[0][i]->Draw("HIST");
            h_dx[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.12, 0.75, 0.56, 0.92);
                leg->AddEntry(h_dx[0][i], cfg_label[0], "L");
                leg->AddEntry(h_dx[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 4: Inter-layer DeltaY per GEM, both configs overlaid (normalized)
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_dy_cmp", "GEM Inter-layer #DeltaY Comparison", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            for (int c = 0; c < 2; c++) {
                double s = h_dy[c][i]->Integral();
                if (s > 0) h_dy[c][i]->Scale(1. / s);
                h_dy[c][i]->SetLineColor(cfg_color[c]);
                h_dy[c][i]->SetLineWidth(2);
            }
            h_dy[0][i]->SetTitle(Form("GEM%d  Inter-layer #DeltaY; #DeltaY (mm); Normalized", i));
            double ymax = std::max(h_dy[0][i]->GetMaximum(), h_dy[1][i]->GetMaximum());
            h_dy[0][i]->GetYaxis()->SetRangeUser(0., ymax * 1.2);
            h_dy[0][i]->Draw("HIST");
            h_dy[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.12, 0.75, 0.56, 0.92);
                leg->AddEntry(h_dy[0][i], cfg_label[0], "L");
                leg->AddEntry(h_dy[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }
}
