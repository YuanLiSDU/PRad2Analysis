
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

static const float cmp_x_lo[4]  = {-380., -30., -380., -30.};
static const float cmp_x_hi[4]  = {  30., 380.,   30., 380.};
static const float cmp_y_lo     = -380., cmp_y_hi = 380.;
static const int   cmp_x_bins   = 120,   cmp_y_bins = 240;
static const float cmp_dxy_lo   = -20., cmp_dxy_hi = 20.;
static const int   cmp_dxy_bins = 300;
static const float Ebeam_cmp    = 3488.43f;

static void TrToGEM(float &x, float &y, float &z, int id) {
    float r = gz_c[id] / z;
    x = x * r; y = y * r; z = gz_c[id];
}

// Project GEM hit (gx, gy at gz_gem) onto each HyCal cluster's z-plane and
// return true if the projected point falls within radius_mm of any cluster.
static bool gemNearHyCal(float gx, float gy, float gz_gem,
                          const ReconEventData &data, float radius_mm = 15.f) {
    float r2 = radius_mm * radius_mm;
    for (int j = 0; j < data.n_clusters; j++) {
        float scale = data.cl_z[j] / gz_gem;
        float dx = gx * scale - data.cl_x[j];
        float dy = gy * scale - data.cl_y[j];
        if (dx*dx + dy*dy < r2) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Fill per-file histograms needed for the comparison:
//   h_2sh[4], h_2mh[4]  – 2-match should-hit / matched-hit (2D)
//   h_dx[4],  h_dy[4]   – inter-layer DeltaX / DeltaY (1D)
// ---------------------------------------------------------------------------
static void fillCompareHists(const TString &fname,
                              TH2F *h_2sh[4], TH2F *h_2mh[4],
                              TH1F *h_dx[4],  TH1F *h_dy[4],
                              TH1F *h_nhits[4],
                              TH1F *h_gemhc[4],
                              TH1F *h_mx[4],  TH1F *h_my[4],
                              TH1F *h_chgx[4], TH1F *h_chgy[4],
                              TH1F *h_szx[4],  TH1F *h_szy[4],
                              TH1F *h_asym[4])
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

    for (int i = 0; i < N/2; i++) {
        chain.GetEntry(i);
        if (i % 10000 == 0) std::cout << "  " << i << "/" << N << "\r" << std::flush;

        // ---- n_gem_hits per chamber for 2-match events ------------------
        // A cluster is "2-match" if it is matched to chambers on both layers
        // (layer A: GEM0/GEM1, layer B: GEM2/GEM3).
        bool evt_has_2match = false;
        for (int j = 0; j < data.n_clusters; j++) {
            if (fabs(data.cl_energy[j] - Ebeam_cmp) > 250.) continue;
            if (fabs(data.cl_x[j]) < 20.75*2.5 && fabs(data.cl_y[j]) < 20.75*2.5) continue;
            bool layerA = (data.matchFlag[j] & (1<<0)) || (data.matchFlag[j] & (1<<1));
            bool layerB = (data.matchFlag[j] & (1<<2)) || (data.matchFlag[j] & (1<<3));
            if (layerA && layerB) { evt_has_2match = true; break; }
        }
        if (evt_has_2match) {
            int cnt_hits[4] = {};
            int cnt_hc  [4] = {};
            for (int h = 0; h < data.n_gem_hits; h++) {
                int id = (int)data.det_id[h];
                if (id < 0 || id >= 4) continue;
                cnt_hits[id]++;
                bool near = gemNearHyCal(data.gem_x[h], data.gem_y[h], gz_c[id], data, 15.f);
                if (near) cnt_hc[id]++;

                float qx = data.gem_x_peak[h], qy = data.gem_y_peak[h];
                int   szx  = (int)data.gem_x_size[h];
                int   szy  = (int)data.gem_y_size[h];
                float qsum = qx + qy;
                float asym = (qsum > 0.f) ? fabs(qx - qy) / qsum : 1.f;

                if (h_chgx && h_chgx[id]) h_chgx[id]->Fill(qx);
                if (h_chgy && h_chgy[id]) h_chgy[id]->Fill(qy);
                if (h_szx  && h_szx[id])  h_szx[id]->Fill(szx);
                if (h_szy  && h_szy[id])  h_szy[id]->Fill(szy);
                if (h_asym && h_asym[id]) h_asym[id]->Fill(asym);
            }
            for (int k = 0; k < 4; k++) {
                if (h_nhits && h_nhits[k]) h_nhits[k]->Fill(cnt_hits[k]);
                if (h_gemhc && h_gemhc[k]) h_gemhc[k]->Fill(cnt_hc[k]);
            }
        }

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
                    if (data.matchFlag[j] & (1<<0)) {
                        h_2mh[0]->Fill(data.matchGEMx[j][0], data.matchGEMy[j][0]);
                    }
                }
            }
            // GEM1: require GEM3 or GEM2
            if (data.matchFlag[j] & (1<<3) || data.matchFlag[j] & (1<<2)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 1);
                h_2sh[1]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<1)) {
                    h_2mh[1]->Fill(data.matchGEMx[j][1], data.matchGEMy[j][1]);
                }
            }
            // GEM2: require GEM0 or GEM1
            if (data.matchFlag[j] & (1<<0) || data.matchFlag[j] & (1<<1)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 2);
                h_2sh[2]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<2)) {
                    h_2mh[2]->Fill(data.matchGEMx[j][2], data.matchGEMy[j][2]);
                }
            }
            // GEM3: require GEM1 or GEM0
            if (data.matchFlag[j] & (1<<1) || data.matchFlag[j] & (1<<0)) {
                float tx=cx, ty=cy, tz=cz; TrToGEM(tx, ty, tz, 3);
                h_2sh[3]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<3)) {
                    h_2mh[3]->Fill(data.matchGEMx[j][3], data.matchGEMy[j][3]);
                }
            }
            // ---- Matched x/y 1D positions (any single match) ----------------
            for (int k = 0; k < 4; k++) {
                if (data.matchFlag[j] & (1<<k)) {
                    if (h_mx && h_mx[k]) h_mx[k]->Fill(data.matchGEMx[j][k]);
                    if (h_my && h_my[k]) h_my[k]->Fill(data.matchGEMy[j][k]);
                }
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
    TH2F *h_2sh  [2][4], *h_2mh [2][4];
    TH2F *h_2eff [2][4];
    TH1F *h_dx   [2][4], *h_dy  [2][4];
    TH1F *h_nhits[2][4];
    TH1F *h_gemhc[2][4];
    TH1F *h_mx   [2][4], *h_my      [2][4];
    TH1F *h_chgx [2][4], *h_chgy    [2][4];
    TH1F *h_szx  [2][4], *h_szy     [2][4];
    TH1F *h_asym    [2][4];

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
            h_dx[c][i]    = new TH1F(Form("hdx_%s_%d",    tag.Data(), i),
                                    Form("%s  GEM%d #DeltaX; #DeltaX (mm); Counts", cfg_label[c], i),
                                    cmp_dxy_bins, cmp_dxy_lo, cmp_dxy_hi);
            h_dy[c][i]    = new TH1F(Form("hdy_%s_%d",    tag.Data(), i),
                                    Form("%s  GEM%d #DeltaY; #DeltaY (mm); Counts", cfg_label[c], i),
                                    cmp_dxy_bins, cmp_dxy_lo, cmp_dxy_hi);
            h_nhits[c][i] = new TH1F(Form("hnhits_%s_%d", tag.Data(), i),
                                    Form("%s  GEM%d # hits (2-match events); # GEM hits; Counts", cfg_label[c], i),
                                    100, -0.5, 200.5);
            h_gemhc[c][i] = new TH1F(Form("hgemhc_%s_%d", tag.Data(), i),
                                    Form("%s  GEM%d # matching multiplicity (r<15 mm, 2-match evts); # matched GEM hits; Counts", cfg_label[c], i),
                                    51, -0.5, 50.5);
            h_mx[c][i]    = new TH1F(Form("hmx_%s_%d",    tag.Data(), i),
                                    Form("%s  GEM%d matched x; x (mm); Counts", cfg_label[c], i),
                                    cmp_x_bins, cmp_x_lo[i], cmp_x_hi[i]);
            h_my[c][i]    = new TH1F(Form("hmy_%s_%d",    tag.Data(), i),
                                    Form("%s  GEM%d matched y; y (mm); Counts", cfg_label[c], i),
                                    cmp_y_bins, cmp_y_lo, cmp_y_hi);
            h_chgx[c][i]  = new TH1F(Form("hchgx_%s_%d",  tag.Data(), i),
                                    Form("%s  GEM%d x charge; Qx (ADC); Counts", cfg_label[c], i),
                                    200, 0., 5000.);
            h_chgy[c][i]  = new TH1F(Form("hchgy_%s_%d",  tag.Data(), i),
                                    Form("%s  GEM%d y charge; Qy (ADC); Counts", cfg_label[c], i),
                                    200, 0., 5000.);
            h_szx[c][i]   = new TH1F(Form("hszx_%s_%d",   tag.Data(), i),
                                    Form("%s  GEM%d x cluster size; x_size; Counts", cfg_label[c], i),
                                    20, -0.5, 19.5);
            h_szy[c][i]   = new TH1F(Form("hszy_%s_%d",   tag.Data(), i),
                                    Form("%s  GEM%d y cluster size; y_size; Counts", cfg_label[c], i),
                                    20, -0.5, 19.5);
            h_asym[c][i]     = new TH1F(Form("hasym_%s_%d",    tag.Data(), i),
                                    Form("%s  GEM%d charge asymmetry; |Qx-Qy|/(Qx+Qy); Counts", cfg_label[c], i),
                                    100, 0., 1.);
        }
    }

    // ---- Fill ---------------------------------------------------------------
    for (int c = 0; c < 2; c++) {
        fillCompareHists(cfg_file[c], h_2sh[c], h_2mh[c], h_dx[c], h_dy[c], h_nhits[c], h_gemhc[c],
                         h_mx[c], h_my[c],
                         h_chgx[c], h_chgy[c], h_szx[c], h_szy[c], h_asym[c]);
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
            for (int c = 0; c < 2; c++) {
                h_dx[c][i]->SetLineColor(cfg_color[c]);
                h_dx[c][i]->SetLineWidth(2);
            }
            h_dx[0][i]->SetTitle(Form("GEM%d #DeltaX of 2 GEMs(2 matching events); #DeltaX (mm); Counts", i));
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
                h_dy[c][i]->SetLineColor(cfg_color[c]);
                h_dy[c][i]->SetLineWidth(2);
            }
            h_dy[0][i]->SetTitle(Form("GEM%d #DeltaY of 2 GEMs(2 matching events); #DeltaY (mm); Counts", i));
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

    // ========================================================================
    // Canvas 5: # GEM hits per chamber in 2-match events, both configs overlaid
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_nhits_cmp", "GEM # Hits per Chamber (2-Match Events)", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            for (int c = 0; c < 2; c++) {
                double s = h_nhits[c][i]->Integral();
                if (s > 0) h_nhits[c][i]->Scale(1. / s);
                h_nhits[c][i]->SetLineColor(cfg_color[c]);
                h_nhits[c][i]->SetLineWidth(2);
            }
            h_nhits[0][i]->GetYaxis()->SetRangeUser(1e-5, 1.);
            gPad->SetLogy();
            h_nhits[0][i]->Draw("HIST");
            h_nhits[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.45, 0.75, 0.92, 0.92);
                leg->AddEntry(h_nhits[0][i], cfg_label[0], "L");
                leg->AddEntry(h_nhits[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 6: # GEM hits per chamber matched to a HyCal cluster (r<15mm)
    //           in 2-match events, both configs overlaid (normalized, log Y)
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_gemhc_cmp",
                                  "GEM hits near HyCal (r<15 mm)", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            for (int c = 0; c < 2; c++) {
                double s = h_gemhc[c][i]->Integral();
                if (s > 0) h_gemhc[c][i]->Scale(1. / s);
                h_gemhc[c][i]->SetLineColor(cfg_color[c]);
                h_gemhc[c][i]->SetLineWidth(2);
            }
            h_gemhc[0][i]->SetTitle(Form("GEM%d matching multiplicity (r<15 mm); # matched GEM hits; Counts", i));
            h_gemhc[0][i]->GetYaxis()->SetRangeUser(1e-5, 1.);
            gPad->SetLogy();
            h_gemhc[0][i]->Draw("HIST");
            h_gemhc[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.45, 0.75, 0.92, 0.92);
                leg->AddEntry(h_gemhc[0][i], cfg_label[0], "L");
                leg->AddEntry(h_gemhc[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 7: Matched GEM x per chamber, both configs overlaid
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_mx_cmp", "GEM Matched x Position", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            for (int c = 0; c < 2; c++) {
                h_mx[c][i]->SetLineColor(cfg_color[c]);
                h_mx[c][i]->SetLineWidth(2);
            }
            h_mx[0][i]->SetTitle(Form("GEM%d matched x; x (mm); Counts", i));
            h_mx[0][i]->Draw("HIST");
            h_mx[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.12, 0.75, 0.56, 0.92);
                leg->AddEntry(h_mx[0][i], cfg_label[0], "L");
                leg->AddEntry(h_mx[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 8: Matched GEM y per chamber, both configs overlaid
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_my_cmp", "GEM Matched y Position", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            for (int c = 0; c < 2; c++) {
                h_my[c][i]->SetLineColor(cfg_color[c]);
                h_my[c][i]->SetLineWidth(2);
            }
            h_my[0][i]->SetTitle(Form("GEM%d matched y; y (mm); Counts", i));
            h_my[0][i]->Draw("HIST");
            h_my[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.12, 0.75, 0.56, 0.92);
                leg->AddEntry(h_my[0][i], cfg_label[0], "L");
                leg->AddEntry(h_my[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 9: Hit charge Qx and Qy per chamber, both configs overlaid
    // ========================================================================
    {
        const int xc[2] = {kBlue, kRed};
        const int yc[2] = {kAzure+7, kOrange+7};
        TCanvas *cv = new TCanvas("c_chg_cmp", "GEM Hit seed Peak Charge (2-match evts)", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetLogy();
            for (int c = 0; c < 2; c++) {
                auto norm = [](TH1F *h){ double s=h->Integral(); if(s>0) h->Scale(1./s); };
                norm(h_chgx[c][i]); norm(h_chgy[c][i]);
                h_chgx[c][i]->SetLineColor(xc[c]); h_chgx[c][i]->SetLineWidth(2); h_chgx[c][i]->SetLineStyle(1);
                h_chgy[c][i]->SetLineColor(yc[c]); h_chgy[c][i]->SetLineWidth(2); h_chgy[c][i]->SetLineStyle(2);
            }
            h_chgx[0][i]->SetTitle(Form("GEM%d hit seed peak charge; Charge (ADC); Normalized", i));
            h_chgx[0][i]->GetYaxis()->SetRangeUser(1e-5, 1.);
            h_chgx[0][i]->Draw("HIST");
            h_chgy[0][i]->Draw("HIST SAME");
            h_chgx[1][i]->Draw("HIST SAME");
            h_chgy[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.38, 0.68, 0.95, 0.92);
                leg->AddEntry(h_chgx[0][i], Form("%s Qx", cfg_label[0]), "L");
                leg->AddEntry(h_chgy[0][i], Form("%s Qy", cfg_label[0]), "L");
                leg->AddEntry(h_chgx[1][i], Form("%s Qx", cfg_label[1]), "L");
                leg->AddEntry(h_chgy[1][i], Form("%s Qy", cfg_label[1]), "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 10: Cluster size x_size and y_size per chamber, both configs overlaid
    // ========================================================================
    {
        const int xc[2] = {kBlue, kRed};
        const int yc[2] = {kAzure+7, kOrange+7};
        TCanvas *cv = new TCanvas("c_sz_cmp", "GEM Cluster Size (2-match evts)", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetLogy();
            for (int c = 0; c < 2; c++) {
                auto norm = [](TH1F *h){ double s=h->Integral(); if(s>0) h->Scale(1./s); };
                norm(h_szx[c][i]); norm(h_szy[c][i]);
                h_szx[c][i]->SetLineColor(xc[c]); h_szx[c][i]->SetLineWidth(2); h_szx[c][i]->SetLineStyle(1);
                h_szy[c][i]->SetLineColor(yc[c]); h_szy[c][i]->SetLineWidth(2); h_szy[c][i]->SetLineStyle(2);
            }
            h_szx[0][i]->SetTitle(Form("GEM%d cluster size; size (strips); Normalized", i));
            h_szx[0][i]->GetYaxis()->SetRangeUser(1e-5, 1.);
            h_szx[0][i]->Draw("HIST");
            h_szy[0][i]->Draw("HIST SAME");
            h_szx[1][i]->Draw("HIST SAME");
            h_szy[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.38, 0.68, 0.95, 0.92);
                leg->AddEntry(h_szx[0][i], Form("%s x_size", cfg_label[0]), "L");
                leg->AddEntry(h_szy[0][i], Form("%s y_size", cfg_label[0]), "L");
                leg->AddEntry(h_szx[1][i], Form("%s x_size", cfg_label[1]), "L");
                leg->AddEntry(h_szy[1][i], Form("%s y_size", cfg_label[1]), "L");
                leg->Draw();
            }
        }
    }

    // ========================================================================
    // Canvas 11: Charge asymmetry |Qx-Qy|/(Qx+Qy) per chamber
    //   Noise on one strip projection → asymmetry → peak near 1
    // ========================================================================
    {
        TCanvas *cv = new TCanvas("c_asym_cmp", "GEM Charge Asymmetry (2-match evts)", 1600, 450);
        cv->Divide(4, 1);
        for (int i = 0; i < 4; i++) {
            cv->cd(i + 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetLogy();
            for (int c = 0; c < 2; c++) {
                double s = h_asym[c][i]->Integral();
                if (s > 0) h_asym[c][i]->Scale(1. / s);
                h_asym[c][i]->SetLineColor(cfg_color[c]);
                h_asym[c][i]->SetLineWidth(2);
            }
            h_asym[0][i]->SetTitle(Form("GEM%d charge asymmetry; |Qx-Qy|/(Qx+Qy); Normalized", i));
            h_asym[0][i]->GetYaxis()->SetRangeUser(1e-5, 1.);
            h_asym[0][i]->Draw("HIST");
            h_asym[1][i]->Draw("HIST SAME");
            if (i == 0) {
                TLegend *leg = new TLegend(0.38, 0.75, 0.92, 0.92);
                leg->AddEntry(h_asym[0][i], cfg_label[0], "L");
                leg->AddEntry(h_asym[1][i], cfg_label[1], "L");
                leg->Draw();
            }
        }
    }
}
