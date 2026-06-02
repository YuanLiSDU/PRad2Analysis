// compare_types.C
// Read background_output.root from four run types, compute C/(A-B) and D/(A-B)
// ratios for mott_yield and moller_yield, and overlay them on one canvas.
//
// Usage:
//   root -l 'emptyTarget/compare_types.C'
//   root -l 'emptyTarget/compare_types.C("pathA.root","pathB.root","pathC.root","pathD.root")'

void compare_types(
    const char *fileA = "../data/empty_target/typeA_may31_gem.root",
    const char *fileB = "../data/empty_target/typeB_may31_gem.root",
    const char *fileC = "../data/empty_target/typeC_may31_gem.root",
    const char *fileD = "../data/empty_target/typeD_may31_gem.root")
{
    // ── Open files and clone histograms ──────────────────────────────────
    const char *paths[4] = { fileA, fileB, fileC, fileD };
    const char *labels[4] = { "A", "B", "C", "D" };

    TH1F *mott  [4];
    TH1F *moller[4];

    for (int i = 0; i < 4; i++) {
        TFile *f = TFile::Open(paths[i]);
        if (!f || f->IsZombie()) {
            std::cerr << "Cannot open " << paths[i] << std::endl;
            return;
        }
        mott  [i] = (TH1F*)f->Get("mott_yield");
        moller[i] = (TH1F*)f->Get("moller_yield");
        if (!mott[i] || !moller[i]) {
            std::cerr << "mott_yield or moller_yield not found in " << paths[i] << std::endl;
            return;
        }
        mott  [i] = (TH1F*)mott  [i]->Clone(Form("mott_%s",   labels[i]));
        moller[i] = (TH1F*)moller[i]->Clone(Form("moller_%s", labels[i]));
        mott  [i]->SetDirectory(nullptr);
        moller[i]->SetDirectory(nullptr);
        f->Close();
        std::cout << "Loaded " << paths[i] << std::endl;
    }

    // ── Compute denominator: A - B ───────────────────────────────────────
    TH1F *mott_denom   = (TH1F*)mott  [0]->Clone("mott_denom");
    TH1F *moller_denom = (TH1F*)moller[0]->Clone("moller_denom");
    mott_denom  ->Add(mott  [1], -1.);
    moller_denom->Add(moller[1], -1.);

    // ── Compute ratios X/(A-B)  (X = C, D) ──────────────────────────────
    // Divide() propagates errors as: e_ratio = ratio * sqrt((eN/N)^2 + (eD/D)^2)
    const int cfg_color[3] = { kRed, kGreen+2, kViolet };
    const char *ratio_label[3] = { "B/(A-B)", "C/(A-B)", "D/(A-B)" };
    // numerator indices in mott[]/moller[]: B=1, C=2, D=3
    const int num_idx[3] = { 1, 2, 3 };

    TH1F *mott_ratio  [3];
    TH1F *moller_ratio[3];

    for (int i = 0; i < 3; i++) {
        mott_ratio[i]   = (TH1F*)mott  [num_idx[i]]->Clone(Form("mott_ratio_%s",   ratio_label[i]));
        moller_ratio[i] = (TH1F*)moller[num_idx[i]]->Clone(Form("moller_ratio_%s", ratio_label[i]));
        mott_ratio  [i]->Divide(mott  [num_idx[i]], mott_denom,   1., 1., "");
        moller_ratio[i]->Divide(moller[num_idx[i]], moller_denom, 1., 1., "");
        mott_ratio  [i]->SetTitle("Background Ratio of e-p Events;Reconstructed Scattering Angle [deg];N_{background} / N_{signal} (%)");
        moller_ratio[i]->SetTitle("Moller Yield Ratio X/(A-B), 2GEMs matching;Reconstructed Scattering Angle [deg];N_{background} / N_{signal} (%)");
    }

    // ── Fix correlated errors for B/(A-B) ────────────────────────────────
    // B appears in both numerator and denominator A-B, so ROOT Divide() over-
    // counts σ_B.  Correct: σ_f = sqrt(B²σ_A² + A²σ_B²) / (A-B)²
    {
        auto fixBoverAmB = [](TH1F *hratio, TH1F *hA, TH1F *hB) {
            for (int b = 1; b <= hratio->GetNbinsX(); b++) {
                double A = hA->GetBinContent(b), sA = hA->GetBinError(b);
                double B = hB->GetBinContent(b), sB = hB->GetBinError(b);
                double d = A - B;
                if (d == 0.) continue;
                hratio->SetBinError(b, std::sqrt(B*B*sA*sA + A*A*sB*sB) / (d*d));
            }
        };
        fixBoverAmB(mott_ratio  [0], mott  [0], mott  [1]);
        fixBoverAmB(moller_ratio[0], moller[0], moller[1]);
    }
    // Scale to percent
    for (int i = 0; i < 3; i++) {
        mott_ratio  [i]->Scale(100.);
        moller_ratio[i]->Scale(100.);
    }

    // ── Draw ─────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);

    TCanvas *cv = new TCanvas("c_compare", "Type Comparison", 1200, 500);
    cv->Divide(2, 1);

    // Helper to find global y-range across three histograms
    auto yRange = [](TH1F **h, int n, double margin = 0.2) -> std::pair<double,double> {
        double lo = 1e9, hi = -1e9;
        for (int i = 0; i < n; i++) {
            for (int b = 1; b <= h[i]->GetNbinsX(); b++) {
                double v = h[i]->GetBinContent(b);
                double e = h[i]->GetBinError(b);
                if (v == 0.) continue;
                lo = std::min(lo, v - e);
                hi = std::max(hi, v + e);
            }
        }
        double span = hi - lo;
        return { lo - margin * span, hi + margin * span };
    };

    // --- Mott ---
    cv->cd(1);
    gPad->SetLogx();
    auto [mott_lo, mott_hi] = yRange(mott_ratio, 3);
    for (int i = 0; i < 3; i++) {
        mott_ratio[i]->SetLineColor(cfg_color[i]);
        mott_ratio[i]->SetMarkerColor(cfg_color[i]);
        mott_ratio[i]->SetMarkerStyle(20 + i);
        mott_ratio[i]->SetLineWidth(2);
        mott_ratio[i]->SetStats(0);
        if (i == 0) {
            mott_ratio[i]->GetYaxis()->SetRangeUser(-0.2, 2.0);
            mott_ratio[i]->GetXaxis()->SetRangeUser(0.5, 3.8);
            mott_ratio[i]->GetXaxis()->CenterTitle();
            mott_ratio[i]->GetYaxis()->CenterTitle();
            mott_ratio[i]->GetXaxis()->SetTitleSize(0.045);
            mott_ratio[i]->GetYaxis()->SetTitleSize(0.045);
            mott_ratio[i]->GetXaxis()->SetTitleOffset(0.9);
            mott_ratio[i]->GetYaxis()->SetTitleOffset(0.9);
            mott_ratio[i]->GetXaxis()->SetMoreLogLabels();
            mott_ratio[i]->GetXaxis()->SetNoExponent();
            mott_ratio[i]->Draw("E");
        } else {
            mott_ratio[i]->Draw("E SAME");
        }
    }
    // Reference line at 1
    TLine *line_mott = new TLine(mott_ratio[0]->GetXaxis()->GetXmin(), 1.,
                                  mott_ratio[0]->GetXaxis()->GetXmax(), 1.);
    line_mott->SetLineStyle(2); line_mott->SetLineColor(kGray+1); line_mott->Draw();
    TLegend *leg1 = new TLegend(0.65, 0.72, 0.92, 0.90);
    for (int i = 0; i < 3; i++) leg1->AddEntry(mott_ratio[i], ratio_label[i], "PE");
    // leg1 drawn later after bmc overlay

    // --- Moller (double-arm range: 1.6 -- 3.0 deg) ---
    cv->cd(2);
    gPad->SetLogx();
    for (int i = 0; i < 3; i++)
        moller_ratio[i]->GetXaxis()->SetRangeUser(1.6, 3.0);
    auto [mol_lo, mol_hi] = yRange(moller_ratio, 3);
    for (int i = 0; i < 3; i++) {
        moller_ratio[i]->SetLineColor(cfg_color[i]);
        moller_ratio[i]->SetMarkerColor(cfg_color[i]);
        moller_ratio[i]->SetMarkerStyle(20 + i);
        moller_ratio[i]->SetLineWidth(2);
        if (i == 0) {
            moller_ratio[i]->GetYaxis()->SetRangeUser(0., mol_hi);
            moller_ratio[i]->GetXaxis()->CenterTitle();
            moller_ratio[i]->GetYaxis()->CenterTitle();
            moller_ratio[i]->GetXaxis()->SetTitleSize(0.045);
            moller_ratio[i]->GetYaxis()->SetTitleSize(0.045);
            moller_ratio[i]->GetXaxis()->SetTitleOffset(0.9);
            moller_ratio[i]->GetYaxis()->SetTitleOffset(0.9);
            moller_ratio[i]->GetXaxis()->SetMoreLogLabels();
            moller_ratio[i]->GetXaxis()->SetNoExponent();
            moller_ratio[i]->Draw("E");
        } else {
            moller_ratio[i]->Draw("E SAME");
        }
    }
    TLine *line_mol = new TLine(moller_ratio[0]->GetXaxis()->GetXmin(), 1.,
                                 moller_ratio[0]->GetXaxis()->GetXmax(), 1.);
    line_mol->SetLineStyle(2); line_mol->SetLineColor(kGray+1); line_mol->Draw();
    TLegend *leg2 = new TLegend(0.65, 0.72, 0.92, 0.90);
    for (int i = 0; i < 3; i++) leg2->AddEntry(moller_ratio[i], ratio_label[i], "PE");
    // leg2 drawn later after bmc overlay

    int Arun = 24941, Brun = 24942, Crun = 24943, Drun = 24945;

    TLatex *tex = new TLatex();
    tex->SetNDC(1);
    tex->SetTextSize(0.038);
    cv->cd(1);
    //tex->DrawLatex(0.13, 0.86, Form("A:%d  B:%d  C:%d  D:%d", Arun, Brun, Crun, Drun));
    cv->cd(2);
    //tex->DrawLatex(0.13, 0.86, Form("A:%d  B:%d  C:%d  D:%d", Arun, Brun, Crun, Drun));

    // ── Overlay (B-C)/(A-B) and (C-D)/(A-B) on the same pads ────────────
    // Numerator: B - C
    TH1F *mott_bmc   = (TH1F*)mott  [1]->Clone("mott_BminusC");
    TH1F *moller_bmc = (TH1F*)moller[1]->Clone("moller_BminusC");
    mott_bmc  ->Add(mott  [2], -1.);
    moller_bmc->Add(moller[2], -1.);

    // (B-C)/(A-B)
    TH1F *mott_bmc_over_amb   = (TH1F*)mott_bmc  ->Clone("mott_BmC_over_AmB");
    TH1F *moller_bmc_over_amb = (TH1F*)moller_bmc->Clone("moller_BmC_over_AmB");
    mott_bmc_over_amb  ->Divide(mott_bmc,   mott_denom,   1., 1., "");
    moller_bmc_over_amb->Divide(moller_bmc, moller_denom, 1., 1., "");

    // Numerator: C - D
    TH1F *mott_cmd   = (TH1F*)mott  [2]->Clone("mott_CminusD");
    TH1F *moller_cmd = (TH1F*)moller[2]->Clone("moller_CminusD");
    mott_cmd  ->Add(mott  [3], -1.);
    moller_cmd->Add(moller[3], -1.);

    // (C-D)/(A-B)
    TH1F *mott_cmd_over_amb   = (TH1F*)mott_cmd  ->Clone("mott_CmD_over_AmB");
    TH1F *moller_cmd_over_amb = (TH1F*)moller_cmd->Clone("moller_CmD_over_AmB");
    mott_cmd_over_amb  ->Divide(mott_cmd,   mott_denom,   1., 1., "");
    moller_cmd_over_amb->Divide(moller_cmd, moller_denom, 1., 1., "");

    // ── Fix correlated errors for (B-C)/(A-B) ────────────────────────────
    {
        auto fixBmCoverAmB = [](TH1F *hratio, TH1F *hA, TH1F *hB, TH1F *hC) {
            for (int b = 1; b <= hratio->GetNbinsX(); b++) {
                double A = hA->GetBinContent(b), sA = hA->GetBinError(b);
                double B = hB->GetBinContent(b), sB = hB->GetBinError(b);
                double C = hC->GetBinContent(b), sC = hC->GetBinError(b);
                double d = A - B;
                if (d == 0.) continue;
                hratio->SetBinError(b, std::sqrt((B-C)*(B-C)*sA*sA +
                                                  (A-C)*(A-C)*sB*sB +
                                                       d*d*sC*sC) / (d*d));
            }
        };
        fixBmCoverAmB(mott_bmc_over_amb,   mott  [0], mott  [1], mott  [2]);
        fixBmCoverAmB(moller_bmc_over_amb, moller[0], moller[1], moller[2]);
    }
    // Scale to percent
    mott_bmc_over_amb  ->Scale(100.);
    moller_bmc_over_amb->Scale(100.);
    mott_cmd_over_amb  ->Scale(100.);
    moller_cmd_over_amb->Scale(100.);

    const int bmc_color[2] = { kBlue+1, kOrange+7 };
    const char *bmc_label[2] = { "(B-C)/(A-B)", "(C-D)/(A-B)" };
    TH1F *mott_bmc_set  [2] = { mott_bmc_over_amb,   mott_cmd_over_amb   };
    TH1F *moller_bmc_set[2] = { moller_bmc_over_amb, moller_cmd_over_amb };

    // Mott pad — overlay bmc series
    cv->cd(1);
    for (int i = 0; i < 2; i++) {
        mott_bmc_set[i]->SetLineColor(bmc_color[i]);
        mott_bmc_set[i]->SetMarkerColor(bmc_color[i]);
        mott_bmc_set[i]->SetMarkerStyle(24 + i);
        mott_bmc_set[i]->SetLineWidth(2);
        mott_bmc_set[i]->SetStats(0);
        mott_bmc_set[i]->Draw("E SAME");
    }
    TLine *line_mott2 = new TLine(mott_ratio[0]->GetXaxis()->GetXmin(), 0.,
                                   mott_ratio[0]->GetXaxis()->GetXmax(), 0.);
    line_mott2->SetLineStyle(3); line_mott2->SetLineColor(kGray+1); line_mott2->Draw();
    for (int i = 0; i < 2; i++) leg1->AddEntry(mott_bmc_set[i], bmc_label[i], "PE");
    leg1->Draw();

    // Moller pad — overlay bmc series
    cv->cd(2);
    for (int i = 0; i < 2; i++) {
        moller_bmc_set[i]->SetLineColor(bmc_color[i]);
        moller_bmc_set[i]->SetMarkerColor(bmc_color[i]);
        moller_bmc_set[i]->SetMarkerStyle(24 + i);
        moller_bmc_set[i]->SetLineWidth(2);
        moller_bmc_set[i]->Draw("E SAME");
    }
    TLine *line_mol2 = new TLine(moller_ratio[0]->GetXaxis()->GetXmin(), 0.,
                                  moller_ratio[0]->GetXaxis()->GetXmax(), 0.);
    line_mol2->SetLineStyle(3); line_mol2->SetLineColor(kGray+1); line_mol2->Draw();
    for (int i = 0; i < 2; i++) leg2->AddEntry(moller_bmc_set[i], bmc_label[i], "PE");
    leg2->Draw();

    cv->SaveAs("compare_types.png");

    // ── Save to ROOT file ─────────────────────────────────────────────────
    TFile *fout = TFile::Open("compare_types_output.root", "RECREATE");
    mott_denom  ->Write();
    moller_denom->Write();
    for (int i = 0; i < 3; i++) {
        mott_ratio  [i]->Write();
        moller_ratio[i]->Write();
    }
    mott_bmc_over_amb  ->Write();
    moller_bmc_over_amb->Write();
    mott_cmd_over_amb  ->Write();
    moller_cmd_over_amb->Write();
    cv->Write();
    fout->Close();
    std::cout << "Saved to compare_types_output.root" << std::endl;

    
}
