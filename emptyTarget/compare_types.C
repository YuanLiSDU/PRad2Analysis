// compare_types.C
// Read background_output.root from four run types, compute C/(A-B) and D/(A-B)
// ratios for mott_yield and moller_yield, and overlay them on one canvas.
//
// Usage:
//   root -l 'emptyTarget/compare_types.C'
//   root -l 'emptyTarget/compare_types.C("pathA.root","pathB.root","pathC.root","pathD.root")'

void compare_types(
    const char *fileA = "../data/0.7GeV/typeA_24653.root",
    const char *fileB = "../data/0.7GeV/typeB.root",
    const char *fileC = "../data/0.7GeV/typeC.root",
    const char *fileD = "../data/0.7GeV/typeD.root")
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
        mott_ratio  [i]->SetTitle("Mott Yield Ratio X/(A-B);#theta (deg);Ratio");
        moller_ratio[i]->SetTitle("Moller Yield Ratio X/(A-B);#theta (deg);Ratio");
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
    auto [mott_lo, mott_hi] = yRange(mott_ratio, 3);
    for (int i = 0; i < 3; i++) {
        mott_ratio[i]->SetLineColor(cfg_color[i]);
        mott_ratio[i]->SetMarkerColor(cfg_color[i]);
        mott_ratio[i]->SetMarkerStyle(20 + i);
        mott_ratio[i]->SetLineWidth(2);
        if (i == 0) {
            mott_ratio[i]->GetYaxis()->SetRangeUser(0., mott_hi);
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
    leg1->Draw();

    // --- Moller (double-arm range: 1.6 -- 3.0 deg) ---
    cv->cd(2);
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
    leg2->Draw();

    cv->SaveAs("compare_types.png");

    // ── Second canvas: (B-C)/A and (B-C)/(A-B) ───────────────────────────
    // Numerator: B - C
    TH1F *mott_bmc   = (TH1F*)mott  [1]->Clone("mott_BminusC");
    TH1F *moller_bmc = (TH1F*)moller[1]->Clone("moller_BminusC");
    mott_bmc  ->Add(mott  [2], -1.);
    moller_bmc->Add(moller[2], -1.);

    // (B-C)/A
    TH1F *mott_bmc_over_a   = (TH1F*)mott_bmc  ->Clone("mott_BmC_over_A");
    TH1F *moller_bmc_over_a = (TH1F*)moller_bmc->Clone("moller_BmC_over_A");
    mott_bmc_over_a  ->Divide(mott_bmc,   mott  [0], 1., 1., "");
    moller_bmc_over_a->Divide(moller_bmc, moller[0], 1., 1., "");

    // (B-C)/(A-B)
    TH1F *mott_bmc_over_amb   = (TH1F*)mott_bmc  ->Clone("mott_BmC_over_AmB");
    TH1F *moller_bmc_over_amb = (TH1F*)moller_bmc->Clone("moller_BmC_over_AmB");
    mott_bmc_over_amb  ->Divide(mott_bmc,   mott_denom,   1., 1., "");
    moller_bmc_over_amb->Divide(moller_bmc, moller_denom, 1., 1., "");

    const int bmc_color[2] = { kBlue+1, kOrange+7 };
    const char *bmc_label[2] = { "(B-C)/A", "(B-C)/(A-B)" };
    TH1F *mott_bmc_set  [2] = { mott_bmc_over_a,   mott_bmc_over_amb   };
    TH1F *moller_bmc_set[2] = { moller_bmc_over_a, moller_bmc_over_amb };

    TCanvas *cv2 = new TCanvas("c_compare_bmc", "(B-C) Ratios", 1200, 500);
    cv2->Divide(2, 1);

    // Mott
    cv2->cd(1);
    auto [bmc_mott_lo, bmc_mott_hi] = yRange(mott_bmc_set, 2);
    for (int i = 0; i < 2; i++) {
        mott_bmc_set[i]->SetTitle("Mott (B-C) Ratios;#theta (deg);Ratio");
        mott_bmc_set[i]->SetLineColor(bmc_color[i]);
        mott_bmc_set[i]->SetMarkerColor(bmc_color[i]);
        mott_bmc_set[i]->SetMarkerStyle(20 + i);
        mott_bmc_set[i]->SetLineWidth(2);
        if (i == 0) {
            mott_bmc_set[i]->GetYaxis()->SetRangeUser(bmc_mott_lo, bmc_mott_hi);
            mott_bmc_set[i]->Draw("E");
        } else {
            mott_bmc_set[i]->Draw("E SAME");
        }
    }
    TLine *line_mott2 = new TLine(mott_bmc_set[0]->GetXaxis()->GetXmin(), 0.,
                                   mott_bmc_set[0]->GetXaxis()->GetXmax(), 0.);
    line_mott2->SetLineStyle(2); line_mott2->SetLineColor(kGray+1); line_mott2->Draw();
    TLegend *leg3 = new TLegend(0.65, 0.72, 0.92, 0.90);
    for (int i = 0; i < 2; i++) leg3->AddEntry(mott_bmc_set[i], bmc_label[i], "PE");
    leg3->Draw();

    // Moller (double-arm range)
    cv2->cd(2);
    for (int i = 0; i < 2; i++)
        moller_bmc_set[i]->GetXaxis()->SetRangeUser(1.6, 3.0);
    auto [bmc_mol_lo, bmc_mol_hi] = yRange(moller_bmc_set, 2);
    for (int i = 0; i < 2; i++) {
        moller_bmc_set[i]->SetTitle("Moller (B-C) Ratios;#theta (deg);Ratio");
        moller_bmc_set[i]->SetLineColor(bmc_color[i]);
        moller_bmc_set[i]->SetMarkerColor(bmc_color[i]);
        moller_bmc_set[i]->SetMarkerStyle(20 + i);
        moller_bmc_set[i]->SetLineWidth(2);
        if (i == 0) {
            moller_bmc_set[i]->GetYaxis()->SetRangeUser(bmc_mol_lo, bmc_mol_hi);
            moller_bmc_set[i]->Draw("E");
        } else {
            moller_bmc_set[i]->Draw("E SAME");
        }
    }
    TLine *line_mol2 = new TLine(moller_bmc_set[0]->GetXaxis()->GetXmin(), 0.,
                                  moller_bmc_set[0]->GetXaxis()->GetXmax(), 0.);
    line_mol2->SetLineStyle(2); line_mol2->SetLineColor(kGray+1); line_mol2->Draw();
    TLegend *leg4 = new TLegend(0.65, 0.72, 0.92, 0.90);
    for (int i = 0; i < 2; i++) leg4->AddEntry(moller_bmc_set[i], bmc_label[i], "PE");
    leg4->Draw();

    cv2->SaveAs("compare_types_bmc.png");

    // ── Save to ROOT file ─────────────────────────────────────────────────
    TFile *fout = TFile::Open("compare_types_output.root", "RECREATE");
    mott_denom  ->Write();
    moller_denom->Write();
    for (int i = 0; i < 3; i++) {
        mott_ratio  [i]->Write();
        moller_ratio[i]->Write();
    }
    mott_bmc_over_a    ->Write();
    moller_bmc_over_a  ->Write();
    mott_bmc_over_amb  ->Write();
    moller_bmc_over_amb->Write();
    cv ->Write();
    cv2->Write();
    fout->Close();
    std::cout << "Saved to compare_types_output.root" << std::endl;

    
}
