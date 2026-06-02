

void plot() {
    const int bin_idx = 9;  // 0-based bin index to plot
    const float binEdge[] = {
        0.500f, 0.550f, 0.600f, 0.650f, 0.700f, 0.750f, 0.775f, 0.800f, 0.825f, 0.850f,
        0.875f, 0.900f, 0.940f, 0.975f, 1.014f, 1.057f, 1.105f, 1.157f, 1.211f, 1.270f,
        1.338f, 1.417f, 1.514f, 1.634f, 1.787f, 2.000f, 2.213f, 2.492f, 2.792f, 3.092f,
        3.392f, 3.692f, 3.992f, 4.292f
    };
    const float Ebeam = 3484.f;

    TFile *fin = TFile::Open("../data/q2_plot_output.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open ../data/q2_plot_output.root" << std::endl;
        return;
    }

    // E_recon histograms were saved inside canvases; extract from canvas
    TCanvas *ctmp = (TCanvas*)fin->Get(Form("c_erecon_bin%d", bin_idx));
    if (!ctmp) {
        std::cerr << "Canvas c_erecon_bin" << bin_idx << " not found" << std::endl;
        fin->Close(); return;
    }
    TH1F *h = nullptr;
    TIter next(ctmp->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(TH1::Class())) { h = (TH1F*)obj; break; }
    }
    if (!h) {
        std::cerr << "No histogram found in canvas" << std::endl;
        fin->Close(); return;
    }
    h = (TH1F*)h->Clone(Form("E_recon_bin%d_plot", bin_idx));
    h->SetDirectory(nullptr);
    fin->Close();

    // ── Normalize Y axis to counts per MeV ───────────────────────────────
    double bw = h->GetBinWidth(1);  // MeV per bin
    h->Scale(1.0 / bw);

    // ── Gaussian re-fit ───────────────────────────────────────────────────
    double resolution = 0.033;
    double sigma_hint = Ebeam * resolution / sqrt(Ebeam / 1000.);
    TF1 *fit = new TF1("gaus_fit", "gaus", Ebeam - 2*sigma_hint, Ebeam + 2*sigma_hint);
    fit->SetParameters(h->GetMaximum(), Ebeam, sigma_hint);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    h->Fit(fit, "RQ");

    double mean      = fit->GetParameter(1);
    double sigma     = fit->GetParameter(2);
    double mean_err  = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    double res_pct   = sigma / mean * sqrt(mean / 1000.) * 100.;
    double res_err   = res_pct * sqrt(pow(sigma_err/sigma, 2) + pow(mean_err/mean, 2) * 0.25);

    // ── Canvas & histogram style ──────────────────────────────────────────
    TCanvas *c = new TCanvas("c_erecon", "E_recon", 800, 600);
    c->SetGrid();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);

    h->SetStats(0);
    h->SetLineColor(kAzure + 1);
    h->SetLineWidth(2);
    h->SetFillColor(kAzure - 9);
    h->SetFillStyle(1001);
    h->SetTitle(Form("Reconstructed Scattering Angle #in [%.3f, %.3f] deg;Energy (MeV);Counts / %.1f MeV",
                     binEdge[bin_idx], binEdge[bin_idx+1], bw));
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->Draw("HIST");
    fit->Draw("SAME");

    // ── Info box ─────────────────────────────────────────────────────────
    // Three columns via a TPaveText with fixed-width text alignment
    TPaveText *pt = new TPaveText(0.13, 0.67, 0.88, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.040);

    // Row 1: headers
    TText *t1 = pt->AddText(Form("Beam energy: %.0f MeV", Ebeam));
    TText *t2 = pt->AddText(Form("Recon. E = %.2f MeV", mean));
    TText *t3 = pt->AddText(Form("Resolution = %.2f %%", res_pct));
    pt->Draw();

    c->Update();
    c->SaveAs(Form("E_recon_bin%d.png", bin_idx));
}

void plot_cross_section() {
    TFile *fin = TFile::Open("../data/q2_plot_output.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open ../data/q2_plot_output.root" << std::endl;
        return;
    }

    TH1F *cross_section = (TH1F*)((TH1F*)fin->Get("cross_section"))->Clone("cross_section_plot");
    TH1F *gen_cross     = (TH1F*)((TH1F*)fin->Get("gen_cross"))    ->Clone("gen_cross_plot");
    TH1F *rel_err       = (TH1F*)((TH1F*)fin->Get("rel_err"))      ->Clone("rel_err_plot");
    cross_section->SetDirectory(nullptr);
    gen_cross    ->SetDirectory(nullptr);
    rel_err      ->SetDirectory(nullptr);
    fin->Close();

    gStyle->SetOptStat(0);

    TCanvas *cv = new TCanvas("c_cross_section", "Cross Section", 1000, 800);

    // ── Upper pad: cross section + gen reference ──────────────────────────
    TPad *pad_top = new TPad("pad_top", "", 0., 0.30, 1., 1.);
    pad_top->SetBottomMargin(0.02);
    pad_top->SetLeftMargin(0.13);
    pad_top->SetRightMargin(0.05);
    pad_top->SetGrid();
    pad_top->SetLogy();
    pad_top->SetLogx();
    pad_top->Draw();
    pad_top->cd();

    cross_section->SetTitle("e-p Cross Section @ 3.5 GeV");
    cross_section->GetXaxis()->SetLabelSize(0.);
    cross_section->GetXaxis()->SetTitleSize(0.);
    cross_section->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
    cross_section->GetYaxis()->SetTitleSize(0.060);
    cross_section->GetYaxis()->SetTitleOffset(0.95);
    cross_section->GetYaxis()->CenterTitle();
    cross_section->SetMarkerStyle(20);
    cross_section->SetMarkerSize(1.0);
    cross_section->SetMarkerColor(kBlack);
    cross_section->SetLineColor(kBlack);
    cross_section->Draw("E1P");

    gen_cross->SetLineColor(kRed + 1);
    gen_cross->SetMarkerStyle(24);
    gen_cross->SetMarkerColor(kRed + 1);
    gen_cross->SetLineColor(kRed + 1);
    gen_cross->Draw("E1P SAME");

    TLegend *leg = new TLegend(0.55, 0.65, 0.90, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.050);
    leg->AddEntry(cross_section, "Data",      "lp");
    leg->AddEntry(gen_cross,     "Generator", "lp");
    leg->Draw();

    // ── Lower pad: relative statistical error ─────────────────────────────
    cv->cd();
    TPad *pad_bot = new TPad("pad_bot", "", 0., 0., 1., 0.30);
    pad_bot->SetTopMargin(0.02);
    pad_bot->SetBottomMargin(0.32);
    pad_bot->SetLeftMargin(0.13);
    pad_bot->SetRightMargin(0.05);
    pad_bot->SetGrid();
    pad_bot->SetLogx();
    pad_bot->SetLogy();
    pad_bot->Draw();
    pad_bot->cd();

    rel_err->SetTitle("");
    rel_err->Scale(1./6.);
    rel_err->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    rel_err->GetYaxis()->SetTitle("Stat. Error (%)");
    rel_err->GetXaxis()->SetTitleSize(0.12);
    rel_err->GetYaxis()->SetTitleSize(0.11);
    rel_err->GetXaxis()->SetLabelSize(0.10);
    rel_err->GetYaxis()->SetLabelSize(0.09);
    rel_err->GetXaxis()->SetTitleOffset(1.0);
    rel_err->GetYaxis()->SetTitleOffset(0.55);
    rel_err->GetXaxis()->CenterTitle();
    rel_err->GetYaxis()->CenterTitle();
    rel_err->GetXaxis()->SetMoreLogLabels();
    rel_err->GetXaxis()->SetNoExponent();
    rel_err->GetYaxis()->SetMoreLogLabels();
    rel_err->GetYaxis()->SetNoExponent();
    rel_err->SetFillColor(kAzure - 9);
    rel_err->SetLineColor(kBlue + 1);
    rel_err->Draw("BAR");

    // Manually add the 0.4 label that ROOT skips on log Y axis
    pad_bot->Update();
    Double_t uymin = pad_bot->GetUymin(); // log10(ymin) in log scale
    Double_t uymax = pad_bot->GetUymax(); // log10(ymax) in log scale
    Double_t yNDC  = pad_bot->GetBottomMargin() +
                     (log10(0.4) - uymin) / (uymax - uymin) *
                     (1. - pad_bot->GetTopMargin() - pad_bot->GetBottomMargin());
    TLatex *lab04 = new TLatex(pad_bot->GetLeftMargin() - 0.01, yNDC, "0.4");
    lab04->SetNDC();
    lab04->SetTextAlign(32);
    lab04->SetTextSize(0.09);
    lab04->Draw();

    cv->cd();
    cv->SaveAs("cross_section.png");
}