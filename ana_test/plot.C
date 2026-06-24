

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
    TFile *fin = TFile::Open("../data/q2_plot_output_0p7.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open ../data/q2_plot_output_0p7.root" << std::endl;
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
    pad_top->SetLogx();
    pad_top->Draw();
    pad_top->cd();

    cross_section->SetTitle("e-p Cross Section @ 0.7 GeV");
    cross_section->GetXaxis()->SetLabelSize(0.);
    cross_section->GetXaxis()->SetTitleSize(0.);
    cross_section->GetYaxis()->SetTitle("d#sigma/d#Omega (arbitrary units)");
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

    TLatex *prelim = new TLatex(0.50, 0.48, "Very Preliminary");
    prelim->SetNDC();
    prelim->SetTextAlign(22);
    prelim->SetTextFont(62);
    prelim->SetTextSize(0.115);
    prelim->SetTextAngle(24);
    prelim->SetTextColorAlpha(kGray - 1, 0.25);
    //prelim->Draw();

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
    rel_err->Scale(sqrt(0.164753/3.18));
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

    TLatex *stat_note1 = new TLatex(0.16, 0.84, "status by June 21st 16:00");
    stat_note1->SetNDC();
    stat_note1->SetTextSize(0.070);
    stat_note1->Draw();

    cv->cd();
    cv->SaveAs("cross_section.png");
}

void plot_cross_section_three_energy() {
    const int nE = 3;
    const char *files[nE] = {
        "../data/q2_plot_output_0p7.root",
        "../data/q2_plot_output_2p2.root",
        "../data/q2_plot_output_3p5.root"
    };
    const char *labels[nE] = {"0.7 GeV", "2.2 GeV", "3.5 GeV"};
    const char *output[nE] = {"cross_section_0p7.png", "cross_section_2p2.png", "cross_section_3p5.png"};
    const double current_charge[nE] = {3.18, 78., 210.}; // mC
    const double run_charge[nE]     = {0.164753, 0.248315, 4.639600}; // mC, charge used to make rel_err

    gStyle->SetOptStat(0);

    for (int i = 0; i < nE; i++) {
        TFile *fin = TFile::Open(files[i]);
        if (!fin || fin->IsZombie()) {
            std::cerr << "Cannot open " << files[i] << std::endl;
            return;
        }

        TH1F *h_cs = (TH1F*)fin->Get("cross_section");
        TH1F *h_gen = (TH1F*)fin->Get("gen_cross");
        TH1F *h_er = (TH1F*)fin->Get("rel_err");
        if (!h_cs || !h_gen || !h_er) {
            std::cerr << "Missing cross_section, gen_cross or rel_err in " << files[i] << std::endl;
            fin->Close();
            return;
        }

        TH1F *cross_section = (TH1F*)h_cs->Clone(Form("cross_section_%d_individual", i));
        TH1F *gen_cross     = (TH1F*)h_gen->Clone(Form("gen_cross_%d_individual", i));
        TH1F *rel_err       = (TH1F*)h_er->Clone(Form("rel_err_%d_individual", i));
        cross_section->SetDirectory(nullptr);
        gen_cross->SetDirectory(nullptr);
        rel_err->SetDirectory(nullptr);
        fin->Close();

        TCanvas *cv = new TCanvas(Form("c_cross_section_%d_individual", i),
                                  Form("Cross Section %s", labels[i]), 1000, 800);

        TPad *pad_top = new TPad(Form("pad_top_%d_individual", i), "", 0., 0.30, 1., 1.);
        pad_top->SetBottomMargin(0.02);
        pad_top->SetLeftMargin(0.13);
        pad_top->SetRightMargin(0.05);
        pad_top->SetGrid();
        pad_top->SetLogx();
        pad_top->Draw();
        pad_top->cd();

        cross_section->SetTitle(Form("e-p Cross Section @ %s", labels[i]));
        cross_section->GetXaxis()->SetLabelSize(0.);
        cross_section->GetXaxis()->SetTitleSize(0.);
        cross_section->GetYaxis()->SetTitle("d#sigma/d#Omega (arbitrary units)");
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

        TLatex *prelim = new TLatex(0.50, 0.48, "Very Preliminary");
        prelim->SetNDC();
        prelim->SetTextAlign(22);
        prelim->SetTextFont(62);
        prelim->SetTextSize(0.115);
        prelim->SetTextAngle(24);
        //prelim->SetTextColorAlpha(kGray - 1, 0.25);
        //prelim->Draw();

        cv->cd();
        TPad *pad_bot = new TPad(Form("pad_bot_%d_individual", i), "", 0., 0., 1., 0.30);
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
        rel_err->Scale(std::sqrt(run_charge[i] / current_charge[i]));
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

        pad_bot->Update();
        Double_t uymin = pad_bot->GetUymin();
        Double_t uymax = pad_bot->GetUymax();
        Double_t yNDC  = pad_bot->GetBottomMargin() +
                         (log10(0.4) - uymin) / (uymax - uymin) *
                         (1. - pad_bot->GetTopMargin() - pad_bot->GetBottomMargin());
        TLatex *lab04 = new TLatex(pad_bot->GetLeftMargin() - 0.01, yNDC, "0.4");
        lab04->SetNDC();
        lab04->SetTextAlign(32);
        lab04->SetTextSize(0.09);
        lab04->Draw();

        TLatex *stat_note1 = new TLatex(0.16, 0.84, "status by June 21st 16:00");
        stat_note1->SetNDC();
        stat_note1->SetTextSize(0.070);
        stat_note1->Draw();

        cv->cd();
        cv->SaveAs(output[i]);
    }
}

void plot_cross_section_all(){
    const int nE = 3;
    const char *files[nE] = {
        "../data/q2_plot_output_0p7.root",
        "../data/q2_plot_output_2p2.root",
        "../data/q2_plot_output_3p5.root"
    };
    const char *labels[nE] = {"0.7 GeV", "2.2 GeV", "3.5 GeV"};
    const int colors[nE] = {kBlue + 1, kGreen + 2, kRed + 1};
    const int markers[nE] = {20, 21, 22};
    const double current_charge[nE] = {3.18, 78., 210.}; // mC
    const double run_charge[nE]     = { 0.164753,  0.248315,  4.639600}; // mC, charge used to make rel_err
    double stat_scale[nE] = {1., 1., 1.};
    for (int i = 0; i < nE; i++) {
        if (current_charge[i] <= 0.) {
            std::cerr << "Invalid current charge for " << labels[i] << std::endl;
            return;
        }
        stat_scale[i] = sqrt(run_charge[i] / current_charge[i]);
        std::cout << labels[i]
                  << ": current charge = " << current_charge[i] << " mC, "
                  << "run charge = " << run_charge[i] << " mC, "
                  << "stat. error scale = " << stat_scale[i] << std::endl;
    }

    TH1F *cross_section[nE] = {nullptr};
    TH1F *gen_cross[nE] = {nullptr};
    TGraph *gen_graph[nE] = {nullptr};
    TH1F *rel_err[nE] = {nullptr};

    double x_min = 1.e30, x_max = -1.e30;
    double cs_min = 1.e30, cs_max = -1.e30;
    double err_min = 1.e30, err_max = -1.e30;

    for (int i = 0; i < nE; i++) {
        TFile *fin = TFile::Open(files[i]);
        if (!fin || fin->IsZombie()) {
            std::cerr << "Cannot open " << files[i] << std::endl;
            return;
        }

        TH1F *h_cs = (TH1F*)fin->Get("cross_section");
        TH1F *h_gen = (TH1F*)fin->Get("gen_cross");
        TH1F *h_er = (TH1F*)fin->Get("rel_err");
        if (!h_cs || !h_gen || !h_er) {
            std::cerr << "Missing cross_section, gen_cross or rel_err in " << files[i] << std::endl;
            fin->Close();
            return;
        }

        cross_section[i] = (TH1F*)h_cs->Clone(Form("cross_section_%d_plot", i));
        gen_cross[i]     = (TH1F*)h_gen->Clone(Form("gen_cross_%d_plot", i));
        rel_err[i]       = (TH1F*)h_er->Clone(Form("rel_err_%d_plot", i));
        cross_section[i]->SetDirectory(nullptr);
        gen_cross[i]->SetDirectory(nullptr);
        rel_err[i]->SetDirectory(nullptr);
        fin->Close();

        rel_err[i]->Scale(stat_scale[i]);

        for (int ib = 1; ib <= cross_section[i]->GetNbinsX(); ib++) {
            double x_low = cross_section[i]->GetXaxis()->GetBinLowEdge(ib);
            double x_up  = cross_section[i]->GetXaxis()->GetBinUpEdge(ib);
            double y     = cross_section[i]->GetBinContent(ib);
            cross_section[i]->SetBinError(ib, (y > 0.) ? 0.01*y : 0.);
            double dy    = cross_section[i]->GetBinError(ib);
            if (x_low > 0. && x_low < x_min) x_min = x_low;
            if (x_up  > 0. && x_up  > x_max) x_max = x_up;
            if (y > 0.) {
                if (y < cs_min) cs_min = y;
                if (y + dy > cs_max) cs_max = y + dy;
            }
        }

        for (int ib = 1; ib <= gen_cross[i]->GetNbinsX(); ib++) {
            double y  = gen_cross[i]->GetBinContent(ib);
            gen_cross[i]->SetBinError(ib, (y > 0.) ? 0.0*y : 0.);
            double dy = gen_cross[i]->GetBinError(ib);
            if (y > 0.) {
                if (y < cs_min) cs_min = y;
                if (y + dy > cs_max) cs_max = y + dy;
            }
        }

        gen_graph[i] = new TGraph();
        gen_graph[i]->SetName(Form("gen_graph_%d_plot", i));
        int ip = 0;
        for (int ib = 1; ib <= gen_cross[i]->GetNbinsX(); ib++) {
            double x = gen_cross[i]->GetXaxis()->GetBinCenter(ib);
            double y = gen_cross[i]->GetBinContent(ib);
            if (x > 0. && y > 0.) {
                gen_graph[i]->SetPoint(ip, x, y);
                ip++;
            }
        }

        for (int ib = 1; ib <= rel_err[i]->GetNbinsX(); ib++) {
            double y = rel_err[i]->GetBinContent(ib);
            if (y > 0.) {
                if (y < err_min) err_min = y;
                if (y > err_max) err_max = y;
            }
        }
    }

    if (x_min >= x_max || cs_min >= cs_max || err_min >= err_max) {
        std::cerr << "Invalid ranges while building all-energy cross-section plot." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    TCanvas *cv = new TCanvas("c_cross_section_all", "Cross Section All Energies", 1000, 800);

    // ── Upper pad: cross sections for all beam energies ───────────────────
    TPad *pad_top = new TPad("pad_top_all", "", 0., 0.30, 1., 1.);
    pad_top->SetBottomMargin(0.02);
    pad_top->SetLeftMargin(0.13);
    pad_top->SetRightMargin(0.05);
    pad_top->SetGrid();
    pad_top->SetLogx();
    pad_top->Draw();
    pad_top->cd();

    TH1F *frame_top = new TH1F("frame_cross_section_all",
                               "e-p Cross Section;Q^{2} (GeV^{2});d#sigma/d#Omega (arbitrary units)",
                               100, x_min*0.9, x_max*1.1);
    frame_top->SetMinimum(cs_min*0.5);
    frame_top->SetMaximum(cs_max*2.0);
    frame_top->GetXaxis()->SetLabelSize(0.);
    frame_top->GetXaxis()->SetTitleSize(0.);
    frame_top->GetYaxis()->SetTitleSize(0.060);
    frame_top->GetYaxis()->SetTitleOffset(0.95);
    frame_top->GetYaxis()->CenterTitle();
    frame_top->Draw();

    TLatex *prelim = new TLatex(0.50, 0.48, "Very Preliminary");
    prelim->SetNDC();
    prelim->SetTextAlign(22);
    prelim->SetTextFont(62);
    prelim->SetTextSize(0.115);
    prelim->SetTextAngle(24);
    prelim->SetTextColorAlpha(kGray + 1, 0.50);
    //prelim->Draw();

    TLegend *leg_top = new TLegend(0.40, 0.68, 0.90, 0.88);
    leg_top->SetBorderSize(0);
    leg_top->SetTextSize(0.038);
    leg_top->SetNColumns(2);

    for (int i = 0; i < nE; i++) {
        cross_section[i]->SetStats(0);
        cross_section[i]->SetMarkerStyle(markers[i]);
        cross_section[i]->SetMarkerSize(1.3);
        cross_section[i]->SetMarkerColor(colors[i]);
        cross_section[i]->SetLineColor(colors[i]);
        cross_section[i]->SetLineWidth(2);
        cross_section[i]->Draw("E1X0P SAME");
        leg_top->AddEntry(cross_section[i], Form("%s Data", labels[i]), "lp");

        gen_graph[i]->SetLineColor(colors[i]);
        gen_graph[i]->SetLineStyle(2);
        gen_graph[i]->SetLineWidth(3);
        gen_graph[i]->Draw("L SAME");
        leg_top->AddEntry(gen_graph[i], "Generator", "l");
    }

    for (int i = 0; i < nE; i++) {
        int n = gen_graph[i]->GetN();
        if (n <= 0) continue;
        double x = 0., y = 0.;
        gen_graph[i]->GetPoint(n - 1, x, y);
        TLatex *e_label = new TLatex(x*1.03, y, labels[i]);
        e_label->SetTextColor(colors[i]);
        e_label->SetTextFont(62);
        e_label->SetTextSize(0.038);
        e_label->SetTextAlign(12);
        e_label->Draw();
    }
    leg_top->Draw();

    // ── Lower pad: relative statistical errors for all beam energies ──────
    cv->cd();
    TPad *pad_bot = new TPad("pad_bot_all", "", 0., 0., 1., 0.30);
    pad_bot->SetTopMargin(0.02);
    pad_bot->SetBottomMargin(0.32);
    pad_bot->SetLeftMargin(0.13);
    pad_bot->SetRightMargin(0.05);
    pad_bot->SetGrid();
    pad_bot->SetLogx();
    pad_bot->SetLogy();
    pad_bot->Draw();
    pad_bot->cd();

    TH1F *frame_bot = new TH1F("frame_rel_err_all",
                               ";Q^{2} (GeV^{2});Stat. Error (%)",
                               100, x_min*0.9, x_max*1.1);
    frame_bot->SetMinimum(err_min*0.5);
    frame_bot->SetMaximum(err_max*2.0);
    frame_bot->GetXaxis()->SetTitleSize(0.12);
    frame_bot->GetYaxis()->SetTitleSize(0.11);
    frame_bot->GetXaxis()->SetLabelSize(0.10);
    frame_bot->GetYaxis()->SetLabelSize(0.09);
    frame_bot->GetXaxis()->SetTitleOffset(1.0);
    frame_bot->GetYaxis()->SetTitleOffset(0.55);
    frame_bot->GetXaxis()->CenterTitle();
    frame_bot->GetYaxis()->CenterTitle();
    frame_bot->GetXaxis()->SetMoreLogLabels();
    frame_bot->GetXaxis()->SetNoExponent();
    frame_bot->GetYaxis()->SetMoreLogLabels();
    frame_bot->GetYaxis()->SetNoExponent();
    frame_bot->Draw();

    for (int i = 0; i < nE; i++) {
        rel_err[i]->SetStats(0);
        rel_err[i]->SetMarkerStyle(markers[i]);
        rel_err[i]->SetMarkerSize(0.9);
        rel_err[i]->SetMarkerColor(colors[i]);
        rel_err[i]->SetLineColor(colors[i]);
        rel_err[i]->SetLineWidth(2);
        rel_err[i]->Draw("E1P SAME");
    }

    TLatex *stat_note1 = new TLatex(0.16, 0.84, "status by June 21st 16:00");
    stat_note1->SetNDC();
    stat_note1->SetTextSize(0.070);
    stat_note1->Draw();

    TLatex *stat_note2 = new TLatex(0.16, 0.74, "(scaled by beam charge)");
    stat_note2->SetNDC();
    stat_note2->SetTextSize(0.070);
    //stat_note2->Draw();

    cv->cd();
    cv->SaveAs("cross_section_all.png");
}

void plot_GE_all(){
    const int nE = 3;
    const char *files[nE] = {
        "../data/q2_plot_output_0p7.root",
        "../data/q2_plot_output_2p2.root",
        "../data/q2_plot_output_3p5.root"
    };
    const char *labels[nE] = {"0.7 GeV", "2.2 GeV", "3.5 GeV"};
    const int colors[nE] = {kBlue + 1, kGreen + 2, kRed + 1};
    const int markers[nE] = {20, 21, 22};
    const double GE_percent_err[nE] = {3., 3., 0.8}; // percent

    TH1F *GE[nE] = {nullptr};
    TH1F *GE_gen[nE] = {nullptr};
    TGraph *GE_gen_graph[nE] = {nullptr};

    double x_min = 1.e30, x_max = -1.e30;
    double ge_min = 1.e30, ge_max = -1.e30;

    for (int i = 0; i < nE; i++) {
        TFile *fin = TFile::Open(files[i]);
        if (!fin || fin->IsZombie()) {
            std::cerr << "Cannot open " << files[i] << std::endl;
            return;
        }

        TH1F *h_ge = (TH1F*)fin->Get("GE");
        TH1F *h_ge_gen = (TH1F*)fin->Get("GE_gen");
        if (!h_ge || !h_ge_gen) {
            std::cerr << "Missing GE or GE_gen in " << files[i] << std::endl;
            fin->Close();
            return;
        }

        GE[i] = (TH1F*)h_ge->Clone(Form("GE_%d_plot", i));
        GE_gen[i] = (TH1F*)h_ge_gen->Clone(Form("GE_gen_%d_plot", i));
        GE[i]->SetDirectory(nullptr);
        GE_gen[i]->SetDirectory(nullptr);
        fin->Close();

        for (int ib = 1; ib <= GE[i]->GetNbinsX(); ib++) {
            double x_low = GE[i]->GetXaxis()->GetBinLowEdge(ib);
            double x_up  = GE[i]->GetXaxis()->GetBinUpEdge(ib);
            double y     = GE[i]->GetBinContent(ib);
            GE[i]->SetBinError(ib, (y > 0.) ? GE_percent_err[i]/100.*y : 0.);
            if (x_low > 0. && x_low < x_min) x_min = x_low;
            if (x_up  > 0. && x_up  > x_max) x_max = x_up;
            if (y > 0.) {
                if (y < ge_min) ge_min = y;
                if (y + GE[i]->GetBinError(ib) > ge_max) ge_max = y + GE[i]->GetBinError(ib);
            }
        }

        GE_gen_graph[i] = new TGraph();
        GE_gen_graph[i]->SetName(Form("GE_gen_graph_%d_plot", i));
        int ip = 0;
        for (int ib = 1; ib <= GE_gen[i]->GetNbinsX(); ib++) {
            double x = GE_gen[i]->GetXaxis()->GetBinCenter(ib);
            double y = GE_gen[i]->GetBinContent(ib);
            if (x > 0. && y > 0.) {
                GE_gen_graph[i]->SetPoint(ip, x, y);
                ip++;
                if (y < ge_min) ge_min = y;
                if (y > ge_max) ge_max = y;
            }
        }
    }

    if (x_min >= x_max || ge_min >= ge_max) {
        std::cerr << "Invalid ranges while building all-energy GE plot." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    TCanvas *cv = new TCanvas("c_GE_all", "Electric Form Factor All Energies", 1000, 800);

    TPad *pad_top = new TPad("pad_GE_top_all", "", 0., 0., 1., 1.);
    pad_top->SetBottomMargin(0.13);
    pad_top->SetLeftMargin(0.13);
    pad_top->SetRightMargin(0.05);
    pad_top->SetGrid();
    pad_top->SetLogx();
    pad_top->Draw();
    pad_top->cd();

    TH1F *frame_top = new TH1F("frame_GE_top_all",
                               "Electric Form Factor G_{E};Q^{2} (GeV^{2});G_{E}",
                               100, x_min*0.9, x_max*1.1);
    frame_top->SetMinimum(0.5);
    frame_top->SetMaximum(1.2);
    frame_top->GetXaxis()->SetTitleSize(0.050);
    frame_top->GetYaxis()->SetTitleSize(0.050);
    frame_top->GetXaxis()->SetLabelSize(0.040);
    frame_top->GetYaxis()->SetLabelSize(0.040);
    frame_top->GetXaxis()->SetTitleOffset(1.0);
    frame_top->GetYaxis()->SetTitleOffset(1.1);
    frame_top->GetXaxis()->CenterTitle();
    frame_top->GetYaxis()->CenterTitle();
    frame_top->GetXaxis()->SetMoreLogLabels();
    frame_top->GetXaxis()->SetNoExponent();
    frame_top->Draw();

    TLegend *leg_top = new TLegend(0.40, 0.68, 0.90, 0.88);
    leg_top->SetBorderSize(0);
    leg_top->SetTextSize(0.038);
    leg_top->SetNColumns(2);

    for (int i = 0; i < nE; i++) {
        GE[i]->SetStats(0);
        GE[i]->SetMarkerStyle(markers[i]);
        GE[i]->SetMarkerSize(1.3);
        GE[i]->SetMarkerColor(colors[i]);
        GE[i]->SetLineColor(colors[i]);
        GE[i]->SetLineWidth(2);
        GE[i]->Draw("E1X0P SAME");
        leg_top->AddEntry(GE[i], Form("%s Data", labels[i]), "lp");

        GE_gen_graph[i]->SetLineColor(colors[i]);
        GE_gen_graph[i]->SetLineStyle(2);
        GE_gen_graph[i]->SetLineWidth(3);
        GE_gen_graph[i]->Draw("L SAME");
        leg_top->AddEntry(GE_gen_graph[i], "Generator", "l");
    }
    leg_top->Draw();

    cv->cd();
    cv->SaveAs("GE_all.png");
}
