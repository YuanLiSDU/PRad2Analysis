
double test_charge = 0.248315/0.8; // mC
double current_charge = 13.58; // mC
double aim_charge = 100.; // mC

void plot_2p2() {
    TFile *fin = TFile::Open("q2_plot_output.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open q2_plot_output.root" << std::endl;
        return;
    }

    TH1F *cross_section = (TH1F*)((TH1F*)fin->Get("cross_section"))->Clone("cross_section_plot");
    TH1F *gen_cross     = (TH1F*)((TH1F*)fin->Get("gen_cross"))    ->Clone("gen_cross_plot");
    TH1F *rel_err       = (TH1F*)((TH1F*)fin->Get("rel_err"))      ->Clone("rel_err_plot");
    TH1F *rel_err_aim = new TH1F("rel_err_aim", "Relative Error Aim", rel_err->GetNbinsX(), rel_err->GetXaxis()->GetXmin(), rel_err->GetXaxis()->GetXmax());
    for(int i = 1; i <= rel_err->GetNbinsX(); i++) {
        rel_err_aim->SetBinContent(i, rel_err->GetBinContent(i));
    }
    cross_section->SetDirectory(nullptr);
    gen_cross    ->SetDirectory(nullptr);
    rel_err      ->SetDirectory(nullptr);
    rel_err_aim  ->SetDirectory(nullptr);
    fin->Close();

    // Scale the rel_err_aim histogram to the aim charge
    rel_err->Scale(sqrt(test_charge/current_charge));
    rel_err_aim->Scale(sqrt(test_charge/ aim_charge));

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

    cross_section->SetTitle("e-p Cross Section @ 2.2 GeV (2 GEMs matching)");
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

    rel_err->SetFillColor(kAzure - 9);
    rel_err->SetLineColor(kBlue + 1);
    rel_err->Draw("BAR");

    //rel_err_aim->SetFillColor(kRed - 9);
    rel_err_aim->SetLineColor(kRed + 1);
    rel_err_aim->SetLineWidth(3);
    rel_err_aim->Draw("HIST SAME");

    rel_err->SetTitle("");
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

    TLatex *charge_current = new TLatex(0.16, 0.84, Form("current charge: %.2f mC", current_charge));
    charge_current->SetNDC();
    charge_current->SetTextColor(kBlue + 1);
    charge_current->SetTextSize(0.075);
    charge_current->Draw();

    TLatex *charge_aim = new TLatex(0.16, 0.73, Form("aim charge: %.0f mC", aim_charge));
    charge_aim->SetNDC();
    charge_aim->SetTextColor(kRed + 1);
    charge_aim->SetTextSize(0.075);
    charge_aim->Draw();

    TLatex *charge_update = new TLatex(0.16, 0.62, "updated at June 16th 07:00");
    charge_update->SetNDC();
    charge_update->SetTextSize(0.075);
    charge_update->Draw();

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
    //lab04->Draw();

    cv->cd();
    cv->SaveAs("cross_section.png");
}
