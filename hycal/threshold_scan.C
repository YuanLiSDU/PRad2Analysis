
void threshold_scan()
{
    const char *base = "/home/liyuan/PRad2Analysis/data/total_sum_thres/";
    struct RunInfo { const char *file; const char *label; int color; int thres; };
    RunInfo runs[] = {
        { "24352_quick.root", "Run 24352", kBlue   , 800},
        { "24357_quick.root", "Run 24357", kRed    , 1200},
        { "24358_quick.root", "Run 24358", kGreen+2, 1500},
        { "24359_quick.root", "Run 24359", kOrange+1, 1800},
    };
    const int nRuns = sizeof(runs) / sizeof(runs[0]);

    // Find peak value near 3500 MeV (search in 3000-4000 MeV range)
    auto getPeakNear3500 = [](TH1F *h) -> double {
        int binLo = h->FindBin(3000);
        int binHi = h->FindBin(4000);
        double peak = 0;
        for (int i = binLo; i <= binHi; ++i) {
            double v = h->GetBinContent(i);
            if (v > peak) peak = v;
        }
        return peak;
    };

    TH1F *hists[nRuns] = {};
    for (int i = 0; i < nRuns; ++i) {
        TString path = TString(base) + runs[i].file;
        TFile *f = TFile::Open(path);
        if (!f || f->IsZombie()) { printf("Cannot open %s\n", path.Data()); return; }
        TH1F *h = (TH1F*)f->Get("energy_plots/total_energy");
        if (!h) { printf("Cannot find total_energy in %s\n", path.Data()); return; }
        TString cloneName = TString("h_") + runs[i].label;
        cloneName.ReplaceAll(" ", "_");
        hists[i] = (TH1F*)h->Clone(cloneName);
        hists[i]->SetDirectory(0);
        f->Close();

        double peak = getPeakNear3500(hists[i]);
        printf("Peak near 3500 MeV: %s = %.1f\n", runs[i].label, peak);
        if (peak <= 0) { printf("Invalid peak value for %s, aborting.\n", runs[i].label); return; }
        hists[i]->Scale(1.0 / peak);

        hists[i]->SetLineColor(runs[i].color);
        hists[i]->SetLineWidth(2);
    }

    hists[0]->SetTitle("Total Energy per Event (normalized to 3500 MeV peak)");
    hists[0]->GetXaxis()->SetTitle("Total Energy (MeV)");
    hists[0]->GetYaxis()->SetTitle("Counts / Peak");

    TCanvas *c = new TCanvas("c_threshold_scan", "Threshold Scan", 800, 600);
    c->SetLogy();

    double ymax = 0;
    for (int i = 0; i < nRuns; ++i)
        ymax = std::max(ymax, hists[i]->GetMaximum());
    hists[0]->SetMaximum(ymax * 1.3);
    hists[0]->SetMinimum(1e-4);

    hists[0]->Draw("HIST");
    for (int i = 1; i < nRuns; ++i)
        hists[i]->Draw("HIST SAME");

    // For first 3 runs: find min bin in 500-2500 MeV as left cut,
    // then find the cumulative 1% point from that cut, and draw a vertical line.
    // The histogram display is unchanged; cut only affects this calculation.
    auto find1PctLine = [&](TH1F *h, int color) -> TLine* {
        int binLo = h->FindBin(500);
        int binHi = h->FindBin(2500);
        // Find minimum bin in [500, 2500] MeV
        double minVal = h->GetBinContent(binLo);
        int cutBin = binLo;
        for (int j = binLo + 1; j <= binHi; ++j) {
            double v = h->GetBinContent(j);
            if (v < minVal) { minVal = v; cutBin = j; }
        }
        // Integral from cutBin to right end
        double total = h->Integral(cutBin, h->GetNbinsX()*3./4.);
        // Find leftmost bin where cumulative sum >= 1% of total
        double cumul = 0;
        double xLine = h->GetBinCenter(h->GetNbinsX());
        for (int j = cutBin; j <= h->GetNbinsX(); ++j) {
            cumul += h->GetBinContent(j);
            if (cumul >= 0.05 * total) { xLine = h->GetBinCenter(j); break; }
        }
        printf("  1%% point after cut at bin %d (%.1f MeV): x = %.1f MeV\n",
               cutBin, h->GetBinCenter(cutBin), xLine);
        TLine *line = new TLine(xLine, h->GetMinimum(), xLine, h->GetMaximum());
        line->SetLineColor(color);
        line->SetLineWidth(2);
        line->SetLineStyle(2); // dashed
        return line;
    };

    for (int i = 0; i < 3; ++i) {
        TLine *l = find1PctLine(hists[i], runs[i].color);
        l->Draw();
    }

    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    for (int i = 0; i < nRuns; ++i)
        leg->AddEntry(hists[i], Form("%s (thres=%d MeV)", runs[i].label, runs[i].thres), "l");
    leg->Draw();

    c->SaveAs(TString(base) + "threshold_scan.pdf");
    printf("Saved: threshold_scan.pdf\n");
}
