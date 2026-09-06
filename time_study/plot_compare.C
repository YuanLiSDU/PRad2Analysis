

void plot_compare() {
    TFile *old_f = TFile::Open("atest_full.root");
    TFile *new_file = TFile::Open("atest_full_new.root");

    TH1F *old_hist = (TH1F*)old_f->Get("h1_time_diff_seed");
    TH1F *new_hist = (TH1F*)new_file->Get("h1_time_diff_seed");

    float count_old = old_hist->GetEntries();
    float count_new = new_hist->GetEntries();

    // count the events outside the seed pulse window 2ns
    int bin_2ns_low = old_hist->FindBin(-2);
    int bin_2ns_high = old_hist->FindBin(2);
    int bin_4ns_low = old_hist->FindBin(-4);
    int bin_4ns_high = old_hist->FindBin(4);

    float count_old_outside_2 = old_hist->Integral(0, bin_2ns_low) + old_hist->Integral(bin_2ns_high, old_hist->GetNbinsX() + 1);
    float count_new_outside_2 = new_hist->Integral(0, bin_2ns_low) + new_hist->Integral(bin_2ns_high, new_hist->GetNbinsX() + 1);
    float fraction_old_outside_2 = count_old_outside_2 / count_old;
    float fraction_new_outside_2 = count_new_outside_2 / count_new;

    float count_old_outside_4 = old_hist->Integral(0, bin_4ns_low) + old_hist->Integral(bin_4ns_high, old_hist->GetNbinsX() + 1);
    float count_new_outside_4 = new_hist->Integral(0, bin_4ns_low) + new_hist->Integral(bin_4ns_high, new_hist->GetNbinsX() + 1);
    float fraction_old_outside_4 = count_old_outside_4 / count_old;
    float fraction_new_outside_4 = count_new_outside_4 / count_new;

    old_hist->Scale(count_new / count_old);

    TCanvas *c = new TCanvas("c", "Comparison", 800, 600);
    c->SetLogy();
    c->SetGrid();
    old_hist->SetStats(0);
    new_hist->SetStats(0);
    old_hist->SetLineColor(kRed);
    new_hist->SetLineColor(kBlue);
    old_hist->SetLineWidth(2);
    new_hist->SetLineWidth(2);
    old_hist->GetXaxis()->SetRangeUser(-6, 25);
    old_hist->SetTitle("#Delta t of Neighbouring Pulses with Seed Pulse");
    old_hist->Draw("HIST");
    new_hist->Draw("HIST SAME");
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(old_hist, Form("Old, outside 2ns: %.3f, outside 4ns: %.3f", fraction_old_outside_2, fraction_old_outside_4), "l");
    legend->AddEntry(new_hist, Form("New, outside 2ns: %.3f, outside 4ns: %.3f", fraction_new_outside_2, fraction_new_outside_4), "l");
    legend->Draw();
}