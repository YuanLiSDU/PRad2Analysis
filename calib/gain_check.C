

void gain_check(){

    TFile *f = TFile::Open("../data/gain_check/prad_024637_gain_corr.root");
    TTree *t = (TTree*)f->Get("gain_corr");

    float refPMT_ratio[3];
    int   batch_id;
    t->SetBranchAddress("refPMT_ratio", refPMT_ratio);
    t->SetBranchAddress("batch_id",     &batch_id);

    // Histogram ranges chosen around the observed value spans
    const char *titles[3] = {
        "refPMT_ratio[0];Ratio;Counts",
        "refPMT_ratio[1];Ratio;Counts",
        "refPMT_ratio[2];Ratio;Counts"
    };
    double lo[3] = { 0.870, 2.355, 0.460 };
    double hi[3] = { 0.880, 2.375, 0.468 };

    TH1F *h[3];
    for(int k = 0; k < 3; k++)
        h[k] = new TH1F(Form("h%d", k), titles[k], 50, lo[k], hi[k]);

    for(Long64_t i = 0; i < t->GetEntries(); i++){
        t->GetEntry(i);
        for(int k = 0; k < 3; k++)
            h[k]->Fill(refPMT_ratio[k]);
    }

    TCanvas *c = new TCanvas("c", "refPMT_ratio", 1200, 400);
    c->Divide(3, 1);
    for(int k = 0; k < 3; k++){
        c->cd(k + 1);
        h[k]->SetLineColor(kBlue + k);
        h[k]->SetLineWidth(2);
        h[k]->Draw("HIST");
        // draw mean line
        double mean = h[k]->GetMean();
        TLine *l = new TLine(mean, 0, mean, h[k]->GetMaximum());
        l->SetLineColor(kRed);
        l->SetLineStyle(2);
        l->Draw();
        TLatex *tx = new TLatex();
        tx->SetNDC();
        tx->SetTextSize(0.045);
        tx->DrawLatex(0.15, 0.85, Form("mean = %.4f", mean));
    }
    //c->SaveAs("refPMT_ratio.pdf");
}