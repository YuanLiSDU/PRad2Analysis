#include "../EventData.h"
#include "../PhysicsTools.h"

void integral_study(){

    TChain *tree = new TChain("events");
    for(int i=0; i<10; i++) 
        tree->Add(Form("../data/integral_study/lms_23915_0.1_raw_%d.root", i+1));

    RawEventData *evp = new RawEventData();
    RawEventData &ev = *evp;
    setupRawBranches(tree, ev);

    TChain *tree2 = new TChain("events");
    for(int i=0; i<10; i++)
        tree2->Add(Form("../data/integral_study/lms_23915_0.06_raw_%d.root", i+1));
    RawEventData ev2;
    setupRawBranches(tree2, ev2);

    //histograms
    TH1F *waveforms[6];
    for (int i=0; i<6; i++){
        waveforms[i] = new TH1F(Form("waveform_%d", i), Form("Waveform_%d;Sample;ADC", i), 100, -0.5, 100-0.5);
    }
    TH1F *smooth_integral = new TH1F("smooth_integral", "Smoothed Integral;Integral;Counts", 300*3, 0, 14000);
    TH1F *raw_integral = new TH1F("raw_integral", "Raw Integral;Integral;Counts", 300*3, 0, 14000);

    TH1F *smooth_integral_0p06 = new TH1F("smooth_integral_0p06", "Smoothed Integral (0.06 threshold);Integral;Counts", 300*3, 0, 14000);
    TH1F *raw_integral_0p06 = new TH1F("raw_integral_0p06", "Raw Integral (0.06 threshold);Integral;Counts", 300*3, 0, 14000);

    int aid_module = 1000 + 727; // W727
    int count = 0;
    float peak_heights[6] = {}, ped_means[6] = {};

    for (Long64_t i=0; i<tree->GetEntries(); i++){
        tree->GetEntry(i);
        cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
        if(ev.nch < 1000) continue;
        for (int ch=0; ch<ev.nch; ch++){
            if (ev.module_id[ch] == aid_module){
                if(count < 6) { // only fill first 6 waveforms for testing
                    for (int s=0; s<ev.nsamples[ch]; s++){
                        waveforms[count]->SetBinContent(s+1, ev.samples[ch][s]);
                    }
                    // store the first peak height
                    if (ev.npeaks[ch] > 0){
                        peak_heights[count] = ev.peak_height[ch][0];
                        ped_means[count] = ev.ped_mean[ch];
                    }
                    count++;
                }
                if(ev.npeaks[ch] == 1){
                    smooth_integral->Fill(ev.peak_integral[ch][0]-3000.);
                    // find peak sample index
                    int peak_idx = 0;
                    float max_val = -1e9;
                    for (int s = 0; s < ev.nsamples[ch]; s++){
                        float val = ev.samples[ch][s] - ev.ped_mean[ch];
                        if (val > max_val){ max_val = val; peak_idx = s; }
                    }
                    // --- 10% threshold ---
                    {
                        float thr = 0.10f * ev.peak_height[ch][0];
                        int left = peak_idx;
                        while (left > 0 && (ev.samples[ch][left-1] - ev.ped_mean[ch]) >= thr) left--;
                        int right = peak_idx;
                        while (right < ev.nsamples[ch]-1 && (ev.samples[ch][right+1] - ev.ped_mean[ch]) >= thr) right++;
                        float raw_int = 0;
                        for (int s = left; s <= right; s++) raw_int += (ev.samples[ch][s] - ev.ped_mean[ch]);
                        raw_integral->Fill(raw_int*3956./3948.-2000.);
                    }
                }
            }
        }
    }

    cout << "Done processing " << tree->GetEntries() << " events." << endl;

    for(Long64_t i=0; i<tree2->GetEntries(); i++){
        tree2->GetEntry(i);
        cout << "Processing event " << i << " / " << tree2->GetEntries() << "\r" << flush;
        if(ev2.nch < 1000) continue;
        for (int ch=0; ch<ev2.nch; ch++){
            if (ev2.module_id[ch] == aid_module){
                if(ev2.npeaks[ch] == 1){
                    smooth_integral_0p06->Fill(ev2.peak_integral[ch][0]*3956./4179.-1000.);
                    int peak_idx = 0;
                    float max_val = -1e9;
                    for (int s = 0; s < ev2.nsamples[ch]; s++){
                        float val = ev2.samples[ch][s] - ev2.ped_mean[ch];
                        if (val > max_val){ max_val = val; peak_idx = s; }
                    }
                    // --- 5% threshold ---
                    {
                        float thr = 0.05f * ev2.peak_height[ch][0];
                        int left = peak_idx;
                        while (left > 0 && (ev2.samples[ch][left-1] - ev2.ped_mean[ch]) >= thr) left--;
                        int right = peak_idx;
                        while (right < ev2.nsamples[ch]-1 && (ev2.samples[ch][right+1] - ev2.ped_mean[ch]) >= thr) right++;
                        float raw_int = 0;
                        for (int s = left; s <= right; s++) raw_int += (ev2.samples[ch][s] - ev2.ped_mean[ch]);
                        raw_integral_0p06->Fill(raw_int*3956./4204.);
                    }
                }
            }
        }
    }

    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    c->Divide(3,2);
    for (int i=0; i<6; i++){
        c->cd(i+1);
        waveforms[i]->Draw("hist");

        double xmin = waveforms[i]->GetXaxis()->GetXmin();
        double xmax = waveforms[i]->GetXaxis()->GetXmax();

        TLine *l10 = new TLine(xmin, 0.10 * peak_heights[i] + ped_means[i], xmax, 0.10 * peak_heights[i] + ped_means[i]);
        l10->SetLineColor(kRed);
        l10->SetLineStyle(2); // dashed
        l10->Draw("same");

        TLine *l2 = new TLine(xmin, 0.06 * peak_heights[i] + ped_means[i], xmax, 0.06 * peak_heights[i] + ped_means[i]);
        l2->SetLineColor(kBlue);
        l2->SetLineStyle(2); // dashed
        l2->Draw("same");
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 500);
    smooth_integral->SetLineColor(kRed);
    smooth_integral->SetLineStyle(1);
    raw_integral->SetLineColor(kBlue);
    raw_integral->SetLineStyle(1);
    smooth_integral_0p06->SetLineColor(kRed);
    smooth_integral_0p06->SetLineStyle(2);
    raw_integral_0p06->SetLineColor(kBlue);
    raw_integral_0p06->SetLineStyle(2);

    smooth_integral->Fit("gaus", "RQ", "", smooth_integral->GetMean() - 2 * smooth_integral->GetRMS(), smooth_integral->GetMean() + 3 * smooth_integral->GetRMS());
    raw_integral->Fit("gaus", "RQ", "", raw_integral->GetMean() - 2 * raw_integral->GetRMS(), raw_integral->GetMean() + 3 * raw_integral->GetRMS());
    smooth_integral_0p06->Fit("gaus", "RQ", "", smooth_integral_0p06->GetMean() - 2 * smooth_integral_0p06->GetRMS(), smooth_integral_0p06->GetMean() + 3 * smooth_integral_0p06->GetRMS());
    raw_integral_0p06->Fit("gaus", "RQ", "", raw_integral_0p06->GetMean() - 2 * raw_integral_0p06->GetRMS(), raw_integral_0p06->GetMean() + 3 * raw_integral_0p06->GetRMS());

    // draw — find global max for axis range
    double ymax = 0;
    for (auto *h : {smooth_integral, raw_integral, smooth_integral_0p06, raw_integral_0p06})
        if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
    smooth_integral->SetMaximum(ymax * 1.15);

    smooth_integral->Draw("hist");
    raw_integral->Draw("hist same");
    smooth_integral_0p06->Draw("hist same");
    raw_integral_0p06->Draw("hist same");

    // draw fit functions with matching color/style
    auto drawFit = [](TH1F *h){
        TF1 *f = h->GetFunction("gaus");
        if(!f) return;
        f->SetLineColor(h->GetLineColor());
        f->SetLineStyle(h->GetLineStyle());
        f->SetLineWidth(2);
        f->Draw("same");
    };
    drawFit(smooth_integral);
    drawFit(raw_integral);
    drawFit(smooth_integral_0p06);
    drawFit(raw_integral_0p06);

    auto getP = [](TH1F *h, int p){ return h->GetFunction("gaus")->GetParameter(p); };
    TLegend *leg = new TLegend(0.55, 0.65, 0.92, 0.92);
    leg->AddEntry(smooth_integral,      Form("Smooth 10%%: %.1f #pm %.1f", getP(smooth_integral,1),      getP(smooth_integral,2)),      "l");
    leg->AddEntry(raw_integral,         Form("Raw    10%%: %.1f #pm %.1f", getP(raw_integral,1),         getP(raw_integral,2)),         "l");
    leg->AddEntry(smooth_integral_0p06, Form("Smooth  5%%: %.1f #pm %.1f", getP(smooth_integral_0p06,1), getP(smooth_integral_0p06,2)), "l");
    leg->AddEntry(raw_integral_0p06,    Form("Raw     5%%: %.1f #pm %.1f", getP(raw_integral_0p06,1),    getP(raw_integral_0p06,2)),    "l");
    leg->Draw();
}