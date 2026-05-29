
int mod = 565+34*6;

void gain_factor(){

    const char *dir = "../data/gain_corr";
    std::vector<TString> files;

    TSystemDirectory sysdir(dir, dir);
    TList *flist = sysdir.GetListOfFiles();
    TIter next(flist);
    TSystemFile *sfile;
    while ((sfile = (TSystemFile*)next())) {
        TString fname = sfile->GetName();
        if (fname.EndsWith("_corr.root"))
            files.push_back(TString(dir) + "/" + fname);
    }
    std::sort(files.begin(), files.end());

    std::vector<double> vals[3];
    std::vector<double> pmt[3];
    std::vector<double> ref_lms, ref_alpha;
    float gain_corr_W[1156][3];
    float refPMT_ratio[3];
    float fit_mean_ref_lms[3];
    float fit_mean_ref_alpha[3];

    for (const TString &fpath : files) {
        TFile *f = TFile::Open(fpath);
        if (!f || f->IsZombie()) { delete f; continue; }
        TTree *t = (TTree*)f->Get("gain_corr");
        if (!t) { delete f; continue; }

        t->SetBranchAddress("gain_corr_W", gain_corr_W);
        t->SetBranchAddress("refPMT_ratio", refPMT_ratio);
        t->SetBranchAddress("fit_mean_ref_lms", fit_mean_ref_lms);
        t->SetBranchAddress("fit_mean_ref_alpha", fit_mean_ref_alpha);
        Long64_t nentries = t->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            for (int k = 0; k < 3; k++) {
                vals[k].push_back(gain_corr_W[mod-1][k]);
                pmt[k].push_back(refPMT_ratio[k]);
            }
            ref_lms.push_back(fit_mean_ref_lms[0]);
            ref_alpha.push_back(fit_mean_ref_alpha[0]);
        }
        delete f;
    }

    int N = (int)vals[0].size();
    Color_t colors[3] = {kBlack, kRed, kBlue};
    TGraph *gr[3];
    for (int k = 0; k < 3; k++) {
        gr[k] = new TGraph(N);
        for (int i = 0; i < N; i++)
            gr[k]->SetPoint(i, i + 1, vals[k][i]);
        gr[k]->SetMarkerStyle(20 + k);
        gr[k]->SetMarkerSize(0.7);
        gr[k]->SetMarkerColor(colors[k]);
        gr[k]->SetLineColor(colors[k]);
    }

    TGraph *gr_avg = new TGraph(N);
    for (int i = 0; i < N; i++) {
        double avg = (vals[1][i] + vals[2][i]) / 2.0;
        gr_avg->SetPoint(i, i + 1, avg);
    }
    gr_avg->SetLineColor(kMagenta+1);
    gr_avg->SetLineWidth(2);

    double ymin_g = 1e9, ymax_g = -1e9;
    for (int k = 1; k < 3; k++)
        for (int i = 0; i < N; i++) {
            ymin_g = std::min(ymin_g, vals[k][i]);
            ymax_g = std::max(ymax_g, vals[k][i]);
        }
    double margin_g = (ymax_g - ymin_g) * 0.1;

    TCanvas *c = new TCanvas("c_gain", Form("Gain correction W[%d]", mod), 1000, 500);
    gr[1]->SetTitle(Form("gain_corr_W%d vs time samples;Entry index;gain_corr_W%d", mod, mod));
    gr[1]->GetYaxis()->SetRangeUser(ymin_g - margin_g, ymax_g + margin_g);
    gr[1]->Draw("AP");
    gr[2]->Draw("P SAME");
    gr_avg->Draw("L SAME");

    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    for (int k = 1; k < 3; k++)
        leg->AddEntry(gr[k], Form("Ref PMT %d", k+1), "p");
    {
        double avg_mean = 0., avg_rms = 0.;
        std::vector<double> avgvals(N);
        for (int i = 0; i < N; i++) avgvals[i] = (vals[1][i] + vals[2][i]) / 2.0;
        for (int i = 0; i < N; i++) avg_mean += avgvals[i];
        avg_mean /= N;
        for (int i = 0; i < N; i++) avg_rms += (avgvals[i] - avg_mean) * (avgvals[i] - avg_mean);
        avg_rms = std::sqrt(avg_rms / N);
        double pct = (avg_mean != 0.) ? avg_rms / std::fabs(avg_mean) * 100. : 0.;
        leg->AddEntry(gr_avg, Form("average  RMS/mean = %.2f%%", pct), "l");
    }
    leg->Draw();

    c->SaveAs(Form("%d_gain_factor.png", mod));
    
    // --- refPMT_ratio plot ---
    TGraph *gr_pmt[3];
    for (int k = 0; k < 3; k++) {
        gr_pmt[k] = new TGraph(N);
        for (int i = 0; i < N; i++)
            gr_pmt[k]->SetPoint(i, i + 1, pmt[k][i]);
        gr_pmt[k]->SetMarkerStyle(20 + k);
        gr_pmt[k]->SetMarkerSize(0.7);
        gr_pmt[k]->SetMarkerColor(colors[k]);
        gr_pmt[k]->SetLineColor(colors[k]);
    }

    TCanvas *c2 = new TCanvas("c_pmt", "refPMT_ratio", 1000, 500);
    gr_pmt[0]->SetTitle("refPMT_ratio(LMS / alpha) vs time samples(every 2000 LMS events);Sample Num;refPMT_ratio");

    double ymin = 1e9, ymax = -1e9;
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < N; i++) {
            ymin = std::min(ymin, pmt[k][i]);
            ymax = std::max(ymax, pmt[k][i]);
        }
    double margin = (ymax - ymin) * 0.1;
    gr_pmt[0]->GetYaxis()->SetRangeUser(ymin - margin, ymax + margin);

    gr_pmt[0]->Draw("AP");
    gr_pmt[1]->Draw("P SAME");
    gr_pmt[2]->Draw("P SAME");

    TLegend *leg2 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg2->SetBorderSize(0);
    for (int k = 0; k < 3; k++) {
        double mean = 0., rms = 0.;
        for (int i = 0; i < N; i++) mean += pmt[k][i];
        mean /= N;
        for (int i = 0; i < N; i++) rms += (pmt[k][i] - mean) * (pmt[k][i] - mean);
        rms = std::sqrt(rms / N);
        double pct = (mean != 0.) ? rms / std::fabs(mean) * 100. : 0.;
        leg2->AddEntry(gr_pmt[k], Form("Ref PMT %d  RMS/mean = %.2f%%", k+1, pct), "p");
    }
    leg2->Draw();
    c2->SaveAs("refPMT_ratio.png");

    // --- Ref PMT1 LMS mean vs alpha mean ---
    TGraph *gr_lms   = new TGraph(N);
    TGraph *gr_alpha = new TGraph(N);
    for (int i = 0; i < N; i++) {
        gr_lms->SetPoint(i,   i + 1, ref_lms[i]);
        gr_alpha->SetPoint(i, i + 1, ref_alpha[i]);
    }
    gr_lms->SetMarkerStyle(20); gr_lms->SetMarkerSize(0.7);
    gr_lms->SetMarkerColor(kBlue);  gr_lms->SetLineColor(kBlue);
    gr_alpha->SetMarkerStyle(21); gr_alpha->SetMarkerSize(0.7);
    gr_alpha->SetMarkerColor(kRed); gr_alpha->SetLineColor(kRed);

    double ymin3 = 1e9, ymax3 = -1e9;
    for (int i = 0; i < N; i++) {
        ymin3 = std::min({ymin3, ref_lms[i], ref_alpha[i]});
        ymax3 = std::max({ymax3, ref_lms[i], ref_alpha[i]});
    }
    double margin3 = (ymax3 - ymin3) * 0.1;

    TCanvas *c3 = new TCanvas("c_ref1", "Ref PMT1 LMS & Alpha mean", 1000, 500);
    gr_lms->SetTitle("Ref PMT1 fit mean vs sample;Sample Num;Fit mean (ADC)");
    gr_lms->GetYaxis()->SetRangeUser(ymin3 - margin3, ymax3 + margin3);
    gr_lms->Draw("AP");
    gr_alpha->Draw("P SAME");

    TLegend *leg3 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg3->SetBorderSize(0);
    leg3->AddEntry(gr_lms,   "LMS mean",   "p");
    leg3->AddEntry(gr_alpha, "Alpha mean", "p");
    leg3->Draw();
    c3->SaveAs("refPMT1_mean.png");

    std::cout << "Total entries: " << N << std::endl;

/*
    // --- peak_time scatter plot from _lms.root files ---
    {
        const char *dir_lms = "../data/gain_corr";
        std::vector<TString> files_lms;
        {
            TSystemDirectory sd(dir_lms, dir_lms);
            TList *fl = sd.GetListOfFiles();
            TIter nx(fl);
            TSystemFile *sf;
            while ((sf = (TSystemFile*)nx())) {
                TString fn = sf->GetName();
                if (fn.EndsWith("_lms.root"))
                    files_lms.push_back(TString(dir_lms) + "/" + fn);
            }
        }
        std::sort(files_lms.begin(), files_lms.end());

        const int W_ID_BASE  = 1000;
        const int BATCH_SIZE = 2000;
        const int target_id  = W_ID_BASE + mod; // module_id in tree

        std::vector<double> pt_batch;   // peak_time values in current batch
        std::vector<double> pt_means;   // mean per batch
        std::vector<double> ph_batch;   // peak_height values in current batch
        std::vector<double> ph_means;   // mean per batch
        int lms_count = 0;

        int   event_type, nch;
        unsigned short module_id[2000];
        unsigned char  module_type[2000], npeaks[2000];
        float peak_time[2000][8];
        float peak_height[2000][8];

        for (const TString &fp : files_lms) {
            TFile *f = TFile::Open(fp);
            if (!f || f->IsZombie()) { delete f; continue; }
            TTree *t = (TTree*)f->Get("lms_gain");
            if (!t) { delete f; continue; }
            t->SetMakeClass(1);
            t->SetBranchAddress("event_type",       &event_type);
            t->SetBranchAddress("hycal.nch",        &nch);
            t->SetBranchAddress("hycal.module_id",   module_id);
            t->SetBranchAddress("hycal.module_type", module_type);
            t->SetBranchAddress("hycal.npeaks",      npeaks);
            t->SetBranchAddress("hycal.peak_time",   peak_time);
            t->SetBranchAddress("hycal.peak_height", peak_height);

            Long64_t nentries = t->GetEntries();
            for (Long64_t i = 0; i < nentries; i++) {
                t->GetEntry(i);
                const bool is_lms = (event_type == 0);
                if (!is_lms) continue;

                for (int ich = 0; ich < nch; ++ich) {
                    if ((int)module_type[ich] != 2) continue;
                    if ((int)module_id[ich]   != target_id) continue;
                    if ((int)npeaks[ich]       != 1) continue;
                    pt_batch.push_back(peak_time[ich][0]);
                    ph_batch.push_back(peak_height[ich][0]);
                    break;
                }

                ++lms_count;
                if (lms_count % BATCH_SIZE == 0) {
                    if (!pt_batch.empty()) {
                        double sum = 0.;
                        for (double v : pt_batch) sum += v;
                        pt_means.push_back(sum / pt_batch.size());
                        pt_batch.clear();
                    }
                    if (!ph_batch.empty()) {
                        double sum = 0.;
                        for (double v : ph_batch) sum += v;
                        ph_means.push_back(sum / ph_batch.size());
                        ph_batch.clear();
                    }
                }
            }
            delete f;
        }
        // flush remaining
        if (!pt_batch.empty()) {
            double sum = 0.;
            for (double v : pt_batch) sum += v;
            pt_means.push_back(sum / pt_batch.size());
        }
        if (!ph_batch.empty()) {
            double sum = 0.;
            for (double v : ph_batch) sum += v;
            ph_means.push_back(sum / ph_batch.size());
        }

        int Npt = (int)pt_means.size();
        TGraph *gr_pt = new TGraph(Npt);
        for (int i = 0; i < Npt; i++)
            gr_pt->SetPoint(i, i + 1, pt_means[i]);
        gr_pt->SetMarkerStyle(20);
        gr_pt->SetMarkerSize(0.7);
        gr_pt->SetMarkerColor(kBlack);

        TCanvas *c4 = new TCanvas("c_pt", Form("Peak time W%d", mod), 1000, 500);
        gr_pt->SetTitle(Form("W%d peak time (per %d LMS events);Batch index;Mean peak time (ns)", mod, BATCH_SIZE));
        gr_pt->Draw("AP");
        c4->SaveAs(Form("%d_peak_time.png", mod));

        int Nph = (int)ph_means.size();
        TGraph *gr_ph = new TGraph(Nph);
        for (int i = 0; i < Nph; i++)
            gr_ph->SetPoint(i, i + 1, ph_means[i]);
        gr_ph->SetMarkerStyle(20);
        gr_ph->SetMarkerSize(0.7);
        gr_ph->SetMarkerColor(kBlack);

        TCanvas *c5 = new TCanvas("c_ph", Form("Peak height W%d", mod), 1000, 500);
        gr_ph->SetTitle(Form("W%d peak height (per %d LMS events);Batch index;Mean peak height (ADC)", mod, BATCH_SIZE));
        gr_ph->Draw("AP");
        c5->SaveAs(Form("%d_peak_height.png", mod));
        std::cout << "Peak time batches: " << Npt << "  total LMS: " << lms_count << std::endl;
    }

    // --- HV plot ---
    {
        const char *dir_hv = "../data/gain_corr";
        std::vector<TString> files_hv;
        {
            TSystemDirectory sd(dir_hv, dir_hv);
            TList *fl = sd.GetListOfFiles();
            TIter nx(fl);
            TSystemFile *sf;
            while ((sf = (TSystemFile*)nx())) {
                TString fn = sf->GetName();
                if (fn.EndsWith("_lms_hv.root"))
                    files_hv.push_back(TString(dir_hv) + "/" + fn);
            }
        }
        std::sort(files_hv.begin(), files_hv.end());

        if (files_hv.empty()) {
            std::cout << "No HV files found." << std::endl;
            return;
        }

        // Find array index for this module in hv_channels
        int hv_cidx = -1;
        {
            TFile *f0 = TFile::Open(files_hv[0]);
            if (f0 && !f0->IsZombie()) {
                TTree *hvc = (TTree*)f0->Get("hv_channels");
                if (hvc) {
                    unsigned short cid;
                    char cname[64];
                    hvc->SetBranchAddress("channel_id", &cid);
                    hvc->SetBranchAddress("name", cname);
                    TString target = Form("W%d", mod);
                    for (Long64_t i = 0; i < hvc->GetEntries(); i++) {
                        hvc->GetEntry(i);
                        if (TString(cname) == target) {
                            hv_cidx = (int)cid;
                            break;
                        }
                    }
                }
                delete f0;
            }
        }
        if (hv_cidx < 0) {
            std::cout << "HV channel for W" << mod << " not found." << std::endl;
            return;
        }
        std::cout << "HV channel index for W" << mod << ": " << hv_cidx << std::endl;

        std::vector<double> hv_t, hv_v;
        double t_unix_s;
        float v0set[1792], dv_arr[1792];

        for (const TString &fp : files_hv) {
            TFile *f = TFile::Open(fp);
            if (!f || f->IsZombie()) { delete f; continue; }
            TTree *hvt = (TTree*)f->Get("hv");
            if (!hvt) { delete f; continue; }
            hvt->SetBranchAddress("t_unix_s", &t_unix_s);
            hvt->SetBranchAddress("v0set",    v0set);
            hvt->SetBranchAddress("dv",       dv_arr);
            Long64_t nentries = hvt->GetEntries();
            for (Long64_t i = 0; i < nentries; i++) {
                hvt->GetEntry(i);
                hv_t.push_back(t_unix_s);
                hv_v.push_back((double)v0set[hv_cidx] + (double)dv_arr[hv_cidx]);
            }
            delete f;
        }

        int Nhv = (int)hv_t.size();
        if (Nhv == 0) { std::cout << "No HV data." << std::endl; return; }

        double t0_hv = hv_t[0];
        TGraph *gr_hv = new TGraph(Nhv);
        for (int i = 0; i < Nhv; i++)
            gr_hv->SetPoint(i, (hv_t[i] - t0_hv) / 60., hv_v[i]);
        gr_hv->SetMarkerStyle(20);
        gr_hv->SetMarkerSize(0.2);
        gr_hv->SetMarkerColor(kBlue+1);
        gr_hv->SetLineColor(kBlue+1);
        gr_hv->SetLineWidth(1);

        TCanvas *c_hv = new TCanvas("c_hv", Form("HV W%d", mod), 1000, 500);
        gr_hv->SetTitle(Form("HV (V0Set + dV) for W%d vs time;Time (min);Voltage (V)", mod));
        gr_hv->Draw("AL");
        c_hv->SaveAs(Form("%d_hv.png", mod));
        std::cout << "HV entries: " << Nhv << std::endl;
    }
*/
}