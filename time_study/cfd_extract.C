#include "../EventData.h"
#include "waveanalyzer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"

namespace fs = std::filesystem;

static std::vector<std::string> collectRootFiles(const std::string &path);
static bool processWaveformDemoAndPlot(const float waveform[100],
                                       float ped_mean,
                                       float ped_rms,
                                       float clk_mhz,
                                       float break_ratio,
                                       int smooth_order,
                                       int module_id,
                                       const std::string &trigger_sel,
                                       const std::string &pdf_path,
                                       Long64_t entry,
                                       int ch);

void cfd_extract(const char *pdf_output = "time_demo.pdf", int max_plots = 10) {

    std::string file_path = "../data/time_study/";
    std::vector<std::string> root_files = collectRootFiles(file_path);

    TChain chain("events");
    for (const auto &file : root_files) {
        chain.Add(file.c_str());
    }

    RawEventData ev;
    setupRawBranches(&chain, ev);

    int plotted = 0;
    const std::string out_pdf = std::string(pdf_output);
    TCanvas pdf_bookend("pdf_bookend", "pdf_bookend", 1100, 700);
    pdf_bookend.SetCanvasSize(1100, 700);
    pdf_bookend.Print((out_pdf + "[").c_str());

    TH1F *time_diff_0p5 = new TH1F("time_diff_0p5", "Time Diff 0.5", 800, -10, 10);
    TH1F *time_diff_0p4 = new TH1F("time_diff_0p4", "Time Diff 0.4", 800, -10, 10);
    TH1F *time_diff_0p3 = new TH1F("time_diff_0p3", "Time Diff 0.3", 800, -10, 10);
    TH1F *time_diff_0p2 = new TH1F("time_diff_0p2", "Time Diff 0.2", 800, -10, 10);

    TH1F *time_diff_0p2_linear = new TH1F("time_diff_0p2_linear", "Time Diff 0.2 Linear", 800, -10, 10);


    int module[2] = {1492-34, 1493-34};
    //int module[2] = {3003, 3002}; //veto

    Long64_t nentries = chain.GetEntries();
    for (Long64_t i = 0; i < nentries/2; ++i) {
        chain.GetEntry(i);
        // Process the event data here

        if (i % 1000 == 0) std::cout << "Processing entry " << i << "/" << nentries << std::endl;

        bool is_sum = ((ev.trigger_bits & 1u << 8) != 0 && ev.nch < 100);
        bool is_lms = ((ev.trigger_bits & 1u << 24) != 0) && ev.nch > 900;

        float cfd_time[2][4] = {}, peak_h[2] = {}, cfd_time_linear[2] = {};

        for (int j = 0; j < ev.nch; j++){
            if (ev.npeaks[j] != 1) continue;

            bool is_pwo = (ev.module_id[j] > 1000 && ev.module_id[j] <= 2156);
            bool is_veto = (ev.module_id[j] >=3001 && ev.module_id[j] <= 3004);

            if (!is_pwo) continue;
            if (!is_sum) continue;
            //if (ev.peak_height[j][0] <= 300.f) continue;

            if (ev.module_id[j] != module[0] && ev.module_id[j] != module[1]) continue;

            float waveform[100];
            for(int k = 0; k < 100; k++){
                waveform[k] = ev.samples[j][k];
            }

            const float ped_mean = ev.ped_mean[j];
            const float ped_rms = ev.ped_rms[j];
            const float break_ratio = 0.1f;
            const int smooth_order = 1;
            const float clk_mhz = 250.0f;

            const SimplePeakWindowResult win = findSinglePeakWindowByRatio(
                waveform, ped_mean, ped_rms, break_ratio, smooth_order);
            if (!win.ok) continue;

            SimpleLogNormalFitResult t_logn[4];
            t_logn[0] = findTimeLogNormalCFDNs(
                waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.5f);
            t_logn[1] = findTimeLogNormalCFDNs(
                waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.4f);
            t_logn[2] = findTimeLogNormalCFDNs(
                waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.3f);
            t_logn[3] = findTimeLogNormalCFDNs(
                waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.2f);
            const SimpleTimeResult t_lin = findTimeLinearCFDNs(
                waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.2f);

            if (ev.module_id[j] == module[0] && ev.peak_height[j][0] > 200.f) {
                for(int t = 0; t < 4; t++){
                    if (t_logn[t].ok) cfd_time[0][t] = t_logn[t].time_ns;
                }
                if (t_lin.ok) cfd_time_linear[0] = t_lin.time_ns;
                peak_h[0] = ev.peak_height[j][0];
            }
            if (ev.module_id[j] == module[1] && ev.peak_height[j][0] > 200.f) {
                for(int t = 0; t < 4; t++){
                    if (t_logn[t].ok) cfd_time[1][t] = t_logn[t].time_ns;
                }
                if (t_lin.ok) cfd_time_linear[1] = t_lin.time_ns;
                peak_h[1] = ev.peak_height[j][0];
            }

            if (peak_h[0] * peak_h[1] > 0.f && cfd_time[0][3] - cfd_time[1][3] < 3.f) {
                if(plotted <= max_plots && ev.module_id[j] == module[0]){
                    const std::string trigger_sel = is_lms ? "LMS" : (is_sum ? "SUM" : "other");
                    const bool made_plot = processWaveformDemoAndPlot(waveform, ped_mean, ped_rms, clk_mhz,
                        break_ratio, smooth_order, ev.module_id[j], trigger_sel, out_pdf, i, j);
                    if (made_plot) {
                        ++plotted;
                    }
                }
            }
        }
        if (peak_h[0] * peak_h[1] > 0.f){
            time_diff_0p5 -> Fill(cfd_time[0][0] - cfd_time[1][0]);
            time_diff_0p4 -> Fill(cfd_time[0][1] - cfd_time[1][1]);
            time_diff_0p3 -> Fill(cfd_time[0][2] - cfd_time[1][2]);
            time_diff_0p2 -> Fill(cfd_time[0][3] - cfd_time[1][3]);
            time_diff_0p2_linear -> Fill(cfd_time_linear[0] - cfd_time_linear[1]);
        }
    }

    pdf_bookend.Print((out_pdf + "]").c_str());
    std::cout << "[time-demo] total plotted pages: " << plotted
              << " -> " << out_pdf << std::endl;

    TCanvas *c1 = new TCanvas("scan_ratio", "scan_ratio", 800, 600);
    c1->cd();
    time_diff_0p5->SetTitle(Form("Time Differences of Mod %d/%d (Diff CFD Ratios)", module[0], module[1]));
    time_diff_0p5->SetStats(0);
    time_diff_0p5->GetYaxis()->SetRangeUser(0, time_diff_0p2->GetMaximum()*1.2);
    time_diff_0p5->SetLineColor(kRed);
    time_diff_0p5->SetLineWidth(1);
    time_diff_0p5->Draw();
    time_diff_0p4->SetLineColor(kBlue);
    time_diff_0p4->SetLineWidth(1);
    time_diff_0p4->Draw("SAME");
    time_diff_0p3->SetLineColor(kGreen);
    time_diff_0p3->SetLineWidth(1);
    time_diff_0p3->Draw("SAME");
    time_diff_0p2->SetLineColor(kMagenta);
    time_diff_0p2->SetLineWidth(1);
    time_diff_0p2->Draw("SAME");

    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(time_diff_0p5, "0.5", "l");
    legend->AddEntry(time_diff_0p4, "0.4", "l");
    legend->AddEntry(time_diff_0p3, "0.3", "l");
    legend->AddEntry(time_diff_0p2, "0.2", "l");
    legend->Draw();

    TCanvas *c2 = new TCanvas("linear_vs_fitting", "linear_vs_fitting", 800, 600);
    c2->cd();
    time_diff_0p2_linear->SetTitle(Form("Time Differences of Mod %d/%d (Linear CFD vs Fitting)", module[0], module[1]));
    time_diff_0p2_linear->SetStats(0);
    time_diff_0p2_linear->GetYaxis()->SetRangeUser(0, time_diff_0p2_linear->GetMaximum()*1.2);
    time_diff_0p2_linear->SetLineColor(kCyan);
    time_diff_0p2_linear->SetLineWidth(1);
    time_diff_0p2_linear->Draw();
    time_diff_0p2->SetLineColor(kMagenta);
    time_diff_0p2->SetLineWidth(1);
    time_diff_0p2->Draw("SAME");
    TLegend *legend2 = new TLegend(0.7,0.7,0.9,0.9);
    legend2->AddEntry(time_diff_0p2_linear, "0.2 Linear", "l");
    legend2->AddEntry(time_diff_0p2, "0.2 Fitting", "l");
    legend2->Draw();
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find(".root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

static bool processWaveformDemoAndPlot(const float waveform[100],
                                       float ped_mean,
                                       float ped_rms,
                                       float clk_mhz,
                                       float break_ratio,
                                       int smooth_order,
                                       int module_id,
                                       const std::string &trigger_sel,
                                       const std::string &pdf_path,
                                       Long64_t entry,
                                       int ch)
{
    const SimplePeakWindowResult win = findSinglePeakWindowByRatio(
        waveform, ped_mean, ped_rms, break_ratio, smooth_order);
    if (!win.ok) return false;

    const SimpleTimeResult t_quad = findTimeQuadraticNs(
        waveform, win.peak_pos, clk_mhz);
    const SimpleTimeResult t_cfd = findTimeLinearCFDNs(
        waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.2f);
    const SimpleLogNormalFitResult t_logn = findTimeLogNormalCFDNs(
        waveform, ped_mean, ped_rms, win.left, win.peak_pos, clk_mhz, 0.2f);
    const int fit_start = std::max(0, win.left - 5);
    const int fit_stop = std::min(99, win.peak_pos + 4);
    const int fit_n = std::max(0, fit_stop - fit_start + 1);

    std::cout
        << "[time-demo] entry=" << entry
        << " ch=" << ch
        << " peak=" << win.peak_pos
        << " left=" << win.left
        << " right=" << win.right
        << " height=" << win.height
        << "\n"
        << "  quadratic_time_ns=" << t_quad.time_ns
        << " (ok=" << t_quad.ok << ")\n"
        << "  linear_cfd_time_ns=" << t_cfd.time_ns
        << " (ok=" << t_cfd.ok
        << ", fallback=" << t_cfd.used_fallback << ")\n"
        << "  lognormal_cfd_time_ns=" << t_logn.time_ns
        << " (ok=" << t_logn.ok
        << ", chi2_ndf=" << t_logn.chi2_per_dof << ")\n"
        << "  lognormal_eq: "
        << (t_logn.ok ? t_logn.equation : (t_logn.fail_reason.empty() ? "fit failed" : t_logn.fail_reason))
        << std::endl;

    double x[100];
    double y[100];
    float ymin = waveform[0];
    float ymax = waveform[0];
    for (int i = 0; i < 100; ++i) {
        x[i] = sampleToNs(static_cast<float>(i), clk_mhz);
        y[i] = waveform[i];
        ymin = std::min(ymin, waveform[i]);
        ymax = std::max(ymax, waveform[i]);
    }

    const float pad = std::max(10.0f, 0.08f * (ymax - ymin));
    float ylow = ymin - pad;
    float yhigh = ymax + pad;

    if (t_logn.ok) {
        for (int i = fit_start; i <= fit_stop; ++i) {
            const float fit_y = evalLogNormalFitAtSample(static_cast<float>(i), ped_mean, t_logn);
            ylow = std::min(ylow, fit_y);
            yhigh = std::max(yhigh, fit_y);
        }
        ylow -= pad;
        yhigh += pad;
    }

    TCanvas c("c_time_demo", "Waveform Time Demo", 1100, 700);
    c.SetCanvasSize(1100, 700);
    c.SetLeftMargin(0.10);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.11);
    c.SetTopMargin(0.08);
    TGraph g(100, x, y);
    g.SetTitle(Form("Waveform Time Demo (entry=%lld, ch=%d, module=%d, trigger=%s)",
                    entry, ch, module_id, trigger_sel.c_str()));
    g.GetXaxis()->SetTitle("time [ns]");
    g.GetYaxis()->SetTitle("ADC");
    g.SetLineColor(kBlack);
    g.SetMarkerStyle(20);
    g.SetMarkerColor(kBlack);
    g.SetMarkerSize(0.7);
    g.GetXaxis()->SetLimits(100.0, 250.0);
    g.GetYaxis()->SetRangeUser(ylow, yhigh);
    g.Draw("AP");

    TLegend leg(0.56, 0.64, 0.89, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&g, "waveform", "p");
    c.Update();

    TLine l_quad;
    TLine l_cfd;
    TLine l_logn;

    if (t_quad.ok) {
        l_quad.SetLineColor(kBlue + 1);
        l_quad.SetLineWidth(2);
        l_quad.SetLineStyle(2);
        l_quad.DrawLine(t_quad.time_ns, ylow, t_quad.time_ns, yhigh);
        leg.AddEntry(&l_quad, Form("quadratic: %.3f ns", t_quad.time_ns), "l");
    }
    if (t_cfd.ok) {
        l_cfd.SetLineColor(kRed + 1);
        l_cfd.SetLineWidth(2);
        l_cfd.SetLineStyle(1);
        l_cfd.DrawLine(t_cfd.time_ns, ylow, t_cfd.time_ns, yhigh);
        leg.AddEntry(&l_cfd,
                     Form("linear CFD: %.3f ns%s", t_cfd.time_ns,
                          t_cfd.used_fallback ? " (fallback)" : ""),
                     "l");
    }
    if (t_logn.ok) {
        std::vector<double> fit_x(fit_n);
        std::vector<double> fit_y(fit_n);
        for (int i = 0; i < fit_n; ++i) {
            const int sample_idx = fit_start + i;
            fit_x[i] = x[sample_idx];
            fit_y[i] = evalLogNormalFitAtSample(static_cast<float>(sample_idx), ped_mean, t_logn);
        }
        auto *fit_graph = new TGraph(fit_n, fit_x.data(), fit_y.data());
        fit_graph->SetLineColor(kGreen + 2);
        fit_graph->SetLineWidth(3);
        fit_graph->SetLineStyle(7);
        fit_graph->Draw("Lsame");

        l_logn.SetLineColor(kGreen + 2);
        l_logn.SetLineWidth(2);
        l_logn.SetLineStyle(7);
        l_logn.DrawLine(t_logn.time_ns, ylow, t_logn.time_ns, yhigh);
        leg.AddEntry(&l_logn,
                     Form("log-normal CFD: %.3f ns", t_logn.time_ns),
                     "l");
    } else {
        leg.AddEntry((TObject *)nullptr,
                     Form("log-normal CFD: %s",
                          t_logn.fail_reason.empty() ? "fit failed" : t_logn.fail_reason.c_str()),
                     "");
    }

    leg.Draw();
    c.Modified();
    c.Update();
    c.Print(pdf_path.c_str());
    std::cout << "[time-demo] appended plot to " << pdf_path << std::endl;
    return true;
}