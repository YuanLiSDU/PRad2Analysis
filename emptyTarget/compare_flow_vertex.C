// compare_flow_vertex.C
// Compare e-p background ratios for three A/B/C run groups, with and without
// vertex cut.  Top row: noVertex files.  Bottom row: vertex-cut files.
//
// For each pad, draw:
//   B/(A-B), C/(A-B), (B-C)/(A-B)
// and fit (B-C)/(A-B) with y = p0 + p1 * ln(theta) on a log-x axis.
//
// Usage:
//   root -l 'emptyTarget/compare_flow_vertex.C'
//   root -l 'emptyTarget/compare_flow_vertex.C("emptyTarget")'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

static TString flowVertexFileName(const char *folder, int run, bool noVertex)
{
    TString path = Form("%s/%d%s.root", folder, run, noVertex ? "_noVertex" : "");
    if (!gSystem->AccessPathName(path)) return path;

    TString local = Form("%d%s.root", run, noVertex ? "_noVertex" : "");
    if (!gSystem->AccessPathName(local)) return local;

    return path;
}

static TH1F *loadMottYield(int run, bool noVertex, const char *groupTag, const char *role, const char *folder)
{
    TString path = flowVertexFileName(folder, run, noVertex);
    TFile *f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << std::endl;
        return nullptr;
    }

    TH1F *h = (TH1F*)f->Get("mott_yield");
    if (!h) {
        std::cerr << "mott_yield not found in " << path << std::endl;
        f->Close();
        return nullptr;
    }

    TH1F *clone = (TH1F*)h->Clone(Form("mott_%s_%s_run%d", groupTag, role, run));
    clone->SetDirectory(nullptr);
    if (clone->GetSumw2N() == 0) clone->Sumw2();
    f->Close();
    std::cout << "Loaded " << path << std::endl;
    return clone;
}

static void fixBoverAmB(TH1F *hratio, TH1F *hA, TH1F *hB)
{
    for (int b = 1; b <= hratio->GetNbinsX(); b++) {
        double A = hA->GetBinContent(b), sA = hA->GetBinError(b);
        double B = hB->GetBinContent(b), sB = hB->GetBinError(b);
        double d = A - B;
        if (d == 0.) continue;
        hratio->SetBinError(b, std::sqrt(B * B * sA * sA + A * A * sB * sB) / (d * d));
    }
}

static void fixBmCoverAmB(TH1F *hratio, TH1F *hA, TH1F *hB, TH1F *hC)
{
    for (int b = 1; b <= hratio->GetNbinsX(); b++) {
        double A = hA->GetBinContent(b), sA = hA->GetBinError(b);
        double B = hB->GetBinContent(b), sB = hB->GetBinError(b);
        double C = hC->GetBinContent(b), sC = hC->GetBinError(b);
        double d = A - B;
        if (d == 0.) continue;
        hratio->SetBinError(b, std::sqrt((B - C) * (B - C) * sA * sA +
                                         (A - C) * (A - C) * sB * sB +
                                          d * d * sC * sC) / (d * d));
    }
}

static std::pair<double, double> flowYRange(TH1F **h, int n, double ymin_default = -0.2, double ymax_default = 2.0)
{
    double lo = 1e9, hi = -1e9;
    for (int i = 0; i < n; i++) {
        for (int b = 1; b <= h[i]->GetNbinsX(); b++) {
            double v = h[i]->GetBinContent(b);
            double e = h[i]->GetBinError(b);
            if (v == 0. && e == 0.) continue;
            lo = std::min(lo, v - e);
            hi = std::max(hi, v + e);
        }
    }
    if (lo > hi) return {ymin_default, ymax_default};

    double span = hi - lo;
    if (span <= 0.) span = 1.;
    lo -= 0.15 * span;
    hi += 0.20 * span;
    return {std::min(lo, ymin_default), std::max(hi, ymax_default)};
}

static bool makeRatiosForGroup(int groupIdx, int runA, int runB, int runC, bool noVertex, const char *folder,
                               TH1F *ratio[3], TF1 *&fit_bmc)
{
    const char *cut_tag = noVertex ? "noVertex" : "vertex";
    TString groupTag = Form("g%d_%s", groupIdx + 1, cut_tag);

    TH1F *hA = loadMottYield(runA, noVertex, groupTag, "A", folder);
    TH1F *hB = loadMottYield(runB, noVertex, groupTag, "B", folder);
    TH1F *hC = loadMottYield(runC, noVertex, groupTag, "C", folder);
    if (!hA || !hB || !hC) return false;

    double a_extra_scale = 1.0;
    if (runA == 24941) {
        a_extra_scale = 105838686. / 30000000.;
    } else if (runA == 24915) {
        a_extra_scale = 84735719. / 30000000.;
    }
    if (a_extra_scale != 1.0) {
        hA->Scale(a_extra_scale);
        std::cout << "Scaled A run " << runA << " by " << a_extra_scale << std::endl;
    }

    TH1F *denom = (TH1F*)hA->Clone(Form("denom_%s", groupTag.Data()));
    denom->SetDirectory(nullptr);
    denom->Add(hB, -1.);

    ratio[0] = (TH1F*)hB->Clone(Form("ratio_B_over_AmB_%s", groupTag.Data()));
    ratio[1] = (TH1F*)hC->Clone(Form("ratio_C_over_AmB_%s", groupTag.Data()));
    TH1F *bmc = (TH1F*)hB->Clone(Form("BminusC_%s", groupTag.Data()));
    bmc->SetDirectory(nullptr);
    bmc->Add(hC, -1.);
    ratio[2] = (TH1F*)bmc->Clone(Form("ratio_BmC_over_AmB_%s", groupTag.Data()));

    ratio[0]->Divide(hB, denom, 1., 1., "");
    ratio[1]->Divide(hC, denom, 1., 1., "");
    ratio[2]->Divide(bmc, denom, 1., 1., "");

    fixBoverAmB(ratio[0], hA, hB);
    fixBmCoverAmB(ratio[2], hA, hB, hC);

    for (int i = 0; i < 3; i++) {
        ratio[i]->SetDirectory(nullptr);
        ratio[i]->Scale(100.);
        ratio[i]->SetStats(0);
        ratio[i]->SetTitle(";Reconstructed Scattering Angle [deg];N_{background} / N_{signal} (%)");
    }

    fit_bmc = new TF1(Form("fit_BmC_over_AmB_%s", groupTag.Data()),
                      "[0]+[1]*log(x)", 0.5, 3.8);
    fit_bmc->SetParameters(0.2, 0.0);
    fit_bmc->SetLineColor(kBlack);
    fit_bmc->SetLineStyle(2);
    fit_bmc->SetLineWidth(2);
    ratio[2]->Fit(fit_bmc, "RQ0");

    return true;
}

static void drawFlowPad(int groupIdx, int runA, int runB, int runC, bool noVertex, const char *folder)
{
    TH1F *ratio[3] = {nullptr, nullptr, nullptr};
    TF1 *fit_bmc = nullptr;
    if (!makeRatiosForGroup(groupIdx, runA, runB, runC, noVertex, folder, ratio, fit_bmc)) return;

    const int colors[3] = {kRed, kGreen + 2, kBlue + 1};
    const int markers[3] = {20, 21, 24};
    const char *labels[3] = {"B/(A-B)", "C/(A-B)", "(B-C)/(A-B)"};

    gPad->SetLogx();
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(noVertex ? 0.08 : 0.13);
    gPad->SetTopMargin(0.09);

    for (int i = 0; i < 3; i++) {
        ratio[i]->SetLineColor(colors[i]);
        ratio[i]->SetMarkerColor(colors[i]);
        ratio[i]->SetMarkerStyle(markers[i]);
        ratio[i]->SetLineWidth(2);
        ratio[i]->GetXaxis()->SetRangeUser(0.5, 3.8);
        if (i == 0) {
            ratio[i]->GetYaxis()->SetRangeUser(-0.2, 2.0);
            ratio[i]->GetXaxis()->CenterTitle();
            ratio[i]->GetYaxis()->CenterTitle();
            ratio[i]->GetXaxis()->SetMoreLogLabels();
            ratio[i]->GetXaxis()->SetNoExponent();
            ratio[i]->Draw("E");
        } else {
            ratio[i]->Draw("E SAME");
        }
    }

    TLine *line0 = new TLine(0.5, 0., 3.8, 0.);
    line0->SetLineColor(kGray + 1);
    line0->SetLineStyle(3);
    line0->Draw("SAME");

    TLine *line1 = new TLine(0.5, 1., 3.8, 1.);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw("SAME");

    fit_bmc->Draw("SAME");

    TLatex tex;
    tex.SetNDC(1);
    tex.SetTextSize(0.050);
    tex.DrawLatex(0.16, 0.84, Form("%s  A:%d  B:%d  C:%d",
                                   noVertex ? "no vertex cut" : "vertex cut",
                                   runA, runB, runC));

    TLegend *leg = new TLegend(0.47, 0.62, 0.91, 0.91);
    leg->SetName(Form("leg_g%d_%s", groupIdx + 1, noVertex ? "noVertex" : "vertex"));
    leg->SetBorderSize(0);
    leg->SetTextSize(0.038);
    for (int i = 0; i < 3; i++) leg->AddEntry(ratio[i], labels[i], "PE");
    leg->AddEntry(fit_bmc,
                  Form("(B-C)/(A-B) fit: %.3f %+.3f ln#theta",
                       fit_bmc->GetParameter(0), fit_bmc->GetParameter(1)),
                  "l");
    leg->Draw();
}

void compare_flow_vertex(const char *folder = "emptyTarget")
{
    gStyle->SetOptStat(0);

    const int runs[3][3] = {
        {24685, 24683, 24874},
        {24941, 24942, 24943},
        {24915, 24910, 24943}
    };

    TCanvas *c = new TCanvas("c_background_ratio_flow_vertex",
                             "Background ratio flow with/without vertex cut",
                             1800, 1000);
    c->Divide(3, 2);

    for (int col = 0; col < 3; col++) {
        c->cd(col + 1);
        drawFlowPad(col, runs[col][0], runs[col][1], runs[col][2], true, folder);

        c->cd(col + 4);
        drawFlowPad(col, runs[col][0], runs[col][1], runs[col][2], false, folder);
    }

    c->SaveAs("compare_flow_vertex.png");
}

static TH1F *loadMollerVertexZ(int run, const char *folder)
{
    TString path = flowVertexFileName(folder, run, false);
    TFile *f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << std::endl;
        return nullptr;
    }

    TH1F *h = (TH1F*)f->Get("moller_vertex_z");
    if (!h) {
        std::cerr << "moller_vertex_z not found in " << path << std::endl;
        f->Close();
        return nullptr;
    }

    TH1F *clone = (TH1F*)h->Clone(Form("moller_vertex_z_run%d_scaled", run));
    clone->SetDirectory(nullptr);
    f->Close();
    return clone;
}

void plot_moller_vertex_z_A_runs(const char *folder = "emptyTarget")
{
    gStyle->SetOptStat(0);

    const int nRuns = 3;
    const int runs[nRuns] = {24685, 24941, 24915};
    const int colors[nRuns] = {kBlue + 1, kRed + 1, kGreen + 2};
    const int markers[nRuns] = {20, 21, 22};

    TH1F *h[nRuns] = {nullptr};
    double ymax = 0.;

    for (int i = 0; i < nRuns; i++) {
        h[i] = loadMollerVertexZ(runs[i], folder);
        if (!h[i]) return;

        const double peak = h[i]->GetMaximum();
        if (peak > 0.) h[i]->Scale(1. / peak);
        ymax = std::max(ymax, h[i]->GetMaximum());

        h[i]->SetStats(0);
        h[i]->SetLineColor(colors[i]);
        h[i]->SetMarkerColor(colors[i]);
        h[i]->SetMarkerStyle(markers[i]);
        h[i]->SetLineWidth(2);
    }

    TCanvas *c = new TCanvas("c_moller_vertex_z_A_runs",
                             "Moller Vertex Z for A Runs", 1000, 700);
    c->SetGrid();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    h[0]->SetTitle("Moller Vertex Z for A Runs;z_{vertex} (mm);Counts / peak");
    h[0]->GetYaxis()->SetRangeUser(0., ymax * 1.25);
    h[0]->GetXaxis()->CenterTitle();
    h[0]->GetYaxis()->CenterTitle();
    h[0]->GetXaxis()->SetTitleSize(0.045);
    h[0]->GetYaxis()->SetTitleSize(0.045);
    h[0]->Draw("HIST");

    for (int i = 1; i < nRuns; i++) {
        h[i]->Draw("HIST SAME");
    }

    TLegend *leg = new TLegend(0.64, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.040);
    for (int i = 0; i < nRuns; i++) {
        leg->AddEntry(h[i], Form("Run %d", runs[i]), "l");
    }
    leg->Draw();

    TLatex label;
    label.SetNDC(1);
    label.SetTextSize(0.040);
    label.DrawLatex(0.15, 0.84, "Each histogram scaled to peak = 1");

    c->SaveAs("moller_vertex_z_A_runs.png");
}

static double flowAExtraScale(int runA)
{
    if (runA == 24941) return 105838686. / 30000000.;
    if (runA == 24915) return 84735719. / 30000000.;
    return 1.0;
}

static TH1F *makeMottSignalAmB(int runA, int runB, bool noVertex,
                               const char *folder, const char *tag)
{
    TH1F *hA = loadMottYield(runA, noVertex, tag, "A", folder);
    TH1F *hB = loadMottYield(runB, noVertex, tag, "B", folder);
    if (!hA || !hB) return nullptr;

    const double a_scale = flowAExtraScale(runA);
    if (a_scale != 1.0) {
        hA->Scale(a_scale);
        std::cout << "Scaled A run " << runA << " by " << a_scale << std::endl;
    }

    TH1F *signal = (TH1F*)hA->Clone(Form("mott_signal_AmB_%s", tag));
    signal->SetDirectory(nullptr);
    signal->Add(hB, -1.);
    return signal;
}

void plot_mott_signal_flow_ratio(const char *folder = "emptyTarget", bool noVertex = false)
{
    gStyle->SetOptStat(0);

    TH1F *sig_800 = makeMottSignalAmB(24685, 24683, noVertex, folder,
                                      noVertex ? "800flow_noVertex" : "800flow_vertex");
    TH1F *sig_1300_1 = makeMottSignalAmB(24941, 24942, noVertex, folder,
                                         noVertex ? "1300flow1_noVertex" : "1300flow1_vertex");
    TH1F *sig_1300_2 = makeMottSignalAmB(24915, 24910, noVertex, folder,
                                         noVertex ? "1300flow2_noVertex" : "1300flow2_vertex");
    if (!sig_800 || !sig_1300_1 || !sig_1300_2) return;

    sig_800->Rebin(3);
    sig_1300_1->Rebin(3);
    sig_1300_2->Rebin(3);

    TH1F *ratio_1 = (TH1F*)sig_800->Clone("ratio_800flow_over_1300flow1");
    TH1F *ratio_2 = (TH1F*)sig_800->Clone("ratio_800flow_over_1300flow2");
    ratio_1->SetDirectory(nullptr);
    ratio_2->SetDirectory(nullptr);
    ratio_1->Divide(sig_800, sig_1300_1, 1., 1., "");
    ratio_2->Divide(sig_800, sig_1300_2, 1., 1., "");

    TH1F *ratios[2] = {ratio_1, ratio_2};
    const int colors[2] = {kBlue + 1, kRed + 1};
    const int markers[2] = {20, 21};
    const char *labels[2] = {"800flow / 1300flow1", "800flow / 1300flow2"};

    double ymin = 1e9, ymax = -1e9;
    for (int i = 0; i < 2; i++) {
        for (int b = 1; b <= ratios[i]->GetNbinsX(); b++) {
            double y = ratios[i]->GetBinContent(b);
            double e = ratios[i]->GetBinError(b);
            if (y == 0. && e == 0.) continue;
            ymin = std::min(ymin, y - e);
            ymax = std::max(ymax, y + e);
        }
    }
    if (ymin > ymax) {
        ymin = 0.;
        ymax = 2.;
    }
    double span = ymax - ymin;
    if (span <= 0.) span = 1.;
    ymin -= 0.15 * span;
    ymax += 0.20 * span;

    TCanvas *c = new TCanvas("c_mott_signal_flow_ratio",
                             "Mott Signal Flow Ratio", 1000, 700);
    c->SetGrid();
    c->SetLogx();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    for (int i = 0; i < 2; i++) {
        ratios[i]->SetStats(0);
        ratios[i]->SetLineColor(colors[i]);
        ratios[i]->SetMarkerColor(colors[i]);
        ratios[i]->SetMarkerStyle(markers[i]);
        ratios[i]->SetLineWidth(2);
        ratios[i]->GetXaxis()->SetRangeUser(0.5, 3.8);
        if (i == 0) {
            ratios[i]->SetTitle(Form("Mott Yield Signal Ratio (%s);Reconstructed Scattering Angle [deg];(A-B) ratio",
                                     noVertex ? "no vertex cut" : "vertex cut"));
            ratios[i]->GetYaxis()->SetRangeUser(ymin, ymax);
            ratios[i]->GetXaxis()->CenterTitle();
            ratios[i]->GetYaxis()->CenterTitle();
            ratios[i]->GetXaxis()->SetMoreLogLabels();
            ratios[i]->GetXaxis()->SetNoExponent();
            ratios[i]->Draw("E");
        } else {
            ratios[i]->Draw("E SAME");
        }
    }

    TLine *line1 = new TLine(0.5, 1., 3.8, 1.);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw("SAME");

    TLegend *leg = new TLegend(0.58, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    for (int i = 0; i < 2; i++) leg->AddEntry(ratios[i], labels[i], "PE");
    leg->Draw();

    TLatex tex;
    tex.SetNDC(1);
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.14, 0.84, "800flow: 24685 - 24683");
    tex.DrawLatex(0.14, 0.79, "1300flow1: 24941 - 24942");
    tex.DrawLatex(0.14, 0.74, "1300flow2: 24915 - 24910");

    c->SaveAs(noVertex ? "mott_signal_flow_ratio_noVertex.png" : "mott_signal_flow_ratio.png");
}
