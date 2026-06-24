// compare_flow_hycal.C
// Compare HyCal-only e-p Mott signal ratios for 800flow and 1300flow runs.
//
// Signals:
//   800flow   = 24685_hycal - 24683_hycal
//   1300flow1 = 24941_hycal - 24942_hycal
//   1300flow2 = 24915_hycal - 24910_hycal
//
// Draw:
//   800flow / 1300flow1
//   800flow / 1300flow2
//
// Usage:
//   root -l 'emptyTarget/compare_flow_hycal.C'
//   root -l 'emptyTarget/compare_flow_hycal.C("emptyTarget")'

#include <algorithm>
#include <cmath>
#include <iostream>

static TString hycalFileName(const char *folder, int run)
{
    TString path = Form("%s/%d_hycal.root", folder, run);
    if (!gSystem->AccessPathName(path)) return path;

    TString local = Form("%d_hycal.root", run);
    if (!gSystem->AccessPathName(local)) return local;

    return path;
}

static TH1F *loadHycalYield(int run, const char *folder, const char *tag, const char *histName)
{
    TString path = hycalFileName(folder, run);
    TFile *f = TFile::Open(path);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << std::endl;
        return nullptr;
    }

    TH1F *h = (TH1F*)f->Get(histName);
    if (!h) {
        std::cerr << histName << " not found in " << path << std::endl;
        f->Close();
        return nullptr;
    }

    TH1F *clone = (TH1F*)h->Clone(Form("%s_%s_run%d", histName, tag, run));
    clone->SetDirectory(nullptr);
    if (clone->GetSumw2N() == 0) clone->Sumw2();
    f->Close();
    std::cout << "Loaded " << path << std::endl;
    return clone;
}

static TH1F *loadHycalMottYield(int run, const char *folder, const char *tag)
{
    return loadHycalYield(run, folder, tag, "mott_yield");
}

static TH1F *loadHycalMollerYield(int run, const char *folder, const char *tag)
{
    return loadHycalYield(run, folder, tag, "moller_yield");
}

static double hycalAExtraScale(int runA)
{
    if (runA == 24941) return 105838686. / 30000000.;
    if (runA == 24915) return 84735719. / 30000000.;
    return 1.0;
}

static TH1F *makeHycalMottSignalAmB(int runA, int runB, const char *folder, const char *tag)
{
    TH1F *hA = loadHycalMottYield(runA, folder, Form("%s_A", tag));
    TH1F *hB = loadHycalMottYield(runB, folder, Form("%s_B", tag));
    if (!hA || !hB) return nullptr;

    /*const double a_scale = hycalAExtraScale(runA);
    if (a_scale != 1.0) {
        hA->Scale(a_scale);
        std::cout << "Scaled A run " << runA << " by " << a_scale << std::endl;
    }*/

    TH1F *signal = (TH1F*)hA->Clone(Form("hycal_mott_signal_AmB_%s", tag));
    signal->SetDirectory(nullptr);
    signal->Add(hB, -1.);
    return signal;
}

static TH1F *makeHycalMollerSignalAmB(int runA, int runB, const char *folder, const char *tag)
{
    TH1F *hA = loadHycalMollerYield(runA, folder, Form("%s_A", tag));
    TH1F *hB = loadHycalMollerYield(runB, folder, Form("%s_B", tag));
    if (!hA || !hB) return nullptr;

    /*const double a_scale = hycalAExtraScale(runA);
    if (a_scale != 1.0) {
        hA->Scale(a_scale);
        std::cout << "Scaled A run " << runA << " by " << a_scale << std::endl;
    }*/

    TH1F *signal = (TH1F*)hA->Clone(Form("hycal_moller_signal_AmB_%s", tag));
    signal->SetDirectory(nullptr);
    signal->Add(hB, -1.);
    return signal;
}

static TH1F *makeHycalEpOverEeAmB(int runA, int runB, const char *folder, const char *tag)
{
    TH1F *mott = makeHycalMottSignalAmB(runA, runB, folder, Form("%s_ep", tag));
    TH1F *moller = makeHycalMollerSignalAmB(runA, runB, folder, Form("%s_ee", tag));
    if (!mott || !moller) return nullptr;

    mott->Rebin(4);
    moller->Rebin(4);

    TH1F *ratio = (TH1F*)mott->Clone(Form("hycal_ep_over_ee_AmB_%s", tag));
    ratio->SetDirectory(nullptr);
    ratio->Divide(mott, moller, 1., 1., "");
    return ratio;
}

static void drawTwoFlowRatios(TH1F *ratio_1, TH1F *ratio_2,
                              const char *title, const char *yTitle,
                              const char *output)
{
    const double thetaMin = 0.6;
    const double thetaMax = 3.8;
    TH1F *ratios[2] = {ratio_1, ratio_2};
    const int colors[2] = {kBlue + 1, kRed + 1};
    const int markers[2] = {20, 21};
    const char *labels[2] = {"800flow / 1300flow1", "800flow / 1300flow2"};

    double ymin = 1e9, ymax = -1e9;
    for (int i = 0; i < 2; i++) {
        for (int b = 1; b <= ratios[i]->GetNbinsX(); b++) {
            const double x = ratios[i]->GetBinCenter(b);
            if (x < thetaMin || x > thetaMax) continue;
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

    TCanvas *c = new TCanvas(Form("c_%s", output), title, 1000, 700);
    c->SetGrid();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    TF1 *fits[2] = {nullptr, nullptr};
    for (int i = 0; i < 2; i++) {
        fits[i] = new TF1(Form("fit_%s_%d", output, i), "pol1", thetaMin, thetaMax);
        fits[i]->SetLineColor(colors[i]);
        fits[i]->SetLineStyle(2);
        fits[i]->SetLineWidth(3);
        ratios[i]->Fit(fits[i], "RQ0");
    }

    for (int i = 0; i < 2; i++) {
        ratios[i]->SetStats(0);
        ratios[i]->SetLineColor(colors[i]);
        ratios[i]->SetMarkerColor(colors[i]);
        ratios[i]->SetMarkerStyle(markers[i]);
        ratios[i]->SetLineWidth(2);
        ratios[i]->GetXaxis()->SetRangeUser(thetaMin, thetaMax);
        if (i == 0) {
            ratios[i]->SetTitle(Form("%s;Reconstructed Scattering Angle [deg];%s", title, yTitle));
            ratios[i]->GetYaxis()->SetRangeUser(ymin, ymax);
            ratios[i]->GetXaxis()->CenterTitle();
            ratios[i]->GetYaxis()->CenterTitle();
            ratios[i]->Draw("E");
        } else {
            ratios[i]->Draw("E SAME");
        }
    }

    for (int i = 0; i < 2; i++) {
        fits[i]->Draw("SAME");
    }

    TLine *line1 = new TLine(thetaMin, 1., thetaMax, 1.);
    line1->SetLineColor(kGray + 1);
    line1->SetLineStyle(2);
    line1->Draw("SAME");

    TLegend *leg = new TLegend(0.58, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    for (int i = 0; i < 2; i++) {
        leg->AddEntry(ratios[i], labels[i], "PE");
        leg->AddEntry(fits[i],
                      Form("fit: %.3f %+.3f#theta",
                           fits[i]->GetParameter(0),
                           fits[i]->GetParameter(1)),
                      "l");
    }
    leg->Draw();

    TLatex tex;
    tex.SetNDC(1);
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.14, 0.84, "800flow: 24685_hycal - 24683_hycal");
    tex.DrawLatex(0.14, 0.79, "1300flow1: 24941_hycal - 24942_hycal");
    tex.DrawLatex(0.14, 0.74, "1300flow2: 24915_hycal - 24910_hycal");

    c->SaveAs(output);
}

void compare_flow_hycal(const char *folder = "emptyTarget")
{
    gStyle->SetOptStat(0);

    TH1F *sig_800 = makeHycalMottSignalAmB(24685, 24683, folder, "800flow");
    TH1F *sig_1300_1 = makeHycalMottSignalAmB(24941, 24942, folder, "1300flow1");
    TH1F *sig_1300_2 = makeHycalMottSignalAmB(24915, 24910, folder, "1300flow2");
    if (!sig_800 || !sig_1300_1 || !sig_1300_2) return;

    sig_800->Rebin(3);
    sig_1300_1->Rebin(3);
    sig_1300_2->Rebin(3);

    TH1F *ratio_1 = (TH1F*)sig_800->Clone("hycal_ratio_800flow_over_1300flow1");
    TH1F *ratio_2 = (TH1F*)sig_800->Clone("hycal_ratio_800flow_over_1300flow2");
    ratio_1->SetDirectory(nullptr);
    ratio_2->SetDirectory(nullptr);
    ratio_1->Divide(sig_800, sig_1300_1, 1., 1., "");
    ratio_2->Divide(sig_800, sig_1300_2, 1., 1., "");

    drawTwoFlowRatios(ratio_1, ratio_2,
                      "HyCal Mott Yield Signal Ratio",
                      "(A-B) ratio",
                      "hycal_mott_signal_flow_ratio.png");

    TH1F *ep_over_ee_800 = makeHycalEpOverEeAmB(24685, 24683, folder, "800flow");
    TH1F *ep_over_ee_1300_1 = makeHycalEpOverEeAmB(24941, 24942, folder, "1300flow1");
    TH1F *ep_over_ee_1300_2 = makeHycalEpOverEeAmB(24915, 24910, folder, "1300flow2");
    if (!ep_over_ee_800 || !ep_over_ee_1300_1 || !ep_over_ee_1300_2) return;

    TH1F *ep_ee_flow_ratio_1 = (TH1F*)ep_over_ee_800->Clone("hycal_ep_ee_ratio_800_over_1300flow1");
    TH1F *ep_ee_flow_ratio_2 = (TH1F*)ep_over_ee_800->Clone("hycal_ep_ee_ratio_800_over_1300flow2");
    ep_ee_flow_ratio_1->SetDirectory(nullptr);
    ep_ee_flow_ratio_2->SetDirectory(nullptr);
    ep_ee_flow_ratio_1->Divide(ep_over_ee_800, ep_over_ee_1300_1, 1., 1., "");
    ep_ee_flow_ratio_2->Divide(ep_over_ee_800, ep_over_ee_1300_2, 1., 1., "");

    drawTwoFlowRatios(ep_ee_flow_ratio_1, ep_ee_flow_ratio_2,
                      "HyCal (e-p)/(e-e) Signal Ratio",
                      "[(A-B)_{e-p}/(A-B)_{e-e}] ratio",
                      "hycal_ep_ee_signal_flow_ratio.png");
}
