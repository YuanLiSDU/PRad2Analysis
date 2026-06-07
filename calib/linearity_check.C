
#include "../PhysicsTools.h"

// Build expected-energy + ±3σ TGraphs for one process ("ep" or "ee")
// sigma = 3.3% / sqrt(E[GeV])  =>  3σ = 0.099 * sqrt(1000 * E_MeV)
void MakeKinematicGraphs(float EBeam,
                         const std::string &type,
                         TGraph *&gCenter,
                         TGraph *&gUp,
                         TGraph *&gDn,
                         int nPts = 200)
{
    double thetaMin = 0.4, thetaMax = 7.6;
    gCenter = new TGraph(nPts);
    gUp     = new TGraph(nPts);
    gDn     = new TGraph(nPts);

    for (int i = 0; i < nPts; ++i) {
        double theta = thetaMin + (thetaMax - thetaMin) * i / (nPts - 1);
        double E     = ExpectedEnergy((float)theta, EBeam, type);
        if (E <= 0) { gCenter->SetPoint(i, theta, 0); gUp->SetPoint(i, theta, 0); gDn->SetPoint(i, theta, 0); continue; }
        double sig3  = 0.099 * std::sqrt(1000.0 * E);  // 3σ in MeV
        gCenter->SetPoint(i, theta, E);
        gUp    ->SetPoint(i, theta, E + sig3);
        gDn    ->SetPoint(i, theta, E - sig3);
    }
}

void StyleGraph(TGraph *g, int color, int lstyle, int lwidth = 2)
{
    g->SetLineColor(color);
    g->SetLineStyle(lstyle);
    g->SetLineWidth(lwidth);
}

void linearity_check()
{
    const float EBeam = 728.9f;  // MeV, run 24655

    // ── open files ──────────────────────────────────────────────────────────
    TFile *fLin    = TFile::Open("prad_024655_quick_check.root",          "READ");
    TFile *fNonLin = TFile::Open("prad_024655_quick_check_nonLinear.root","READ");
    if (!fLin || fLin->IsZombie() || !fNonLin || fNonLin->IsZombie()) {
        std::cerr << "Cannot open input files" << std::endl; return;
    }

    TH2F *hLin    = (TH2F*)fLin   ->Get("energy_plots/h2_energy_theta");
    TH2F *hNonLin = (TH2F*)fNonLin->Get("energy_plots/h2_energy_theta");
    hLin   ->SetDirectory(nullptr);
    hNonLin->SetDirectory(nullptr);
    fLin->Close(); fNonLin->Close();

    // merge 2 x-bins into 1  (0.05 deg/bin → 0.10 deg/bin)
    hLin   ->RebinX(2);
    hNonLin->RebinX(2);

    // ── kinematic curves ────────────────────────────────────────────────────
    // ep: black, ee: red (kRed+1)
    const int COL_EP = kBlack;
    const int COL_EE = kRed+1;

    TGraph *gEP_cen, *gEP_up, *gEP_dn;
    TGraph *gEE_cen, *gEE_up, *gEE_dn;
    MakeKinematicGraphs(EBeam, "ep", gEP_cen, gEP_up, gEP_dn);
    MakeKinematicGraphs(EBeam, "ee", gEE_cen, gEE_up, gEE_dn);

    // center lines: dashed (style 2); sigma bands: solid (style 1)
    StyleGraph(gEP_cen, COL_EP, 2, 2);
    StyleGraph(gEP_up,  COL_EP, 1, 2);
    StyleGraph(gEP_dn,  COL_EP, 1, 2);
    StyleGraph(gEE_cen, COL_EE, 2, 2);
    StyleGraph(gEE_up,  COL_EE, 1, 2);
    StyleGraph(gEE_dn,  COL_EE, 1, 2);

    // label positions: evaluate expected energy at chosen theta
    const double labelTheta_ep = 3.0;
    const double labelTheta_ee = 3.0;
    const double Eep_label = ExpectedEnergy((float)labelTheta_ep, EBeam, "ep");
    const double Eee_label = ExpectedEnergy((float)labelTheta_ee, EBeam, "ee");

    // ── canvas ──────────────────────────────────────────────────────────────
    TCanvas *c = new TCanvas("c_linearity","Linearity Check", 1600, 700);
    c->Divide(2, 1, 0.005, 0.01);

    auto DrawPad = [&](TVirtualPad *pad, TH2F *h, const char *title) {
        pad->cd();
        pad->SetLogz();
        pad->SetLeftMargin(0.12);
        pad->SetRightMargin(0.13);
        pad->SetBottomMargin(0.12);
        h->SetTitle(Form("%s;#theta (deg);Energy (MeV)", title));
        h->SetStats(0);
        h->GetXaxis()->SetRangeUser(0.0, 4.0);
        h->GetYaxis()->SetRangeUser(50, 1000);
        h->Draw("COLZ");
        gEP_cen->Draw("L same");
        gEP_up ->Draw("L same");
        gEP_dn ->Draw("L same");
        gEE_cen->Draw("L same");
        gEE_up ->Draw("L same");
        gEE_dn ->Draw("L same");

        // TLatex labels next to the curves
        TLatex *latEP = new TLatex(labelTheta_ep + 0.05, Eep_label + 25, "e-p");
        latEP->SetTextColor(COL_EP);
        latEP->SetTextSize(0.042);
        latEP->SetTextFont(42);
        latEP->Draw();

        TLatex *latEE = new TLatex(labelTheta_ee + 0.05, Eee_label + 12, "e-e");
        latEE->SetTextColor(COL_EE);
        latEE->SetTextSize(0.042);
        latEE->SetTextFont(42);
        latEE->Draw();
    };

    DrawPad(c->GetPad(1), hLin,    "0.7 GeV Reconstructed Energy vs. Angle");
    DrawPad(c->GetPad(2), hNonLin, "0.7 GeV (non-linearity correction)");

    c->SaveAs("linearity_check.pdf");
    c->SaveAs("linearity_check.png");
    std::cout << "Saved linearity_check.pdf / .png" << std::endl;
}