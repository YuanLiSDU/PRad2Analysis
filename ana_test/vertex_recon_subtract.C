#include "../EventData.h"
#include "../PhysicsTools.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// Usage:
//   root -l 'ana_test/vertex_recon_subtract.C("A.root","B.root","C.root","D.root")'
//   root -l 'ana_test/vertex_recon_subtract.C("A.root","B.root","C.root","D.root",QA,QB,QC,QD)'
// Outputs:
//   vertex_recon_types.png       : A, B, C, D charge-normalized vertex spectra
//   vertex_recon_type_diffs.png  : A-B, B-C, C-D charge-normalized differences

float Ebeam_subtract = 2243.5f; // MeV
float resolution_subtract = 0.037f;

static TH1F *makeVertexHist(const char *name, const char *label)
{
    TH1F *h = new TH1F(
        name,
        Form("%s Vertex Z;z_{vertex} (mm);Counts / charge", label),
        200, -3500, 6500);
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

static bool fillVertexHist(const char *path, TH1F *h_vz)
{
    TFile *file = TFile::Open(path);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << path << std::endl;
        return false;
    }

    TTree *tree = (TTree*)file->Get("recon");
    if (!tree) {
        std::cerr << "Cannot find tree 'recon' in " << path << std::endl;
        file->Close();
        return false;
    }

    ReconEventData ev;
    setupReconBranches(tree, ev);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (i % 10000 == 0) {
            std::cout << "  " << path << ": " << i << " / " << nentries << "\r" << std::flush;
        }

        if (ev.n_clusters != 1) continue;
        if (ev.matchNum != 1) continue;
        if (!isMott(ev.cl_energy[0], Ebeam_subtract, resolution_subtract)) continue;
        if (fabs(ev.cl_x[0]) < 20.25 * 2.5 && fabs(ev.cl_y[0]) < 20.25 * 2.5) continue;
        if (fabs(ev.cl_x[0]) > 20.25 * 16. || fabs(ev.cl_y[0]) > 20.25 * 16.) continue;

        const float x1 = ev.mHit_gx[0][1];
        const float y1 = ev.mHit_gy[0][1];
        const float z1 = ev.mHit_gz[0][1];
        const float dx = ev.mHit_gx[0][0] - ev.mHit_gx[0][1];
        const float dy = ev.mHit_gy[0][0] - ev.mHit_gy[0][1];
        const float dz = ev.mHit_gz[0][0] - ev.mHit_gz[0][1];

        const float den = dx * dx + dy * dy;
        if (den < 1e-12f) continue;

        const float t = -(x1 * dx + y1 * dy) / den;
        const float vz = z1 + t * dz;
        h_vz->Fill(vz);
    }

    std::cout << std::endl;
    file->Close();
    return true;
}

static void setHistStyle(TH1F *h, Color_t color, Style_t marker)
{
    h->SetStats(0);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetLineWidth(2);
}

static void setOverlayRange(TH1F *h_frame, TH1F **hist, int nhist)
{
    double ymin = 0.;
    double ymax = 0.;

    for (int ih = 0; ih < nhist; ih++) {
        for (int bin = 1; bin <= hist[ih]->GetNbinsX(); bin++) {
            const double y = hist[ih]->GetBinContent(bin);
            const double e = hist[ih]->GetBinError(bin);
            ymin = std::min(ymin, y - e);
            ymax = std::max(ymax, y + e);
        }
    }

    const double span = ymax - ymin;
    if (span <= 0.) {
        ymin = -1.;
        ymax = 1.;
    } else {
        ymin -= 0.12 * span;
        ymax += 0.18 * span;
    }

    h_frame->GetYaxis()->SetRangeUser(ymin, ymax);
}

void vertex_recon_subtract(
    const char *fileA = "/lustre24/expphy/volatile/hallb/prad/prad2_replays/2.2GeV/A/prad_025292/prad_025292_filtered.root",
    const char *fileB = "/lustre24/expphy/volatile/hallb/prad/prad2_replays/2.2GeV/B/prad_025297/prad_025297_filtered.root",
    const char *fileC = "/lustre24/expphy/volatile/hallb/prad/prad2_replays/2.2GeV/C/prad_025291/prad_025291_filtered.root",
    const char *fileD = "/lustre24/expphy/volatile/hallb/prad/prad2_replays/2.2GeV/D/prad_025309/prad_025309_filtered.root",
    double chargeA = 201080.,
    double chargeB = 1062770.,
    double chargeC = 876277.,
    double chargeD = 786689.)
{
    const char *files[4] = { fileA, fileB, fileC, fileD };
    const char *labels[4] = { "A", "B", "C", "D" };
    double charges[4] = { chargeA, chargeB, chargeC, chargeD };
    Color_t colors[4] = { kBlue + 1, kRed + 1, kGreen + 2, kViolet + 1 };
    Style_t markers[4] = { 20, 21, 22, 23 };

    for (int i = 0; i < 4; i++) {
        if (charges[i] <= 0.) {
            std::cerr << "Beam charge for type " << labels[i] << " must be positive." << std::endl;
            return;
        }
    }

    TH1F *h_vz[4];
    for (int i = 0; i < 4; i++) {
        h_vz[i] = makeVertexHist(Form("h_vz_%s", labels[i]), labels[i]);
        if (!fillVertexHist(files[i], h_vz[i])) return;
        h_vz[i]->Scale(1. / charges[i]);
        setHistStyle(h_vz[i], colors[i], markers[i]);
    }

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_vertex_types", "Charge-normalized vertex comparison", 1000, 700);
    c->cd();
    h_vz[0]->SetTitle("Vertex Z by Type;z_{vertex} (mm);Counts / charge");
    setOverlayRange(h_vz[0], h_vz, 4);

    h_vz[0]->Draw("HIST");
    for (int i = 1; i < 4; i++) {
        h_vz[i]->Draw("HIST SAME");
    }

    TLegend *leg = new TLegend(0.62, 0.66, 0.88, 0.88);
    leg->SetBorderSize(0);
    for (int i = 0; i < 4; i++) {
        leg->AddEntry(h_vz[i], Form("Type %s / Q, Q = %.6g", labels[i], charges[i]), "l");
    }
    leg->Draw();

    c->SaveAs("vertex_recon_types.png");

    TH1F *h_diff[3];
    const char *diff_labels[3] = { "A - B", "B - C", "C - D" };
    Color_t diff_colors[3] = { kBlack, kOrange + 7, kCyan + 2 };
    Style_t diff_markers[3] = { 20, 21, 22 };

    for (int i = 0; i < 3; i++) {
        h_diff[i] = (TH1F*)h_vz[i]->Clone(Form("h_vz_%sminus%s", labels[i], labels[i + 1]));
        h_diff[i]->SetTitle("Vertex Z Type Differences;z_{vertex} (mm);#Delta Counts / charge");
        h_diff[i]->Add(h_vz[i + 1], -1.);
        setHistStyle(h_diff[i], diff_colors[i], diff_markers[i]);
    }
    h_diff[0]->SetLineWidth(3);

    TCanvas *c_diff = new TCanvas("c_vertex_type_diffs", "Charge-normalized vertex differences", 1000, 700);
    c_diff->cd();
    setOverlayRange(h_diff[0], h_diff, 3);

    h_diff[0]->Draw("E1");
    for (int i = 1; i < 3; i++) {
        h_diff[i]->Draw("E1 SAME");
    }

    TLine *zero = new TLine(h_diff[0]->GetXaxis()->GetXmin(), 0.,
                            h_diff[0]->GetXaxis()->GetXmax(), 0.);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray + 1);
    zero->Draw("SAME");

    TLegend *leg_diff = new TLegend(0.62, 0.70, 0.88, 0.88);
    leg_diff->SetBorderSize(0);
    for (int i = 0; i < 3; i++) {
        leg_diff->AddEntry(h_diff[i], diff_labels[i], "lep");
    }
    leg_diff->Draw();

    c_diff->SaveAs("vertex_recon_type_diffs.png");
}

void vertex_recon_subtract_old_ab(
    const char *fileA,
    double chargeA,
    const char *fileB,
    double chargeB,
    const char *labelA = "A",
    const char *labelB = "B")
{
    if (chargeA <= 0. || chargeB <= 0.) {
        std::cerr << "Beam charge must be positive." << std::endl;
        return;
    }

    TH1F *h_A = makeVertexHist("h_vz_A", labelA);
    TH1F *h_B = makeVertexHist("h_vz_B", labelB);

    if (!fillVertexHist(fileA, h_A)) return;
    if (!fillVertexHist(fileB, h_B)) return;

    h_A->Scale(1. / chargeA);
    h_B->Scale(1. / chargeB);

    TH1F *h_sub = (TH1F*)h_A->Clone("h_vz_AminusB");
    h_sub->SetTitle(Form("Vertex Z: %s / Q_{%s} - %s / Q_{%s};z_{vertex} (mm);Counts / charge",
                         labelA, labelA, labelB, labelB));
    h_sub->Add(h_B, -1.);

    setHistStyle(h_A, kBlue + 1, 20);
    setHistStyle(h_B, kRed + 1, 21);
    setHistStyle(h_sub, kBlack, 24);
    h_sub->SetLineWidth(3);

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_vertex_subtract", "Charge-normalized vertex subtraction", 1000, 700);
    c->cd();
    TH1F *h_range[3] = { h_A, h_B, h_sub };
    setOverlayRange(h_sub, h_range, 3);

    h_sub->Draw("E1");
    h_A->Draw("HIST SAME");
    h_B->Draw("HIST SAME");
    h_sub->Draw("E1 SAME");

    TLine *zero = new TLine(h_sub->GetXaxis()->GetXmin(), 0.,
                            h_sub->GetXaxis()->GetXmax(), 0.);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray + 1);
    zero->Draw("SAME");

    TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h_A, Form("%s / Q, Q = %.6g", labelA, chargeA), "l");
    leg->AddEntry(h_B, Form("%s / Q, Q = %.6g", labelB, chargeB), "l");
    leg->AddEntry(h_sub, Form("%s / Q - %s / Q", labelA, labelB), "lep");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC(1);
    tex->SetTextSize(0.032);
    tex->DrawLatex(0.14, 0.86, Form("%s: %s", labelA, fileA));
    tex->DrawLatex(0.14, 0.81, Form("%s: %s", labelB, fileB));

    c->SaveAs("vertex_recon_subtract.png");
}
