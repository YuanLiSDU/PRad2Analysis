// summary_pdf.C
// Read histograms from four background_output.root files and save them
// to a single multi-page PDF.  One page per file, 3x4 pad layout.
//
// 2D hists        : COLZ + logz
// E vs angle 2D  : x=[0,4.5 deg], y=[0,1000 MeV], logz
// mott_yield      : logy
// moller_yield    : linear
//
// Usage:
//   root -l 'emptyTarget/summary_pdf.C'
//   root -l 'emptyTarget/summary_pdf.C("A.root","B.root","C.root","D.root","out.pdf")'

void summary_pdf(
    const char *fileA  = "../data/0.7GeV/typeA_24653.root",
    const char *fileB  = "../data/0.7GeV/typeB.root",
    const char *fileC  = "../data/0.7GeV/typeC.root",
    const char *fileD  = "../data/0.7GeV/typeD.root",
    const char *outpdf = "summary.pdf")
{
    const char *paths [4] = { fileA, fileB, fileC, fileD };
    const char *labels[4] = { "typeA", "typeB", "typeC", "typeD" };

    // ── Histogram configuration ───────────────────────────────────────────
    struct HCfg {
        const char *name;
        bool is2D;
        bool logy;
        bool eAngle;
        float xlo, xhi;
        const char *drawOpt; // 1D draw option, e.g. "E" or "HIST"
    };

    const HCfg cfg[] = {
        { "hit_all",         true,  false, false, 0, 0,          ""     },
        { "E_angle",         true,  false, true,  0, 0,          ""     },
        { "hits_mott",       true,  false, false, 0, 0,          ""     },
        { "E_angle_mott",    true,  false, true,  0, 0,          ""     },
        { "hits_moller",     true,  false, false, 0, 0,          ""     },
        { "E_angle_moller",  true,  false, true,  0, 0,          ""     },
        { "mott_yield",      false, true,  false, 0, 0,          "E"    },
        { "moller_yield",    false, false, false, 0, 0,          "E"    },
        { "yield_ratio",     false, false, false, 1.6f, 3.0f,   "HIST" },
        { "moller_center",   false, false, false, 0, 0,          "E"    },
        { "moller_center_y", false, false, false, 0, 0,          "E"    },
        { "moller_vertex_z", false, false, false, -1500.f, 6000.f, "E" },
    };

    // ── Page groups ───────────────────────────────────────────────────────
    // Each group: hist indices (up to 3), n hists, canvas W, canvas H, png tag
    const int   gIdx [5][3] = { {0,1,-1}, {2,3,-1}, {4,5,-1}, {6,7,8}, {9,10,11} };
    const int   gN   [5]    = { 2, 2, 2, 3, 3 };
    const int   gW   [5]    = { 1600, 1600, 1600, 1500, 1500 };
    const int   gH   [5]    = { 700,  700,  700,  500,  500  };
    const char *gTag [5]    = { "hit_E_angle", "mott_hit_E_angle", "moller_hit_E_angle",
                                "yield", "moller_center_vertex" };
    const int nGroups = 5;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Helper: draw one histogram on the current pad and add a title label
    auto drawHist = [&](TObject *obj, const HCfg &h, const char *lbl) {
        gPad->SetLogz(0); gPad->SetLogy(0);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(cfg->is2D ? 0.15 : 0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.12);
        if (h.is2D) {
            TH2 *hh = (TH2*)obj;
            if (h.eAngle) { hh->GetXaxis()->SetRangeUser(0.,4.5); hh->GetYaxis()->SetRangeUser(0.,1000.); }
            gPad->SetLogz();
            hh->Draw("COLZ");
        } else {
            TH1 *hh = (TH1*)obj;
            hh->SetLineColor(kBlue); hh->SetLineWidth(2);
            if (h.logy) gPad->SetLogy();
            if (h.xlo < h.xhi) hh->GetXaxis()->SetRangeUser(h.xlo, h.xhi);
            hh->Draw(h.drawOpt);
            if (hh->GetListOfFunctions()->GetSize() > 0) hh->GetListOfFunctions()->Draw("SAME");
        }
        TString title = h.is2D ? ((TH2*)obj)->GetTitle() : ((TH1*)obj)->GetTitle();
        if (title.Contains(";")) title = title(0, title.Index(";"));
        TLatex tex; tex.SetNDC(); tex.SetTextAlign(23); tex.SetTextSize(0.055); tex.SetTextFont(42);
        tex.DrawLatex(0.5, 0.955, Form("[%s]  %s", lbl, title.Data()));
    };

    // Use a resizable canvas — recreate per group size
    TCanvas *cv = nullptr;
    cv = new TCanvas("cv_summary", "Summary", gW[0], gH[0]);
    cv->Print(Form("%s[", outpdf));

    for (int iF = 0; iF < 4; iF++) {
        TFile *f = TFile::Open(paths[iF]);
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << paths[iF] << std::endl; continue; }
        std::cout << "Processing " << paths[iF] << std::endl;

        for (int ig = 0; ig < nGroups; ig++) {
            cv->SetCanvasSize(gW[ig], gH[ig]);
            cv->Clear();
            cv->Divide(gN[ig], 1, 0.003, 0.003);

            for (int ip = 0; ip < gN[ig]; ip++) {
                int ih = gIdx[ig][ip];
                cv->cd(ip + 1);
                TObject *obj = f->Get(cfg[ih].name);
                if (!obj) { std::cerr << "  Not found: " << cfg[ih].name << std::endl; continue; }
                drawHist(obj, cfg[ih], labels[iF]);
            }

            cv->Print(outpdf);
            cv->SaveAs(Form("../data/0.7GeV/%s_%s.png", labels[iF], gTag[ig]));
        }

        f->Close();
    }

    cv->Print(Form("%s]", outpdf));
    std::cout << "Saved to " << outpdf << std::endl;
}
