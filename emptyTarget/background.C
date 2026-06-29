#include "../EventData.h"
#include "../PhysicsTools.h"

// ── Fill histograms for one or more input files ───────────────────────────
static void fillHists(const std::vector<TString> &fnames,
                      TH2F *hit_all, TH2F *E_angle,
                      TH2F *hits_mott, TH2F *hits_moller,
                      TH2F *E_angle_mott, TH2F *E_angle_moller,
                      TH1F *mott_yield, TH1F *moller_yield, TH1F *yield_ratio,
                      TH1F *moller_center_x, TH1F *moller_center_y, TH1F *moller_vertex_z,
                      float Ebeam)
{
    float resolution = 0.037f;

    TChain *tree = new TChain("recon");
    for (const auto &fname : fnames)
        tree->Add(fname);
    std::cout << "Processing " << fnames.size() << " file(s), total entries: " << tree->GetEntries() << std::endl;

    ReconEventData ev;
    setupReconBranches(tree, ev);

    MollerData HC_moller_events;

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        if (i % 10000 == 0)
            std::cout << "  " << i << " / " << tree->GetEntries() << "\r" << std::flush;
        tree->GetEntry(i);
        std::vector<HCHit> moller_Hits_candidate;
        for (int j = 0; j < ev.n_clusters; j++) {
            float x = ev.cl_x[j];
            float y = ev.cl_y[j];
            float z = ev.cl_z[j];
            float E = ev.cl_energy[j];
            float theta = atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;

            if( fabs(x) < 20.25 * 2.5 && fabs(y) < 20.25 * 2.5 ) continue;
            if( fabs(x) > 20.25 * 16. || fabs(y) > 20.25 * 16. ) continue;

            if( ev.cl_nblocks[j] <= 3 ) continue; // remove single-block clusters which are mostly noise

            hit_all->Fill(x, y);
            E_angle->Fill(theta, E);

            if (isMott(E, Ebeam, resolution)) {
                hits_mott->Fill(x, y);
                E_angle_mott->Fill(theta, E);
                mott_yield->Fill(theta);
            }

            if (E > 80. && E < Ebeam - resolution * Ebeam / sqrt(Ebeam/1000.))
                moller_Hits_candidate.emplace_back(x, y, z, E);
        }

        std::sort(moller_Hits_candidate.begin(), moller_Hits_candidate.end(),
                  [](const HCHit &a, const HCHit &b){ return a.energy > b.energy; });

        MollerData mollerData_event;
        int nCand = moller_Hits_candidate.size();
        for (int ii = 0; ii < nCand; ++ii) {
            const HCHit &hi = moller_Hits_candidate[ii];
            float theta_i = atan2(std::sqrt(hi.x*hi.x + hi.y*hi.y), hi.z) * 180.f / M_PI;
            for (int jj = nCand - 1; jj > ii; --jj) {
                const HCHit &hj = moller_Hits_candidate[jj];
                float theta_j = atan2(std::sqrt(hj.x*hj.x + hj.y*hj.y), hj.z) * 180.f / M_PI;
                if (isMoller_kinematic(theta_i, hi.energy, theta_j, hj.energy, Ebeam, resolution)) {
                    MollerEvent candidate = {DataPoint(hi.x, hi.y, hi.z, hi.energy), DataPoint(hj.x, hj.y, hj.z, hj.energy)};
                    if(fabs(GetMollerPhiDiff(candidate)) > 10.f) continue;
                    mollerData_event.emplace_back(candidate);
                }
            }
        }

        if(mollerData_event.size() == 0) continue;

        if(mollerData_event.size() > 1) {
            auto getPt = [](const MollerEvent &mev) -> float {
                float sin_t1 = sqrt(mev.first.x*mev.first.x + mev.first.y*mev.first.y) / std::sqrt(mev.first.z*mev.first.z + mev.first.x*mev.first.x + mev.first.y*mev.first.y);
                float sin_t2 = sqrt(mev.second.x*mev.second.x + mev.second.y*mev.second.y) / std::sqrt(mev.second.z*mev.second.z + mev.second.x*mev.second.x + mev.second.y*mev.second.y);
                return fabs(mev.first.E * sin_t1 - mev.second.E * sin_t2);
            };
            auto best = std::min_element(mollerData_event.begin(), mollerData_event.end(),
                                         [&](const MollerEvent &a, const MollerEvent &b){ return getPt(a) < getPt(b); });
            MollerEvent bestPair = *best;
            mollerData_event.clear();
            mollerData_event.push_back(bestPair);
        }

        const MollerEvent &mev = mollerData_event.front();
        HC_moller_events.push_back(mev);
        float t1 = atan2(std::sqrt(mev.first.x*mev.first.x   + mev.first.y*mev.first.y),   mev.first.z)  * 180.f / M_PI;
        float t2 = atan2(std::sqrt(mev.second.x*mev.second.x + mev.second.y*mev.second.y), mev.second.z) * 180.f / M_PI;
        hits_moller->Fill(mev.first.x,  mev.first.y);
        hits_moller->Fill(mev.second.x, mev.second.y);
        E_angle_moller->Fill(t1, mev.first.E);
        E_angle_moller->Fill(t2, mev.second.E);
        moller_yield->Fill(t1);
        moller_yield->Fill(t2);
    }
    std::cout << std::endl;
    std::cout << "  Total double arm Mollers : " << HC_moller_events.size() << std::endl;

    for(int i = 1; i <= moller_yield->GetNbinsX(); i++){
        double mott   = mott_yield->GetBinContent(i);
        double moller = moller_yield->GetBinContent(i);
        yield_ratio->SetBinContent(i, (moller > 0) ? mott / moller : 0.);
    }
    for(int i = 1; i <= HC_moller_events.size(); i++){
        auto center = GetMollerCenter(HC_moller_events[i-1], HC_moller_events[i]);
        moller_center_x->Fill(center[0]);
        moller_center_y->Fill(center[1]);
        float x1 = HC_moller_events[i-1].first.x, y1 = HC_moller_events[i-1].first.y, z1 = HC_moller_events[i-1].first.z;
        float x2 = HC_moller_events[i-1].second.x, y2 = HC_moller_events[i-1].second.y, z2 = HC_moller_events[i-1].second.z;
        float vertex = sqrt( (Ebeam + M_ELECTRON) * sqrt(x1*x1 + y1*y1) * sqrt(x2*x2 + y2*y2) / (2. * M_ELECTRON) );
        moller_vertex_z->Fill(6270.f - vertex); // Subtract nominal target center to get vertex relative to target center
    }

    delete tree;
}

// ── Main function ─────────────────────────────────────────────────────────
// Usage:
//   root -l background.C                                          (default file)
//   root -l 'background.C("file1.root")'
//   root -l 'background.C("file1.root", 1234)'
//   root -l 'background.C("file1.root,file2.root,file3.root", 1234)'
void background(const char *files = "../data/empty_target/prad_24386.filtered.root",
                double lc = 1.0)
{
    float Ebeam = 2239.51f; // MeV

    // Parse comma-separated file list
    std::vector<TString> fileList;
    TString allFiles(files);
    TObjArray *tokens = allFiles.Tokenize(",");
    for (int i = 0; i < tokens->GetEntries(); i++) {
        TString f = ((TObjString *)tokens->At(i))->GetString();
        f.Strip(TString::kBoth);
        if (!f.IsNull()) fileList.push_back(f);
    }
    delete tokens;

    std::cout << "Files (" << fileList.size() << "):" << std::endl;
    for (const auto &f : fileList) std::cout << "  " << f << std::endl;
    std::cout << "livecharge = " << lc << std::endl;

    // Allocate histograms
    TH2F *hit_all       = new TH2F("hit_all",        "Hit Position Distribution(all clusters);X (mm);Y (mm)",                    600, -360, 360, 600, -360, 360);
    TH2F *E_angle       = new TH2F("E_angle",         "E vs Scattering Angle;Scattering Angle (deg);Energy (MeV)",                120, 0, 6, 5000, 0, 5000);
    TH2F *hits_mott     = new TH2F("hits_mott",       "Hit Position Distribution (e-p);X (mm);Y (mm)",                           600, -360, 360, 600, -360, 360);
    TH2F *hits_moller   = new TH2F("hits_moller",     "Hit Position Distribution (2 arm Moller);X (mm);Y (mm)",                  600, -360, 360, 600, -360, 360);
    TH2F *E_angle_mott  = new TH2F("E_angle_mott",    "E vs Scattering Angle (e-p);Scattering Angle (deg);Energy (MeV)",         120, 0, 6, 5000, 0, 5000);
    TH2F *E_angle_moller= new TH2F("E_angle_moller",  "E vs Scattering Angle (2 arm Moller);Scattering Angle (deg);Energy (MeV)",120, 0, 6, 5000, 0, 5000);
    TH1F *mott_yield    = new TH1F("mott_yield",      "e-p Yield;Theta (deg);Yield(arbitrary units)",     60*2, 0, 6);
    TH1F *moller_yield  = new TH1F("moller_yield",    "2 arm Moller Yield;Theta (deg);Yield(arbitrary units)",  60*2, 0, 6);
    TH1F *yield_ratio   = new TH1F("yield_ratio",     "e-p/Moller Yield Ratio;Theta (deg);Yield Ratio",   60*2, 0, 6);
    TH1F *moller_center_x = new TH1F("moller_center",   "Moller Center X ; X(mm);Counts", 100, -10, 10);
    TH1F *moller_center_y = new TH1F("moller_center_y", "Moller Center Y ; Y(mm);Counts", 100, -10, 10);
    TH1F *moller_vertex_z = new TH1F("moller_vertex_z", "Moller Vertex Z ; Z(mm);Counts", 250*10, -4000, 6000);

    // Fill histograms
    fillHists(fileList, hit_all, E_angle, hits_mott, hits_moller,
              E_angle_mott, E_angle_moller,
              mott_yield, moller_yield, yield_ratio,
              moller_center_x, moller_center_y, moller_vertex_z,
              Ebeam);

    // Scale yields by 1/livecharge so y-axis is yield per unit charge
    if (lc != 1.0) {
        mott_yield  ->Scale(1.0 / lc);
        moller_yield->Scale(1.0 / lc);
        // Scale() already propagates bin errors correctly (sigma = sqrt(N)/lc).
        // No manual SetBinError needed.

        // recompute ratio from scaled yields
        for (int i = 1; i <= yield_ratio->GetNbinsX(); i++) {
            double mott   = mott_yield  ->GetBinContent(i);
            double moller = moller_yield->GetBinContent(i);
            yield_ratio->SetBinContent(i, (moller > 0) ? mott / moller : 0.);
        }
        mott_yield  ->GetYaxis()->SetTitle("Yield / livecharge");
        moller_yield->GetYaxis()->SetTitle("Yield / livecharge");
    }

    // Draw: overview canvas
    TCanvas *cv1 = new TCanvas("c_overview", "Background Analysis", 1200, 800);
    cv1->Divide(3, 2);
    cv1->cd(1); hit_all      ->Draw("COLZ");
    cv1->cd(2); E_angle      ->Draw("COLZ");
    cv1->cd(3); hits_mott    ->Draw("COLZ");
    cv1->cd(4); E_angle_mott ->Draw("COLZ");
    cv1->cd(5); hits_moller  ->Draw("COLZ");
    cv1->cd(6); E_angle_moller->Draw("COLZ");
    cv1->SaveAs("background_overview.png");

    // Draw: yield histograms
    TCanvas *cv2 = new TCanvas("c_yield", "e-p and Moller Yield", 1200, 400);
    cv2->Divide(3, 1);
    cv2->cd(1); mott_yield  ->SetLineColor(kBlue); mott_yield  ->SetLineWidth(2); mott_yield  ->Draw("E");
    cv2->cd(2); moller_yield->SetLineColor(kBlue); moller_yield->SetLineWidth(2); moller_yield->Draw("E");
    cv2->cd(3); yield_ratio ->SetLineColor(kBlue); yield_ratio ->SetLineWidth(2); yield_ratio ->Draw("E");
    cv2->SaveAs("background_yield.png");

    // Draw: Moller center and vertex
    gStyle->SetOptFit(1111); // show chi2/ndf, prob, params and errors on plot
    TCanvas *cv3 = new TCanvas("c_moller_center_vertex", "Moller Center and Vertex", 1200, 400);
    moller_center_x->Fit("gaus", "r", "", -3, 3);
    moller_center_y->Fit("gaus", "r", "", -3, 3);
    moller_vertex_z->Fit("gaus", "r", "", -100, 100);
    cv3->Divide(3, 1);
    cv3->cd(1); moller_center_x->SetLineColor(kBlue); moller_center_x->SetLineWidth(2); moller_center_x->Draw("E"); gPad->Update();
    cv3->cd(2); moller_center_y->SetLineColor(kBlue); moller_center_y->SetLineWidth(2); moller_center_y->Draw("E"); gPad->Update();
    cv3->cd(3); moller_vertex_z->SetLineColor(kBlue); moller_vertex_z->SetLineWidth(2); moller_vertex_z->Draw("E"); gPad->Update();
    cv3->SaveAs("background_moller_center_vertex.png");

    // ── Save everything to a ROOT file ────────────────────────────────────
    TFile *fout = TFile::Open("background_output.root", "RECREATE");
    // Histograms: TF1 fit functions are stored in GetListOfFunctions()
    // and are written automatically together with each histogram.
    hit_all       ->Write();
    E_angle       ->Write();
    hits_mott     ->Write();
    hits_moller   ->Write();
    E_angle_mott  ->Write();
    E_angle_moller->Write();
    mott_yield    ->Write();
    moller_yield  ->Write();
    yield_ratio   ->Write();
    moller_center_x->Write();
    moller_center_y->Write();
    moller_vertex_z->Write();
    // Canvases: preserve TPaveStats (fit parameter boxes) as drawn
    cv1->Write();
    cv2->Write();
    cv3->Write();
    fout->Close();
    std::cout << "Saved to background_output.root" << std::endl;
}