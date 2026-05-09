#include "../EventData.h"
#include "../PhysicsTools.h"

// ── Fill histograms for one input file ───────────────────────────────────
static void fillHists(const TString &fname,
                      TH2F *hit_all, TH2F *E_angle,
                      TH2F *hits_mott, TH2F *hits_moller,
                      TH2F *E_angle_mott, TH2F *E_angle_moller,
                      TH1F *mott_yield, TH1F *moller_yield, TH1F *yield_ratio,
                      float Ebeam)
{
    float resolution = 0.026 * 1.5;

    TChain *tree = new TChain("recon");
    tree->Add(fname);
    std::cout << "Processing: " << fname << " (" << tree->GetEntries() << " entries)" << std::endl;

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

            hit_all->Fill(x, y);
            E_angle->Fill(theta, E);

            if (isMott(E, Ebeam, resolution)) {
                hits_mott->Fill(x, y);
                E_angle_mott->Fill(theta, E);
                mott_yield->Fill(theta);
            }

            if (E > 150. && E < Ebeam - resolution * Ebeam / sqrt(Ebeam/1000.))
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
                float sin_t1 = sqrt(mev.first.x*mev.first.x + mev.first.y*mev.first.y) / mev.first.z;
                float sin_t2 = sqrt(mev.second.x*mev.second.x + mev.second.y*mev.second.y) / mev.second.z;
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
    std::cout << "  Total Moller candidates: " << HC_moller_events.size() << std::endl;

    for(int i = 1; i <= moller_yield->GetNbinsX(); i++){
        double mott   = mott_yield->GetBinContent(i);
        double moller = moller_yield->GetBinContent(i);
        yield_ratio->SetBinContent(i, (moller > 0) ? mott / moller : 0.);
    }

    delete tree;
}

// ── Main function ─────────────────────────────────────────────────────────
// Usage:
//   root -l background.C                                          (default file)
//   root -l 'background.C("file1.root")'
//   root -l 'background.C("file1.root", 1234)'
//   root -l 'background.C("file1.root", 1234, "file2.root", 5678)'
void background(const char *file1 = "../data/empty_target/prad_24386.filtered.root",
                double lc1 = 1.0,
                const char *file2 = "",
                double lc2 = 1.0)
{
    float Ebeam = 3488.43f; // MeV

    TString cfg_file[2]  = { file1, file2 };
    double  livecharge[2]= { lc1,   lc2   };
    int nCfg = (TString(file2) == "") ? 1 : 2;

    for (int c = 0; c < nCfg; c++)
        std::cout << "Config" << c << ": " << cfg_file[c] << "  livecharge=" << livecharge[c] << std::endl;

    const char *label[2] = {"type A", "type B"};

    // Allocate histograms for each config
    TH2F *hit_all      [2], *E_angle      [2];
    TH2F *hits_mott    [2], *hits_moller  [2];
    TH2F *E_angle_mott [2], *E_angle_moller[2];
    TH1F *mott_yield   [2], *moller_yield [2], *yield_ratio[2];

    for (int c = 0; c < nCfg; c++) {
        TString t = Form("_c%d", c);
        hit_all      [c] = new TH2F("hit_all"+t,       "Hit Position Distribution(all clusters);X (mm);Y (mm)",                   600, -360, 360, 600, -360, 360);
        E_angle      [c] = new TH2F("E_angle"+t,       "E vs Scattering Angle;Scattering Angle (deg);Energy (MeV)",               120, 0, 6, 5000, 0, 5000);
        hits_mott    [c] = new TH2F("hits_mott"+t,     "Hit Position Distribution (e-p);X (mm);Y (mm)",                          600, -360, 360, 600, -360, 360);
        hits_moller  [c] = new TH2F("hits_moller"+t,   "Hit Position Distribution (2 arm Moller);X (mm);Y (mm)",                 600, -360, 360, 600, -360, 360);
        E_angle_mott [c] = new TH2F("E_angle_mott"+t,  "E vs Scattering Angle (e-p);Scattering Angle (deg);Energy (MeV)",        120, 0, 6, 5000, 0, 5000);
        E_angle_moller[c]= new TH2F("E_angle_moller"+t,"E vs Scattering Angle (2 arm Moller);Scattering Angle (deg);Energy (MeV)",120, 0, 6, 5000, 0, 5000);
        mott_yield   [c] = new TH1F("mott_yield"+t,    "e-p Yield;Theta (deg);Yield(arbitrary units)",    60, 0, 6);
        moller_yield [c] = new TH1F("moller_yield"+t,  "2 arm Moller Yield;Theta (deg);Yield(arbitrary units)", 60, 0, 6);
        yield_ratio  [c] = new TH1F("yield_ratio"+t,   "e-p/Moller Yield Ratio;Theta (deg);Yield Ratio",  60, 0, 6);
        // prepend label to all histogram titles
        TString lbl = Form("[%s] ", label[c]);
        hit_all      [c]->SetTitle(lbl + hit_all      [c]->GetTitle());
        E_angle      [c]->SetTitle(lbl + E_angle      [c]->GetTitle());
        hits_mott    [c]->SetTitle(lbl + hits_mott    [c]->GetTitle());
        hits_moller  [c]->SetTitle(lbl + hits_moller  [c]->GetTitle());
        E_angle_mott [c]->SetTitle(lbl + E_angle_mott [c]->GetTitle());
        E_angle_moller[c]->SetTitle(lbl + E_angle_moller[c]->GetTitle());
        mott_yield   [c]->SetTitle(lbl + mott_yield   [c]->GetTitle());
        moller_yield [c]->SetTitle(lbl + moller_yield [c]->GetTitle());
        yield_ratio  [c]->SetTitle(lbl + yield_ratio  [c]->GetTitle());
    }

    // Fill histograms
    for (int c = 0; c < nCfg; c++) {
        fillHists(cfg_file[c], hit_all[c], E_angle[c], hits_mott[c], hits_moller[c],
                  E_angle_mott[c], E_angle_moller[c],
                  mott_yield[c], moller_yield[c], yield_ratio[c], Ebeam);
        // scale yields by 1/livecharge so y-axis is yield per unit charge
        if (livecharge[c] != 1.0) {
            mott_yield  [c]->Scale(1.0 / livecharge[c]);
            moller_yield[c]->Scale(1.0 / livecharge[c]);
            // recompute ratio from scaled yields
            for (int i = 1; i <= yield_ratio[c]->GetNbinsX(); i++) {
                double mott   = mott_yield  [c]->GetBinContent(i);
                double moller = moller_yield[c]->GetBinContent(i);
                yield_ratio[c]->SetBinContent(i, (moller > 0) ? mott / moller : 0.);
            }
            mott_yield  [c]->GetYaxis()->SetTitle("Yield / livecharge");
            moller_yield[c]->GetYaxis()->SetTitle("Yield / livecharge");
        }
    }

    // Draw: one canvas row per config
    for (int c = 0; c < nCfg; c++) {
        TString tag = Form("_c%d", c);
        TString lbl = (nCfg > 1) ? Form(" [%s]", label[c]) : "";

        TCanvas *cv1 = new TCanvas("c_overview"+tag, "Background Analysis"+lbl, 1200, 800);
        cv1->Divide(3, 2);
        cv1->cd(1); hit_all      [c]->Draw("COLZ");
        cv1->cd(2); E_angle      [c]->Draw("COLZ");
        cv1->cd(3); hits_mott    [c]->Draw("COLZ");
        cv1->cd(4); E_angle_mott [c]->Draw("COLZ");
        cv1->cd(5); hits_moller  [c]->Draw("COLZ");
        cv1->cd(6); E_angle_moller[c]->Draw("COLZ");
        cv1->SaveAs(Form("background_overview_c%d.png", c));

        TCanvas *cv2 = new TCanvas("c_yield"+tag, "e-p and Moller Yield"+lbl, 1200, 400);
        cv2->Divide(3, 1);
        cv2->cd(1); mott_yield  [c]->Draw("HIST");
        cv2->cd(2); moller_yield[c]->Draw("HIST");
        cv2->cd(3); yield_ratio [c]->Draw("HIST");
        cv2->SaveAs(Form("background_yield_c%d.png", c));
    }
}