#include "../EventData.h"
#include "../PhysicsTools.h"

namespace fs = std::filesystem;

static std::vector<std::string> collectRootFiles(const std::string &path);

void background(){

    float resolution = 0.026 * 1.5;
    float Ebeam = 3488.43; // MeV

    TFile *f = TFile::Open("../data/empty_target/prad_24386.filtered.root");
    TTree *tree = (TTree*)f->Get("recon");
    ReconEventData ev;
    setupReconBranches(tree, ev);

    // output histograms
        //overall distribution
    TH2F *hit_all = new TH2F("hit_all", "Hit Position Distribution(all clusters);X (mm);Y (mm)", 600, -360, 360, 600, -360, 360);
    TH2F *E_angle = new TH2F("E_angle", "E vs Scattering Angle;#Scattering Angle (deg);Energy (MeV)", 120, 0, 6, 5000, 0, 5000);
        
        //event selection
    TH2F *hits_mott = new TH2F("hits_mott", "Hit Position Distribution (e-p);X (mm);Y (mm)", 600, -360, 360, 600, -360, 360);
    TH2F *hits_moller = new TH2F("hits_moller", "Hit Position Distribution (2 arm Moller);X (mm);Y (mm)", 600, -360, 360, 600, -360, 360);
    TH2F *E_angle_mott = new TH2F("E_angle_mott", "E vs Scattering Angle (e-p);#Scattering Angle (deg);Energy (MeV)", 120, 0, 6, 5000, 0, 5000);
    TH2F *E_angle_moller = new TH2F("E_angle_moller", "E vs Scattering Angle (2 arm Moller);#Scattering Angle (deg);Energy (MeV)", 120, 0, 6, 5000, 0, 5000);
    
        //check matching quality
    TH1D *match_deltaX = new TH1D("match_deltaX", "Delta X of Matched Hits;#DeltaX (mm);Energy (MeV)", 300, -15, 15);
    TH1D *match_deltaY = new TH1D("match_deltaY", "Delta Y of Matched Hits;#DeltaY (mm);Energy (MeV)", 300, -15, 15);
    
        //mott and moller yield
    TH1F *mott_yield = new TH1F("e-p_yield", "e-p Yield;Theta (deg);Yield(arbitrary units)", 60, 0, 6);
    TH1F *moller_yield = new TH1F("moller_yield", "2 arm Moller Yield;Theta (deg);Yield(arbitrary units)", 60, 0, 6);
    TH1F *yield_ratio = new TH1F("yield_ratio", "e-p/Moller Yield Ratio;Theta (deg);Yield Ratio", 60, 0, 6);

    MollerData HC_moller_events;

    // Loop over events and fill histograms
    //two options:
    //1. no matching, only use hycal data
    //2. apply matching, and check matching quality, only fill the matching clusters in hist
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        cout << "Processing event " << i << " / " << tree->GetEntries() << "\r" << flush;
        tree->GetEntry(i);
        std::vector<HCHit> moller_Hits_candidate; // store potential Moller hits for this event
        for (int j = 0; j < ev.n_clusters; j++) {
            float x = ev.cl_x[j];
            float y = ev.cl_y[j];
            float z = ev.cl_z[j];
            float E = ev.cl_energy[j];
            float theta = atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;

            if( fabs(x) < 20.25 * 2.5 && fabs(y) < 20.25 * 2.5 ) continue;
            if( fabs(x) > 20.25 * 16. || fabs(y) > 20.25 * 16. ) continue;

            //if(theta < 0.6 || theta > 3.0) continue;

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
        // Sort by energy descending (highest first)
        std::sort(moller_Hits_candidate.begin(), moller_Hits_candidate.end(),
                  [](const HCHit &a, const HCHit &b){ return a.energy > b.energy; });

        // Find all Moller pairs: outer loop from highest energy, inner from lowest
        MollerData mollerData_event;
        int nCand = moller_Hits_candidate.size();
        for (int ii = 0; ii < nCand; ++ii) {
            const HCHit &hi = moller_Hits_candidate[ii];
            float theta_i = atan2(std::sqrt(hi.x*hi.x + hi.y*hi.y), hi.z) * 180.f / M_PI;
            for (int jj = nCand - 1; jj > ii; --jj) {
                const HCHit &hj = moller_Hits_candidate[jj];
                float theta_j = atan2(std::sqrt(hj.x*hj.x + hj.y*hj.y), hj.z) * 180.f / M_PI;
                if (isMoller_kinematic(theta_i, hi.energy, theta_j, hj.energy, Ebeam, resolution)) { 
                    // Require back-to-back in phi within 10 degrees
                    MollerEvent candidate = {DataPoint(hi.x, hi.y, hi.z, hi.energy), DataPoint(hj.x, hj.y, hj.z, hj.energy)};
                    if(fabs(GetMollerPhiDiff(candidate)) > 10.f) continue;
                    mollerData_event.emplace_back(candidate);
                }
            }
        }

        if(mollerData_event.size() == 0) continue;

        if(mollerData_event.size() > 1) {
            // Keep the pair with smallest total transverse momentum magnitude
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

        // Fill Moller histograms from found pairs
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
    cout << "Event processing completed. Total Moller candidates found: " << HC_moller_events.size() << endl;

    for(int i = 1; i < moller_yield->GetNbinsX(); i++){
        double mott = mott_yield->GetBinContent(i);
        double moller = moller_yield->GetBinContent(i);
        if(moller > 0) yield_ratio->SetBinContent(i, mott / moller);
        else yield_ratio->SetBinContent(i, 0);
    }

    //Draw histograms
    TCanvas *c = new TCanvas("c", "Background Analysis", 1200, 800);
    c->Divide(3, 2);
    c->cd(1); hit_all->Draw("COLZ");
    c->cd(2); E_angle->Draw("COLZ");
    c->cd(3); hits_mott->Draw("COLZ");
    c->cd(4); E_angle_mott->Draw("COLZ");
    c->cd(5); hits_moller->Draw("COLZ");
    c->cd(6); E_angle_moller->Draw("COLZ");
    c->SaveAs("background_overview.png");

    TCanvas *c2 = new TCanvas("c2", "e-p and Moller Yield", 1200, 800);
    c2->Divide(3, 1);
    c2->cd(1); mott_yield->Draw("HIST");
    c2->cd(2); moller_yield->Draw("HIST");
    c2->cd(3); yield_ratio->Draw("HIST");

    c2->SaveAs("background_yield.png");
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