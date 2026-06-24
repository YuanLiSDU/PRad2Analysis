#include "EventData.h"
#include "PhysicsTools.h"

#include <cmath>
#include <iostream>

float Ebeam = 3488.43f; // MeV
float resolution = 0.035f;

static bool getLineLineVertex(const float p1[3], const float p2[3],
                              const float q1[3], const float q2[3],
                              float &vx, float &vy, float &vz, float &dca)
{
    const float ux = p2[0] - p1[0];
    const float uy = p2[1] - p1[1];
    const float uz = p2[2] - p1[2];
    const float vx2 = q2[0] - q1[0];
    const float vy2 = q2[1] - q1[1];
    const float vz2 = q2[2] - q1[2];
    const float wx = p1[0] - q1[0];
    const float wy = p1[1] - q1[1];
    const float wz = p1[2] - q1[2];

    const float a = ux * ux + uy * uy + uz * uz;
    const float b = ux * vx2 + uy * vy2 + uz * vz2;
    const float c = vx2 * vx2 + vy2 * vy2 + vz2 * vz2;
    const float d = ux * wx + uy * wy + uz * wz;
    const float e = vx2 * wx + vy2 * wy + vz2 * wz;
    const float den = a * c - b * b;

    if (a < 1e-12f || c < 1e-12f || fabs(den) < 1e-12f) return false;

    const float s = (b * e - c * d) / den;
    const float t = (a * e - b * d) / den;

    const float cx1 = p1[0] + s * ux;
    const float cy1 = p1[1] + s * uy;
    const float cz1 = p1[2] + s * uz;
    const float cx2 = q1[0] + t * vx2;
    const float cy2 = q1[1] + t * vy2;
    const float cz2 = q1[2] + t * vz2;

    vx = 0.5f * (cx1 + cx2);
    vy = 0.5f * (cy1 + cy2);
    vz = 0.5f * (cz1 + cz2);
    dca = std::sqrt((cx1 - cx2) * (cx1 - cx2) +
                    (cy1 - cy2) * (cy1 - cy2) +
                    (cz1 - cz2) * (cz1 - cz2));
    return true;
}

void vertex_recon_moller()
{
    TFile *file = TFile::Open("../data/recon/3.5GeV/prad_024917_recon.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open input file" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("recon");
    if (!tree) {
        std::cerr << "Cannot find tree 'recon'" << std::endl;
        file->Close();
        return;
    }

    ReconEventData ev;
    setupReconBranches(tree, ev);

    TH1F *h_vz = new TH1F("h_moller_vz",
                          "Moller Vertex Z from Two GEM Track Lines;z_{vertex} (mm);Counts",
                          500, -3500, 6500);
    TH1F *h_dca = new TH1F("h_moller_track_dca",
                           "Closest Distance between Two Moller Track Lines;DCA (mm);Counts",
                           400, 0, 40);
    TH2F *h_vxy = new TH2F("h_moller_vxy",
                           "Moller Vertex X-Y from Two GEM Track Lines;x_{vertex} (mm);y_{vertex} (mm)",
                           300, -30, 30, 300, -30, 30);

    Long64_t n_moller = 0;
    Long64_t n_vertex = 0;

    for (Long64_t i = 0; i < tree->GetEntries()/100; i++) {
        tree->GetEntry(i);

        if(i%10000 == 0)
            std::cout << "  " << i << " / " << tree->GetEntries() << "\r" << std::flush;

        if(ev.matchNum != 2) continue;
        
        if( fabs(ev.cl_x[0]) < 20.25 * 2.5 && fabs(ev.cl_y[0]) < 20.25 * 2.5 ) continue;
        if( fabs(ev.cl_x[0]) > 20.25 * 16. || fabs(ev.cl_y[0]) > 20.25 * 16. ) continue;
        if( fabs(ev.cl_x[1]) < 20.25 * 2.5 && fabs(ev.cl_y[1]) < 20.25 * 2.5 ) continue;
        if( fabs(ev.cl_x[1]) > 20.25 * 16. || fabs(ev.cl_y[1]) > 20.25 * 16. ) continue;

        float theta1 = atan2(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]), ev.cl_z[0]) * 180.f / M_PI;
        float theta2 = atan2(std::sqrt(ev.cl_x[1]*ev.cl_x[1] + ev.cl_y[1]*ev.cl_y[1]), ev.cl_z[1]) * 180.f / M_PI;
        float x[2][4], y[2][4], z[2][4], E[2];

        x[0][0] = 0.f; y[0][0] = 0.f; z[0][0] = 0.f; // target center
        x[1][0] = 0.f; y[1][0] = 0.f; z[1][0] = 0.f;

        x[0][1] = ev.mHit_gx[0][1]; y[0][1] = ev.mHit_gy[0][1]; z[0][1] = ev.mHit_gz[0][1];
        x[0][2] = ev.mHit_gx[0][0]; y[0][2] = ev.mHit_gy[0][0]; z[0][2] = ev.mHit_gz[0][0];
        x[0][3] = ev.cl_x[0]; y[0][3] = ev.cl_y[0]; z[0][3] = ev.cl_z[0];

        x[1][1] = ev.mHit_gx[1][1]; y[1][1] = ev.mHit_gy[1][1]; z[1][1] = ev.mHit_gz[1][1];
        x[1][2] = ev.mHit_gx[1][0]; y[1][2] = ev.mHit_gy[1][0]; z[1][2] = ev.mHit_gz[1][0];
        x[1][3] = ev.cl_x[1]; y[1][3] = ev.cl_y[1]; z[1][3] = ev.cl_z[1];

        E[0] = ev.cl_energy[0];
        E[1] = ev.cl_energy[1];

        if(!isMoller_kinematic(theta1, E[0], theta2, E[1], Ebeam, resolution))
            continue;

        MollerEvent moller = {DataPoint(x[0][1], y[0][1], z[0][1], E[0]), DataPoint(x[1][1], y[1][1], z[1][1], E[1])};

        if(fabs(GetMollerPhiDiff(moller)) > 10.f)
            continue;

        if(theta1 < 0.8 || theta2 < 0.8)
            continue;

        n_moller++;

        const float p1[3] = {x[0][1], y[0][1], z[0][1]};
        const float p2[3] = {x[0][2], y[0][2], z[0][2]};
        const float q1[3] = {x[1][1], y[1][1], z[1][1]};
        const float q2[3] = {x[1][2], y[1][2], z[1][2]};

        float vx = 0.f, vy = 0.f, vz = 0.f, dca = 0.f;
        if (!getLineLineVertex(p1, p2, q1, q2, vx, vy, vz, dca)) continue;

        h_vz->Fill(vz);
        h_dca->Fill(dca);
        h_vxy->Fill(vx, vy);
        n_vertex++;

    }

    std::cout << std::endl;
    std::cout << "  Selected Moller events : " << n_moller << std::endl;
    std::cout << "  Reconstructed vertices : " << n_vertex << std::endl;

    TCanvas *c = new TCanvas("c_moller_vertex", "Moller vertex reconstruction", 1200, 450);
    c->Divide(3, 1);
    c->cd(1);
    h_vz->SetLineColor(kBlue + 1);
    h_vz->SetLineWidth(2);
    h_vz->Draw();
    c->cd(2);
    h_dca->SetLineColor(kRed + 1);
    h_dca->SetLineWidth(2);
    h_dca->Draw();
    c->cd(3);
    h_vxy->Draw("COLZ");
    c->SaveAs("vertex_recon_moller.png");

    TFile *fout = TFile::Open("vertex_recon_moller_output.root", "RECREATE");
    h_vz->Write();
    h_dca->Write();
    h_vxy->Write();
    c->Write();

    std::cout << "Saved to vertex_recon_moller_output.root" << std::endl;
}
