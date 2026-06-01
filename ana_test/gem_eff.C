#include "../EventData.h"
#include "../PhysicsTools.h"

const int Nbins = 33;
const float binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

float Ebeam = 3484.f; // MeV
float resolution = 0.033f;

// Returns true if (x, y) [mm] lies inside the HyCal active acceptance:
// outside the beam hole (2.5 module widths) and inside the outer edge (16 module widths).
bool inHyCal(double xmm, double ymm) {
    const double module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
}

void gem_eff(){

    TFile *f = TFile::Open("../../A/24917_recon_filter.root");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file" << std::endl;
        return;
    }
    TTree *tree = new TTree();
    f->GetObject("recon", tree);

    ReconEventData ev;
    setupReconBranches(tree, ev);

    float gem_up_should[Nbins]= {};
    float gem_down_should[Nbins]= {};
    float gem_up_count[Nbins] = {};
    float gem_down_count[Nbins] = {};

    TH1F *gem_up_eff = new TH1F("gem_up_eff", "Upstream GEM Efficiency;Theta (deg);Efficiency", Nbins, binEdge);
    TH1F *gem_down_eff = new TH1F("gem_down_eff", "Downstream GEM Efficiency;Theta (deg);Efficiency", Nbins, binEdge);
    TH1F *gem_eff = new TH1F("gem_eff", "Overall GEM Efficiency;Theta (deg);Efficiency", Nbins, binEdge);

    for (Long64_t i = 0; i < tree->GetEntries()/100; i++) {
        if (i % 10000 == 0)
            std::cout << "  " << i << " / " << tree->GetEntries() << "\r" << std::flush;
        tree->GetEntry(i);

        if(ev.n_clusters != 1) continue;
        if(ev.matchFlag[0] == 0) continue; 

        //one cluster with upstream GEM(2,3) match
        //calculate downstream GEM(0,1) efficiency
        if(ev.matchFlag[0] & (1<<2 | 1<<3)){ // matched with GEM2 or GEM3
            int id = -1;
            if(ev.matchFlag[0] & (1<<2)) id = 2;
            else if(ev.matchFlag[0] & (1<<3)) id = 3;

            float z = ev.matchGEMz[0][id];
            float scale = 6270.f / z;
            float x = ev.matchGEMx[0][id] * scale;
            float y = ev.matchGEMy[0][id] * scale;
            z = 6270.f;

            if(!inHyCal(x, y)) continue;

            float theta = atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            int bin = gem_eff->FindBin(theta);
            gem_down_should[bin-1]++;

            if(ev.matchFlag[0] & (1<<0 | 1<<1)) // also matched with GEM0 or GEM1
                gem_down_count[bin-1]++;
        }

        //one cluster with downstream GEM(0,1) match
        //calculate downstream GEM(2,3) efficiency
        if(ev.matchFlag[0] & (1<<0 | 1<<1)){ // matched with GEM0 or GEM1
            int id = -1;
            if(ev.matchFlag[0] & (1<<0)) id = 0;
            else if(ev.matchFlag[0] & (1<<1)) id = 1;

            float z = ev.matchGEMz[0][id];
            float scale = 6270.f / z;
            float x = ev.matchGEMx[0][id] * scale;
            float y = ev.matchGEMy[0][id] * scale;
            z = 6270.f;

            if(!inHyCal(x, y)) continue;

            float theta = atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            int bin = gem_eff->FindBin(theta);
            gem_up_should[bin-1]++;

            if(ev.matchFlag[0] & (1<<2 | 1<<3)) // also matched with GEM2 or GEM3 (upstream)
                gem_up_count[bin-1]++;
        }

    }

    for(int i=0; i<Nbins; i++){
        float up_eff = (gem_up_should[i] > 0) ? gem_up_count[i] / gem_up_should[i] : 0.f;
        float down_eff = (gem_down_should[i] > 0) ? gem_down_count[i] / gem_down_should[i] : 0.f;
        float overall_eff = up_eff * down_eff; // assuming independent efficiencies for upstream and downstream GEMs

        //error bars using binomial
        float up_err = (gem_up_should[i] > 0) ? sqrt(up_eff * (1 - up_eff) / gem_up_should[i]) : 0.f;
        float down_err = (gem_down_should[i] > 0) ? sqrt(down_eff * (1 - down_eff) / gem_down_should[i]) : 0.f;
        float overall_err = sqrt(pow(up_err * down_eff, 2) + pow(down_err * up_eff, 2)); // propagate errors for product of two independent variables

        gem_up_eff->SetBinContent(i+1, up_eff);
        gem_up_eff->SetBinError(i+1, up_err);
        gem_down_eff->SetBinContent(i+1, down_eff);
        gem_down_eff->SetBinError(i+1, down_err);
        gem_eff->SetBinContent(i+1, overall_eff);
        gem_eff->SetBinError(i+1, overall_err);
    }

    TCanvas *c = new TCanvas("c", "GEM Efficiency", 1200, 600);
    c->SetGrid();
    gem_up_eff->SetLineColor(kRed);
    gem_up_eff->SetMarkerColor(kRed);
    gem_up_eff->SetMarkerStyle(20);
    gem_up_eff->SetTitle("GEM Efficiency;Theta (deg);Efficiency");
    gem_up_eff->SetStats(0);
    gem_up_eff->GetYaxis()->SetRangeUser(0.7, 1.05);
    gem_up_eff->Draw("E1P");
    gem_down_eff->SetLineColor(kBlue);
    gem_down_eff->SetMarkerColor(kBlue);
    gem_down_eff->SetMarkerStyle(20);
    gem_down_eff->SetStats(0);
    gem_down_eff->Draw("E1P SAME");
    gem_eff->SetLineColor(kBlack);
    gem_eff->SetMarkerColor(kBlack);
    gem_eff->SetMarkerStyle(20);
    gem_eff->SetStats(0);
    gem_eff->Draw("E1P SAME");

    TLegend *leg = new TLegend(0.70, 0.15, 0.92, 0.40);
    leg->SetBorderSize(0);
    leg->AddEntry(gem_up_eff,   "Upstream GEM",   "lp");
    leg->AddEntry(gem_down_eff, "Downstream GEM", "lp");
    leg->AddEntry(gem_eff,      "Overall",        "lp");
    leg->Draw();

    c->SaveAs("gem_efficiency.png");

    cout << "GEM Efficiency Results:" << endl;
    for(int i=0; i<Nbins; i++){
        cout << gem_eff->GetBinCenter(i+1) << ", ";
    }

}