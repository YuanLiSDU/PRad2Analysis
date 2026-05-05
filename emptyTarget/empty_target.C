// Physical constants
static constexpr float M_PROTON  = 938.272f;   // MeV
static constexpr float M_ELECTRON = 0.511f;    // MeV
static constexpr float DEG2RAD = 3.14159265f / 180.f;

static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

namespace fs = std::filesystem;

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;
    float     EBeam = 0.f;

    // HyCal clusters
    float total_energy = 0.f;
    int     n_clusters = 0;
    float cl_x[kMaxClusters]       = {};
    float cl_y[kMaxClusters]       = {};
    float cl_z[kMaxClusters]       = {};
    float cl_energy[kMaxClusters]  = {};
    uint8_t cl_nblocks[kMaxClusters] = {};
    uint16_t cl_center[kMaxClusters]  = {};
    uint32_t cl_flag[kMaxClusters]    = {};

    // GEM reconstructed hits
    int        n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    uint8_t gem_x_size[kMaxGemHits]   = {};
    uint8_t gem_y_size[kMaxGemHits]   = {};

    // Detector matching results
    int      match_num = 0;
    float matchHC_x[kMaxClusters] = {};
    float matchHC_y[kMaxClusters] = {};
    float matchHC_z[kMaxClusters] = {};
    float matchHC_energy[kMaxClusters] = {};
    uint16_t matchHC_center[kMaxClusters] = {};
    uint32_t matchHC_flag[kMaxClusters] = {};
    float matchG_x[kMaxClusters][2] = {}; // up/down GEM for each cluster
    float matchG_y[kMaxClusters][2] = {};
    float matchG_z[kMaxClusters][2] = {};
    uint8_t matchG_det_id[kMaxClusters][2] = {};

    //veto informations
    int      veto_nch = 0;
    uint8_t veto_id[4]   = {}; // 0,1,2,3 for veto1-4
    int veto_npeaks[4] = {};
    float veto_peak_time[4][8]     = {};
    float veto_peak_integral[4][8] = {};

    //LMS reference PMT information
    int      lms_nch = 0;
    uint8_t lms_id[4]   = {}; // 0,1,2,3 for lms1-4
    int lms_npeaks[4] = {};
    float lms_peak_time[4][8]     = {};
    float lms_peak_integral[4][8] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;
};

// Aliases for the shared replay data structures
using EventVars_Recon = ReconEventData;

// ── Tree branch struct ───────────────────────────────────────────────────
void setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("total_energy", &ev.total_energy);
    tree->SetBranchAddress("EBeam",        &ev.EBeam);
    // HyCal cluster branches
    tree->SetBranchAddress("n_clusters",   &ev.n_clusters);
    tree->SetBranchAddress("cl_x",         ev.cl_x);
    tree->SetBranchAddress("cl_y",         ev.cl_y);
    tree->SetBranchAddress("cl_z",         ev.cl_z);
    tree->SetBranchAddress("cl_energy",    ev.cl_energy);
    tree->SetBranchAddress("cl_nblocks",   ev.cl_nblocks);
    tree->SetBranchAddress("cl_center",    ev.cl_center);
    tree->SetBranchAddress("cl_flag",      ev.cl_flag);
    // GEM part
    tree->SetBranchAddress("n_gem_hits",   &ev.n_gem_hits);
    tree->SetBranchAddress("det_id",       ev.det_id);
    tree->SetBranchAddress("gem_x",        ev.gem_x);
    tree->SetBranchAddress("gem_y",        ev.gem_y);
    tree->SetBranchAddress("gem_x_charge", ev.gem_x_charge);
    tree->SetBranchAddress("gem_y_charge", ev.gem_y_charge);
    tree->SetBranchAddress("gem_x_peak",   ev.gem_x_peak);
    tree->SetBranchAddress("gem_y_peak",   ev.gem_y_peak);
    tree->SetBranchAddress("gem_x_size",   ev.gem_x_size);
    tree->SetBranchAddress("gem_y_size",   ev.gem_y_size);
    // Matching results
    tree->SetBranchAddress("match_num",       &ev.match_num);
    tree->SetBranchAddress("matchHC_x",       ev.matchHC_x);
    tree->SetBranchAddress("matchHC_y",       ev.matchHC_y);
    tree->SetBranchAddress("matchHC_z",       ev.matchHC_z);
    tree->SetBranchAddress("matchHC_energy",  ev.matchHC_energy);
    tree->SetBranchAddress("matchHC_center",  ev.matchHC_center);
    tree->SetBranchAddress("matchHC_flag",    ev.matchHC_flag);
    tree->SetBranchAddress("matchG_x",        ev.matchG_x);
    tree->SetBranchAddress("matchG_y",        ev.matchG_y);
    tree->SetBranchAddress("matchG_z",        ev.matchG_z);
    tree->SetBranchAddress("matchG_det_id",   ev.matchG_det_id);
}

struct DataPoint
{
    float x;
    float y;
    float z;
    float E;

    DataPoint() {};
    DataPoint(float xi, float yi, float zi, float Ei) : x(xi), y(yi), z(zi), E(Ei) {};
};
typedef std::pair<DataPoint, DataPoint> MollerEvent;
typedef std::vector<MollerEvent> MollerData;

static std::vector<std::string> collectRootFiles(const std::string &path);
float ExpectedEnergy(float theta_deg, float EBeam, const std::string &type);
std::array<float, 2> GetMollerCenter(MollerEvent &event1, MollerEvent &event2);

float GetPhiAngle(float x, float y)
{
    float phi = std::atan2(y, x) * 180.f / static_cast<float>(TMath::Pi());
    if (phi < 0) phi += 360.f;
    return phi;
}
float GetMollerPhiDiff(MollerEvent &event1)
{
    // Calculate the azimuthal angle difference (phi) for a Moller event
    float x1 = event1.first.x, y1 = event1.first.y;
    float x2 = event1.second.x, y2 = event1.second.y;
    float phi1 = GetPhiAngle(x1, y1);
    float phi2 = GetPhiAngle(x2, y2);
    float phi_diff = fabs(phi1 - phi2) - 180.f; // Expecting back-to-back, so difference should be around 180 degrees
    return phi_diff;
}

void processType(TTree *tree, EventVars_Recon &ev,
                 TH2F *E_theta, TH2F *E_theta_mott, TH2F *E_theta_moller,
                 TH2F *all_hits, TH2F *mott_hits, TH2F *moller_hits,
                 TH2F *moller_center, TH1F *moller_vertex,
                 TH1F *mott_yield, TH1F *moller_yield,
                 MollerData &mollerData, const std::string &typeName)
{
    Long64_t nEntries = tree->GetEntries();
    cout << "Processing Type " << typeName << " events: " << nEntries << " entries found." << endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if(i % 1000 == 0) cout << "Processed " << i << " / " << nEntries << " events." << "\r" << flush;
        for (int j = 0; j < ev.n_clusters; ++j) {
            float x = ev.cl_x[j];
            float y = ev.cl_y[j];
            float z = ev.cl_z[j];
            float energy = ev.cl_energy[j];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180 / M_PI;
            E_theta->Fill(theta, energy);
            all_hits->Fill(x, y);
            if (fabs(energy - ev.EBeam) < 3. * ev.EBeam * 0.026 / sqrt(ev.EBeam/1000.)) { // Mott selection
                E_theta_mott->Fill(theta, energy);
                mott_hits->Fill(x, y);
                mott_yield->Fill(theta);
            }
        }
        if(ev.n_clusters == 2){ // Moller events should have 2 clusters
            float x1 = ev.cl_x[0], x2 = ev.cl_x[1];
            float y1 = ev.cl_y[0], y2 = ev.cl_y[1];
            float z1 = ev.cl_z[0], z2 = ev.cl_z[1];
            float energy1 = ev.cl_energy[0], energy2 = ev.cl_energy[1];
            float theta1 = std::atan(std::sqrt(x1*x1 + y1*y1) / z1) * 180 / M_PI;
            float theta2 = std::atan(std::sqrt(x2*x2 + y2*y2) / z2) * 180 / M_PI;
            float expectedMollerE1 = ExpectedEnergy(theta1, ev.EBeam, "ee");
            float expectedMollerE2 = ExpectedEnergy(theta2, ev.EBeam, "ee");
            if(fabs(energy1 - expectedMollerE1) < 5. * expectedMollerE1 * 0.026 / std::sqrt(expectedMollerE1/1000.) &&
               fabs(energy2 - expectedMollerE2) < 5. * expectedMollerE2 * 0.026 / std::sqrt(expectedMollerE2/1000.) &&
               fabs(energy1 + energy2 - ev.EBeam) < 3. * ev.EBeam * 0.026 / std::sqrt(ev.EBeam/1000.)) { // Moller selection
                //cout << ev.EBeam << " " << theta1 << " " << energy1 << " " << expectedMollerE1 << " | "
                  //   << theta2 << " " << energy2 << " " << expectedMollerE2 << endl;
                MollerEvent event(DataPoint(x1, y1, z1, energy1), DataPoint(x2, y2, z2, energy2));
                float phi_diff = GetMollerPhiDiff(event);
                if(fabs(phi_diff) > 10.f) continue; // Require back-to-back in phi within 10 degrees

                E_theta_moller->Fill(theta1, energy1);
                E_theta_moller->Fill(theta2, energy2);
                moller_hits->Fill(x1, y1);
                moller_hits->Fill(x2, y2);
                moller_yield->Fill(theta1);
                moller_yield->Fill(theta2);
                float vertex = sqrt( (ev.EBeam + M_ELECTRON) * sqrt(x1*x1 + y1*y1) * sqrt(x2*x2 + y2*y2) / (2. * M_ELECTRON) );
                moller_vertex->Fill(vertex);
                mollerData.push_back(event);
            }
        }
    }
    for (int i = 1; i < (int)mollerData.size(); i++) {
        auto center = GetMollerCenter(mollerData[i-1], mollerData[i]);
        moller_center->Fill(center[0], center[1]);
    }
    cout << "Type " << typeName << " processing completed." << endl;
}

void empty_target(){

    std::vector<std::string> A_files, B_files, C_files, D_files;
    auto f = collectRootFiles("typeA/");
    A_files.insert(A_files.end(), f.begin(), f.end());
        f = collectRootFiles("typeB/");
    B_files.insert(B_files.end(), f.begin(), f.end());
        f = collectRootFiles("typeC/");
    C_files.insert(C_files.end(), f.begin(), f.end());
        f = collectRootFiles("typeD/");
    D_files.insert(D_files.end(), f.begin(), f.end());

    TChain *Achain = new TChain("recon");
    for (const auto &file : A_files) Achain->Add(file.c_str());
    TChain *Bchain = new TChain("recon");
    for (const auto &file : B_files) Bchain->Add(file.c_str());
    TChain *Cchain = new TChain("recon");
    for (const auto &file : C_files) Cchain->Add(file.c_str());
    TChain *Dchain = new TChain("recon");
    for (const auto &file : D_files) Dchain->Add(file.c_str());

    TTree *Atree = Achain;
    TTree *Btree = Bchain;
    TTree *Ctree = Cchain;
    TTree *Dtree = Dchain;

    EventVars_Recon A_ev, B_ev, C_ev, D_ev;
    setupReconBranches(Atree, A_ev);
    setupReconBranches(Btree, B_ev);
    setupReconBranches(Ctree, C_ev);
    setupReconBranches(Dtree, D_ev);

    // output histograms
    //no selection, all clusters
    TH2F *E_theta_A = new TH2F("E_theta_A", "Type A: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_B = new TH2F("E_theta_B", "Type B: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_C = new TH2F("E_theta_C", "Type C: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_D = new TH2F("E_theta_D", "Type D: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    //after selection, mott and moller
    TH2F *E_theta_A_mott = new TH2F("E_theta_A_mott", "Type A Mott: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_B_mott = new TH2F("E_theta_B_mott", "Type B Mott: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_C_mott = new TH2F("E_theta_C_mott", "Type C Mott: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_D_mott = new TH2F("E_theta_D_mott", "Type D Mott: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_A_moller = new TH2F("E_theta_A_moller", "Type A Moller: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_B_moller = new TH2F("E_theta_B_moller", "Type B Moller: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_C_moller = new TH2F("E_theta_C_moller", "Type C Moller: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    TH2F *E_theta_D_moller = new TH2F("E_theta_D_moller", "Type D Moller: E vs #theta;#theta (deg);Energy (GeV)", 57, 0.3, 6, 1000, 0, 4000);
    //hit position distribution
    TH2F *all_hits_A = new TH2F("all_hits_A", "Type A: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *all_hits_B = new TH2F("all_hits_B", "Type B: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *all_hits_C = new TH2F("all_hits_C", "Type C: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *all_hits_D = new TH2F("all_hits_D", "Type D: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *mott_hits_A = new TH2F("mott_hits_A", "Type A Mott: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *mott_hits_B = new TH2F("mott_hits_B", "Type B Mott: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *mott_hits_C = new TH2F("mott_hits_C", "Type C Mott: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *mott_hits_D = new TH2F("mott_hits_D", "Type D Mott: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *moller_hits_A = new TH2F("moller_hits_A", "Type A Moller: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *moller_hits_B = new TH2F("moller_hits_B", "Type B Moller: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *moller_hits_C = new TH2F("moller_hits_C", "Type C Moller: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    TH2F *moller_hits_D = new TH2F("moller_hits_D", "Type D Moller: hit position distribution;X (mm);Y (mm)", 300, -360, 360, 300, -360, 360);
    //Moller center and moller vertex distribution
    TH2F *moller_center_A = new TH2F("moller_center_A", "Type A Moller: cluster center distribution;X (mm);Y (mm)", 100, -50, 50, 100, -50, 50);
    TH2F *moller_center_B = new TH2F("moller_center_B", "Type B Moller: cluster center distribution;X (mm);Y (mm)", 100, -50, 50, 100, -50, 50);
    TH2F *moller_center_C = new TH2F("moller_center_C", "Type C Moller: cluster center distribution;X (mm);Y (mm)", 100, -50, 50, 100, -50, 50);
    TH2F *moller_center_D = new TH2F("moller_center_D", "Type D Moller: cluster center distribution;X (mm);Y (mm)", 100, -50, 50, 100, -50, 50);
    TH1F *moller_vertex_A = new TH1F("moller_vertex_A", "Type A Moller: vertex distribution;X (mm);Y (mm)", 500, 5000,6000);
    TH1F *moller_vertex_B = new TH1F("moller_vertex_B", "Type B Moller: vertex distribution;X (mm);Y (mm)", 500, 5000,6000);
    TH1F *moller_vertex_C = new TH1F("moller_vertex_C", "Type C Moller: vertex distribution;X (mm);Y (mm)", 500, 5000,6000);
    TH1F *moller_vertex_D = new TH1F("moller_vertex_D", "Type D Moller: vertex distribution;X (mm);Y (mm)", 500, 5000,6000);
    //mott and moller yields
    TH1F *mott_yield_A = new TH1F("mott_yield_A", "Type A Mott: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *mott_yield_B = new TH1F("mott_yield_B", "Type B Mott: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *mott_yield_C = new TH1F("mott_yield_C", "Type C Mott: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *mott_yield_D = new TH1F("mott_yield_D", "Type D Mott: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *moller_yield_A = new TH1F("moller_yield_A", "Type A Moller: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *moller_yield_B = new TH1F("moller_yield_B", "Type B Moller: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *moller_yield_C = new TH1F("moller_yield_C", "Type C Moller: yield;Theta (deg);Yield", 57, 0.3, 6);
    TH1F *moller_yield_D = new TH1F("moller_yield_D", "Type D Moller: yield;Theta (deg);Yield", 57, 0.3, 6);

    MollerData mollerData_A, mollerData_B, mollerData_C, mollerData_D;

    // Loop over events and fill histograms
    processType(Atree, A_ev, E_theta_A, E_theta_A_mott, E_theta_A_moller, all_hits_A, mott_hits_A, moller_hits_A, moller_center_A, moller_vertex_A, mott_yield_A, moller_yield_A, mollerData_A, "A");
    processType(Btree, B_ev, E_theta_B, E_theta_B_mott, E_theta_B_moller, all_hits_B, mott_hits_B, moller_hits_B, moller_center_B, moller_vertex_B, mott_yield_B, moller_yield_B, mollerData_B, "B");
    processType(Ctree, C_ev, E_theta_C, E_theta_C_mott, E_theta_C_moller, all_hits_C, mott_hits_C, moller_hits_C, moller_center_C, moller_vertex_C, mott_yield_C, moller_yield_C, mollerData_C, "C");
    processType(Dtree, D_ev, E_theta_D, E_theta_D_mott, E_theta_D_moller, all_hits_D, mott_hits_D, moller_hits_D, moller_center_D, moller_vertex_D, mott_yield_D, moller_yield_D, mollerData_D, "D");

    //Draw histograms
    TCanvas *c1 = new TCanvas("c1", "Energy vs Theta", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); E_theta_A->Draw("COLZ");
    c1->cd(2); E_theta_B->Draw("COLZ");
    c1->cd(3); E_theta_C->Draw("COLZ");
    c1->cd(4); E_theta_D->Draw("COLZ");
    c1->SaveAs("E_theta_all.png");

    TCanvas *c2 = new TCanvas("c2", "Energy vs Theta (Mott)", 1200, 800);
    c2->Divide(2,2);
    c2->cd(1); E_theta_A_mott->Draw("COLZ");
    c2->cd(2); E_theta_B_mott->Draw("COLZ");
    c2->cd(3); E_theta_C_mott->Draw("COLZ");
    c2->cd(4); E_theta_D_mott->Draw("COLZ");
    c2->SaveAs("E_theta_mott.png");

    TCanvas *c3 = new TCanvas("c3", "Energy vs Theta (Moller)", 1200, 800);
    c3->Divide(2,2);
    c3->cd(1); E_theta_A_moller->Draw("COLZ");
    c3->cd(2); E_theta_B_moller->Draw("COLZ");
    c3->cd(3); E_theta_C_moller->Draw("COLZ");
    c3->cd(4); E_theta_D_moller->Draw("COLZ");
    c3->SaveAs("E_theta_moller.png");

    TCanvas *c4 = new TCanvas("c4", "Hit Position Distribution", 1200, 800);
    c4->Divide(2,2);
    c4->cd(1); all_hits_A->Draw("COLZ");
    c4->cd(2); all_hits_B->Draw("COLZ");
    c4->cd(3); all_hits_C->Draw("COLZ");
    c4->cd(4); all_hits_D->Draw("COLZ");
    c4->SaveAs("hit_position_all.png");

    TCanvas *c5 = new TCanvas("c5", "Hit Position Distribution (Mott)", 1200, 800);
    c5->Divide(2,2);
    c5->cd(1); mott_hits_A->Draw("COLZ");
    c5->cd(2); mott_hits_B->Draw("COLZ");
    c5->cd(3); mott_hits_C->Draw("COLZ");
    c5->cd(4); mott_hits_D->Draw("COLZ");
    c5->SaveAs("hit_position_mott.png");

    TCanvas *c6 = new TCanvas("c6", "Hit Position Distribution (Moller)", 1200, 800);
    c6->Divide(2,2);
    c6->cd(1); moller_hits_A->Draw("COLZ");
    c6->cd(2); moller_hits_B->Draw("COLZ");
    c6->cd(3); moller_hits_C->Draw("COLZ");
    c6->cd(4); moller_hits_D->Draw("COLZ");
    c6->SaveAs("hit_position_moller.png");

    TCanvas *c7 = new TCanvas("c7", "Moller Cluster Center Distribution", 1200, 800);
    c7->Divide(2,2);
    c7->cd(1); moller_center_A->Draw("COLZ");
    c7->cd(2); moller_center_B->Draw("COLZ");
    c7->cd(3); moller_center_C->Draw("COLZ");
    c7->cd(4); moller_center_D->Draw("COLZ");
    c7->SaveAs("moller_cluster_center.png");

    TCanvas *c8 = new TCanvas("c8", "Moller Vertex Distribution", 1200, 800);
    c8->Divide(2,2);
    c8->cd(1); moller_vertex_A->Draw();
    c8->cd(2); moller_vertex_B->Draw();
    c8->cd(3); moller_vertex_C->Draw();
    c8->cd(4); moller_vertex_D->Draw();
    c8->SaveAs("moller_vertex.png");

    TCanvas *c9 = new TCanvas("c9", "Mott Yield", 1200, 800);
    c9->Divide(2,2);
    c9->cd(1); mott_yield_A->Draw();
    c9->cd(2); mott_yield_B->Draw();
    c9->cd(3); mott_yield_C->Draw();
    c9->cd(4); mott_yield_D->Draw();
    c9->SaveAs("mott_yield.png");

    TCanvas *c10 = new TCanvas("c10", "Moller Yield", 1200, 800);
    c10->Divide(2,2);
    c10->cd(1); moller_yield_A->Draw();
    c10->cd(2); moller_yield_B->Draw();
    c10->cd(3); moller_yield_C->Draw();
    c10->cd(4); moller_yield_D->Draw();
    c10->SaveAs("moller_yield.png");



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

float ExpectedEnergy(float theta_deg, float EBeam, const std::string &type)
{
    float theta = theta_deg * DEG2RAD;
    float cos_t = std::cos(theta);
    float sin_t = std::sin(theta);

    if (type == "ep") {
        // elastic e-p: E' = E * M / (M + E*(1 - cos_t))
        // where M = proton mass
        float expectE = EBeam * M_PROTON / (M_PROTON + EBeam * (1.f - cos_t));
        return expectE;
    }
    if (type == "ee") {
        // Moller scattering: exact lab-frame formula from 4-momentum conservation
        // E' = m * [(gamma+1) + (gamma-1)*cos^2(theta)] / [(gamma+1) - (gamma-1)*cos^2(theta)]
        float gamma = EBeam / M_ELECTRON;
        float num = (gamma + 1.f) + (gamma - 1.f) * cos_t * cos_t;
        float den = (gamma + 1.f) - (gamma - 1.f) * cos_t * cos_t;
        if (den <= 0) return 0.f;
        float expectE = M_ELECTRON * num / den;
        return expectE;
    }
    return 0.f;
}

std::array<float, 2> GetMollerCenter(MollerEvent &event1, MollerEvent &event2)
{
    float x1[2], y1[2];
    float x2[2], y2[2];

    x1[0] = event1.first.x; y1[0] = event1.first.y;
    x1[1] = event1.second.x; y1[1] = event1.second.y;
    x2[0] = event2.first.x; y2[0] = event2.first.y;
    x2[1] = event2.second.x; y2[1] = event2.second.y;

    //two lines: y = ax + b, y = cx + d
    float dx1 = x1[0] - x1[1];
    float dx2 = x2[0] - x2[1];
    if (std::abs(dx1) < 1e-6f || std::abs(dx2) < 1e-6f)
        return {0.f, 0.f};  // vertical line — degenerate

    float a = (y1[0] - y1[1]) / dx1;
    float b = y1[0] - a * x1[0];
    float c = (y2[0] - y2[1]) / dx2;
    float d = y2[0] - c * x2[0];

    if (std::abs(a - c) < 1e-6f)
        return {0.f, 0.f};  // parallel lines — no intersection

    float x_cross = (d - b) / (a - c);
    float y_cross = a * x_cross + b;

    return {x_cross, y_cross};

}