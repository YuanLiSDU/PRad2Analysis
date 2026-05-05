static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    Float_t     EBeam;

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

    //veto information
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

void setupReconBranches(TTree *tree, ReconEventData &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("total_energy", &ev.total_energy, "total_energy/F");
    tree->Branch("EBeam",        &ev.EBeam,        "EBeam/F");
    // HyCal cluster branches
    // detector coordinate system (crystal surface)
    tree->Branch("n_clusters",   &ev.n_clusters,   "n_clusters/I");
    tree->Branch("cl_x",         ev.cl_x,          "cl_x[n_clusters]/F");
    tree->Branch("cl_y",         ev.cl_y,          "cl_y[n_clusters]/F");
    tree->Branch("cl_z",         ev.cl_z,          "cl_z[n_clusters]/F");
    tree->Branch("cl_energy",    ev.cl_energy,     "cl_energy[n_clusters]/F");
    tree->Branch("cl_nblocks",   ev.cl_nblocks,    "cl_nblocks[n_clusters]/b");
    tree->Branch("cl_center",    ev.cl_center,     "cl_center[n_clusters]/s");
    tree->Branch("cl_flag",      ev.cl_flag,       "cl_flag[n_clusters]/i");
    // GEM part
    //detector coordinate system (GEM plane)
    tree->Branch("n_gem_hits",   &ev.n_gem_hits,   "n_gem_hits/I");
    tree->Branch("det_id",       ev.det_id,        "det_id[n_gem_hits]/b");
    tree->Branch("gem_x",        ev.gem_x,         "gem_x[n_gem_hits]/F");
    tree->Branch("gem_y",        ev.gem_y,         "gem_y[n_gem_hits]/F");
    tree->Branch("gem_x_charge", ev.gem_x_charge,  "gem_x_charge[n_gem_hits]/F");
    tree->Branch("gem_y_charge", ev.gem_y_charge,  "gem_y_charge[n_gem_hits]/F");
    tree->Branch("gem_x_peak",   ev.gem_x_peak,    "gem_x_peak[n_gem_hits]/F");
    tree->Branch("gem_y_peak",   ev.gem_y_peak,    "gem_y_peak[n_gem_hits]/F");
    tree->Branch("gem_x_size",   ev.gem_x_size,    "gem_x_size[n_gem_hits]/b");
    tree->Branch("gem_y_size",   ev.gem_y_size,    "gem_y_size[n_gem_hits]/b");
    // Matching results
    // the x,y,z positions in this part are in the target and beam center coordinate system
    tree->Branch("match_num",       &ev.match_num,       "match_num/I");
    tree->Branch("matchHC_x",       ev.matchHC_x,        "matchHC_x[match_num]/F");
    tree->Branch("matchHC_y",       ev.matchHC_y,        "matchHC_y[match_num]/F");
    tree->Branch("matchHC_z",       ev.matchHC_z,        "matchHC_z[match_num]/F");
    tree->Branch("matchHC_energy",  ev.matchHC_energy,   "matchHC_energy[match_num]/F");
    tree->Branch("matchHC_center",  ev.matchHC_center,   "matchHC_center[match_num]/s");
    tree->Branch("matchHC_flag",    ev.matchHC_flag,     "matchHC_flag[match_num]/i");
    tree->Branch("matchG_x",        ev.matchG_x,         "matchG_x[match_num][2]/F");
    tree->Branch("matchG_y",        ev.matchG_y,         "matchG_y[match_num][2]/F");
    tree->Branch("matchG_z",        ev.matchG_z,         "matchG_z[match_num][2]/F");
    tree->Branch("matchG_det_id",   ev.matchG_det_id,    "matchG_det_id[match_num][2]/b");
    //veto information
    tree->Branch("veto_nch",       &ev.veto_nch,       "veto_nch/I");
    tree->Branch("veto_id",        ev.veto_id,        "veto_id[veto_nch]/b");
    tree->Branch("veto_npeaks",       &ev.veto_npeaks,       "veto_npeaks[veto_nch]/I");
    tree->Branch("veto_peak_time",    ev.veto_peak_time,     Form("veto_peak_time[veto_nch][%d]/F", 8));
    tree->Branch("veto_peak_integral",ev.veto_peak_integral, Form("veto_peak_integral[veto_nch][%d]/F", 8));
    //LMS information
    tree->Branch("lms_nch",       &ev.lms_nch,       "lms_nch/I");
    tree->Branch("lms_id",        ev.lms_id,         "lms_id[lms_nch]/b");
    tree->Branch("lms_npeaks",       &ev.lms_npeaks,       "lms_npeaks[lms_nch]/I");
    tree->Branch("lms_peak_time",    ev.lms_peak_time,     Form("lms_peak_time[lms_nch][%d]/F", 8));
    tree->Branch("lms_peak_integral",ev.lms_peak_integral, Form("lms_peak_integral[lms_nch][%d]/F", 8));
    // Raw 0xE10C SSP trigger bank words (variable-length per event)
    tree->Branch("ssp_raw", &ev.ssp_raw);
};

struct iputEventData {
     // ── Input tree variables ─────────────────────────────────────────────
    // Event-level
    Int_t       in_EventNumber;
    Int_t       in_nHit;
    Float_t     in_EBeam;
    Float_t     in_TotalE;
    UShort_t    in_trgType;
    UShort_t    in_trgTime;

    // Hit-level arrays (variable-length, indexed by nHit)
    static constexpr int kMaxHits = 500;
    UInt_t      in_Hit_Flag[kMaxHits];
    Float_t     in_Hit_X[kMaxHits];
    Float_t     in_Hit_Y[kMaxHits];
    Float_t     in_Hit_Z[kMaxHits];
    Float_t     in_Hit_E[kMaxHits];
    Short_t     in_Hit_NModule[kMaxHits];
    Short_t     in_Hit_CID[kMaxHits];
};

void setupInputBranches(TTree *intree, iputEventData &inputData)
{
// ── SetBranchAddress ─────────────────────────────────────────────────
    intree->SetBranchAddress("EventNumber",      &inputData.in_EventNumber);
    intree->SetBranchAddress("nHit",             &inputData.in_nHit);
    intree->SetBranchAddress("EBeam",            &inputData.in_EBeam);
    intree->SetBranchAddress("TotalE",           &inputData.in_TotalE);
    intree->SetBranchAddress("TriggerType",      &inputData.in_trgType);
    intree->SetBranchAddress("TriggerTime",      &inputData.in_trgTime);

    intree->SetBranchAddress("Hit.Flag",         inputData.in_Hit_Flag);
    intree->SetBranchAddress("Hit.X",            inputData.in_Hit_X);
    intree->SetBranchAddress("Hit.Y",            inputData.in_Hit_Y);
    intree->SetBranchAddress("Hit.Z",            inputData.in_Hit_Z);
    intree->SetBranchAddress("Hit.E",            inputData.in_Hit_E);
    intree->SetBranchAddress("Hit.NModule",      inputData.in_Hit_NModule);
    intree->SetBranchAddress("Hit.CID",          inputData.in_Hit_CID);
}

void convert(string inputFile, string outputFile);

void convert2prad2()
{   
    for(int run = 1314; run <=1314; run++){
        string inputFile = Form("typeA/physReplay_newAna_%d.root", run);
        string outputFile = Form("prad2Replay_%d.root", run);
        convert(inputFile, outputFile);
    }
}

void convert(string inputFile, string outputFile)
{   
    TFile *infile = TFile::Open(inputFile.c_str(), "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }
    TTree *intree = (TTree*)infile->Get("T");
    if (!intree) {
        std::cerr << "Error: TTree 'events' not found in input file!" << std::endl;
        infile->Close();
        return;
    }
    iputEventData inputData;
    setupInputBranches(intree, inputData);

    TFile *outfile = TFile::Open(outputFile.c_str(), "RECREATE");
    TTree *outtree = new TTree("recon", "PRad2 replay reconstruction");
    ReconEventData ev;
    setupReconBranches(outtree, ev);

    Long64_t nEntries = intree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        intree->GetEntry(i);
        ev.event_num = inputData.in_EventNumber;
        ev.EBeam = inputData.in_EBeam;
        ev.trigger_type = inputData.in_trgType;
        ev.timestamp = inputData.in_trgTime;
        ev.total_energy = inputData.in_TotalE;
        ev.n_clusters = inputData.in_nHit;
        ev.trigger_bits = 0;
        for (int j = 0; j < inputData.in_nHit; ++j) {
            ev.cl_x[j] = inputData.in_Hit_X[j];
            ev.cl_y[j] = inputData.in_Hit_Y[j];
            ev.cl_z[j] = inputData.in_Hit_Z[j];
            ev.cl_energy[j] = inputData.in_Hit_E[j];
            ev.cl_nblocks[j] = inputData.in_Hit_NModule[j];
            ev.cl_center[j] = inputData.in_Hit_CID[j];
            ev.cl_flag[j] = inputData.in_Hit_Flag[j];
        }
        outtree->Fill();
    }
    outtree->Write();
    outfile->Close();
    infile->Close();
    cout << "Conversion completed successfully! file: " << outputFile << std::endl;
}