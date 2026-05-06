static constexpr int MAX_SAMPLES  = 200;   // samples per channel per event
static constexpr int MAX_SLOTS    = 32;    // slot IDs 0..31 (VME 3-20, Fastbus 0-25)
static constexpr int MAX_ROCS     = 10;    // number of ROC crates
static constexpr int MAX_PEAKS    = 8;  

static constexpr int APV_STRIP_SIZE     = 128;   // channels per APV25 chip
static constexpr int SSP_TIME_SAMPLES   = 6;     // fixed by SSP firmware
static constexpr int MAX_APVS_PER_MPD   = 16;    // APV slots per MPD
static constexpr int MAX_MPDS           = 64;    // MPDs across all crates

static constexpr int kMaxChannels  = MAX_ROCS * MAX_SLOTS * 16;
static constexpr int kMaxGemStrips = MAX_MPDS * MAX_APVS_PER_MPD * APV_STRIP_SIZE;
static constexpr int kMaxClusters  = 100;
static constexpr int kMaxGemHits   = 400;

// ── Raw replay ("events" tree) ───────────────────────────────────────────

struct RawEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits      = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;

    // HyCal per-channel data
    int          nch = 0;
    uint16_t     module_id[kMaxChannels] = {};
    int nsamples[kMaxChannels] = {};
    uint16_t     samples[kMaxChannels][MAX_SAMPLES] = {};
    float   ped_mean[kMaxChannels] = {};
    float   ped_rms[kMaxChannels]  = {};
    float   integral[kMaxChannels] = {};
    float   gain_factor[kMaxChannels] = {};

    //Veto per-channel data
    int          veto_nch = 0;
    uint8_t veto_id[4]   = {}; // 1,2,3,4 for veto1-4
    int veto_nsamples[4] = {};
    uint16_t     veto_samples[4][MAX_SAMPLES] = {};
    float   veto_ped_mean[4] = {};
    float   veto_ped_rms[4]  = {};
    float   veto_integral[4] = {};

    //LMS reference PMT data
    int lms_nch = 0;
    uint8_t lms_id[4] = {}; // 1,2,3 for lms1-3, 0 for Pin
    int lms_nsamples[4] = {};
    uint16_t lms_samples[4][MAX_SAMPLES] = {};
    float   lms_ped_mean[4] = {};
    float   lms_ped_rms[4]  = {};
    float   lms_integral[4] = {};

    // Optional peak data
    int npeaks[kMaxChannels] = {};
    float   peak_height[kMaxChannels][MAX_PEAKS]   = {};
    float   peak_time[kMaxChannels][MAX_PEAKS]     = {};
    float   peak_integral[kMaxChannels][MAX_PEAKS] = {};

    //optional veto peak data
    int veto_npeaks[4] = {};
    float   veto_peak_height[4][MAX_PEAKS]   = {};
    float   veto_peak_time[4][MAX_PEAKS]     = {};
    float   veto_peak_integral[4][MAX_PEAKS] = {};

    //optional LMS peak data
    int lms_npeaks[4] = {};
    float   lms_peak_height[4][MAX_PEAKS]   = {};
    float   lms_peak_time[4][MAX_PEAKS]     = {};
    float   lms_peak_integral[4][MAX_PEAKS] = {};

    // GEM per-strip data
    int        gem_nch = 0;
    uint8_t mpd_crate[kMaxGemStrips]  = {};
    uint8_t mpd_fiber[kMaxGemStrips]  = {};
    uint8_t apv[kMaxGemStrips]        = {};
    uint8_t strip[kMaxGemStrips]      = {};
    int16_t ssp_samples[kMaxGemStrips][SSP_TIME_SAMPLES] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;
};

// ── Reconstructed replay ("recon" tree) ──────────────────────────────────

struct ReconEventData {
    int      event_num    = 0;
    uint8_t  trigger_type = 0;   // main trigger (from event tag: tag - 0x80)
    uint32_t trigger_bits = 0;   // FP trigger bits (multi-bit, from TI master d[5])
    long long  timestamp    = 0;
    float EBeam = 0.f;

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
    // Matching results
    uint32_t matchFlag[kMaxClusters] = {};
    float    matchGEMx[kMaxClusters][4] = {};
    float    matchGEMy[kMaxClusters][4] = {};
    float    matchGEMz[kMaxClusters][4] = {};
    int      matchNum = 0; // number of clusters with matches (for quick access, can be derived from matchFlag)
    //for quick simple access to each matched hit on HC and GEM planes
    // HC_Energy, HC_x/y/z, GEM_x/y/z (in mm, beam center and target center coordinate)
    float    mHit_E[kMaxClusters] = {};
    float    mHit_x[kMaxClusters] = {};
    float    mHit_y[kMaxClusters] = {};
    float    mHit_z[kMaxClusters] = {};
    float    mHit_gx[kMaxClusters][2] = {};
    float    mHit_gy[kMaxClusters][2] = {};
    float    mHit_gz[kMaxClusters][2] = {};
    float    mHit_gid[kMaxClusters][2] = {}; //det_id for matched GEM hits

    // GEM reconstructed hits
    int        n_gem_hits = 0;
    uint8_t det_id[kMaxGemHits]       = {};
    float   gem_x[kMaxGemHits]        = {};
    float   gem_y[kMaxGemHits]        = {};
    float   gem_z[kMaxGemHits]        = {};
    float   gem_x_charge[kMaxGemHits] = {};
    float   gem_y_charge[kMaxGemHits] = {};
    float   gem_x_peak[kMaxGemHits]   = {};
    float   gem_y_peak[kMaxGemHits]   = {};
    uint8_t gem_x_size[kMaxGemHits]   = {};
    uint8_t gem_y_size[kMaxGemHits]   = {};
    uint8_t gem_x_mTbin[kMaxGemHits]   = {};
    uint8_t gem_y_mTbin[kMaxGemHits]   = {};

    //veto information
    int      veto_nch = 0;
    uint8_t veto_id[4]   = {}; // 0,1,2,3 for veto1-4
    int veto_npeaks[4] = {};
    float veto_peak_time[4][MAX_PEAKS]     = {};
    float veto_peak_integral[4][MAX_PEAKS] = {};

    //LMS reference PMT information
    int      lms_nch = 0;
    uint8_t lms_id[4]   = {}; // 0,1,2,3 for lms1-4
    int lms_npeaks[4] = {};
    float lms_peak_time[4][MAX_PEAKS]     = {};
    float lms_peak_integral[4][MAX_PEAKS] = {};

    // Raw 0xE10C SSP trigger bank words (one variable-length entry per event)
    std::vector<uint32_t> ssp_raw;
};

void setupReconBranches(TTree *tree, ReconEventData &ev){
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
    // Matching results
    tree->SetBranchAddress("matchFlag",    ev.matchFlag);
    tree->SetBranchAddress("matchGEMx",    ev.matchGEMx);
    tree->SetBranchAddress("matchGEMy",    ev.matchGEMy);
    tree->SetBranchAddress("matchGEMz",    ev.matchGEMz);
    tree->SetBranchAddress("match_num",     &ev.matchNum);
    // Quick and simple matching results for fast checks
    tree->SetBranchAddress("mHit_E", ev.mHit_E);
    tree->SetBranchAddress("mHit_x", ev.mHit_x);
    tree->SetBranchAddress("mHit_y", ev.mHit_y);
    tree->SetBranchAddress("mHit_z", ev.mHit_z);
    tree->SetBranchAddress("mHit_gx", ev.mHit_gx);
    tree->SetBranchAddress("mHit_gy", ev.mHit_gy);
    tree->SetBranchAddress("mHit_gz", ev.mHit_gz);
    tree->SetBranchAddress("mHit_gid", ev.mHit_gid);
    // GEM branches
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
    tree->SetBranchAddress("gem_x_mTbin",   ev.gem_x_mTbin);
    tree->SetBranchAddress("gem_y_mTbin",   ev.gem_y_mTbin);
};