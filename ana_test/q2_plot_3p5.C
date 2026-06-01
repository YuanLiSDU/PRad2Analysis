#include "../EventData.h"
#include "../PhysicsTools.h"

float livecharge[200] = {0.};

const int Nbins = 33;
const float binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

float gem_efficiency[Nbins] = {0.883383, 0.920705, 0.921194, 0.907423, 0.916365, 
    0.928948, 0.921558, 0.918408, 0.917756, 0.908999, 0.900851, 0.898294, 0.901366,
     0.884914, 0.881327, 0.881446, 0.875971, 0.864394, 0.850832, 0.855859, 0.867368,
      0.865487, 0.825393, 0.791823, 0.831812, 0.851094, 0.853007, 0.85037, 0.827107,
       0.822984, 0.852204, 0.840495, 0.822109};

float Ebeam = 3484.f; // MeV
float resolution = 0.033f;

bool useGEMs = true;

// Returns true if (x, y) [mm] lies inside the HyCal active acceptance:
// outside the beam hole (2.5 module widths) and inside the outer edge (16 module widths).
bool inHyCal(double xmm, double ymm) {
    const double module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
}

double GetEprime( double EBeam, double_t theta)//MeV
{
    /*Given the beam energy in MeV, and scattering angle theta in rad, calculate the
      expected scattering electron energy E'*/
    return EBeam * M_PROTON / ( EBeam * (1. - cos(theta)) + M_PROTON );
}
double GetQ2(double theta){//GeV^2
    return 2. * Ebeam * GetEprime(Ebeam, theta)*(1. - cos(theta))/1e6;
}
double GetMottCS(double theta){//in mb/sr
    double E_e = GetEprime(Ebeam, theta);
    double alpha = 1./137.;
    // Natural-unit result is in 1/MeV^2; convert to mb using (hbar*c)^2.
    const double hbarc_MeV_fm = 197.3269804;
    const double fm2_to_mb = 10.0; // 1 fm^2 = 10 mb
    double mott_nat = pow(alpha*cos(theta/2.)/(2*Ebeam*sin(theta/2.)*sin(theta/2.)), 2) * E_e / Ebeam;
    double mottcs = mott_nat * hbarc_MeV_fm * hbarc_MeV_fm * fm2_to_mb;
    return mottcs;
}
double GetTau(double theta){ //dimensionless
    double Q2 = GetQ2(theta);
    return Q2/(4.*pow(M_PROTON*1.e-3, 2));
}
double GetEpsilon(double theta){//unitless
    double tau = GetTau(theta);
    return 1./(1.+2.*(1.+tau)*pow(tan(theta/2.), 2));
}
double GetKellyGM(double Q2){
    const double_t a21_K = 0.12; // Magnetic form factor
    const double_t b21_K = 10.97;
    const double_t b22_K = 18.86;
    const double_t b23_K = 6.55;
    const double_t mu_p = 2.79284734463; // proton magnetic moment
    //return the Kelly magnetic form factor value at Q2 (GeV^2)
    double_t t = fabs(Q2) / (4. * (M_PROTON* M_PROTON/1.e6)); // Tau
    return mu_p*(1. + a21_K*t)/(1. + b21_K*t + b22_K*pow(t, 2) + b23_K*pow(t, 3));
}
Double_t GetKellyGE(double Q2)
{   
    const double a11_K = -0.24; // Electric form factor
    const double b11_K = 10.98;
    const double b12_K = 12.82;
    const double b13_K = 21.97;
    //return the Kelly electric form factor value at Q2 (GeV^2)
    double_t t = fabs(Q2) / (4. * (M_PROTON* M_PROTON/1.e6)); // Tau
    return (1. + a11_K*t)/(1. + b11_K*t + b12_K*pow(t,2) + b13_K*pow(t,3));
}
double GetGE(double theta, double CS_mb){
    double tau = GetTau(theta);
    double mottcs = GetMottCS(theta); // in mb/sr
    double eps = GetEpsilon(theta);
    double GM = GetKellyGM(GetQ2(theta));
    return sqrt((1.+tau)*CS_mb/mottcs-tau/eps*GM*GM);
}

// Calculate hydrogen areal density from cell pressure via ideal gas law.
//   pressure_mTorr : cell gas pressure [mTorr]
//   Returns        : H atom areal density [atoms/cm^2]
// Assumptions: H2 gas, T = 17.8 K, cylindrical cell L = 4 cm,
//              beam travels along the cylinder axis.
double GetArealDensity(double pressure_mTorr) {
    const double kB  = 1.380649e-16; // Boltzmann constant [erg/K = dyn·cm/K]
    const double T   = 19.0;         // cell temperature [K]
    const double L   = 4.0;          // cell length along beam axis [cm]
    // 1 mTorr = 1.33322 dyn/cm^2
    const double mTorr_to_dyncm2 = 1.33322;
    double P      = pressure_mTorr * mTorr_to_dyncm2; // [dyn/cm^2]
    double n_mol  = P / (kB * T);                      // H2 number density [molecules/cm^3]
    return 2.0 * n_mol * L;                            // H atoms/cm^2 (factor 2: H2 -> 2H)
}

bool Vetoed(float cl_time, float sci_time, float sci_int){
    // Simple veto logic: if the cluster time is within a certain window of the scintillator time, and the scintillator signal is above a threshold, we consider it a vetoed event.
    const float time_shift = 35.f; // ns
    const float time_window = 7.f; // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window) && (sci_int > int_threshold);
}

// ── Fill histograms for one or more input files ───────────────────────────
static void fillHists(const std::vector<TString> &fnames,
                      TH2F *hit_all, TH2F *E_angle,
                      TH2F *hits_mott, TH2F *hits_moller,
                      TH2F *E_angle_mott, TH2F *E_angle_moller,
                      TH1F *mott_yield, TH1F *moller_yield,
                      TH1F *veto_time, TH1F *veto_int,
                      TH1F *E_recon[Nbins],
                      float Ebeam)
{

    TChain *tree = new TChain("recon");
    for (const auto &fname : fnames)
        tree->Add(fname);
    std::cout << "Processing " << fnames.size() << " file(s), total entries: " << tree->GetEntries() << std::endl;

    ReconEventData ev;
    setupReconBranches(tree, ev);

    for (Long64_t i = 0; i < tree->GetEntries()/10; i++) {
        if (i % 10000 == 0)
            std::cout << "  " << i << " / " << tree->GetEntries() << "\r" << std::flush;
        tree->GetEntry(i);
        if(!useGEMs){
            for (int j = 0; j < ev.n_clusters; j++) {
                float x = ev.cl_x[j];
                float y = ev.cl_y[j];
                float z = ev.cl_z[j];
                float E = ev.cl_energy[j];
                float theta = atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;

                if( !inHyCal(x, y) ) continue;

                if( ev.cl_nblocks[j] < 4 ) continue; // remove single-block clusters which are mostly noise

                hit_all->Fill(x, y);
                E_angle->Fill(theta, E);

                if (isMott(E, Ebeam, resolution)) {
                    hits_mott->Fill(x, y);
                    E_angle_mott->Fill(theta, E);
                    mott_yield->Fill(theta);
                    break; // only consider the first Mott-like cluster in each event to avoid double counting
                }
            }
            if (ev.n_clusters == 2){
                float x1 = ev.cl_x[0], y1 = ev.cl_y[0], z1 = ev.cl_z[0], E1 = ev.cl_energy[0];
                float x2 = ev.cl_x[1], y2 = ev.cl_y[1], z2 = ev.cl_z[1], E2 = ev.cl_energy[1];
                if (!inHyCal(x1, y1) || !inHyCal(x2, y2)) continue;
                float theta1 = atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
                float theta2 = atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
                if (isMoller_kinematic(theta1, E1, theta2, E2, Ebeam, resolution)) {
                    if(fabs(GetMollerPhiDiff({DataPoint(x1, y1, z1, E1), DataPoint(x2, y2, z2, E2)})) > 10.f) continue;
                    hits_moller->Fill(x1, y1);
                    hits_moller->Fill(x2, y2);
                    E_angle_moller->Fill(theta1, E1);
                    E_angle_moller->Fill(theta2, E2);
                    moller_yield->Fill(theta1);
                    moller_yield->Fill(theta2);
                }
            }
        }else{
            for (int j = 0; j < ev.matchNum; j++) {
                GEMHit hit;
                hit.x = ev.mHit_gx[j][1];
                hit.y = ev.mHit_gy[j][1];
                hit.z = ev.mHit_gz[j][1];
                projectToHyCal(hit, 6270.f);
                float E = ev.mHit_E[j];
                float theta = atan2(std::sqrt(hit.x*hit.x + hit.y*hit.y), hit.z) * 180.f / M_PI;

                if( !inHyCal(hit.x, hit.y) ) continue;

                hit_all->Fill(hit.x, hit.y);
                E_angle->Fill(theta, E);

                if (isMott(E, Ebeam, resolution)) {
                    hits_mott->Fill(hit.x, hit.y);
                    E_angle_mott->Fill(theta, E);
                    mott_yield->Fill(theta);
                    int bin = mott_yield->FindBin(theta);
                    if (bin >= 1 && bin <= Nbins)
                        E_recon[bin-1]->Fill(E);
                }
            }
        }
    }
    std::cout << std::endl;

    delete tree;
};

void q2_plot_3p5(const char *files = "../../A/24917_recon_filter.root",
                double lc = 4639600./10.){

    double lc_B = 2292700;

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

    double density = GetArealDensity(1032.); // atoms/cm^2 per mTorr
    std::cout << "Hydrogen areal density = " << density << " atoms/cm^2" << std::endl;
    double luminosity = density * lc * 1e-9 / 1.6e-19; // livecharge in nC -> C -> N_beam; Luminosity [cm^-2]
    std::cout << "Luminosity = " << luminosity << " cm^-2" << std::endl;

    // Allocate histograms
    TH2F *hit_all       = new TH2F("hit_all",        "Hit Position Distribution(all clusters);X (mm);Y (mm)",                    600, -360, 360, 600, -360, 360);
    TH2F *E_angle       = new TH2F("E_angle",         "E vs Scattering Angle;Scattering Angle (deg);Energy (MeV)",                120, 0, 6, 5000, 0, 5000);
    TH2F *hits_mott     = new TH2F("hits_mott",       "Hit Position Distribution (e-p);X (mm);Y (mm)",                           600, -360, 360, 600, -360, 360);
    TH2F *hits_moller     = new TH2F("hits_moller",     "Hit Position Distribution (e-e);X (mm);Y (mm)",                           600, -360, 360, 600, -360, 360);
    TH2F *E_angle_mott  = new TH2F("E_angle_mott",    "E vs Scattering Angle (e-p);Scattering Angle (deg);Energy (MeV)",         120, 0, 6, 5000, 0, 5000);
    TH2F *E_angle_moller = new TH2F("E_angle_moller", "E vs Scattering Angle (e-e);Scattering Angle (deg);Energy (MeV)",         120, 0, 6, 5000, 0, 5000);
    TH1F *mott_yield    = new TH1F("mott_yield",      "e-p Yield;Theta (deg);Yield(arbitrary units)", Nbins, binEdge);
    TH1F *moller_yield  = new TH1F("moller_yield",    "2 arm Moller Yield;Theta (deg);Yield(arbitrary units)",  Nbins, binEdge);
    TH1F *yield_ratio   = new TH1F("yield_ratio",     "e-p/Moller Yield Ratio;Theta (deg);Yield Ratio",   Nbins, binEdge);
    TH1F *veto_time = new TH1F("veto_time", "Cluster Time - Sci Time;Time Difference (ns);Counts", 100, -50, 50);
    TH1F *veto_int = new TH1F("veto_int", "Scintillator Integral;Integral (arb. units);Counts", 500, 0, 10000);
    TH1F *E_recon[Nbins];
    for (int i = 0; i < Nbins; i++) {
        E_recon[i] = new TH1F(Form("E_recon_bin%d", i), Form("Reconstructed Energy (bin %d);Energy (MeV);Counts", i), 350, 0, Ebeam*1.1);
    }

    TH2F *hit_all_B       = new TH2F("hit_all_B",        "Hit Position Distribution(all clusters, type B);X (mm);Y (mm)",                    600, -360, 360, 600, -360, 360);
    TH2F *E_angle_B       = new TH2F("E_angle_B",         "E vs Scattering Angle (type B);Scattering Angle (deg);Energy (MeV)",                120, 0, 6, 5000, 0, 5000);
    TH2F *hits_mott_B     = new TH2F("hits_mott_B",       "Hit Position Distribution (e-p, type B);X (mm);Y (mm)",                           600, -360, 360, 600, -360, 360);
    TH2F *hits_moller_B     = new TH2F("hits_moller_B",     "Hit Position Distribution (e-e, type B);X (mm);Y (mm)",                           600, -360, 360, 600, -360, 360);
    TH2F *E_angle_mott_B  = new TH2F("E_angle_mott_B",    "E vs Scattering Angle (e-p, type B);Scattering Angle (deg);Energy (MeV)",         120, 0, 6, 5000, 0, 5000);
    TH2F *E_angle_moller_B = new TH2F("E_angle_moller_B", "E vs Scattering Angle (e-e, type B);Scattering Angle (deg);Energy (MeV)",         120, 0, 6, 5000, 0, 5000);
    TH1F *mott_yield_B = new TH1F("mott_yield_B",      "e-p Yield (type B);Theta (deg);Yield(arbitrary units)", Nbins, binEdge);
    TH1F *moller_yield_B = new TH1F("moller_yield_B",    "2 arm Moller Yield (type B);Theta (deg);Yield(arbitrary units)",  Nbins, binEdge);
    TH1F *yield_ratio_B = new TH1F("yield_ratio_B",     "e-p/Moller Yield Ratio (type B);Theta (deg);Yield Ratio",   Nbins, binEdge);
    TH1F *veto_time_B = new TH1F("veto_time_B", "Cluster Time - Sci Time (type B);Time Difference (ns);Counts", 100, -50, 50);
    TH1F *veto_int_B = new TH1F("veto_int_B", "Scintillator Integral (type B);Integral (arb. units);Counts", 500, 0, 10000);
    TH1F *E_recon_B[Nbins];
    for (int i = 0; i < Nbins; i++) {
        E_recon_B[i] = new TH1F(Form("E_recon_B_bin%d", i), Form("Reconstructed Energy (type B, bin %d);Energy (MeV);Counts", i), 350, 0, Ebeam*1.1);
    }

    fillHists(fileList, hit_all, E_angle, hits_mott, hits_moller, E_angle_mott, 
        E_angle_moller, mott_yield, moller_yield, veto_time, veto_int, E_recon, Ebeam);
    std::vector<TString> fileListB = {"../../B/24916_recon_filter.root"};
    fillHists(fileListB, hit_all_B, E_angle_B, hits_mott_B, hits_moller_B, E_angle_mott_B, 
        E_angle_moller_B, mott_yield_B, moller_yield_B, veto_time_B, veto_int_B, E_recon_B, Ebeam);

    // Scale yields by 1/livecharge so y-axis is yield per unit charge
    mott_yield  ->Scale(1.0 / lc);
    moller_yield ->Scale(1.0 / lc);
    mott_yield_B->Scale(1.0 / lc_B);
    moller_yield_B->Scale(1.0 / lc_B);
        // Scale() already propagates bin errors correctly (sigma = sqrt(N)/lc).
        // No manual SetBinError needed.

    mott_yield  ->GetYaxis()->SetTitle("Yield / livecharge");
    moller_yield ->GetYaxis()->SetTitle("Yield / livecharge");
    mott_yield_B->GetYaxis()->SetTitle("Yield / livecharge");
    moller_yield_B->GetYaxis()->SetTitle("Yield / livecharge");

    // Draw: overview canvas
    TCanvas *cv1 = new TCanvas("c_overview", "Background Analysis", 1200, 800);
    cv1->Divide(3, 2);
    cv1->cd(1); hit_all      ->Draw("COLZ");
    cv1->cd(2); E_angle      ->Draw("COLZ");
    cv1->cd(3); hits_mott    ->Draw("COLZ");
    cv1->cd(4); E_angle_mott ->Draw("COLZ");
    cv1->cd(5); hits_moller ->Draw("COLZ");
    cv1->cd(6); E_angle_moller->Draw("COLZ");
    cv1->SaveAs("overview.png");
    
    //calculate corresponding Q2 values for angle bin edges
    double q2_edges[Nbins+1];
    for (int i = 0; i <= Nbins; i++) {
        q2_edges[i] = GetQ2(binEdge[i] * M_PI / 180.f); // bin edges in Q2 space
    }

    //Generator reference — read generator ntuple, same cuts, normalize by L_gen
    //calculate acceptance correction for each bin: acceptance = N_accepted / N_generated in that bin
    double geo_acceptance[Nbins];
    TH1F *gen_cross = nullptr;
    TH1F *acceptance_q2 = nullptr;
    {
        const double L_gen_pb  = 2.2018e-02;       // generator integrated luminosity [pb^-1]
        const double L_gen_cm2 = L_gen_pb * 1e36;  // [cm^-2]  (1 pb = 1e-36 cm^2)
        const double cm2_to_mb = 1e27;
        const double theta_min = binEdge[0]     * M_PI / 180.;  // rad
        const double theta_max = binEdge[Nbins] * M_PI / 180.;  // rad

        TH1F *gen_raw = new TH1F("gen_raw", "", Nbins, binEdge); // theta in degrees
        double in_acceptance[Nbins] = {}, out_acceptance[Nbins] = {};
        {   

            TFile *gf = TFile::Open("../../../3.5GeV_e-.root");
            if (!gf || gf->IsZombie()) {
                std::cerr << "ERROR: cannot open 3.5GeV_e-.root" << std::endl;
                return;
            }
            TNtuple *nt = (TNtuple*)gf->Get("ntp");
            if (!nt) {
                std::cerr << "ERROR: ntuple 'ntp' not found in 3.5GeV_e-.root" << std::endl;
                gf->Close(); return;
            }
            Float_t E_l, theta_l, phi_l;
            nt->SetBranchAddress("E_l",     &E_l);
            nt->SetBranchAddress("theta_l", &theta_l);
            nt->SetBranchAddress("phi_l",   &phi_l);
            for (Long64_t i = 0; i < nt->GetEntries(); i++) {
                nt->GetEntry(i);
                if (theta_l < theta_min || theta_l > theta_max) continue;
                E_l += gRandom->Gaus(0., Ebeam * resolution / sqrt(Ebeam/1000.)); // apply energy smearing
                if (!isMott(E_l, Ebeam, resolution))            continue;
                gen_raw->Fill(theta_l * 180. / M_PI);  // store in degrees

                double x = 6270.f * tan(theta_l) * cos(phi_l);
                double y = 6270.f * tan(theta_l) * sin(phi_l);
                float theta = theta_l * 180.f / M_PI;
                if( inHyCal(x, y) && theta > 0.3f )
                    in_acceptance[gen_raw->FindBin(theta_l * 180. / M_PI) - 1]++;
                else out_acceptance[gen_raw->FindBin(theta_l * 180. / M_PI) - 1]++;
            }
            gf->Close();
        }
        std::cout << "Generator: " << gen_raw->GetEntries() << " events pass cuts" << std::endl;

        gen_cross = new TH1F("gen_cross", "", Nbins, q2_edges);
        for (int i = 1; i <= gen_raw->GetNbinsX(); i++) {
            double theta_low  = binEdge[i-1] * M_PI / 180.;
            double theta_high = binEdge[i]   * M_PI / 180.;
            double dOmega = 2.*M_PI*(cos(theta_low) - cos(theta_high));
            double N = gen_raw->GetBinContent(i);
            if (dOmega > 0.) {
                gen_cross->SetBinContent(i, N / (L_gen_cm2 * dOmega) * cm2_to_mb);
                gen_cross->SetBinError  (i, sqrt(N) / (L_gen_cm2 * dOmega) * cm2_to_mb);
            }
        }
        delete gen_raw;
        
        acceptance_q2 = new TH1F("acceptance_q2", "Detector Acceptance;Q^{2} (GeV^{2});Acceptance(%)", Nbins, q2_edges);
        //calculate geometric acceptance for each bin
        for (int i = 1; i <= gen_cross->GetNbinsX(); i++) {
            double theta_low  = binEdge[i-1] * M_PI / 180.;
            double theta_high = binEdge[i]   * M_PI / 180.;
            double dOmega = 2.*M_PI*(cos(theta_low) - cos(theta_high));
            double total_gen = in_acceptance[i-1] + out_acceptance[i-1];
            if (dOmega > 0.) {
                geo_acceptance[i-1] = (total_gen > 0.) ? in_acceptance[i-1] / total_gen : 0.;
                acceptance_q2->SetBinContent(i, geo_acceptance[i-1] * 100.);
                std::cout << "Bin " << i << ": theta [" << binEdge[i-1] << ", " << binEdge[i] << "] deg, "
                          << "Q2 [" << GetQ2(theta_high) << ", " << GetQ2(theta_low) << "] GeV^2, "
                          << "Acceptance = " << geo_acceptance[i-1] * 100. << " %" << std::endl;
            }
        }
    }

    // Calculate cross section from yield: dσ/dΩ = yield / (luminosity × dΩ)
    TH1F *cross_section = new TH1F("cross_section", "Cross Section;Q^{2} (GeV^{2});d#sigma/d#Omega (mb/sr)", Nbins, q2_edges);
    for (int i = 1; i <= mott_yield->GetNbinsX(); i++) {
        // solid angle of this theta bin [sr]
        double theta_low  = binEdge[i-1] * M_PI / 180.;
        double theta_high = binEdge[i]   * M_PI / 180.;
        double dOmega = 2. * M_PI * (cos(theta_low) - cos(theta_high));

        double yield = mott_yield->GetBinContent(i);  // N_i / lc  [/nC]
        yield -= mott_yield_B->GetBinContent(i); // subtract type B background

        //geo acceptance correction
        yield /= geo_acceptance[i-1];
        //gem efficiency correction
        yield /= gem_efficiency[i-1];

        double err   = sqrt(pow(mott_yield->GetBinError(i), 2) + pow(mott_yield_B->GetBinError(i), 2)); // combine errors in quadrature
        err /= geo_acceptance[i-1]; // propagate acceptance correction to error
        err /= gem_efficiency[i-1]; // propagate efficiency correction to error

        // dσ/dΩ = yield × (e/nC→1) / (density × dΩ),  e = 1.6e-19 C, 1 nC = 1e-9 C
        // lc cancels:  N_i/(L·dΩ) = (yield·lc) / (density·lc/1.6e-10 · dΩ)
        //            = yield · 1.6e-10 / (density · dΩ)   [cm²/sr]
        // × 1e27 → mb/sr
        const double e_over_nC  = 1.6e-19 / 1e-9;  // = 1.6e-10
        const double cm2_to_mb  = 1e27;
        if (dOmega > 0.) {
            double dsig = yield * e_over_nC / (density * dOmega) * cm2_to_mb;
            double derr = err   * e_over_nC / (density * dOmega) * cm2_to_mb;
            cross_section->SetBinContent(i, dsig);
            cross_section->SetBinError(i, derr);
        }
    }

    TCanvas *cv2 = new TCanvas("c_mott", "Mott Yield", 900, 700);
    cv2->SetGrid();
    mott_yield->Draw("EP");

    TCanvas *cv2_moller = new TCanvas("c_moller", "Moller Yield", 900, 700);
    cv2_moller->SetGrid();
    moller_yield->Draw("EP");

    TCanvas *cv2_ratio = new TCanvas("c_ratio", "Yield Ratio", 900, 700);
    cv2_ratio->SetGrid();
    yield_ratio->Draw("EP");

    // Build relative-error histogram (percentage)
    TH1F *rel_err = new TH1F("rel_err", ";Q^{2} (GeV^{2});Statistical Error (%)", Nbins, q2_edges);
    for (int i = 1; i <= cross_section->GetNbinsX(); i++) {
        double val = cross_section->GetBinContent(i);
        double err = cross_section->GetBinError(i);
        if (val != 0.) {
            rel_err->SetBinContent(i, err / val * 100.);
            rel_err->SetBinError(i, 0.);
        }
    }

    TCanvas *cv3 = new TCanvas("c_cross", "Mott Cross Section", 1400, 900);
    // Upper pad: cross section
    TPad *pad_top = new TPad("pad_top", "", 0., 0.3, 1., 1.);
    pad_top->SetBottomMargin(0.02);
    pad_top->SetGrid();
    pad_top->SetLogy();
    pad_top->SetLogx();
    pad_top->Draw();
    pad_top->cd();
    cross_section->SetStats(0);
    cross_section->GetXaxis()->SetLabelSize(0.);
    cross_section->GetXaxis()->SetTitleSize(0.);
    cross_section->SetMarkerStyle(20);
    cross_section->SetMarkerSize(1.0);
    cross_section->SetTitle("e-p Cross Section @ 3.5GeV (corrected for acceptance)");
    cross_section->Draw("E1P");

    // overlay generator reference
    gen_cross->SetStats(0);
    gen_cross->SetLineColor(kRed + 1);
    gen_cross->SetMarkerStyle(24);
    gen_cross->SetMarkerColor(kRed + 1);
    gen_cross->Draw("E1P SAME");

    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(cross_section, "Data (measured lumi)", "lp");
    leg->AddEntry(gen_cross,     "Generator (L_{gen})",  "lp");
    leg->Draw();

    cv3->cd();
    // Lower pad: relative error
    TPad *pad_bot = new TPad("pad_bot", "", 0., 0., 1., 0.3);
    pad_bot->SetTopMargin(0.02);
    pad_bot->SetBottomMargin(0.35);
    pad_bot->SetGrid();
    pad_bot->SetLogx();
    pad_bot->SetLogy();
    pad_bot->SetTicky();
    pad_bot->Draw();
    pad_bot->cd();
    rel_err->SetStats(0);
    rel_err->GetXaxis()->SetLabelSize(0.10);
    rel_err->GetXaxis()->SetTitleSize(0.12);
    rel_err->GetYaxis()->SetLabelSize(0.09);
    rel_err->GetYaxis()->SetTitleSize(0.10);
    rel_err->GetYaxis()->SetTitleOffset(0.5);
    rel_err->SetFillColor(kAzure - 9);
    rel_err->SetLineColor(kBlue + 1);
    rel_err->Draw("BAR");

    cv3->SaveAs("cross_section.png");

    TCanvas *cv4 = new TCanvas("c_acceptance", "Geometric Acceptance", 900, 700);
    cv4->SetGrid();
    cv4->SetLogx();
    acceptance_q2->SetStats(0);
    acceptance_q2->SetMarkerStyle(20);
    acceptance_q2->SetMarkerSize(1.0);
    acceptance_q2->SetTitle("Detector Geometric Acceptance");
    acceptance_q2->GetYaxis()->SetRangeUser(0., 105);
    acceptance_q2->Draw("P");
    cv4->SaveAs("acceptance.png");

    TCanvas *cv5 = new TCanvas("c_yield", "Mott Yield Comparison", 900, 700);
    cv5->SetGrid();
    cv5->SetLogx();

    TH1F *bg_ratio = (TH1F *)mott_yield_B->Clone("bg_ratio");
    bg_ratio->SetTitle("Type-B / Type-A Mott Yield;#theta (deg);mott_{yield_B} / mott_{yield_A}");
    bg_ratio->SetStats(0);

    // Protect against division by zero while keeping proper error propagation.
    for (int i = 1; i <= bg_ratio->GetNbinsX(); i++) {
        double a = mott_yield->GetBinContent(i);
        if (a <= 0.) {
            bg_ratio->SetBinContent(i, 0.);
            bg_ratio->SetBinError(i, 0.);
        }
    }
    bg_ratio->Divide(mott_yield);
    bg_ratio->SetMarkerStyle(20);
    bg_ratio->SetMarkerSize(1.0);
    bg_ratio->SetMarkerColor(kBlue + 1);
    bg_ratio->SetLineColor(kBlue + 1);
    bg_ratio->Draw("E1P");
    cv5->SaveAs("bg_ratio_B_over_A.png");

    //calculate GE and plot
    TH1F *GE = new TH1F("GE", "Electric Form Factor;Q^{2} (GeV^{2});G_{E}", Nbins, q2_edges);
    TH1F *GE_gen = new TH1F("GE_gen", "Kelly Electric Form Factor;Q^{2} (GeV^{2});G_{E}", Nbins, q2_edges);
    TH1F *GE_theory = new TH1F("GE_theory", "Kelly Electric Form Factor;Q^{2} (GeV^{2});G_{E}", Nbins, q2_edges);
    for (int i = 1; i <= GE_theory->GetNbinsX(); i++) {
        double theta = (binEdge[i-1] + binEdge[i]) / 2. * M_PI / 180.; // bin center in radians
        double Q2 = GetQ2(theta);
        GE_theory->SetBinContent(i, GetKellyGE(Q2));
    }
    for(int i =1; i <= gen_cross->GetNbinsX(); i++){
        double theta = (binEdge[i-1] + binEdge[i]) / 2. * M_PI / 180.; // bin center in radians
       // if(theta < 0.7f * M_PI / 180.) continue;
        double gen_cs = gen_cross->GetBinContent(i);
        double GM  = GetKellyGM(GetQ2(theta));
        double eps = GetEpsilon(theta);
        double tau = GetTau(theta);
        double mott = GetMottCS(theta);
        double arg = (1. + tau) * gen_cs / mott - tau / eps * GM * GM;
        if (gen_cs <= 0. || arg <= 0.) { GE_gen->SetBinContent(i, 0.); GE_gen->SetBinError(i, 0.); continue; }
        GE_gen->SetBinContent(i, sqrt(arg));
        double dCS = gen_cross->GetBinError(i);
        double dGE_dCS = 0.5 * (1. + tau) / sqrt(arg) / mott;
        GE_gen->SetBinError(i, dCS * fabs(dGE_dCS));
    }
    for (int i = 1; i <= cross_section->GetNbinsX(); i++) {
        double theta = (binEdge[i-1] + binEdge[i]) / 2. * M_PI / 180.; // bin center in radians
        double CS_mb = cross_section->GetBinContent(i); // in mb/sr
        double GM  = GetKellyGM(GetQ2(theta));
        double eps = GetEpsilon(theta);
        double tau = GetTau(theta);
        double mott = GetMottCS(theta);
        double arg = (1. + tau) * CS_mb / mott - tau / eps * GM * GM;
        if (CS_mb <= 0. || arg <= 0.) { GE->SetBinContent(i, 0.); GE->SetBinError(i, 0.); continue; }
        GE->SetBinContent(i, sqrt(arg));
        double dCS = cross_section->GetBinError(i);
        double dGE_dCS = 0.5 * (1. + tau) / sqrt(arg) / mott;
        GE->SetBinError(i, dCS * fabs(dGE_dCS));
    }
    TCanvas *cv6 = new TCanvas("c_GE", "Electric Form Factor G_E", 900, 700);
    cv6->SetGrid();
    cv6->SetLogx();
    GE->SetStats(0);
    GE->SetMarkerStyle(20);
    GE->SetMarkerSize(1.0);
    GE->SetTitle("Electric Form Factor G_{E}");
    GE->Draw("E1P");
    GE_gen->SetStats(0);
    GE_gen->SetLineColor(kRed + 1);
    GE_gen->SetMarkerStyle(24);
    GE_gen->SetMarkerColor(kRed + 1);
    GE_gen->Draw("P SAME");
    GE_theory->SetLineColor(kGreen + 2);
    GE_theory->SetLineStyle(2);
    GE_theory->Draw("L SAME");
   
    TLegend *leg_GE = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg_GE->SetBorderSize(0);
    leg_GE->AddEntry(GE, "Data (measured lumi)", "lp");
    leg_GE->AddEntry(GE_gen, "Kelly Parametrization", "lp");
    leg_GE->AddEntry(GE_theory, "Kelly Parametrization (theory)", "l");
    leg_GE->Draw();
    cv6->SaveAs("GE.png");

    //fit E_recon in each bin to check energy reconstruction, save to ROOT file
    TFile *fout_erecon = TFile::Open("E_recon_fits.root", "RECREATE");
    for (int i = 0; i < Nbins; i++) {
        if (E_recon[i]->GetEntries() < 100) continue; // skip bins with too few entries for a meaningful fit
        double sigma_hint = Ebeam * resolution / sqrt(Ebeam / 1000.);
        TF1 *fit = new TF1(Form("fit_bin%d", i), "gaus", Ebeam - 1*sigma_hint, Ebeam + 1*sigma_hint);
        // gaus parameters: [0]=amplitude, [1]=mean, [2]=sigma
        fit->SetParameters(E_recon[i]->GetMaximum(), Ebeam, sigma_hint);
        fit->SetLineColor(kRed + 1);
        fit->SetLineWidth(2);
        E_recon[i]->Fit(fit, "RQ"); // R=range, Q=quiet (no delete: TF1 owned by histogram)

        double mean      = fit->GetParameter(1);
        double sigma     = fit->GetParameter(2);
        double mean_err  = fit->GetParError(1);
        double sigma_err = fit->GetParError(2);
        std::cout << "Bin " << i << " [" << binEdge[i] << "-" << binEdge[i+1] << " deg]"
                  << ": mean = " << mean << " +/- " << mean_err << " MeV"
                  << ", sigma = " << sigma << " +/- " << sigma_err << " MeV"
                  << ", sigma/E = " << sigma/mean*100. << "%" << std::endl;

        TCanvas *ctmp = new TCanvas(Form("c_erecon_bin%d", i),
                                    Form("E_recon [%.3f-%.3f deg]", binEdge[i], binEdge[i+1]),
                                    700, 500);
        ctmp->SetGrid();
        E_recon[i]->SetStats(0);
        E_recon[i]->Draw();  // fit line drawn automatically (still in histogram's function list)

        TPaveText *pt = new TPaveText(0.13, 0.62, 0.52, 0.88, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(1);
        pt->SetTextAlign(12);
        pt->AddText(Form("#theta #in [%.3f, %.3f] deg", binEdge[i], binEdge[i+1]));
        pt->AddText(Form("Mean = %.2f #pm %.2f MeV", mean, mean_err));
        pt->AddText(Form("#sigma = %.2f #pm %.2f MeV", sigma, sigma_err));
        pt->AddText(Form("#sigma/E = %.3f%%", sigma / mean * 100.));
        pt->Draw();

        fout_erecon->cd();
        ctmp->Write();
    }
    fout_erecon->Close();
    std::cout << "Saved E_recon fits to E_recon_fits.root" << std::endl;


