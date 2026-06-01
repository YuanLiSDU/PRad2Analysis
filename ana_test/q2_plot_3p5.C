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

float Ebeam = 3484.f; // MeV
float resolution = 0.03f;

bool useGEMs = true;

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

                if( fabs(x) < 20.75 * 2.5 && fabs(y) < 20.75 * 2.5 ) continue;
                if( fabs(x) > 20.75 * 16. || fabs(y) > 20.75 * 16. ) continue;

                if( ev.cl_nblocks[j] < 5 ) continue; // remove single-block clusters which are mostly noise

                hit_all->Fill(x, y);
                E_angle->Fill(theta, E);

                if (isMott(E, Ebeam, resolution)) {
                    //find veto info for this cluster
                    float cl_time = ev.cl_time[j];
                    float sci_time = -999.f;
                    float sci_int = -999.f;
                    bool veto = false;
                    //if(theta < 1.1f){
                    /*for(int k = 0; k < ev.veto_nch; k++){
                        for(int p = 0; p < ev.veto_npeaks[k]; p++){
                            sci_time = ev.veto_peak_time[k][p];
                            sci_int = ev.veto_peak_integral[k][p];
                            if(cl_time-sci_time < 28.f || cl_time-sci_time > 42.f) continue;
                            veto_time->Fill(cl_time - sci_time);
                            veto_int->Fill(sci_int);
                            if(Vetoed(cl_time, sci_time, sci_int)){
                                veto = true;
                                break;
                            }
                        }
                    }*/
                    //}
                    //if(theta < 1.1f && ((fabs(x) > 20.75 * 2.1 && fabs(y) > 20.) || (fabs(y) > 20.75 * 2.1 && fabs(x) > 20.))) continue; // apply veto for low-theta clusters near the center
                    //if(theta < 1.1f && veto) break; // apply veto for low-theta clusters
                    //if(theta < 0.5f) break; 
                    hits_mott->Fill(x, y);
                    E_angle_mott->Fill(theta, E);
                    mott_yield->Fill(theta);
                    break; // only consider the first Mott-like cluster in each event to avoid double counting
                }
            }
            if (ev.n_clusters == 2){
                float x1 = ev.cl_x[0], y1 = ev.cl_y[0], z1 = ev.cl_z[0], E1 = ev.cl_energy[0];
                float x2 = ev.cl_x[1], y2 = ev.cl_y[1], z2 = ev.cl_z[1], E2 = ev.cl_energy[1];
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

                if( fabs(hit.x) < 20.75 * 2.5 && fabs(hit.y) < 20.75 * 2.5 ) continue;
                if( fabs(hit.x) > 20.75 * 16. || fabs(hit.y) > 20.75 * 16. ) continue;

                //if( ev.cl_nblocks[j] < 3 ) continue; // remove single-block clusters which are mostly noise

                hit_all->Fill(hit.x, hit.y);
                E_angle->Fill(theta, E);

                if (isMott(E, Ebeam, resolution)) {
                    hits_mott->Fill(hit.x, hit.y);
                    E_angle_mott->Fill(theta, E);
                    mott_yield->Fill(theta, 1./0.85);
                }
            }
        }
    }
    std::cout << std::endl;

    delete tree;
};

void q2_plot(const char *files = "../../A/24917_recon_filter.root",
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

    fillHists(fileList, hit_all, E_angle, hits_mott, hits_moller, E_angle_mott, 
        E_angle_moller, mott_yield, moller_yield, veto_time, veto_int, Ebeam);
    std::vector<TString> fileListB = {"../../B/24916_recon_filter.root"};
    fillHists(fileListB, hit_all_B, E_angle_B, hits_mott_B, hits_moller_B, E_angle_mott_B, 
        E_angle_moller_B, mott_yield_B, moller_yield_B, veto_time_B, veto_int_B, Ebeam);

    TCanvas *cv_veto = new TCanvas("c_veto", "Veto Distributions", 1200, 500);
    cv_veto->Divide(2, 2);
    cv_veto->cd(1); veto_time->Draw();
    cv_veto->cd(2); veto_int->Draw();
    cv_veto->cd(3); veto_time_B->Draw();
    cv_veto->cd(4); veto_int_B->Draw();
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

    double q2_edges[Nbins+1];
    double geo_acceptance[Nbins];
    for (int i = 0; i <= Nbins; i++) {
        q2_edges[i] = GetQ2(binEdge[i] * M_PI / 180.f); // bin edges in Q2 space
    }

    TH1F *cross_section = new TH1F("cross_section", "Cross Section;Q^{2} (GeV^{2});d#sigma/d#Omega (mb/sr)", Nbins, q2_edges);
    for (int i = 1; i <= mott_yield->GetNbinsX(); i++) {
        // solid angle of this theta bin [sr]
        double theta_low  = binEdge[i-1] * M_PI / 180.;
        double theta_high = binEdge[i]   * M_PI / 180.;
        double dOmega = 2. * M_PI * (cos(theta_low) - cos(theta_high));

        double yield = mott_yield->GetBinContent(i);  // N_i / lc  [/nC]
        yield -= mott_yield_B->GetBinContent(i); // subtract type B background
        double err   = sqrt(pow(mott_yield->GetBinError(i), 2) + pow(mott_yield_B->GetBinError(i), 2)); // combine errors in quadrature

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

    //Generator reference — read generator ntuple, same cuts, normalize by L_gen
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
            TNtuple *nt = (TNtuple*)gf->Get("ntp");
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
                bool cuts = false;
                float theta = theta_l * 180.f / M_PI;
                //if(theta < 1.1f && ((fabs(x) > 20.75 * 2.1 && fabs(y) > 20.) || (fabs(y) > 20.75 * 2.1 && fabs(x) > 20.))) cuts = true;
                if( (fabs(x) > 20.75 * 2.5 || fabs(y) > 20.75 * 2.5)
                    && (fabs(x) < 20.75 * 16. && fabs(y) < 20.75 * 16.) && theta > 0.3f)
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
            if (dOmega > 0.) {
                geo_acceptance[i-1] = (in_acceptance[i-1]) / (in_acceptance[i-1] + out_acceptance[i-1]);
                acceptance_q2->SetBinContent(i, geo_acceptance[i-1] * 100.);
                std::cout << "Bin " << i << ": theta [" << binEdge[i-1] << ", " << binEdge[i] << "] deg, "
                          << "Q2 [" << GetQ2(theta_high) << ", " << GetQ2(theta_low) << "] GeV^2, "
                          << "Acceptance = " << geo_acceptance[i-1] * 100. << " %" << std::endl;
            }
        }

    }

    // Moller generator reference — read text file (9 columns per line):
    //   E1  theta1  phi1  E2  theta2  phi2  gamaE  gamaTheta  gamaPhi
    //   (energies MeV, angles rad)
    TH1F *gen_moller_cs       = new TH1F("gen_moller_cs",       "Moller dσ/dΩ (gen);Q^{2} (GeV^{2});d#sigma/d#Omega (mb/sr)", Nbins, q2_edges);
    TH1F *double_arm_geo_acc  = new TH1F("double_arm_geo_acc",  "Moller Double-Arm Geometric Acceptance;Q^{2} (GeV^{2});Acceptance (%)", Nbins, q2_edges);
    {
        const double L_gen_nb  = 0.440658;            // [nb^-1]
        const double L_gen_cm2 = L_gen_nb * 1e33;     // 1 nb = 1e-33 cm^2  ⇒  nb^-1 → cm^-2
        const double cm2_to_mb = 1e27;
        const double theta_min = 1.45     * M_PI / 180.;  // rad
        const double theta_max = binEdge[Nbins] * M_PI / 180.;  // rad

        const char *fname_moller =
            "../data/0.7GeV/ee_pass6_728.90MeV_5000000_theta_0.5_5.0_COMIN_50_COMAX_650_total_e-._0.440658nb-1dat";
        std::ifstream fin(fname_moller);
        if (!fin.is_open()) {
            std::cerr << "ERROR: cannot open " << fname_moller << std::endl;
        }

        TH1F *gen_moller_raw = new TH1F("gen_moller_raw", "", Nbins, q2_edges);
        double m_in [Nbins] = {};
        double m_out[Nbins] = {};

        auto inHyCal = [](double xmm, double ymm) {
            return (fabs(xmm) > 20.75 * 2.5 || fabs(ymm) > 20.75 * 2.5)
                && (fabs(xmm) < 20.75 * 16. && fabs(ymm) < 20.75 * 16.);
        };

        float E1, theta1, phi1, E2, theta2, phi2, gE, gTh, gPh;
        Long64_t nread = 0;
        while (fin >> E1 >> theta1 >> phi1 >> E2 >> theta2 >> phi2 >> gE >> gTh >> gPh) {
            nread++;
            if (theta1 < theta_min || theta1 > theta_max) continue;
            if (theta2 < theta_min || theta2 > theta_max) continue;

            float E1s = E1 + gRandom->Gaus(0., E1 * resolution / sqrt(E1/1000.));
            float E2s = E2 + gRandom->Gaus(0., E2 * resolution / sqrt(E2/1000.));
            if (!isMoller_kinematic(theta1 * 180./M_PI, E1s,
                                     theta2 * 180./M_PI, E2s, Ebeam, resolution)) continue;
            
            // Project to HyCal plane (z = 6270 mm)
            double x1 = 6270. * tan(theta1) * cos(phi1);
            double y1 = 6270. * tan(theta1) * sin(phi1);
            double x2 = 6270. * tan(theta2) * cos(phi2);
            double y2 = 6270. * tan(theta2) * sin(phi2);
            bool both_in = inHyCal(x1, y1) && inHyCal(x2, y2);

            // Fill cross-section histogram for both arms (per-arm dσ/dΩ)
            double q2_1 = GetQ2(theta1);
            double q2_2 = GetQ2(theta2);
            gen_moller_raw->Fill(q2_1);
            gen_moller_raw->Fill(q2_2);

            int b1 = gen_moller_raw->FindBin(q2_1);
            int b2 = gen_moller_raw->FindBin(q2_2);
            if (b1 >= 1 && b1 <= Nbins) {
                if(both_in) m_in[b1-1]++;
                else        m_out[b1-1]++;
            }
            if (b2 >= 1 && b2 <= Nbins) {
                if(both_in) m_in[b2-1]++;
                else        m_out[b2-1]++;
            }
        }
        fin.close();
        std::cout << "Moller generator: read " << nread
                  << " events, " << gen_moller_raw->GetEntries()/2.
                  << " pass cuts" << std::endl;

        // dσ/dΩ per bin (single-arm, summed over both electrons)
        for (int i = 1; i <= Nbins; i++) {
            double theta_low  = binEdge[i-1] * M_PI / 180.;
            double theta_high = binEdge[i]   * M_PI / 180.;
            double dOmega = 2. * M_PI * (cos(theta_low) - cos(theta_high));
            double N = gen_moller_raw->GetBinContent(i);
            if (dOmega > 0.) {
                gen_moller_cs->SetBinContent(i, N / (L_gen_cm2 * dOmega) * cm2_to_mb);
                gen_moller_cs->SetBinError  (i, sqrt(N) / (L_gen_cm2 * dOmega) * cm2_to_mb);
            }
            double tot = m_in[i-1] + m_out[i-1];
            double eff = (tot > 0.) ? m_in[i-1] / tot : 0.;
            double eff_err = (tot > 0.) ? sqrt(eff * (1. - eff) / tot) : 0.;
            double_arm_geo_acc->SetBinContent(i, eff * 100.);
            double_arm_geo_acc->SetBinError  (i, eff_err * 100.);
            std::cout << "Moller bin " << i
                      << ": Q2 [" << GetQ2(theta_high) << ", " << GetQ2(theta_low) << "] GeV^2,"
                      << " dσ/dΩ = " << gen_moller_cs->GetBinContent(i) << " mb/sr,"
                      << " double-arm acc = " << eff * 100. << " %" << std::endl;
        }
        delete gen_moller_raw;
    }

    // Calculate yield ratio and its error
    for (int i = 1; i <= Nbins; i++) {
        double ep_acc = geo_acceptance[i-1];
        double ee_acc = double_arm_geo_acc->GetBinContent(i) / 100.;
        double mott = mott_yield->GetBinContent(i) - mott_yield_B->GetBinContent(i); // subtract type B background
        double moller = moller_yield->GetBinContent(i) - moller_yield_B->GetBinContent(i); // subtract type B background
        mott /= ep_acc; // correct for acceptance
        moller /= ee_acc; // correct for acceptance
        double ratio = (moller > 0) ? mott / moller : 0;
        //if(yield_ratio->GetBinCenter(i) < 1.6 || yield_ratio->GetBinCenter(i) > 2.9) ratio = 0; // only keep ratio in the 1.6-2.9 degree region where Moller yield is significant
        yield_ratio->SetBinContent(i, ratio);
        double mott_err = sqrt(pow(mott_yield->GetBinError(i), 2) + pow(mott_yield_B->GetBinError(i), 2)) / ep_acc; // combine errors in quadrature and correct for acceptance
        double moller_err = sqrt(pow(moller_yield->GetBinError(i), 2) + pow(moller_yield_B->GetBinError(i), 2)) / ee_acc; // combine errors in quadrature and correct for acceptance
        double ratio_err = (moller > 0) ? ratio * sqrt(pow(mott_err/mott, 2) + pow(moller_err/moller, 2)) : 0;
        yield_ratio->SetBinError(i, ratio_err);
    }

    // Calculate cross section from ep/ee ratio and generator Moller cross section
    TH1F *cross_section_from_ratio = new TH1F("cross_section_from_ratio", "Cross Section from Yield Ratio;Q^{2} (GeV^{2});d#sigma/d#Omega (mb/sr)", Nbins, q2_edges);
    for (int i = 1; i <= Nbins; i++) {
        double ratio = yield_ratio->GetBinContent(i);
        double ratio_err = yield_ratio->GetBinError(i);
        double moller_cs = gen_moller_cs->GetBinContent(i);
        double moller_cs_err = gen_moller_cs->GetBinError(i);
        double cs = ratio * moller_cs;
        double cs_err = cs * sqrt(pow(ratio_err/ratio, 2) + pow(moller_cs_err/moller_cs, 2));
        cross_section_from_ratio->SetBinContent(i, cs);
        cross_section_from_ratio->SetBinError(i, cs_err);
    }

    TCanvas *moller_cv = new TCanvas("c_moller_cs", "Moller Cross Section (Generator)", 900, 700);
    gen_moller_cs->SetStats(0);
    gen_moller_cs->SetMarkerStyle(20);
    gen_moller_cs->SetMarkerSize(1.0);
    gen_moller_cs->SetTitle("Moller Scattering Cross Section from Generator");
    gen_moller_cs->Draw("E1P");

    TCanvas *acc_cv = new TCanvas("c_double_arm_acc", "Moller Double-Arm Geometric Acceptance", 900, 700);
    double_arm_geo_acc->SetStats(0);
    double_arm_geo_acc->SetMarkerStyle(20);
    double_arm_geo_acc->SetMarkerSize(1.0);
    double_arm_geo_acc->SetTitle("Moller Double-Arm Geometric Acceptance");
    double_arm_geo_acc->GetYaxis()->SetRangeUser(0., 105);
    double_arm_geo_acc->Draw("E1P");

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
    for(int i = 1; i <= cross_section->GetNbinsX(); i++) {
        double acc = geo_acceptance[i-1];
        if (acc > 0.) {
            cross_section->SetBinContent(i, cross_section->GetBinContent(i) / acc);
            cross_section->SetBinError(i, cross_section->GetBinError(i) / acc);
        }
    }
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

    cross_section_from_ratio->SetStats(0);
    cross_section_from_ratio->SetLineColor(kGreen + 1);
    cross_section_from_ratio->SetMarkerStyle(25);
    cross_section_from_ratio->SetMarkerColor(kGreen + 1);
    cross_section_from_ratio->Draw("E1P SAME");

    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(cross_section, "Data (measured lumi)", "lp");
    leg->AddEntry(gen_cross,     "Generator (L_{gen})",  "lp");
    leg->AddEntry(cross_section_from_ratio, "From Yield Ratio", "lp");
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
    TH1F *GE_from_ratio = new TH1F("GE_from_ratio", "Electric Form Factor from Yield Ratio;Q^{2} (GeV^{2});G_{E}", Nbins, q2_edges);
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
    for (int i = 1; i <= cross_section_from_ratio->GetNbinsX(); i++) {
        double theta = (binEdge[i-1] + binEdge[i]) / 2. * M_PI / 180.; // bin center in radians
        double CS_mb = cross_section_from_ratio->GetBinContent(i); // in mb/sr
        double GM  = GetKellyGM(GetQ2(theta));
        double eps = GetEpsilon(theta);
        double tau = GetTau(theta);
        double mott = GetMottCS(theta);
        double arg = (1. + tau) * CS_mb / mott - tau / eps * GM * GM;
        if (CS_mb <= 0. || arg <= 0.) { GE_from_ratio->SetBinContent(i, 0.); GE_from_ratio->SetBinError(i, 0.); continue; }
        GE_from_ratio->SetBinContent(i, sqrt(arg));
        double dCS = cross_section_from_ratio->GetBinError(i);
        double dGE_dCS = 0.5 * (1. + tau) / sqrt(arg) / mott;
        GE_from_ratio->SetBinError(i, dCS * fabs(dGE_dCS));
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
    GE_from_ratio->SetStats(0);
    GE_from_ratio->SetLineColor(kBlue + 1);
    GE_from_ratio->SetMarkerStyle(25);
    GE_from_ratio->SetMarkerColor(kBlue + 1);
    //GE_from_ratio->Draw("P SAME");
    TLegend *leg_GE = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg_GE->SetBorderSize(0);
    leg_GE->AddEntry(GE, "Data (measured lumi)", "lp");
    leg_GE->AddEntry(GE_gen, "Kelly Parametrization", "lp");
    leg_GE->AddEntry(GE_theory, "Kelly Parametrization (theory)", "l");
    //leg_GE->AddEntry(GE_from_ratio, "From Yield Ratio", "lp");
    leg_GE->Draw();
    cv6->SaveAs("GE.png");


}


