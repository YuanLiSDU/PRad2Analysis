#include "../EventData.h"
#include "../PhysicsTools.h"

// Input settings: each file uses the beam energy at the same array index.
const int nFiles = 3;
const char *inputFiles[nFiles] = {
    "../data/recon/0.7GeV/prad_024655_recon.root",
    "../data/recon/2.2GeV/prad_025202_filtered.root",
    "../data/recon/3.5GeV/prad_024917_recon.root"
};
float beamEnergies[nFiles] = {728.9f, 2239.5f, 3484.f}; // MeV

float resolution = 0.035f;
float mollerEnergyFraction = 0.75f;
float mollerEnergyNSigma = 3.f;

void drawCombinedSigmaEnergyFit(TGraph *mottGraph, TGraph *mollerGraph)
{
    TGraph *combinedGraph = new TGraph();
    combinedGraph->SetName("combinedSigmaVsEnergy");

    for (int pointIndex = 0; pointIndex < mottGraph->GetN(); ++pointIndex) {
        combinedGraph->SetPoint(combinedGraph->GetN(),
                                mottGraph->GetX()[pointIndex],
                                mottGraph->GetY()[pointIndex]);
    }
    for (int pointIndex = 0; pointIndex < mollerGraph->GetN(); ++pointIndex) {
        combinedGraph->SetPoint(combinedGraph->GetN(),
                                mollerGraph->GetX()[pointIndex],
                                mollerGraph->GetY()[pointIndex]);
    }

    if (combinedGraph->GetN() < 2) {
        std::cerr << "Not enough valid points for the sigma-vs-energy fit." << std::endl;
        return;
    }

    combinedGraph->Sort();
    TCanvas *canvas = new TCanvas("c_sigma_vs_energy", "Position sigma vs energy", 900, 650);
    combinedGraph->SetTitle("Position resolution vs electron energy;Electron energy E (GeV);Average #sigma (mm)");
    combinedGraph->SetMarkerStyle(20);
    combinedGraph->SetMarkerSize(1.2);
    combinedGraph->SetMarkerColor(kBlue + 1);
    combinedGraph->Draw("AP");

    TF1 *energyFit = new TF1("combinedEnergyFit",
                              "[0]/sqrt(x) + [1]/x + [2]", 0.1, 10.);
    energyFit->SetParNames("k1", "k2", "k3");
    energyFit->SetParameters(combinedGraph->GetY()[0]
                             * std::sqrt(combinedGraph->GetX()[0]), 0., 0.);
    energyFit->SetLineColor(kRed + 1);
    energyFit->SetLineWidth(2);
    combinedGraph->Fit(energyFit, "Q");

    TLatex *fitLabel = new TLatex();
    fitLabel->SetNDC();
    fitLabel->SetTextSize(0.034);
    fitLabel->DrawLatex(0.17, 0.85,
                        "#sigma = k_{1}/#sqrt{E [GeV]} + k_{2}/E [GeV] + k_{3}");
    fitLabel->DrawLatex(0.17, 0.79,
                        Form("k_{1} = %.4f #pm %.4f mm #sqrt{GeV}",
                             energyFit->GetParameter(0), energyFit->GetParError(0)));
    fitLabel->DrawLatex(0.17, 0.73,
                        Form("k_{2} = %.4f #pm %.4f mm GeV",
                             energyFit->GetParameter(1), energyFit->GetParError(1)));
    fitLabel->DrawLatex(0.17, 0.67,
                        Form("k_{3} = %.4f #pm %.4f mm",
                             energyFit->GetParameter(2), energyFit->GetParError(2)));
    canvas->Update();
}

void position_sigma(){

    gStyle->SetOptStat(1111);

    // Filled once per input file and plotted after the file loop.
    TGraph *sigmaVsEnergy = new TGraph();
    sigmaVsEnergy->SetName("sigmaVsEnergy");
    TGraph *mollerSigmaVsEnergy = new TGraph();
    mollerSigmaVsEnergy->SetName("mollerSigmaVsEnergy");

    for (int fileIndex = 0; fileIndex < nFiles; ++fileIndex) {
        TFile *file = TFile::Open(inputFiles[fileIndex]);
        if (!file || file->IsZombie()) {
            std::cerr << "Cannot open input file: " << inputFiles[fileIndex] << std::endl;
            if (file) file->Close();
            continue;
        }

        TTree *tree = (TTree*)file->Get("recon");
        if (!tree) {
            std::cerr << "Cannot find tree 'recon' in: " << inputFiles[fileIndex] << std::endl;
            file->Close();
            continue;
        }

        ReconEventData ev;
        setupReconBranches(tree, ev);

        TH1F *deltaX_gem1 = new TH1F(Form("deltaX_gem1_%d", fileIndex), "Difference in X between two GEM hits;#DeltaX (mm);Counts", 400, -50, 50);
        TH1F *deltaY_gem1 = new TH1F(Form("deltaY_gem1_%d", fileIndex), "Difference in Y between GEM1 hit and cluster;#DeltaY (mm);Counts", 400, -50, 50);

        TH1F *deltaX_gem2 = new TH1F(Form("deltaX_gem2_%d", fileIndex), "Difference in X between GEM2 hit and cluster;#DeltaX (mm);Counts", 400, -50, 50);
        TH1F *deltaY_gem2 = new TH1F(Form("deltaY_gem2_%d", fileIndex), "Difference in Y between GEM2 hit and cluster;#DeltaY (mm);Counts", 400, -50, 50);

        TH1F *mollerDeltaX_gem1 = new TH1F(Form("mollerDeltaX_gem1_%d", fileIndex), "Moller: difference in X between GEM1 hit and cluster;#DeltaX (mm);Counts", 400, -50, 50);
        TH1F *mollerDeltaY_gem1 = new TH1F(Form("mollerDeltaY_gem1_%d", fileIndex), "Moller: difference in Y between GEM1 hit and cluster;#DeltaY (mm);Counts", 400, -50, 50);
        TH1F *mollerDeltaX_gem2 = new TH1F(Form("mollerDeltaX_gem2_%d", fileIndex), "Moller: difference in X between GEM2 hit and cluster;#DeltaX (mm);Counts", 400, -50, 50);
        TH1F *mollerDeltaY_gem2 = new TH1F(Form("mollerDeltaY_gem2_%d", fileIndex), "Moller: difference in Y between GEM2 hit and cluster;#DeltaY (mm);Counts", 400, -50, 50);

        // Keep the histograms alive after their input file is closed.
        deltaX_gem1->SetDirectory(nullptr);
        deltaY_gem1->SetDirectory(nullptr);
        deltaX_gem2->SetDirectory(nullptr);
        deltaY_gem2->SetDirectory(nullptr);
        mollerDeltaX_gem1->SetDirectory(nullptr);
        mollerDeltaY_gem1->SetDirectory(nullptr);
        mollerDeltaX_gem2->SetDirectory(nullptr);
        mollerDeltaY_gem2->SetDirectory(nullptr);

        const Long64_t entriesToProcess = tree->GetEntries() / 5;
        std::cout << "Processing " << inputFiles[fileIndex]
                  << " (Ebeam = " << beamEnergies[fileIndex] << " MeV)" << std::endl;

        for (Long64_t i = 0; i < entriesToProcess; i++) {
            tree->GetEntry(i);

            if(i%10000 == 0)
                std::cout << "  " << i << " / " << entriesToProcess << "\r" << std::flush;

            if(ev.n_clusters == 1 && ev.matchNum == 1){
                if(!isMott(ev.cl_energy[0], beamEnergies[fileIndex], resolution)) continue;
                if( fabs(ev.cl_x[0]) < 20.25 * 2.5 && fabs(ev.cl_y[0]) < 20.25 * 2.5 ) continue;
                if( fabs(ev.cl_x[0]) > 20.25 * 16. || fabs(ev.cl_y[0]) > 20.25 * 16. ) continue;

                float x[4], y[4], z[4], E;

                x[0] = 0.f; y[0] = 0.f; z[0] = 0.f; // target center
                x[1] = ev.mHit_gx[0][1]; y[1] = ev.mHit_gy[0][1]; z[1] = ev.mHit_gz[0][1];
                x[2] = ev.mHit_gx[0][0]; y[2] = ev.mHit_gy[0][0]; z[2] = ev.mHit_gz[0][0];
                x[3] = ev.cl_x[0]; y[3] = ev.cl_y[0]; z[3] = ev.cl_z[0];
                float scale1 = z[3] / z[1], scale2 = z[3] / z[2];
                x[1] *= scale1; y[1] *= scale1; z[1] *= scale1;
                x[2] *= scale2; y[2] *= scale2; z[2] *= scale2;
                E = ev.cl_energy[0];
                float theta = atan2(std::sqrt(x[3]*x[3] + y[3]*y[3]), z[3]) * 180.f / M_PI;

                deltaX_gem1->Fill(x[1] - x[3]);
                deltaY_gem1->Fill(y[1] - y[3]);
                deltaX_gem2->Fill(x[2] - x[3]);
                deltaY_gem2->Fill(y[2] - y[3]);
            }
            if(ev.n_clusters == 2 && ev.matchNum == 2) {
                float theta[2];
                for (int hitIndex = 0; hitIndex < 2; ++hitIndex) {
                    theta[hitIndex] = atan2(std::sqrt(ev.cl_x[hitIndex]*ev.cl_x[hitIndex]
                                                     + ev.cl_y[hitIndex]*ev.cl_y[hitIndex]),
                                            ev.cl_z[hitIndex]) * 180.f / M_PI;
                }

                if (!isMoller_kinematic(theta[0], ev.cl_energy[0],
                                        theta[1], ev.cl_energy[1],
                                        beamEnergies[fileIndex], resolution)) continue;

                MollerEvent moller = {
                    DataPoint(ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]),
                    DataPoint(ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1])
                };
                if (std::fabs(GetMollerPhiDiff(moller)) > 10.f) continue;

                const float selectedEnergy = mollerEnergyFraction * beamEnergies[fileIndex];
                const float selectedEnergySigma = resolution * selectedEnergy
                                                / std::sqrt(selectedEnergy / 1000.f);

                for (int hitIndex = 0; hitIndex < 2; ++hitIndex) {
                    if (std::fabs(ev.cl_energy[hitIndex] - selectedEnergy)
                        > mollerEnergyNSigma * selectedEnergySigma) continue;
                    if (std::fabs(ev.cl_x[hitIndex]) < 20.25 * 2.5
                        && std::fabs(ev.cl_y[hitIndex]) < 20.25 * 2.5) continue;
                    if (std::fabs(ev.cl_x[hitIndex]) > 20.25 * 16.
                        || std::fabs(ev.cl_y[hitIndex]) > 20.25 * 16.) continue;

                    const float gem1Scale = ev.cl_z[hitIndex] / ev.mHit_gz[hitIndex][1];
                    const float gem2Scale = ev.cl_z[hitIndex] / ev.mHit_gz[hitIndex][0];
                    const float gem1X = ev.mHit_gx[hitIndex][1] * gem1Scale;
                    const float gem1Y = ev.mHit_gy[hitIndex][1] * gem1Scale;
                    const float gem2X = ev.mHit_gx[hitIndex][0] * gem2Scale;
                    const float gem2Y = ev.mHit_gy[hitIndex][0] * gem2Scale;

                    mollerDeltaX_gem1->Fill(gem1X - ev.cl_x[hitIndex]);
                    mollerDeltaY_gem1->Fill(gem1Y - ev.cl_y[hitIndex]);
                    mollerDeltaX_gem2->Fill(gem2X - ev.cl_x[hitIndex]);
                    mollerDeltaY_gem2->Fill(gem2Y - ev.cl_y[hitIndex]);
                }
            }
        }
        std::cout << std::endl;

        TCanvas *c = new TCanvas(Form("c_%d", fileIndex),
                                 Form("Ebeam = %.1f MeV", beamEnergies[fileIndex]),
                                 800, 600);
        c->Divide(2,2);
        c->cd(1); deltaX_gem1->Draw(); deltaX_gem1->Fit("gaus","Q","", -4, 4);
        c->cd(2); deltaY_gem1->Draw(); deltaY_gem1->Fit("gaus","Q","", -4, 4);
        c->cd(3); deltaX_gem2->Draw(); deltaX_gem2->Fit("gaus","Q","", -4, 4);
        c->cd(4); deltaY_gem2->Draw(); deltaY_gem2->Fit("gaus","Q","", -4, 4);
        c->Update();

        TH1F *histograms[4] = {
            deltaX_gem1, deltaY_gem1, deltaX_gem2, deltaY_gem2
        };
        double sigmaSum = 0.;
        int validSigmaCount = 0;
        for (int histogramIndex = 0; histogramIndex < 4; ++histogramIndex) {
            TF1 *gaussianFit = histograms[histogramIndex]->GetFunction("gaus");
            if (gaussianFit) {
                sigmaSum += std::fabs(gaussianFit->GetParameter(2));
                ++validSigmaCount;
            }
        }

        if (validSigmaCount == 4) {
            const double averageSigma = sigmaSum / 4.;
            const double energyGeV = beamEnergies[fileIndex] / 1000.;
            sigmaVsEnergy->SetPoint(sigmaVsEnergy->GetN(), energyGeV, averageSigma);
            std::cout << "Average sigma at " << energyGeV << " GeV: "
                      << averageSigma << " mm" << std::endl;
        } else {
            std::cerr << "Only " << validSigmaCount
                      << " of 4 Gaussian fits succeeded for "
                      << inputFiles[fileIndex] << std::endl;
        }

        TCanvas *cMoller = new TCanvas(Form("c_moller_%d", fileIndex),
                                       Form("Moller, Ebeam = %.1f MeV", beamEnergies[fileIndex]),
                                       800, 600);
        cMoller->Divide(2,2);
        cMoller->cd(1); mollerDeltaX_gem1->Draw(); mollerDeltaX_gem1->Fit("gaus","Q","", -4, 4);
        cMoller->cd(2); mollerDeltaY_gem1->Draw(); mollerDeltaY_gem1->Fit("gaus","Q","", -4, 4);
        cMoller->cd(3); mollerDeltaX_gem2->Draw(); mollerDeltaX_gem2->Fit("gaus","Q","", -4, 4);
        cMoller->cd(4); mollerDeltaY_gem2->Draw(); mollerDeltaY_gem2->Fit("gaus","Q","", -4, 4);
        cMoller->Update();

        TH1F *mollerHistograms[4] = {
            mollerDeltaX_gem1, mollerDeltaY_gem1,
            mollerDeltaX_gem2, mollerDeltaY_gem2
        };
        double mollerSigmaSum = 0.;
        int validMollerSigmaCount = 0;
        for (int histogramIndex = 0; histogramIndex < 4; ++histogramIndex) {
            TF1 *gaussianFit = mollerHistograms[histogramIndex]->GetFunction("gaus");
            if (gaussianFit) {
                mollerSigmaSum += std::fabs(gaussianFit->GetParameter(2));
                ++validMollerSigmaCount;
            }
        }

        if (validMollerSigmaCount == 4) {
            const double averageMollerSigma = mollerSigmaSum / 4.;
            const double energyGeV = mollerEnergyFraction * beamEnergies[fileIndex] / 1000.;
            mollerSigmaVsEnergy->SetPoint(mollerSigmaVsEnergy->GetN(),
                                          energyGeV, averageMollerSigma);
            std::cout << "Average Moller sigma at " << energyGeV << " GeV: "
                      << averageMollerSigma << " mm" << std::endl;
        } else {
            std::cerr << "Only " << validMollerSigmaCount
                      << " of 4 Moller Gaussian fits succeeded for "
                      << inputFiles[fileIndex] << std::endl;
        }

        file->Close();
    }

    drawCombinedSigmaEnergyFit(sigmaVsEnergy, mollerSigmaVsEnergy);
}
