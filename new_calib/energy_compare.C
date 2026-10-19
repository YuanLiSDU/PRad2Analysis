
#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>

#include <memory>
#include <stdexcept>

namespace EnergyCompare {

struct AllGraphs {
    std::unique_ptr<TGraphErrors> energy;
    std::unique_ptr<TGraphErrors> resolution;
};

std::unique_ptr<TGraphErrors> cloneAll(TCanvas &canvas, const char *padName,
                                      const char *newName, const char *path)
{
    auto *pad = dynamic_cast<TPad *>(canvas.GetPrimitive(padName));
    if (!pad)
        throw std::runtime_error(TString::Format("Missing %s in %s", padName, path).Data());

    TGraphErrors *all = nullptr;
    for (auto *object : *pad->GetListOfPrimitives()) {
        auto *graph = dynamic_cast<TGraphErrors *>(object);
        // Both source canvases use blue filled circles for the "All" sample.
        if (!graph || graph->GetMarkerColor() != kBlue + 1 ||
            graph->GetMarkerStyle() != 20)
            continue;
        if (all)
            throw std::runtime_error(TString::Format("Ambiguous All graph in %s/%s", path, padName).Data());
        all = graph;
    }
    if (!all || all->GetN() == 0)
        throw std::runtime_error(TString::Format("Missing or empty All graph in %s/%s", path, padName).Data());

    auto clone = std::unique_ptr<TGraphErrors>(
        static_cast<TGraphErrors *>(all->Clone(newName)));
    clone->GetListOfFunctions()->Clear();
    return clone;
}

AllGraphs loadAll(const char *path, const char *tag)
{
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    if (!file || file->IsZombie())
        throw std::runtime_error(TString::Format("Cannot open %s", path).Data());

    auto *canvas = file->Get<TCanvas>("recon_energy_resolution_vs_angle");
    if (!canvas)
        throw std::runtime_error(TString::Format("Missing recon_energy_resolution_vs_angle in %s", path).Data());

    return {
        cloneAll(*canvas, "upper_pad", TString::Format("energy_all_%s", tag), path),
        cloneAll(*canvas, "lower_pad", TString::Format("resolution_all_%s", tag), path)
    };
}

void drawPair(TPad &pad, std::unique_ptr<TGraphErrors> &first,
              std::unique_ptr<TGraphErrors> &second, const char *title,
              const char *label1, const char *label2, const char *name)
{
    pad.cd();
    pad.SetLeftMargin(0.15);
    pad.SetRightMargin(0.04);
    pad.SetBottomMargin(0.14);
    pad.SetTopMargin(0.10);
    pad.SetTicks(1, 1);

    first->SetMarkerColor(kBlue + 1);
    first->SetLineColor(kBlue + 1);
    first->SetMarkerStyle(20);
    second->SetMarkerColor(kRed + 1);
    second->SetLineColor(kRed + 1);
    second->SetMarkerStyle(21);
    first->SetMarkerSize(0.85);
    second->SetMarkerSize(0.85);

    auto *legend = new TLegend(0.57, 0.74, 0.94, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(first.get(), label1, "pe");
    legend->AddEntry(second.get(), label2, "pe");

    auto *graphs = new TMultiGraph(name, title);
    graphs->Add(first.release(), "P");
    graphs->Add(second.release(), "P");
    graphs->Draw("A");
    legend->Draw();
}

} // namespace EnergyCompare

// Defaults are relative to this macro, not the ROOT working directory.
// Run: root -l -b -q 'new_calib/energy_compare.C()'
void energy_compare(const char *file1 = nullptr, const char *file2 = nullptr,
                    const char *outputPrefix = nullptr)
{
    const TString directory = gSystem->DirName(__FILE__);
    const TString path1 = file1 ? TString(file1) : directory + "/round1_recon_25211.root";
    const TString path2 = file2 ? TString(file2) : directory + "/round1_recon_test.root";
    const TString output = outputPrefix ? TString(outputPrefix) : directory + "/energy_compare";
    if (output + ".root" == path1 || output + ".root" == path2)
        throw std::runtime_error("Output ROOT file must not overwrite either input");

    auto first = EnergyCompare::loadAll(path1, "25211");
    auto second = EnergyCompare::loadAll(path2, "test");
    const TString label1 = "Test run(not calibrated)";
    const TString label2 = "Calibrated Runs";

    auto *canvas = new TCanvas("c_energy_compare", "All events: energy comparison", 1000, 900);
    canvas->Divide(1, 2, 0.0, 0.01);
    EnergyCompare::drawPair(*static_cast<TPad *>(canvas->cd(1)),
        first.energy, second.energy,
        "Reconstructed energy vs angle (All);Scattering angle [deg];E_{recon} [MeV]",
        label1, label2, "energy_all_compare");
    EnergyCompare::drawPair(*static_cast<TPad *>(canvas->cd(2)),
        first.resolution, second.resolution,
        "Energy resolution vs angle (All);Scattering angle [deg];#sigma / E #times #sqrt{E [GeV]} [%]",
        label1, label2, "resolution_all_compare");
    canvas->Update();

    TFile outputFile(output + ".root", "RECREATE");
    if (outputFile.IsZombie())
        throw std::runtime_error(TString::Format("Cannot create %s.root", output.Data()).Data());
    if (canvas->Write("energy_compare") <= 0)
        throw std::runtime_error("Failed to write comparison canvas");
    outputFile.Close();
    canvas->SaveAs(output + ".png");
    canvas->SaveAs(output + ".pdf");
    Info("energy_compare", "Saved %s.png, %s.pdf and %s.root",
         output.Data(), output.Data(), output.Data());
}