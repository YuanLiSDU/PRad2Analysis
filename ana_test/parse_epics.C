// parse_epics.C
// Read the "epics" TTree from a PRad2 ROOT file and write all EPICS values
// to a plain text file, one entry per block.
// Usage: root -l -q 'parse_epics.C("input.root","output.txt")'

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>

void parse_epics(const char *input  = "../data/0.7GeV/prad_024639_filter.root",
                 const char *output = "epics_values.txt")
{
    TFile *fin = TFile::Open(input);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input file: " << input << std::endl;
        return;
    }

    TTree *tree = (TTree*)fin->Get("epics");
    if (!tree) {
        std::cerr << "No 'epics' tree found in " << input << std::endl;
        fin->Close();
        return;
    }

    // branch variables
    Int_t    event_num   = 0;
    UInt_t   unix_time   = 0;
    UInt_t   sync_cnt    = 0;
    UInt_t   run_num     = 0;
    Bool_t   good        = false;
    std::vector<std::string> *channel = nullptr;
    std::vector<double>      *value   = nullptr;

    tree->SetBranchAddress("event_number_at_arrival", &event_num);
    tree->SetBranchAddress("unix_time",               &unix_time);
    tree->SetBranchAddress("sync_counter",            &sync_cnt);
    tree->SetBranchAddress("run_number",              &run_num);
    tree->SetBranchAddress("good",                    &good);
    tree->SetBranchAddress("channel",                 &channel);
    tree->SetBranchAddress("value",                   &value);

    std::ofstream fout(output);
    if (!fout.is_open()) {
        std::cerr << "Cannot open output file: " << output << std::endl;
        fin->Close();
        return;
    }

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // format unix_time as human-readable string
        time_t t = (time_t)unix_time;
        char tbuf[32];
        strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", localtime(&t));

        fout << "--- Entry " << i << " ---\n";
        fout << "  " << std::left << std::setw(45) << "run_number"
             << " = " << run_num  << "\n";
        fout << "  " << std::left << std::setw(45) << "event_number_at_arrival"
             << " = " << event_num << "\n";
        fout << "  " << std::left << std::setw(45) << "sync_counter"
             << " = " << sync_cnt  << "\n";
        fout << "  " << std::left << std::setw(45) << "unix_time"
             << " = " << unix_time << "  (" << tbuf << ")\n";
        fout << "  " << std::left << std::setw(45) << "good"
             << " = " << (good ? "true" : "false") << "\n";

        if (channel && value) {
            size_t n = std::min(channel->size(), value->size());
            for (size_t j = 0; j < n; ++j)
                fout << "  " << std::left << std::setw(45) << (*channel)[j]
                     << " = " << (*value)[j] << "\n";
        }
        fout << "\n";
    }
    fout.close();
    fin->Close();

    std::cout << "Wrote " << nentries << " EPICS entries to: " << output << "\n";
}
