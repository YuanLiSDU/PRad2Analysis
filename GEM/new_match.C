#include "../EventData.h"
#include "../PhysicsTools.h"

float Ebeam = 2239.5f; // MeV

string inputDir_old = "";
string inputDir_new = "";

void new_match() {
    TChain chain("recon");
    chain.Add(inputDir_old.c_str());

    ReconEventData data;
    setupReconBranches(&chain, data);
}