#include "../EventData.h"
#include "../PhysicsTools.h"

std::string file_path = "../data/time_study/";

void moller_time_diff(){

    TChain chain("recon");
    chain.Add(file_path.c_str());

    ReconEventData ev;
    setupReconBranches(&chain, ev);

    for (Long64_t i = 0; i < chain.GetEntries(); i++) {
        if (i % 10000 == 0)
            std::cout << "  " << i << " / " << chain.GetEntries() << "\r" << std::flush;
        chain.GetEntry(i);

        bool is_sum = ((ev.trigger_bits & 1u << 8) != 0 && ev.nch < 100);
        
        if(ev.matchNum == 2) continue;

        for(int j = 0; j < ev.matchNum; j++){
            float t1 = ev.mHit_gz[j][1];
            float t2 = ev.mHit_gz[j][0];
            float dt = t1 - t2;
            if(fabs(dt) > 10.) continue;
            std::cout << "Event: " << ev.event_num << ", GEM hit time difference: " << dt << std::endl;
        }
    }
}

