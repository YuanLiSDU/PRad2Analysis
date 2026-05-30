

void rand_trigger(){

    TFile *f = new TFile("../data/0.7GeV/prad_024640_filter.root");
    TTree *t = (TTree*)f->Get("events");

    RawEventData ev;
    setupRawBranches(t, ev);

    for (Long64_t i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        if(i % 10000 == 0) std::cout << "Entry " << i << " / " << t->GetEntries() << "\r" << std::flush;

        bool isPulser = (ev.trigger_bits & 1u << 15) != 0;
        bool isRawsum = (ev.trigger_bits & 1u << 8) != 0;
        
        if(!isPulser && !isRawsum) continue; // only look at pulser and rawsum triggers

    }

}