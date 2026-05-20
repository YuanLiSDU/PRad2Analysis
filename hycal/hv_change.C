
//"description": "Gain = exp(k * ln(Vset) + b), r2: fitting quality(1 is perfect)",
// "lms": {
//        "k": 5.65102,
//        "b": -39.9941,
//        "r2": 0.999919
//      }
void hv_change(){

        double Vset[300];
        double gain[300][2];
        double ratio[300];

        for(int i = 0; i < 300; i++){
            Vset[i] = 1000. + i;
            gain[i][0] = exp(5.85102 * log(Vset[i]+0.5) - 39.9941);
            gain[i][1] = exp(5.85102 * log(Vset[i]-0.5) - 39.9941);
            ratio[i] = (1 - gain[i][1] / gain[i][0]) * 100.;
        }

        TCanvas *c_gain = new TCanvas("c_gain", "Gain vs HV", 800, 600);
        TGraph *gr_gain = new TGraph(300);
        gr_gain->SetTitle("Gain change vs HV; Vset (V); Gain change(%)");
        for (int i = 0; i < 300; i++) {
            gr_gain->SetPoint(i, Vset[i], ratio[i]);
        }
        gr_gain->Draw("AP");
}