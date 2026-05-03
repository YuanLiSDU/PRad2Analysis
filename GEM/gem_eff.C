
#include "../EventData.h"
#include "../PhysicsTools.h"
#include <algorithm>
#include <numeric>

float gx[4] = {0, 0, 0, 0};
float gy[4] = {0, 0, 0, 0};
float gz[4] = {5805.0,  5850.0, 5413.0, 5458.9};

float Ebeam = 3900.f; // MeV, can adjust as needed for different beam energies

void setupReconBranches(TTree *tree, ReconEventData &ev);
void TransformToGEMFrame(float &x, float &y, float &z, int gem_id);

// 从文件名中提取第一段连续数字作为 run number，例如 24190_all.root -> 24190
int extractRunNumber(const TString &fname) {
    TString base = gSystem->BaseName(fname.Data());
    TString num_str = "";
    for (int i = 0; i < base.Length(); i++) {
        char c = base[i];
        if (isdigit(c)) {
            num_str += c;
        } else if (num_str.Length() > 0) {
            break; // 遇到第一段数字结束就停止
        }
    }
    return (num_str.Length() > 0) ? num_str.Atoi() : -1;
}

// 处理单个文件，计算各 GEM 的效率和误差（Binomial 误差）
// 若 out_gem_eff / out_gem_2match_eff 不为 nullptr，则填充 2D 效率图
void processOneFile(const TString &fname,
                    double eff[4],       double eff_err[4],
                    double eff2[4],      double eff2_err[4],
                    TH2F *out_gem_eff[4],
                    TH2F *out_gem_2match_eff[4])
{
    TChain chain("recon");
    chain.Add(fname);

    ReconEventData data;
    setupReconBranches(&chain, data);

    float gem_x_range_lo[4] = {-250., 0., -250., 0.};
    float gem_x_range_hi[4] = {0., 250., 0., 250.};
    const float gem_y_lo = -250., gem_y_hi = 250.;

    // 用文件名做直方图名称前缀，避免多次调用时名称冲突
    TString tag = gSystem->BaseName(fname.Data());
    tag.ReplaceAll(".root", "");
    tag.ReplaceAll(".", "_");

    TH2F *h_should[4], *h_match[4], *h_2sh[4], *h_2mh[4];
    for (int i = 0; i < 4; i++) {
        h_should[i] = new TH2F(Form("h_sh_%s_%d",  tag.Data(), i), "", 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_lo, gem_y_hi);
        h_match[i]  = new TH2F(Form("h_mh_%s_%d",  tag.Data(), i), "", 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_lo, gem_y_hi);
        h_2sh[i]    = new TH2F(Form("h_2sh_%s_%d", tag.Data(), i), "", 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_lo, gem_y_hi);
        h_2mh[i]    = new TH2F(Form("h_2mh_%s_%d", tag.Data(), i), "", 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, gem_y_lo, gem_y_hi);
    }

    int nEntries = chain.GetEntries();
    std::cout << "Processing " << fname << " (" << nEntries << " entries) ..." << std::endl;

    for (int i = 0; i < nEntries; i++) {
        chain.GetEntry(i);
        if (i % 10000 == 0)
            std::cout << "  " << i << "/" << nEntries << "\r" << std::flush;

        for (int j = 0; j < data.n_clusters; j++) {
            if (fabs(data.cl_energy[j] - Ebeam) > 400.) continue;
            if (fabs(data.cl_x[j]) < 20.75*2.5 && fabs(data.cl_y[j]) < 20.75*2.5) continue;

            // 每次使用局部副本，避免坐标被累积修改
            float cx = data.cl_x[j], cy = data.cl_y[j], cz = data.cl_z[j];

            for (int k = 0; k < 4; k++) {
                float tx = cx, ty = cy, tz = cz;
                TransformToGEMFrame(tx, ty, tz, k);
                h_should[k]->Fill(tx, ty);
                if ((data.matchFlag[j] & (1 << k)) != 0)
                    h_match[k]->Fill(data.matchGEMx[j][k], data.matchGEMy[j][k]);
            }

            // gem 0: 要求 gem2 或 gem3 也匹配
            if (data.matchFlag[j] & (1<<2) || data.matchFlag[j] & (1<<3)) {
                float tx=cx, ty=cy, tz=cz; TransformToGEMFrame(tx, ty, tz, 0);
                h_2sh[0]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<0)) h_2mh[0]->Fill(data.matchGEMx[j][0], data.matchGEMy[j][0]);
            }
            // gem 2: 要求 gem0 或 gem1 也匹配
            if (data.matchFlag[j] & (1<<0) || data.matchFlag[j] & (1<<1)) {
                float tx=cx, ty=cy, tz=cz; TransformToGEMFrame(tx, ty, tz, 2);
                h_2sh[2]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<2)) h_2mh[2]->Fill(data.matchGEMx[j][2], data.matchGEMy[j][2]);
            }
            // gem 1: 要求 gem3 或 gem2 也匹配
            if (data.matchFlag[j] & (1<<3) || data.matchFlag[j] & (1<<2)) {
                float tx=cx, ty=cy, tz=cz; TransformToGEMFrame(tx, ty, tz, 1);
                h_2sh[1]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<1)) h_2mh[1]->Fill(data.matchGEMx[j][1], data.matchGEMy[j][1]);
            }
            // gem 3: 要求 gem1 或 gem0 也匹配
            if (data.matchFlag[j] & (1<<1) || data.matchFlag[j] & (1<<0)) {
                float tx=cx, ty=cy, tz=cz; TransformToGEMFrame(tx, ty, tz, 3);
                h_2sh[3]->Fill(tx, ty);
                if (data.matchFlag[j] & (1<<3)) h_2mh[3]->Fill(data.matchGEMx[j][3], data.matchGEMy[j][3]);
            }
        }
    }
    std::cout << std::endl;

    for (int i = 0; i < 4; i++) {
        double ns  = h_should[i]->Integral(), nm  = h_match[i]->Integral();
        double n2s = h_2sh[i]->Integral(),    n2m = h_2mh[i]->Integral();
        eff[i]      = (ns  > 0) ? nm  / ns  : 0.;
        eff_err[i]  = (ns  > 0) ? sqrt(eff[i]  * (1 - eff[i])  / ns)  : 0.;
        eff2[i]     = (n2s > 0) ? n2m / n2s : 0.;
        eff2_err[i] = (n2s > 0) ? sqrt(eff2[i] * (1 - eff2[i]) / n2s) : 0.;
        std::cout << Form("  GEM%d: eff=%.2f%% (%.0f/%.0f),  2match_eff=%.2f%% (%.0f/%.0f)",
                          i, eff[i]*100, nm, ns, eff2[i]*100, n2m, n2s) << std::endl;

        if (out_gem_eff && out_gem_eff[i])
            out_gem_eff[i]->Divide(h_match[i], h_should[i], 1, 1, "B");
        if (out_gem_2match_eff && out_gem_2match_eff[i])
            out_gem_2match_eff[i]->Divide(h_2mh[i], h_2sh[i], 1, 1, "B");
    }

    for (int i = 0; i < 4; i++) {
        delete h_should[i]; delete h_match[i];
        delete h_2sh[i];    delete h_2mh[i];
    }
}

void gem_eff(){
    // 收集命令行输入的 root 文件
    std::vector<TString> input_files;
    int argc = gApplication->Argc();
    char **argv = gApplication->Argv();
    for (int i = 1; i < argc; i++) {
        TString arg(argv[i]);
        if (arg.EndsWith(".root") && !arg.BeginsWith("-")) {
            input_files.push_back(arg);
            std::cout << "Adding file: " << arg << std::endl;
        }
    }
    if (input_files.empty()) {
        input_files.push_back("../data/24017.root");
        std::cout << "No input files specified, using default: ../data/24017.root" << std::endl;
    }

    int nFiles = input_files.size();
    bool single_file = (nFiles == 1);

    float gem_x_range_lo[4] = {-250., 0., -250., 0.};
    float gem_x_range_hi[4] = {0., 250., 0., 250.};

    // 单文件时创建 2D 效率图
    TH2F *gem_eff_2d[4]       = {};
    TH2F *gem_2match_eff_2d[4]= {};
    if (single_file) {
        for (int i = 0; i < 4; i++) {
            gem_eff_2d[i]        = new TH2F(Form("h2_eff_%d",       i), Form("GEM%d Efficiency; x (mm); y (mm)",        i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, -250, 250);
            gem_2match_eff_2d[i] = new TH2F(Form("h2_2match_eff_%d",i), Form("GEM%d 2-Match Efficiency; x (mm); y (mm)",i), 25, gem_x_range_lo[i], gem_x_range_hi[i], 50, -250, 250);
        }
    }

    // 收集每个文件的结果
    std::vector<double> run_nums;
    std::vector<double> v_eff[4],  v_eff_err[4];
    std::vector<double> v_eff2[4], v_eff2_err[4];
    std::vector<double> v_xerr;

    for (auto &fname : input_files) {
        double eff[4], eff_err[4], eff2[4], eff2_err[4];
        processOneFile(fname, eff, eff_err, eff2, eff2_err,
                       single_file ? gem_eff_2d        : nullptr,
                       single_file ? gem_2match_eff_2d : nullptr);
        run_nums.push_back(extractRunNumber(fname));
        v_xerr.push_back(0.);
        for (int k = 0; k < 4; k++) {
            v_eff[k].push_back(eff[k]);     v_eff_err[k].push_back(eff_err[k]);
            v_eff2[k].push_back(eff2[k]);   v_eff2_err[k].push_back(eff2_err[k]);
        }
    }

    // 按 run number 排序
    std::vector<int> order(nFiles);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b){ return run_nums[a] < run_nums[b]; });

    auto reorder = [&](const std::vector<double> &v) {
        std::vector<double> r(nFiles);
        for (int i = 0; i < nFiles; i++) r[i] = v[order[i]];
        return r;
    };
    auto sruns = reorder(run_nums);
    auto sxerr = reorder(v_xerr);

    // 颜色和点型
    int colors[4]  = {kRed, kBlue, kGreen+2, kMagenta};
    int markers[4] = {20, 21, 22, 23};

    // 多文件时画效率随 run 变化的趋势图
    if (!single_file) {
        // 简单效率趋势图
        TCanvas *c_trend = new TCanvas("c_trend", "GEM Overall Efficiency vs Run", 900, 600);
        c_trend->SetGrid();
        TMultiGraph *mg = new TMultiGraph("mg_eff", "GEM Overall Efficiency vs Run; Run Number; Efficiency");
        TLegend *leg = new TLegend(0.12, 0.15, 0.38, 0.42);
        for (int k = 0; k < 4; k++) {
            auto se = reorder(v_eff[k]), see = reorder(v_eff_err[k]);
            TGraphErrors *gr = new TGraphErrors(nFiles, sruns.data(), se.data(), sxerr.data(), see.data());
            gr->SetName(Form("gr_eff_%d", k));
            gr->SetMarkerStyle(markers[k]); gr->SetMarkerSize(1.2);
            gr->SetMarkerColor(colors[k]);  gr->SetLineColor(colors[k]);
            mg->Add(gr, "PL");
            leg->AddEntry(gr, Form("GEM%d", k), "PL");
        }
        mg->Draw("A");
        mg->GetYaxis()->SetRangeUser(0., 1.05);
        leg->Draw();
        c_trend->Update();

        // 2-match 效率趋势图
        TCanvas *c_trend2 = new TCanvas("c_trend2", "GEM 2-Match Efficiency vs Run", 900, 600);
        c_trend2->SetGrid();
        TMultiGraph *mg2 = new TMultiGraph("mg_2match_eff", "GEM 2-Match Efficiency vs Run; Run Number; Efficiency");
        TLegend *leg2 = new TLegend(0.12, 0.15, 0.38, 0.42);
        for (int k = 0; k < 4; k++) {
            auto se = reorder(v_eff2[k]), see = reorder(v_eff2_err[k]);
            TGraphErrors *gr = new TGraphErrors(nFiles, sruns.data(), se.data(), sxerr.data(), see.data());
            gr->SetName(Form("gr_2match_eff_%d", k));
            gr->SetMarkerStyle(markers[k]); gr->SetMarkerSize(1.2);
            gr->SetMarkerColor(colors[k]);  gr->SetLineColor(colors[k]);
            mg2->Add(gr, "PL");
            leg2->AddEntry(gr, Form("GEM%d", k), "PL");
        }
        mg2->Draw("A");
        mg2->GetYaxis()->SetRangeUser(0., 1.05);
        leg2->Draw();
        c_trend2->Update();
    }

    // 单文件时画 2D 效率分布图
    if (single_file) {
        TCanvas *c_eff = new TCanvas("c_eff", "GEM Efficiency", 1200, 800);
        c_eff->Divide(2, 2);
        for (int i = 0; i < 4; i++) {
            c_eff->cd(i+1);
            gem_eff_2d[i]->SetStats(0);
            gem_eff_2d[i]->Draw("COLZ");
        }

        TCanvas *c_2match_eff = new TCanvas("c_2match_eff", "GEM 2-Match Efficiency", 1200, 800);
        c_2match_eff->Divide(2, 2);
        for (int i = 0; i < 4; i++) {
            c_2match_eff->cd(i+1);
            gem_2match_eff_2d[i]->SetStats(0);
            gem_2match_eff_2d[i]->Draw("COLZ");
        }
    }
}

void setupReconBranches(TTree *tree, ReconEventData &ev){
    tree->SetBranchAddress("event_num",    &ev.event_num);
    tree->SetBranchAddress("trigger_bits", &ev.trigger_bits);
    tree->SetBranchAddress("timestamp",    &ev.timestamp);
    tree->SetBranchAddress("total_energy", &ev.total_energy);
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
    //quick and simple matching results for quick check
    tree->SetBranchAddress("mHit_E", ev.mHit_E);
    tree->SetBranchAddress("mHit_x", ev.mHit_x);
    tree->SetBranchAddress("mHit_y", ev.mHit_y);
    tree->SetBranchAddress("mHit_z", ev.mHit_z);
    tree->SetBranchAddress("mHit_gx", ev.mHit_gx);
    tree->SetBranchAddress("mHit_gy", ev.mHit_gy);
    tree->SetBranchAddress("mHit_gz", ev.mHit_gz);
    tree->SetBranchAddress("mHit_gid", ev.mHit_gid);
    // GEM part
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

void TransformToGEMFrame(float &x, float &y, float &z, int gem_id) {
    float ratio = gz[gem_id] / z;
    x = x * ratio;
    y = y * ratio;
    z = gz[gem_id];
}