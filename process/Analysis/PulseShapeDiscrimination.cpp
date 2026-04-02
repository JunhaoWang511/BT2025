// 通过电子+强子的双模板函数拟合波形，得到强子份额来鉴别粒子
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TLegend.h"
#include "../DataModel2025.hh"

Double_t ff(double *x, double *par)
{
    Double_t val = par[0] * exp(-(x[0] - par[1]) / par[2]) + par[3] * exp(-(x[0] - par[1]) / par[4]) + par[5] * exp(-(x[0] - par[1]) / par[6]) + par[7] * exp(-(x[0] - par[1]) / par[8]);

    if (x[0] >= par[1] && x[0] <= 3000)
        return val * TMath::Power((x[0] - par[1]), par[9]);
    else
        return 0;
}

// 画出电子、强子和合并后的模板函数
TF1 *DrawFunction()
{
    // 从电子、强子的平均波形拟合得到模板函数
    TF1 *fEtemp = new TF1("f_electronT", ff, -100, 3200, 10);
    TF1 *fHtemp = new TF1("f_hadronT", ff, -100, 3200, 10);
    // const double parE[10] = {6.88426, 1214.95, 67.8797, 31.5932, 49.7781, -19.9323, 43.9053, -18.5503, 59.7757, 2};
    // const double parH[10] = {0.276155, 1216.13, 85.9418, -9.43851, 43.5488, 11.8893, 47.6825, -2.72658, 56.2289, 2};
    const double parE[10] = {1.16602, 1207.99, 60.9593, -0.0858108, 24.2351, 2.07133, 56.8116, -3.15428, 58.5705, 2.7416};
    const double parH[10] = {1.9561, 1208.96, 61.5725, -0.0345951, 23.5987, 2.26606, 59.2444, -4.18932, 60.4074, 2.7612};
    fEtemp->SetParameters(parE);
    fHtemp->SetParameters(parH);
    // 将模板函数归一化
    double normalizeE, normalizeH;
    normalizeE = fEtemp->GetMaximum();
    normalizeH = fHtemp->GetMaximum();
    TF1 *fE = new TF1("f_electron", [fEtemp, normalizeE](double *x, double *par)
                      { fEtemp->SetParameter(1,par[0]); return fEtemp->Eval(x[0]) / normalizeE; }, -100, 3200, 1);
    TF1 *fH = new TF1("f_hadron", [fHtemp, normalizeH](double *x, double *par)
                      {fHtemp->SetParameter(1,par[0]); return fHtemp->Eval(x[0]) / normalizeH; }, -100, 3200, 1);
    fE->SetParameter(0, 1000);
    fH->SetParameter(0, 1000);
    // 电子和强子模板合并得到新的模板函数
    TF1 *ftotal = new TF1("f_total", [fE, fH](double *x, double *par)
                          { fE->SetParameter(0,par[0]); fH->SetParameter(0,par[0]) ; return par[1]*((1-par[2])*fE->Eval(x[0]) + par[2] * fH->Eval(x[0])); }, -100, 3200, 3);
    double par[3] = {1000, 2, 1};
    ftotal->SetParameters(par);
    // TCanvas *can = new TCanvas("can", "can", 900, 600);
    // ftotal->SetLineColor(kGreen);
    // fE->SetLineColor(kRed);
    // fH->SetLineColor(kBlue);
    // ftotal->SetNpx(1000);
    // fE->SetNpx(1000);
    // fH->SetNpx(1000);
    // ftotal->SetTitle("Template Function;Time[ns];Amplitude");
    // ftotal->GetXaxis()->SetRangeUser(1000, 1500);
    // ftotal->Draw("l");
    // fE->Draw("l same");
    // fH->Draw("l same");
    // TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg->AddEntry(ftotal, "f_{Electron}+f_{Hadron}", "l");
    // leg->AddEntry(fE, "f_{Electron}", "l");
    // leg->AddEntry(fH, "f_{Hadron}", "l");
    // leg->Draw();
    return ftotal;
}

// 模板拟合
bool TemplateFit(DataModel2025 *hit, TF1 *fun)
{
    // 画出波形和拟合结果
    bool fitresult_output = false;
    static TLatex *tex = new TLatex();
    // 模板函数参数初始化
    double init_par[3] = {0, hit->LowGainPeak, 0};
    fun->SetParameters(init_par);
    // 拟合范围：峰前12个点开始，总点数40
    int cc = 40, dd = 12;
    static TGraph *gr_wave = new TGraph(cc);
    gr_wave->SetTitle("waveform;time[ns];ADC value");
    gr_wave->SetMarkerStyle(8);
    gr_wave->SetMarkerSize(0.5);
    static TGraph *gr_wave_wrong = new TGraph(256);
    gr_wave_wrong->SetTitle("waveform;time[ns];ADC value");
    gr_wave_wrong->SetMarkerStyle(8);
    gr_wave_wrong->SetMarkerSize(0.5);

    // 寻找峰位置
    int HGMaxID = 0;
    double HGMaxAmp = 0;
    int LGMaxID = 0;
    double LGMaxAmp = 0;
    for (int j = 100; j < 120; j++)
    {
        if (hit->HAmplitude[j] >= HGMaxAmp)
        {
            HGMaxAmp = hit->HAmplitude[j];
            HGMaxID = j;
        }
        if (hit->LAmplitude[j] >= LGMaxAmp)
        {
            LGMaxAmp = hit->LAmplitude[j];
            LGMaxID = j;
        }
    }

    // 筛去震荡事例
    bool IsOsc = false;
    if ((HGMaxID < 100 || HGMaxID > 116) && HGMaxAmp < 16000)
        IsOsc = true;
    // 判断高低增益通道切换
    bool choose_HG = false, choose_LG = false;
    double HGpeak = hit->HighGainPeak, LGpeak = hit->LowGainPeak;
    double HGped = 2356.82, LGped = 2414.79, HGnoise = 23.9481, LGnoise = 3.93446;
    if (HGpeak > HGped + 6 * HGnoise && HGpeak < 16000 && !IsOsc)
        choose_HG = true;
    else if (LGpeak > LGped + 6 * LGnoise && LGpeak < 16000 && !IsOsc)
        choose_LG = true;
    if (!IsOsc)
    {
        // 默认选用LG通道
        choose_HG = false, choose_LG = true;
        // 选用HG通道
        if (choose_HG)
        {
            for (int j = 0; j < cc; j++)
            {
                gr_wave->SetPoint(j, j * 12.5, hit->HAmplitude[HGMaxID - dd + j] - HGped);
                gr_wave->SetTitle("waveform_HG");
            }
            if (!fitresult_output)
                gr_wave->Fit(fun, "Q");
            else
            {
                gr_wave->Fit(fun);
                gr_wave->Draw("ap");
                fun->Draw("same");
                tex->DrawLatexNDC(0.45, 0.8, Form("PSD value=%.2lf", fun->GetParameter(2)));
                gPad->Update();
                sleep(1);
            }
        }
        // 选用LG通道
        else if (choose_LG)
        {
            for (int j = 0; j < cc; j++)
            {
                gr_wave->SetPoint(j, j * 12.5, hit->LAmplitude[LGMaxID - dd + j] - LGped);
                gr_wave->SetTitle("waveform_LG");
            }
            if (!fitresult_output)
                gr_wave->Fit(fun, "Q");
            else
            {
                gr_wave->Fit(fun);
                gr_wave->Draw("ap");
                fun->Draw("same");
                tex->DrawLatexNDC(0.45, 0.8, Form("PSD value=%.2lf", fun->GetParameter(2)));
                gPad->Update();
                sleep(1);
            }
            fun->SetParameter(0, (LGMaxID - dd) * 12.5 + fun->GetParameter(0));
            return true;
        }
        // 波形过大或过小
        else
        {
            std::cout << "waveform too small or too large!" << std::endl;
            for (int j = 0; j < 256; j++)
            {
                gr_wave_wrong->SetPoint(j, j * 12.5, hit->HAmplitude[j]);
                gr_wave_wrong->SetTitle("waveform_HG");
                gr_wave_wrong->GetYaxis()->SetRangeUser(2000, 1.3 * HGMaxAmp);
            }
            gr_wave_wrong->Draw("ap");
            gPad->Update();
            sleep(1);
        }
    }
    // 震荡波形
    else
    {
        std::cout << "oscillation waveform!" << std::endl;
        for (int j = 0; j < 256; j++)
        {
            gr_wave_wrong->SetPoint(j, j * 12.5, hit->HAmplitude[j]);
            gr_wave_wrong->SetTitle("waveform_HG");
            gr_wave_wrong->GetYaxis()->SetRangeUser(2000, 1.3 * HGMaxAmp);
        }
        gr_wave_wrong->Draw("ap");
        gPad->Update();
        sleep(1);
    }
    return false;
}

// 粒子鉴别和时间响应
void PulseShapeDiscrimination(TString decodefile, TString recfile, TString t0file)
{
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "xy");
    // gStyle->SetPadGridX(true);
    // gStyle->SetPadGridY(true);
    // 定义合并的模板函数
    TF1 *fun = DrawFunction();
    // 设置波形、重建、T0数据存储变量
    TFile *infile_decode = new TFile(decodefile.Data(), "READ");
    TTree *tr_decode = (TTree *)infile_decode->Get("decode_data");
    DataModel2025 *Hit[25];
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            Hit[5 * i + j] = new DataModel2025();
            tr_decode->SetBranchAddress(Form("Hit_%d_%d", i + 1, j + 1), Hit[5 * i + j]);
        }

    TFile *infile_rec = new TFile(recfile.Data(), "READ");
    TTree *tr_rec = (TTree *)infile_rec->Get("rec_data");
    std::vector<int> *SeedID = 0;
    std::vector<int> *HitID = 0;
    std::vector<double> *EnergyShower = 0;
    std::vector<double> *HitEnergy = 0;
    std::vector<double> *ShowerX = 0;
    std::vector<double> *ShowerY = 0;
    std::vector<double> *HitTime = 0;
    int triggerID;
    tr_rec->SetBranchAddress("EventID", &triggerID);
    tr_rec->SetBranchAddress("HitID", &HitID);
    tr_rec->SetBranchAddress("HitEnergy", &HitEnergy);
    tr_rec->SetBranchAddress("HitTime", &HitTime);
    tr_rec->SetBranchAddress("ShowerID", &SeedID);
    tr_rec->SetBranchAddress("ShowerE5x5", &EnergyShower);
    tr_rec->SetBranchAddress("ShowerPosX5x5", &ShowerX);
    tr_rec->SetBranchAddress("ShowerPosY5x5", &ShowerY);

    TFile *infile_t0 = new TFile(t0file.Data(), "READ");
    TTree *tr_t0 = (TTree *)infile_t0->Get("CoarseCaliTree");
    Long64_t T0ID;
    double T0time[2];
    double T0Stamp;
    tr_t0->SetBranchAddress("triggerID", &T0ID);
    tr_t0->SetBranchAddress("T0time", &T0time);
    tr_t0->SetBranchAddress("triggerTS", &T0Stamp);

    // 存储拟合的能量、时间、强子份额变量
    TFile *outfile = new TFile("PSDresult.root", "recreate");
    TTree *tr_out = new TTree("PSDtree", "PSDtree");
    double time_fit, energy_fit, PSDvalue_fit;
    tr_out->Branch("Time", &time_fit);
    tr_out->Branch("Energy", &energy_fit);
    tr_out->Branch("HadronFraction", &PSDvalue_fit);
    double rec_time;
    tr_out->Branch("RecTime", &rec_time);
    // 循环所有重建事例
    TH1F *his_Hfrac = new TH1F("his", ";hadron fraction;counts", 150, -5, 10);
    his_Hfrac->SetDirectory(nullptr);
    TH1F *his_T = new TH1F("hisT", "Time Resolution;Time[ns];counts", 2000, 2450, 2550);
    his_T->SetDirectory(nullptr);
    TH1F *his_E = new TH1F("hisE", ";Energy diff[MeV];counts", 200, -100, 100);
    his_E->SetDirectory(nullptr);
    TH2F *his_ET = new TH2F("his_ET", ";Seed Energy[MeV];Time[ns]", 5000, 0, 5000, 1000, 2450, 2550);
    his_ET->SetDirectory(nullptr);
    TH2F *his_EHfrac = new TH2F("his_EHfrac", ";Seed Energy[MeV];hadron fraction", 5000, 0, 5000, 150, -5, 10);
    his_EHfrac->SetDirectory(nullptr);
    double time_flight, time_middle, time_hit, time_incident;
    double time_seed, energy_seed;
    int EventNum = std::min(tr_decode->GetEntries(), tr_t0->GetEntries());
    for (int i = 0; i < EventNum; i++)
    {
        if (i % (EventNum / 100) == 0 && i != 0)
        {
            int process = 100 * i / (EventNum);
            std::cout << "\rProcessing: >>>>>>>>>>>>>>> " << process << "% <<<<<<<<<<<<<<<";
            std::cout.flush();
        }

        tr_decode->GetEntry(i);
        tr_rec->GetEntry(i);
        tr_t0->GetEntry(i);
        if (triggerID != T0ID)
        {
            continue;
            std::cerr << "triggerID does not match." << std::endl;
        }
        // 筛选seed晶体是3-3通道,且T0测量时间有效的事例
        if (SeedID->size() != 1 || SeedID->at(0) != 326034 || (T0time[0] == 0 || T0time[1] == 0))
            continue;
        // 寻找Seed晶体并读出时间和能量信息
        for (int j = 0; j < HitID->size(); j++)
        {
            if (HitID->at(j) != 326034)
                continue;
            else
            {
                time_seed = HitTime->at(j);
                energy_seed = HitEnergy->at(j);
            }
        }
        // 两个T0之间间距5.15m，T0中点到ECAL间距4.34m
        time_flight = (T0time[1] - T0time[0]) / 1000;
        time_middle = (T0time[1] + T0time[0]) / 2 / 1000;
        time_incident = time_middle + time_flight / 5.15 * 4.08;
        // 波形拟合
        // 筛选2GeV电子簇射事例
        // if (EnergyShower->at(0) < 1800)
        //     continue;
        // 筛选2GeV质子簇射事例
        // if (EnergyShower->at(0) < 200 || time_flight < 17.5)
        //     continue;
        // 筛选2GeV MIP事例
        // if (EnergyShower->at(0) > 200)
        //     continue;
        // 筛选3GeV电子簇射事例
        // if (EnergyShower->at(0) < 2500)
        //     continue;
        // 筛选3GeV质子簇射事例
        // if (EnergyShower->at(0) < 300 || time_flight < 17.6 || time_flight > 50)
        //     continue;
        // 筛选3GeV Pion簇射事例
        // if (EnergyShower->at(0) < 200 || time_flight > 17.2 || EnergyShower->at(0) > 2000)
        //     continue;

        // if (EnergyShower->at(0) < 300)
        //     continue;
        bool valid = false;
        valid = TemplateFit(Hit[12], fun);
        if (valid)
        {
            time_fit = fun->GetParameter(0) - time_incident;
            energy_fit = fun->GetParameter(1) / 849.61 * 181.8;
            PSDvalue_fit = fun->GetParameter(2);
            rec_time = time_seed - time_incident;
            tr_out->Fill();
            his_ET->Fill(EnergyShower->at(0), time_fit);
            his_EHfrac->Fill(energy_fit, PSDvalue_fit);
            his_T->Fill(time_fit);
            his_E->Fill(energy_fit - energy_seed);
            his_Hfrac->Fill(PSDvalue_fit);
        }
    }
    his_Hfrac->Draw();
    new TCanvas();
    his_T->Draw();
    new TCanvas();
    his_E->Draw();
    new TCanvas();
    his_ET->Draw("colz");
    new TCanvas();
    his_EHfrac->Draw("colz");
    infile_decode->Close();
    infile_rec->Close();
    infile_t0->Close();
    outfile->cd();
    tr_out->Write();
    outfile->Close();
}