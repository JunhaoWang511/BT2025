// CMS时间测量方法：通过波形相邻点的比值曲线计算达峰时间
#include <iostream>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
// 获得3-3通道幅度比例-时间模板函数
void TempFun()
{
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "xy");
    double LGamp[256] = {0.0711472, 0.087709, 0.0671214, 0.0314754, 0.0211306, 0.0559104, 0.0283669, 0.0300995, 0.00153678, 0.0127988, -0.00386491, -0.00121502, -0.0207324, -0.00541917, -0.0335487, 0.00319296, -0.0505691, -0.052531, -0.0371923, -0.032504, -0.0310517, -0.0364789, -0.0355361, -0.0361222, -0.0337271, -0.0391542, -0.0259558, -0.0249366, -0.00967427, -0.011254, -0.0248601, -0.0375745, -0.075794, -0.0542637, -0.0814759, -0.0766858, -0.0729403, -0.0587481, -0.0459573, -0.0427214, -0.0457789, -0.0366827, -0.0159677, -0.0203502, 0.011907, 0.0180476, 0.0213854, 0.0685483, 0.0785363, 0.106054, 0.152351, 0.139, 0.11803, 0.148656, 0.112118, 0.148988, 0.137089, 0.121622, 0.119151, 0.115431, 0.0765489, 0.0614394, 0.0809314, 0.0782815, 0.0584074, 0.0624586, 0.0225575, 0.0260737, -0.00549561, -0.00187749, 0.0146078, 0.00811053, 0.0304307, 0.0521903, 0.067631, 0.142159, 0.146643, 0.200533, 0.167741, 0.222573, 0.23516, 0.261276, 0.209833, 0.238956, 0.226726, 0.215821, 0.192379, 0.201119, 0.171486, 0.143153, 0.13533, 0.120934, 0.108755, 0.10343, 0.0878109, 0.115966, 0.463993, 9.01377, 63.0748, 223.631, 523.665, 940.351, 1413.95, 1878.63, 2283.25, 2598.12, 2813.1, 2932.45, 2968.29, 2936.42, 2852.93, 2732.52, 2587.66, 2428.35, 2262.32, 2095.34, 1931.61, 1774.01, 1624.32, 1483.77, 1353.01, 1232.31, 1121.49, 1020.34, 928.265, 844.709, 769.004, 700.555, 638.677, 582.751, 532.074, 486.152, 444.38, 406.371, 371.801, 340.36, 311.667, 285.655, 262.012, 240.612, 221.227, 203.704, 187.917, 173.758, 161.033, 149.666, 139.505, 130.423, 122.283, 114.981, 108.662, 102.827, 97.8131, 93.1137, 89.0941, 85.3383, 82.0301, 78.7959, 75.9443, 73.1032, 70.5895, 68.0568, 65.8784, 63.6323, 61.7714, 59.8614, 58.3256, 56.7182, 55.4914, 54.1322, 53.1599, 52.0681, 51.1643, 50.263, 49.4191, 48.6157, 47.8165, 47.0956, 46.3491, 45.6829, 45.0362, 44.4149, 43.8378, 43.3093, 42.7835, 42.2868, 41.797, 41.37, 40.9049, 40.5443, 40.0825, 39.7831, 39.3857, 39.1054, 38.7173, 38.5296, 38.2223, 38.0019, 37.6844, 37.441, 37.1027, 36.8531, 36.4765, 36.165, 35.7734, 35.4786, 35.0445, 34.7277, 34.282, 34.0196, 33.6516, 33.4039, 33.0916, 32.938, 32.7029, 32.5876, 32.4502, 32.4039, 32.33, 32.3237, 32.2338, 32.2235, 32.1681, 32.0852, 32.0008, 31.8736, 31.6655, 31.4727, 31.1998, 30.9788, 30.6756, 30.5059, 30.2926, 30.1847, 30.0256, 29.9513, 29.7811, 29.735, 29.6298, 29.6139, 29.5214, 29.4721, 29.4032, 29.3872, 29.2531, 29.196, 29.0644, 28.9716, 28.7672, 28.6294, 28.3929, 28.2112, 27.9184, 27.7153, 27.3815, 27.1695};
    double HGamp[256] = {1.28935, 1.31153, 1.34297, 1.30247, 1.21638, 1.06013, 0.920056, 0.740186, 0.56983, 0.355583, 0.222914, 0.0156584, -0.160944, -0.291751, -0.460079, -0.594692, -0.687316, -0.811752, -0.891056, -0.955384, -1.03531, -1.09331, -1.13352, -1.15991, -1.07204, -1.11585, -1.16152, -1.22498, -1.34127, -1.37842, -1.39749, -1.47857, -1.50124, -1.5882, -1.66717, -1.71884, -1.75818, -1.8034, -1.82888, -1.77063, -1.75094, -1.73733, -1.73642, -1.78772, -1.81701, -1.84323, -1.81771, -1.84638, -1.78734, -1.80005, -1.82027, -1.79429, -1.72819, -1.75934, -1.73621, -1.78292, -1.77572, -1.86264, -1.89089, -1.85358, -1.83463, -1.81622, -1.79703, -1.81916, -1.77481, -1.80104, -1.81047, -1.89383, -1.88485, -1.93015, -1.85399, -1.78697, -1.65372, -1.58596, -1.4643, -1.35997, -1.1463, -1.0835, -1.10696, -1.09901, -1.1017, -1.09868, -1.11316, -1.14038, -1.17786, -1.21936, -1.23723, -1.23396, -1.18705, -1.24447, -1.26395, -1.42918, -1.54633, -1.67532, -1.70246, -1.59089, -0.331837, 17.4338, 109.214, 358.492, 811.276, 1446.88, 2196.02, 2971.81, 3699.17, 4324.21, 4816.49, 5165.72, 5376.54, 5462.64, 5442.68, 5336.61, 5164.22, 4943.01, 4688.42, 4413.13, 4127.57, 3839.8, 3556.15, 3281.11, 3018.05, 2769.17, 2535.87, 2318.69, 2117.86, 1932.94, 1763.43, 1608.43, 1467.07, 1338.14, 1220.77, 1113.85, 1016.39, 927.538, 846.569, 772.797, 705.54, 644.234, 588.544, 538, 492.216, 450.726, 413.381, 379.638, 349.265, 321.996, 297.568, 275.716, 256.198, 238.668, 223.079, 209.154, 196.812, 185.639, 175.667, 166.612, 158.44, 150.867, 143.916, 137.42, 131.54, 125.954, 120.817, 115.907, 111.403, 107.192, 103.408, 99.9075, 96.888, 94.0652, 91.6393, 89.4511, 87.4228, 85.6055, 83.9103, 82.3232, 80.7787, 79.3622, 77.9606, 76.577, 75.247, 73.9691, 72.7575, 71.5848, 70.3895, 69.2524, 68.1461, 67.1629, 66.2474, 65.3912, 64.535, 63.7176, 62.986, 62.2646, 61.5558, 60.9672, 60.3872, 59.8672, 59.3438, 58.9337, 58.4374, 57.9339, 57.3298, 56.7682, 56.0063, 55.2554, 54.4396, 53.6008, 52.8134, 51.9987, 51.1745, 50.4271, 49.7673, 49.2303, 48.7603, 48.3649, 48.063, 47.8683, 47.7565, 47.7873, 47.8447, 47.906, 47.9238, 47.9026, 47.8493, 47.7331, 47.5938, 47.3481, 47.0986, 46.7927, 46.4848, 46.0716, 45.6867, 45.2625, 44.9098, 44.5463, 44.2474, 44.0508, 43.8396, 43.6057, 43.4902, 43.3608, 43.232, 43.0208, 42.7709, 42.5461, 42.3078, 42.087, 41.8113, 41.5514, 41.2356, 40.8725, 40.4802, 40.0331, 39.7087, 39.7692};
    double LGratio[43], HGratio[43];
    TGraph *gr_LG = new TGraph(43);
    TGraph *gr_HG = new TGraph(43);
    gr_LG->SetMarkerStyle(20);
    gr_LG->SetMarkerColor(kBlue);
    gr_HG->SetMarkerStyle(20);
    gr_HG->SetMarkerColor(kRed);
    for (int i = 0; i < 43; i++)
    {
        LGratio[i] = LGamp[96 + i] / LGamp[97 + i];
        HGratio[i] = HGamp[96 + i] / HGamp[97 + i];
        gr_LG->SetPoint(i, LGratio[i], (i + 96) * 12.5 - 1348.35);
        gr_HG->SetPoint(i, HGratio[i], (i + 96) * 12.5 - 1363.38);
    }
    TF1 *fLG = new TF1("fLG", "pol5", 0, 1.2);
    TF1 *fHG = new TF1("fHG", "pol5", 0, 1.2);
    fLG->SetLineColor(kBlue);
    fHG->SetLineColor(kRed);
    gr_HG->Fit(fHG, "R", "", 0.2, 1.05);
    gr_LG->Fit(fLG, "R", "", 0.2, 1.05);
    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    gr_HG->GetYaxis()->SetRangeUser(-200, 200);
    gr_HG->SetTitle(";R=A(T)/A(T+25 ns);T-T_{max} [ns]");
    gr_HG->Draw("ap");
    gr_LG->Draw("p same");

    double waveHG[256], waveLG[256], scaleLG, scaleHG;
    scaleHG = *std::max_element(HGamp, HGamp + 256);
    scaleLG = *std::max_element(LGamp, LGamp + 256);

    TGraph *gr_waveHG = new TGraph(256);
    TGraph *gr_waveLG = new TGraph(256);
    gr_waveHG->SetMarkerStyle(21);
    gr_waveLG->SetMarkerStyle(21);
    gr_waveHG->SetMarkerColor(kRed);
    gr_waveLG->SetMarkerColor(kBlue);
    for (int i = 0; i < 256; i++)
    {
        waveHG[i] = HGamp[i] / scaleHG * 0.8;
        waveLG[i] = LGamp[i] / scaleLG * 0.8;
        gr_waveHG->SetPoint(i, waveHG[i], 12.5 * i - 1363.38);
        gr_waveLG->SetPoint(i, waveLG[i], 12.5 * i - 1348.35);
    }
    gr_waveHG->Draw("p same");
    gr_waveLG->Draw("p same");

    TLegend *leg = new TLegend(0.2, 0.4, 0.4, 0.7);
    leg->AddEntry(gr_HG, "HG ratio curve", "p");
    leg->AddEntry(gr_LG, "LG ratio curve", "p");
    leg->AddEntry(fHG, "HG ratio function", "l");
    leg->AddEntry(fLG, "LG ratio function", "l");
    leg->AddEntry(gr_waveHG, "HG waveform", "p");
    leg->AddEntry(gr_waveLG, "LG waveform", "p");
    leg->Draw();

    double HGpars[6], LGpars[6];
    fHG->GetParameters(HGpars);
    fLG->GetParameters(LGpars);
    std::cout << "HG parameters:" << std::endl;
    for (int i = 0; i < 6; i++)
    {
        if (i != 5)
            std::cout << HGpars[i] << ',';
        else
            std::cout << HGpars[i] << std::endl;
    }
    std::cout << "LG parameters: " << std::endl;
    for (int i = 0; i < 6; i++)
    {
        if (i != 5)
            std::cout << LGpars[i] << ',';
        else
            std::cout << LGpars[i] << std::endl;
    }
}

// 计算3-3通道高低通道波形的时间
double CalculateTime(std::array<double, 256> HGwave, std::array<double, 256> LGwave)
{
    // 波形前后两个采样点幅度和它们的比例及误差
    double Amp, AmpNext;
    double ratio[13], sigma2[13], time[13];
    // 计算得到的时间
    double timeSum = 0, sigmaSum = 0, timeAve;
    // 3-3通道的基线、噪声
    static const double HGped = 2356.82, LGped = 2414.79, HGNoise = 23.9481, LGNoise = 3.93446;
    // 时间取值函数
    static TF1 *fHG = new TF1("fHG", "pol5", 0.4, 1.05);
    static TF1 *fLG = new TF1("fLG", "pol5", 0.4, 1.05);
    double parHG[6] = {-518.454, 3390.34, -11791.6, 20203.2, -16713.5, 5427.55};
    double parLG[6] = {-361.488, 2168.32, -7677.99, 13489.5, -11391.4, 3769.26};
    fHG->SetParameters(parHG);
    fLG->SetParameters(parLG);
    // TCanvas *can = new TCanvas("can", "can", 900, 600);
    // can->Divide(2, 1);
    // can->cd(1);
    // fHG->Draw();
    // can->cd(2);
    // fLG->Draw();

    double HGmax = 0, LGmax = 0;
    int HGpos, LGpos;
    for (int i = 90; i < 120; i++)
    {
        if (HGwave[i] > HGmax)
        {
            HGmax = HGwave[i];
            HGpos = i;
        }
        if (LGwave[i] > LGmax)
        {
            LGmax = LGwave[i];
            LGpos = i;
        }
    }
    // 判断振荡波形
    bool osc = false;
    if (HGmax < 16000)
    {
        if (HGpos < 103 || HGpos > 113)
            osc = true;
    }
    // 高增益通道
    if (HGmax > (HGped + 6 * HGNoise) && HGmax < 16000 && !osc)
    {
        // 扣除基线后算出比例及对应误差(峰前10个点开始，总共13个点)
        // 计算每个采样点的外推时间
        for (int i = 0; i < 13; i++)
        {
            Amp = HGwave[HGpos - 10 + i] - HGped;
            AmpNext = HGwave[HGpos - 9 + i] - HGped;
            ratio[i] = Amp / AmpNext;
            // sigma2[i] = TMath::Power(HGNoise, 2) * (1 / TMath::Power(Amp, 2) + TMath::Power(Amp, 2) / TMath::Power(AmpNext, 4)) * TMath::Power(fHG->Eval(ratio[i]), 2);
            // sigma2[i] = TMath::Power(HGNoise, 2) * (1 / TMath::Power(Amp, 2) + TMath::Power(Amp, 2) / TMath::Power(AmpNext, 4));
            sigma2[i] = 1;
            if (ratio[i] > 0.4 && ratio[i] < 1.05)
            {
                time[i] = (HGpos - 10 + i) * 12.5 - fHG->Eval(ratio[i]);
            }
            else
                time[i] = 0;
        }
        // 对可用采样点的外推时间取平均
        for (int i = 3; i < 8; i++)
        {
            if (time[i] != 0)
            {
                timeSum += time[i] / sigma2[i];
                sigmaSum += 1 / sigma2[i];
            }
        }
        if (sigmaSum != 0)
            timeAve = timeSum / sigmaSum;
        else
            timeAve = 777777;
    }
    // 低增益通道
    else if (LGmax > (LGped + 6 * LGNoise) && LGmax < 16000 && !osc)
    {
        for (int i = 0; i < 13; i++)
        {
            Amp = LGwave[LGpos - 10 + i] - LGped;
            AmpNext = LGwave[LGpos - 9 + i] - LGped;
            ratio[i] = Amp / AmpNext;
            // sigma2[i] = TMath::Power(LGNoise, 2) * (1 / TMath::Power(Amp, 2) + TMath::Power(Amp, 2) / TMath::Power(AmpNext, 4)) * TMath::Power(fLG->Eval(ratio[i]), 2);
            // sigma2[i] = TMath::Power(LGNoise, 2) * (1 / TMath::Power(Amp, 2) + TMath::Power(Amp, 2) / TMath::Power(AmpNext, 4)) ;
            sigma2[i] = 1;
            if (ratio[i] > 0.4 && ratio[i] < 1.05)
            {
                time[i] = (LGpos - 10 + i) * 12.5 - fLG->Eval(ratio[i]);
            }
            else
                time[i] = 0;
        }
        for (int i = 3; i < 8; i++)
        {
            if (time[i] != 0)
            {
                timeSum += time[i] / sigma2[i];
                sigmaSum += 1 / sigma2[i];
            }
        }
        if (sigmaSum != 0)
            timeAve = timeSum / sigmaSum;
        else
            timeAve = 777777;
    }
    // 震荡波形和幅度过高或过低事例
    else
    {
        if (osc)
            timeAve = 888888;
        else
            timeAve = 999999;
    }
    return timeAve;
}

// 输入decode文件，处理3-3通道
void TimeRatioMethod(TString filename)
{
    // 读入事例波形
    TFile *infile = new TFile(filename.Data(), "READ");
    TTreeReader reader("decode_data", infile);
    TTreeReaderArray<double> HGamp(reader, "Hit_3_3.HAmplitude");
    TTreeReaderArray<double> LGamp(reader, "Hit_3_3.LAmplitude");
    TTreeReaderValue<int> TriggerID(reader, "TriggerID");
    std::array<double, 256> HGwave, LGwave;
    // 保存计算得到的时间信息
    TFile *outfile = new TFile("time_ratio.root", "RECREATE");
    TTree *tr_out = new TTree("Tr_time", "Tr_time");
    double time;
    Long64_t Tid;
    tr_out->Branch("time", &time, "time/D");
    tr_out->Branch("triggerID", &Tid, "triggerID/L");
    // 处理事例波形
    int nEntries = reader.GetEntries();
    int interval = nEntries / 20;
    for (int i = 0; i < reader.GetEntries(); i++)
    {
        int progress = static_cast<float>(i + 1) / nEntries * 100;
        if ((i + 1) % interval == 0)
        {
            std::cout << "Progress: " << progress + 1 << "%\r" << std::endl;
            std::cout.flush();
        }
        reader.Next();
        Tid = *TriggerID;
        for (int j = 0; j < 256; j++)
        {
            HGwave[j] = HGamp[j];
            LGwave[j] = LGamp[j];
        }
        time = CalculateTime(HGwave, LGwave);
        tr_out->Fill();
    }
    infile->Close();
    outfile->cd();
    tr_out->Write();
    outfile->Close();
}

// 计算时间分辨
void TimeResolution(std::string timefile, std::string ecalfile, std::string t0file)
{
    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42);
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLegendTextSize(0.06);
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    TFile *InfileTime = new TFile(timefile.data(), "READ");
    TTree *TrTime = (TTree *)InfileTime->Get("Tr_time");
    double time;
    Long64_t Tid;
    TrTime->SetBranchAddress("time", &time);
    TrTime->SetBranchAddress("triggerID", &Tid);

    TFile *InfileECAL = new TFile(ecalfile.data(), "READ");
    TTree *TrECAL = (TTree *)InfileECAL->Get("rec_data");
    std::vector<int> *SeedID = 0;
    std::vector<int> *HitID = 0;
    std::vector<double> *EnergyShower = 0;
    std::vector<double> *HitEnergy = 0;
    std::vector<double> *ShowerX = 0;
    std::vector<double> *ShowerY = 0;
    std::vector<double> *HitTime = 0;
    int triggerID;
    TrECAL->SetBranchAddress("EventID", &triggerID);
    TrECAL->SetBranchAddress("ShowerID", &SeedID);
    TrECAL->SetBranchAddress("ShowerE5x5", &EnergyShower);
    TrECAL->SetBranchAddress("HitID", &HitID);
    TrECAL->SetBranchAddress("HitEnergy", &HitEnergy);
    TrECAL->SetBranchAddress("ShowerPosX5x5", &ShowerX);
    TrECAL->SetBranchAddress("ShowerPosY5x5", &ShowerY);
    TrECAL->SetBranchAddress("HitTime", &HitTime);

    TFile *InfileT0 = new TFile(t0file.data(), "READ");
    TTree *TrT0 = (TTree *)InfileT0->Get("CoarseCaliTree");
    Long64_t T0ID;
    double T0time[2];
    double T0Stamp;
    TrT0->SetBranchAddress("triggerID", &T0ID);
    TrT0->SetBranchAddress("T0time", &T0time);
    TrT0->SetBranchAddress("triggerTS", &T0Stamp);
    double time_flight, time_middle, time_hit, time_incident;
    double posT01, posT02, posECAL;

    // determine electron beam energy
    TH1F *energy_test = new TH1F("energy", "energy", 1000, 0, 10000);
    double maxenergy = 0, maxheight = 0;
    for (int i = 0; i < TrECAL->GetEntries(); i++)
    {
        TrECAL->GetEntry(i);
        for (size_t j = 0; j < EnergyShower->size(); j++)
        {
            if (EnergyShower->at(j) > 10000)
                continue;
            energy_test->Fill(EnergyShower->at(j));
            if ((EnergyShower->at(j) > maxenergy) && (energy_test->GetBinContent(static_cast<int>(EnergyShower->at(j) / energy_test->GetBinWidth(0))) > 50))
                maxenergy = EnergyShower->at(j);
        }
    }
    int startbin = 0, maxbin = 0;
    startbin = maxenergy / 2 / 10;
    for (int i = startbin; i < energy_test->GetNbinsX(); i++)
    {
        if (energy_test->GetBinContent(i) > maxheight)
        {
            maxheight = energy_test->GetBinContent(i);
            maxbin = i;
        }
    }
    double energy = std::round(maxbin * energy_test->GetBinWidth(0) / 100) * 100;
    double SeedCut = energy / 2;

    double TimeSeed, EnergySeed;
    TH2F *his_ET = new TH2F("his_ET", ";Seed Energy[MeV];Time[ns]", 3500, 0, 3500, 1000, 2550, 2750);
    TH1F *his_T = new TH1F("hisT", "Time Resolution;Time[ns];counts", 2000, 2550, 2750);
    TH1F *his_TMIP = new TH1F("hisTMIP", "Time Resolution;Time[ns];counts", 2000, 2550, 2750);
    TH1F *his_E = new TH1F("his", ";Energy[MeV];counts", 3000, 0, 3000);
    his_ET->SetDirectory(nullptr);
    his_T->SetDirectory(nullptr);
    his_E->SetDirectory(nullptr);
    int EventSelected = 0;
    for (int i = 0; i < std::min(TrT0->GetEntries(), TrECAL->GetEntries()); i++)
    {
        TrT0->GetEntry(i);
        TrECAL->GetEntry(i);
        TrTime->GetEntry(i);
        if (triggerID != T0ID)
        {
            continue;
            std::cerr << "triggerID not match." << std::endl;
        }
        if (SeedID->size() != 1 || SeedID->at(0) != 326034 || (T0time[0] == 0 || T0time[1] == 0))
            continue;
        for (int j = 0; j < HitID->size(); j++)
        {
            if (HitID->at(j) != 326034)
                continue;
            else
            {
                TimeSeed = HitTime->at(j);
                EnergySeed = HitEnergy->at(j);
            }
        }
        // 两个T0之间间距5.15m，T0中点到ECAL间距4.34m
        time_flight = (T0time[1] - T0time[0]) / 1000;
        time_middle = (T0time[1] + T0time[0]) / 2 / 1000;
        time_incident = time_middle + time_flight / 5.15 * 4.34;

        // 筛选MIP pion事例
        if (time_flight < 18 && EnergyShower->at(0) < 250)
        // if (time_flight < 17.2 && EnergyShower->at(0) < 250)
        {
            his_ET->Fill(EnergyShower->at(0), time - time_incident);
            his_TMIP->Fill(time - time_incident);
        }

        // 去除低能强子事例
        if (EnergySeed < SeedCut)
            continue;

        his_ET->Fill(EnergySeed, time - time_incident);
        // his_ET->Fill(EnergyShower->at(0), TimeSeed - time_incident);
        // his_T->Fill(time - time_incident);
        his_T->Fill(time - time_middle);
        if (time < 10000)
            EventSelected++;
        his_E->Fill(EnergyShower->at(0));
    }

    TCanvas *can_ET = new TCanvas();
    his_ET->Draw("colz");

    TCanvas *can_T = new TCanvas();
    can_T->cd();

    // TLatex *tex = new TLatex();
    // tex->SetTextSize(0.04);
    // tex->SetTextFont(42);
    // tex->SetTextAlign(22);

    // double time_peak = his_TMIP->GetMean();
    // double time_error = his_TMIP->GetStdDev();
    // double scale_factor = his_TMIP->GetMaximum();
    // his_TMIP->GetXaxis()->SetRangeUser(time_peak - 40, time_peak + 10);
    // his_TMIP->SetLineColor(kBlue);
    // TH1F *his1 = (TH1F *)his_TMIP->Clone("hisTMIP_norm");
    // his1->SetDirectory(nullptr);
    // his1->Scale(1 / scale_factor);
    // his1->Draw("hist");
    // TF1 *ftimeMIP = new TF1("f1", "gaus", 2550, 2750);
    // ftimeMIP->SetLineColor(kBlue);
    // ftimeMIP->SetNpx(1000);
    // ftimeMIP->SetParameter(1, time_peak);
    // ftimeMIP->SetParameter(2, time_error);
    // his1->Fit(ftimeMIP, "QR", "", time_peak - 6, time_peak + 6);
    // ftimeMIP->Draw("same");
    // tex->SetTextColor(kBlue);
    // tex->DrawLatexNDC(0.3, 0.6, Form("Mean T_{MIP}=%.2f ns", ftimeMIP->GetParameter(1)));

    // time_peak = his_T->GetMean();
    // time_error = his_T->GetStdDev();
    // his_T->GetXaxis()->SetRangeUser(time_peak - 6, time_peak + 6);
    // his_T->SetLineColor(kRed);
    // scale_factor = his_T->GetMaximum();
    // TH1F *his2 = (TH1F *)his_T->Clone("hisT_norm");
    // his2->SetDirectory(nullptr);
    // his2->Scale(1 / scale_factor);
    // his2->Draw("hist same");
    // TF1 *ftime = new TF1("f2", "gaus", 2550, 2750);
    // ftime->SetLineColor(kRed);
    // ftime->SetNpx(1000);
    // ftime->SetParameter(1, time_peak);
    // ftime->SetParameter(2, time_error);
    // his2->Fit(ftime, "QR", "", time_peak - 3, time_peak + 3);
    // ftime->Draw("same");
    // tex->SetTextColor(kRed);
    // tex->DrawLatexNDC(0.3, 0.5, Form("Mean T_{Eshower}=%.2f ns", ftime->GetParameter(1)));

    // tex->SetTextColor(kBlack);
    // tex->SetTextSize(0.05);
    // tex->DrawLatexNDC(0.3, 0.4, Form("T_{diff}=%.2f ps", (ftimeMIP->GetParameter(1) - ftime->GetParameter(1)) * 1000));
    // tex->DrawLatexNDC(0.3, 0.3, Form("#sigma_{Tdiff}=%.2f ps", TMath::Sqrt(TMath::Power(ftime->GetParError(1) * 1000, 2) + TMath::Power(ftimeMIP->GetParError(1) * 1000, 2))));
    // TLegend *leg = new TLegend(0.15, 0.65, 0.35, 0.85);
    // leg->SetTextSize(0.05);
    // leg->AddEntry(his1, "MIP #pi", "l");
    // leg->AddEntry(his2, "Eshower", "l");
    // leg->Draw();
    // can_T->SaveAs("time_differenceCMS.png");

    double time_peak = his_T->GetMean();
    double time_error = his_T->GetStdDev();
    his_T->GetXaxis()->SetRangeUser(time_peak - 3, time_peak + 3);
    his_T->Draw();
    TF1 *ftime = new TF1("f2", "gaus", 2550, 2750);
    ftime->SetNpx(1000);
    ftime->SetParameter(1, time_peak);
    ftime->SetParameter(2, time_error);
    his_T->Fit(ftime, "QSR", "", time_peak - 3, time_peak + 3);
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.08);
    tex->SetTextFont(42);
    tex->SetTextColor(kRed);
    tex->SetTextAlign(22);
    tex->DrawLatexNDC(0.25, 0.6, Form("#sigma=%.0lf ps", ftime->GetParameter(2) * 1000));
    can_T->SaveAs("time_resolution.png");
    std::cout << "beam energy: " << energy << " MeV." << std::endl;
    std::cout << "total Event number: " << EventSelected << std::endl;
    std::cout << "time resolution: " << ftime->GetParameter(2) * 1000 << " ps" << std::endl;
    new TCanvas();
    his_E->Draw();
    InfileECAL->Close();
    InfileT0->Close();
}
