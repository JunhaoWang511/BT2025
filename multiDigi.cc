#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <TString.h>
#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>
#include <TChain.h>
#include <TLatex.h>
#include <TStyle.h>
#include "DataModel2025.hh"
#include "data_model.hh"
#include "Parameter.hh"
#include "Decode2025.hh"
#include <eigen3/Eigen/Dense>
using namespace Eigen;

using namespace std;
using namespace TMath;
#define _Nsample 40
// 判断输入文件是不是txt
bool isTextFile(const string &fileName)
{
    size_t dotPos = fileName.find_last_of('.');
    if (dotPos == std::string::npos)
        return false;

    std::string extension = fileName.substr(dotPos + 1);
    if (extension == "txt")
        return true;
    else
        return false;
}
// 将3-3通道的低增益拼接到高增益
void CombineWave(double *LGamp, double *HGamp, double *CombineAmp)
{
    // 滤波向量
    static const double FilterVec_H2L[8] = {0.817645978068048, 1.21694377160809, -1.92417373282404, 2.21537086557430, -2.54878567302396, 1.94923847264921, -0.991979185381939, 0.197343964786675};
    static const double FilterVec_L2H[8] = {1.13899632255057, -0.992992540406039, 1.34191143401557, -0.464510047618861, -1.15199803190155, 2.52488621629364, -2.07852584018316, 0.754333693894717};
    // 滤波向量，不用乘上高低增益比
    static const double FilterVec_L2H1[8] = {10.4656121588991, -5.72952676356295, 6.00640823733476, 0.487130016570159, -7.99434621285014, 15.0185528533186, -12.5446979502488, 5.07359496609259};
    static const double HGped = 2356.82, LGped = 2414.79;
    static const double GainRatio = 10.0951;
    // 扣除高低增益基线
    for (int i = 0; i < _Npoints; i++)
    {
        LGamp[i] -= LGped;
        HGamp[i] -= HGped;
    }
    double FilterAmp[_Npoints];
    double temp;
    // 对低增益做滤波,并乘上比例系数
    // for (int i = 7; i < _Npoints; i++)
    // {
    //     temp = 0;
    //     for (int k = 0; k < 8; k++)
    //     {
    //         temp += FilterVec_L2H[k] * LGamp[i - k];
    //     }
    //     FilterAmp[i] = temp * GainRatio;
    // }
    // 对低增益做滤波,不用乘上比例系数
    for (int i = 7; i < _Npoints; i++)
    {
        temp = 0;
        for (int k = 0; k < 8; k++)
        {
            temp += FilterVec_L2H1[k] * LGamp[i - k];
        }
        FilterAmp[i] = temp;
    }
    // 拼接高增益和滤波波形
    for (int i = 0; i < _Npoints; i++)
    {
        if (HGamp[i] < 16000 - HGped)
            CombineAmp[i] = HGamp[i];
        else
        {
            CombineAmp[i] = FilterAmp[i];
        }
    }
}

// 模板函数形式是4个指数乘上一个幂函数
Double_t f4pow(double *x, double *par)
{
    Double_t val = par[0] * exp(-(x[0] - par[1]) / par[2]) + par[3] * exp(-(x[0] - par[1]) / par[4]) + par[5] * exp(-(x[0] - par[1]) / par[6]) + par[7] * exp(-(x[0] - par[1]) / par[8]);

    if (x[0] >= par[1] && x[0] <= 3000)
        return val * TMath::Power((x[0] - par[1]), 2);
    else
        return 0;
}

// 采样1250个模板函数
void TemplateSample(double sample[][_Nsample])
{
    // 3-3LG通道的模板函数参数
    double pars[9] = {2.93485, 867.419, 66.7128, -9.86, 47.2225, 16.8887, 51.7385, -9.96349, 59.2686};
    TF1 *ff = new TF1("ff", f4pow, 0, 3200, 9);
    ff->SetParameters(pars);
    ff->SetNpx(320000);
    double max = ff->GetMaximum();
    TF1 *fun = new TF1("fun", [ff, max](double *x, double *p)
                       { return ff->Eval(x[0]) / max; }, 0, 3200, 9);
    // fun->SetLineColor(kRed);
    // fun->Draw();
    // 采样起始点从70个点（862.5 ns）到71个点（875 ns），采样1250份（间隔10 ps），采样长度为40个点
    double sampletime;
    for (int i = 0; i < 1250; i++)
        for (int j = 0; j < _Nsample; j++)
        {
            sampletime = 862.5 + 0.01 * i + 12.5 * j;
            sample[i][j] = fun->Eval(sampletime);
        }
}
void test()
{
    double sampletemp[1250][_Nsample];
    TemplateSample(sampletemp);
    // 分别是模板向量，幅度向量和时幅向量
    const double (&tempVec)[1250][_Nsample] = sampletemp;
    const double ampVec[20] = {0.021883354, 0.035369323, 0.048912107, 0.061108757, 0.071082743, 0.078427294, 0.083095147, 0.085283467, 0.085331614, 0.083644033, 0.080631583, 0.076677201, 0.072116958, 0.067228187, 0.062230281, 0.057288467, 0.052520386, 0.048003958, 0.043785325, 0.039886313};
    const double amptimeVec[20] = {-15.689054, -15.757584, -14.150650, -11.493371, -8.3474713, -5.1451263, -2.1800620, 0.37561494, 2.4447146, 4.0175842, 5.1282024, 5.8351812, 6.2078214, 6.3165472, 6.2269922, 5.9966289, 5.6734446, 5.2958130, 4.8931847, 4.4871620};
    TGraph *gr = new TGraph(1250);
    double amp, time;
    double wavetmp[20];
    for (int i = 0; i < 1250; i++)
    {
        amp = std::inner_product(ampVec, ampVec + 20, tempVec[i], 0.);
        time = std::inner_product(amptimeVec, amptimeVec + 20, tempVec[i], 0.);
        time /= amp;
        gr->SetPoint(i, time * 100 + 625, amp);
    }
    gr->Draw("ap");
    gPad->SaveAs("template_ATcor.png");
}
bool OnePulseFit(int cc, VectorXd timingvec, double *time, double *amp, double *pedestal, double *ochi2, TF1 *f1)
{

    VectorXd vec_ones(cc);
    MatrixXd mat_noi = MatrixXd::Zero(cc, cc);

    for (int i = 0; i < cc; i++)
    {
        vec_ones(i) = 1;
        mat_noi(i, i) = 1;
    }

    double t_trg = 0.0;
    // double t_trg = -400;
    VectorXd x(3);
    double fit_a, chi2, pes;
    double deltachi2, tchi2;
    chi2 = 9999;

    for (int l = 0; l < 6; l++)
    {
        VectorXd est_data(cc);
        VectorXd d_est_data(cc);
        for (int i = 0; i < cc; i++)
        {
            // double i1 = (i * 25 + t_trg);
            double i1 = (i * 12.5 + t_trg);
            est_data(i) = f1->Eval(i1);
            d_est_data(i) = f1->Derivative(i1);
        }
        MatrixXd A(3, 3);
        A(0, 0) = est_data.transpose() * mat_noi * est_data;
        A(0, 1) = est_data.transpose() * mat_noi * d_est_data;
        A(0, 2) = est_data.transpose() * mat_noi * vec_ones;
        A(1, 0) = d_est_data.transpose() * mat_noi * est_data;
        A(1, 1) = d_est_data.transpose() * mat_noi * d_est_data;
        A(1, 2) = d_est_data.transpose() * mat_noi * vec_ones;
        A(2, 0) = vec_ones.transpose() * mat_noi * est_data;
        A(2, 1) = vec_ones.transpose() * mat_noi * d_est_data;
        A(2, 2) = vec_ones.transpose() * mat_noi * vec_ones;

        VectorXd b(3);
        b(0) = est_data.transpose() * mat_noi * timingvec;
        b(1) = d_est_data.transpose() * mat_noi * timingvec;
        b(2) = vec_ones.transpose() * mat_noi * timingvec;
        x = A.colPivHouseholderQr().solve(b);
        double dt;
        dt = x(1) / x(0);
        t_trg = t_trg + dt;
        fit_a = x(0);
        pes = x(2);
        tchi2 = (timingvec - x(0) * est_data - x(1) * d_est_data - x(2) * vec_ones)
                    .transpose() *
                mat_noi *
                (timingvec - x(0) * est_data - x(1) * d_est_data - x(2) * vec_ones);
        deltachi2 = chi2 - tchi2;
        chi2 = tchi2;
    }

    *time = -t_trg;
    *amp = fit_a;
    *pedestal = pes;
    *ochi2 = chi2;

    return 1;
}

double TemplateFit(double *LGamp)
{
    // 3-3LG通道的模板函数参数
    double pars[9] = {2.93485, 867.419, 66.7128, -9.86, 47.2225, 16.8887, 51.7385, -9.96349, 59.2686};
    TF1 *ff = new TF1("ff", f4pow, 0, 3200, 9);
    ff->SetParameters(pars);
    int MaxID = 0;
    double MaxAmp = 0;
    for (int i = 90; i < 125; i++)
    {
        if (LGamp[i] > MaxAmp)
        {
            MaxAmp = LGamp[i];
            MaxID = i;
        }
    }
    int cc = 10;
    int dd = 12;
    VectorXd lgtimingvec(cc);
    for (int j = 0; j < cc; j++)
    {
        lgtimingvec(j) = LGamp[MaxID - dd + j];
    }
    ff->SetParameter(1, 0);
    double time, amp, pedestal, chi2;
    OnePulseFit(cc, lgtimingvec, &time, &amp, &pedestal, &chi2, ff);
    return amp * ff->GetMaximum();
}
// int main(int argc, char const *argv[])
// {
//     gStyle->SetPadGridX(true);
//     gStyle->SetPadGridY(true);
//     vector<string> datafiles;

//     if (isTextFile(argv[1]))
//     {
//         ifstream filestream(argv[1]);
//         if (filestream.is_open())
//         {
//             string filename;
//             while (getline(filestream, filename))
//                 datafiles.push_back(filename);
//         }
//         else
//             cout << "error in reading filelist" << endl;
//     }
//     else
//         datafiles.push_back(argv[1]);

//     TChain *tree = new TChain("decode_data");
//     for (size_t i = 0; i < datafiles.size(); i++)
//         tree->AddFile(datafiles.at(i).data());
//     int TriggerID;
//     Long64_t TimeCode;
//     tree->SetBranchAddress("TriggerID", &TriggerID);
//     tree->SetBranchAddress("TimeCode", &TimeCode);
//     vector<decode_data_col *> Hit(25);
//     string Name[25];
//     for (int i = 0; i < 25; i++)
//     {
//         Name[i] = to_string(i / 5 + 1) + "_" + to_string(i % 5 + 1);
//     }
//     for (int i = 0; i < 25; i++)
//         Hit[i] = new decode_data_col();
//     for (int i = 0; i < 25; i++)
//     {
//         string name = "Hit_" + Name[i];
//         tree->SetBranchAddress(name.data(), Hit[i]);
//     }

//     TString outputfile;
//     if (argc >= 3)
//         outputfile = argv[2];
//     else
//     {
//         outputfile = "ECALDigiOnline.root";
//         cout << "Auto save file as ECALDigiOnline.root..." << endl;
//     }

//     // 输出文件
//     TFile *fout = new TFile(outputfile.Data(), "recreate");
//     TTree *mTree = new TTree("decode_data", "decode_data");
//     DataModel2025 *mHit[25];
//     for (int i = 0; i < 25; i++)
//         mHit[i] = new DataModel2025();
//     Long64_t mEventID;
//     int mTriggerID;
//     Long64_t mTimeCode;
//     float mTime[_Npoints];
//     float mTemperature[10];
//     mTree->Branch("TriggerID", &mTriggerID, "TriggerID/I");
//     mTree->Branch("TimeCode", &mTimeCode, "TimeCode/L");
//     mTree->Branch("Time", &mTime, Form("Time[%d]/F", _Npoints));

//     for (int i = 0; i < 25; i++)
//     {
//         Name[i] = std::to_string(i / 5 + 1) + "_" + std::to_string(i % 5 + 1);
//         std::string name = "Hit_" + Name[i];
//         std::string leaf_list = "";

//         leaf_list += "CrystalID/L";
//         leaf_list += ":Temperature1/D";
//         leaf_list += ":Temperature2/D";
//         leaf_list += Form(":LAmplitude[%d]/D", _Npoints);
//         leaf_list += Form(":HAmplitude[%d]/D", _Npoints);
//         leaf_list += Form(":LNoise[%d]/D", _Nnoise);
//         leaf_list += Form(":HNoise[%d]/D", _Nnoise);
//         leaf_list += ":LowGainPedestal/D";
//         leaf_list += ":HighGainPedestal/D";
//         leaf_list += ":LowGainPeak/D";
//         leaf_list += ":HighGainPeak/D";

//         mTree->Branch(name.data(), &mHit[i]->CrystalID, leaf_list.data());
//         mTree->Branch(Form("%s_TimeStamp", name.data()), &mHit[i]->TimeStamp, "TimeStamp/D");
//         mTree->Branch(Form("%s_CoarseTime", name.data()), &mHit[i]->CoarseTime);
//         mTree->Branch(Form("%s_FineTime", name.data()), &mHit[i]->FineTime);
//         mTree->Branch(Form("%s_Amplitude", name.data()), &mHit[i]->Amplitude);
//     }
//     double timestamp, coarsetime, finetime, amplitude;
//     int nEntries = tree->GetEntries();
//     int interval = nEntries / 20;
//     double sampletemp[1250][_Nsample];
//     TemplateSample(sampletemp);
//     const double LGped[25] = {2297.49, 2305.03, 2338.2, 2377.62, 2354.75, 2365.34, 2347.87, 2365.32, 2349.4, 2326.33, 2347.65, 2354.71, 2414.79, 2379.85, 2381.56, 2308.86, 2323.2, 2365.7, 2365.65, 2318.67, 2337, 2363.79, 2352.87, 2328.34, 2329.16};
//     // 分别是模板向量，幅度向量和时幅向量
//     const double (&tempVec)[1250][_Nsample] = sampletemp;
//     // 延迟50ns的提取向量
//     // const double ampVec[20] = {0.021883354, 0.035369323, 0.048912107, 0.061108757, 0.071082743, 0.078427294, 0.083095147, 0.085283467, 0.085331614, 0.083644033, 0.080631583, 0.076677201, 0.072116958, 0.067228187, 0.062230281, 0.057288467, 0.052520386, 0.048003958, 0.043785325, 0.039886313};
//     // const double amptimeVec[20] = {-15.689054, -15.757584, -14.150650, -11.493371, -8.3474713, -5.1451263, -2.1800620, 0.37561494, 2.4447146, 4.0175842, 5.1282024, 5.8351812, 6.2078214, 6.3165472, 6.2269922, 5.9966289, 5.6734446, 5.2958130, 4.8931847, 4.4871620};
//     const double ampVec[20] = {-6.38357040248042e-05, -0.00330697007822949, -0.00398696499554231, 0.00179119025323171, 0.0134049745184865, 0.0285970474714096, 0.0449272265503205, 0.0604065482502062, 0.0736885867430082, 0.0840412605666282, 0.0912268879037425, 0.0953595590576070, 0.0967738233148716, 0.0959185309201458, 0.0932787330236419, 0.0893231471006163, 0.0844725072062953, 0.0790836856326403, 0.0734449149428907, 0.0677782432443050};
//     const double amptimeVec[20] = {-0.0472543391531407, -3.43817924184712, -8.37953139807132, -12.0715530570965, -13.7001056339471, -13.3994007185142, -11.6731104324930, -9.09225696548941, -6.15481131398440, -3.23549247413817, -0.582418053622125, 1.66481780411593, 3.44986840988616, 4.77456618906330, 5.67744119483497, 6.21693646975686, 6.45923915076050, 6.47009666969433, 6.30980287939698, 6.03054868872617};
//     double wavetmp[20], timetmp, amptmp;
//     // 原始波形
//     TGraph *gr = new TGraph(_Npoints);
//     gr->SetTitle(";;ADC value");
//     gr->GetYaxis()->SetTitleSize(0.06);
//     gr->GetYaxis()->SetTitleOffset(0.5);
//     gr->GetYaxis()->SetLabelSize(0.04);
//     gr->SetMarkerStyle(8);
//     gr->SetMarkerSize(0.4);
//     // 剩余波形
//     TGraph *gr1 = new TGraph(_Npoints);
//     gr1->SetTitle(";Time[ns];");
//     gr1->GetXaxis()->SetTitleSize(0.15);
//     gr1->GetXaxis()->SetLabelSize(0.1);
//     gr1->GetYaxis()->SetLabelSize(0.1);
//     gr1->SetMarkerStyle(8);
//     gr1->SetMarkerSize(0.4);
//     gr1->SetMarkerColor(kRed);
//     // 模板波形
//     TGraph *gr2 = new TGraph(_Nsample);
//     gr2->SetMarkerStyle(8);
//     gr2->SetMarkerSize(0.4);
//     gr2->SetMarkerColor(kGreen);
//     TCanvas *can = new TCanvas();
//     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1);
//     pad1->SetBottomMargin(0.02); // 减小空白
//     pad1->Draw();
//     can->cd();
//     TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
//     pad2->SetTopMargin(0.02);
//     pad2->SetBottomMargin(0.30);
//     pad2->Draw();
//     TLatex *tex = new TLatex();
//     tex->SetTextSize(0.08);
//     tex->SetTextAlign(10);
//     tex->SetTextFont(42);
//     // 作图开关选项
//     bool DrawGraph = false;
//     TH2F *his_correctAT[25];
//     for (int i = 0; i < 25; i++)
//     {
//         his_correctAT[i] = new TH2F(Form("his_%d", i), ";finetime;ratio", 25, 0, 1250, 100, 0.9, 1.1);
//     }
//     TH1F *his_Aratio = new TH1F("his1", ";amplitude ratio;counts", 100, 0.95, 1.05);
//     TH2F *his_ATcor = new TH2F("his2", ";finetime[ps];amplitude ratio", 25, 0, 1250, 100, 0.95, 1.05);
//     TProfile *pro_ATcor;
//     for (int i = 0; i < nEntries; i++)
//     {
//         int progress = static_cast<float>(i + 1) / nEntries * 100;
//         if ((i + 1) % interval == 0)
//         {
//             cout << "Progress: " << progress + 1 << "%\r" << endl;
//             std::cout.flush();
//         }
//         tree->GetEntry(i);
//         mTriggerID = TriggerID;
//         mTimeCode = TimeCode;
//         for (int j = 0; j < 256; j++)
//             mTime[j] = 12.5 * j;
//         double LGAmp[_Npoints];
//         for (int j = 0; j < 25; j++)
//         {
//             mHit[j]->clear();
//             mHit[j]->Set(Hit[j]->CrystalID, Hit[j]->Temperature1, Hit[j]->Temperature2, Hit[j]->LAmplitude, Hit[j]->HAmplitude, Hit[j]->LNoise, Hit[j]->HNoise, Hit[j]->LowGainPedestal, Hit[j]->HighGainPedestal, Hit[j]->LowGainPeak, Hit[j]->HighGainPeak);
//             // 只用LG通道的波形做拟合,扣除基线
//             std::copy(Hit[j]->LAmplitude, Hit[j]->LAmplitude + _Npoints, LGAmp);
//             for (int m = 0; m < _Npoints; m++)
//                 LGAmp[m] -= LGped[j];
//             for (int k = 0; k < 150; k++)
//             {
//                 std::copy(LGAmp + k, LGAmp + k + 20, wavetmp);
//                 // 采样起始点是粗时间，拟合得到幅度和细时间(ns)
//                 amptmp = std::inner_product(ampVec, ampVec + 20, wavetmp, 0.);
//                 // 如果幅度小于20或者细时间超过一个采样点间隔，认为不是一个本底/信号
//                 if (amptmp < 20)
//                     continue;
//                 timetmp = std::inner_product(amptimeVec, amptimeVec + 20, wavetmp, 0.);
//                 timetmp /= amptmp;
//                 // 这里要判断timetmp的取值范围！！
//                 // if (fabs(timetmp) > 6.25)
//                 if (fabs(timetmp) > 20)
//                     continue;
//                 timestamp = 65 * 12.5;
//                 coarsetime = k * 12.5;
//                 // finetime = -timetmp * 100 + 625;
//                 finetime = -timetmp * 100;
//                 amplitude = amptmp;
//                 mHit[j]->AddHit(timestamp, coarsetime, finetime, amplitude);
//                 // if (amplitude > 200 && (coarsetime - timestamp) > 300 && (coarsetime - timestamp) < 500)
//                 // {
//                 //     his_correctAT[j]->Fill(finetime, amplitude / (Hit[j]->LowGainPeak - LGped[j]));
//                 // }
//                 // if (j == 13 && amplitude > 200 && (coarsetime - timestamp) > 300 && (coarsetime - timestamp) < 500)
//                 // {
//                 //     his_ATcor->Fill(finetime, amplitude / (Hit[j]->LowGainPeak - LGped[j]));
//                 //     his_Aratio->Fill(amplitude / (Hit[j]->LowGainPeak - Hit[j]->LowGainPedestal));
//                 // }
//                 if (DrawGraph)
//                 {
//                     pad1->cd();
//                     for (int m = 0; m < _Npoints; m++)
//                         gr->SetPoint(m, 12.5 * m, LGAmp[m]);
//                     gr->Draw("ap");
//                     gr->GetXaxis()->SetLabelSize(0);
//                 }
//                 // 在波形上扣除这个本底/信号
//                 int tempN = static_cast<int>(finetime);
//                 // std::cout << "tempN=" << tempN << std::endl;
//                 for (int m = 0; m < _Nsample; m++)
//                 {
//                     // LGAmp[k + m] -= tempVec[tempN % 1250][m] * amptmp;
//                     if (DrawGraph)
//                     {
//                         gr2->SetPoint(m, (k + m) * 12.5, tempVec[tempN][m] * amptmp);
//                         gr2->Draw("p same");
//                         tex->DrawLatexNDC(0.6, 0.7, Form("peak value=%.0lf", Hit[j]->LowGainPeak - LGped[j]));
//                         tex->DrawLatexNDC(0.6, 0.6, Form("fit result=%.0lf", amplitude));
//                         // tex->DrawLatexNDC(0.6, 0.4, Form("fine time=%.2lf ns", timetmp));
//                     }
//                 }
//                 if (DrawGraph)
//                 {
//                     pad2->cd();
//                     for (int m = 0; m < _Npoints; m++)
//                         gr1->SetPoint(m, 12.5 * m, LGAmp[m]);
//                     gr1->Draw("ap");
//                     if (amplitude > 200)
//                         can->SaveAs(Form("subtract_Event%i_Hit%i_Point%i.png", i, j, k));
//                 }
//             }
//         }
//         mTree->Fill();
//     }
//     double factor[25][25];
//     for (int i = 0; i < 25; i++)
//     {
//         TProfile *prof = his_correctAT[i]->ProfileX();
//         for (int j = 0; j < 25; j++)
//         {
//             factor[i][j] = prof->GetBinContent(j + 1);
//         }
//     }
//     // std::cout << "calibration factor:" << std::endl;
//     // for (int i = 0; i < 25; i++)
//     // {
//     //     std::cout << '{';
//     //     for (int j = 0; j < 25; j++)
//     //     {
//     //         if (factor[i][j] == 0)
//     //             factor[i][j] = 1;
//     //         if (j != 24)
//     //             std::cout << 1 / factor[i][j] << ',';
//     //         else
//     //             std::cout << 1 / factor[i][j] << "}";
//     //     }
//     //     std::cout << ',';
//     // }
//     // his_ATcor->Draw("colz");
//     // gPad->SaveAs("hisATcor.png");
//     // pro_ATcor = his_ATcor->ProfileX("fpx");
//     // pro_ATcor->GetYaxis()->SetTitle("amplitude ratio");
//     // pro_ATcor->SetErrorOption("i");
//     // pro_ATcor->Draw();
//     // gPad->SaveAs("proATcor.png");
//     // his_Aratio->Draw();
//     // gPad->SaveAs("hisAratio.png");
//     fout->cd();
//     mTree->Write();
//     fout->Close();

//     return 0;
// }

// 两个提取向量，间距0.5个采样点
int main(int argc, char const *argv[])
{
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    vector<string> datafiles;

    if (isTextFile(argv[1]))
    {
        ifstream filestream(argv[1]);
        if (filestream.is_open())
        {
            string filename;
            while (getline(filestream, filename))
                datafiles.push_back(filename);
        }
        else
            cout << "error in reading filelist" << endl;
    }
    else
        datafiles.push_back(argv[1]);

    TChain *tree = new TChain("decode_data");
    for (size_t i = 0; i < datafiles.size(); i++)
        tree->AddFile(datafiles.at(i).data());
    int TriggerID;
    Long64_t TimeCode;
    tree->SetBranchAddress("TriggerID", &TriggerID);
    tree->SetBranchAddress("TimeCode", &TimeCode);
    vector<decode_data_col *> Hit(25);
    string Name[25];
    for (int i = 0; i < 25; i++)
    {
        Name[i] = to_string(i / 5 + 1) + "_" + to_string(i % 5 + 1);
    }
    for (int i = 0; i < 25; i++)
        Hit[i] = new decode_data_col();
    for (int i = 0; i < 25; i++)
    {
        string name = "Hit_" + Name[i];
        tree->SetBranchAddress(name.data(), Hit[i]);
    }

    TString outputfile;
    if (argc >= 3)
        outputfile = argv[2];
    else
    {
        outputfile = "ECALDigiOnline.root";
        cout << "Auto save file as ECALDigiOnline.root..." << endl;
    }

    // 输出文件
    TFile *fout = new TFile(outputfile.Data(), "recreate");
    TTree *mTree = new TTree("decode_data", "decode_data");
    DataModel2025 *mHit[25];
    for (int i = 0; i < 25; i++)
        mHit[i] = new DataModel2025();
    Long64_t mEventID;
    int mTriggerID;
    Long64_t mTimeCode;
    float mTime[_Npoints];
    float mTemperature[10];
    mTree->Branch("TriggerID", &mTriggerID, "TriggerID/I");
    mTree->Branch("TimeCode", &mTimeCode, "TimeCode/L");
    mTree->Branch("Time", &mTime, Form("Time[%d]/F", _Npoints));

    for (int i = 0; i < 25; i++)
    {
        Name[i] = std::to_string(i / 5 + 1) + "_" + std::to_string(i % 5 + 1);
        std::string name = "Hit_" + Name[i];
        std::string leaf_list = "";

        leaf_list += "CrystalID/L";
        leaf_list += ":Temperature1/D";
        leaf_list += ":Temperature2/D";
        leaf_list += Form(":LAmplitude[%d]/D", _Npoints);
        leaf_list += Form(":HAmplitude[%d]/D", _Npoints);
        leaf_list += Form(":LNoise[%d]/D", _Nnoise);
        leaf_list += Form(":HNoise[%d]/D", _Nnoise);
        leaf_list += ":LowGainPedestal/D";
        leaf_list += ":HighGainPedestal/D";
        leaf_list += ":LowGainPeak/D";
        leaf_list += ":HighGainPeak/D";

        mTree->Branch(name.data(), &mHit[i]->CrystalID, leaf_list.data());
        mTree->Branch(Form("%s_TimeStamp", name.data()), &mHit[i]->TimeStamp, "TimeStamp/D");
        mTree->Branch(Form("%s_CoarseTime", name.data()), &mHit[i]->CoarseTime);
        mTree->Branch(Form("%s_FineTime", name.data()), &mHit[i]->FineTime);
        mTree->Branch(Form("%s_Amplitude", name.data()), &mHit[i]->Amplitude);
    }
    double timestamp, coarsetime, finetime, amplitude;
    int nEntries = tree->GetEntries();
    int interval = nEntries / 20;
    double sampletemp[1250][_Nsample];
    TemplateSample(sampletemp);
    const double LGped[25] = {2297.49, 2305.03, 2338.2, 2377.62, 2354.75, 2365.34, 2347.87, 2365.32, 2349.4, 2326.33, 2347.65, 2354.71, 2414.79, 2379.85, 2381.56, 2308.86, 2323.2, 2365.7, 2365.65, 2318.67, 2337, 2363.79, 2352.87, 2328.34, 2329.16};
    // 分别是模板向量，幅度向量和时幅向量
    const double (&tempVec)[1250][_Nsample] = sampletemp;
    const double ampVec1[10] = {-1.8562517e-05, -0.094403470, -0.22365875, -0.27526923, -0.22859872, -0.10309496, 0.069436969, 0.25839261, 0.43991157, 0.59815506};
    const double amptimeVec1[10] = {0.0015056472, 7.9746465, 20.061841, 27.231616, 27.339226, 21.369777, 11.392180, -0.46948477, -12.482087, -23.433152};
    const double ampVec2[10] = {-0.023987174, -0.13178224, -0.21082495, -0.21469844, -0.14624425, -0.027363541, 0.11655549, 0.26374905, 0.39851646, 0.51128031};
    const double amptimeVec2[10] = {2.2196555, 12.852344, 22.274154, 26.019414, 23.796107, 17.042891, 7.6387993, -2.7110364, -12.711270, -21.512936};
    double wavetmp[10], timetmp, amptmp, timetmp1, amptmp1, timetmp2, amptmp2;
    bool useTemp1, useTemp2;
    // 原始波形
    TGraph *gr = new TGraph(_Npoints);
    gr->SetTitle(";;ADC value");
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetTitleOffset(0.5);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(0.4);
    // 剩余波形
    TGraph *gr1 = new TGraph(_Npoints);
    gr1->SetTitle(";Time[ns];");
    gr1->GetXaxis()->SetTitleSize(0.15);
    gr1->GetXaxis()->SetLabelSize(0.1);
    gr1->GetYaxis()->SetLabelSize(0.1);
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerSize(0.4);
    gr1->SetMarkerColor(kRed);
    // 模板波形
    TGraph *gr2 = new TGraph(_Nsample);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerSize(0.4);
    gr2->SetMarkerColor(kGreen);
    TCanvas *can = new TCanvas();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1);
    pad1->SetBottomMargin(0.02); // 减小空白
    pad1->Draw();
    can->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.30);
    pad2->Draw();
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.08);
    tex->SetTextAlign(10);
    tex->SetTextFont(42);
    // 作图开关选项
    bool DrawGraph = false;
    TH1F *his_Aratio = new TH1F("his1", ";amplitude ratio;counts", 100, 0.95, 1.05);
    TH2F *his_ATcor = new TH2F("his2", ";finetime[ns];amplitude ratio", 25, 0, 1250, 100, 0.95, 1.05);
    TProfile *pro_ATcor;
    for (int i = 0; i < nEntries; i++)
    {
        int progress = static_cast<float>(i + 1) / nEntries * 100;
        if ((i + 1) % interval == 0)
        {
            cout << "Progress: " << progress + 1 << "%\r" << endl;
            std::cout.flush();
        }
        tree->GetEntry(i);
        mTriggerID = TriggerID;
        mTimeCode = TimeCode;
        for (int j = 0; j < 256; j++)
            mTime[j] = 12.5 * j;
        double LGAmp[_Npoints];
        for (int j = 0; j < 25; j++)
        {
            mHit[j]->clear();
            mHit[j]->Set(Hit[j]->CrystalID, Hit[j]->Temperature1, Hit[j]->Temperature2, Hit[j]->LAmplitude, Hit[j]->HAmplitude, Hit[j]->LNoise, Hit[j]->HNoise, Hit[j]->LowGainPedestal, Hit[j]->HighGainPedestal, Hit[j]->LowGainPeak, Hit[j]->HighGainPeak);
            // 只用LG通道的波形做拟合,扣除基线
            std::copy(Hit[j]->LAmplitude, Hit[j]->LAmplitude + _Npoints, LGAmp);
            for (int m = 0; m < _Npoints; m++)
                LGAmp[m] -= LGped[j];
            for (int k = 0; k < 150; k++)
            {
                std::copy(LGAmp + k, LGAmp + k + 10, wavetmp);
                // 采样起始点是粗时间，拟合得到幅度和细时间(ns)
                amptmp1 = std::inner_product(ampVec1, ampVec1 + 10, wavetmp, 0.);
                timetmp1 = std::inner_product(amptimeVec1, amptimeVec1 + 10, wavetmp, 0.);
                timetmp1 /= amptmp1;
                amptmp2 = std::inner_product(ampVec2, ampVec2 + 10, wavetmp, 0.);
                timetmp2 = std::inner_product(amptimeVec2, amptimeVec2 + 10, wavetmp, 0.);
                timetmp2 /= amptmp2;
                // 判断使用提取向量1或者提取向量2
                // useTemp1 = (amptmp1 > 20 && fabs(timetmp1) < 3.125);
                // useTemp2 = (amptmp2 > 20 && fabs(timetmp2) < 3.125);
                useTemp1 = (amptmp1 > 20 && fabs(timetmp1) < 4);
                useTemp2 = (amptmp2 > 20 && fabs(timetmp2) < 4);

                if (useTemp1 && !useTemp2)
                {
                    timetmp = timetmp1;
                    amptmp = amptmp1;
                }
                else if (!useTemp1 && useTemp2)
                {
                    timetmp = timetmp2 + 6.25;
                    amptmp = amptmp2;
                }
                else if (useTemp1 && useTemp2)
                {
                    if (fabs(timetmp1) < fabs(timetmp2))
                    {
                        timetmp = timetmp1;
                        amptmp = amptmp1;
                    }
                    else
                    {
                        timetmp = timetmp2 + 6.25;
                        amptmp = amptmp2;
                    }
                }
                else
                    continue;
                timestamp = 65 * 12.5;
                coarsetime = k * 12.5;
                finetime = timetmp * 100 + 625;
                amplitude = amptmp;
                mHit[j]->AddHit(timestamp, coarsetime, finetime, amplitude);
                if (j == 12 && amplitude > 200 && (coarsetime - timestamp) > 300 && (coarsetime - timestamp) < 500)
                {
                    his_ATcor->Fill(finetime, amplitude / (Hit[j]->LowGainPeak - LGped[j]));
                    his_Aratio->Fill(amplitude / (Hit[j]->LowGainPeak - Hit[j]->LowGainPedestal));
                }
                if (DrawGraph)
                {
                    pad1->cd();
                    for (int m = 0; m < _Npoints; m++)
                        gr->SetPoint(m, 12.5 * m, LGAmp[m]);
                    gr->Draw("ap");
                    gr->GetXaxis()->SetLabelSize(0);
                }
                // 在波形上扣除这个本底/信号
                int tempN = static_cast<int>(finetime);
                // std::cout << "tempN=" << tempN << std::endl;
                for (int m = 0; m < _Nsample; m++)
                {
                    // LGAmp[k + m] -= tempVec[tempN % 1250][m] * amptmp;
                    if (DrawGraph)
                    {
                        gr2->SetPoint(m, (k + m) * 12.5, tempVec[tempN][m] * amptmp);
                        gr2->Draw("p same");
                        tex->DrawLatexNDC(0.6, 0.7, Form("peak value=%.0lf", Hit[j]->LowGainPeak - LGped[j]));
                        tex->DrawLatexNDC(0.6, 0.6, Form("fit result=%.0lf", amplitude));
                        // tex->DrawLatexNDC(0.6, 0.4, Form("fine time=%.2lf ns", timetmp));
                    }
                }
                if (DrawGraph)
                {
                    pad2->cd();
                    for (int m = 0; m < _Npoints; m++)
                        gr1->SetPoint(m, 12.5 * m, LGAmp[m]);
                    gr1->Draw("ap");
                    if (amplitude > 200)
                        can->SaveAs(Form("subtract_Event%i_Hit%i_Point%i.png", i, j, k));
                }
            }
        }
        mTree->Fill();
    }
    // his_ATcor->Draw("colz");
    // gPad->SaveAs("hisATcor.png");
    // pro_ATcor = his_ATcor->ProfileX("fpx");
    // pro_ATcor->GetYaxis()->SetTitle("amplitude ratio");
    // pro_ATcor->SetErrorOption("i");
    // pro_ATcor->Draw();
    // gPad->SaveAs("proATcor.png");
    // his_Aratio->Draw();
    // gPad->SaveAs("hisAratio.png");
    fout->cd();
    mTree->Write();
    fout->Close();

    return 0;
}

// 拼接高低增益通道波形
// int main(int argc, char const *argv[])
// {
//     gStyle->SetPadGridX(true);
//     gStyle->SetPadGridY(true);
//     vector<string> datafiles;

//     if (isTextFile(argv[1]))
//     {
//         ifstream filestream(argv[1]);
//         if (filestream.is_open())
//         {
//             string filename;
//             while (getline(filestream, filename))
//                 datafiles.push_back(filename);
//         }
//         else
//             cout << "error in reading filelist" << endl;
//     }
//     else
//         datafiles.push_back(argv[1]);

//     TChain *tree = new TChain("decode_data");
//     for (size_t i = 0; i < datafiles.size(); i++)
//         tree->AddFile(datafiles.at(i).data());
//     int TriggerID;
//     Long64_t TimeCode;
//     tree->SetBranchAddress("TriggerID", &TriggerID);
//     tree->SetBranchAddress("TimeCode", &TimeCode);
//     vector<decode_data_col *> Hit(25);
//     string Name[25];
//     for (int i = 0; i < 25; i++)
//     {
//         Name[i] = to_string(i / 5 + 1) + "_" + to_string(i % 5 + 1);
//     }
//     for (int i = 0; i < 25; i++)
//         Hit[i] = new decode_data_col();
//     for (int i = 0; i < 25; i++)
//     {
//         string name = "Hit_" + Name[i];
//         tree->SetBranchAddress(name.data(), Hit[i]);
//     }

//     TString outputfile;
//     if (argc >= 3)
//         outputfile = argv[2];
//     else
//     {
//         outputfile = "ECALDigiOnline.root";
//         cout << "Auto save file as ECALDigiOnline.root..." << endl;
//     }

//     // 输出文件
//     TFile *fout = new TFile(outputfile.Data(), "recreate");
//     TTree *mTree = new TTree("decode_data", "decode_data");
//     DataModel2025 *mHit[25];
//     for (int i = 0; i < 25; i++)
//         mHit[i] = new DataModel2025();
//     Long64_t mEventID;
//     int mTriggerID;
//     Long64_t mTimeCode;
//     float mTime[_Npoints];
//     float mTemperature[10];
//     mTree->Branch("TriggerID", &mTriggerID, "TriggerID/I");
//     mTree->Branch("TimeCode", &mTimeCode, "TimeCode/L");
//     mTree->Branch("Time", &mTime, Form("Time[%d]/F", _Npoints));

//     for (int i = 0; i < 25; i++)
//     {
//         Name[i] = std::to_string(i / 5 + 1) + "_" + std::to_string(i % 5 + 1);
//         std::string name = "Hit_" + Name[i];
//         std::string leaf_list = "";

//         leaf_list += "CrystalID/L";
//         leaf_list += ":Temperature1/D";
//         leaf_list += ":Temperature2/D";
//         leaf_list += Form(":LAmplitude[%d]/D", _Npoints);
//         leaf_list += Form(":HAmplitude[%d]/D", _Npoints);
//         leaf_list += Form(":LNoise[%d]/D", _Nnoise);
//         leaf_list += Form(":HNoise[%d]/D", _Nnoise);
//         leaf_list += ":LowGainPedestal/D";
//         leaf_list += ":HighGainPedestal/D";
//         leaf_list += ":LowGainPeak/D";
//         leaf_list += ":HighGainPeak/D";

//         mTree->Branch(name.data(), &mHit[i]->CrystalID, leaf_list.data());
//         mTree->Branch(Form("%s_TimeStamp", name.data()), &mHit[i]->TimeStamp, "TimeStamp/D");
//         mTree->Branch(Form("%s_CoarseTime", name.data()), &mHit[i]->CoarseTime);
//         mTree->Branch(Form("%s_FineTime", name.data()), &mHit[i]->FineTime);
//         mTree->Branch(Form("%s_Amplitude", name.data()), &mHit[i]->Amplitude);
//     }
//     double timestamp, coarsetime, finetime, amplitude;
//     int nEntries = tree->GetEntries();
//     int interval = nEntries / 20;
//     double sampletemp[1250][_Nsample];
//     TemplateSample(sampletemp);
//     const double LGped[25] = {2297.49, 2305.03, 2338.2, 2377.62, 2354.75, 2365.34, 2347.87, 2365.32, 2349.4, 2326.33, 2347.65, 2354.71, 2414.79, 2379.85, 2381.56, 2308.86, 2323.2, 2365.7, 2365.65, 2318.67, 2337, 2363.79, 2352.87, 2328.34, 2329.16};
//     const double HGped[25] = {2320.09, 2302.05, 2320.62, 2395.44, 2397.38, 2330.23, 2310.94, 2317.42, 2298.1, 2326.95, 2359.09, 2434.24, 2356.82, 2367.07, 2303.62, 2342.62, 2358.67, 2336.86, 2385.94, 2315.43, 2326.41, 2329.69, 2301.34, 2355.02, 2299.13};
//     const double GainRatio[25] = {10.0069, 10.0467, 10.1415, 10.0949, 10.0537, 10.149, 10.1223, 10.1895, 10.409, 10.1114, 9.99876, 10.2592, 10.0951, 10.0875, 10.1467, 9.98977, 10.002, 10.169, 10.0875, 10.0776, 10.0623, 10.1243, 10.1987, 10.1202, 10.0715};
//     // 分别是模板向量，幅度向量和时幅向量
//     const double (&tempVec)[1250][_Nsample] = sampletemp;
//     const double ampVecLG[20] = {-6.38357040248042e-05, -0.00330697007822949, -0.00398696499554231, 0.00179119025323171, 0.0134049745184865, 0.0285970474714096, 0.0449272265503205, 0.0604065482502062, 0.0736885867430082, 0.0840412605666282, 0.0912268879037425, 0.0953595590576070, 0.0967738233148716, 0.0959185309201458, 0.0932787330236419, 0.0893231471006163, 0.0844725072062953, 0.0790836856326403, 0.0734449149428907, 0.0677782432443050};
//     const double amptimeVecLG[20] = {-0.0472543391531407, -3.43817924184712, -8.37953139807132, -12.0715530570965, -13.7001056339471, -13.3994007185142, -11.6731104324930, -9.09225696548941, -6.15481131398440, -3.23549247413817, -0.582418053622125, 1.66481780411593, 3.44986840988616, 4.77456618906330, 5.67744119483497, 6.21693646975686, 6.45923915076050, 6.47009666969433, 6.30980287939698, 6.03054868872617};
//     const double ampVecHG[20] = {-0.00172075776600942, -0.00864246126568028, -0.0121287113716010, -0.00879953116721358, 0.000935957651454447, 0.0151517712372137, 0.0316124912449845, 0.0483499572188746, 0.0638817779548702, 0.0772390195651733, 0.0879028492750695, 0.0957072280162231, 0.100738352221455, 0.103245678705214, 0.103570192419650, 0.102090572799003, 0.0991854413438790, 0.0952088814952126, 0.0904762526322076};
//     const double amptimeVecHG[20] = {-0.870388648985807, -5.29163393443170, -10.0652370036382, -13.3735967418439, -14.7431955352232, -14.3531597645547, -12.6414123909776, -10.0936295059682, -7.14128317740175, -4.12221677481193, -1.27496995523588, 1.25045432601793, 3.37539600635657, 5.07439039640397, 6.35878408794224, 7.26321816595203, 7.83518682682308, 8.12742328760399, 8.19266932763705};
//     double wavetmp[20], timetmp, amptmp;
//     // 原始波形
//     TGraph *gr = new TGraph(_Npoints);
//     gr->SetTitle(";;ADC value");
//     gr->GetYaxis()->SetTitleSize(0.06);
//     gr->GetYaxis()->SetTitleOffset(0.5);
//     gr->GetYaxis()->SetLabelSize(0.04);
//     gr->SetMarkerStyle(8);
//     gr->SetMarkerSize(0.4);
//     // 剩余波形
//     TGraph *gr1 = new TGraph(_Npoints);
//     gr1->SetTitle(";Time[ns];");
//     gr1->GetXaxis()->SetTitleSize(0.15);
//     gr1->GetXaxis()->SetLabelSize(0.1);
//     gr1->GetYaxis()->SetLabelSize(0.1);
//     gr1->SetMarkerStyle(8);
//     gr1->SetMarkerSize(0.4);
//     gr1->SetMarkerColor(kRed);
//     // 模板波形
//     TGraph *gr2 = new TGraph(_Nsample);
//     gr2->SetMarkerStyle(8);
//     gr2->SetMarkerSize(0.4);
//     gr2->SetMarkerColor(kGreen);
//     TCanvas *can = new TCanvas();
//     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1);
//     pad1->SetBottomMargin(0.02); // 减小空白
//     pad1->Draw();
//     can->cd();
//     TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);
//     pad2->SetTopMargin(0.02);
//     pad2->SetBottomMargin(0.30);
//     pad2->Draw();
//     TLatex *tex = new TLatex();
//     tex->SetTextSize(0.08);
//     tex->SetTextAlign(10);
//     tex->SetTextFont(42);
//     // 作图开关选项
//     bool DrawGraph = false;
//     TH1F *his_Aratio = new TH1F("his1", ";amplitude ratio;counts", 100, 0.95, 1.05);
//     TH2F *his_ATcor = new TH2F("his2", ";finetime[ns];amplitude ratio", 25, 0, 1250, 100, 0.95, 1.05);
//     TProfile *pro_ATcor;
//     for (int i = 0; i < nEntries; i++)
//     {
//         int progress = static_cast<float>(i + 1) / nEntries * 100;
//         if ((i + 1) % interval == 0)
//         {
//             cout << "Progress: " << progress + 1 << "%\r" << endl;
//             std::cout.flush();
//         }
//         tree->GetEntry(i);
//         mTriggerID = TriggerID;
//         mTimeCode = TimeCode;
//         for (int j = 0; j < 256; j++)
//             mTime[j] = 12.5 * j;
//         double LGAmp[_Npoints];
//         double HGAmp[_Npoints];
//         double FilterAmp[_Npoints];
//         for (int j = 0; j < 25; j++)
//         {
//             if (j != 12)
//             {
//                 mHit[j]->clear();
//                 mHit[j]->Set(Hit[j]->CrystalID, Hit[j]->Temperature1, Hit[j]->Temperature2, Hit[j]->LAmplitude, Hit[j]->HAmplitude, Hit[j]->LNoise, Hit[j]->HNoise, Hit[j]->LowGainPedestal, Hit[j]->HighGainPedestal, Hit[j]->LowGainPeak, Hit[j]->HighGainPeak);
//                 // 只用LG通道的波形做拟合,扣除基线
//                 std::copy(Hit[j]->LAmplitude, Hit[j]->LAmplitude + _Npoints, LGAmp);
//                 for (int m = 0; m < _Npoints; m++)
//                     LGAmp[m] -= LGped[j];
//                 for (int k = 0; k < 150; k++)
//                 {
//                     std::copy(LGAmp + k, LGAmp + k + 20, wavetmp);
//                     // 采样起始点是粗时间，拟合得到幅度和细时间(ns)
//                     amptmp = std::inner_product(ampVecLG, ampVecLG + 20, wavetmp, 0.);
//                     // 如果幅度小于20或者细时间超过一个采样点间隔，认为不是一个本底/信号
//                     if (amptmp < 20)
//                         continue;
//                     timetmp = std::inner_product(amptimeVecLG, amptimeVecLG + 20, wavetmp, 0.);
//                     timetmp /= amptmp;
//                     // 这里要判断timetmp的取值范围！！
//                     if (fabs(timetmp) > 6.25)
//                         continue;
//                     timestamp = 65 * 12.5;
//                     coarsetime = k * 12.5;
//                     finetime = -timetmp * 100 + 625;
//                     amplitude = amptmp;
//                     mHit[j]->AddHit(timestamp, coarsetime, finetime, amplitude);
//                     if (j == 12 && amplitude > 200 && (coarsetime - timestamp) > 300 && (coarsetime - timestamp) < 500)
//                     {
//                         his_ATcor->Fill(finetime, amplitude / (Hit[j]->LowGainPeak - LGped[j]));
//                         his_Aratio->Fill(amplitude / (Hit[j]->LowGainPeak - Hit[j]->LowGainPedestal));
//                     }
//                     if (DrawGraph)
//                     {
//                         pad1->cd();
//                         for (int m = 0; m < _Npoints; m++)
//                             gr->SetPoint(m, 12.5 * m, LGAmp[m]);
//                         gr->Draw("ap");
//                         gr->GetXaxis()->SetLabelSize(0);
//                     }
//                     // 在波形上扣除这个本底/信号
//                     int tempN = static_cast<int>(finetime);
//                     // std::cout << "tempN=" << tempN << std::endl;
//                     for (int m = 0; m < _Nsample; m++)
//                     {
//                         LGAmp[k + m] -= tempVec[tempN][m] * amptmp;
//                         if (DrawGraph)
//                         {
//                             gr2->SetPoint(m, (k + m) * 12.5, tempVec[tempN][m] * amptmp);
//                             gr2->Draw("p same");
//                             tex->DrawLatexNDC(0.6, 0.7, Form("peak value=%.0lf", Hit[j]->LowGainPeak - LGped[j]));
//                             tex->DrawLatexNDC(0.6, 0.6, Form("fit result=%.0lf", amplitude));
//                             // tex->DrawLatexNDC(0.6, 0.4, Form("fine time=%.2lf ns", timetmp));
//                         }
//                     }
//                     if (DrawGraph)
//                     {
//                         pad2->cd();
//                         for (int m = 0; m < _Npoints; m++)
//                             gr1->SetPoint(m, 12.5 * m, LGAmp[m]);
//                         gr1->Draw("ap");
//                         if (amplitude > 200)
//                             can->SaveAs(Form("subtract_Event%i_Hit%i_Point%i.png", i, j, k));
//                     }
//                 }
//             }
//             else
//             {
//                 mHit[j]->clear();
//                 mHit[j]->Set(Hit[j]->CrystalID, Hit[j]->Temperature1, Hit[j]->Temperature2, Hit[j]->LAmplitude, Hit[j]->HAmplitude, Hit[j]->LNoise, Hit[j]->HNoise, Hit[j]->LowGainPedestal, Hit[j]->HighGainPedestal, Hit[j]->LowGainPeak, Hit[j]->HighGainPeak);
//                 // HG没饱和，不用拼接
//                 if (Hit[j]->HighGainPeak < 16000)
//                 {
//                     std::copy(Hit[j]->HAmplitude, Hit[j]->HAmplitude + _Npoints, HGAmp);
//                     for (int m = 0; m < _Npoints; m++)
//                         HGAmp[m] -= HGped[j];
//                 }
//                 // HG饱和需要拼接
//                 else
//                 {
//                     CombineWave(Hit[j]->LAmplitude, Hit[j]->HAmplitude, HGAmp);
//                 }
//                 for (int k = 0; k < 150; k++)
//                 {
//                     std::copy(HGAmp + k, HGAmp + k + 20, wavetmp);
//                     // 采样起始点是粗时间，拟合得到幅度和细时间(ns)
//                     amptmp = std::inner_product(ampVecHG, ampVecHG + 20, wavetmp, 0.);
//                     // 如果幅度小于100或者细时间超过一个采样点间隔，认为不是一个本底/信号
//                     if (amptmp < 100)
//                         continue;
//                     timetmp = std::inner_product(amptimeVecHG, amptimeVecHG + 20, wavetmp, 0.);
//                     timetmp /= amptmp;
//                     // 这里要判断timetmp的取值范围！！
//                     if (fabs(timetmp) > 6.25)
//                         continue;
//                     timestamp = 65 * 12.5;
//                     coarsetime = k * 12.5;
//                     finetime = -timetmp * 100 + 625;
//                     amplitude = amptmp / GainRatio[j];
//                     mHit[j]->AddHit(timestamp, coarsetime, finetime, amplitude);
//                     if (j == 12 && amplitude > 200 && (coarsetime - timestamp) > 300 && (coarsetime - timestamp) < 500)
//                     {
//                         his_ATcor->Fill(finetime, amplitude / (Hit[j]->LowGainPeak - LGped[j]));
//                         his_Aratio->Fill(amplitude / (Hit[j]->LowGainPeak - Hit[j]->LowGainPedestal));
//                     }
//                     if (DrawGraph)
//                     {
//                         pad1->cd();
//                         for (int m = 0; m < _Npoints; m++)
//                             gr->SetPoint(m, 12.5 * m, HGAmp[m]);
//                         gr->Draw("ap");
//                         gr->GetXaxis()->SetLabelSize(0);
//                     }
//                     // 在波形上扣除这个本底/信号
//                     int tempN = static_cast<int>(finetime);
//                     // std::cout << "tempN=" << tempN << std::endl;
//                     for (int m = 0; m < _Nsample; m++)
//                     {
//                         HGAmp[k + m] -= tempVec[tempN][m] * amptmp;
//                         if (DrawGraph)
//                         {
//                             gr2->SetPoint(m, (k + m) * 12.5, tempVec[tempN][m] * amptmp);
//                             gr2->Draw("p same");
//                             tex->DrawLatexNDC(0.6, 0.7, Form("peak value=%.0lf", Hit[j]->LowGainPeak - LGped[j]));
//                             tex->DrawLatexNDC(0.6, 0.6, Form("fit result=%.0lf", amplitude));
//                             // tex->DrawLatexNDC(0.6, 0.4, Form("fine time=%.2lf ns", timetmp));
//                         }
//                     }
//                     if (DrawGraph)
//                     {
//                         pad2->cd();
//                         for (int m = 0; m < _Npoints; m++)
//                             gr1->SetPoint(m, 12.5 * m, HGAmp[m]);
//                         gr1->Draw("ap");
//                         if (amplitude > 200)
//                             can->SaveAs(Form("subtract_Event%i_Hit%i_Point%i.png", i, j, k));
//                     }
//                 }
//             }
//         }
//         mTree->Fill();
//     }
//     fout->cd();
//     mTree->Write();
//     fout->Close();

//     return 0;
// }