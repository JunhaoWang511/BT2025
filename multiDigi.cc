#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <TString.h>
#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>
#include <TH1.h>
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
    // 采样起始点从69.5个点（856.25 ns）到70.5个点（868.75 ns），采样1250份（间隔10 ps），采样长度为40个点
    double sampletime;
    for (int i = 0; i < 1250; i++)
        for (int j = 0; j < _Nsample; j++)
        {
            sampletime = 856.25 + 0.01 * i + 12.5 * j;
            sample[i][j] = fun->Eval(sampletime);
        }
}

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
    const double ampVec[20] = {-6.38357040248042e-05, -0.00330697007822949, -0.00398696499554231, 0.00179119025323171, 0.0134049745184865, 0.0285970474714096, 0.0449272265503205, 0.0604065482502062, 0.0736885867430082, 0.0840412605666282, 0.0912268879037425, 0.0953595590576070, 0.0967738233148716, 0.0959185309201458, 0.0932787330236419, 0.0893231471006163, 0.0844725072062953, 0.0790836856326403, 0.0734449149428907, 0.0677782432443050};
    const double amptimeVec[20] = {-0.0472543391531407, -3.43817924184712, -8.37953139807132, -12.0715530570965, -13.7001056339471, -13.3994007185142, -11.6731104324930, -9.09225696548941, -6.15481131398440, -3.23549247413817, -0.582418053622125, 1.66481780411593, 3.44986840988616, 4.77456618906330, 5.67744119483497, 6.21693646975686, 6.45923915076050, 6.47009666969433, 6.30980287939698, 6.03054868872617};
    double wavetmp[20], timetmp, amptmp;
    TGraph *gr = new TGraph(_Npoints);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(0.4);
    TGraph *gr1 = new TGraph(_Npoints);
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerSize(0.4);
    gr1->SetMarkerColor(kRed);
    TGraph *gr2 = new TGraph(_Nsample);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerSize(0.4);
    gr2->SetMarkerColor(kGreen);
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.06);
    tex->SetTextAlign(10);
    bool DrawGraph = false;
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
            for (int k = 0; k < 160; k++)
            {
                std::copy(LGAmp + k, LGAmp + k + 20, wavetmp);
                // 采样起始点是粗时间，拟合得到幅度和细时间(ns)
                amptmp = std::inner_product(ampVec, ampVec + 20, wavetmp, 0);
                // 修正拟合幅度
                amptmp += 10;
                // 如果幅度小于20或者细时间超过一个采样点间隔，认为不是一个本底/信号
                if (amptmp < 20)
                    continue;
                timetmp = std::inner_product(amptimeVec, amptimeVec + 20, wavetmp, 0);
                timetmp /= amptmp;
                // 这里要判断timetmp的取值范围！！
                if (fabs(timetmp) > 6.25)
                    continue;
                // std::cout << " amplitude=" << amptmp << std::endl;
                // std::cout << "fine time=" << timetmp << std::endl;
                timestamp = 64 * 12.5;
                coarsetime = k * 12.5;
                finetime = -timetmp * 100 + 625;
                amplitude = amptmp;
                mHit[j]->AddHit(timestamp, coarsetime, finetime, amplitude);
                if (DrawGraph)
                {
                    for (int m = 0; m < _Npoints; m++)
                        gr->SetPoint(m, 12.5 * m, LGAmp[m]);
                    gr->Draw("ap");
                }
                // 在波形上扣除这个本底/信号
                int tempN = static_cast<int>(finetime);
                // std::cout << "tempN=" << tempN << std::endl;
                for (int m = 0; m < _Nsample; m++)
                {
                    LGAmp[k + m] -= tempVec[tempN][m] * amptmp;
                    if (DrawGraph)
                    {
                        gr2->SetPoint(m, (k + m) * 12.5, tempVec[tempN][m] * amptmp);
                        gr2->Draw("p same");
                    }
                }
                if (DrawGraph)
                {
                    for (int m = 0; m < _Npoints; m++)
                        gr1->SetPoint(m, 12.5 * m, LGAmp[m]);
                    gr1->Draw("p same");
                    tex->DrawLatexNDC(0.5, 0.7, Form("amp peak=%.0lf", Hit[j]->LowGainPeak - LGped[j]));
                    tex->DrawLatexNDC(0.5, 0.6, Form("fit peak=%.0lf", amplitude));
                    tex->DrawLatexNDC(0.5, 0.4, Form("fine time=%.2lf", timetmp));
                    if (amplitude > 500)
                        gPad->SaveAs(Form("subtract_Event%i_Hit%i_Point%i.png", i, j, k));
                }
            }
        }
        mTree->Fill();
    }

    fout->cd();
    mTree->Write();
    fout->Close();

    return 0;
}