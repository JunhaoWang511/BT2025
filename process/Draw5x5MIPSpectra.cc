#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>

#include "DataModel2025.hh"

using namespace std;
using namespace TMath;

double langaufun(Double_t *x, Double_t *par)
{

    // Fit parameters:
    // par[0]=Width (scale) parameter of Landau density
    // par[1]=Most Probable (MP, location) parameter of Landau density
    // par[2]=Total area (integral -inf to inf, normalization constant)
    // par[3]=Width (sigma) of convoluted Gaussian function
    //
    // In the Landau distribution (represented by the CERNLIB approximation),
    // the maximum is located at x=-0.22278298 with the location parameter=0.
    // This shift is corrected within this function, so that the actual
    // maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
    Double_t mpshift = -0.22278298;      // Landau maximum location

    // Control constants
    Double_t np = 100.0; // number of convolution steps
    Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow, xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp - xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for (i = 1.0; i <= np / 2; i++)
    {
        xx = xlow + (i - .5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);

        xx = xupp - (i - .5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

void Draw5x5MIPSpectra(TString rootname)
{
    gStyle->SetOptStat(0000);
    // gStyle->SetOptFit(1110);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.3);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42);
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLegendTextSize(0.04);
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    size_t dotpos = rootname.Last('.');
    size_t slashpos = rootname.Last('/');
    TString filename = rootname(slashpos + 1, dotpos - slashpos - 1);
    TFile *infile = new TFile(rootname, "READ");
    TTree *tr = (TTree *)infile->Get("decode_data");
    Int_t Nentries = (Int_t)tr->GetEntries();

    DataModel2025 hit[25];
    TH1F *HGMip_his[25], *LGMip_his[25];

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            tr->SetBranchAddress(Form("Hit_%d_%d", i + 1, j + 1), &hit[5 * i + j]);
            HGMip_his[5 * i + j] = new TH1F(Form("HGMIPSpectrum_%d_%d", i + 1, j + 1), Form("his%d;ADC value;counts", 5 * i + j), 200, 4000, 14000);
            LGMip_his[5 * i + j] = new TH1F(Form("LGMIPSpectrum_%d_%d", i + 1, j + 1), Form("his%d;ADC value;counts", 5 * i + j), 200, 400, 1400);
        }
    double energycut = 2000, peakHG, peakLG, tempHG, tempLG;
    int channel, chcount, process;
    for (int n = 0; n < Nentries; n++)
    {
        if (n % (Nentries / 100) == 0 && n != 0)
        {
            process = 100 * n / (Nentries);
            std::cout << "\rProcessing: >>>>>>>>>>>>>>> " << process << "% <<<<<<<<<<<<<<<";
            std::cout.flush();
        }
        tr->GetEntry(n);
        chcount = 0;
        channel = -1;
        tempHG = 0;
        tempLG = 0;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
            {
                int ch = 5 * i + j;
                tempLG = hit[ch].LowGainPeak - hit[ch].LowGainPedestal;
                tempHG = hit[ch].HighGainPeak - hit[ch].HighGainPedestal;
                if (tempHG > energycut)
                {
                    chcount++;
                    channel = ch;
                    peakLG = tempLG;
                    peakHG = tempHG;
                }
            }
        if (chcount != 1)
            continue;
        else
        {
            if (peakHG > HGMip_his[channel]->GetBinCenter(0) && peakHG < HGMip_his[channel]->GetBinCenter(HGMip_his[channel]->GetNbinsX()))
                HGMip_his[channel]->Fill(peakHG);
            if (peakLG > LGMip_his[channel]->GetBinCenter(0) && peakLG < LGMip_his[channel]->GetBinCenter(LGMip_his[channel]->GetNbinsX()))
                LGMip_his[channel]->Fill(peakLG);
        }
    }
    std::cout << std::endl;

    double MPVHG[25], PeakHG[25];
    double MPVLG[25], PeakLG[25];
    TCanvas *canHG = new TCanvas("HGMipSpectra", "canvas", 1600, 900);
    canHG->Divide(5, 5);
    TCanvas *canLG = new TCanvas("LGMipSpectra", "canvas", 1600, 900);
    canLG->Divide(5, 5);
    double parHG[4] = {1000, 8000, 10000, 50};
    double parLG[4] = {100, 800, 1000, 5};
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            TF1 *ff = new TF1("lan_gaus_conv", langaufun, 0, 20000, 4);
            canHG->cd(21 + i - 5 * j);
            HGMip_his[5 * i + j]->Draw();
            parHG[1] = HGMip_his[5 * i + j]->GetMaximumBin() * HGMip_his[5 * i + j]->GetBinWidth(0) + HGMip_his[5 * i + j]->GetBinLowEdge(1);
            ff->SetParameters(parHG);
            HGMip_his[5 * i + j]->Fit(ff, "QR", "", HGMip_his[5 * i + j]->GetMaximumBin() * HGMip_his[5 * i + j]->GetBinWidth(0) + HGMip_his[5 * i + j]->GetBinLowEdge(1) - 1000, HGMip_his[5 * i + j]->GetMaximumBin() * HGMip_his[5 * i + j]->GetBinWidth(0) + HGMip_his[5 * i + j]->GetBinLowEdge(1) + 4000);
            if (HGMip_his[5 * i + j]->GetFunction("lan_gaus_conv"))
                MPVHG[5 * i + j] = HGMip_his[5 * i + j]->GetFunction("lan_gaus_conv")->GetMaximumX();
            PeakHG[5 * i + j] = HGMip_his[5 * i + j]->GetMaximumBin() * HGMip_his[5 * i + j]->GetBinWidth(0) + HGMip_his[5 * i + j]->GetBinLowEdge(1);
            HGMip_his[5 * i + j]->GetXaxis()->SetRangeUser(PeakHG[5 * i + j] - 2 * HGMip_his[5 * i + j]->GetRMS(), PeakHG[5 * i + j] + 5 * HGMip_his[5 * i + j]->GetRMS());
            HGMip_his[5 * i + j]->GetYaxis()->SetRangeUser(0, HGMip_his[5 * i + j]->GetMaximum() * 1.2);

            canLG->cd(21 + i - 5 * j);
            LGMip_his[5 * i + j]->Draw();
            parLG[1] = LGMip_his[5 * i + j]->GetMaximumBin() * LGMip_his[5 * i + j]->GetBinWidth(0) + LGMip_his[5 * i + j]->GetBinLowEdge(1);
            ff->SetParameters(parLG);
            LGMip_his[5 * i + j]->Fit(ff, "QR", "", LGMip_his[5 * i + j]->GetMaximumBin() * LGMip_his[5 * i + j]->GetBinWidth(0) + LGMip_his[5 * i + j]->GetBinLowEdge(1) - 100, LGMip_his[5 * i + j]->GetMaximumBin() * LGMip_his[5 * i + j]->GetBinWidth(0) + LGMip_his[5 * i + j]->GetBinLowEdge(1) + 400);
            if (LGMip_his[5 * i + j]->GetFunction("lan_gaus_conv"))
                MPVLG[5 * i + j] = LGMip_his[5 * i + j]->GetFunction("lan_gaus_conv")->GetMaximumX();
            PeakLG[5 * i + j] = LGMip_his[5 * i + j]->GetMaximumBin() * LGMip_his[5 * i + j]->GetBinWidth(0) + LGMip_his[5 * i + j]->GetBinLowEdge(1);
            LGMip_his[5 * i + j]->GetXaxis()->SetRangeUser(PeakLG[5 * i + j] - 2 * LGMip_his[5 * i + j]->GetRMS(), PeakLG[5 * i + j] + 5 * LGMip_his[5 * i + j]->GetRMS());
            LGMip_his[5 * i + j]->GetYaxis()->SetRangeUser(0, LGMip_his[5 * i + j]->GetMaximum() * 1.2);
        }
    TFile *fout = new TFile("./5x5MIP-" + filename + ".root", "RECREATE");
    fout->cd();
    canHG->Write();
    canLG->Write();
    for (int i = 0; i < 25; i++)
    {
        HGMip_his[i]->Write();
        LGMip_his[i]->Write();
    }
    fout->Close();
    std::ofstream outfile("MIPpeak.txt");
    outfile << "HG MIP peak value:\n";
    for (int i = 0; i < 25; i++)
        outfile << PeakHG[i] << ", ";
    outfile << "\nHG MIP peak fit value:\n";
    for (int i = 0; i < 25; i++)
        outfile << MPVHG[i] << ", ";
    outfile << "\nLG MIP peak value:\n";
    for (int i = 0; i < 25; i++)
        outfile << PeakLG[i] << ", ";
    outfile << "\nLG MIP peak fit value:\n";
    for (int i = 0; i < 25; i++)
        outfile << MPVLG[i] << ", ";
    outfile.close();
    std::cout << "generate root file: " << "./5x5MIP-" + filename + ".root" << ", and parameter file: MIPpeak.txt" << std::endl;
}

#if !defined(__CINT__) && !defined(__CLING__)
int main(int argv, char *argc[])
{
    if (argv < 2)
    {
        cout << "Usage: ./Draw5x5MIPSpectra decode file" << endl;
        return 1;
    }
    TString rootname = argc[1];
    Draw5x5MIPSpectra(rootname);
    return 0;
}
#endif