#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TAxis.h>
#include <TPaveStats.h>
#include <TChain.h>
#include <TLatex.h>
#include <iomanip>
// par1: parameters.txt(H/LG pedestal,noise,gain ratio)
// par2: MIPpeak.txt(H/LG MIP spectra max value)
void generate_calipars(TString channelfile, TString mipfile)
{
    std::fstream infile(channelfile.Data(), std::ios::in);
    std::fstream infile1(mipfile.Data(), std::ios::in);
    double HGped[25], LGped[25], HGnoise[25], LGnoise[25], HLratio[25], HGswitch[25];
    double HGpeak[25], LGpeak[25], HGpeakfit[25], LGpeakfit[25];
    std::string temp;
    char tempc;

    std::getline(infile, temp);
    for (int i = 0; i < 25; i++)
        infile >> HGped[i] >> tempc;
    std::getline(infile, temp);

    std::getline(infile, temp);
    for (int i = 0; i < 25; i++)
        infile >> LGped[i] >> tempc;
    std::getline(infile, temp);

    std::getline(infile, temp);
    for (int i = 0; i < 25; i++)
        infile >> HGnoise[i] >> tempc;
    std::getline(infile, temp);

    std::getline(infile, temp);
    for (int i = 0; i < 25; i++)
        infile >> LGnoise[i] >> tempc;
    std::getline(infile, temp);

    std::getline(infile, temp);
    for (int i = 0; i < 25; i++)
        infile >> HLratio[i] >> tempc;
    std::getline(infile, temp);

    for (int i = 0; i < 25; i++)
        HGswitch[i] = 16000 - HGped[i];

    std::getline(infile1, temp);
    for (int i = 0; i < 25; i++)
        infile1 >> HGpeak[i] >> tempc;
    std::getline(infile1, temp);

    std::getline(infile1, temp);
    for (int i = 0; i < 25; i++)
        infile1 >> HGpeakfit[i] >> tempc;
    std::getline(infile1, temp);

    std::getline(infile1, temp);
    for (int i = 0; i < 25; i++)
        infile1 >> LGpeak[i] >> tempc;
    std::getline(infile1, temp);

    std::getline(infile1, temp);
    for (int i = 0; i < 25; i++)
        infile1 >> LGpeakfit[i] >> tempc;
    std::getline(infile1, temp);

    infile.close();
    infile1.close();

    std::fstream outfile("CaliPara.dat", std::ios::out);
    for (int i = 0; i < 25; i++)
    {
        outfile << std::left << std::setprecision(6) << std::setw(10) << LGped[i] << ' ' << std::setw(10) << LGnoise[i] << ' ' << std::setw(10) << HGped[i] << ' ' << std::setw(10) << HGnoise[i] << ' ' << std::setw(10) << HLratio[i] << ' ' << std::setw(10) << HGswitch[i] << std::endl;
    }
    outfile.close();

    std::fstream outfile1("HGMipPeak.dat", std::ios::out);
    for (int i = 0; i < 25; i++)
    {
        outfile1 << HGpeakfit[i] + HGped[i] << std::endl;
    }
    outfile1.close();

    std::fstream outfile2("LGMipPeak.dat", std::ios::out);
    for (int i = 0; i < 25; i++)
    {
        outfile2 << LGpeakfit[i] + LGped[i] << std::endl;
    }
    outfile2.close();
    std::cout << "Generate files: " << "CaliPara.dat" << ", HGMipPeak.dat" << ", LGMipPeak.dat" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "parameters too few, parameter1: calibration file(pedestal, noise, gain ratio); parameter2: mip file (mip peak)" << std::endl;
        return 1;
    }
    generate_calipars(argv[1], argv[2]);
    return 0;
}