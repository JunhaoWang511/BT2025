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
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPaveStats.h>
// input 5x5Noise-decode.root
void DrawState(TString filename)
{
    gStyle->SetOptStat(0000);
    gStyle->SetOptTitle(0);
    // gStyle->SetOptFit(1110);
    gStyle->SetPadLeftMargin(0.12);
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
    gStyle->SetTitleOffset(1, "y");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    TFile *infile = new TFile(filename.Data(), "READ");
    TH1F *his_Lnoise, *his_Hnoise;
    TGraph *gr_ratio;
    TLatex *tex = new TLatex();
    tex->SetTextFont(42);
    tex->SetTextSize(0.08);
    tex->SetTextColor(kBlue);
    tex->SetTextAlign(20);
    his_Lnoise = (TH1F *)infile->Get("LGnoise_3_3");
    his_Hnoise = (TH1F *)infile->Get("HGnoise_3_3");
    gr_ratio = (TGraph *)infile->Get("GainRatio_3_3");
    double lpeak = 3264.4, lped = 2414.79;
    double hpeak = 10911.4, hped = 2356.82;
    double electric_factorLG = 1.5;
    double electric_factorHG = 14.81;
    double gain_ratio = gr_ratio->GetFunction("fratio")->GetParameter(0);
    double hnoise = his_Hnoise->GetFunction("gaus")->GetParameter(2);
    double lnoise = his_Lnoise->GetFunction("gaus")->GetParameter(2);
    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    his_Lnoise->Draw();
    tex->DrawLatexNDC(0.3, 0.7, Form("#sigma=%.2f fC", lnoise / electric_factorHG * gain_ratio));
    tex->DrawLatexNDC(0.3, 0.6, Form("~%.2f MeV", lnoise / (lpeak - lped) * 178));
    can1->SaveAs("LGnoise.png");
    can1->SaveAs("LGnoise.pdf");

    TCanvas *can2 = new TCanvas("can2", "can2", 900, 600);
    his_Hnoise->Draw();
    tex->DrawLatexNDC(0.3, 0.7, Form("#sigma=%.2f fC", hnoise / electric_factorHG));
    tex->DrawLatexNDC(0.3, 0.6, Form("~%.2f MeV", hnoise / (hpeak - hped) * 178));
    can2->SaveAs("HGnoise.png");
    can2->SaveAs("HGnoise.pdf");

    TCanvas *can3 = new TCanvas("can3", "can3", 900, 600);
    gr_ratio->Draw("AP");
    tex->DrawLatexNDC(0.6, 0.4, Form("Slope=%.2f", gain_ratio));
    TLegend *leg = new TLegend(0.7, 0.6, 0.8, 0.8);
    leg->AddEntry(gr_ratio, Form("#frac{high gain}{low gain}"), "p");
    leg->SetFillStyle(0);
    leg->SetTextColor(kBlue);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.08);
    // leg->Draw();
    can3->SaveAs("gainratio.png");
    can3->SaveAs("gainratio.pdf");
    std::cout << Form("equivalent noise energy: HG= %.2f, LG= %.2f", hnoise / (hpeak - hped) * 178, lnoise / (lpeak - lped) * 178) << std::endl;
}
// draw noise, pedestal, temperature and gain ratio stability
void DrawStateStability()
{
    gStyle->SetOptStat(0000);
    gStyle->SetOptTitle(0);
    // gStyle->SetOptFit(1110);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.3);
    gStyle->SetPadTickX(1);
    // gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gStyle->SetTextFont(42);
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLegendTextSize(0.04);
    gStyle->SetTitleOffset(1, "y");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetMarkerStyle(8);
    gStyle->SetMarkerSize(1);

    std::vector<double> pedestal = {2320.09, 2323.89, 2322.1, 2323.52, 2321.15, 2320.9, 2322.22};
    std::vector<double> noise = {26.5333, 26.3687, 26.0047, 25.5379, 25.1639, 24.8113, 25.6121};
    std::vector<double> ratio = {10.0069, 10.0369, 10.0423, 10.0191, 10.0406, 10.0182, 9.88754};
    std::vector<double> temperature = {28.6981, 28.7027, 28.5784, 28.545, 28.485, 28.905, 28.8465};
    std::vector<TString> runID = {"1202_25", "1203_04", "1204_22", "1205_18", "1206_29", "1207_25", "1208_02"};
    int size = runID.size();
    if (pedestal.size() != size || noise.size() != size || ratio.size() != size || temperature.size() != size)
        std::cerr << "size does not match!" << std::endl;
    TGraph *gr_ped = new TGraph(size);
    TGraph *gr_noise = new TGraph(size);
    TGraph *gr_ratio = new TGraph(size);
    TGraph *gr_temp = new TGraph(size);
    gr_ped->SetTitle(";Date_Run;HG pedestal[ADC value]");
    gr_noise->SetTitle(";Date_Run;HG noise[ADC value]");
    gr_ratio->SetTitle(";Date_Run;HG/LG ratio");
    gr_temp->SetTitle(";Date_Run;Temperature/#circC");
    for (int i = 0; i < size; i++)
    {
        gr_ped->SetPoint(i, i + 1, pedestal[i]);
        gr_noise->SetPoint(i, i + 1, noise[i]);
        gr_ratio->SetPoint(i, i + 1, ratio[i]);
        gr_temp->SetPoint(i, i + 1, temperature[i]);
        gr_ped->GetXaxis()->SetLabelSize(0);
        gr_noise->GetXaxis()->SetLabelSize(0);
        gr_ratio->GetXaxis()->SetLabelSize(0);
        gr_temp->GetXaxis()->SetLabelSize(0);
        // gr_ped->GetXaxis()->SetBinLabel(i + 1, runID.at(i).Data());
        // gr_noise->GetXaxis()->SetBinLabel(i + 1, runID.at(i).Data());
        // gr_ratio->GetXaxis()->SetBinLabel(i + 1, runID.at(i).Data());
        // gr_temp->GetXaxis()->SetBinLabel(i + 1, runID.at(i).Data());
    }

    // TH1F *hframe = new TH1F("hframe", ";;Y Value", size, 0.5, size + 0.5);
    // hframe->SetStats(0);
    // for (int i = 0; i < size; i++)
    // {
    //     hframe->GetXaxis()->SetBinLabel(i + 1, runID.at(i).Data());
    // }

    TText *t = new TText();
    t->SetTextAlign(12);
    t->SetTextSize(0.04);
    t->SetTextAngle(90);

    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    gr_ped->GetYaxis()->SetRangeUser(2200, 2400);
    gr_ped->GetXaxis()->SetLimits(0,9);
    gr_ped->Draw("ap");
    gPad->Update();
    for (int i = 0; i < size; i++)
        t->DrawText(i + 1, 2200 - 0.23 * (2400 - 2200), runID.at(i).Data());
    can1->SaveAs("ped_stable.png");
    can1->SaveAs("ped_stable.pdf");

    TCanvas *can2 = new TCanvas("can2", "can2", 900, 600);
    gr_noise->GetYaxis()->SetRangeUser(0, 50);
    gr_noise->GetXaxis()->SetLimits(0,9);
    gr_noise->Draw("ap");
    gPad->Update();
    for (int i = 0; i < size; i++)
        t->DrawText(i + 1, 0 - 0.23 * (50 - 0), runID.at(i).Data());
    can2->SaveAs("noise_stable.png");
    can2->SaveAs("noise_stable.pdf");

    TCanvas *can3 = new TCanvas("can3", "can3", 900, 600);
    gr_ratio->GetYaxis()->SetRangeUser(0, 20);
    gr_ratio->GetXaxis()->SetLimits(0,9);
    gr_ratio->Draw("ap");
    gPad->Update();
    for (int i = 0; i < size; i++)
        t->DrawText(i + 1, 0 - 0.23 * (20 - 0), runID.at(i).Data());
    can3->SaveAs("gain_stable.png");
    can3->SaveAs("gain_stable.pdf");

    TCanvas *can4 = new TCanvas("can4", "can4", 900, 600);
    gr_temp->GetYaxis()->SetRangeUser(20, 40);
    gr_temp->GetXaxis()->SetLimits(0,9);
    gr_temp->Draw("ap");
    gPad->Update();
    for (int i = 0; i < size; i++)
        t->DrawText(i + 1, 20 - 0.23 * (40 - 20), runID.at(i).Data());
    can4->SaveAs("temp_stable.png");
    can4->SaveAs("temp_stable.pdf");

    TCanvas *can5 = new TCanvas("can5", "can5", 900, 600);
    TGraph *gr_ped1 = (TGraph *)gr_ped->Clone();
    gr_ped1->SetMarkerColor(kRed);
    gr_ped1->GetYaxis()->SetLabelColor(kRed);
    gr_ped1->GetYaxis()->SetTitleColor(kRed);
    gr_ped1->GetXaxis()->SetLimits(0,9);
    gr_ped1->Draw("ap");
    TGraph *gr_noise1 = (TGraph *)gr_noise->Clone();
    gr_noise1->SetMarkerColor(kBlue);
    for (int i = 0; i < size; i++)
        t->DrawText(i + 1, 2200 - 0.23 * (2400 - 2200), runID.at(i).Data());
    for (int i = 0; i < gr_noise1->GetN(); i++)
    {
        gr_noise1->SetPoint(i, i + 1, 2200 + (2400 - 2200) / (50 - 0) * gr_noise1->GetPointY(i));
    }
    gr_noise1->Draw("p same");
    TGaxis *Axis = new TGaxis(gr_noise1->GetN() + 2, 2200, gr_noise1->GetN() + 2, 2400, 0, 50, 510, "+L");
    Axis->SetTitle("HG noise[ADC value]");
    Axis->SetTitleColor(kBlue);
    Axis->SetTitleFont(42);
    Axis->SetTitleSize(0.05);
    Axis->SetTitleOffset(0.8);
    Axis->SetLabelColor(kBlue);
    Axis->SetLabelFont(42);
    Axis->Draw();
    TLegend *leg = new TLegend(0.6, 0.25, 0.8, 0.45);
    leg->SetTextSize(0.05);
    leg->AddEntry(gr_ped1, "Baseline", "p");
    leg->AddEntry(gr_noise1, "Noise", "p");
    leg->Draw();
    can5->SaveAs("PedNoise_stable.png");
    can5->SaveAs("PedNoise_stable.pdf");
}