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

void resolution_E()
{
    gStyle->SetOptStat(0);
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
    gStyle->SetTitleOffset(0.8, "xy");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    double energy2510[11] = {0.2, 0.3, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double energyres2510[11] = {3.29, 2.88, 2.59, 2.48, 2.36, 2.25, 2.22, 2.22, 2.14, 2.14, 2.25};
    double energy2512[5] = {0.4, 0.5, 0.6, 0.8, 1};
    double energyres2512[5] = {3.26, 2.65, 2.72, 2.4, 2.09};
    // double energy25[11] = {0.2, 0.3, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    // double energyRec[11] = {0.2077, 0.3111, 0.5079, 0.609, 0.8075, 1.0197, 1.4963, 1.9872, 2.4724, 2.956, 3.3757};
    double energy2511[11] = {0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double energyRec2511[11] = {0.2040, 0.4034, 0.5001, 0.6039, 0.8064, 1.0075, 1.5005, 2.0093, 2.4932, 2.9866, 3.4665};
    double energyres2511[11] = {4.47, 3.26, 2.87, 2.7, 2.35, 2.22, 1.98, 1.83, 1.83, 1.77, 1.74};
    double energy2511online[9] = {0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3.5};
    double energyRec2511online[9] = {0.1904, 0.3843, 0.5819, 0.7763, 0.9787, 1.47, 1.9654, 2.4559, 3.4285};
    double energyres2511online[9] = {6.08, 3.62, 3.05, 2.7, 2.51, 2.23, 2.06, 2.03, 1.96};
    // double energyresMC[11] = {3.24, 2.86, 2.45, 2.36, 2.27, 2.25, 2.22, 2.2, 2.2, 2.18, 2.27};
    // double energyMC[11] = {0.2, 0.3, 0.5, 0.6, 0.8, 1, 1.4, 2, 2.5, 3, 3.5};
    double energyres24[3] = {2.43, 2.33, 2.3};
    double energy24[3] = {1, 1.5, 2};
    TGraph *gr_en2511 = new TGraph(11, energy2511, energyres2511);
    TGraph *gr_en2511online = new TGraph(9, energy2511online, energyres2511online);
    TGraph *gr_en2512 = new TGraph(5, energy2512, energyres2512);
    TGraph *gr_en2510 = new TGraph(11, energy2510, energyres2510);
    TGraph *gr_en24 = new TGraph(3, energy24, energyres24);
    // TGraph *gr_enMC = new TGraph(11, energyMC, energyresMC);
    gr_en24->SetMarkerStyle(8);
    gr_en24->SetMarkerSize(1);
    gr_en24->SetMarkerColor(kBlack);
    gr_en2512->SetMarkerStyle(8);
    gr_en2512->SetMarkerSize(1);
    gr_en2512->SetMarkerColor(kGreen);
    // gr_enMC->SetMarkerStyle(8);
    // gr_enMC->SetMarkerSize(1);
    // gr_enMC->SetMarkerColor(kGreen);
    gr_en2510->SetMarkerStyle(8);
    gr_en2510->SetMarkerSize(1);
    gr_en2510->SetMarkerColor(kBlue);
    gr_en2511online->SetMarkerStyle(8);
    gr_en2511online->SetMarkerSize(1);
    gr_en2511online->SetMarkerColor(kBlue);
    gr_en2511->SetMarkerStyle(8);
    gr_en2511->SetMarkerSize(1);
    gr_en2511->SetMarkerColor(kRed);
    gr_en2511->SetTitle(";E_{beam}[GeV];Energy Resolution[%]");
    gr_en2511->GetYaxis()->SetRangeUser(1.5, 5);

    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    gr_en2511->Draw("ap");
    // gr_en2511online->Draw("p same");
    // gr_en2512->Draw("p same");
    // gr_en2510->Draw("p same");
    // gr_en24->Draw("p same");
    // gr_enMC->Draw("p same");
    TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg->SetTextFont(42);
    leg->SetTextSize(0.06);
    // leg->AddEntry(gr_en24, "2024", "p");
    // leg->AddEntry(gr_en2510, "2025 autumn", "p");
    // leg->AddEntry(gr_en2511, "2025 winter", "p");
    leg->AddEntry(gr_en2511, "Single fitting", "p");
    leg->AddEntry(gr_en2511online, "Pipeline fitting", "p");
    // leg->AddEntry(gr_en2512, "winter (front)", "p");
    // leg->AddEntry(gr_enMC, "MC result", "p");
    // leg->Draw();
    can1->SaveAs("EnergyResolution.png");
    can1->SaveAs("EnergyResolution.pdf");

    TCanvas *can2 = new TCanvas("can2", "can2", 900, 600);
    // TGraph *gr_linear = new TGraph(11, energy2511, energyRec2511);
    TGraph *gr_linear = new TGraph(9, energy2511online, energyRec2511online);
    gr_linear->GetXaxis()->SetRangeUser(0, 3.8);
    gr_linear->GetYaxis()->SetRangeUser(0, 3.8);
    gr_linear->SetTitle(";E_{beam}[GeV];E_{rec}[GeV]");
    gr_linear->Draw("ap");
    gr_linear->SetMarkerStyle(8);
    gr_linear->SetMarkerSize(1);
    gr_linear->SetMarkerColor(kBlue);
    gr_linear->Fit("pol1");
    gr_linear->GetFunction("pol1")->SetLineColor(kRed);
    can2->SaveAs("EnergyLinearity.png");
    can2->SaveAs("EnergyLinearity.pdf");
}
void resolution_P()
{
    gStyle->SetOptStat(0);
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

    double energy25[9] = {0.5, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    // double positionres25[9] = {7.08, 6.32, 5.82, 5.11, 4.56, 4.14, 3.94, 3.65, 3.48};
    double positionres25[9] = {7.11, 6.46, 5.97, 5.48, 4.88, 4.38, 4.03, 3.69, 3.56};
    double energy24[6] = {1, 1.5, 2, 2.5, 3, 3.5};
    double positionres24[6] = {5.4, 4.5, 3.98, 3.7, 3.5, 3.3};
    double positionresMC[6] = {5.3, 4.46, 3.84, 3.6, 3.44, 3.12};

    TGraph *gr_pos25 = new TGraph(9, energy25, positionres25);
    gr_pos25->SetMarkerStyle(8);
    gr_pos25->SetMarkerSize(1);
    gr_pos25->SetMarkerColor(kRed);
    gr_pos25->SetTitle(";energy[GeV];position resolution[mm]");
    gr_pos25->GetYaxis()->SetRangeUser(2, 8);

    TGraph *gr_pos24 = new TGraph(6, energy24, positionres24);
    gr_pos24->SetMarkerStyle(8);
    gr_pos24->SetMarkerSize(1);
    gr_pos24->SetMarkerColor(kBlue);

    TGraph *gr_MC = new TGraph(6, energy24, positionresMC);
    gr_MC->SetMarkerStyle(8);
    gr_MC->SetMarkerSize(1);
    gr_MC->SetMarkerColor(kGreen);

    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    gr_pos25->Draw("ap");
    // gr_pos24->Draw("p same");
    // gr_MC->Draw("p same");

    // TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.7);
    // leg->SetTextFont(42);
    // leg->SetTextSize(0.06);
    // leg->AddEntry(gr_pos24, "2024 result", "p");
    // leg->AddEntry(gr_pos25, "2025 result", "p");
    // leg->AddEntry(gr_MC, "MC result", "p");
    // leg->Draw();
    can1->SaveAs("PositionResolution.png");
    can1->SaveAs("PositionResolution.pdf");
}
// 在线提取和拟合点数以及CMS时间提取方法的结果
void resolution_ETP()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
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

    double energy[10] = {0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double time_fit[10] = {1207, 289, 294, 273, 258, 233, 213, 206, 203, 203};
    double time_fit10p[10] = {404, 237, 204, 168, 148, 128, 119, 114, 114, 117};
    double time_offline[10] = {486, 283, 242, 223, 210, 190, 177, 170, 165, 160};
    double time_online[10] = {577, 285, 243, 223, 209, 190, 176, 170, 0, 0};
    double time_pipeline_improvement[10] = {543, 273, 218, 197, 188, 195, 185, 208, 176, 165};
    double time_pipeline_improvement_LG[10] = {538, 275, 210, 184, 167, 153, 146, 141, 140, 137};
    double time_CMS[10] = {571, 295, 250, 211, 188, 164, 148, 141, 138, 139};

    double energy_fit[10] = {4.47, 3.26, 2.7, 2.35, 2.22, 1.98, 1.83, 1.83, 1.77, 1.74};
    double energy_fit10p[10] = {4.45, 3.33, 2.75, 2.37, 2.25, 1.98, 1.85, 1.82, 1.75, 1.69};
    double energy_offline[10] = {0, 0, 3.06, 2.67, 2.53, 2.23, 2.05, 2.03, 1.97, 2};
    double energy_online[10] = {6.08, 3.62, 3.05, 2.7, 2.51, 2.23, 2.06, 2.03, 0, 1.96};
    // 每个通道各自模板，模板插值，波形拼接，拟合上升沿
    double energy_pipeline_improvement[10] = {0, 3.33, 2.75, 2.42, 2.27, 2.05, 1.89, 1.89, 1.86, 1.85};
    double energy_pipeline_improvement_LG[10] = {0, 3.49, 2.89, 2.52, 2.37, 2.06, 1.92, 1.89, 1.86, 1.86};

    double energyresponse_fit20p[10] = {0.204, 0.4034, 0.6039, 0.8064, 1.0075, 1.5005, 2.0093, 2.4932, 2.9866, 3.4665};
    double energyresponse_fit10p[10] = {0.2034, 0.407, 0.6106, 0.8157, 1.0193, 1.5272, 2.0357, 2.54, 3.0411, 3.5141};

    double postion_fit10p[10] = {0, 8.5, 7.26, 6.34, 5.67, 4.72, 4.28, 3.78, 3.51, 3.24};

    TGraph *gr_timeCMS = new TGraph(10, energy, time_CMS);
    gr_timeCMS->SetMarkerStyle(8);
    gr_timeCMS->SetMarkerSize(1);
    gr_timeCMS->SetMarkerColor(kRed);
    gr_timeCMS->SetTitle(";Energy[GeV];Time resolution[ps]");
    gr_timeCMS->GetYaxis()->SetRangeUser(150, 600);

    TGraph *gr_timeonline = new TGraph(10, energy, time_online);
    gr_timeonline->SetMarkerStyle(8);
    gr_timeonline->SetMarkerSize(1);
    gr_timeonline->SetMarkerColor(kRed);
    gr_timeonline->SetTitle(";Energy[GeV];Time resolution[ps]");
    gr_timeonline->GetYaxis()->SetRangeUser(150, 600);

    TGraph *gr_timeoffline = new TGraph(10, energy, time_offline);
    gr_timeoffline->SetMarkerStyle(kFullTriangleDown);
    gr_timeoffline->SetMarkerSize(1);
    gr_timeoffline->SetMarkerColor(kBlue);
    gr_timeoffline->SetTitle(";Energy[GeV];Time resolution[ps]");

    TGraph *gr_timefit = new TGraph(10, energy, time_fit);
    gr_timefit->SetMarkerStyle(8);
    gr_timefit->SetMarkerSize(1);
    gr_timefit->SetMarkerColor(kBlack);
    gr_timefit->SetTitle(";Energy[GeV];Time resolution[ps]");

    TGraph *gr_timefit10p = new TGraph(10, energy, time_fit10p);
    gr_timefit10p->SetMarkerStyle(8);
    gr_timefit10p->SetMarkerSize(1);
    gr_timefit10p->SetMarkerColor(kRed);
    gr_timefit10p->SetTitle(";E_{beam}[GeV];Time resolution[ps]");

    TGraph* gr_timepipelineimprovement = new TGraph(10, energy, time_pipeline_improvement);
    gr_timepipelineimprovement->SetMarkerStyle(kFullTriangleUp);
    gr_timepipelineimprovement->SetMarkerSize(1.5);
    gr_timepipelineimprovement->SetMarkerColor(kBlue);
    gr_timepipelineimprovement->SetTitle(";E_{beam}[GeV];Time resolution[ps]");

    TGraph *gr_timepipelineimprovement_LG = new TGraph(10, energy, time_pipeline_improvement_LG);
    gr_timepipelineimprovement_LG->SetMarkerStyle(kFullTriangleDown);
    gr_timepipelineimprovement_LG->SetMarkerSize(1.5);
    gr_timepipelineimprovement_LG->SetMarkerColor(kRed);
    gr_timepipelineimprovement_LG->SetTitle(";E_{beam}[GeV];Time resolution[ps]");

    TGraph *gr_energyonline = new TGraph(10, energy, energy_online);
    gr_energyonline->SetMarkerStyle(8);
    gr_energyonline->SetMarkerSize(1);
    gr_energyonline->SetMarkerColor(kRed);
    gr_energyonline->SetTitle(";Energy[GeV];Energy resolution[%]");
    gr_energyonline->GetYaxis()->SetRangeUser(1.5, 6.1);

    TGraph *gr_energyoffline = new TGraph(10, energy, energy_offline);
    gr_energyoffline->SetMarkerStyle(kFullTriangleDown);
    gr_energyoffline->SetMarkerSize(1);
    gr_energyoffline->SetMarkerColor(kBlue);
    gr_energyoffline->SetTitle(";Energy[GeV];Energy resolution[%]");

    TGraph *gr_energyfit = new TGraph(10, energy, energy_fit);
    gr_energyfit->SetMarkerStyle(8);
    gr_energyfit->SetMarkerSize(1);
    gr_energyfit->SetMarkerColor(kBlack);
    gr_energyfit->SetTitle(";Energy[GeV];Energy resolution[%]");

    TGraph *gr_energyfit10p = new TGraph(10, energy, energy_fit10p);
    gr_energyfit10p->SetMarkerStyle(8);
    gr_energyfit10p->SetMarkerSize(1);
    gr_energyfit10p->SetMarkerColor(kRed);
    gr_energyfit10p->SetTitle(";Energy[GeV];Energy resolution[%]");

    TGraph* gr_energypipelineimprovement = new TGraph(10, energy, energy_pipeline_improvement);
    gr_energypipelineimprovement->SetMarkerStyle(kFullTriangleUp);
    gr_energypipelineimprovement->SetMarkerSize(1.5);
    gr_energypipelineimprovement->SetMarkerColor(kBlue);
    gr_energypipelineimprovement->SetTitle(";E_{beam};Energy resolution[%]");

    TGraph* gr_energypipelineimprovement_LG = new TGraph(10, energy, energy_pipeline_improvement_LG);
    gr_energypipelineimprovement_LG->SetMarkerStyle(kFullTriangleDown);
    gr_energypipelineimprovement_LG->SetMarkerSize(1.5);
    gr_energypipelineimprovement_LG->SetMarkerColor(kRed);
    gr_energypipelineimprovement_LG->SetTitle(";E_{beam};Energy resolution[%]");

    TGraph *gr_positionfit10p = new TGraph(10, energy, postion_fit10p);
    gr_positionfit10p->SetMarkerStyle(8);
    gr_positionfit10p->SetMarkerSize(1);
    gr_positionfit10p->SetMarkerColor(kRed);
    gr_positionfit10p->SetTitle(";E_{beam};Position resolution[mm]");

    TCanvas *can0 = new TCanvas("can0", "can0", 900, 600);
    gr_positionfit10p->GetYaxis()->SetRangeUser(2, 10);
    gr_positionfit10p->Draw("ap");
    can0->SaveAs("PositionResolution.png");

    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    gr_timefit10p->Draw("ap");
    // gr_timeoffline->Draw("p same");
    gr_timepipelineimprovement->Draw("p same");
    gr_timepipelineimprovement_LG->Draw("p same");
    // gr_timeoffline->Draw("p same");
    TLegend *leg1 = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.06);
    leg1->AddEntry(gr_timefit10p, "offline fit", "p");
    // leg1->AddEntry(gr_timeoffline, "pipeline fit", "p");
    leg1->AddEntry(gr_timepipelineimprovement, "pipeline fit with HG", "p");
    leg1->AddEntry(gr_timepipelineimprovement_LG, "pipeline fit with LG", "p");
    // leg1->AddEntry(gr_timeoffline, "offline TQ", "p");
    leg1->Draw();
    can1->SaveAs("TimeResolution.png");

    TCanvas *can2 = new TCanvas("can2", "can2", 900, 600);
    gr_energyfit->Draw("ap");
    // gr_energyoffline->Draw("p same");
    // gr_energyonline->Draw("p same");
    gr_energypipelineimprovement->Draw("p same");
    gr_energypipelineimprovement_LG->Draw("p same");
    TLegend *leg2 = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.06);
    leg2->AddEntry(gr_energyfit, "offline fit", "p");
    // leg2->AddEntry(gr_energyonline, "pipeline fit", "p");
    leg2->AddEntry(gr_energypipelineimprovement, "pipeline fit with HG", "p");
    leg2->AddEntry(gr_energypipelineimprovement_LG, "pipeline fit with LG", "p");
    // leg2->AddEntry(gr_energyoffline, "offline TQ", "p");
    leg2->Draw();
    can2->SaveAs("EnergyResolution.png");

    TCanvas *can3 = new TCanvas("can3", "can3", 900, 600);
    gr_energyfit->Draw("ap");
    gr_energyfit10p->Draw("p same");
    TLegend *leg3 = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.06);
    leg3->AddEntry(gr_energyfit, "20 point", "p");
    leg3->AddEntry(gr_energyfit10p, "10 point", "p");
    leg3->Draw();
    can3->SaveAs("EnergyImprovement.png");

    TCanvas *can4 = new TCanvas("can4", "can4", 900, 600);
    gr_timefit->Draw("ap");
    gr_timefit10p->Draw("p same");
    TLegend *leg4 = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg4->SetTextFont(42);
    leg4->SetTextSize(0.06);
    leg4->AddEntry(gr_timefit, "20 point", "p");
    leg4->AddEntry(gr_timefit10p, "10 point", "p");
    // leg4->Draw();
    can4->SaveAs("TimeImprovement.png");

    TCanvas *can5 = new TCanvas("can5", "can5", 900, 600);
    gr_timeCMS->GetYaxis()->SetRangeUser(100, 600);
    gr_timeCMS->Draw("ap");
    gr_timefit10p->Draw("p same");
    gr_timeoffline->Draw("p same");
    gr_timefit10p->SetMarkerColor(kBlack);
    gr_timeoffline->SetMarkerColor(kBlue);
    gr_timeCMS->SetMarkerColor(kRed);
    TLegend *leg5 = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg5->SetTextFont(42);
    leg5->SetTextSize(0.06);
    leg5->AddEntry(gr_timefit10p, "offline fit", "p");
    leg5->AddEntry(gr_timeoffline, "online TQ", "p");
    leg5->AddEntry(gr_timeCMS, "CMS method", "p");
    leg5->Draw();
    can5->SaveAs("TimeCMSmethod.png");

    double point[6] = {5, 8, 10, 15, 20, 40};
    double time_point[6] = {643, 184, 145, 179, 217, 257};
    double energy_point[6] = {5.27638, 2.26757, 2.25285, 2.25506, 2.24175, 2.25116};

    TGraph *gr_pointT = new TGraph(6, point, time_point);
    gr_pointT->SetTitle(";Point number;Time resolution[ps]");
    gr_pointT->SetMarkerStyle(20);
    TGraph *gr_pointE = new TGraph(6, point, energy_point);
    gr_pointE->SetTitle(";Point number;Energy resolution[%]");
    gr_pointE->SetMarkerStyle(20);

    TCanvas *can6 = new TCanvas("can6", "can6", 900, 600);
    gr_pointE->Draw("ap");
    can6->SaveAs("point-time_curve.png");
    TCanvas *can7 = new TCanvas("can7", "can7", 900, 600);
    gr_pointT->Draw("ap");
    can7->SaveAs("point-energy_curve.png");

    TCanvas *can8 = new TCanvas("can8", "can8", 900, 600);
    TGraph *gr_energylinear20p = new TGraph(10, energy, energyresponse_fit20p);
    gr_energylinear20p->GetYaxis()->SetRangeUser(0, 3.8);
    gr_energylinear20p->GetXaxis()->SetRangeUser(0, 3.8);
    gr_energylinear20p->GetYaxis()->SetTitleOffset(0.8);
    gr_energylinear20p->SetMarkerStyle(20);
    gr_energylinear20p->SetTitle(";Beam momentum[GeV/c];E_{rec}[GeV]");
    gr_energylinear20p->Draw("ap");
    gr_energylinear20p->Fit("pol1", "R", "", *std::min_element(energy, energy + sizeof(energy) / sizeof(double)), *std::max_element(energy, energy + sizeof(energy) / sizeof(double)));
    gr_energylinear20p->GetFunction("pol1")->SetLineColor(kBlack);
    gr_energylinear20p->GetFunction("pol1")->SetLineStyle(1);
    TGraph *gr_energylinear10p = new TGraph(10, energy, energyresponse_fit10p);
    gr_energylinear10p->SetMarkerStyle(20);
    gr_energylinear10p->SetMarkerColor(kRed);
    gr_energylinear10p->SetTitle(";Beam momentum[GeV/c];E_{rec}[GeV]");
    // 画第二坐标轴
    double xmin, xmax, ymin1, ymax1, ymin2, ymax2;
    xmin = 0, xmax = 3.8, ymin1 = 0, ymax1 = 3.8, ymin2 = -0.2, ymax2 = 3.6;
    double scale = (ymax1 - ymin1) / (ymax2 - ymin2);
    for (int i = 0; i < gr_energylinear10p->GetN(); i++)
    {
        double x, y;
        gr_energylinear10p->GetPoint(i, x, y);
        gr_energylinear10p->SetPoint(i, x, ymin1 + (y - ymin2) * scale);
    }
    gr_energylinear10p->Draw("p same");
    gr_energylinear10p->Fit("pol1", "R", "", *std::min_element(energy, energy + sizeof(energy) / sizeof(double)), *std::max_element(energy, energy + sizeof(energy) / sizeof(double)));
    gr_energylinear10p->GetFunction("pol1")->SetLineColor(kRed);
    gr_energylinear10p->GetFunction("pol1")->SetLineStyle(1);

    double uxmax = gPad->GetUxmax();
    double uymin = gPad->GetUymin();
    double uymax = gPad->GetUymax();
    // TGaxis *axis = new TGaxis(uxmax, uymin, uxmax, uymax, ymin2, ymax2, 510, "+L");
    TGaxis *axis = new TGaxis(3.83, 0, 3.83, 3.8, ymin2, ymax2, 510, "+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetLabelFont(42);
    axis->SetTitleColor(kRed);
    axis->SetTitle("E_{rec}[GeV]");
    axis->SetTitleFont(42);
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(0.8);
    axis->Draw();
    gPad->SetTicky(0);

    TLegend *leg6 = new TLegend(0.15, 0.6, 0.35, 0.8);
    leg6->SetTextSize(0.06);
    leg6->SetTextFont(42);
    leg6->AddEntry(gr_energylinear20p, "20 point", "p");
    leg6->AddEntry(gr_energylinear10p, "10 point", "p");
    leg6->Draw();

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.05);
    tex->SetTextColor(kBlack);
    tex->DrawLatexNDC(0.4, 0.73, Form("#chi^{2}=%.2e", gr_energylinear20p->GetFunction("pol1")->GetChisquare()));
    tex->SetTextColor(kRed);
    tex->DrawLatexNDC(0.4, 0.63, Form("#chi^{2}=%.2e", gr_energylinear10p->GetFunction("pol1")->GetChisquare()));
    can8->SaveAs("EnergyResponse.png");
}
void resolution_T()
{
    gStyle->SetOptStat(0);
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

    double energy[5] = {1, 2, 2.5, 3, 3.5};
    double time_offline[5] = {261, 216, 206, 203, 203};
    double time_online[5] = {220, 177, 170, 175, 170};
    double energy24[7] = {0.5, 1, 1.5, 2, 2.5, 3, 3.5};
    double time_offline24[7] = {363, 289, 289, 285, 276, 261, 253};

    TGraph *gr_timeoffline = new TGraph(5, energy, time_offline);
    gr_timeoffline->SetMarkerStyle(8);
    gr_timeoffline->SetMarkerSize(1);
    gr_timeoffline->SetMarkerColor(kRed);
    gr_timeoffline->SetTitle(";Energy[GeV];Time resolution[ps]");
    gr_timeoffline->GetYaxis()->SetRangeUser(180, 280);

    TGraph *gr_timeoffline24 = new TGraph(7, energy24, time_offline24);
    gr_timeoffline24->SetMarkerStyle(8);
    gr_timeoffline24->SetMarkerSize(1);
    gr_timeoffline24->SetMarkerColor(kBlack);
    gr_timeoffline24->SetTitle(";Energy[GeV];Time resolution[ps]");
    gr_timeoffline24->GetYaxis()->SetRangeUser(180, 380);

    TGraph *gr_timeonline = new TGraph(5, energy, time_online);
    gr_timeonline->SetMarkerStyle(8);
    gr_timeonline->SetMarkerSize(1);
    gr_timeonline->SetMarkerColor(kBlue);
    gr_timeonline->SetTitle(";Energy[GeV];Time resolution[ps]");

    TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
    // gr_timeoffline->Draw("ap");
    gr_timeoffline24->Draw("ap");
    gr_timeoffline->Draw("p same");
    // gr_timeonline->Draw("p same");

    TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.7);
    leg->SetTextFont(42);
    leg->SetTextSize(0.06);
    // leg->AddEntry(gr_timeoffline, "Single fitting", "p");
    // leg->AddEntry(gr_timeonline, "Pipeline fitting", "p");
    leg->AddEntry(gr_timeoffline24, "2024", "p");
    leg->AddEntry(gr_timeoffline, "2025", "p");
    leg->Draw();
    can1->SaveAs("TimeResolution.png");
    can1->SaveAs("TimeResolution.pdf");
}

// draw QT online process performance
void LED_resolution()
{
    gStyle->SetOptStat(0);
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
    // winter
    double template_rec[6] = {2.25, 2.29, 2.34, 2.36, 2.48, 2.57};
    double peak_rec[6] = {2.25, 2.4, 2.42, 2.51, 2.77, 2.92};
    double TQ_rec[6] = {2.52, 2.6, 2.62, 2.66, 2.76, 2.9};
    // 10个点拟合+模板插值
    double TQ_rec2temp[6] = {2.33, 2.37, 2.43, 2.48, 2.57, 2.62};
    double frec[6] = {0, 1, 3, 5, 9, 14};
    // autumn
    // double template_rec[5] = {2.24, 2.37, 2.69, 3.51, 3.98};
    // double TQ_rec[5] = {2.36, 2.43, 2.52, 2.97, 3.23};
    // double frec[5] = {0, 1, 3, 9, 14};
    TGraph *gr_template = new TGraph(6, frec, template_rec);
    gr_template->SetMarkerStyle(kFullCircle);
    gr_template->SetMarkerSize(1.2);
    gr_template->SetMarkerColorAlpha(kRed, 0.5);
    TGraph *gr_peak = new TGraph(6, frec, peak_rec);
    gr_peak->SetMarkerStyle(kFullTriangleUp);
    gr_peak->SetMarkerSize(1.2);
    gr_peak->SetMarkerColorAlpha(kGreen, 0.5);
    TGraph *gr_pipeline = new TGraph(6, frec, TQ_rec);
    gr_pipeline->SetMarkerStyle(kFullTriangleDown);
    gr_pipeline->SetMarkerSize(1.2);
    gr_pipeline->SetMarkerColorAlpha(kBlue, 0.5);
    TGraph *gr_pipeline1 = new TGraph(6, frec, TQ_rec2temp);
    gr_pipeline1->SetMarkerStyle(kFullSquare);
    gr_pipeline1->SetMarkerSize(1.2);
    gr_pipeline1->SetMarkerColorAlpha(kMagenta, 0.5);
    TCanvas *can = new TCanvas("can", "can", 900, 600);
    gr_template->SetTitle(";LED Frequency[MHz];Energy Resolution[%]");
    gr_template->GetYaxis()->SetRangeUser(2, 3.5);
    gr_template->Draw("ap");
    gr_peak->Draw("p same");
    gr_pipeline->Draw("p same");
    gr_pipeline1->Draw("p same");
    TLegend *leg = new TLegend(0.15, 0.6, 0.6, 0.88);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(gr_template, "Template fitting", "p");
    leg->AddEntry(gr_pipeline, "Pipeline fitting", "p");
    leg->AddEntry(gr_peak, "Peak searching", "p");
    leg->AddEntry(gr_pipeline1, "Improved pipeline fitting", "p");
    leg->Draw();
    can->SaveAs("TQresolution.png");
    can->SaveAs("TQresolution.pdf");
}

// 10/20个点拟合的能量响应线性
void Fit_Eres()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "xy");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    double energy[10] = {0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double energy_fit10p[10] = {4.45, 3.33, 2.75, 2.37, 2.25, 1.98, 1.85, 1.82, 1.75, 1.69};
    TGraph *gr_energy = new TGraph(10, energy, energy_fit10p);
    gr_energy->GetYaxis()->SetRangeUser(1.5, 5);
    gr_energy->SetTitle(";E_{beam}[GeV];Energy resolution[%]");
    gr_energy->SetMarkerStyle(20);
    gr_energy->SetMarkerColor(kBlue);
    TF1 *f_energy = new TF1("fres", "TMath::Sqrt(TMath::Power([0]/x,2)+TMath::Power([1]/TMath::Sqrt(x),2)+TMath::Power([2],2))", 0, 4);
    f_energy->SetParameter(0, 0.002);
    f_energy->SetParameter(1, 0.01);
    f_energy->SetParameter(2, 0.02);

    gr_energy->Fit(f_energy);
    gr_energy->Draw("ap");

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.06);
    tex->SetTextFont(132);
    tex->SetTextColor(kRed);
    double pars[3];
    f_energy->GetParameters(pars);
    TString str = Form("#frac{#sigma_{E}}{E} = #frac{%.2lf%}{E(GeV)} #oplus #frac{%.2lf%}{#sqrt{E(GeV)}} #oplus %.2lf%", pars[0], pars[1], pars[2]);
    std::cout << str << std::endl;
    tex->DrawLatexNDC(0.3, 0.6, str.Data());
}

// 10/20/40个点拟合的时间响应一致性(用MIP修正)
void FitT_unifrom()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetTitleOffset(1, "y");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    double Energy[8] = {0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double deltaT10p[8] = {2390.83, 2396.49, 2414.14, 2394.61, 2378.83, 2348.08, 2293.26, 2235.63};
    double deltaT20p[8] = {2601.1, 2896.71, 2649.52, 2876.69, 3029.61, 2932.5, 2920.77, 2882.84};
    double deltaT40p[8] = {2708.68, 2781.51, 2796.34, 2892.27, 2955.8, 2949.81, 2943.72, 2902.15};

    TGraph *gr_time10p = new TGraph(8, Energy, deltaT10p);
    TGraph *gr_time20p = new TGraph(8, Energy, deltaT20p);
    TGraph *gr_time40p = new TGraph(8, Energy, deltaT40p);
    gr_time10p->SetTitle(";Beam Momentum[GeV];#Delta_{T}[ps]");
    gr_time10p->SetMarkerStyle(20);
    gr_time10p->SetMarkerColor(kRed);
    gr_time20p->SetMarkerStyle(20);
    gr_time20p->SetMarkerColor(kBlue);
    gr_time40p->SetMarkerStyle(20);
    gr_time40p->SetMarkerColor(kGreen);

    gr_time10p->GetYaxis()->SetRangeUser(2000, 3500);
    gr_time10p->Draw("ap");
    // gr_time20p->Draw("p same");
    // gr_time40p->Draw("p same");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->AddEntry(gr_time10p, "10 point", "p");
    leg->AddEntry(gr_time20p, "20 point", "p");
    leg->AddEntry(gr_time40p, "40 point", "p");
    // leg->Draw();
}

// CMS方法提取的时间响应一致性（用MIP修正）
void FitT_CMSunifrom()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetTitleOffset(1, "y");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    double Energy[8] = {0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5};
    double deltaTCMS[8] = {13497.42, 13476.1, 13491.3, 13542.94, 13574.34, 13608.87, 13696.25, 13781.71};

    TGraph *gr_timeCMS = new TGraph(8, Energy, deltaTCMS);
    gr_timeCMS->SetTitle(";Beam Momentum[GeV];#Delta_{T}[ps]");
    gr_timeCMS->SetMarkerStyle(20);
    gr_timeCMS->SetMarkerColor(kBlue);
    gr_timeCMS->GetYaxis()->SetRangeUser(13000, 14500);
    gr_timeCMS->Draw("ap");
}

// pipeline时间判选宽度对提取效率（能量分辨）的影响
void pipeline_Tcut()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetTitleOffset(1, "y");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    double timecut[11] = {6.25, 6.50, 7, 8, 9, 10, 12.5, 15, 17.5, 20, 30};
    double energyres[11] = {2.53, 2.45, 2.43, 2.41, 2.41, 2.41, 2.42, 2.43, 2.45, 2.45, 2.57};

    TGraph *gr = new TGraph(11, timecut, energyres);
    gr->SetTitle(";Time cut width[ns];Energy resolution[%]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->GetYaxis()->SetRangeUser(2.35, 2.6);
    gr->Draw("ap");
}