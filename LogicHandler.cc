#pragma once

#include <QFileDialog>

#include "DNARepairSimulation.hh"

#include "LogicHandler.hh"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TProfile.h>

#include <iostream>
#include <cstdio>

int LogicHandler::argc_main = 0;
char** LogicHandler::argv_main = nullptr;


LogicHandler::LogicHandler(QWidget* parent)
    : QObject(parent)
{
}

void LogicHandler::onSimulationRequested(double a, double b, double Nirrep ,double Dz)
{
    DNARepairSimulation sim;
    sim.runSimulation(m_macroFilePath.toStdString(), a, b, Nirrep, Dz);
}

void LogicHandler::onMacSelectionRequested()
{
    m_macroFilePath = QFileDialog::getOpenFileName(dynamic_cast<QWidget*>(parent()), tr("Open Mac"), "C:", tr("Mac Files (*.in *.mac)"));
    emit fileSelected(m_macroFilePath);
}

void LogicHandler::onROOTRequested(const QString& simPath, const QString& expPath, const QString& legendTitle, const int timeScale)
{
    gROOT->Reset();
    gStyle->SetOptStat(0);

    auto* c1 = new TCanvas("c1", "", 200, 20, 800, 600);

    FILE* fsim = fopen(simPath.toLocal8Bit(), "r");
    if (!fsim) {
        std::cerr << "Failed to open simulation data file" << std::endl;
        return;
    }

    Float_t tsim, ysim;
    Float_t Tsim[1000], Ysim[1000];
    Float_t maxY = -1e9;
    int nlines = 0;

    while (fscanf(fsim, " %f %f", &tsim, &ysim) == 2) {
        if (maxY < ysim) maxY = ysim;
        Tsim[nlines] = tsim;
        Ysim[nlines] = ysim;
        nlines++;
    }
    fclose(fsim);

    std::cout << "Maximum value of function = " << maxY << std::endl;

    for (int i = 0; i < nlines; ++i)
        Ysim[i] = 100. * Ysim[i] / maxY;

    TGraph* grSim = new TGraph(nlines, Tsim, Ysim);
    grSim->SetLineWidth(4);
    grSim->SetLineColor(kBlue);
    grSim->SetLineStyle(9);
    grSim->GetYaxis()->SetTitle(legendTitle.toLocal8Bit());
    grSim->GetXaxis()->SetTitle("Time (h)");
    grSim->GetYaxis()->SetRangeUser(0., 105.);
    grSim->GetXaxis()->SetRangeUser(0.0, double(timeScale));
    grSim->Draw();
    /////////////////////////////////////////////////////////////////
    auto* hexp = new TProfile("hexp", "selected YourEnzyme", 100, 0.0, 25, 0, 105);

    FILE* fexp = fopen(expPath.toLocal8Bit(), "r");
    if (!fexp) 
    {
        std::cerr << "Failed to open experimental data file!" << std::endl;
    }

    if (fexp)
    {
        float texp, yexp;
        int ncols = 0;

        while (true) {
            ncols = fscanf(fexp, " %f %f", &texp, &yexp);
            if (ncols < 0) break;
            hexp->Fill(texp, yexp);
        }

        fclose(fexp);

        // Настройки графика
        hexp->SetMarkerSize(3);
        hexp->SetMarkerColor(kRed);
        hexp->SetMarkerStyle(27);
        hexp->SetFillStyle(3005);

        // Важно: отрисовать на том же холсте
        hexp->Draw("P SAME");
    }

    auto* leg = new TLegend(0.87, 0.75, 0.60, 0.85);
    leg->SetTextSize(0.035);
    leg->SetFillColor(0);
    if (fexp)
        leg->AddEntry(hexp, "Experiment", "P");

    leg->AddEntry(grSim, "Simulation", "l");
    leg->Draw();

    c1->Update();

    QString plotPath = simPath.left(simPath.lastIndexOf(".")) + "plot.png";
    
    c1->SaveAs(plotPath.toLocal8Bit());
    emit rootGenerated(plotPath);
}

void LogicHandler::setArgs(int argc, char **argv)
{
    argc_main = argc;
    argv_main = argv;
}
