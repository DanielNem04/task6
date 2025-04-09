#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

void plot_histos_kt_averages() {
    const int numGroups = 10;
    const int numSubDirs = 4;

    const std::vector<std::string> subdirNames = {"1", "2", "3", "4"};
    const std::vector<std::string> subdirTitles = {
        "#it{k}_{T} = 0.0-0.2 GeV",
        "#it{k}_{T} = 0.2-0.4 GeV",
        "#it{k}_{T} = 0.4-0.6 GeV",
        "#it{k}_{T} = 0.6-0.8 GeV"
    };

    const std::vector<std::string> subdirShort = {"0_02", "02_04", "04_06", "06_08"};
    const std::vector<std::string> paramNames = {"N", "R", "alpha"};
    const std::vector<std::string> yAxisLabels = {"#LTN#GT", "R [fm]", "#alpha"};
    const std::string basePath = "histos/";

    std::vector<double> xValues = {351.4, 299, 253.9, 215.3, 181.6, 151.5, 125.7, 102.7, 82.9, 65.9};
    std::vector<double> xErrors(numGroups, 0.0);

    gSystem->Exec("mkdir -p plots_kt");

    std::vector<std::vector<TGraphErrors*>> graphs(3, std::vector<TGraphErrors*>(numGroups));

    for (int p = 0; p < 3; ++p) {
        for (int g = 0; g < numGroups; ++g) {
            graphs[p][g] = new TGraphErrors();
            graphs[p][g]->SetTitle(Form("%s vs. k_{T} (#LTNpart#GT = %.1f);k_{T} [GeV];%s",
                                        paramNames[p].c_str(), xValues[g], paramNames[p].c_str()));
            graphs[p][g]->SetMarkerStyle(20);
            graphs[p][g]->SetMarkerSize(1.2);
            graphs[p][g]->SetLineWidth(2);
        }
    }

    for (int groupIdx = 0; groupIdx < numGroups; ++groupIdx) {
        std::string groupPath = basePath + "group_" + std::to_string(groupIdx) + "/";

        for (int subIdx = 0; subIdx < numSubDirs; ++subIdx) {
            std::string filePath = groupPath + subdirNames[subIdx] + "/fit_results.txt";
            std::ifstream infile(filePath);
            if (!infile.is_open()) {
                std::cerr << "Nem olvasható: " << filePath << std::endl;
                continue;
            }

            std::string header;
            std::getline(infile, header);

            double Nsum = 0, Nerr2sum = 0;
            double Rsum = 0, Rerr2sum = 0;
            double Asum = 0, Aerr2sum = 0;
            int count = 0;

            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                std::string tag;
                double N, Nerr, R, Rerr, A, Aerr;

                if (!(iss >> tag >> N >> Nerr >> R >> Rerr >> A >> Aerr)) continue;

                Nsum += N;      Nerr2sum += Nerr * Nerr;
                Rsum += R;      Rerr2sum += Rerr * Rerr;
                Asum += A;      Aerr2sum += Aerr * Aerr;
                count++;
            }

            if (count == 0) continue;

            double Navg = Nsum / count;
            double Ravg = Rsum / count;
            double Aavg = Asum / count;

            double NerrAvg = std::sqrt(Nerr2sum) / count;
            double RerrAvg = std::sqrt(Rerr2sum) / count;
            double AerrAvg = std::sqrt(Aerr2sum) / count;

            double x = (subIdx + 0.5) * 0.2;

            graphs[0][groupIdx]->SetPoint(subIdx, x, Navg);
            graphs[0][groupIdx]->SetPointError(subIdx, 0, NerrAvg);

            graphs[1][groupIdx]->SetPoint(subIdx, x, Ravg);
            graphs[1][groupIdx]->SetPointError(subIdx, 0, RerrAvg);

            graphs[2][groupIdx]->SetPoint(subIdx, x, Aavg);
            graphs[2][groupIdx]->SetPointError(subIdx, 0, AerrAvg);
        }
    }

    for (int p = 0; p < 3; ++p) {
        for (int g = 0; g < numGroups; ++g) {
            TCanvas* c = new TCanvas(Form("kT_c%d_%d", p, g), graphs[p][g]->GetTitle(), 800, 600);
            graphs[p][g]->GetXaxis()->SetTitle("k_{T} [GeV]");
            graphs[p][g]->GetYaxis()->SetTitle(yAxisLabels[p].c_str());

            if (p == 2) {
                graphs[p][g]->GetYaxis()->SetRangeUser(1.0, 2.0);
            }

            graphs[p][g]->Draw("AP");

            std::string filename = Form("plots_kt/%s_vs_kT_Npart%.1f.png",
                                         paramNames[p].c_str(), xValues[g]);
            c->SaveAs(filename.c_str());
        }
    }
}
