#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <exception>
#include <chrono>
#include <cmath>
#include "levy_calc.h"

using namespace std;

double rmin = 1.0, rmax = 200.0;
int NDF = 0;

vector<double> levy_fit(const string& path) {
    auto start = chrono::high_resolution_clock::now();
    cout << "[TRACE] Entered levy_fit: " << path << endl;
    TFile f(path.c_str());
    if (!f.IsOpen() || f.IsZombie()) {
        cerr << "[ERROR] Cannot open ROOT file: " << path << endl;
        return {};
    }

    TH1* hist = nullptr;
    for (auto&& key : *f.GetListOfKeys()) {
        if (string(key->GetName()).find("rho_hist") != string::npos) {
            hist = dynamic_cast<TH1*>(f.Get(key->GetName()));
            break;
        }
    }

    if (!hist || hist->GetEntries() == 0) {
        cerr << "[WARNING] No valid 'rho_hist' in: " << path << endl;
        return {};
    }

    int N = hist->GetNbinsX();
    TGraphErrors gr;
    for (int i = 1; i <= N; ++i) {
        double y = hist->GetBinContent(i);
        double err = hist->GetBinError(i);
        if (y == 0 && err == 0) continue;
        gr.SetPoint(gr.GetN(), hist->GetBinCenter(i), y);
        gr.SetPointError(gr.GetN() - 1, 0, err);
    }

    if (gr.GetN() == 0) {
        cerr << "[WARNING] No valid data points in graph for: " << path << endl;
        return {};
    }

    auto fitfn = [](const double* x, const double* p) {
        if (p[1] <= 0 || p[2] <= 0) return 0.0;
        double Rcc = p[1] * pow(2., 1. / p[2]);
        if (Rcc <= 0 || !std::isfinite(Rcc)) return 0.0;

        if (x[0] > 1e4 || x[0]/Rcc > 1e4) {
            cerr << "[WARNING] Skipping too large x or x/Rcc: x = " << x[0] << ", Rcc = " << Rcc << endl;
            return 0.0;
        }

        double val = levy_calc(x[0] / Rcc, Rcc, p[2]);
        if (!std::isfinite(val)) {
            cerr << "[WARNING] Non-finite result from levy_calc: x = " << x[0]
                 << ", Rcc = " << Rcc << ", alpha = " << p[2] << endl;
            return 0.0;
        }
        return p[0] * val;
    };

    auto chi2 = [&](const double* p) {
        double sum = 0; NDF = 0;
        for (int i = 0; i < gr.GetN(); ++i) {
            double x, y;
            gr.GetPoint(i, x, y);
            if (x < rmin || x > rmax) continue;
            double err = gr.GetErrorY(i);
            if (err == 0) continue;
            double model = fitfn(&x, p);
            if (!std::isfinite(model)) continue;
            double diff = (y - model) / err;
            sum += diff * diff;
            NDF++;
        }
        return sum - 3;
    };

    ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kCombined);
    min.SetFunction(ROOT::Math::Functor(chi2, 3));
    min.SetMaxFunctionCalls(10000);  // korlátozottabb
    min.SetMaxIterations(1000);     // gyorsabb kilépés
    min.SetTolerance(1e-3);        // kevésbé szigorú

		min.SetLimitedVariable(0, "N",     1.0,  0.01,  -1,  50.0);
    min.SetLimitedVariable(1, "R",     8.0,  0.01,  -1,    100.0);
    min.SetLimitedVariable(2, "alpha", 1.2,  0.01,  -1,      5.0);


    cout << "[TRACE] Starting minimization..." << endl;
    min.Minimize();
    cout << "[TRACE] Minimization done. Status = " << min.Status() << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "[TIMER] Fit time: " << elapsed.count() << " seconds" << endl;

    if (min.Status() != 0) {
        cerr << "[ERROR] Fit failed for: " << path << endl;
        return {};
    }

    const double *p = min.X(), *e = min.Errors();
    double chi = chi2(p);

    if (NDF < 5 || chi > 1e6) {
        cerr << "[ERROR] Unreliable fit: NDF too small or chi2 too large in " << path << endl;
        return {};
    }

    TCanvas c("c", "", 800, 600);
    gStyle->SetOptStat(0);
    c.SetLogx(); c.SetLogy();
    gr.SetTitle((path + ";#rho [fm];D(#rho)").c_str());
    gr.SetMarkerStyle(20);
    gr.Draw("AP");

    TF1 fit("fit", fitfn, rmin, rmax, 3);
    fit.SetParameters(p);
    fit.SetLineColor(kRed);
    fit.Draw("same");

    TLatex latex;
    latex.SetNDC(); latex.SetTextSize(0.03);
    latex.DrawLatex(0.15, 0.85, Form("N = %.3f #pm %.3f", p[0], e[0]));
    latex.DrawLatex(0.15, 0.80, Form("R = %.3f #pm %.3f", p[1], e[1]));
    latex.DrawLatex(0.15, 0.75, Form("#alpha = %.3f #pm %.3f", p[2], e[2]));
    latex.DrawLatex(0.15, 0.70, Form("#chi^{2}/NDF = %.2f / %d", chi, NDF));

    string pngpath = path;
    pngpath.replace(pngpath.find(".root"), 5, ".png");
    c.SaveAs(pngpath.c_str());

    return {p[0], e[0], p[1], e[1], p[2], e[2]};
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <directory>" << endl;
        return 1;
    }

    for (auto& file : filesystem::recursive_directory_iterator(argv[1])) {
        if (file.path().extension() != ".root") continue;
        cout << "[INFO] Processing: " << file.path() << endl;
        try {
            auto res = levy_fit(file.path().string());
            if (!res.empty()) {
                string out = file.path().parent_path().string() + "/fit_results.txt";

                if (!filesystem::exists(out)) {
                    ofstream header(out);
                    header << "filename N pmN R pmR alpha pmAlpha\n";
                    header.close();
                }

                ofstream(out, ios::app) << file.path().filename().string() << " "
                                        << res[0] << " " << res[1] << " "
                                        << res[2] << " " << res[3] << " "
                                        << res[4] << " " << res[5] << "\n";
                cout << "[OK] Fit succeeded." << endl;
            } else {
                cout << "[WARN] Fit skipped or failed." << endl;
            }
        } catch (const exception& e) {
            cerr << "[EXCEPTION] While processing " << file.path() << ": " << e.what() << endl;
        }
    }
    return 0;
}
