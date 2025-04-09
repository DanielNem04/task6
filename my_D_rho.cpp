#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <utility>
#include <regex>
#include <filesystem>
#include <sys/stat.h>

const int npMax = 84000;

struct FourVector {
    float t, x, y, z;
};

struct Particle {
    float m, px, py, pz, x, y, z, t;
    int id, ist;
};

FourVector compute_K(const Particle& p1, const Particle& p2) {
    return {0.5f * ((std::sqrt(p1.m * p1.m + p1.px * p1.px + p1.py * p1.py + p1.pz * p1.pz) +
                    std::sqrt(p2.m * p2.m + p2.px * p2.px + p2.py * p2.py + p2.pz * p2.pz))),
            0.5f * (p1.px + p2.px),
            0.5f * (p1.py + p2.py),
            0.5f * (p1.pz + p2.pz)};
}

std::pair<float, int> compute_rho_LCMS_kT(const Particle& p1, const Particle& p2, int index) {
    FourVector K = compute_K(p1, p2);
    float r_x = p1.x - p2.x;
    float r_y = p1.y - p2.y;
    float r_z = p1.z - p2.z;
    float t = p1.t - p2.t;

    float K_perp = std::sqrt(K.x * K.x + K.y * K.y);
    float K_long2 = K.t * K.t - K.z * K.z;
    if (K_perp < 1e-5 || K_long2 < 1e-5) return std::make_pair(-1.0f, -1);
    float K_long = std::sqrt(K_long2);

    float rho_out = (r_x * K.x / K_perp) + (r_y * K.y / K_perp) - (K_perp / K_long) * (K.t * t - K.z * r_z);
    float rho_side = (-r_x * K.y / K_perp) + (r_y * K.x / K_perp);
    float rho_long = (K.t * r_z - K.z * t) / K_long;

    float rho = std::sqrt(rho_out * rho_out + rho_side * rho_side + rho_long * rho_long);
    int ikT = floor(K_perp / 0.2f);

    return std::make_pair(rho, ikT);
}

std::string getGroupOutputPath(const std::string& filename, int index_value) {
    std::regex pattern("z-gg2my59a-(\\d+)\\.root");
    std::smatch match;

    if (std::regex_search(filename, match, pattern)) {
        int index = std::stoi(match[1]);
        int group = index % 10;

        std::string baseName = match[0].str().substr(0, match[0].str().find_last_of("."));
        std::string dir = "histos/group_" + std::to_string(group) + "/" + std::to_string(index_value);

        std::filesystem::create_directories(dir);

        return dir + "/" + baseName + "-DrhoHist.root";
    } else {
        throw std::runtime_error("Hibás fájlnév: " + filename);
    }
}

void histomaker(int index_value, std::vector<float>& v, int nE, const std::string& input_filename) {
    if (v.empty()) {
        std::cout << "No data for index " << index_value << ", skipping histogram." << std::endl;
        return;
    }

    const int nbins = 60;
    float rho_min = 0.8f;
    float rho_max = 1000.0f;
    float binfactor = TMath::Power(rho_max / rho_min, 1.0f / nbins);
    float rho_bins[nbins + 1];
    rho_bins[0] = rho_min;
    for (int ibin = 1; ibin <= nbins; ibin++)
        rho_bins[ibin] = rho_min * TMath::Power(binfactor, ibin);

    std::string base_name = std::filesystem::path(input_filename).stem().string();
    std::string hname = "rho_hist_" + std::to_string(index_value) + "_" + base_name;

    TH1D* h_rho = new TH1D(hname.c_str(), "D(#rho); #rho (fm); D(#rho)", nbins, rho_bins);
    h_rho->Sumw2();

    for (float rho : v) {
        if (rho > 0) {
            float weight = 1.0f / (4 * TMath::Pi() * rho * rho) / ( v.size() );
            h_rho->Fill(rho, weight);
        }
    }

    for (int i = 1; i <= h_rho->GetNbinsX(); i++) {
        float bin_content = h_rho->GetBinContent(i);
        float bin_width = h_rho->GetBinWidth(i);
        h_rho->SetBinContent(i, bin_content / bin_width);
        h_rho->SetBinError(i, h_rho->GetBinError(i) / bin_width);
    }

    std::string output_filename = getGroupOutputPath(input_filename, index_value);

    TFile* outfile = new TFile(output_filename.c_str(), "RECREATE");
    h_rho->Write();
    outfile->Close();
    delete outfile;
    delete h_rho;

    std::cout << "Histogram saved: " << output_filename << std::endl;
}

int process_single_file(const std::string& input_filename) {
    TFile* inputfile = new TFile(input_filename.c_str(), "READ");
    if (!inputfile || inputfile->IsZombie()) {
        std::cerr << "Cannot open ROOT file: " << input_filename << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)inputfile->Get("teposevent");
    if (!tree) {
        std::cerr << "Could not find tree 'teposevent' in file: " << input_filename << std::endl;
        return 1;
    }

    int nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        std::cerr << "No entries found in tree in file: " << input_filename << std::endl;
        return 1;
    }

    int np;
    std::vector<int> id(npMax), ist(npMax);
    std::vector<float> mass(npMax), px(npMax), py(npMax), pz(npMax), x(npMax), y(npMax), z(npMax), t(npMax);
    std::vector<float> kT_02, kT_04, kT_06, kT_08, kT_10;
    kT_02.reserve(npMax);
    kT_04.reserve(npMax);
    kT_06.reserve(npMax);
    kT_08.reserve(npMax);
    kT_10.reserve(npMax);

    tree->SetBranchAddress("np", &np);
    tree->SetBranchAddress("id", id.data());
    tree->SetBranchAddress("ist", ist.data());
    tree->SetBranchAddress("e", mass.data());
    tree->SetBranchAddress("px", px.data());
    tree->SetBranchAddress("py", py.data());
    tree->SetBranchAddress("pz", pz.data());
    tree->SetBranchAddress("x", x.data());
    tree->SetBranchAddress("y", y.data());
    tree->SetBranchAddress("z", z.data());
    tree->SetBranchAddress("t", t.data());

    for (int iEvent = 0; iEvent < nEntries; iEvent++) {
        if (tree->GetEntry(iEvent) <= 0) continue;

        std::vector<Particle> selectedParticles;
        for (int i = 0; i < np; i++) {
            if (id[i] == 120 && ist[i] == 0) {
                float pT = std::sqrt(px[i]*px[i] + py[i]*py[i]);
                float p = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
                float eta = 0.5f * std::log((p + pz[i]) / (p - pz[i]));

                if (pT > 0.15f && pT < 1.0f && std::abs(eta) < 1.0f) {
                    selectedParticles.push_back({mass[i], px[i], py[i], pz[i], x[i], y[i], z[i], t[i], id[i], ist[i]});
                }
            }
        }

        size_t nSelected = selectedParticles.size();
        for (size_t i = 0; i < nSelected; i++) {
            for (size_t j = i + 1; j < nSelected; j++) {
                auto [rho, kTindex] = compute_rho_LCMS_kT(selectedParticles[i], selectedParticles[j], 1);
                switch (kTindex) {
                    case 1: kT_02.push_back(rho); break;
                    case 2: kT_04.push_back(rho); break;
                    case 3: kT_06.push_back(rho); break;
                    case 4: kT_08.push_back(rho); break;
                    case 5: kT_10.push_back(rho); break;
                }
            }
        }
    }

    histomaker(1, kT_02, nEntries, input_filename);
    histomaker(2, kT_04, nEntries, input_filename);
    histomaker(3, kT_06, nEntries, input_filename);
    histomaker(4, kT_08, nEntries, input_filename);
    histomaker(5, kT_10, nEntries, input_filename);

    inputfile->Close();
    delete inputfile;
    return 0;
}

void process_all_root_files_in_directory(const std::string& directory_path) {
    for (const auto& entry : std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == ".root") {
            std::string file_path = entry.path().string();
            std::cout << "\n--- Processing file: " << file_path << " ---" << std::endl;
            process_single_file(file_path);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.root> OR <directory_with_root_files>" << std::endl;
        return 1;
    }

    std::string input = argv[1];
    if (std::filesystem::is_directory(input)) {
        process_all_root_files_in_directory(input);
    } else {
        process_single_file(input);
    }

    return 0;
}
