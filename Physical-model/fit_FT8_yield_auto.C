#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kTol = 1e-8;

std::vector<std::string> discoverRateFiles(const std::string &suffix = "_rate.txt") {
    std::vector<std::string> files_out;
    TSystemDirectory dir(".", ".");
    if (TList *files = dir.GetListOfFiles()) {
        TIter next(files);
        while (TSystemFile *f = static_cast<TSystemFile *>(next())) {
            std::string name = f->GetName();
            if (!f->IsDirectory() && name.size() > suffix.size() &&
                name.rfind(suffix) == name.size() - suffix.size()) {
                files_out.push_back(name);
            }
        }
        delete files;
    }
    std::sort(files_out.begin(), files_out.end());
    return files_out;
}

double parseCouplingLabel(const std::string &label) {
    if (label == "sm") {
        return 0.0;
    }
    const std::string prefix = "FT8_";
    if (label.rfind(prefix, 0) != 0) {
        throw std::runtime_error("Unrecognised label in rate file: " + label);
    }
    std::string suffix = label.substr(prefix.size());
    double sign = 1.0;
    if (!suffix.empty() && suffix[0] == 'm') {
        sign = -1.0;
        suffix.erase(suffix.begin());
    }
    std::replace(suffix.begin(), suffix.end(), 'p', '.');
    if (suffix.empty()) {
        throw std::runtime_error("Empty coupling value for label " + label);
    }
    return sign * std::stod(suffix);
}

std::map<double, double> loadRates(const std::string &path) {
    std::ifstream src(path);
    if (!src.is_open()) {
        throw std::runtime_error("Failed to open rate file: " + path);
    }
    std::map<double, double> rates;
    std::string label;
    double value = 0.0;
    while (src >> label >> value) {
        rates[parseCouplingLabel(label)] = value;
    }
    if (rates.empty()) {
        throw std::runtime_error("No entries were read from " + path);
    }
    return rates;
}

bool approximatelyEq(double a, double b) {
    return std::abs(a - b) < kTol;
}

double rateForCoupling(const std::map<double, double> &rates, double coupling) {
    auto direct = rates.find(coupling);
    if (direct != rates.end()) {
        return direct->second;
    }
    for (const auto &kv : rates) {
        if (approximatelyEq(kv.first, coupling)) {
            return kv.second;
        }
    }
    std::ostringstream msg;
    msg << "Missing rate for coupling " << coupling;
    throw std::runtime_error(msg.str());
}

std::vector<double> makeYVals(const std::map<double, double> &rates,
                              const std::vector<double> &x_vals) {
    std::vector<double> y_vals;
    y_vals.reserve(x_vals.size());
    for (double coupling : x_vals) {
        y_vals.push_back(rateForCoupling(rates, coupling));
    }
    return y_vals;
}

bool containsCoupling(const std::vector<double> &targets, double value) {
    return std::any_of(targets.begin(), targets.end(),
                       [&](double v) { return approximatelyEq(v, value); });
}

int colorForIndex(size_t idx) {
    static const int palette[] = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
        kOrange + 1, kCyan + 2, kViolet + 6, kAzure + 4
    };
    return palette[idx % (sizeof(palette) / sizeof(int))];
}

struct YearGraphs {
    std::string label;
    TGraph *fitGraph = nullptr;
    TGraph *extraGraph = nullptr;
    TGraph *extrapolatedGraph = nullptr;
    TF1 *fitFunc = nullptr;
    int color = kBlack;
};

struct UserConfig {
    // Set to false to hide the |FT8|>8 points and skip extrapolation helper graph.
    bool drawExtrapolated = true;
    // Leave empty to include all discovered years, or list the years you want explicitly.
    std::vector<std::string> selectedYears = {
        "2016",
        "2016APV",
        "2017",
        "2018",
    };
};

UserConfig &config() {
    static UserConfig cfg;
    return cfg;
}

}  // namespace

void fit_FT8_yield_auto() {
    gStyle->SetOptFit(0);

    const auto &cfg = config();

    const auto rate_files = discoverRateFiles();
    if (rate_files.empty()) {
        Error("fit_FT8_yield_auto", "No *_rate.txt files were found in the current directory.");
        return;
    }

    const bool applyFilter = !cfg.selectedYears.empty();

    const std::vector<double> fit_couplings = {
        -4, -2, -1, -0.5, -0.3, -0.2, -0.1,
         0.0,
         0.1,  0.2,  0.3,  0.5,  1,   2,   4
    };
    // Only two extrapolated points: one at -16 and one at +16
    const std::vector<double> extrapolate_targets = { -16, 16 };
    const double displayMaxX = 16.0;  // Only draw up to |FT8|=16

    double globalXMin = std::numeric_limits<double>::max();
    double globalXMax = std::numeric_limits<double>::lowest();
    double globalYMin = std::numeric_limits<double>::max();
    double globalYMax = std::numeric_limits<double>::lowest();

    auto updateRange = [&](double x, double y) {
        if (x < globalXMin) globalXMin = x;
        if (x > globalXMax) globalXMax = x;
        if (y < globalYMin) globalYMin = y;
        if (y > globalYMax) globalYMax = y;
    };

    std::vector<YearGraphs> year_graphs;
    year_graphs.reserve(rate_files.size());

    size_t processedIndex = 0;
    for (const auto &file : rate_files) {
        const std::string label = file.substr(0, file.size() - std::string("_rate.txt").size());

        if (applyFilter) {
            if (std::find(cfg.selectedYears.begin(), cfg.selectedYears.end(), label) ==
                cfg.selectedYears.end()) {
                continue;
            }
        }

        std::map<double, double> rates;
        try {
            rates = loadRates(file);
        } catch (const std::exception &ex) {
            Error("fit_FT8_yield_auto", "Skipping %s: %s", file.c_str(), ex.what());
            continue;
        }

        std::vector<double> fit_x;
        std::vector<double> fit_y;
        std::vector<double> extra_x;
        std::vector<double> extra_y;

        for (const auto &kv : rates) {
            if (containsCoupling(fit_couplings, kv.first)) {
                fit_x.push_back(kv.first);
                fit_y.push_back(kv.second);
                updateRange(kv.first, kv.second);
            } else if (cfg.drawExtrapolated && std::abs(kv.first) <= displayMaxX) {
                extra_x.push_back(kv.first);
                extra_y.push_back(kv.second);
                updateRange(kv.first, kv.second);
            }
        }

        if (fit_x.empty()) {
            Warning("fit_FT8_yield_auto",
                    "No fit points found for %s, skipping.", label.c_str());
            continue;
        }

        const int color = colorForIndex(processedIndex++);

        const double fitMinX = *std::min_element(fit_x.begin(), fit_x.end());
        const double fitMaxX = *std::max_element(fit_x.begin(), fit_x.end());
        double funcMinX = -displayMaxX;
        double funcMaxX = displayMaxX;

        TGraph *gr_fit = new TGraph(static_cast<int>(fit_x.size()));
        for (int i = 0; i < gr_fit->GetN(); ++i) {
            gr_fit->SetPoint(i, fit_x[i], fit_y[i]);
        }
        gr_fit->SetMarkerStyle(20);
        gr_fit->SetMarkerColor(color);
        gr_fit->SetLineColor(0);
        gr_fit->SetLineStyle(0);

        TF1 *poly2 = new TF1(Form("poly2_%s", label.c_str()), "[0] + [1]*x + [2]*x*x",
                              funcMinX, funcMaxX);
        poly2->SetNpx(1000);
        gr_fit->Fit(poly2, "RQ");

        printf("==== %s ====\n", label.c_str());
        for (int p = 0; p < 3; ++p) {
            printf("  c%d = %.10f\n", p, poly2->GetParameter(p));
        }

        YearGraphs yg;
        yg.label = label;
        yg.fitGraph = gr_fit;
        yg.fitFunc = poly2;
        yg.color = color;

        if (cfg.drawExtrapolated && !extra_x.empty()) {
            TGraph *gr_extra = new TGraph(static_cast<int>(extra_x.size()));
            for (int i = 0; i < gr_extra->GetN(); ++i) {
                gr_extra->SetPoint(i, extra_x[i], extra_y[i]);
            }
            gr_extra->SetMarkerStyle(24);
            gr_extra->SetMarkerColor(color);
            gr_extra->SetLineColor(0);
            gr_extra->SetLineStyle(0);
            yg.extraGraph = gr_extra;
        }

        // Extrapolated points from fit function (positive and negative two points)
        if (cfg.drawExtrapolated) {
            TGraph *gr_ext = new TGraph(static_cast<int>(extrapolate_targets.size()));
            for (int i = 0; i < gr_ext->GetN(); ++i) {
                double x = extrapolate_targets[i];
                gr_ext->SetPoint(i, x, poly2->Eval(x));
                updateRange(x, poly2->Eval(x));
            }
            gr_ext->SetMarkerStyle(24);  // open circle
            gr_ext->SetMarkerColor(color);
            gr_ext->SetLineColor(color);
            gr_ext->SetLineStyle(2);
            yg.extrapolatedGraph = gr_ext;
        }

        year_graphs.push_back(yg);
    }

    if (year_graphs.empty()) {
        Error("fit_FT8_yield_auto", "No valid rate files could be processed.");
        return;
    }

    double xSpan = globalXMax - globalXMin;
    double ySpan = globalYMax - globalYMin;
    double xPad = xSpan > 0 ? 0.1 * xSpan : 1.0;
    double yPad = ySpan > 0 ? 0.1 * ySpan : 1.0;

    // Clamp x-range to displayMaxX
    double xMin = -displayMaxX;
    double xMax = displayMaxX;

    TCanvas *c = new TCanvas("c_fit_FT8_all", "FT8 Yield", 950, 650);
    TH2F *frame = new TH2F("frame_FT8_auto", "FT8 Yield",
                           1, xMin - xPad, xMax + xPad,
                           1, globalYMin - yPad, globalYMax + yPad);
    frame->SetStats(false);
    frame->SetTitle(";FT8;Yield");
    frame->GetXaxis()->SetTitle("FT8");
    frame->GetYaxis()->SetTitle("Yield");
    frame->Draw();

    for (const auto &yg : year_graphs) {
        yg.fitGraph->Draw("P SAME");
        if (yg.extraGraph) {
            yg.extraGraph->Draw("P SAME");
        }
        if (yg.extrapolatedGraph) {
            yg.extrapolatedGraph->Draw("P SAME");
        }
    }

    for (const auto &yg : year_graphs) {
        yg.fitFunc->SetLineColor(yg.color);
        yg.fitFunc->SetLineWidth(2);
        yg.fitFunc->Draw("SAME");
    }

    // Legend placed near top center (slightly lower)
    TLegend *leg = new TLegend(0.38, 0.66, 0.68, 0.84);
    leg->SetBorderSize(0);
    for (const auto &yg : year_graphs) {
        leg->AddEntry(yg.fitGraph, yg.label.c_str(), "p");
    }
    bool extra_shown = std::any_of(year_graphs.begin(), year_graphs.end(),
                                   [](const YearGraphs &yg) { return yg.extraGraph != nullptr; });
    bool extrap_shown = std::any_of(year_graphs.begin(), year_graphs.end(),
                                    [](const YearGraphs &yg) { return yg.extrapolatedGraph != nullptr; });
    if (extra_shown) {
        // Use the first available extra graph for legend entry
        for (const auto &yg : year_graphs) {
            if (yg.extraGraph) {
                leg->AddEntry(yg.extraGraph, "Not used for fit (|FT8| > 4)", "p");
                break;
            }
        }
    }
    // if (extrap_shown) {
    //     for (const auto &yg : year_graphs) {
    //         if (yg.extrapolatedGraph) {
    //             leg->AddEntry(yg.extrapolatedGraph, "Extrapolated (fit)", "p");
    //             break;
    //         }
    //     }
    // }
    std::string fitLabel = "Quadratic fit";
    if (!year_graphs.empty() && year_graphs.front().fitFunc) {
        auto f = year_graphs.front().fitFunc;
        fitLabel = Form("Quadratic fit: y = %.3g %+.3g x %+.3g x^{2}",
                        f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
    }
    leg->AddEntry(year_graphs.front().fitFunc, fitLabel.c_str(), "l");
    leg->Draw();

    c->SaveAs("fit_FT8_yield_auto.pdf");
    c->SaveAs("fit_FT8_yield_auto.png");
}
