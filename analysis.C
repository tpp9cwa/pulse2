#include "TMath.h"
#include <array>
#include "PulseFitter.h"
#include "CamelFitter.h"

constexpr unsigned int N_Channels = 8;
constexpr unsigned int N_Samples = 1024;

int binCalc(double x_min, double x_max, double vgain) {
    int n = (int) ((x_max - x_min) / (1000 * vgain / 0x1000));
    std::cout << "n bins: " << n << "\n";
    return n;
}

class TreeReader {
    public:
        TreeReader(TString filename) {
            file_ = new TFile(filename);
            tree_ = (TTree*) file_->Get("tree");

            tree_->SetBranchAddress("vertical_gain", vertical_gain_.data());
            tree_->SetBranchAddress("vertical_offset", vertical_offset_.data());
            tree_->SetBranchAddress("horizontal_interval", &horizontal_interval_);
            tree_->SetBranchAddress("channels", channels_[0].data());
        }

        UInt_t num_entries() const {
            return tree_->GetEntries();
        }

        void get_entry(UInt_t i) {
            tree_->GetEntry(i);
        }

        std::array<Double_t, N_Samples> time() const {
            std::array<Double_t, N_Samples> t;
            for (int i = 0; i < N_Samples; i++) {
                t[i] = static_cast<Double_t>(i) * horizontal_interval_;
            }
            return t;
        }

        std::array<Double_t, N_Samples> voltages(UInt_t channel) const {
            std::array<Double_t, N_Samples> volts;
            const Double_t gain = vertical_gain_[channel];
            const Double_t offset = vertical_offset_[channel];
            for (int i = 0; i < N_Samples; i++) {
                volts[i] = gain * (static_cast<Double_t>(channels_[channel][i]) - offset);
            }
            return volts;
        }

        Float_t horizontal_interval() const { return horizontal_interval_; }
        Float_t vertical_gain(UInt_t channel) const { return vertical_gain_[channel]; }
        Float_t vertical_offset(UInt_t channel) const { return vertical_offset_[channel]; }

        UInt_t get_bin(Double_t t) const {
            return static_cast<UInt_t>(t / horizontal_interval_);
        }

    private:
        TFile* file_;
        TTree* tree_;

        Float_t horizontal_interval_;
        std::array<Float_t, N_Channels> vertical_offset_, vertical_gain_;
        std::array<std::array<Float_t, N_Samples>, N_Channels> channels_;
};

void plot_samples(TreeReader& tree, int channel) {

    std::cout << "Plotting samples\n";

    constexpr int num_to_plot = 3;
    auto mg = new TMultiGraph();
    mg->SetTitle(TString::Format("Sample buffer channel %d;Time [ns];Voltage [mV]", channel));
    for (int i = 0; i < num_to_plot; i++) {
        tree.get_entry(i);
        auto g = new TGraph(N_Samples, tree.time().data(), tree.voltages(channel).data());
        g->SetLineColor(kGray + i);
        mg->Add(g);
    }
    mg->Draw("ac");
}


void analysis(TString filename, int channel) {

    gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto canvas = new TCanvas("tcsum", "summary");
    canvas->Divide(2,2);
    canvas->cd(1);

    plot_samples(tree, channel);


    canvas->cd(2);

    // Average histogram plots

    auto hprof = new TProfile("hprof", "Average waveform;Time [ns];Voltage [mV]", N_Samples,
            -0.5 * tree.horizontal_interval(), 
            (static_cast<Double_t>(N_Samples) - 0.5) * tree.horizontal_interval());

    for (int i = 0; i < tree.num_entries(); i++) {
        tree.get_entry(i);
        const auto times = tree.time();
        const auto volts = tree.voltages(channel);
        for (int s = 0; s < N_Samples; s++) {
            hprof->Fill(times[s], volts[s]);
        }
    }
    hprof->Draw();

    
    PulseFitter pulse_fit(hprof->ProjectionX());
    pulse_fit.pulse_fit->Draw("same");

    //auto tfout=new TFile("avepulse.root","recreate");
    //auto h=hprof->ProjectionX();
    //h->Write();
    //tfout->Write();
    //tfout->Close();
    
    // Integration

    const PulseParams* params = pulse_fit.parameters();


    const double baseline = params->pedestal;

    const double start = params->x0;
    const int i_start = hprof->GetBin(start);

    const double peak = pulse_fit.x_peak;
    const int i_peak = pulse_fit.i_peak;

    const double stop = pulse_fit.x_peak + 60;
    const int i_stop = hprof->GetBin(stop);

    // try a shorter integration range around peak of profile hist
    // i_start = i_peak - 40;
    // i_stop = i_peak + 20;

    std::cout << "Sample peak at: " << peak << "\n"
        << "Integrating over range: " << start << " : " << stop << "\n"
        << "Integrating over bins: " << i_start << " : " << i_stop << "\n";

    const double y_min = 0;
    const double y_max = pulse_fit.y_peak * 1.1;

    for (double x : {start, stop, peak}) {
        auto line = new TLine(x, y_min, x, y_max);
        line->SetLineStyle(2);
        if (x == peak) line->SetLineColor(kRed);
        line->Draw();
    }

    double i_scale = pulse_fit.y_peak / hprof->Integral(i_start, i_stop);

    const double x_min = -baseline;
    const double x_max = y_max * 3;
    int nx = binCalc(x_min, x_max, tree.vertical_gain(channel));
    auto phd = new TH1D("phd", "Pulse Heights;Voltage [mV];Count", nx, x_min, x_max);
    auto pid = new TH1D("pid", "Pulse Integral around peak / #Deltat;Voltage [mV];Count", nx, x_min, x_max);
    
    
    for (int i = 0; i < tree.num_entries(); i++) {
        tree.get_entry(i);
        const auto times = tree.time();
        const auto volts = tree.voltages(channel);

        phd->Fill(volts[i_peak] - baseline);

        double sum = 0;
        for (int n = i_start; n < i_stop; n++) {
            sum += volts[n] - baseline;
        }
        pid->Fill(sum * i_scale);
    }
    

    canvas->cd(3);
    phd->SetLineColor(kBlack);
    phd->Draw();
    //find_peaks(phd);

    canvas->cd(4);
    pid->SetLineColor(kBlack);
    pid->Draw();
    //find_peaks(pid);

    canvas->Update();
    canvas->SaveAs(TString::Format("peaks_%d.png", channel));
}
