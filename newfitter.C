#include "TMath.h"
#include <array>
#include "PulseFitter.h"
#include "CamelFitter.h"
#include <TSpectrum.h>
#include <iostream>
#include <vector>
#include <algorithm>

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

void plot_samples(TreeReader& tree, int channel, int entry=0) {

    tree.get_entry(entry);
    auto g = new TGraph(N_Samples, tree.time().data(), tree.voltages(channel).data());
    g->SetTitle(TString::Format("Sample buffer channel %d;Time [ns];Voltage [mV]", channel));
    g->Draw("ac");
}

int CountLocalMaxima(TProfile *profile, int windowSize) {
    int count = 0;
    int nbins = profile->GetNbinsX();

    // Smoothed data vector
    std::vector<double> smoothedData(nbins);

    // Apply moving average smoothing to the data
    for (int i = 1; i <= nbins; ++i) {
        double sum = 0.0;
        int start = std::max(1, i - windowSize / 2);
        int end = std::min(nbins, i + windowSize / 2);
        for (int j = start; j <= end; ++j) {
            sum += profile->GetBinContent(j);
        }
        smoothedData[i - 1] = sum / (end - start + 1);
    }

    // Find local maxima
    for (int i = 1; i <= nbins - 2; ++i) {
        if (smoothedData[i] > smoothedData[i - 1] && smoothedData[i] > smoothedData[i + 1]) {
            ++count;
        }
    }

    return count;
}


void newfitter(TString filename, int channel, int entry=0) {

    gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto canvas = new TCanvas("tcsum", "summary");
    canvas->Divide(2,1);
    canvas->cd(1);

    plot_samples(tree, channel, entry);

    auto tf=new TFile("avepulse.root");
    TH1D* hprof_px = (TH1D*) ( tf->Get("hprof_px") );
    hprof_px->Fit("pol0","R0","",0,400);
    double offset = hprof_px->GetFunction("pol0")->GetParameter(0);
    for (int i=1; i<=hprof_px->GetNbinsX(); ++i) hprof_px->SetBinContent(i,hprof_px->GetBinContent(i)-offset);
    hprof_px->Scale(1.0/hprof_px->GetMaximum());



    TProfile* trialprof = new TProfile("trialprof", "Trial", N_Samples, -0.5*tree.horizontal_interval(), (static_cast<Double_t>(N_Samples)-0.5)*tree.horizontal_interval());
    
    tree.get_entry(entry);
    const auto times = tree.time();
    const auto volts = tree.voltages(channel);

    auto vmax = volts[0];
    auto tmax = times[0];
    int imax = 0;
    for (int s = 0; s < N_Samples; s++) {
      trialprof->Fill(times[s], volts[s]);
      if(volts[s]>vmax){
	vmax=volts[s];
	tmax=times[s];
	imax=s;
      }
    }

    int windowSize = 40;
    int maximaCount = CountLocalMaxima(trialprof, windowSize);
    int minCount = maximaCount;

    for(int w = 1; w < 75 ; w++){
      maximaCount = CountLocalMaxima(trialprof, w);
      if(maximaCount < minCount){
	minCount = maximaCount;
	windowSize = w;
      }
    }
    printf("Peaks: %d\t Window: %d\n", minCount, windowSize);
    
    // histogram to function!
    TF1 *hfn = new TF1("hfn",[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
		       hprof_px->GetBinLowEdge(1),
		       hprof_px->GetBinLowEdge(hprof_px->GetNbinsX()+1),3);
    hfn->SetLineColor(kRed);

    hfn->SetParameter(0,vmax-12);
    hfn->SetParLimits(0, vmax/3, 28);
    hfn->SetParameter(1,tmax-450);
    hfn->SetParameter(2,7);
    //hfn->SetParLimits(2, 7, vmax/4);

    printf("%d %d\n",vmax-12,tmax-450);

    trialprof->Fit(hfn, "WR", "", 450, 650);
    //trialprof->Fit(hfn, "W");

    hfn->Draw("same");
    

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
    //tf->Close();
}
