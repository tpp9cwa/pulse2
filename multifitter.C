#include "TMath.h"
#include <array>
#include "PulseFitter.h"
#include "CamelFitter.h"
#include <TSpectrum.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <algorithm>

constexpr unsigned int N_Channels = 8;
constexpr unsigned int N_Samples = 1024;

static TH1D* hpro;
static int numberofpeaks = 0;

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


void smoothHistogram(TH1D* original, TH1D* smoothed, int windowSize) {
    int nbins = original->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        double sum = 0.0;
        int count = 0;
        for (int j = i - windowSize; j <= i + windowSize; ++j) {
            if (j > 0 && j <= nbins) {
                sum += original->GetBinContent(j);
                count++;
            }
        }
        smoothed->SetBinContent(i, sum / count);
    }
}

std::vector<double> find_peaks2(TH1D* hist1, TH1D* dummy=0) {

    std::vector<double> peaks;
    std::vector<double> troughs;

    TH1D* hist = (TH1D*)hist1->Clone("hist");
    hist->Rebin(2);
    //hist->Scale(0.5);
    //hist->Smooth(1);
    //smoothHistogram(hist1, hist, 50);
    //smoothHistogram(hist1, hist, 50);
    //hist->Draw("same");
    const double bins_per_mv = (hist->FindBin(100) - hist->FindBin(0)) / 100.0;
    //std::cout << "bins per mV: " << bins_per_mv << "\n";
    
    // Find extrema in a 1.5 mV window
    const int window = static_cast<int>(bins_per_mv * 30);
    //std::cout << "window: " << window << "\n";

    int i_prev = 0;
    bool find_max = true;
    for (int i = 0; i < hist->GetNbinsX() + 1 - window; i++) {
        // Find window extrema (skip empty bins in case of over binning)
        int i_window_ext = i;
        for (int j = 1; j < window; j++) {
            int index = i + j;
            if (hist->GetBinContent(index) == 0) continue;
            if (extremum(find_max, hist->GetBinContent(index), hist->GetBinContent(i_window_ext))) {
                i_window_ext = index;
            }
        }

        if (hist->GetBinContent(i_window_ext) == 0) {
            continue;
        }

	int backer = i -50;
	if(backer < 0){
	  backer = 0;
	}
	int here = hist->GetBinContent(i);
	int then = hist->GetBinContent(backer);
        // If derivative swaps sign, then we're at an extrema
        if (!extremum(find_max, hist->GetBinContent(i_window_ext), hist->GetBinContent(i_prev))) {

            if (find_max) {
                
                auto line = new TLine(hist->GetBinCenter(i_prev), 0, hist->GetBinCenter(i_prev), hist->GetBinContent(i_prev));
                
            }

            if (find_max) {if(hist->GetBinContent(i_prev)>10)peaks.push_back(hist->GetBinCenter(i_prev));}
            else troughs.push_back(hist->GetBinCenter(i_prev));

            // Search for the other extrema
            find_max = !find_max;
	    //if(!find_max && hist->GetBinContent(i_prev)>10)printf("Peak at %lf\n", hist->GetBinCenter(i_prev));
        }

        i_prev = i_window_ext;
    }
    if(dummy){
    *dummy = *hist;
    }
    return peaks;

}



void multifitter(TString filename, int channel, int entry=0) {

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
    std::vector<double> peaklist = find_peaks2(trialprof->ProjectionX("trialhist","e"));
    int maximaCount = peaklist.size();

    printf("Peaks: %d\n", maximaCount);

    // histogram to function!
    /*TF1 *hfn = new TF1("hfn",[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
		       hprof_px->GetBinLowEdge(1),
		       hprof_px->GetBinLowEdge(hprof_px->GetNbinsX()+1),3);
		       hfn->SetLineColor(kRed);*/

    int low = 1;
    int high = peaklist[0];
    int index = 0;
    int val = 0;
    for(int m = 0; m < peaklist.size(); m++){
      if(trialprof->GetBinContent(peaklist[m]) > trialprof->GetBinContent(peaklist[index])){
	val = trialprof->GetBinContent(peaklist[m]);
	index = m;
      }
    }

    if(index < maximaCount-1){
      high = peaklist[index+1]-15;
    }
    else{
      high = hprof_px->GetNbinsX()+1;
    }

    if(index > 1){
      low = peaklist[index-1] + 50;
    }
    else{
      low = 1;
    }
    
    TF1 *hfn = new TF1("hfn",[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
    hprof_px->GetBinLowEdge(low),
    hprof_px->GetBinLowEdge(high+1),3);
    hfn->SetLineColor(kRed);

    int up2 = trialprof->GetBinContent(peaklist[index]);
      if(up2 <= 5){up2 = 60;}
      hfn->SetParameter(0, trialprof->GetBinContent(peaklist[index]));
      hfn->SetParLimits(0, 5, up2);
      hfn->SetParameter(1, peaklist[index]-450);
      hfn->SetParLimits(1, peaklist[index]-450-50, peaklist[index]-450+(high-peaklist[index]));
      hfn->SetParameter(2, 12);
      hfn->SetParLimits(2, 0, 30);
    
    trialprof->Fit(hfn, "WR");
    
      /*std::vector<TF1*> funcs;
    int low = 1;
    int high = peaklist[0];
    for(int f = 0; f < maximaCount; f++){
      if(f < maximaCount - 1){
	high = peaklist[f+1]-25;
      }
      else{
	high = hprof_px->GetNbinsX();
      }
      TString base = "func";
      Int_t witch = f;
      TString yes = base + witch;
      TF1* func = new TF1(yes.Data(),[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
		       hprof_px->GetBinLowEdge(low),
		       hprof_px->GetBinLowEdge(1+high),3);
      low = peaklist[f];

      funcs.push_back(func);

      int up2 = trialprof->GetBinContent(peaklist[f])-3;
      if(up2 <= 5){up2 = 30;}
      funcs[f]->SetParameter(0, trialprof->GetBinContent(peaklist[f])-10);
      funcs[f]->SetParLimits(0, 5, up2);
      funcs[f]->SetParameter(1, peaklist[f]-450);
      funcs[f]->SetParLimits(1, peaklist[f]-450-50, peaklist[f]-450+(high-peaklist[f]));
      funcs[f]->SetParameter(2, 12);
      funcs[f]->SetParLimits(2, 0, 30);

      trialprof->Fit(funcs[f], "WR");

      funcs[f]->Draw("same");
      }

      */

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

double pulse(double *x, double *par, int num){

  return par[num*3+2]+par[num*3+0]*hpro->Interpolate(x[0]-par[num*3+1]);
  
}

Double_t multipeak(double *x, double *par){

  Double_t fout = 0;
  int npeaks = numberofpeaks;
  for(int i = 0; i < npeaks; i++){
    if(x[0] - par[1 + i*3] >= 444){
      fout += pulse(x, par, i);
    }
  }
  return fout;
}

void multifitter2lim(TString filename, int channel, int entry=0) {

    gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto canvas = new TCanvas("tcsum", "summary");
    canvas->Divide(2,1);
    canvas->cd(1);

    plot_samples(tree, channel, entry);

    //Creating average pulse histogram

    auto tf=new TFile("avepulse.root");
    TH1D* hprof_px = (TH1D*) ( tf->Get("hprof_px") );
    hprof_px->Fit("pol0","R0","",0,400);
    double offset = hprof_px->GetFunction("pol0")->GetParameter(0);
    for (int i=1; i<=hprof_px->GetNbinsX(); ++i) hprof_px->SetBinContent(i,hprof_px->GetBinContent(i)-offset);
    hprof_px->Scale(1.0/hprof_px->GetMaximum());
    hpro = hprof_px;

    //Filling profile with single buffer data

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

    //Peak finding; peaklist contains peak locations
    
    std::vector<double> peaklist = find_peaks2(trialprof->ProjectionX("trialhist","e"));
    int maximaCount = peaklist.size();

    printf("Peaks: %d\n", maximaCount);

    //Function to fit

    numberofpeaks = 2;

    int indexmax = 0;
    for(int i = 0; i < peaklist.size(); i++){
      if(peaklist[i]>peaklist[indexmax]){
	indexmax=i;
      }
    }
    std::vector<double> peaklist2;
    if(peaklist.size() > 2){
      if(indexmax == 0){
	peaklist2.push_back(peaklist[0]);
	peaklist2.push_back(peaklist[1]);
      }
      else if(indexmax == peaklist.size()-1){
	peaklist2.push_back(peaklist[peaklist.size()-2]);
	peaklist2.push_back(peaklist[peaklist.size()-1]);
      }
      else{
	peaklist2.push_back(peaklist[indexmax-1]);
	peaklist2.push_back(peaklist[indexmax]);
	peaklist2.push_back(peaklist[indexmax+1]);
	numberofpeaks = 3;
      }
    }
    TF1 *hfn = new TF1("hfn", multipeak, hprof_px->GetBinLowEdge(1),
		       hprof_px->GetBinLowEdge(hprof_px->GetNbinsX()+1), peaklist2.size()*3);
    hfn->SetParameter(0, peaklist2.size());
    cout << "Reduced Peak Number: " << peaklist2.size() << endl;
    cout << "Max Peak: " << peaklist[indexmax]-444 << endl;
    for(int i = 0; i < peaklist2.size() ; i++){
      hfn->SetParameter(3*i + 0, trialprof->GetBinContent(peaklist2[i]-444)-5);
      hfn->SetParameter(3*i + 1, peaklist2[i]-444);
      hfn->SetParameter(3*i + 2, 12);

      hfn->SetParLimits(3*i + 0, trialprof->GetBinContent(peaklist2[i]-444)-10, 50);
      hfn->SetParLimits(3*i + 1, peaklist2[i]-444-50, peaklist2[i]-444+50);
      hfn->SetParLimits(3*i + 2, 0, 20);
      TLine *line = new TLine(peaklist2[i]-444, 0, peaklist2[i]-444, 50);
      line->Draw("same");
      canvas->Update();
    }

    //Fit
    
    trialprof->Fit(hfn, "WR");
    
    //tf->Close();
}

void multifitter2(TString filename, int channel, int entry=0) {

    gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto canvas = new TCanvas("tcsum", "summary");
    canvas->Divide(2,1);
    canvas->cd(1);

    plot_samples(tree, channel, entry);

    //Creating average pulse histogram

    auto tf=new TFile("avepulse.root");
    TH1D* hprof_px = (TH1D*) ( tf->Get("hprof_px") );
    hprof_px->Fit("pol0","R0","",0,400);
    double offset = hprof_px->GetFunction("pol0")->GetParameter(0);
    for (int i=1; i<=hprof_px->GetNbinsX(); ++i) hprof_px->SetBinContent(i,hprof_px->GetBinContent(i)-offset);
    hprof_px->Scale(1.0/hprof_px->GetMaximum());
    hpro = hprof_px;

    //Filling profile with single buffer data

    TProfile* trialprof = new TProfile("trialprof", "Trial", N_Samples, -0.5*tree.horizontal_interval(), (static_cast<Double_t>(N_Samples)-0.5)*tree.horizontal_interval());
    
    tree.get_entry(entry);
    const auto times = tree.time();
    const auto volts = tree.voltages(channel);

    for (int s = 0; s < N_Samples; s++) {
      trialprof->Fill(times[s], volts[s]);
    }

    //Peak finding; peaklist contains peak locations
    TH1D* dummy = new TH1D();
    std::vector<double> peaklist = find_peaks2(trialprof->ProjectionX("trialhist","e"), dummy);
    int maximaCount = peaklist.size();

    printf("Peaks: %d\n", maximaCount);

    //Function to fit

    numberofpeaks = maximaCount;
    
    TF1 *hfn = new TF1("hfn", multipeak, hprof_px->GetBinLowEdge(1),
		       hprof_px->GetBinLowEdge(hprof_px->GetNbinsX()+1), peaklist.size()*3);

    hfn->SetNpx(500);

    hfn->SetParameter(0, peaklist.size());
    for(int i = 0; i < peaklist.size() ; i++){
      hfn->SetParameter(3*i + 0, trialprof->GetBinContent(peaklist[i]-444));
      hfn->SetParameter(3*i + 1, peaklist[i]-444);
      hfn->SetParameter(3*i + 2, 13);

      hfn->SetParLimits(3*i + 0, trialprof->GetBinContent(peaklist[i]-444)-10, 50);
      hfn->SetParLimits(3*i + 1, peaklist[i]-444-50, peaklist[i]-444+50);
      hfn->SetParLimits(3*i + 2, 0, 20);
      TLine *line = new TLine(peaklist[i]-444, 0, peaklist[i]-444, 50);
      line->Draw("same");
      canvas->Update();
    }

    //Fit
    
    trialprof->Fit(hfn, "WRQ");
    trialprof->Fit(hfn, "WRQ");
    trialprof->Fit(hfn, "WRQ");
    trialprof->Fit(hfn, "WRQ");
    trialprof->Fit(hfn, "WR");
    dummy->SetLineColor(kGreen);
    dummy->Draw("same");
    //tf->Close();
}

std::vector<double> getHeights(TString filename, int channel, int entry){
  gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto tf=new TFile("avepulse.root");
    TH1D* hprof_px = (TH1D*) ( tf->Get("hprof_px") );
    hprof_px->Fit("pol0","QR0","",0,400);
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
    std::vector<double> peaklist = find_peaks2(trialprof->ProjectionX("trialhist","e"));
    int maximaCount = peaklist.size();

    // histogram to function!
    //TF1 *hfn = new TF1("hfn",[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
    //hprof_px->GetBinLowEdge(1),
    //hprof_px->GetBinLowEdge(hprof_px->GetNbinsX()+1),3);
    //hfn->SetLineColor(kRed);

    std::vector<double> heights;
    std::vector<TF1*> funcs;
    int low = 1;
    int high = peaklist[0];    
    
    for(int f = 0; f < maximaCount; f++){
      if(f < maximaCount - 1){
	high = peaklist[f+1]-25;
      }
      else{
	high = hprof_px->GetNbinsX();
      }
      TString base = "func";
      Int_t witch = f;
      TString yes = base + witch;
      TF1* func = new TF1(yes.Data(),[=](double *x, double *par){return par[2]+par[0]*hprof_px->Interpolate(x[0]-par[1]);},
		       hprof_px->GetBinLowEdge(low),
		       hprof_px->GetBinLowEdge(1+high),3);
      low = peaklist[f];

      funcs.push_back(func);

      int up2 = trialprof->GetBinContent(peaklist[f])-3;
      if(up2 <= 5){up2 = 30;}
      funcs[f]->SetParameter(0, trialprof->GetBinContent(peaklist[f])-10);
      funcs[f]->SetParLimits(0, 5, up2);
      funcs[f]->SetParameter(1, peaklist[f]-450);
      funcs[f]->SetParLimits(1, peaklist[f]-450-50, peaklist[f]-450+(high-peaklist[f]));
      funcs[f]->SetParameter(2, 12);
      funcs[f]->SetParLimits(2, 0, 30);

      trialprof->Fit(funcs[f], "QWR0");

      heights.push_back(funcs[f]->GetParameter(0));

      delete func;
      func = nullptr;

      }

    delete trialprof;

    return heights;
    
}

void pulsedist(TString filename, int channel){

  gStyle->SetOptStat(0);

    TreeReader tree(filename);

    auto canvas = new TCanvas("tcsum", "summary");
    canvas->Divide(1,1);
    canvas->cd(1);

    TH1D* dist = new TH1D("dist", "Pulse Heights;Height;Frequency",1001, 0, 101);
    
    for(int i = 0; i < 500;i++){
      std::vector<double> heights = getHeights(filename, channel, i);
      for(int h = 0; h < heights.size(); h++){
	dist->Fill(heights[h]);
	system("clear");
	printf("%lf%%\n", 100*(double)i/tree.num_entries());
	//printf("entry: %d, height: %lf\n",i,heights[h]);
      }
    }

    dist->Draw();
  
}
