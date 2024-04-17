// 

#include "TMath.h"

int binCalc(double &xmin, double &xmax, double vgain){
  return (int) ((xmax-xmin) / fabs(vgain)); // roughly # of ADC counts for range
}

double pulse(double *t, double *par){
  double x=t[0];
  double x0=par[0];
  double A=par[1];
  double a=par[2];
  double tau=par[3];
  double ped=par[4];
  double pulse = ped;
  if (x>x0) pulse += A*pow((x-x0),a)*exp(-(x-x0)/tau);
  return pulse;
}



// to do: improve init parameter evaluation
class pulseFitter {
public:
  TF1 *fn;
  pulseFitter(TH1* h, bool test=false){
    double xmin = h->GetBinCenter(1);
    double xmax = h->GetBinCenter(h->GetNbinsX());
    /*
    // histogram to function!
    TF1 *hfn = new TF1("hfn",[=](double *x, double *){return h->Interpolate(x[0]);},
		       xmin,xmax,0);
    TGraph *g = new TGraph(hfn,"d");  // derivative of histogram "function"
    auto y = g->GetY();
    auto x = g->GetX();
    int n  = g->GetN();
    int i  = std::max_element(y,y+n)-y;
    double x0=x[i];  // max value of derivative
    */

    fn = new TF1("pulse",pulse,xmin,xmax,5);
    fn->SetLineColor(kRed-4);

    int imax = h->GetMaximumBin();
    double xM = h->GetBinCenter(imax);
    double ymax = h->GetBinContent(imax);
    int x0 = xM;
    for (int i=imax; i>0 ; --i){
      if (h->GetBinContent(i)<ymax*0.1){
	x0 = h->GetBinCenter(i);
	break;
      }
    }
    fn->SetParameter(0, x0);

    double baseline=h->GetBinContent(1);
    fn->FixParameter(4, baseline);

    fn->SetParameter(1,(ymax-baseline)/(xM-x0)); 
    double a = 1;
    fn->SetParameter(2, a);
    double tau = (xM-x0)/a/log(xM-x0);
    fn->SetParameter(3, tau);
    
    fn->SetRange(xmin,xM+100);
    if (test) {
      cout << "A " <<  ymax-baseline << endl;
      cout << "x0 " << x0 << endl; 
      cout << "xM " << xM << endl;
      cout << "a " << a << endl;
      cout << "tau " << tau << endl;
      cout << "ped " << baseline << endl;
      return;
    }
    h->Fit("pulse","N");
    //fn->SetRange(0,95); // (x0+xM)/2);
    h->Fit("pulse","0");  // final fit only includes the main pulse area after the peak
  }
};


void pulse_height(TString file="C2--CV-40_54V-partslaser--5k--00000.root",int ch=3,
		  bool savePlot=false){

  const double dTbaseline = 3;  // time before x0 to start integral in ns
  const double nTau =0;         // end of integration region in units of tail decay parameter

  auto tf=new TFile(file);
  auto tree=(TTree*)tf->Get("tree");


  TString chan = TString::Format("ch%d",ch);
  // output file
  TString fnout=file;
  fnout.ReplaceAll(".root","_"+chan+".pdf");
  int pos=fnout.Last('/');
  if ( pos>-1 ) fnout.Remove(0,pos+1);

  gStyle->SetOptStat(0);

  // find the size of the sample buffer
  Int_t samples; 
  Float_t horiz_interval;
  Float_t vertical_offset[8];
  Float_t vgains[8];

  TString m; 

  tree->SetBranchAddress("vertical_gain", vgains); 
  tree->SetBranchAddress("samples", &samples);
  tree->SetBranchAddress("horizontal_interval", &horiz_interval);
  tree->SetBranchAddress("vertical_offset", vertical_offset);

  tree->GetEntry(0);

  Float_t vgain = vgains[ch];
    
  cout << "Sample count: " << samples << endl; 
  cout << "Horizontal interval: " << horiz_interval << endl; 
  cout << "Ch:" << ch << " vertical gain:" << vgain << endl; 
  int LEN = samples; 
  cout << "Processing buffers of length: " << LEN << endl;
  cout << "Number of buffers: " << tree->GetEntries() << endl;

  auto tcsum = new TCanvas("tcsum","summary");
  tcsum->Divide(2,2);
  tcsum->cd(1);

  // Draw a sample buffer
  m.Form("(channels[%d]-vertical_offset[%d])*1000*vertical_gain[%d]:time*1e9 >> htrace", ch, ch, ch);
  //tree->Draw(m,"event==1","GOFF",1);
  //std::cout << m << std::endl;
  int nplot=3;
  int n;
  auto mg = new TMultiGraph();
  mg->SetTitle(TString::Format("Sample buffer ch%d;ns;mV",ch));
  for (int i=0; i<nplot; ++i){
    n = tree->Draw(m,"","GOFF",1,i);
    TGraph *g = new TGraph(n,tree->GetV2(),tree->GetV1()); 
    g->SetLineColor(kGray+i);
    mg->Add(g);
  }
  mg->Draw("ac"); 

  // find polarity
  //tree->Draw(m,"event<10","",10);
  tree->Draw(m,"event<10","GOFF",10);
  TH2F *htrace = (TH2F*)gDirectory->Get("htrace");
  auto projY =  htrace->ProjectionY();
  double min = projY->GetBinCenter(1);
  double max = projY->GetBinCenter(projY->GetNbinsX());
  cout << "min,max " << min << " " << max << endl;
  //htrace->ProjectionY()->Draw();

  double polarity=1.0;
  if (fabs(min)>fabs(max)){  // careful this might not be smart enough
    polarity = -1.0;  
  }
  cout << "signal polarity is: " << polarity << endl;

  // convert calibration to mV and ns
  vgain*=1000*polarity;
  double dt = horiz_interval * 1e9;


  int16_t *pool = new int16_t[LEN*8]; 
  Double_t *time = new Double_t[LEN]; 
  // int16_t **volts = new int16_t *[8];

  int16_t volts[8][1024];

  // for (int i = 0; i < 8; i++) {
  //   volts[i] = &(pool[LEN*i]);
  // }

  Double_t startx;
  Long_t event;
  
  tree->SetBranchAddress("time",time);
  // tree->SetBranchAddress("channels",pool);
  tree->SetBranchAddress("channels", volts);

  // determine the average pulse shape
  auto hprof= new TProfile("hprof","Average waveform;time [ns];mV",LEN,-dt/2,
			   LEN*dt-dt/2);
  double vmin=10000, vmax=-10000;
  for (int i=0; i<tree->GetEntries(); ++i){
    tree->GetEntry(i);
    for (int n=0; n<LEN; ++n) {
      double v = (volts[ch][n]-vertical_offset[ch])*vgain; 
      vmin = std::min(vmin, v);
      vmax = std::max(vmax, v);
      hprof->Fill(n*dt,v);   
    }
  }

  min = hprof->GetMinimum();
  max = hprof->GetMaximum();
  cout << "min/max samples: " << vmin << "/" << vmax << " mV" << endl;

  tcsum->cd(2);

  // get a rough fit to the pulse shape
  pulseFitter fP(hprof->ProjectionX(),true);
  hprof->Draw();
  fP.fn->Draw("same");
  tcsum->Update();
 

  // set up integration window and baseline
  double start = fP.fn->GetParameter(0)-dTbaseline;
  
  double baseline = fP.fn->GetParameter(4);
  int iBLS = hprof->FindBin(start);  // use 1 sample before rise as baseline
  cout << "sample baseline at: " << hprof->GetBinCenter(iBLS) << endl;

  // note: there is a bias in the peak loaction, probably due to averaging over pileup
  // the average pulse shape is less peaked and the maximum value is later than for single pulses
  int ipeak = hprof->GetMaximumBin();
  double peak=hprof->GetBinCenter(ipeak);

  int istart = iBLS;
  double stop = peak + nTau * fP.fn->GetParameter(3);
  int istop = hprof->FindBin(stop);
  // try a shorter integration range around peak of profile hist
  istart = ipeak-40;
  //istart = hprof->FindBin(92);
  istop = ipeak+20;

  cout << "Sample peak at: " << peak << endl;
  cout << "Integrating over range: " << start << " : " << stop << endl;
  cout << "Integrating over bins: " << istart << " : " << istop << endl;

  double ymin=0;
  double ymax=hprof->GetMaximum()*1.1;
  auto l1= new TLine(start,ymin,start,ymax);
  l1->SetLineStyle(2);
  l1->Draw();
  auto l2= new TLine (stop,ymin,stop,ymax);
  l2->SetLineStyle(2);
  l2->Draw();
  auto l3= new TLine (peak,ymin,peak,ymax);
  l3->SetLineStyle(2);
  l3->SetLineColor(kRed);
  l3->Draw();
  
  
  // find scale factor for the integral to normalize it to the pHd
  double iScale = hprof->GetMaximum() / hprof->Integral(istart,istop);


  TString plotfile=file;
  int loc=plotfile.Last('/');
  if ( loc>-1 ) plotfile.Remove(0,loc+1);
  TString outfile=plotfile;
  outfile.ReplaceAll(".root","_out.root");
  auto tfout = new TFile(outfile,"RECREATE");

  //auto phd = new TH1F("phd","PulseHeights;mV;frequency",160,-0.5,ymax*1000*2);
  double xmin=-20.0;
  double xmax=ymax*3;  // hack seems to work
  int nx = binCalc(xmin,xmax,vgain);
  auto phd = new TH1F("phd","PulseHeights;mV;frequency",nx,xmin,xmax);
  auto pid = (TH1F*)(phd->Clone("pid"));
  pid->SetTitle("PulseIntegral around peak / #Deltat;mV (equiv);frequency");

  for (int ievt=0; ievt<tree->GetEntries(); ++ievt){
    tree->GetEntry(ievt);
    baseline = volts[ch][ipeak-iBLS]*vgain;
    double height=(volts[ch][ipeak]-volts[ch][iBLS])*vgain;
    phd->Fill(height);
    double sum=0;
    // integrate over fixed region around nominal peak
    for (int n=istart; n<istop; ++n) {
      double V=(volts[ch][n]-volts[ch][iBLS])*vgain;
      sum+=(V); 
    }
    pid->Fill(sum * iScale);   // new
  }

  tcsum->cd(3);
  phd->DrawCopy();
  tcsum->cd(4);
  pid->DrawCopy();

  if (savePlot) tcsum->SaveAs(fnout);
  
  delete[] time;
  delete[] pool; 

  tfout->Write();
  tfout->Close();
  
  return;

  /*  
  // DCR analysis
  // loop over N buffers at a time to read a long sample window
  // Use TSpectrum to extract peaks
  int NBUF=50;
  auto buffer = new TH1F("buffer","buffer",LEN*NBUF,0,LEN*NBUF);
  int bbins=buffer->GetNbinsX();
  int offset=0;
  int ievt=0;
  while (ievt<NBUF){
    tree->GetEntry(ievt);
    for (int n=0; n<LEN; ++n){
       double V=-(volts[ch][n]*1000*vgain-baseline);
       buffer->SetBinContent(n+offset,V);
    }
    offset+=LEN;
    ievt++;
  }

  auto tcbuf=new TCanvas("rawbuffer","raw buffer");

  double *source = new double[bbins];
  for (int n=0; n<bbins; ++n) source[n]=buffer->GetBinContent(n+1);
  auto back=(TH1F*)buffer->Clone();
  auto TS=new TSpectrum(100*NBUF);

  TS->Background(source,bbins,150,TSpectrum::kBackDecreasingWindow,
		TSpectrum::kBackOrder2,kTRUE,
		TSpectrum::kBackSmoothing15,kFALSE);
  
   for (int n = 0; n < bbins; n++) back->SetBinContent(n + 1,source[n]+baseline);  //HACK
   back->SetLineColor(kRed);
   buffer->Draw();
   back->Draw("SAME L");
   delete[] source;
   auto tcbuf2=new TCanvas("Backgroundsubtracted","Background subtracted");
   auto buffer2 = (TH1F*)buffer->Clone("buffer2");
   buffer2->Add(back,-1.0);
   buffer2->Draw();

   int nrebin=20;
   buffer2->Rebin(nrebin);
   int npeaks = TS->Search(buffer2,2,"nobackground noMarkov",0.05);
   cout << npeaks << " peaks found" << endl;
   Double_t *xpeaks;
   xpeaks = TS->GetPositionX();
   Double_t *ypeaks;
   ypeaks = TS->GetPositionY();
   /// hacks ///
   double ymin = 6; 
   double t0=0;
   double dt=1;   // ns
   ///
   /// sort in time
   auto idx = new int[npeaks];
   TMath::Sort(npeaks,xpeaks,idx,false);

   auto hdt=new TH1I("hdt","hdt",100,0,3000);
   double last=xpeaks[idx[0]];
   double sum=0;
   int count=0;
   for (int n=1; n<npeaks; ++n){
     int i=idx[n];
     if (ypeaks[i]<ymin) continue;
     double dt=(xpeaks[i]-last)*nrebin*0.05;
     //cout << last << " " << xpeaks[i] << endl;
     hdt->Fill(dt);
     sum+=dt;
     count++;
     last=xpeaks[i];
   }
   
   auto tdt=new TCanvas("dt","dt");
   hdt->Fit("expo","","",300,3000);

   double meanDC=sum/count * nrebin * 0.05;
   cout << " Mean dt = " << sum/count * nrebin * 0.05 << " ns" << endl;
   cout << " DCR = " << 1/(meanDC*1e-9) << " Hz" << endl;
  */
}

