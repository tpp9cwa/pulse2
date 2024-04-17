struct PulseParams {
    constexpr static size_t size() { return sizeof(PulseParams) / sizeof(double); }

    PulseParams()
        : x0(0), A(0), a(0), tau(0), pedestal(0)
    {
        static_assert(sizeof(PulseParams) == size() * sizeof(double));
    }

    static PulseParams* from_array(double* arr) {
        return reinterpret_cast<PulseParams*>(arr);
    }

    double* to_array() {
        return reinterpret_cast<double*>(this);
    }

    double x0, A, a, tau, pedestal;

    double pulse(double x) const {
        double value = pedestal;
        if (x > x0) {
            value += A * pow((x - x0), a) * exp(-(x - x0) / tau);
        }
        return value;
    }
};

struct PulseFitter {
    PulseFitter(TH1* hist)
        : hist(hist)
    {
        x_min = hist->GetBinCenter(1);
        x_max = hist->GetBinCenter(hist->GetNbinsX());
        
        const double initial_pedestal = hist->GetBinContent(1);

        i_peak = hist->GetMaximumBin();
        x_peak = hist->GetBinCenter(i_peak);
        y_peak = hist->GetBinContent(i_peak);

        const double delay = 10; // ns delay for start of rise

        for (int i = i_peak; i > 0; i--) {
            if (std::fabs(hist->GetBinContent(i) - initial_pedestal) < 0.1 * (y_peak - initial_pedestal)) {
                pedestal_end = hist->GetBinCenter(i) - delay;
                i_pedestal_end = hist->FindBin(pedestal_end);
                break;
            }
        }

        // Fit the pedestal

        TF1* pedestal_fit = new TF1("pedestal",
                [] (double* t, double* par) {
                    return par[0];
                },
                x_min, pedestal_end, 1);

        pedestal_fit->SetParameter(0, initial_pedestal);
        pedestal_fit->SetRange(x_min, pedestal_end);
        hist->Fit("pedestal", "NR");

        pedestal = pedestal_fit->GetParameter(0);

        fit_pulse();
    }
    
    void fit_pulse() {
        pulse_fit = new TF1("pulse",
                [] (double* t, double* par) {
                    return PulseParams::from_array(par)->pulse(t[0]);
                },
                x_min, x_max, PulseParams::size());

        pulse_fit->SetLineColor(kRed - 4);
     

        PulseParams initial;
        initial.x0 = pedestal_end;
        initial.pedestal = pedestal;
        initial.A = (y_peak - initial.pedestal) / (x_peak - initial.x0);
        initial.a = 1;
        initial.tau = (x_peak - initial.x0) / (initial.a * log(x_peak - initial.x0));

        pulse_fit->SetParameter(0, initial.x0);
        pulse_fit->SetParameter(1, initial.A);
        pulse_fit->SetParameter(2, initial.a);
        pulse_fit->SetParameter(3, initial.tau);
        pulse_fit->FixParameter(4, initial.pedestal);

        pulse_fit->SetRange(pedestal_end, x_peak);

        hist->Fit("pulse", "NR");
    }

    PulseParams* parameters() { return PulseParams::from_array(pulse_fit->GetParameters()); }

    TH1* hist;
    TF1* pulse_fit;

    // Peaks
    double x_peak;
    double y_peak;
    int i_peak;

    // Bounds
    double x_min;
    double x_max;

    // Pedestal
    int i_pedestal_end;
    double pedestal_end;
    double pedestal;
};
