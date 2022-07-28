#include "base_algorithm.h"

vector<complex<double>> range_chirp_for_correlation(int size_range, double fs, double K_r,double tau_p, bool approx){
    //(D.) Define correlation chirp in range
    // Basic definitions
    vector<complex<double>> range_chirp(size_range); //empty vector to be filled with chirp values
    vector<double> tau = fill_up(-tau_p/2, tau_p/2, 1/fs);  // time axis in range
    //Define chirp in range
    //Get size of chirp

    size_t size_chirp_r = tau.size();
    size_t start = (size_range - size_chirp_r)/2 - 1;
    double temp;
    int count;

    for(size_t i = 0;i < size_chirp_r;i++){
        if(approx){
            temp = tau[i]*tau[i] * K_r;//в единицах PI
            count = static_cast<int>(temp) / 2;
            temp -= count * 2.0;
            if(temp <= 0.5) {
                range_chirp[i+start] = complex<double>(1.0,1.0);
            }else if(temp > 0.5 && temp <= 1.0){
                range_chirp[i+start] = complex<double>(-1.0,1.0);
            }else if(temp > 1.0 && temp <= 1.5){
                range_chirp[i+start] = complex<double>(-1.0,-1.0);
            }else{
                range_chirp[i+start] = complex<double>(1.0,-1.0);
            }
            range_chirp[i+start] /= sqrt(2.0);
        }else{
            range_chirp[i+start] = exp(tau[i]*tau[i]*M_PI*K_r*static_cast<complex<double>>(I));
        }


    }

    //Position chirp in range vector (centered)
    // --> used for correlation procedure

    // Transform vector in frequency domain (Fourier transform)
    vector<complex<double>> RANGE_CHIRP(size_range);
    fftw_plan plan_f = fftw_plan_dft_1d(size_range, (fftw_complex*) &range_chirp[0],
                                        (fftw_complex*) &RANGE_CHIRP[0], FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan

    fftw_execute(plan_f); // Fourier Transform
    fftw_destroy_plan(plan_f);

    // Define chirp for correlation
    // --> conjugate complex of signal chirp
    vector<complex<double>> CON_RANGE_CHIRP(size_range);
    conjugate(CON_RANGE_CHIRP, RANGE_CHIRP);
    return CON_RANGE_CHIRP;
}

vector<complex<double>> azimuth_chirp_for_correlation(int size_azimuth,  double V, double Lambda, double R_0, double ta, double prf, bool approx){
    // E.) Define correlation chirp in azimuth

    // Basic definitions
    vector<complex<double>> azimuth_chirp(size_azimuth);  // empty vector to be filled with chirp values
    vector<double> t = fill_up(-ta/2, ta/2, 1/prf);// time axis in azimuth


    // FM Rate Azimuth Chirp
    double K_a = (-2*V*V)/(Lambda*R_0);

    // Define chirp in azimuth
    // Get size of chirp
    size_t size_chirp_a = t.size();
    size_t start2 = (size_azimuth - size_chirp_a)/2 - 1;

    double temp;
    int count;

    for(size_t i = 0;i < size_chirp_a;i++){
        if(approx){
            temp = t[i]*t[i] * K_a;//в единицах PI
            count = static_cast<int>(temp) / 2;
            temp -= count * 2.0;
            if(temp >= -0.5) {
                azimuth_chirp[i+start2] = complex<double>(1.0,-1.0);
            }else if(temp < -0.5 && temp >= -1.0){
                azimuth_chirp[i+start2] = complex<double>(-1.0,-1.0);
            }else if(temp < -1.0 && temp >= -1.5){
                azimuth_chirp[i+start2] = complex<double>(-1.0,1.0);
            }else{
                azimuth_chirp[i+start2] = complex<double>(1.0,1.0);
            }
            azimuth_chirp[i+start2] /= sqrt(2.0);
        }else{
            azimuth_chirp[i+start2] = exp(t[i]*t[i]*M_PI*K_a*static_cast<complex<double>>(I));
        }


    }




    // Position chirp in azimuth vector (centered)
    // --> used for correlation procedure


    // Transform vector in frequency domain (Fourier transform)
    fftw_plan plan;
    vector<complex<double>> AZIMUTH_CHIRP(size_azimuth);


    plan = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &azimuth_chirp[0],
                            (fftw_complex*)&AZIMUTH_CHIRP[0], FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
    fftw_execute(plan); // Fourier Transform
    fftw_destroy_plan(plan);
    // Define chirp for correlation
    // --> conjugate complex of signal chirp
    vector<complex<double>> CON_AZIMUTH_CHIRP(size_azimuth);
    conjugate(CON_AZIMUTH_CHIRP, AZIMUTH_CHIRP);
    return CON_AZIMUTH_CHIRP;
}

void SAR(){
    vector<vector<std::complex<double>>> Raw_data = read_file();

    // Get image size
    int size_azimuth = static_cast<int>(Raw_data.size());
    int size_range = static_cast<int>(Raw_data[0].size());

    vector<complex<double>> vek;
    // B.) Sensor parameters (ERS satellite)
    double fs = 18.962468e+6;  //Range Sampling Frequency [Hz]
    double K_r = 4.18989015e+11;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 37.12e-6;       // Chirp duration [s]
    double V = 7098.0194;                // Effective satellite velocity [m/s]
    double Lambda = 0.05656;            // Length of carrier wave [m]
    double R_0 = 852358.15;              // Range to center of antenna footprint [m]
    double ta = 0.6;                     // Aperture time [s]
    double prf = 1679.902;               // Pulse Repitition Frequency [Hz]

    //RANGE DOMAIN

    vector<complex<double>> CON_RANGE_CHIRP(size_range);
    CON_RANGE_CHIRP = range_chirp_for_correlation(size_range, fs, K_r, tau_p, true);

    //-------------------------------------------------------------------------------------------------------

    //AZIMUTH DOMAIN

    vector<complex<double>> CON_AZIMUTH_CHIRP(size_azimuth);
    CON_AZIMUTH_CHIRP = azimuth_chirp_for_correlation(size_azimuth, V, Lambda, R_0, ta, prf, true);

    //--------------------------------------------------------------------------------------------------------

    vector<vector<complex<double>>>processed(size_azimuth, vector<complex<double>>(size_range));// vector with finished image

    //Focusing
    //1)Range compression
    vector<complex<double>> vek_out(size_range);
    vector<complex<double>> vek_in(size_range);
    vector<complex<double>> if_vek(size_range);
    vector<complex<double>> CORR(size_range);
    vector<complex<double>> norm_r(size_range, complex<double>(1.0 / size_range));
    fftw_plan plan_r;
    for(size_t k1 = 0; k1 < size_azimuth;k1++) {
        vek_in = Raw_data[k1];
        plan_r = fftw_plan_dft_1d(size_range, (fftw_complex*) &vek_in[0],
                                  (fftw_complex*) &vek_out[0], FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
        fftw_execute(plan_r); // Fourier Transform
        CORR = vek_out * CON_RANGE_CHIRP;  // Multiply vectors
        plan_r = fftw_plan_dft_1d(size_range, (fftw_complex*) &CORR[0],
                                  (fftw_complex*) &if_vek[0], FFTW_BACKWARD, FFTW_ESTIMATE);// Inverse Fourier Transform
        fftw_execute(plan_r);
        if_vek = if_vek * norm_r;
        iFFTshift(if_vek);// Sorting after FFT
        processed [k1] = if_vek;  // Store row in matrix
    }
    fftw_destroy_plan(plan_r);
    //2)Azimuth compression
    vector<complex<double>> vek2(size_azimuth);// = processed[:, k2];  // Select row in azimuth
    vector<complex<double>> VEK2(size_azimuth);
    vector<complex<double>> if_vek2(size_azimuth);
    vector<complex<double>> CORR2(size_azimuth);
    vector<complex<double>> norm_a(size_azimuth, complex<double>(1.0 / size_azimuth));
    fftw_plan plan_a;
    for(size_t k2 = 0; k2 < size_range;k2++) {
        for(size_t i = 0; i < size_azimuth;i++){
            vek2[i] = processed[i][k2];
        }
        plan_a = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &vek2[0],
                                  (fftw_complex*) &VEK2[0], FFTW_FORWARD, FFTW_ESTIMATE);//making draft plan
        fftw_execute(plan_a); // Fourier Transform
        CORR2 = VEK2*CON_AZIMUTH_CHIRP;
        plan_a = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &CORR2[0],
                                  (fftw_complex*) &if_vek2[0], FFTW_BACKWARD, FFTW_ESTIMATE);// Inverse Fourier Transform
        fftw_execute(plan_a);
        if_vek2 = if_vek2 * norm_a;
        iFFTshift(if_vek2);// Sorting after FFT
        for(size_t i = 0; i < size_azimuth;i++){
            processed[i][k2] = if_vek2[i];
        }

    }
    fftw_destroy_plan(plan_a);
    Write_in_file(processed);
}