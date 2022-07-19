#include <iostream>
#include <fftw.h>
#include <stdio.h>
#include <complex>
#include <tgmath.h>
#include <vector>
#include <cmath>
#include "extra_tools.h"


void read_file(){
    //Test_Image.mat
}
void SAR(){
    vector<vector<std::complex<double>>> Raw_data;
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();

    vector<complex<double>> vek;
    double fs = 18.962468*10e6;  //Range Sampling Frequency [Hz]
    double K_r = 4.18989015*10e11;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 37.12*10e-6;       // Chirp duration [s]
    double V = 7098.0194;                // Effective satellite velocity [m/s]
    double Lambda = 0.05656;            // Length of carrier wave [m]
    double R_0 = 852358.15;              // Range to center of antenna footprint [m]
    double ta = 0.6;                     // Aperture time [s]
    double prf = 1679.902;               // Pulse Repitition Frequency [Hz]

    //RANGE DOMAIN
    // Basic definitions
    vector<std::complex<double>> range_chirp(size_range); //empty vector to be filled with chirp values
    vector<std::complex<double>> RANGE_CHIRP(size_range);
    vector<std::complex<double>> tau = fill_up(-tau_p/2, tau_p/2, 1/fs);  // time axis in range
    //vector<std::complex<double>> omega = fill_up(-fs/2, fs/2, 1/tau_p);  // frequency axis in range

    //Define chirp in range
    vector<std::complex<double>> ra_chirp_temp (tau.size());
    for(size_t i = 0;i < tau.size();i++){
        ra_chirp_temp[i] = exp(tau[i]*tau[i]*M_PI*K_r*static_cast<complex<double>>(I));
    }

    //Get size of chirp
    size_t size_chirp_r = tau.size();

    //Position chirp in range vector (centered)
    // --> used for correlation procedure
    size_t index_start = ceil(((size_range-size_chirp_r)/2.0)-1.0);
    size_t index_end = size_chirp_r+ceil(((size_range-size_chirp_r)/2.0)-2.0);
    //range_chirp[0, index_start:index_end+1] = ra_chirp_temp;

    // Transform vector in frequency domain (Fourier transform)
    fftw_plan plan_f;
    plan_f = fftw_create_plan(size_range, FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
    fftw_one(plan_f, range_chirp, RANGE_CHIRP); // Fourier Transform
    fftw_destroy_plan(plan_f);


    // Define chirp for correlation
    // --> conjugate complex of signal chirp
    //CON_RANGE_CHIRP = np.conjugate(RANGE_CHIRP)


    for(size_t k1 = 0; k1 < size_azimuth;k1++) {
        fftw_complex vek_in = Raw_data[k1];  // Select row in range
        fftw_complex vek_out[size_range];
        fftw_plan plan = fftw_create_plan(size_range, FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
        fftw_one(plan, vek_in, vek_out); // Fourier Transform
        auto CORR = vek_out * CON_RANGE_CHIRP;  // Multiply vectors
        if_vek = np.fft.ifft(CORR);  // Inverse Fourier Transform
        if_vek_sort = np.fft.ifftshift(if_vek);// Sorting after FFT
        processed [k1] = if_vek_sort;  // Store row in matrix
        fftw_destroy_plan(plan);
    }
}

int main() {
    std::string file_name;

    return 0;
}
