#include <iostream>
#include <fftw3.h>
#include <stdio.h>
#include <complex.h>
#include <tgmath.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <array>
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
    vector<complex<double>> range_chirp(size_range); //empty vector to be filled with chirp values
    //vector<std::complex<double>> RANGE_CHIRP(size_range);
    vector<complex<double>> tau = fill_up(-tau_p/2, tau_p/2, 1/fs);  // time axis in range
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
    ra_chirp_temp.insert(ra_chirp_temp.begin(),size_range - size_chirp_r,0);
    range_chirp = ra_chirp_temp;

    // Transform vector in frequency domain (Fourier transform)
    fftw_plan plan_f;
    vector<complex<double>> RANGE_CHIRP(size_range);


    plan_f = fftw_plan_dft_1d(size_range, reinterpret_cast<fftw_complex*>(&range_chirp),
                              reinterpret_cast<fftw_complex*>(&RANGE_CHIRP), FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
    fftw_execute(plan_f); // Fourier Transform
    fftw_destroy_plan(plan_f);


    // Define chirp for correlation
    // --> conjugate complex of signal chirp
    vector<complex<double>> CON_RANGE_CHIRP(size_range);
    conjugate(CON_RANGE_CHIRP, RANGE_CHIRP);

    vector<vector<complex<double>>>processed(size_azimuth, vector<complex<double>>(size_range));


    for(size_t k1 = 0; k1 < size_azimuth;k1++) {
        size_t size = Raw_data[k1].size();
        vector<complex<double>> vek_out(size_range);
        fftw_plan plan = fftw_plan_dft_1d(size_range, reinterpret_cast<fftw_complex*>(&Raw_data[k1]),
                                           reinterpret_cast<fftw_complex*>(&vek_out), FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan
        fftw_execute(plan); // Fourier Transform
        vector<complex<double>> CORR = vek_out * CON_RANGE_CHIRP;  // Multiply vectors

        vector<complex<double>> if_vek(size_range);
        plan = fftw_plan_dft_1d(size_range, reinterpret_cast<fftw_complex*>(&CORR),
                                           reinterpret_cast<fftw_complex*>(&if_vek), FFTW_BACKWARD, FFTW_ESTIMATE);// Inverse Fourier Transform

        iFFTshift(if_vek);// Sorting after FFT
        processed [k1] = if_vek;  // Store row in matrix
        fftw_destroy_plan(plan);
    }
}

int main() {
    std::string file_name;

    return 0;
}
