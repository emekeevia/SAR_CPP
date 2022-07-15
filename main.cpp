#include <iostream>
#include <fftw.h>
#include <complex.h>
#include <vector>


void read_file(){
    //Test_Image.mat
}
void SAR(){
    vector<vector<complex<double>>> Raw_data;
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();
    vector<complex<double>> vek;
    auto fs = 18.962468*10e6;  //Range Sampling Frequency [Hz]
    auto K_r = 4.18989015*10e11;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    auto tau_p = 37.12*10e-6;       // Chirp duration [s]
    auto V = 7098.0194;                // Effective satellite velocity [m/s]
    auto Lambda = 0.05656;            // Length of carrier wave [m]
    auto R_0 = 852358.15;              // Range to center of antenna footprint [m]
    auto ta = 0.6;                     // Aperture time [s]
    auto prf = 1679.902;               // Pulse Repitition Frequency [Hz]

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
