#pragma once

#include <complex.h>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include "extra_tools.h"

vector<complex<double>> range_chirp_for_correlation(int size_range, double fs, double K_r,double tau_p,bool approx);
vector<complex<double>> azimuth_chirp_for_correlation(int size_azimuth,  double V, double Lambda, double R_0, double ta, double prf, bool approx);

void SAR();

