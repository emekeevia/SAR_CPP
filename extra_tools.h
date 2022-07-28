#pragma once

#include <iostream>
#include <cmath>
#include <array>
#include <sstream>
#include <vector>
#include <fftw3.h>
#include <complex>
#include <algorithm>
#include <fstream>

using namespace std;



vector<double> fill_up(const double& start,const double& finish, const double& step);


void conjugate(vector<complex<double>>& sopr, const vector<complex<double>>& orig);

vector<complex<double>> operator*(const vector<complex<double>>& v1, const vector<complex<double>>& v2);

void iFFTshift(vector<complex<double>>& v);

std::vector<std::vector<std::complex<double>>> read_file();

void Write_in_file(vector<vector<complex<double>>>& v);


template<typename T>
T operator+(T v1, T v2){
	if(v1.size() != v2.size()){
		throw -1;
	}
	size_t count = v2.size();
	T temp;
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]+v2[i];
	}
	return temp;
}

template<typename T>
vector<T> operator-(vector<T>& v1, vector<T>& v2){
	if(v1.size() != v2.size()){
		throw -1;
	}
	size_t count = v1.size();
	T temp(count);
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}

template<typename T, int l>
array<T,l> operator-(array<T,l>& v1, array<T,l>& v2){
	size_t count = v1.size();
	array<T,l> temp;
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}



//template<typename T, typename D>
//T operator*(const T& v, D a){
//	T temp;
//	for(size_t i = 0; i < 3;i++){
//		temp[i] = v[i] * a;
//	}
//	return temp;
//}

//template<typename T>
//T operator*(const T& v1, const T& v2){
//    T temp;
//    if(v1.size() != v2.size()){
//        throw("Not equal size!");
//    }
//    for(size_t i = 0; i < v1.size();i++){
//        temp[i] = v1[i]*v2[i];
//    }
//    return temp;
//}


template<typename T>
double abs_vect(const T& v){
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}




template<typename T>
double scal_mult(const T& v1,const T& v2){
	return static_cast<double>(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
	os << "{";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}
	}
	os << "}";
	return os;
}

template<typename T, long long unsigned int l>
ostream& operator<<(ostream& os, const array<T,l>& v){
	os << "[";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}

	}
	os << "]";
	return os;
}


