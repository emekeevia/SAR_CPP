#include "extra_tools.h"


vector<complex<double>> fill_up(const double& start,const double& finish, const double& step){
    vector<complex<double>> mass;
    for(double temp = start; temp < finish;temp+=step){
        mass.push_back(temp);
    }
    return mass;
}

void conjugate(vector<complex<double>>& sopr, const vector<complex<double>>& orig){
    size_t size = orig.size();
    for(size_t i = 0; i < size;i++){
        sopr[i] = conj(orig[i]);
    }
}

vector<complex<double>> operator*(const vector<complex<double>>& v1, const vector<complex<double>>& v2){
    size_t size = v1.size();
    if(size != v2.size()){
        throw("Not equal size!");
    }
    vector<complex<double>> temp(size);
    for(size_t i = 0; i < size;i++){
        temp[i] = v1[i]*v2[i];
    }
    return temp;
}

void iFFTshift(vector<complex<double>>& v){
    size_t N = v.size();
    rotate(v.begin(),v.begin() + N/2,v.end());
}






