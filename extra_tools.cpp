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
std::vector<std::vector<std::complex<double>>> read_file(){
    //Test_Image.mat
    std::vector<std::vector<std::complex<double>>> v;
    std::string line;
    vector<complex<double>> temp;
    stringstream ss;
    complex<double> ch;


    std::ifstream in("Test_Image.txt"); // окрываем файл для чтения
    if (in.is_open()){
        while (getline(in, line))
        {
            ss << line;
            while(ss>>ch){
                temp.push_back(ch);
            }
            v.push_back(temp);
            temp.clear();
        }
    }
    return v;
}

void Write_in_file(vector<vector<complex<double>>>& v){
    std::ofstream out;          // поток для записи
    out.open("SAR_out.txt");
    size_t azimuth_size = v.size();
    size_t range_size = v[0].size();

    for(size_t i = azimuth_size - 1; i >= 0;i--){
        for(size_t j = 0; j < range_size;j++){
            out << abs(v[i][j]) << ' ';
        }
        out << '\n';
    }
    out.close();
}






