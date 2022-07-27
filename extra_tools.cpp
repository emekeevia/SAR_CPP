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
    std::string line, elem;
    std::vector<std::complex<double>> temp;
    std::stringstream ss, micro_ss;
    std::complex<double> ch;
    double real, im;
    char sign;
    //int count = 0;




    std::ifstream in_; // окрываем файл для чтения
    in_.open("/home/ivan/CLionProjects/SAR_CPP/Refreshed_test.txt", std::ios::in);
    if (in_.is_open()){
        while (getline(in_, line))
        {
            ss << line;
            while(ss>>elem){
                micro_ss << elem;
                micro_ss >> real;
                micro_ss.get(sign);
                micro_ss >> im;
                if(sign == '-') {
                    im *= -1;
                }
                micro_ss.str(std::string());
                ch = complex<double>(real, im);
//                if(count < 2){
//                    count++;
//                    std::cout << elem << endl;
//                    std::cout << real << endl;
//                    std::cout << sign << endl;
//                    std::cout << im << endl;
//                    std::cout << "'" <<micro_ss.str() << "'" <<endl;
//                }
                temp.push_back(ch);
            }
            v.push_back(temp);
            temp.clear();
        }
    }
    in_.close();
    return v;
}

void Write_in_file(vector<vector<complex<double>>>& v){
    std::ofstream out;          // поток для записи
    out.open("SAR_out.txt");
    size_t azimuth_size = v.size();
    size_t range_size = v[0].size();

    for(size_t i = 0; i < azimuth_size;i++){
        for(size_t j = 0; j < range_size;j++){
            if(j != range_size - 1){
                out << abs(v[i][j]) << ' ';
            }else {
                out << abs(v[i][j]);
            }
        }
        out << '\n';
    }
    out.close();
}






