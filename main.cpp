#include <complex>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "demod.cpp"
#include "demodSIMD.cpp"

int main(){
    std::complex<double> num{0.9, 0.9};
    double d1 = 0;
    double d2 = 0;
    double d3 = 0;
    double d4 = 0;
    double noise_var = 1;
    std::vector<std::complex<double>> test2 = {{0.2, 0.9}, {-0.2,0.9}, {-0.2,-0.9}, {0.2,-0.9}};
    std::vector<float> test1_for_simd_i(16, -0.2);
    std::vector<float> test1_for_simd_q(16, 0.9);
    //int N = 4;
    std::vector<float> demodulated((test1_for_simd_i.size()), 0);

    for(const auto& n : test2){
        std::cout<<n<<std::endl;
        Demodulate(n, d1, d2, d3, d4, noise_var);
        std::cout << "Demodulated" << std::endl;
        std::cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<std::endl<<std::endl;
        std::cout << "Demodulate_no_if" << std::endl;
        Demodulate_no_if(n, d1, d2, d3, d4, noise_var);
        std::cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<std::endl<<std::endl;
    }
    Demodulate16QAMSIMD(test1_for_simd_i, test1_for_simd_q, demodulated);
    std::cout << "DEMOD len " << demodulated.size() << std::endl;
    for (auto elem : demodulated){
        std::cout << elem << " ";
    }
    
    std::cout<<std::endl;
}