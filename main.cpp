#include <complex>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "demod.cpp"
#include "demodSIMD.cpp"
#include <cstdlib>
#include <ctime>
#include <chrono>

const uint BITS_PER_SYMB = 4;
const double RAND_MAX_DOUBLE = RAND_MAX;
const int NUM_OF_SIMULATION = 100;

int main(){
    std::vector<int> ns_power = {6,7,8,9,10,11,12,13,14,15,16};

    for (auto p : ns_power) {
        int NUMBER_OF_SYMBOLS = std::pow(2,p);
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds duration_demod, duration_demod_if, duration_demod_simd;


        //std::chrono::duration duration_demod, duration_demod_if, duration_demod_simd;
        for (int i = 0; i < NUM_OF_SIMULATION; ++i){

            std::srand((unsigned)std::time(0));
            double noise_var = 1;
            std::vector<float> test_for_simd_i, test_for_simd_q;
            test_for_simd_i.reserve(NUMBER_OF_SYMBOLS * BITS_PER_SYMB);
            test_for_simd_q.reserve(NUMBER_OF_SYMBOLS * BITS_PER_SYMB);
            std::vector<std::complex<float>> test;
            test.reserve(NUMBER_OF_SYMBOLS * BITS_PER_SYMB);
            std::vector<float> demodulated(NUMBER_OF_SYMBOLS * BITS_PER_SYMB, 0);
            std::vector<float> demodulated_no_if(NUMBER_OF_SYMBOLS * BITS_PER_SYMB, 0);
            std::vector<float> demodulated_SIMD(NUMBER_OF_SYMBOLS * BITS_PER_SYMB, 0);
            
            for(int n = 0; n < NUMBER_OF_SYMBOLS; ++n){
                float I = (std::rand() - (RAND_MAX_DOUBLE/2))/(RAND_MAX_DOUBLE/2);
                float Q = (std::rand() - (RAND_MAX_DOUBLE/2))/(RAND_MAX_DOUBLE/2);
                for (int i = 0; i < BITS_PER_SYMB; ++i){
                    test_for_simd_i.push_back(I);
                    test_for_simd_q.push_back(Q);
                }
                test.push_back({I,Q});
            }
            {
                auto start = std::chrono::high_resolution_clock::now();
                Demodulate(test, demodulated, noise_var);
                auto end = std::chrono::high_resolution_clock::now();
                duration_demod += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            }
            
            {
                auto start = std::chrono::high_resolution_clock::now();
                Demodulate_no_if(test, demodulated_no_if, noise_var);
                auto end = std::chrono::high_resolution_clock::now();
                duration_demod_if += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            }
            
            {
                auto start = std::chrono::high_resolution_clock::now();
                Demodulate16QAMSIMD(test_for_simd_i, test_for_simd_q, demodulated_SIMD);
                auto end = std::chrono::high_resolution_clock::now();
                duration_demod_simd += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            }
        }
        std::cout << "NUMBER_OF_SYMBOLS: 2^" << p <<std::endl;
        std::cout << "Execution time: " << (duration_demod.count())/NUM_OF_SIMULATION << " microseconds" << std::endl;
        std::cout << "Execution time: " << (duration_demod_if.count())/NUM_OF_SIMULATION << " microseconds" << std::endl;
        std::cout << "Execution time: " << (duration_demod_simd.count())/NUM_OF_SIMULATION << " microseconds" << std::endl;
        std::cout << "____________" << std::endl << std::endl;

    }

    /*
    for (int i = 0; i < 4; ++i){
        std::cout << "Demod: " << demodulated[i * BITS_PER_SYMB] << ", " << demodulated[i * BITS_PER_SYMB + 1] << ", " 
        << demodulated[i * BITS_PER_SYMB + 2] << ", " << demodulated[i * BITS_PER_SYMB + 3] << std::endl;
        std::cout << "Demod_no_if: " << demodulated_no_if[i * BITS_PER_SYMB] << ", " << demodulated_no_if[i * BITS_PER_SYMB + 1] << ", " 
        << demodulated_no_if[i * BITS_PER_SYMB + 2] << ", " << demodulated_no_if[i * BITS_PER_SYMB + 3] << std::endl;
        std::cout << "Demod_SIMD: " << demodulated_SIMD[i * BITS_PER_SYMB] << ", " << demodulated_SIMD[i * BITS_PER_SYMB + 1] << ", " 
        << demodulated_SIMD[i * BITS_PER_SYMB + 2] << ", " << demodulated_SIMD[i * BITS_PER_SYMB + 3] << std::endl;
    }
    */
}