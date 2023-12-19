#include <cmath>

const double STEPFORQAM16 = 1/sqrt(10);

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0)); 
}

template <typename T>
bool GetSvAndSign(const T& rx_symb, T& ref_sym, T& sv0, T& sv1, T& sv2, T& sv3){
    if (rx_symb.real <= -2*STEPFORQAM16 && rx_symb.im >= 2*STEPFORQAM16 ){
        ref_sym = {-3*STEPFORQAM16, 3 *STEPFORQAM16};
        sv0 = {3*STEPFORQAM16, 3 *STEPFORQAM16};
        sv1 = {-3*STEPFORQAM16, -3 *STEPFORQAM16};
        sv2 = {-1*STEPFORQAM16, 3 *STEPFORQAM16};
        sv3 = {-3*STEPFORQAM16, 1 *STEPFORQAM16};
    }else if (rx_symb.real <= 0*STEPFORQAM16 && rx_symb.im >= 2*STEPFORQAM16 ){
        ref_sym = {-1*STEPFORQAM16, 3 *STEPFORQAM16};
        sv0 = {1*STEPFORQAM16, 3 *STEPFORQAM16};
        sv1 = {-3*STEPFORQAM16, 3 *STEPFORQAM16};
        sv2 = {-1*STEPFORQAM16, -3 *STEPFORQAM16};
        sv3 = {-1*STEPFORQAM16, 1 *STEPFORQAM16};
    }else if (rx_symb.real >= 2*STEPFORQAM16 && rx_symb.im >= 2*STEPFORQAM16 ){
        ref_sym = {3*STEPFORQAM16, 3 *STEPFORQAM16};
        sv0 = {1*STEPFORQAM16, 3 *STEPFORQAM16};
        sv1 = {-3*STEPFORQAM16, 3 *STEPFORQAM16};
        sv2 = {-1*STEPFORQAM16, -3 *STEPFORQAM16};
        sv3 = {-1*STEPFORQAM16, 1 *STEPFORQAM16};
    }
}

template <typename T>
void Demodulate(const std::vector<std::complex<T>>& rx_symbols, std::vector<T>& demodulated, double noise_var){
    const uint BITS_PER_SYMB = 4;
    T d1, d2, d3, d4;
    for (uint i = 0; i < rx_symbols.size(); ++i){
        T real = rx_symbols[i].real();
        T imag = rx_symbols[i].imag();
        T real_abs = std::abs(real);
        T imag_abs = std::abs(imag);
        if (real_abs>= 2 * STEPFORQAM16 && imag_abs >= 2 * STEPFORQAM16){
            auto d = - pow((real_abs - 3*STEPFORQAM16),2) - pow((imag_abs - 3*STEPFORQAM16),2);
            d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((imag_abs + 1*STEPFORQAM16),2) + pow((real_abs - 3*STEPFORQAM16),2)); 
            d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((imag_abs - 3*STEPFORQAM16),2) + pow((real_abs + 1*STEPFORQAM16),2)); 
            d3 = -1 * (1/(noise_var))*( d + pow((imag_abs - 1*STEPFORQAM16),2) + pow((real_abs - 3*STEPFORQAM16),2));  
            d4 = -1 * (1/(noise_var))*( d + pow((imag_abs - 3*STEPFORQAM16),2) + pow((real_abs - 1*STEPFORQAM16),2));
        } 
        else if (real_abs <= 2 * STEPFORQAM16 && imag_abs >= 2 * STEPFORQAM16){
            auto d = - pow((real_abs - 1*STEPFORQAM16),2) - pow((imag_abs - 3*STEPFORQAM16),2);
            d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2)); 
            d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
            d3 = 1 * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));  
            d4 = -1 * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));
        }
        else if (real_abs >= 2 * STEPFORQAM16 && imag_abs <= 2 * STEPFORQAM16){
            auto d = - pow((real_abs - 3*STEPFORQAM16),2) - pow((imag_abs - 1*STEPFORQAM16),2);
            d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2)); 
            d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
            d3 = -1 * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));  
            d4 = 1 * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));
        }else {
            auto d = - pow((real_abs - 1*STEPFORQAM16),2) - pow((imag_abs - 1*STEPFORQAM16),2);
            d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2)); 
            d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
            d3 =  1 * (1/(noise_var)) * ( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));  
            d4 = 1 * (1/(noise_var)) * ( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));
        }
        demodulated[i*BITS_PER_SYMB] = std::move(d1);
        demodulated[i*BITS_PER_SYMB+1] = std::move(d2);
        demodulated[i*BITS_PER_SYMB+2] = std::move(d3);
        demodulated[i*BITS_PER_SYMB+3] = std::move(d4);

    }
}

template <typename T>
void Demodulate_no_if(const std::vector<std::complex<T>>& rx_symbols, std::vector<T>& demodulated, double noise_var){
    const uint BITS_PER_SYMB = 4;
    T re, im;
    T d1, d2, d3, d4;
    for (uint i = 0; i < rx_symbols.size(); ++i){
        re = rx_symbols[i].real();
        im = rx_symbols[i].imag();
        float dist_to_symb = std::pow(std::abs(std::abs(re) - 2 * STEPFORQAM16) - STEPFORQAM16, 2) + std::pow(std::abs(std::abs(im) - 2 * STEPFORQAM16) - STEPFORQAM16, 2);
        auto step_for_real = (static_cast<int>(std::abs(re)/(2*STEPFORQAM16)) * 2 + 1) * STEPFORQAM16;
        auto step_for_imag = (static_cast<int>(std::abs(im)/(2*STEPFORQAM16)) * 2 + 1) * STEPFORQAM16;
        d1 = sign(re) * (std::pow((std::abs(re) + STEPFORQAM16), 2) + std:: pow(std::abs(im) - step_for_imag, 2) - dist_to_symb);
        d2 = sign(im) * (std::pow((std::abs(re) - step_for_real), 2) + std:: pow(std::abs(im) + STEPFORQAM16, 2) - dist_to_symb);
        d3 = (-1) * sign(std::abs(re) - 2 * STEPFORQAM16) * (std::pow(std::abs(std::abs(re) - 2* STEPFORQAM16) + STEPFORQAM16, 2) + std::pow(std::abs(std::abs(im) - 2* STEPFORQAM16) - STEPFORQAM16, 2) - dist_to_symb);
        d4 = (-1) * sign(std::abs(im) - 2 * STEPFORQAM16) * (std::pow(std::abs(std::abs(re) - 2* STEPFORQAM16) - STEPFORQAM16, 2) + std::pow(std::abs(std::abs(im) - 2* STEPFORQAM16) + STEPFORQAM16, 2) - dist_to_symb);
        demodulated[i*BITS_PER_SYMB] = std::move(d1);
        demodulated[i*BITS_PER_SYMB+1] = std::move(d2);
        demodulated[i*BITS_PER_SYMB+2] = std::move(d3);
        demodulated[i*BITS_PER_SYMB+3] = std::move(d4);
    }
}