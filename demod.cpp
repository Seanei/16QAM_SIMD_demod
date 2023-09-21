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

template <typename T, typename Y>
void Demodulate(const T& rx_symb, Y& d1, Y& d2, Y& d3, Y& d4, double noise_var){
    double real = rx_symb.real();
    double imag = rx_symb.imag();
    double real_abs = std::abs(real);
    double imag_abs = std::abs(imag);
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
}

template <typename T, typename Y>
void Demodulate_no_if(const T& rx_symb, Y& d1, Y& d2, Y& d3, Y& d4, double noise_var){
    float dist_to_symb = std::pow(std::abs(std::abs(rx_symb.real()) - 2 * STEPFORQAM16) - STEPFORQAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2 * STEPFORQAM16) - STEPFORQAM16, 2);
    auto step_for_real = (static_cast<int>(std::abs(rx_symb.real())/(2*STEPFORQAM16)) * 2 + 1) * STEPFORQAM16;
    auto step_for_imag = (static_cast<int>(std::abs(rx_symb.imag())/(2*STEPFORQAM16)) * 2 + 1) * STEPFORQAM16;
    d1 = sign(rx_symb.real()) * (std::pow((std::abs(rx_symb.real()) + STEPFORQAM16), 2) + std:: pow(std::abs(rx_symb.imag()) - step_for_imag, 2) - dist_to_symb);
    d2 = sign(rx_symb.imag()) * (std::pow((std::abs(rx_symb.real()) - step_for_real), 2) + std:: pow(std::abs(rx_symb.imag()) + STEPFORQAM16, 2) - dist_to_symb);
    d3 = (-1) * sign(std::abs(rx_symb.real()) - 2 * STEPFORQAM16) * (std::pow(std::abs(std::abs(rx_symb.real()) - 2* STEPFORQAM16) + STEPFORQAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2* STEPFORQAM16) - STEPFORQAM16, 2) - dist_to_symb);
    d4 = (-1) * sign(std::abs(rx_symb.imag()) - 2 * STEPFORQAM16) * (std::pow(std::abs(std::abs(rx_symb.real()) - 2* STEPFORQAM16) - STEPFORQAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2* STEPFORQAM16) + STEPFORQAM16, 2) - dist_to_symb);
}