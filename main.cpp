#include <complex>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <immintrin.h>

const float STEP_FOR_QAM16 = 1/sqrt(10);
const float HUGE_NUMB_FOR_COMP = 100;
const double STEPFORQAM16 = 1/sqrt(10);

template <typename T> int sign(T val) {
    //std::cout << "SIGN " << (T(0) < val) - (val < T(0)) << std::endl;
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
        std::cout << "d  " << d <<std::endl;
        d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((imag_abs + 1*STEPFORQAM16),2) + pow((real_abs - 3*STEPFORQAM16),2)); 
        d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((imag_abs - 3*STEPFORQAM16),2) + pow((real_abs + 1*STEPFORQAM16),2)); 
        d3 = -1 * (1/(noise_var))*( d + pow((imag_abs - 1*STEPFORQAM16),2) + pow((real_abs - 3*STEPFORQAM16),2));  
        d4 = -1 * (1/(noise_var))*( d + pow((imag_abs - 3*STEPFORQAM16),2) + pow((real_abs - 1*STEPFORQAM16),2));
    } 
    else if (real_abs <= 2 * STEPFORQAM16 && imag_abs >= 2 * STEPFORQAM16){
        auto d = - pow((real_abs - 1*STEPFORQAM16),2) - pow((imag_abs - 3*STEPFORQAM16),2);
        std::cout << "d  " << d <<std::endl;
        d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2)); 
        d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
        d3 = 1 * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));  
        d4 = -1 * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));
    }
    else if (real_abs >= 2 * STEPFORQAM16 && imag_abs <= 2 * STEPFORQAM16){
        auto d = - pow((real_abs - 3*STEPFORQAM16),2) - pow((imag_abs - 1*STEPFORQAM16),2);
        std::cout << "d  " << d <<std::endl;
        d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2)); 
        d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
        d3 = -1 * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));  
        d4 = 1 * (1/(noise_var))*( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));
    }else {
        auto d = - pow((real_abs - 1*STEPFORQAM16),2) - pow((imag_abs - 1*STEPFORQAM16),2);
        std::cout << "d  " << d <<std::endl;
        d1 = (real >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs + 1*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2)); 
        d2 = (imag >=0 ? 1 : -1) * (1/(noise_var))*( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs + 1*STEPFORQAM16),2)); 
        d3 =  1 * (1/(noise_var)) * ( d + pow((real_abs - 3*STEPFORQAM16),2) + pow((imag_abs - 1*STEPFORQAM16),2));  
        d4 = 1 * (1/(noise_var)) * ( d + pow((real_abs - 1*STEPFORQAM16),2) + pow((imag_abs - 3*STEPFORQAM16),2));
    }
}

template <typename T, typename Y>
void Demodulate_no_if(const T& rx_symb, Y& d1, Y& d2, Y& d3, Y& d4, double noise_var){
    float dist_to_symb = std::pow(std::abs(std::abs(rx_symb.real()) - 2 * STEP_FOR_QAM16) - STEP_FOR_QAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2 * STEP_FOR_QAM16) - STEP_FOR_QAM16, 2);
    std::cout << "dist_to_symb  " << dist_to_symb <<std::endl;
    std::cout << static_cast<int>(std::abs(rx_symb.real())/(2*STEP_FOR_QAM16)) <<std::endl;

    auto step_for_real = (static_cast<int>(std::abs(rx_symb.real())/(2*STEP_FOR_QAM16)) * 2 + 1) * STEP_FOR_QAM16;
    std::cout << std::abs(rx_symb.imag())/(2*STEP_FOR_QAM16)/(2*STEP_FOR_QAM16) <<std::endl;
    auto step_for_imag = (static_cast<int>(std::abs(rx_symb.imag())/(2*STEP_FOR_QAM16)) * 2 + 1) * STEP_FOR_QAM16;
    d1 = sign(rx_symb.real()) * (std::pow((std::abs(rx_symb.real()) + STEP_FOR_QAM16), 2) + std:: pow(std::abs(rx_symb.imag()) - step_for_imag, 2) - dist_to_symb);
    d2 = sign(rx_symb.imag()) * (std::pow((std::abs(rx_symb.real()) - step_for_real), 2) + std:: pow(std::abs(rx_symb.imag()) + STEP_FOR_QAM16, 2) - dist_to_symb);
    d3 = (-1) * sign(std::abs(rx_symb.real()) - 2 * STEP_FOR_QAM16) * (std::pow(std::abs(std::abs(rx_symb.real()) - 2* STEP_FOR_QAM16) + STEP_FOR_QAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2* STEP_FOR_QAM16) - STEP_FOR_QAM16, 2) - dist_to_symb);
    d4 = (-1) * sign(std::abs(rx_symb.imag()) - 2 * STEP_FOR_QAM16) * (std::pow(std::abs(std::abs(rx_symb.real()) - 2* STEP_FOR_QAM16) - STEP_FOR_QAM16, 2) + std::pow(std::abs(std::abs(rx_symb.imag()) - 2* STEP_FOR_QAM16) + STEP_FOR_QAM16, 2) - dist_to_symb);
}


void DemodulateSIMD(const std::vector<float>& inputI, const std::vector<float>& inputQ, std::vector<float>& output){
    constexpr auto FLOATS_IN_AVX_REGISTER = 16u;
    const auto vectorizableSamples = (inputI.size() / FLOATS_IN_AVX_REGISTER) 
                                     * FLOATS_IN_AVX_REGISTER;
    auto i = 0u;
    const float mask1[] = {0, 0, -2*STEP_FOR_QAM16, -2*STEP_FOR_QAM16};
    const float mask2_i[] = {STEP_FOR_QAM16, -STEP_FOR_QAM16, STEP_FOR_QAM16, -STEP_FOR_QAM16};
    const float mask2_q[] = {-STEP_FOR_QAM16, STEP_FOR_QAM16, -STEP_FOR_QAM16, STEP_FOR_QAM16};
    const float flags_for_mask2_i[] = {2*STEP_FOR_QAM16, HUGE_NUMB_FOR_COMP,  HUGE_NUMB_FOR_COMP, HUGE_NUMB_FOR_COMP};
    const float flags_for_mask2_q[] = {HUGE_NUMB_FOR_COMP, 2*STEP_FOR_QAM16, HUGE_NUMB_FOR_COMP, HUGE_NUMB_FOR_COMP};
    const uint16_t mask_for_sign_i1 = 4369; //0001 0001 0001 0001 <- begin
    const uint16_t mask_for_sign_q1 = 8738; // 0010 0010 0010 0010
    const uint16_t mask_for_sign_i2 = 17476; // 0100 0100 0100 0100
    const uint16_t mask_for_sign_q2 = 34952; //1000 1000 1000 1000


    auto mask1_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask1));
    auto mask2_i_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask2_i));
    auto mask2_q_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask2_q));
    auto mask1_dist = _mm512_set1_ps(-2*STEP_FOR_QAM16);
    auto flags_for_mask2_i_reg =  _mm512_broadcast_f32x4(_mm_load_ps(flags_for_mask2_i));
    auto flags_for_mask2_q_reg =  _mm512_broadcast_f32x4(_mm_load_ps(flags_for_mask2_i));


    auto mask2_dist = _mm512_set1_ps(-STEP_FOR_QAM16);

    for (; i < vectorizableSamples; i += FLOATS_IN_AVX_REGISTER) {
        // Distance to the closest
        auto dist_i = _mm512_loadu_ps(inputI.data() + i);
        auto dist_q = _mm512_loadu_ps(inputQ.data() + i);

        dist_i = _mm512_abs_ps(dist_i);
        dist_q = _mm512_abs_ps(dist_q);

        dist_i = _mm512_add_ps(dist_i, mask1_dist);
        dist_q = _mm512_add_ps(dist_q, mask1_dist);

        dist_i = _mm512_abs_ps(dist_i);
        dist_q = _mm512_abs_ps(dist_q);

        dist_i = _mm512_add_ps(dist_i, mask2_dist);
        dist_q = _mm512_add_ps(dist_q, mask2_dist);

        dist_i = _mm512_mul_ps (dist_i, dist_i); //Power of 2
        dist_q = _mm512_mul_ps (dist_q, dist_q);   
        auto dist =   _mm512_add_ps(dist_i, dist_q); 
        
        // Distance to neghbour
        auto iRegister = _mm512_loadu_ps(inputI.data() + i);
        auto qRegister = _mm512_loadu_ps(inputQ.data() + i);

        auto flags_i_sign_1 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_for_sign_i1), iRegister, _mm512_set1_ps(0) , _CMP_LT_OS);
        auto flags_q_sign_1 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_for_sign_q1), qRegister, _mm512_set1_ps(0) , _CMP_LT_OS);

        iRegister = _mm512_abs_ps (iRegister); //First abs
        qRegister = _mm512_abs_ps (qRegister);

        auto flags_i_mask = _mm512_cmp_ps_mask(iRegister, flags_for_mask2_i_reg, _CMP_GT_OS);
        auto flags_q_mask = _mm512_cmp_ps_mask(qRegister, flags_for_mask2_q_reg, _CMP_GT_OS);

        auto flags_i_sign_2 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_for_sign_i2), iRegister, _mm512_set1_ps(2*STEP_FOR_QAM16) , _CMP_GT_OS);
        auto flags_q_sign_2 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_for_sign_q2), qRegister, _mm512_set1_ps(2*STEP_FOR_QAM16) , _CMP_GT_OS);
        auto sign = _kor_mask16(flags_i_sign_1, flags_q_sign_1);
        sign = _kor_mask16(sign, flags_i_sign_2);
        sign = _kor_mask16(sign, flags_q_sign_2);

        iRegister = _mm512_add_ps (iRegister, mask1_reg); //First Mask
        qRegister = _mm512_add_ps (qRegister, mask1_reg);

        iRegister = _mm512_abs_ps (iRegister); //Second abs
        qRegister = _mm512_abs_ps (qRegister);
        // perform the addition
        iRegister = _mm512_add_ps (iRegister, mask2_i_reg); //Second Mask
        qRegister = _mm512_add_ps (qRegister, mask2_q_reg);

        iRegister = _mm512_mask_sub_ps(iRegister, flags_i_mask, iRegister, flags_for_mask2_i_reg);
        qRegister = _mm512_mask_sub_ps(qRegister, flags_q_mask, qRegister, flags_for_mask2_q_reg);

        iRegister = _mm512_mul_ps (iRegister, iRegister); //Power of 2
        qRegister = _mm512_mul_ps (qRegister, qRegister);
        
        auto result = _mm512_add_ps(iRegister, qRegister);
        result = _mm512_sub_ps(result, dist);
        result = _mm512_mask_sub_ps(result, sign, _mm512_set1_ps(0.0), result);
        // store data back in the data vector
        _mm512_storeu_ps(output.data() + i, result);
    }
}


int main(){
    std::complex<double> num{0.9, 0.9};
    double d1 = 0;
    double d2 = 0;
    double d3 = 0;
    double d4 = 0;
    double noise_var = 1;
    std::vector<std::complex<double>> test2 = {{0.2, 0.9}, {-0.2,0.9}, {-0.2,-0.9}, {0.2,-0.9}};
    //std::vector<std::complex<double>> test3 = {{0.9, 0.2}, {-0.9,0.2}, {-0.9,-0.2}, {0.9,-0.2}};
    //std::vector<std::complex<double>> test1 = {{0.2, 0.2}, {-0.2,0.2}, {-0.2,-0.2}, {0.2,-0.2}};
    std::vector<float> test1_for_simd_i(16, -0.2);
    std::vector<float> test1_for_simd_q(16, -0.9);
    //int N = 4;
    std::vector<float> demodulated((test1_for_simd_i.size()), 0);
/*
    std::ifstream in("data.bin", std::ios_base::in | std::ios_base::binary);
    int nsamps = 900;
    std::vector<std::complex<double> > sx(nsamps);
    in.read(reinterpret_cast<char*>(sx.data()), 2*nsamps * sizeof(double));
    in.close();
    for (int i = 0; i<4;++i){
        std::cout<<sx[i]<<std::endl;
    }
    */
   //Demodulate(num, d1, d2, d3, d4, noise_var);

    for(const auto& n : test2){
        std::cout<<n<<std::endl;
        Demodulate(n, d1, d2, d3, d4, noise_var);
        std::cout << "Demodulated" << std::endl;
        std::cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<std::endl<<std::endl;
        std::cout << "Demodulate_no_if" << std::endl;
        Demodulate_no_if(n, d1, d2, d3, d4, noise_var);
        std::cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<std::endl<<std::endl;
    }
    DemodulateSIMD(test1_for_simd_i, test1_for_simd_q, demodulated);
    std::cout << "DEMOD len " << demodulated.size() << std::endl;
    for (auto elem : demodulated){
        std::cout << elem << " ";
    }
    
    std::cout<<std::endl;
}