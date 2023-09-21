#include <immintrin.h>
#include <vector>

const float STEP_FOR_QAM16 = 1/sqrt(10);

void Demodulate16QAMSIMD(const std::vector<float>& inputI, const std::vector<float>& inputQ, std::vector<float>& output){
    // Input vector is 4 times repeated I and Q components
    const float mask1[] = {0, 0, -2*STEP_FOR_QAM16, -2*STEP_FOR_QAM16};
    const float mask2_i[] = {STEP_FOR_QAM16, -STEP_FOR_QAM16, STEP_FOR_QAM16, -STEP_FOR_QAM16};
    const float mask2_q[] = {-STEP_FOR_QAM16, STEP_FOR_QAM16, -STEP_FOR_QAM16, STEP_FOR_QAM16};
    const uint16_t mask_1 = 4369; //0001 0001 0001 0001 <- begin // mask for each 1st byte
    const uint16_t mask_2 = 8738; // 0010 0010 0010 0010 // mask for each 2nd byte
    const uint16_t mask_3 = 17476; // 0100 0100 0100 0100 // mask for each 3rd byte
    const uint16_t mask_4 = 34952; //1000 1000 1000 1000 // mask for each 4th byte

    auto mask1_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask1));
    auto mask2_i_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask2_i));
    auto mask2_q_reg =  _mm512_broadcast_f32x4(_mm_load_ps(mask2_q));

    constexpr auto FLOATS_IN_AVX_REGISTER = 16u; 
    const auto vectorizableSamples = (inputI.size() / FLOATS_IN_AVX_REGISTER) * FLOATS_IN_AVX_REGISTER;
    auto i = 0u;
    for (; i < vectorizableSamples; i += FLOATS_IN_AVX_REGISTER) {
        // Distance to the closest symbol
        auto dist_i = _mm512_loadu_ps(inputI.data() + i); //Load data
        auto dist_q = _mm512_loadu_ps(inputQ.data() + i);

        dist_i = _mm512_abs_ps(dist_i); // First abs
        dist_q = _mm512_abs_ps(dist_q);

        dist_i = _mm512_sub_ps(dist_i, _mm512_set1_ps(2*STEP_FOR_QAM16)); // real - 2*STEP_FOR_QAM16
        dist_q = _mm512_sub_ps(dist_q, _mm512_set1_ps(2*STEP_FOR_QAM16)); // imag - 2*STEP_FOR_QAM16

        dist_i = _mm512_abs_ps(dist_i); // Second abs
        dist_q = _mm512_abs_ps(dist_q);

        dist_i = _mm512_sub_ps(dist_i, _mm512_set1_ps(STEP_FOR_QAM16)); // real - STEP_FOR_QAM16
        dist_q = _mm512_sub_ps(dist_q, _mm512_set1_ps(STEP_FOR_QAM16)); // imag - STEP_FOR_QAM16

        dist_i = _mm512_mul_ps (dist_i, dist_i); //Power of 2
        dist_q = _mm512_mul_ps (dist_q, dist_q);  

        auto dist =   _mm512_add_ps(dist_i, dist_q); // real + imag
        
        // Distance to neghbour
        auto iRegister = _mm512_loadu_ps(inputI.data() + i);  //Load data
        auto qRegister = _mm512_loadu_ps(inputQ.data() + i);

        auto flags_i_sign_1 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_1), iRegister, _mm512_set1_ps(0) , _CMP_LT_OS); //Find flags for sign. Compare real with 0 according to mask_1 for sign of 1st bit
        auto flags_q_sign_1 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_2), qRegister, _mm512_set1_ps(0) , _CMP_LT_OS); //Find flags for sign. Compare imag with 0 according to mask_2 for sign of 2nd bit

        iRegister = _mm512_abs_ps (iRegister); //First abs
        qRegister = _mm512_abs_ps (qRegister);

        auto flags_i_mask = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_2), iRegister, _mm512_set1_ps(2*STEP_FOR_QAM16), _CMP_GT_OS); //Find flags for resize if im or real are greater then 2*STEP_FOR_QAM16. It needs only for 1st bit for real and 2nd for imag. 
        auto flags_q_mask = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_1), qRegister, _mm512_set1_ps(2*STEP_FOR_QAM16), _CMP_GT_OS);

        auto flags_i_sign_2 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_3), iRegister, _mm512_set1_ps(2*STEP_FOR_QAM16) , _CMP_GT_OS); //Find flags for sign. Compare real with 2*STEP_FOR_QAM16 according to mask_3 for sign of 3rd bit
        auto flags_q_sign_2 = _mm512_mask_cmp_ps_mask(_cvtu32_mask16(mask_4), qRegister, _mm512_set1_ps(2*STEP_FOR_QAM16) , _CMP_GT_OS); //Find flags for sign. Compare real with 2*STEP_FOR_QAM16 according to mask_4 for sign of 4th bit
        auto sign = _kor_mask16(flags_i_sign_1, flags_q_sign_1); // OR for all sign masks
        sign = _kor_mask16(sign, flags_i_sign_2);
        sign = _kor_mask16(sign, flags_q_sign_2);

        iRegister = _mm512_add_ps (iRegister, mask1_reg); // real[3,4] - 2 * STEP_FOR_QAM16
        qRegister = _mm512_add_ps (qRegister, mask1_reg); // imag[3,4] - 2 * STEP_FOR_QAM16

        iRegister = _mm512_abs_ps (iRegister); //Second abs
        qRegister = _mm512_abs_ps (qRegister);
        // perform the addition
        iRegister = _mm512_add_ps (iRegister, mask2_i_reg); // real[1,3] + STEP_FOR_QAM16, real[2,4] - STEP_FOR_QAM16
        qRegister = _mm512_add_ps (qRegister, mask2_q_reg); // imag[2,4] + STEP_FOR_QAM16, real[1,3] - STEP_FOR_QAM16

        iRegister = _mm512_mask_sub_ps(iRegister, flags_i_mask, iRegister, _mm512_set1_ps(2*STEP_FOR_QAM16)); // real - 2*STEP_FOR_QAM16 according to mask flags_i_mask
        qRegister = _mm512_mask_sub_ps(qRegister, flags_q_mask, qRegister, _mm512_set1_ps(2*STEP_FOR_QAM16)); // imag - 2*STEP_FOR_QAM16 according to mask flags_q_mask

        iRegister = _mm512_mul_ps (iRegister, iRegister); //Power of 2
        qRegister = _mm512_mul_ps (qRegister, qRegister);
        
        auto result = _mm512_add_ps(iRegister, qRegister); // real + imag

        // Distance_to_neighbour - Distance_to_closest
        result = _mm512_sub_ps(result, dist); 
        result = _mm512_mask_sub_ps(result, sign, _mm512_set1_ps(0.0), result);

        // store result in the output vector
        _mm512_storeu_ps(output.data() + i, result);
    }
}