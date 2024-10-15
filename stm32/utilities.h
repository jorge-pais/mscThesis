#ifndef UTILITIES_H
#define UTILITIES_H

#include "arm_math.h"
#include "math.h"

#define Lharm 200
#define MAX_PULSE_BUFFER 512

// Glottal pulse magnitude and normalized relative delay
float32_t mag[] = {31391239.2593, 25509661.0703, 13796316.3286, 8922293.6888, 6308154.8657, 4404891.5172, 3472284.0635, 2669235.4559, 2123884.261, 1794168.7554, 1438579.5183, 1256119.7964, 1064617.8879, 909351.494, 820899.342, 696679.9187, 636641.6077, 564688.1864, 499617.7436, 466703.7485, 408347.812, 382744.1622, 348445.4386, 314785.245, 300195.5905, 267675.447, 254960.7555, 236048.6748, 216192.8902, 209044.4959, 188808.9749, 181845.4439, 170343.3578, 157536.0452, 153840.5031, 140244.8755, 136170.3901, 128664.5681, 119853.1088, 117912.34, 108249.0459, 105754.3466, 100589.7114, 94225.0894, 93234.9735, 86065.3583, 84490.951, 80787.3204, 76012.2687, 75560.7831, 70058.6868, 69047.0618, 66301.718, 62608.7523, 62471.7061, 58132.8765, 57478.7198, 55387.757, 52459.8529, 52509.4422, 49010.61, 48590.4151, 46961.4123, 44591.6236, 44752.2983, 41877.6718, 41614.2702, 40320.6185, 38368.8194, 38594.6482, 36195.2235, 36038.8583, 34994.503, 33362.8311, 33625.2958, 31595.1877, 31512.92, 30657.7238, 29275.9977, 29557.1314, 27819.2147, 27788.7495, 27079.664, 25896.3862, 26184.7763, 24681.6564, 24687.6169, 24093.173, 23069.79, 23358.179, 22046.3345, 22077.9458, 21574.7086, 20681.8431, 20965.6514, 19811.528, 19861.1383, 19431.3668, 18646.2795, 18922.6756, 17900.0185, 17962.1333, 17592.1981, 16897.0391, 17164.3427, 16252.3354, 16322.9735, 16002.2642, 15382.8443, 15640.1371, 14822.0774, 14898.3335, 14618.491, 14063.4012, 14310.2704, 13572.6118, 13652.3553, 13406.7216, 12906.6936, 13143.0599, 12474.7075, 12556.3706, 12339.5903, 11887.0229, 12113.0248, 11504.81, 11587.2384, 11394.9642, 10983.5697, 11199.4851, 10643.7686, 10726.1137, 10554.7871, 10179.3228, 10385.5151, 9875.8842, 9957.5249, 9804.2105, 9460.2729, 9657.155, 9188.1886, 9268.6745, 9130.9335, 8814.7991, 9002.8081, 8569.8942, 8648.9029, 8524.6957, 8233.1959, 8412.7766, 8011.9678, 8089.2741, 7976.8843, 7707.3069, 7878.9007, 7506.7996, 7582.2513, 7480.2262, 7230.2373, 7394.2743, 7047.9419, 7121.4424, 7028.5457, 6796.1266, 6953.0214, 6629.9028, 6701.3973, 6616.5719, 6399.968, 6550.1168, 6247.9808, 6317.4464, 6239.7838, 6037.4633, 6181.2424, 5898.1325, 5965.5699, 5894.2858, 5704.9062, 5842.6713, 5576.8645, 5642.2925, 5576.7071, 5399.0865, 5531.1731, 5281.1468, 5344.597, 5284.1187, 5117.2133, 5243.9365, 5008.3402, 5069.8538, 5013.9659, 4856.8502, 4978.5061, 4756.1377, 4815.7629, 4764.0124, 4615.8631};
float32_t nrd[] = {0.0636, 0.5883, 0.0059, 0.4135, -0.1979, 0.1985, 0.5963, -0.0137, 0.3889, -0.2182, 0.1793, 0.5807, -0.0271, 0.3769, -0.2279, 0.1702, 0.573, -0.0339, 0.3706, -0.2332, 0.1652, 0.5686, -0.0379, 0.3668, -0.2365, 0.162, 0.5658, -0.0405, 0.3643, -0.2387, 0.1598, 0.5638, -0.0424, 0.3625, -0.2403, 0.1582, 0.5624, -0.0438, 0.3611, -0.2415, 0.157, 0.5612, -0.0448, 0.3601, -0.2424, 0.156, 0.5604, -0.0457, 0.3592, -0.2432, 0.1552, 0.5597, -0.0464, 0.3586, -0.2438, 0.1546, 0.5591, -0.0469, 0.358, -0.2443, 0.1541, 0.5586, -0.0474, 0.3575, -0.2447, 0.1536, 0.5582, -0.0478, 0.3571, -0.2451, 0.1532, 0.5578, -0.0482, 0.3568, -0.2454, 0.1529, 0.5575, -0.0485, 0.3565, -0.2457, 0.1526, 0.5572, -0.0487, 0.3562, -0.246, 0.1524, 0.557, -0.049, 0.3559, -0.2462, 0.1521, 0.5567, -0.0492, 0.3557, -0.2464, 0.1519, 0.5565, -0.0494, 0.3555, -0.2466, 0.1517, 0.5564, -0.0495, 0.3553, -0.2467, 0.1516, 0.5562, -0.0497, 0.3552, -0.2469, 0.1514, 0.5561, -0.0498, 0.355, -0.247, 0.1513, 0.5559, -0.05, 0.3549, -0.2471, 0.1511, 0.5558, -0.0501, 0.3548, -0.2473, 0.151, 0.5557, -0.0502, 0.3547, -0.2474, 0.1509, 0.5556, -0.0503, 0.3546, -0.2475, 0.1508, 0.5555, -0.0504, 0.3545, -0.2476, 0.1507, 0.5554, -0.0505, 0.3544, -0.2476, 0.1506, 0.5553, -0.0506, 0.3543, -0.2477, 0.1505, 0.5552, -0.0507, 0.3542, -0.2478, 0.1505, 0.5552, -0.0507, 0.3541, -0.2479, 0.1504, 0.5551, -0.0508, 0.3541, -0.2479, 0.1503, 0.555, -0.0509, 0.354, -0.248, 0.1503, 0.555, -0.0509, 0.3539, -0.248, 0.1502, 0.5549, -0.051, 0.3539, -0.2481, 0.1501, 0.5549, -0.051, 0.3538, -0.2481, 0.1501, 0.5548, -0.0511, 0.3538, -0.2482, 0.15, 0.5548, -0.0511, 0.3537, -0.2482, 0.15, 0.5547, -0.0512, 0.3537};

// F0 contour template function (just for Fs = 22050 Hz)
#define contourLen 162
float32_t contour[] = {100.0, 108.025, 113.55, 117.337, 119.806, 121.253, 121.879, 121.833, 121.236, 120.186, 118.767, 117.059, 115.129, 113.046, 110.865, 108.658, 106.461, 104.34, 102.341, 100.506, 98.87, 97.459, 96.297, 95.392, 94.748, 94.356, 94.201, 94.261, 94.51, 94.919, 95.457, 96.092, 96.8, 97.552, 98.324, 99.098, 99.858, 100.588, 101.275, 101.914, 102.495, 103.015, 103.471, 103.861, 104.185, 104.442, 104.635, 104.767, 104.838, 104.855, 104.819, 104.735, 104.607, 104.44, 104.238, 104.007, 103.749, 103.472, 103.177, 102.87, 102.556, 102.238, 101.919, 101.604, 101.295, 100.995, 100.708, 100.434, 100.175, 99.934, 99.71, 99.507, 99.322, 99.158, 99.013, 98.888, 98.781, 98.693, 98.622, 98.567, 98.527, 98.502, 98.489, 98.487, 98.495, 98.513, 98.537, 98.569, 98.605, 98.645, 98.689, 98.735, 98.781, 98.828, 98.875, 98.921, 98.966, 99.009, 99.05, 99.088, 99.123, 99.155, 99.184, 99.21, 99.232, 99.252, 99.268, 99.281, 99.292, 99.299, 99.304, 99.306, 99.306, 99.304, 99.3, 99.294, 99.287, 99.278, 99.268, 99.257, 99.245, 99.232, 99.219, 99.205, 99.191, 99.177, 99.163, 99.149, 99.135, 99.121, 99.107, 99.094, 99.08, 99.068, 99.055, 99.043, 99.031, 99.02, 99.009, 98.999, 98.989, 98.979, 98.97, 98.961, 98.952, 98.943, 98.935, 98.927, 98.92, 98.912, 98.905, 98.898, 98.891, 98.884, 98.877, 98.87, 98.863, 98.857, 98.85, 98.843, 98.837, 98.83};

// Complex number type
typedef struct{
    float32_t real;
    float32_t imag;
} complex_t;
//} __attribute__(__packed) complex_t; // packed to avoid memory padding (allowing this to be cast to and from float32_t)

#define getPulseLen(fs, f0) ((uint16_t) floorf(1/f0 * fs + 0.5))

/// @brief Generate a single LF pulse and store it in the output buffer out
/// @param fs Sample rate
/// @param f0 Fundamental Frequency
/// @param pulseBuff Which pulse buffer to write to
/// @return Size of output vector (actual pulse len)
uint16_t genPulse(float32_t fs, float32_t f0, float32_t * pulseBuff, float32_t gain){
    float32_t T0 = 1/f0;
    float32_t Ts = 1/fs;

    // compute the output length
    uint16_t outLen = getPulseLen(fs, f0); 
    float32_t max = 0;

    // zero out the pulse buffer before writing to it
    memset(pulseBuff, 0.0f, MAX_PULSE_BUFFER);

    float32_t magn, phase, argsin;

    // Calculate and accumulate each harmonic
    for (int32_t i = 0; i < Lharm; i++){
        // avoid computing these three every cycle
        magn = mag[i]; 
        phase = 2 * PI * nrd[i];
        argsin = (i+1) * 2 * PI * f0 * Ts;

        // generate entire waveform
        for (int j = 0; j < outLen; j++)
            pulseBuff[j] += magn * sin(argsin * j + phase); 
    }
    
    // Determine max sample for normalization
    float32_t sample;
    for (int32_t i = 0; i < outLen; i++){
        sample = (float32_t) fabs((double) pulseBuff[i]);
        
        if(sample > max) max = sample; 
    }
    
    // normalize
    float32_t scale = gain * (INT16_MAX) / max;
		//float32_t scale = gain / max; // debug
    for (int32_t i = 0; i < outLen; i++)
        pulseBuff[i] = scale * pulseBuff[i];
    
    return outLen;
}

/// @brief Calculate the tilt gain compensation for a given frequency
/// @param t out - gain array
/// @param n in - twice the length of the gain array
/// @param fs in - sampling frequency
/// @param f0 in - fundamental frequency
void tiltGain(float32_t* t, uint16_t n, float32_t fs, float32_t f0){
    // Tilt gain model approximation

    const float32_t m1_IP = -5.968431149938542f, b1_IP = 0.0f, m2_IP = -32.9522014311249f, b2_IP = 8.122924250743393f;
    float32_t n2 = n/2;
    float32_t ratio;

    for (int32_t i = 0; i < n2; i++){
        ratio = i * (fs / f0) / n;

        if (ratio < 1.0f)
            t[i] = 1.0f;
        else if (ratio < 2.0f)
            t[i] = (pow(ratio, m1_IP/20)) * pow(10, b1_IP/20);
        else
            t[i] = (pow(ratio, m2_IP/20)) * pow(10, b2_IP/20);
				
        t[i] = 1/t[i]; // invert the gain, as arm only has the function to multiply to numbers
        t[n - 1 - i] = t[i]; // copy the flipped second half of the filter
    }

    return;
}

/// @brief Implement a direct form II IIR filter with no zeros (all-pole)
/// @param in input array
/// @param out 
/// @param len length of the input array
/// @param b0 gain coefficient
/// @param a denominator coefficients
/// @param coeffLen number of denominator coefficients
void iirDirectTypeII(float32_t* in, float32_t* out, int len, float32_t b0, float32_t* a, int coeffLen){
    // missing out of bounds check

    // set output to all zeros
    memset(out, 0.0f, len);
	
    // variable length arrays are allowed in C99
    float32_t w[coeffLen]; // delay line initialized to all zeros
    memset(w, 0.0f, coeffLen);

    for (int32_t i = 0; i < len; i++){
        // First difference equation is to determine w[n] from the denominator
        w[0] = in[i];
        for(int32_t j = 1; j < coeffLen; j++) // compute the output delay line
            w[0] -= a[j] * w[j];

        // Second difference equation is simplified as there is only b0
        out[i] = b0 * w[0];

        // Shift delay line
        for (int32_t j = coeffLen - 1; j > 0; j--)
            w[j] = w[j - 1];
    }

    return;
}

/// @brief Levinson-Durbin algorithm for LPC
/// @param r IN - Auto-correlation coefficients
/// @param n_len IN - Length of the auto-correlation sequence
/// @param a OUT - LPC coefficients
/// @param e OUT - Prediction Error
void levinson(float32_t* r, uint16_t n_len, float32_t* a, float32_t* e){
    float32_t lambda[n_len + 1];
    a[0] = 1.0;
    *e = r[0];

    float32_t lambdaSum, temp;
    for (int32_t i = 1; i <= n_len; i++){
        lambdaSum = 0.0f;

        for (int32_t j = 1; j <= i; j++)
            lambdaSum += r[j] * a[i - j];

        lambda[i] = -lambdaSum / (*e);

        for (int32_t j = 1; j <= i / 2; j++){
            temp = a[j] + lambda[i] * a[i - j];
            a[i - j] += lambda[i] * a[j];
            a[j] = temp;
        }

        a[i] = lambda[i];
        (*e) = (*e) * (1 - lambda[i] * lambda[i]);
    }

    return;
}

/// @brief Compute the spectral centre of gravity from a frame's PSD
/// @param psd Power Spectral Density
/// @param n Total number of samples (1024)
/// @param fs Sampling Frequency
/// @return Spectral Centre of Gravity
float32_t centreOfGravity(float32_t* psd, int n, float32_t fs){
    float32_t num = 0, dem = 0;

    for (size_t i = 0; i < n / 2; i++){
        dem += psd[i];
        num += (i * fs / n) * psd[i];
    }

    return num / dem;
}

/// @brief Mapping sigmoid for computing the centre of gravity
/// @param cog Centre of gravity frequency
/// @param Fvuv Vocalic whispered threshold
/// @param Fuv Unvoiced frequency threshold
/// @param tau Sharpness coefficient (0.99)
/// @return mapped mixing coefficient
float32_t mapping(float32_t cog, float32_t Fvuv, float32_t Fuv, float32_t tau){
    float32_t k = (cog - sqrtf(Fvuv*Fuv)) * 2/(Fuv - Fvuv) * log10((1+tau)/(1-tau));    
    
    return 1/(1 + exp(-k));
}
////
//// Random number shenanigans
////
static uint32_t seed = 1337;
#define LCG_A 1664525       // LCG multiplier
#define LCG_C 1013904223    // LCG increment
#define LCG_M 0xFFFFFFFF    // LCG modulus 2^32 - 1

/// @brief Generate a pseudo-random integer using a Linear Congruential Generator
/// @return Random integer
uint32_t randInt(){
    seed = (LCG_A * seed + LCG_C) % LCG_M;
    return seed;
}

/// @brief Generate a pseudo-random float
/// @return Random 32bit float
float32_t randFloat(){
    return 2.0f * ((float32_t) randInt() / (float32_t) LCG_M) - 1.0f;
}

/// @brief Computes the Weighted Average Decision (WAD) based on energy and zero-crossing rate.
/// @param frame Pointer to an array of 512 float32_t elements representing the signal frame.
/// @param z1 Energy threshold below which the function returns 0.0f.
/// @param z2 Energy threshold above which the function returns 1.0f.
/// @param zCross Zero-crossing rate threshold to decide the WAD when energy is between z1 and z2.
/// @return Returns 1.0f if the energy is above z2 or if the zero-crossing rate is above zCross when energy is between z1 and z2. Returns 0.0f otherwise.
float32_t WAD(float32_t* frame, float32_t z1, float32_t z2, float32_t zCross) {
    float energy = 0.0f;
    int zeroCrossings = 0;
    
    for (int i = 0; i < 512; i++)
        energy += frame[i] * frame[i];
    
    // Check energy against thresholds
    if (energy >= z2)
        return 1.0f;
    else if (energy < z1)
        return 0.0f;

    // Calculate the zero-crossing rate
    for (int i = 0; i < 511; i++)
        if (frame[i] * frame[i + 1] < 0)
            zeroCrossings++;
    
    float zeroCrossingRate = (float)zeroCrossings / 512.0f;
    
    // Return WAD based on zero-crossing rate
    return (zeroCrossingRate > zCross) ? 1.0f : 0.0f;
}


#endif // UTILITIES_H