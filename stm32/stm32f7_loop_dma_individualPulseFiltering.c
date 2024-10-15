// stm32f7_loop_dma_FunSP.c
#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "lfCoeffs.h"
#include "math.h"
#include "arm_math.h"
#include "arm_const_structs.h"

#define SOURCE_FILE_NAME "stm32f7_loop_dma_FunSP.c"

#define PULSE_GAIN 0.75f
#define N 1024 //1024
#define N2 512
#define FS 22050.0f // AUDIO_FREQUENCY_22K
#define F0 110.0f   // 110 male - 220 female
#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE

extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

#define LPCorder 22

// TODO - figure out how to make the compiler pack the values of the struct to ensure correct casting to float32_t
typedef struct{
    float32_t real;
    float32_t imag;
} complex_t;

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

// Window and odft exponential constants
float32_t sinW[N];
complex_t dirExp[N];
complex_t invExp[N];

complex_t currBuff[N];
complex_t tempBuff[N]; 
float32_t oldInBuff[N2] = {0.0f};
float32_t oldOutBuff[N2] = {0.0f};

/// @brief Processes one DMA transfer block worth of data
/// @param mainBuf Pointer to the main buffer
void process_buffer(float32_t *mainBuf) {
    int16_t *rx_buf, *tx_buf;

    if (rx_buffer_proc == PING)
        rx_buf = (int16_t *)PING_IN;
    else
        rx_buf = (int16_t *)PONG_IN;

    if (tx_buffer_proc == PING)
        tx_buf = (int16_t *)PING_OUT;
    else
        tx_buf = (int16_t *)PONG_OUT;

    int16_t sampleL, sampleR, out;
    for (int32_t i = 0; i < (PING_PONG_BUFFER_SIZE); i++){
        // read both input samples
        sampleL = *(rx_buf++); sampleR = *(rx_buf++);
			
        // overlapAdd the output buffer (possible to save one multiplation here!)
        *(tx_buf++) = (int16_t) floor(mainBuf[i] + 0.5f);
        //*(tx_buf++) = (int16_t) floor(oldOutBuff[i] * sinW[i] + mainBuf[i] * sinW[i + N2] + 0.5); // this is giving chopping audio idk why
        *(tx_buf++) = sampleR; // RIGHT channel remains untouched

        oldOutBuff[i] = (float32_t) mainBuf[i];
		
        // overlapAdd the input buffer for the currbuff
        currBuff[i].real = sinW[i] * oldInBuff[i];
        currBuff[i].imag = 0.0f;
        currBuff[i + N2].real = ((float32_t) sampleL) * sinW[i + N2];
        currBuff[i + N2].imag = 0.0f;

        oldInBuff[i] = (float32_t) sampleL;
    }

    RX_buffer_full = 0;
    TX_buffer_empty = 0;
}

int main(void) {
    // initialize the board and the audio codec
    stm32f7_wm8994_init(AUDIO_FREQUENCY_22K, // 22050 Hz
                    IO_METHOD_DMA,
                    INPUT_DEVICE_INPUT_LINE_1,
                    OUTPUT_DEVICE_HEADPHONE,
                    WM8994_HP_OUT_ANALOG_GAIN_0DB,
                    WM8994_LINE_IN_GAIN_0DB,
                    WM8994_DMIC_GAIN_9DB,
                    SOURCE_FILE_NAME,
                    NOGRAPH);

    uint16_t idxMainBuf = 0, vibrato = 0, contourIndex = 0; // reset main buffer index
    uint16_t gPulseLen, lastPulseLen, filteredLength;

    float32_t mainBuf[2 * PING_PONG_BUFFER_SIZE] = {0.0f};
    float32_t gPulseBuf[MAX_PULSE_BUFFER] = {0.0f}; // 220 ~ 22050/100
    float32_t colorPulse[MAX_PULSE_BUFFER] = {0.0f};
		
    float32_t tilt[N] = {0.0f};

    // * generate single prototype pulse 10% shorter than that of the fundamental period
    //gPulseLen = genPulse(FS, F0 * 1.1f, gPulseBuf, PULSE_GAIN);
    gPulseLen = genPulse(FS, F0, gPulseBuf, PULSE_GAIN);
		
    // pre compute the pulse spectral tilt of the filter (1/lfgain)
    tiltGain(tilt, N, FS, F0);
    
		// * precompute sin window and odft exponentials
    for (uint32_t i = 0; i < N; i++) {
        sinW[i] = (float32_t) sin(PI / N * (i + 0.5f));
			
        dirExp[i].real = (float32_t) cos(-PI * i / N);
        dirExp[i].imag = (float32_t) sin(-PI * i / N);

        invExp[i].real = (float32_t) cos(PI * i / N);
        invExp[i].imag = (float32_t) sin(PI * i / N);
    }

    float32_t autocorrelation[LPCorder + 1] = {0.0f};
    float32_t aCoeff[LPCorder + 1] = {0.0f};
    float32_t error = 0.0f, b0 = 0.0f;

    float32_t currPSD[N] = {0.0f};
    float32_t lastPSD[N] = {0.0f};
    float32_t pAlpha = 0.0f;
    
    for(;;){
        while (!(RX_buffer_full && TX_buffer_empty)){}

        //* 1. Copy first half to the output during the DMA processing
        process_buffer(mainBuf);

        //* 1.5. Compute filter from currBuff
        arm_cmplx_mult_cmplx_f32((float32_t *)(currBuff), (float32_t *)(dirExp), (float32_t *)(tempBuff), N); //odft exponential //working fine
        arm_cfft_f32(&arm_cfft_sR_f32_len1024, (float32_t *)(tempBuff), 0, 1); // ifft flag = 0, doBitReverse = 1
        
        arm_cmplx_mag_squared_f32((float32_t *)(tempBuff), (float32_t *)(currPSD), N); // Compute PSD 

        //* 2. Copy second half to first half
        for (int32_t i = 0; i < N2; i++){
            mainBuf[i] = mainBuf[i + N2];
            mainBuf[i + N2] = 0.0f;
        }

        // Place pulses inside the buffer
        while (idxMainBuf < N2){
            // follow f0 contourIndex
            //! we can save one division by instead of saving the contour freq, saving the period
            //lastPulseLen = getPulseLen(FS, contour[contourIndex]); // dummy pulse length just for placement
            lastPulseLen = getPulseLen(FS, F0); // dummy pulse length just for placement
            pAlpha = (idxMainBuf != 0) ? idxMainBuf/N2 : 1/N2 ; // pulse filtering coefficient (check for zero)
            
            // add vibrato
            //lastPulseLen += (uint16_t) floor((float32_t) 0.5f + (float32_t) 7.0f * sin((float32_t)2.0f * PI * vibrato / (float32_t)28.0f));

            //* 2.5. Generate a new filter based on the current pulse position
            for (uint32_t i = 0; i < N; i++){
                //tempBuff[i].real = pAlpha * currPSD[i] + (1 - pAlpha) * lastPSD[i];
                tempBuff[i].real = (pAlpha * currPSD[i] + (1 - pAlpha) * lastPSD[i]) * tilt[i]; // apply tilt gain
                tempBuff[i].imag = 0.0f;
            }

            // Use weiner khinchine (apply IODFT) to compute the autocorrelation sequence from the PSD
            arm_cfft_f32(&arm_cfft_sR_f32_len1024, (float32_t *)(tempBuff), 1, 1); // ifft flag = 1, doBitReverse = 1
            arm_cmplx_mult_cmplx_f32((float32_t *)(tempBuff), (float32_t *)(invExp), (float32_t *)(currBuff), N); // iodft

            // only first (LPCorder + 1) values are needed for levinson recursion
            for (uint32_t i = 0; i <= LPCorder; i++)
                autocorrelation[i] = tempBuff[i].real;    // assuming here the complex part will be zero

            levinson(autocorrelation, LPCorder, aCoeff, &error); // LPC analysis
            b0 = sqrtf(error); // gain factor

            //* 3. Filter and accumulate the pulse into the mainbuffer
            // for now just a dummy filtering operation as the pulse is filled with zeros at the end
            filteredLength = MAX_PULSE_BUFFER;
            //iirDirectTypeII(gPulseBuf, colorPulse, MAX_PULSE_BUFFER, 1, vowel, LPCorder + 1);
            iirDirectTypeII(gPulseBuf, colorPulse, MAX_PULSE_BUFFER, b0, aCoeff, LPCorder + 1);

            // 4. Copy *filtered* pulse to main buffer
            for (int32_t i = 0; i < filteredLength; i++)
                mainBuf[idxMainBuf++] += 0.00002f * colorPulse[i];

            // update main buffer index (this eliminates the need for a second condition?)
            idxMainBuf -= (filteredLength - lastPulseLen);
        }

        for (uint32_t i = 0; i < N; i++)
            lastPSD[i] = currPSD[i]; // copy values for next frame

        idxMainBuf -= N2;
    }
}
