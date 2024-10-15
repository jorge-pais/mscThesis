// stm32f7_loop_dma_FunSP.c
#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "lfCoeffs.h"

#define SOURCE_FILE_NAME "stm32f7_loop_dma_FunSP.c"

#define PULSE_GAIN 0.75
#define N 2*PING_PONG_BUFFER_SIZE //1024
#define N2 PING_PONG_BUFFER_SIZE
#define FS 22050.0 // AUDIO_FREQUENCY_22K
#define F0 250.0   // 110 male - 220 female
#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE

extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

#define LPCorder 22
// /a/ vowel template with tilt compensation (F0 = 258Hz)
float32_t vowel[23] = {1.0, -3.56018445, 7.97021317, -12.71042867, 16.03850806, 
-15.68614206, 11.31567925, -3.71767768, -3.93117113, 9.17894871, 
-9.93343347, 6.97790725, -1.82948168, -2.76171435, 5.46348437, 
-5.65824728, 4.37828272, -2.45242487, 0.98733188, -0.15435453, 
-0.02705114, 0.01881448, 0.01900302};

// TODO - figure out how to make the compiler pack the values of the struct to ensure correct casting to float32_t
typedef struct{
    float32_t real;
    float32_t imag;
} complex_t;

#define getPulseLen(fs, f0) ((uint16_t) floor(1/f0 * fs + 0.5))

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
    for (size_t i = 0; i < Lharm; i++){
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
    for (size_t i = 0; i < outLen; i++){
        sample = (float32_t) fabs((double) pulseBuff[i]);
        
        if(sample > max) max = sample; 
    }
    
    // normalize
    float32_t scale = gain * (INT16_MAX) / max;
		//float32_t scale = gain / max; // debug
    for (size_t i = 0; i < outLen; i++)
        pulseBuff[i] = scale * pulseBuff[i];
    
    return outLen;
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

    for (int i = 0; i < len; i++){
        // First difference equation is to determine w[n] from the denominator
        w[0] = in[i];
        for(int j = 1; j < coeffLen; j++) // compute the output delay line
            w[0] -= a[j] * w[j];

        // Second difference equation is simplified as there is only b0
        out[i] = b0 * w[0];

        // Shift delay line
        for (int j = coeffLen - 1; j > 0; j--)
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
    for (size_t i = 1; i <= N; i++){
        lambdaSum = 0.0f;

        for (size_t j = 1; j <= i; j++)
            lambdaSum += a[j] * r[i - j];

        lambda[i] = -lambdaSum / (*e);

        for (size_t j = 1; j <= i / 2; j++){
            temp = a[j] + lambda[i] * a[i - j];
            a[i - j] += lambda[i] * a[j];
            a[j] = temp;
        }

        a[i] = lambda[i];
        (*e) = (*e) * (1 - lambda[i] * lambda[i]);
    }

    return;
}

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

    int16_t sampleL, sampleR;
    for (int i = 0; i < (PING_PONG_BUFFER_SIZE); i++){
        // read both input samples
        sampleL = *(rx_buf++); sampleR = *(rx_buf++);

        // copy the input buffer to the output buffer
        *(tx_buf++) = (int16_t) floor((float32_t) 0.5 + mainBuf[i]);
        *(tx_buf++) = sampleR; // RIGHT channel remains untouched
    }

    RX_buffer_full = 0;
    TX_buffer_empty = 0;
}

int main(void) {

    uint16_t idxMainBuf = 0, vibrato = 0, contourIndex = 0; // reset main buffer index
    uint16_t gPulseLen, lastPulseLen, filteredLength;

    float32_t mainBuf[2 * PING_PONG_BUFFER_SIZE] = {0.0f};
    float32_t gPulseBuf[MAX_PULSE_BUFFER] = {0.0f}; // 220 ~ 22050/100
    float32_t colorPulse[MAX_PULSE_BUFFER] = {0.0f};

    stm32f7_wm8994_init(AUDIO_FREQUENCY_22K, // 22050 Hz
                        IO_METHOD_DMA,
                        INPUT_DEVICE_INPUT_LINE_1,
                        OUTPUT_DEVICE_HEADPHONE,
                        WM8994_HP_OUT_ANALOG_GAIN_0DB,
                        WM8994_LINE_IN_GAIN_0DB,
                        WM8994_DMIC_GAIN_9DB,
                        SOURCE_FILE_NAME,
                        NOGRAPH);

    // generate single prototype pulse 10% shorter than that of the fundamental period
    gPulseLen = genPulse(FS, F0 * 1.1, gPulseBuf, PULSE_GAIN);
    //gPulseLen = genPulse(FS, F0, gPulseBuf, PULSE_GAIN);

    for(;;){
        while (!(RX_buffer_full && TX_buffer_empty)){}

        // 1. Copy first half to the output during the DMA processing
        process_buffer(mainBuf);

        // 2. Copy second half to first half
        for (size_t i = 0; i < N2; i++){
            mainBuf[i] = mainBuf[i + N2];
            mainBuf[i + N2] = 0.0f;
        }

        // Place pulses inside the buffer
        while (idxMainBuf < N2){
            // follow f0 contourIndex
            // we can save one division by instead of saving the contour freq, saving the period
            //lastPulseLen = getPulseLen(FS, contour[contourIndex]); // dummy pulse length just for placement
            lastPulseLen = getPulseLen(FS, F0); // dummy pulse length just for placement

            // add vibrato
            //lastPulseLen += (uint16_t) floor((float32_t) 0.5f + (float32_t) 7.0f * sin((float32_t)2.0f * PI * vibrato / (float32_t)28.0f));

            // 3. Filter and accumulate the pulse into the mainbuffer
            // for now just a dummy filtering operation as the pulse is filled with zeros at the end
            filteredLength = MAX_PULSE_BUFFER;
            iirDirectTypeII(gPulseBuf, colorPulse, MAX_PULSE_BUFFER, 1, vowel, LPCorder + 1);

            // 4. Copy *filtered* pulse to main buffer
            for (size_t i = 0; i < filteredLength; i++)
                mainBuf[idxMainBuf++] += 0.1 * colorPulse[i];

            // update main buffer index (this eliminates the need for a second condition?)
            idxMainBuf -= (filteredLength - lastPulseLen);
            
            vibrato = (vibrato+1) % 28;
            //contourIndex = (contourIndex < contourLen) ? contourIndex + 1 : contourIndex;
        }

        idxMainBuf -= N2;
    }
}
