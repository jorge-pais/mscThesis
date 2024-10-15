// stm32f7_template.c

#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "lfCoeffs.h"

#define SOURCE_FILE_NAME "stm32f7_DMA_template.c"

#define FS AUDIO_FREQUENCY_22K
#define N PING_PONG_BUFFER_SIZE // 512 samples
#define PULSE_GAIN 0.75
#define MAX_PULSE_BUFFER 200 // 200 ~= (2*Fs / 220)

#define F0 110.0 // 110 Hz male - 220 Hz female 

/// @brief Generate a single LF pulse and store it in the output buffer out
/// @param fs Sample rate
/// @param f0 Fundamental Frequency
/// @param pulseBuff Which pulse buffer to write to
/// @return Size of output vector (actual pulse len)
uint16_t genPulse(float32_t fs, float32_t f0, float32_t * pulseBuff, float32_t gain){
    float32_t T0 = 1/f0;
    float32_t Ts = 1/fs;

    // compute the output length
    uint16_t outLen = floor(T0/Ts + 0.5); 
    float32_t max = 0;

    // zero out the pulse buffer before writting to it
    memset(pulseBuff, 0.0f, MAX_PULSE_BUFFER);

    float32_t magn, phase, argsin;

    // Calculate and accumulate each harmonic
    for (size_t i = 0; i < Lmag; i++){
        // avoid computing these three every cycle
        magn = mag[i]; 
        phase = 2 * PI * nrd[i];
        argsin = (i+1) * 2 * PI * f0 * Ts;

        // generate entire waveform
        for (int j = 0; j < outLen; j++)
            pulseBuff[j] += magn * sin(argsin * j + phase); 
    }
    
    // Determine max sample for normalization >()
    float32_t sample;
    for (size_t i = 0; i < outLen; i++){
        sample = (float32_t) fabs((double) pulseBuff[i]);
        
        if(sample > max) max = sample; 
    }
    
    // normalize
    float32_t scale = gain * (INT16_MAX) / max;
    for (size_t i = 0; i < outLen; i++)
        pulseBuff[i] = scale * pulseBuff[i];
    
    return outLen;
}

// DMA buffer processing variables
extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

/// @brief This function processes one DMA transfer block worth of data
/// @return Pointer to the output buffer
int16_t* process_buffer(void){
    // in / out buffers, PING_PONG_BUFFER_SIZE * (L + R) samples long
    int16_t *rx_buf, *tx_buf;

    // Select input and output buffers
    if (rx_buffer_proc == PING)
        rx_buf = (int16_t *)PING_IN;
    else
        rx_buf = (int16_t *)PONG_IN;
    
    if (tx_buffer_proc == PING) 
        tx_buf = (int16_t *)PING_OUT;
    else 
        tx_buf = (int16_t *)PONG_OUT;

    // copy the input buffer to the output buffer
    for (int i = 0; i < PING_PONG_BUFFER_SIZE; i++){
        *tx_buf++ = *rx_buf++;
        *tx_buf++ = *rx_buf++;
    }
    
    RX_buffer_full = 0;
    TX_buffer_empty = 0;

    // pointer to the start of the buffer
    // simpler than 
    return (tx_buf - 2 * PING_PONG_BUFFER_SIZE);
}

int main(void){
    // setup the audio codec and enable DMA processing 
    stm32f7_wm8994_init(AUDIO_FREQUENCY_22K, //22050 Hz
                        IO_METHOD_DMA,
                        INPUT_DEVICE_INPUT_LINE_1,
                        OUTPUT_DEVICE_HEADPHONE,
                        WM8994_HP_OUT_ANALOG_GAIN_0DB,
                        //WM8994_HP_OUT_ANALOG_ATTEN_9DB,
                        WM8994_LINE_IN_GAIN_0DB,
                        WM8994_DMIC_GAIN_9DB,
                        SOURCE_FILE_NAME,
                        NOGRAPH);
		
    int16_t *tx_buf;
    int16_t idxMainBuf = 0; // buffer indices
    float32_t mainBuf[2 * N] = {0.0f}; // 200 ~ 2*22050/220

    float32_t gPulseBuf[MAX_PULSE_BUFFER] = {0.0f}; 
    int16_t lastPulseLen = genPulse(FS, F0, gPulseBuf, PULSE_GAIN);
    for(;;){

        while (!(RX_buffer_full && TX_buffer_empty)){}
        tx_buf = process_buffer();
    
        // 1. Copy first half to output
        // second half to first half
        for (size_t i = 0; i < N; i++){
            *(tx_buf) = mainBuf[i];
            tx_buf += 2; // increment twice to skip right channel
            mainBuf[i] = mainBuf[i + N];
            mainBuf[i + N] = 0.0f;
        }
        
        // Fill first half of the buffer
        //while(idxMainBuf < N && !(idxMainBuf + lastPulseLen >= N)){
        while(idxMainBuf < N){

            for (size_t i = 0; i < lastPulseLen; i++) // this works under the presumption that the pulse isn't longer than 200 samples long!
                mainBuf[idxMainBuf++] = gPulseBuf[i];
        }
        idxMainBuf -= N; // point to first half of the main buffer
    }   
}

/* int main(void)
{
    int16_t idxMainBuf = 0, vibrato = 0; // reset main buffer index
    int16_t k, idxptr = 0, GpulseLen, GpulseLen_local;
    float32_t MainBuffer[2 * PING_PONG_BUFFER_SIZE] = {0.0};
    float32_t GpulseBuf[MAX_PULSE_BUFFER] = {0.0}; // 220 ~ 22050/100

    stm32f7_wm8994_init(AUDIO_FREQUENCY_22K, // 22050 Hz
                        IO_METHOD_DMA,
                        INPUT_DEVICE_INPUT_LINE_1,
                        OUTPUT_DEVICE_HEADPHONE,
                        WM8994_HP_OUT_ANALOG_GAIN_0DB,
                        WM8994_LINE_IN_GAIN_0DB,
                        WM8994_DMIC_GAIN_9DB,
                        SOURCE_FILE_NAME,
                        NOGRAPH);
    // generate single prototype of glottal pulse
    GpulseLen = genPulse(FS, F0, GpulseBuf, PULSE_GAIN);

    while (1){

        while (!(RX_buffer_full && TX_buffer_empty)){}
        process_buffer(MainBuffer);

        // left-shift second half of MainBuffer and reset second-half
        for (k = 0; k < N; k++){
            *(MainBuffer + k) = *(MainBuffer + N + k);
            *(MainBuffer + N + k) = 0.0;
        }

        while ((idxMainBuf < N) && (1)) // "1" is to include condition of next glottal pulse still start before index N
        {
            idxptr = 0; // reset glottal pulse buffer pointer
            GpulseLen_local = GpulseLen + (int16_t)floor((float32_t)0.5 + (float32_t)7.0 * sin((float32_t)2.0 * PI * vibrato / (float32_t)28.0));
            while (idxptr < GpulseLen_local)
            {
                MainBuffer[idxMainBuf++] = GpulseBuf[idxptr++]; // 0.0025 is to be adapted, later
            }
            vibrato = (vibrato++) % 28;
        }
        idxMainBuf -= N;
    }
} */