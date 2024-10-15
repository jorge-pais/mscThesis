// stm32f7_loop_dma_FunSP.c
#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
//#include "lfCoeffs.h"
#include "utilities.h"
#include "math.h"
#include "arm_math.h"
#include "arm_const_structs.h"

#define SOURCE_FILE_NAME "stm32f7_loop_dma_FunSP.c"

#define PULSE_GAIN 0.75f
#define N 1024 //1024
#define N2 512
#define FS 22050.0f // AUDIO_FREQUENCY_22K
#define F0 110.0f   // 110 male - 220 female
//#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE

extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

#define LPCorder 22

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
    float32_t gNoiseBuf[MAX_PULSE_BUFFER] = {0.0f};
    float32_t colorPulse[MAX_PULSE_BUFFER] = {0.0f};
		
    float32_t tilt[N] = {0.0f};

    // * generate single prototype pulse 10% shorter than that of the fundamental period
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

        // compute the PSD's centre of gravity
        float32_t cog = 0.0f, voiced = 0.0f;
        cog = centreOfGravity(currPSD, N, FS);
        voiced = mapping(cog, 2500, 4200, 0.99);

        // Place pulses inside the buffer
        while (idxMainBuf < N2){
            lastPulseLen = getPulseLen(FS, F0); // dummy pulse length just for placement
            pAlpha = (idxMainBuf != 0) ? idxMainBuf/N2 : 1/N2 ; // pulse filtering coefficient (check for zero)

            //* 2.5. Generate a new filter based on the current pulse position
            for (uint32_t i = 0; i < N; i++){
                tempBuff[i].real = (pAlpha * currPSD[i] + (1 - pAlpha) * lastPSD[i]) * tilt[i]; // apply tilt gain
                tempBuff[i].imag = 0.0f;

                lastPSD[i] = currPSD[i]; // save last psd
            }

            // Use weiner khinchine (apply IODFT) to compute the autocorrelation sequence from the PSD
            arm_cfft_f32(&arm_cfft_sR_f32_len1024, (float32_t *)(tempBuff), 1, 1); // ifft flag = 1, doBitReverse = 1
            arm_cmplx_mult_cmplx_f32((float32_t *)(tempBuff), (float32_t *)(invExp), (float32_t *)(currBuff), N); // iodft

            // only first (LPCorder + 1) values are needed for levinson recursion
            for (uint32_t i = 0; i <= LPCorder; i++)
                autocorrelation[i] = tempBuff[i].real;    // assuming here the complex part will be zero

            levinson(autocorrelation, LPCorder, aCoeff, &error); // LPC analysis
            b0 = sqrtf(error); // gain factor

            //* Add noise to the pulse before filtering
            float32_t uvMix = 0.5f + voiced * 0.5f;

            for (uint32_t i = 0; i < MAX_PULSE_BUFFER; i++)
                gNoiseBuf[i] = uvMix * gPulseBuf[i] + (1.0f - uvMix) * randFloat(); // add noise to the pulse
 
            //* 3. Filter and accumulate the pulse into the mainbuffer
            filteredLength = MAX_PULSE_BUFFER;
            iirDirectTypeII(gPulseBuf, colorPulse, MAX_PULSE_BUFFER, b0, aCoeff, LPCorder + 1);

            //* 4. Copy filtered pulse to main buffer
            for (int32_t i = 0; i < filteredLength; i++)
                mainBuf[idxMainBuf++] += 0.00002f * colorPulse[i];

            // update main buffer index
            idxMainBuf -= (filteredLength - lastPulseLen);
        }

        for (uint32_t i = 0; i < N; i++)
            lastPSD[i] = currPSD[i];
        
        idxMainBuf -= N2;
    }
}
