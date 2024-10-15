// stm32f7_loop_dma_FunSP.c
#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "utilities.h"

#define SOURCE_FILE_NAME "stm32f7_loop_dma_FunSP.c"

#define N 1024 // DFT SIZE
#define N2 512 // PING_PONG_BUFFER_SIZE
#define FS 22050.0f // AUDIO_FREQUENCY_22K

#define Nlpc 22

#define F0 258.0f   // 110 male - 220 female
#define PULSE_GAIN 0.75
#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE

extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG


// for overlap add and odft calculation
float32_t sinW[N];
complex_t dirExp[N];
complex_t invExp[N];

// input and output buffers for overlap add
float32_t oldInBuff[N2]  = {0.0f};
float32_t oldOutBuff[N2] = {0.0f};

// current input buffer so, the layout is as follows:
// [oldInBuff | new DMA block]
complex_t currBuff[N];

/// @brief Processes one DMA transfer block worth of data
/// @param mainBuf Pointer to the main buffer
void process_buffer(float32_t *mainBuf) {
    int i;
    int16_t *rx_buf, *tx_buf;

    if (rx_buffer_proc == PING)
        rx_buf = (int16_t *)PING_IN;
    else
        rx_buf = (int16_t *)PONG_IN;

    if (tx_buffer_proc == PING)
        tx_buf = (int16_t *)PING_OUT;
    else
        tx_buf = (int16_t *)PONG_OUT;

    int16_t sampleInputL, sampleInputR;

    for (i = 0; i < (PING_PONG_BUFFER_SIZE); i++){
        // read both input channels
        sampleInputL = *(rx_buf++); // read input
        sampleInputR = *(rx_buf++);

        // copy processed input to the output buffer
        *(tx_buf++) = (int16_t) floor((float32_t) 0.5 + mainBuf[i]);
        *(tx_buf++) = sampleInputR; // right channel remains untouched

        // copy old buff to first half of currBuff
        currBuff[i].real = (float32_t) oldInBuff[i] * sinW[i];
        currBuff[i].imag = 0.0f;

        // copy new buff to second half of currBuff
        currBuff[i + N2].real = (float32_t) sampleInputL * sinW[i + N2];
        currBuff[i + N2].imag = 0.0f;

        // update old buffer with current samples
        oldInBuff[i] = sampleInputL;
    }

    RX_buffer_full = 0;
    TX_buffer_empty = 0;
}

int main(void) {
    uint16_t idxMainBuf = 0, vibrato = 0, contourIndex = 0; // reset main buffer index
    uint16_t gPulseLen, lastPulseLen, filteredLength;

    float32_t mainBuf[N] = {0.0f};
    float32_t gPulseBuf[MAX_PULSE_BUFFER] = {0.0f}; // 220 ~ 22050/100
    float32_t colorPulse[MAX_PULSE_BUFFER] = {0.0f};
    
    stm32f7_wm8994_init(AUDIO_FREQUENCY_22K, // 22050 Hz
                        IO_METHOD_DMA,
                        INPUT_DEVICE_INPUT_LINE_1,
                        //INPUT_DEVICE_DIGITAL_MICROPHONE_2, // select microphone input (digital microphone doesn't do anything???)
                        OUTPUT_DEVICE_HEADPHONE,
                        WM8994_HP_OUT_ANALOG_GAIN_0DB,
                        WM8994_LINE_IN_GAIN_0DB,
                        WM8994_DMIC_GAIN_9DB,
                        SOURCE_FILE_NAME,
                        NOGRAPH);

    // generate single prototype pulse 10% shorter than that of the fundamental period
    gPulseLen = genPulse(FS, F0 * 1.1, gPulseBuf, PULSE_GAIN);
    //gPulseLen = genPulse(FS, F0, gPulseBuf, PULSE_GAIN);

    // precompute the sine window and the complex exponentials for odft
    for (uint16_t i = 0; i < N; i++){
        sinW[i] = (float32_t) sin(PI / N * (i + 0.5)); // sine window

        // odft exponential
        dirExp[i].real = (float32_t) cos(PI * i / N);
        dirExp[i].imag = (float32_t) -sin(PI * i / N);

        // iodft exponential
        invExp[i].real = (float32_t) cos(PI * i / N);
        invExp[i].imag = (float32_t) sin(PI * i / N);
    }

		complex_t X[N]; // where to store the odft and psd results
		complex_t psd[N];
		
    for(;;){
        while (!(RX_buffer_full && TX_buffer_empty)){}

        //* 1. Copy first half to the output during the DMA processing
        //* also read the input buffer and window it for the odft
        process_buffer(mainBuf);


        // First compute the odft for the current frame
        arm_cmplx_mult_cmplx_f32((float32_t *)(currBuff), (float32_t *)(dirExp), (float32_t *)(currBuff), N);
        arm_cfft_f32(&arm_cfft_sR_f32_len1024, (float32_t *)(currBuff), 0, 1); // ifftFlag = 0; doBitReverse = 1

        // Then get the autocorrelation sequence through the wiener khinchin theorem, then the LPC coefficients
        arm_cmplx_mag_squared_f32((float32_t*) currBuff, (float32_t*) currBuff, N);
        arm_cfft_f32(&arm_cfft_sR_f32_len1024, (float32_t*) currBuff, 1, 1); // ifftFlag = 1; doBitReverse = 1
        arm_cmplx_mult_cmplx_f32((float32_t*) (currBuff), (float32_t*) (invExp), (float32_t*)(currBuff), N); // IODFT

		// autocorrelation sequence
        float32_t r[Nlpc + 1];
        for (int i = 0; i < Nlpc; i++)
            r[i] = currBuff[i].real;
        
        // Compute the LPC coefficients
		float32_t currLPC[N_LPC + 1];
        float32_t b0;
        levinson(r, N_LPC, currLPC, &b0, NULL);
        b0 = sqrt(b0);

        //* 2. Copy second half to first half
        for (size_t i = 0; i < N2; i++){
            mainBuf[i] = mainBuf[i + N];
            mainBuf[i + N] = 0.0f;
        }
        
        //* Place pulses inside the buffer
        while (idxMainBuf < N2){
            
            //..... follow f0 contourIndex
            //..... we can save one division by instead of saving
            //..... the contour freq, saving the period
            //lastPulseLen = getPulseLen(FS, contour[contourIndex]);
            lastPulseLen = getPulseLen(FS, F0); // dummy pulse length just for placement

            //..... add vibrato
            //lastPulseLen += (uint16_t) floor((float32_t) 0.5f + (float32_t) 7.0f * sin((float32_t)2.0f * PI * vibrato / (float32_t)28.0f));

            //* 3. Filter and accumulate the pulse into the mainbuffer
            // for now just a dummy filtering operation as the pulse is filled with zeros at the end
            filteredLength = MAX_PULSE_BUFFER;
            iirDirectTypeII(gPulseBuf, colorPulse, MAX_PULSE_BUFFER, b0, currLPC, N_LPC);

            //* 4. Copy *filtered* pulse to main buffer
            for (size_t i = 0; i < filteredLength; i++)
                mainBuf[idxMainBuf++] += 0.7 * colorPulse[i];

            //* update main buffer index with the difference in size between the original and filtered pulse
            idxMainBuf -= (filteredLength - lastPulseLen);
            
            vibrato = (vibrato+1) % 28;
            //contourIndex = (contourIndex < contourLen) ? contourIndex + 1 : contourIndex;
        }

        idxMainBuf -= N;
    }
}
