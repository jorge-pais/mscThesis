// stm32f7_loop_dma_FunSP.c
#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "utilities.h"

#define SOURCE_FILE_NAME "stm32f7_loop_dma_FunSP.c"

#define PULSE_GAIN 0.75
#define N PING_PONG_BUFFER_SIZE
#define FS 22050.0 // AUDIO_FREQUENCY_22K
#define F0 110.0   // 110 male - 220 female
//#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE

//#define getPulseLen(fs, f0) ((uint16_t) floor(1/f0 * fs + 0.5))

extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

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

    for (i = 0; i < (PING_PONG_BUFFER_SIZE); i++){
        // copy the input buffer to the output buffer
        *(tx_buf++) = (int16_t) floor((float32_t) 0.5 + mainBuf[i]);
        
        rx_buf++; // LEFT channel unused
        *(tx_buf++) = *(rx_buf++); // RIGHT channel remains untouched
    }

    RX_buffer_full = 0;
    TX_buffer_empty = 0;
}

int main(void) {

    uint16_t idxMainBuf = 0, vibrato = 0, contourIndex = 0; // reset main buffer index
    uint16_t gPulseLen, lastPulseLen, filteredLength;

    float32_t mainBuf[2 * PING_PONG_BUFFER_SIZE] = {0.0f};
    float32_t gPulseBuf[MAX_PULSE_BUFFER] = {0.0f}; // 220 ~ 22050/100

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

    // precompute complex exponential vectors for the ODFT
    

    for(;;){
        while (!(RX_buffer_full && TX_buffer_empty)){}

        // 1. Copy first half to the output during the DMA processing
        process_buffer(mainBuf); 

        // 2. Copy second half to first half
        for (size_t i = 0; i < N; i++){
            mainBuf[i] = mainBuf[i + N];
            mainBuf[i + N] = 0.0f;
        }

        // Place pulses inside the buffer
        while (idxMainBuf < N){
            // follow f0 contourIndex // TODO precompute the contour
            // we can save one division by instead of saving the contour freq, saving the period
            lastPulseLen = getPulseLen(FS, contour[contourIndex]); // dummy pulse length just for placement

            // add vibrato
            lastPulseLen += (uint16_t) floor((float32_t) 0.5f + (float32_t) 7.0f * sin((float32_t)2.0f * PI * vibrato / (float32_t)28.0f));

            // 3. Filter and accumulate the pulse into the mainbuffer
            // for now just a dummy filtering operation as the pulse is filled with zeros at the end
            filteredLength = lastPulseLen + 20;

            // 4. Copy *filtered* pulse to main buffer
            for (size_t i = 0; i < filteredLength; i++)
                mainBuf[idxMainBuf++] += gPulseBuf[i];

            // update main buffer index (this eliminates the need for a second condition?)
            idxMainBuf -= (filteredLength - lastPulseLen);
            
            vibrato = (vibrato+1) % 28;
            contourIndex = (contourIndex < contourLen) ? contourIndex + 1 : contourIndex;
        }

        idxMainBuf -= N;
    }
}
