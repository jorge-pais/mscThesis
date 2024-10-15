// stm32f7_template.c

#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "lcd_log.h"
#include <math.h>

#define SOURCE_FILE_NAME "stm32f7_DMA_template.c"

#define FS AUDIO_FREQUENCY_22K
//#define PI 3.141592653589793
#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE //256
#define N_PULSE_BUFFER 1
#define PULSE_GAIN 0.75

/// @brief Individual pulse buffers
/// in this state it takes up 12k of memory
static float32_t pulse[N_PULSE_BUFFER][MAX_PULSE_BUFFER];

/// TODO: Pass these const arrays to another header file
// Hard coded glottal pulse parameters, magnitude and nrd
// (maybe precision got messed up copying from matlab :/ )
#define Lmag 20
const float32_t mag[] = { 
    2442622.5627, 1986197.9268, 1076016.1427, 697031.8704,
    494252.9668,  346305.8243, 274035.1474,  211492.4250,
    169194.7322,  143811.7299, 116065.7454,  101804.5773,
    87080.0532,   74941.8874,  68346.0143,   58389.1305,
    53970.0127,   48463.5746,  43285.9485,   40796.0093};
const float32_t nrd[] = { 
    0.0000, 0.4607, 0.8146, 0.1582, 
    0.4829, 0.8155, 0.1496, 0.4756, 
    0.8143, 0.1435, 0.4772, 0.8146,
    0.1429, 0.4830, 0.8146, 0.1487, 
    0.4877, 0.8169, 0.1575, 0.4901};

// 641 * 32bit => 2564 bytes
// we could change this to uint8 => 641 bytes 
#define contourSize 641

/// @brief Generate a single LF pulse and store it in the output buffer out
/// @param fs Sample rate
/// @param f0 Fundamental Frequency
/// @param pulseBuff pointer to buffer to write to
/// @return Size of output vector (actual pulse len)
int genPulse(float32_t fs, float32_t f0, float32_t * pulseBuff, float32_t gain){
    float32_t T0 = 1/f0;
    float32_t Ts = 1/fs;

    int outLen = floor(T0/Ts + 0.5);

    float32_t max = 0, float32_t sample;

    /// TODO: change nested loop order (?)
    for (uint16_t i = 0; i < outLen; i++){ // per sample calculation
        sample = 0;

        for (uint8_t j = 0; j < Lmag; j++) // per harmonic calculation
            sample += mag[j] * sin((j+1)*2*PI*f0*(i*Ts) + 2*PI*nrd[j]);

        if (fabs(sample) > max) // Update maximum absolute value
            max = fabs(sample);

        pulseBuff[i] = sample;
    }

    // I'd like to believe this wound't be necessary with compiler optimizations
    float32_t scale = gain * (INT16_MAX) / max;
		
    // set remaining buffer to 0 and normalize the samples within outLen
    for (uint16_t i = outLen; i < MAX_PULSE_BUFFER; i++)
        pulseBuff[i] = 0;
		
    for (uint16_t i = 0; i < outLen; i++)
        pulseBuff[i] = scale * pulseBuff[i];
    
    // in this case we just return the raw pulse length, this way we can know 
    return outLen;
}

// DMA buffer processing variables
extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

// Pulse position and length variables (for last pulse)
volatile int lastPulseLen = 0;
volatile int lastPulsePos = 0;
// What buffer was last read to
volatile int8_t currBuf = 0, lastBuf = 0;
volatile int8_t currF0 = 0;

/// @brief This function processes one DMA transfer block worth of data
/// @param void 
void process_buffer(void){

    //int16_t *rx_buf, *tx_buf; 
    int16_t *tx_buf; // in this example we're not gonna be reading anything for now

    // Select input and output buffers
    /* if (rx_buffer_proc == PING) {
        rx_buf = (int16_t *)PING_IN;
    }
    else {
        rx_buf = (int16_t *)PONG_IN;
    } */
    
    if (tx_buffer_proc == PING) 
        tx_buf = (int16_t *)PING_OUT;
    else 
        tx_buf = (int16_t *)PONG_OUT;

    //////// Generate the output buffer
    float sample;

    for (size_t i = 0; i < PING_PONG_BUFFER_SIZE; i++){
        if(lastPulsePos < lastPulseLen){ // Copy the previous buffer, if there is one yet to copy
            sample = pulse[1][lastPulsePos++];
        }
        else{
            lastPulsePos = 0;
			lastPulseLen = genPulse(FS, 210, pulse[0], PULSE_GAIN);
            currF0 = (currF0 + 1) % contourSize;
            lastBuf = currBuf;
            currBuf = (currBuf + 1) % N_PULSE_BUFFER;
            sample = 0.0f;
        }

        *tx_buf++ = sample;
        *tx_buf++ = 0.0f; // silent channel
    }

    RX_buffer_full = 0;
    TX_buffer_empty = 0;
}

int main(void){
    // first setup the audio codec and enable DMA processing 
    // todo: CHECK BUFFER SIZE  
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
		
    for(;;){
        // First wait for the DMA operation to fill the RX buffer
        // and empty out the TX_buffer
        while (!(RX_buffer_full && TX_buffer_empty)) {}
				
        // then run the processing routine
        process_buffer();
    }
}
