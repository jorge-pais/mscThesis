// stm32f7_template.c

#include "stm32f7_wm8994_init.h"
#include "stm32f7_display.h"
#include "lcd_log.h"
#include <math.h>

#define SOURCE_FILE_NAME "stm32f7_DMA_template.c"

#define FS AUDIO_FREQUENCY_22K
//#define PI 3.141592653589793
#define MAX_PULSE_BUFFER PING_PONG_BUFFER_SIZE //256
#define N_PULSE_BUFFER 3
#define PULSE_GAIN 0.25

#define N 1024 // DFT size
#define N2 512

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

// Hard-coded not parameterized F0 contour
#define CONTOUR_SIZE 162
const float f0contour[] = {100.0, 108.025, 113.55, 117.337, 119.806, 121.253, 121.879, 121.833, 121.236, 120.186, 118.767, 117.059, 115.129, 113.046, 110.865, 108.658, 106.461, 104.34, 102.341, 100.506, 98.87, 97.459, 96.297, 95.392, 94.748, 94.356, 94.201, 94.261, 94.51, 94.919, 95.457, 96.092, 96.8, 97.552, 98.324, 99.098, 99.858, 100.588, 101.275, 101.914, 102.495, 103.015, 103.471, 103.861, 104.185, 104.442, 104.635, 104.767, 104.838, 104.855, 104.819, 104.735, 104.607, 104.44, 104.238, 104.007, 103.749, 103.472, 103.177, 102.87, 102.556, 102.238, 101.919, 101.604, 101.295, 100.995, 100.708, 100.434, 100.175, 99.934, 99.71, 99.507, 99.322, 99.158, 99.013, 98.888, 98.781, 98.693, 98.622, 98.567, 98.527, 98.502, 98.489, 98.487, 98.495, 98.513, 98.537, 98.569, 98.605, 98.645, 98.689, 98.735, 98.781, 98.828, 98.875, 98.921, 98.966, 99.009, 99.05, 99.088, 99.123, 99.155, 99.184, 99.21, 99.232, 99.252, 99.268, 99.281, 99.292, 99.299, 99.304, 99.306, 99.306, 99.304, 99.3, 99.294, 99.287, 99.278, 99.268, 99.257, 99.245, 99.232, 99.219, 99.205, 99.191, 99.177, 99.163, 99.149, 99.135, 99.121, 99.107, 99.094, 99.08, 99.068, 99.055, 99.043, 99.031, 99.02, 99.009, 98.999, 98.989, 98.979, 98.97, 98.961, 98.952, 98.943, 98.935, 98.927, 98.92, 98.912, 98.905, 98.898, 98.891, 98.884, 98.877, 98.87, 98.863, 98.857, 98.85, 98.843, 98.837, 98.83};

/// @brief Generate a single LF pulse and store it in the output buffer out
/// @param fs Sample rate
/// @param f0 Fundamental Frequency
/// @param pulseBuff Which pulse buffer to write to
/// @return Size of output vector (actual pulse len)
int genPulse(float32_t fs, float32_t f0, float32_t * pulseBuff, float32_t gain){
    float32_t T0 = 1/f0;
    float32_t Ts = 1/fs;

    int outLen = floor(T0/Ts + 0.5); // Round this up (?) i think C might implicitily floor the result when casting from an floating point to an integer type
    float32_t max = 0, sample;

    /// TODO: change nested loop order (?)
    for (int i = 0; i < outLen; i++){ // per sample calculation
        sample = 0.0;

        for (int j = 0; j < Lmag; j++) // per harmonic calculation
            sample += mag[j] * sin((j+1)*2*PI*f0*(i*Ts) + 2*PI*nrd[j]);
			
				if (fabs(sample) > max) // Update maximum absolute value
						max = fabs(sample);
				
				pulseBuff[i] = sample;
    }
    // I'd like to believe this wound't be necessary with compiler optimizations
    float32_t scale = gain * (INT16_MAX) / max;
		
    // set remaining buffer to 0 and normalize the samples within outLen
    /*for (uint16_t i = outLen; i < MAX_PULSE_BUFFER; i++)
        pulseBuff[i] = 0;*/
		
    for (uint16_t i = 0; i < outLen; i++)
        pulseBuff[i] = scale * pulseBuff[i];
    
    return outLen;
}

// DMA buffer processing variables
extern volatile int32_t TX_buffer_empty;       // these may not need to be int32_t
extern volatile int32_t RX_buffer_full;        // they were extern volatile int16_t in F4 version
extern int16_t rx_buffer_proc, tx_buffer_proc; // will be assigned token values PING or PONG

// Pulse position and length variables (for last pulse)
int lastPulseLen = 0;
int lastPulsePos = 0;
// What buffer was last read to
int currBuf = 0, lastBuf = 0;
int currF0 = 0;

float32_t old_outBuffer[PING_PONG_BUFFER_SIZE];

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

    for (int i = 0; i < PING_PONG_BUFFER_SIZE; i++)
    {
        *tx_buf++ = old_outBuffer[i];
        *tx_buf++ = 0;
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
				
        // then run the DMA processing routine
        process_buffer();
        
        float32_t sample;
        int overlap = 0, i = 0;
        while (lastPulsePos < lastPulseLen){
            old_outBuffer[i++] = pulse[0][lastPulsePos++];
            overlap++;
        }
        lastPulseLen = 0; //reset

        /// 2. Generate buffers and copy them
        for (i = overlap; i < PING_PONG_BUFFER_SIZE; i++){
            if(lastPulsePos >= lastPulseLen){
                lastPulseLen = genPulse(FS, 80, pulse[0], PULSE_GAIN);
                lastPulsePos = 0;
            }

            old_outBuffer[i] = pulse[0][lastPulsePos++];
        }

        /* //////// Generate the output buffer
        float32_t sample;

        for (int i = 0; i < PING_PONG_BUFFER_SIZE; i++){
            if(lastPulsePos < lastPulseLen){
                sample = pulse[0][lastPulsePos++];
            }
            else{
                lastPulseLen = genPulse(FS, 200, pulse[0], PULSE_GAIN);
                lastPulsePos = 0;
								if(currF0 < CONTOUR_SIZE - 1)
                    currF0++;

                lastBuf = currBuf;
                currBuf = (currBuf + 1) % N_PULSE_BUFFER;
                sample = pulse[0][lastPulsePos++];
            }
            old_outBuffer[i] = sample;
        } */
    }
}
