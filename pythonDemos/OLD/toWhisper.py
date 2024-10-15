#!/bin/python3

import numpy as np
from scipy.io import wavfile
from scipy.signal import lfilter
from numpy.random import default_rng
from librosa import lpc

GAIN = 10
DEFAULT_EMPHASIS = 0.6
INPUT_FILE = '/home/jorgep/OneDrive_JPais/Faculdade/Mestrado/Tese/data/4572-112383-0003.wav'
OUTPUT_FILE = 'output.wav'
1
def pre_emphasis(signal, coeff = DEFAULT_EMPHASIS):
    return np.append(signal[0], signal[1:] - coeff * signal[:-1])

def de_emphasis(signal, coeff = DEFAULT_EMPHASIS):
    return lfilter([1], [1, -coeff], signal)

def gen_white_noise(length):
    rng = default_rng()
    return rng.standard_normal(length)

def apply_window(signal, window_type='hamming'):
    if window_type == 'hamming':
        window = np.hamming(len(signal))
    elif window_type == 'hanning':
        window = np.hanning(len(signal))
    elif window_type == 'hamming':
        window = np.blackman(len(signal))
    else:
        window = np.ones(len(signal))
        
    return signal*window

def whisper_synthesis(input_file, output_file, lpc_order = 20, frame_length=10, sampling_rate = 44100):
    sr, signal = wavfile.read(input_file)
    if signal.ndim > 1: # Stereo to mono by averaging both channels (?)
        signal = signal.mean(axis=1)

    signal = pre_emphasis(signal)
    frames = int(np.ceil(len(signal) / (frame_length * 1e-3 * sr)))
    synthesized_signal = np.zeros_like(signal)

    for i in range(frames):
        start = int(i * frame_length * 1e-3 * sr)
        end = int((i + 1) * frame_length * 1e-3 * sr)
        frame = signal[start:end]
        frame = apply_window(frame)

        # LPC analysis
        a = lpc(frame, order=lpc_order)
        error_signal = lfilter(a, [1], frame) # inverse filtering
        
        # Synthesize white noise for excitation
        white_noise = gen_white_noise(len(frame))

        # Apply filter to noise
        whisper_frame = lfilter([1], a, white_noise)

        # Blend into synthesized signal
        synthesized_signal[start:end] += whisper_frame

    # de-emphasis filtering and amplitude normalization
    synthesized_signal = de_emphasis(synthesized_signal)
    synthesized_signal = GAIN * synthesized_signal

    wavfile.write(output_file, sr, synthesized_signal.astype(np.int16))


if __name__ == "__main__":
    whisper_synthesis(INPUT_FILE, OUTPUT_FILE)

    import matplotlib.pyplot as plt

    # Load the input and output signals
    input_sr, input_signal = wavfile.read(INPUT_FILE)
    output_sr, output_signal = wavfile.read(OUTPUT_FILE)

    input_signal = input_signal[:][1]

    # Compute the spectrogram for the input signal
    input_spec, input_freqs, input_times, _ = plt.specgram(input_signal, Fs=input_sr)

    # Compute the spectrogram for the output signal
    output_spec, output_freqs, output_times, _ = plt.specgram(output_signal, Fs=output_sr)

    # Display the spectrograms
    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.title('Input Signal Spectrogram')
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.imshow(input_spec, aspect='auto', origin='lower', extent=[input_times.min(), input_times.max(), input_freqs.min(), input_freqs.max()])
    plt.colorbar(label='Magnitude')  # Add this line to create a colorbar

    plt.subplot(2, 1, 2)
    plt.title('Output Signal Spectrogram')
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.imshow(output_spec, aspect='auto', origin='lower', extent=[output_times.min(), output_times.max(), output_freqs.min(), output_freqs.max()])
    plt.colorbar(label='Magnitude')  # Add this line to create a colorbar

    plt.tight_layout()
    plt.savefig('spectrograms.png')  # Save the figure as an image instead of displaying it