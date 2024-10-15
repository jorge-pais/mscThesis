import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the magnitudes (in dB) of the Liljencrants-Fant model harmonics
LFmag = loadmat('../other/DAFx18_AJF_JMT/LFmag.mat')['LFmag'].flatten()
LFNRD = loadmat('../other/DAFx18_AJF_JMT/LFNRD.mat')['LFNRD'].flatten()

# Create a linear model for the NRD coefficients
P = np.polyfit(np.arange(2, 21), LFNRD[1:20], 1)
tmp_modelo = np.polyval(P, np.arange(2, 121))
shift_NRD = LFNRD[19] - tmp_modelo[18]
tmp_modelo += shift_NRD
tmp_modelo[:19] = LFNRD[:19]

# Make reference to ZERO
LFmag -= LFmag[0]
eixox = np.arange(1, 121)

logeixox = np.log10(eixox)
P2 = np.polyfit(logeixox[1:47], LFmag[1:47], 1)
modelo2 = np.polyval(P2, logeixox[:120])
offset = LFmag[1] - np.polyval(P2, logeixox[1])
modelo2 += offset
modelo2[:2] = LFmag[:2]

# Fit a line for the first two values
P1 = np.polyfit(logeixox[:2], LFmag[:2], 1)
modelo_IP = np.polyval(P1, logeixox[:2])

offset = np.polyval(P1, logeixox[1]) - np.polyval(P2, logeixox[1])
P2[1] += offset
modelo1 = np.polyval(P2, logeixox[1:47])

m1_IP, b1_IP = P1
m2_IP, b2_IP = P2

# Calculate the spectral tilt gain
N = 1024
N2 = N // 2
Fs = 22050

f0 = 250
currperiodIP = Fs / f0

ganho = np.zeros(N2)
for kk in range(N2):
    ratio = kk * currperiodIP / N
    if ratio < 1.0:
        ganho[kk] = 1.0
    elif ratio < 2.0:
        ganho[kk] = (ratio ** (m1_IP / 20)) * 10 ** (b1_IP / 20)
    else:
        ganho[kk] = (ratio ** (m2_IP / 20)) * 10 ** (b2_IP / 20)

# Plot for correctness check
plt.figure()
ratio = np.arange(N2) * currperiodIP / N
indice = ratio <= 50
plt.semilogx(np.arange(len(LFmag)) + 1, LFmag, label='LF Spectral Tilt')
plt.semilogx(ratio[indice], 20 * np.log10(ganho[indice]), '--', label='Model Approximation')
plt.legend()
plt.xlabel('Frequency/F0')
plt.ylabel('Gain (dB)')
plt.show()

# Load the spectral template
sofia_a_hlpc = loadmat('../other/templatesAJF/sofia_a_hlpc.mat')['sumhLPC'].flatten()
currIP_spec_model = sofia_a_hlpc

prototype = currIP_spec_model - 20 * np.log10(ganho)
plt.figure()
plt.plot(np.arange(N2) * Fs / N / 1000.0, currIP_spec_model, 'k', label='Original PSD Model')
plt.plot(np.arange(N2) * Fs / N / 1000.0, prototype, 'b--', label='Tilt Compensated Model')
plt.ylabel('PSD (dB)')
plt.xlabel('Frequency (kHz)')
plt.legend()
plt.show()