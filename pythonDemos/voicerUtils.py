import numpy as np
from scipy.io import loadmat


def odft(x, N):
    """
    Computes the Odd Discrete Fourier Transform (ODFT) of a given signal.

    Parameters:
    x (array-like): The input signal.
    N (int): The number of samples in the DFT.

    Returns:
    array-like: The ODFT of the input signal.
    """
    f = np.arange(N) # from 0 to N-1
    dirExp = np.exp(-1j * np.pi * f / N)

    X = np.fft.fft(x*dirExp)

    #plt.plot(abs(X)) # plot the DFT magnitude

    return X

def iodft(X, N):
    """
    Performs the inverse odd discrete Fourier transform (IODFT) on the input signal.

    Parameters:
    X (array-like): The input signal in the frequency domain.
    N (int): The length of the input signal.

    Returns:
    array-like: The output signal in the time domain after applying the IODFT.
    """
    x = np.fft.ifft(X)

    f = np.arange(N) # from 0 to N-1
    x = x * np.exp(1j * np.pi * f / N)

    return x

a = np.random.normal(0, 1, 1024)
assert iodft(odft(a, 1024), 1024).all() == a.all()

def levinson(r, N):
    """
    Apply the Levinson-Durbin recursion to compute the prediction filter coefficients.

    Parameters:
    r (array-like): Autocorrelation sequence of length N+1.
    N (int): Order of the prediction filter.

    Returns:
    tuple: A tuple containing the prediction filter coefficients (A), the prediction error (E),
           and the reflection coefficients (k).

    """
    k = np.zeros(N+1) # Reflection coefficients
    A = np.zeros(N+1) # Prediction filter coefficients
    A[0] = 1.0
    E = np.real(r[0]) # Prediction error

    # Levinson-Durbin recursion
    for i in range(1, N+1):
        lambda_ = -np.dot(r[1:i+1], A[i-1::-1]) / E
        A[1:i+1] += lambda_ * A[i-1::-1]
        E *= (1.0 - lambda_ * lambda_)
        k[i] = -lambda_
    return A, E, k[1:]

x = np.array([1.0, 0.7, 0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0])

a1 = levinson(x, 6)[0]
#print("a1 (levinson):", a1)

# matlab result
result = np.array([1.000000000000000, -0.158730158730140, -0.603174603174617, -0.158730158730125, -0.587301587301598, -0.317460317460331, 0.984126984126981])
assert np.allclose(a1, result)

FS = 22050

LFmag = np.array(
[127.75712729, 125.96045048, 120.63637574, 116.86505272, 113.8789857,
 110.78919591, 108.75612537, 106.50589634, 104.56773675, 103.15588621,
 101.2940813, 100.1553461, 98.79837372, 97.49449253, 96.69426387,
  95.32664017, 94.64305041, 93.70830887, 92.72693877, 92.21235364,
  91.19303829, 90.74026643, 90.05456809, 89.24164179, 89.00770799,
  88.16411508, 87.86905499, 87.36129279, 86.73759912, 86.62027454,
  85.91529439, 85.73638821, 85.37634018, 84.8546795, 84.86114343,
  84.25057655, 84.19294799, 83.92136407, 83.50342335, 83.60809649,
  83.10197466, 83.12072764, 82.93968198, 82.62918072, 82.81524698,
  82.38319663, 82.49187896, 82.40771467, 82.15257198, 82.4169382])
tmpmag = 10 ** (LFmag/20) 

plf = [0.335465, -0.207431]
nrd_approx = np.polyval(plf, np.arange(2, len(LFmag)))
nrd_approx = np.insert(nrd_approx, 0, 0.0)
nrd_approx = np.mod(nrd_approx, 1)

def singleLFpulse(F0, Fs, mag, nrd, gain = 0.95):
    """
    Generates a single low-frequency pulse using the given parameters.
    """
    T = 1/F0 # fundamental period
    t = np.arange(0, T, 1/Fs)
    
    L = len(nrd) # number of harmonics to use (how is nrd lower than mag in this case?)
    dgf = np.zeros(len(t)) # output glottal flow derivative

    for i in range(L):
        dgf = dgf + mag[i] * np.sin((i+1)*2*np.pi*F0*t + 2*np.pi*nrd[i])

    dgf = dgf/max(abs(dgf)) * gain # normalize output signal amplitude

    return dgf

# Load delay coefficients and magnitudes from matlab files
mat = loadmat('LFNRD.mat')
LFnrd = mat["LFNRD"][0]
mat = loadmat('LFmag.mat')
LFmag = mat["LFmag"][0]

# NRD linear model + adjustments
# linear model of the first 20 nrd 
P = np.polyfit(np.arange(1,19), LFnrd[1:19], 1)
tmp_modelo = np.polyval(P, np.arange(1,120)) # extrapolate to 120
# shift the model to match the 20th value
shift_NRD = LFnrd[19] - tmp_modelo[19]
tmp_modelo = tmp_modelo + shift_NRD
tmp_modelo[:18] = LFnrd[:18] 

LFmag = LFmag - LFmag[0]
eixox = np.arange(1, 121)
logeixox = np.log10(eixox)

P2 = np.polyfit(logeixox[1:48], LFmag[1:48], 1)
modelo2 = np.polyval(P2, logeixox[:120])
offset = LFmag[1] - np.polyval(P2, logeixox[1])
modelo2 = modelo2 + offset
modelo2[:2] = LFmag[:2]

P1 = np.polyfit(logeixox[:2], LFmag[:2], 1)
modelo_IP = np.polyval(P1, logeixox[:2])

offset = np.polyval(P1, logeixox[1]) - np.polyval(P2, logeixox[1])
P2[1] = P2[1] + offset
modelo1 = np.polyval(P2, logeixox[1:48])

# three rules exist for compensating the magnitude curve
#      f/F0 < 1         -> G = 1
#      1 <= f/F0 < 2    -> G = m1*(f/F0) + b1
#      f/F0 >= 2        -> G = m2*(f/F0) + b2

m1_IP = P1[0]; b1_IP = P1[1]
m2_IP = P2[0]; b2_IP = P2[1]

# we can then define a function that calculates the gain
# that should be applied to the magnitude curve at a specific frequency

def tiltGain(F0, Fs, N):
    """
    Calculates the gain that should be applied to the magnitude curve at a specific frequency.
    """
    T0 = Fs/F0
    N2 = N//2

    gain = np.zeros(N2)

    for kk in range(0,N2):
        ratio = kk * (Fs/F0) / N
        if ratio < 1.0:
            gain[kk] = 1.0
        elif ratio < 2.0:
            gain[kk] = (ratio**(m1_IP/20)) * 10**(b1_IP/20)
        else:
            gain[kk] = (ratio**(m2_IP/20)) * 10**(b2_IP/20)

    return gain


def v_glotlf(d=0, t=None, p=None):
    if p is None:
        p = np.array([0.6, 0.1, 0.2])
    elif len(p) < 3:
        default_p = np.array([0.6, 0.1, 0.2])
        p = np.concatenate([p, default_p[len(p):]])

    if t is None:
        tt = np.linspace(0, 0.99, 100)
    else:
        tt = t - np.floor(t)  # only interested in the fractional part of t

    te = p[0]
    mtc = te - 1
    e0 = 1
    wa = np.pi / (te * (1 - p[2]))
    a = -np.log(-p[1] * np.sin(wa * te)) / te
    inta = e0 * ((wa / np.tan(wa * te) - a) / p[1] + wa) / (a**2 + wa**2)

    rb0 = p[1] * inta
    rb = rb0
    for _ in range(4):
        kk = 1 - np.exp(mtc / rb)
        err = rb + mtc * (1 / kk - 1) - rb0
        derr = 1 - (1 - kk) * (mtc / rb / kk)**2
        rb = rb - err / derr

    e1 = 1 / (p[1] * (1 - np.exp(mtc / rb)))
    ta = tt < te
    tb = ~ta

    if d == 0:
        u = np.zeros_like(tt)
        u[ta] = e0 * (np.exp(a * tt[ta]) * (a * np.sin(wa * tt[ta]) - wa * np.cos(wa * tt[ta])) + wa) / (a**2 + wa**2)
        u[tb] = e1 * (np.exp(mtc / rb) * (tt[tb] - 1 - rb) + np.exp((te - tt[tb]) / rb) * rb)
    elif d == 1:
        u = np.zeros_like(tt)
        u[ta] = e0 * np.exp(a * tt[ta]) * np.sin(wa * tt[ta])
        u[tb] = e1 * (np.exp(mtc / rb) - np.exp((te - tt[tb]) / rb))
    elif d == 2:
        u = np.zeros_like(tt)
        u[ta] = e0 * np.exp(a * tt[ta]) * (a * np.sin(wa * tt[ta]) + wa * np.cos(wa * tt[ta]))
        u[tb] = e1 * np.exp((te - tt[tb]) / rb) / rb
    else:
        raise ValueError('Derivative must be 0, 1, or 2')

    return u