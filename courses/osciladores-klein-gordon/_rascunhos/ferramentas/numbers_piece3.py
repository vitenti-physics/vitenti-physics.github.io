import numpy as np
def passo_verlet(F, x, p, dt, m):
    p = p + 0.5 * dt * F(x); x = x + dt * p / m; p = p + 0.5 * dt * F(x); return x, p
def forca_cadeia(h, K, K0):
    return K * (np.roll(h, -1) - 2 * h + np.roll(h, 1)) - K0 * h
def energia_cadeia(h, p, m, K, K0):
    return (np.sum(p**2, axis=-1) / (2 * m) + 0.5 * K * np.sum((np.roll(h, -1, axis=-1) - h) ** 2, axis=-1) + 0.5 * K0 * np.sum(h**2, axis=-1))
def simula_cadeia(h0, p0, m, K, K0, dt, n_passos, a_cada=1):
    F = lambda h: forca_cadeia(h, K, K0)
    h, p = h0.copy(), p0.copy(); hs, ps = [h.copy()], [p.copy()]
    for n in range(1, n_passos + 1):
        h, p = passo_verlet(F, h, p, dt, m)
        if n % a_cada == 0: hs.append(h.copy()); ps.append(p.copy())
    return dt * a_cada * np.arange(len(hs)), np.array(hs), np.array(ps)
def frequencias(N, m, K, K0):
    a = np.arange(N)
    return np.sqrt(K0 / m + 4 * K / m * np.sin(np.pi * a / N) ** 2)
def solucao_exata(h0, p0, t, m, K, K0):
    N = len(h0); w = frequencias(N, m, K, K0)
    ht0 = np.fft.fft(h0, norm="ortho"); pt0 = np.fft.fft(p0, norm="ortho")
    t = np.asarray(t)[:, None]
    ht = ht0 * np.cos(w * t) + pt0 / m * t * np.sinc(w * t / np.pi)
    pt = pt0 * np.cos(w * t) - m * w * ht0 * np.sin(w * t)
    return np.fft.ifft(ht, axis=1, norm="ortho").real, np.fft.ifft(pt, axis=1, norm="ortho").real
N, m, K = 15, 1.0, 1.0
h0 = np.zeros(N); p0 = np.zeros(N); h0[0] = 1.0; p0[1] = p0[2] = 0.5
for K0 in [0.0, 0.2]:
    t, hs, ps = simula_cadeia(h0, p0, m, K, K0, 0.05, 1200, 5)
    he, pe = solucao_exata(h0, p0, t, m, K, K0)
    err = np.abs(hs - he).max(axis=1)
    E_ex = energia_cadeia(he, pe, m, K, K0)
    print(f"K0={K0}: Verlet vs exact: max err t<=10 {err[t<=10].max():.2e}, t<=30 {err[t<=30].max():.2e}, t<=60 {err.max():.2e}; exact energy drift {np.abs(E_ex/E_ex[0]-1).max():.1e}; P exact {pe.sum(1)[[0,-1]]}")
    # Parseval
    ht = np.fft.fft(h0, norm="ortho"); pt = np.fft.fft(p0, norm="ortho"); w = frequencias(N, m, K, K0)
    Hm = np.sum(np.abs(pt)**2 / (2*m) + m * w**2 * np.abs(ht)**2 / 2)
    print(f"   H direct {energia_cadeia(h0, p0, m, K, K0):.12f}  H modes {Hm:.12f}")
    errs = []
    for dt in [0.1, 0.05, 0.025]:
        n = int(round(20 / dt))
        t, hs, ps = simula_cadeia(h0, p0, m, K, K0, dt, n)
        he, pe = solucao_exata(h0, p0, t, m, K, K0)
        errs.append(np.abs(hs - he).max())
    print(f"   convergence t<=20: {['%.3e' % e for e in errs]}, ratios {errs[0]/errs[1]:.2f} {errs[1]/errs[2]:.2f}")
    # mode amplitudes of initial condition
    print("   |ht0| by a:", np.round(np.abs(np.fft.fft(h0, norm='ortho')), 3))
# CFL demo
K0 = 0.2
wmax = frequencias(N, m, K, K0).max(); dtc = 2 / wmax
print(f"omega_max {wmax:.4f}, dt_c {dtc:.4f}")
for f in [0.98, 1.02]:
    dt = f * dtc
    t, hs, ps = simula_cadeia(h0, p0, m, K, K0, dt, 300)
    E = energia_cadeia(hs, ps, m, K, K0)
    print(f"factor {f}: dt={dt:.4f}, t_end={t[-1]:.1f}, E_end/E0={E[-1]/E[0]:.3e}, max E/E0={E.max()/E[0]:.3e}")
# unstable growth rate predicted: Verlet eigenvalue for eps = omega_max dt
for f in [1.02]:
    eps = wmax * f * dtc
    tr = 2 - eps**2
    lam = (-tr + np.sqrt(tr**2 - 4)) / 2
    print(f"predicted |lambda| per step at factor {f}: {lam:.4f}; E growth per step lambda^2 = {lam**2:.4f}; over 300 steps ~ {lam**600:.2e} (amplitude^2)")
