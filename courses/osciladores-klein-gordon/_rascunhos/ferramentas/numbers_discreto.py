import numpy as np
c, rho, L = 1.0, 1.0, 40.0
def numeros_de_onda(N, L): return 2 * np.pi * np.fft.fftfreq(N, d=L / N)
def omega_cadeia(k, dx, mu): return np.sqrt(mu**2 + 4 * c**2 / dx**2 * np.sin(k * dx / 2) ** 2)
def evolui_com_p(h0, p0, t, dx, mu):
    k = numeros_de_onda(len(h0), L); w = omega_cadeia(k, dx, mu); m = rho * dx
    ft0, pt0 = np.fft.fft(h0), np.fft.fft(p0)
    t = np.asarray(t)[:, None]
    ft = ft0 * np.cos(w * t) + pt0 / m * t * np.sinc(w * t / np.pi)
    pt = pt0 * np.cos(w * t) - m * w * ft0 * np.sin(w * t)
    return np.fft.ifft(ft, axis=1).real, np.fft.ifft(pt, axis=1).real
def energia_e_fluxo(h, p, dx, mu):
    m, K, K0 = rho * dx, rho * c**2 / dx, rho * mu**2 * dx
    hn, pn = np.roll(h, -1, axis=-1), np.roll(p, -1, axis=-1)
    e_mola = K / 2 * (hn - h) ** 2 + (p**2 + pn**2) / (4 * m) + K0 * (h**2 + hn**2) / 4
    S = -K / 2 * (hn - h) * (p + pn) / m
    return e_mola, S
x0, sig = L / 2, 0.5
N = 512; dx = L / N; x = np.arange(N) * dx; xm = x + dx / 2
g = np.exp(-((x - x0) ** 2) / (2 * sig**2))
mu = 1.0
tt = np.linspace(0, 15, 301)
h, p = evolui_com_p(g, 0 * g, tt, dx, mu)
e, S = energia_e_fluxo(h, p, dx, mu)
dens = e / dx
msk = dens > 1e-3 * dens.max()
print("max |S| dx/(c e_mola) where dens > 1e-3 max:", np.abs(S[msk] * dx / (c * e[msk])).max(), " overall (e>1e-14):", np.abs(S[e > 1e-14] * dx / (c * e[e > 1e-14])).max())
i = -1
j = np.argmax(dens[i] * (xm > x0))
inside = np.abs(xm - x0) < c * tt[i] / 2
print(f"t=15: max density at x - x0 = {xm[j] - x0:.2f}, S dx/e there = {S[i, j] * dx / e[i, j]:.3f}; fraction within |x-x0|<ct/2: {e[i, inside].sum() / e[i].sum():.3f}")
# window check
tf = np.linspace(0, 12, 4801)
for mu_ in [0.0, 1.0]:
    h, p = evolui_com_p(g, 0 * g, tf, dx, mu_)
    m, K, K0 = rho * dx, rho * c**2 / dx, rho * mu_**2 * dx
    hn, hp = np.roll(h, -1, 1), np.roll(h, 1, 1)
    ej = p**2 / (2 * m) + K0 * h**2 / 2 + K / 4 * ((hn - h) ** 2 + (h - hp) ** 2)
    _, S = energia_e_fluxo(h, p, dx, mu_)
    jA, jB = np.searchsorted(x, x0 - 5), np.searchsorted(x, x0 + 5)
    Ew = ej[:, jA:jB].sum(1); out = S[:, jB - 1] - S[:, jA - 1]
    acc = np.concatenate([[0], np.cumsum((out[1:] + out[:-1]) / 2 * np.diff(tf))])
    print(f"mu={mu_}: window E {Ew[0]:.4f} -> {Ew[-1]:.4f}; max err {np.abs(Ew - Ew[0] + acc).max():.1e}")
# packet on the chain: energy velocity = c^2 P_rede / H
mu = 1.0; k0, sp_, xp = 2.0, 2.0, 10.0
k = numeros_de_onda(N, L); w = omega_cadeia(k, dx, mu); m = rho * dx
ca = np.where(k > 0, np.exp(-((k - k0) ** 2) * sp_**2 / 2) * np.exp(-1j * k * xp), 0)
tp = np.linspace(0, 20, 81)
psi = np.fft.ifft(ca[None] * np.exp(-1j * w[None] * tp[:, None]), axis=1)
dpsi = np.fft.ifft(ca[None] * (-1j * w[None]) * np.exp(-1j * w[None] * tp[:, None]), axis=1)
h, p = psi.real, m * dpsi.real
H = (p**2).sum(1) / (2 * m) + (rho * c**2 / dx / 2 * (np.roll(h, -1, 1) - h) ** 2).sum(1) + (rho * mu**2 * dx / 2 * h**2).sum(1)
Pr = -(p * (np.roll(h, -1, 1) - np.roll(h, 1, 1))).sum(1) / (2 * dx)
print("packet: c^2 P_rede/H:", (c**2 * Pr / H)[[0, -1]], "; P_rede drift", np.abs(Pr / Pr[0] - 1).max())
# spectrum figure on the chain
a = np.arange(-10, 11); kk = 2 * np.pi * a / L
kap = np.sin(kk * dx) / dx
print("max |kappa/k - 1| for 0<|a|<=10:", np.abs(kap[a != 0] / kk[a != 0] - 1).max(), "; max |omega_chain/omega_cont -1|:", np.abs(omega_cadeia(kk, dx, 1.0) / np.sqrt(1 + kk**2) - 1).max())
from math import comb
print("dim H_n for N=15:", [comb(n + 14, n) for n in range(5)], " N=512 n=2:", comb(513, 2))
