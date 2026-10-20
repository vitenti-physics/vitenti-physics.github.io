import numpy as np
c, rho, L = 1.0, 1.0, 40.0
def numeros_de_onda(N, L):
    return 2 * np.pi * np.fft.fftfreq(N, d=L / N)          # k_a for a in numpy order
def omega_rede(k, dx, mu):
    return np.sqrt(mu**2 + 4 * c**2 / dx**2 * np.sin(k * dx / 2) ** 2)
def omega_cont(k, mu):
    return np.sqrt(mu**2 + c**2 * k**2)
def evolui(phi0, Pi0, t, w):
    ft0, pt0 = np.fft.fft(phi0), np.fft.fft(Pi0)
    t = np.asarray(t)[:, None]
    ft = ft0 * np.cos(w * t) + pt0 / rho * t * np.sinc(w * t / np.pi)
    return np.fft.ifft(ft, axis=1).real
x0, sig = L / 2, 0.5
print("=== dispersion relative error at k = 2, mu = 0 and 1 ===")
for dx in [1.0, 0.5, 0.25, 0.125]:
    for mu in [0.0, 1.0]:
        r = omega_rede(2.0, dx, mu) / omega_cont(2.0, mu) - 1
        print(f"dx={dx}: mu={mu}: rel err {r:.3e}" + (f"  (pred -(k dx)^2/24 = {-(2*dx)**2/24:.3e})" if mu == 0 else ""))
print("=== continuum mu=0 vs d'Alembert ===")
N = 4096; xs = np.arange(N) * L / N; k = numeros_de_onda(N, L)
g = lambda y: np.exp(-((y - x0) ** 2) / (2 * sig**2))
ts = np.array([0.0, 5.0, 10.0, 15.0])
ref = {mu: evolui(g(xs), 0 * xs, ts, omega_cont(k, mu)) for mu in [0.0, 1.0]}
dal = 0.5 * (g(xs[None, :] - c * ts[:, None]) + g(xs[None, :] + c * ts[:, None]))
print("max |spectral - dAlembert| =", np.abs(ref[0.0] - dal).max())
print("=== lattice -> continuum convergence at t = 15 ===")
for mu in [0.0, 1.0]:
    errs = []
    for Nn in [128, 256, 512, 1024]:
        xn = np.arange(Nn) * L / Nn; kn = numeros_de_onda(Nn, L); dx = L / Nn
        sol = evolui(g(xn), 0 * xn, ts, omega_rede(kn, dx, mu))
        errs.append(np.abs(sol[-1] - ref[mu][-1, :: N // Nn]).max())
    print(f"mu={mu}: errors {['%.2e' % e for e in errs]}, ratios {[round(errs[i]/errs[i+1], 2) for i in range(3)]}")
print("=== outside the light cone, N = 512, t = 15 ===")
Nn = 512; xn = np.arange(Nn) * L / Nn; kn = numeros_de_onda(Nn, L); dx = L / Nn
tt = np.linspace(0, 15, 301)
for mu in [0.0, 1.0]:
    sol = evolui(g(xn), 0 * xn, tt, omega_rede(kn, dx, mu))
    fora = np.abs(xn[None, :] - x0) > c * tt[:, None] + 5 * sig
    sol_c = evolui(g(xs), 0 * xs, tt, omega_cont(k, mu))
    fora_c = np.abs(xs[None, :] - x0) > c * tt[:, None] + 5 * sig
    print(f"mu={mu}: max |phi| outside |x-x0| > ct + 5 sigma: lattice {np.abs(sol[fora]).max():.1e}, continuum {np.abs(sol_c[fora_c]).max():.1e}; max|phi| inside at t=15: {np.abs(sol[-1]).max():.3f}; at x0, t=15: {sol[-1, Nn//2]:.3f}")
print("=== moving packet, mu = 1, k0 = 2 ===")
mu, k0, sp_, xp = 1.0, 2.0, 2.0, 10.0
Nn = 512; xn = np.arange(Nn) * L / Nn; kn = numeros_de_onda(Nn, L); dx = L / Nn
for nome, w in [("rede", omega_rede(kn, dx, mu)), ("contínuo", omega_cont(kn, mu))]:
    ca = np.where(kn > 0, np.exp(-((kn - k0) ** 2) * sp_**2 / 2) * np.exp(-1j * kn * xp), 0)
    tt = np.linspace(0, 20, 81)
    psi = np.fft.ifft(ca[None, :] * np.exp(-1j * w[None, :] * tt[:, None]), axis=1)
    dens = np.abs(psi) ** 2
    centro = (dens * xn).sum(1) / dens.sum(1)
    vel = np.polyfit(tt, centro, 1)[0]
    print(f"{nome}: fitted center velocity {vel:.4f}; neg-k weight {np.abs(np.exp(-((-np.abs(kn[kn<=0]) - k0) ** 2) * sp_**2 / 2)).max():.1e}")
print("group velocity continuum c^2 k0/omega:", c**2 * k0 / np.sqrt(mu**2 + c**2 * k0**2), "; lattice:", c**2 * np.sin(k0 * dx) / (dx * omega_rede(k0, dx, mu)), "; phase velocity omega/k:", np.sqrt(mu**2 + k0**2) / k0)
print("dx for N=512:", L / 512, " sigma/dx:", 0.5 / (L / 512), " kmax=pi/dx:", np.pi / (L / 512))
