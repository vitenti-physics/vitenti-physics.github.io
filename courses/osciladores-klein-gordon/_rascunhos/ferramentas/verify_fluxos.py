import sympy as sp
import numpy as np
print("=== chain: discrete energy balance ===")
K, K0, m = sp.symbols("K K_0 m", positive=True)
hm, h0, hp, hpp, hmm = sp.symbols("h_m h_0 h_p h_pp h_mm")
vm, v0, vp, vpp, vmm = sp.symbols("v_m v_0 v_p v_pp v_mm")
acc0 = (K * (hp - 2 * h0 + hm) - K0 * h0) / m
e0 = m * v0**2 / 2 + K0 * h0**2 / 2 + K / 4 * ((hp - h0) ** 2 + (h0 - hm) ** 2)
de0 = sp.diff(e0, h0) * v0 + sp.diff(e0, hp) * vp + sp.diff(e0, hm) * vm + sp.diff(e0, v0) * acc0
S = lambda ha, hb, va, vb: -K / 2 * (hb - ha) * (va + vb)
print("de_j/dt + S_{j+1/2} - S_{j-1/2} =", sp.simplify(de0 + S(h0, hp, v0, vp) - S(hm, h0, vm, v0)))
print("=== continuum: energy and momentum conservation ===")
t, x = sp.symbols("t x")
rho, tau, mu = sp.symbols("rho tau mu", positive=True)
phi = sp.Function("phi")(t, x)
Pi = rho * sp.diff(phi, t)
E = Pi**2 / (2 * rho) + tau * sp.diff(phi, x) ** 2 / 2 + rho * mu**2 * phi**2 / 2
Sf = -tau / rho * Pi * sp.diff(phi, x)
Pd = -Pi * sp.diff(phi, x)
Tf = Pi**2 / (2 * rho) + tau * sp.diff(phi, x) ** 2 / 2 - rho * mu**2 * phi**2 / 2
kg = {sp.diff(phi, t, 2): tau / rho * sp.diff(phi, x, 2) - mu**2 * phi}
print("dE/dt + dS/dx =", sp.simplify(sp.expand((sp.diff(E, t) + sp.diff(Sf, x)).subs(kg))))
print("dP/dt + dT/dx =", sp.simplify(sp.expand((sp.diff(Pd, t) + sp.diff(Tf, x)).subs(kg))))
c = sp.sqrt(tau / rho)
A, B = sp.symbols("A B", real=True)  # A = Pi, B = phi'
diff_ = c * (A**2 / (2 * rho) + tau * B**2 / 2) - sp.Abs(tau / rho * A * B)
print("c*E_kin+grad - |S| = (sqrt-form):", sp.simplify(c * (A**2 / (2 * rho) + tau * B**2 / 2) - tau / rho * A * B - c / 2 * (A / sp.sqrt(rho) - sp.sqrt(tau) * B) ** 2))
print("=== numerics (lattice N = 512, L = 40) ===")
cc, rr, L = 1.0, 1.0, 40.0
N = 512; dx = L / N; xj = np.arange(N) * dx
k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
def sol(phi0, Pi0, t, mu_):
    w = np.sqrt(mu_**2 + 4 * cc**2 / dx**2 * np.sin(k * dx / 2) ** 2)
    a0, b0 = np.fft.fft(phi0), np.fft.fft(Pi0)
    t = np.atleast_1d(t)[:, None]
    ph = np.fft.ifft(a0 * np.cos(w * t) + b0 / rr * t * np.sinc(w * t / np.pi), axis=1).real
    Pi = np.fft.ifft(b0 * np.cos(w * t) - rr * w * a0 * np.sin(w * t), axis=1).real
    return ph, Pi
Kc, mm = rr * cc**2 / dx, rr * dx
def e_S(ph, Pi, mu_):
    v = Pi / rr                     # hdot = p/m = Pi dx/(rho dx)
    K0c = rr * mu_**2 * dx
    dh = np.roll(ph, -1, axis=-1) - ph
    e = mm * v**2 / 2 + K0c * ph**2 / 2 + Kc / 4 * (dh**2 + np.roll(dh, 1, axis=-1) ** 2)
    S = -Kc / 2 * dh * (v + np.roll(v, -1, axis=-1))   # S_{j+1/2}
    return e, S
x0, s = L / 2, 0.5
g = np.exp(-((xj - x0) ** 2) / (2 * s**2))
for mu_ in [0.0, 1.0]:
    tt = np.linspace(0, 12, 4801)
    ph, Pi = sol(g, 0 * g, tt, mu_)
    e, S = e_S(ph, Pi, mu_)
    jA, jB = np.searchsorted(xj, x0 - 5), np.searchsorted(xj, x0 + 5)   # window j in [jA, jB)
    Ew = e[:, jA:jB].sum(1)
    flux_out = S[:, jB - 1] - S[:, jA - 1]                                # S_{jB-1/2} - S_{jA-1/2}
    integ = np.concatenate([[0], np.cumsum((flux_out[1:] + flux_out[:-1]) / 2 * np.diff(tt))])
    print(f"mu={mu_}: E_window(0)={Ew[0]:.6f}, E_window(12)={Ew[-1]:.6f}, max|E_w(t)-E_w(0)+int flux| = {np.abs(Ew - Ew[0] + integ).max():.2e}, total E drift {np.abs(e.sum(1)/e.sum(1)[0]-1).max():.1e}")
    dens = e / dx
    m_ = dens > 1e-3 * dens.max()
    Smid = (S + np.roll(S, 1, axis=-1)) / 2
    vE = Smid[m_] / dens[m_]
    print(f"   local energy velocity S/E (where E > 1e-3 max): min {vE.min():.4f}, max {vE.max():.4f}")
# packet: energy velocity and momentum, continuum spectral on same grid
mu_ = 1.0; k0, sp_, xp = 2.0, 2.0, 10.0
w = np.sqrt(mu_**2 + 4 * cc**2 / dx**2 * np.sin(k * dx / 2) ** 2)
ca = np.where(k > 0, np.exp(-((k - k0) ** 2) * sp_**2 / 2) * np.exp(-1j * k * xp), 0)
tt = np.linspace(0, 20, 81)
psi = np.fft.ifft(ca[None] * np.exp(-1j * w[None] * tt[:, None]), axis=1)
dpsi = np.fft.ifft(ca[None] * (-1j * w[None]) * np.exp(-1j * w[None] * tt[:, None]), axis=1)
ph, Pi = psi.real, rr * dpsi.real
e, S = e_S(ph, Pi, mu_)
Etot = e.sum(1); Stot = S.sum(1)
print("packet: total E drift", np.abs(Etot / Etot[0] - 1).max(), "; energy velocity sum S / sum e:", (Stot / Etot)[[0, 40, 80]])
weights_E = np.abs(ca) ** 2 * 1.0          # energy weights ~ omega^2 |c|^2 / ... compute exactly below
vg = cc**2 * np.sin(k * dx) / (dx * w)
alpha_sq = np.abs(ca) ** 2 * w             # |alpha_a|^2 propto rho w |c_a|^2 /2 (since psi_+ = phi/2... ) see text
print("weighted mean vg, weights w^2|c|^2 (energy):", (w**2 * np.abs(ca) ** 2 * vg).sum() / (w**2 * np.abs(ca) ** 2).sum(), " weights |c|^2:", (np.abs(ca)**2 * vg).sum() / (np.abs(ca)**2).sum())
# lattice field momentum -sum Pi (phi_{j+1}-phi_{j-1})/2
Pf = -(Pi * (np.roll(ph, -1, axis=1) - np.roll(ph, 1, axis=1)) / 2).sum(1)
print("lattice field momentum: relative variation", np.abs(Pf / Pf[0] - 1).max(), "; P/E =", (Pf / (Etot / 1))[0] , " k0/omega(k0):", k0 / np.sqrt(mu_**2 + k0**2))
