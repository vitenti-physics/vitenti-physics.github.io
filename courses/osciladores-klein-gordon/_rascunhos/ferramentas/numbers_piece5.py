import numpy as np
c, rho, L, mu = 1.0, 1.0, 40.0, 1.0
def malha(N):
    dx = L / N; x = np.arange(N) * dx
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    w = np.sqrt(mu**2 + 4 * c**2 / dx**2 * np.sin(k * dx / 2) ** 2)
    return dx, x, k, w
def alfas(phi, Pi, dx, k, w):
    ht = np.fft.fft(phi) * dx / np.sqrt(L); pt = np.fft.fft(Pi) * dx / np.sqrt(L)
    return (rho * w * ht + 1j * pt) / np.sqrt(2 * rho * w)
def campo(coef, t, x, k, w):   # sum_a coef_a f_a(t, x)
    t = np.atleast_1d(t)
    return (np.exp(1j * np.outer(x, k))[None] * (coef * np.exp(-1j * np.outer(t, w)) / np.sqrt(2 * rho * w))[:, None, :]).sum(-1) / np.sqrt(L)
N = 512; dx, x, k, w = malha(N)
x0, s = L / 2, 0.5
phi = np.exp(-((x - x0) ** 2) / (2 * s**2))
al = alfas(phi, 0 * x, dx, k, w)
psi = campo(al, 0.0, x, k, w)[0]
print("max|2 Re psi_+ - phi| =", np.abs(2 * psi.real - phi).max(), "; norm <psi_+,psi_+> = sum|alpha|^2 =", np.sum(np.abs(al) ** 2), "; H = sum w|alpha|^2 =", np.sum(w * np.abs(al) ** 2))
for d in [1, 2, 3, 5, 8]:
    i = np.argmin(np.abs(x - (x0 + d)))
    print(f"|x-x0|={d}: phi={phi[i]:.2e}  |psi_+|={np.abs(psi[i]):.2e}  Im psi_+={psi[i].imag:.2e}")
print("max |Im psi_+| =", np.abs(psi.imag).max(), " at x-x0 =", x[np.argmax(np.abs(psi.imag))] - x0)
# two particles
N2 = 256; dx2, x2, k2, w2 = malha(N2)
def pacote(k0, xp, sp_=2.0):
    co = np.where(np.sign(k2) == np.sign(k0), np.exp(-((k2 - k0) ** 2) * sp_**2 / 2) * np.exp(-1j * k2 * xp), 0)
    return co / np.sqrt(np.sum(np.abs(co) ** 2))
cL, cR = pacote(2.0, 10.0), pacote(-2.0, 30.0)
print("<psi_L, psi_R> =", np.vdot(cL, cR))
Psi_ab = (np.outer(cL, cR) + np.outer(cR, cL)) / np.sqrt(2)
print("norm^2 of Psi =", np.sum(np.abs(Psi_ab) ** 2))
vg = (np.abs(cL) ** 2 * c**2 * np.sin(k2 * dx2) / (dx2 * w2)).sum()
print("mean group velocity of psi_L:", vg, " crossing time 10/vg =", 10 / vg)
for tt in [0.0, 10 / vg, 2 * 10 / vg]:
    pL, pR = campo(cL, tt, x2, k2, w2)[0], campo(cR, tt, x2, k2, w2)[0]
    Psi = (np.outer(pL, pR) + np.outer(pR, pL)) / np.sqrt(2)
    D = np.abs(Psi) ** 2
    i, j = np.unravel_index(np.argmax(D), D.shape)
    print(f"t={tt:.2f}: max|Psi - Psi^T| = {np.abs(Psi - Psi.T).max():.1e}, peak at (x1,x2)=({x2[i]:.2f},{x2[j]:.2f})")
