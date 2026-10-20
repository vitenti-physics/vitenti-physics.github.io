import numpy as np
from scipy.linalg import sqrtm
rng = np.random.default_rng(3)
c, rho, L = 1.0, 1.0, 40.0
def cadeia(N, mu):
    dx = L / N; m = rho * dx; K = rho * c**2 / dx; K0 = rho * mu**2 * dx
    T = np.roll(np.eye(N), 1, axis=1); D = 2 * np.eye(N) - T - T.T
    a = np.fft.fftfreq(N, d=1 / N); k = 2 * np.pi * a / L
    w = np.sqrt(K0 / m + 4 * K / m * np.sin(np.pi * a / N) ** 2)
    V = np.exp(2j * np.pi * np.outer(np.arange(N), a) / N) / np.sqrt(N)   # column a = v_a
    return dict(N=N, dx=dx, m=m, K=K, K0=K0, D=D, T=T, k=k, w=w, V=V, S=K * D + K0 * np.eye(N))
def evolve(C, h0, p0, t):
    ht, pt = C["V"].conj().T @ h0, C["V"].conj().T @ p0
    w, m = C["w"], C["m"]
    t = np.atleast_1d(t)[:, None]
    H = ht * np.cos(w * t) + pt / m * t * np.sinc(w * t / np.pi)
    P = pt * np.cos(w * t) - m * w * ht * np.sin(w * t)
    return (H @ C["V"].T).real, (P @ C["V"].T).real
print("=== bond energy bound and lattice momentum ===")
for mu in [0.0, 1.0]:
    C = cadeia(512, mu)
    x = np.arange(C["N"]) * C["dx"]
    h0 = np.exp(-((x - 20) ** 2) / (2 * 0.5**2)); p0 = 0 * x
    tt = np.linspace(0, 15, 301)
    hs, ps = evolve(C, h0, p0, tt)
    v = ps / C["m"]; dh = np.roll(hs, -1, 1) - hs
    S = -C["K"] / 2 * dh * (v + np.roll(v, -1, 1))
    eb = C["K"] / 2 * dh**2 + (ps**2 + np.roll(ps, -1, 1) ** 2) / (4 * C["m"]) + C["K0"] * (hs**2 + np.roll(hs, -1, 1) ** 2) / 4
    H = eb.sum(1)
    E_direct = (ps**2).sum(1) / (2 * C["m"]) + 0.5 * np.einsum("ti,ij,tj->t", hs, C["S"], hs)
    msk = eb > 1e-12
    print(f"mu={mu}: sum e_bond == H: {np.abs(H - E_direct).max():.1e}; max |S| dx /(c e_bond) = {np.abs(S[msk] * C['dx'] / (c * eb[msk])).max():.4f}")
    Plat = -(ps * (np.roll(hs, -1, 1) - np.roll(hs, 1, 1))).sum(1) / (2 * C["dx"])
    print(f"   P_rede drift {np.abs(Plat - Plat[0]).max():.1e} (P_rede(0)={Plat[0]:.3e}); max|dx sum S - c^2 P_rede| = {np.abs(C['dx'] * S.sum(1) - c**2 * Plat).max():.1e}")
C = cadeia(64, 1.0)
h0 = rng.normal(size=64); p0 = rng.normal(size=64)
ht, pt = C["V"].conj().T @ h0, C["V"].conj().T @ p0
al = (C["m"] * C["w"] * ht + 1j * pt) / np.sqrt(2 * C["m"] * C["w"])
kap = np.sin(C["k"] * C["dx"]) / C["dx"]
Plat = -(p0 * (np.roll(h0, -1) - np.roll(h0, 1))).sum() / (2 * C["dx"])
print("P_rede vs sum kappa |alpha|^2:", Plat, (kap * np.abs(al) ** 2).sum(), "; H vs sum w|alpha|^2:", p0 @ p0 / (2 * C["m"]) + 0.5 * h0 @ C["S"] @ h0, (C["w"] * np.abs(al) ** 2).sum())
print("P_rede random state drift:", end=" ")
hs, ps = evolve(C, h0, p0, np.linspace(0, 30, 50))
PP = -(ps * (np.roll(hs, -1, 1) - np.roll(hs, 1, 1))).sum(1) / (2 * C["dx"])
print(np.abs(PP - PP[0]).max())
print("=== chain product and modes ===")
m, w, V = C["m"], C["w"], C["V"]
def modo_h(t): return V * np.exp(-1j * w * t) / np.sqrt(2 * m * w)       # f_a^j(t), column a
def prod(h1, hd1, h2, hd2): return 1j * m * (h1.conj().T @ hd2 - hd1.conj().T @ h2)
f = modo_h(2.1); fd = -1j * w * f
print("max|<f_a,f_b>-d|", np.abs(prod(f, fd, f, fd) - np.eye(64)).max(), " max|<f*,f*>+d|", np.abs(prod(f.conj(), fd.conj(), f.conj(), fd.conj()) + np.eye(64)).max(), " max|<f,f*>|", np.abs(prod(f, fd, f.conj(), fd.conj())).max())
f0 = modo_h(0.0); fd0 = -1j * w * f0
al2 = prod(f0, fd0, h0[:, None], (p0 / m)[:, None])[:, 0]
print("alpha = <f_a, h> vs formula:", np.abs(al2 - al).max())
print("=== one-particle Schroedinger with matrix sqrt ===")
Om = sqrtm(C["S"] / m)
psi_coef = rng.normal(size=64) + 1j * rng.normal(size=64)
psi = lambda t: modo_h(t) @ psi_coef
dt = 1e-6
lhs = 1j * (psi(1.0 + dt) - psi(1.0 - dt)) / (2 * dt)
print("max |i dpsi/dt - sqrt(S/m) psi| =", np.abs(lhs - Om @ psi(1.0)).max(), " (scale", np.abs(Om @ psi(1.0)).max(), ")")
print("=== quantum chain, N = 3, occupation basis ===")
hbar = 1.0
C3 = cadeia(3, 1.0)
nc = 7
a1 = np.diag(np.sqrt(np.arange(1, nc)), 1)
I1 = np.eye(nc)
def op(single, i):
    mats = [I1] * 3; mats[i] = single
    out = mats[0]
    for M in mats[1:]: out = np.kron(out, M)
    return out
A = [op(a1, i) for i in range(3)]
f3 = C3["V"] / np.sqrt(2 * C3["m"] * C3["w"])
hat_h = [np.sqrt(hbar) * sum(A[a] * f3[j, a] + A[a].conj().T * np.conj(f3[j, a]) for a in range(3)) for j in range(3)]
hat_p = [np.sqrt(hbar) * sum(A[a] * (-1j * C3["m"] * C3["w"][a] * f3[j, a]) + A[a].conj().T * np.conj(-1j * C3["m"] * C3["w"][a] * f3[j, a]) for a in range(3)) for j in range(3)]
occ = np.array(np.meshgrid(*[np.arange(nc)] * 3, indexing="ij")).reshape(3, -1).T
low = np.where(occ.sum(1) <= nc - 3)[0]
comm_err = max(np.abs((hat_h[j] @ hat_p[k] - hat_p[k] @ hat_h[j])[np.ix_(low, low)] - 1j * hbar * (j == k) * np.eye(len(low))).max() for j in range(3) for k in range(3))
print("max |[h_j, p_k] - i hbar delta| on low states:", comm_err)
Hch = sum(P @ P for P in hat_p) / (2 * C3["m"]) + 0.5 * sum(C3["S"][j, k] * hat_h[j] @ hat_h[k] for j in range(3) for k in range(3))
Hn = sum(hbar * C3["w"][a] * (A[a].conj().T @ A[a] + 0.5 * np.eye(nc**3)) for a in range(3))
print("max |H_chain - sum hbar w (n + 1/2)| on low states:", np.abs((Hch - Hn)[np.ix_(low, low)]).max(), "; zero-point sum hbar w/2 =", hbar * C3["w"].sum() / 2)
print("=== zero-point energy vs N (L = 40, mu = 1) ===")
for N in [16, 64, 256, 1024]:
    CN = cadeia(N, 1.0)
    print(f"N={N}: sum hbar w_a / 2 = {CN['w'].sum() / 2:.4e}")
