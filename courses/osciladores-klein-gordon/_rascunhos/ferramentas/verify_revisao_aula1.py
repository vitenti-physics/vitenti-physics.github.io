"""Checks for the formulas added in the Aula 1 review (2026-10-20)."""
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(7)
Om = np.array([[0.0, 1.0], [-1.0, 0.0]])
S_eta = lambda a: np.array([[1.0, 0.0], [-a, 1.0]])
S_xi = lambda a: np.array([[1.0, a], [0.0, 1.0]])

# gradient and Omega: chi_dot = Omega grad H, energy conserved
m, kappa = 1.3, 2.1
Hm = np.diag([kappa, 1 / m])
A = Om @ Hm
assert np.allclose(A, [[0, 1 / m], [-kappa, 0]])
assert np.allclose(A @ A, -(kappa / m) * np.eye(2))
for _ in range(5):
    chi = rng.normal(size=2)
    grad = Hm @ chi
    assert abs(grad @ Om @ grad) < 1e-15
# derivatives = powers of A along a solution
chi0 = rng.normal(size=2)
t, h = 0.7, 1e-3
sol = lambda tt: expm(tt * A) @ chi0
d2 = (sol(t + h) - 2 * sol(t) + sol(t - h)) / h**2
assert np.allclose(d2, A @ A @ sol(t), atol=1e-5)
# A^T Omega + Omega A = 0 for any symmetric H (2N = 6)
Hs = rng.normal(size=(6, 6)); Hs = Hs + Hs.T
I3 = np.eye(3); Om6 = np.block([[0 * I3, I3], [-I3, 0 * I3]])
A6 = Om6 @ Hs
assert np.allclose(A6.T @ Om6 + Om6 @ A6, 0)
assert np.allclose(expm(0.4 * A6).T @ Om6 @ expm(0.4 * A6), Om6)
# area figure value
c1, c2 = np.array([1.0, 0.4]), np.array([0.2, 1.2])
A4 = np.array([[0, 1.0], [-4.0, 0]])
for tt in (0, 0.3, 0.7):
    print("area", tt, (expm(tt * A4) @ c1) @ Om @ (expm(tt * A4) @ c2))
# M^T Omega M = det M Omega; eigenvalues from trace and det
for _ in range(5):
    M = rng.normal(size=(2, 2))
    assert np.allclose(M.T @ Om @ M, np.linalg.det(M) * Om)
    lam = np.linalg.eigvals(M)
    assert np.isclose(lam.sum(), np.trace(M)) and np.isclose(lam.prod(), np.linalg.det(M))
# Euler eigenvalues 1 +- i eps
eps = 0.3
assert np.allclose(sorted(np.linalg.eigvals(np.eye(2) + eps * Om), key=np.imag), [1 - 1j * eps, 1 + 1j * eps])
# Euler balance numbers
e1 = np.log(1.01) / (40 * np.pi)
print("Euler 1%: eps =", e1, "steps =", 40 * np.pi / e1)
n = int(20 * 2 * np.pi / 0.05)
print("eps=0.05:", n, "steps, H_n/H_0 =", (1 + 0.05**2) ** n)
# Sp(2,R) as three shears; rotation as shears
for _ in range(5):
    M = rng.normal(size=(2, 2)); M /= np.sqrt(abs(np.linalg.det(M)))
    if np.linalg.det(M) < 0: M[:, 0] *= -1
    F = S_eta((1 - M[1, 1]) / M[0, 1]) @ S_xi(M[0, 1]) @ S_eta((1 - M[0, 0]) / M[0, 1])
    assert np.allclose(F, M)
Mex = np.array([[np.cos(eps), np.sin(eps)], [-np.sin(eps), np.cos(eps)]])
assert np.allclose(S_eta(np.tan(eps / 2)) @ S_xi(np.sin(eps)) @ S_eta(np.tan(eps / 2)), Mex)
Mlow = np.array([[2.0, 0.0], [0.7, 0.5]])  # M12 = 0: four shears
Mp = S_xi(-1.0) @ Mlow
F = S_xi(1.0) @ S_eta((1 - Mp[1, 1]) / Mp[0, 1]) @ S_xi(Mp[0, 1]) @ S_eta((1 - Mp[0, 0]) / Mp[0, 1])
assert np.allclose(F, Mlow)
# force diagram numbers
K, K0, dx = 1.0, 0.5, 1.0
hm, hj, hp = 0.3, 0.8, 0.5
fr, fl = K * np.array([dx, hp - hj]), K * np.array([-dx, hm - hj])
print("forces", fr, fl, fr + fl, -K0 * hj, K * (hp - 2 * hj + hm))
# chain: F = K(Th-h)+K(T^-1h-h)-K0h = -(KD+K0I)h ; mean height and momentum per particle
N = 9
T = np.roll(np.eye(N), 1, axis=1)
D = 2 * np.eye(N) - T - T.T
hv, pv = rng.normal(size=N), rng.normal(size=N)
Fv = K * (T @ hv - hv) + K * (T.T @ hv - hv) - K0 * hv
assert np.allclose(Fv, -(K * D + K0 * np.eye(N)) @ hv)
assert np.allclose(Fv, K * (np.roll(hv, -1) - 2 * hv + np.roll(hv, 1)) - K0 * hv)
assert np.isclose(Fv.mean(), -K0 * hv.mean())  # pbar_dot = -K0 hbar
# Omega and A for the chain, Hamilton eqs
I = np.eye(N); Z = 0 * I
OmN = np.block([[Z, I], [-I, Z]])
HN = np.block([[K * D + K0 * I, Z], [Z, I / m]])
chiN = np.concatenate([hv, pv])
assert np.allclose(OmN @ HN @ chiN, np.concatenate([pv / m, Fv]))
# Fourier eigenvectors: T^N = I, lambda^N = 1, v^j = lambda^j v^0, normalized
for a in range(N):
    lam = np.exp(2j * np.pi * a / N)
    v = lam ** np.arange(N) / np.sqrt(N)
    assert np.allclose(T @ v, lam * v) and np.isclose(np.vdot(v, v), 1)
# dispersion without omega_0
for K0_ in (0.0, 0.2):
    w2 = np.sort(np.linalg.eigvalsh((K * D + K0_ * I) / m))
    a = np.arange(N)
    assert np.allclose(w2, np.sort((K0_ + 4 * K * np.sin(np.pi * a / N) ** 2) / m))
print("all checks passed")
