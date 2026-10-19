import sympy as sp
import numpy as np

print("=== symbolic checks, N = 5 ===")
N = 5
h = sp.symbols(f"h0:{N}", real=True)
p = sp.symbols(f"p0:{N}", real=True)
K, K0, m = sp.symbols("K K_0 m", positive=True)
T = sp.zeros(N, N)
for j in range(N):
    T[j, (j + 1) % N] = 1
hv = sp.Matrix(h)
print("T h =", list(T * hv))
print("T^-1 == T^T:", T.inv() == T.T, "; T^N == I:", T**N == sp.eye(N))
D = 2 * sp.eye(N) - T - T.inv()
print("D =", D.tolist())
print("D symmetric:", D == D.T, "; D 1 =", list(D * sp.ones(N, 1)))
V = K / 2 * sum((h[(j + 1) % N] - h[j]) ** 2 for j in range(N)) + K0 / 2 * sum(hj**2 for hj in h)
print("V - 1/2 h^T (KD+K0 I) h =", sp.expand(V - (hv.T * (K * D + K0 * sp.eye(N)) * hv)[0] / 2))
F = [-sp.diff(V, hj) for hj in h]
newton = [K * (h[(j + 1) % N] - 2 * h[j] + h[(j - 1) % N]) - K0 * h[j] for j in range(N)]
print("-dV/dh^j - [K(h^{j+1}-2h^j+h^{j-1}) - K0 h^j] =", [sp.expand(a - b) for a, b in zip(F, newton)])
print("sum_j F_j =", sp.expand(sum(F)), " (= -K0 N hbar)")
Om = sp.BlockMatrix([[sp.zeros(N), sp.eye(N)], [-sp.eye(N), sp.zeros(N)]]).as_explicit()
Hm = sp.BlockMatrix([[K * D + K0 * sp.eye(N), sp.zeros(N)], [sp.zeros(N), sp.eye(N) / m]]).as_explicit()
A = Om * Hm
print("A blocks: upper-right = I/m:", A[:N, N:] == sp.eye(N) / m, "; lower-left = -(KD+K0 I):", A[N:, :N] == -(K * D + K0 * sp.eye(N)))
print("eigs of D (N=5):", sorted(set(sp.nsimplify(e) for e in D.eigenvals())))
print("alternating vector (N=6): D a =", end=" ")
N6 = 6
T6 = sp.zeros(N6, N6)
for j in range(N6):
    T6[j, (j + 1) % N6] = 1
D6 = 2 * sp.eye(N6) - T6 - T6.T
a = sp.Matrix([(-1) ** j for j in range(N6)])
print(list(D6 * a), "= 4 a:", D6 * a == 4 * a)

print("\n=== natural-length spring: vertical force to linear order ===")
dx, l0, dh = sp.symbols("Delta_x ell_0 delta_h", positive=True)
d = sp.sqrt(dx**2 + dh**2)
Fv = K * (1 - l0 / d) * dh
print("F_vertical series in dh:", sp.series(Fv, dh, 0, 3).removeO(), " ; tau/dx * dh with tau = K(dx - l0):", sp.simplify(K * (dx - l0) / dx * dh))
Fh = K * (1 - l0 / d) * dx
print("F_horizontal series:", sp.series(Fh, dh, 0, 3).removeO())

print("\n=== numerics of the simulation ===")
def passo_verlet(F, x, p, dt, m):
    p = p + 0.5 * dt * F(x)
    x = x + dt * p / m
    p = p + 0.5 * dt * F(x)
    return x, p
def forca_cadeia(h, K, K0):
    return K * (np.roll(h, -1) - 2 * h + np.roll(h, 1)) - K0 * h
def energia_cadeia(h, p, m, K, K0):
    return (np.sum(p**2, axis=-1) / (2 * m) + 0.5 * K * np.sum((np.roll(h, -1, axis=-1) - h) ** 2, axis=-1) + 0.5 * K0 * np.sum(h**2, axis=-1))
def simula(h0, p0, m, K, K0, dt, n_passos, a_cada):
    F = lambda h: forca_cadeia(h, K, K0)
    h, pp = h0.copy(), p0.copy()
    hs, ps = [h.copy()], [pp.copy()]
    for n in range(1, n_passos + 1):
        h, pp = passo_verlet(F, h, pp, dt, m)
        if n % a_cada == 0:
            hs.append(h.copy()); ps.append(pp.copy())
    return dt * a_cada * np.arange(len(hs)), np.array(hs), np.array(ps)
Nn, mm, KK = 15, 1.0, 1.0
h0 = np.zeros(Nn); p0 = np.zeros(Nn)
h0[0] = 1.0; p0[1] = p0[2] = 0.5
dt, n_passos = 0.05, 1200
rng = np.random.default_rng(1)
hr = rng.normal(size=Nn)
Tn = np.roll(np.eye(Nn), 1, axis=1)
print("np.roll(h,-1) == T h:", np.allclose(np.roll(hr, -1), Tn @ hr), "; np.roll(h,1) == T^-1 h:", np.allclose(np.roll(hr, 1), Tn.T @ hr))
Dn = 2 * np.eye(Nn) - Tn - Tn.T
print("force via roll vs via matrix, max diff:", np.abs(forca_cadeia(hr, KK, 0.2) - (-(KK * Dn + 0.2 * np.eye(Nn)) @ hr)).max())
for K0 in [0.0, 0.2]:
    t, hs, ps = simula(h0, p0, mm, KK, K0, dt, n_passos, 5)
    E = energia_cadeia(hs, ps, mm, KK, K0)
    P = ps.sum(axis=1)
    hbar = hs.mean(axis=1)
    w0 = np.sqrt(K0 / mm)
    if K0 == 0:
        pred = h0.mean() + P[0] * t / (Nn * mm)
    else:
        pred = h0.mean() * np.cos(w0 * t) + P[0] / (Nn * mm * w0) * np.sin(w0 * t)
    wmax = np.sqrt((4 * KK + K0) / mm)
    print(f"K0={K0}: E0={E[0]:.4f}, max|E-E0|/E0={np.abs(E/E[0]-1).max():.2e}, P0={P[0]}, max|P-P0|={np.abs(P-P[0]).max():.2e}, "
          f"hbar(0)={hbar[0]:.4f}, hbar(end)={hbar[-1]:.4f}, max|hbar-pred|={np.abs(hbar-pred).max():.2e}, "
          f"omega0={w0:.4f}, omega_max={wmax:.4f}, 2/omega_max={2/wmax:.4f}, range h=({hs.min():.2f},{hs.max():.2f}), "
          f"range h-hbar=({(hs-hbar[:,None]).min():.2f},{(hs-hbar[:,None]).max():.2f}), n frames={len(t)}")

print("\n=== propagation and A^2 ===")
res = {K0: simula(h0, p0, mm, KK, K0, dt, n_passos, 1) for K0 in [0.0, 0.2]}
for K0, (t, hs, ps) in res.items():
    dev = hs - hs.mean(axis=1)[:, None]
    for j in [3, 5, 7]:
        k = np.argmax(np.abs(dev[:, j] - dev[0, j]) > 0.05)
        print(f"K0={K0}: |h^{j}-hbar| changes by 0.05 first at t={t[k]:.2f}")
d0 = res[0.0][1] - res[0.0][1].mean(axis=1)[:, None]
d2 = res[0.2][1] - res[0.2][1].mean(axis=1)[:, None]
for tmax in [5, 10, 60]:
    sel = res[0.0][0] <= tmax
    print(f"max |dev(K0=0.2) - dev(K0=0)| for t <= {tmax}: {np.abs(d2[sel] - d0[sel]).max():.3f}")
S = sp.Symbol("S", commutative=False)
Ablk = sp.Matrix([[0, 1 / m], [-S, 0]])
print("A^2 (blocks, S = KD + K0 I):", (Ablk * Ablk).applyfunc(sp.expand))
