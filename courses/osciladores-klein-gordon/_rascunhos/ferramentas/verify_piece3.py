import sympy as sp
import numpy as np

print("=== symbolic, N = 6 ===")
N = 6
b = sp.symbols(f"b0:{N}")
T = sp.zeros(N, N)
for j in range(N):
    T[j, (j + 1) % N] = 1
C = sp.zeros(N, N)
for r in range(N):
    C += b[r] * T**r
print("C^j_k == b_{k-j}:", all(C[j, k] == b[(k - j) % N] for j in range(N) for k in range(N)))
print("first rows:", [list(C[j, :]) for j in range(3)])
Cg = sp.Matrix(N, N, lambda j, k: sp.Symbol(f"c{j}{k}"))
eqs = sp.solve(list(Cg * T - T * Cg), list(Cg.free_symbols), dict=True)[0]
Csol = Cg.subs(eqs)
print("CT = TC forces circulant:", all(sp.simplify(Csol[j, k] - Csol[0, (k - j) % N]) == 0 for j in range(N) for k in range(N)))
w = sp.exp(2 * sp.pi * sp.I / N)
for a in [1, 2]:
    va = sp.Matrix([w ** (a * j) for j in range(N)]) / sp.sqrt(N)
    lam = sum(b[r] * w ** (a * r) for r in range(N))
    print(f"a={a}: C v_a - lambda_a v_a == 0:", sp.simplify(C * va - lam * va) == sp.zeros(N, 1), "; T v_a = e^{2pi i a/N} v_a:", sp.simplify(T * va - w**a * va) == sp.zeros(N, 1))
x = sp.symbols("x", real=True)
print("lambda_a(D): 2 - 2cos(x) - 4 sin^2(x/2) =", sp.simplify(2 - 2 * sp.cos(x) - 4 * sp.sin(x / 2) ** 2))
print("Hartley: (1-i)e^{ix}/2 + (1+i)e^{-ix}/2 - (cos x + sin x) =", sp.simplify(sp.expand_complex((1 - sp.I) * sp.exp(sp.I * x) / 2 + (1 + sp.I) * sp.exp(-sp.I * x) / 2) - (sp.cos(x) + sp.sin(x))))

print("\n=== numerics ===")
rng = np.random.default_rng(42)
for Nn in [8, 15]:
    j = np.arange(Nn)
    V = np.exp(2j * np.pi * np.outer(j, j) / Nn) / np.sqrt(Nn)   # column a = v_a
    Tn = np.roll(np.eye(Nn), 1, axis=1)
    bb = rng.normal(size=Nn)
    Cn = sum(bb[r] * np.linalg.matrix_power(Tn, r) for r in range(Nn))
    L = V.conj().T @ Cn @ V
    print(f"N={Nn}: max|V^dag V - I| = {np.abs(V.conj().T @ V - np.eye(Nn)).max():.1e}, "
          f"max off-diag of V^dag C V = {np.abs(L - np.diag(np.diag(L))).max():.1e}, "
          f"max|diag - N ifft(b)| = {np.abs(np.diag(L) - Nn * np.fft.ifft(bb)).max():.1e}")
    h = rng.normal(size=Nn)
    print(f"   v_a^dag h == fft(h, ortho)[a]: {np.allclose(V.conj().T @ h, np.fft.fft(h, norm='ortho'))};  h == ifft(ht, ortho): {np.allclose(np.fft.ifft(np.fft.fft(h, norm='ortho'), norm='ortho').real, h)}")
    chi = (np.cos(2 * np.pi * np.outer(j, j) / Nn) + np.sin(2 * np.pi * np.outer(j, j) / Nn)) / np.sqrt(Nn)  # column a = chi_a
    Dn = 2 * np.eye(Nn) - Tn - Tn.T
    lamD = 4 * np.sin(np.pi * j / Nn) ** 2
    print(f"   Hartley: max|chi^T chi - I| = {np.abs(chi.T @ chi - np.eye(Nn)).max():.1e}, max|D chi_a - lambda_a chi_a| = {np.abs(Dn @ chi - chi * lamD).max():.1e}")
    for K0 in [0.0, 0.2]:
        wnum = np.sqrt(np.clip(np.linalg.eigvalsh(Dn + K0 * np.eye(Nn)), 0, None))
        wf = np.sqrt(K0 + lamD)
        print(f"   K0={K0}: max|omega numeric - formula| = {np.abs(np.sort(wnum) - np.sort(wf)).max():.1e}, omega_max = {wf.max():.4f}, 2/omega_max = {2 / wf.max():.4f}, bound sqrt(4+K0) = {np.sqrt(4 + K0):.4f}")
