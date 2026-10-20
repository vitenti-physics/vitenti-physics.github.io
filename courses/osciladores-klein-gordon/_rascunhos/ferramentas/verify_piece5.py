import itertools, math
import numpy as np
import sympy as sp

print("=== conservation of the KG product (symbolic) ===")
t, x, c, mu, rho = sp.symbols("t x c mu rho", positive=True)
p1, p2 = sp.Function("phi1")(t, x), sp.Function("phi2")(t, x)
dens = sp.I * rho * (sp.conjugate(p1) * sp.diff(p2, t) - sp.conjugate(sp.diff(p1, t)) * p2)
ddt = sp.diff(dens, t)
kg = lambda f: c**2 * sp.diff(f, x, 2) - mu**2 * f      # phi_tt from KG
ddt = ddt.subs({sp.diff(p2, t, 2): kg(p2), sp.conjugate(sp.diff(p1, t, 2)): sp.conjugate(kg(p1))})
flux = sp.I * rho * c**2 * (sp.conjugate(p1) * sp.diff(p2, x) - sp.conjugate(sp.diff(p1, x)) * p2)
print("d/dt density - d/dx flux =", sp.simplify(sp.expand(ddt - sp.diff(flux, x))))

print("=== lattice checks (N = 64, L = 40, mu = 1) ===")
N, L, cc, mm, rr = 64, 40.0, 1.0, 1.0, 1.0
dx = L / N
xj = np.arange(N) * dx
a = np.fft.fftfreq(N, d=1 / N).astype(int)
k = 2 * np.pi * a / L
w = np.sqrt(mm**2 + 4 * cc**2 / dx**2 * np.sin(k * dx / 2) ** 2)
u = np.exp(1j * np.outer(xj, k)) / np.sqrt(L)                  # column a = u_a(x_j)
def modo(tt):                                                 # f_a(t, x_j) and Pi_a = rho d_t f_a
    f = u * np.exp(-1j * w * tt) / np.sqrt(2 * rr * w)
    return f, rr * (-1j * w) * f
def prod(f1, P1, f2, P2):                                     # <phi1, phi2> = i sum dx (phi1^* Pi2 - Pi1^* phi2)
    return 1j * dx * (f1.conj().T @ P2 - P1.conj().T @ f2)
for tt in [0.0, 3.7]:
    f, P = modo(tt)
    G = prod(f, P, f, P); Gc = prod(f.conj(), P.conj(), f.conj(), P.conj()); Gm = prod(f, P, f.conj(), P.conj())
    print(f"t={tt}: max|<f_a,f_b> - delta| = {np.abs(G - np.eye(N)).max():.1e}, max|<f*_a,f*_b> + delta| = {np.abs(Gc + np.eye(N)).max():.1e}, max|<f_a,f*_b>| = {np.abs(Gm).max():.1e}")
x0, s = L / 2, 0.5
phi0 = np.exp(-((xj - x0) ** 2) / (2 * s**2)); Pi0 = 0.3 * np.sin(2 * np.pi * xj / L)
f, P = modo(0.0)
alpha = prod(f, P, phi0[:, None], Pi0[:, None])[:, 0]
ht = np.fft.fft(phi0) * dx / np.sqrt(L) ; pt = np.fft.fft(Pi0) * dx / np.sqrt(L)   # tilde phi^a, tilde Pi^a (Riemann sums)
alpha2 = (rr * w * ht + 1j * pt) / np.sqrt(2 * rr * w)
print("alpha = <f_a, phi> vs (rho w phi~ + i Pi~)/sqrt(2 rho w):", np.abs(alpha - alpha2).max())
rec = (f @ alpha + f.conj() @ alpha.conj())
recP = (P @ alpha + P.conj() @ alpha.conj())
print("reconstruction max err phi, Pi:", np.abs(rec - phi0).max(), np.abs(recP - Pi0).max())
H_direct = dx * np.sum(Pi0**2 / (2 * rr) + rr * cc**2 / 2 * ((np.roll(phi0, -1) - phi0) / dx) ** 2 + rr * mm**2 / 2 * phi0**2)
print("H direct", H_direct, " sum w|alpha|^2", np.sum(w * np.abs(alpha) ** 2))
psi_plus = f @ alpha
print("psi_+ real part = phi/2? max|2 Re psi_+ - phi| =", np.abs(2 * psi_plus.real - phi0).max())

print("=== classical Poisson brackets of alpha (lattice, exact linear algebra) ===")
# phase space coordinates phi_j, Pi_j with {phi_j, Pi_k} = delta_jk/dx ; alpha linear: alpha_a = sum_j A_aj phi_j + B_aj Pi_j
A_ = 1j * dx * (-(P.conj().T))   # coefficient of phi_j in alpha_a : i dx (-Pi_a^*)
B_ = 1j * dx * (f.conj().T)      # coefficient of Pi_j : i dx f_a^*
# {alpha_a, conj(alpha_b)} = sum_j (A_aj conj(B_bj) - B_aj conj(A_bj)) / dx
PB = (A_ @ B_.conj().T - B_ @ A_.conj().T) / dx
PB2 = (A_ @ B_.T - B_ @ A_.T) / dx
print("max|{alpha_a, alpha_b^*} + i delta| =", np.abs(PB + 1j * np.eye(N)).max(), "; max|{alpha_a, alpha_b}| =", np.abs(PB2).max())

print("=== field commutator sum ===")
yj = xj
comm = np.zeros((N, N), complex)
fx, Px = modo(1.3)
comm = fx @ Px.conj().T - fx.conj() @ Px.T                    # sum_a f_a(x) Pi_a^*(y) - f_a^*(x) Pi_a(y)  (hbar = 1)
print("max|[phi(x_j), Pi(x_k)]/i - delta_jk/dx| =", np.abs(comm / 1j - np.eye(N) / dx).max())

print("=== Fock space, d = 3 modes, n <= 4 ===")
d = 3
rng = np.random.default_rng(0)
def sym(T):
    n = T.ndim
    if n <= 1: return T
    return sum(np.transpose(T, perm) for perm in itertools.permutations(range(n))) / math.factorial(n)
def adag(chi, T):                       # n -> n+1
    n = T.ndim
    out = np.multiply.outer(chi, T) if n > 0 else chi * T
    return math.sqrt(n + 1) * sym(out)
def ann(psi, T):                        # n -> n-1
    n = T.ndim
    if n == 0: return np.zeros(())
    return math.sqrt(n) * np.tensordot(psi.conj(), T, axes=(0, 0))
psi = rng.normal(size=d) + 1j * rng.normal(size=d); chi = rng.normal(size=d) + 1j * rng.normal(size=d)
for n in range(0, 4):
    T = sym(rng.normal(size=(d,) * n) + 1j * rng.normal(size=(d,) * n)) if n > 0 else np.array(1.0 + 0.5j)
    lhs = ann(psi, adag(chi, T)) - (adag(chi, ann(psi, T)) if n > 0 else 0)
    print(f"n={n}: max|[a(psi), a+(chi)] T - <psi,chi> T| = {np.abs(lhs - np.vdot(psi, chi) * T).max():.1e};  [a+(psi),a+(chi)]T: {np.abs(adag(psi, adag(chi, T)) - adag(chi, adag(psi, T))).max():.1e}")
e = np.eye(d)
T = np.array(1.0)
for _ in range(3): T = adag(e[0], T)
print("a_0^dag^3 |0>: norm^2 =", np.vdot(T, T).real, "(expected 3! = 6)")
num = sum(adag(e[b], ann(e[b], T)) for b in range(d))
print("number operator on it:", np.abs(num - 3 * T).max())
