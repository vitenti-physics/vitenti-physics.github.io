import sympy as sp
print("=== Taylor of the second difference ===")
x, dx, k, c, mu, rho, tau = sp.symbols("x Delta_x k c mu rho tau", positive=True)
f = sp.Function("phi")
expr = sum(sp.diff(f(x), x, n) * (dx**n + (-dx) ** n) / sp.factorial(n) for n in range(0, 8)) - 2 * f(x)
print(sp.simplify(expr))
print("=== Newton with the scalings -> KG ===")
m = rho * dx; K = tau / dx; K0 = rho * mu**2 * dx
phi2, phi4, phi0 = sp.symbols("phi2 phi4 phi0")
rhs = (K * (dx**2 * phi2 + dx**4 * phi4 / 12) - K0 * phi0) / m
print("phi_tt =", sp.expand(rhs), "; with c^2 = tau/rho")
print("omega0^2 = K0/m =", sp.simplify(K0 / m), "; K/m =", sp.simplify(K / m))
print("=== Hamiltonian density ===")
Pi, dphi = sp.symbols("Pi dphi")
pj = Pi * dx
print("p^2/2m =", sp.simplify(pj**2 / (2 * m)), "; K/2 (dh)^2 with dh = dphi*dx:", sp.simplify(K / 2 * (dphi * dx) ** 2), "; K0/2 phi^2 =", sp.simplify(K0 / 2 * phi0**2))
print("=== dispersion ===")
a, N, L = sp.symbols("a N L", positive=True)
w2_lat = K0 / m + 4 * K / m * sp.sin(sp.pi * a / N) ** 2
w2_k = w2_lat.subs(a, k * L / (2 * sp.pi)).subs(N, L / dx)
print("omega^2(k) lattice =", sp.simplify(w2_k))
print("series in dx:", sp.series(sp.simplify(w2_k).subs(tau, rho * c**2), dx, 0, 4))
w = sp.sqrt(mu**2 + 4 * c**2 / dx**2 * sp.sin(k * dx / 2) ** 2)
vg = sp.diff(w, k)
print("group velocity lattice:", sp.simplify(vg), " == c^2 sin(k dx)/(dx omega):", sp.simplify(vg - c**2 * sp.sin(k * dx) / (dx * w)) == 0)
print("vg^2 <= c^2 cos^2(k dx/2): c^2cos^2 - vg^2 =", sp.simplify(c**2 * sp.cos(k * dx / 2) ** 2 - vg**2))
wc = sp.sqrt(mu**2 + c**2 * k**2)
print("continuum group velocity:", sp.simplify(sp.diff(wc, k)), "; phase velocity:", wc / k)
