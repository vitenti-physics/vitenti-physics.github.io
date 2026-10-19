import sympy as sp
e, xi, eta = sp.symbols('varepsilon xi eta', real=True)
# symplectic Euler step in (xi, eta)
M_se = sp.Matrix([[1 - e**2, e], [-e, 1]])
print("det SE =", sp.simplify(M_se.det()))
v = sp.Matrix([xi, eta])
v1 = M_se * v
Hse = lambda a, b: (a**2 + b**2 - e*a*b)/2
print("SE conserves Htilde:", sp.simplify(Hse(v1[0], v1[1]) - Hse(xi, eta)))
# Verlet: p half, x, p half with F=-m w^2 x, in xi/eta: eta -= e/2 xi; xi += e eta; eta -= e/2 xi
xi1 = xi; eta1 = eta
eta1 = eta1 - e/2*xi1
xi1 = xi1 + e*eta1
eta1 = eta1 - e/2*xi1
M_v = sp.Matrix([[sp.diff(xi1, xi), sp.diff(xi1, eta)], [sp.diff(eta1, xi), sp.diff(eta1, eta)]])
print("Verlet matrix:", sp.simplify(M_v), "det =", sp.simplify(M_v.det()), "trace =", sp.simplify(M_v.trace()))
Hv = lambda a, b: b**2/2 + (1 - e**2/4)*a**2/2
print("Verlet conserves Htilde_V:", sp.simplify(Hv(xi1, eta1) - Hv(xi, eta)))
# RK4
w = sp.symbols('w')
R = 1 + w + w**2/2 + w**3/6 + w**4/24
R2 = sp.expand(R.subs(w, -sp.I*e) * R.subs(w, sp.I*e))
print("|R|^2 =", R2)
# Euler implicit
print("implicit Euler |1/(1+i e)|^2 =", sp.simplify(1/((1+sp.I*e)*(1-sp.I*e))))
# A^2
m, om = sp.symbols('m omega', positive=True)
A = sp.Matrix([[0, 1/m], [-m*om**2, 0]])
print("A^2 =", A*A)
t = sp.symbols('t', real=True)
E = sp.cos(om*t)*sp.eye(2) + sp.sin(om*t)/om*A
Om = sp.Matrix([[0, 1], [-1, 0]])
print("E^T Omega E - Omega =", sp.simplify(E.T*Om*E - Om))
print("exp(tA) check:", sp.simplify(sp.exp(t*A) - E))

# ---- second pass (2026-10-19): powers of A, step matrices as polynomials in eps ----
print("\n=== powers of A ===")
for n in range(2, 7):
    print(f"A^{n} =", sp.simplify(A**n))
N = 14
serie = sum(((t*A)**n / sp.factorial(n) for n in range(N)), sp.zeros(2, 2))
# truncated series vs closed form, numerically at omega=1, m=1, t=2.7
num = {m: 1, om: 1, t: 2.7}
print("series(14 terms) - closed form, max abs =", max(abs(float(v)) for v in (serie - E).subs(num)))
Om = sp.Matrix([[0, 1], [-1, 0]])
Mex = sp.cos(e)*sp.eye(2) + sp.sin(e)*Om
print("exact step in (xi,eta):", Mex, " = rotation", sp.simplify(Mex - sp.Matrix([[sp.cos(e), sp.sin(e)], [-sp.sin(e), sp.cos(e)]])))
ME = sp.eye(2) + e*Om
print("Euler = sqrt(1+e^2) * rotation(arctan e):", sp.simplify(ME - sp.sqrt(1+e**2)*Mex.subs(e, sp.atan(e))))
Seta = lambda a: sp.Matrix([[1, 0], [-a, 1]])
Sxi = lambda a: sp.Matrix([[1, a], [0, 1]])
MSE = Sxi(e)*Seta(e)
print("S_xi S_eta =", MSE, "== SE matrix:", sp.simplify(MSE - M_se))
MV = Seta(e/2)*Sxi(e)*Seta(e/2)
print("Verlet = S_eta(e/2) S_xi(e) S_eta(e/2):", sp.simplify(MV - M_v))
print("Verlet - rotation, series:", sp.simplify((MV - Mex).applyfunc(lambda x: sp.series(x, e, 0, 5).removeO())))
print("SE - rotation, series:", sp.simplify((MSE - Mex).applyfunc(lambda x: sp.series(x, e, 0, 4).removeO())))
th = sp.symbols('theta')
print("Verlet trace = 2 - e^2:", sp.simplify(MV.trace() - (2 - e**2)), "; SE trace:", sp.simplify(MSE.trace() - (2 - e**2)))
print("cos(2 asin(e/2)) =", sp.simplify(sp.cos(2*sp.asin(e/2))))
print("theta_V series:", sp.series(2*sp.asin(e/2), e, 0, 6))
QSE = sp.Matrix([[1, -e/2], [-e/2, 1]]); QV = sp.Matrix([[1 - e**2/4, 0], [0, 1]])
print("M_SE^T Q M_SE - Q =", sp.simplify(MSE.T*QSE*MSE - QSE), "; M_V^T Q M_V - Q =", sp.simplify(MV.T*QV*MV - QV))
# RK4 for linear system: k's
dt = sp.symbols('Delta_t', positive=True)
Ag = sp.MatrixSymbol('A', 2, 2)
Xg = sp.MatrixSymbol('X', 2, 1)
k1 = Ag*Xg
k2 = Ag*(Xg + dt/2*k1)
k3 = Ag*(Xg + dt/2*k2)
k4 = Ag*(Xg + dt*k3)
rk = Xg + dt/6*(k1 + 2*k2 + 2*k3 + k4)
taylor4 = (sp.Identity(2) + dt*Ag + dt**2/2*Ag**2 + dt**3/6*Ag**3 + dt**4/24*Ag**4)*Xg
print("RK4 - 4th-order Taylor (symbolic matrices):", sp.expand(rk - taylor4))
MRK = sum(((e*Om)**k / sp.factorial(k) for k in range(5)), sp.zeros(2, 2))
print("M_RK4 =", sp.simplify(MRK), "; det =", sp.expand(MRK.det()))
print("RK4 eigenvalues |lambda|^2:", [sp.expand(sp.simplify(l*sp.conjugate(l))) for l in MRK.eigenvals()])
