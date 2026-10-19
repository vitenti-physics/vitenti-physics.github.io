# Minicurso "Osciladores Acoplados e a Equação de Klein–Gordon": state and plan

Record of the work started on 2026-10-19. Keep this file current; it lives
in `_rascunhos/`, which Quarto does not render. Old notes, logos and earlier
reveal.js decks are under `/home/sandro/Documents/UEL/Lectures/` (the Insync
copy of the same folder is a mirror).

## 1. What the course is

Ementa (from Sandro): starts from a discrete 1D system of coupled harmonic
oscillators and studies the passage from discrete to continuous: circulant
matrices, Fourier series, the continuum limit giving the Klein–Gordon equation,
with computational implementations and numerical solutions throughout.

- **Format:** two sessions of 1h30. Audience: mid-late undergrad and graduate.
- **Material:** Quarto notes (full derivations, every step) + reveal.js slides,
  same section structure. Slides are projected; students follow the notes on
  the site. Each slide carries an aside linking to the matching notes section.
- **Language:** talk to Sandro in English; all course material in pt-BR.
- **Writing rules:** `~/Projects/NumCosmoDevNotes/memory/writing-rules.md` and
  `qmd-figure-conventions.md` in the same folder. Read them before writing.
  Key points: no em-dashes; `$$` alone on its lines with no blank line around
  the display; `\equiv` only for definitions; one symbol one meaning across all
  pages; upper index for vector components; figure cells `#| code-fold: true`
  with measured numbers in the caption; no metaphors, no agent verbs on
  quantities; verify every formula numerically before writing it.
- **Workflow:** piece by piece, in order. Write one section in notes AND
  slides, render, verify, report, stop for Sandro's review. Do not commit or
  push unless asked.

## 2. Course plan (agreed)

| Session | Piece | Content | ~min |
|---|---|---|---|
| Aula 1 | 1 | Single HO: Hamilton, matrix form $\dot X=\Omega\mathsf H X$, exact $e^{tA}$, Euler / symplectic Euler / Verlet / RK4 | 30 |
| Aula 1 | 2 | Periodic chain with neighbor springs $K$ and on-site springs $K_0$; `np.roll` = shift operator; zero mode demo | 30 |
| Aula 1 | 3 | Circulant matrices: polynomial in $T$, Fourier vectors, eigenvalues = DFT of first row, dispersion relation, normal modes, real bases (Hartley), FFT exact solution vs Verlet, CFL $\Delta t<2/\omega_{\max}$ | 30 |
| Aula 2 | 4 | Continuum limit to Klein–Gordon, discrete vs continuum dispersion, Fourier series as limit of DFT, wave packet inside the light cone | 50 |
| Aula 2 | 5 | Classical creation/annihilation variables: complexified phase space, positive-frequency modes, $\alpha_a$, $H=\sum\omega_a|\alpha_a|^2$, $\{\alpha_a,\alpha_b^*\}=-i\delta_{ab}$, complex structure $J$, one-line quantization preview | 40 |

Part 5 stays classical, ending with only a pointer to quantization (Sandro
did not object to this default; confirm when reaching it).

## 3. Current state

- **Piece 1: written and reviewed (second session, 2026-10-19); slides
  rebuilt in the format of section 7a, awaiting Sandro's look.** `aula1.qmd`
  (sections `#sec-ho*`) and `slides-aula1.qmd` (15 slides). Renders with no
  warnings. Review pass: every formula re-verified with sympy
  (`_rascunhos/ferramentas/verify_piece1.py`); one fix, the energy
  derivation multiplies by $\dot x$, not $m\dot x$. Three figures added
  to the notes for the slides: `fig-ho-animacao` (mass on a spring plus the
  phase-space point, jshtml animation), `fig-ho-rotacao` (exact step as a
  rotation by $\varepsilon$), `fig-ho-um-passo` (eight steps of each method
  with $\varepsilon = 0.8$; measured $H_8/H_0$: Euler 52.3, RK4 0.974,
  symplectic Euler $H/H_0 \in [0.77, 1.57]$, Verlet $[0.84, 1.00]$).
  Third pass (same day), after Sandro's review "matrix exponential too
  abstract, stepping methods need more detail but intuitive, A needs a
  clear expression and a few powers": `#sec-ho-exata` now derives $e^{tA}$
  from the series with $A^2 \dots A^6$ written out and the even/odd split;
  `#sec-ho-numerico` presents every method as a polynomial matrix
  $M(\varepsilon)$ replacing the exact rotation
  $M_{\rm ex}(\varepsilon) = \cos\varepsilon\,I + \sin\varepsilon\,\Omega$
  (order, determinant, eigenvalues), with the symplectic ones as products
  of shears $S_\xi$, $S_\eta$ (new `fig-ho-cisalhamento`), the Verlet
  phase $\theta_{\rm V} = 2\arcsin(\varepsilon/2)$, RK4 as the fourth-order
  truncation of the series, and a summary table (`#sec-ho-resumo`). The
  real symplectic product $\Omega(X, Y) = X^{\mathsf T}\Omega Y$ of two
  solutions is introduced in `#sec-ho-exata` (conserved, equals
  $\det = 1$ in 2D) with a numerical check per method at the end of
  `#sec-ho-implementacao`; the complexified product stays in Piece 5. All
  identities verified in `verify_piece1.py` (second block).
  The $z = \xi + i\eta$ paragraph in `#sec-ho-adimensionais` stays as it
  is (Sandro, 2026-10-19: his idea for avoiding complex variables was too
  advanced for this course).
  Measured numbers (2513 steps, $\varepsilon=0.05$): Euler $5.30\times10^2 =
  (1+\varepsilon^2)^n-1$; RK4 $5.45\times10^{-7}=n\varepsilon^6/72$; Verlet max
  $6.25\times10^{-4}=\varepsilon^2/4$; symplectic Euler max
  $2.56\times10^{-2}\approx\varepsilon/2$; `expm` vs closed form $2\times10^{-15}$.
  Formulas for the integrators were derived with sympy: symplectic Euler
  conserves $(\xi^2+\eta^2-\varepsilon\xi\eta)/2$, Verlet conserves
  $\eta^2/2+(1-\varepsilon^2/4)\xi^2/2$, Verlet stable iff $\varepsilon<2$, RK4
  $|R|^2=1-\varepsilon^6/72+\varepsilon^8/576$.
  Slide overflow is now measured by `_rascunhos/ferramentas/mede_slides.sh`
  (section 7a); every slide of the rebuilt deck reports 0 px overflow.
- **Piece 2: written (2026-10-19), awaiting Sandro's review.** Section
  `#sec-cadeia` of `aula1.qmd` (model, forces, translation operator $T$
  and $D$, Hamiltonian and block form, zero mode with $P$, $\bar h$ and the
  Poisson bracket, Verlet simulation, 6 exercises) and 11 slides (deck now
  31 slides, all at 0 px overflow). Figures: `fig-cadeia-esquema` (chain +
  ring), `fig-cadeia-matrizes` ($T$, $T^{-1}$, $D$ for $N=8$),
  `fig-modo-zero`, `fig-cadeia-xt` ($h^j-\bar h$ vs $j,t$),
  `fig-cadeia-animacao` (duplicated in the deck with the simulation, see
  7a). Sandro's initial condition kept in shape ($h^0 = 1$, particles 1
  and 2 with velocity 0.5; his old script used $h^0 = 10$), $N=15$,
  $m=K=1$, $K_0 \in \{0, 0.2\}$, $\Delta t=0.05$, 1200 steps. Measured:
  energy error $\le 1.5\times10^{-3}$; $K_0=0$: $P$ conserved to
  $4.7\times10^{-15}$, $\bar h$ on the free-particle line to
  $1.3\times10^{-14}$, $\bar h(60) = 4.07$; $K_0=0.2$: $\bar h$ follows
  the $\omega_0$ oscillator to $8.1\times10^{-5}$; particle 7 first moves
  by more than 0.05 at $t=4.6$. Symbols added to `#sec-notacao`: $N$,
  $\Delta x$, $h^j$/$h$, $K$, $K_0$, $F_j$, $\omega_0$, $\ell_0$, $\tau$,
  $T$, $D$, $V$, $p_j$/$p$, $\mathbf 1$, $\beta$ (uniform shift; $c$ was
  avoided because of the speed in Piece 4), $P$, $\bar h$, $\{f,g\}$.
  Checks in `_rascunhos/ferramentas/verify_piece2.py`.
- **Piece 3: written (2026-10-19), awaiting Sandro's review.** Section
  `#sec-circulantes` of `aula1.qmd`: circulants as polynomials in $T$ and
  the proof that $CT = TC$ iff circulant; eigenvectors of $T$ (Fourier
  basis $v_a$, orthonormality by the geometric sum, symmetric range and
  aliasing); DFT with numpy `norm="ortho"`; $\lambda_a(C)$ = $N$ ifft of
  the first row, $\lambda_a(D) = 4\sin^2(\pi a/N)$; dispersion relation
  $\omega_a$ ($\omega_{a=0} = \omega_0$, consistent with Piece 2); normal
  modes and $H$ as a sum over modes; Hartley basis $\chi_a$ (standing vs
  travelling waves); exact solution by FFT vs Verlet; stability
  $\Delta t < 2/\omega_{\max}$; 7 exercises. 11 slides (deck: 42 slides,
  all at 0 px overflow). Measured: dispersion numeric vs formula
  $2.3\times10^{-15}$; $\omega_{\max} = 1.989$ ($K_0=0$) and $2.039$
  ($K_0=0.2$) for $N=15$, below the even-$N$ bound; Verlet error vs FFT at
  $t=60$: $1.15\times10^{-2}$ and $1.22\times10^{-2}$; halving $\Delta t$
  divides it by 4.00; at $1.02\times2/\omega_{\max}$ only $a=\pm7$ are
  unstable and $E$ grows by 2.2226 per step, equal to $\lambda^2$ of the
  Verlet matrix ($E/E_0 = 7\times10^{102}$ after 300 steps); at
  $0.98\times$ it stays in $[0.51, 1.04]$. Piece 2 sentence corrected:
  for odd $N$ the top frequency is below $\sqrt{(4K+K_0)/m}$. Checks in
  `_rascunhos/ferramentas/verify_piece3.py` and `numbers_piece3.py`.
  Aula 1 is now complete (pieces 1 to 3); `aula1.qmd` ends with the site
  footer include. The old draft `aula1-completa-v0.qmd` is fully
  superseded. Written
  BEFORE the writing rules were known and BEFORE the symbol renames below: it
  still uses $k$ for the HO spring, $h=\omega\Delta t$, $P$ for scaled momentum,
  `F` for the DFT matrix and $H^{(a)}$ for the Hartley basis, and lower indices
  $h_j$. Rework it under the rules when its turn comes; its code (np.roll
  force, Verlet chain, FFT exact solution, CFL demo, circulant check) ran
  correctly.
- **Pieces 4 and 5: not written.** Plan in section 5 below.
- `index.qmd`: course page with material table and the notation table
  (`#sec-notacao`), the single arbiter for symbols. Extend it with each piece.

## 4. Conventions fixed so far

- Symbols in use (see `index.qmd#sec-notacao`): $m$, $\kappa$ (HO spring;
  $k$ reserved for wavenumber), $\omega\equiv\sqrt{\kappa/m}$, $x,p$, $H$, $X$,
  $\mathsf H$, $\Omega=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ (same as
  `courses/mathematical-physics/propagators-symplectic.qmd`), $A\equiv\Omega\mathsf H$,
  $\Delta t$, $t_n$, $\xi\equiv\sqrt m\,\omega x$, $\eta\equiv p/\sqrt m$,
  $\varepsilon\equiv\omega\Delta t$, $z\equiv\xi+i\eta$, $\tilde H_{\rm SE}$,
  $\tilde H_{\rm V}$, $F$ (force), $R$ (RK4 amplification).
- Planned for pieces 2–5 (not yet in the table): $N$, $\Delta x$, $h^j$
  (heights, upper index), $p_j$, $K$, $K_0$, $\omega_0^2\equiv K_0/m$, shift
  $T$ with $(Th)^j=h^{j+1}$, $D\equiv 2I-T-T^{-1}$, $P\equiv\sum_j p_j$ (total
  momentum), Fourier vectors $v_a$ with components $v_a^{\,j}=e^{2\pi iaj/N}/\sqrt N$,
  circulant eigenvalues $\lambda_a(C)$, Hartley vectors (not $H$; e.g.
  $\chi_a$), $L\equiv N\Delta x$, $k_a\equiv 2\pi a/L$, continuum: $\rho$, $\tau$
  (tension), $\mu$ (mass term), $c^2\equiv\tau/\rho$, $\phi(t,x)$, $\pi(t,x)$.
  Check collisions before fixing: $c$ (speed) vs circulant first-row
  coefficients (use another letter for those), $T$ (shift) vs tension ($\tau$).
- **Hermitian product:** $\langle V,U\rangle = +i\,V^\dagger\Omega U$, as in the
  mathematical-physics course (see that project's memory
  `mathematical-physics-conventions.md`). With
  $u(t)=(f,\,m\dot f)$, $f=e^{-i\omega t}/\sqrt{2m\omega}$: $\langle u,u\rangle=1$,
  $\langle u^*,u^*\rangle=-1$, $\langle u,u^*\rangle=0$, and for real $X$,
  $\alpha=\langle u,X\rangle=(m\omega x+ip)/\sqrt{2m\omega}$, $H=\omega|\alpha|^2$,
  $\{\alpha,\alpha^*\}=-i$. Note $z=\sqrt{2\omega}\,\alpha$ (Piece 1 already
  points to this).

## 5. Plan for pieces 4 and 5 (math checked in the original session, re-verify)

**Piece 4 (continuum limit).** $x_j=j\Delta x$, $h^j(t)=\phi(t,x_j)$; scalings
$m=\rho\Delta x$, $K=\tau/\Delta x$, $K_0=\rho\mu^2\Delta x$; Taylor
$h^{j+1}-2h^j+h^{j-1}=\Delta x^2\phi''+\Delta x^4\phi''''/12+\dots$; result
$\phi_{tt}-c^2\phi_{xx}+\mu^2\phi=0$, $c^2=\tau/\rho$. Hamiltonian density with
$\pi=p_j/\Delta x$, $\{\phi(x),\pi(y)\}=\delta(x-y)$. Dispersion
$\omega^2=\mu^2+(4c^2/\Delta x^2)\sin^2(k\Delta x/2)\to\mu^2+c^2k^2$; group
velocity $<c$. $v_a^{\,j}=\sqrt{\Delta x}\,e^{ik_ax_j}/\sqrt L$; completeness
needs $a\in\mathbb Z$ (symmetric range before the limit). Numerics:
dispersion convergence plot; Gaussian packet / single displaced mass evolved
with FFT, $x$–$t$ image with light-cone lines, $\mu=0$ vs $\mu>0$ (lattice
tails beyond the cone are exponentially small, mention).

**Piece 5 (classical $a,a^*$).** Single HO first (formulas above), complex
structure $J=-A/\omega$, $g=\mathsf H/\omega=\Omega(\cdot,J\cdot)$ positive.
General $J$ with $J^2=-I$ compatible with $\Omega$:
$J=\begin{pmatrix}\sinh\zeta & -e^{-\gamma}\cosh\zeta\\ e^{\gamma}\cosh\zeta & -\sinh\zeta\end{pmatrix}$
(positivity needs the lower-left entry $>0$); only $\zeta=0$,
$e^\gamma=m\omega$ commutes with $A$. Chain:
$u_a=(v_af_a,\,m\dot f_av_a)$, $f_a=e^{-i\omega_at}/\sqrt{2m\omega_a}$;
$\langle u_a,u_b\rangle=\delta_{ab}$, $\langle u_a,u_b^*\rangle=0$;
$\alpha_a=(m\omega_a\tilde h^a+i\tilde p^a)/\sqrt{2m\omega_a}$ with $\alpha_a$
and $\alpha_{-a}$ independent (unlike $\tilde h^{-a}=\tilde h^{a*}$);
$H=\sum\omega_a|\alpha_a|^2$, $\{\alpha_a,\alpha_b^*\}=-i\delta_{ab}$ (verified).
Zero mode with $K_0=0$ has no $\alpha_0$: the mass term regulates it (ties to
the zero-mode demo of Piece 2). $|\alpha_a|^2$ are action variables. Numerics:
check $\alpha_a(t)=\alpha_a(0)e^{-i\omega_at}$ from a simulation; build initial
data with only $a>0$ to get a one-directional packet.

## 6. Sandro's own points to keep

- The zero-mode drift in his old script is deliberate: it shows that a
  massless field has no preferred value (shift symmetry $h^j\to h^j+c$,
  conserved $P$). Keep as the Piece 2 demo with his initial conditions
  (one mass displaced, masses 2 and 3 with velocity 0.5), contrast with
  $K_0>0$. I suggested "valor preferido / origem" instead of "escala".
- The old student notes contain many errors: never reuse their derivations.
  Errata: `~/Insync/vitenti@uel.br/Google Drive/UEL/Lectures/qft_in_curved/notas_antigas/REVIEW.md`
  (old `main.tex` and `plot_oscilators.py` are in the same folder).

## 7. Build and site infrastructure (added in this work, uncommitted)

- `_quarto.yml`: `project.post-render: ./ipynb_to_py.sh`; sidebar section
  "Minicourses" with Overview, Aula 1, Slides Aula 1 (add Aula 2 pages when
  they exist; listing a missing file breaks the render).
- `ipynb-strip-xref.lua` (copied from NumCosmo docs), `ipynb_to_py.sh`
  (uses `nbconvert --to python`; `--to script` produced `.txt`),
  `add-py-link.html` (adds "Python" download and "Abrir no Colab" pointing to
  `gh-pages`; Colab works only after publishing).
- `_quarto-drafts.yml`: profile with `draft-mode: visible` and
  `cookie-consent: false`. Plain `quarto render` writes draft pages as empty
  90-byte files; review with `quarto preview --profile drafts`.
- Notes pages: formats `html` + `ipynb` (`eval: true`), `format-links: [ipynb]`.
  `eval: true` is required: slides use `{{< embed aula1.qmd#fig-... >}}`,
  which fails when the notebook is not executed.
- The cookie banner cannot be disabled per page; in class, click it once.
- CI (`.github/workflows/quarto.yml`) renders on push to `main` and installs
  `jupyter matplotlib scipy astropy`; code must stay fast and use only
  numpy/scipy/matplotlib. Use `anim.to_jshtml()` for animations (no ffmpeg,
  also works in JupyterLite later).
- Pre-existing uncommitted changes by Sandro, not ours: `.gitignore`,
  `EMEnergy.pdf`, `MovingChargeRadiation.pdf`. Leave them alone.

## 7a. Slide format (fixed 2026-10-19 after Sandro's review)

Sandro's review of the first deck: equations cut off at the right, font too
big, "much more visual approach", use our logos (not DESC). The deck now
follows his earlier lectures (`~/Documents/UEL/Lectures/*/presentation.qmd`):

- `format: revealjs`, 1920x1080, `margin: 0.04`, `center: true`,
  `theme: [serif, slides.scss]` (Crimson Pro, root font 36px),
  `include-after-body: slides-titulo.html` (hides footer and slide number on
  the title slide), title slide black with `logos.png` bottom-left
  (`logos.png` = Grupo de Cosmologia white-text logo + NumCosmo logo,
  composed with PIL from `~/Documents/UEL/Lectures/data_modelling/
  Logo_Grupo_de_Cosmologia_UEL_Oficial_Branco.png` and
  `.../primordial_to_cosmic_web/Logo_NumCosmo_Black_BackGround.png`).
- Notes figures are rendered at `figure.dpi = 150` (setup cell) so they fill
  the 1920 px slide; `slides.scss` caps embedded images at 100% width and
  640 px height. Animations keep `dpi=100` (jshtml size: 80 frames = 3.8 MB).
  Figures are exported with a tight bounding box, so a row of equal-aspect
  panels has its height fixed by the panel count, not by `figsize[1]`;
  shrink the whole `figsize` to make such a figure shorter on a slide
  (`fig-ho-cisalhamento` is 10.5 x 3.2 in for that reason).
- One idea per slide: at most three display lines of math, two or three
  bullets, and a figure wherever the notes have one.
- **Embedded figure cells must contain only the figure.** `{{< embed >}}`
  copies the whole cell output, so a `print` in a figure cell shows up on
  the slide; put numerical checks in a separate cell next to it. Captions
  of embedded figures must not contain cross-references (`@sec-...`,
  `@eq-...`): they are rendered inside the deck, where the targets do not
  exist, and Quarto warns "Unable to resolve crossref".
- **Stale intermediates.** If a render fails with "File
  'aula1_files/figure-ipynb/...png' not found in resource path" (seen
  twice), remove `courses/osciladores-klein-gordon/aula1_files`,
  `slides-aula1_files`, any `*.quarto_ipynb*` there, and `.quarto/embed`,
  then render again. Figures come from the
  notes via `{{< embed aula1.qmd#fig-... >}}`; the notes figure is therefore
  designed with the slide in mind (wide, 1x4 panels rather than 2x2).
- **HTML outputs cannot be embedded across documents.** Quarto renders the
  embedded notebook through its ipynb path, which drops every `text/html`
  output (verified with a 3-cell test: even `HTML("<b>oi</b>")` is lost, and
  `ipynb-output: all` fails to parse). Animation cells (`to_jshtml`) are
  therefore duplicated in the slides with `#| echo: false` and a one-line
  comment saying so. The downloadable `aula1.ipynb` also lacks those outputs;
  students regenerate them by running the notebook.
- Overflow check: from the repo root, after `quarto render ... --profile
  drafts`, run `courses/osciladores-klein-gordon/_rascunhos/ferramentas/
  mede_slides.sh courses/osciladores-klein-gordon/slides-aula1.html`. It
  serves `_site` with node, injects a script that waits for MathJax 2
  (`MathJax.Hub.Queue`), forces each section visible, and reports right and
  bottom overflow in slide pixels; headless Firefox is held on the page by a
  slow image so the report is posted before the screenshot. Acceptance:
  every slide `right+0 bottom+0`. The script unhides each slide and its
  parent vertical stack, loads every lazy `data-src` image first, and
  ignores MathJax's hidden assistive MathML (`.MJX_Assistive_MathML`),
  which otherwise reports up to 200 px of bottom overflow that is not
  visible. Screenshots: `slides-aula1-check.html?slide=N` (flat index,
  navigated by the injected script; the URL hash does not work headless). Screenshots of slide N: same server,
  `slides-aula1-check.html#/N` with `firefox --headless --screenshot`.

## 8. Next action

Sandro reviews Piece 3 (notes `#sec-circulantes`, slides 32 to 42; Piece
2 review may come with it). Then Aula 2: create `aula2.qmd` and
`slides-aula2.qmd` with the same front matter as Aula 1 (notes: html +
ipynb with `eval: true`; deck: section 7a format), add both to the
`_quarto.yml` sidebar only once the files exist, update the material
table in `index.qmd`, and write Piece 4 (continuum limit, section 5)
under the symbols of section 4 ($x_j = j\Delta x$, $L \equiv N\Delta x$,
$k_a \equiv 2\pi a/L$, $\rho$, $\tau$, $\mu$, $c^2 \equiv \tau/\rho$,
$\phi$, $\pi$). Aula 2 must be self-contained for code: it is a separate
notebook, so the functions it needs from Aula 1 are redefined there.
Render with `--profile drafts`, run `mede_slides.sh`, report, stop.
