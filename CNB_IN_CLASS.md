# Cosmic Neutrino Background anisotropies in CLASS

This note documents the modifications made to `class_public` (CLASS v3.4.0) to
compute the Cosmic Neutrino Background (CNB) angular power spectrum of
`Tesi_Finale.pdf`, natively inside the code.

Equation numbers refer to the Thesis.

---

## 1. What is computed

For each requested neutrino mass, the code computes the CNB anisotropy transfer
function and the resulting spectra

- `C_l^{nunu}` — the CNB auto-correlation,
- `C_l^{Tnu}`  — its cross-correlation with the CMB temperature (when `tCl` is
  also requested),

from eq. (3.2.24) in the **momentum approach** (the default; see Sec. 2 bis) or
eq. (3.2.32) in the **energy approach**, momentum-averaged following
eq. (3.2.38) and integrated over `k` following eq. (3.2.37).

The Sachs-Wolfe and integrated Sachs-Wolfe terms are kept as two separate
transfer types internally, so they can be switched on and off independently
(`cnb_contributions`) in order to reproduce the component-by-component figures
of Sec. 4.3.

### Usage

```ini
output = tCl,cnbCl

cnb_approach      = momentum          # default; 'energy' is not recommended
cnb_contributions = sw, isw           # optional, default is both
cnb_m_nu          = 1.e-5, 1.e-4, 1.e-3, 4.e-3   # eV; default: the ncdm masses
cnb_z_dec         = 6.e9              # neutrino decoupling redshift
cnb_q_sampling    = laguerre          # default; 'chi' is accurate for heavy masses (Sec. 5)
cnb_q_size        = 15                # laguerre only: Gauss-Laguerre nodes in q
cnb_q_phase_step  = 4                 # chi only: Bessel phase per momentum panel
cnb_q_panel_order = 4                 # chi only: Gauss-Legendre nodes per panel
cnb_q_weight_tol  = 1.e-3             # chi only: phase resolved where weight > tol*max
cnb_q_max         = 25                # chi only: largest sampled q/T_nu0
k_max_tau0_over_l_max = 12            # needed for any CNB run: the default 1.5 truncates C_l^nunu
```

A ready-made input file using the cosmological parameters of Table 4.1 is
provided in `cnb.ini`. The output columns are labelled
`nu[1.000e-04eV]-nu[1.000e-04eV]` and `T-nu[1.000e-04eV]`.

If `cnb_m_nu` is omitted, the masses of the `ncdm` species declared in the
background are used. Giving `cnb_m_nu` explicitly is what allows a scan of the
mass dependence of the CNB signal *without* changing the cosmology, which is what
Figs. 4.5-4.9 do.

---

## 2. The two physical difficulties, and how they were solved

### 2.1 The radial coordinate is momentum dependent

For photons the Bessel functions are always evaluated at `k(tau_0 - tau)`. Massive
neutrinos travel at a momentum-dependent speed `v = q/eps`, so the correct radial
coordinate is the *neutrino* comoving distance

```
chi_nu(tau) = int_tau^tau_0 (q/eps) dtau' ,     eps = sqrt(q^2 + a^2 M^2)
```

which depends on both the neutrino mass and its momentum. This is the origin of
the broad, mass-dependent last scattering surface of Sec. 4.2 / Fig. 4.4.

CLASS already stores the radial coordinate *per transfer type*, in
`ptw->tau0_minus_tau[]`, from which `transfer_radial_coordinates()` builds the
Bessel argument. The CNB therefore fits the existing machinery by simply writing
`chi_nu(tau)` into that array instead of `tau_0 - tau`. Two consequences:

- the trapezoidal weights built from that array are weights in `dchi_nu`, not in
  `dtau`, so the source is multiplied by the Jacobian `dtau/dchi_nu = eps/q`;
- the array stays monotonically decreasing, as the rest of the transfer module
  assumes.

### 2.2 The momentum average must happen after the Bessel convolution

Because each momentum has its own `chi_nu`, the average of eq. (3.2.38) cannot be
applied to the `C_l`, nor to `Delta_l` at fixed `l` in the harmonic module — the
different momenta must first be convolved with *different* Bessel functions.

The average is therefore accumulated on the fly in
`transfer_cnb_compute_for_each_q()`, which loops over neutrino momenta inside the
loop over `l`. This also avoids storing one transfer function per momentum: at
CMB resolution that would have cost several hundred MB.

### 2.3 The potentials at neutrino decoupling

The Sachs-Wolfe term needs `phi` at neutrino decoupling (`z ~ 6e9`,
`tau_dec ~ 8e-5 Mpc`), where every wavenumber of interest is deep outside the
horizon, so `phi` sits at its constant super-horizon **radiation-era** value.

The subtlety is *where* that value can be read. CLASS's standard source sampling
starts when the universe stops being opaque, around `z ~ 1400`, which is already
**matter dominated** — and the super-horizon `phi` drops by the familiar factor
`~9/10` across the radiation-to-matter transition. Reading it there gives 0.6481
instead of the radiation-era 0.6996, an 8% misnormalisation of the Sachs-Wolfe
amplitude (17% in `C_l`), and a value that drifts with the precision settings.

The code therefore prepends a coarse, logarithmically spaced set of sampling
points reaching back into the radiation era
(`perturbations_prepend_radiation_era_sampling()`, precision parameter
`cnb_radiation_era_tau_size`, default 40). The earliest point sits at
`z = min(cnb_z_dec, thermo_z_initial/5)`, i.e. `z = 1e6` with default precision,
where `phi` is within 0.03% of its plateau. These points are cheap: the
large-k modes are integrated from those times anyway, and the sources are smooth
there. `phi` is then read at the first sampled time and the smallest sampled
wavenumber, guarded by a `k*tau < 1e-2` test.

This restricts the implementation to **adiabatic initial conditions**, for which
super-horizon `phi` is constant in radiation domination; isocurvature modes are
rejected with an explicit error in `perturbations_cnb_indices()`.

Side effect: because the source time sampling changes when `cnbCl` is requested,
the CMB spectra shift at the `1e-4` relative level compared with a run without
`cnbCl`. That is a precision-level effect, not a bias, but it means the two runs
are not bit-identical.

---

## 2 bis. Which phase-space assumption is correct

`cnb_approach = momentum` is the **default and the physically correct choice**.
The two options are not two defensible modelling assumptions; only one solves
the equation it is substituted into.

Neutrinos decouple while ultra-relativistic at `T ~ 1 MeV` and are collisionless
thereafter. In the homogeneous background the comoving momentum is conserved
(`dq/dtau = 0`), so the background Boltzmann equation reduces to

    d f0 / d tau |_q  =  0

- `momentum`: `f0(q) = 1/(exp(q)+1)` — manifestly time independent, satisfies
  this identically. It is the same `f0` CLASS uses for its ncdm species, and the
  same one this code uses for the momentum-average weights `q^2 eps f0(q)`.
- `energy`: `f0 = 1/(exp(sqrt(q^2+a^2 M^2))+1)` — at **fixed q** this drifts with
  the scale factor (0.0411 -> 0.0130 -> 0 for `q = 3.15` as `a` goes from 1e-3 to
  1e-1 with `m_nu = 0.05 eV`). It is not a solution of the collisionless
  Boltzmann equation. It describes a species that stayed in kinetic equilibrium
  to the present day; relic neutrinos did not, and never re-thermalise, so they
  keep their ultra-relativistic Fermi-Dirac shape in momentum forever.

There is also an internal inconsistency: in the `energy` branch the code uses
`f0(eps)` for `Gamma` (`transfer_cnb_Gamma()`) but `f0(q)` for the momentum
average weights — two different distributions in one expression.

**Why it nonetheless looks reasonable.** There is a legitimate factor `v^2`
nearby. Keeping the correct `f0(q)` and merely re-expressing the same physical
`Psi_f` in an energy-shift variable rather than a momentum-shift one gives,
since `dln q/dln eps = eps^2/q^2`,

    Gamma_eps = (eps^2/q^2) Gamma_q    =>    N_eps = v^2 N_q    exactly

i.e. a relabelling, not a second prescription. The `energy` option carries that
`p^2/E^2` but *also* replaces `f0(q)` by `f0(eps)`, and that is the step that
breaks.

The `energy` code path is retained because it costs nothing and allows
comparison with the literature, but it should not be used for results.

---

## 3. The two potentials, and the Sachs-Wolfe coefficient

**This section replaces an earlier one that was wrong.** An earlier version of
this implementation claimed the Thesis had omitted a boundary term and used
`C_SW = 3/2`. That was a misreading, and it produced `C_l^SW` twelve times too
large. Both errors are corrected below; keep the massless-limit test of Sec. 5
in place if you touch this code again.

### 3.1 The naming convention is opposite to CLASS

The Thesis writes the metric with `Phi` in the time-time component and `Psi` in
the spatial one (its eq. 3.1.6 and 3.1.7). CLASS and Ma & Bertschinger do the
**opposite**: `psi` is the Newtonian (time-time) potential, `phi` the spatial
curvature perturbation. So

    Thesis Phi  ==  CLASS psi
    Thesis Psi  ==  CLASS phi

The Thesis equations are correct in their own convention. Translating them
symbol-by-symbol into CLASS names — which is what the first version of this code
did — interchanges the two potentials everywhere.

That is not cosmetic. Neutrino anisotropic stress splits them by 14% at neutrino
decoupling (`psi/phi = 0.8596` measured at `k*tau = 4.6e-6`). And the error
hides in the regime where you are most likely to test: for a massless neutrino
`(eps/q)^2 = 1`, the ISW integrand collapses to the symmetric `phi' + psi'` and
the swap changes nothing at all. It only bites for massive neutrinos, where
`(eps/q)^2` can reach `1e2` and it matters which potential it multiplies.

The Boltzmann equation fixes the assignment: the source is
`[phi' - (eps/q) ikmu psi]`, so the factor generated by the integration by parts
belongs to `psi`.

### 3.2 The Sachs-Wolfe coefficient is 1/2, not 3/2

Two terms carry the angular dependence `j_l(k chi_nu(tau_dec))`:

- the adiabatic initial condition, Thesis eq. **(1.4.16)**, `N_0 = -Phi/2` —
  *negative*. In CLASS names that is `Theta_0 = -psi/2`, the standard
  super-horizon relation `delta_nu = -2 psi`. Verified directly:
  `delta_ur/4 = -0.49958 * psi` at `k*tau = 4.6e-6`.
- the lower boundary term from the by-parts step, `+(eps_dec/q)^2 psi(tau_dec)`,
  which is `+psi` to within `a_dec^2 M^2/q^2 ~ 1e-15`.

Their sum is `+psi/2`. The `+1/2` in Thesis eq. (3.2.18)/(3.2.24) is therefore
the **already-combined result**, not a bare initial condition — nothing was
omitted. Reading it as an initial condition and adding the boundary term a
second time double-counts it.

| approach | `C_SW` |
|---|---|
| momentum | `-1/2 + (eps_dec/q)^2` (= 1/2 in practice) |
| energy   | `(1/2) * Gamma_dec/Gamma_0` |

In the massless limit `C_SW * psi = Theta_0 + psi` is exactly the textbook
Sachs-Wolfe amplitude — the same thing CLASS uses for the CMB itself, in the
`switch_sw * g * (delta_g/4 + psi)` term of `perturbations_sources()`.

## 3 bis. Testing

`notebooks/test_cnb_spectra.ipynb` runs twelve checks against the implementation
and prints a pass/fail summary. It uses a classy built from this modified tree and
kept in `classy_cnb/` (git-ignored), deliberately **outside** site-packages so that
any `classy` already installed in the environment is left untouched:

```python
sys.path.insert(0, os.path.abspath("../classy_cnb"))
import classy_cnb as classy
```

That directory is not committed, because the compiled extension is large and
specific to one platform and Python version. To build it from the repository root:

```bash
make libclass.a

mkdir -p classy_cnb && cd classy_cnb
cp ../python/classy.pyx classy_cnb.pyx
cp ../python/cclassy.pxd .
cython --cplus -3 classy_cnb.pyx -I ../python
g++ --std=c++11 -O3 -pthread -fPIC -shared \
    -I$(python -c "import sysconfig;print(sysconfig.get_paths()['include'])") \
    -I$(python -c "import numpy;print(numpy.get_include())") \
    -I../include -I../external/RecfastCLASS -I../external/heating \
    -I../external/HyRec2020 -I../external/Halofit -I../external/HMcode \
    classy_cnb.cpp -o classy_cnb.so ../libclass.a -lm
```

Rerun this after any change to the C sources, or the notebook will silently test
stale code. Alternatively, install the wrapper normally with
`pip install .` from the repository root — the CNB additions are in
`python/classy.pyx`, so a standard `classy` build picks them up too.

The notebook found two real defects that the command-line checks had missed: the
`phi` anchoring described in Sec. 2.3, and the momentum-average convergence
documented in Sec. 5.

---

## 4. Files modified

All new code is gated behind `ppt->has_cl_cnb` / `ptr->has_cnb`; when `cnbCl` is
not requested, no CNB code runs.

### `include/perturbations.h`
- `enum cnb_approaches {cnb_momentum_approach, cnb_energy_approach}`
- `#define _CNB_NUM_MAX_ 16`
- `has_cl_cnb`, `has_cnb_sw`, `has_cnb_isw`
- `cnb_approach`, `cnb_num`, `cnb_M[]`, `cnb_m_in_eV[]`, `cnb_z_dec`, `cnb_q_size`
- prototype of `perturbations_cnb_indices()`

### `source/perturbations.c`
- `perturbations_cnb_indices()` — resolves the neutrino masses (explicit list, or
  the `ncdm` species) into the dimensionless `M = m_nu/T_nu0`, and rejects
  isocurvature initial conditions.
- `perturbations_prepend_radiation_era_sampling()` — extends `ppt->tau_sampling`
  back into the radiation era (Sec. 2.3).
- sets `has_cmb = _TRUE_` when the CNB is requested. This does **not** compute any
  photon source; it selects the time-sampling strategy, which the CNB needs for the
  same reason the CMB does (early start, step set by the metric variation rate).
- requests the `phi`, `phi_prime` and `psi` source functions.

### `source/input.c`
- `cnbCl` added to the `output` keywords.
- reads `cnb_approach`, `cnb_contributions`, `cnb_m_nu`, `cnb_z_dec`, `cnb_q_sampling`,
  `cnb_q_size`, `cnb_q_phase_step`, `cnb_q_panel_order`, `cnb_q_weight_tol`, `cnb_q_max`,
  with defaults set in `input_default_params()`.

### `include/transfer.h`
- transfer types `index_tt_cnb_sw`, `index_tt_cnb_isw` (one per neutrino mass).
- `cnb_psi_dec` alongside `cnb_phi_dec`: the SW term needs **psi** (Sec. 3.1).
- the free-streaming tables, flattened as `cnb_q_offset[index_cnb] + index_q_nu`
  (each species has its own `cnb_q_num[index_cnb]` momenta under the chi sampling):
  `cnb_q`, `cnb_weight`, `cnb_chi`, `cnb_eps`, `cnb_deps`, `cnb_Gamma`,
  `cnb_dGamma`, `cnb_chi_dec`, `cnb_eps_dec`, `cnb_Gamma_dec`, `cnb_Gamma_0`,
  plus `cnb_tau_dec` and `cnb_phi_dec`.
- workspace arrays `cnb_phi`, `cnb_phi_prime`, `cnb_psi`, `cnb_psi_prime`,
  `cnb_ddpsi`.
- `transfer_compute_for_each_l()` gained an optional output pointer.

### `source/transfer.c`
New routines:
- `transfer_cnb_distribution()` — the Fermi-Dirac `f0(q)`.
- `transfer_cnb_Gamma()` — `Gamma = dln f0/dln eps` and its derivative, written
  with the logistic function to avoid overflowing `exp()` at large `eps`.
- `transfer_cnb_precompute()` — Gauss-Laguerre momentum grid, normalised weights
  `q^2 eps f0(q)` of eq. (3.2.38), and the tables of `eps(tau)`, `chi_nu(tau)`,
  `Gamma(tau)`. `chi_nu` is integrated backwards from today on
  `ppt->tau_sampling`, then continued on a logarithmic grid down to `tau_dec` to
  obtain the distance to the neutrino last scattering surface.
- `transfer_cnb_free()`.
- `transfer_cnb_interpolate_potentials()` — interpolates `phi`, `phi'`, `psi` at
  the current `k`. `phi'` is taken exactly from the perturbation module; `psi'` is
  not available there and is obtained by splining `psi` in conformal time on the
  ISW-resolving sampling.
- `transfer_cnb_source()` — the sources of eq. (3.2.24) and (3.2.32), with the
  potentials mapped into CLASS names and the SW coefficient of Sec. 3 above. The SW term is returned as a
  single-point (Dirac) source at `tau_dec`, which `transfer_integrate()` already
  handles.
- `transfer_cnb_compute_for_each_q()` — the loop over neutrino momenta with the
  on-the-fly momentum average.

Existing routines touched:
- `transfer_init()` — sets `has_cnb` and the CNB sizes before `transfer_indices()`,
  and calls `transfer_cnb_precompute()`.
- `transfer_indices()`, `transfer_get_l_list()`, `transfer_source_tau_size()`,
  `transfer_get_source_correspondence()` — CNB branches.
- `transfer_compute_for_each_q()` — dispatches the CNB types to the dedicated
  routine and skips them in the normal type loop.
- `transfer_compute_for_each_l()` — returns the value through the new output
  pointer when it is non-NULL, instead of storing it (the single pre-existing call
  site passes NULL and is unaffected).
- `transfer_workspace_init()` / `_free()`, `transfer_free()` — allocation.

### `include/harmonic.h`, `source/harmonic.c`
- `has_nn`, `has_tn`, `index_ct_nn`, `index_ct_tn`, `n_size`, `cnb_m_in_eV[]`.
- in `harmonic_compute_cl()`, the total CNB transfer function is the sum of its SW
  and ISW parts (as the temperature one is the sum of `t0+t1+t2`), and the
  integrands of eq. (3.2.37) are filled for `nunu` and `Tnu`.

### `source/output.c`
- column titles for the CNB spectra, labelled by neutrino mass.

### `include/precisions.h`, `include/quadrature.h`
- `cnb_radiation_era_tau_size` (default 40) and
  `_MAX_NUMBER_OF_LAGUERRE_POINTS_` (80).

### `python/cclassy.pxd`, `python/classy.pyx`
- the new harmonic fields, and a `cnb_cl()` method returning a dict with `ell`,
  `m_nu`, `nn` and `tn` (the last two of shape `(n_mass, lmax+1)`).

### `explanatory.ini`
- documentation of `cnbCl` and of the `cnb_*` parameters.

---

## 5. Validation

Run with the parameters of Table 4.1 (`cnb.ini`).

**Massless limit against the textbook Sachs-Wolfe amplitude.** This is the test
that anchors everything, because both of its ingredients come from the
perturbation module and none from the CNB code:

    Theta_0 + psi = -0.300429 + 0.601368 = 0.300939
    => l(l+1)C_l/2pi = A_s * 0.300939^2 = 1.9064e-10

At `m_nu = 1e-5 eV`, `l = 2`, `n_s = 1` the code returns **1.9021e-10**
(momentum) and **1.9010e-10** (energy), i.e. 0.2%. The residual is the small
`chi_nu` spread that survives even at this mass, not a normalisation error. The
interchanged-potential version of the code fails this test by a factor 3.5 in
amplitude; nothing else we tried caught it.

**Sachs-Wolfe coefficient at finite mass.** Implied
`C_SW = sqrt(Dl/A_s)/psi_dec` at `l = 2`, against eq. (C_SW) above:

| | implied | expected |
|---|---|---|
| momentum, all masses | 0.4999 / 0.4996 / 0.4700 / 0.4415 | 0.5 |
| energy, `m = 1e-5 eV` | 0.4998 | 0.4998 |
| energy, `m = 1e-4 eV` | 0.4859 | 0.4860 |

for `m_nu = 1e-5, 1e-4, 1e-3, 4e-3 eV`. Agreement is 0.02% wherever the plateau
is genuinely flat; the energy-approach number is the stronger check, since the
expected value moves with mass and the code tracks it. Where the plateau drifts
(`m_nu >= 1e-3 eV`, drift 30-54% between `l = 2` and `l = 10`) the closed form
itself no longer applies and the disagreement measures that, not a code error.

**Neutrino last scattering surface (Fig. 4.4).** Range of `chi_nu(tau_dec)` over
the sampled momenta, against the CMB value of 14159 Mpc:

| `m_nu` | `chi_nu(tau_dec)` |
|---|---|
| 1e-5 eV | 13744 – 14159 Mpc |
| 1e-4 eV | 8582 – 14159 Mpc |
| 1e-3 eV | 3149 – 14141 Mpc |
| 4e-3 eV | 1532 – 13896 Mpc |

The spread also shows up directly in the SW spectrum, which stays flat to
`l ~ 100` for the lightest mass but starts falling around `l ~ 5` for the heaviest,
as the different momenta decohere.

**Momentum vs energy approach.** They agree to 3e-4 at `m_nu = 1e-5 eV`
(Fig. 4.6 caption) and diverge as the mass grows (Sec. 4.1).

**Sign of the cross-correlation.** At `m_nu >= 1e-3 eV` the momentum approach gives
a negative `C_l^{Tnu}` at low `l`, suppressed above `l ~ 20`, while the energy
approach is positive and strongly suppressed — as described in Sec. 4.3.

**Convergence of the momentum average.** With the default Gauss-Laguerre sampling
this is the weak point of the implementation. For a massive neutrino `chi_nu(q)`
varies steeply at small `q`, so `Delta_l(q)` oscillates rapidly in `q`, and a rule
tuned to the smooth weight `exp(-q)` cannot follow it. Comparing Laguerre runs with
each other (30 vs 60 nodes) badly *understates* the problem, because the 60-node
reference is itself far from converged.

`cnb_q_sampling = chi` fixes it. `transfer_cnb_q_sampling_chi()` places composite
Gauss-Legendre panels in `q` so that the Bessel phase `k chi_nu(tau,q)` advances by
`cnb_q_phase_step` per panel. The rate is evaluated at `k = min(l_max/chi_nu, k_max)`
and maximised over all sampled times; late-time `chi_nu` of a non-relativistic
neutrino goes like `q`, twice as steep in `ln q` as `chi_nu(tau_dec)`. It is scaled
down where the effective weight `q * q^2 eps f0 * (1+M^2/q^2)` is below
`cnb_q_weight_tol` of its maximum. Extra edges keep `dq <= 1`, `dlnq <= 0.5`.

Checks made on 2026-10-16, `l_max = 200`, `k_max_tau0_over_l_max = 12`:

- **Laguerre path unchanged:** bit-identical to the previous code, full and ISW-only.
- **Self-convergence** (relative to `cnb_q_phase_step = 0.5`/`1`), worst case over `l <= 200`:

  | `m_nu` | term | Laguerre 30 | Laguerre 80 | chi, step 8 | chi, step 4 |
  |---|---|---|---|---|---|
  | 1e-3 eV | SW  | 2.7 | 1.3 | 1.9e-3 | 3e-7 |
  | 1e-3 eV | ISW | 0.71 | 0.26 | 6e-5 | 4e-7 |
  | 4e-3 eV | SW  | 6.9 | 3.8 | 4.6e-3 | 2e-6 |
  | 4e-3 eV | ISW | 2.1 | 0.94 | 1e-4 | 4e-7 |

- **Independent reference for SW:** with `psi(tau_dec)` k-independent,
  `C_l^SW = 4pi int dk/k P(k) [psi_dec <S(q) j_l(k chi_dec(q))>_W]^2`, computed in numpy on
  4000-16000 `sqrt(q)` points with `chi_dec` from the CLASS background table (script kept
  outside the repo). The chi-sampled CLASS SW spectrum agrees with it to `<= 0.33%` at
  1e-3 eV and `<= 0.09%` at 4e-3 eV for all `l <= 200` (the same floor as at 1e-5 eV,
  i.e. CLASS's own k/l sampling). Laguerre-30 is off by +29% / +275% (1e-3 eV, `l` = 50 / 200)
  and +108% / +694% (4e-3 eV).
- **Light masses:** chi agrees with Laguerre-80 to `< 5e-4` for 1e-4 eV and `< 1e-6` for 1e-5 eV.

Cost: 180-200 momenta for a light species, ~1300 for 1e-3 or 4e-3 eV at step 4. Each
ISW momentum costs ~0.06 s with 643 sampled times, so the four-mass run above takes
222 s against 13 s for Laguerre-30. ISW alone is converged to 0.1% already at step 16,
so separate SW/ISW momentum grids would cut that by about three.

**Non-interference.** Enabling `cnbCl` changes the `TT` spectrum by at most
`1e-4` in relative terms, through the extended source sampling of Sec. 2.3.

---

## 6. Open points

**Resolved: the cross-correlation maximum near `m_nu ~ 1e-4 eV`.** Sec. 4.3 of
the Thesis is right, but the statement is multipole dependent. At `l = 2`
`C_l^Tnu` falls monotonically with mass; the maximum only appears above
`l ~ 25`, growing from +2% at `l = 24` to +13% at `l = 95`. This code puts it at
`m_nu ~ 9e-5 eV`. The Thesis quotes `1e-4 eV` because `single_mass_prova.py`
line 94 evaluates the Fermi-Dirac distribution at `T = 2e-4 eV` rather than
`T_nu0 = 1.682e-4 eV`; since everything depends on `M = m_nu/T_nu0`, that
rescales every mass by 1.19, and `1e-4/1.19 = 8.4e-5 eV`. Cross-checked against
that script's own saved outputs in `~/neutrini/cross/`, which show the same
maximum at the same multipoles.

**Still open: the momentum-approach ISW diverges as `q -> 0` for heavy masses**,
since its source carries `eps^2/q^2 = 1/v^2`. This is what the equations say
rather than a bug, and is presumably why Sec. 4.1 of the Thesis presents the
energy approach as complementary. The integral converges, and the chi sampling
of Sec. 5 resolves it; under the Laguerre sampling it makes heavy-mass results
sensitive to `cnb_q_size`.

**Resolved: the momentum quadrature** for `m_nu >= 1e-3 eV`. The Gauss-Laguerre
rule is not converged above `l ~ 20` for heavy masses, and `cnb_q_sampling = chi`
(Sec. 5) is the fix: a rule adapted to `chi_nu(q)` rather than to `exp(-q)`.
Laguerre stays the default because it is ~15x cheaper and exact enough for
`m_nu <= 1e-4 eV`; use `chi` for anything heavier.
