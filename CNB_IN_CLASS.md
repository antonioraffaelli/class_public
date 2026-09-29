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

from eq. (3.2.24) in the **momentum approach** and eq. (3.2.32) in the **energy
approach**, momentum-averaged following eq. (3.2.38) and integrated over `k`
following eq. (3.2.37).

The Sachs-Wolfe and integrated Sachs-Wolfe terms are kept as two separate
transfer types internally, so they can be switched on and off independently
(`cnb_contributions`) in order to reproduce the component-by-component figures
of Sec. 4.3.

### Usage

```ini
output = tCl,cnbCl

cnb_approach      = momentum          # or 'energy'
cnb_contributions = sw, isw           # optional, default is both
cnb_m_nu          = 1.e-5, 1.e-4, 1.e-3, 4.e-3   # eV; default: the ncdm masses
cnb_z_dec         = 6.e9              # neutrino decoupling redshift
cnb_q_size        = 15                # Gauss-Laguerre nodes in q
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

### 2.3 phi at neutrino decoupling

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

## 3. Correction to eq. (3.2.18) / (3.2.30)

The integration by parts of eq. (3.2.17) and (3.2.29) produces a boundary term at
the lower limit `eta_in`, which is written out in those equations but does **not**
appear in eq. (3.2.18) and (3.2.30), nor in the Sachs-Wolfe parts of (3.2.24) and
(3.2.32). It has been restored here.

Since that term is proportional to `phi(eta_in) exp(ikmu chi(eta_in,eta_in))` with
`chi(eta_in,eta_in) = 0`, it carries exactly the same angular dependence
`j_l(k chi_nu(tau_dec))` as the Sachs-Wolfe term coming from the initial condition
`N(eta_in) = phi(eta_in)/2`. It therefore renormalises the SW coefficient.

The two approaches get **different** coefficients, because the energy-approach
Boltzmann equation (3.2.28) carries an extra factor `p^2/E^2` which exactly cancels
the `E^2/p^2` generated by the by-parts step:

| approach | boundary term | resulting SW coefficient |
|---|---|---|
| momentum, (3.2.17) | `+(E^2/p^2) phi(eta_in) Gamma` | `1/2 + (eps_dec/q)^2` |
| energy, (3.2.29)   | `+phi(eta_in) Gamma(eta_in)`   | `(3/2) * Gamma_dec/Gamma_0` |

In the momentum approach `(eps_dec/q)^2` differs from one only by
`a(tau_dec)^2 M^2 / q^2 ~ 1e-15`, so the coefficient is `3/2` for all practical
purposes; it is nevertheless kept exact in the code, since it costs nothing.

Compared with the Thesis expressions (coefficient `1/2`), this raises the SW
contribution to `C_l^{nunu}` by a factor 9 in the momentum approach, and changes
the SW/ISW balance — hence the sign structure of `C_l^{Tnu}` at low multipoles.

---

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
- reads `cnb_approach`, `cnb_contributions`, `cnb_m_nu`, `cnb_z_dec`, `cnb_q_size`,
  with defaults set in `input_default_params()`.

### `include/transfer.h`
- transfer types `index_tt_cnb_sw`, `index_tt_cnb_isw` (one per neutrino mass).
- the free-streaming tables, flattened as `index_cnb*cnb_q_size + index_q_nu`:
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
- `transfer_cnb_source()` — the sources of eq. (3.2.24) and (3.2.32), including
  the corrected SW coefficient of Sec. 3 above. The SW term is returned as a
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
- documentation of `cnbCl` and of the five `cnb_*` parameters.

---

## 5. Validation

Run with the parameters of Table 4.1 (`cnb.ini`).

**Sachs-Wolfe term against its closed form.** With `n_s = 1` and a single
last-scattering distance, `int dk/k j_l(k chi)^2 = 1/(2l(l+1))` gives exactly

```
l(l+1) C_l / 2pi = A_s (C_SW phi_dec)^2 ,   independent of l
```

The code returns a flat plateau, with the implied `C_SW` equal to

| | implied `C_SW` | expected |
|---|---|---|
| momentum | 1.4998 | `1/2 + (eps_dec/q)^2 = 1.5` |
| energy, `m = 1e-5 eV` | 1.4993 | `(3/2) Gamma_dec/Gamma_0 -> 1.5` |
| energy, `m = 1e-3 eV` | 0.631 | decreasing with mass |

i.e. the full chain — Bessel convolution against the neutrino radial coordinate,
Dirac source handling, momentum average, primordial integration and harmonic
normalisation — is exact to ~0.02%.

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

**Convergence of the momentum average.** This is the weak point of the
implementation. Relative to a 60-node run, the default `cnb_q_size = 30` gives:

| `m_nu` | `l = 2` | `l = 10` | `l = 50` | `l = 150` |
|---|---|---|---|---|
| 1e-5 eV | 2e-8 | 7e-8 | 3e-6 | 3e-7 |
| 1e-4 eV | 8e-7 | 2e-5 | 4e-4 | 1e-3 |
| 1e-3 eV | 6e-4 | 4e-3 | 2e-1 | 6e-1 |
| 4e-3 eV | 7e-5 | 5e-4 | 4e-1 | 7e-1 |

i.e. excellent for `m_nu <= 1e-4 eV` at all `l`, but **not converged** for
`m_nu >= 1e-3 eV` above `l ~ 20`, and raising `cnb_q_size` helps only slowly. The
cause is that for a massive neutrino `chi_nu(q)` varies steeply at small `q`, so
`Delta_l(q)` oscillates rapidly in `q` and a Gauss-Laguerre rule tuned to the
smooth weight `exp(-q)` resolves it poorly. Fixing this properly would mean
replacing the momentum quadrature with one adapted to `chi_nu(q)` (for instance
sampling uniformly in `chi_nu` rather than in `q`). Until then, heavy-mass,
high-`l` results should be treated as indicative. `cnb_q_size` is capped at 80,
beyond which the Gauss-Laguerre weights underflow.

**Non-interference.** Enabling `cnbCl` changes the `TT` spectrum by at most
`1e-4` in relative terms, through the extended source sampling of Sec. 2.3.

---

## 6. Open point

I could not reproduce the statement in Sec. 4.3 that the cross-correlation has
**maximum amplitude around `m_nu ~ 1e-4 eV`**. The code gives an essentially flat
`C_l^{Tnu}` between `1e-5` and `1e-4 eV` (4.75e-10 vs 4.70e-10 at `l = 2`), with no
peak. The `chi_nu` ranges above are consistent with that: at `1e-5 eV` *all*
neutrinos already share the CMB last scattering surface, so the overlap should be
maximal there rather than at `1e-4 eV`. Worth cross-checking against the original
Python script.

Note also that the momentum approach diverges as `q -> 0` for heavy masses, since
its ISW source carries a factor `E^2/p^2 = 1/v^2`. This is what the equations say,
and is presumably why Sec. 4.1 presents the energy approach as complementary, but
it does make the heavy-mass momentum-approach results sensitive to `cnb_q_size`.
