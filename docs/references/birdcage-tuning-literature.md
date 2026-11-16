# Birdcage tuning — how the field does it (literature note)

Collected 2026-09-28 (interactive session, operator request) while deciding how
to choose the ring capacitor for the human-scale birdcage (`WF-10`, §10 chain
H8′). Public literature only — tracked, citable from the plan. **Read depth is
stated per entry**: "read" means the page itself was fetched and read;
"abstract / snippet" means only a search-result summary was seen, so check the
paper before leaning on a detail.

## The finding in one paragraph

The standard workflow for tuning a birdcage (or any multi-capacitor MRI coil)
is **EM / RF-circuit co-simulation**: every capacitor site becomes a 50 Ω port
in *one* driven full-wave solve (HFSS, CST, …); capacitors are then chosen in
the circuit domain, and fields for any tuning are rebuilt by superposition.
The capacitor is chosen **for a mode** — either by construction (a circuit
model that targets the k = 1 current pattern) or by an optimiser whose cost
includes field homogeneity, with the mode then confirmed from the
capacitor-current pattern (≈ sinusoidal around the end ring for k = 1).
**Full-wave eigenmode solves are used to analyse the mode spectrum**, less to
tune a loaded, lossy coil (the problem is complex and, with lumped
capacitors, frequency-dependent). Lumped-element models alone break down at
high field, which is why co-simulation replaced them.

**What this repo does with it** (operator ruling 2026-09-28, §10 H8′): our
`PORT-15` circuit layer *is* co-simulation. On a C16-symmetric coil the port
matrix splits exactly by azimuthal harmonic, so the k = 1 mode is tuned in
its own block (`C = 1 / (ω Im λ₁)`), checked against an in-model solve and a
k = ±1 current share ≥ 0.9 — the "choose for a mode, confirm by current
pattern" step of the literature, done with linear algebra instead of an
optimiser. `TH-17` T5a/T5b did the same on F-small and confirmed it against a
true eigensolve.

## Sources

| Source | What it contributes | Read depth |
|---|---|---|
| Kozlov M., Turner R., "Fast MRI coil analysis based on 3-D electromagnetic and RF circuit co-simulation", *J. Magn. Reson.* 200(1):147–152, 2009 | Introduces co-simulation for MRI coils: capacitors and feeds replaced by 50 Ω ports in one 3-D EM run, tuning done in the circuit domain, fields rebuilt per tuning. | Snippet only — seen as cited by the next entry; bibliographic details from memory, **verify volume/pages/DOI** |
| [Effects of tuning conditions on near field of MRI transmit birdcage coil at 64 MHz (ResearchGate)](https://www.researchgate.net/publication/309341951_Effects_of_tuning_conditions_on_near_field_of_MRI_transmit_birdcage_coil_at_64_MHz) · [PubMed 28269678](https://pubmed.ncbi.nlm.nih.gov/28269678/) | 64 MHz birdcage near field vs tuning condition, using co-simulation; describes the ports-at-capacitors workflow. Relevant to `WF-7` / `WF-10` (same frequency, same question of what tuning does to the field). | Snippet |
| [dezanche/Elliptical_Birdcage_Coil (GitHub)](https://github.com/dezanche/Elliptical_Birdcage_Coil/) | A worked, public example: HFSS 2022R1, 36 lumped 50 Ω ports at the capacitor gaps, discrete sweep 100–150 MHz / 2 MHz; CoSimPy circuit layer; Nelder–Mead over 7 symmetry-reduced capacitor values minimising B₁⁺ standard deviation + 100 × the imaginary admittance at 128 MHz; mode confirmed by (co)sinusoidal end-ring capacitor currents and by re-running HFSS with the optimal capacitors in place. The closest analogue to our `PORT-15` + `WF-10` pipeline. | Read |
| Gurler N., Ider Y. Z., "Numerical methods and software tools for simulation, design, and resonant mode analysis of radio frequency birdcage coils used in MRI", *Concepts Magn. Reson. B* 45(1), 2015 — [Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/cmr.b.21279) · [ResearchGate](https://www.researchgate.net/publication/277724665_Numerical_methods_and_software_tools_for_simulation_design_and_resonant_mode_analysis_of_radio_frequency_birdcage_coils_used_in_MRI_SIMULATION_DESIGN_RESONANT_MODE_ANALYSIS_OF_RF_COILS) | FEM time-harmonic simulation and **resonant-mode analysis** of birdcages; capacitor estimation by maximising port-impedance magnitude at the target frequency while minimising B₁ variance; HFSS re-runs with the optimal capacitors agree with co-simulation. The eigen/mode-analysis side of the field. | Abstract / snippet |
| Chin C.-L., Collins C. M., Li S., Dardzinski B. J., Smith M. B., "BirdcageBuilder: design of specified-geometry birdcage coils with desired current pattern and resonant frequency", *Concepts Magn. Reson.* 15(2):156–163, 2002 — [Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/cmr.10030) | Analytic circuit design: capacitance computed from a **chosen current pattern** (k = 1), calculated inductances and the target frequency — high-pass, low-pass and band-pass. The "select the mode first, then solve for C" principle H8′ adopts. Lumped model; not accurate at high field. | Abstract / snippet |
| [OpenBirdcageBuilder (GitHub)](https://github.com/dbnyu/OpenBirdcageBuilder) · [Open Source Imaging page](https://www.opensourceimaging.org/project/open-birdcagebuilder/) | Open-source re-implementation of BirdcageBuilder (low-/high-/band-pass, tubular or flat legs and rings). A cheap cross-check on a first-guess `C` for F-human or the 32-leg `WF-11` coil. | Snippet |
| Kazemivalipour E. et al., "Enhancing fine-tuning efficiency and design optimization of an eight-channel 3T transmit array via equivalent circuit modeling and Eigenmode analysis", *Med. Phys.* 2025 — [PMC11972053](https://pmc.ncbi.nlm.nih.gov/articles/PMC11972053/) · [Wiley](https://aapm.onlinelibrary.wiley.com/doi/full/10.1002/mp.17612) | 3 T transmit array tuned with an equivalent-circuit model plus eigenmode analysis. **Whether "eigenmode" there means a full-wave eigensolver or an eigen-decomposition of a coupling/impedance matrix was not established** (the page returned a CAPTCHA). Worth reading before `WF-11` (3 T). | Title only |

## Pointers into this repo

- `PORT-15` (§7) — the circuit layer (termination reduction, `tuning_sweep`).
- `TH-17` (§7; §10 chain T5a/T5b) — circulant reduction `z_k(f)` and tuning by
  the eigen route on F-small, confirmed by an eigensolve.
- `WF-10` (§7; §10 chain H8 / H8′ / H9) — the human-scale tuning this note was
  written for.
- `WF-11` (§7; §10 H11) — the 3 T, 32-leg coil; the 128 MHz literature above
  (elliptical birdcage, Kazemivalipour) is its nearest prior work.
