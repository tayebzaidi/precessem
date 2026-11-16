# Anatomical human body models — availability and licences (reference note)

Collected 2026-09-28 (interactive session, operator request) for the
human-scale and implant work (Phase 7 implants, `WF-9`, `WF-10` / `WF-11`,
`POST-9` and the worst-case implant-heating research direction). Public
information only; the licence terms were read from each provider's page or a
search summary on that date — **re-check the terms before adopting a model**,
they change (the Visible Human Project's did in 2019).

## Summary

Several anatomical models are free for research; **almost none is open source
in the strict sense** — nearly all carry non-commercial (NC) and/or
no-derivatives (ND) terms. XCAT and the IT'IS Virtual Population are the
field's standards and are licensed / paid.

| Model | Coverage | Format | Licence (as read 2026-09-28) | Fit for this platform |
|---|---|---|---|---|
| [AustinMan / AustinWoman](https://web.corral.tacc.utexas.edu/AustinManEMVoxels/AustinMan/index.html) (UT Austin, from the Visible Human data) | Whole body, male and female | Voxels, built for EM | [CC BY-NC-ND 3.0](http://web.corral.tacc.utexas.edu/AustinManEMVoxels/AustinMan/citing_the_model/index.html) — free for research, teaching, non-commercial use | Best free whole-body option for a 1.5 T body-coil study. Needs voxel → tetrahedral conversion. **ND: do not redistribute meshes derived from it**; using them in-house is fine. |
| [MIDA](https://itis.swiss/virtual-population/regional-human-models/mida-model) (IT'IS / FDA) | Head and neck, 153 structures, 0.5 mm | Segmentation / surfaces | [Free of charge](https://cdrh-rst.fda.gov/mida-multimodal-imaging-based-model-human-head-and-neck) under a privacy licence agreement (handling fee only) | The natural choice for DBS lead heating: deep brain structures; already used in DBS / neurostimulation safety studies ([PLOS One](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0124126)). |
| [SimNIBS "Ernie" / "Ernie Extended"](https://simnibs.github.io/simnibs/build/html/dataset.html) | Head; head + upper shoulders ([paper](https://direct.mit.edu/imag/article/doi/10.1162/imag_a_00379/125235/Head-and-shoulders-The-impact-of-an-extended-head)) | Tetrahedral FEM mesh (Gmsh) | CC BY-NC 4.0 | Least effort into FEniCSx — already our mesh type. Tissue labels carry low-frequency conductivities; RF properties at 64 / 128 MHz must be re-assigned. Shoulders help coil loading. |
| [Visible Human Project](https://www.nlm.nih.gov/research/visible/getting_data.html) (NLM) | Whole body; raw cryosection, CT, MR | Images, unsegmented | Since July 2019 [no licence agreement](https://en.wikipedia.org/wiki/Visible_Human_Project), only terms requiring acknowledgement of NLM | Most permissive, but raw — segmentation is the work (what AustinMan did). |
| [ICRP 145 mesh-type reference phantoms](https://www.icrp.org/publication.asp?id=ICRP+Publication+145) | Whole body, adult reference male / female | Polygon / tetrahedral mesh | Publication [free to access since 2023](https://www.icrp.org/page.asp?id=491); **licence of the data files themselves not confirmed** | Mesh-native, well documented, built for radiation dosimetry; micrometre-scale target layers would need simplifying for EM. |
| [IT'IS Virtual Population](https://itis.swiss/virtual-population/virtual-population/overview) (Duke, Ella, …) | Whole body, many subjects | Voxels; CAD versions for ANSYS / CST | Paid, per model ([Duke](https://itis.swiss/virtual-population/virtual-population/vip2/duke/), [Ella](https://itis.swiss/virtual-population/virtual-population/vip2/ella/ella-v2-0)) | The field's standard; cost and terms are the barrier. |
| XCAT (Duke) | Whole body, parameterised anatomy | NURBS | Restrictive licence (operator, 2026-09-28) | Its patient-variability parameterisation is the attraction and the reason for the licence. |

## Recommendation (2026-09-28)

- **Head / DBS lead (the `POST-9` pitch, aims 1–2):** SimNIBS **Ernie
  Extended** first (already tetrahedral, free for research, extends past the
  head), **MIDA** when anatomical detail matters.
- **Whole body at 1.5 T:** **AustinMan / AustinWoman**, respecting ND.
- **Tissue properties** at 64 / 128 MHz come from elsewhere for every model —
  the usual source is the IT'IS tissue-property database (**terms not checked
  here; verify**). Gelled-saline phantom values stay in §11.

## Repository rule (proposed, to be enacted when a model is adopted)

Every model above carries NC and/or ND terms, or a privacy agreement. Treat
model files and any mesh derived from them like the AED results: **gitignored,
never in a tracked file or a commit, never pushed**; only the loader code, the
tissue-property mapping and qualitative or relative results go public. The
existing NDA rule for the research coil's sizing (`GEO-33`) is the precedent.

## Sources

- [AustinMan Electromagnetic Voxels](https://web.corral.tacc.utexas.edu/AustinManEMVoxels/AustinMan/index.html) · [Citing the model / licence](http://web.corral.tacc.utexas.edu/AustinManEMVoxels/AustinMan/citing_the_model/index.html) · [paper (ResearchGate)](https://www.researchgate.net/publication/239938054_AustinMan_and_AustinWoman_High_fidelity_reproducible_and_open-source_electromagnetic_voxel_models)
- [MIDA — IT'IS Foundation](https://itis.swiss/virtual-population/regional-human-models/mida-model) · [MIDA — FDA CDRH](https://cdrh-rst.fda.gov/mida-multimodal-imaging-based-model-human-head-and-neck) · [MIDA paper — PLOS One](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0124126)
- [SimNIBS datasets](https://simnibs.github.io/simnibs/build/html/dataset.html) · [Ernie Extended — Imaging Neuroscience](https://direct.mit.edu/imag/article/doi/10.1162/imag_a_00379/125235/Head-and-shoulders-The-impact-of-an-extended-head)
- [Visible Human Project — getting the data (NLM)](https://www.nlm.nih.gov/research/visible/getting_data.html) · [Wikipedia (2019 licence change)](https://en.wikipedia.org/wiki/Visible_Human_Project)
- [ICRP Publication 145](https://www.icrp.org/publication.asp?id=ICRP+Publication+145) · [ICRP — Publication 145 now available](https://www.icrp.org/page.asp?id=491)
- [IT'IS Virtual Population overview](https://itis.swiss/virtual-population/virtual-population/overview) · [Duke](https://itis.swiss/virtual-population/virtual-population/vip2/duke/) · [Ella](https://itis.swiss/virtual-population/virtual-population/vip2/ella/ella-v2-0)
