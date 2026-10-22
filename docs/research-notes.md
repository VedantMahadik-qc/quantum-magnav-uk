# Phase 0 research notes (2026-10-22)

Background research done before writing any code. Numbers marked **[derived]** are my own
arithmetic from the cited figures; everything else is quoted from the linked source.

## 1. Prior work

| Work | What it did | Overlap with this project |
|---|---|---|
| Canciani & Raquet, *Magnetic Anomaly Navigation Accuracy with Respect to Map Quality and Altitude*, ION ITM 2016 ([link](https://www.ion.org/publications/pdf.cfm?articleID=13415)) | 1,500 km flight test: ~3.2 km error using existing North American maps; simulation suggests ~160 m with an improved map | Map quality + altitude as dominant error sources. Same question, North America |
| Canciani & Raquet 2016 (as cited in Sengupta 2025) | Up to 13 m accuracy using a marginalized particle filter, high-quality map, low altitude | Particle-filter MagNav is established |
| MagNav.jl, MIT–USAF AI Accelerator ([GitHub](https://github.com/MIT-AI-Accelerator/MagNav.jl)); Gnadt et al., *Signal Enhancement for MagNav Challenge Problem* ([arXiv:2007.12158](https://arxiv.org/abs/2007.12158)); Gnadt PhD thesis 2022 ([MIT](https://dspace.mit.edu/handle/1721.1/145137)) | Julia package: compensation (Tolles–Lawson, NN), maps, EKF/MPF navigation. Challenge data: Cessna 208B near Ottawa, 5 scalar (Cs) + vector fluxgates, 400 m / 800 m altitudes | Reference implementation and dataset. Focus is calibration/compensation |
| Gupta, Sengupta, Phernetton, Sosanya, *Lower Bounds on MagNav Performance as a Function of Magnetic Anomaly Map Quality*, DASC 2024 ([IEEE](https://ieeexplore.ieee.org/document/10749342/)) | CRLB heatmap over a 250 × 400 km region of the North American map (NAMAM); CRLB correlates with gradient magnitude | **Closest to the "navigability map" idea**, but uses CRLB (not a full filter) and covers North America |
| Sengupta, *A Horizontal Accuracy Metric for Magnetic Navigation*, NAVIGATION 72(4), 2025 ([doi:10.33012/navi.717](https://navi.ion.org/content/72/4/navi.717)) | Error budget = sensor variance + map error + calibration error; calibration "can often dominate"; ~5 nT achievable (citing Gnadt 2022b), 25 nT used as conservative (MagNav challenge) | Source for calibration range. Supports the hypothesis qualitatively |
| *Global Uncertainty-Aware Planning for MagNav* (2024, [ResearchGate](https://www.researchgate.net/publication/384075217)); Kang et al., *Optimal Trajectory Generation for Improved MagNav*, 2026 ([arXiv:2609.27553](https://arxiv.org/abs/2609.27553)) | Entropy/CRLB-based path planning, particle filter | Information maps used for routing, not national surveys |
| Muradoglu et al. (Q-CTRL), *Quantum-assured MagNav…*, 2025 ([arXiv:2504.08167](https://arxiv.org/abs/2504.08167)) | Field trials with quantum magnetometers, best 22 m (0.006% of distance), up to 19,000 ft | State of the art with real hardware |
| SandboxAQ AQNav + Acubed ([GPS World](https://www.gpsworld.com/sandboxaq-and-acubed-advance-magnetic-navigation/)) | 150+ flight-hours over NAMAM; RNP2 100% of time, RNP1 95% | Operational evidence, North America |
| Angappan et al., *Towards unified Geophysical Data Requirements for MagNav*, 2026 ([arXiv:2606.14623](https://arxiv.org/abs/2606.14623)) | Position paper on map standards/uncertainty | Map-quality framing |
| UK activity: Royal Navy + Univ. Sussex quantum magnetometer trials ([Royal Navy, Oct 2025](https://www.royalnavy.mod.uk/news/2025/october/15/20251015-rn-university-of-sussex-quantum-trials)) | Sea trials of quantum magnetometers for nav | UK interest, no published navigability map found |

### Novelty assessment (honest)
- **Not novel:** particle-filter MagNav, the error-source taxonomy, the conclusion that map and
  calibration matter more than sensor noise (Canciani 2016, Sengupta 2025), or regional
  information/CRLB heatmaps (Gupta 2024).
- **Not found anywhere public:** a UK-wide navigability map built on the BGS aeromagnetic data, or a
  UK error budget. (Absence of evidence: searched arXiv, ION/NAVI, IEEE and general web.
  Defence work may exist unpublished.)
- **Defensible contribution:** an open, reproducible UK case study that (a) uses a full particle
  filter rather than CRLB alone, so it captures multimodal/divergent failures that a CRLB misses;
  (b) sweeps each error source independently with cited sensor parameters; (c) deals honestly with
  a 1955–65 analogue-era map. Pitch it as "a reproducible UK case study", not "first ever".

## 2. Map data

### BGS Great Britain Aeromagnetic Survey (primary)
- Endpoint: <https://ogcapi.bgs.ac.uk/collections/aeromag>, CSV via `items?f=csv`, GeoJSON via
  `items?f=json`. **No login.** Verified 2026-10-22 with curl: `numberMatched = 541508`.
- **Format: scattered points, not a grid.** These are flight-line/contour-line intercepts digitised from
  hand-contoured maps (10 nT contour interval). Columns: `survey_area, LINE_NUMB-SEG, order_by,
  map_sheet, longitude, latitude, grid_east, grid_north, ground_elev, mtc, aod, mag_bul,
  mag_igrf55, mag_igrf90, l1, l2, l3, feature_id`. `grid_east/north` = British National Grid (m),
  `mtc` = mean terrain clearance (m), `aod` = altitude above Ordnance Datum (m),
  `mag_igrf90` = anomaly (nT) w.r.t. a variant of IGRF-90.
- Survey: 1955–1965, Hunting Geology & Geophysics + Canadian Aeroservices; mostly 2 km line
  spacing, 10 km tie lines; "most of the country was flown at 1000 ft (305 m) mean terrain
  clearance" (collection description).
- Clearance heights, counted via the API (`mtc=` filter): 305 m: 430,892 pts; 152 m: 68,722
  (survey HG57, ~49.96–50.61°N, i.e. SW England); 457 m: 8,233 (mostly HG56). About 33,661 points
  have other values that I haven't enumerated yet. Phase 2 will check this, including the Highlands.
- Coverage bbox (API extent): −9 to 3°E, 49.5 to 61°N.
- Anomaly range in a 13.5k-point sample: −986 to +1524 nT.
- **Gridded version:** BGS shows a 1 km × 1 km grid ("variable tension technique, and smoothed")
  only as WMS images in GeoIndex ([data.gov.uk](https://www.data.gov.uk/dataset/05fd168e-a5ab-4e32-acf0-20a5821b122f/map-based-index-geoindex-magnetic-anomalies1)).
  **I found no downloadable grid of values**, so we grid it ourselves (e.g. Python `verde`, one
  offline step).
- **Licence: needs confirming.** The OGC API says "mostly Open Government Licence, please see
  individual collections", and the aeromag collection doesn't state one. The older data.gov.uk
  record has restrictive NERC wording. Fatiando's Harmonica/Ensaio redistributes this dataset as
  **OGL v3** with attribution "Contains NERC materials ©NERC [year]"
  ([notebook](https://www.fatiando.org/harmonica/v0.2.0/_downloads/d637f9133e6357827833b3a71e40928b/britain_magnetic.ipynb)).
  Action: email BGS enquiries to confirm, and don't commit raw data to the repo in the meantime
  (download script only).
- Quality caveats for the error budget:
  - 2 km line spacing at 305 m clearance means lines are ~6.5× the flight height apart
    **[derived]**, so short-wavelength anomalies between lines are undersampled. By contrast the
    MagNav challenge "Perth mini-survey" deliberately used line spacing < AGL height so the map is
    fully sampled (Gnadt et al. 2020).
  - Contour digitisation at a 10 nT interval: if the error is uniform within one interval, σ ≈
    10/√12 ≈ 2.9 nT **[derived, my assumption, not a published figure]**.
  - Mixed clearance heights (152/305/457 m) mean the map isn't at one altitude. Upward
    continuation to the flight altitude must handle that.

### Fallback: EMAG2v3 (Meyer, Saltus & Chulliat 2017)
- Global 2 arc-minute grid at **4 km above the geoid** (sea level over oceans), CSV/GeoTIFF, public
  download, no login, includes a per-cell error estimate
  ([NCEI](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.geophysical_models:EMAG2_V3/html),
  [doi:10.7289/V5H70CVX](https://doi.org/10.7289/V5H70CVX)). Licence: US public domain, but with
  the explicit constraint "Not to be used for navigation". That's fine for a simulation study, but
  the README must say so.
- Trade-off: 2′ ≈ 3.7 km N–S and ≈ 2.1 km E–W at 55°N **[derived]**. At 4 km altitude,
  short-wavelength anomalies are already attenuated, so it only represents high-altitude flight.
  Over GB it is largely built from the same BGS survey, so it is **not an independent map**.
- Useful role: a cross-check. Comparing BGS upward-continued to 4 km against EMAG2v3 gives an
  empirical map-error proxy.

## 3. Sensor parameters

| Class | Example | Noise floor | Drift / systematic | Source |
|---|---|---|---|---|
| Fluxgate (vector, 3-axis) | Bartington Mag-03 (±70/100 µT range) | <6 pT rms/√Hz @1 Hz (low-noise); 6–10 (standard); 10–20 (basic) | Offset error ±5 nT; offset tempco ±0.1 nT/°C; scaling tempco ±15 ppm/°C (±70 µT). In a ~50 µT field the scale term alone gives ≈0.75 nT/°C **[derived]**, so it dominates the offset tempco | [Bartington DS0013 datasheet](https://www.bartington.com/products/mag-03/) |
| Optically pumped, Cs scalar | Geometrics G-823A (flown in MagNav challenge) | <0.004 nT/√Hz rms (4 pT/√Hz) | Heading error ±0.15 nT over 360°; absolute accuracy <3 nT; no drift rate specified | [Geometrics G-823A](https://www.geometrics.com/product/g-823a/), [manual](https://www.geometrics.com/wp-content/uploads/2018/10/MAN822AO_revB5.pdf) |
| Optically pumped, Rb scalar (compact) | QuSpin QTFM Gen-2 | <5 pT/√Hz (typ. 3) | Heading error <3 nT uncompensated; no drift rate specified | [QuSpin QTFM Gen-2](https://quspin.com/qtfm-gen-2/) |
| NV diamond (lab, best case) | Fescenko et al., PRR 2, 023394 (2020) | ~0.9 pT/√Hz, 10–1000 Hz; within ×2 down to ~10 Hz, ~1/f below | Zero-field splitting shift χ ≈ −0.1 MHz/K → −24 nT/K in single-resonance mode (with 250× flux concentrators); dual-resonance cancels it to first order | [PMC7591154](https://pmc.ncbi.nlm.nih.gov/articles/PMC7591154/) |
| NV diamond (field-portable) | Dai et al., Adv. Quantum Technol. 2025 ([arXiv:2508.03237](https://arxiv.org/abs/2508.03237)) | ~2.14 nT/√Hz (optimal) | Not reported in abstract | as linked |
| NV diamond (vehicle-borne) | *On the road with a diamond magnetometer*, Diamond Relat. Mater. 2025 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0925963525000020)) | 0.3 ± 0.2 nT/√Hz, 10–150 Hz | Not found (paywalled; authors to confirm) | as linked |
| NV, intrinsic thermal shift | Acosta et al., PRL 104, 070801 (2010) | – | dD/dT ≈ −74 kHz/K → ≈2.6 µT/K equivalent if uncompensated (÷ 28 GHz/T) **[derived]** | Acosta 2010 (I cited this from memory, so check the PDF before using it in the README) |

Observations relevant to the hypothesis:
- Earth-field scalar OPMs are **pT/√Hz**, not fT/√Hz. fT/√Hz figures come from SERF-regime
  devices, which need near-zero ambient field and do not operate unshielded in Earth's ~50 µT field.
  This *strengthens* the hypothesis: even pT-level noise is far below 5–25 nT calibration error.
- Field-portable NV sensors are currently ~0.3–2 nT/√Hz, *noisier* than a fluxgate. NV's draw is
  vector + absolute + solid-state, not raw sensitivity.
- For OPMs, the published "drift" is really an orientation-dependent heading error, not a random
  walk. I'll model it that way rather than invent a drift rate.

### Not-yet-sourced numbers (do NOT use until sourced)
- INS drift rates for tactical/navigation-grade IMUs (needed in Phase 2; e.g. velocity random walk,
  bias instability from a datasheet such as a tactical-grade MEMS IMU).
- Long-term (hours) bias stability for fluxgate and NV in flight conditions.
