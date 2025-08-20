# Metamaterials Simulation Setups (High-Iteration Editions)

This package contains **ready-to-run** FDTD/Meep-style simulation scaffolds for two projects:

- **Project A (RF):** Programmable Huygens metasurface unit cell (2.4–6 GHz), periodic model with Bloch boundaries.
- **Project B (Roof):** 1D multilayer radiative-cooling film (visible–mid-IR reflectance/emissivity).

These scripts are configured for **10× heavier iterations** than a quick demo:
- finer **resolution**
- longer **simulation time**
- tighter **convergence** thresholds

> If you hit performance limits, set `FAST_MODE: true` in the YAML to auto-scale down.

## Requirements

- Python 3.9+
- meep (Python bindings), numpy, matplotlib, pyyaml
- Optional: mpirun for parallel meep

## Structure

```
metamaterials_sim_setups/
├── rf/
│   ├── ucell_periodic_meep.py
│   ├── param_sweep_rf.yaml
│   └── batch_runner_rf.py
├── roof/
│   ├── stack_reflectance_meep.py
│   ├── param_sweep_roof.yaml
│   └── batch_runner_roof.py
└── shared/
    ├── materials_library.py
    └── analysis_templates.ipynb
```

## Quick Start (RF)

```bash
cd rf
python3 batch_runner_rf.py --config param_sweep_rf.yaml --out out_rf
python3 ucell_periodic_meep.py --config out_rf/resolved_config.json
```

## Quick Start (Roof)

```bash
cd roof
python3 batch_runner_roof.py --config param_sweep_roof.yaml --out out_roof
python3 stack_reflectance_meep.py --config out_roof/resolved_config.json
```

Results will include spectra (S-parameters, reflectance/emissivity) and CSV logs.

## License

CERN-OHL-P v2 for hardware-like content; Code under MIT. See headers.
