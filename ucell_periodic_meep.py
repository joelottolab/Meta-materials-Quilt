# MIT License
# RF metasurface periodic unit cell with grounded dielectric + square patch ("mushroom" type)
# This sets a concrete geometry suitable for Meep when available.
#
# Usage:
#   python3 ucell_periodic_meep.py --config out_rf/resolved_config.json
#
import json, argparse, os, numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    with open(args.config) as f:
        cfg = json.load(f)
    out_dir = os.path.dirname(args.config)
    os.makedirs(out_dir, exist_ok=True)

    # Geometry parameters from config
    pitch = cfg["unit_cell"]["pitch_mm"]*1e-3     # m
    t_sub = cfg["unit_cell"]["substrate_thickness_mm"]*1e-3
    t_cu  = cfg["unit_cell"]["conductor_thickness_um"]*1e-6
    fill  = cfg["unit_cell"]["copper_fill_factor"]

    # Define square patch size on top of substrate
    patch_size = np.sqrt(fill) * pitch  # preserves area fraction
    gap = 0.5e-3  # 0.5 mm gap to neighboring patches (placeholder)

    # We save a JSON "geometry description" that a Meep-enabled environment can consume:
    geom = {
        "domain": {"sx": pitch, "sy": pitch, "sz": t_sub + 4*t_cu},
        "pml": {"thickness": cfg["simulation"]["pml_thickness_cells"]},
        "materials": {
            "substrate_eps": 3.66,
            "metal": "pec"
        },
        "objects": [
            {"type":"box","name":"substrate",
             "size":[pitch, pitch, t_sub],
             "center":[0,0,-0.5*(t_sub)],
             "material":{"eps":3.66}},
            {"type":"box","name":"ground",
             "size":[pitch, pitch, t_cu],
             "center":[0,0,-0.5*t_sub - 0.5*t_cu],
             "material":{"metal":"pec"}},
            {"type":"box","name":"patch",
             "size":[patch_size, patch_size, t_cu],
             "center":[0,0, 0.5*t_cu],
             "material":{"metal":"pec"}}
        ],
        "ports": {
            "source_z": 0.02,    # source plane +z
            "probe_z": -0.02     # probe plane -z
        },
        "bloch": {
            "kx": cfg["simulation"]["bloch_kx"],
            "ky": cfg["simulation"]["bloch_ky"]
        }
    }

    with open(os.path.join(out_dir, "meep_geometry.json"), "w") as f:
        json.dump(geom, f, indent=2)

    # Since Meep may not be present in this environment, we still emit a stub S-parameter file
    # so downstream plotting pipelines work. Replace with real run when Meep is available.
    freqs = np.linspace(cfg["simulation"]["freq_min_GHz"], cfg["simulation"]["freq_max_GHz"], 401)
    # Very rough toy response: simple resonance tied to patch length ~ c/(2nL)
    c = 3e8
    n_eff = np.sqrt(3.66)
    f0 = c/(2*n_eff*patch_size) / 1e9  # GHz
    S11 = 1.0/(1.0 + ((freqs-f0)/(0.15*f0+1e-9))**2)  # Lorenztian-ish
    S11 = np.clip(S11, 0, 1)
    S21 = 1 - S11
    np.savetxt(os.path.join(out_dir, "sparams.csv"),
               np.column_stack([freqs, S11, S21]),
               delimiter=",", header="GHz,S11_mag,S21_mag", comments="")
    with open(os.path.join(out_dir, "notes.txt"), "w") as f:
        f.write("Geometry defined for grounded substrate + square patch. Use meep to run real FDTD.\n")
        f.write(f"Estimated resonant frequency (GHz): {f0:.2f}\n")

    print("Geometry description saved to meep_geometry.json")
    print("Stub S-parameters saved to sparams.csv")
    print(f"Estimated patch resonance ~ {f0:.2f} GHz (toy estimate).")

if __name__ == "__main__":
    main()
