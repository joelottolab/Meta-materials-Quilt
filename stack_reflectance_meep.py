# MIT License
# Transfer Matrix Method (TMM) solver for multilayer reflectance & emissivity
# Usage:
#   python3 stack_reflectance_meep.py --config out_roof/resolved_config.json
#
# Notes:
# - Normal incidence
# - Materials specified by (n,k) per wavelength; defaults provided if not in config
# - For opaque stacks (no transmission), emissivity ε ≈ 1 - R by Kirchhoff's law
#
import json, argparse, os, numpy as np

def nk_const(n, k=0.0):
    return lambda lam_um: (n, k)

# Simple default material models (replace with real nk tables as needed)
DEFAULT_MATERIALS = {
    "Air": nk_const(1.0, 0.0),
    "PDMS": nk_const(1.41, 0.0),
    "SiO2_np": nk_const(1.45, 0.0),   # effective medium proxy for nanoparticle layer
    "TiO2_np": nk_const(2.40, 0.02),  # slight loss to mimic realistic TiO2
    "Adhesive": nk_const(1.40, 0.0),
    "Roof": nk_const(1.6, 0.02),      # lossy substrate (e.g., asphalt/painted metal)
}

def load_material_fn(name):
    if name in DEFAULT_MATERIALS:
        return DEFAULT_MATERIALS[name]
    raise ValueError(f"Unknown material: {name}")

def tmm_reflectance_normal(lam_um, layers):
    """
    Compute reflectance at normal incidence for stack:
    layers = [(n_complex, thickness_um), ...] with first and last semi-infinite.
    """
    # Convert to complex refractive indices at each wavelength
    lam = lam_um
    k0 = 2.0*np.pi/lam  # in 1/um

    # Build characteristic matrix
    # First medium (semi-infinite)
    n0 = layers[0][0]
    ns = layers[-1][0]
    eta0 = n0
    etas = ns

    M11 = np.ones_like(lam, dtype=complex)
    M12 = np.zeros_like(lam, dtype=complex)
    M21 = np.zeros_like(lam, dtype=complex)
    M22 = np.ones_like(lam, dtype=complex)

    for j in range(1, len(layers)-1):
        nj = layers[j][0]
        dj_um = layers[j][1]
        beta = k0*nj*dj_um  # phase thickness
        cj = np.cos(beta)
        sj = 1j*np.sin(beta)
        # Layer admittance (normal incidence): Y = n
        Yj = nj
        # Interface matrix for layer j
        M11, M12, M21, M22 = (M11*cj + (sj/Yj)*M12,
                               M11*sj*Yj + M12*cj,
                               M21*cj + (sj/Yj)*M22,
                               M21*sj*Yj + M22*cj)

    # Input admittance looking into stack
    Yin = (M11*etas + M12)/(M21*etas + M22)
    r = (Yin - eta0)/(Yin + eta0)
    R = np.real(r*np.conj(r))
    R = np.clip(R, 0.0, 1.0)
    return R

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    with open(args.config) as f:
        cfg = json.load(f)

    out_dir = os.path.dirname(args.config)
    os.makedirs(out_dir, exist_ok=True)

    lam = np.linspace(cfg["simulation"]["lambda_min_um"],
                      cfg["simulation"]["lambda_max_um"],
                      cfg["simulation"]["points"])

    # Build layer list: [(n_complex(lam), thickness_um), ...]
    layers_desc = cfg["layers"]
    layers = []
    for i, layer in enumerate(layers_desc):
        name = layer["material"]
        t_um = float(layer["thickness_um"])
        nk_fn = load_material_fn(name)
        n_list = []
        for L in lam:
            n, k = nk_fn(L)
            n_list.append(n + 1j*k)
        n_arr = np.array(n_list)
        layers.append([n_arr, t_um])

    R = tmm_reflectance_normal(lam, layers)
    eps = 1.0 - R  # opaque assumption

    np.savetxt(os.path.join(out_dir, "spectrum.csv"),
               np.column_stack([lam, R, eps]),
               delimiter=",", header="lambda_um,Reflectance,Emissivity", comments="")

    mask = (lam>=8.0)&(lam<=13.0)
    eps_win = float(np.mean(eps[mask]))
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump({"emissivity_8_13um": eps_win}, f, indent=2)

    print("Done. Saved:", os.path.join(out_dir, "spectrum.csv"))
    print("ε(8–13µm) =", eps_win)

if __name__ == "__main__":
    main()
